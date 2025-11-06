#!/usr/bin/env python3
"""
Compare particle trajectories across multiple simulation runs.

This script computes trajectory deviation metrics to quantify stochasticity
in DFC slump test results when running on GPUs. It uses run 1 as the baseline
and compares all other runs against it.

The main metrics computed are:
1. Trajectory Deviation: Average point-wise distance between matched trajectories
2. Relative Trajectory Deviation: Deviation normalized by average path length
3. Path Length Difference: Difference in total distance traveled

Usage:
    python compare_trajectory_runs.py --baseline OUT_slump_V1/particle_trajectories.csv \\
                                      --runs OUT_slump_V2/particle_trajectories.csv \\
                                            OUT_slump_V3/particle_trajectories.csv
    
    # Or specify runs using a pattern
    python compare_trajectory_runs.py --baseline run1/particle_trajectories.csv \\
                                      --runs run*/particle_trajectories.csv
"""

import argparse
import csv
import math
import sys
import glob
import os
from typing import Dict, List, Tuple

# Check for required dependencies
try:
    import numpy as np
except ImportError:
    print("Error: numpy is required. Install it with: pip install numpy")
    sys.exit(1)

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:
    print("Error: scipy is required. Install it with: pip install scipy")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is required. Install it with: pip install matplotlib")
    sys.exit(1)

def load_trajectories(filename: str) -> Dict[int, List[Tuple[float, float, float, float]]]:
    """
    Load trajectory data from CSV file.
    
    Args:
        filename: Path to the trajectory CSV file
        
    Returns:
        Dictionary mapping particle_id to list of (time, x, y, z) tuples
    """
    trajectories = {}
    
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            particle_id = int(row['particle_id'])
            time = float(row['time'])
            x = float(row['x'])
            y = float(row['y'])
            z = float(row['z'])
            
            if particle_id not in trajectories:
                trajectories[particle_id] = []
            trajectories[particle_id].append((time, x, y, z))
    
    # Sort each trajectory by time
    for particle_id in trajectories:
        trajectories[particle_id].sort(key=lambda point: point[0])
    
    return trajectories

def calculate_path_length(trajectory: List[Tuple[float, float, float, float]]) -> float:
    """
    Calculate total path length traveled by a particle.
    
    Args:
        trajectory: List of (time, x, y, z) tuples
        
    Returns:
        Total path length
    """
    if len(trajectory) < 2:
        return 0.0
    
    total_length = 0.0
    for i in range(1, len(trajectory)):
        prev_pos = trajectory[i-1][1:4]  # (x, y, z)
        curr_pos = trajectory[i][1:4]
        
        dx = curr_pos[0] - prev_pos[0]
        dy = curr_pos[1] - prev_pos[1]
        dz = curr_pos[2] - prev_pos[2]
        
        total_length += math.sqrt(dx**2 + dy**2 + dz**2)
    
    return total_length

# Constant for time epsilon to avoid division by zero in interpolation
TIME_EPSILON = 1e-10

def interpolate_trajectory_at_time(trajectory: List[Tuple[float, float, float, float]], 
                                  target_time: float) -> Tuple[float, float, float]:
    """
    Interpolate particle position at a specific time.
    
    Args:
        trajectory: List of (time, x, y, z) tuples
        target_time: Time at which to interpolate position
        
    Returns:
        Tuple of (x, y, z) at target_time
    """
    if not trajectory:
        raise ValueError("Empty trajectory")
    
    # If target time is before first point, return first position
    if target_time <= trajectory[0][0]:
        return trajectory[0][1:4]
    
    # If target time is after last point, return last position
    if target_time >= trajectory[-1][0]:
        return trajectory[-1][1:4]
    
    # Find bracketing time points
    for i in range(len(trajectory) - 1):
        t1, x1, y1, z1 = trajectory[i]
        t2, x2, y2, z2 = trajectory[i + 1]
        
        if t1 <= target_time <= t2:
            # Linear interpolation
            if t2 - t1 < TIME_EPSILON:  # Avoid division by zero
                return (x1, y1, z1)
            
            alpha = (target_time - t1) / (t2 - t1)
            x = x1 + alpha * (x2 - x1)
            y = y1 + alpha * (y2 - y1)
            z = z1 + alpha * (z2 - z1)
            return (x, y, z)
    
    # Should not reach here
    return trajectory[-1][1:4]

def compute_trajectory_deviation(baseline_traj: List[Tuple[float, float, float, float]], 
                                 comparison_traj: List[Tuple[float, float, float, float]]) -> Tuple[float, float, int]:
    """
    Compute deviation between two trajectories.
    
    This function interpolates both trajectories at common time points and
    computes the average Euclidean distance between corresponding positions.
    
    Args:
        baseline_traj: Baseline trajectory (time, x, y, z)
        comparison_traj: Comparison trajectory (time, x, y, z)
        
    Returns:
        Tuple of (average_deviation, max_deviation, num_points)
    """
    if not baseline_traj or not comparison_traj:
        return 0.0, 0.0, 0
    
    # Get common time range
    t_start = max(baseline_traj[0][0], comparison_traj[0][0])
    t_end = min(baseline_traj[-1][0], comparison_traj[-1][0])
    
    if t_start >= t_end:
        return 0.0, 0.0, 0
    
    # Collect all unique time points from both trajectories
    time_points = set()
    for t, _, _, _ in baseline_traj:
        if t_start <= t <= t_end:
            time_points.add(t)
    for t, _, _, _ in comparison_traj:
        if t_start <= t <= t_end:
            time_points.add(t)
    
    time_points = sorted(time_points)
    
    if not time_points:
        return 0.0, 0.0, 0
    
    # Compute deviations at each time point
    deviations = []
    for t in time_points:
        pos1 = interpolate_trajectory_at_time(baseline_traj, t)
        pos2 = interpolate_trajectory_at_time(comparison_traj, t)
        
        dx = pos2[0] - pos1[0]
        dy = pos2[1] - pos1[1]
        dz = pos2[2] - pos1[2]
        
        deviation = math.sqrt(dx**2 + dy**2 + dz**2)
        deviations.append(deviation)
    
    avg_deviation = sum(deviations) / len(deviations)
    max_deviation = max(deviations)
    
    return avg_deviation, max_deviation, len(time_points)

def compute_cumulative_path_length_at_time(trajectory: List[Tuple[float, float, float, float]], 
                                          target_time: float) -> float:
    """
    Compute cumulative path length traveled up to a specific time.
    
    Args:
        trajectory: List of (time, x, y, z) tuples
        target_time: Time up to which to compute path length
        
    Returns:
        Cumulative path length traveled up to target_time
    """
    if not trajectory or len(trajectory) < 2:
        return 0.0
    
    cumulative_length = 0.0
    
    # Find all trajectory segments up to target_time
    for i in range(1, len(trajectory)):
        t_prev = trajectory[i-1][0]
        t_curr = trajectory[i][0]
        
        if t_curr <= target_time:
            # Full segment is before target_time
            prev_pos = trajectory[i-1][1:4]
            curr_pos = trajectory[i][1:4]
            
            dx = curr_pos[0] - prev_pos[0]
            dy = curr_pos[1] - prev_pos[1]
            dz = curr_pos[2] - prev_pos[2]
            
            cumulative_length += math.sqrt(dx**2 + dy**2 + dz**2)
        elif t_prev < target_time < t_curr:
            # Partial segment - interpolate position at target_time
            prev_pos = trajectory[i-1][1:4]
            curr_pos = trajectory[i][1:4]
            
            # Interpolate position at target_time
            alpha = (target_time - t_prev) / (t_curr - t_prev) if (t_curr - t_prev) > TIME_EPSILON else 0.0
            interp_x = prev_pos[0] + alpha * (curr_pos[0] - prev_pos[0])
            interp_y = prev_pos[1] + alpha * (curr_pos[1] - prev_pos[1])
            interp_z = prev_pos[2] + alpha * (curr_pos[2] - prev_pos[2])
            
            dx = interp_x - prev_pos[0]
            dy = interp_y - prev_pos[1]
            dz = interp_z - prev_pos[2]
            
            cumulative_length += math.sqrt(dx**2 + dy**2 + dz**2)
            break
        else:
            # Segment is after target_time, stop
            break
    
    return cumulative_length

def compute_trajectory_deviation_evolution(baseline_traj: List[Tuple[float, float, float, float]], 
                                          comparison_traj: List[Tuple[float, float, float, float]]) -> Tuple[List[float], List[float]]:
    """
    Compute relative deviation evolution over time between two trajectories.
    
    This function returns the relative deviation (deviation / cumulative path length)
    at each time point, allowing visualization of how relative deviation changes
    over the simulation. This normalizes the deviation by how much the particle
    has traveled, making it easier to compare particles with different path lengths.
    
    Args:
        baseline_traj: Baseline trajectory (time, x, y, z)
        comparison_traj: Comparison trajectory (time, x, y, z)
        
    Returns:
        Tuple of (time_points, relative_deviations) where both are lists of same length
    """
    if not baseline_traj or not comparison_traj:
        return [], []
    
    # Get common time range
    t_start = max(baseline_traj[0][0], comparison_traj[0][0])
    t_end = min(baseline_traj[-1][0], comparison_traj[-1][0])
    
    if t_start >= t_end:
        return [], []
    
    # Collect all unique time points from both trajectories
    time_points = set()
    for t, _, _, _ in baseline_traj:
        if t_start <= t <= t_end:
            time_points.add(t)
    for t, _, _, _ in comparison_traj:
        if t_start <= t <= t_end:
            time_points.add(t)
    
    time_points = sorted(time_points)
    
    if not time_points:
        return [], []
    
    # Compute relative deviations at each time point
    relative_deviations = []
    for t in time_points:
        pos1 = interpolate_trajectory_at_time(baseline_traj, t)
        pos2 = interpolate_trajectory_at_time(comparison_traj, t)
        
        # Compute absolute deviation
        dx = pos2[0] - pos1[0]
        dy = pos2[1] - pos1[1]
        dz = pos2[2] - pos1[2]
        
        deviation = math.sqrt(dx**2 + dy**2 + dz**2)
        
        # Compute cumulative path length in baseline trajectory up to time t
        cumulative_length = compute_cumulative_path_length_at_time(baseline_traj, t)
        
        # Compute relative deviation (avoid division by zero)
        if cumulative_length > TIME_EPSILON:
            relative_deviation = deviation / cumulative_length
        else:
            # At the very start, use a small value or skip
            relative_deviation = 0.0
        
        relative_deviations.append(relative_deviation)
    
    return time_points, relative_deviations

def plot_deviation_evolution(baseline_trajs: Dict[int, List], 
                             comparison_trajs: Dict[int, List],
                             matching: Dict[int, int],
                             run_name: str,
                             output_file: str) -> None:
    """
    Plot relative deviation evolution over time for all matched particles.
    
    Creates a plot with one curve per particle showing how relative deviation
    (deviation / cumulative path length) from baseline changes throughout the simulation.
    This allows comparison of particles with different trajectory lengths by normalizing
    the deviation by the distance traveled.
    
    Args:
        baseline_trajs: Baseline trajectories {particle_id: [(time, x, y, z), ...]}
        comparison_trajs: Comparison trajectories
        matching: Dictionary mapping baseline particle IDs to comparison particle IDs
        run_name: Name of the comparison run for plot title
        output_file: Path to save the plot
    """
    plt.figure(figsize=(12, 8))
    
    # Use a colormap for different particles
    num_particles = len(matching)
    colors = plt.cm.rainbow(np.linspace(0, 1, num_particles))
    
    max_relative_deviation = 0.0
    max_time = 0.0
    
    # Plot deviation evolution for each matched particle
    for idx, (baseline_pid, comp_pid) in enumerate(sorted(matching.items())):
        baseline_traj = baseline_trajs[baseline_pid]
        comp_traj = comparison_trajs[comp_pid]
        
        time_points, relative_deviations = compute_trajectory_deviation_evolution(baseline_traj, comp_traj)
        
        if time_points and relative_deviations:
            # Plot with alpha for transparency when many curves
            alpha = 0.3 if num_particles > 20 else 0.6
            plt.plot(time_points, relative_deviations, color=colors[idx], alpha=alpha, linewidth=1.0)
            
            # Track maximum values for axis scaling
            if relative_deviations:
                max_relative_deviation = max(max_relative_deviation, max(relative_deviations))
            if time_points:
                max_time = max(max_time, max(time_points))
    
    plt.xlabel('Simulation Time (s)', fontsize=12)
    plt.ylabel('Relative Deviation (deviation / cumulative path length)', fontsize=12)
    plt.title(f'Particle Trajectory Relative Deviation Evolution - {run_name}\n({num_particles} particles)', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Add some padding to the y-axis
    if max_relative_deviation > 0:
        plt.ylim(0, max_relative_deviation * 1.1)
    
    # Add text annotation
    plt.text(0.02, 0.98, f'Each curve represents one particle\nShowing relative deviation (normalized by cumulative path length)', 
             transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"    Saved deviation evolution plot: {output_file}")

def match_particles_by_initial_position(baseline_trajs: Dict[int, List], 
                                       comparison_trajs: Dict[int, List]) -> Dict[int, int]:
    """
    Match particles between two runs based on their initial positions.
    
    Uses the Hungarian algorithm (linear_sum_assignment) to find the optimal
    matching that minimizes total initial position distance.
    
    Args:
        baseline_trajs: Baseline trajectories {particle_id: [(time, x, y, z), ...]}
        comparison_trajs: Comparison trajectories
        
    Returns:
        Dictionary mapping baseline particle IDs to comparison particle IDs
    """
    baseline_ids = sorted(baseline_trajs.keys())
    comparison_ids = sorted(comparison_trajs.keys())
    
    if not baseline_ids or not comparison_ids:
        return {}
    
    # Build cost matrix based on initial positions
    n = len(baseline_ids)
    m = len(comparison_ids)
    cost_matrix = np.zeros((n, m))
    
    for i, bid in enumerate(baseline_ids):
        baseline_init_pos = baseline_trajs[bid][0][1:4]  # (x, y, z) at t=0
        
        for j, cid in enumerate(comparison_ids):
            comparison_init_pos = comparison_trajs[cid][0][1:4]
            
            # Euclidean distance between initial positions
            dx = baseline_init_pos[0] - comparison_init_pos[0]
            dy = baseline_init_pos[1] - comparison_init_pos[1]
            dz = baseline_init_pos[2] - comparison_init_pos[2]
            cost_matrix[i, j] = math.sqrt(dx**2 + dy**2 + dz**2)
    
    # Find optimal assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    # Create mapping
    matching = {}
    for i, j in zip(row_ind, col_ind):
        matching[baseline_ids[i]] = comparison_ids[j]
    
    return matching

def compare_runs(baseline_file: str, comparison_files: List[str]) -> None:
    """
    Compare multiple runs against a baseline run.
    
    Args:
        baseline_file: Path to baseline trajectory CSV file
        comparison_files: List of paths to comparison trajectory CSV files
    """
    print("=" * 80)
    print("TRAJECTORY COMPARISON ACROSS RUNS")
    print("=" * 80)
    print()
    
    # Load baseline trajectories
    print(f"Loading baseline from: {baseline_file}")
    try:
        baseline_trajs = load_trajectories(baseline_file)
        print(f"  Loaded {len(baseline_trajs)} particles")
    except Exception as e:
        print(f"Error loading baseline: {e}")
        sys.exit(1)
    
    # Compute baseline statistics
    baseline_path_lengths = {pid: calculate_path_length(traj) 
                            for pid, traj in baseline_trajs.items()}
    avg_baseline_path_length = sum(baseline_path_lengths.values()) / len(baseline_path_lengths)
    print(f"  Average path length: {avg_baseline_path_length:.2f} units")
    print()
    
    # Compare each run
    all_results = []
    
    for comp_idx, comp_file in enumerate(comparison_files, start=2):
        print(f"Comparing Run {comp_idx}: {comp_file}")
        
        try:
            comp_trajs = load_trajectories(comp_file)
            print(f"  Loaded {len(comp_trajs)} particles")
        except Exception as e:
            print(f"  Error loading file: {e}")
            print()
            continue
        
        # Match particles by initial position
        matching = match_particles_by_initial_position(baseline_trajs, comp_trajs)
        print(f"  Matched {len(matching)} particle pairs")
        
        # Compute deviations for each matched particle
        deviations = []
        max_deviations = []
        relative_deviations = []
        path_length_diffs = []
        
        for baseline_pid, comp_pid in matching.items():
            baseline_traj = baseline_trajs[baseline_pid]
            comp_traj = comp_trajs[comp_pid]
            
            # Compute trajectory deviation
            avg_dev, max_dev, num_points = compute_trajectory_deviation(baseline_traj, comp_traj)
            
            if num_points > 0:
                deviations.append(avg_dev)
                max_deviations.append(max_dev)
                
                # Compute relative deviation (normalized by path length)
                baseline_path_len = baseline_path_lengths[baseline_pid]
                if baseline_path_len > 0:
                    relative_dev = avg_dev / baseline_path_len
                    relative_deviations.append(relative_dev)
                
                # Compute path length difference
                comp_path_len = calculate_path_length(comp_traj)
                path_length_diff = abs(comp_path_len - baseline_path_len)
                path_length_diffs.append(path_length_diff)
        
        # Compute summary statistics
        if deviations:
            avg_deviation = sum(deviations) / len(deviations)
            max_deviation_overall = max(max_deviations)
            avg_relative_deviation = sum(relative_deviations) / len(relative_deviations) if relative_deviations else 0.0
            avg_path_length_diff = sum(path_length_diffs) / len(path_length_diffs)
            
            print(f"  Average trajectory deviation: {avg_deviation:.4f} units")
            print(f"  Maximum trajectory deviation: {max_deviation_overall:.4f} units")
            print(f"  Average relative deviation: {avg_relative_deviation:.6f} ({avg_relative_deviation*100:.4f}%)")
            print(f"  Average path length difference: {avg_path_length_diff:.4f} units")
            
            # Generate deviation evolution plot for this run
            # Extract filename without extension for plot naming
            comp_basename = os.path.splitext(os.path.basename(comp_file))[0]
            plot_filename = f"deviation_evolution_run{comp_idx}_{comp_basename}.png"
            print(f"  Generating deviation evolution plot...")
            try:
                plot_deviation_evolution(baseline_trajs, comp_trajs, matching, 
                                       f"Run {comp_idx}", plot_filename)
            except Exception as e:
                print(f"    Warning: Could not generate plot: {e}")
            
            all_results.append({
                'run_idx': comp_idx,
                'file': comp_file,
                'num_matched': len(matching),
                'avg_deviation': avg_deviation,
                'max_deviation': max_deviation_overall,
                'avg_relative_deviation': avg_relative_deviation,
                'avg_path_length_diff': avg_path_length_diff
            })
        else:
            print("  No deviations computed (no matched trajectories)")
        
        print()
    
    # Print summary table
    if all_results:
        print("=" * 80)
        print("SUMMARY TABLE")
        print("=" * 80)
        print()
        print(f"{'Run':<6} {'Matched':<10} {'Avg Dev':<12} {'Max Dev':<12} {'Rel Dev %':<12} {'Path Diff':<12}")
        print("-" * 80)
        
        for result in all_results:
            print(f"{result['run_idx']:<6} {result['num_matched']:<10} "
                  f"{result['avg_deviation']:<12.4f} {result['max_deviation']:<12.4f} "
                  f"{result['avg_relative_deviation']*100:<12.4f} {result['avg_path_length_diff']:<12.4f}")
        
        print("-" * 80)
        
        # Compute statistics across all runs
        avg_deviations = [r['avg_deviation'] for r in all_results]
        rel_deviations = [r['avg_relative_deviation'] for r in all_results]
        
        print(f"{'Mean':<6} {'':<10} {sum(avg_deviations)/len(avg_deviations):<12.4f} "
              f"{'':<12} {sum(rel_deviations)/len(rel_deviations)*100:<12.4f} {'':<12}")
        print(f"{'StdDev':<6} {'':<10} {np.std(avg_deviations):<12.4f} "
              f"{'':<12} {np.std(rel_deviations)*100:<12.4f} {'':<12}")
        
        print()
        print("=" * 80)
        print("INTERPRETATION")
        print("=" * 80)
        print()
        print("Metrics:")
        print("  - Avg Dev: Average point-wise distance between matched trajectories")
        print("  - Max Dev: Maximum point-wise distance observed")
        print("  - Rel Dev %: Deviation as percentage of path length (measures relative error)")
        print("  - Path Diff: Average difference in total distance traveled")
        print()
        print("The Relative Deviation percentage shows the trajectory stochasticity")
        print("as a fraction of the total motion. Lower values indicate more reproducible results.")
        print()

def main():
    parser = argparse.ArgumentParser(
        description='Compare particle trajectories across multiple simulation runs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--baseline', '-b', 
                       required=True,
                       help='Baseline trajectory CSV file (e.g., run 1)')
    parser.add_argument('--runs', '-r',
                       nargs='+',
                       required=True,
                       help='Comparison trajectory CSV files (e.g., runs 2-10)')
    
    args = parser.parse_args()
    
    # Expand glob patterns in comparison files
    comparison_files = []
    for pattern in args.runs:
        matched_files = glob.glob(pattern)
        if matched_files:
            comparison_files.extend(matched_files)
        else:
            comparison_files.append(pattern)  # Keep as-is if no glob match
    
    # Remove baseline from comparison files if it was included
    comparison_files = [f for f in comparison_files if f != args.baseline]
    
    if not comparison_files:
        print("Error: No comparison files specified")
        sys.exit(1)
    
    print(f"Baseline: {args.baseline}")
    print(f"Comparing against {len(comparison_files)} runs")
    print()
    
    compare_runs(args.baseline, comparison_files)

if __name__ == '__main__':
    main()
