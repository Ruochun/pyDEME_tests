#!/usr/bin/env python3
"""
Visualization script for particle trajectories from pyDEME_DFCSlump.py or pyDEME_DFCLBox.py

This script reads the particle_trajectories.csv file and creates visualizations
of the particle trajectories during the simulation.

Usage:
    python visualize_trajectories.py [--input INPUT_FILE] [--output OUTPUT_FILE] [--particles PARTICLE_IDS]

Examples:
    # Visualize all tracked particles
    python visualize_trajectories.py --input OUT_slump_V1/particle_trajectories.csv
    
    # Visualize specific particles (e.g., particles 0, 5, 10)
    python visualize_trajectories.py --input OUT_slump_V1/particle_trajectories.csv --particles 0,5,10
    
    # Save to a specific output file
    python visualize_trajectories.py --input OUT_slump_V1/particle_trajectories.csv --output my_trajectories.png
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Configuration constants
MAX_LEGEND_PARTICLES = 20  # Maximum particles to show in legend before simplifying
DEFAULT_TIME_PLOT_PARTICLES = 10  # Default number of particles for time evolution plots

def load_trajectory_data(filename):
    """Load trajectory data from CSV file."""
    df = pd.read_csv(filename)
    return df

def plot_3d_trajectories(df, particle_ids=None, output_file=None):
    """
    Create a 3D visualization of particle trajectories.
    
    Args:
        df: DataFrame with columns [time, particle_id, x, y, z]
        particle_ids: List of particle IDs to plot. If None, plot all particles.
        output_file: If provided, save the plot to this file instead of displaying.
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get unique particle IDs
    all_particle_ids = df['particle_id'].unique()
    
    # Filter to requested particles
    if particle_ids is not None:
        plot_particle_ids = [pid for pid in particle_ids if pid in all_particle_ids]
        if not plot_particle_ids:
            print(f"Warning: None of the requested particles {particle_ids} found in data.")
            print(f"Available particles: {sorted(all_particle_ids)}")
            return
    else:
        plot_particle_ids = all_particle_ids
    
    # Use a colormap for different particles
    colors = plt.cm.rainbow(np.linspace(0, 1, len(plot_particle_ids)))
    
    # Plot each particle's trajectory
    for idx, particle_id in enumerate(sorted(plot_particle_ids)):
        particle_data = df[df['particle_id'] == particle_id].sort_values('time')
        
        # Plot the trajectory line
        ax.plot(particle_data['x'], particle_data['y'], particle_data['z'], 
                color=colors[idx], alpha=0.6, linewidth=1.5, 
                label=f'Particle {int(particle_id)}')
        
        # Mark start position with a larger marker
        ax.scatter(particle_data['x'].iloc[0], 
                  particle_data['y'].iloc[0], 
                  particle_data['z'].iloc[0], 
                  color=colors[idx], marker='o', s=50, alpha=0.8)
        
        # Mark end position with a different marker
        ax.scatter(particle_data['x'].iloc[-1], 
                  particle_data['y'].iloc[-1], 
                  particle_data['z'].iloc[-1], 
                  color=colors[idx], marker='s', s=50, alpha=0.8)
    
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_title('Particle Trajectories During Slump Test\n(circles=start, squares=end)')
    
    # Add legend (limit to avoid clutter)
    if len(plot_particle_ids) <= MAX_LEGEND_PARTICLES:
        ax.legend(loc='upper left', fontsize=8)
    else:
        ax.text2D(0.02, 0.98, f'Showing {len(plot_particle_ids)} particles', 
                  transform=ax.transAxes, fontsize=10, verticalalignment='top')
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved trajectory plot to {output_file}")
    else:
        plt.show()

def plot_2d_projections(df, particle_ids=None, output_file=None):
    """
    Create 2D projection plots (XY, XZ, YZ views).
    
    Args:
        df: DataFrame with columns [time, particle_id, x, y, z]
        particle_ids: List of particle IDs to plot. If None, plot all particles.
        output_file: If provided, save the plot to this file instead of displaying.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Get unique particle IDs
    all_particle_ids = df['particle_id'].unique()
    
    # Filter to requested particles
    if particle_ids is not None:
        plot_particle_ids = [pid for pid in particle_ids if pid in all_particle_ids]
    else:
        plot_particle_ids = all_particle_ids
    
    # Use a colormap for different particles
    colors = plt.cm.rainbow(np.linspace(0, 1, len(plot_particle_ids)))
    
    # Plot each particle's trajectory in different projections
    for idx, particle_id in enumerate(sorted(plot_particle_ids)):
        particle_data = df[df['particle_id'] == particle_id].sort_values('time')
        
        # XY projection
        axes[0].plot(particle_data['x'], particle_data['y'], 
                    color=colors[idx], alpha=0.6, linewidth=1.5)
        axes[0].scatter(particle_data['x'].iloc[0], particle_data['y'].iloc[0], 
                       color=colors[idx], marker='o', s=30)
        
        # XZ projection
        axes[1].plot(particle_data['x'], particle_data['z'], 
                    color=colors[idx], alpha=0.6, linewidth=1.5)
        axes[1].scatter(particle_data['x'].iloc[0], particle_data['z'].iloc[0], 
                       color=colors[idx], marker='o', s=30)
        
        # YZ projection
        axes[2].plot(particle_data['y'], particle_data['z'], 
                    color=colors[idx], alpha=0.6, linewidth=1.5,
                    label=f'P{int(particle_id)}' if len(plot_particle_ids) <= 10 else None)
        axes[2].scatter(particle_data['y'].iloc[0], particle_data['z'].iloc[0], 
                       color=colors[idx], marker='o', s=30)
    
    axes[0].set_xlabel('X Position')
    axes[0].set_ylabel('Y Position')
    axes[0].set_title('XY Projection (Top View)')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('X Position')
    axes[1].set_ylabel('Z Position')
    axes[1].set_title('XZ Projection (Front View)')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].set_xlabel('Y Position')
    axes[2].set_ylabel('Z Position')
    axes[2].set_title('YZ Projection (Side View)')
    axes[2].grid(True, alpha=0.3)
    if len(plot_particle_ids) <= DEFAULT_TIME_PLOT_PARTICLES:
        axes[2].legend(fontsize=8)
    
    plt.suptitle(f'Particle Trajectory Projections ({len(plot_particle_ids)} particles)')
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved projection plots to {output_file}")
    else:
        plt.show()

def plot_time_evolution(df, particle_ids=None, output_file=None):
    """
    Plot the evolution of positions over time.
    
    Args:
        df: DataFrame with columns [time, particle_id, x, y, z]
        particle_ids: List of particle IDs to plot. If None, plot a subset.
        output_file: If provided, save the plot to this file instead of displaying.
    """
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    # Get unique particle IDs
    all_particle_ids = df['particle_id'].unique()
    
    # Filter to requested particles, or use first 10 if none specified
    if particle_ids is not None:
        plot_particle_ids = [pid for pid in particle_ids if pid in all_particle_ids]
    else:
        plot_particle_ids = sorted(all_particle_ids)[:DEFAULT_TIME_PLOT_PARTICLES]  # Limit for clarity
    
    # Use a colormap for different particles
    colors = plt.cm.rainbow(np.linspace(0, 1, len(plot_particle_ids)))
    
    # Plot each particle's position components over time
    for idx, particle_id in enumerate(sorted(plot_particle_ids)):
        particle_data = df[df['particle_id'] == particle_id].sort_values('time')
        
        axes[0].plot(particle_data['time'], particle_data['x'], 
                    color=colors[idx], alpha=0.7, linewidth=1.5, 
                    label=f'Particle {int(particle_id)}')
        axes[1].plot(particle_data['time'], particle_data['y'], 
                    color=colors[idx], alpha=0.7, linewidth=1.5)
        axes[2].plot(particle_data['time'], particle_data['z'], 
                    color=colors[idx], alpha=0.7, linewidth=1.5)
    
    axes[0].set_ylabel('X Position')
    axes[0].set_title('X Position vs Time')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8, ncol=2)
    
    axes[1].set_ylabel('Y Position')
    axes[1].set_title('Y Position vs Time')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].set_ylabel('Z Position')
    axes[2].set_xlabel('Time')
    axes[2].set_title('Z Position vs Time')
    axes[2].grid(True, alpha=0.3)
    
    plt.suptitle(f'Position Evolution Over Time ({len(plot_particle_ids)} particles)')
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved time evolution plot to {output_file}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(
        description='Visualize particle trajectories from pyDEME slump test simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--input', '-i', 
                       default='OUT_slump_V1/particle_trajectories.csv',
                       help='Input CSV file with trajectory data (default: OUT_slump_V1/particle_trajectories.csv)')
    parser.add_argument('--output', '-o', 
                       help='Output file for saving plots (optional, will display if not provided)')
    parser.add_argument('--particles', '-p', 
                       help='Comma-separated list of particle IDs to visualize (e.g., 0,5,10). If not provided, all particles are shown.')
    parser.add_argument('--plot-type', '-t', 
                       choices=['3d', '2d', 'time', 'all'],
                       default='all',
                       help='Type of plot to generate (default: all)')
    
    args = parser.parse_args()
    
    # Parse particle IDs if provided
    particle_ids = None
    if args.particles:
        try:
            particle_ids = [int(pid.strip()) for pid in args.particles.split(',')]
        except ValueError:
            print("Error: Particle IDs must be comma-separated integers")
            return 1
    
    # Load data
    try:
        df = load_trajectory_data(args.input)
        print(f"Loaded trajectory data from {args.input}")
        print(f"Data contains {len(df)} data points for {df['particle_id'].nunique()} particles")
        print(f"Time range: {df['time'].min():.3f} to {df['time'].max():.3f}")
    except FileNotFoundError:
        print(f"Error: File {args.input} not found")
        return 1
    except Exception as e:
        print(f"Error loading data: {e}")
        return 1
    
    # Generate plots
    if args.plot_type in ['3d', 'all']:
        output_3d = args.output if args.output and args.plot_type == '3d' else None
        plot_3d_trajectories(df, particle_ids, output_3d)
    
    if args.plot_type in ['2d', 'all']:
        output_2d = args.output if args.output and args.plot_type == '2d' else None
        plot_2d_projections(df, particle_ids, output_2d)
    
    if args.plot_type in ['time', 'all']:
        output_time = args.output if args.output and args.plot_type == 'time' else None
        plot_time_evolution(df, particle_ids, output_time)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
