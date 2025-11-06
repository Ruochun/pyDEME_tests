#!/usr/bin/env python3
"""
Test script for compare_trajectory_runs.py functionality.

This script tests trajectory comparison logic without requiring actual
simulation data.
"""

import os
import sys
import tempfile
import math
import traceback

# Try to import compare_trajectory_runs module
try:
    from compare_trajectory_runs import (
        load_trajectories,
        calculate_path_length,
        interpolate_trajectory_at_time,
        compute_trajectory_deviation,
        compute_trajectory_deviation_evolution,
        match_particles_by_initial_position
    )
except ImportError as e:
    print(f"Error: Could not import compare_trajectory_runs module: {e}")
    print("Make sure compare_trajectory_runs.py is in the same directory.")
    sys.exit(1)

def create_test_trajectory_file(filename, trajectories):
    """
    Create a test trajectory CSV file.
    
    Args:
        filename: Output file path
        trajectories: Dict of {particle_id: [(time, x, y, z), ...]}
    """
    with open(filename, 'w') as f:
        f.write("time,particle_id,x,y,z\n")
        for particle_id, points in trajectories.items():
            for time, x, y, z in points:
                f.write(f"{time},{particle_id},{x},{y},{z}\n")

def test_load_trajectories():
    """Test loading trajectories from CSV file."""
    print("Testing trajectory loading...")
    
    # Create test data
    test_trajs = {
        0: [(0.0, 1.0, 2.0, 3.0), (0.1, 1.1, 2.1, 3.1)],
        1: [(0.0, 4.0, 5.0, 6.0), (0.1, 4.1, 5.1, 6.1)]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        temp_file = f.name
    
    try:
        create_test_trajectory_file(temp_file, test_trajs)
        
        # Load and verify
        loaded = load_trajectories(temp_file)
        
        assert len(loaded) == 2, f"Expected 2 particles, got {len(loaded)}"
        assert 0 in loaded and 1 in loaded, "Missing particle IDs"
        assert len(loaded[0]) == 2, f"Expected 2 points for particle 0, got {len(loaded[0])}"
        
        # Check first point of first particle
        t, x, y, z = loaded[0][0]
        assert abs(t - 0.0) < 1e-6, f"Time mismatch: {t}"
        assert abs(x - 1.0) < 1e-6, f"X mismatch: {x}"
        assert abs(y - 2.0) < 1e-6, f"Y mismatch: {y}"
        assert abs(z - 3.0) < 1e-6, f"Z mismatch: {z}"
        
        print("✓ Trajectory loading test passed")
        return True
    finally:
        os.unlink(temp_file)

def test_calculate_path_length():
    """Test path length calculation."""
    print("Testing path length calculation...")
    
    # Simple straight-line trajectory
    traj1 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0, 0.0),
        (2.0, 2.0, 0.0, 0.0)
    ]
    length1 = calculate_path_length(traj1)
    assert abs(length1 - 2.0) < 1e-6, f"Expected length 2.0, got {length1}"
    
    # 3D trajectory
    traj2 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0, 1.0)
    ]
    length2 = calculate_path_length(traj2)
    expected = math.sqrt(3)  # sqrt(1^2 + 1^2 + 1^2)
    assert abs(length2 - expected) < 1e-6, f"Expected length {expected}, got {length2}"
    
    # Empty trajectory
    traj3 = []
    length3 = calculate_path_length(traj3)
    assert length3 == 0.0, f"Expected length 0.0 for empty trajectory, got {length3}"
    
    print("✓ Path length calculation test passed")
    return True

def test_interpolate_trajectory():
    """Test trajectory interpolation."""
    print("Testing trajectory interpolation...")
    
    # Linear trajectory from (0,0,0) to (10,10,10) over time 0 to 10
    traj = [
        (0.0, 0.0, 0.0, 0.0),
        (10.0, 10.0, 10.0, 10.0)
    ]
    
    # Test interpolation at midpoint
    x, y, z = interpolate_trajectory_at_time(traj, 5.0)
    assert abs(x - 5.0) < 1e-6, f"Expected x=5.0, got {x}"
    assert abs(y - 5.0) < 1e-6, f"Expected y=5.0, got {y}"
    assert abs(z - 5.0) < 1e-6, f"Expected z=5.0, got {z}"
    
    # Test at start
    x, y, z = interpolate_trajectory_at_time(traj, 0.0)
    assert abs(x - 0.0) < 1e-6, f"Expected x=0.0, got {x}"
    
    # Test at end
    x, y, z = interpolate_trajectory_at_time(traj, 10.0)
    assert abs(x - 10.0) < 1e-6, f"Expected x=10.0, got {x}"
    
    # Test before start (should return start position)
    x, y, z = interpolate_trajectory_at_time(traj, -1.0)
    assert abs(x - 0.0) < 1e-6, f"Expected x=0.0 for time before start, got {x}"
    
    # Test after end (should return end position)
    x, y, z = interpolate_trajectory_at_time(traj, 15.0)
    assert abs(x - 10.0) < 1e-6, f"Expected x=10.0 for time after end, got {x}"
    
    print("✓ Trajectory interpolation test passed")
    return True

def test_compute_trajectory_deviation():
    """Test trajectory deviation computation."""
    print("Testing trajectory deviation computation...")
    
    # Identical trajectories
    traj1 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0, 1.0),
        (2.0, 2.0, 2.0, 2.0)
    ]
    traj2 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0, 1.0),
        (2.0, 2.0, 2.0, 2.0)
    ]
    avg_dev, max_dev, num_points = compute_trajectory_deviation(traj1, traj2)
    assert abs(avg_dev) < 1e-6, f"Expected 0 deviation for identical trajectories, got {avg_dev}"
    assert abs(max_dev) < 1e-6, f"Expected 0 max deviation, got {max_dev}"
    assert num_points > 0, f"Expected positive number of points, got {num_points}"
    
    # Parallel trajectories with constant offset
    traj3 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0, 0.0),
        (2.0, 2.0, 0.0, 0.0)
    ]
    traj4 = [
        (0.0, 0.0, 1.0, 0.0),  # Offset by 1 in y direction
        (1.0, 1.0, 1.0, 0.0),
        (2.0, 2.0, 1.0, 0.0)
    ]
    avg_dev, max_dev, num_points = compute_trajectory_deviation(traj3, traj4)
    assert abs(avg_dev - 1.0) < 1e-6, f"Expected deviation 1.0, got {avg_dev}"
    assert abs(max_dev - 1.0) < 1e-6, f"Expected max deviation 1.0, got {max_dev}"
    
    print("✓ Trajectory deviation computation test passed")
    return True

def test_match_particles():
    """Test particle matching by initial position."""
    print("Testing particle matching...")
    
    # Create two sets of trajectories with known correspondence
    baseline = {
        0: [(0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)],
        1: [(0.0, 5.0, 5.0, 5.0), (1.0, 6.0, 6.0, 6.0)]
    }
    
    # Comparison has same initial positions but different particle IDs
    comparison = {
        10: [(0.0, 5.0, 5.0, 5.0), (1.0, 6.1, 6.1, 6.1)],  # Matches baseline[1]
        20: [(0.0, 0.0, 0.0, 0.0), (1.0, 1.1, 1.1, 1.1)]   # Matches baseline[0]
    }
    
    matching = match_particles_by_initial_position(baseline, comparison)
    
    assert len(matching) == 2, f"Expected 2 matches, got {len(matching)}"
    assert matching[0] == 20, f"Expected baseline[0] -> comparison[20], got {matching[0]}"
    assert matching[1] == 10, f"Expected baseline[1] -> comparison[10], got {matching[1]}"
    
    print("✓ Particle matching test passed")
    return True

def test_integration():
    """Integration test with complete workflow."""
    print("Testing complete workflow...")
    
    # Create baseline trajectory file
    baseline_trajs = {
        0: [(0.0, 0.0, 0.0, 0.0), (0.5, 0.5, 0.5, 0.5), (1.0, 1.0, 1.0, 1.0)],
        1: [(0.0, 5.0, 0.0, 0.0), (0.5, 5.5, 0.5, 0.5), (1.0, 6.0, 1.0, 1.0)]
    }
    
    # Create comparison trajectory file with slight deviations
    comparison_trajs = {
        0: [(0.0, 0.0, 0.0, 0.0), (0.5, 0.5, 0.51, 0.5), (1.0, 1.0, 1.02, 1.0)],
        1: [(0.0, 5.0, 0.0, 0.0), (0.5, 5.5, 0.49, 0.5), (1.0, 6.0, 0.98, 1.0)]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        baseline_file = f.name
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        comparison_file = f.name
    
    try:
        create_test_trajectory_file(baseline_file, baseline_trajs)
        create_test_trajectory_file(comparison_file, comparison_trajs)
        
        # Load both files
        baseline_loaded = load_trajectories(baseline_file)
        comparison_loaded = load_trajectories(comparison_file)
        
        # Match particles
        matching = match_particles_by_initial_position(baseline_loaded, comparison_loaded)
        assert len(matching) == 2, f"Expected 2 matches, got {len(matching)}"
        
        # Compute deviations
        deviations = []
        for base_id, comp_id in matching.items():
            avg_dev, max_dev, num_points = compute_trajectory_deviation(
                baseline_loaded[base_id],
                comparison_loaded[comp_id]
            )
            deviations.append(avg_dev)
            assert num_points > 0, "Expected positive number of comparison points"
        
        # Check that deviations are small but non-zero
        assert all(d > 0 for d in deviations), "Expected some deviation"
        assert all(d < 0.1 for d in deviations), "Deviations should be small for this test data"
        
        print("✓ Integration test passed")
        return True
    finally:
        os.unlink(baseline_file)
        os.unlink(comparison_file)

def test_relative_deviation():
    """Test that relative deviation is computed correctly."""
    print("Testing relative deviation computation...")
    
    # Create trajectory with known path length
    traj1 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 10.0, 0.0, 0.0)  # Path length = 10
    ]
    
    # Create similar trajectory with small deviation
    traj2 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 10.0, 0.1, 0.0)  # Deviates by 0.1 in y
    ]
    
    path_length = calculate_path_length(traj1)
    avg_dev, _, _ = compute_trajectory_deviation(traj1, traj2)
    relative_dev = avg_dev / path_length
    
    # Relative deviation should be around 0.5% (0.05 average deviation / 10 path length)
    assert 0.0 < relative_dev < 0.02, f"Expected relative deviation ~0.5%, got {relative_dev*100}%"
    
    print("✓ Relative deviation test passed")
    return True

def test_deviation_evolution():
    """Test deviation evolution computation."""
    print("Testing deviation evolution computation...")
    
    # Create two trajectories with known deviations
    traj1 = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0, 0.0),
        (2.0, 2.0, 0.0, 0.0)
    ]
    
    # Trajectory 2 has constant offset of 0.5 in y
    traj2 = [
        (0.0, 0.0, 0.5, 0.0),
        (1.0, 1.0, 0.5, 0.0),
        (2.0, 2.0, 0.5, 0.0)
    ]
    
    time_points, deviations = compute_trajectory_deviation_evolution(traj1, traj2)
    
    # Should have 3 time points
    assert len(time_points) == 3, f"Expected 3 time points, got {len(time_points)}"
    assert len(deviations) == 3, f"Expected 3 deviations, got {len(deviations)}"
    
    # All deviations should be approximately 0.5 (constant offset)
    for dev in deviations:
        assert abs(dev - 0.5) < 1e-6, f"Expected deviation ~0.5, got {dev}"
    
    # Time points should match input times
    assert abs(time_points[0] - 0.0) < 1e-6, f"First time point incorrect: {time_points[0]}"
    assert abs(time_points[1] - 1.0) < 1e-6, f"Second time point incorrect: {time_points[1]}"
    assert abs(time_points[2] - 2.0) < 1e-6, f"Third time point incorrect: {time_points[2]}"
    
    print("✓ Deviation evolution test passed")
    return True

def run_all_tests():
    """Run all tests."""
    print("=" * 80)
    print("Running Trajectory Comparison Tests")
    print("=" * 80)
    print()
    
    tests = [
        test_load_trajectories,
        test_calculate_path_length,
        test_interpolate_trajectory,
        test_compute_trajectory_deviation,
        test_match_particles,
        test_integration,
        test_relative_deviation,
        test_deviation_evolution
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except AssertionError as e:
            print(f"✗ Test failed: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ Test error: {e}")
            traceback.print_exc()
            failed += 1
            failed += 1
        print()
    
    print("=" * 80)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 80)
    
    return failed == 0

if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
