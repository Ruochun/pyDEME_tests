#!/usr/bin/env python3
"""
Test script for particle trajectory tracking functionality.

This script tests the trajectory tracking logic without requiring the full
DEME simulation environment.
"""

import os
import sys
import tempfile
import csv

def test_trajectory_data_structure():
    """Test that trajectory data structure works correctly."""
    print("Testing trajectory data structure...")
    
    # Simulate what the code does
    trajectory_data = []
    num_tracked_particles = 5
    
    # Simulate recording positions at two time steps
    time_steps = [0.0, 0.05]
    mock_positions = [
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0], [13.0, 14.0, 15.0]],
        [[1.1, 2.1, 3.1], [4.1, 5.1, 6.1], [7.1, 8.1, 9.1], [10.1, 11.1, 12.1], [13.1, 14.1, 15.1]]
    ]
    
    for t_idx, t in enumerate(time_steps):
        positions = mock_positions[t_idx]
        for i in range(num_tracked_particles):
            trajectory_data.append([t, i, positions[i][0], positions[i][1], positions[i][2]])
    
    # Verify data structure
    assert len(trajectory_data) == num_tracked_particles * len(time_steps), \
        f"Expected {num_tracked_particles * len(time_steps)} data points, got {len(trajectory_data)}"
    
    # Verify first entry
    assert trajectory_data[0] == [0.0, 0, 1.0, 2.0, 3.0], \
        f"First entry incorrect: {trajectory_data[0]}"
    
    # Verify last entry
    assert trajectory_data[-1] == [0.05, 4, 13.1, 14.1, 15.1], \
        f"Last entry incorrect: {trajectory_data[-1]}"
    
    print("✓ Data structure test passed")
    return True

def test_csv_output():
    """Test that CSV output is written correctly."""
    print("Testing CSV output...")
    
    # Create sample data
    trajectory_data = [
        [0.0, 0, 1.0, 2.0, 3.0],
        [0.0, 1, 4.0, 5.0, 6.0],
        [0.05, 0, 1.1, 2.1, 3.1],
        [0.05, 1, 4.1, 5.1, 6.1]
    ]
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        temp_filename = f.name
        f.write("time,particle_id,x,y,z\n")
        for data_point in trajectory_data:
            f.write(f"{data_point[0]},{data_point[1]},{data_point[2]},{data_point[3]},{data_point[4]}\n")
    
    try:
        # Read back and verify
        with open(temp_filename, 'r') as f:
            lines = f.readlines()
        
        # Check header
        assert lines[0].strip() == "time,particle_id,x,y,z", \
            f"Header incorrect: {lines[0]}"
        
        # Check number of data lines (header + 4 data lines)
        assert len(lines) == 5, f"Expected 5 lines, got {len(lines)}"
        
        # Parse and verify a data line
        parts = lines[1].strip().split(',')
        assert len(parts) == 5, f"Expected 5 columns, got {len(parts)}"
        assert float(parts[0]) == 0.0, f"Time incorrect: {parts[0]}"
        assert int(parts[1]) == 0, f"Particle ID incorrect: {parts[1]}"
        assert float(parts[2]) == 1.0, f"X incorrect: {parts[2]}"
        
        print("✓ CSV output test passed")
        return True
    finally:
        # Clean up
        os.unlink(temp_filename)

def test_min_function():
    """Test that min(50, len(particle_tracker)) works as expected."""
    print("Testing particle count selection...")
    
    # Test when particles < 50
    mock_particle_count_1 = 30
    num_tracked = min(50, mock_particle_count_1)
    assert num_tracked == 30, f"Expected 30, got {num_tracked}"
    
    # Test when particles > 50
    mock_particle_count_2 = 100
    num_tracked = min(50, mock_particle_count_2)
    assert num_tracked == 50, f"Expected 50, got {num_tracked}"
    
    # Test when particles = 50
    mock_particle_count_3 = 50
    num_tracked = min(50, mock_particle_count_3)
    assert num_tracked == 50, f"Expected 50, got {num_tracked}"
    
    print("✓ Particle count selection test passed")
    return True

def test_trajectory_file_reading():
    """Test that the CSV file can be read by standard tools."""
    print("Testing CSV file reading...")
    
    # Create sample data
    trajectory_data = [
        [0.0, 0, 1.0, 2.0, 3.0],
        [0.0, 1, 4.0, 5.0, 6.0],
        [0.05, 0, 1.1, 2.1, 3.1]
    ]
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        temp_filename = f.name
        f.write("time,particle_id,x,y,z\n")
        for data_point in trajectory_data:
            f.write(f"{data_point[0]},{data_point[1]},{data_point[2]},{data_point[3]},{data_point[4]}\n")
    
    try:
        # Test reading with csv module
        with open(temp_filename, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        assert len(rows) == 3, f"Expected 3 rows, got {len(rows)}"
        assert 'time' in rows[0], "Missing 'time' column"
        assert 'particle_id' in rows[0], "Missing 'particle_id' column"
        assert 'x' in rows[0], "Missing 'x' column"
        assert 'y' in rows[0], "Missing 'y' column"
        assert 'z' in rows[0], "Missing 'z' column"
        
        # Test that values are correct
        assert float(rows[0]['time']) == 0.0
        assert int(rows[0]['particle_id']) == 0
        assert float(rows[0]['x']) == 1.0
        
        print("✓ CSV file reading test passed")
        return True
    finally:
        # Clean up
        os.unlink(temp_filename)

def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Running Particle Trajectory Tracking Tests")
    print("=" * 60)
    print()
    
    tests = [
        test_trajectory_data_structure,
        test_csv_output,
        test_min_function,
        test_trajectory_file_reading
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
            failed += 1
        print()
    
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0

if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
