# Particle Trajectory Tracking

## Overview

The `pyDEME_DFCSlump.py` script now tracks the xyz positions of approximately 50 particles throughout the simulation and saves them stepwise for later visualization and analysis.

## What Changed

### Modified Files
- **pyDEME_DFCSlump.py**: 
  - Added tracking of ~50 particles (the first 50 particles in the simulation)
  - Records particle positions at each output frame (20 frames per second of simulation time)
  - Saves trajectory data to `OUT_slump_V1/particle_trajectories.csv`

### New Files
- **visualize_trajectories.py**: A visualization script to create plots from the trajectory data

## Output Format

The trajectory data is saved in CSV format with the following columns:
- `time`: Simulation time in seconds
- `particle_id`: ID of the tracked particle (0-49)
- `x`: X position
- `y`: Y position  
- `z`: Z position

Example data:
```csv
time,particle_id,x,y,z
0.0,0,-45.2341,125.678,23.4567
0.0,1,-38.1234,130.234,19.8765
0.05,0,-45.2298,125.654,23.4512
0.05,1,-38.1201,130.211,19.8721
...
```

## Usage

### Running the Simulation

Run the simulation as before:
```bash
python pyDEME_DFCSlump.py --out_name output.txt -s 
```

This will create the output directory `OUT_slump_V1/` containing:
- Full simulation snapshots: `slump_va_M000000.csv`, `slump_va_M000001.csv`, etc.
- Mesh files: `slump_mesh000000.vtk`, `slump_mesh000001.vtk`, etc.
- Final positions of all particles: `output.txt` (or your specified output name)
- **NEW**: Particle trajectories: `particle_trajectories.csv`

### Visualizing Trajectories

The `visualize_trajectories.py` script provides several visualization options:

#### Basic Usage (requires matplotlib and pandas)
```bash
# Install dependencies if needed
pip install matplotlib pandas

# Visualize all tracked particles (creates 3D plot, 2D projections, and time evolution)
python visualize_trajectories.py --input OUT_slump_V1/particle_trajectories.csv

# Visualize specific particles
python visualize_trajectories.py --input OUT_slump_V1/particle_trajectories.csv --particles 0,5,10,15,20

# Save to a file instead of displaying
python visualize_trajectories.py --input OUT_slump_V1/particle_trajectories.csv --output trajectories.png --plot-type 3d
```

#### Visualization Types

1. **3D Trajectory Plot** (`--plot-type 3d`):
   - Shows the full 3D path of particles
   - Start positions marked with circles
   - End positions marked with squares

2. **2D Projections** (`--plot-type 2d`):
   - Three views: XY (top), XZ (front), YZ (side)
   - Useful for understanding motion in specific planes

3. **Time Evolution** (`--plot-type time`):
   - Shows X, Y, Z coordinates vs time
   - Helps understand the dynamics of particle motion

4. **All plots** (`--plot-type all`, default):
   - Generates all three plot types

## Custom Analysis

You can also write your own analysis scripts by loading the CSV file:

```python
import pandas as pd

# Load trajectory data
df = pd.read_csv('OUT_slump_V1/particle_trajectories.csv')

# Get data for a specific particle
particle_0 = df[df['particle_id'] == 0]

# Calculate displacement over time
initial_pos = particle_0.iloc[0][['x', 'y', 'z']].values
final_pos = particle_0.iloc[-1][['x', 'y', 'z']].values
displacement = final_pos - initial_pos

# Calculate velocity (approximate)
particle_0['vx'] = particle_0['x'].diff() / particle_0['time'].diff()
particle_0['vy'] = particle_0['y'].diff() / particle_0['time'].diff()
particle_0['vz'] = particle_0['z'].diff() / particle_0['time'].diff()

# And more...
```

## Technical Details

### Number of Tracked Particles
The code tracks exactly 50 particles (or fewer if the total number of particles is less than 50). The first 50 particles from the input data are selected for tracking.

### Tracking Frequency
Positions are recorded at the same frequency as the full simulation output, which is 20 frames per second of simulation time (configurable via `frame_time` in the code).

### Performance Impact
The trajectory tracking has minimal performance impact as it only reads particle positions at output frames (which are already being computed) and appends them to a list in memory.

## Extending the Code

To track different particles or modify the tracking behavior:

1. **Change number of tracked particles**: Modify the line `num_tracked_particles = min(50, len(particle_tracker))` in the initialization section before the main simulation loop

2. **Track specific particles**: Instead of tracking the first N particles, you can modify the loop to track specific particle IDs:
   ```python
   # Track specific particles (e.g., particles 10, 20, 30, ...)
   tracked_particle_ids = [10, 20, 30, 40, 50]
   for particle_id in tracked_particle_ids:
       if particle_id < len(positions):
           trajectory_data.append([t, particle_id, positions[particle_id][0], 
                                  positions[particle_id][1], positions[particle_id][2]])
   ```

3. **Change tracking frequency**: Modify `frame_time = 1. / 20.` to change how often positions are recorded (e.g., `frame_time = 1. / 10.` for 10 frames per second)

4. **Add additional data**: You can track velocity, acceleration, or other properties by querying the appropriate methods on `particle_tracker`
