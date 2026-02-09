import ovito
from ovito import data
from ovito.pipeline import Pipeline
from ovito.io import import_file, export_file
from ovito.modifiers import ConstructSurfaceModifier

import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

import ovito._extensions.pyscript
# import argparse
import sys
def surface_reconstrution(
        input_file, 
        output_file, 
        radius, 
        smoothing_level,
        compute_distances,       
        select_surface_particles, 
    ):
    """
    Analyze the surface mesh of a given particle dataset using OVITO's ConstructSurfaceModifier.

    Parameters:
    - input_file (str): Path to the input `.xyz` file.
    - output_file (str): Path for saving the output `.vtk` file.
    - radius (float): Alpha shape radius for the surface modifier.
    - smoothing_level (int): Level of smoothing applied to the surface mesh.
    - select_surface_particles (bool): If True, select the surface particles.
    
    Returns:
    - None
    """

    # Import the data file into the OVITO pipeline
    # pipeline = import_file(input_file, columns=['Particle Identifier', 'Position.X', 'Position.Y', 'Position.Z'])
    pipeline = import_file(input_file, columns=['Position.X', 'Position.Y', 'Position.Z'])

    # Create the ConstructSurfaceModifier with given parameters
    modifier = ConstructSurfaceModifier()
    modifier.method = ConstructSurfaceModifier.Method.AlphaShape
    modifier.radius = radius
    modifier.smoothing_level = smoothing_level
    modifier.compute_distances = compute_distances
    modifier.select_surface_particles = select_surface_particles
    modifier.identify_regions = True
    modifier.grid_resolution=int(1e6)

    # Insert the modifier into the pipeline
    pipeline.modifiers.append(modifier)

    # Compute the pipeline to generate the surface mesh
    data = pipeline.compute()

    # Access and print data size information
    num_particles = data.particles.count
    print(f"Number of particles (atoms): {num_particles}")

    print("Aggregate quantities:")
    print(f"Surface area: {data.attributes['ConstructSurfaceMesh.surface_area']}")
    print(f"Solid volume: {data.attributes['ConstructSurfaceMesh.filled_volume']}")

    print("Volumetric regions:")
    regions = data.surfaces['surface'].regions
    # for vol, area, filled, exterior in zip(regions['Volume'], regions['Surface Area'], regions['Filled'], regions['Exterior']):
    #     print(f"  volume={vol}, surface area={area}, filled={filled}, exterior={exterior}")

    # Check for the surface mesh in the computed data
    if 'surface' in data.surfaces:
        surface_mesh = data.surfaces['surface']
        
        # Print information about volumetric regions
        print("Volumetric regions:")
        regions = surface_mesh.regions
        for vol, area, filled, exterior in zip(regions['Volume'], regions['Surface Area'], regions['Filled'], regions['Exterior']):
            print(f"  Volume={vol}, Surface Area={area}, Filled={filled}, Exterior={exterior}")
    else:
        print("No surface mesh was generated.")
    
    # Export the surface mesh to a file
    export_file(pipeline, output_file+".vtk", "vtk/trimesh")
    if compute_distances == False or select_surface_particles == False:
        export_file(pipeline, output_file+".xyz", "xyz", columns = ["Particle Identifier", "Position.X", "Position.Y", "Position.Z"])
    else:
        export_file(pipeline, output_file+".xyz", "xyz", columns = ["Particle Identifier", "Position.X", "Position.Y", "Position.Z","Surface Distance","Selection"])
    print(f"Surface mesh exported to {output_file}")
    

# Check if the script is run directly (instead of being imported as a module)
if __name__ == "__main__":
    print(sys.argv[0])
    if len(sys.argv) < 5 or len(sys.argv) > 7:
        print("Please enter input_file, output_file, radius, smoothing_level and optional(compute_distances,select_surface_particles)")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    radius = float(sys.argv[3])
    smoothing_level = int(sys.argv[4])
    if len(sys.argv) < 6:
        compute_distances = True
    else: compute_distances = bool(int(sys.argv[5]))
    if len(sys.argv) < 7:
        select_surface_particles = True
    else: select_surface_particles = bool(int(sys.argv[6]))

    print(f"inputs are radius={radius},smoothing_level={smoothing_level},compute_distances = {compute_distances},select_surface_particles={select_surface_particles}")
    surface_reconstrution(
        input_file,
        output_file,
        radius,
        smoothing_level,
        compute_distances,
        select_surface_particles,
    )
    