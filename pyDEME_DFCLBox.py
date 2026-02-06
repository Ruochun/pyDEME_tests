import DEME
from DEME import HCPSampler

import numpy
import os
import time

import pychrono as chrono

# parse in a string
import argparse
parser = argparse.ArgumentParser()
# parser.add_argument("--in_name", "-i", help="input file name(s) for particles", type=str, nargs='+', required=True)
parser.add_argument("--out_name", "-o", help="output file name", type=str, default="output.txt")
parser.add_argument("--settle", "-s", help="before slumping, do NOT settle the particles to low KE", default=True, action='store_false')
parser.add_argument("--particle_file", "-p", help="input particle files", type=str, default="./DFCgeo000-data-particles-LBOX_SCC1.dat")
args = parser.parse_args()
# in_name = args.in_name
out_name = args.out_name
settle = args.settle
particle_file = args.particle_file

# Concrete properties
h_layer=4.0
minD=5.0
maxD=10.0
cement=571
WtoC=0.45 
AtoC=2.62 
rho_c=3150
rho_w=1000
vair=0.03
nF=0.5	

# volume of specimen
specimenVol=600*200*100
Vol=specimenVol

va=1.0-cement/rho_c-(WtoC*cement)/rho_w-vair
va0=(1-pow((minD/maxD),nF))*va
Va0=va0*Vol
# rho=cement*(1.0+WtoC+AtoC)*1E-12
rho=1.47919e-9 

# Define material parameters of fresh concrete
mortar_layer=h_layer
ENm=2.0  # 4.0E-2;
ENa=100 
h=mortar_layer	
alpha=0.25
beta=0.5
np=1.0
sgmTmax=2.0E-4
sgmTau0=1E-7
kappa0=100.
eta_inf=2.0E-5 # can be used to slow down the pile collapse

# Define some families
FIXED_FAM = 10
TEMP_FAM = 11

# These are placeholders
flocbeta = 0.01
flocm = 1.0
flocTcr = 100.
lambda_init = 2.0

def AddLBox(sys, mat, mesh_dir):
	Box = sys.AddWavefrontMeshObject(mesh_dir, mat)
	Box.SetMass(1)
	# q = chrono.QuatFromAngleAxis(-(chrono.CH_PI_2), chrono.VECT_X)
	# Box.SetInitQuat([q.e1, q.e2, q.e3, q.e0])
	# Box.SetInitPos([0,133,0])
	Box.SetFamily(FIXED_FAM)

	# Also add the three cylinders
	cyls = sys.AddExternalObject()
	pos1 = [56., 0., 0.]
	rad = 6.
	cyls.AddCylinder(pos1, [0, 1, 0], rad, mat, -1)

	pos2 = [56., 0., -53.]
	cyls.AddCylinder(pos2, [0, 1, 0], rad, mat, -1)

	pos3 = [56., 0., 53.]
	cyls.AddCylinder(pos3, [0, 1, 0], rad, mat, -1)
	cyls.SetFamily(FIXED_FAM)
	return sys.Track(Box)

def ReadDFCparticles(sys, file_name, rho, mat, h_layer):

	particles = []

	with open(file_name, "r") as f:
		# read in a file (ignoring the lines with // comments)
		for line in f:
			if line.startswith("//"):
				continue
			parts = line.split()
			if len(parts) == 5:
				n = int(parts[0])
				x = float(parts[1])
				y = float(parts[2])
				z = float(parts[3])
				d = float(parts[4])
				particles.append([n, x, y, z, d])

	q = chrono.QuatFromAngleAxis(-(chrono.CH_PI_2), chrono.VECT_X)
	# ChVector3d t(0,0,0)
	csys = chrono.ChCoordsysd(chrono.ChVector3d(0, 0, 0), q)
	
	my_templates = []
	body_points_new = []
	totVol = 0.0
	Vols = []
	radii = []
	for particle in particles:	
		# Our input file loads r directly
		radius=particle[4]+h_layer
		radii.append(radius)
		Vol = 4/3*chrono.CH_PI*pow(radius, 3)
		# mass = rho/2*Vol
		totVol += Vol
		Vols.append(Vol)

	# rho_new = rho*specimenVol/totVol
	# print(f"Original rho: {rho}, new rho: {rho_new}")
	print(f"rho is {rho}, mass is around {rho * Vols[0]}")

	i = 0
	for particle in particles:
		x_pos=particle[1]
		y_pos=particle[2]
		z_pos=particle[3]
		radius=radii[i]
       
		pos = chrono.ChVector3d(x_pos, y_pos, z_pos)
		# pos=csys.TransformPointLocalToParent(pos)
		# print("pos: ", pos, "\t", "rotated pos: ", pos)

		mass = rho * Vols[i]

		# create templates for each sphere
		my_templates.append(sys.LoadSphereType(
            mass,
            radius,
            mat
        ))
		body_points_new.append([pos.x, pos.y, pos.z])

		i += 1

	# Add the particles to the system
	added = sys.AddClumps(my_templates, body_points_new)
	return sys.Track(added)
		      

if __name__ == "__main__":

	# get the current directory
	current_dir = os.getcwd()

	out_dir = current_dir + "/OUT_LBox/"; 
	# Create oputput directories
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	else:
		os.system("rm -r "+out_dir+"/*")
	print("Output directory: ", out_dir)

	DEMSim = DEME.DEMSolver()
	DEMSim.SetVerbosity("INFO")
	DEMSim.SetOutputFormat("CSV")
	DEMSim.SetOutputContent(["ABSV", "XYZ"])
	DEMSim.SetMeshOutputFormat("VTK")
	# This is needed, as contact wildcards are extracted later on, so we need to notify the solver to include them
	DEMSim.SetContactOutputContent(["OWNER", "CNT_WILDCARD"])
	DEMSim.SetJitifyClumpTemplates(False)
	DEMSim.SetJitifyMassProperties(False)
	DEMSim.UseCubForceCollection()
	DEMSim.InstructBoxDomainDimension(2000.0, 2000.0, 2000.0)
	DEMSim.SetErrorOutAvgContacts(250)

	DEMSim.SetGravitationalAcceleration([0, -9810, 0])

	floor_mat = DEMSim.LoadMaterial({"ENm_mat":ENm, "ENa_mat":ENa, "mortar_layer_mat":mortar_layer, 
										"alpha_mat":alpha, "beta_mat":beta, "np_mat":np, "sgmTmax_mat":sgmTmax, 
										"sgmTau0_mat":sgmTau0, "kappa0_mat":kappa0, "eta_inf_mat":eta_inf, 
										"flocbeta_mat":flocbeta, "flocTcr_mat":flocTcr, "flocm_mat":flocm, "lambda_init_mat":lambda_init})
	
	AddLBox(DEMSim, floor_mat, current_dir + "/LBox_NoCyl.obj")
	# AddLBox(DEMSim, floor_mat, current_dir + "/LBox_Mesh.obj")
	
	# Optionally we can add a gate plane here
	# gate = DEMSim.AddPlane([56.0, 0.0, 0.0], [-1.0, 0.0, 0.0], floor_mat)
	# gate.SetFamily(TEMP_FAM)

	# Let family 10,11 be fixed
	DEMSim.SetFamilyFixed(TEMP_FAM)
	DEMSim.SetFamilyFixed(FIXED_FAM)

	# Change default contact force algorithm to fresh concrete ones
	my_force_model = DEMSim.ReadContactForceModel("DFCModel.cu")
	# Those following lines are needed. We must let the solver know that those var names are history variable etc.
	my_force_model.SetMustHaveMatProp(set([]))
	my_force_model.SetMustPairwiseMatProp(set([]))
	my_force_model.SetPerContactWildcards(
		set(["contact_info_step_time", "contact_info_lambda", "contact_info_strain_x", "contact_info_strain_y", "contact_info_strain_z"])
	)
	my_force_model.SetPerOwnerWildcards(set(["ParticleThixoAvg"]))
	# Add a utility function for the force model to use
	my_force_model.DefineCustomModelPrerequisites("""
													void __device__ XdirToDxDyDz(const float3 &Vxdir,
																				const float3 &Vsingular, 
																				float3 &Vx,
																				float3 &Vy, 
																				float3 &Vz) {
														const float3 mVnull = make_float3(0);
														float zlen;
														if (length(Vxdir) < 1e-8) {
														Vx = make_float3(1, 0, 0); // default X direction
														} else {
														Vx = normalize(Vxdir);
														}
														Vz = cross(Vx, Vsingular);
														zlen = length(Vz);

														if (zlen < 1e-6) {
														float3 mVsingualr;
														if (abs(Vsingular.z) < 0.9)
															mVsingualr = make_float3(0, 0, 1);
														if (abs(Vsingular.y) < 0.9)
															mVsingualr = make_float3(0, 1, 0);
														if (abs(Vsingular.x) < 0.9)
															mVsingualr = make_float3(1, 0, 0);
														Vz = cross(Vx, mVsingualr);
														zlen = length(Vz);
														}

														Vz = normalize(Vz);
														Vy = cross(Vz, Vx);
													};
													"""
													)

	# Read and load predefined particle list
	particle_tracker = ReadDFCparticles(DEMSim, particle_file, rho, floor_mat, mortar_layer)

	KE_finder = DEMSim.CreateInspector("clump_kinetic_energy")
	# maxZ_finder = DEMSim.CreateInspector("clump_max_z")

	timestep = 1.0e-8
	DEMSim.SetInitTimeStep(timestep)
	DEMSim.SetMaxVelocity(400000.)
	DEMSim.SetErrorOutVelocity(50000.)

	# DEMSim.SetCDUpdateFreq(5)
	# DEMSim.DisableAdaptiveUpdateFreq()
	DEMSim.SetInitBinNumTarget(int(1e7))
	DEMSim.DisableAdaptiveBinSize()
	DEMSim.Initialize()

	max_step_size = 1e-5

	stepsize_increase_mult = 1.01
	currframe = 0
	t = 0.0
	
	step_num = 0

	frame_time = 1. / 20.
	out_time = 0.
	
	# Track positions of ~50 particles throughout simulation
	num_tracked_particles = 50
	trajectory_data = []  # Will store [time, particle_id, x, y, z] for each frame

	if settle:

		# 1. phase: waiting with a plate as the gate
		print("\n 1. phase: waiting with a plate as the gate \n")

		while t < 5.0:
			if (step_num % 20 == 0) and (timestep < max_step_size):
				timestep = min(timestep*stepsize_increase_mult,max_step_size)
				DEMSim.UpdateStepSize(timestep)

			if (t >= out_time):
				print(f"============ Time {t}, settling =============")
				print(f"Current step size: {timestep}")
				print(f"Outputting frame: {currframe}")
				filename = os.path.join(out_dir, f"slump_va_M{currframe:06d}.csv")
				meshname = os.path.join(out_dir, f"slump_mesh{currframe:06d}.vtk")
				DEMSim.WriteSphereFile(filename)
				# DEMSim.WriteMeshFile(meshname)
				print(f"Current number of contacts: {DEMSim.GetNumContacts()}")
				# print(f"Bin num: {DEMSim.GetBinNum()}")
				# print(f"Bin size: {DEMSim.GetBinSize()}")
				# DEMSim.ShowThreadCollaborationStats()
				# DEMSim.ShowMemStats()
				# DEMSim.PrintKinematicScratchSpaceUsage()
				
				# Record positions of tracked particles
				positions = particle_tracker.Positions()
				for i in range(num_tracked_particles):
					trajectory_data.append([t, i, positions[i][0], positions[i][1], positions[i][2]])

				currframe += 1
				out_time += frame_time

				ke = KE_finder.GetValue()
				print(f"Kinetic energy: {ke}")
				if (ke < 0.5) and (step_num > 4000):
					print("Kinetic energy is low, finishing up the settling phase.")
					break

			DEMSim.DoDynamics(timestep)
			t += timestep
			step_num += 1
    

	# 2. phase: remove containers and allow flow 
	
	print("\n 2. phase: flow  \n")

	phase1_end = t
	timestep=1e-5
	stepsize_increase_mult = 1.01

	DEMSim.DoDynamicsThenSync(0)
	DEMSim.UpdateStepSize(timestep)
	
	# Disable contact between some entities and the particles
	DEMSim.DisableContactBetweenFamilies(TEMP_FAM, 0)

	start_time = time.time()
	while t < phase1_end + 20.:
		if (step_num % 10 == 0) and (timestep < max_step_size):
			timestep = min(timestep*stepsize_increase_mult,max_step_size)
			DEMSim.UpdateStepSize(timestep)

		if (t >= out_time):
			print(f"============ Time {t}, flowing =============")
			print(f"Current step size: {timestep}")
			print(f"Outputting frame: {currframe}")
			filename = os.path.join(out_dir, f"slump_va_M{currframe:06d}.csv")
			meshname = os.path.join(out_dir, f"slump_mesh{currframe:06d}.vtk")
			DEMSim.WriteSphereFile(filename)
			# DEMSim.WriteMeshFile(meshname)
			print(f"Current number of contacts: {DEMSim.GetNumContacts()}")
			# print(f"Bin num: {DEMSim.GetBinNum()}")
			# print(f"Bin size: {DEMSim.GetBinSize()}")
			# DEMSim.ShowThreadCollaborationStats()
			# DEMSim.ShowMemStats()
			# DEMSim.PrintKinematicScratchSpaceUsage()
			
			# Record positions of tracked particles
			positions = particle_tracker.Positions()
			for i in range(num_tracked_particles):
				trajectory_data.append([t, i, positions[i][0], positions[i][1], positions[i][2]])

			currframe += 1
			out_time += frame_time

			# maxZ = maxZ_finder.GetValue()
			# print(f"Highest point: {maxZ}")

		DEMSim.DoDynamics(timestep)
		t += timestep
		step_num += 1

	end_time = time.time()
	print(f"Total simulation time: {end_time - start_time} seconds")

	# output the xyz of all particles to file
	filename = os.path.join(out_dir, out_name)
	xyz = particle_tracker.Positions()
	with open(filename, "w") as f:
		f.write("# X,Y,Z\n")
		for pos in xyz:
			f.write(f"{pos[0]},{pos[1]},{pos[2]}\n")
	
	# Save trajectory data for tracked particles
	trajectory_filename = os.path.join(out_dir, "particle_trajectories.csv")
	print(f"Saving trajectory data for {num_tracked_particles} particles to {trajectory_filename}")
	with open(trajectory_filename, "w") as f:
		f.write("time,particle_id,x,y,z\n")
		for data_point in trajectory_data:
			f.write(f"{data_point[0]},{data_point[1]},{data_point[2]},{data_point[3]},{data_point[4]}\n")
	print(f"Trajectory data saved. Total data points: {len(trajectory_data)}")
