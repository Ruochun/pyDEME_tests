import DEME
from DEME import HCPSampler

import numpy
import os
import time

import pychrono as chrono

# parse in a string
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--out_name", "-o", help="output file name", type=str, default="output.txt")
parser.add_argument("--permutation", "-p", help="amount of random permutation to apply to the initial particles", type=float, default=-1.0)
args = parser.parse_args()
out_name = args.out_name
permutation = args.permutation

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
specimenVol=( chrono.CH_PI*pow(100,2)*600/3- chrono.CH_PI*pow(50,2)*300/3)
Vol=specimenVol

va=1.0-cement/rho_c-(WtoC*cement)/rho_w-vair
va0=(1-pow((minD/maxD),nF))*va
Va0=va0*Vol
rho=cement*(1.0+WtoC+AtoC)*1E-12

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

# Define fixed family
FIXED_FAM = 10

# These are placeholders
flocbeta = 0.01
flocm = 1.0
flocTcr = 100.
lambda_init = 2.0

def AddSlumpCone(sys, mat, density, radius, height, current_dir):
	obj_file = current_dir + "/coni_100_50_300.obj"
	SlumpCone = sys.AddWavefrontMeshObject(obj_file, mat)
	SlumpCone.SetMass(1)
	q = chrono.QuatFromAngleAxis(-(chrono.CH_PI_2), chrono.VECT_X)
	SlumpCone.SetInitQuat([q.e1, q.e2, q.e3, q.e0])
	# SlumpCone.SetInitPos([0,133,0])
	SlumpCone.SetFamily(FIXED_FAM)
	return sys.Track(SlumpCone)

def ReadDFCparticles(sys, file_name, rho, mat, h_layer, permutation=-1.0):

	particles = []

	if len(file_name) == 1:

		with open(file_name[0], "r") as f:
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


	elif len(file_name) == 2:
		i = 0
		with open(file_name[0], "r") as f:
			# read in a file (ignoring the lines with // comments)
			for line in f:
				if line.startswith("//"):
					continue
				parts = line.split()
				if len(parts) == 3:
					x = float(parts[0])
					y = float(parts[1])
					z = float(parts[2])
					d = float(0)
					particles.append([i, x, y, z, d])
					i += 1

		i = 0
		with open(file_name[1], "r") as f:
			# read in a file (ignoring the lines with // comments)
			for line in f:
				if line.startswith("//"):
					continue
				parts = line.split()
				if len(parts) == 5:
					d = float(parts[4])
					particles[i][4] = d
					i += 1
		if i != len(particles):
			raise ValueError("Number of particles in the two files do not match")

	else:
		raise ValueError("file_name must be a string or a list of two strings")

	q = chrono.QuatFromAngleAxis(-(chrono.CH_PI_2), chrono.VECT_X)
	# ChVector3d t(0,0,0)
	csys = chrono.ChCoordsysd(chrono.ChVector3d(0, 0, 0), q)
	
	my_templates = []
	body_points_new = []
	totVol = 0.0
	Vols = []
	radii = []
	for particle in particles:	
		radius=particle[4]/2+h_layer
		radii.append(radius)
		Vol = 4/3*chrono.CH_PI*pow(radius, 3)
		# mass = rho/2*Vol
		totVol += Vol
		Vols.append(Vol)

	rho_new = rho*specimenVol/totVol

	print(f"Original rho: {rho}, new rho: {rho_new}")

	# permute particles if needed
	if permutation > 0.0:
		each_permutation = permutation / len(particles)
		# generate random normalized vector in 3D
		for i in range(len(particles)):
			random_dir_3D = numpy.random.rand(3)
			random_dir_3D /= numpy.linalg.norm(random_dir_3D)
			# move each particle along this direction by a distance of each_permutation
			particles[i][1] += random_dir_3D[0] * each_permutation 
			particles[i][2] += random_dir_3D[1] * each_permutation
			particles[i][3] += random_dir_3D[2] * each_permutation

	i = 0
	for particle in particles:
		x_pos=particle[1]
		y_pos=particle[2]
		z_pos=particle[3]
		radius=radii[i]
       
		pos = chrono.ChVector3d(x_pos, y_pos, z_pos)
		if len(file_name) == 1:
			pos=csys.TransformPointLocalToParent(pos)
		# print("pos: ", pos, "\t", "rotated pos: ", pos)
		mass = rho_new * Vols[i]

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

	out_dir = current_dir + "/OUT_slump_V1/"; 
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
	DEMSim.InstructBoxDomainDimension(1000.0, 1000.0, 1000.0)
	DEMSim.SetErrorOutAvgContacts(250)

	DEMSim.SetGravitationalAcceleration([0, -9810, 0])

	floor_mat = DEMSim.LoadMaterial({"ENm_mat":ENm, "ENa_mat":ENa, "mortar_layer_mat":mortar_layer, 
										"alpha_mat":alpha, "beta_mat":beta, "np_mat":np, "sgmTmax_mat":sgmTmax, 
										"sgmTau0_mat":sgmTau0, "kappa0_mat":kappa0, "eta_inf_mat":eta_inf, 
										"flocbeta_mat":flocbeta, "flocTcr_mat":flocTcr, "flocm_mat":flocm, "lambda_init_mat":lambda_init})

	# Floor is at y=0, pointing towards +y
	bottom_plane = DEMSim.AddBCPlane([0,0,0], [0,1,0], floor_mat)
	# bottom_plane.SetFamily(FIXED_FAM)

	cyl_radius=100
	cyl_height=150
	density=7.8E-9

	SlumpConeBody = AddSlumpCone(DEMSim, floor_mat, density*1000, cyl_radius, cyl_height, current_dir)
	# Let family 10 be fixed
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
	particle_tracker = ReadDFCparticles(DEMSim, [current_dir + "/slump_va_A076800.dat", current_dir + "/DFCgeo000-data-particles-SCC-2.dat"], rho, floor_mat, mortar_layer, permutation)
	# particle_tracker = ReadDFCparticles(DEMSim, [current_dir + "/DFCgeo000-data-particles-SCC-2.dat"], rho, floor_mat, mortar_layer)

	# Top plate is at y=300, pointing towards -y
	topPlate = DEMSim.AddBCPlane([0,300,0], [0,-1,0], floor_mat)
	topPlate.SetFamily(FIXED_FAM)

	KE_finder = DEMSim.CreateInspector("clump_kinetic_energy")
	# maxZ_finder = DEMSim.CreateInspector("clump_max_z")

	timestep = 1.0e-8
	DEMSim.SetInitTimeStep(timestep)
	DEMSim.SetMaxVelocity(400000.)
	DEMSim.SetErrorOutVelocity(500000.)

	# DEMSim.SetCDUpdateFreq(5)
	# DEMSim.DisableAdaptiveUpdateFreq()
	DEMSim.Initialize()

	# 1. phase: waiting with a plate on top of the SlumpCone
	print("\n 1. phase: waiting with a plate on top of the SlumpCone \n")

	max_step_size = 1e-5

	stepsize_increase_mult = 1.01
	currframe = 0
	t = 0.0
	
	step_num = 0

	frame_time = 1. / 20.
	out_time = 0.


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
			DEMSim.WriteMeshFile(meshname)
			print(f"Current number of contacts: {DEMSim.GetNumContacts()}")
			# print(f"Bin num: {DEMSim.GetBinNum()}")
			# print(f"Bin size: {DEMSim.GetBinSize()}")
			# DEMSim.ShowThreadCollaborationStats()
			# DEMSim.ShowMemStats()
			# DEMSim.PrintKinematicScratchSpaceUsage()

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
	timestep=1e-8
	stepsize_increase_mult = 1.01

	DEMSim.DoDynamicsThenSync(0)
	DEMSim.UpdateStepSize(timestep)
	
	# Disable contact between some entities and the particles
	DEMSim.DisableContactBetweenFamilies(FIXED_FAM, 0)

	start_time = time.time()
	while t < phase1_end + 60.:
		if (step_num % 10 == 0) and (timestep < max_step_size):
			timestep = min(timestep*stepsize_increase_mult,max_step_size)
			DEMSim.UpdateStepSize(timestep)
			# ke = KE_finder.GetValue()
			# print(f"Kinetic energy: {ke}")


		if (t >= out_time):
			print(f"============ Time {t}, slumping =============")
			print(f"Current step size: {timestep}")
			print(f"Outputting frame: {currframe}")
			filename = os.path.join(out_dir, f"slump_va_M{currframe:06d}.csv")
			meshname = os.path.join(out_dir, f"slump_mesh{currframe:06d}.vtk")
			DEMSim.WriteSphereFile(filename)
			DEMSim.WriteMeshFile(meshname)
			print(f"Current number of contacts: {DEMSim.GetNumContacts()}")
			# print(f"Bin num: {DEMSim.GetBinNum()}")
			# print(f"Bin size: {DEMSim.GetBinSize()}")
			# DEMSim.ShowThreadCollaborationStats()
			# DEMSim.ShowMemStats()
			# DEMSim.PrintKinematicScratchSpaceUsage()

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
