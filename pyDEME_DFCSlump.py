import os
import DEME
import pychrono as chrono

# Concrete properties
h_layer=4.0; 
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
ENm=2.0   # 4.0E-2;
ENa=100
h=mortar_layer	
alpha=0.25
beta=0.5
np=1.0
sgmTmax=2.0E-4
sgmTau0=1E-7
kappa0=100.
eta_inf=2.0E-5

def AddSlumpCone(sys, mat, density, radius, height, current_dir):
    obj_file = current_dir + "/coni_100_50_300.obj"
    SlumpCone = sys.AddWavefrontMeshObject(obj_file, mat)
    SlumpCone.SetMass(1e8)
    q = chrono.QuatFromAngleAxis(chrono.CH_PI_2, chrono.VECT_Z)
    SlumpCone.SetInitQuat([q.e1, q.e2, q.e3, q.e0])
    SlumpCone.SetInitPos([0,133,0])
    return sys.Track(SlumpCone)

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
				particles.append((n, x, y, z, d))

	q = chrono.QuatFromAngleAxis(-(chrono.CH_PI_2), chrono.VECT_X)
	# ChVector3d t(0,0,0)
	csys = chrono.ChCoordsysd(chrono.ChVector3d(0, 0, 0), q)
	
	my_templates = []
	body_points_new = []
	for particle in particles:
		x_pos=particle[1]
		y_pos=particle[2]
		z_pos=particle[3]
       
		pos = chrono.ChVector3d(x_pos, y_pos, z_pos)
		rotatedpos=csys.TransformPointLocalToParent(pos)
		print("pos: ", pos, "\t", "rotated pos: ", rotatedpos)
		
		radius=particle[4]/2+h_layer
		mass = rho/2*4/3*pow(radius, 3)*chrono.CH_PI

        # create templates for each sphere
		my_templates.append(sys.LoadSphereType(
            mass,
            radius,
            mat
        ))
		body_points_new.append([rotatedpos.x(), rotatedpos.y(), rotatedpos.z()])
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
		
    DEMSim = DEME.DEMSolver()
    DEMSim.SetVerbosity("INFO")
    DEMSim.SetOutputFormat("CSV")
    DEMSim.SetOutputContent(["ABSV", "XYZ"])
    DEMSim.SetMeshOutputFormat("VTK")
    # This is needed, as contact wildcards are extracted later on, so we need to notify the solver to include them
    DEMSim.SetContactOutputContent(["OWNER", "CNT_WILDCARD"])
    DEMSim.SetJitifyClumpTemplates(False)
    DEMSim.SetJitifyMassProperties(False)
    DEMSim.InstructBoxDomainDimension(1000.0, 1000.0, 1000.0)
    DEMSim.SetErrorOutAvgContacts(250)
	
    DEMSim.SetGravitationalAcceleration([0, -9810, 0])
	
    floor_mat = DEMSim.LoadMaterial({"ENm_mat":ENm, "ENa_mat":ENa, "mortar_layer_mat":mortar_layer, 
                                       "alpha_mat":alpha, "beta_mat":beta, "np_mat":np, "sgmTmax_mat":sgmTmax, 
                                       "sgmTau0_mat":sgmTau0, "kappa0_mat":kappa0, "eta_inf_mat":eta_inf})
	
    # Floor is at y=0, pointing towards +y
    bottom_plane = DEMSim.AddBCPlane([0,0,0], [0,1,0], floor_mat)
    
    cyl_radius=100
    cyl_height=150
    density=7.8E-9

    SlumpConeBody = AddSlumpCone(DEMSim, floor_mat, density*1000, cyl_radius, cyl_height, current_dir)

    # Change default contact force algorithm to fresh concrete ones
    my_force_model = DEMSim.ReadContactForceModel("DFCModel.mod")
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
    ReadDFCparticles(DEMSim, current_dir + "/DFCgeo000-data-particles-SCC-2/", rho, floor_mat, mortar_layer)
    
	# Top plate is at y=300, pointing towards -y
    topPlate = DEMSim.AddBCPlane([0,300,0], [0,-1,0], floor_mat)

	//
	////////////////////////////////////////////
    // Bind all existing visual shapes to the visualization system
    vis->AttachSystem(&sys);
	///////////////////////////////////////////////////////////////////////////////////////////////////////    	     
	//
    // SOLVER & STEPPER ---  Modify some setting of the physical system for the simulation
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////////// 
    // Modify some setting of the physical system for the simulation, if you want
	std::cout<<"\nSOLVER: Pardiso MKL - Stepper: HHT SetAlpha(-0.2)\n"<<std::endl;
	auto solver = chrono_types::make_shared<ChSolverPardisoMKL>();    
    sys.SetSolver(solver); 
    solver->UseSparsityPatternLearner(true);
    //solver->LockSparsityPattern(true);   
    solver->SetVerbose(false);
	//
	auto mystepper=chrono_types::make_shared<ChTimestepperHHT>(&sys);    
	mystepper->SetAlpha(-0.2); // alpha=-0.2 default value
	mystepper->SetMaxIters(50);
	mystepper->SetAbsTolerances(1e-06, 1e-06);  
	mystepper->SetVerbose(false);	    
	mystepper->SetModifiedNewton(true);
	mystepper->SetStepControl(false);
    sys.SetTimestepper(mystepper);

    // ===== Simulation loop =============
    int step_num=0;
    char filename[200];
    //
    double timestep = 1.0e-5;
    double current_time;
    //
    auto body_list = sys.GetBodies();
	//
	
    // 1. phase: remove particles above the prism height ////////////////////////////////////////////////
    //
    std::cout<<"\n 1. phase: waiting with a plate on top of the SlumpCone\n"<<std::endl;	
    /////////////////////////////////////////////////////////////////////////////////////////////////////
	//
    //
    int kk=0;
    //
	//
    timestep = 1.0e-8;
    while (  vis->Run() &  sys.GetChTime()<0.75) {
		//
		//vis->Run() & 
    	vis->BeginScene();
		vis->Render();  
		vis->EndScene();
		//	
		//
		if (step_num%20==0){	
			timestep=std::min(timestep*=1.05, 1.0e-5);
		}
		//
        sys.DoStepDynamics(timestep);
					
		if (std::fmod(step_num, 800) ==0) {
		    //
		    sprintf(filename, "%s/slump_va_A%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, true);
		    sprintf(filename, "%s/slump_va_M%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, false);		    
		    std::cout<< "filename: "<<filename<<"\t";
			//
			double rho_new=materialFCM->CalculateEffDensity(sys, specimenVol, rho);
			if (rho_new<rho/4.)
			    rho_new=rho/4.;
			materialFCM->ModifyDensity(sys, rho_new);
			//
			std::vector<std::shared_ptr<ChBody>> bodylist=sys.GetBodies();
			double aggVolFrac=calc_aggVolFrac(bodylist, mortar_layer, specimenVol);
			double current_time= sys.GetChTime();
			double IE;
			container->IterOnContactList(current_time,IE);
			double KE=calculateKE(sys);
			std::cout<<" time: "<< current_time << " va0: " << va0 <<  " aggVolFrac: "<< aggVolFrac << " IE: "<< IE << " KE: "<< KE << "\n";
			
			kk++;
		}

		step_num++;
		 
		}

	//
	// 2. phase: remove containers and allow flow  ////////////////////////////////////////////////////
	//
	std::cout<<"\n 2. phase: flow  \n"<<std::endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////
	SlumpConeBody->EnableCollision(false);
	topPlate->EnableCollision(false);
	//
	timestep=1e-8;
	//
	auto t1=sys.GetChTime();
	//
	while ( vis->Run() & sys.GetChTime()<t1+9.5) {
		//vis->Run() & 
		vis->BeginScene();
		vis->Render();  
		vis->EndScene(); 
		//
		sys.DoStepDynamics(timestep);
       	
		if (step_num%10)
		timestep=std::min(timestep*=1.01, 1.0e-5);
		// 
		if (std::fmod(step_num, 800) ==0) {
		    sprintf(filename, "%s/slump_va_A%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, true);
		    sprintf(filename, "%s/slump_va_M%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, false);		    
		    std::cout<< "filename: "<<filename<<"\t";
		    // 
		    //
		    double current_time= sys.GetChTime();
		    std::vector<std::shared_ptr<ChBody>> bodylist=sys.GetBodies();
		    double aggVolFrac=calc_aggVolFrac(bodylist, mortar_layer, specimenVol);
		    double IE;
		    container->IterOnContactList(current_time,IE);
		    double KE=calculateKE(sys);
		    std::cout<<" time: "<< current_time << " va0: " << va0 <<  " aggVolFrac: "<< aggVolFrac << " IE: "<< IE << " KE: "<< KE << "\n";
     
		 }
		 
		 step_num++;
    } 

    return 0;
}