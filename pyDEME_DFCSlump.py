import os
import DEME

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
    DEMSim.InstructBoxDomainDimension(500.0, 500.0, 500.0)
    DEMSim.SetErrorOutAvgContacts(250)
	
    DEMSim.SetGravitationalAcceleration([0, -9810, 0])
	
    # Floor is at y=0, pointing towards +y
	DEMSim.Add
	
    floor_mat = DEMSim.LoadMaterial({"ENm_mat":ENm, "ENa_mat":ENa, "mortar_layer_mat":mortar_layer, 
                                       "alpha_mat":alpha, "beta_mat":beta, "np_mat":np, "sgmTmax_mat":sgmTmax, 
                                       "sgmTau0_mat":sgmTau0, "kappa0_mat":kappa0, "eta_inf_mat":eta_inf,
                                       "flocbeta_mat":flocbeta, "flocTcr_mat":flocTcr, "flocm_mat":flocm, "lambda_init_mat":lambda_init})
    
    # Create the floor body
    auto floor_mat = chrono_types::make_shared<ChContactMaterialSMC>();
    floor_mat->SetYoungModulus(1e2); //2e3
    floor_mat->SetFriction(0.5);
    floor_mat->SetRestitution(0.0f);
    floor_mat->SetAdhesion(0.0);
    //floor_mat->SetCohesion(0);
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(900, 10, 900, 2.4e-9, true, true, floor_mat);
    floorBody->SetPos(ChVector3d(0, -5, 0));
    floorBody->SetFixed(true);
    sys.Add(floorBody);
    //	
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetYoungModulus(1e2); //2e3
    mat->SetFriction(0.5);
    mat->SetRestitution(0.0f);
    mat->SetAdhesion(0.0);
    //mat->SetCohesion(0);	
    double cyl_radius=100;
    double cyl_height=150;
    double density=7.8E-9;
    //
    auto SlumpConeBody =AddSlumpCone(sys, mat, density*1000, cyl_radius, cyl_height, current_dir);
    SlumpConeBody->GetVisualShape(0)->SetOpacity(0.4f);
    SlumpConeBody->SetFixed(true);
    //
    //
    //volume of specimen
    double specimenVol=( CH_PI*pow(100,2)*600/3- CH_PI*pow(50,2)*300/3);
    double Vol=specimenVol;
    //
    // Concrete properties
    double h_layer=4.0; 
    double minD=5.0;
    double maxD=10.0;
    double cement=571; 
    double WtoC=0.45;  
    double AtoC=2.62;  
    double rho_c=3150;
    double rho_w=1000;
    double vair=0.03;
    double nF=0.5;	
    //
    double va=1.0-cement/rho_c-(WtoC*cement)/rho_w-vair;
    double va0=(1-pow((minD/maxD),nF))*va;	
    double Va0=va0*Vol;
    double rho=cement*(1.0+WtoC+AtoC)*1E-12;
    //
    //
    // Define material parameters of fresh concrete
    //
    float mortar_layer=h_layer;
    float ENm=2.0;   //4.0E-2;
    float ENa=100;
    float h=mortar_layer;		
    float alpha=0.25;
    float beta=0.5;
    float np=1.0;
    float sgmTmax=2.0E-4;
    float sgmTau0=1E-7;
    float kappa0=100.;
    float eta_inf=2.0E-5;    
    //
    //
	//////////////////////////////////////////////////
    //	Change default contact force algorithm to fresh concrete ones
    //////////////////////////////////////////////////
    auto Contact_force_algorithm = chrono_types::make_unique<FCContactForce>();
    auto materialFCM=chrono_types::make_shared<ChMaterialFCM>(ENm, ENa, mortar_layer, 
									alpha, beta, np, sgmTmax, sgmTau0, kappa0, eta_inf);
    Contact_force_algorithm->Set_Material(materialFCM);
    sys.SetContactForceTorqueAlgorithm(std::move(Contact_force_algorithm));
    //
    auto container = chrono_types::make_shared<MyContactContainer>();
    sys.SetContactContainer(container);	
    //
    //////////////////////////////////////////////////
    // Read and load predefined particle list
    //////////////////////////////////////////////////
    std::string data_path=current_dir;
    std::string file_name="DFCgeo000-data-particles-SCC-2";
    ReadDFCparticles(sys, vis, data_path, file_name, rho, mortar_layer);
    // 
    ////////////////////////////////////////////////
    //
	auto topPlate = AddTopWall(sys);
	auto top_plate_collision = topPlate->GetCollisionModel()->GetShapeInstance(0).shape;
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