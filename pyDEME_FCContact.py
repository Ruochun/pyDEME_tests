# =============================================================================
#
# =============================================================================

import DEME
from DEME import HCPSampler

import numpy
import os

math_PI = 3.1415927

# For clarity, we attach the force model file here.
# This model can be found in file ForceModelWithElectrostatic.cu, too.

import pychrono as chrono
# import pychrono.irrlicht as chronoirr


def force_model():
    model = '''
    do {
      // Get material properties
      //
      float mortar_layer = mortar_layer_mat[bodyAMatType];
      float ENm = ENm_mat[bodyAMatType];
      float ENa = ENa_mat[bodyAMatType];
      float alpha = alpha_mat[bodyAMatType];
      float beta = beta_mat[bodyAMatType];
      float np = np_mat[bodyAMatType];
      float sgmTmax = sgmTmax_mat[bodyAMatType];
      float sgmTau0 = sgmTau0_mat[bodyAMatType];
      float kappa0 = kappa0_mat[bodyAMatType];
      float eta_inf = eta_inf_mat[bodyAMatType];
      float flocbeta = flocbeta_mat[bodyAMatType];
      float flocm = flocm_mat[bodyAMatType];
      float flocTcr = flocTcr_mat[bodyAMatType];

      //
      float h = mortar_layer;
      float R1 = 5000.0;
      float R2 = 5000.0;
      double Rmin, Rmax, delta_rmin;
      double Lij, L0ij;
      // Material parameters of viscous part which should be defined inside
      // contactinfo
      double eta0 = kappa0 * eta_inf;
      double Deps0 = sgmTau0 / eta0;

      //
      // Get state variables relevant to this contact
      //
      float3 statevar;
      statevar.x = contact_info_strain_x;
      statevar.y = contact_info_strain_y;
      statevar.z = contact_info_strain_z;
      float lambda0 = contact_info_lambda;
      //
      // ==== SIGMA-t influenced by Flocculation ======
      sgmTmax = sgmTmax * (1. + lambda0);
      //
      if (lambda0 == 0.) {
        lambda0 = lambda_init_mat[bodyAMatType];
        float val0 = ParticleThixoAvg[AOwner];
        float val1 = ParticleThixoAvg[BOwner];
        //
        lambda0 = DEME_MIN(val0, val1);
        contact_info_lambda = lambda0;
        // std::cout << "lambda-4 " << map_contact_info[mykey].lambda <<  "key "
        // << mykey << std::endl;
        contact_info_step_time = time;
      }
      //
      // initialize force values
      //
      force = make_float3(0, 0, 0);

      if (overlapDepth < 0) {
        break;
      }

      float3 Vx, Vy, Vz;
      float3 normal_dir = -B2A; // DEME uses B2A as normal direction, opposite
                                // of Chrono as it seems

      XdirToDxDyDz(normal_dir, make_float3(0, 1, 0), Vx, Vy, Vz);
      // std::cout<<" Vy: "<<Vy<<" Vz: "<<Vz<<"\t";

      R1 = (ARadius < 50.) ? ARadius : 50.;
      R2 = (BRadius < 50.) ? BRadius : 50.;

      double delta = overlapDepth;

      if (ContactType == deme::SPHERE_SPHERE_CONTACT) {
        /// Concrete - Concrete Interaction
        // Contact of two DFC spheres
        // center to center distance between 2 objects
        L0ij = R1 + R2 - h;
        Lij = abs(R1 + R2 - delta);
        Rmin = DEME_MIN(R1, R2);
        Rmax = DEME_MAX(R1, R2);
        delta_rmin = delta / 2 * (2.0 * Rmax - delta) / Lij;

      } else if (ContactType ==
                 deme::SPHERE_MESH_CONTACT) { /// COncrete-container
                                              /// interaction
        h = mortar_layer / 2.;
        // sgmTmax=this->material->Get_sgmTmax()/2;
        L0ij = R1 + R2 - h;
        Lij = abs(R1 + R2 - delta);
        Rmin = R1;
        // Rmin=std::min(R1, R2);
        Rmax = DEME_MAX(R1, R2);
        delta_rmin = delta;

      } else { /// Fiber-container or another object- container interaction
        h = mortar_layer / 2.;
        // sgmTmax=this->material->Get_sgmTmax()/2;
        L0ij = R1 + R2 - h;
        Lij = abs(R1 + R2 - delta);
        Rmin = DEME_MIN(R1, R2);
        Rmax = DEME_MAX(R1, R2);
        delta_rmin = delta / 2.;
      }

      double ai = abs(Rmin - delta_rmin);
      double radius2 = Rmin * Rmin - ai * ai;
      double contact_area = deme::PI * radius2;

      //
      // modify penetration according to initial overlap
      //
      double delta_new = delta - h;
      //
      //
      // Relative velocities at contact at midplane (p0=(p1+p2)/2)
      //
      // In DEME, contact point in global is just called contactPnt

      // In DEME, local contact point is called locCPA and locCPB
      // DEME uses its own way to calculate rel vel
      float3 rotVelCPA, rotVelCPB;
      {
        // We also need the relative velocity between A and B in global frame to
        // use in the damping terms To get that, we need contact points'
        // rotational velocity in GLOBAL frame This is local rotational velocity
        // (the portion of linear vel contributed by rotation)
        rotVelCPA = cross(ARotVel, locCPA);
        rotVelCPB = cross(BRotVel, locCPB);
        // This is mapping from local rotational velocity to global
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPA.x, rotVelCPA.y,
                                                rotVelCPA.z, AOriQ.w, AOriQ.x,
                                                AOriQ.y, AOriQ.z);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPB.x, rotVelCPB.y,
                                                rotVelCPB.z, BOriQ.w, BOriQ.x,
                                                BOriQ.y, BOriQ.z);
      }
      float time_scale_val = 1.;
      float3 relvel = (BLinVel + rotVelCPB) - (ALinVel + rotVelCPA); // A2B vel
      relvel = relvel / time_scale_val;

      float relvel_n_mag = dot(relvel, normal_dir);
      float3 relvel_n = relvel_n_mag * normal_dir;
      float3 relvel_t = relvel - relvel_n;
      float relvel_t_mag = length(relvel_t);
      //
      // Calculate displacement increment in normal and tangential direction
      //
      float dT = ts;
      float delta_t = relvel_t_mag * dT;
      float delta_n = relvel_n_mag * dT;
      float3 v_delta_t = relvel_t * dT;
      //
      //  Calculate the strain increment in each local direction
      //
      float depsN = delta_n / Lij;
      float delta_M = dot(relvel, Vy) * dT; // v_delta_t.Dot(Vy);
      float delta_L = dot(relvel, Vz) * dT; // v_delta_t.Dot(Vz);
      float depsM = delta_M / Lij;
      float depsL = delta_L / Lij;
      //
      //
      //
      float epsA = log(1 - h / (L0ij));
      float epsN = log(Lij / L0ij);


      
      // float epsN=statevar.x+depsN;
      float epsM = statevar.y + depsM;
      float epsL = statevar.z + depsL;
      //
      float epsT = powf((epsM * epsM + epsL * epsL), 0.5);
      float epsQ = powf(epsN * epsN + alpha * epsT * epsT, 0.5);
      //
      //
      //
      statevar.x = epsN;
      statevar.y = epsM;
      statevar.z = epsL;
      contact_info_strain_x = statevar.x;
      contact_info_strain_y = statevar.y;
      contact_info_strain_z = statevar.z;
      contact_info_step_time = time;
      //
      //
      //
      if (epsN < 0) {
        ////////////////////////////////////////////////////////////
        // Compressive contact;
        ////////////////////////////////////////////////////////////
        float sgmN = 0;
        float sgmM = 0;
        float sgmL = 0;
        float sgmT = 0;
        //

        float stot;
        if (epsN >= epsA) {
          // stot=epsQ*ENm;
          sgmN = epsN * ENm;
          contact_info_strain_y = 0;
          contact_info_strain_z = 0;
        } else {
          sgmN = (epsA)*ENm + (epsN - epsA) * ENa;
          sgmM = 0;
          sgmL = 0;
        }

        // std::cout<<"delta_new: "<<delta_new<<" epsN: "<<epsN<< " epsQ:
        // "<<epsQ<<" stot: "<<stot<<" sgmN: "<<sgmN<<std::endl;
        //////////////////////////////////////////////////////////
        //
        // Viscous stresses
        //
        //////////////////////////////////////////////////////////
        //
        // depsN=relvel_n_mag/Lij*dT;
        // depsT=relvel_t_mag/Lij*dT;
        // calculate strain rates
        float vdepsN = depsN / dT;
        float vdepsM = depsM / dT;
        float vdepsL = depsL / dT;
        float vDeps = powf(
            (beta * vdepsN * vdepsN + vdepsM * vdepsM + vdepsL * vdepsL), 0.5);
        //
        // std::cout << "Lij : " << Lij << " vdeps_Comp : " << vDeps <<
        // std::endl;
        float eta;
        float lambda = lambda0;

        {
          // double lambda = rungeKutta4<double>(dfloc, lambda0, current_time,
          // flocbeta, vDeps, flocm, flocTcr, dT, 1); lambda=0;
          // map_contact_info[mykey].lambda=lambda;
          if (vDeps <= Deps0) {
            eta = eta0;
          } else {
            eta = eta_inf * powf(abs(vDeps), np - 1.) +
                  sgmTau0 * (1.0 + lambda) / vDeps;
          }
        }

        float sgmN_vis = beta * eta * vdepsN;
        float sgmM_vis = eta * vdepsM;
        float sgmL_vis = eta * vdepsL;
        float sgmT_vis = powf(sgmM_vis * sgmM_vis + sgmL_vis * sgmL_vis, 0.5);
        //
        //
        //
        float forceN = contact_area * (sgmN + sgmN_vis);
        float forceM = contact_area * (sgmM + sgmM_vis);
        float forceL = contact_area * (sgmL + sgmL_vis);
        // std::cout<<"epsN "<<epsN<<" sgmN "<<sgmN<<std::endl;
        // force[0]=forceN;force[1]=contact_area*(sgmM-sgmM_vis);force[2]=contact_area*(sgmL-sgmL_vis);
        force = forceN * normal_dir + forceM * Vy + forceL * Vz;
        // std::cout<<" force "<<force<<std::endl;
        /*
        force = -forceN * normal_dir;
        if (relvel_t_mag >= sys.GetSlipVelocityThreshold())
                force -= (forceT / relvel_t_mag) * relvel_t;
*/
        // return -force;

      } else {
        ////////////////////////////////////////////////////////////
        // Tensile contact
        ////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////
        //
        // Calculate Stress from material stiffness
        //
        //////////////////////////////////////////////////////////
        float sgmN = ENm * epsN;
        float sgmM = 0;
        float sgmL = 0;
        if (sgmN > sgmTmax)
          sgmN = sgmTmax;
        //
        //////////////////////////////////////////////////////////
        //
        // Viscous stresses
        //
        //////////////////////////////////////////////////////////
        //
        // calculate strain rates
        float vdepsN = depsN / dT;
        float vdepsM = depsM / dT;
        float vdepsL = depsL / dT;
        float vDeps = powf(
            (beta * vdepsN * vdepsN + vdepsM * vdepsM + vdepsL * vdepsL), 0.5);
        //
        // std::cout << "Lij : " << Lij << " vdeps_Tension : " << vDeps <<
        // std::endl;
        float eta;
        float lambda = lambda0;

        {
          // double lambda = rungeKutta4<double>(dfloc, lambda0, current_time,
          // flocbeta, vDeps, flocm, flocTcr, dT, 1); lambda=0;
          // map_contact_info[mykey].lambda=lambda;
          if (vDeps <= Deps0) {
            eta = eta0;
          } else {
            eta = eta_inf * powf(abs(vDeps), np - 1.) +
                  sgmTau0 * (1.0 + lambda) / vDeps;
          }
        }
        // std::cout<<"relvel"<<relvel<<" Deps0: "<<Deps0<<" vDeps"<<vDeps<<"
        // eta: "<<eta<<std::endl;
        float sgmN_vis = beta * eta * vdepsN;
        float sgmM_vis = eta * vdepsM;
        float sgmL_vis = eta * vdepsL;
        // double sgmT=powf(sgmM*sgmM+sgmL*sgmL,0.5);
        //////////////////////////////////////////////////////////
        //
        // Combine Viscous stresses and stiffness stresses and calculate forces
        //
        //////////////////////////////////////////////////////////
        // exit(0);
        float forceN = contact_area * (sgmN + sgmN_vis);
        float forceM = contact_area * (sgmM + sgmM_vis);
        float forceL = contact_area * (sgmL + sgmL_vis);
        // std::cout<<"epsN "<<epsN<<" sgmN "<<sgmN<<std::endl;
        // double forceT = sgmT * contact_area;
        // force[0]=forceN;force[1]=contact_area*(sgmM-sgmM_vis);force[2]=contact_area*(sgmL-sgmL_vis);
        force = forceN * normal_dir + forceM * Vy + forceL * Vz;
        /*
        force = -forceN * normal_dir;
        if (relvel_t_mag >= sys.GetSlipVelocityThreshold())
                force -= (forceT / relvel_t_mag) * relvel_t;
        */
        // return -force;
      }
      //   ChVector3d torque = Vcross((p1 - p0), -force);
      // return std::make_pair(-force, torque);
      //   return {-force, torque};

      if (isnan(length(force))) {
        printf("Force NaN!!!!\\n");
      }
    } while (0);
    '''

    return model


density_steel = 7.8e-9


def AddBMCContainer(DEMSim, mat, density, obj_file):
    bmcContainer = DEMSim.AddWavefrontMeshObject(obj_file, mat)
    q = chrono.QuatFromAngleAxis(chrono.CH_PI_2, chrono.VECT_Z)
    # bmcContainer.SetInitQuat([q.e1, q.e2, q.e3, q.e0])
    # MOve the container down a bit
    bmcContainer.InformCentroidPrincipal([0, 36, 0], [0, 0, 0, 1])
    bmcContainer.SetFamily(2)
    # Fix family 2
    DEMSim.SetFamilyFixed(2)
    DEMSim.DisableContactBetweenFamilies(1, 2)
    tracked = DEMSim.Track(bmcContainer)
    return tracked


def AddBMCLid(DEMSim, mat, density, obj_file):
    lid = DEMSim.AddWavefrontMeshObject(obj_file, mat)
    lid.SetInitPos([0, 60, 0])
    lid.SetFamily(4)
    DEMSim.SetFamilyFixed(4)
    DEMSim.DisableContactBetweenFamilies(1, 4)
    DEMSim.DisableFamilyOutput(4)
    tracked = DEMSim.Track(lid)
    return tracked


def CreateBMCrod(DEMSim, mat, obj_file):
    # trimesh = chrono.ChTriangleMeshConnected()
    # trimesh.LoadWavefrontMesh(obj_file, False, True)

    # chrono.ChInertiaUtils.PrincipalInertia(minertia, mcog, mmass, minertia)
    bmcRod = DEMSim.AddWavefrontMeshObject(obj_file, mat)
    bmcRod.InformCentroidPrincipal([0, 0, 0], [0, 0, 0, 1])
    bmcRod.SetMass(3.42506e-05)
    bmcRod.SetMOI([0.00296771, 0.0931488, 0.09585])
    bmcRod.SetInitPos([0, -35, 0])
    # Family 3 is not rotational, 10 is fixed
    bmcRod.SetFamily(10)
    DEMSim.DisableContactBetweenFamilies(1, 3)

    # For now, give this family a prescribed motion
    DEMSim.SetFamilyPrescribedAngVel(3, "0", "0.5", "0")
    DEMSim.SetFamilyFixed(10)
    tracked = DEMSim.Track(bmcRod)
    return tracked

def AddFallingItems(DEMSim, mat, coords_file, rad_file, h_layer, rho):
    body_points = []
    body_radius = []

    with open(coords_file, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                x, y, z = map(float, line.split())
                body_points.append([x, y, z])

    body_points = numpy.array(body_points)
    # print(body_points)

    with open(rad_file, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                r = float(line)
                body_radius.append(r)

    body_radius = numpy.array(body_radius)

    numSpheres = len(body_radius)
    body_points_new = [] # numpy.empty((numSpheres, 3), dtype=float)
    my_templates = []
    for i in range(numSpheres):
        x_pos = body_points[i, 0]
        y_pos = body_points[i, 1]
        z_pos = body_points[i, 2]
        y_level = y_pos-0
        if (y_pos>60 or y_pos<-27) or (x_pos*x_pos+z_pos*z_pos>30.*30.):
            continue
        body_points_new.append([x_pos, y_pos, z_pos])
        radius = body_radius[i]
        mass = rho/2*4/3*(radius**3)*chrono.CH_PI

        # create templates for each sphere
        my_templates.append(DEMSim.LoadSphereType(
            mass,
            radius,
            mat
        ))
        # print(f"Adding sphere {i} at {pos} with radius {radius} and mass {mass}")

    particles = DEMSim.AddClumps(my_templates, body_points_new)
    return DEMSim.Track(particles)


if __name__ == "__main__":
    out_dir = "DemoOutput_FCContact/"
    out_dir = os.path.join(os.getcwd(), out_dir)
    os.makedirs(out_dir, exist_ok=True)

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

    # # E, nu, CoR, mu, Crr...
    # floor_mat = DEMSim.LoadMaterial(
    #     {"E": 1e2, "nu": 0.3, "CoR": 0.0, "mu": 0.0, "Crr": 0.0}
    # )

    specimenVol = 40.0 * 40.0 * 250.0
    Vol = specimenVol
    h_layer = 2.5
    minD = 2
    maxD = 4
    cement = 845
    WtoC = 0.21
    AtoC = 1.24
    rho_c = 3150.0
    rho_w = 1000.0
    vair = 0.03
    nF = 0.5
    va = 1.0 - cement / rho_c - (WtoC * cement) / rho_w - vair
    va0 = (1.0 - pow((minD / maxD), nF)) * va
    Va0 = va0 * Vol
    rho = cement * (1.0 + WtoC + AtoC) * 1e-12
    targetmass = Va0 * rho
    targetVol = specimenVol * 2.5

    mortar_layer=h_layer
    ENm=1.5E-3
    ENa=100.0*0.001
    h=mortar_layer
    alpha=0.25
    beta=0.5
    np=1.0
    sgmTmax=3.0E-4
    sgmTau0=2.0E-5
    kappa0=100.
    eta_inf=1.25E-5
    flocbeta = 0.01
    flocm = 1.0
    flocTcr = 100.
    lambda_init = 2.0

    materialFCM = DEMSim.LoadMaterial({"ENm_mat":ENm, "ENa_mat":ENa, "mortar_layer_mat":mortar_layer, 
                                       "alpha_mat":alpha, "beta_mat":beta, "np_mat":np, "sgmTmax_mat":sgmTmax, 
                                       "sgmTau0_mat":sgmTau0, "kappa0_mat":kappa0, "eta_inf_mat":eta_inf,
                                       "flocbeta_mat":flocbeta, "flocTcr_mat":flocTcr, "flocm_mat":flocm, "lambda_init_mat":lambda_init})

    density = 7.8e-9
    bmcContainer = AddBMCContainer(DEMSim, materialFCM, density, "bmcContainer.obj")
    bmcLid = AddBMCLid(DEMSim, materialFCM, density, "bmc_lid.obj")
    bmcRod = CreateBMCrod(DEMSim,materialFCM, "bmcRod.obj")

    particles = AddFallingItems(DEMSim, materialFCM, "bmc_particle_New4_coords.dat", "bmc_particle_New4_radius.dat", h_layer, rho*0.7051*0.966)

    
    # We can specify the force model using a string.
    # This force model is the standard Hertzian--Mindlin model, plus an electrostatic force, based on something
    # called a `geometry wildcard'. This is an extra property that we can associate with each geometry entity,
    # such as triangle and sphere. We use this value to derive the electrostatic force. But first, we need to
    # declare in the force model that a geometry wildcard is in use...
    # my_force_model = DEMSim.DefineContactForceModel(force_model())
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
    
    DEMSim.SetForceCalcThreadsPerBlock(256)

    # Some inspectors
    KE_finder = DEMSim.CreateInspector("clump_kinetic_energy")

    timestep = 1e-8
    DEMSim.SetInitTimeStep(timestep)
    DEMSim.SetGravitationalAcceleration([0, -9810, 0])
    DEMSim.SetMaxVelocity(5000.)
    DEMSim.SetErrorOutVelocity(5000.)

    # DEMSim.SetCDUpdateFreq(5)
    # DEMSim.DisableAdaptiveUpdateFreq()
    DEMSim.Initialize()

    num_owners = DEMSim.GetNumOwners()
    print(f"Numer of all entities: {num_owners}")
    particles_IDs = particles.GetOwnerIDs()
    # print(f"Particles IDs: {particles_IDs}")
    bmcContainer_ID = bmcContainer.GetOwnerID()
    bmcLid_ID = bmcLid.GetOwnerID()
    bmcRod_ID = bmcRod.GetOwnerID()
    # print(f"Container ID: {bmcContainer_ID}, Lid ID: {bmcLid_ID}, Rod ID: {bmcRod_ID}")
    # All entities are these arrays concatenated, which is simply 0 to max_number_of_entities - 1.
    IDs = particles_IDs + [bmcContainer_ID, bmcLid_ID, bmcRod_ID]
    assert IDs[-1] == num_owners - 1, "The number of owners does not match the IDs."
    ParticleThixoInfo = numpy.zeros(num_owners, dtype=float)

    max_step_size = 5e-6

    stepsize_increase_mult = 1.0001
    currframe = 0
    t = 0.0
    step_num = 0

    frame_time = 1. / 100.
    out_time = 0.
    while t < 0.5:
        if (step_num % 10 == 0) and (timestep < max_step_size):
            timestep = min(timestep*stepsize_increase_mult,max_step_size)
            DEMSim.UpdateStepSize(timestep)

        if (t >= out_time):
            """
            # we do not need to get ParticleThixoAvg, as it's cleared every time anyway
            # ParticleThixoAvg = DEMSim.GetAllOwnerWildcardValue("ParticleThixoAvg")
            ParticleThixoAvg = numpy.zeros(num_owners, dtype=float)
            ParticleThixoCount = numpy.zeros(num_owners, dtype=int)

            # We can use a tiny number to discourage output of non-force-producing (aka preemptive) contacts.
            contact_info = DEMSim.GetContactDetailedInfo(1e-20)
            # Contact pairs are already sorted based on ID 1
            contact_type = contact_info.GetContactType()
            As = contact_info.GetAOwner()
            Bs = contact_info.GetBOwner()
            lambda_values = contact_info.Get("contact_info_lambda")
            for i in range(len(contact_type)):
                lambda_value = lambda_values[i]
                bodyA = As[i]
                bodyB = Bs[i]
                ParticleThixoAvg[bodyA] += lambda_value
                ParticleThixoAvg[bodyB] += lambda_value
                ParticleThixoCount[bodyA] += 1
                ParticleThixoCount[bodyB] += 1
            # Now we can average it, using numpy's vectorized operations
            ParticleThixoAvg = numpy.divide(
                ParticleThixoAvg, ParticleThixoCount, out=numpy.zeros_like(ParticleThixoAvg), where=ParticleThixoCount != 0
            )
            # print(f"ParticleThixoAvg: {ParticleThixoAvg}")   
            # Feed the updated values to the solver  
            DEMSim.SetOwnerWildcardValue(0, "ParticleThixoAvg", ParticleThixoAvg)
            """

            print(f"============ Time {t} =============")
            print(f"Outputting frame: {currframe}")
            filename = os.path.join(out_dir, f"DEMdemo_output_{currframe:04d}.csv")
            meshname = os.path.join(out_dir, f"DEMdemo_mesh_{currframe:04d}.vtk")
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
            if (ke < 1e-3) and (step_num > 4000):
                print("Kinetic energy is low, finishing up the settling phase.")
                break
            
        
        DEMSim.DoDynamics(timestep)
        t += timestep
        step_num += 1

    DEMSim.ShowTimingStats()
    DEMSim.ShowThreadCollaborationStats()
    DEMSim.ShowMemStats()

    # now rotate the rod
    DEMSim.DoDynamicsThenSync(0)
    DEMSim.ChangeFamily(10, 3)  # Change family of the rod to 3, which is not fixed

    print(f"============================================")
    print(f"============ Settling Finished =============")
    print(f"============================================")


    # t = 0.0
    # timestep = 1e-8
    # DEMSim.UpdateStepSize(timestep)

    out_time = t
    while t < 60.0:
        if (step_num % 20 == 0) and (timestep < max_step_size):
            timestep = min(timestep*stepsize_increase_mult,max_step_size)
            DEMSim.UpdateStepSize(timestep)

        if (t >= out_time):
            """
            # we do not need to get ParticleThixoAvg, as it's cleared every time anyway
            # ParticleThixoAvg = DEMSim.GetAllOwnerWildcardValue("ParticleThixoAvg")
            ParticleThixoAvg = numpy.zeros(num_owners, dtype=float)
            ParticleThixoCount = numpy.zeros(num_owners, dtype=int)

            # We can use a tiny number to discourage output of non-force-producing (aka preemptive) contacts.
            contact_info = DEMSim.GetContactDetailedInfo(1e-20)
            # Contact pairs are already sorted based on ID 1
            contact_type = contact_info.GetContactType()
            As = contact_info.GetAOwner()
            Bs = contact_info.GetBOwner()
            lambda_values = contact_info.Get("contact_info_lambda")
            for i in range(len(contact_type)):
                lambda_value = lambda_values[i]
                bodyA = As[i]
                bodyB = Bs[i]
                ParticleThixoAvg[bodyA] += lambda_value
                ParticleThixoAvg[bodyB] += lambda_value
                ParticleThixoCount[bodyA] += 1
                ParticleThixoCount[bodyB] += 1
            # Now we can average it, using numpy's vectorized operations
            ParticleThixoAvg = numpy.divide(
                ParticleThixoAvg, ParticleThixoCount, out=numpy.zeros_like(ParticleThixoAvg), where=ParticleThixoCount != 0
            ) 
            # Feed the updated values to the solver  
            DEMSim.SetOwnerWildcardValue(0, "ParticleThixoAvg", ParticleThixoAvg)
            """
                
            print(f"============ Time {t} =============")
            print(f"Outputting frame: {currframe}")
            filename = os.path.join(out_dir, f"DEMdemo_output_{currframe:04d}.csv")
            meshname = os.path.join(out_dir, f"DEMdemo_mesh_{currframe:04d}.vtk")
            DEMSim.WriteSphereFile(filename)
            DEMSim.WriteMeshFile(meshname)
            print(f"Current number of contacts: {DEMSim.GetNumContacts()}")
            # DEMSim.ShowThreadCollaborationStats()
            
            currframe += 1
            out_time += frame_time

        DEMSim.DoDynamics(timestep)
        t += timestep
        step_num += 1

    DEMSim.ShowTimingStats()
    DEMSim.ShowThreadCollaborationStats()
    DEMSim.ShowMemStats()

    print("FCContact demo exiting...")


