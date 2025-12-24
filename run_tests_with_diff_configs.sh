# Run the simulation with different parameters

python3 pyDEME_DFCSlump.py --out_name "slump_va_final_A.csv" -s --particle_files "DFC_gen/DFC-A_slump_va_A076800.vtk" "DFC_gen/DFCgeo000-data-particles-SCC-A.dat"
cp OUT_slump_V1/particle_trajectories.csv particle_trajectories_A.csv

python3 pyDEME_DFCSlump.py --out_name "slump_va_final_B.csv" -s --particle_files "DFC_gen/DFC-B_slump_va_A076800.vtk" "DFC_gen/DFCgeo000-data-particles-SCC-B.dat"
cp OUT_slump_V1/particle_trajectories.csv particle_trajectories_B.csv

python3 pyDEME_DFCSlump.py --out_name "slump_va_final_C.csv" -s --particle_files "DFC_gen/DFC-C_slump_va_A076800.vtk" "DFC_gen/DFCgeo000-data-particles-SCC-C.dat"
cp OUT_slump_V1/particle_trajectories.csv particle_trajectories_C.csv

python3 pyDEME_DFCSlump.py --out_name "slump_va_final_D.csv" -s --particle_files "DFC_gen/DFC-D_slump_va_A076800.vtk" "DFC_gen/DFCgeo000-data-particles-SCC-D.dat"
cp OUT_slump_V1/particle_trajectories.csv particle_trajectories_D.csv

python3 pyDEME_DFCSlump.py --out_name "slump_va_final_E.csv" -s --particle_files "DFC_gen/DFC-E_slump_va_A076800.vtk" "DFC_gen/DFCgeo000-data-particles-SCC-E.dat"
cp OUT_slump_V1/particle_trajectories.csv particle_trajectories_E.csv
