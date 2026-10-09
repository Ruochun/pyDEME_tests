
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_1.dat -d ./1-4 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_1.dat -d ./1-5 -s

# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_2.dat -d ./2-2 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_2.dat -d ./2-3 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_2.dat -d ./2-4 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_2.dat -d ./2-5 -s

# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_3.dat -d ./3-2 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_3.dat -d ./3-3 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_3.dat -d ./3-4 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_3.dat -d ./3-5 -s

# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_4.dat -d ./4-2 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_4.dat -d ./4-3 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_4.dat -d ./4-4 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_4.dat -d ./4-5 -s

# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_5.dat -d ./5-2 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_5.dat -d ./5-3 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_5.dat -d ./5-4 -s
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_5.dat -d ./5-5 -s

# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_3.dat -d ./3-9 -s > ./time3-9.txt 2>&1
# python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_3.dat -d ./3-10 -s > ./time3-10.txt 2>&1

for i in {1..10}
do
    for j in {1..10}
    do
        rm -rf ./${i}-${j}
        mkdir -p ./${i}-${j}
        python pyDEME_DFCSlump.py -p ./Slump_init_cond/init_${i}.dat -d ./${i}-${j} -s > ./time${i}-${j}.txt 2>&1
    done
done 

# for i in {1..5}
# do
#     for j in {1..10}
#     do
#         rm -rf ./${i}-${j}
#         mkdir -p ./${i}-${j}
#         python pyDEME_DFCLBox.py -p ./Slump_init_cond/init_${i}.dat -d ./${i}-${j} -s > ./time${i}-${j}.txt 2>&1
#     done
# done 



# run csv_to_xyz for all the output files, in folders (1..5)-(1..5)
# for i in {1..5}
# do
#     for j in {1..10}
#     do
#         python csv_to_xyz.py ./DFC_Slump_different_init_runs/${i}-${j}/slump_va_M000599.csv ./DFC_Slump_different_init_runs/rs${i}-${j}.xyz
#     done
# done    

# # reconstruct surface out of these xyz files
# for i in {1..5}
# do
#     for j in {1..10}
#     do
#         echo "Processing run ${i}-${j}"
#         python surface_construction.py ./DFC_Slump_different_init_runs/rs${i}-${j}.xyz ./DFC_Slump_different_init_runs/rs${i}-${j}-res 70. 0 1 1
#     done
# done    