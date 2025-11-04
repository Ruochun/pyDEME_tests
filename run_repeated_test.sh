# Run the simulation with different parameters

for i in {1..10}
do
    python pyDEME_DFCSlump.py --out_name "slump_va_final_${i}.csv" --permutation 0.0
    cp OUT_slump_V1/slump_va_final_${i}.csv slump_va_final_${i}.csv
done


for i in {1..10}
do
    python pyDEME_DFCSlump.py --out_name "slump_va_final_permu_${i}.csv" --permutation 20000.0
    cp OUT_slump_V1/slump_va_final_permu_${i}.csv slump_va_final_permu_${i}.csv
done
