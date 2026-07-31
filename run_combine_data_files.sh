for i in {1..10}
do
    # convert 1 to 10 to A to K
    num2letter=$(echo "$i" | awk '{printf "%c", $1+64}')
    # echo "Combining data files for simulation $num2letter"
    ln -sf ./Slump_init_cond/Slump-DFC-${num2letter}_sim01/DFC-${num2letter}_slump_va_A076800.vtk ./Slump_init_cond/
    python combine_dat_files.py ./Slump_init_cond/Slump-DFC-${num2letter}_sim01/DFC-${num2letter}_slump_va_A076800_coords.dat \
        ./Slump_init_cond/Slump-DFC-${num2letter}_sim01/DFC-${num2letter}_slump_va_A076800_radius.dat ./Slump_init_cond/init_${i}.dat
done




#
