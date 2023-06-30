ntasks=3
export OMP_NUM_THREADS=1
#source ../env.sh
source  /home/rschade/tmp/eri/cp2k/tools/toolchain/install/setup
#rm *.wfn

nx=1
ny=1
nz=1

sed "s/_NREPX_/$nx/g" input_bulk_B88_3.inp | sed "s/_NREPY_/$ny/g" | sed "s/_NREPZ_/$nz/g" > input_bulk_B88_${nx}_${ny}_${nz}.inp
sed "s/_NREPX_/$nx/g" input_bulk_HFX_3.inp | sed "s/_NREPY_/$ny/g" | sed "s/_NREPZ_/$nz/g" > input_bulk_HFX_${nx}_${ny}_${nz}.inp

#mpirun -n $ntasks ../exe/local/cp2k.psmp input_bulk_B88_${nx}_${ny}_${nz}.inp > input_bulk_B88_${nx}_${ny}_${nz}.out 2> input_bulk_B88_${nx}_${ny}_${nz}.err
#cp LiH_bulk_3-RESTART.wfn B88.wfn
mpirun -n $ntasks -output-filename  out ../exe/local/cp2k.psmp input_bulk_HFX_${nx}_${ny}_${nz}.inp > input_bulk_HFX_${nx}_${ny}_${nz}.out 2> input_bulk_HFX_${nx}_${ny}_${nz}.err
