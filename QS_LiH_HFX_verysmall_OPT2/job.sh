#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --cpus-per-task=64
#SBATCH -A pc2-mitarbeiter
#SBATCH -t 30:00

export OMP_NUM_THREADS=64
source ../env.sh
rm *.wfn

nx=3
ny=3
nz=3

sed "s/_NREPX_/$nx/g" input_bulk_B88_3.inp | sed "s/_NREPY_/$ny/g" | sed "s/_NREPZ_/$nz/g" > input_bulk_B88_${nx}_${ny}_${nz}.inp
sed "s/_NREPX_/$nx/g" input_bulk_HFX_3.inp | sed "s/_NREPY_/$ny/g" | sed "s/_NREPZ_/$nz/g" > input_bulk_HFX_${nx}_${ny}_${nz}.inp

srun ../exe/local/cp2k.psmp input_bulk_B88_${nx}_${ny}_${nz}.inp > input_bulk_B88_${nx}_${ny}_${nz}.out 2> input_bulk_B88_${nx}_${ny}_${nz}.err
cp LiH_bulk_3-RESTART.wfn B88.wfn
srun ../exe/local/cp2k.psmp input_bulk_HFX_${nx}_${ny}_${nz}.inp > input_bulk_HFX_${nx}_${ny}_${nz}.out 2> input_bulk_HFX_${nx}_${ny}_${nz}.err
