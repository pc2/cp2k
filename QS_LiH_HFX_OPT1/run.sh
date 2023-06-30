source env.sh
mpirun -np 16 --output-filename out ../exe/local/cp2k.psmp input_bulk_HFX_3.inp > input_bulk_HFX_3.out
