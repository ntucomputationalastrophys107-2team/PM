# PM
Final Project-IV: Particle Mesh


Compilation:

mpicxx -o PM PM.cpp -L$(FFTW_PATH)/lib -I$(FFTW_PATH)/include -lfftw3 -lrt -fopenmp

Execution:

mpirun -np 2 PM > log

Plot:

python plot.py
