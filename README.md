# PM
Final Project-IV: Particle Mesh


Compilation:

(CPU)
mpicxx -o PM PM.cpp -L$(FFTW_PATH)/lib -I$(FFTW_PATH)/include -lfftw3 -lrt -fopenmp

(GPU)
nvcc -arch sm_13 -o PM PM.cu -L$(OPENMPI_PATH)/lib64 -I$(OPEN_MPI)/include -lmpi -lrt -Xcompiler -fopenmp

Execution:

mpirun -np 2 ./PM > log

Plot:

python plot.py
