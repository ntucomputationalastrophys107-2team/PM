# PM
Final Project-IV: Particle Mesh


Compilation:

g++ -o PM PM.cpp -L$(FFTW_PATH)/lib -I$(FFTW_PATH)/include -lfftw3 -lrt -fopenmp

Execution:

./PM > log

Plot:

python plot.py
