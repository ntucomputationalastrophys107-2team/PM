#include <cstdio>
#include <cstdlib>
#include <cmath>

//-----------------------------------
// A Partilce Mesh Code
// ----------------------------------
/*  Parallelizaion to be finished */

// constants
const double L = 1.0;           // length of the 3-D domain box
const int    N = 16;            // number of grid in each direction
const double dx = L/N;          // spatial resolution
const double dt = 0.1;          // time step
const int    ParN  = 2;         // number of particles
const double G = 1.0;           // gravitational constant
const double end_time = 10.0;   // end time of the evolution
const double ParM[ParN]={ 1.0, 1.0};  // mass of each particle

// schemes
const int BC = 1;               // boundary condition ( 1=periodic, 2=isolated )
const int Scheme_MD = 1;        // scheme of mass deposition ( 1=NGP, 2=CIC, 3=TSC )
const int Scheme_PS = 1;        // scheme of poisson solver ( 1=FFT, 2=SOR )
const int Scheme_OI = 1;        // shceme of orbit integration ( 1=KDK, 2=DKD, 3=RK4 )


// FUNCTION Init: Set the initial condition
void Init( double x[ParN][3], double v[ParN][3] ){
    /* To be modified for the test problem */
    x[0][0] = 0.25*L;
    x[0][1] = 0.5*L;
    x[0][2] = 0.5*L;
    v[0][0] = 0.0;
    v[0][1] = 0.0;
    v[0][2] = 1.0;
    x[1][0] = 0.75*L;
    x[1][1] = 0.5*L;
    x[1][2] = 0.5*L;
    v[1][0] = 0.0;
    v[1][1] = 0.0;
    v[1][2] = -1.0;
    return;
}// FUNCTION Init


// FUNCTION MassDeposition: Deposit particle mass onto grids
void MassDeposition( double x[ParN][3], double rho[N][N][N] ){
    /* To be finished */

    return;
}// FUNCTION MassDeposition


// FUNCTION PoissonSolver: Solve the Poisson equation to get the potential
void PoissonSolver( double rho[N][N][N], double phi[N][N][N] ){
    /* To be finished */

    return;
}// FUNCTION PoissonSolver


// FUNCTION Update: Update the system by dt
void Update( double x[ParN][3], double v[ParN][3], double phi[N][N][N] ){
    /* To be finished */

    return;
}// FUNCTION Update


// FUNCTION Energy: Calculate the total energy of the system
double Energy( double x[ParN][3], double v[ParN][3] ){
    double e;
    /* To be finished */
    e = 1.0;

    return e;
}//FUNCTION Energy


// FUNCTION Momentum: Calculate the total momentum of the system
double Momentum( double v[ParN][3] ){
    double p;
    /* To be finished */
    p = 1.0;

    return p;
}// FUNCTION Momentum


// FUNCTION CheckBoundary: Deal with the particles outside the box
void CheckBoundary( double x[ParN][3], double v[ParN][3] ){
    /* To be finished */

    return;
}// FUNCTION CheckBoundary


// FUNCTION main
int main( int argc, char *argv[] ){
    // variables
    double t = 0.0;           // time
    double ParX[ParN][3];     // position of the particle
    double ParV[ParN][3];     // velocity of the particle
    double Rho[N][N][N]={0};  // density
    double Phi[N][N][N]={0};  // potential

    double E0;                // initial energy
    double E;                 // energy
    double Eerr;              // energy conservation error
    double P0;                // initial momentum |P0|
    double P;                 // moentum |P|
    double Perr;              // momentum conservation error

    // initialization
    Init( ParX, ParV );
    E0 = Energy( ParX, ParV );
    P0 = Momentum( ParV );

    // main loop
    printf(" Start!\n");
    while( t<end_time ){
        printf("Time: %f -> %f, dt=%f\n", t, t+dt, dt );

        MassDeposition( ParX, Rho );
        PoissonSolver( Rho, Phi );
        Update( ParX, ParV, Phi );
        CheckBoundary( ParX, ParV );

        E = Energy( ParX, ParV );
        Eerr = (E-E0)/E0;
        P = Momentum( ParV );
        Perr = (P-P0)/P0;
        t = t+dt;
        /* Visualizaion() to be finished */
    }
    printf(" End!\n\n" );

    // output the results
    printf("Energy   Error: %10.3e\n", Eerr );
    printf("Momentum Error: %10.3e\n", Perr );

    return EXIT_SUCCESS;
}// FUNCTION main
