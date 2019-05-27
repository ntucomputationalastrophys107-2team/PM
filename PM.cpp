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


// FUNCTION Acceleration: Calculate the acceleration of each particle
void Acceleration( double x[ParN][3], double a[ParN][3] ){
    double Rho[N][N][N];  // density
    double Phi[N][N][N];  // potential
    MassDeposition( x, Rho );
    PoissonSolver( Rho, Phi );
    /* To be finished */

    return;
}// FUNCTION Acceleration


// FUNCTION Update: Update the system by dt
void Update( double x[ParN][3], double v[ParN][3] ){
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

// FUNCTION OutputData: Output the results into files
void OutputData( double t, double x[ParN][3], double e, double eerror, double p, double perror ){
    FILE *fp;
    fp = fopen( "Data_ParticlePosition", "a" );

    if(t==0.0){
        fprintf( fp, "#        Time");
        for(int i=0;i<ParN;i++) fprintf( fp,"            x_%d            y_%d            z_%d",i+1,i+1,i+1);
        fprintf( fp, "\n");
    }

    fprintf( fp, "%13.6e", t );
    for(int i=0;i<ParN;i++) fprintf( fp, "  %13.6e  %13.6e  %13.6e",x[i][0],x[i][1],x[i][2] );
    fprintf( fp, "\n");

    fclose(fp);

    FILE *fc;
    fc = fopen( "Data_ConservedQuantity", "a" );

    if(t==0.0)
        fprintf( fc, "#        Time         Energy    EnergyError       Momentum  MomentumError\n");
    fprintf( fc, "%13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n", t, e, eerror, p, perror );

    fclose(fc);

    return;
}// FUNCTION OutputData

// FUNCTION main
int main( int argc, char *argv[] ){
    // variables
    double t = 0.0;           // time
    double ParX[ParN][3];     // position of the particle
    double ParV[ParN][3];     // velocity of the particle

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
    OutputData( t, ParX, E0, 0.0, P0, 0.0 );

    // main loop
    printf(" Start!\n");
    while( t<end_time ){
        printf("Time: %13.6e -> %13.6e, dt=%f\n", t, t+dt, dt );

        Update( ParX, ParV );
        CheckBoundary( ParX, ParV );

        E = Energy( ParX, ParV );
        Eerr = (E-E0)/fabs(E0);
        P = Momentum( ParV );
        Perr = (P-P0)/fabs(P0);
        t = t+dt;
        OutputData( t, ParX, E, Eerr, P, Perr );
    }
    printf(" End!\n\n" );

    // output the results
    printf("Info:\n");
    printf("Number of Particles = %d\n", ParN);
    printf("Number of Grids     = %dx%dx%d\n",N,N,N);
    switch( BC ){
        case 1:
            printf("Boundary Condition  = Periodic\n");
            break;
        case 2:
            printf("Boundary Condition  = Isolated\n");
            break;
    }
    switch( Scheme_MD ){
        case 1:
            printf("Mass Deposition     = NGP\n");
            break;
        case 2:
            printf("Mass Deposition     = CIC\n");
            break;
        case 3:
            printf("Mass Deposition     = TSC\n");
            break;
    }
    switch( Scheme_PS ){
        case 1:
            printf("Poisson Solver      = FFT\n");
            break;
    }
    switch( Scheme_OI ){
        case 1:
            printf("Orbit Integration   = KDK\n");
            break;
        case 2:
            printf("Orbit Integration   = DKD\n");
            break;
        case 3:
            printf("Orbit Integration   = Rk4\n");
            break;
    }
    printf("\n");
    printf("Energy   Error: %13.6e\n", Eerr );
    printf("Momentum Error: %13.6e\n", Perr );

    return EXIT_SUCCESS;
}// FUNCTION main
