#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>

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
    if(BC==1){
       
       double Ksquare[N][N][N/2+1];
       int ni, nj;
       for(int i=0;i<N;i++)
          for(int j=0;j<N;j++)
             for(int k=0;k<N/2+1;k++){
                if( i>N/2) ni=i-N; else ni=i;
                if( j>N/2) nj=j-N; else nj=j;
                Ksquare[i][j][k] = (ni*2*M_PI/L)*(ni*2*M_PI/L) + (nj*2*M_PI/L)*(nj*2*M_PI/L) + (k*2*M_PI/L)*(k*2*M_PI/L);
             }
       Ksquare[0][0][0] = 1.0;
    
       fftw_complex out[N][N][N/2+1];
       fftw_plan p, q;
       p = fftw_plan_dft_r2c_3d( N, N, N, &rho[0][0][0], &out[0][0][0], FFTW_ESTIMATE );
       q = fftw_plan_dft_c2r_3d( N, N, N, &out[0][0][0], &phi[0][0][0], FFTW_ESTIMATE );

       fftw_execute(p);
       for(int i=0;i<N;i++)
          for(int j=0;j<N;j++)
             for(int k=0;k<N/2+1;k++){
                out[i][j][k][0] = out[i][j][k][0]/Ksquare[i][j][k]*4.0*M_PI*G*(-1.0/(N*N*N));
                out[i][j][k][1] = out[i][j][k][1]/Ksquare[i][j][k]*4.0*M_PI*G*(-1.0/(N*N*N));
             }
       fftw_execute(q);

       fftw_destroy_plan(p);
       fftw_destroy_plan(q);
    }
    else if(BC==2){
       if(Scheme_PS==1){
          double rho_zeropadding[2*N][2*N][2*N];
          double Greens_R[2*N][2*N][2*N];
          double phi_zeropadding[2*N][2*N][2*N];

          for(int i=0;i<N;i++)
             for(int j=0;j<N;j++)
                for(int k=0;k<N;k++){
                   rho_zeropadding[i][j][k] = rho[i][j][k];
                   if( i==0 && j==0 && k==0) Greens_R[i][j][k] = 0.0;
                   else Greens_R[i][j][k] = -1.0/sqrt(i*i+j*j+k*k);
                }

          for(int i=N;i<2*N;i++)
             for(int j=N;j<2*N;j++)
                for(int k=N;k<2*N;k++){
                   rho_zeropadding[i][j][k] = 0.0;
                   Greens_R[i][j][k] = -1.0/sqrt((2*N-i)*(2*N-i)+(2*N-j)*(2*N-j)+(2*N-k)*(2*N-k));
                }

       
          fftw_complex rho_zeropadding_K[2*N][2*N][N+1];
          fftw_complex Greens_R_K[2*N][2*N][N+1];
          fftw_complex out[2*N][2*N][N+1];
          fftw_plan p_rho, p_greens, q;
          p_rho = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, &rho_zeropadding[0][0][0], &rho_zeropadding_K[0][0][0], FFTW_ESTIMATE );
          p_greens = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, &Greens_R[0][0][0], &Greens_R_K[0][0][0], FFTW_ESTIMATE );
          q = fftw_plan_dft_c2r_3d( N, N, N, &out[0][0][0], &phi_zeropadding[0][0][0], FFTW_ESTIMATE );
          
          fftw_execute(p_rho);
          fftw_execute(p_greens);
          for(int i=0;i<2*N;i++)
             for(int j=0;j<2*N;j++)
                for(int k=0;k<N+1;k++){
                   out[i][j][k][0] = ( rho_zeropadding_K[i][j][k][0]*Greens_R_K[i][j][k][0]-rho_zeropadding_K[i][j][k][1]*Greens_R_K[i][j][k][1])*(G*dx*dx/(2*N*2*N*2*N));
                   out[i][j][k][1] = ( rho_zeropadding_K[i][j][k][0]*Greens_R_K[i][j][k][1]+rho_zeropadding_K[i][j][k][1]*Greens_R_K[i][j][k][0])*(G*dx*dx/(2*N*2*N*2*N));
                }
          fftw_execute(q);

          fftw_destroy_plan(p_rho);
          fftw_destroy_plan(p_greens);
          fftw_destroy_plan(q);

          for(int i=0;i<N;i++)
             for(int j=0;j<N;j++)
                for(int k=0;k<N;k++)
                   phi[i][j][k] = phi_zeropadding[i][j][k];


       }
       else if(Scheme_PS==2){
          ;
       }
       else
          printf("ERROR: Unsupported Scheme_PS!\n");
    }
    else
       printf("ERROR: Unsuppoted Boundary Condition!\n");


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

    // main loop
    printf(" Start!\n");
    while( t<end_time ){
        printf("Time: %f -> %f, dt=%f\n", t, t+dt, dt );

        Update( ParX, ParV );
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
