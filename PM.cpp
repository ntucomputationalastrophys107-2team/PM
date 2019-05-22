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
const int Scheme_PS = 1;        // scheme of poisson solver ( 1=FFT )
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
    if( BC==1 ){             // Periodic Boundary Condition
       if( Scheme_PS == 1){  // FFT
       
       double KK;     // K*K = Kx*Kx + Ky*Ky + Kz*Kz
       int nx, ny;
       fftw_complex rho_K[N][N][N/2+1];  // density   in K space
       fftw_complex phi_K[N][N][N/2+1];  // potential in K space
       fftw_plan rhoXtorhoK, phiKtophiX;
       rhoXtorhoK = fftw_plan_dft_r2c_3d( N, N, N, &rho[0][0][0], &rho_K[0][0][0], FFTW_ESTIMATE ); // Fourier Transform from rho(x) to rho(k)
       phiKtophiX = fftw_plan_dft_c2r_3d( N, N, N, &phi_K[0][0][0], &phi[0][0][0], FFTW_ESTIMATE ); // Inverse Fourier Transform from phi(k) to phi(x)

       fftw_execute( rhoXtorhoK ); // Fourier Transform from rho(x) to rho(k)
       for(int i=0;i<N;i++){
       for(int j=0;j<N;j++){
       for(int k=0;k<N/2+1;k++){
           if( i>N/2) nx=i-N; else nx=i;
           if( j>N/2) ny=j-N; else ny=j;
           KK = (nx*2*M_PI/L)*(nx*2*M_PI/L) + (ny*2*M_PI/L)*(ny*2*M_PI/L) + (k*2*M_PI/L)*(k*2*M_PI/L); // K*K = Kx*Kx + Ky*Ky + Kz*Kz
           phi_K[i][j][k][0] = rho_K[i][j][k][0]/KK*4.0*M_PI*G*(-1.0/(N*N*N)); // real part phi(k) = -4*Pi*G*rho(k)/k^2 and normalize by 1/N^3
           phi_K[i][j][k][1] = rho_K[i][j][k][1]/KK*4.0*M_PI*G*(-1.0/(N*N*N)); // imag part phi(k) = -4*Pi*G*rho(k)/k^2 and normalize by 1/N^3
       }}}
       phi_K[0][0][0][0] = 0.0;       // set DC to 0
       phi_K[0][0][0][1] = 0.0;       // set DC to 0
       fftw_execute( phiKtophiX ); // Inverse Fourier Transform from phi(k) to phi(x)

       fftw_destroy_plan( rhoXtorhoK );
       fftw_destroy_plan( phiKtophiX );

       }
       else printf("ERROR: Unsupported Scheme_PS!\n");
    }
    else if( BC==2 ){        // Isolated Boundary Condition
       if( Scheme_PS == 1){  // FFT

       double mas_0pad[2*N][2*N][2*N]; // zero padding on the mass array
       double greensfn[2*N][2*N][2*N]; // symmetric discrete Green's function -1/r
       double phi_0pad[2*N][2*N][2*N]; // zero padding on the potential array
       int nx,ny,nz;                   // number of distance interval
       for(int i=0;i<2*N;i++){
       for(int j=0;j<2*N;j++){
       for(int k=0;k<2*N;k++){
           if( i<N && j<N && k<N ) mas_0pad[i][j][k] = rho[i][j][k]*dx*dx*dx;  // mass of cell = density * cell volume
           else mas_0pad[i][j][k] = 0.0;                                       // zero padding

           if( i>=N ) nx=2*N-i; else nx=i;  // symmetrization 
           if( j>=N ) ny=2*N-j; else ny=j;  // symmetrization
           if( k>=N ) nz=2*N-k; else nz=k;  // symmetrization

           if( i==0 && j==0 && k==0) greensfn[i][j][k] = 0.0;  // ignore self force
           else greensfn[i][j][k] = -1.0/(dx*sqrt( nx*nx + ny*ny + nz*nz )); // -1/r
       }}}
       fftw_complex mas_0pad_K[2*N][2*N][N+1];  // mass in K space
       fftw_complex greensfn_K[2*N][2*N][N+1];  // Green's function in K space
       fftw_complex phi_0pad_K[2*N][2*N][N+1];  // potential in K space
       fftw_plan masXtomasK, greXtogreK, phiKtophiX;
       masXtomasK = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, &mas_0pad[0][0][0], &mas_0pad_K[0][0][0], FFTW_ESTIMATE ); // Fourier Transform from mas(x) to mas(k)
       greXtogreK = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, &greensfn[0][0][0], &greensfn_K[0][0][0], FFTW_ESTIMATE ); // Fourier Transform from gre(x) to gre(k)
       phiKtophiX = fftw_plan_dft_c2r_3d( 2*N, 2*N, 2*N, &phi_0pad_K[0][0][0], &phi_0pad[0][0][0], FFTW_ESTIMATE ); // Inverse Fourier Transform from phi(k) to phi(x)
          
       fftw_execute( masXtomasK ); // Fourier Transform from mas(x) to mas(k)
       fftw_execute( greXtogreK ); // Fourier Transform from gre(x) to gre(k)
       for(int i=0;i<2*N;i++){
       for(int j=0;j<2*N;j++){
       for(int k=0;k<N+1;k++){
           phi_0pad_K[i][j][k][0] = G*( mas_0pad_K[i][j][k][0]*greensfn_K[i][j][k][0]-mas_0pad_K[i][j][k][1]*greensfn_K[i][j][k][1])*(1.0/(2*N*2*N*2*N));// real part phi(k) = G*mas(k)*gre(k) and normalize by 1/(2N)^3
           phi_0pad_K[i][j][k][1] = G*( mas_0pad_K[i][j][k][0]*greensfn_K[i][j][k][1]+mas_0pad_K[i][j][k][1]*greensfn_K[i][j][k][0])*(1.0/(2*N*2*N*2*N));// imag part phi(k) = G*mas(k)*gre(k) and normalize by 1/(2N)^3
       }}}
       fftw_execute( phiKtophiX ); // Inverse Fourier Transform from phi(k) to phi(x)

       fftw_destroy_plan( masXtomasK );
       fftw_destroy_plan( greXtogreK );
       fftw_destroy_plan( phiKtophiX );

       for(int i=0;i<N;i++)
       for(int j=0;j<N;j++)
       for(int k=0;k<N;k++)
           phi[i][j][k] = phi_0pad[i][j][k]; // remove the padding
       
       }
       else printf("ERROR: Unsupported Scheme_PS!\n");
    }
    else printf("ERROR: Unsuppoted Boundary Condition!\n");

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
