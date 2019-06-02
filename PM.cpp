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
const double dt = 0.01;          // time step
const int    ParN  = 2;         // number of particles
const double G = 1.0;           // gravitational constant
const double end_time = 10.0;   // end time of the evolution
const double ParM[ParN]={ 1.0, 1.0};  // mass of each particle

// schemes
const int BC = 2;               // boundary condition ( 1=periodic, 2=isolated )
const int Scheme_MD = 1;        // scheme of mass deposition ( 1=NGP, 2=CIC, 3=TSC )
const int Scheme_PS = 1;        // scheme of poisson solver ( 1=FFT )
const int Scheme_OI = 1;        // shceme of orbit integration ( 1=KDK, 2=DKD, 3=RK4 )
void CheckBoundary( double x[ParN][3], double v[ParN][3] );


// FUNCTION Init: Set the initial condition
void Init( double x[ParN][3], double v[ParN][3] ){
    /* To be modified for the test problem */
    x[0][0] = 0.25*L;
    x[0][1] = 0.5*L;
    x[0][2] = 0.5*L;
    v[0][0] = 0.0;
    v[0][1] = 0.01;
    v[0][2] = 1.0;
    x[1][0] = 0.75*L;
    x[1][1] = 0.5*L;
    x[1][2] = 0.5*L;
    v[1][0] = 0.0;
    v[1][1] = 0.01;
    v[1][2] = -1.0;

    return;
}// FUNCTION Init


// FUNCTION MassDeposition: Deposit particle mass onto grids
void MassDeposition( double x[ParN][3], double rho[N][N][N] ){
    for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
    for(int k=0;k<N;k++)
        rho[i][j][k]=0.0;   // initialization as zero

    if(Scheme_MD==1){ //NGP
        for(int n=0;n<ParN;n++){
            if( x[n][0]>0 && x[n][0]<L && x[n][1]>0 && x[n][1]<L && x[n][2]>0 && x[n][2]<L )
            rho[int(x[n][0]/dx)][int(x[n][1]/dx)][int(x[n][2]/dx)] += ParM[n]/(dx*dx*dx);
        } 
    }
    if(Scheme_MD==2){  //CIC
        int index_x[8];
        int index_y[8];
        int index_z[8];
        double weighting[3][2];
        for (int n=0;n<ParN;n++){
            if( x[n][0]>0.5*dx && x[n][0]<L-0.5*dx && x[n][1]>0.5*dx && x[n][1]<L-0.5*dx && x[n][2]>0.5*dx && x[n][2]<L-0.5*dx ){
            for(int i=0;i<8;i++){
                if(i<4){
                    index_x[i]=int(x[n][0]/dx);
                    index_y[i]=int(x[n][1]/dx);
                    index_z[i]=int(x[n][2]/dx);
                }
                if(i>=4){
                    index_x[i]=int(x[n][0]/dx)+1;
                    index_y[i]=int(x[n][1]/dx)+1;
                    index_z[i]=int(x[n][2]/dx)+1;
                }
            }
            for(int dim=0;dim<3;dim++){
                if(x[n][dim]-int(x[n][dim]/dx)*dx<=dx*0.5){
                    weighting[dim][0]=double((int(x[n][dim]/dx)+0.5)*dx-x[n][dim]);
                    weighting[dim][1]=double(-(int(x[n][dim]/dx)-0.5)*dx+x[n][dim]);
                }
                if(x[n][dim]-int(x[n][dim]/dx)*dx>dx*0.5){
                    weighting[dim][0]=double((int(x[n][dim]/dx)+1.5)*dx-x[n][dim]);
                    weighting[dim][1]=double(-(int(x[n][dim]/dx)+0.5)*dx+x[n][dim]);
                }
            }
            for(int i=0;i<8;i++){
                rho[index_x[i]][index_y[i]][index_z[i]]+=weighting[0][i%2]*weighting[1][i%2]*weighting[2][i%2]*ParM[n]/pow(dx,3);
            }
            }
        }            
    }
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
    double acc[N+1][N+1][N+1][3];
    if (BC == 1){
        for ( int i = 0; i < N; i = i + 1 ){
            for ( int j = 0; j < N; j = j + 1 ){
                for (int k = 0; k < N; k = k + 1){
                    acc[i][j][k][0] = (0.25/dx)*(Phi[i][j][k]+Phi[i][(j-1)%N][k]+Phi[i][j][(k-1)%N]+Phi[i][(j-1)%N][(k-1)%N]-Phi[(i-1)%N][j][k]-Phi[(i-1)%N][(j-1)%N][k]-Phi[(i-1)%N][j][(k-1)%N]-Phi[(i-1)%N][(j-1)%N][(k-1)%N]);
                    acc[i][j][k][1] = (0.25/dx)*(Phi[i][j][k]+Phi[(i-1)%N][j][k]+Phi[i][j][(k-1)%N]+Phi[(i-1)%N][j][(k-1)%N]-Phi[i][(j-1)%N][k]-Phi[(i-1)%N][(j-1)%N][k]-Phi[i][(j-1)%N][(k-1)%N]-Phi[(i-1)%N][(j-1)%N][(k-1)%N]);
                    acc[i][j][k][2] = (0.25/dx)*(Phi[i][j][k]+Phi[(i-1)%N][j][k]+Phi[i][(j-1)%N][k]+Phi[(i-1)%N][(j-1)%N][k]-Phi[i][j][(k-1)%N]-Phi[(i-1)%N][j][(k-1)%N]-Phi[i][(j-1)%N][(k-1)%N]-Phi[(i-1)%N][(j-1)%N][(k-1)%N]);
                }
            }
        }
        for ( int j = 0; j < N; j = j + 1 ){
            for ( int k = 0; k < N; k = k + 1 ){
                for (int d = 0; d < 3; d = d + 1){
                    acc[N][j][k][d] = acc[1][j][k][d];
                }
            }
        }
        for ( int i = 0; i < N+1; i = i + 1 ){
            for ( int k = 0; k < N; k = k + 1 ){
                for (int d = 0; d < 3; d = d + 1){
                    acc[i][N][k][d] = acc[i][0][k][d];
                }
            }
        }
        for ( int i = 0; i < N+1; i = i + 1 ){
            for ( int j = 0; j < N+1; j = j + 1 ){
                for (int d = 0; d < 3; d = d + 1){
                    acc[i][j][N][d] = acc[i][j][0][d];
                }
            }
        }
    }////acc for periodic BC
    if (BC == 2){
        for ( int i = 1; i < N; i = i + 1 ){
            for ( int j = 1; j < N; j = j + 1 ){
                for (int k = 1; k < N; k = k + 1){
                    acc[i][j][k][0] = (0.25/dx)*(Phi[i][j][k]+Phi[i][j-1][k]+Phi[i][j][k-1]+Phi[i][j-1][k-1]-Phi[i-1][j][k]-Phi[i-1][j-1][k]-Phi[i-1][j][k-1]-Phi[i-1][j-1][k-1]);
                    acc[i][j][k][1] = (0.25/dx)*(Phi[i][j][k]+Phi[i-1][j][k]+Phi[i][j][k-1]+Phi[i-1][j][k-1]-Phi[i][j-1][k]-Phi[i-1][j-1][k]-Phi[i][j-1][k-1]-Phi[i-1][j-1][k-1]);
                    acc[i][j][k][2] = (0.25/dx)*(Phi[i][j][k]+Phi[i-1][j][k]+Phi[i][j-1][k]+Phi[i-1][j-1][k]-Phi[i][j][k-1]-Phi[i-1][j][k-1]-Phi[i][j-1][k-1]-Phi[i-1][j-1][k-1]);
                }
            }
        }
        for ( int i = 1; i < N; i = i + 1 ){
            for ( int j = 1; j < N; j = j + 1 ){
                for ( int d = 0; d < 3; d = d + 1 ){
                    acc[0][i][j][d]   = 2*acc[1][i][j][d]-acc[2][i][j][d];
                    acc[N+1][i][j][d] = 2*acc[N][i][j][d]-acc[N-1][i][j][d];
                    acc[i][j][0][d]   = 2*acc[i][j][1][d]-acc[i][j][2][d];
                    acc[i][j][N+1][d] = 2*acc[i][j][N][d]-acc[i][j][N-1][d];
                    acc[i][0][j][d]   = 2*acc[i][1][j][d]-acc[i][2][j][d];
                    acc[i][N+1][j][d] = 2*acc[i][N][j][d]-acc[i][N-1][j][d];
                }
            }
        }
        for ( int i = 1; i < N; i = i + 1 ){
            for ( int d = 0; d < 3; d = d + 1 ){
                acc[i][0][0][d]      = 0.5*(2*acc[i][1][0][d]  -acc[i][2][0][d]    +2*acc[i][0][1][d]  -acc[i][0][2][d]);
                acc[i][N+1][N+1][d]  = 0.5*(2*acc[i][N][N+1][d]-acc[i][N-1][N+1][d]+2*acc[i][N+1][N][d]-acc[i][N+1][N-1][d]);
                acc[0][i][0][d]      = 0.5*(2*acc[1][i][0][d]  -acc[2][i][0][d]    +2*acc[0][i][1][d]  -acc[0][i][2][d]);
                acc[N+1][i][N+1][d]  = 0.5*(2*acc[N][i][N+1][d]-acc[N-1][i][N+1][d]+2*acc[N+1][i][N][d]-acc[N+1][i][N-1][d]);
                acc[0][0][i][d]      = 0.5*(2*acc[1][0][i][d]  -acc[2][0][i][d]    +2*acc[0][1][i][d]  -acc[0][2][i][d]);
                acc[N+1][N+1][i][d]  = 0.5*(2*acc[N][N+1][i][d]-acc[N-1][N+1][i][d]+2*acc[N+1][N][i][d]-acc[N+1][N-1][i][d]);
            }
        }
        for ( int d = 0; d < 3; d = d + 1 ){
            acc[0][0][0][d]       = (2*acc[1][0][0][d]    -acc[2][0][0][d]      +2*acc[0][1][0][d]    -acc[0][2][0][d]      +2*acc[0][0][1][d]    -acc[0][0][2][d]      )/3;
            acc[0][0][N+1][d]     = (2*acc[1][0][N+1][d]  -acc[2][0][N+1][d]    +2*acc[0][1][N+1][d]  -acc[0][2][N+1][d]    +2*acc[0][0][N][d]    -acc[0][0][N-1][d]    )/3;
            acc[0][N+1][0][d]     = (2*acc[1][N+1][0][d]  -acc[2][N+1][0][d]    +2*acc[0][N][0][d]    -acc[0][N-1][0][d]    +2*acc[0][N+1][1][d]  -acc[0][N+1][2][d]    )/3;
            acc[0][N+1][N+1][d]   = (2*acc[1][N+1][N+1][d]-acc[2][N+1][N+1][d]  +2*acc[0][N][N+1][d]  -acc[0][N-1][N+1][d]  +2*acc[0][N+1][N][d]  -acc[0][N+1][N-1][d]  )/3;
            acc[N+1][0][0][d]     = (2*acc[N][0][0][d]    -acc[N-1][0][0][d]    +2*acc[N+1][1][0][d]  -acc[N+1][2][0][d]    +2*acc[N+1][0][1][d]  -acc[N+1][0][2][d]    )/3;
            acc[N+1][0][N+1][d]   = (2*acc[N][0][N+1][d]  -acc[N-1][0][N+1][d]  +2*acc[N+1][1][N+1][d]-acc[N+1][2][N+1][d]  +2*acc[N+1][0][N][d]  -acc[N+1][0][N-1][d]  )/3;
            acc[N+1][N+1][0][d]   = (2*acc[N][N+1][0][d]  -acc[N-1][N+1][0][d]  +2*acc[N+1][N][0][d]  -acc[N+1][N-1][0][d]  +2*acc[N+1][N+1][1][d]-acc[N+1][N+1][2][d]  )/3;
            acc[N+1][N+1][N+1][d] = (2*acc[N][N+1][N+1][d]-acc[N-1][N+1][N+1][d]+2*acc[N+1][N][N+1][d]-acc[N+1][N-1][N+1][d]+2*acc[N+1][N+1][N][d]-acc[N+1][N+1][N-1][d])/3;
        }
        for ( int i = 1; i < N; i = i + 1 ){
            for ( int d = 0; d < 3; d = d + 1 ){
                acc[i][0][0][d]      = 0.5*(2*acc[i][1][0][d]  -acc[i][2][0][d]    +2*acc[i][0][1][d]  -acc[i][0][2][d]);
                acc[i][N+1][N+1][d]  = 0.5*(2*acc[i][N][N+1][d]-acc[i][N-1][N+1][d]+2*acc[i][N+1][N][d]-acc[i][N+1][N-1][d]);
                acc[0][i][0][d]      = 0.5*(2*acc[1][i][0][d]  -acc[2][i][0][d]    +2*acc[0][i][1][d]  -acc[0][i][2][d]);
                acc[N+1][i][N+1][d]  = 0.5*(2*acc[N][i][N+1][d]-acc[N-1][i][N+1][d]+2*acc[N+1][i][N][d]-acc[N+1][i][N-1][d]);
                acc[0][0][i][d]      = 0.5*(2*acc[1][0][i][d]  -acc[2][0][i][d]    +2*acc[0][1][i][d]  -acc[0][2][i][d]);
                acc[N+1][N+1][i][d]  = 0.5*(2*acc[N][N+1][i][d]-acc[N-1][N+1][i][d]+2*acc[N+1][N][i][d]-acc[N+1][N-1][i][d]);
            }
        }
    }////acc for isolated BC

    for ( int i = 0; i < ParN; i = i + 1 ){
        int ParIndex[3];
        double ParFrax[3];
        if (x[i][0] > L || x[i][0] < 0 || x[i][1] > L || x[i][1] < 0 || x[i][2] > L || x[i][2] < 0 ){
            for ( int d = 0; d < 3; d = d+1){
                a[i][d] = 0;
            }
        }
        else{
            for ( int d = 0; d < 3; d = d+1){
                ParIndex[d] = int(x[i][d]/dx);
                ParFrax[d] = x[i][d]/dx - int(x[i][d]/dx);
            }
            a[i][0] = 0.25*(1-ParFrax[0])*(acc[ParIndex[0]  ][ParIndex[1]+1][ParIndex[2]+1][0]+acc[ParIndex[0]  ][ParIndex[1]+1][ParIndex[2]  ][0]+acc[ParIndex[0]  ][ParIndex[1]  ][ParIndex[2]+1][0]+acc[ParIndex[0]  ][ParIndex[1]  ][ParIndex[2]  ][0])
                     +0.25*    ParFrax[0]*(acc[ParIndex[0]+1][ParIndex[1]+1][ParIndex[2]+1][0]+acc[ParIndex[0]+1][ParIndex[1]+1][ParIndex[2]  ][0]+acc[ParIndex[0]+1][ParIndex[1]  ][ParIndex[2]+1][0]+acc[ParIndex[0]+1][ParIndex[1]  ][ParIndex[2]  ][0]);
            a[i][1] = 0.25*(1-ParFrax[1])*(acc[ParIndex[0]+1][ParIndex[1]  ][ParIndex[2]+1][1]+acc[ParIndex[0]+1][ParIndex[1]  ][ParIndex[2]  ][1]+acc[ParIndex[0]  ][ParIndex[1]  ][ParIndex[2]+1][1]+acc[ParIndex[0]  ][ParIndex[1]  ][ParIndex[2]  ][1])
                     +0.25*    ParFrax[1]*(acc[ParIndex[0]+1][ParIndex[1]+1][ParIndex[2]+1][1]+acc[ParIndex[0]+1][ParIndex[1]+1][ParIndex[2]  ][1]+acc[ParIndex[0]  ][ParIndex[1]+1][ParIndex[2]+1][1]+acc[ParIndex[0]  ][ParIndex[1]+1][ParIndex[2]  ][1]);
            a[i][2] = 0.25*(1-ParFrax[2])*(acc[ParIndex[0]+1][ParIndex[1]+1][ParIndex[2]  ][2]+acc[ParIndex[0]+1][ParIndex[1]  ][ParIndex[2]  ][2]+acc[ParIndex[0]  ][ParIndex[1]+1][ParIndex[2]  ][2]+acc[ParIndex[0]  ][ParIndex[1]  ][ParIndex[2]  ][2])
                     +0.25*    ParFrax[2]*(acc[ParIndex[0]+1][ParIndex[1]+1][ParIndex[2]+1][2]+acc[ParIndex[0]+1][ParIndex[1]  ][ParIndex[2]+1][2]+acc[ParIndex[0]  ][ParIndex[1]+1][ParIndex[2]+1][2]+acc[ParIndex[0]  ][ParIndex[1]  ][ParIndex[2]+1][2]);
       }
    }
    return;
}// FUNCTION Acceleration


// FUNCTION Update: Update the system by dt
void Update( double x[ParN][3], double v[ParN][3] ){
   double a[ParN][3];     // acceleration of the particle
   if (Scheme_OI == 1){
     Acceleration( x, a );
     for ( int i = 0; i < ParN; i = i + 1 ){
        for ( int d = 0; d < 3; d = d + 1 ){
           v[i][d] += a[i][d]*dt/2;
        }
     } //////kick
     for ( int i = 0; i < ParN; i = i + 1 ){
        for ( int d = 0; d < 3; d = d + 1 ){
           x[i][d] += v[i][d]*dt;
        }
     } /////Drift
     CheckBoundary( x, v );
     Acceleration( x,  a );
     for ( int i = 0; i < ParN; i = i + 1 ){
        for ( int d = 0; d < 3; d = d + 1 ){
           v[i][d] += a[i][d]*dt/2;
        }
     } /////Kick
   } ///KDK
   else if (Scheme_OI == 2){
      for ( int i = 0; i < ParN; i = i + 1 ){
         for ( int d = 0; d < 3; d = d + 1 ){
            x[i][d] += v[i][d]*dt/2;
         }
      } /////Drift
      CheckBoundary( x, v );
      Acceleration( x, a );
      for ( int i = 0; i < ParN; i = i + 1 ){
         for ( int d = 0; d < 3; d = d + 1 ){
            v[i][d] += a[i][d]*dt;
         }
      } /////Kick
      for ( int i = 0; i < ParN; i = i + 1 ){
         for ( int d = 0; d < 3; d = d + 1 ){
            x[i][d] += v[i][d]*dt/2;
         }
      } /////Drift
   } ///DKD
   else if (Scheme_OI == 3){
      double k[4][ParN][3][2];   ////k_n,Par_ID,dim,pos/vel
      Acceleration( x, a );
      for ( int i = 0; i < ParN; i = i + 1 ){
         for (int d = 0; d < 3; d = d+1){
           k[0][i][d][0] = v[i][d];
           k[0][i][d][1] = a[i][d];
           x[i][d] += (dt/2)* k[0][i][d][0];
           v[i][d] += (dt/2)* k[0][i][d][1];
         }
      } /////////// Step1
      CheckBoundary( x, v );
      Acceleration( x, a );
      for ( int i = 0; i < ParN; i = i + 1 ){
         for (int d = 0; d < 3; d = d+1){
            k[1][i][d][0] = v[i][d];
            k[1][i][d][1] = a[i][d];
            x[i][d] += (dt/2)* k[1][i][d][0] -(dt/2)* k[0][i][d][0];
            v[i][d] += (dt/2)* k[1][i][d][1] -(dt/2)* k[0][i][d][1];
         }
      }/////////// Step2
      CheckBoundary( x, v );
      Acceleration( x, a );
      for ( int i = 0; i < ParN; i = i + 1 ){
         for (int d = 0; d < 3; d = d+1){
            k[2][i][d][0] = v[i][d];
            k[2][i][d][1] = a[i][d];
            x[i][d] += (dt)* k[2][i][d][0] -(dt/2)* k[1][i][d][0];
            v[i][d] += (dt)* k[2][i][d][1] -(dt/2)* k[1][i][d][1];
         }
      }/////////// Step3
      CheckBoundary( x, v );
      Acceleration( x, a );
      for ( int i = 0; i < ParN; i = i + 1 ){
         for (int d = 0; d < 3; d = d+1){
            k[3][i][d][0] = v[i][d];
            k[3][i][d][1] = a[i][d];
            x[i][d] += (dt/6)*(k[0][i][d][0] + 2*k[1][i][d][0] +2*k[2][i][d][0] + k[3][i][d][0]) -(dt)* k[2][i][d][0];
            v[i][d] += (dt/6)*(k[0][i][d][1] + 2*k[1][i][d][1] +2*k[2][i][d][1] + k[3][i][d][1]) -(dt)* k[2][i][d][1];
         }
      }/////////// Step4
   } ///RK4
   else {
     printf("Wrong Scheme_OI Mode");
   }
    /* To be finished */

    return;
}// FUNCTION Update


// FUNCTION Energy: Calculate the total energy of the system
double Energy( double x[ParN][3], double v[ParN][3] ){
    double e;                // total energy
    double ekin = 0.0;       // kinetic energy
    double epot = 0.0;       // potentail energy

    // kineteic energy
    for(int i=0;i<ParN;i++){
       ekin += 0.5*ParM[i]*( v[i][0]*v[i][0] + v[i][1]+v[i][1] + v[i][2]*v[i][2] );
    }

    // potential energy
    double Rho[N][N][N];      // density
    double Phi[N][N][N];      // potential
    MassDeposition( x, Rho ); // get density
    PoissonSolver( Rho, Phi );// get potential
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
    for(int k=0;k<N;k++){
       epot += 0.5*Rho[i][j][k]*Phi[i][j][k]*dx*dx*dx;
    }}}

    // total energy
    e = ekin + epot;

    return e;
}//FUNCTION Energy


// FUNCTION Momentum: Calculate the total momentum of the system
double Momentum( double v[ParN][3] ){
    double p;                // total momentum
    double px = 0.0;         // momentum in x-direction
    double py = 0.0;         // momentum in y-direction
    double pz = 0.0;         // momentum in z-direction

    for(int i=0;i<ParN;i++){ // sum over all particles
       px += ParM[i]*v[i][0];
       py += ParM[i]*v[i][1];
       pz += ParM[i]*v[i][2];
    }

    p = sqrt( px*px + py*py + pz*pz );

    return p;
}// FUNCTION Momentum


// FUNCTION CheckBoundary: Deal with the particles outside the box
void CheckBoundary( double x[ParN][3], double v[ParN][3] ){
    for ( int i = 0; i < ParN; i = i + 1 ){
       if (BC == 1){
          for ( int d = 0; d < 3; d = d + 1 ){
             if (x[i][d] < 0 ){
                x[i][d] += L;
             }
             if (x[i][d] > L ){
                x[i][d] -= L;
             }
          }
       }
       else if (BC == 2){
          if (x[i][0] < 0 || x[i][1] < 0 || x[i][2] < 0 || x[i][0] > L || x[i][1] > L || x[i][2] > L ){
             x[i][0] = 1000*L;
             x[i][1] = 1000*L;
             x[i][2] = 1000*L;
             v[i][0] = 0;
             v[i][1] = 0;
             v[i][2] = 0;
          }
       }
       else {
          printf("Wrong BC");
       }
    }


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
