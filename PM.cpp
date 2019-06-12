#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>
#include <ctime>
#include <omp.h>
#include <mpi.h>


//-----------------------------------
// A Partilce Mesh Code
// ----------------------------------
/*  Parallelizaion to be finished */

// constants
const double L = 1.0;           // length of the 3-D domain box
const int    N = 64;            // number of grid in each direction
const double dx = L/N;          // spatial resolution
const double dt = 0.001;        // time step
const int    ParN  = 10;        // number of particles
const double G = 1.0;           // gravitational constant
const double end_time = 10.0;   // end time of the evolution
static double *ParM = NULL;     // mass of each particle
const int NThread = 4;          // numer of threads

// schemes
const int BC = 1;               // boundary condition ( 1=Periodic, 2=Isolated )
const int Scheme_SG = 2;        // scheme of self-gravity ( 1=PM, 2=DN )
const int Scheme_MD = 2;        // scheme of mass deposition ( 1=NGP, 2=CIC, 3=TSC )
const int Scheme_PS = 1;        // scheme of poisson solver ( 1=FFT )
const int Scheme_OI = 2;        // scheme of orbit integration ( 1=KDK, 2=DKD, 3=RK4 )


// FUNCTION Init: Set the initial condition
void Init( double *x, double *v, const int NRank, const int MyRank ){
    /* To be modified for the test problem */
    ParM = new double[ParN];  // all particles
    srand(100);
    for(int i=0;i<ParN/NRank;i++){
        x[i*3+0] = (double)rand()/RAND_MAX*0.5*L+0.25*L;
        x[i*3+1] = (double)rand()/RAND_MAX*0.5*L+0.25*L;
        x[i*3+2] = (double)rand()/RAND_MAX*0.5*L+0.25*L;
        v[i*3+0] = ((double)rand()/RAND_MAX*-0.5)*0.2*L;
        v[i*3+1] = ((double)rand()/RAND_MAX*-0.5)*0.2*L;
        v[i*3+2] = ((double)rand()/RAND_MAX*-0.5)*0.2*L;
    }
    for(int i=0;i<ParN;i++)
        ParM[i] = 1.0;

    return;
}// FUNCTION Init


// FUNCTION timespec diff: Caltulate the time difference
struct timespec diff( struct timespec start, struct timespec end ){
    struct timespec temp;
    if( (end.tv_nsec-start.tv_nsec)<0 ){
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    }
    else{
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }

    return temp;
}// FUNCTION timespec diff


// FUNCTION CheckBoundary: Deal with the particles outside the box
void CheckBoundary( double *x, double *v, const int NRank ){
#   pragma omp parallel for num_threads(NThread)
    for(int i=0;i<ParN/NRank;i=i+1){
        if( BC==1 ){      // Periodic Boundary Condition
            int offset;
            for(int d=0;d<3;d=d+1){
                if( x[i*3+d]<0 ){
                    offset = (-x[i*3+d])/L;
                    x[i*3+d] += (offset+1)*L;
                }
                if( x[i*3+d]>L ){
                    offset = x[i*3+d]/L;
                    x[i*3+d] -= offset*L;
                }
            }
        }
        else if( BC==2 ){ // Isolated Boundary Condition
            if( x[i*3+0]<0 || x[i*3+1]<0 || x[i*3+2]<0 || x[i*3+0]>L || x[i*3+1]>L || x[i*3+2]>L ){
                x[i*3+0] = 1000*L;
                x[i*3+1] = 1000*L;
                x[i*3+2] = 1000*L;
                v[i*3+0] = 0;
                v[i*3+1] = 0;
                v[i*3+2] = 0;
            }
       }
       else{
          printf("Wrong BC\n");
       }
    }

    return;
}// FUNCTION CheckBoundary


// FUNCTION Index: Convert 3D index in N*N*N box into 1D index
int Index( int i, int j, int k ){
    return k+N*(j+N*i);
}// FUNCTION Index


// FUNCTION MassDeposition: Deposit particle mass onto grids
void MassDeposition( double *x, double *rho, const int NRank, const int MyRank ){
#   pragma omp parallel for collapse( 3 ) num_threads(NThread)
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
    for(int k=0;k<N;k++){
        rho[Index( i, j, k )] = 0.0;   // initialization as zero
    }}}

    if( Scheme_MD==1 ){       // Nearest-Grid-Point
#       pragma omp parallel for num_threads(NThread)
        for(int n=0;n<ParN/NRank;n++){
            if( x[n*3+0]>dx && x[n*3+0]<L-dx && x[n*3+1]>dx && x[n*3+1]<L-dx && x[n*3+2]>dx && x[n*3+2]<L-dx )
#           pragma omp critical
            rho[Index( (int)(x[n*3+0]/dx), (int)(x[n*3+1]/dx), (int)(x[n*3+2]/dx) )] += ParM[n+MyRank*ParN/NRank]/(dx*dx*dx);
        } 
    }
    else if( Scheme_MD==2 ){  // Cloud-In-Cell
#       pragma omp parallel for num_threads(NThread)
        for(int n=0;n<ParN/NRank;n++){
            double weighting[3][2];
            if( x[n*3+0]>1.5*dx && x[n*3+0]<L-1.5*dx && x[n*3+1]>1.5*dx && x[n*3+1]<L-1.5*dx && x[n*3+2]>1.5*dx && x[n*3+2]<L-1.5*dx ){

            for(int i=0;i<3;i++){  // weighting of distribution
                weighting[i][1] = (x[n*3+i]/dx-0.5) - (int)(x[n*3+i]/dx-0.5);
                weighting[i][0] = 1 - weighting[i][1];
            }

            for(int i=0;i<2;i++)
            for(int j=0;j<2;j++)
            for(int k=0;k<2;k++) // deposit density into 8 cells
#               pragma omp critical
                rho[Index( (int)(x[n*3+0]/dx-0.5)+i, (int)(x[n*3+1]/dx-0.5)+j, (int)(x[n*3+2]/dx-0.5)+k )] += weighting[0][i]*weighting[1][j]*weighting[2][k]*ParM[n+MyRank*ParN/NRank]/(dx*dx*dx);
            }
        }
    }
    else if( Scheme_MD==3 ){  // Triangular-Shaped-Cloud
#       pragma omp parallel for num_threads(NThread)
        for(int n=0;n<ParN/NRank;n++){
            double devide[3];
            double weighting[3][3];
            if( x[n*3+0]>2.0*dx && x[n*3+0]<L-2.0*dx && x[n*3+1]>2.0*dx && x[n*3+1]<L-2.0*dx && x[n*3+2]>2.0*dx && x[n*3+2]<L-2.0*dx ){
            for(int i=0;i<3;i++){  // weighting of distribution
                devide[i]= x[n*3+i]/dx - (int)(x[n*3+i]/dx);//find the distance between n grid and center of particle
                weighting[i][0] = pow((1.0-devide[i]),2)*0.5;
                weighting[i][1] = 0.5*(1.0 + 2.0*devide[i]-2.0*pow(devide[i],2));
                weighting[i][2] = pow(devide[i],2)*0.5;
            }

            for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            for(int k=0;k<3;k++) // deposit density into 27 cells
#               pragma omp critical
                rho[Index( (int)(x[n*3+0]/dx-1.0)+i, (int)(x[n*3+1]/dx-1.0)+j, (int)(x[n*3+2]/dx-1.0)+k )] += weighting[0][i]*weighting[1][j]*weighting[2][k]*ParM[n+MyRank*ParN/NRank]/(dx*dx*dx);
            }
        }
    }

    return;
}// FUNCTION MassDeposition


// FUNCTION PoissonSolver: Solve the Poisson equation to get the potential
void PoissonSolver( double *rho, double *phi ){
    if( BC==1 ){            // Periodic Boundary Condition
        if( Scheme_PS==1 ){ // Fast Fourier Transform

        fftw_complex *rho_K;  // density   in K space
        rho_K = (fftw_complex*) fftw_malloc( N*N*(N/2+1) * sizeof(fftw_complex) );
        fftw_complex *phi_K;  // potential in K space
        phi_K = (fftw_complex*) fftw_malloc( N*N*(N/2+1) * sizeof(fftw_complex) );

        fftw_plan rhoXtorhoK, phiKtophiX;
        rhoXtorhoK = fftw_plan_dft_r2c_3d( N, N, N, rho, rho_K, FFTW_ESTIMATE ); // Fourier Transform from rho(x) to rho(k)
        phiKtophiX = fftw_plan_dft_c2r_3d( N, N, N, phi_K, phi, FFTW_ESTIMATE ); // Inverse Fourier Transform from phi(k) to phi(x)

        fftw_execute( rhoXtorhoK ); // Fourier Transform from rho(x) to rho(k)

#       pragma omp parallel for collapse( 3 ) num_threads(NThread)
        for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
        for(int k=0;k<N/2+1;k++){
            double KK;            // K*K = Kx*Kx + Ky*Ky + Kz*Kz
            int nx, ny;
            if( i>N/2 ) nx=i-N; else nx=i;
            if( j>N/2 ) ny=j-N; else ny=j;
            KK = (nx*2*M_PI/L)*(nx*2*M_PI/L) + (ny*2*M_PI/L)*(ny*2*M_PI/L) + (k*2*M_PI/L)*(k*2*M_PI/L); // K*K = Kx*Kx + Ky*Ky + Kz*Kz
            phi_K[k+(N/2+1)*(j+N*i)][0] = rho_K[k+(N/2+1)*(j+N*i)][0]/KK*4.0*M_PI*G*(-1.0/(N*N*N)); // real part phi(k) = -4*Pi*G*rho(k)/k^2 and normalize by 1/N^3
            phi_K[k+(N/2+1)*(j+N*i)][1] = rho_K[k+(N/2+1)*(j+N*i)][1]/KK*4.0*M_PI*G*(-1.0/(N*N*N)); // imag part phi(k) = -4*Pi*G*rho(k)/k^2 and normalize by 1/N^3
        }}}
        phi_K[0][0] = 0.0;       // set DC to 0
        phi_K[0][1] = 0.0;       // set DC to 0
        fftw_execute( phiKtophiX ); // Inverse Fourier Transform from phi(k) to phi(x)

        fftw_destroy_plan( rhoXtorhoK );
        fftw_destroy_plan( phiKtophiX );

        fftw_free( rho_K );
        fftw_free( phi_K );

        }
        else printf("ERROR: Unsupported Scheme_PS!\n");
    }
    else if( BC==2 ){       // Isolated Boundary Condition
        if( Scheme_PS==1 ){ // Fast Fourier Transform

        double *mas_0pad;          // zero padding on the mass array
        mas_0pad = (double*) fftw_malloc( (2*N)*(2*N)*(2*N) * sizeof(double) );
        double *greensfn;          // symmetric discrete Green's function -1/r
        greensfn = (double*) fftw_malloc( (2*N)*(2*N)*(2*N) * sizeof(double) );
        double *phi_0pad;          // zero padding on the potential array
        phi_0pad = (double*) fftw_malloc( (2*N)*(2*N)*(2*N) * sizeof(double) );
        fftw_complex *mas_0pad_K;  // mass in K space
        mas_0pad_K = (fftw_complex*) fftw_malloc( (2*N)*(2*N)*(N+1) * sizeof(fftw_complex) );
        fftw_complex *greensfn_K;  // Green's function in K space
        greensfn_K = (fftw_complex*) fftw_malloc( (2*N)*(2*N)*(N+1) * sizeof(fftw_complex) );
        fftw_complex *phi_0pad_K;  // potential in K space
        phi_0pad_K = (fftw_complex*) fftw_malloc( (2*N)*(2*N)*(N+1) * sizeof(fftw_complex) );

        fftw_plan masXtomasK, greXtogreK, phiKtophiX;
        masXtomasK = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, mas_0pad, mas_0pad_K, FFTW_ESTIMATE ); // Fourier Transform from mas(x) to mas(k)
        greXtogreK = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, greensfn, greensfn_K, FFTW_ESTIMATE ); // Fourier Transform from gre(x) to gre(k)
        phiKtophiX = fftw_plan_dft_c2r_3d( 2*N, 2*N, 2*N, phi_0pad_K, phi_0pad, FFTW_ESTIMATE ); // Inverse Fourier Transform from phi(k) to phi(x)

#      pragma omp parallel for collapse( 3 ) num_threads(NThread)
        for(int i=0;i<2*N;i++){
        for(int j=0;j<2*N;j++){
        for(int k=0;k<2*N;k++){
            int nx,ny,nz;              // number of distance interval
            if( i<N && j<N && k<N ) mas_0pad[k+(2*N)*(j+(2*N)*i)] = rho[Index( i, j, k )]*dx*dx*dx;  // mass of cell = density * cell volume
            else mas_0pad[k+(2*N)*(j+(2*N)*i)] = 0.0;                                       // zero padding

            if( i>=N ) nx=2*N-i; else nx=i;  // symmetrization
            if( j>=N ) ny=2*N-j; else ny=j;  // symmetrization
            if( k>=N ) nz=2*N-k; else nz=k;  // symmetrization

            if( i==0 && j==0 && k==0 ) greensfn[0] = 0.0;  // ignore self force
            else greensfn[k+(2*N)*(j+(2*N)*i)] = -1.0/(dx*sqrt( nx*nx + ny*ny + nz*nz )); // -1/r
        }}}

        fftw_execute( masXtomasK ); // Fourier Transform from mas(x) to mas(k)
        fftw_execute( greXtogreK ); // Fourier Transform from gre(x) to gre(k)

#       pragma omp parallel for num_threads(NThread)
        for(int i=0;i<2*N;i++){
        for(int j=0;j<2*N;j++){
        for(int k=0;k<N+1;k++){
            phi_0pad_K[k+(N+1)*(j+(2*N)*i)][0] = G*( mas_0pad_K[k+(N+1)*(j+(2*N)*i)][0]*greensfn_K[k+(N+1)*(j+(2*N)*i)][0]-mas_0pad_K[k+(N+1)*(j+(2*N)*i)][1]*greensfn_K[k+(N+1)*(j+(2*N)*i)][1])*(1.0/(2*N*2*N*2*N));// real part phi(k) = G*mas(k)*gre(k) and normalize by 1/(2N)^3
            phi_0pad_K[k+(N+1)*(j+(2*N)*i)][1] = G*( mas_0pad_K[k+(N+1)*(j+(2*N)*i)][0]*greensfn_K[k+(N+1)*(j+(2*N)*i)][1]+mas_0pad_K[k+(N+1)*(j+(2*N)*i)][1]*greensfn_K[k+(N+1)*(j+(2*N)*i)][0])*(1.0/(2*N*2*N*2*N));// imag part phi(k) = G*mas(k)*gre(k) and normalize by 1/(2N)^3
        }}}

        fftw_execute( phiKtophiX ); // Inverse Fourier Transform from phi(k) to phi(x)

#       pragma omp parallel for collapse( 3 ) num_threads(NThread)
        for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
        for(int k=0;k<N;k++){
            phi[Index( i, j, k )] = phi_0pad[k+(2*N)*(j+(2*N)*i)]; // remove the padding
        }}}

        fftw_destroy_plan( masXtomasK );
        fftw_destroy_plan( greXtogreK );
        fftw_destroy_plan( phiKtophiX );

        fftw_free( mas_0pad_K );
        fftw_free( greensfn_K );
        fftw_free( phi_0pad_K );
        fftw_free( mas_0pad );
        fftw_free( greensfn );
        fftw_free( phi_0pad );

        }
        else printf("ERROR: Unsupported Scheme_PS!\n");
    }
    else printf("ERROR: Unsuppoted Boundary Condition!\n");

    return;
}// FUNCTION PoissonSolver


// FUNCTION Acceleration: Calculate the acceleration of each particle
void Acceleration( double *x, double *a, const int NRank, const int MyRank ){
#   pragma omp parallel for collapse( 2 ) num_threads(NThread)
    for(int i=0;i<ParN/NRank;i++)
    for(int j=0;j<3;j++)
        a[i*3+j] = 0.0;   // initialization as zero

    if( Scheme_SG==1 ){ // Particle Mesh
        double *Rho = new double[N*N*N]; // density
        double *Phi = new double[N*N*N]; // potential
        struct timespec MD_start, MD_end, PS_start, PS_end;
        double MD_time_used, PS_time_used;

        clock_gettime(CLOCK_MONOTONIC, &MD_start);
        /////////measure start
        MassDeposition( x, Rho, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &MD_end);
        struct timespec MD_temp = diff(MD_start, MD_end);
        MD_time_used = MD_temp.tv_sec + (double) MD_temp.tv_nsec / 1000000000.0;
        printf("MD_Time = %f\n", MD_time_used);

        clock_gettime(CLOCK_MONOTONIC, &PS_start);
        /////////measure start
        PoissonSolver( Rho, Phi );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &PS_end);
        struct timespec PS_temp = diff(PS_start, PS_end);
        PS_time_used = PS_temp.tv_sec + (double) PS_temp.tv_nsec / 1000000000.0;
        printf("PS_Time = %f\n", PS_time_used);

        double *XSendBuf = new double[ParN/NRank*3];  // Buffer for sending x between MPI ranks
        double *ASendBuf = new double[ParN/NRank*3];  // Buffer for sending a between MPI ranks
        MPI_Request Request[2];
        for(int rank=0;rank<NRank;rank++){

        if( Scheme_MD==1 ){       // Nearest-Grid-Point
#           pragma omp parallel for num_threads(NThread)
            for(int n=0;n<ParN/NRank;n++){
                if( x[n*3+0]>dx && x[n*3+0]<L-dx && x[n*3+1]>dx && x[n*3+1]<L-dx && x[n*3+2]>dx && x[n*3+2]<L-dx ){

                // x-direction
                    a[n*3+0] += -( Phi[Index( (int)(x[n*3+0]/dx)+1, (int)(x[n*3+1]/dx), (int)(x[n*3+2]/dx) )] - Phi[Index( (int)(x[n*3+0]/dx)-1, (int)(x[n*3+1]/dx), (int)(x[n*3+2]/dx) )] ) /(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank]);

                // y-direction
                    a[n*3+1] += -( Phi[Index( (int)(x[n*3+0]/dx), (int)(x[n*3+1]/dx)+1, (int)(x[n*3+2]/dx) )] - Phi[Index( (int)(x[n*3+0]/dx), (int)(x[n*3+1]/dx-1), (int)(x[n*3+2]/dx) )] ) /(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank]);

                // z-direction
                    a[n*3+2] += -( Phi[Index( (int)(x[n*3+0]/dx), (int)(x[n*3+1]/dx), (int)(x[n*3+2]/dx)+1 )] - Phi[Index( (int)(x[n*3+0]/dx), (int)(x[n*3+1]/dx), (int)(x[n*3+2]/dx)-1 )] ) /(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank]);
                }
            }
        }
        else if( Scheme_MD==2 ){  // Cloud-In-Cell
#           pragma omp parallel for num_threads(NThread)
            for(int n=0;n<ParN/NRank;n++){
                double weighting[3][2];
                if( x[n*3+0]>1.5*dx && x[n*3+0]<L-1.5*dx && x[n*3+1]>1.5*dx && x[n*3+1]<L-1.5*dx && x[n*3+2]>1.5*dx && x[n*3+2]<L-1.5*dx ){

                for(int i=0;i<3;i++){ // weighting of distribution
                    weighting[i][1] = (x[n*3+i]/dx-0.5) - (int)(x[n*3+i]/dx-0.5);
                    weighting[i][0] = 1 - weighting[i][1];
                }

                for(int i=0;i<2;i++){
                for(int j=0;j<2;j++){
                for(int k=0;k<2;k++){ // get acceleration from 8 cells

                // x-direction
                    a[n*3+0] += -( Phi[Index( (int)(x[n*3+0]/dx-0.5)+i+1, (int)(x[n*3+1]/dx-0.5)+j, (int)(x[n*3+2]/dx-0.5)+k )] - Phi[Index( (int)(x[n*3+0]/dx-0.5)+i-1, (int)(x[n*3+1]/dx-0.5)+j, (int)(x[n*3+2]/dx-0.5)+k) ] )/(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank])*weighting[0][i]*weighting[1][j]*weighting[2][k];

                // y-direction
                    a[n*3+1] += -( Phi[Index( (int)(x[n*3+0]/dx-0.5)+i, (int)(x[n*3+1]/dx-0.5)+j+1, (int)(x[n*3+2]/dx-0.5)+k )] - Phi[Index( (int)(x[n*3+0]/dx-0.5)+i, (int)(x[n*3+1]/dx-0.5)+j-1, (int)(x[n*3+2]/dx-0.5)+k )] )/(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank])*weighting[0][i]*weighting[1][j]*weighting[2][k];

                // z-direction
                    a[n*3+2] += -( Phi[Index( (int)(x[n*3+0]/dx-0.5)+i, (int)(x[n*3+1]/dx-0.5)+j, (int)(x[n*3+2]/dx-0.5)+k+1 )] - Phi[Index( (int)(x[n*3+0]/dx-0.5)+i, (int)(x[n*3+1]/dx-0.5)+j, (int)(x[n*3+2]/dx-0.5)+k-1 )] )/(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank])*weighting[0][i]*weighting[1][j]*weighting[2][k];

                }}}
              }
           }
        }
        else if( Scheme_MD==3 ){  // Triangular-Shaped-Cloud
#           pragma omp parallel for num_threads(NThread)
            for(int n=0;n<ParN/NRank;n++){
                double weighting[3][3];
                double devide[3];
                if( x[n*3+0]>2.0*dx && x[n*3+0]<L-2.0*dx && x[n*3+1]>2.0*dx && x[n*3+1]<L-2.0*dx && x[n*3+2]>2.0*dx && x[n*3+2]<L-2.0*dx ){

                for(int i=0;i<3;i++){ // weighting of distribution
                    devide[i]= x[n*3+i]/dx - (int)(x[n*3+i]/dx);//find the distance between n grid and center of particle
                    weighting[i][0] = pow((1.0-devide[i]),2)*0.5;
                    weighting[i][1] = 0.5*(1.0 + 2.0*devide[i]-2.0*pow(devide[i],2));
                    weighting[i][2] = pow(devide[i],2)*0.5;
                }

                for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){ // get acceleration from 27 cells

                // x-direction
                    a[n*3+0] += -( Phi[Index( (int)(x[n*3+0]/dx-1.0)+i+1, (int)(x[n*3+1]/dx-1.0)+j, (int)(x[n*3+2]/dx-1.0)+k )] - Phi[Index( (int)(x[n*3+0]/dx-1.0)+i-1, (int)(x[n*3+1]/dx-1.0)+j, (int)(x[n*3+2]/dx-1.0)+k )] )/(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank])*weighting[0][i]*weighting[1][j]*weighting[2][k];

                // y-direction
                    a[n*3+1] += -( Phi[Index( (int)(x[n*3+0]/dx-1.0)+i, (int)(x[n*3+1]/dx-1.0)+j+1, (int)(x[n*3+2]/dx-1.0)+k )] - Phi[Index( (int)(x[n*3+0]/dx-1.0)+i, (int)(x[n*3+1]/dx-1.0)+j-1, (int)(x[n*3+2]/dx-1.0)+k )] )/(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank])*weighting[0][i]*weighting[1][j]*weighting[2][k];

                // z-direction
                    a[n*3+2] += -( Phi[Index( (int)(x[n*3+0]/dx-1.0)+i, (int)(x[n*3+1]/dx-1.0)+j, (int)(x[n*3+2]/dx-1.0)+k+1 )] - Phi[Index( (int)(x[n*3+0]/dx-1.0)+i, (int)(x[n*3+1]/dx-1.0)+j, (int)(x[n*3+2]/dx-1.0)+k-1 )] )/(2.0*dx*ParM[n+((MyRank+rank)%NRank)*ParN/NRank])*weighting[0][i]*weighting[1][j]*weighting[2][k];

                }}}
                }
            }
        }

        for(int i=0;i<ParN/NRank*3;i++){ // prepare buffers for MPI
            XSendBuf[i] = x[i];
            ASendBuf[i] = a[i];
        }
        MPI_Isend( XSendBuf, ParN/NRank*3, MPI_DOUBLE, (MyRank+NRank+1)%NRank, 123, MPI_COMM_WORLD, &Request[0] ); // send x to RHS Rank
        MPI_Isend( ASendBuf, ParN/NRank*3, MPI_DOUBLE, (MyRank+NRank+1)%NRank, 456, MPI_COMM_WORLD, &Request[0] ); // send a to RHS Rank
        MPI_Recv( x, ParN/NRank*3, MPI_DOUBLE, (MyRank+NRank-1)%NRank, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );   // recv x from LHS Rank
        MPI_Recv( a, ParN/NRank*3, MPI_DOUBLE, (MyRank+NRank-1)%NRank, 456, MPI_COMM_WORLD, MPI_STATUS_IGNORE );   // recv a from LHS Rank

        }
        delete [] XSendBuf;
        delete [] ASendBuf;

        delete [] Rho;
        delete [] Phi;
    }
    else if( Scheme_SG==2 && NRank==1 ){  // Direct N-body
#   pragma omp parallel for num_threads(NThread)
    for(int i=0;i<ParN;i++){
        if( x[i*3+0]>0 && x[i*3+0]<L && x[i*3+1]>0 && x[i*3+1]<L && x[i*3+2]>0 && x[i*3+2]<L ){
        for(int j=0;j<ParN;j++){
            double r;  // distance between two particles
            if( i!=j && x[j*3+0]>0 && x[j*3+0]<L && x[j*3+1]>0 && x[j*3+1]<L && x[j*3+2]>0 && x[j*3+2]<L ){
                r = sqrt( (x[j*3+0]-x[i*3+0])*(x[j*3+0]-x[i*3+0]) + (x[j*3+1]-x[i*3+1])*(x[j*3+1]-x[i*3+1]) + (x[j*3+2]-x[i*3+2])*(x[j*3+2]-x[i*3+2]) );
#               pragma omp critical
                a[i*3+0] += G*ParM[j]/(r*r)*(x[j*3+0]-x[i*3+0])/r;
#               pragma omp critical
                a[i*3+1] += G*ParM[j]/(r*r)*(x[j*3+1]-x[i*3+1])/r;
#               pragma omp critical
                a[i*3+2] += G*ParM[j]/(r*r)*(x[j*3+2]-x[i*3+2])/r;
            }
        }
        }
    }
    }

    return;
}// FUNCTION Acceleration


// FUNCTION Update: Update the system by dt
void Update( double *x, double *v, const int NRank, const int MyRank ){
    double *a= new double[ParN/NRank*3];     // acceleration of the particle
    struct timespec AC1_start, AC1_end, AC2_start, AC2_end, AC3_start, AC3_end, AC4_start, AC4_end;
    double AC1_time_used, AC2_time_used, AC3_time_used, AC4_time_used;
    if( Scheme_OI==1 ){      // Kick-Drift-Kick
        clock_gettime(CLOCK_MONOTONIC, &AC1_start);
        /////////measure start
        CheckBoundary( x, v );
        Acceleration( x, a, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &AC1_end);
        struct timespec AC1_temp = diff(AC1_start, AC1_end);
        AC1_time_used = AC1_temp.tv_sec + (double) AC1_temp.tv_nsec / 1000000000.0;
        printf("AC_Time = %f\n", AC1_time_used);

#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                v[i*3+d] += a[i*3+d]*dt/2;
            }
        } // Kick
#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                x[i*3+d] += v[i*3+d]*dt;
            }
        } // Drift

        clock_gettime(CLOCK_MONOTONIC, &AC2_start);
        /////////measure start
        CheckBoundary( x, v, NRank );
        Acceleration( x, a, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &AC2_end);
        struct timespec AC2_temp = diff(AC2_start, AC2_end);
        AC2_time_used = AC2_temp.tv_sec + (double) AC2_temp.tv_nsec / 1000000000.0;
        printf("AC2_Time = %f\n", AC2_time_used);

#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                v[i*3+d] += a[i*3+d]*dt/2;
            }
        } // Kick
    }
    else if( Scheme_OI==2 ){ // Drift-Kick-Drift
#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                x[i*3+d] += v[i*3+d]*dt/2;
            }
        } // Drift
        clock_gettime(CLOCK_MONOTONIC, &AC1_start);
        /////////measure start
        CheckBoundary( x, v, NRank );
        Acceleration( x, a, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &AC1_end);
        struct timespec AC1_temp = diff(AC1_start, AC1_end);
        AC1_time_used = AC1_temp.tv_sec + (double) AC1_temp.tv_nsec / 1000000000.0;
        printf("AC_Time = %f\n", AC1_time_used);

#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                v[i*3+d] += a[i*3+d]*dt;
            }
        } // Kick
#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                x[i*3+d] += v[i*3+d]*dt/2;
            }
        } // Drift
    }
    else if( Scheme_OI==3 ){ // Runge-Kutta 4th
        double *rk = new double[4*ParN/NRank*3*2];   // k_n, Par_ID, dim, pos/vel
        clock_gettime(CLOCK_MONOTONIC, &AC1_start);
        /////////measure start
        CheckBoundary( x, v );
        Acceleration( x, a, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &AC1_end);
        struct timespec AC1_temp = diff(AC1_start, AC1_end);
        AC1_time_used = AC1_temp.tv_sec + (double) AC1_temp.tv_nsec / 1000000000.0;
        printf("AC_Time = %f\n", AC1_time_used);

#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                rk[((0*ParN/NRank+i)*3+d)*2+0] = v[i*3+d];
                rk[((0*ParN/NRank+i)*3+d)*2+1] = a[i*3+d];
                x[i*3+d] += (dt/2)* rk[((0*ParN/NRank+i)*3+d)*2+0];
                v[i*3+d] += (dt/2)* rk[((0*ParN/NRank+i)*3+d)*2+1];
            }
        } // Step1
        clock_gettime(CLOCK_MONOTONIC, &AC2_start);
        /////////measure start
        CheckBoundary( x, v, NRank );
        Acceleration( x, a, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &AC2_end);
        struct timespec AC2_temp = diff(AC2_start, AC2_end);
        AC2_time_used = AC2_temp.tv_sec + (double) AC2_temp.tv_nsec / 1000000000.0;
        printf("AC2_Time = %f\n", AC2_time_used);

#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                rk[((1*ParN/NRank+i)*3+d)*2+0] = v[i*3+d];
                rk[((1*ParN/NRank+i)*3+d)*2+1] = a[i*3+d];
                x[i*3+d] += (dt/2)* rk[((1*ParN/NRank+i)*3+d)*2+0] -(dt/2)* rk[((0*ParN/NRank+i)*3+d)*2+0];
                v[i*3+d] += (dt/2)* rk[((1*ParN/NRank+i)*3+d)*2+1] -(dt/2)* rk[((0*ParN/NRank+i)*3+d)*2+1];
            }
        } // Step2
        clock_gettime(CLOCK_MONOTONIC, &AC3_start);
        /////////measure start
        CheckBoundary( x, v, NRank );
        Acceleration( x, a, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &AC3_end);
        struct timespec AC3_temp = diff(AC3_start, AC3_end);
        AC3_time_used = AC3_temp.tv_sec + (double) AC3_temp.tv_nsec / 1000000000.0;
        printf("AC3_Time = %f\n", AC3_time_used);

#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                rk[((2*ParN/NRank+i)*3+d)*2+0] = v[i*3+d];
                rk[((2*ParN/NRank+i)*3+d)*2+1] = a[i*3+d];
                x[i*3+d] += (dt)* rk[((2*ParN/NRank+i)*3+d)*2+0] -(dt/2)* rk[((1*ParN/NRank+i)*3+d)*2+0];
                v[i*3+d] += (dt)* rk[((2*ParN/NRank+i)*3+d)*2+1] -(dt/2)* rk[((1*ParN/NRank+i)*3+d)*2+1];
            }
        } // Step3
        clock_gettime(CLOCK_MONOTONIC, &AC4_start);
        /////////measure start
        CheckBoundary( x, v, NRank );
        Acceleration( x, a, NRank, MyRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &AC4_end);
        struct timespec AC4_temp = diff(AC4_start, AC4_end);
        AC4_time_used = AC4_temp.tv_sec + (double) AC4_temp.tv_nsec / 1000000000.0;
        printf("AC4_Time = %f\n", AC4_time_used);

#       pragma omp parallel for collapse( 2 ) num_threads(NThread)
        for(int i=0;i<ParN/NRank;i=i+1){
            for(int d=0;d<3;d=d+1){
                rk[((3*ParN/NRank+i)*3+d)*2+0] = v[i*3+d];
                rk[((3*ParN/NRank+i)*3+d)*2+1] = a[i*3+d];
                x[i*3+d] += (dt/6)*(rk[((0*ParN/NRank+i)*3+d)*2+0] + 2*rk[((1*ParN/NRank+i)*3+d)*2+0] +2*rk[((2*ParN/NRank+i)*3+d)*2+0] + rk[((3*ParN/NRank+i)*3+d)*2+0]) -(dt)* rk[((2*ParN/NRank+i)*3+d)*2+0];
                v[i*3+d] += (dt/6)*(rk[((0*ParN/NRank+i)*3+d)*2+1] + 2*rk[((1*ParN/NRank+i)*3+d)*2+1] +2*rk[((2*ParN/NRank+i)*3+d)*2+1] + rk[((3*ParN/NRank+i)*3+d)*2+1]) -(dt)* rk[((2*ParN/NRank+i)*3+d)*2+1];
            }
        } // Step4
        delete [] rk;
    }
    else{
        printf("Wrong Scheme_OI Mode\n");
    }
    delete [] a;

    return;
}// FUNCTION Update


// FUNCTION Energy: Calculate the total energy of the system
double Energy( double *x, double *v ){
    double e;                // total energy
    double ekin = 0.0;       // kinetic energy
    double epot = 0.0;       // potentail energy

    // kineteic energy
#   pragma omp parallel for reduction( +:ekin ) num_threads(NThread)
    for(int i=0;i<ParN;i++){
        ekin += 0.5*ParM[i]*( v[i*3+0]*v[i*3+0] + v[i*3+1]+v[i*3+1] + v[i*3+2]*v[i*3+2] );
    }

    // potential energy
    double *Rho = new double[N*N*N]; // density
    double *Phi = new double[N*N*N]; // potential
    MassDeposition( x, Rho, 1, 0 );  // get density
    PoissonSolver( Rho, Phi );       // get potential

#   pragma omp parallel for collapse( 2 ) num_threads(NThread)
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
    for(int k=0;k<N;k++){
#       pragma omp critical
        epot += 0.5*Rho[Index( i, j, k )]*Phi[Index( i, j, k )]*dx*dx*dx;
    }}}
    delete [] Rho;
    delete [] Phi;

    // total energy
    e = ekin + epot;

    return e;
}//FUNCTION Energy


// FUNCTION Momentum: Calculate the total momentum of the system
double Momentum( double *v ){
    double p;                // total momentum
    double px = 0.0;         // momentum in x-direction
    double py = 0.0;         // momentum in y-direction
    double pz = 0.0;         // momentum in z-direction

#   pragma omp parallel for reduction( +:px,py,pz ) num_threads(NThread)
    for(int i=0;i<ParN;i++){ // sum over all particles
        px += ParM[i]*v[i*3+0];
        py += ParM[i]*v[i*3+1];
        pz += ParM[i]*v[i*3+2];
    }

    p = sqrt( px*px + py*py + pz*pz );

    return p;
}// FUNCTION Momentum


// FUNCTION OutputData: Output the results into files
void OutputData( double t, double *x, double e, double eerror, double p, double perror ){
    FILE *fp;
    fp = fopen( "Data_ParticlePosition", "a" );

    if( t==0.0 ){
        fprintf( fp, "#        Time");
        for(int i=0;i<ParN;i++) fprintf( fp,"            x_%d            y_%d            z_%d",i+1,i+1,i+1);
        fprintf( fp, "\n");
    }

    fprintf( fp, "%13.6e", t );
    for(int i=0;i<ParN;i++) fprintf( fp, "  %13.6e  %13.6e  %13.6e",x[i*3+0],x[i*3+1],x[i*3+2] );
    fprintf( fp, "\n");

    fclose(fp);

    FILE *fc;
    fc = fopen( "Data_ConservedQuantity", "a" );

    if( t==0.0 )
        fprintf( fc, "#        Time         Energy    EnergyError       Momentum  MomentumError\n");
    fprintf( fc, "%13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n", t, e, eerror, p, perror );

    fclose(fc);

    return;
}// FUNCTION OutputData


// FUNCTION main
int main( int argc, char *argv[] ){
    // initialize MPI
    int NRank, MyRank;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );    // My MPI rank ID
    MPI_Comm_size( MPI_COMM_WORLD, &NRank );     // Total number of MPI ranks
    if( ParN%NRank!=0 )  // ensure the rank number can divide the number of particles
    {
       fprintf( stderr, "ERROR: NRank (%d) cannot divide ParN=(%d)\n", NRank,ParN );
       MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    double *OutputXSendBuf = new double [ParN/NRank*3]; // Buffer for sending x between MPI ranks
    double *OutputVSendBuf = new double [ParN/NRank*3]; // Buffer for sending v between MPI ranks
    double *OutputXRecvBuf = new double [ParN*3];       // Buffer for gathering x between MPI ranks
    double *OutputVRecvBuf = new double [ParN*3];       // Buffer for gathering v between MPI ranks

    // variables
    double t = 0.0;           // time
    double *ParX = new double[ParN/NRank*3];     // position of the particle
    double *ParV = new double[ParN/NRank*3];     // velocity of the particle
    double E0;                // initial energy
    double E;                 // energy
    double Eerr;              // energy conservation error
    double P0;                // initial momentum |P0|
    double P;                 // moentum |P|
    double Perr;              // momentum conservation error
    struct timespec UD_start, UD_end; //timer
    double UD_time_used;

    // initialization
    Init( ParX, ParV, NRank, MyRank );
    for(int i=0;i<ParN/NRank*3;i++){ // prepare buffers
        OutputXSendBuf[i] = ParX[i];
        OutputVSendBuf[i] = ParV[i];
    }
    MPI_Gather( OutputXSendBuf, ParN/NRank*3, MPI_DOUBLE, OutputXRecvBuf, ParN/NRank*3,MPI_DOUBLE, 0, MPI_COMM_WORLD ); // gather x from all ranks
    MPI_Gather( OutputVSendBuf, ParN/NRank*3, MPI_DOUBLE, OutputVRecvBuf, ParN/NRank*3,MPI_DOUBLE, 0, MPI_COMM_WORLD ); // gather v from all ranks
    if( MyRank==0 ){
        E0 = Energy( OutputXRecvBuf, OutputVRecvBuf );
        P0 = Momentum( OutputVRecvBuf );
        OutputData( t, OutputXRecvBuf, E0, 0.0, P0, 0.0 );
    }

    // main loop
    if( MyRank==0 )printf(" Start!\n");
    while( t<end_time ){
        if( MyRank==0 )printf("Time: %13.6e -> %13.6e, dt=%f\n", t, t+dt, dt );

        clock_gettime(CLOCK_MONOTONIC, &UD_start);
        /////////measure start
        Update( ParX, ParV, NRank, MyRank );
        CheckBoundary( ParX, ParV, NRank );
        /////////measure end
        clock_gettime(CLOCK_MONOTONIC, &UD_end);
        struct timespec UD_temp = diff(UD_start, UD_end);
        UD_time_used = UD_temp.tv_sec + (double) UD_temp.tv_nsec / 1000000000.0;
        printf("UD_Time = %f\n", UD_time_used);

        t = t+dt;
        if( (int)(t/dt)%NRank==0 ){
        for(int i=0;i<ParN/NRank*3;i++){
            OutputXSendBuf[i] = ParX[i];
            OutputVSendBuf[i] = ParV[i];
        }
        MPI_Gather( OutputXSendBuf, ParN/NRank*3, MPI_DOUBLE, OutputXRecvBuf, ParN/NRank*3,MPI_DOUBLE, 0, MPI_COMM_WORLD ); // gather x from all ranks
        MPI_Gather( OutputVSendBuf, ParN/NRank*3, MPI_DOUBLE, OutputVRecvBuf, ParN/NRank*3,MPI_DOUBLE, 0, MPI_COMM_WORLD ); // gather v from all ranks
        if( MyRank==0 ){
            E = Energy( OutputXRecvBuf, OutputVRecvBuf );
            Eerr = (E-E0)/fabs(E0);
            P = Momentum( OutputVRecvBuf );
            Perr = (P-P0)/fabs(P0);
            OutputData( t, OutputXRecvBuf, E, Eerr, P, Perr );
        }
        }
    }
    if( MyRank==0 )printf(" End!\n\n" );
    delete [] ParX;
    delete [] ParV;
    delete [] ParM;

    // output the results
    if( MyRank==0 ){
    printf("Info:\n");
    printf("Number of Particles = %d\n", ParN);
    printf("Number of Grids     = %dx%dx%d\n",N,N,N);
    printf("Number of MPI Ranks = %d\n", NRank);
    printf("Number of Threads   = %d\n", NThread);
    switch( BC ){
        case 1:
            printf("Boundary Condition  = Periodic\n");
            break;
        case 2:
            printf("Boundary Condition  = Isolated\n");
            break;
    }
    switch( Scheme_SG ){
        case 1:
            printf("Self-Gravity        = PM\n");
            break;
        case 2:
            printf("Self-Gravity        = DN\n");
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
            printf("Orbit Integration   = RK4\n");
            break;
    }
    printf("\n");
    printf("Energy   Error: %13.6e\n", Eerr );
    printf("Momentum Error: %13.6e\n", Perr );
    }

    // terminate MPI
    MPI_Finalize();
    delete [] OutputXSendBuf;
    delete [] OutputVSendBuf;
    delete [] OutputXRecvBuf;
    delete [] OutputVRecvBuf;

    return EXIT_SUCCESS;
}// FUNCTION main
