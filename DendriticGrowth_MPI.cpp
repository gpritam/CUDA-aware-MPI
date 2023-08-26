//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Parallel program for dendritic solidification using Finite Difference method.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 7.03.2023
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=DendriticGrowth_MPI.cpp
// make run Nproc=1
//________________________________________________________________________________________________

#include "General.h"
#include "General_MPI_2D.h"

using namespace std::chrono;

#define TECPLOT
//#define VISIT

//________________________________________________________________________________________________
// Npx : Number of processors along x-direction
// Npy : Number of processors along y-direction
// 
// Nx : Total number of control volumes in x direction
// Ny : Total number of control volumes in y direction
//________________________________________________________________________________________________

int Nx = 512, Ny = 512, offset = 2, Npx, Npy;

const int SaveInterval = 50;

double Lx = 9.0, Ly = 9.0, dt = 5.0E-5, Tf = 0.4;

double tau = 3.0E-4, epsilonb = 1.0E-2, kappa = 1.8, delta = 2.0E-2;
double beta = 6.0, alpha = 0.9, Gamma = 10.0, Teq = 1.0, Theta0 = 0.2, r0 = 5.0*Lx*Lx/250000.0;

double **phi, **phiold, **T, **Told, **Epsilon, **Epsilon_theta, *Temporary, dx, dy, t = 0.0;
int TimeStep = 0;

char *s;

double *phi_gpu, *phiold_gpu, *T_gpu, *Told_gpu, *Epsilon_gpu, *Epsilon_theta_gpu;
double *Laplacian, *term0, *term1, *term2, *term3, *m, *Temporary_gpu;

//________________________________________________________________________________________________
// extern function prototypes
//________________________________________________________________________________________________
extern "C++" void Allocate_CUDA ( );
extern "C++" void Deallocate_CUDA ( );
extern "C++" void CPUtoGPU_CUDA ( );
extern "C++" void GPUtoCPU_CUDA ( );
extern "C++" void GPUtoCPUEfficient_CUDA ( double **&A, 
                                           double *A_gpu, 
                                           const int Nx, 
                                           const int Ny, 
                                           const int offset );
extern "C++" void CPUtoGPUEfficient_CUDA ( double **A, 
                                           double *A_gpu, 
                                           const int Nx, 
                                           const int Ny, 
                                           const int offset );
extern "C++" void UpdateOldValues_CUDA ( );
extern "C++" void ComputeEpsilon_CUDA ( );
extern "C++" void ComputePhiT_CUDA ( const double PI );

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SetSubDomains ( int &Npx, 
                     int &Npy, 
                     int &Nx, 
                     int &Ny )
{
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    bool reverse = false;
    double ratio = Lx/Ly;
    
    Npx = size;
    Npy = 1;
    
    if (Lx < Ly)
    {
        reverse = true;
        ratio = Ly/Lx;
    }
    
    int imax = floor(sqrt(double(size)));
    double distance = Absolute(ratio-size);
    
    for (int i{2}; i <= imax; ++i)
    {
        if (size % i == 0)
        {
            if (Absolute(ratio-size/i) < distance)
            {
                distance = Absolute(ratio-size/i);
                Npx = size/i;
                Npy = i;
            }
        }
    }
    
    if (reverse)
        Swap(Npx,Npy);
    
    Nx = ceil(double(Nx)/Npx);
    Ny = ceil(double(Ny)/Npy);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Allocate ()
{
    Allocate(phi,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(T,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(phiold,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(Told,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(Epsilon,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(Epsilon_theta,Nx+2*offset,Ny+2*offset,offset,offset);
    
    Allocate(Temporary,(2*Nx*offset+2*Ny*offset));
    
    Allocate(s,200);
    
    Allocate_CUDA();
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate ()
{
    Deallocate(phi,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(T,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(phiold,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(Told,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(Epsilon,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(Epsilon_theta,Nx+2*offset,Ny+2*offset,offset,offset);
    
    Deallocate(Temporary,(2*Nx*offset+2*Ny*offset));
    
    Deallocate(s,200);
    
    Deallocate_CUDA();
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CheckNaN ()
{
    bool ErrorOccured = false;
    
    for (int j{}; j < offset; ++j)
    {
        for (int i{}; i < Nx; ++i)
        {
            if (phi[i][j] != phi[i][j])
            {
                Deallocate();
                
                std::cout << "NaN encountered!" << std::endl;
                
                ErrorOccured = true;
                break;
            }
        }
        
        if (ErrorOccured)
            break;
    }
    
    if (!ErrorOccured)
    {
        for (int j{}; j < offset; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                if (phi[i][Ny-1-j] != phi[i][Ny-1-j])
                {
                    Deallocate();
                    
                    std::cout << "NaN encountered!" << std::endl;
                    
                    ErrorOccured = true;
                    break;
                }
            }
            
            if (ErrorOccured)
                break;
        }
    }
    
    if (!ErrorOccured)
    {
        for (int i{}; i < offset; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                if (phi[i][j] != phi[i][j])
                {
                    Deallocate();
                    
                    std::cout << "NaN encountered!" << std::endl;
                    
                    ErrorOccured = true;
                    break;
                }
            }
            
            if (ErrorOccured)
                break;
        }
    }
    
    if (!ErrorOccured)
    {
        for (int i{}; i < offset; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                if (phi[Nx-1-i][j] != phi[Nx-1-i][j])
                {
                    Deallocate();
                    
                    std::cout << "NaN encountered!" << std::endl;
                    
                    ErrorOccured = true;
                    break;
                }
            }
            
            if (ErrorOccured)
                break;
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (ErrorOccured)
        MPI_Abort(MPI_COMM_WORLD,1);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Initialize ()
{
    double x, y, r;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            x = (I*Nx+i)*dx;
            y = (J*Ny+j)*dy;
            
            r = (x-Lx/2.0)*(x-Lx/2.0) + (y-Ly/2.0)*(y-Ly/2.0);
            
            phi[i][j] = (r < r0 ? 1.0 : 0.0);
            
            T[i][j] = 0.0;
        }
    }
    
    ExchangeData(phi,Nx,Ny,offset);
    ExchangeData(T,Nx,Ny,offset);
}

//________________________________________________________________________________________________
// Main code
//________________________________________________________________________________________________
int main ( int argc, char *argv[] )
{
    MPI_Init(&argc, &argv);
    
    std::cout.flags( std::ios::dec | std::ios::fixed );
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SetSubDomains(Npx,Npy,Nx,Ny);
    
    SetMPIEnvironment MPI_ReadyToUse(Nx,Ny,Npx,Npy,offset);
    
    dx = Lx/(Npx*Nx);
    dy = Ly/(Npy*Ny);
    
    Allocate();
    
    Initialize();
    
    CPUtoGPU_CUDA();
    
    #ifdef VISIT
    if (Rank == 0)
    {
        std::ofstream Output("Output/Field.visit", std::ios::out);
        
        if ( !Output )
            ErrorMessage("Output file couldnot be opened!");
        
        Output << "!NBLOCKS " << size << std::endl;
        
        Output.close();
    }
    #endif
    
    // Update concentration
    while (t < Tf)
    {
        if (Rank == 0)
            std::cout << "Time step = " << TimeStep << ", Current time = " << t << std::endl;
        
        // Step 1
        UpdateOldValues_CUDA();
        
        // Step 2: Compute Epsilon and (d Epsilon/d theta)
        ComputeEpsilon_CUDA();
        
        GPUtoCPUEfficient_CUDA(Epsilon,Epsilon_gpu,Nx,Ny,offset);
        GPUtoCPUEfficient_CUDA(Epsilon_theta,Epsilon_theta_gpu,Nx,Ny,offset);
        
        ExchangeData(Epsilon,Nx,Ny,offset);
        ExchangeData(Epsilon_theta,Nx,Ny,offset);
        
        CPUtoGPUEfficient_CUDA(Epsilon,Epsilon_gpu,Nx,Ny,offset);
        CPUtoGPUEfficient_CUDA(Epsilon_theta,Epsilon_theta_gpu,Nx,Ny,offset);
        
        // Step 3: Update phase-field parameter and temperature
        ComputePhiT_CUDA(PI);
        
        GPUtoCPUEfficient_CUDA(phi,phi_gpu,Nx,Ny,offset);
        GPUtoCPUEfficient_CUDA(T,T_gpu,Nx,Ny,offset);
        
        ExchangeData(phi,Nx,Ny,offset);
        ExchangeData(T,Nx,Ny,offset);
        
        CPUtoGPUEfficient_CUDA(phi,phi_gpu,Nx,Ny,offset);
        CPUtoGPUEfficient_CUDA(T,T_gpu,Nx,Ny,offset);
        
        // Step 4: Write file
        if (TimeStep % SaveInterval == 0)
        {
            GPUtoCPU_CUDA();
            
            #ifdef TECPLOT
            sprintf(s,"Output/Field-%d-%04d.tec",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output( s, std::ios::out );
            Output.flags( std::ios::dec );
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "TITLE = Flow" << std::endl << "VARIABLES = \"X\", \"Y\", \"phi\", \"T\" " << std::endl;
            Output << "Zone T = U I = " << Ny << " J = " << Nx << std::endl;
            
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    Output << (I*Nx+i)*dx << "\t" << (J*Ny+j)*dy << "\t" << phi[i][j] << "\t" << T[i][j] << std::endl;
            
            Output.close();
            #endif
            
            #ifdef VISIT
            sprintf(s,"Output/Field-%d-%04d.vtk",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output( s, std::ios :: out );
            Output.flags( std::ios::dec );
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "# vtk DataFile Version 3.1" << std::endl;
            Output << "Dendritic growth" << std::endl;
            Output << "ASCII" << std::endl;
            Output << "DATASET STRUCTURED_GRID" << std::endl;
            Output << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl; 
            Output << "POINTS " << Nx*Ny << " FLOAT" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    Output << (I*Nx+i)*dx << "\t" << (J*Ny+j)*dy << "\t" << 0.0 << std::endl;
            
            Output << std::endl << "POINT_DATA " << Nx*Ny << std::endl;
            Output << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    Output << phi[i][j] << std::endl;
            
            Output.close();
            
            if (Rank == 0)
            {
                std::ofstream Output("Output/Field.visit", std::ios::app);
                
                if ( !Output )
                    ErrorMessage("Output file couldnot be opened!");
                
                for (int i{}; i < size; ++i)
                {
                    sprintf(s,"Field-%d-%04d.vtk",TimeStep/SaveInterval,i);
                    
                    Output << s << std::endl;
                }
                
                Output.close();
            }
            #endif
        }
        
        TimeStep++;
        
        t += dt;
        
        CheckNaN();
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    Deallocate();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (Rank == 0)
    {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        
        std::cout << std::endl << "Total time elapsed : " << duration.count() << " milliseconds." << std::endl << std::endl;
    }
    
    // Finalize the MPI environment.
    MPI_Finalize();
    
    return 0;
}
