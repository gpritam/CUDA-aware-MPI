//________________________________________________________________________________________________
//________________________________________________________________________________________________
// Parallel program for spinoidal decomposition using Finite Difference method.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 17.02.2023
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=SpinoidalDecomposition_MPI.cpp
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

const int SaveInterval = 2000;

double dt = 1.0E-3, Tf = 250.0, M = 1.0, kappa = 0.5, alpha = 1.0;

double **C, **D, *Temporary, dx, dy, Lx = Nx, Ly = Ny, t = 0.0;
double *C_gpu, *D_gpu, *Laplacian, *mu, *dCdx, *dCdy, *Temporary_gpu;

int TimeStep = 0;

char *s;

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
extern "C++" void ComputeD_CUDA ( );
extern "C++" void ComputeC_CUDA ( );

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
    Allocate(C,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(D,Nx+2*offset,Ny+2*offset,offset,offset);
    
    Allocate(Temporary,(2*Nx*offset+2*Ny*offset));
    
    Allocate(s,200);
    
    Allocate_CUDA();
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate ()
{
    Deallocate(C,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(D,Nx+2*offset,Ny+2*offset,offset,offset);
    
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
            if (C[i][j] != C[i][j])
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
                if (C[i][Ny-1-j] != C[i][Ny-1-j])
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
                if (C[i][j] != C[i][j])
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
                if (C[Nx-1-i][j] != C[Nx-1-i][j])
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
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            C[i][j] = 0.4 + 0.02*(0.5 - RandomNumber());
    
    ExchangeData(C,Nx,Ny,offset);
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
    
    Lx = Npx*Nx;
    Ly = Npy*Ny;
    
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
        
        // Step 1: Compute D
        ComputeD_CUDA();
        
        GPUtoCPUEfficient_CUDA(D,D_gpu,Nx,Ny,offset);
        ExchangeData(D,Nx,Ny,offset);
        CPUtoGPUEfficient_CUDA(D,D_gpu,Nx,Ny,offset);
        
        // Step 2: Update C
        ComputeC_CUDA();
        
        GPUtoCPUEfficient_CUDA(C,C_gpu,Nx,Ny,offset);
        ExchangeData(C,Nx,Ny,offset);
        CPUtoGPUEfficient_CUDA(C,C_gpu,Nx,Ny,offset);
        
        // Step 3: Write file
        if (TimeStep % SaveInterval == 0)
        {
            GPUtoCPU_CUDA();
            
            #ifdef TECPLOT
            sprintf(s,"Output/Field-%d-%04d.tec",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output(s, std::ios::out);
            Output.flags(std::ios::dec);
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "TITLE = Flow" << std::endl << "VARIABLES = \"X\", \"Y\", \"C\" " << std::endl;
            Output << "Zone T = U I = " << Ny << " J = " << Nx << std::endl;
            
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    Output << (I*Nx+i)*dx/Lx << "\t" << (J*Ny+j)*dy/Ly << "\t" << C[i][j] << std::endl;
            
            Output.close();
            #endif
            
            #ifdef VISIT
            sprintf(s,"Output/Field-%d-%04d.vtk",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output(s,std::ios::out);
            Output.flags(std::ios::dec);
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "# vtk DataFile Version 3.1" << std::endl;
            Output << "Spinoidal decomposition" << std::endl;
            Output << "ASCII" << std::endl;
            Output << "DATASET STRUCTURED_GRID" << std::endl;
            Output << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl; 
            Output << "POINTS " << Nx*Ny << " FLOAT" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    Output << (I*Nx+i)*dx/Lx << "\t" << (J*Ny+j)*dy/Ly << "\t" << 0.0 << std::endl;
            
            Output << std::endl << "POINT_DATA " << Nx*Ny << std::endl;
            Output << "SCALARS C float" << std::endl << "LOOKUP_TABLE default" << std::endl;
            
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    Output << C[i][j] << std::endl;
            
            Output.close();
            
            if (Rank == 0)
            {
                std::ofstream Output("Output/Field.visit",std::ios::app);
                
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
