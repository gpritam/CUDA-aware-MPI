#include "General_CUDA_2D.cuh"

#define IJ(i,j) (((i)+offset)*(Ny+2*offset)+((j)+offset))
#define ij(i,j) ((i)*Ny+(j))
#define TotalSize ((Nx+2*offset)*(Ny+2*offset))
#define ExtraSize (2*Nx*offset+2*Ny*offset)

//__________________________________________________________________
// External variables
//__________________________________________________________________
extern int Nx, Ny, Npx, Npy, offset, Rank, size;
extern double Lx, Ly, dx, dy, dt, t;
extern double tau, epsilonb, kappa, delta, beta, alpha, Gamma, Teq, Theta0;

extern double **phi, **phiold, **T, **Told, **Epsilon, **Epsilon_theta, *Temporary;
extern double *phi_gpu, *phiold_gpu, *T_gpu, *Told_gpu, *Epsilon_gpu, *Epsilon_theta_gpu;
extern double *Laplacian, *term0, *term1, *term2, *term3, *m, *Temporary_gpu;

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void RowCopy ( double *A_gpu, 
                          double *Temporary_gpu, 
                          const int Nx, 
                          const int Ny, 
                          const int offset, 
                          const int j, 
                          const int Start )
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    
    if (i < Nx)
        Temporary_gpu[Start+i] = A_gpu[IJ(i,j)];
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void ColumnCopy ( double *A_gpu, 
                             double *Temporary_gpu, 
                             const int Nx, 
                             const int Ny, 
                             const int offset, 
                             const int i, 
                             const int Start )
{
    int j = blockDim.x*blockIdx.x + threadIdx.x;
    
    if (j < Ny)
        Temporary_gpu[Start+j] = A_gpu[IJ(i,j)];
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void RowWrite ( double *A_gpu, 
                           double *Temporary_gpu, 
                           const int Nx, 
                           const int Ny, 
                           const int offset, 
                           const int j, 
                           const int Start )
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    
    if (i < Nx)
        A_gpu[IJ(i,j)] = Temporary_gpu[Start+i];
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void ColumnWrite ( double *A_gpu, 
                              double *Temporary_gpu, 
                              const int Nx, 
                              const int Ny, 
                              const int offset, 
                              const int i, 
                              const int Start )
{
    int j = blockDim.x*blockIdx.x + threadIdx.x;
    
    if (j < Ny)
        A_gpu[IJ(i,j)] = Temporary_gpu[Start+j];
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void UpdateOldValues ( double *phi_gpu, 
                                  double *T_gpu, 
                                  double *phiold_gpu, 
                                  double *Told_gpu, 
                                  const int N )
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    
    if (i < N)
    {
        phiold_gpu[i] = phi_gpu[i];
        Told_gpu[i] = T_gpu[i];
    }
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void ComputeEpsilon ( double *Epsilon_gpu, 
                                 double *Epsilon_theta_gpu, 
                                 double *phiold_gpu, 
                                 double *term0, 
                                 double *term1, 
                                 double *term2, 
                                 const double dx, 
                                 const double dy, 
                                 const double epsilonb, 
                                 const double Theta0, 
                                 const double beta, 
                                 const double delta, 
                                 const int Nx, 
                                 const int Ny, 
                                 const int offset )
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    int j = blockDim.y*blockIdx.y + threadIdx.y;
    
    if ((i < Nx) && (j < Ny))
    {
        // phix
        term0[ij(i,j)] = 0.5*(phiold_gpu[IJ(i+1,j)] - phiold_gpu[IJ(i-1,j)])/dx;
        
        // phiy
        term1[ij(i,j)] = 0.5*(phiold_gpu[IJ(i,j+1)] - phiold_gpu[IJ(i,j-1)])/dy;
        
        // Theta
        term2[ij(i,j)] = atan2(term1[ij(i,j)],term0[ij(i,j)]);
        
        Epsilon_gpu[IJ(i,j)] = epsilonb + epsilonb*delta*cos(beta*(term2[ij(i,j)]-Theta0));
        Epsilon_theta_gpu[IJ(i,j)] = -epsilonb*beta*delta*sin(beta*(term2[ij(i,j)]-Theta0));
    }
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
__device__ double PhiX ( double *phiold_gpu, 
                         const double dx, 
                         const int Nx, 
                         const int Ny, 
                         const int offset, 
                         const int i, 
                         const int j )
{
    return 0.5*(phiold_gpu[IJ(i+1,j)]-phiold_gpu[IJ(i-1,j)])/dx;
}

//__________________________________________________________________
// 
//__________________________________________________________________
__device__ double PhiY ( double *phiold_gpu, 
                         const double dy, 
                         const int Nx, 
                         const int Ny, 
                         const int offset, 
                         const int i, 
                         const int j )
{
    return 0.5*(phiold_gpu[IJ(i,j+1)]-phiold_gpu[IJ(i,j-1)])/dy;
}

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void ComputePhiT ( double *phi_gpu, 
                              double *T_gpu, 
                              double *Epsilon_gpu, 
                              double *Epsilon_theta_gpu, 
                              double *phiold_gpu, 
                              double *Told_gpu, 
                              double *Laplacian, 
                              double *term0, 
                              double *term1, 
                              double *term2, 
                              double *term3, 
                              double *m, 
                              const double dt, 
                              const double dx, 
                              const double dy, 
                              const double alpha, 
                              const double kappa, 
                              const double Teq, 
                              const double Gamma, 
                              const double tau, 
                              const double PI, 
                              const int Nx, 
                              const int Ny, 
                              const int offset )
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    int j = blockDim.y*blockIdx.y + threadIdx.y;
    
    if ((i < Nx) && (j < Ny))
    {
        // Update phase-field parameter
        Laplacian[ij(i,j)] = (phiold_gpu[IJ(i+1,j)] - 2.0*phiold_gpu[IJ(i,j)] + phiold_gpu[IJ(i-1,j)])/(dx*dx) + (phiold_gpu[IJ(i,j+1)] - 2.0*phiold_gpu[IJ(i,j)] + phiold_gpu[IJ(i,j-1)])/(dy*dy);
        
        m[ij(i,j)] = alpha*atan(Gamma*(Teq-Told_gpu[IJ(i,j)]))/PI;
        
        term0[ij(i,j)] = 0.5*(Epsilon_gpu[IJ(i,j+1)]*Epsilon_theta_gpu[IJ(i,j+1)]*PhiX(phiold_gpu,dx,Nx,Ny,offset,i,j+1) - Epsilon_gpu[IJ(i,j-1)]*Epsilon_theta_gpu[IJ(i,j-1)]*PhiX(phiold_gpu,dx,Nx,Ny,offset,i,j-1))/dy;
        
        term1[ij(i,j)] = 0.5*(Epsilon_gpu[IJ(i+1,j)]*Epsilon_theta_gpu[IJ(i+1,j)]*PhiY(phiold_gpu,dy,Nx,Ny,offset,i+1,j) - Epsilon_gpu[IJ(i-1,j)]*Epsilon_theta_gpu[IJ(i-1,j)]*PhiY(phiold_gpu,dy,Nx,Ny,offset,i-1,j))/dx;
        
        term2[ij(i,j)] = Epsilon_gpu[IJ(i,j)]*Epsilon_gpu[IJ(i,j)]*Laplacian[ij(i,j)];
        
        term3[ij(i,j)] = phiold_gpu[IJ(i,j)]*(1.0-phiold_gpu[IJ(i,j)])*(phiold_gpu[IJ(i,j)]-0.5+m[ij(i,j)]);
        
        phi_gpu[IJ(i,j)] += dt*(term0[ij(i,j)] - term1[ij(i,j)] + term2[ij(i,j)] + term3[ij(i,j)])/tau;
        
        // Update temperature
        Laplacian[ij(i,j)] = (Told_gpu[IJ(i+1,j)] - 2.0*Told_gpu[IJ(i,j)] + Told_gpu[IJ(i-1,j)])/(dx*dx) + (Told_gpu[IJ(i,j+1)] - 2.0*Told_gpu[IJ(i,j)] + Told_gpu[IJ(i,j-1)])/(dy*dy);
        
        T_gpu[IJ(i,j)] += dt*Laplacian[ij(i,j)] + kappa*(phi_gpu[IJ(i,j)]-phiold_gpu[IJ(i,j)]);
    }
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
void Allocate_CUDA ( )
{
    InitializeCUDA();
    
    int deviceCount;
    
    cudaGetDeviceCount(&deviceCount);
    
    int device_id = Rank % deviceCount;
    
    cudaSetDevice(device_id);
    
    cudaMalloc(&phi_gpu,TotalSize*sizeof(double));
    cudaMalloc(&phiold_gpu,TotalSize*sizeof(double));
    cudaMalloc(&T_gpu,TotalSize*sizeof(double));
    cudaMalloc(&Told_gpu,TotalSize*sizeof(double));
    cudaMalloc(&Epsilon_gpu,TotalSize*sizeof(double));
    cudaMalloc(&Epsilon_theta_gpu,TotalSize*sizeof(double));
    
    cudaMalloc(&Laplacian,Nx*Ny*sizeof(double));
    cudaMalloc(&term0,Nx*Ny*sizeof(double));
    cudaMalloc(&term1,Nx*Ny*sizeof(double));
    cudaMalloc(&term2,Nx*Ny*sizeof(double));
    cudaMalloc(&term3,Nx*Ny*sizeof(double));
    cudaMalloc(&m,Nx*Ny*sizeof(double));
    
    cudaMalloc(&Temporary_gpu,ExtraSize*sizeof(double));
}

//__________________________________________________________________
// 
//__________________________________________________________________
void Deallocate_CUDA ( )
{
    cudaFree(phi_gpu);
    cudaFree(phiold_gpu);
    cudaFree(T_gpu);
    cudaFree(Told_gpu);
    cudaFree(Epsilon_gpu);
    cudaFree(Epsilon_theta_gpu);
    
    cudaFree(Laplacian);
    cudaFree(term0);
    cudaFree(term1);
    cudaFree(term2);
    cudaFree(term3);
    cudaFree(m);
    
    cudaFree(Temporary_gpu);
}

//__________________________________________________________________
// 
//__________________________________________________________________
void CPUtoGPU_CUDA ( )
{
    cudaMemcpy(phi_gpu, &phi[-offset][-offset], TotalSize*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(T_gpu, &T[-offset][-offset], TotalSize*sizeof(double), cudaMemcpyHostToDevice);
}

//__________________________________________________________________
// 
//__________________________________________________________________
void GPUtoCPU_CUDA ( )
{
    cudaMemcpy(&phi[-offset][-offset], phi_gpu, TotalSize*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&T[-offset][-offset], T_gpu, TotalSize*sizeof(double), cudaMemcpyDeviceToHost);
}

//__________________________________________________________________
// 
//__________________________________________________________________
void GPUtoCPUEfficient_CUDA ( double **&A, 
                              double *A_gpu, 
                              const int Nx, 
                              const int Ny, 
                              const int offset )
{
    int Start;
    
    SetKernelLaunching(Nx);
    
    // Bottom
    Start = 0;
    
    for (int j = 0; j < offset; ++j)
    {
        RowCopy <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,j,(Start+j*Nx));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
    
    // Top
    Start = Nx*offset;
    
    for (int j = 0; j < offset; ++j)
    {
        RowCopy <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,(Ny-1-j),(Start+j*Nx));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
    
    SetKernelLaunching(Ny);
    
    // Left
    Start = 2*Nx*offset;
    
    for (int i = 0; i < offset; ++i)
    {
        ColumnCopy <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,i,(Start+i*Ny));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
    
    // Right
    Start = 2*Nx*offset + Ny*offset;
    
    for (int i = 0; i < offset; ++i)
    {
        ColumnCopy <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,(Nx-1-i),(Start+i*Ny));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
    
    cudaMemcpy(Temporary, Temporary_gpu, ExtraSize*sizeof(double), cudaMemcpyDeviceToHost);
    
    // Bottom
    Start = 0;
    
    for (int j = 0; j < offset; ++j)
        for (int i = 0; i < Nx; ++i)
            A[i][j] = Temporary[Start+j*Nx+i];
    
    // Top
    Start = Nx*offset;
    
    for (int j = 0; j < offset; ++j)
        for (int i = 0; i < Nx; ++i)
            A[i][Ny-1-j] = Temporary[Start+j*Nx+i];
    
    // Left
    Start = 2*Nx*offset;
    
    for (int i = 0; i < offset; ++i)
        for (int j = 0; j < Ny; ++j)
            A[i][j] = Temporary[Start+i*Ny+j];
    
    // Right
    Start = 2*Nx*offset + Ny*offset;
    
    for (int i = 0; i < offset; ++i)
        for (int j = 0; j < Ny; ++j)
            A[Nx-1-i][j] = Temporary[Start+i*Ny+j];
}

//__________________________________________________________________
// 
//__________________________________________________________________
void CPUtoGPUEfficient_CUDA ( double **A, 
                              double *A_gpu, 
                              const int Nx, 
                              const int Ny, 
                              const int offset )
{
    int Start;
    
    // Bottom
    Start = 0;
    
    for (int j = 0; j < offset; ++j)
        for (int i = 0; i < Nx; ++i)
            Temporary[Start+j*Nx+i] = A[i][-j-1];
    
    // Top
    Start = Nx*offset;
    
    for (int j = 0; j < offset; ++j)
        for (int i = 0; i < Nx; ++i)
            Temporary[Start+j*Nx+i] = A[i][Ny+j];
    
    // Left
    Start = 2*Nx*offset;
    
    for (int i = 0; i < offset; ++i)
        for (int j = 0; j < Ny; ++j)
            Temporary[Start+i*Ny+j] = A[-i-1][j];
    
    // Right
    Start = 2*Nx*offset + Ny*offset;
    
    for (int i = 0; i < offset; ++i)
        for (int j = 0; j < Ny; ++j)
            Temporary[Start+i*Ny+j] = A[Nx+i][j];
    
    cudaMemcpy(Temporary_gpu, Temporary, ExtraSize*sizeof(double), cudaMemcpyHostToDevice);
    
    SetKernelLaunching(Nx);
    
    // Bottom
    Start = 0;
    
    for (int j = 0; j < offset; ++j)
    {
        RowWrite <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,(-j-1),(Start+j*Nx));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
    
    // Top
    Start = Nx*offset;
    
    for (int j = 0; j < offset; ++j)
    {
        RowWrite <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,(Ny+j),(Start+j*Nx));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
    
    SetKernelLaunching(Ny);
    
    // Left
    Start = 2*Nx*offset;
    
    for (int i = 0; i < offset; ++i)
    {
        ColumnWrite <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,(-i-1),(Start+i*Ny));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
    
    // Right
    Start = 2*Nx*offset + Ny*offset;
    
    for (int i = 0; i < offset; ++i)
    {
        ColumnWrite <<< CUDA_GridDimension,CUDA_BlockDimension >>> (A_gpu,Temporary_gpu,Nx,Ny,offset,(Nx+i),(Start+i*Ny));
        
        cudaDeviceSynchronize();
        CUDA_ErrorMessage();
    }
}

//__________________________________________________________________
// 
//__________________________________________________________________
void UpdateOldValues_CUDA ( )
{
    SetKernelLaunching(TotalSize);
    
    UpdateOldValues <<< CUDA_GridDimension,CUDA_BlockDimension >>> (phi_gpu,T_gpu,phiold_gpu,Told_gpu,TotalSize);
    
    cudaDeviceSynchronize();
    CUDA_ErrorMessage();
}

//__________________________________________________________________
// 
//__________________________________________________________________
void ComputeEpsilon_CUDA ( )
{
    SetKernelLaunching(Nx,Ny);
    
    ComputeEpsilon <<< CUDA_GridDimension,CUDA_BlockDimension >>> (Epsilon_gpu,Epsilon_theta_gpu,phiold_gpu,term0,term1,term2,dx,dy, epsilonb,Theta0,beta,delta,Nx,Ny,offset);
    
    cudaDeviceSynchronize();
    CUDA_ErrorMessage();
}

//__________________________________________________________________
// 
//__________________________________________________________________
void ComputePhiT_CUDA ( const double PI )
{
    SetKernelLaunching(Nx,Ny);
    
    ComputePhiT <<< CUDA_GridDimension,CUDA_BlockDimension >>> (phi_gpu,T_gpu,Epsilon_gpu,Epsilon_theta_gpu,phiold_gpu,Told_gpu,Laplacian, 
    term0,term1,term2,term3,m,dt,dx,dy,alpha,kappa,Teq,Gamma,tau,PI,Nx,Ny,offset);
    
    cudaDeviceSynchronize();
    CUDA_ErrorMessage();
}
