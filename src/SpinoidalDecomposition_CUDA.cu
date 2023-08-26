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
extern double M, kappa, alpha;

extern double **C, **D, *Temporary;
extern double *C_gpu, *D_gpu;
extern double *Laplacian, *mu, *dCdx, *dCdy, *Temporary_gpu;

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
__global__ void ComputeD ( double *D_gpu, 
                           double *C_gpu, 
                           double *Laplacian, 
                           double *mu, 
                           const double dx, 
                           const double dy, 
                           const double alpha, 
                           const double kappa, 
                           const int Nx, 
                           const int Ny, 
                           const int offset )
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    int j = blockDim.y*blockIdx.y + threadIdx.y;
    
    if ((i < Nx) && (j < Ny))
    {
        Laplacian[ij(i,j)] = (C_gpu[IJ(i+1,j)] - 2.0*C_gpu[IJ(i,j)] + C_gpu[IJ(i-1,j)])/(dx*dx) 
                           + (C_gpu[IJ(i,j+1)] - 2.0*C_gpu[IJ(i,j)] + C_gpu[IJ(i,j-1)])/(dy*dy);
        
        mu[ij(i,j)] = 2.0*alpha*C_gpu[IJ(i,j)]*(1.0 - C_gpu[IJ(i,j)])*(1.0 - 2.0*C_gpu[IJ(i,j)]);
        
        D_gpu[IJ(i,j)] = mu[ij(i,j)] - kappa*Laplacian[ij(i,j)];
    }
    
    __syncthreads();
}

//__________________________________________________________________
// 
//__________________________________________________________________
__global__ void ComputeC ( double *D_gpu, 
                           double *C_gpu, 
                           double *Laplacian, 
                           const double dx, 
                           const double dy, 
                           const double dt, 
                           const double M, 
                           const int Nx, 
                           const int Ny, 
                           const int offset )
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    int j = blockDim.y*blockIdx.y + threadIdx.y;
    
    if ((i < Nx) && (j < Ny))
    {
        Laplacian[ij(i,j)] = (D_gpu[IJ(i+1,j)] - 2.0*D_gpu[IJ(i,j)] + D_gpu[IJ(i-1,j)])/(dx*dx) 
                           + (D_gpu[IJ(i,j+1)] - 2.0*D_gpu[IJ(i,j)] + D_gpu[IJ(i,j-1)])/(dy*dy);
        
        C_gpu[IJ(i,j)] += dt*M*Laplacian[ij(i,j)];
        
        if (C_gpu[IJ(i,j)] >= 0.99999) C_gpu[IJ(i,j)] = 0.99999;
        if (C_gpu[IJ(i,j)] <= 0.00001) C_gpu[IJ(i,j)] = 0.00001;
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
    
    cudaMalloc(&C_gpu,TotalSize*sizeof(double));
    cudaMalloc(&D_gpu,TotalSize*sizeof(double));
    
    cudaMalloc(&Laplacian,Nx*Ny*sizeof(double));
    cudaMalloc(&mu,Nx*Ny*sizeof(double));
    cudaMalloc(&dCdx,Nx*Ny*sizeof(double));
    cudaMalloc(&dCdy,Nx*Ny*sizeof(double));
    
    cudaMalloc(&Temporary_gpu,ExtraSize*sizeof(double));
}

//__________________________________________________________________
// 
//__________________________________________________________________
void Deallocate_CUDA ( )
{
    cudaFree(C_gpu);
    cudaFree(D_gpu);
    
    cudaFree(Laplacian);
    cudaFree(mu);
    cudaFree(dCdx);
    cudaFree(dCdy);
    
    cudaFree(Temporary_gpu);
}

//__________________________________________________________________
// 
//__________________________________________________________________
void CPUtoGPU_CUDA ( )
{
    cudaMemcpy(C_gpu, &C[-offset][-offset], TotalSize*sizeof(double), cudaMemcpyHostToDevice);
}

//__________________________________________________________________
// 
//__________________________________________________________________
void GPUtoCPU_CUDA ( )
{
    cudaMemcpy(&C[-offset][-offset], C_gpu, TotalSize*sizeof(double), cudaMemcpyDeviceToHost);
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
void ComputeD_CUDA ( )
{
    SetKernelLaunching(Nx,Ny);
    
    ComputeD <<< CUDA_GridDimension,CUDA_BlockDimension >>> (D_gpu,C_gpu,Laplacian,mu,dx,dy,alpha,kappa,Nx,Ny,offset);
    
    cudaDeviceSynchronize();
    CUDA_ErrorMessage();
}

//__________________________________________________________________
// 
//__________________________________________________________________
void ComputeC_CUDA ( )
{
    SetKernelLaunching(Nx,Ny);
    
    ComputeC <<< CUDA_GridDimension,CUDA_BlockDimension >>> (D_gpu,C_gpu,Laplacian,dx,dy,dt,M,Nx,Ny,offset);
    
    cudaDeviceSynchronize();
    CUDA_ErrorMessage();
}
