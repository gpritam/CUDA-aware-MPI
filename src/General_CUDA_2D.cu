#include "General_CUDA_2D.cuh"

//_______________________________________________________________________________
// Set constants on the CUDA-enabled GPU. The code will run on the first CUDA-enabled gpu
//_______________________________________________________________________________
int MaximumThreadsPerBlock;
int MaximumBlocksX, MaximumBlocksY, MaximumBlocksZ;

int MaximumThreads1DX;
int MaximumThreads2DX, MaximumThreads2DY;
int MaximumThreads3DX, MaximumThreads3DY, MaximumThreads3DZ;

dim3 CUDA_GridDimension;
dim3 CUDA_BlockDimension;

//_______________________________________________________________________________
// Print an error message and exit the program
//_______________________________________________________________________________
void ErrorMessageCUDA ( const char* message )
{
    std::cout << message << std::endl;
    
    exit(1);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void InitializeCUDA ()
{
    int DeviceNumber;
    
    cudaGetDeviceCount(&DeviceNumber);
    
    if (DeviceNumber == 0)
        ErrorMessageCUDA("No CUDA capable devices detected.");
    else
    {
        cudaDeviceProp DeviceProperty;
        
        cudaGetDeviceProperties(&DeviceProperty, 0);
        
        MaximumThreadsPerBlock = DeviceProperty.maxThreadsPerBlock;
        
        MaximumBlocksX = DeviceProperty.maxGridSize[0];
        MaximumBlocksY = DeviceProperty.maxGridSize[1];
        MaximumBlocksZ = DeviceProperty.maxGridSize[2];
        
        MaximumThreads1DX = MaximumThreadsPerBlock;
        
        MaximumThreads2DX = floor(sqrt((double)(MaximumThreadsPerBlock)));
        MaximumThreads2DY = MaximumThreads2DX;
        
        MaximumThreads3DX = floor(cbrt((double)(MaximumThreadsPerBlock)));;
        MaximumThreads3DY = MaximumThreads3DX;
        MaximumThreads3DZ = MaximumThreads3DX;
    }
}

//_______________________________________________________________________________
// Number of Blocks in one dimension
//_______________________________________________________________________________
void SetKernelLaunching ( const int d1 )
{
    CUDA_BlockDimension.x = MaximumThreads1DX;
    CUDA_BlockDimension.y = 1;
    CUDA_BlockDimension.z = 1;
    
    CUDA_GridDimension.x = ceil(d1/((double)(CUDA_BlockDimension.x)));
    CUDA_GridDimension.y = 1;
    CUDA_GridDimension.z = 1;
    
    if (CUDA_GridDimension.x > MaximumBlocksX)
        ErrorMessageCUDA("Number of blocks along x direction is more than maximum limit!");
}

//_______________________________________________________________________________
// Number of Blocks in two dimensions
//_______________________________________________________________________________
void SetKernelLaunching ( const int d1, 
                          const int d2 )
{
    CUDA_BlockDimension.x = MaximumThreads2DX;
    CUDA_BlockDimension.y = MaximumThreads2DY;
    CUDA_BlockDimension.z = 1;
    
    CUDA_GridDimension.x = ceil(d1/((double)(CUDA_BlockDimension.x)));
    CUDA_GridDimension.y = ceil(d2/((double)(CUDA_BlockDimension.y)));
    CUDA_GridDimension.z = 1;
    
    if (CUDA_GridDimension.x > MaximumBlocksX)
        ErrorMessageCUDA("Number of blocks along x direction is more than maximum limit!");
    
    if (CUDA_GridDimension.y > MaximumBlocksY)
        ErrorMessageCUDA("Number of blocks along y direction is more than maximum limit!");
}

//_______________________________________________________________________________
// Number of Blocks in three dimensions
//_______________________________________________________________________________
void SetKernelLaunching ( const int d1, 
                          const int d2, 
                          const int d3 )
{
    CUDA_BlockDimension.x = MaximumThreads3DX;
    CUDA_BlockDimension.y = MaximumThreads3DY;
    CUDA_BlockDimension.z = MaximumThreads3DZ;
    
    CUDA_GridDimension.x = ceil(d1/((double)(CUDA_BlockDimension.x)));
    CUDA_GridDimension.y = ceil(d2/((double)(CUDA_BlockDimension.y)));
    CUDA_GridDimension.z = ceil(d3/((double)(CUDA_BlockDimension.z)));
    
    if (CUDA_GridDimension.x > MaximumBlocksX)
        ErrorMessageCUDA("Number of blocks along x direction is more than maximum limit!");
    
    if (CUDA_GridDimension.y > MaximumBlocksY)
        ErrorMessageCUDA("Number of blocks along y direction is more than maximum limit!");
    
    if (CUDA_GridDimension.z > MaximumBlocksZ)
        ErrorMessageCUDA("Number of blocks along z direction is more than maximum limit!");
}

//_______________________________________________________________________________
// Check for error in CUDA
//_______________________________________________________________________________
void CUDA_ErrorMessage ( )
{
    cudaError_t error = cudaGetLastError();
    
    if ( error != cudaSuccess )
    {
        std::cout << "CUDA Error : " << cudaGetErrorString(error) << std::endl;
        
        exit(1);
    }
}
