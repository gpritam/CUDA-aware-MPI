#include "General_MPI_2D.h"

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
int IJtoRank ( Int I0, 
               Int J0, 
               Int Npy )
{
    return (I0*Npy + J0);
}

//________________________________________________________________________________________________
// Constructor
//________________________________________________________________________________________________
SetMPIEnvironment::SetMPIEnvironment ( Int nx, 
                                       Int ny, 
                                       Int npx, 
                                       Int npy, 
                                       Int offset )
{
    Nx = nx;
    Ny = ny;
    Npx = npx;
    Npy = npy;
    Offset = offset;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (Rank == 0)
    {
        // Check if sufficient processors are available or not
        if (size != Npx*Npy)
        {
            std::cout << "Required number of processors : " << Npx*Npy << std::endl;
            
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        
        // Check if sufficient grid points are there or not to conduct second order computation
        if ( (Nx < 3) || (Ny < 3) )
        {
            std::cout << "Minimum number of points along each direction should be " << 3 << std::endl;
            
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    
    // Get two indices for every Rank
    // Note that, Rank = I*Npy + J;
    J = Rank % Npy;
    I = Rank / Npy;
    
    // Find neighboring ranks
    LeftID   = (I == 0       ? IJtoRank(Npx-1,J,Npy) : IJtoRank(I-1,J,Npy));
    RightID  = (I == (Npx-1) ? IJtoRank(0,J,Npy)     : IJtoRank(I+1,J,Npy));
    BottomID = (J == 0       ? IJtoRank(I,Npy-1,Npy) : IJtoRank(I,J-1,Npy));
    TopID    = (J == (Npy-1) ? IJtoRank(I,0,Npy)     : IJtoRank(I,J+1,Npy));
    
    // Variables for MPI communication
    // Note that P[i][j] = P[i*(Ny+2*Offset) + j];
    MPI_Type_vector(Nx, Offset, (Ny+2*Offset), MPI_DOUBLE, &BottomTop);
    MPI_Type_commit(&BottomTop);
    
    MPI_Type_vector(Offset, Ny, (Ny+2*Offset), MPI_DOUBLE, &LeftRight);
    MPI_Type_commit(&LeftRight);
    
    request = new (std::nothrow) MPI_Request [8];
    status  = new (std::nothrow) MPI_Status  [8];
}

//_______________________________________________________________________________
// Destructor
//_______________________________________________________________________________
SetMPIEnvironment::~SetMPIEnvironment ()
{
    delete [] request;
    delete [] status;
    
    if (Rank == 0)
        std::cout << "Destroying MPI environment..." << std::endl;
}

//_______________________________________________________________________________
// MPI communication
//_______________________________________________________________________________
void ExchangeData ( double **&ScalarField, 
                    Int Nx, 
                    Int Ny, 
                    Int Offset )
{
    int communications = 0;
    
    // Bottom face
    if (BottomID == Rank)
    {
        for (int i0{}; i0 < Nx; ++i0)
            for (int i1{}; i1 < Offset; ++i1)
                ScalarField[i0][-i1-1] = ScalarField[i0][Ny-i1-1];
    }
    else
    {
        MPI_Isend(&ScalarField[0][0],       1, BottomTop, BottomID, 0, MPI_COMM_WORLD, &request[communications]);
        MPI_Irecv(&ScalarField[0][-Offset], 1, BottomTop, BottomID, 2, MPI_COMM_WORLD, &request[communications+1]);
        
        communications += 2;
    }
    
    // Top face
    if (TopID == Rank)
    {
        for (int i0{}; i0 < Nx; ++i0)
            for (int i1{}; i1 < Offset; ++i1)
                ScalarField[i0][Ny+i1] = ScalarField[i0][i1];
    }
    else
    {
        MPI_Isend(&ScalarField[0][Ny-Offset], 1, BottomTop, TopID, 2, MPI_COMM_WORLD, &request[communications]);
        MPI_Irecv(&ScalarField[0][Ny],        1, BottomTop, TopID, 0, MPI_COMM_WORLD, &request[communications+1]);
        
        communications += 2;
    }
    
    // Left face
    if (LeftID == Rank)
    {
        for (int i0{}; i0 < Offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                ScalarField[-i0-1][i1] = ScalarField[Nx-i0-1][i1];
    }
    else
    {
        MPI_Isend(&ScalarField[0][0],       1, LeftRight, LeftID, 3, MPI_COMM_WORLD, &request[communications]);
        MPI_Irecv(&ScalarField[-Offset][0], 1, LeftRight, LeftID, 1, MPI_COMM_WORLD, &request[communications+1]);
        
        communications += 2;
    }
    
    // Right face
    if (RightID == Rank)
    {
        for (int i0{}; i0 < Offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                ScalarField[Nx+i0][i1] = ScalarField[i0][i1];
    }
    else
    {
        MPI_Isend(&ScalarField[Nx-Offset][0], 1, LeftRight, RightID, 1, MPI_COMM_WORLD, &request[communications]);
        MPI_Irecv(&ScalarField[Nx][0],        1, LeftRight, RightID, 3, MPI_COMM_WORLD, &request[communications+1]);
        
        communications += 2;
    }
    
    for (int i{}; i < communications; ++i)
        MPI_Wait(&request[i],&status[i]);
    
    MPI_Barrier(MPI_COMM_WORLD);
}
