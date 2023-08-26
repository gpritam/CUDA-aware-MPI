#pragma once

#ifndef General_MPI_2D_H
#define General_MPI_2D_H

#include <mpi.h>
#include "General.h"

inline MPI_Request *request;
inline MPI_Status *status;

// A pair of unique indices for every Rank
inline int I, J;

// Neighbors' ranks
inline int LeftID, RightID, TopID, BottomID;

// Variables for MPI communication
inline MPI_Datatype BottomTop, LeftRight;

inline int Rank, size;

int IJtoRank(Int I0, Int J0, Int Npy);
void ExchangeData(double **&ScalarField, Int Nx, Int Ny, Int Offset);

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
class SetMPIEnvironment
{
    private:
    
    int Nx, Ny, Npx, Npy, Offset;
    
    public:
    
    SetMPIEnvironment(Int, Int, Int, Int, Int);
    ~SetMPIEnvironment();
};
#endif
