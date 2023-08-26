#pragma once

#ifndef General_CUDA_2D_H
#define General_CUDA_2D_H

#include <cuda_runtime.h>
#include <bits/stdc++.h>

//_______________________________________________________________________________
// Set constants on the CUDA-enabled GPU. The code will run on the first CUDA-enabled gpu
//_______________________________________________________________________________
extern int MaximumThreadsPerBlock;
extern int MaximumBlocksX, MaximumBlocksY, MaximumBlocksZ;

extern int MaximumThreads1DX;
extern int MaximumThreads2DX, MaximumThreads2DY;
extern int MaximumThreads3DX, MaximumThreads3DY, MaximumThreads3DZ;

extern dim3 CUDA_GridDimension;
extern dim3 CUDA_BlockDimension;

void ErrorMessageCUDA(const char* message);
void InitializeCUDA();

void SetKernelLaunching(const int d1);
void SetKernelLaunching(const int d1, const int d2);
void SetKernelLaunching(const int d1, const int d2, const int d3);

void CUDA_ErrorMessage();
#endif
