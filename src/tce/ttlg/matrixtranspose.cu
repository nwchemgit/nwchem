#define BLOCK_DIM 16
#include <stdio.h>
#include <cublas_v2.h>
#include <cstdlib>
#include <cuda.h>
#include <omp.h>
#define type double


__global__ void transpose(type *odata, type *idata, int width, int height)
{
	__shared__ type block[BLOCK_DIM][BLOCK_DIM+1];

	// read the matrix tile into shared memory

	unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
	if((xIndex < width) && (yIndex < height))
	{
		unsigned int index_in = yIndex * width + xIndex;
		block[threadIdx.y][threadIdx.x] = idata[index_in];
	}
	__syncthreads();

	// write the transposed matrix tile to global memory
	xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
	yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;

	if((xIndex < height) && (yIndex < width))
	{
		unsigned int index_out = yIndex * height + xIndex;
		odata[index_out] = block[threadIdx.x][threadIdx.y];
	}

}
extern "C" void matrix_transpose(type *odata, type *idata, int width, int height)
{
#ifdef printd
	printf("In matrix transpose..");
#endif

	int ndim = 2;
	int lda[2];
	lda[0] = width, lda[1] = height;
	cublasHandle_t handle;
	cublasCreate(&handle);
	cublasOperation_t aflag = CUBLAS_OP_T, bflag = CUBLAS_OP_T;
	double alpha = 1, beta = 0;
	//cublasDgeam(handle,aflag, bflag,width, height, &alpha, idata, width, &beta, idata, width, odata, height);
#ifdef NOHTIME
#include "includes/nohtimestart.h"
#endif

	cublasDgeam(handle,aflag, bflag,height, width, &alpha, idata, width, &beta, idata, width, odata, height);
#ifdef NOHTIME
#include "includes/nohtimestop.h"
#endif

	return;

}
