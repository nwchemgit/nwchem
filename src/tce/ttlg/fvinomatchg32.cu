#include <cuda_runtime.h> 
#include <cuComplex.h> 
#include <complex.h> 
#include <stdio.h> 
#include <omp.h>
#define type double

#define STR1(X) #X
#define STR(X) STR1(X) 
#define STRINGIFY(X,Y) X ## Y
#define CON(X,Y) STRINGIFY(X,Y)
#define KDir kernels


#include "includes/ourmacros.h"

extern __shared__ type tile[];
__device__ __forceinline__ void fvinomatchg32_main(const type * __restrict__ Atmp,  type * __restrict__ Btmp,const  int lda1,const  int ldb1) 
{ 
#define TILE_DIM_X  32 
#define TILE_DIM_Y  32 
#define THREADS_PER_ROW  8
#define SHM2 33 
	int id = threadIdx.x;
	int rowId = id%32;
	int colId = id/TILE_DIM_X;
#ifdef printd
	if(blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 && threadIdx.x == 0)
		printf("lda1 = %d, ldb1 = %d\n", lda1, ldb1);
#endif
	for(int j = colId; j < TILE_DIM_Y; j+= THREADS_PER_ROW)
	{
		tile[j * SHM2 + rowId] = Atmp[rowId + j * lda1];
	}
	__syncthreads();
	for(int j = colId; j < TILE_DIM_Y; j+= THREADS_PER_ROW)
	{
		//Btmp[j * ldb1 + rowId] =tile[rowId * SHM2 + j];
		Btmp[j * ldb1 + rowId] =tile[rowId * SHM2 + j];
	}
#undef SHM2
} 

__device__ __forceinline__ 
void fvinomatchg32_rem(const type * __restrict__ Atmp,  type * __restrict__ Btmp, const int lda1, const int ldb1, const  int remainderx,const  int remaindery) 
{
#define TILE_DIM_X  32 
#define TILE_DIM_Y  32 
#define THREADS_PER_ROW  8
	int id = threadIdx.x;
#define SHM2 33 
	int rowId = id%32;
	int colId = id/TILE_DIM_X;

	int limit = remaindery > TILE_DIM_Y ? TILE_DIM_Y :  remaindery;	
	if(rowId < remainderx)
		for(int j = colId; j < limit; j+= THREADS_PER_ROW)
		{
			tile[j * SHM2 + rowId] = Atmp[rowId + j * lda1];
		}
	__syncthreads();
	if(rowId >= remaindery)
		return;
	limit = remainderx > TILE_DIM_Y ? TILE_DIM_Y :  remainderx;	
	for(int j = colId; j < limit; j+= THREADS_PER_ROW)
	{
		Btmp[j * ldb1 + rowId] =tile[rowId * SHM2 + j];
	}
#undef SHM2
}   
__device__ __forceinline__ void fvinomatchg32_main_coars(const type * __restrict__ Atmp,  type * __restrict__ Btmp,const  int lda1,const  int ldb1, const int acoars, const int bcoars, const int size) 
{ 
#define TILE_DIM_X  32 
#define TILE_DIM_Y  32 
#define THREADS_PER_ROW  8
#define SHM2 33 
	int id = threadIdx.x;
	int rowId = id%32;
	int colId = id/TILE_DIM_X;
#ifdef printd
	if(blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 && threadIdx.x == 0)
		printf("lda1 = %d, ldb1 = %d\n", lda1, ldb1);
#endif
	for(int i = 0; i < size; i++)
	{
		for(int j = colId; j < TILE_DIM_Y; j+= THREADS_PER_ROW)
		{
			tile[j * SHM2 + rowId] = Atmp[i*acoars+rowId + j * lda1];
		}
		__syncthreads();
		for(int j = colId; j < TILE_DIM_Y; j+= THREADS_PER_ROW)
		{
			Btmp[i*bcoars +j * ldb1 + rowId] =tile[rowId * SHM2 + j];
		}
		__syncthreads();
	}
#undef SHM2
} 

	__device__ __forceinline__ 
void fvinomatchg32_rem_coars(const type * __restrict__ Atmp,  type * __restrict__ Btmp, const int lda1, const int ldb1, const  int remainderx,const  int remaindery, const int acoars, const int bcoars, const int size) 
{
#define TILE_DIM_X  32 
#define TILE_DIM_Y  32 
#define THREADS_PER_ROW  8
	int id = threadIdx.x;
#define SHM2 33 
	int rowId = id%32;
	int colId = id/TILE_DIM_X;
	
#ifdef printd
//if(threadIdx.x == 0 && blockIdx.x%5 == 0)
//printf("lda1 = %d, ldb1 = %d, remainderx = %d, remaindery = %d, acoars = %d, bcoars = %d, size = %d\n", lda1, ldb1, remainderx, remaindery, acoars, bcoars, size);
#endif
	for(int i = 0; i < size; i++)
	{
		if(rowId < remainderx)
			for(int j = colId; j < remaindery; j+= THREADS_PER_ROW)
			{
				tile[j * SHM2 + rowId] = Atmp[i*acoars+rowId + j * lda1];
			}
		__syncthreads();
		if(rowId < remaindery)
		for(int j = colId; j < remainderx; j+= THREADS_PER_ROW)
		{
			Btmp[i*bcoars+j * ldb1 + rowId] =tile[rowId * SHM2 + j];
		}
		__syncthreads();
	}
#undef SHM2
}   
#define FNAME fvinomatchg32.h
#include "includes/macro.h"
#undef FNAME
#define FNAME fvinomatchg32_coars.h
#include "includes/macro.h"
#undef FNAME



void fvinomatchg32_CallerWrapper(int ndim, const type * __restrict__ A, type * B,  const int numblocks, const int numthreads, const int shmsize
		, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
		, const int coarsa, const int coarsb, const int lda_kernel1, const int ldb_kernel1,  const int inputrem, const int outputrem, const int remainderx, const int remaindery, const int size)
{
		#ifdef printd
			printf("ndim  = %d, size = %d, lda1 = %d, ldb1 = %d, irem = %d, orem = %d\n", ndim, size, lda_kernel1, ldb_kernel1, inputrem, outputrem);
		#endif

	if(size > 0)
	{
		#ifdef printd
			printf("Coarsening... No. of blocks = %d\n", numblocks/size);
		#endif
		dim3 thread_blocks(numblocks/size, 1, 1);
		switch(ndim)
		{
			EXPANDDIMS(fvinomatchg32_coars_kernel_, thread_blocks, numthreads, shmsize, (A,  B, lda_s,ldb_s, idx_s, coarsa,coarsb, inputrem, outputrem, lda_kernel1, ldb_kernel1, remainderx, remaindery, size))
			default:
			{
			}
		}
	}
	else
	{
		dim3 thread_blocks(numblocks, 1, 1);
		switch(ndim)
		{
			EXPANDDIMS(fvinomatchg32_kernel_, thread_blocks, numthreads, shmsize, (A, B, lda_s,ldb_s, idx_s, inputrem, outputrem, lda_kernel1, ldb_kernel1, remainderx, remaindery))
			default:
			{
			}
		}
	}
}
void swap(int array[], int ind1, int ind2);

int  cancoarsen(int *lda, int newndim);
	extern "C" 
void  fvinomatchg32_transpose_kernel(int ndim, const type *A, type *B,  const int *lda, const int *ldb, const int* params, const int * rperm) 
{ 
	// int numBlocks = computeNumBlocksCode ; 
#ifdef printd
	printf("\nA Dims: %d \t %d \t %d\t %d\t %d\n", lda[0], lda[1], lda[2], lda[3], lda[4]);
	printf("\nParams: %d \t %d \t %d\t %d\t %d\t %d\t %d \t %d \t %d \t %d \t %d\n", params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10]);
	printf("\nB Dims: %d \t %d \t %d\t %d\t %d\n", ldb[0], ldb[1], ldb[2], ldb[3], ldb[4]);
	printf("\nR perm: %d \t %d \t %d\t %d\t %d\n", rperm[0], rperm[1], rperm[2], rperm[3], rperm[4]);
#endif	

	//#ifndef orig1

	int numBlocks = params[6];//((size[1] + 8 -1)/8) * size[2] * ((size[3] + 8 -1)/8) * size[4] ; 

	//	printf("\n%d \t %d \t %d\t %d\t %d\n", ldb[0], ldb[1], ldb[2], ldb[3], ldb[4]);
	const int blockA = params[0];
	const int blockB = params[0];
	int irem, orem;	
	//const int shm = params[5] * params[10]; 
	int *d_lda_s, *d_ldb_s,  *d_idx_s; 
	const int remainder1 = lda[0] % 32;
	const int remainder2 = lda[params[7]] % 32;
	if(remainder1 == 0) irem = lda[0];
	else irem = (lda[0] - remainder1)/blockA;
	if(remainder2 == 0) orem = ldb[0];
	else orem = (ldb[0] - remainder2)/blockB;

	int lda_s[20], ldb_s[20], idx_s[20], temp[20];
	lda_s[0] = 1;
	ldb_s[0] = 1;

	int i;
	idx_s[0] = (lda[0] + blockA - 1) / blockA;
	for(i = 1; i < ndim; i++)
	{
		if( i == params[7])
		{
			idx_s[i] = (lda[i] + blockB - 1)/blockB;
		}
		else
		{
			idx_s[i] = lda[i];
		}
		lda_s[i] = lda_s[i-1] * lda[i-1];
		ldb_s[i] = ldb_s[i-1] * ldb[i-1];
	}
	for(i = 0; i < ndim; i++)
	{
		temp[i] = ldb_s[rperm[i]];
	}

	const int lda_kernel1 =  lda_s[params[7]]; //lda,0
	const  int ldb_kernel1 = ldb_s[params[8]]; //ldb,0
	lda_s[0] *= params[0];
        temp[0] *= params[0];
        lda_s[params[7]] *= params[0];
        temp[params[7]] *= params[0];

	if(params[7] != 1)
	{
		swap(idx_s, 1, params[7]);
		swap(lda_s, 1, params[7]);
		swap(temp, 1, params[7]);
	}
	int newndim = ndim;

	int acoars = -1, bcoars = -1, size = -1;
#ifndef NOCOARSEN
	int cd = cancoarsen(idx_s+2, ndim-2);
	if(cd >= 0)
	{
#ifdef printd
		printf("cd = %d\n", cd);
#endif
		acoars = lda_s[2+cd];
		bcoars = temp[2+cd];
		size = idx_s[2+cd];
		for(int j = cd+1; j < newndim; j++)
		{
			idx_s[2+j-1] = idx_s[2+j];
			lda_s[2+j-1] = lda_s[2+j];
			temp[2+j-1] = temp[2+j];

		}
		newndim--;
	}
#endif

	SAFECUDAMALLOC(&d_lda_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_ldb_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_idx_s,newndim*sizeof(int)); 
	SAFECUDAMEMCPY(d_idx_s, idx_s,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_lda_s, lda_s,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_ldb_s, temp,newndim*sizeof(int), cudaMemcpyHostToDevice); 

#ifdef MODEL
printf("\t%d\t%d\t", lda_kernel1, ldb_kernel1);
printf("\t%d\t%d\t%d\t%d\t", lda[0]/32,lda[0]%32, ldb[0]/32,ldb[0]%32 );
printf("\t%lf\t", ((lda[0]/32) * (ldb[0]/32) + (double)(lda[0]/32) * (ldb[0]%32) /32+ (double)(lda[0]%32) * (ldb[0]/32) /32 + (double)(lda[0]%32) * (ldb[0]%32) /(32*32) )/ (int)(((lda[0]+31)/32) * ((ldb[0]+31)/32)));
//printf("\t%lf\t%d %d \t", ((lda[0]/32) * (ldb[0]/32) + (double)(lda[0]/32) * (ldb[0]%32) /32+ (double)(lda[0]%32) * (ldb[0]/32) /32 + (double)(lda[0]%32) * (ldb[0]%32) /(32*32) ), (int)((lda[0]+31)/32) * ((ldb[0]+31)/32), ldb[0]);
#endif

#ifdef NOHTIME
#include "includes/nohtimestart.h"
#endif
	fvinomatchg32_CallerWrapper(newndim, A,  B,
			numBlocks, params[2],
			params[10] * params[5]*sizeof(type)
			, d_lda_s,d_ldb_s,d_idx_s
			, acoars, bcoars
			, lda_kernel1, ldb_kernel1,irem, orem, remainder1, remainder2, size);


#ifdef NOHTIME
#include "includes/nohtimestop.h"
#endif
	cudaDeviceSynchronize();

	{cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess){
			printf("\nKernel ERROR in fvi_nomatch_g32: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);
			//exit(-1);
		}}
	cudaFree(d_lda_s);
	cudaFree(d_ldb_s);
	cudaFree(d_idx_s);
}



