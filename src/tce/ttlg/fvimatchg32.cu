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

#define FNAME fvimatchg32.h
#include "includes/macro.h"
#undef FNAME



void fvimatchg32CallerWrapper(int ndim, const type * A, type * B,const int size0, const int numblocks, const int numthreads
		, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s, type alpha, type beta)
{

	//	 printf("\n***%d %d %d \n", idx_ss[1],idx_ss[2], numblocks/(idx_ss[1]*idx_ss[2]));
	//dim3 param3(idx_ss[1],idx_ss[2], numblocks/(idx_ss[1]*idx_ss[2]));
	dim3 thread_blocks(numblocks/1, 1, 1);
	switch(ndim)
	{
		EXPANDDIMS(fvimg32_kernel_, thread_blocks, numthreads,0, ( A,  B, size0, lda_s,ldb_s,idx_s, alpha, beta))
		default:
		{
		}

	}
}



	extern "C" 
void  fvimatchg32_transpose_kernel(int ndim, const type *A, type *B,  const int *lda, const int *ldb, const int* params, const int * perm, type alpha, type beta) 
{ 
	// int numBlocks = computeNumBlocksCode ; 
#ifdef printd
	printf("\nA Dims: %d \t %d \t %d\t %d\t %d\n", lda[0], lda[1], lda[2], lda[3], lda[4]);
	printf("\nParams: %d \t %d \t %d\t %d\t %d\t %d\t %d\n", params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
	printf("\nB Dims: %d \t %d \t %d\t %d\t %d\n", ldb[0], ldb[1], ldb[2], ldb[3], ldb[4]);
	printf("\n Perm: %d \t %d \t %d\t %d\t %d\n", perm[0], perm[1], perm[2], perm[3], perm[4]);
#endif	


	int numBlocks = params[6];//((size[1] + 8 -1)/8) * size[2] * ((size[3] + 8 -1)/8) * size[4] ; 

	int *d_lda_s, *d_ldb_s,  *d_idx_s; 
	int lda_s[20], ldb_s[20], idx_s[20], temp[20];
	lda_s[0] = 1;
	ldb_s[0] = 1;
	int i;
	lda_s[1] = lda_s[0] * lda[0];
	ldb_s[1] = ldb_s[0] * ldb[0];
	for(i = 1; i < ndim; i++)
	{
		idx_s[i] = ldb[i];

		lda_s[i] = lda_s[i-1] * lda[i-1];
		ldb_s[i] = ldb_s[i-1] * ldb[i-1];
	}
	for(i = 1; i < ndim; i++)
	{
#ifdef printd
		printf("%d ", idx_s[i]);
#endif
		temp[i] = lda_s[perm[i]];
	}
#ifdef printd
	printf("\n");
#endif


	SAFECUDAMALLOC(&d_lda_s,ndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_ldb_s,ndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_idx_s,ndim*sizeof(int)); 
	SAFECUDAMEMCPY(d_idx_s, idx_s,ndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_lda_s, temp,ndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_ldb_s, ldb_s,ndim*sizeof(int), cudaMemcpyHostToDevice); 


#ifdef NOHTIME
#include "includes/nohtimestart.h"
#endif

	fvimatchg32CallerWrapper(ndim, A,  B, lda[0],
			numBlocks, params[2]
			, d_lda_s,d_ldb_s,d_idx_s, alpha, beta);

#ifdef NOHTIME
#include "includes/nohtimestop.h"
#endif

	cudaDeviceSynchronize();

	{cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess){
			printf("\nKernel ERROR in dCuKernel %s (line: %d)\n", cudaGetErrorString(err), __LINE__);
			//exit(-1);
		}}
	cudaFree(d_lda_s);
	cudaFree(d_ldb_s);
	cudaFree(d_idx_s);
}



