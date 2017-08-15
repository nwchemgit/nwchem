#include <cuda_runtime.h> 
#include <cuComplex.h> 
#include <complex.h> 
#include <stdio.h> 
#include<omp.h>
#define type double

#define STR1(X) #X
#define STR(X) STR1(X) 
#define STRINGIFY(X,Y) X ## Y
#define CON(X,Y) STRINGIFY(X,Y)
#define KDir kernels

#include "includes/ourmacros.h"


/*
   __device__ __inline__ double ld_gbl_cg(const double *addr) {
   double return_value;
   asm("ld.global.cg.f64 %0, [%1];" : "=d"(return_value) : "l"(addr));
   return return_value;
   }*/

extern __shared__ type tile[];

__device__ __forceinline__ void fvimatchl32_main(const  type *   Atmp, type *  Btmp,const  int lda1,const  int ldb1,const int size0, const int  plain, const int sbp, const unsigned  short * __restrict__ offset, type alpha, type beta) 
{ 

	const int id = threadIdx.x; 
	const int i2 = id / 32; 
	const int rest = id %32; 
	const type *Adisp = Atmp +i2 * lda1; 
	const int tile_disp = i2 * (sbp);
	int regs[8];int j = rest;
#pragma unroll 8
	for(int i = 0; i < 8; i++)
	{
		regs[i] = offset[j];
		j+=32;
		if(j >= plain)
			break;
	}
#pragma unroll
	for( int j = rest; j < plain; j+=32){
		//tile[tile_disp+j] =  ld_gbl_cg(&Adisp[j]); 
		tile[tile_disp+j] =  Adisp[j]; 
	}
	__syncthreads();
	type *Bdisp = &Btmp[i2 * ldb1]; 
	const  int tile_disp1 = i2 * size0;
	j = rest;
#pragma unroll 8
	for(int i = 0; i < 8; i++)
	{
		Bdisp[j] = beta*Bdisp[j] + alpha* tile[regs[i]+tile_disp1]; 
		j+=32;
		if(j >= plain)
			break;
	}
} 

	__device__ __forceinline__ 
void fvimatchl32_rem(const type *  Atmp, type *  Btmp, const int lda1, const int ldb1,const  int remainderx,const  int remaindery,const  int size0, const int plain, const int sbp, const int ilimit, const int olimit, const unsigned short int* __restrict__ offset, type alpha, type beta) 
{

	int id = threadIdx.x; 
	int i2 = id / 32; 
	int rest = id %32; 
	if(i2 < remaindery){ 
		const type *Adisp = &Atmp[i2 * lda1]; 
		const int tile_disp = i2 * (sbp); 
#pragma unroll		
		for( int  j = rest; j < ilimit; j+=32){ 
			tile[tile_disp+j] =  Adisp[j]; 
		} 
	} 
	__syncthreads(); 
	if(i2 >= remainderx)
		return;
	type *Bdisp = &Btmp[i2 * ldb1]; 
	int regs[8];int j = rest;
#pragma unroll 8
	for(int i = 0; i < 8; i++)
	{
		if(j >= olimit)
			break;
		regs[i] = offset[j];
		j+=32;
	}
	const  int tile_disp1 = i2 * size0; 
	j = rest;
#pragma unroll 8
	for(int i = 0; i < 8; i++)
	{
		if(j >= olimit)
			break;
		Bdisp[j] = beta*Bdisp[j]+ alpha* tile[regs[i]+tile_disp1]; 
		j+=32;
	}
}   
__device__ __forceinline__ void fvimatchl32_main_coars(const  type *   Atmp, type *  Btmp,const  int lda1,const  int ldb1,const int size0, const int  plain, const int sbp, const unsigned  short * __restrict__ offset, const int acoars, const int bcoars, const int size, type alpha, type beta) 
{ 

	const int id = threadIdx.x; 
	const int i2 = id / 32; 
	const int rest = id %32; 
	const int tile_disp = i2 * (sbp);
	int regs[8];
	for(int c = 0; c < size; c++)
	{
		int j = rest;
		const type *Adisp = Atmp +i2 * lda1 + c*acoars; 
#pragma unroll 8
		for(int i = 0; i < 8; i++)
		{
			regs[i] = offset[j];
			j+=32;
			if(j >= plain)
				break;
		}
#pragma unroll
		for(j = rest; j < plain; j+=32){
			tile[tile_disp+j] =  Adisp[j]; 
		}
		__syncthreads();
		type *Bdisp = Btmp + i2 * ldb1 + c*bcoars; 
		const  int tile_disp1 = i2 * size0;
		j = rest;
#pragma unroll 8
		for(int i = 0; i < 8; i++)
		{
			Bdisp[j] = alpha* tile[regs[i]+tile_disp1] + beta*Bdisp[j]; 
			j+=32;
			if(j >= plain)
				break;
		}
		__syncthreads();
	}
} 

	__device__ __forceinline__ 
void fvimatchl32_rem_coars(const type *  Atmp, type *  Btmp, const int lda1, const int ldb1,const  int remainderx,const  int remaindery,const  int size0, const int plain, const int sbp, const int ilimit, const int olimit, const unsigned short int* __restrict__ offset, const int acoars, const int bcoars, const int size, type alpha, type beta) 
{

	int id = threadIdx.x; 
	int i2 = id / 32; 
	int rest = id %32; 
	for(int c = 0; c< size; c++)
	{
		__syncthreads(); 
		if(i2 < remaindery){ 
			const type *Adisp = Atmp + i2 * lda1 + c*acoars; 
			const int tile_disp = i2 * (sbp); 
#pragma unroll		
			for( int  j = rest; j < ilimit; j+=32){ 
				tile[tile_disp+j] =  Adisp[j]; 
			} 
		} 
		__syncthreads(); 
		if(i2 >= remainderx)
			continue;
		type *Bdisp = Btmp + i2 * ldb1 + c*bcoars; 
		int regs[8];
		int j = rest;
#pragma unroll 8
		for(int i = 0; i < 8; i++)
		{
			if(j >= olimit)
				break;
			regs[i] = offset[j];
			j+=32;
		}
		const  int tile_disp1 = i2 * size0; 
		j = rest;
#pragma unroll 8
		for(int i = 0; i < 8; i++)
		{
			if(j >= olimit)
				break;
			Bdisp[j] = alpha* tile[regs[i]+tile_disp1] + beta*Bdisp[j]; 
			j+=32;
		}
	}
}

#define FNAME fvimatchl32.h
#include "includes/macro.h"
#undef FNAME
#define FNAME fvimatchl32_coars.h
#include "includes/macro.h"
#undef FNAME


void fvimatchl32CallerWrapper(const int ndim,const  type * __restrict__  A, type * B,const int size0, const  int param0, const int param1, const int numthreads1, const int numblocks, const int numthreads, const int shmsize
		, const int * __restrict__ lda_s, const int* __restrict__  ldb_s, const int* __restrict__ idx_s
		, const int remainder1,		const int remainder2, const int lda_kernel1, const int ldb_kernel1, unsigned short int* offset,
		const int ilimit, const int olimit, const int plain, const int sbp, const int ldal, const int ldap2l, const int acoars, const int bcoars, const int size,type alpha,type beta
		)
{
#ifdef printd
	printf("ndim = %d, size = %d, numblocks = %d, numthreads = %d, acoars = %d, bcoars = %d\n", ndim, size, numblocks, numthreads, acoars, bcoars);
#endif

	if(size > 0)
	{
#ifdef printd
		printf("Coarsening... No. of blocks = %d\n", numblocks/size);
#endif
		dim3 thread_blocks(numblocks/size, 1, 1);
		switch(ndim)
		{
			EXPANDDIMS(fvimatchl32_coars_kernel_, thread_blocks, numthreads, shmsize, (A,  B,size0, ldal, ldap2l,param0, param1,plain,sbp, lda_s,ldb_s,idx_s,remainder1,remainder2,ilimit,olimit,lda_kernel1, ldb_kernel1,offset, acoars, bcoars, size, alpha, beta))
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
			EXPANDDIMS(fvimatchl32_kernel_, thread_blocks, numthreads, shmsize, (A,  B,size0, ldal, ldap2l,param0, param1,plain,sbp, lda_s,ldb_s,idx_s,remainder1,remainder2,ilimit,olimit,lda_kernel1, ldb_kernel1,offset, alpha, beta))
			default:
			{
			}
		}
	}
}

void swap(int array[], int ind1, int ind2);

int  cancoarsen(int *lda, int newndim);


	extern "C" 
void  fvimatchl32_transpose_kernel(int ndim,   type *A, type *B,  const int *lda, const int *ldb, const int* params, const int * rperm, type alpha, type beta) 
{

	//printf("l32 normal\n");

	// int numBlocks = computeNumBlocksCode ; 
#ifdef printd
	printf("\nA Dims: %d \t %d \t %d\t %d\t %d\n", lda[0], lda[1], lda[2], lda[3], lda[4]);
	printf("\nParams: %d \t %d \t %d\t %d\t %d\t %d\t %d\n", params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
	printf("\nB Dims: %d \t %d \t %d\t %d\t %d\n", ldb[0], ldb[1], ldb[2], ldb[3], ldb[4]);
	printf("\n         R perm: %d \t %d \t %d\t %d\t %d\n", rperm[0], rperm[1], rperm[2], rperm[3], rperm[4]);
#endif	


	int numBlocks = params[6];//((size[1] + 8 -1)/8) * size[2] * ((size[3] + 8 -1)/8) * size[4] ; 

	int *d_lda_s, *d_ldb_s,  *d_idx_s; 
	const int size0 = lda[0];
	const int remainder1 = lda[1] % params[0];
	const int remainder2 = lda[params[3]] % params[0];
	int lda_s[20], ldb_s[20], idx_s[20], temp[20];
	lda_s[0] = 1;
	ldb_s[0] = 1;
	int i, blockA=params[0];
	int blockB = blockA;
	idx_s[1] = (lda[1] + blockA - 1) / blockA;
	lda_s[1] = lda_s[0] * lda[0];
	ldb_s[1] = ldb_s[0] * ldb[0];
	for(i = 2; i < ndim; i++)
	{
		if( i == params[3])
		{
			idx_s[i] = (lda[i] + blockA - 1)/blockA;
		}
		else
		{
			idx_s[i] = lda[i];

		}
		lda_s[i] = lda_s[i-1] * lda[i-1];
		ldb_s[i] = ldb_s[i-1] * ldb[i-1];
	}


	for(i = 1; i < ndim; i++)
	{
		temp[i] = ldb_s[rperm[i]];
	}
	const int lda_kernel1 = lda_s[params[3]];
	const  int ldb_kernel1 = ldb_s[params[4]];

	lda_s[1] *= blockA;
	lda_s[params[3]] *= blockA;
	temp[1] *= blockA;
	temp[params[3]] *= blockA;

	unsigned short int offset[9000];
	unsigned short int *d_offset;
	unsigned short limit = lda[0] * params[0];
	int tlimit = -1;
	for(i = 0; i < limit ; i++)
	{
		offset[i] = (i/lda[0]) * (lda[0] * params[0] + params[1])     +(i%lda[0]);
		if(i / lda[0] >= remainder2 && tlimit == -1) tlimit = i;
	}
#ifdef printd
	printf("Offset memory size = %d \n", limit);
#endif

	if(params[3] != 2)
	{
		swap(idx_s, 2, params[3]);
		swap(lda_s, 2, params[3]);
		swap(temp, 2, params[3]);
	}
	int newndim = ndim;

	int acoars = -1, bcoars = -1, size = -1;
#ifndef NOCOARSEN
	int noblock = 3;// (params[3] != 2);
	int cd = cancoarsen(idx_s+ noblock, ndim- noblock);

	if(cd >= 0)
	{
#ifdef printd
		printf("cd = %d, noblock = %d, ndim-noblock = %d\n", cd, noblock, ndim-noblock);
#endif
		acoars = lda_s[noblock+cd];
		bcoars = temp[noblock+cd];
		size = idx_s[noblock+cd];
		for(int j = cd+1; j < newndim; j++)
		{
			idx_s[noblock+j-1] = idx_s[noblock+j];
			lda_s[noblock+j-1] = lda_s[noblock+j];
			temp[noblock+j-1] = temp[noblock+j];

		}
		newndim--;
	}
#endif

	newndim--;

	SAFECUDAMALLOC(&d_offset,limit*sizeof(short)); 
	SAFECUDAMEMCPY(d_offset, offset,limit*sizeof(short), cudaMemcpyHostToDevice); 
	SAFECUDAMALLOC(&d_lda_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_ldb_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_idx_s,newndim*sizeof(int)); 
	SAFECUDAMEMCPY(d_idx_s, idx_s+1,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_lda_s, lda_s+1,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_ldb_s, temp+1,newndim*sizeof(int), cudaMemcpyHostToDevice); 

	const int ilimit = remainder1 * size0;
	const int olimit = remainder2 * size0;
	const int plain = params[0] * size0;
	const int sbp = plain+params[1];
	const int ldal = (lda[1] - remainder1)/blockA;
	const int ldp2l = (lda[params[3]] - remainder2)/blockB;
/*
#ifdef MODEL
printf("\t%d\t%d\t", plain, blockA);
printf("\t%d\t%d\t", plain/32, plain%32);
double f1, f2, f3, f4, f;
int minlimit = min(ilimit, olimit);
printf("\tf1=%lf\t", f1 =  ((plain/32)  + (double)(plain%32) /32)/ (int)((plain+31)/32));
printf("\tf2=%lf\t", f2 =  ((ilimit/32)  + (double)(ilimit%32) /32)/ (int)(max(1,(ilimit+31)/32)));
printf("\tf3=%lf\t", f3 =  ((olimit/32)  + (double)(olimit%32) /32)/ (int)(max(1,(olimit+31)/32)));
printf("\tf4=%lf\t", f4 =  ((minlimit/32)  + (double)(minlimit%32) /32)/ (int)(max(1,(minlimit+31)/32)));
//printf("\t%d\t%d\t", lda[1], ldb[1]);
int asize = lda[1];
int bsize = lda[1];
//printf("\t%d\t%d\t%d\t%d\t", asize/blockA, asize%blockA, bsize/blockB,bsize%blockB );
//int amax = min(blockA, 32);
//int bmax = min(blockB, 32);
int amax = blockA;
int bmax = blockB;
printf("\tf=%lf\t", f = ((asize/amax) * (bsize/bmax) *f1 + (double)(asize/amax) * (bsize%bmax > 0) *f3+ (double)(asize%amax>0) * (bsize/bmax)*f2 + (double)(asize%amax > 0) * (bsize%bmax > 0) *f4 )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax)));
printf("\t%lf\t", f);
#endif
*/

#ifdef NOHTIME
#include "includes/nohtimestart.h"
#endif
	fvimatchl32CallerWrapper(newndim, A,  B,lda[0],params[0], params[1], params[3]-1,
			numBlocks, params[2], params[5]*sizeof(type)
			, d_lda_s,d_ldb_s,d_idx_s
			,remainder1,remainder2,lda_kernel1, ldb_kernel1, d_offset, ilimit, olimit, plain, sbp, ldal, ldp2l, acoars, bcoars, size, alpha, beta);


#ifdef NOHTIME
#include "includes/nohtimestop.h"
#endif


	{cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess){
			printf("\nKernel ERROR in fvimatchl32: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);
			//exit(-1);
		}}
	cudaFree(d_lda_s);
	cudaFree(d_ldb_s);
	cudaFree(d_idx_s);
	cudaFree(d_offset);
}



