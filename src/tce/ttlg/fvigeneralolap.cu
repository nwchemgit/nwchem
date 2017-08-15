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

extern __shared__ type tile[];


__device__ __forceinline__  void fvinomgeneralolap_main_coars(const type * __restrict__ Atmp, type *  Btmp,  const int tb_size,  const int* __restrict__ aexpr, const int* __restrict__ bexpr, const int* __restrict__ texpr1, const int * __restrict__  texpr2, const int ilimit, const int olimit, const int rowinc, const int shm2, const int numelems_blk, const int acoars, const int bcoars, const int size,
		type alpha, type beta) 
{
	//if(blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 && threadIdx.x == 1)
	//printf("ilimit = %d, olimit = %d, tbsize = %d, Atmp = %p, Btmp = %p, txpr2[10] = %d, numblocks = %d \n",ilimit, olimit, tb_size, Atmp, Btmp, texpr2[10], gridDim.x*gridDim.y*gridDim.z); 
	for(int i = 0; i < size; i++)
	{
		for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
		{
			tile[texpr1[Id]] = Atmp[aexpr[Id] + i*acoars];
		}
		__syncthreads();
		for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
		{
			Btmp[bexpr[Id] + i *bcoars] = alpha* tile[texpr2[Id]] + beta* Btmp[bexpr[Id] + i *bcoars];
		}
		__syncthreads();
	}

}
__device__ __forceinline__ void fvinomgeneralolap_rem_coars(const type * __restrict__ Atmp, type * __restrict__ Btmp,  const int tb_size,  const int* __restrict__ aexpr, const int* __restrict__ bexpr, const int* __restrict__ texpr1, const int * __restrict__  texpr2, const int ilimit, const int olimit, const int rowinc, const int shm2, const int ilimitr, const int olimitr, const int numelems_blk, const int acoars, const int bcoars, const int size, type alpha, type beta) 
{
	for(int i = 0; i < size; i++)
	{
		for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
		{
			tile[texpr1[Id]] = Atmp[aexpr[Id]+i * acoars];
		}
		__syncthreads();
		for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
		{
			int toffset2 = texpr2[Id];
			if(toffset2 % (shm2) < ilimitr && toffset2/shm2 < olimitr)
				Btmp[bexpr[Id]+i*bcoars] =  alpha* tile[toffset2] + beta *  Btmp[bexpr[Id]+i*bcoars];
		}
		__syncthreads();
	}
}
__device__ __forceinline__  void fvinomgeneralolap_main(const type * __restrict__ Atmp, type *  Btmp,  const int tb_size,  const int* __restrict__ aexpr, const int* __restrict__ bexpr, const int* __restrict__ texpr1, const int * __restrict__  texpr2, const int ilimit, const int olimit, const int rowinc, const int shm2, const int numelems_blk, type alpha, type beta) 
{
	//if(blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 && threadIdx.x == 1)
//	printf("ilimit = %d, olimit = %d, tbsize = %d, Atmp = %p, Btmp = %p, txpr2[10] = %d, numblocks = %d \n",ilimit, olimit, tb_size, Atmp, Btmp, texpr2[10], gridDim.x*gridDim.y*gridDim.z); 

	for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
	{
		tile[texpr1[Id]] = Atmp[aexpr[Id]];
	}
	__syncthreads();
	for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
	{
		int toffset2 = texpr2[Id];
			Btmp[bexpr[Id]] =alpha* tile[toffset2] + beta*Btmp[bexpr[Id]];
	}

}
__device__ __forceinline__ void fvinomgeneralolap_rem(const type * __restrict__ Atmp, type * __restrict__ Btmp,  const int tb_size,  const int* __restrict__ aexpr, const int* __restrict__ bexpr, const int* __restrict__ texpr1, const int * __restrict__  texpr2, const int ilimit, const int olimit, const int rowinc, const int shm2, const int ilimitr, const int olimitr, const int numelems_blk, type alpha, type beta) 
{
	for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
	{
		tile[texpr1[Id]] = Atmp[aexpr[Id]];
	}
	__syncthreads();
	for(int Id=threadIdx.x; Id < numelems_blk; Id+= tb_size)
	{
		int toffset2 = texpr2[Id];
		if(toffset2 % (shm2) < ilimitr && toffset2/shm2 < olimitr)
			Btmp[bexpr[Id]] = alpha* tile[toffset2] + beta*Btmp[bexpr[Id]];
	}
}

#define FNAME fvigeneralolap.h
#include "includes/macro.h"
#undef FNAME
#define FNAME fvigeneralolap_coars.h
#include "includes/macro.h"
#undef FNAME
/*
//general
__global__ void fvinomgeneralolap_kernel (const int ndim, const type *  A, type * B, const int ilimit, const int olimit,const int param2, const int param3, const int param4
, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
, const int remainder1,         const int remainder2,
const int* __restrict__ offseti, const int* __restrict__ offseto, const int* __restrict__ tile1, const int* __restrict__ tile2,
const int ilimitr, const int olimitr, const int inputrem, const int outputrem, const int rowinc, const int shm2, const int numelems_blk
)
{

int tmp;
int val0 = blockIdx.x;
int val1 = blockIdx.y;
int aexpr =0, bexpr = 0;
if(ndim > 1)
aexpr = val0 * lda_s[0] + val1*lda_s[1], bexpr = val0 * ldb_s[0] + val1 * ldb_s[1];
else
{
aexpr = val0 * lda_s[0], bexpr = val0 * ldb_s[0];
}
int idx;
idx = blockIdx.z;
int ii1 = -1 , iip2 = -1;

if(param2 == 0)
iip2 = blockIdx.x;
else if(param2 == 1)
iip2 = blockIdx.y;
if(param3 == 0)
ii1 = blockIdx.x;
else if(param3 == 1)
ii1 = blockIdx.y;
#pragma unroll
for(int i = 2; i < ndim; i++)
{
tmp = idx/idx_s[i];
int index = idx - tmp * idx_s[i];
aexpr += index * lda_s[i];
bexpr += index * ldb_s[i];
idx = tmp;
if(i == param2) ii1 = index;
else if(i == param3) iip2 = index;
}
const double *Atmp = A + aexpr;
double *Btmp = B + bexpr;
if(ii1 < inputrem && iip2 < outputrem)
{
fvinomgeneralolap_cuSharedMemTranspose_vec256(Atmp,Btmp,  param4,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, numelems_blk );
}
else if(ii1 >= inputrem && iip2 < outputrem)
{ //remainder in size1

fvinomgeneralolap_cuSharedMemRemTranspose_vec256(Atmp,Btmp,  param4,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, ilimitr, olimit, numelems_blk );
}
else if(iip2 >= outputrem && ii1 < inputrem)
{ //remainder in size2 
fvinomgeneralolap_cuSharedMemRemTranspose_vec256(Atmp,Btmp,  param4,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, ilimit, olimitr, numelems_blk );
}
else
{
fvinomgeneralolap_cuSharedMemRemTranspose_vec256(Atmp,Btmp,  param4,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, ilimitr, olimitr, numelems_blk );
}


return;
//#undef ndim
}
 */



void fvinomgeneralolap_CallerWrapper(int ndim, type *  A, type * B,const int ilimit, const int olimit, const int blockAI, const int blockBI, const int numblocks, const int numthreads, const  int shm
		, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s,
		const int coarsa, const int coarsb, const int * idx_ss, const int shm2,
		const int* __restrict__ aexpr, const int* __restrict__ bexpr, const int* __restrict__ texpr1, const int* __restrict__ texpr2, const int ilimitr, const int olimitr,
		const int inputrem, const int outputrem, const int numelems_blk, const int size, type alpha, type beta
		)

{

	/*	int second, third;
		if(ndim > 2)
		{
		second = idx_ss[1]; third = numblocks/(idx_ss[0]*idx_ss[1]);
		}
		else if(ndim > 1)
		{
		second = idx_ss[1]; third = 1;// numblocks/(idx_ss[0]*idx_ss[1]);
		}
		else
		{
		second = third = 1;
		}
		dim3 thread_blocks(idx_ss[0], second, third);*/
	const int rowinc = (numthreads+ilimit-1)/ilimit;

#ifdef printd
	printf("thread_blocks = %d, numthreads = %d, shm = %d\n", numblocks, numthreads, shm);
	printf("size = %d, ndim = %d, shm = %d\n", size, ndim, shm);
#endif
	if(size > 0)
	{
		dim3 thread_blocks(numblocks/size, 1, 1);
		switch(ndim)
		{
			EXPANDDIMS(fvinomgeneralolap_coars_kernel_, thread_blocks, numthreads, shm, (A,  B, ilimit, olimit, blockAI, blockBI, numthreads, lda_s,ldb_s, idx_s, coarsa,coarsb, aexpr, bexpr, texpr1, texpr2, ilimitr, olimitr, inputrem, outputrem, rowinc, shm2, numelems_blk, size, alpha, beta)) 
			default:
			{
				//			fvinomgeneralolap_coars_kernel<<<thread_blocks, numthreads, shm>>>(ndim, A,  B, ilimit, olimit, blockAI, blockBI, numthreads, lda_s,ldb_s, idx_s, coarsa, coarsb, aexpr, bexpr, texpr1, texpr2, ilimitr, olimitr, inputrem, outputrem, rowinc, shm2, numelems_blk);
			}
		}
	}
	else
	{
		dim3 thread_blocks(numblocks, 1, 1);
		switch(ndim)
		{
			EXPANDDIMS(fvinomgeneralolap_kernel_, thread_blocks, numthreads, shm, (A,  B, ilimit, olimit, blockAI, blockBI, numthreads, lda_s,ldb_s, idx_s, aexpr, bexpr, texpr1, texpr2, ilimitr, olimitr, inputrem, outputrem, rowinc, shm2, numelems_blk, alpha, beta)) 
			default:
			{
				//			fvinomgeneralolap_kernel<<<thread_blocks, numthreads, shm>>>(ndim, A,  B, ilimit, olimit, blockAI, blockBI, numthreads, lda_s,ldb_s, idx_s, aexpr, bexpr, texpr1, texpr2, ilimitr, olimitr, inputrem, outputrem, rowinc, shm2, numelems_blk);
			}
		}
	}

}

int ispresent(int a, int *array, int n)
{
	for(int i = 0; i < n; i++)
	{
		if(array[i] == a) return 1;
	}
	return 0;
}

int getoff(int index,const int * dims,const int * stride, int n)
{
	int ret = 0;
	for(int i = 0; i < n; i++)
	{
		int ii = index % dims[i];
		ret += ii * stride[i];
		index/= dims[i];
	}

	return ret;
} 
void makeconsecutive(int *tmp, int k)
{
	int tmp2[20], permi[20], j,i;
	//printf("\n perm before: ");
	//for(i = 0; i < k; i++) printf("%d ", tmp[i]);//permi[i];
	for(i = 0; i < k; i ++) tmp2[i] = tmp[i];
	for(i = 0; i < k; i ++)
	{
		for(j = 0; j < k-i-1; j++)
		{
			if(tmp2[j] > tmp2[j+1])
			{
				int tmp = tmp2[j];
				tmp2[j] = tmp2[j+1];
				tmp2[j+1] = tmp;
			}

		}
	}
	for(i = 0; i < k; i++)
	{
		for(j=0; j <k; j++)
		{
			if(tmp[i] == tmp2[j])
			{
				permi[i] = j; break;
			}
		}
	}
	for(i = 0; i < k; i++) tmp[i] = permi[i];
}

void swap(int array[], int ind1, int ind2);

int  cancoarsen(int *lda, int newndim)
{
	if(newndim < 1) return -1;
	unsigned long vol = 1;
	for(int i = 0; i < newndim; i++)
	{
		vol *= lda[i];
	}
	if(vol < 32*100) return -1;
	for(int i = 0; i < newndim; i++)
	{
		if(lda[i] >= 4 && lda[i] <= 31)
			return i;
	}
	/*	for(int i = 0; i < newndim; i++)
		{
		if(lda[i] >= 2 && lda[i] <= 300)
		return i;
		}*/
	return -1;
}

	extern "C" 
void  fvigeneralolap_transpose_kernel(int ndim, type *A, type *B,  int *lda, const int *ldb, const int* params, const int * perm, const int* rperm, type alpha, type beta) 
{
	// int numBlocks = computeNumBlocksCode ; 
#ifdef printd
	printf("\nA Dims: %d \t %d \t %d\t %d\t %d\n", lda[0], lda[1], lda[2], lda[3], lda[4]);
	printf("\nAll diff Params: %d \t %d \t %d\t %d\t %d\t %d\t %d\t %d\t %d \t%d\t %d\t %d\n", params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10], params[11]);
	printf("\nB Dims: %d \t %d \t %d\t %d\t %d\n", ldb[0], ldb[1], ldb[2], ldb[3], ldb[4]);
	printf("\nR perm: %d \t %d \t %d\t %d\t %d\n", rperm[0], rperm[1], rperm[2], rperm[3], rperm[4]);
#endif
	int alimit = params[3];
	int blimit = params[4];
	int blockA=params[0];
	int blockB = params[11];
	int ilimit = params[7];
	//int olimit = params[8];
	int i = 0, j = 0;
	int size = 1;

	for(i = 0; i < blimit; i++)
	{
		if(perm[i] > alimit)
		{
			size *= ldb[i];
		}
	}
	if(perm[i] > alimit)
	{
		if(blockB == 1)
			size *= ldb[i];
		else
			size*= blockB;
	}

#ifdef printd
	printf("In .cu, alimit = %d, blimit = %d, bsize = %d, blockB = %d, blockA = %d\n", alimit, blimit, size, blockB, blockA);
#endif
	//for(int y = 0; y < j; y++)
	//printf("bo[%d] = %d ",y, bo[y]);
	int olimit = size;
#ifdef SLICE
printf("\t%d\t%d\t", ilimit, olimit);
#endif
	//exit(0);
	int numBlocks = params[6];//((size[1] + 8 -1)/8) * size[2] * ((size[3] + 8 -1)/8) * size[4] ; 
	const int pad =  ((ilimit %2)+1)%2;

	int *d_lda_s, *d_ldb_s,  *d_idx_s; 
	const int remainder1 = lda[params[3]] % blockA;
	int remainder2;
	remainder2 = lda[perm[params[4]]] % blockB;
	const int ilimitr = (ilimit * remainder1) / blockA; 
	int olimitr = (olimit * remainder2) / blockB; 
	int irem, orem;

	if(remainder1 == 0) irem = lda[alimit];
	else irem = (lda[alimit] - remainder1)/blockA;
	if(remainder2 == 0) orem = lda[perm[blimit]];
	else orem = (ldb[blimit] - remainder2)/blockB;
	if(perm[params[4]] == params[3])
	{
		olimitr = olimit;
		remainder2 = 0;
		orem = lda[perm[blimit]];
	}
	else
	{
		//remainder2 = lda[perm[params[4]]] % blockB;
	}
#ifdef printd
	printf("\nrem1 = %d, rem2 = %d\n", remainder1, remainder2);
	printf("\nilimit = %d, olimit = %d, ilimitr = %d, olimitr = %d\n", ilimit, olimit, ilimitr, olimitr);
#endif

	int *input_base, *output_base, *tile_base1, *tile_base2;
	int *aexpr, *bexpr, *texpr1, *texpr2;
	//	int *ablock, *bblock, *d_ablock, *d_bblock;;

	int lda_s[20], ldb_s[20], idx_s[20], temp[20];
	lda_s[0] = 1;
	ldb_s[0] = 1;
	idx_s[0] = 1;
	for(i = 1; i < alimit; i++)
	{
		idx_s[i] = 1;
		lda_s[i] = lda_s[i-1] * lda[i-1];
		ldb_s[i] = ldb_s[i-1] * ldb[i-1];
	}
	if(rperm[alimit] < blimit || rperm[alimit] == blimit && blockB == 1)
		idx_s[alimit] = 1;
	else{
		idx_s[alimit] = (lda[alimit] + blockA - 1) / blockA;
	}
	lda_s[i] = lda_s[i-1] * lda[i-1];
	ldb_s[i] = ldb_s[i-1] * ldb[i-1];

	for(i = alimit+1; i < ndim; i++)
	{
		lda_s[i] = lda_s[i-1] * lda[i-1];
		ldb_s[i] = ldb_s[i-1] * ldb[i-1];
		if(rperm[i] < blimit)
		{
			idx_s[i] = 1;// (lda[i] + blockA - 1) / blockA;
		}
		else if(rperm[i] == blimit)
		{
			idx_s[i] = (lda[i] + blockB - 1) / blockB;
		}
		else
		{
			idx_s[i] = lda[i];
		}
	}

	for(i = 0; i < ndim; i++)
	{
		temp[i] = ldb_s[rperm[i]];
#ifdef printd
		printf("Idx[%d] = %d\n", i, idx_s[i]);
#endif
	}

	aexpr = (int*)malloc(ilimit* olimit * sizeof(int));
	bexpr = (int*)malloc(ilimit * olimit * sizeof(int));
	texpr1 = (int*)malloc(ilimit* olimit * sizeof(int));
	texpr2 = (int*)malloc(ilimit * olimit* sizeof(int));

	SAFECUDAMALLOC(&input_base,ilimit*olimit*sizeof(int)); 
	SAFECUDAMALLOC(&output_base,ilimit*olimit*sizeof(int)); 
	SAFECUDAMALLOC(&tile_base1, ilimit*olimit *sizeof(int)); 
	SAFECUDAMALLOC(&tile_base2, ilimit*olimit*sizeof(int)); 

	int outD[20],  outD_s[20], B_s[20];
	outD_s[0] = 1;
	outD[0] = ldb[0];

	int inD[20], inD_s[20];
	inD_s[0] = 1;
	inD[0] = lda[0];
	//B_s[0] = ldb_s[rperm[0]];
	B_s[0] = 1;//ldb_s[rperm[0]];

	int permD_s[20], permD[20];
	//	permD_s[alimit+1] = outD_s[0];
	//permD[alimit+1] = outD[0];
	int OO_C = 0, C_C = 0;
	int onlyOut[20];int onlyOutI[20];

	for(i = 0; i <= alimit; i++)
	{
		if(rperm[i] <= blimit)
			C_C++;	
	}
	const int OI_C = alimit + 1;
	if(perm[0] > alimit)
	{
		OO_C++;
		onlyOut[0] = ldb[0];
		onlyOutI[0] = 0; 
		permD[0] = OI_C;
	}
	else 
	{
	//	C_C++;
	permD[0] = perm[0];
	}
	for(i = 1; i < blimit; i++)
	{
		outD[i] = ldb[i];
		outD_s[i] = outD_s[i-1] * outD[i-1];
		//	B_s[i] = ldb_s[rperm[i]];
		B_s[i] = ldb_s[i];
		if(perm[i] > alimit)
		{
			onlyOut[OO_C] = ldb[i];
			onlyOutI[OO_C++] = i;
			permD[i] = OI_C + i;	
		}
		else
		{
	//		C_C++;
		permD[i] = perm[i];	
		}
	}
	if(blimit == 0) {i = 0; OO_C = 0;}
	if(blockB == 1)
	{
		outD[i] = ldb[i];
	}
	else
	{
		outD[i] = blockB;
	}
	if(i > 0)
	{
		outD_s[i] = outD_s[i-1] * outD[i-1];
	}
	else
	{
		outD_s[i] =  1;
	}
	B_s[i] = ldb_s[i];
	if(perm[i] <= alimit)
	{
	//	C_C++;
	permD[i] = perm[i];	
	}
	else
	{
		//if(blockB == 1)
		onlyOut[OO_C] = ldb[i];
		//else
		//onlyOut[OO_C] = blockB;
		onlyOutI[OO_C++] = i;
		permD[i] = OI_C + i;	
	}
	i++;
	for(j  = 0; j < alimit; j++)
	{
		if(rperm[j] > blimit)
		{
			outD[i] = lda[j];
			permD[i] = j;	
			outD_s[i] = outD_s[i-1] * outD[i-1];
			B_s[i] = ldb_s[rperm[j]];
			i++;
		}
		else
		{ 
			//tmp [k++] = j;
			//C_C++;
		}

	}
	if(rperm[j] > blimit)
	{
		//printf("BI = %d, rperm[%d] = %d, blimit = %d\n", i, j, rperm[j], blimit);
		if(blockA == 1)
		{
			outD[i] = lda[j];
		}
		else
		{
			outD[i] = blockA;
		}
		permD[i] = j;	
		B_s[i] = ldb_s[rperm[j]];
		outD_s[i] = outD_s[i-1] * outD[i-1];
		i++;
	}
	else
	{ 
		//tmp [k++] = rperm[j];
		//C_C++;
	}
	int BI = i;

	for(i  = 1; i < alimit; i++)
	{
		inD[i] = lda[i];
		inD_s[i] = inD_s[i-1] * inD[i-1] ;
	}
	if(alimit == 0) i = 0;
	if(blockA == 1)
	{
		inD[i] = lda[i];
	}
	else
	{
		inD[i] = blockA;
	}
	if(i > 0)
		inD_s[i] = inD_s[i-1] * inD[i-1];
	i++;

	for(j  = 0; j < blimit; j++)
	{
		if(perm[j] > alimit)
		{
			inD[i] = ldb[j];
			inD_s[i] = inD_s[i-1] * inD[i-1];
			i++;
		} 
		else{
			//	C_C++;
		}

	}
	if(perm[j] > alimit)
	{
		if(blockB == 1)
		{
			inD[i] = ldb[j];
		}
		else
		{
			inD[i] = blockB;
		}
		inD_s[i] = inD_s[i-1] * inD[i-1];
		i++;
	}
	else{
		//		C_C++;
	}
	int AI = i;
	makeconsecutive(permD, AI);
	//permD[0] = 1, permD[1] = 2, permD[2] = 0;
	//inD[0] = 32, inD[1] = 2, inD[2] = 30;
	//inD_s[0] = 1, inD_s[1] = 32, inD_s[2] = 64;
	for(i = 0; i < AI; i++)
	{
		permD_s[i] = inD_s[permD[i]];
	}
	if(BI != AI)
	{
		printf("No. of dimensions in I and O non-matching...\n");
		//return;
	}
#ifdef printd
	printf("\nOO_C = %d, C_C = %d\n ", OO_C, C_C);
	printf("\nOUT_D: ");
	for(int i = 0; i < BI; i++)
	{
		printf("%d ",outD[i]);
	}
	printf("\nOUT_D_S: ");
	for(int i = 0; i < BI; i++)
	{
		printf("%d ",outD_s[i]);
	}
	printf("\nIn_D: ");
	for(int i = 0; i < AI; i++)
	{
		printf("%d ",inD[i]);
	}
	printf("\nIN_D_S: ");
	for(int i = 0; i < AI; i++)
	{
		printf("%d ", inD_s[i]);
	}
	printf("\n");
	printf("\nB_S: ");
	for(int i = 0; i < BI; i++)
	{
		printf("%d ", B_s[i]);
	}
	printf("\nPerm_D: ");
	for(int i = 0; i < BI; i++)
	{
		printf("%d ",permD[i]);
	}
	printf("\n");
	printf("\nPerm_S: ");
	for(int i = 0; i < BI; i++)
	{
		printf("%d ", permD_s[i]);
	}
	printf("\n");
	printf("\n");
#endif

	for(int rowId=0; rowId < olimit; rowId++)
	{
		int tmp = rowId;
		int aoff=0,j;

		for(j = 0; j < OO_C; j++)
		{
			int dval = onlyOut[j];
			int val = tmp%dval;
			tmp /= dval;
			aoff += val * lda_s[perm[onlyOutI[j]]];
		}
		for(int colId=0; colId < ilimit; colId++)
		{
			aexpr[rowId*ilimit + colId] = aoff + colId;
			texpr1[rowId * (ilimit) + colId] =  rowId * (ilimit+pad) + colId;

			int off = getoff(rowId * ilimit + colId, outD, permD_s, BI);	
			texpr2[rowId * (ilimit) + colId] =  off + pad * (off/ilimit);
			off = getoff(rowId * ilimit + colId, outD, B_s, AI);	
			bexpr[rowId* ilimit + colId] = 	off;

		}
	}

#ifdef printd
	printf("\n...A...\n");
	for(int rowId=0; rowId < olimit; rowId++)
	{
		printf("%d ", aexpr[rowId]);
		//printf("\n");
		//for(int colId=0; colId < ilimit; colId++)
		{
			//		printf("%d ", bexpr[rowId * ilimit + colId]);
			//	printf("%d ", texpr2[rowId * ilimit + colId]);
		}
		//printf("\n");
	}
	printf("\n...B...\n");
	for(int rowId=0; rowId < olimit; rowId++)
	{
		//	printf("%d ", bexpr[rowId]);
		//printf("\n");
		for(int colId=0; colId < ilimit; colId++)
		{
			printf("%d ", bexpr[rowId * ilimit + colId]);
			//	printf("%d ", texpr2[rowId * ilimit + colId]);
		}
		printf("\n");
	}
	printf("\n...T...\n");
	for(int rowId=0; rowId < olimit; rowId++)
	{
		//		printf("%d ", bexpr[rowId]);
		//	printf("\n");
		for(int colId=0; colId < ilimit; colId++)
		{
			//	printf("%d ", bexpr[rowId * ilimit + colId]);
			printf("%d ", texpr2[rowId * ilimit + colId]);
		}
		printf("\n");
	}
#endif

	lda_s[params[3]] *= params[0];///lda_s[i-1] * lda[i-1];
	temp[params[3]] *= params[0];// ldb_s[i-1] * ldb[i-1];

	if(params[3] != perm[params[4]])//no double blocking
	{
		lda_s[perm[params[4]]] *= params[11];///lda_s[i-1] * lda[i-1];
		temp[perm[params[4]]] *= params[11];// ldb_s[i-1] * ldb[i-1];
	}

	int c = 0, d = 0;
	c = alimit + 1;//c = No. of dimensions to be removed from input for thread blocking, b = same for output but only for those which are not in input
	if(blockA > 1) c--;
	int ablockI, bblockI;
	//int dims[20];
	ablockI = alimit-c;
	bblockI = perm[blimit]-c;
	int tempbblockI = bblockI;
	#ifdef printd
		printf("\nablockI = %d, bblockI = %d\n", ablockI, bblockI);
	#endif

 for(int i = c; i < ndim; i++)
        {
                if(((rperm[i] < blimit) || ((rperm[i] == blimit) && (blockB ==1))))
                {
                        idx_s[i] = 1;// idx_s[j];
                        /*for(int j = i+1; j < ndim-d; j ++)
                        {
                                idx_s[j-1] = idx_s[j];
                                lda_s[j-1] = lda_s[j];
                                temp[j-1] = temp[j];
                        }*/
                        d++;
                        if((i < bblockI + c)  || (i == bblockI + c) && (blockB == 1))
                                tempbblockI--;
                }
        }
 bblockI = tempbblockI;
        int cnt = 0;
        for(int i = c; i < ndim; i++)
        {
                if(idx_s[i] == 1)
                {
                        for(int j = i+1; j < ndim; j ++)
                        {
                                idx_s[j-1] = idx_s[j];
                                lda_s[j-1] = lda_s[j];
                                temp[j-1] = temp[j];
                        }
                        cnt++;
                        i--;

                }
                if(cnt > ndim) break;
        }

/*

	for(int i = c; i < ndim; i++)
	{
		//	dims[i] = lda[i];
		if((rperm[i] < blimit || (rperm[i] == blimit && blockB ==1)))
		{ 
			for(j = i+1; j < ndim-d; j ++)
			{
				idx_s[j-1] = idx_s[j];
				lda_s[j-1] = lda_s[j];
				temp[j-1] = temp[j]; 
				//	dims[j-1] = lda[j];
			}
			d++;
			if((i < bblockI + c)  || (i == bblockI + c) && (blockB == 1))
				bblockI--;
			if(i < alimit)
				ablockI--;
		}
	}*/
	int newndim = ndim - (c + d);

	#ifdef printd
		printf("\nChanged ablockI = %d, bblockI = %d\n", ablockI, bblockI);
	#endif
	//Find the largest dimension and make it the first as only Dimx can have > 65k size
	/*int max = 0;
	for(int i = 1; i < newndim; i++)
	{
		if(idx_s[c+i] > idx_s[max+c]) max = i;
	}
	//printf("\nmax: %d ", max);
	if(max > c)
	{
		swap(idx_s, c, max+c);
		swap(lda_s, c, max+c);
		swap(temp, c, max+c);

		if(max == ablockI) ablockI = 0;
		else if(ablockI == 0) ablockI = max;
		if(max == bblockI) bblockI = 0;
		else if(bblockI ==0) bblockI = max;
	}*/
	if(ablockI > 0)//move it to start, junk part
	{
		swap(idx_s, ablockI+c, c);
		swap(lda_s, ablockI+c, c);
		swap(temp, ablockI+c, c);
	}
	int bi = 0;
	if(bblockI == 0 && ablockI > 0)
		bi = ablockI+c;
	else
		bi = bblockI+c;
	if(bblockI >= 0 && bblockI != ablockI)//move it to start
	{
		swap(idx_s, bi, c+ (ablockI >=0));
		swap(lda_s, bi, c+ (ablockI >=0));
		swap(temp, bi, c+ (ablockI >=0));
	}
	if(bblockI >= 0) { 
		if(ablockI == bblockI || ablockI < 0) bblockI = 0;
		else bblockI = 1;
	}

	if(ablockI >=0) ablockI = 0;

	int nblkdims = 0;
	if(ablockI >= 0) nblkdims++;
	if((bblockI >= 0) && (bblockI != ablockI)) nblkdims++;
#ifdef printd
	printf("\nIDx: ");
	for(int i = 0; i < newndim; i++)
	{
		printf("%d ",idx_s[i+c]); 
	}
	printf("ndim = %d, c = %d, d = %d, newndim = %d, nblkdims = %d\n", ndim, c, d, newndim, nblkdims);
#endif
	int acoars = 0, bcoars = 0;
	size = -1; 
#ifdef printd
	printf("\nirem = %d, orem = %d, alimit = %d, blimit = %d, ablockI = %d, bblockI = %d, newdim = %d, c = %d, d = %d, olddim = %d\n\n", irem, orem, alimit, blimit, ablockI, bblockI, newndim, c, d, ndim);
#endif
#ifndef NOCOARSEN
	int cd = cancoarsen(idx_s+c+nblkdims, newndim-nblkdims);
	if(cd >= 0)
	{
		int offset = c + cd + nblkdims;
		acoars = lda_s[offset];
		bcoars = temp[offset];
		size = idx_s[offset];
		for(int j = cd+1+nblkdims; j < newndim; j++)
		{
			idx_s[c+j-1] = idx_s[c+j];
			lda_s[c+j-1] = lda_s[c+j];
			temp[c+j-1] = temp[c+j];

		}	
	//	ablockI--;
	//	bblockI--;
		newndim--;
	}
#ifdef printd
	printf("\nirem = %d, orem = %d, alimit = %d, blimit = %d, ablockI = %d, bblockI = %d, newdim = %d, c = %d, d = %d, olddim = %d, acoars = %d, bcoars = %d, cd = %d\n\n", irem, orem, alimit, blimit, ablockI, bblockI, newndim, c, d, ndim, acoars, bcoars, cd);
#endif
#endif

	SAFECUDAMALLOC(&d_lda_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_ldb_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_idx_s,newndim*sizeof(int)); 
	SAFECUDAMEMCPY(d_idx_s, idx_s+c,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_lda_s, lda_s+c,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_ldb_s, temp+c,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(input_base, aexpr, ilimit*olimit*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(output_base, bexpr, ilimit*olimit*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(tile_base1, texpr1, ilimit*olimit*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(tile_base2, texpr2,ilimit* olimit*sizeof(int), cudaMemcpyHostToDevice); 

#ifdef MODEL
	{
const int olimit = params[8]; 
int olimitr = (olimit * remainder2) / blockB; 
printf("\t%d\t%d\t", ilimit, olimit);
printf("\t%d\t%d\t%d\t%d\t", ilimit/32, ilimit%32, olimit/32,olimit%32 );
double f1, f2, f3, f4, f;
printf("\tf1=%lf\t", f1 =  ((ilimit/32) * (olimit/32) + (double)(ilimit/32) * (olimit%32) /32+ (double)(ilimit%32) * (olimit/32) /32 + (double)(ilimit%32) * (olimit%32) /(32*32) )/ (int)(((ilimit+31)/32) * ((olimit+31)/32)));
printf("\tf2=%lf\t", f2 =  ((ilimitr/32) * (olimit/32) + (double)(ilimitr/32) * (olimit%32) /32+ (double)(ilimitr%32) * (olimit/32) /32 + (double)(ilimitr%32) * (olimit%32) /(32*32) )/ max(1,(int)(((ilimitr+31)/32) * ((olimit+31)/32))));
printf("\tf3=%lf\t", f3 =  ((ilimit/32) * (olimitr/32) + (double)(ilimit/32) * (olimitr%32) /32+ (double)(ilimit%32) * (olimitr/32) /32 + (double)(ilimit%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((ilimit+31)/32) * ((olimitr+31)/32))));
printf("\tf4=%lf\t", f4 =  ((ilimitr/32) * (olimitr/32) + (double)(ilimitr/32) * (olimitr%32) /32+ (double)(ilimitr%32) * (olimitr/32) /32 + (double)(ilimitr%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((olimitr+31)/32) * ((olimitr+31)/32))));
printf("\t%d\t%d\t", lda[alimit], ldb[blimit]);
int asize = lda[alimit];
int bsize = ldb[blimit];
printf("MKL \t%d\t%d\t%d\t%d\t", asize/blockA, asize%blockA, bsize/blockB,bsize%blockB );
//int amax = min(blockA, 32);
//int bmax = min(blockB, 32);
int amax = blockA;
int bmax = blockB;
printf("\tf=%lf\t", f = ((asize/amax) * (bsize/bmax) *f1 + (double)((asize/amax) * (bsize%bmax > 0) *f3)+ (double)((asize%amax > 0) * (bsize/bmax)*f2)  + (double)((asize%amax>0) * (bsize%bmax > 0) *f4) )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax)));

//printf("\tf=%lf\t", f = ((asize/amax) * (bsize/bmax) *f1 + (double)(asize/amax) * (bsize%bmax) *f3/bmax+ (double)(asize%amax) * (bsize/bmax)*f2 /amax + (double)(asize%amax) * (bsize%bmax) *f4/(amax*bmax) )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax)));
printf("\t%lf\t", f);
	}
#endif

#ifdef NOHTIME
#include "includes/nohtimestart.h"
#endif

	fvinomgeneralolap_CallerWrapper(newndim, A,  B,ilimit,olimit,  ablockI,bblockI
			,numBlocks, params[2], (ilimit+pad) * olimit *sizeof(type)
			, d_lda_s,d_ldb_s,d_idx_s
			,acoars,bcoars,idx_s+c, (ilimit+pad), input_base, output_base, tile_base1, tile_base2, ilimitr, olimitr, irem, orem, ilimit*olimit, size, alpha, beta);

#ifdef NOHTIME
#include "includes/nohtimestop.h"
#endif


	{cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess){
			printf("\nKernel ERROR in fvi_nomatch_generalolap: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);
			//exit(-1);
		}}

	free(aexpr);
	free(bexpr);
	free(texpr1);
	free(texpr2);

	cudaFree(d_lda_s);
	cudaFree(d_ldb_s);
	cudaFree(d_idx_s);
	cudaFree(input_base);
	cudaFree(output_base);
	cudaFree(tile_base1);
	cudaFree(tile_base2);
}



