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
__device__  void fvinomgeneral_main(const type * __restrict__ Atmp, type * __restrict__ Btmp,  const int tb_size,  const int* __restrict__ aexpr, const int* __restrict__ bexpr, const int ilimit, const int olimit, type alpha, type beta) 
{

	const int TPR = tb_size/32;
	int in_s_colId = threadIdx.x % 32;
	int in_s_rowId = threadIdx.x / 32;
#ifdef printd	
	if(blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 && threadIdx.x == 0)
	{
		printf("ilimit = %d, olimit = %d, TPR = %d\n", ilimit, olimit, TPR);
		printf("\n%d %d\n\n\n", aexpr[10], bexpr[10]);
		printf("\nAtmp = %p, Btmp = %p\n\n\n", Atmp, Btmp);
	}
#endif
	for(int rowBatchId=0; rowBatchId < (olimit+31)/32; rowBatchId++)
	{
		int in_g_rowId = rowBatchId*32 + threadIdx.x / 32;
		int out_g_colId = rowBatchId * 32 + threadIdx.x % 32;
		for(int colBatchId=0; colBatchId < (ilimit+31)/32; colBatchId++)
		{
			int in_g_colId = colBatchId * 32 + threadIdx.x % 32;
			if(in_g_colId < ilimit)
			{
				for(int cur_row = in_g_rowId, local_r=0; local_r < 32 && cur_row < olimit && cur_row < (rowBatchId+1)* 32; local_r+=TPR, cur_row+=TPR)
					//for(int cur_row = in_g_rowId, local_r=0; local_r < 32 && cur_row < olimit; local_r+=TPR, cur_row+=TPR)
				{
					/*	if(blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 && threadIdx.x == 0)
						{
						printf("%d %d\n",cur_row, aexpr[cur_row] + in_g_colId);
						}*/
					tile[ (in_s_rowId + local_r) * 33 + in_s_colId] = Atmp[aexpr[cur_row] + in_g_colId];
				}
			}

			__syncthreads();

			int out_s_colId = in_s_rowId;
			int out_g_rowId = colBatchId*32 + threadIdx.x / 32;
			int out_s_rowId = in_s_colId;

			if(out_g_colId < olimit)
			{
				for(int cur_row = out_g_rowId, local_c=0; local_c < 32 && cur_row <ilimit && cur_row < (colBatchId+1)*32; local_c+=TPR, cur_row += TPR)
					//for(int cur_row = out_g_rowId, local_c=0; local_c < 32 && cur_row <ilimit; local_c+=TPR, cur_row += TPR)
				{
					//	if(blockIdx.x == 1 && blockIdx.y == 1 && blockIdx.z == 1 && threadIdx.x == 1)
					{
						//		printf("%p %p %d %d\n",Atmp, Btmp, cur_row, bexpr[cur_row] + out_g_colId);
					}
					Btmp[bexpr[cur_row] + out_g_colId ]  = alpha* tile[out_s_rowId * 33 + out_s_colId+ local_c] + beta* Btmp[bexpr[cur_row] + out_g_colId ];
				}
			}       
			__syncthreads();
			//	break;
		}
	}

}

#define FNAME fvinomatchgeneral.h
#include "includes/macro.h"
#undef FNAME





void fvinomatchgeneral_kernel_CallerWrapper(int ndim, type *  A, type * B, const int ilimit, const int olimit, const int blockAI, const int blockBI, const int numblocks,  int numthreads,  int shm
		, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
		, const int remainder1,		const int remainder2, const int * idx_ss,
		const int* __restrict__ aexpr, const int* __restrict__ bexpr, const int ilimitr, const int olimitr,
		const int inputrem, const int outputrem,
		type alpha, type beta
		)

{

	/*int second, third;
	  if(ndim > 2)
	  {
	  second = idx_ss[1]; third = numblocks/(idx_ss[0]*idx_ss[1]);
	  }
	  else if(ndim > 1)
	  {
	  second = idx_ss[1]; third = 1;
	  }
	  else
	  {
	  second = third = 1;
	  }
	  dim3 thread_blocks(idx_ss[0], second, third);*/
	dim3 thread_blocks(numblocks, 1, 1);

	switch(ndim)
	{
		EXPANDDIMS(fvinomgeneral_kernel_, thread_blocks, numthreads, shm, (A,  B, ilimit, olimit, blockAI, blockBI, numthreads, lda_s,ldb_s, idx_s, remainder1,remainder2, aexpr, bexpr, ilimitr, olimitr, inputrem, outputrem, alpha, beta))
		default: {}
			 //                       fvinomgeneralolap_kernel<<<thread_blocks, numthreads, shm>>>(ndim, A,  B, ilimit, olimit, i_blkindex, b_blkindex, numthreads, lda_s,ldb_s, idx_s, remainder1,remainder2, aexpr, bexpr, texpr1, texpr2, ilimitr, olimitr, inputrem, outputrem);
	}


}

void swap(int array[], int ind1, int ind2)
{
	if(ind1 == ind2) return;
	int tmp = array[ind1];
	array[ind1] = array[ind2];
	array[ind2] = tmp;
}
	extern "C" 
void  fvinomatchgeneral_transpose_kernel(int ndim, type *A, type *B,  const int *lda, const int *ldb, const int* params, const int * perm, const int* rperm, type alpha, type beta) 
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
	int olimit = params[8];
#ifdef SLICE
printf("\t%d\t%d\t", ilimit, olimit);
#endif
	int i = 0;
	//printf("blockA = %d, blockB = %d\n",blockA, blockB);
	//for(int y = 0; y < j; y++)
	//printf("bo[%d] = %d ",y, bo[y]);

	//exit(0);
	int numBlocks = params[6];//((size[1] + 8 -1)/8) * size[2] * ((size[3] + 8 -1)/8) * size[4] ; 

	int *d_lda_s, *d_ldb_s,  *d_idx_s; 
	const int remainder1 = lda[params[3]] % blockA;
	const int remainder2 = lda[perm[params[4]]] % blockB;
	const int ilimitr = ilimit * remainder1 / blockA; 
	const int olimitr = olimit * remainder2 / blockB; 
#ifdef MODEL
	printf("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", ilimit, olimit, remainder1, remainder2,ilimitr, olimitr, lda[params[3]] / blockA, lda[perm[params[4]]] / blockB);
#endif
#ifdef printd
	printf("\nrem1 = %d, rem2 = %d\n", remainder1, remainder2);
	printf("\nilimit = %d, olimit = %d", ilimit, olimit);
#endif

	int *input_base, *output_base;
	int *aexpr, *bexpr;

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
	if(blockA == 1)
	{
		idx_s[alimit] = 1;//(lda[i] + blockA - 1) / blockA;
	}
	else
	{
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
			if(blockB == 1)
			{
				idx_s[i] = 1;
			}
			else
			{
				idx_s[i] = (lda[i] + blockB - 1) / blockB;
			}
		}
		else
		{
			idx_s[i] = lda[i];
		}
	}

	for(i = 0; i < ndim; i++)
	{
		temp[i] = ldb_s[rperm[i]];
	//			printf("temp[%d] = %d\n", i, temp[i]);
#ifdef printd
				printf("idx[%d] = %d\n", i, idx_s[i]);
#endif
	}


	int irem, orem;
	if(remainder1 == 0) irem = lda[alimit];
	else irem = (lda[alimit] - remainder1)/blockA;
	if(remainder2 == 0) orem = ldb[blimit];
	else orem = (ldb[blimit] - remainder2)/blockB;


	aexpr = (int*)malloc(olimit * sizeof(int));
	bexpr = (int*)malloc(ilimit * sizeof(int));


	SAFECUDAMALLOC(&input_base, olimit*sizeof(int)); 
	SAFECUDAMALLOC(&output_base, ilimit*sizeof(int)); 


	const int TPR = params[2]/32;
	for(int i = 0; i < params[2]; i++)
	{
		for(int rowBatchId=0; rowBatchId < (olimit+31)/32; rowBatchId++)
		{
			for(int colBatchId=0; colBatchId < (ilimit + 31)/32; colBatchId++)
			{
				int in_g_colId = colBatchId * 32 + i % 32;
				int in_g_rowId = rowBatchId*32 + i / 32;
				int cur_row = in_g_rowId;
				if(in_g_colId < ilimit)
				{
					for(int local_r=0; local_r < 32 && cur_row < olimit &&cur_row < in_g_rowId + 32; local_r++, cur_row+=TPR)
					{
						int tmp = cur_row;
						int ii[20];int aoff=0,j;
						for(j = 0; j < blimit; j++)
						{
							ii[j] = tmp%ldb[j];
							tmp /= ldb[j];//tmp/bo[j];
							aoff += ii[j]* lda_s[perm[j]];
							//printf(" j = %d\t aoff = %d\n", j, aoff);	
						}
						aoff += (tmp)* lda_s[perm[j]];
						aexpr[cur_row] = aoff;// i2*lda0 + i1 * lda1;
					}
				}
				//exit(0);
				int out_g_colId = rowBatchId * 32 + i % 32;
				int out_g_rowId = colBatchId*32 + i / 32;
				if(out_g_colId < olimit)
				{
					for(int cur_row = out_g_rowId, local_c=0; local_c < 32 && cur_row < ilimit && cur_row < out_g_rowId + 32; local_c++, cur_row +=TPR)
					{
						int tmp = cur_row;
						int ii[20];int boff=0,j;
						for(j = 0; j < alimit; j++)
						{
							ii[j] = tmp%lda[j];
							tmp = tmp/lda[j];
							boff += ii[j]* ldb_s[rperm[j]];
							//printf(" j = %d\t boff = %d\n", j, boff);	
						}
						boff += tmp* ldb_s[rperm[j]];
						bexpr[cur_row] = boff;
						//						printf("currow = %d, b = %d\n", cur_row, boff);

					}
				}
			}
		}
	}

#ifdef printd
	printf("\nA..\n");
	for(int i = 0; i < olimit; i++)
	{
		printf("%d ", aexpr[i]);
	}
	printf("\n");
	printf("\nB..\n");

	for(int i = 0; i < ilimit; i++)
	{
		printf("%d ", bexpr[i]);
	}
	printf("\n");
#endif

	lda_s[params[3]] *= params[0];
	temp[params[3]] *= params[0];
	lda_s[perm[params[4]]] *= params[11];
	temp[perm[params[4]]] *= params[11];

	//Lets remove unwanted dimensions for thread block indexing

	int c = 0, d = 0;
	c = alimit + 1;//c = No. of dimensions to be removed from input for thread blocking
	if(blockA > 1) c--;
	int ablockI, bblockI;
	ablockI = alimit-c;
	bblockI = perm[blimit]-c;
	int tempbblockI = bblockI;
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
		//	printf("\n i = %d\n", i);
			if((i < bblockI + c)  || (i == bblockI + c) && (blockB == 1))
			{
				tempbblockI--;
			}
		}
	}
	bblockI = tempbblockI;
#ifdef printd
	printf("\nd = %d, bblockI_changed = %d\n", d, bblockI);
#endif
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
	const int newndim = ndim - (c + d);
#ifdef printd
	for(i = c; i < ndim-d; i++)
	{
				printf("idx[%d] = %d\n", i, idx_s[i]);
	}


	printf("ndim = %d, c = %d, d = %d, newndim = %d, ablockI = %d, bblockI = %d\n", ndim, c, d, newndim, ablockI, bblockI);
#endif
	//Find the largest dimension and make it the first as only Dimx can have > 65k size
	/*int max = 0;
	  for(int i = 1; i < newndim; i++)
	  {
	  if(idx_s[c+i] > idx_s[c+max]) max = i;
	  }
	//printf("\nmax: %d ", max);
	swap(idx_s, c, max+c);
	swap(lda_s, c, max+c);
	swap(temp, c, max+c);

	if(max == ablockI) ablockI = 0;
	else if(ablockI == 0) ablockI = max;
	if(max == bblockI) bblockI = 0;
	else if(bblockI ==0) bblockI = max;
	 */

	if(blockB == 1) bblockI = -1;
	if(ablockI > 0)//move it to first//shouldnt happen
	{
		swap(idx_s, ablockI+c, c);
		swap(lda_s, ablockI+c, c);
		swap(temp, ablockI+c, c);
	}
	
	if(bblockI >= 0)//move it to second
	{
		if(bblockI != 0 || ablockI < 0)
		{
			swap(idx_s, bblockI+c, c+ 1);
			swap(lda_s, bblockI+c,  c+ 1);
			swap(temp, bblockI+c,  c+1);
		}
		else
		{
			swap(idx_s, ablockI+c, c+ 1);
			swap(lda_s, ablockI+c, c+ 1);
			swap(temp, ablockI+c, c+ 1);
		}
	}


	if(bblockI >= 0) {
		// if(ablockI < 0) bblockI = 0;
		bblockI = 1;
	}

	if(ablockI >=0) ablockI = 0; 
#ifdef printd
	printf("\nIDx: ");
	for(int i = 0; i < newndim; i++)
	{
		printf("%d ",idx_s[i+c]);
	}
	printf("ndim = %d, c = %d, d = %d, newndim = %d, ablockI = %d, bblockI = %d\n", ndim, c, d, newndim, ablockI, bblockI);
#endif


	SAFECUDAMALLOC(&d_lda_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_ldb_s,newndim*sizeof(int)); 
	SAFECUDAMALLOC(&d_idx_s,newndim*sizeof(int)); 

	SAFECUDAMEMCPY(d_idx_s, idx_s+c,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_lda_s, lda_s+c,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(d_ldb_s, temp+c,newndim*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(input_base, aexpr, olimit*sizeof(int), cudaMemcpyHostToDevice); 
	SAFECUDAMEMCPY(output_base, bexpr, ilimit*sizeof(int), cudaMemcpyHostToDevice); 

#ifdef MODEL
	printf("\tilimit=%d\tolimit=%d\t", ilimit, olimit);
	printf("\t%d\t%d\t%d\t%d\t", ilimit/32, ilimit%32, olimit/32,olimit%32 );
	double f1, f2, f3, f4, f;
	printf("\tf1=%lf\t", f1 =  ((ilimit/32) * (olimit/32) + (double)(ilimit/32) * (olimit%32) /32+ (double)(ilimit%32) * (olimit/32) /32 + (double)(ilimit%32) * (olimit%32) /(32*32) )/ (int)(((ilimit+31)/32) * ((olimit+31)/32)));
	printf("\tf2=%lf\t", f2 =  ((ilimitr/32) * (olimit/32) + (double)(ilimitr/32) * (olimit%32) /32+ (double)(ilimitr%32) * (olimit/32) /32 + (double)(ilimitr%32) * (olimit%32) /(32*32) )/ max(1,(int)(((ilimitr+31)/32) * ((olimit+31)/32))));
	printf("\tf3=%lf\t", f3 =  ((ilimit/32) * (olimitr/32) + (double)(ilimit/32) * (olimitr%32) /32+ (double)(ilimit%32) * (olimitr/32) /32 + (double)(ilimit%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((ilimit+31)/32) * ((olimitr+31)/32))));
	printf("\tf4=%lf\t", f4 =  ((ilimitr/32) * (olimitr/32) + (double)(ilimitr/32) * (olimitr%32) /32+ (double)(ilimitr%32) * (olimitr/32) /32 + (double)(ilimitr%32) * (olimitr%32) /(32*32) )/ max(1,(int)(((ilimitr+31)/32) * ((olimitr+31)/32))));
	printf("\t%d\t%d\t%d\t%d\t", lda[alimit], ldb[blimit],blockA,blockB);
	int asize = lda[alimit];
	int bsize = ldb[blimit];
	printf("\t%d\t%d\t%d\t%d\t", asize/blockA, asize%blockA, bsize/blockB,bsize%blockB );
	//int amax = min(blockA, 32);
	//int bmax = min(blockB, 32);
	int amax = blockA;
	int bmax = blockB;
	printf("\tf=%lf\t", f = ((asize/amax) * (bsize/bmax) *f1 + (double)((asize/amax) * (bsize%bmax > 0) *f3)+ (double)((asize%amax > 0) * (bsize/bmax)*f2)  + (double)((asize%amax>0) * (bsize%bmax > 0) *f4) )/ (int)(((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax)));
	//printf("\tg=%lf\t%lf\t%lf\t%lf, den=%d\t", (asize/amax) * (bsize/bmax) *f1 , (double)((asize/amax) * (bsize%bmax) *f3)/bmax, (double)((asize%amax) * (bsize/bmax)*f2)/amax , (double)((asize%amax) * (bsize%bmax) *f4)/(amax*bmax),  (((asize+amax-1)/amax) * ((bsize+bmax-1)/bmax)));
	printf("\t%lf\t", f);
#endif

#ifdef NOHTIME
#include "includes/nohtimestart.h"
#endif

	//fvinomgeneral_kernel_CallerWrapper(newndim, A,  B,ilimit,olimit, ablockI, bblockI
	fvinomatchgeneral_kernel_CallerWrapper(newndim, A,  B,ilimit,olimit, ablockI, bblockI
			,numBlocks, params[2],params[10]* params[5]*sizeof(type)
			, d_lda_s,d_ldb_s,d_idx_s
			,remainder1,remainder2,idx_s+c,  input_base, output_base, ilimitr, olimitr, irem, orem, alpha, beta);

#ifdef NOHTIME
#include "includes/nohtimestop.h"
#endif

	{cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess){
			printf("\nKernel ERROR in fvi_nomatch_general: %s (line: %d)\n", cudaGetErrorString(err), __LINE__);
			//exit(-1);
		}}

	free(aexpr);
	free(bexpr);

	cudaFree(d_lda_s);
	cudaFree(d_ldb_s);
	cudaFree(d_idx_s);
	cudaFree(input_base);
	cudaFree(output_base);
	//cudaFree(d_ablock);
	//cudaFree(d_bblock);
	//#endif
}



