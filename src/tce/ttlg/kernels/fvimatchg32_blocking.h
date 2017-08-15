__global__ void CON(fvimatchg32_blocking_kernel_, ndim)  (const type * A, type * B,const int size0, const int size1, const int size2, const  int param0
		, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
		,const int lda_kernel1, const int ldb_kernel1, type alpha, type beta)
{
	const int block=param0;

	 int i=1;
        int aexpr = 0, bexpr = 0;
        int idx = blockIdx.x;
        int iip2=-1, ii1=-1;
        int  tmp = idx/idx_s[i];
        int index = idx - tmp * idx_s[i];
      //  aexpr += index * lda_s[i];
       // bexpr += index * ldb_s[i];
        idx = tmp;
        ii1 = index;
        i++;
        tmp = idx/idx_s[i];
        index = idx - tmp * idx_s[i];
       // aexpr += index * lda_s[i];
        //bexpr += index * ldb_s[i];
        idx = tmp;
        iip2 = index;

#pragma unroll
	//printf("%d\n ", idx);
	for(i = 3; i < ndim; i++)
	{
		tmp = idx/idx_s[i];
		int index = idx - tmp * idx_s[i];
		idx = tmp;
		aexpr += index * lda_s[i];
		bexpr += index * ldb_s[i];
	}


	int j=0; 
	int itemp = ii1*block; 
	int m0 = threadIdx.x/size0; 
	int i0 = threadIdx.x % size0; 
	int jump = 128/size0;
#pragma unroll 4 
	for(int i1 = iip2*block; i1<size1 && j<block; i1++,j++) { 
		int m=m0; 
		for(int i2 = (itemp+m); i2<size2 && m<block; m+=jump,i2+=jump){ 
			B[bexpr + i0 + i2*size0 + i1*ldb_kernel1] = alpha*A[aexpr + i0 + i1*size0 + i2*lda_kernel1] + beta*B[bexpr + i0 + i2*size0 + i1*ldb_kernel1]; 
		} 
	}
}
#undef ndim
