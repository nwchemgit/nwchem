__global__ void CON(fvinomatchg32_kernel_, ndim)  (const type *   A, type * B,
		const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
		, const int inputrem, const int outputrem, const int lda_kernel1, const int ldb_kernel1, const int remainder1, const int remainder2)
{
	int i=0;
	int aexpr = 0, bexpr = 0;	
	int idx = blockIdx.x;
	int iip2=-1, ii1=-1;
	int  tmp = idx/idx_s[i];
	int index = idx - tmp * idx_s[i];
	aexpr += index * lda_s[i];
	bexpr += index * ldb_s[i];
	idx = tmp;
	ii1 = index;
	i++;
	tmp = idx/idx_s[i];
	index = idx - tmp * idx_s[i];
	aexpr += index * lda_s[i];
	bexpr += index * ldb_s[i];
	idx = tmp;
	iip2 = index;
	for(i = 2; i < ndim; i++)
	{

		int  tmp = idx/idx_s[i];
		int index = idx - tmp * idx_s[i];
		aexpr += index * lda_s[i];
		bexpr += index * ldb_s[i];
		idx = tmp;
	}
	const double *Atmp = A + aexpr;
	double *Btmp = B + bexpr;
	if(ii1 < inputrem && iip2 < outputrem)
//	if(1 ||(ii1 < inputrem && iip2 < outputrem))
	{
		fvinomatchg32_main(Atmp,Btmp, lda_kernel1, ldb_kernel1);
	}
	else if(ii1 >= inputrem && iip2 < outputrem)
	{ 
		fvinomatchg32_rem(Atmp,Btmp, lda_kernel1, ldb_kernel1 ,remainder1, 32);
	}
	else if(iip2 >= outputrem && ii1 < inputrem)
	{ 
		fvinomatchg32_rem(Atmp,Btmp, lda_kernel1, ldb_kernel1, 32, remainder2);
	}
	else
	{
		fvinomatchg32_rem(Atmp,Btmp,lda_kernel1, ldb_kernel1,remainder1,remainder2);
	}

	return;

}
#undef ndim

