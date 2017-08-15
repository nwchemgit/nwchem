__global__ void CON(fvinomatchg32_coars_kernel_, ndim)  (const type *   A, type * B,
		const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
		,const int acoars, const int bcoars
		, const int inputrem, const int outputrem, const int lda_kernel1, const int ldb_kernel1, const int remainder1, const int remainder2, const int size)
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
	{
		fvinomatchg32_main_coars(Atmp,Btmp, lda_kernel1, ldb_kernel1, acoars, bcoars, size );
	}
	else if(ii1 >= inputrem && iip2 < outputrem)
	{ 
		fvinomatchg32_rem_coars(Atmp,Btmp, lda_kernel1, ldb_kernel1 ,remainder1, 32, acoars, bcoars, size);
	}
	else if(iip2 >= outputrem && ii1 < inputrem)
	{ 
		fvinomatchg32_rem_coars(Atmp,Btmp, lda_kernel1, ldb_kernel1, 32, remainder2, acoars, bcoars, size);
	}
	else
	{
		fvinomatchg32_rem_coars(Atmp,Btmp,lda_kernel1, ldb_kernel1,remainder1,remainder2,acoars, bcoars, size);
	}

	return;

}
#undef ndim

