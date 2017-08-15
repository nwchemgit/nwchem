__global__ void CON(fvinomgeneral_kernel_,  ndim)  (const type *  A, type * B, const int ilimit, const int olimit, const int blockAI, const int blockBI, const int num_threads
		, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
		, const int remainder1,         const int remainder2, 
		const int* __restrict__ offseti, const int* __restrict__ offseto,
		const int ilimitr, const int olimitr, const int inputrem, const int outputrem,
		type alpha, type beta
		)
{
	//int ii[ndim];
	//ii[0] = blockIdx.x;
	//ii[1] = blockA*blockIdx.x;
	//ii[1] = blockIdx.y;
	int aexpr = 0, bexpr = 0;// = ii[0] * lda_s[0] + ii[1]*lda_s[1], bexpr = ii[0] * ldb_s[0] + ii[1] * ldb_s[1];
	int idx = blockIdx.x;
	int iip2=-1, ii1=-1, i=0;
	{
		int tmp = idx/idx_s[0];
		int index = idx - tmp * idx_s[0];
		aexpr += index * lda_s[0];
		bexpr += index * ldb_s[0];
		idx = tmp;
		if(blockAI == 0)
			ii1 = index;
		//	else if(blockBI == 0)
		//		iip2 = index;
	}
	if(ndim > 1)
	{
		int tmp = idx/idx_s[1];
		int index = idx - tmp * idx_s[1];
		aexpr += index * lda_s[1];
		bexpr += index * ldb_s[1];
		//if(blockAI == 1)
		//	ii1 = index;
		//else
		if(blockBI == 1)
			iip2 = index;
		idx = tmp;
	}
	
	for(i = 2; i < ndim; i++)
	{
		int tmp = idx/idx_s[i];
		//ii[i] = idx - tmp * idx_s[i];
		int index = idx - tmp * idx_s[i];
		//int index = ii[i];
		aexpr += index * lda_s[i];
		bexpr += index * ldb_s[i];
		idx = tmp;
	}
	//int iip2 = ii[2], ii1 = ii[1];
	const double *Atmp = A + aexpr;
	double *Btmp = B + bexpr;
	if((ii1 < inputrem && iip2 < outputrem))
	{
		fvinomgeneral_main(Atmp,Btmp,  num_threads,  offseti, offseto, ilimit, olimit, alpha, beta );
	}
	else if(ii1 >= inputrem && iip2 < outputrem)
	{ //remainder in size1
		fvinomgeneral_main(Atmp,Btmp,  num_threads,  offseti, offseto, ilimitr, olimit, alpha, beta );
	}
	else if(iip2 >= outputrem && ii1 < inputrem)
	{ //remainder in size2 
		fvinomgeneral_main(Atmp,Btmp,  num_threads,  offseti, offseto, ilimit, olimitr, alpha, beta );
	}
	else
	{       
		fvinomgeneral_main(Atmp,Btmp,  num_threads,  offseti, offseto, ilimitr, olimitr, alpha, beta );
	}


	return;
}
#undef ndim

