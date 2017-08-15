__global__ void CON(fvinomgeneralolap_coars_kernel_,ndim) (const type *  A, type * B, const int ilimit, const int olimit,const int blockAI,const int blockBI, const int numthreads
		, const int * __restrict__ lda_s, const int* __restrict__ ldb_s, const int* __restrict__ idx_s
		, const int acoars,         const int bcoars, 
		const int* __restrict__ offseti, const int* __restrict__ offseto, const int* __restrict__ tile1, const int* __restrict__ tile2,
	const int ilimitr, const int olimitr, const int inputbdmv, const int outputbdmv, const int rowinc, const int shm2, const int numelems_blk, const int size,
	type alpha, type beta
)
{
	int aexpr =0, bexpr = 0;
	int idx;
	idx = blockIdx.x;
	int inputbdv = -1 , outputbdv = -1; //ensures correctness when no remainder and some dimensions are missing from thread blocks
	
/*	if(blockAI == 0)
		inputbdv = val0;
	else if(blockAI == 1)
		inputbdv = val1;
	if(blockBI == 0)
		outputbdv = val0;
	else if(blockBI == 1)
		outputbdv = val1;*/
	int i=0;
	if(i<ndim){
		int tmp = idx/idx_s[i];
		int index = idx - tmp * idx_s[i];
		aexpr += index * lda_s[i];
		bexpr += index * ldb_s[i];
		idx = tmp;
		if(i == blockBI) outputbdv = index;
		if(i == blockAI) inputbdv = index;
		i++;
	}
	if(i<ndim){
		int tmp = idx/idx_s[i];
		int index = idx - tmp * idx_s[i];
		aexpr += index * lda_s[i];
		bexpr += index * ldb_s[i];
		idx = tmp;
		if(i == blockAI) inputbdv = index;
		if(i == blockBI) outputbdv = index;
		i++;
	}	
	#pragma unroll
	for(i = 2; i < ndim; i++)
	{
		int tmp = idx/idx_s[i];
		int index = idx - tmp * idx_s[i];
		aexpr += index * lda_s[i];
		bexpr += index * ldb_s[i];
		idx = tmp;
	}
	const double *Abase = A + aexpr;
	double *Bbase = B + bexpr;
	//if(blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 && threadIdx.x == 5)
	//	printf("inputbdv = %d, outputbdv = %d \n", inputbdv, outputbdv); 
	//if(1 ||(inputbdv < inputbdmv && outputbdv < outputbdmv))
	if( (inputbdv < inputbdmv && outputbdv < outputbdmv))
	{
		fvinomgeneralolap_main_coars(Abase,Bbase,  numthreads,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, numelems_blk, acoars, bcoars, size, alpha, beta );
	}
	 else if(inputbdv >= inputbdmv && outputbdv < outputbdmv)
        { //remainder in size1

		fvinomgeneralolap_rem_coars(Abase,Bbase,  numthreads,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, ilimitr, olimit, numelems_blk, acoars, bcoars, size, alpha, beta );
        }
        else if(outputbdv >= outputbdmv && inputbdv < inputbdmv)
        { //remainder in size2 
		fvinomgeneralolap_rem_coars(Abase,Bbase,  numthreads,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, ilimit, olimitr, numelems_blk, acoars, bcoars, size, alpha, beta );
        }
        else
        {       
		fvinomgeneralolap_rem_coars(Abase,Bbase,  numthreads,  offseti, offseto,tile1, tile2, ilimit, olimit, rowinc, shm2, ilimitr, olimitr, numelems_blk, acoars, bcoars, size, alpha, beta );
        }
	return;
}
#undef ndim
