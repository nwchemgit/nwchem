__global__ void CON(fvimatchl32_coars_kernel_,  ndim)  (const type *  A, type * B,const int  size0, const int ldal,const int ldap2l ,const  int param0, const int param1, const int plain, const int sbp
, const int * __restrict__   lda_s, const int * __restrict__   ldb_s, const int * __restrict__    idx_s
, const int remainder1,         const int remainder2, const int ilimit, const int olimit, const int lda_kernel1, const int ldb_kernel1, const unsigned short* __restrict__ offset, const int acoars, const int bcoars, const int size, type alpha, type beta)
{
	const int blockA = param0;
	const int blockB = param0;
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


	if(ii1 < ldal && iip2 < ldap2l)
        {
                fvimatchl32_main_coars(Atmp,Btmp, lda_kernel1, ldb_kernel1,size0, plain, sbp, offset, acoars, bcoars, size, alpha, beta);
        }
        else if(ii1 >= ldal && iip2 < ldap2l)
        { 
                fvimatchl32_rem_coars(Atmp,Btmp, lda_kernel1, ldb_kernel1,remainder1, blockB, size0,  plain, sbp, ilimit, plain, offset, acoars, bcoars, size, alpha, beta);
         }
        else if(iip2 >= ldap2l && ii1 < ldal)
        { 
                fvimatchl32_rem_coars(Atmp,Btmp, lda_kernel1, ldb_kernel1,blockA, remainder2,size0, plain, sbp, plain, olimit, offset, acoars, bcoars, size, alpha, beta);
        }
        else
        {
                fvimatchl32_rem_coars(Atmp,Btmp,lda_kernel1, ldb_kernel1,remainder1,remainder2, size0, plain, sbp, ilimit,olimit, offset, acoars, bcoars, size, alpha, beta);
        }

        return;
}
#undef ndim
                         
