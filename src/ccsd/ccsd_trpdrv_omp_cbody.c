#include <stdio.h>

#if defined(MKL)
# include <mkl.h>
# ifndef MKL_INT
#  error MKL_INT not defined!
# endif
  typedef MKL_INT cblas_int;
#elif defined(ACCELERATE)
  /* The location of cblas.h is not in the system include path when -framework Accelerate is provided. */
# include <Accelerate/Accelerate.h>
  typedef int cblas_int;
#else
# include <cblas.h>
  typedef int cblas_int;
#endif

void ccsd_trpdrv_omp_cbody_(double * restrict f1n, double * restrict f1t,
                            double * restrict f2n, double * restrict f2t,
                            double * restrict f3n, double * restrict f3t,
                            double * restrict f4n, double * restrict f4t,
                            double * restrict eorb,
                            int    * restrict ncor_, int * restrict nocc_, int * restrict nvir_,
                            double * restrict emp4_, double * restrict emp5_,
                            int    * restrict a_, int * restrict i_, int * restrict j_, int * restrict k_, int * restrict klo_,
                            double * restrict tij, double * restrict tkj, double * restrict tia, double * restrict tka,
                            double * restrict xia, double * restrict xka, double * restrict jia, double * restrict jka,
                            double * restrict kia, double * restrict kka, double * restrict jij, double * restrict jkj,
                            double * restrict kij, double * restrict kkj,
                            double * restrict dintc1, double * restrict dintx1, double * restrict t1v1,
                            double * restrict dintc2, double * restrict dintx2, double * restrict t1v2)
{
    double emp4 = *emp4_;
    double emp5 = *emp5_;

    double emp4i = 0.0;
    double emp5i = 0.0;
    double emp4k = 0.0;
    double emp5k = 0.0;

#pragma omp parallel shared(eorb,f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t) \
                     shared(t1v1,dintc1,dintx1,t1v2,dintc2,dintx2)
   {

        const int ncor = *ncor_;
        const int nocc = *nocc_;
        const int nvir = *nvir_;

        const int lnov = nocc * nvir;
        const int lnvv = nvir * nvir;

        /* convert from Fortran to C offset convention... */
        const int k   = *k_ - 1;
        const int klo = *klo_ - 1;

        /* Performance Note:
           By definition, the following does not scale to more than 8 threads
           unless nested parallelism (i.e. inside of DGEMM) is used.
           It may be prudent to write a manually threaded wrapper for the
           cases where single-threaded BLAS is used. */

        const cblas_int nv = nvir;
        const cblas_int no = nocc;

        #pragma omp sections
        {
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                            nv, nv, nv, 1.0, jia, nv, &tkj[(k-klo)*lnvv], nv, 0.0, f1n, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, tia, nv, &kkj[(k-klo)*lnov], no, 1.0, f1n, nv);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                            nv, nv, nv, 1.0, kia, nv, &tkj[(k-klo)*lnvv], nv, 0.0, f2n, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, xia, nv, &kkj[(k-klo)*lnov], no, 1.0, f2n, nv);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, nv, 1.0, jia, nv, &tkj[(k-klo)*lnvv], nv, 0.0, f3n, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, tia, nv, &jkj[(k-klo)*lnov], no, 1.0, f3n, nv);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, nv, 1.0, kia, nv, &tkj[(k-klo)*lnvv], nv, 0.0, f4n, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, xia, nv, &jkj[(k-klo)*lnov], no, 1.0, f4n, nv);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                            nv, nv, nv, 1.0, &jka[(k-klo)*lnvv], nv, tij, nv, 0.0, f1t, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, &tka[(k-klo)*lnov], nv, kij, no, 1.0, f1t, nv);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                            nv, nv, nv, 1.0, &kka[(k-klo)*lnvv], nv, tij, nv, 0.0, f2t, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, &xka[(k-klo)*lnov], nv, kij, no, 1.0, f2t, nv);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, nv, 1.0, &jka[(k-klo)*lnvv], nv, tij, nv, 0.0, f3t, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, &tka[(k-klo)*lnov], nv, jij, no, 1.0, f3t, nv);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, nv, 1.0, &kka[(k-klo)*lnvv], nv, tij, nv, 0.0, f4t, nv);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            nv, nv, no, -1.0, &xka[(k-klo)*lnov], nv, jij, no, 1.0, f4t, nv);
            }
        }

        /* convert from Fortran to C offset convention... */
        const int a   = *a_ - 1;
        const int i   = *i_ - 1;
        const int j   = *j_ - 1;

        const double eaijk = eorb[a] - (eorb[ncor+i] + eorb[ncor+j] + eorb[ncor+k]);

        /* b,c loop over [0,nvir) to eliminate the need for offset corrections... */
        #pragma omp for collapse(2) schedule(static) reduction(+:emp5i,emp4i) reduction(+:emp5k,emp4k)
        for (int b = 0; b < nvir; ++b) {
            for (int c = 0; c < nvir; ++c) {
                double denom = -1.0 / (eorb[ncor+nocc+b] + eorb[ncor+nocc+c] + eaijk);
                emp4i += denom * (f1t[b+c*nvir]+f1n[c+b*nvir]+f2t[c+b*nvir]+f3n[b+c*nvir]+f4n[c+b*nvir])
                               * (f1t[b+c*nvir]-f2t[b+c*nvir]*2-f3t[b+c*nvir]*2+f4t[b+c*nvir])
                       - denom * (f1n[b+c*nvir]+f1t[c+b*nvir]+f2n[c+b*nvir]+f3n[c+b*nvir])
                               * (f1t[b+c*nvir]*2-f2t[b+c*nvir]-f3t[b+c*nvir]+f4t[b+c*nvir]*2)
                       + denom * 3 * (f1n[b+c*nvir]*(f1n[b+c*nvir]+f3n[c+b*nvir]+f4t[c+b*nvir]*2)
                                      +f2n[b+c*nvir]*f2t[c+b*nvir]+f3n[b+c*nvir]*f4t[b+c*nvir]);
                emp4k += denom * (f1n[b+c*nvir]+f1t[c+b*nvir]+f2n[c+b*nvir]+f3t[b+c*nvir]+f4t[c+b*nvir])
                               * (f1n[b+c*nvir]-f2n[b+c*nvir]*2-f3n[b+c*nvir]*2+f4n[b+c*nvir])
                       - denom * (f1t[b+c*nvir]+f1n[c+b*nvir]+f2t[c+b*nvir]+f3t[c+b*nvir])
                               * (f1n[b+c*nvir]*2-f2n[b+c*nvir]-f3n[b+c*nvir]+f4n[b+c*nvir]*2)
                       + denom * 3 * (f1t[b+c*nvir]*(f1t[b+c*nvir]+f3t[c+b*nvir]+f4n[c+b*nvir]*2)
                                      +f2t[b+c*nvir]*f2n[c+b*nvir]+f3t[b+c*nvir]*f4n[b+c*nvir]);
                emp5i += denom * t1v1[b] * dintx1[c] * (f1t[b+c*nvir]+f2n[b+c*nvir]+f4n[c+b*nvir]
                                                        -(f3t[b+c*nvir]+f4n[b+c*nvir]+f2n[c+b*nvir]
                                                          +f1n[b+c*nvir]+f2t[b+c*nvir]+f3n[c+b*nvir])*2
                                                        +(f3n[b+c*nvir]+f4t[b+c*nvir]+f1n[c+b*nvir])*4)
                       + denom * t1v1[b] * dintc1[c] * (f1n[b+c*nvir]+f4n[b+c*nvir]+f1t[c+b*nvir]
                                                        -(f2n[b+c*nvir]+f3n[b+c*nvir]+f2t[c+b*nvir])*2);
                emp5k += denom * t1v2[b] * dintx2[c] * (f1n[b+c*nvir]+f2t[b+c*nvir]+f4t[c+b*nvir]
                                                        -(f3n[b+c*nvir]+f4t[b+c*nvir]+f2t[c+b*nvir]
                                                          +f1t[b+c*nvir]+f2n[b+c*nvir]+f3t[c+b*nvir])*2
                                                        +(f3t[b+c*nvir]+f4n[b+c*nvir]+f1t[c+b*nvir])*4)
                       + denom * t1v2[b] * dintc2[c] * (f1t[b+c*nvir]+f4t[b+c*nvir]+f1n[c+b*nvir]
                                                        -(f2t[b+c*nvir]+f3t[b+c*nvir]+f2n[c+b*nvir])*2);
            }
        }
    }

    emp4 += emp4i;
    emp5 += emp5i;

    if (*i_ != *k_) {
        emp4 += emp4k;
        emp5 += emp5k;
    }

    *emp4_ = emp4;
    *emp5_ = emp5;

    return;
}

