#if defined(MKL)
#include <mkl.h>
#ifdef MKL_ILP64
#error Use the MKL library for 32-bit integers!
#endif
#elif defined(ACCELERATE)
/* The location of cblas.h is not in the system include path when -framework Accelerate is provided. */
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

void ccsd_trpdrv_omp_body_(double * restrict f1n, double * restrict f1t,
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

    const int ncor = *ncor_;
    const int nocc = *nocc_;
    const int nvir = *nvir_;

    /* convert from Fortran to C offset convention... */
    const int a   = *a_ - 1;
    const int i   = *i_ - 1;
    const int j   = *j_ - 1;
    const int k   = *k_ - 1;
    const int klo = *klo_ - 1;

    double emp4i = 0.0;
    double emp5i = 0.0;
    double emp4k = 0.0;
    double emp5k = 0.0;

#pragma omp parallel shared(eorb) \
            shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t) \
            shared(t1v1,dintc1,dintx1) \
            shared(t1v2,dintc2,dintx2) \
            firstprivate(ncor,nocc,nvir,i,j,k,klo)
   {

        /* Performance Note:
           By definition, the following does not scale to more than 8 threads
           unless nested parallelism (i.e. inside of DGEMM) is used.
           It may be prudent to write a manually threaded wrapper for the
           cases where single-threaded BLAS is used. */

        const int lnov = nocc * nvir;
        const int lnvv = nvir * nvir;

        #pragma omp sections
        {
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasTrans,
                            nvir, nvir, nvir, 1.0, jia, nvir, &tkj[(k-klo)*lnvv], nvir, 0.0, f1n, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, tia, nvir, &kkj[(k-klo)*lnov], nocc, 1.0, f1n, nvir);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasTrans,
                            nvir, nvir, nvir, 1.0, kia, nvir, &tkj[(k-klo)*lnvv], nvir, 0.0, f2n, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, xia, nvir, &kkj[(k-klo)*lnov], nocc, 1.0, f2n, nvir);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nvir, 1.0, jia, nvir, &tkj[(k-klo)*lnvv], nvir, 0.0, f3n, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, tia, nvir, &jkj[(k-klo)*lnov], nocc, 1.0, f3n, nvir);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nvir, 1.0, kia, nvir, &tkj[(k-klo)*lnvv], nvir, 0.0, f4n, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, xia, nvir, &jkj[(k-klo)*lnov], nocc, 1.0, f4n, nvir);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasTrans,
                            nvir, nvir, nvir, 1.0, &jka[(k-klo)*lnvv], nvir, tij, nvir, 0.0, f1t, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, &tka[(k-klo)*lnov], nvir, kij, nocc, 1.0, f1t, nvir);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasTrans,
                            nvir, nvir, nvir, 1.0, &kka[(k-klo)*lnvv], nvir, tij, nvir, 0.0, f2t, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, &xka[(k-klo)*lnov], nvir, kij, nocc, 1.0, f2t, nvir);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nvir, 1.0, &jka[(k-klo)*lnvv], nvir, tij, nvir, 0.0, f3t, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, &tka[(k-klo)*lnov], nvir, jij, nocc, 1.0, f3t, nvir);
            }
            #pragma omp section
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nvir, 1.0, &kka[(k-klo)*lnvv], nvir, tij, nvir, 0.0, f4t, nvir);
                cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
                            nvir, nvir, nocc, -1.0, &xka[(k-klo)*lnov], nvir, jij, nocc, 1.0, f4t, nvir);
            }
        }

        const double eaijk = eorb[a] - (eorb[ncor+i] + eorb[ncor+j] + eorb[ncor+k]);

        /* b,c loop over [0,nvir) to eliminate the need for offset corrections... */
        #pragma omp for collapse(2) schedule(static) reduction(+:emp5i,emp4i) reduction(+:emp5k,emp4k)
        for (int b = 0; b < nvir; ++b) {
            for (int c = 0; c < nvir; ++c) {
                /* WARNING: Do not add IVDEP here.  Code will be incorrect. */
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
    if (i != k) {
        emp4 += emp4k;
        emp5 += emp5k;
    }

    *emp4_ = emp4;
    *emp5_ = emp5;

    return;
}

