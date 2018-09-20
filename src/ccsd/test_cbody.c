#include <stdio.h>
#include <stdlib.h>

/* Do not allow the test to allocate more than MAX_MEM gigabytes. */
#ifndef MAX_MEM
#define MAX_MEM 4
#endif

#define MIN(x,y) (x<y ? x : y)
#define MAX(x,y) (x>y ? x : y)

#ifdef _OPENMP
#include <omp.h>
#else
#warning No timer!
double omp_get_wtime() { return 0.0; }
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
                            double * restrict dintc2, double * restrict dintx2, double * restrict t1v2);

double * make_array(int n)
{
    double * a = malloc(n*sizeof(double));
    for (int i=0; i<n; i++) {
        a[i] = 1.0/(100.0+i);
    }
    return a;
}

int main(int argc, char* argv[])
{
    int ncor, nocc, nvir;
    int maxiter = 100;
    int nkpass = 1;

    if (argc<3) {
        printf("Usage: ./test_cbody nocc nvir [maxiter] [nkpass]\n");
        return argc;
    } else {
        ncor = 0;
        nocc = atoi(argv[1]);
        nvir = atoi(argv[2]);
        if (argc>3) {
            maxiter = atoi(argv[3]);
            /* if negative, treat as "infinite" */
            if (maxiter<0) maxiter = 1<<30;
        }
        if (argc>4) {
            nkpass = atoi(argv[4]);
        }
    }

    if (nocc<1 || nvir<1) {
        printf("Arguments must be non-negative!\n");
        return 1;
    }

    printf("Test driver for cbody with nocc=%d, nvir=%d, maxiter=%d, nkpass=%d\n", nocc, nvir, maxiter, nkpass);

    const int nbf = ncor + nocc + nvir;
    const int lnvv = nvir * nvir;
    const int lnov = nocc * nvir;
    const int kchunk = (nocc - 1)/nkpass + 1;

    const double memory = (nbf+8.0*lnvv+
                           lnvv+kchunk*lnvv+lnov*nocc+kchunk*lnov+lnov*nocc+kchunk*lnov+lnvv+
                           kchunk*lnvv+lnvv+kchunk*lnvv+lnov*nocc+kchunk*lnov+lnov*nocc+
                           kchunk*lnov+lnov+nvir*kchunk+nvir*nocc+
                           6.0*lnvv)*sizeof(double);
    printf("This test requires %f GB of memory.\n", 1.0e-9*memory);

    if (1.0e-9*memory > MAX_MEM) {
        printf("You need to increase MAX_MEM (%d)\n", MAX_MEM);
        printf("or set nkpass (%d) to a larger number.\n", nkpass);
        return MAX_MEM;
    }

    double * eorb = make_array(nbf);

    double * f1n = make_array(lnvv);
    double * f2n = make_array(lnvv);
    double * f3n = make_array(lnvv);
    double * f4n = make_array(lnvv);
    double * f1t = make_array(lnvv);
    double * f2t = make_array(lnvv);
    double * f3t = make_array(lnvv);
    double * f4t = make_array(lnvv);

    double * Tij  = make_array(lnvv);
    double * Tkj  = make_array(kchunk*lnvv);
    double * Tia  = make_array(lnov*nocc);
    double * Tka  = make_array(kchunk*lnov);
    double * Xia  = make_array(lnov*nocc);
    double * Xka  = make_array(kchunk*lnov);
    double * Jia  = make_array(lnvv);
    double * Jka  = make_array(kchunk*lnvv);
    double * Kia  = make_array(lnvv);
    double * Kka  = make_array(kchunk*lnvv);
    double * Jij  = make_array(lnov*nocc);
    double * Jkj  = make_array(kchunk*lnov);
    double * Kij  = make_array(lnov*nocc);
    double * Kkj  = make_array(kchunk*lnov);
    double * Dja  = make_array(lnov);
    double * Djka = make_array(nvir*kchunk);
    double * Djia = make_array(nvir*nocc);

    double * dintc1 = make_array(lnvv);
    double * dintc2 = make_array(lnvv);
    double * dintx1 = make_array(lnvv);
    double * dintx2 = make_array(lnvv);
    double * t1v1   = make_array(lnvv);
    double * t1v2   = make_array(lnvv);

    int ntimers = MIN(maxiter,nocc*nocc*nocc*nocc);
    double * timers = calloc(ntimers,sizeof(double));

    double emp4=0.0, emp5=0.0;
    int a=1, i=1, j=1, k=1, klo=1;

    int iter = 0;

    for (int klo=1; klo<=nocc; klo+=kchunk) {
        const int khi = MIN(nocc, klo+kchunk-1);
        int a=1;
        for (int j=1; j<=nocc; j++) {
            for (int i=1; i<=nocc; i++) {
                for (int k=klo; k<=MIN(khi,i); k++) {
                    double t0 = omp_get_wtime();
                    ccsd_trpdrv_omp_cbody_(f1n, f1t, f2n, f2t, f3n, f3t, f4n, f4t, eorb,
                                           &ncor, &nocc, &nvir, &emp4, &emp5, &a, &i, &j, &k, &klo,
                                           Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka, Jij, Jkj, Kij, Kkj,
                                           dintc1, dintx1, t1v1, dintc2, dintx2, t1v2);
                    double t1 = omp_get_wtime();
                    timers[iter] = (t1-t0);

                    iter++;
                    if (iter==maxiter) {
                        printf("Stopping after %d iterations...\n", iter);
                        goto maxed_out;
                    }

                    /* prevent NAN for large maxiter... */
                    if (emp4 >  1000.0) emp4 -= 1000.0;
                    if (emp4 < -1000.0) emp4 += 1000.0;
                    if (emp5 >  1000.0) emp5 -= 1000.0;
                    if (emp5 < -1000.0) emp5 += 1000.0;
                }
            }
        }
    }

maxed_out:
    printf("");

    double tsum =  0.0;
    double tmax = -1.0e10;
    double tmin =  1.0e10;
    for (int i=0; i<iter; i++) {
        //printf("timers[%d] = %f\n", i, timers[i]);
        tsum += timers[i];
        tmax  = MAX(tmax,timers[i]);
        tmin  = MIN(tmin,timers[i]);
    }
    double tavg = tsum / iter;
    printf("TIMING: min=%f, max=%f, avg=%f\n", tmin, tmax, tavg);

    double dgemm_flops = ((8.0*nvir)*nvir)*(nvir+nocc);
    double dgemm_mops  = 8.0*(4.0*nvir*nvir + 2.0*nvir*nocc);

    /* The inner loop of tengy touches 86 f[1234][nt] elements and 8 other arrays...
     * We will just assume flops=mops even though flops>mops */
    double tengy_ops = ((1.0*nvir)*nvir)*(86+8);

    printf("OPS: dgemm_flops=%10.3e dgemm_mops=%10.3e tengy_ops=%10.3e\n",
            dgemm_flops, dgemm_mops, tengy_ops);
    printf("PERF: GF/s=%10.3e GB/s=%10.3e\n",
            1.0e-9*(dgemm_flops+tengy_ops)/tavg, 8.0e-9*(dgemm_mops+tengy_ops)/tavg);

    printf("These are meaningless but should not vary for a particular input:\n");
    printf("emp4=%f emp5=%f\n", emp4, emp5);

    printf("SUCCESS\n");

    return 0;
}
