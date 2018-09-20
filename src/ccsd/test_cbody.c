#include <stdio.h>
#include <stdlib.h>

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
        a[i] = 1.0/(10.0+i);
    }
    return a;
}

int main(int argc, char* argv[])
{
    int ncor, nocc, nvir;

    if (argc!=3) {
        printf("Usage: ./test_cbody nocc nvir\n");
        return argc;
    } else {
        ncor = 0;
        nocc = atoi(argv[1]);
        nvir = atoi(argv[2]);
    }

    if (nocc<1 || nvir<1) {
        printf("Arguments must be non-negative!\n");
        return 1;
    }

    printf("Test driver for cbody with nocc=%d, nvir=%d\n", nocc, nvir);

    const int nbf = ncor + nocc + nvir;
    const int lnvv = nvir * nvir;
    const int lnoo = nvir * nocc;

    const int nkpass = 1; /* assume unlimited memory */
    const int kchunk = (nocc - 1)/nkpass + 1;

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

    double emp4=0.0, double emp5=0.0;
    int a=1, i=1, j=1, k=1, klo=1;

    return 0;
}
