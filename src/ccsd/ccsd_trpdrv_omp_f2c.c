/* ccsd_trpdrv_omp_f2c.F -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static doublereal c_b8 = -1.;

/* Subroutine */ int ccsd_trpdrv_omp_body__(doublereal *f1n, doublereal *f1t, 
	doublereal *f2n, doublereal *f2t, doublereal *f3n, doublereal *f3t, 
	doublereal *f4n, doublereal *f4t, doublereal *eorb, integer *ncor, 
	integer *nocc, integer *nvir, doublereal *emp4, doublereal *emp5, 
	integer *a, integer *i__, integer *j, integer *k, integer *klo, 
	doublereal *tij, doublereal *tkj, doublereal *tia, doublereal *tka, 
	doublereal *xia, doublereal *xka, doublereal *jia, doublereal *jka, 
	doublereal *kia, doublereal *kka, doublereal *jij, doublereal *jkj, 
	doublereal *kij, doublereal *kkj, doublereal *dintc1, doublereal *
	dintx1, doublereal *t1v1, doublereal *dintc2, doublereal *dintx2, 
	doublereal *t1v2)
{
    /* System generated locals */
    integer f1n_dim1, f1n_offset, f1t_dim1, f1t_offset, f2n_dim1, f2n_offset, 
	    f2t_dim1, f2t_offset, f3n_dim1, f3n_offset, f3t_dim1, f3t_offset, 
	    f4n_dim1, f4n_offset, f4t_dim1, f4t_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6, i__7;

    /* Local variables */
    static integer chunking, b, c__, bb, cc, lnov, lnvv;
    static doublereal emp4i, emp5i, emp4k, emp5k, eaijk;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal denom;

/* chunking is the loop blocking size in the loop nest */
/* formerly associated with the tengy routine. */
/* we have not explored this paramater space but 32 is */
/* optimal for TLB blocking in matrix transpose on most */
/* architectures (especially x86). */
    /* Parameter adjustments */
    --eorb;
    --t1v2;
    --dintx2;
    --dintc2;
    --t1v1;
    --dintx1;
    --dintc1;
    f4t_dim1 = *nvir;
    f4t_offset = 1 + f4t_dim1;
    f4t -= f4t_offset;
    f4n_dim1 = *nvir;
    f4n_offset = 1 + f4n_dim1;
    f4n -= f4n_offset;
    f3t_dim1 = *nvir;
    f3t_offset = 1 + f3t_dim1;
    f3t -= f3t_offset;
    f3n_dim1 = *nvir;
    f3n_offset = 1 + f3n_dim1;
    f3n -= f3n_offset;
    f2t_dim1 = *nvir;
    f2t_offset = 1 + f2t_dim1;
    f2t -= f2t_offset;
    f2n_dim1 = *nvir;
    f2n_offset = 1 + f2n_dim1;
    f2n -= f2n_offset;
    f1t_dim1 = *nvir;
    f1t_offset = 1 + f1t_dim1;
    f1t -= f1t_offset;
    f1n_dim1 = *nvir;
    f1n_offset = 1 + f1n_dim1;
    f1n -= f1n_offset;
    --tij;
    --tkj;
    --tia;
    --tka;
    --xia;
    --xka;
    --jia;
    --jka;
    --kia;
    --kka;
    --jij;
    --jkj;
    --kij;
    --kkj;

    /* Function Body */
    chunking = 32;
    lnov = *nocc * *nvir;
    lnvv = *nvir * *nvir;
    emp4i = 0.;
    emp5i = 0.;
    emp4k = 0.;
    emp5k = 0.;
/* $omp parallel */
/* $omp& shared(eorb) */
/* $omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t) */
/* $omp& shared(t1v1,dintc1,dintx1) */
/* $omp& shared(t1v2,dintc2,dintx2) */
/* $omp& private(eaijk,denom) */
/* $omp& firstprivate(ncor,nocc,nvir,lnov,lnvv,i,j,k,klo) */

/* Performance Note: */

/* By definition, the following does not scale to more than 8 threads */
/* unless nested parallelism (i.e. inside of DGEMM) is used. */
/* It may be prudent to write a manually threaded wrapper for the */
/* cases where single-threaded BLAS is used. */

/* $omp sections */
/* $omp section */
    dgemm_("n", "t", nvir, nvir, nvir, &c_b4, &jia[1], nvir, &tkj[(*k - *klo) 
	    * lnvv + 1], nvir, &c_b5, &f1n[f1n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &tia[1], nvir, &kkj[(*k - *klo) 
	    * lnov + 1], nocc, &c_b4, &f1n[f1n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp section */
    dgemm_("n", "t", nvir, nvir, nvir, &c_b4, &kia[1], nvir, &tkj[(*k - *klo) 
	    * lnvv + 1], nvir, &c_b5, &f2n[f2n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &xia[1], nvir, &kkj[(*k - *klo) 
	    * lnov + 1], nocc, &c_b4, &f2n[f2n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp section */
    dgemm_("n", "n", nvir, nvir, nvir, &c_b4, &jia[1], nvir, &tkj[(*k - *klo) 
	    * lnvv + 1], nvir, &c_b5, &f3n[f3n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &tia[1], nvir, &jkj[(*k - *klo) 
	    * lnov + 1], nocc, &c_b4, &f3n[f3n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp section */
    dgemm_("n", "n", nvir, nvir, nvir, &c_b4, &kia[1], nvir, &tkj[(*k - *klo) 
	    * lnvv + 1], nvir, &c_b5, &f4n[f4n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &xia[1], nvir, &jkj[(*k - *klo) 
	    * lnov + 1], nocc, &c_b4, &f4n[f4n_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp section */
    dgemm_("n", "t", nvir, nvir, nvir, &c_b4, &jka[(*k - *klo) * lnvv + 1], 
	    nvir, &tij[1], nvir, &c_b5, &f1t[f1t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &tka[(*k - *klo) * lnov + 1], 
	    nvir, &kij[1], nocc, &c_b4, &f1t[f1t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp section */
    dgemm_("n", "t", nvir, nvir, nvir, &c_b4, &kka[(*k - *klo) * lnvv + 1], 
	    nvir, &tij[1], nvir, &c_b5, &f2t[f2t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &xka[(*k - *klo) * lnov + 1], 
	    nvir, &kij[1], nocc, &c_b4, &f2t[f2t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp section */
    dgemm_("n", "n", nvir, nvir, nvir, &c_b4, &jka[(*k - *klo) * lnvv + 1], 
	    nvir, &tij[1], nvir, &c_b5, &f3t[f3t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &tka[(*k - *klo) * lnov + 1], 
	    nvir, &jij[1], nocc, &c_b4, &f3t[f3t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp section */
    dgemm_("n", "n", nvir, nvir, nvir, &c_b4, &kka[(*k - *klo) * lnvv + 1], 
	    nvir, &tij[1], nvir, &c_b5, &f4t[f4t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
    dgemm_("n", "n", nvir, nvir, nocc, &c_b8, &xka[(*k - *klo) * lnov + 1], 
	    nvir, &jij[1], nocc, &c_b4, &f4t[f4t_offset], nvir, (ftnlen)1, (
	    ftnlen)1);
/* $omp end sections */
    eaijk = eorb[*a] - (eorb[*ncor + *i__] + eorb[*ncor + *j] + eorb[*ncor + *
	    k]);
/* $omp do collapse(2) */
/* $omp& schedule(static) */
/* $omp& reduction(+:emp5i,emp4i) */
/* $omp& reduction(+:emp5k,emp4k) */
/* WARNING: Do not add IVDEP here.  Code will be incorrect. */
    i__1 = *nvir;
    i__2 = chunking;
    for (bb = 1; i__2 < 0 ? bb >= i__1 : bb <= i__1; bb += i__2) {
	i__3 = *nvir;
	i__4 = chunking;
	for (cc = 1; i__4 < 0 ? cc >= i__3 : cc <= i__3; cc += i__4) {
/* Computing MIN */
	    i__6 = bb + chunking - 1;
	    i__5 = min(i__6,*nvir);
	    for (b = bb; b <= i__5; ++b) {
/* Computing MIN */
		i__7 = cc + chunking - 1;
		i__6 = min(i__7,*nvir);
		for (c__ = cc; c__ <= i__6; ++c__) {
		    denom = -1. / (eorb[*ncor + *nocc + b] + eorb[*ncor + *
			    nocc + c__] + eaijk);
/* fusing emp[45][ki] accumulates may help vectorization... */
		    emp4i = emp4i + denom * (f1t[b + c__ * f1t_dim1] + f1n[
			    c__ + b * f1n_dim1] + f2t[c__ + b * f2t_dim1] + 
			    f3n[b + c__ * f3n_dim1] + f4n[c__ + b * f4n_dim1])
			     * (f1t[b + c__ * f1t_dim1] - f2t[b + c__ * 
			    f2t_dim1] * 2 - f3t[b + c__ * f3t_dim1] * 2 + f4t[
			    b + c__ * f4t_dim1]) - denom * (f1n[b + c__ * 
			    f1n_dim1] + f1t[c__ + b * f1t_dim1] + f2n[c__ + b 
			    * f2n_dim1] + f3n[c__ + b * f3n_dim1]) * (f1t[b + 
			    c__ * f1t_dim1] * 2 - f2t[b + c__ * f2t_dim1] - 
			    f3t[b + c__ * f3t_dim1] + f4t[b + c__ * f4t_dim1] 
			    * 2) + denom * 3 * (f1n[b + c__ * f1n_dim1] * (
			    f1n[b + c__ * f1n_dim1] + f3n[c__ + b * f3n_dim1] 
			    + f4t[c__ + b * f4t_dim1] * 2) + f2n[b + c__ * 
			    f2n_dim1] * f2t[c__ + b * f2t_dim1] + f3n[b + c__ 
			    * f3n_dim1] * f4t[b + c__ * f4t_dim1]);
		    emp4k = emp4k + denom * (f1n[b + c__ * f1n_dim1] + f1t[
			    c__ + b * f1t_dim1] + f2n[c__ + b * f2n_dim1] + 
			    f3t[b + c__ * f3t_dim1] + f4t[c__ + b * f4t_dim1])
			     * (f1n[b + c__ * f1n_dim1] - f2n[b + c__ * 
			    f2n_dim1] * 2 - f3n[b + c__ * f3n_dim1] * 2 + f4n[
			    b + c__ * f4n_dim1]) - denom * (f1t[b + c__ * 
			    f1t_dim1] + f1n[c__ + b * f1n_dim1] + f2t[c__ + b 
			    * f2t_dim1] + f3t[c__ + b * f3t_dim1]) * (f1n[b + 
			    c__ * f1n_dim1] * 2 - f2n[b + c__ * f2n_dim1] - 
			    f3n[b + c__ * f3n_dim1] + f4n[b + c__ * f4n_dim1] 
			    * 2) + denom * 3 * (f1t[b + c__ * f1t_dim1] * (
			    f1t[b + c__ * f1t_dim1] + f3t[c__ + b * f3t_dim1] 
			    + f4n[c__ + b * f4n_dim1] * 2) + f2t[b + c__ * 
			    f2t_dim1] * f2n[c__ + b * f2n_dim1] + f3t[b + c__ 
			    * f3t_dim1] * f4n[b + c__ * f4n_dim1]);
		    emp5i = emp5i + denom * t1v1[b] * dintx1[c__] * (f1t[b + 
			    c__ * f1t_dim1] + f2n[b + c__ * f2n_dim1] + f4n[
			    c__ + b * f4n_dim1] - (f3t[b + c__ * f3t_dim1] + 
			    f4n[b + c__ * f4n_dim1] + f2n[c__ + b * f2n_dim1] 
			    + f1n[b + c__ * f1n_dim1] + f2t[b + c__ * 
			    f2t_dim1] + f3n[c__ + b * f3n_dim1]) * 2 + (f3n[b 
			    + c__ * f3n_dim1] + f4t[b + c__ * f4t_dim1] + f1n[
			    c__ + b * f1n_dim1]) * 4) + denom * t1v1[b] * 
			    dintc1[c__] * (f1n[b + c__ * f1n_dim1] + f4n[b + 
			    c__ * f4n_dim1] + f1t[c__ + b * f1t_dim1] - (f2n[
			    b + c__ * f2n_dim1] + f3n[b + c__ * f3n_dim1] + 
			    f2t[c__ + b * f2t_dim1]) * 2);
		    emp5k = emp5k + denom * t1v2[b] * dintx2[c__] * (f1n[b + 
			    c__ * f1n_dim1] + f2t[b + c__ * f2t_dim1] + f4t[
			    c__ + b * f4t_dim1] - (f3n[b + c__ * f3n_dim1] + 
			    f4t[b + c__ * f4t_dim1] + f2t[c__ + b * f2t_dim1] 
			    + f1t[b + c__ * f1t_dim1] + f2n[b + c__ * 
			    f2n_dim1] + f3t[c__ + b * f3t_dim1]) * 2 + (f3t[b 
			    + c__ * f3t_dim1] + f4n[b + c__ * f4n_dim1] + f1t[
			    c__ + b * f1t_dim1]) * 4) + denom * t1v2[b] * 
			    dintc2[c__] * (f1t[b + c__ * f1t_dim1] + f4t[b + 
			    c__ * f4t_dim1] + f1n[c__ + b * f1n_dim1] - (f2t[
			    b + c__ * f2t_dim1] + f3t[b + c__ * f3t_dim1] + 
			    f2n[c__ + b * f2n_dim1]) * 2);
		}
	    }
	}
    }
/* $omp end do */
/* $omp end parallel */
    *emp4 += emp4i;
    *emp5 += emp5i;
    if (*i__ != *k) {
	*emp4 += emp4k;
	*emp5 += emp5k;
    }
/* (i.ne.k) */
    return 0;
} /* ccsd_trpdrv_omp_body__ */

