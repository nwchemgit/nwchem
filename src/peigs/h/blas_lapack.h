/*
 $Id: blas_lapack.h,v 1.25 2005-07-26 06:06:52 edo Exp $
  c defines */

#define NO_EVEC 0

/*
  CPU definitions and machine precision definitions
  */

#ifndef CRAY_T3D
#ifdef ALPHA

#define DLAMCHE        2.22044604925031308e-16
#define DLAMCHP        2.22044604925031308e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.22507385850720138e-308
#define DLAMCHU        2.22507385850720138e-308
#define SLAMCHE        2.22044604925031308e-16
#define SLAMCHP        2.22044604925031308e-16
#define SLAMCHB        2.e0
#define SLAMCHS        2.22507385850720138e-308
#define SLAMCHU        2.22507385850720138e-308


/*
 Single Precision results
 depsilon      5.9604644775390625e-08
 dbase      2.0000000000000000e+00
 dsafeulp      1.1754943508222875e-38
 */
#endif
#endif

#ifdef HPPA
/*
 Double Precision results
 depsilon   1.1102230246251565e-16
 dbase      2.0000000000000000e+00
 dsafeulp   2.2250738585072013e-308
 dlamch(u)  2.2250738585072013e-308
 */
#define DLAMCHE 1.1102230246251565e-16
#define DLAMCHP 1.1102230246251565e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072013e-308
#define DLAMCHU 2.2250738585072013e-308

#endif

#ifdef IEEE
/*
 Double Precision results
 depsilon   1.1102230246251565e-16
 dbase      2.0000000000000000e+00
 dsafeulp   2.2250738585072013-308
 dlamch(u)  2.2250738585072013-308
 */
#define DLAMCHE 1.1102230246251565e-16
#define DLAMCHP 1.1102230246251565e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072013e-308
#define DLAMCHU 2.2250738585072013e-308

#endif

#ifdef SPARC

/*
  sparc
  */
/*
 Double Precision results
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014e-308

 Single Precision results
   depsilon      5.9604644775390625e-08
      dbase      2.0000000000000000e+00
   dsafeulp      1.1754943508222875e-38
*/

#define DLAMCHE        1.1102230246251565e-16
#define DLAMCHP        1.1102230246251565e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.2250738585072014e-308
#define DLAMCHU        2.2250738585072014e-308

#endif

#ifdef SPARC64
/*./teslamch
 
  Double Precision results
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014-308
 dlamch(u)       2.2250738585072014-308
 
  Single Precision results
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014-308
 slamch(u)       2.2250738585072014-308
*/
#define DLAMCHE        2.22044604925031308e-16
#define DLAMCHP        2.22044604925031308e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.22507385850720138e-308
#define DLAMCHU        2.22507385850720138e-308
/* this values were OK for WS 5.0, break for WS 6.0
#define DLAMCHE        1.1102230246251565e-16
#define DLAMCHP        1.1102230246251565e-16
#define DLAMCHB        2.e0
#define DLAMCHS        2.2250738585072014e-308
#define DLAMCHU        2.2250738585072014e-308
*/

#endif
#ifdef __crayx1
#define DLAMCHE 2.2204460492503130e-16
#define DLAMCHP 2.2204460492503130e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHU 2.2250738585072014e-308

#define dscal_  sscal_
#define ddot_   sdot_
#define daxpy_  saxpy_
#define dnrm2_  snrm2_
#define dasum_  sasum_
#define dcopy_  scopy_
#define idamax_ isamax_

/*
  lapack
  */

#define dlagtf_ slagtf_
#define dlagts_ slagts_
#define dlamch_ slamch_
#define dlarnv_ slarnv_


#endif
#ifdef PENTIUM
/* wild ass guess; same as sparc */
#define DLAMCHE 2.2204460492503131e-16
#define DLAMCHP 2.2204460492503131e-16
#define DLAMCHB 2.0000000000000000e+00
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHU 2.2250738585072014e-308

#endif
#ifdef MIPS
/*
  SGI PowerChallenge R4010 fpu ... R8000 chip
  seems to have the same reps .... probably have
  to print out hex to get real differences
  */

#define DLAMCHE 1.1102230246251565e-16
#define DLAMCHP 1.1102230246251565e-16
#define DLAMCHB 2.e0
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHU 2.2250738585072014e-308


#endif
#ifdef RS6000

/* rs6000 */

#define DLAMCHE 0.111022302462515654e-15
#define DLAMCHP 0.111022302462515654e-15
#define DLAMCHB 2.e0
#define DLAMCHS 0.22250738585072013e-307
#define DLAMCHU 0.22250738585072013e-307


/*
  depsilon  0.111022302462515654e-15 
  dbase   2.00000000000000000     
  dsafeulp  0.22250738585072013e-307 
  depsilon  0.5960464478e-07 
  dbase   2.000000000     
  dsafeulp  0.1175494351e-37 
  */
#endif
#ifdef RS600064

/* rs6000 64 -bit*/

#define DLAMCHE 1.1102230246251565e-16
#define DLAMCHP 1.1102230246251565e-16
#define DLAMCHB 2.e0
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHU 2.2250738585072014e-308


/*
   depsilon      1.1102230246251565e-16
      dbase      2.0000000000000000e+00
   dsafeulp      2.2250738585072014e-308
 dlamch(u)       2.2250738585072014e-308
  */
#endif
#ifdef i860

/*
  with -Knoieee
  */ 

#define DLAMCHE 1.1102230246251565e-016
#define DLAMCHB 2.e0
#define DLAMCHS 2.2250738585072014e-308
#define DLAMCHE 1.1102230246251565e-016
#define DLAMCHP 1.1102230246251565e-016
#define DLAMCHB  2.e0
#define DLAMCHS  2.2250738585072014e-308
#define DLAMCHU  2.2250738585072014e-308


#endif
#ifdef FUJITSU_VPP

#define DLAMCHE  1.1102230246251565e-016
#define DLAMCHP  1.1102230246251565e-016
#define DLAMCHB  2.e0
#define DLAMCHS  2.2250738585072014e-308
#define DLAMCHU  2.2250738585072014e-308

#endif

#ifdef KSR

#define DLAMCHE  1.1102230246251565e-16
#define DLAMCHB  2.
#define DLAMCHS  2.2250738585072014e-308
#define DLAMCHE  1.1102230246251565e-16
#define DLAMCHP  1.1102230246251565e-16
#define DLAMCHB  2.
#define DLAMCHS  2.2250738585072014e-308
#define DLAMCHU  2.2250738585072014e-308



#endif
#ifdef KSR8

#define DLAMCHE  1.1102230246251565e-16
#define DLAMCHB  2.e0
#define DLAMCHS  2.2250738585072014e-308
#define DLAMCHE  1.1102230246251565e-16
#define DLAMCHP  1.1102230246251565e-16
#define DLAMCHB  2.e0
#define DLAMCHS  2.2250738585072014e-308
#define DLAMCHU  2.2250738585072014e-308


#define dscal_  sscal_
#define ddot_   sdot_
#define daxpy_  saxpy_
#define dnrm2_  snrm2_
#define dasum_  sasum_
#define dcopy_  scopy_

/*
  lapack
  */

#define dlagtf_ slagtf_
#define dlagts_ slagts_
#define dlamch_ slamch_
#define dlarnv_ slarnv_

/*
  peigs
  */

#define heapsort_ sheapsort_
#define neblw2_   sneblw2_ 
#define dstebz3_  sstebz3_
#define dstebz1_  sstebz1_
#define dlaebz2_  slaebz2_
#define damax_    samax_
#endif

#ifdef CRAY_T3D
#undef DLAMCHE  2.22044604925031308e-16
#undef DLAMCHP  2.22044604925031308e-16
#undef DLAMCHB  2.e0
#undef DLAMCHS  2.22507385850720138e-308
#undef DLAMCHU  2.22507385850720138e-308
#undef SLAMCHE  2.22044604925031308e-16
#undef SLAMCHP  2.22044604925031308e-16
#undef SLAMCHB  2.e0
#undef SLAMCHS  2.22507385850720138e-308
#undef SLAMCHU  2.22507385850720138e-308

#define DLAMCHE  2.22044604925031308e-16
#define DLAMCHP  2.22044604925031308e-16
#define DLAMCHB  2.e0
#define DLAMCHS  2.22507385850720138e-308
#define DLAMCHU  2.22507385850720138e-308
#define SLAMCHE  2.22044604925031308e-16
#define SLAMCHP  2.22044604925031308e-16
#define SLAMCHB  2.e0
#define SLAMCHS  2.22507385850720138e-308
#define SLAMCHU  2.22507385850720138e-308
 

#define dscal_  SSCAL
#define ddot_   SDOT
#define daxpy_  SAXPY
#define dnrm2_  SNRM2
#define dasum_  SASUM
#define dcopy_  SCOPY
#define idamax_ ISAMAX
#define xerbla_ XERBLA

/*
  lapack
  */

#define dlagtf_ SLAGTF
#define dlagts_ SLAGTS
#define dlamch_ SLAMCH
#define dlarnv_ SLARNV
#define xerbl2_ XERBL2



/*
peigs
*/

#define heapsort_ SHEAPSORT
#define neblw1_   SNEBLW1
#define neblw2_   SNEBLW2
#define dstebz3_  SSTEBZ3
#define dstebz1_  SSTEBZ1
#define dlaebz2_  SLAEBZ2
#define pairup_   PAIRUP
#define peigs_cmod_ PEIGS_CMOD

#define sumdc_    SUMDC
#define sumd_     SUMD
#define sumdv_    SUMDV
#define sumi_     SUMI
#define sumiv_    SUMIV
#define damax_    SAMAX

/*
  mx
*/

#define mxpara_   MXPARA
#define mxmynd_   MXMYND
#define mxtick_   MXTICK
#define mxread_   MXREAD
#define mxwrit_   MXWRIT
#define mxsync_   MXSYNC
#define mxmynd_   MXMYND
#define mxnprc_   MXNPRC
#define mxclock_  MXCLOCK
#define mxinit_   MXINIT
#define mxlbuf_   MXLBUF
#define mxpend_   MXPEND
#define maxdv_    MAXDV
#define menode_   MENODE
#define mxbrod_   MXBROD
#define mxcombv1_ MXCOMBV1
#define mxinit_   MXINIT
#define mxend_    MXEND
#define mxpara_   MXPARA
#define mxtime_   MXTIME

/*
  peigs ctof 
*/

#define  choleski_     CHOLESKI
#define  inversel_     INVERSEL
#define  fmemreq_      FMEMREQ
#define  pdspev_       PDSPEV
#define  pdspgv_       PDSPGV
#define  tresid_       TRESID
#define  sonenrm_      SONENRM
#define  bortho_       BORTHO
#define  mxm35_        MXM35
#define  mxm2_         MXM2
#define  mxm4_         MXM4
#define  mxm5x_        MXM5X
#define  mxm88_                    MXM88
#define  mxm_                      MXM
#define ortho_                     ORTHO
#define pdspevx_                   PDSPEVX
#define pdspgvx_                   PDSPGVX
#define pdsptri_                   PDSPTRI
#define pstein_                    PSTEIN
#define resid_                     RESID
#define xstop_                     XSTOP
#define dgetavec_                  DGETAVEC
#define dlasq1_                    DLASQ1
#define dshellsort2_               DSHELLSORT2
#define dshellsort_             DSHELLSORT
#define maxd_                   MAXD
#define maxi_                   MAXI
#define dgetavec_               DGETAVEC
#define dlas2_                  SLAS2
#define dlascl_                 SLASCL
#define dsterf_                 SSTERF
#define ga_sync_                 GA_SYNC

#endif
