#ifdef KSR8

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
#define  mxm88_                 MXM88
#define  mxm_                   MXM
#define ortho_                  ORTHO
#define pdspevx_                PDSPEVX
#define pdspgvx_                PDSPGVX
#define pdsptri_                PDSPTRI
#define pstein_                 PSTEIN
#define resid_                  RESID
#define xstop_                  XSTOP
#define dgetavec_               DGETAVEC
#define dlasq1_                 DLASQ1
#define dshellsort2_            DSHELLSORT2
#define dshellsort_             DSHELLSORT
#define maxd_                   MAXD
#define maxi_                   MAXI
#define dgetavec_               DGETAVEC
#define dlas2_                  SLAS2
#define dlascl_                 SLASCL
#define DNRM2 SNRM2
#define DLASCL SLASCL
#define DLAMCH SLAMCH
#define DLAS2  SLAS2
#define DCOPY  SCOPY

#endif






