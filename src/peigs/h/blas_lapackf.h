*
* $Id: blas_lapackf.h,v 1.9 2002-04-23 18:37:22 edo Exp $
*

#ifdef IEEE

#define DLAMCHE 1.1102230246251565D-16
#define DLAMCHP 1.1102230246251565D-16
#define DLAMCHB 2.0000000000000000D+00
#define DLAMCHS 2.2250738585072013D-308
#define DLAMCHU 2.2250738585072013D-308

#endif

#ifdef KSR8

#define dscal  sscal
#define ddot   sdot
#define daxpy  saxpy
#define dnrm2  snrm2
#define dasum  sasum
#define dcopy  scopy

#define dlagtf slagtf
#define dlagts slagts
#define dlamch slamch
#define dlarnv slarnv

#define heapsort sheapsort
#define neblw2   sneblw2 
#define dstebz3  sstebz3
#define dstebz1  sstebz1
#define dlaebz2  slaebz2
#define damax    samax

#endif

#ifdef CRAY_T3D

#define dscal  SSCAL
#define ddot   SDOT
#define daxpy  SAXPY
#define dnrm2  SNRM2
#define dasum  SASUM
#define dcopy  SCOPY
#define idamax ISAMAX
#define xerbla XERBLA

#define dlagtf SLAGTF
#define dlagts SLAGTS
#define dlamch SLAMCH
#define dlarnv SLARNV
#define xerbl2 XERBL2




#define heapsort SHEAPSORT
#define neblw1   SNEBLW1
#define neblw2   SNEBLW2
#define dstebz3  SSTEBZ3
#define dstebz1  SSTEBZ1
#define dlaebz2  SLAEBZ2
#define pairup   PAIRUP
#define peigscmod PEIGSCMOD

#define sumdc    SUMDC
#define sumd     SUMD
#define sumdv    SUMDV
#define sumi     SUMI
#define sumiv    SUMIV
#define damax    SAMAX


#define mxpara   MXPARA
#define mxmynd   MXMYND
#define mxtick   MXTICK
#define mxread   MXREAD
#define mxwrit   MXWRIT
#define mxsync   MXSYNC
#define mxmynd   MXMYND
#define mxnprc   MXNPRC
#define mxclock  MXCLOCK
#define mxinit   MXINIT
#define mxlbuf   MXLBUF
#define mxpend   MXPEND
#define maxdv    MAXDV
#define menode   MENODE
#define mxbrod   MXBROD
#define mxcombv1 MXCOMBV1
#define mxinit   MXINIT
#define mxend    MXEND
#define mxpara   MXPARA
#define mxtime   MXTIME

#define  choleski     CHOLESKI
#define  inversel     INVERSEL
#define  fmemreq      FMEMREQ
#define  pdspev       PDSPEV
#define  pdspgv       PDSPGV
#define  tresid       TRESID
#define  sonenrm      SONENRM
#define  bortho       BORTHO
#define  mxm35        MXM35
#define  mxm2         MXM2
#define  mxm4         MXM4
#define  mxm5x        MXM5X
#define  mxm88                 MXM88
#define  mxm                   MXM
#define ortho                  ORTHO
#define pdspevx                PDSPEVX
#define pdspgvx                PDSPGVX
#define pdsptri                PDSPTRI
#define pstein                 PSTEIN
#define resid                  RESID
#define xstop                  XSTOP
#define dgetavec               DGETAVEC
#define dlasq1                 DLASQ1
#define dshellsort2            DSHELLSORT2
#define dshellsort             DSHELLSORT
#define maxd                   MAXD
#define maxi                   MAXI
#define dgetavec               DGETAVEC
#define dlas2                  SLAS2
#define dlascl                 SLASCL
#define DNRM2 SNRM2
#define DLASCL SLASCL
#define DLAMCH SLAMCH
#define DLAS2  SLAS2
#define DCOPY  SCOPY

#endif






