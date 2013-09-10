      SUBROUTINE DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR,EXPS2)
*
* Density matrices between L and R
*
* I12 = 1 => only one-body density
* I12 = 2 => one- and two-body-density matrices
*
* Jeppe Olsen,      Oct 94
* GAS modifications Aug 95
* Two body density added, '96
*
* Table-Block driven, June 97
*
* Two-body density is stored as 
*     rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
*                ijkl = ij*(ij-1)/2+kl, ij.ge.kl
*
* If the twobody density matrix is calculated, then also the
* expectation value of the spin is evaluated.
* The latter is realized as
* S**2
*      = S+S- + Sz(Sz-1)
*      = -Sum(ij) a+i alpha a+j beta a i beta a j alpha + Nalpha +
*        1/2(N alpha - N beta))(1/2(N alpha - Nbeta) - 1)
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KSTSTS, KSTSTD, KCONSPA, KCONSPB, KCB, KSB
      INTEGER*8 KINSCR, KSIOIO, KCIOIO, KC2, KSSCR, KCSCR
      INTEGER*8 KI1, KI2, KI3, KI4, KXI1S, KXI2S, KXI3S, KXI4S
      INTEGER*8 KSBLTP, KCBLTP, KSVST
      INTEGER*8 KRHO1S, KRHO1P, KXNATO, KRHO1SM, KXNATSM, KOCCSM
      INTEGER*8 KLOCSTR, KLREO, KLZ, KLZSCR
      INTEGER*8 KLLBTL, KLLEBTL, KLI1BTL, KLIBTL, KLSCLFCL
      INTEGER*8 KLLBTR, KLLEBTR, KLI1BTR, KLIBTR, KLSCLFCR
      ! for addressing of WORK
      REAL*8 INPRDD
*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "spinfo.inc"
#include "glbbas.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
#include "lucinp.inc"
#include "clunit.inc"
*. Scratch for string information
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. Specific input
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*)
      IF( I12 .eq. 0 ) RETURN
*     before I forget it :
      CALL QENTER('DENSI')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
      ZERO = 0.0D0
      CALL DZERO(RHO1,NACOB ** 2 )
      IF(I12.EQ.2)
     &CALL DZERO(RHO2,NACOB ** 2 *(NACOB**2+1)/2)
*
C?     WRITE(6,*) ' ISSPC ICSPC in DENSI2 ',ISSPC,ICSPC
*
* Info for this internal space
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
      JATP = 1
      JBTP = 2
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offsets for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      ILSM = ISSM
      IRSM = ICSM
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. connection matrices for supergroups
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(WORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF (IPRDEN.GE.10) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
*. replace with MXINKA !!!
CTF
* MAXI and MAXK should not be zero, which might happen in very small cases.
* Set to one then.
      MAXI = max(1,MIN(MXINKA,MXSTBL))
      MAXK = max(1,MIN(MXINKA,MXSTBL))
*Largest active orbital block belonging to given type and symmetry
CTF
* Using local MXTSOB_L (MXTSOB is now a list parameter!)
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
      MAXIJ = MXTSOB_L ** 2
*.Local scratch arrays for blocks of C and sigma
      IF (IPRDEN.GE.10) write(6,*) ' DENSI2 : MXSB MXTSOB_L MXSOOB ',
     &       MXSB,MXTSOB_L,MXSOOB
      IF (IPRDEN.GE.10) WRITE(6,*) ' ICISTR,LBLOCK ',ICISTR,LBLOCK
      IF(ICISTR.EQ.1) THEN
        CALL MEMMAN(KCB,LBLOCK,'ADDL  ',2,'KCB   ')
        CALL MEMMAN(KSB,LBLOCK,'ADDL  ',2,'KSB   ')
      END IF
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF (IPRDEN.GE.10)
     &   WRITE(6,*) ' Density scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KSIOIO),IPRDIA)
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KCIOIO),IPRDIA)
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESCPH(WORK(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSTFSMSPGP,MXPNSMST,
     &     NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &     NELFSPGP,
     &     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &     IPHGAS,NHLFSPGP,MNHL,IADVICE)
      IF (IPRDEN.GE.10) THEN
        WRITE(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF (IPRDEN.GE.10)
     &    WRITE(6,*) ' Space for resolution matrices ',LSCR2
      LSCR12 = MAX(LBLOCK,2*LSCR2)
*. It is assumed that the third block already has been allocated, so
      KC2 = KVEC3
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
      CALL MEMMAN(KI1,  LSCR3       ,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KI2,  LSCR3       ,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KI3,  LSCR3       ,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KI4,  LSCR3       ,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
      CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
      CALL MEMMAN(KXI3S,LSCR3       ,'ADDL  ',2,'XI3S  ')
      CALL MEMMAN(KXI4S,LSCR3       ,'ADDL  ',2,'XI4S  ')
*. Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.0 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
* scratch space containing active one body
      CALL MEMMAN(KRHO1S,NACOB ** 2,'ADDL  ',2,'RHO1S ')
*. For natural orbitals
      CALL MEMMAN(KRHO1P,NACOB*(NACOB+1)/2,'ADDL  ',2,'RHO1P ')
      CALL MEMMAN(KXNATO,NACOB **2,'ADDL  ',2,'XNATO ')
*. Natural orbitals in symmetry blocks
      CALL MEMMAN(KRHO1SM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KXNATSM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KOCCSM,NACOB ,'ADDL  ',2,'RHO1S ')
*
*. Space for one block of string occupations and two arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
      DO K12 = 1, 1
        CALL MEMMAN(KLOCSTR(K12),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      END DO
      DO I1234 = 1, 2
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
*. Arrays for partitioning of Left vector = sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTL ,NTTS  ,'ADDL  ',1,'LBT_L  ')
      CALL MEMMAN(KLLEBTL,NTTS  ,'ADDL  ',1,'LEBT_L ')
      CALL MEMMAN(KLI1BTL,NTTS  ,'ADDL  ',1,'I1BT_L ')
      CALL MEMMAN(KLIBTL ,8*NTTS,'ADDL  ',1,'IBT_L  ')
      CALL MEMMAN(KLSCLFCL,NTTS, 'ADDL  ',2,'SCLF_L')
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),
     &     WORK(KLI1BTL),WORK(KLIBTL),0,ITTSS_ORD)
*. Number of BLOCKS
        NBLOCKL = IFRMR(WORK(KLI1BTL),1,NBATCHL)
     &         + IFRMR(WORK(KLLBTL),1,NBATCHL) - 1
*. Arrays for partitioning of Right  vector = C
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTR ,NTTS  ,'ADDL  ',1,'LBT_R  ')
      CALL MEMMAN(KLLEBTR,NTTS  ,'ADDL  ',1,'LEBT_R ')
      CALL MEMMAN(KLI1BTR,NTTS  ,'ADDL  ',1,'I1BT_R ')
      CALL MEMMAN(KLIBTR ,8*NTTS,'ADDL  ',1,'IBT_R  ')
      CALL MEMMAN(KLSCLFCR,NTTS, 'ADDL  ',2,'SCLF_R')
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KCIOIO),ISMOST(1,ICSM),
     &     NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),
     &     WORK(KLI1BTR),WORK(KLIBTR),0,ITTSS_ORD)
*. Number of BLOCKS
        NBLOCKR = IFRMR(WORK(KLI1BTR),1,NBATCHR)
     &         + IFRMR(WORK(KLLBTR),1,NBATCHR) - 1
C?      WRITE(6,*) ' DENSI2T :NBLOCKR =',NBLOCKR

      IF(ICISTR.EQ.1) THEN
         WRITE(6,*) ' Sorry, ICISTR = 1 is out of fashion'
         WRITE(6,*) ' Switch to ICISTR = 2 - or reprogram '
         Call Abend2( ' DENSI2T : ICISTR = 1 in use ' )
      ELSE IF(ICISTR.GE.2) THEN
        S2_TERM1 = 0.0D0
        CALL GASDN2(I12,RHO1,RHO2,L,R,L,R,WORK(KC2),
     &       WORK(KCIOIO),WORK(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),WORK(KCBLTP),WORK(KSBLTP),NACOB,
     &       WORK(KNSTSO(IATP)),WORK(KISTSO(IATP)),
     &       WORK(KNSTSO(IBTP)),WORK(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,
     &       MAXK,MAXI,LBLOCK,LBLOCK,WORK(KCSCR),WORK(KSSCR),
     &       SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &       WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &       WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),WORK(KINSCR),
     &       MXPOBS,IPRDEN,WORK(KRHO1S),LUL,LUR,
     &       PSSIGN,PSSIGN,WORK(KRHO1P),WORK(KXNATO),
     &       NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),WORK(KLI1BTL),
     &       WORK(KLIBTL),
     &       NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),WORK(KLI1BTR),
     &       WORK(KLIBTR),WORK(KCONSPA),WORK(KCONSPB),
     &       WORK(KLSCLFCL),WORK(KLSCLFCR),S2_TERM1,IUSE_PH,IPHGAS,
     &       MXTSOB)
C     KLLBTR  KLLEBTR KLI1BTR KLIBTR
      END IF
C?    WRITE(6,*) ' Memcheck in densi2 after GASDN2'
C?    CALL MEMCHK
*
*
*. Add terms from hole-hole commutator
      IF(IUSE_PH.EQ.1) THEN
*. Overlap between left and right vector
       XLR = INPRDD(L,R,LUR,LUL,1,-1)
       CALL RHO1_HH(RHO1,XLR)
      END IF
* Natural Orbitals
      CALL LNATORB(RHO1,NSMOB,NTOOBS,NACOBS,NINOBS,
     &             IREOST,WORK(KXNATO),
     &             WORK(KRHO1SM),WORK(KOCCSM),NACOB,
     &             WORK(KRHO1P),IPRDEN)
*
      IF (IPRDEN.GE.5) THEN
        WRITE(6,*) ' One-electron density matrix '
        WRITE(6,*) ' ============================'
        CALL WRTMT_LU(RHO1,NTOOB,NTOOB,NTOOB,NTOOB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-electron density '
          CALL PRSYM(RHO2,NACOB**2)
        END IF
      END IF
*
      IF (I12.EQ.2) THEN
* <L!S**2|R>
        EXPS2 = S2_TERM1 + NAEL +
     &          0.5*(NAEL-NBEL)*(0.5*(NAEL-NBEL)-1)
        IF(IPRDEN.GT.0) THEN
          WRITE(6,*) ' Term 1 to S2 ', S2_TERM1
          WRITE(6,*) ' Expectation value of S2 ', EXPS2
        END IF
      ELSE
        EXPS2 = 0.0D0
      END IF

*. Eliminate local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'DENSI ')
      CALL QEXIT('DENSI')
C     WRITE(6,*) ' Leaving DENSI '
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDN2(I12,RHO1,RHO2,
     &           L,R,CB,SB,C2,ICOCOC,ISOCOC,ICSMOS,ISSMOS,
     &           ICBLTP,ISBLTP,NACOB,NSSOA,ISSOA,NSSOB,ISSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,
     &           MXPNGAS,NOBPTS,IOBPTS,MAXK,MAXI,LC,LS,
     &           CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &           MXPOBS,IPRNT,RHO1S,LUL,LUR,PSL,PSR,RHO1P,XNATO ,
     &           NBATCHL,LBATL,LEBATL,I1BATL,IBLOCKL,
     &           NBATCHR,LBATR,LEBATR,I1BATR,IBLOCKR,
     &           ICONSPA,ICONSPB,SCLFAC_L,SCLFAC_R,S2_TERM1,
     &           IUSE_PH,IPHGAS,MXTSOB)
*
*
* Jeppe Olsen , Winter of 1991
* GAS modificatios, August 1995
*
* Table driven, June 97
*
* Last revision : Jan. 98 (IUSE_PH,IPHGAS added)
*
* =====
* Input
* =====
*
* I12    : = 1 => calculate one-electrondensity matrix
*          = 2 => calculate one-and two- electrondensity matrix
* RHO1   : Initial one-electron density matrix
* RHO2   : Initial two-electron density matrix
*
* ICOCOC : Allowed type combinations for C
* ISOCOC : Allowed type combinations for S(igma)
* ICSMOS : Symmetry array for C
* ISSMOS : Symmetry array for S
* ICBLTP : Block types for C
* ISBLTP : Block types for S
*
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* ISSOA : Offset for strings if given type and symmetry, alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* ISSOB : Offset for strings if given type and symmetry, beta strings
* NBEL  : Number of active beta electrons
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
*
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* RHO1S: Scratch array for one body
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
* The L and R vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB),ISOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST),ISSMOS(NSMST)
      INTEGER ICBLTP(*),ISBLTP(*)
      INTEGER NSSOA(NSMST,NOCTPA),ISSOA(NSMST,NOCTPA)
      INTEGER NSSOB(NSMST,NOCTPB),ISSOB(NSMST,NOCTPB)
      INTEGER SXSTSM(NSMSX,NSMST)
      INTEGER STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER ADSXA(MXPOBS,2*MXPOBS),ASXAD(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,NSMOB),IOBPTS(MXPNGAS,NSMOB)
      INTEGER NELFSPGP(MXPNGAS,*)
*. Info on batches and blocks
      INTEGER  LBATL(NBATCHL),LEBATL(NBATCHL),I1BATL(NBATCHL),
     &         IBLOCKL(8,*)
      INTEGER  LBATR(NBATCHR),LEBATR(NBATCHR),I1BATR(NBATCHR),
     &         IBLOCKR(8,*)
*. Interaction between supergroups
      INTEGER ICONSPA(NOCTPA,NOCTPA),ICONSPB(NOCTPB,NOCTPB)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),XI1S(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION RHO1S(*)
      DIMENSION SCLFAC_L(*),SCLFAC_R(*)
*.
      INTEGER LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      INTEGER RASM(4),RBSM(4),RATP(4),RBTP(4),RSGN(5),RTRP(5)
      REAL * 8 INPROD,L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*)
      DIMENSION RHO1P(*),XNATO(*)
*
      CALL QENTER('GASDN')
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GASDN2 speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' NACOB,MAXK,NGAS,IDC,MXPOBS',
     &             NACOB,MAXK,NGAS,IDC,MXPOBS
        WRITE(6,*) ' LUL, LUR ', LUL,LUR
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial L vector '
        IF(LUL.EQ.0) THEN
          CALL WRTRS2(L,ISSMOS,ISBLTP,ISOCOC,NOCTPA,NOCTPB,
     &                NSSOA,NSSOB,NSMST)
        ELSE
          CALL WRTVCD(L,LUL,1,-1)
        END IF
        WRITE(6,*) ' Initial R vector '
        IF(LUR.EQ.0) THEN
          CALL WRTRS2(R,ICSMOS,ICBLTP,ICOCOC,NOCTPA,NOCTPB,
     &                NSSOA,NSSOB,NSMST)
        ELSE
          CALL WRTVCD(R,LUR,1,-1)
        END IF
      END IF
* Loop over batches over L blocks
      IF(LUL.NE.0) CALL REWINO(LUL)
      DO 10001 IBATCHL = 1, NBATCHL
*. Obtain L blocks
        NBLKL = LBATL(IBATCHL)
        IF(NTEST.GE.200)
     &    WRITE(6,*) ' Left batch, number of blocks',IBATCHL,NBLKL
        DO IIL  = 1,NBLKL
          IL  = I1BATL(IBATCHL)-1+IIL
          IATP = IBLOCKL(1,IL)
          IBTP = IBLOCKL(2,IL)
          IASM = IBLOCKL(3,IL)
          IBSM = IBLOCKL(4,IL)
          IOFF = IBLOCKL(5,IL)
          IF(NTEST.GE.200)
     &    WRITE(6,*) 'IATP IBTP IASM IBSM',IATP,IBTP,IASM,IBSM
          ISCALE = 0
          IF(NTEST.GE.200)
     &    WRITE(6,*) 'IOFF ',IOFF
          CALL GSTTBL(L,SB(IOFF),IATP,IASM,IBTP,IBSM,ISOCOC,
     &                NOCTPA,NOCTPB,NSSOA,NSSOB,PSL,ISOOSC,IDC,
     &                PSL,LUL,C2,NSMST,ISCALE,SCLFAC_L(IL))
        END DO
*. Loop over batches  of R vector
        IF(LUR.NE.0) CALL REWINO(LUR)
        DO 9001 IBATCHR = 1, NBATCHR
*. Read R blocks into core
        NBLKR = LBATR(IBATCHR)
        IF(NTEST.GE.200)
     &    WRITE(6,*) ' Right batch, number of blocks',IBATCHR,NBLKR
        DO IIR  = 1,NBLKR
          IR  = I1BATR(IBATCHR)-1+IIR
          JATP = IBLOCKR(1,IR)
          JBTP = IBLOCKR(2,IR)
          JASM = IBLOCKR(3,IR)
          JBSM = IBLOCKR(4,IR)
          JOFF = IBLOCKR(5,IR)
          IF(NTEST.GE.200)
     &    WRITE(6,*) ' JATP JBTP JASM JBSM ',JATP,JBTP,JASM,JBSM
*. Read R blocks into core
*
*. Only blocks interacting with current batch of L are read in
*. Loop over L  blocks in batch
          DO IIL = 1, NBLKL
            IL  = I1BATL(IBATCHL)-1+IIL
            IATP = IBLOCKL(1,IL)
            IBTP = IBLOCKL(2,IL)
            IASM = IBLOCKL(3,IL)
            IBSM = IBLOCKL(4,IL)
*. Well, permutations of L blocks
            CALL PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PS,PL,
     &              LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
            DO IPERM = 1, NPERM
              IIASM = LASM(IPERM)
              IIBSM = LBSM(IPERM)
              IIATP = LATP(IPERM)
              IIBTP = LBTP(IPERM)

              IAEXC = ICONSPA(IIATP,JATP)
              IBEXC = ICONSPB(IIBTP,JBTP)
              IF(IAEXC.EQ.0.AND.IIASM.NE.JASM) IAEXC = 1
              IF(IBEXC.EQ.0.AND.IIBSM.NE.JBSM) IBEXC = 1
              IABEXC = IAEXC + IBEXC
              IF(IABEXC.LE.I12) THEN
                INTERACT = 1
              END IF
            END DO
          END DO
*.          ^ End of checking whether C-block is needed
          ISCALE = 0
          IF(INTERACT.EQ.1) THEN
            ISCALE = 0
            CALL GSTTBL(R,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                  NOCTPA,NOCTPB,NSSOA,NSSOB,PSR,ICOOSC,IDC,
     &                  PCL,LUR,C2,NSMST,ISCALE,SCLFAC_R(IR))
          ELSE
C             WRITE(6,*) ' TTSS for C block skipped  '
C             CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
            CALL IFRMDS(LBL,-1,1,LUR)
            CALL SKPRCD2(LBL,-1,LUR)
            SCLFAC_R(IR) = 0.0D0
          END IF
*
*
          IF(NTEST.GE.100) THEN
            IF(INTERACT.EQ.1) THEN
              WRITE(6,*) ' TTSS for C block read in  '
              CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
            ELSE
              WRITE(6,*) ' TTSS for C block skipped  '
              CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
            END IF
          END IF
        END DO
*. Loop over L and R blocks in batches and obtain  contribution from
* given L and R blocks
          DO 10000 IIL = 1, NBLKL
            IL  = I1BATL(IBATCHL)-1+IIL
          IF(SCLFAC_L(IL).NE.0.0D0) THEN
            IATP = IBLOCKL(1,IL)
            IBTP = IBLOCKL(2,IL)
            IASM = IBLOCKL(3,IL)
            IBSM = IBLOCKL(4,IL)
            IOFF = IBLOCKL(5,IL)
*
            NIA = NSSOA(IASM,IATP)
            NIB = NSSOB(IBSM,IBTP)
*. Possible permutations of L blocks
            CALL PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PLR,
     &           LATP,LBTP,LASM,LBSM,LSGN,LTRP,NLPERM)
            DO 9999 ILPERM = 1, NLPERM
C             write(6,*) ' Loop 9999 ILPERM = ', ILPERM
              IIASM = LASM(ILPERM)
              IIBSM = LBSM(ILPERM)
              IIATP = LATP(ILPERM)
              IIBTP = LBTP(ILPERM)
              NIIA = NSSOA(IIASM,IIATP)
              NIIB = NSSOB(IIBSM,IIBTP)
*
              IF(LTRP(ILPERM).EQ.1) THEN
                LROW = NSSOA(LASM(ILPERM-1),LATP(ILPERM-1))
                LCOL = NSSOB(LBSM(ILPERM-1),LBTP(ILPERM-1))
                CALL TRPMT3(SB(IOFF),LROW,LCOL,C2)
                CALL COPVEC(C2,SB(IOFF),LROW*LCOL)
               END IF
              IF(LSGN(ILPERM).EQ.-1)
     &        CALL SCALVE(SB(IOFF),-1.0D0,NIA*NIB)

              DO 9000 IIR = 1, NBLKR
                IR  = I1BATR(IBATCHR)-1+IIR
              IF(SCLFAC_R(IR).NE.0.0D0) THEN
                JATP = IBLOCKR(1,IR)
                JBTP = IBLOCKR(2,IR)
                JASM = IBLOCKR(3,IR)
                JBSM = IBLOCKR(4,IR)
                JOFF = IBLOCKR(5,IR)
*
                NJA = NSSOA(JASM,JATP)
                NJB = NSSOB(JBSM,JBTP)
*
                IAEXC = ICONSPA(JATP,IIATP)
                IBEXC = ICONSPB(JBTP,IIBTP)
*
                IF(IAEXC.EQ.0.AND.JASM.NE.IIASM) IAEXC = 1
                IF(IBEXC.EQ.0.AND.JBSM.NE.IIBSM) IBEXC = 1
                IABEXC = IAEXC + IBEXC
*
                IF(IABEXC.LE.I12) THEN
                  INTERACT = 1
                ELSE
                  INTERACT = 0
                END IF
*
                IF(INTERACT.EQ.1) THEN
*. Possible permutations of this block
                   CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,
     &                  PSR,PLR,RATP,RBTP,RASM,RBSM,RSGN,RTRP,
     &                  NRPERM)
*. Well, spin permutations are simple to handle
* if there are two terms just calculate and and multiply with
* 1+PSL*PSR
                     IF(NRPERM.EQ.1) THEN
                       FACTOR = 1.0D0
                     ELSE
                       FACTOR = 1.0D0 +PSL*PSR
                     END IF
                     SCLFAC = FACTOR*SCLFAC_L(IL)*SCLFAC_R(IR)
                     IF(INTERACT.EQ.1.AND.SCLFAC.NE.0.0D0) THEN
                     IF(NTEST.GE.20) THEN
                       WRITE(6,*) ' RSDNBB will be called for '
                       WRITE(6,*) ' L block : '
                       WRITE(6,'(A,5I5)')
     &                 ' IIASM IIBSM IIATP IIBTP',
     &                   IIASM,IIBSM,IIATP,IIBTP
                       WRITE(6,*) ' R  block : '
                       WRITE(6,'(A,5I5)')
     &                 ' JASM JBSM JATP JBTP',
     &                   JASM,JBSM,JATP,JBTP
                       WRITE(6,*) ' IOFF,JOFF ', IOFF,JOFF
                       WRITE(6,*) ' SCLFAC = ', SCLFAC
                     END IF
                     CALL GSDNBB2(I12,RHO1,RHO2,
     &                    IIASM,IIATP,IIBSM,IIBTP,
     &                    JASM,JATP,JBSM,JBTP,NGAS,
     &                    NELFSPGP(1,IOCTPA-1+IIATP),
     &                    NELFSPGP(1,IOCTPB-1+IIBTP),
     &                    NELFSPGP(1,IOCTPA-1+JATP),
     &                    NELFSPGP(1,IOCTPB-1+JBTP),
     &                    NAEL,NBEL,IAGRP,IBGRP,
     &                    SB(IOFF),CB(JOFF),C2,
     &                    ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                    MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,
     &                    SSCR,CSCR,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                    X,NSMOB,NSMST,NSMSX,NSMDX,
     &                    NIIA,NIIB,NJA,NJB,MXPOBS,
     &                    IPRNT,NACOB,RHO1S,SCLFAC,
     &                    S2_TERM1,IUSE_PH,IPHGAS,
     &                    MXTSOB)
                          IF(NTEST.GE.500) THEN
                            write(6,*) ' Updated rho1 '
                            call WRTMT_LU(rho1,nacob,nacob,nacob,nacob)
                          END IF
*
                     END IF
                   END IF
                END IF
 9000         CONTINUE
*. End of loop over R blocks in Batch
 9999     CONTINUE
*. Transpose or scale L block to restore order ??
          IF(LTRP(NLPERM+1).EQ.1) THEN
            CALL TRPMT3(SB(IOFF),NIB,NIA,C2)
            CALL COPVEC(C2,SB(IOFF),NIA*NIB)
          END IF
          IF(LSGN(NLPERM+1).EQ.-1)
     &    CALL SCALVE(SB(IOFF),-1.0D0,NIA*NIB)
*
          END IF
10000     CONTINUE
*. End of loop over L blocks in batch
 9001   CONTINUE
*.      ^ End of loop over batches of R blocks
10001 CONTINUE
*.    ^ End of loop over batches of L blocks
      CALL QEXIT('GASDN')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETD1(RHO1B,ISM,IGAS,JSM,JGAS)
*
* Extract TS block of one-electron density matrix
*
      IMPLICIT REAL*8(A,H,O-Z)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "orbinp.inc"
*. output
      DIMENSION RHO1B(*)
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
*
      II = IOBPTS(IGAS,ISM)
      IJ = IOBPTS(JGAS,JSM)

*
      DO I = 1, NI
        DO J = 1, NJ
          IABS = I-1+II
          JABS = J-1+IJ
          IJABS = (JABS-1)*NTOOB + IABS
          RHO1B((J-1)*NI+I) = WORK(KRHO1-1+IJABS)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Block of one-electron density matrix'
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
        WRITE(6,*) 'IGAS,ISM,JGAS,JSM',IGAS,ISM,JGAS,JSM
        CALL WRTMT_LU(RHO1B,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETD2(RHO2B,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,
     &                 ICOUL)
*. Extract given TS block from the 2e-density matrix
*.
*. Jeppe Olsen, Some day in Hfors CITY, winter 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Initial implementation, KISS to the MAX !!!
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "orbinp.inc"
*. Output
      DIMENSION RHO2B(*)
*. output is in the form i,j,k,l
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
*
      IELMNT = 0
      DO L = 1, NL
        DO K = 1, NK
          DO J = 1, NJ
            DO I = 1, NI
              IELMNT = IELMNT + 1
              RHO2B(IELMNT) = GETD2E(I,IGAS,ISM,J,JGAS,JSM,
     &                               K,KGAS,KSM,L,LGAS,LSM)
            END DO
          END DO
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' BLock of two-electron density'
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Type and symmetry of the 4 orbitals (i j k l )'
        WRITE(6,'(1H ,8(1X,I4) )')
     &  IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
*
        NIJ = NI*NJ
        NKL = NK*NL
        CALL WRTMT_LU(RHO2B,NIJ,NKL,NIJ,NKL)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION GETD2E(I,IGAS,ISM,J,JGAS,JSM,K,KGAS,KSM,L,LGAS,LSM)
*
* Obtain element of two-electron density matrix
* Currently stored without symmetry
*
*. 2-electron density is assumed stored in wotk(krho2)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "orbinp.inc"
*
      IABS = I + IOBPTS(IGAS,ISM)-1
      JABS = J + IOBPTS(JGAS,JSM)-1
      KABS = K + IOBPTS(KGAS,KSM)-1
      LABS = L + IOBPTS(LGAS,LSM)-1
*
      IJ = (JABS-1)*NTOOB+IABS
      KL = (LABS-1)*NTOOB+KABS
      IF(IJ.GE.KL) THEN
        IJKL = IJ*(IJ-1)/2+KL
      ELSE
        IJKL = KL*(KL-1)/2+IJ
      END IF
*
      X = WORK(KRHO2-1+IJKL)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Element of two electron density matrix'
        WRITE(6,*) 'I,IGAS,ISM,J,JGAS,JSM,K,KGAS,KSM,L,LGAS,LSM'
        WRITE(6,'( 1H ,(12I4) )')
     &  I,IGAS,ISM,J,JGAS,JSM,K,KGAS,KSM,L,LGAS,LSM
        WRITE(6,*) ' IJ, KL ', IJ,KL
        WRITE(6,*) ' IABS JABS KABS LABS',IABS,JABS,KABS,LABS
        WRITE(6,*) 'Address and value', IJKL, X
      END IF
*
      GETD2E = X
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GSBBD1(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &                  IUSE_PH,IPHGAS)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen
* Particle-Hole version of Jan. 98
*
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1 : Updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), ITSOB(*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(*),XI1S(*)
      DIMENSION I2(*),XI2S(*)
*.Local arrays ( assume MPNGAS = 16 ) !!!
      DIMENSION ITP(16*16),JTP(16*16)
*
      DIMENSION IJ_REO(2),IJ_DIM(2),IJ_SM(2),IJ_TP(2),IJ_AC(2)
      DIMENSION IJ_OFF(2)
      DIMENSION ISCR(2)
      DIMENSION ICGRP(16),ISGRP(16)
*
*.Local arrays
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GSBBD1 in action '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' GSBBD1, sclfac ',SCLFAC
      END IF
*
      IFRST = 1
      JFRST = 1
*. Number of partitionings over column strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Groups defining supergroups
C          GET_SPGP_INF(ISPGP,ITP,IGRP)
      CALL GET_SPGP_INF(ICCTP,IGRP,ICGRP)
      CALL GET_SPGP_INF(ISCTP,IGRP,ISGRP)

* Type of single excitations that connects the two column strings
      CALL SXTYP2_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL,IPHGAS)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(NTEST.GE.1000)
     &WRITE(6,*) ' ISCSM,ICCSM IJSM ', ISCSM,ICCSM,IJSM
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
*. Hvilken vej skal vi valge,
        NOP = 2
        IJ_AC(1) = 2
        IJ_AC(2) = 1
        IJ_TP(1) = ITYP
        IJ_TP(2) = JTYP
        IF(IUSE_PH.EQ.1) THEN
          CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TP,IJ_AC,IJ_REO,SIGNIJ)
        ELSE
          IJ_REO(1) = 1
          IJ_REO(2) = 2
          SIGNIJ = 1.0D0
        END IF
*
        ISCR(1) = IJ_AC(1)
        ISCR(2) = IJ_AC(2)
        IJ_AC(1) = ISCR(IJ_REO(1))
        IJ_AC(2) = ISCR(IJ_REO(2))
*
        ISCR(1) = ITYP
        ISCR(2) = JTYP
        IJ_TP(1) = ISCR(IJ_REO(1))
        IJ_TP(2) = ISCR(IJ_REO(2))

        DO 800 ISM = 1, NSMOB
*. new i and j so new intermediate strings
          KFRST = 1
*
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.1000) write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
          IBIORB = IOBPTS(ITYP,ISM)
          IBJORB = IOBPTS(JTYP,JSM)
*. Reorder
*
          ISCR(1) = ISM
          ISCR(2) = JSM
          IJ_SM(1) = ISCR(IJ_REO(1))
          IJ_SM(2) = ISCR(IJ_REO(2))
*
          ISCR(1) = NIORB
          ISCR(2) = NJORB
          IJ_DIM(1) = ISCR(IJ_REO(1))
          IJ_DIM(2) = ISCR(IJ_REO(2))
*
          ISCR(1) = IBIORB
          ISCR(2) = IBJORB
          IJ_OFF(1) = ISCR(IJ_REO(1))
          IJ_OFF(2) = ISCR(IJ_REO(2))
*

          IF(NTEST.GE.2000)
     &    WRITE(6,*) ' NIORB NJORB ', NIORB,NJORB
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          SCLFACS = SCLFAC*SIGNIJ
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
          CALL ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC(1))
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          ONE = 1.0D0
          CALL ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,
     &         I2,XI2S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC(1))
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I2,XI2S,IJ_DIM(1),
     &             NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF
C         WRITE(6,*) ' NKAEFF NKASTR', NKAEFF,NKASTR

*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = MIN(KTOP + MAXK,NKAEFF)
              IF(KTOP.EQ.NKAEFF) THEN
                KEND = 1
              ELSE
                KEND = 0
              END IF
              LKABTC = KTOP - KBOT +1

*. This is the place to start over partitioning of I strings
              DO 701 IPART = 1, NIPART
                IBOT = (IPART-1)*MAXI+1
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP - IBOT + 1
* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
                DO JJORB = 1,IJ_DIM(2)
                  ICGOFF = 1 + (JJORB-1)*LKABTC*NIBTC
                  CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                 LKABTC,I1(KBOT+(JJORB-1)*NKASTR),
     &                 XI1S(KBOT+(JJORB-1)*NKASTR) )
                END DO
* Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
                DO IIORB = 1,IJ_DIM(1)
*.Gather S Block
                  ISGOFF = 1 + (IIORB-1)*LKABTC*NIBTC
                  CALL MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,
     &                   LKABTC,I2(KBOT+(IIORB-1)*NKASTR),
     &                   XI2S(KBOT+(IIORB-1)*NKASTR) )
                END DO
*
                IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' CSCR and SSCR '
                 CALL WRTMT_LU(CSCR,IJ_DIM(2),NKI,IJ_DIM(2),NKI)
                 CALL WRTMT_LU(SSCR,IJ_DIM(1),NKI,IJ_DIM(1),NKI)
                END IF
*
*. And then the hard  work
                NKI = LKABTC*NIBTC
                FACTORC = 0.0D0
                FACTORAB = 1.0D0
                CALL MATML7(RHO1S,SSCR,CSCR,IJ_DIM(1),IJ_DIM(2),NKI,
     &               IJ_DIM(1),NKI,IJ_DIM(2),FACTORC,FACTORAB,1)
*
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' Block to one-body density '
                  CALL WRTMT_LU(RHO1S,IJ_DIM(1),IJ_DIM(2),
     &                              IJ_DIM(1),IJ_DIM(2))
                END IF
*. Scatter out to complete matrix
                DO JJORB = 1, IJ_DIM(2)
                  JORB = IJ_OFF(2)-1+JJORB
                  DO IIORB = 1, IJ_DIM(1)
                    IORB = IJ_OFF(1)-1+IIORB
                    RHO1((JORB-1)*NACOB+IORB) =
     &              RHO1((JORB-1)*NACOB+IORB) +
     &              RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)
                  END DO
                END DO
*               /\ End of hard work

  701     CONTINUE
*. /\ end of this I partitioning
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
*. End of loop over I partitioninigs
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
C!    Call Abend2( ' enforrced stop in RSBBD1 ' )
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,SB,CB,
     &                  ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC)
*
* Contributions to two-electron density matrix from column excitations
*
* GAS version, '96 , Jeppe Olsen
*
* =====
* Input
* =====
* RHO2  : two body density matrix to be updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO2 : Updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S, I2,XI2S : For holding creations/annihilations
*              type and symmetry
*
* Jeppe Olsen, Fall of 96
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST), SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO2(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
C     INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSBBD2A in action '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
      IFRST = 1
      JFRST = 1
*
* Type of single excitations that connects the two column strings
      CALL DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,NGAS,ISEL,ICEL)
*.Symmetry of Double excitation that connects IBSM and JBSM
*. For general use : STSTSX => STSTDX
      IDXSM = STSTSX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0) GOTO 2001
      IF(NTEST.GE.1000)
     &WRITE(6,*) ' ISCSM,ICCSM ', ISCSM,ICCSM
      DO 2000 IDXTP =  1, NDXTP
        ITYP = ITP(IDXTP)
        JTYP = JTP(IDXTP)
        KTYP = KTP(IDXTP)
        LTYP = LTP(IDXTP)
        IF(NTEST.GE.1000)
     &  write(6,*) ' ITYP JTYP KTYP LTYP ', ITYP,JTYP,KTYP,LTYP
        DO 1950 IKOBSM = 1, NSMOB
          JLOBSM = SXDXSX(IKOBSM,IDXSM)
          IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
          KFRST = 1
*. Loop over of symmetry of i orbitals
          DO 1940 ISM = 1, NSMOB
          KSM = ADSXA(ISM,IKOBSM)
          NI = NOBPTS(ITYP,ISM)
          NK = NOBPTS(KTYP,KSM)
          IF(NI.EQ.0.OR.NK.EQ.0) GOTO 1940
*. Loop over batches of j orbitals
          DO 1930 JSM = 1, NSMOB
          LSM = ADSXA(JSM,JLOBSM)
          NJ = NOBPTS(JTYP,JSM)
          NL = NOBPTS(LTYP,LSM)
          IF(NJ.EQ.0.OR.NL.EQ.0) GOTO 1930
*
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          KOFF = IOBPTS(KTYP,KSM)
          LOFF = IOBPTS(LTYP,LSM)
*
          IF(IOFF.LT.KOFF) GOTO 1930
          IF(JOFF.LT.LOFF) GOTO 1930
*
*
* =========================================================================
*                    Use N-2 projection method
* =========================================================================
*
              IFIRST = 1
*. Loop over batches of I strings
              NPART = NROW/MAXI
              IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
              IF(NTEST.GE.2000)
     &        write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
              DO 1801 IPART = 1, NPART
                IBOT = 1+(IPART-1)*MAXI
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                KBOT = 1- MAXK
                KTOP = 0
 1800           CONTINUE
                  KBOT = KBOT + MAXK
                  KTOP = KTOP + MAXK
*
* =========================================================
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
* =========================================================
*
                  IONE = 1
                  JLBOFF = 1
                  IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                    NJL = NJ*(NJ+1)/2
                    JLSM = 1
                  ELSE
                    NJL = NJ * NL
                    JLSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
                  CALL QENTER('ADADS')
                  I12 = 1
                  K12 = 1
                  IONE = 1
                  CALL ADADST_GAS(IONE,JSM,JTYP,NJ,
     &                            IONE,LSM,LTYP,NL,
     &                        ICCTP,ICCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        JFRST,KFRST,I12,K12,SCLFAC)
                  JFRST = 0
                  KFRST = 0
*
                  CALL QEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
                  CALL QENTER('MATCG')
                  DO  IJL = 1, NJL
                    CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                    I1JL = (L-1)*NJ+J
*. JAN28
                    IF(JLSM.NE.0) THEN
                      IJLE = J*(J-1)/2+L
                    ELSE
                      IJLE = IJL
                    END IF
*. JAN28
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
C                   JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                    JLOFF = (JLBOFF-1+IJLE-1)*NKBTC*NIBTC+1
                    IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,I1JL),XI1S(1,I1JL))
                    END IF
                  END DO
                  CALL QEXIT ('MATCG')
*
*
* =========================================================
*
*. obtain sb(KB,IA,ik) = sum(IB)<KB!a kb a ib !IB>S(IA,IB)
*
* =========================================================
*
                  IONE = 1
                  IKBOFF = 1
                  IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                    NIK = NI*(NI+1)/2
                    IKSM = 1
                  ELSE
                    NIK = NI * NK
                    IKSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
CT                CALL QENTER('ADADS')
                  I12 = 2
                  K12 = 1
                  IONE = 1
                  IF(IFRST.EQ.1) KFRST = 1
                  ONE = 1.0D0
                  CALL ADADST_GAS(IONE,ISM,ITYP,NI,
     &                            IONE,KSM,KTYP,NK,
     &                        ISCTP,ISCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        IFRST,KFRST,I12,K12,ONE   )
                  IFRST = 0
                  KFRST = 0
*
CT                CALL QEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  I = 0
                  K = 1
*
CT                CALL QENTER('MATCG')
                  DO  IIK = 1, NIK
                    CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                    I1IK = (K-1)*NI+I
*. JAN28
                    IF(IKSM.NE.0) THEN
                      IIKE = I*(I-1)/2+K
                    ELSE
                      IIKE = IIK
                    END IF
*. JAN28
*.SB(IA,KB,ik) = +/-S(IA,a+ka+iIA)
C                   IKOFF = (IKBOFF-1+IIK-1)*NKBTC*NIBTC+1
                    IKOFF = (IKBOFF-1+IIKE-1)*NKBTC*NIBTC+1
                    IF(IKSM.EQ.1.AND.I.EQ.K) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(SSCR(IKOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,I1IK),XI1S(1,I1IK))
                    END IF
                  END DO
CT                CALL QEXIT ('MATCG')
*
*
* =================================================================
*
* RHO2C(ik,jl)  = RHO2C(ik,jl) - sum(Ia,Kb)SB(Ia,Kb,ik)*CB(Ia,Kb,jl)
*
* =================================================================
*
* The minus ??
*
* Well, the density matrices are constructed as

* <I!a+i a+k aj al!> = -sum(K) <I!a+ia+k!K><J!aj al!K>, and
* the latter matrices are the ones we are constructing
*
              IOFF = IOBPTS(ITYP,ISM)
              JOFF = IOBPTS(JTYP,JSM)
              KOFF = IOBPTS(KTYP,KSM)
              LOFF = IOBPTS(LTYP,LSM)
              NTESTO = NTEST
C?            IF(IOFF.EQ.3.AND.JOFF.EQ.3.AND.KOFF.EQ.4.AND.LOFF.EQ.4)
C?   &            NTEST = 5000
                  LDUMMY = NKBTC*NIBTC
                  IF(NTEST.GE.2000) THEN
                    WRITE(6,*) ' CSCR matrix '
                    CALL WRTMT_LU(CSCR,LDUMMY,NJL,LDUMMY,NJL)
                    WRITE(6,*) ' SSCR matrix '
                    CALL WRTMT_LU(SSCR,LDUMMY,NIK,LDUMMY,NIK)
                  END IF

                  IF(IFIRST.EQ.1) THEN
                    FACTOR = 0.0D0
                  ELSE
                    FACTOR = 1.0D0
                  END IF
C                 MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
                  LDUMMY = NKBTC*NIBTC
                  ONEM = -1.0D0
                  CALL MATML7(X,SSCR,CSCR,NIK,NJL,
     &                        LDUMMY,NIK,LDUMMY,NJL,
     &                        FACTOR,ONEM,1)
                  IFIRST = 0
                  IF(NTEST.GE.2000) THEN
                    WRITE(6,*) ' Updated X matrix'
                    CALL WRTMT_LU(X,NIK,NJL,NIK,NJL)
                  END IF

*
                IF(KEND.EQ.0) GOTO 1800
*. End of loop over partitionings of resolution strings
 1801         CONTINUE
*. Rho2(ik,jl) has been constructed for ik,jl belonging to
*. Scatter out to density matrix
              IOFF = IOBPTS(ITYP,ISM)
              JOFF = IOBPTS(JTYP,JSM)
              KOFF = IOBPTS(KTYP,KSM)
              LOFF = IOBPTS(LTYP,LSM)
              CALL ADTOR2(RHO2,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,
     &                    NACOB)
C                  ADTOR2(RHO2,RHO2T,ITYPE,
C    &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)

 1930       CONTINUE
 1940     CONTINUE
 1950   CONTINUE
 2000 CONTINUE
 2001 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXK,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  CJRES,SIRES,NORB,NTESTG,SCLFAC,S2_TERM1,
     &                  MXTSOB)
*
* alpha-beta contribution to two-particle density matrix
* from given c-block and s-block.
*
* S2_TERM1 = - <L!a+i alpha a+jbeta a i beta a j alpha !R>
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* X : Space for block of two-electron integrals
*
* Jeppe Olsen, Fall of 1996
*
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
*.Input
      DIMENSION CB(*),SB(*)
*. Output
      DIMENSION RHO2(*)
*.Scratch
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION CJRES(*),SIRES(*)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
*
      CALL QENTER('GSD2B')
      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' GSBBD2B speaking '
        WRITE(6,*) ' ================== '
      END IF
C?    WRITE(6,*) ' NJAS NJB = ',NJA,NJB
C?    WRITE(6,*) ' IAGRP IBGRP = ', IAGRP,IBGRP
C?    WRITE(6,*) ' MXPNGAS = ', MXPNGAS
C?    WRITE(6,*) ' NSMOB = ', NSMOB
      IROUTE = 3
*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC)
      CALL SXTYP_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all Ka strings
*. a+j!ka> = +/-/0 * !Ja>
          CALL ADSTN_GAS(JSM,JTYP,JATP,JASM,IAGRP,
     &                   I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,
     &                   SCLFAC)
*. a+i!ka> = +/-/0 * !Ia>
          ONE    = 1.0D0
          CALL ADSTN_GAS(ISM,ITYP,IATP,IASM,IAGRP,
     &                   I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,
     &                   ONE   )
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
C             COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
              CALL COMPRS2LST(I1,XI1S,NJ,I3,XI3S,NI,NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF

*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
          DO 1801 IKABTC = 1, NKABTC
C?          write(6,*) ' Batch over kstrings ', IKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, NJ
              CALL GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*. Obtain S(ka,i,Ib) for Ka in batch
            DO II = 1, NI
              CALL GET_CKAJJB(SB,NI,NIA,SIRES,LKABTC,NIB,
     &             II,I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*
              DO 1930 KSM = 1, NSMOB
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
C?              WRITE(6,*) ' Loop 1930, KSM LSM ',KSM,LSM
                KOFF = IOBPTS(KTYP,KSM)
                LOFF = IOBPTS(LTYP,LSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP) IKORD=1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
                ONE = 1.0D0
                CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE)
                IF(NKBSTR.EQ.0) GOTO 1930
*
*. Update two-electron density matrix
*  Rho2b(ij,kl) =  Sum(ka)S(Ka,i,Ib)<Ib!Eb(kl)!Jb>C(Ka,j,Jb)
*
                ZERO = 0.0D0
                CALL SETVEC(X,ZERO,NI*NJ*NK*NL)
*
C               WRITE(6,*) ' Before call to ABTOR2'
                CALL ABTOR2(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,X,NI,NJ,NK,NL,NKBSTR,
     &               I4,XI4S,I2,XI2S,IKORD,MXTSOB)
*. contributions to Rho2(ij,kl) has been obtained, scatter out
C?              WRITE(6,*) ' Before call to ADTOR2'
C?              WRITE(6,*) ' RHO2B (X) matrix '
C?              call WRTMT_LU(x,ni*nj,nk*nl,ni*nj,nk*nl)
*. Contribution to S2
                IF(KTYP.EQ.JTYP.AND.KSM.EQ.JSM.AND.
     &            ITYP.EQ.LTYP.AND.ISM.EQ.LSM) THEN
                  DO I = 1, NI
                    DO J = 1, NJ
                      IJ = (J-1)*NI+I
                      JI = (I-1)*NJ+J
                      NIJ = NI*NJ
                      S2_TERM1 = S2_TERM1-X((JI-1)*NIJ+IJ)
                    END DO
                  END DO
                END IF

     &
                CALL ADTOR2(RHO2,X,2,
     &                NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
C?              write(6,*) ' updated density matrix '
C?              call prsym(rho2,NORB*NORB)

 1930         CONTINUE
 2000       CONTINUE
 1801     CONTINUE
*. End of loop over partitioning of alpha strings
 1940   CONTINUE
 2001 CONTINUE
*
 9999 CONTINUE
*
*
      CALL QEXIT('GSD2B')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GSDNBB2(I12,RHO1,RHO2,
     &                  IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                  NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  NAEL,NBEL,
     &                  IJAGRP,IJBGRP,
     &                  SB,CB,C2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  X,NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,RHO1S,
     &                  SCLFAC,S2_TERM1,IUSE_PH,IPHGAS,MXTSOB)
*
* Contributions to density matrix from sigma block (iasm iatp, ibsm ibtp ) and
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS : Number of As'es
* IAOC : Occpation of each AS for alpha strings in L
* IBOC : Occpation of each AS for beta  strings in L
* JAOC : Occpation of each AS for alpha strings in R
* JBOC : Occpation of each AS for beta  strings in R
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* MXPNGAS : Largest number of As'es allowed by program
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* Rho1, RHo2 : Updated density blocks
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*. Input
      DIMENSION CB(*),SB(*)
*. Output
      DIMENSION RHO1(*),RHO2(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION  I1(*),XI1S(*),I2(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION C2(*)
*
      CALL QENTER('GSDNB')
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      NTESTO= NTEST
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSDNBB2 :  R block '
        WRITE(6,*) ' ==================='
        CALL WRTMT_LU(CB,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' GSDNBB2 :  L block '
        WRITE(6,*) ' ==================='
        CALL WRTMT_LU(SB,NIA,NIB,NIA,NIB)
*
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in L '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in L '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in R '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in R '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB
*
        WRITE(6,*) 'SCLFAC =',SCLFAC
      END IF
      IACTIVE = 0
*
      IF(IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
*  beta contribution to RHO1
* =============================
*
C?      WRITE(6,*) ' GSBBD1 will be called (beta)'
        CALL GSBBD1(RHO1,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &       NGAS,IBOC,JBOC,
     &       SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS)
C?      WRITE(6,*) ' GSBBD1 was called '
C?      WRITE(6,*) ' Memory check '
C?      CALL MEMCHK
*
* ================================
* beta-beta contribution to RHO2
* ================================
*
        IF(I12.EQ.2.AND.NBEL.GE.2) THEN
C?        WRITE(6,*) ' GSBBD2A will be called (beta)'
          CALL GSBBD2A(RHO2,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &         NGAS,IBOC,JBOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC)
C?        WRITE(6,*) ' GSBBD2A was called '
*
C              GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
C    &         NGAS,ISEL,ICEL,SB,CB,
C    &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
C    &         NOBPTS,IOBPTS,MAXI,MAXK,
C    &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
C    &         NSMOB,NSMST,NSMSX,MXPOBS)
        END IF
      END IF
*
      IF(IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* =============================
*  alpha contribution to RHO1
* =============================
*
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
C?        WRITE(6,*) ' GSBBD1 will be called (alpha)'
        CALL GSBBD1(RHO1,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &       NGAS,IAOC,JAOC,SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS)
C?        WRITE(6,*) ' GSBBD1 was called '
        IF(I12.EQ.2.AND.NAEL.GE.2) THEN
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
C?        WRITE(6,*) ' GSBBD2A will be called (alpha)'
          CALL GSBBD2A(RHO2,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &         NGAS,IAOC,JAOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC)
C?        WRITE(6,*) ' GSBBD2A was called '
        END IF
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMAT(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
* ===================================
*  alpha-beta contribution to RHO2
* ===================================
*
      IF(I12.EQ.2.AND.NAEL.GE.1.AND.NBEL.GE.1) THEN
*. Routine uses transposed blocks
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
C?      WRITE(6,*) ' GSBBD2B will be called '
        IUSEAB = 0
        CALL GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                    JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                    IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                    SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                    NOBPTS,IOBPTS,MAXK,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                    NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                    SSCR,CSCR,NACOB,NTEST,SCLFAC,S2_TERM1,
     &                    MXTSOB)
C?      WRITE(6,*) ' GSBBD2B was called '
     &
C     GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
C    &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
C    &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
C    &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
C    &                  NOBPTS,IOBPTS,MAXK,
C    &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
C    &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C    &                  CJRES,SIRES,NORB,NTEST,
C    &                  MXTSOB)
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMAT(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
      CALL QEXIT('GSDNB')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PERTDN
     &(N,LU0,LUN,ISM,ISPC,VEC1,VEC2,RHO1N,RHO2N,LUSC1,LUSC2)
*
* Construct one body density matrix of order N
*
*      Jeppe + Dage, Nov. 11 1995
*                    Debugged Jan 31 '97
*
* Note : I12 added, April 98
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 KLDEN1, KLDEN2
      ! for addressing of WORK
*
*. Should not be called with ICISTR = 1
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
C     INCLUDE 'clunit.inc'
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Output
      DIMENSION RHO1N(*),RHO2N(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'PERTDN')
*
      LRHO1 = NTOOB**2
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
      CALL MEMMAN(KLDEN1,LRHO1,'ADDL  ',2,'KLDEN1')
      IF(I12.EQ.2) THEN
        CALL MEMMAN(KLDEN2,LRHO2,'ADDL  ',2,'KLDEN2')
      END IF
*
      LBLK = -1
      ZERO = 0.0D0
      CALL SETVEC(RHO1N,ZERO,LRHO1)
      IF (I12.EQ.2) THEN
        CALL SETVEC(RHO2N,ZERO,LRHO2)
      END IF
*
      DO L = 0, N
C?      write(6,*) ' Will load next pair of vectors '
        NMINL = N - L
CTOBE   IF(L.LE.NMINL) THEN
*. put correction vector L and NMINL on LUSC1 and LUSC2, respectively
          IF(L.EQ.0) THEN
             CALL COPVCD(LU0,LUSC1,VEC1,1,LBLK)
          ELSE
             CALL SKPVCD(LUN,L-1,VEC1,1,LBLK)
             CALL REWINO(LUSC1)
             CALL COPVCD(LUN,LUSC1,VEC1,0,LBLK)
          END IF
*
          IF(NMINL.EQ.0) THEN
             CALL COPVCD(LU0,LUSC2,VEC1,1,LBLK)
          ELSE
             CALL SKPVCD(LUN,NMINL-1,VEC1,1,LBLK)
             CALL REWINO(LUSC2)
             CALL COPVCD(LUN,LUSC2,VEC1,0,LBLK)
          END IF
C?      write(6,*) ' next pair of vectors loaded '
* Do the densi
          LEQR = 0
C?        WRITE(6,*) ' Calling DENSI2 '
          CALL DENSI2(I12,WORK(KLDEN1),WORK(KLDEN2),VEC1,VEC2,
     &                LUSC1,LUSC2,EXPS2)
          WRITE(6,*) ' Home from DENSI2 '
C         WRITE(6,*) ' Densities fresh from DENSI2 '
C         CALL WRTMT_LU(WORK(KLDEN1),NTOOB,NTOOB,NTOOB,NTOOB)
C         IF(I12.EQ.2) THEN
C         CALL PRSYM(WORK(KLDEN2),NTOOB**2)
C         END IF
*
CTOBE     IF(L.NE.NMINL) THEN
*. The matrix <L! E !NMINL> was calculated, add <NMINL! E ! L>
*. as simple transposition
CTOBE        CALL TRPAD(WORK(KLDEN),ONE,NTOOB)
CTOBE     END IF
          ONE = 1.0D0
          CALL VECSUM(RHO1N,RHO1N,WORK(KLDEN1),ONE,ONE,LRHO1)
          IF(I12.EQ.2) THEN
            CALL VECSUM(RHO2N,RHO2N,WORK(KLDEN2),ONE,ONE,LRHO2)
          END IF
CTOBE   END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Density matrix of order in perturbation ', N
        WRITE(6,*) ' ==========================================='
        WRITE(6,*)
        WRITE(6,*) ' One-body density '
        WRITE(6,*) ' ================ '
        CALL WRTMT_LU(RHO1N,NTOOB,NTOOB,NTOOB,NTOOB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-body density '
          WRITE(6,*) ' ================ '
          CALL PRSYM(RHO2N,NTOOB**2)
        END IF
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ', IDUM,'PERTDN')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE REORHO1(RHO1I,RHO1O,IRHO1SM)
*
* The density matric rho1 is given in complete form
* Extract symmetry blocks with symmetry IRHO1SM
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "multd2h.inc"
*. Input
      DIMENSION RHO1I(NTOOB,NTOOB)
*. Output
      DIMENSION RHO1O(*)
*
      IMOFF = 1
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IRHO1SM)
        IOFF = IBSO(ISM)
        JOFF = IBSO(JSM)
*
        NI  = NOCOBS(ISM)
        NJ =  NOCOBS(JSM)
        DO I = 1, NI
          DO J = 1, NJ
            IP = IREOST(IOFF-1+I)
            JP = IREOST(JOFF-1+J)
            RHO1O(IMOFF-1+(J-1)*NI+I) = RHO1I(IP,JP)
          END DO
        END DO
        IMOFF = IMOFF + NI*NJ
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' REORHO1 in action '
        WRITE(6,*) ' Symmetry of blocks extracted ',IRHO1SM
        WRITE(6,*) ' Input density '
        CALL WRTMT_LU(RHO1I,NTOOB,NTOOB,NTOOB,NTOOB)
        WRITE(6,*)
        WRITE(6,*) ' extracted blocks : '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
        CALL PRHONE(RHO1O,NOCOBS,IRHO1SM,NSMOB,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE RHO1_HH(RHO1,XLR)
*
* Add terms from hole-hole commutator to one-particle density matrix :
*
*  2*<L!R> to diagonal for Hole orbitals
*
* Jeppe Olsen, Jan. 1998 (<= Just to show we are geared for the millenium
*                            change)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
*Input and output
      DIMENSION RHO1(NACOB,NACOB)
*
      DO IORB = 1, NACOB
       IF(IPHGAS(ITPFTO(IORB)).EQ.2)
     & RHO1(IORB,IORB) = RHO1(IORB,IORB) + 2.0D0*XLR
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' One-body density matrix with hole-commutator term '
        WRITE(6,*) ' (See manual)'
        CALL WRTMT_LU(RHO1,NACOB,NACOB,NACOB,NACOB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRADEN(I12,RHO1,RHO2,NL,NR,LUL,LUR)
*
* Transition density matrices between the NL states stored on LUL
* and the NR states stored on LUR
* Density matrices between L and R stored on LU
*
* I12 = 1 => only one-bodydensity
* I12 = 2 => one- and two-body-density matrices
*
* Jeppe Olsen, July 97 (from densi2)
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KSTSTS, KSTSTD, KCONSPA, KCONSPB, KCB, KSB
      INTEGER*8 KINSCR, KSIOIO, KCIOIO, KLL, KLR, KC2, KSSCR, KCSCR
      INTEGER*8 KI1, KI2, KI3, KI4, KXI1S, KXI2S, KXI3S, KXI4S
      INTEGER*8 KSBLTP, KCBLTP, KSVST
      INTEGER*8 KRHO1S, KRHO1P, KXNATO, KRHO1SM, KXNATSM, KOCCSM
      INTEGER*8 KLOCSTR, KLREO, KLZ, KLZSCR
      INTEGER*8 KLLBTL, KLLEBTL, KLI1BTL, KLIBTL, KLSCLFCL
      INTEGER*8 KLLBTR, KLLEBTR, KLI1BTR, KLIBTR, KLSCLFCR
      ! for addressing of WORK
*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "spinfo.inc"
#include "glbbas.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*./LUCINP/
*. Used : NSMOB
#include "lucinp.inc"
#include "clunit.inc"
*. Scratch for string information
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*.Output
      DIMENSION RHO1(*),RHO2(*)
*. Before I forget it :
      CALL QENTER('TRADN')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'TRADEN')
      ZERO = 0.0D0
      CALL SETVEC(RHO1,ZERO ,NL*NR*NACOB ** 2 )
      IF(I12.EQ.2)
     &CALL SETVEC(RHO2,ZERO ,NL*NR*NACOB ** 2 *(NACOB**2+1)/2)
*
* Info for this internal space
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
      JATP = 1
      JBTP = 2
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offsets for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      ILSM = ISSM
      IRSM = ICSM
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. connection matrices for supergroups
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(WORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRDEN.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
*. replace with MXINKA !!!
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
      MAXIJ = MXTSOB_L ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(IPRDEN.GE.2) write(6,*) ' DENSI2 : MXSB MXTSOB_L MXSOOB ',
     &       MXSB,MXTSOB_L,MXSOOB
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' ICISTR,LBLOCK ',ICISTR,LBLOCK
      IF(ICISTR.EQ.1) THEN
        CALL MEMMAN(KCB,LBLOCK,'ADDL  ',2,'KCB   ')
        CALL MEMMAN(KSB,LBLOCK,'ADDL  ',2,'KSB   ')
      END IF
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      write(6,*) ' TRADEN : ISSPC,ICSPC : ',ISSPC,ICSPC
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KSIOIO),IPRDIA)
      CALL IAIBCM_GAS(LCMBSPC(ICSPC),ICMBSPC(1,ICSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KCIOIO),IPRDIA)
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESC(WORK(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
      IF(IPRDEN.GE.2) THEN
        WRITE(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*. Allocation of blocks of vector
      CALL MEMMAN(KLL,LBLOCK,'ADDL  ',2,'KLL   ')
      CALL MEMMAN(KLR,LBLOCK,'ADDL  ',2,'KLR   ')
      LSCR12 = MAX(LBLOCK,2*LSCR2)
      KC2 = KVEC3
      CALL MEMMAN(KC2,LSCR12,'ADDL  ',2,'KC2   ')
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
      CALL MEMMAN(KI1,  LSCR3       ,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KI2,  LSCR3       ,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KI3,  LSCR3       ,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KI4,  LSCR3       ,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
      CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
      CALL MEMMAN(KXI3S,LSCR3       ,'ADDL  ',2,'XI3S  ')
      CALL MEMMAN(KXI4S,LSCR3       ,'ADDL  ',2,'XI4S  ')
*. Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
      write(6,*) ' ISSM ICSM IDC ',ISSM,ICSM,IDC
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.0 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
* scratch space containing active one body
      CALL MEMMAN(KRHO1S,NACOB ** 2,'ADDL  ',2,'RHO1S ')
*. For natural orbitals
      CALL MEMMAN(KRHO1P,NACOB*(NACOB+1)/2,'ADDL  ',2,'RHO1P ')
      CALL MEMMAN(KXNATO,NACOB **2,'ADDL  ',2,'XNATO ')
*. Natural orbitals in symmetry blocks
      CALL MEMMAN(KRHO1SM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KXNATSM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KOCCSM,NACOB ,'ADDL  ',2,'RHO1S ')
*
*. Space for two  blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*. Space for two blocks of string occupations and four arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      DO K12 = 1, 2
        CALL MEMMAN(KLOCSTR(K12),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      END DO
      DO I1234 = 1, 4
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
*. Arrays for partitioning of Left vector = sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTL ,NTTS  ,'ADDL  ',1,'LBT_L  ')
      CALL MEMMAN(KLLEBTL,NTTS  ,'ADDL  ',1,'LEBT_L ')
      CALL MEMMAN(KLI1BTL,NTTS  ,'ADDL  ',1,'I1BT_L ')
      CALL MEMMAN(KLIBTL ,8*NTTS,'ADDL  ',1,'IBT_L  ')
      CALL MEMMAN(KLSCLFCL,NTTS, 'ADDL  ',2,'SCLF_L')
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),
     &     WORK(KLI1BTL),WORK(KLIBTL),0,ITTSS_ORD)
*. Number of BLOCKS
        NBLOCKL = IFRMR(WORK(KLI1BTL),1,NBATCHL)
     &         + IFRMR(WORK(KLLBTL),1,NBATCHL) - 1
*. Arrays for partitioning of Right  vector = C
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTR ,NTTS  ,'ADDL  ',1,'LBT_R  ')
      CALL MEMMAN(KLLEBTR,NTTS  ,'ADDL  ',1,'LEBT_R ')
      CALL MEMMAN(KLI1BTR,NTTS  ,'ADDL  ',1,'I1BT_R ')
      CALL MEMMAN(KLIBTR ,8*NTTS,'ADDL  ',1,'IBT_R  ')
      CALL MEMMAN(KLSCLFCR,NTTS, 'ADDL  ',2,'SCLF_R')
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KCIOIO),ISMOST(1,ICSM),
     &     NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),
     &     WORK(KLI1BTR),WORK(KLIBTR),0,ITTSS_ORD)
*. Number of BLOCKS
        NBLOCKR = IFRMR(WORK(KLI1BTR),1,NBATCHR)
     &         + IFRMR(WORK(KLLBTR),1,NBATCHR) - 1
        WRITE(6,*) ' DENSI2T :NBLOCKR =',NBLOCKR
      WRITE(6,*) ' Memcheck in densi2 before GASDN2'
      CALL LMEMCHK('densi2 before GASDN2')
*
      WRITE(6,*) ' LUL LUR = ',LUL,LUR
      WRITE(6,*) ' LUSC1 LUSC2 ',LUSC1,LUSC2
      CALL REWINO(LUL)
      LENRHO1 = NTOOB ** 2
      ITRADEN = 0
      DO IL = 1, NL
        IF(IL.EQ.1) THEN
          LULEFF = LUL
        ELSE
          LULEFF = LUSC1
          CALL REWINO(LULEFF)
          CALL COPVCD(LUL,LULEFF,WORK(KLL),0,LBLK)
        END IF
        CALL REWINO(LUR)
        DO IR = 1, NR
          ITRADEN = ITRADEN + 1
          IF(IR.EQ.1) THEN
            LUREFF = LUR
          ELSE
            LUREFF = LUSC2
            CALL REWINO(LUREFF)
            CALL COPVCD(LUR,LUREFF,WORK(KLL),0,LBLK)
          END IF
*
          WRITE(6,*) ' Transition density will be called for : '
          WRITE(6,*) '             Left  state : ', IL
          WRITE(6,*) '             Right state : ', IR
*
      IF(ICISTR.EQ.1) THEN
         WRITE(6,*) ' Sorry, ICISTR = 1 is out of fashion'
         WRITE(6,*) ' Switch to ICISTR = 2 - or reprogram '
         Call Abend2( ' DENSI2T : ICISTR = 1 in use ' )
      ELSE IF(ICISTR.GE.2) THEN
        CALL GASDN2(I12,RHO1(1+(ITRADEN-1)*LENRHO1),RHO2,
     &       WORK(KLL),WORK(KLR),
     &       WORK(KLL),WORK(KLR),WORK(KC2),
     &       WORK(KCIOIO),WORK(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),WORK(KCBLTP),WORK(KSBLTP),NACOB,
     &       WORK(KNSTSO(IATP)),WORK(KISTSO(IATP)),
     &       WORK(KNSTSO(IBTP)),WORK(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,
     &       MAXK,MAXI,LBLOCK,LBLOCK,WORK(KCSCR),WORK(KSSCR),
     &       SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &       WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &       WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),WORK(KINSCR),
     &       MXPOBS,IPRDEN,WORK(KRHO1S),LULEFF,LUREFF,
     &       PSSIGN,PSSIGN,WORK(KRHO1P),WORK(KXNATO),
     &       NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),WORK(KLI1BTL),
     &       WORK(KLIBTL),
     &       NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),WORK(KLI1BTR),
     &       WORK(KLIBTR),WORK(KCONSPA),WORK(KCONSPB),
     &       WORK(KLSCLFCL),WORK(KLSCLFCR),S2_TERM1,IUSE_PH,IPHGAS,
     &       MXTSOB)
      END IF
*
      END DO
*     ^ End of loop over right states
      END DO
*     ^ End of loop over left states
*
      CALL LMEMCHK('TRADEN')
*
*. Eliminate local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'TRADEN')
      CALL QEXIT('TRADN')
      WRITE(6,*) ' Leaving TRADEN'
      RETURN
      END
