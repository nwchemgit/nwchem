*/ a number of clones of DENSI2 and two slave-routines which deliver
*/ the one and two-particle densities in spin-blocks
      SUBROUTINE DENSI2_AB
     &     (I12,RHO1,RHO2,L,R,LUL,LUR,EXPS2,IDOSRHO1,SRHO1,ISEPAB)
*
* Density matrices between L and R
*
* I12 = 1 => only one-bodydensity
* I12 = 2 => one- and two-body-density matrices
*
* Jeppe Olsen,      Oct 94
* GAS modifications Aug 95
* Two body density added, '96
*
* Table-Block driven, June 97
* Spin density added, Jan. 99
*
* Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
* ijkl = ij*(ij-1)/2+kl, ij.ge.kl
*
* If the twobody density matrix is calculated, then also the
* expectation value of the spin is evalueated.
* The latter is realixed as
* S**2 
*      = S+S- + Sz(Sz-1)
*      = -Sum(ij) a+i alpha a+j beta a i beta a j alpha + Nalpha +
*        1/2(N alpha - N beta))(1/2(N alpha - Nbeta) - 1)
* If IDOSRHO1 = 1, spin density is also calculated
c      IMPLICIT REAL*8(A-H,O-Z) 
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc' 
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
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
*
      INCLUDE 'csmprd.inc'
c      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
c      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
c     &              ADSXA(MXPOBS,2*MXPOBS),
c     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      INCLUDE 'lucinp.inc'
      INCLUDE 'clunit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Scratch for string information
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. Specific input 
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*),SRHO1(*)
*. Before I forget it :
      CALL LUCIAQENTER('DENSI')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
      ZERO = 0.0D0
      CALL SETVEC(RHO1,ZERO ,NACOB ** 2 )
      IF (ISEPAB.EQ.1)
     &  CALL SETVEC(RHO1(NOCOB**2+1),ZERO ,NACOB ** 2 )
      IF(I12.EQ.2) 
     &CALL SETVEC(RHO2,ZERO ,NACOB**2 *(NACOB**2+1)/2)
      IF (I12.EQ.2.AND.ISEPAB.EQ.1)
     &CALL SETVEC(RHO2(NACOB**2 *(NACOB**2+1)/2 + 1),ZERO ,
     &     NACOB**4 + NACOB**2 *(NACOB**2+1)/2)
*
      IF(IDOSRHO1.EQ.1) THEN  
        CALL SETVEC(SRHO1,ZERO,NACOB ** 2)
      END IF
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
      MAXA = MAX(MAXA,MAXA0)   
      MAXB = MAX(MAXB,MAXB0)   
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRDEN.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
*. replace with MXINKA !!!
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
C?    WRITE(6,*) ' DENSI2 : MAXI MAXK ', MAXI,MAXK
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(IPRDEN.GE.2) write(6,*) ' DENSI2 : MXSB MXTSOB MXSOOB ',
     &       MXSB,MXTSOB,MXSOOB 
      IF(ISIMSYM.NE.1) THEN
        LSCR1 = MXSOOB
      ELSE
        LSCR1 = MXSOOB_AS
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
      IF(ICISTR.EQ.1) THEN
        CALL MEMMAN(KCB,LSCR1,'ADDL  ',2,'KCB   ')
        CALL MEMMAN(KSB,LSCR1,'ADDL  ',2,'KSB   ')
      END IF
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
*
      CALL IAIBCM(ISSPC,WORK(KSIOIO))
      CALL IAIBCM(ICSPC,WORK(KCIOIO))
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESCPH(int_mb(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSTFSMSPGP,MXPNSMST,
     &     NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &     NELFSPGP,
     &     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &     IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,
     &     MX_NSPII)
      IF(IPRDEN.GE.2) THEN
        WRITE(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
      LSCR12 = MAX(LSCR1,2*LSCR2)
*. It is assumed that the third block already has been allocated, so
      KC2 = KVEC3
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL MEMMAN(KI1,  LSCR3       ,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KI2,  LSCR3       ,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KI3,  LSCR3       ,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KI4,  LSCR3       ,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
      CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
      CALL MEMMAN(KXI3S,LSCR3       ,'ADDL  ',2,'XI3S  ')
      CALL MEMMAN(KXI4S,LSCR3       ,'ADDL  ',2,'XI4S  ')
*. Arrays giving block type
COLD  CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
COLD  CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
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
      write(*,*) 'VOG MALAKIA? 2'
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
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),
     &     WORK(KLI1BTL),WORK(KLIBTL),0,ISIMSYM)
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
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),
     &     NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),
     &     WORK(KLI1BTR),WORK(KLIBTR),0,ISIMSYM)
*. Number of BLOCKS
        NBLOCKR = IFRMR(WORK(KLI1BTR),1,NBATCHR)
     &         + IFRMR(WORK(KLLBTR),1,NBATCHR) - 1
C?      WRITE(6,*) ' DENSI2T :NBLOCKR =',NBLOCKR

      IF(ICISTR.EQ.1) THEN
         WRITE(6,*) ' Sorry, ICISTR = 1 is out of fashion'
         WRITE(6,*) ' Switch to ICISTR = 2 - or reprogram '
         STOP' DENSI2T : ICISTR = 1 in use '
      ELSE IF(ICISTR.GE.2) THEN
        S2_TERM1 = 0.0D0
        CALL GASDN2_AB(I12,RHO1,RHO2,L,R,L,R,WORK(KC2),
     &       WORK(KCIOIO),WORK(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),WORK(KCBLTP),WORK(KSBLTP),NACOB,
     &       WORK(KNSTSO(IATP)),WORK(KISTSO(IATP)),
     &       WORK(KNSTSO(IBTP)),WORK(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,      
     &       MAXK,MAXI,LSCR1,LSCR1,WORK(KCSCR),WORK(KSSCR),
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
     &       IDOSRHO1,SRHO1,ISEPAB)
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
      CALL NATORB(RHO1,NSMOB,NTOOBS,NACOBS,NINOBS,
     &            IREOST,WORK(KXNATO),
     &            WORK(KRHO1SM),WORK(KOCCSM),NACOB,
     &            WORK(KRHO1P),IPRDEN)
*
      IF(IPRDEN.GE.5) THEN
        WRITE(6,*) ' One-electron density matrix '
        WRITE(6,*) ' ============================'
        CALL WRTMAT(RHO1,NTOOB,NTOOB,NTOOB,NTOOB) 
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-electron density '
          CALL PRSYM(RHO2,NACOB**2)
        END IF
      END IF
*
      IF(I12.EQ.2) THEN
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
*
      IF(IDOSRHO1.EQ.1.AND.IPRDEN.GE.2) THEN
        WRITE(6,*) ' One-electron spindensity <0!E(aa) - E(bb)!0> '
        CALL WRTMAT(SRHO1,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF

*. Eliminate local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'DENSI ')
      CALL LUCIAQEXIT('DENSI')
C     WRITE(6,*) ' Leaving DENSI '
      RETURN
      END
      SUBROUTINE GASDN2_AB(I12,RHO1,RHO2,
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
     &           IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,ISEPAB)
*
*
* Jeppe Olsen , Winter of 1991
* GAS modificatios, August 1995
*
* Table driven, June 97
*
* Last revision : Jan. 98 (IUSE_PH,IPHGAS added)
*                 Jan. 99 (IDOSRHO1,SRHO1 added)
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
      CALL LUCIAQENTER('GASDN')
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
      IF(NTEST.GE.1000) THEN
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
          ISCALE = 1
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
          ISCALE = 1
          IF(INTERACT.EQ.1) THEN
            ISCALE = 1
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
                     CALL GSDNBB2_AB(I12,RHO1,RHO2,
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
     &                    S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,ISEPAB)
                          IF(NTEST.GE.500) THEN
                            write(6,*) ' Updated rho1 '
                            call wrtmat(rho1,nacob,nacob,nacob,nacob)
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
      CALL LUCIAQEXIT('GASDN')
      RETURN
      END
      SUBROUTINE GSDNBB2_AB(I12,RHO1,RHO2,
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
     &                  SCLFAC,S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &                  ISEPAB)
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
      CALL LUCIAQENTER('GSDNB')
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      NTESTO= NTEST

      IRHO1A = 1
      IRHO1B = 1
      IRHO2AA = 1
      IRHO2AB = 1
      IRHO2BB = 1

      IF (ISEPAB.EQ.1) THEN
        IRHO1A = 1
        IRHO1B = NACOB**2 + 1
        IRHO2AA = 1
        IRHO2BB =    NACOB**2 *(NACOB**2+1)/2  + 1
        IRHO2AB = 2*(NACOB**2 *(NACOB**2+1)/2) + 1
      END IF

      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSDNBB2 :  R block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' GSDNBB2 :  L block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
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
        IF (ISEPAB.EQ.1) THEN
          WRITE(6,*) ' I will sepate spin-orbital contributions'
          WRITE(6,*) ' offsets: '
          WRITE(6,*) ' RHO1 :  ', IRHO1A, IRHO1B
          IF (I12.EQ.2)
     &      WRITE(6,*) ' RHO2 :  ', IRHO2AA, IRHO2BB, IRHO2AB
        END IF

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
        IAB = 2
        CALL GSBBD1(RHO1(IRHO1B),NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &       NGAS,IBOC,JBOC,
     &       SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB)
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
          CALL GSBBD2A(RHO2(IRHO2BB),
     &         NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &         NGAS,IBOC,JBOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,0,XDUM)
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
        IAB = 1
        CALL GSBBD1(RHO1(IRHO1A),NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &       NGAS,IAOC,JAOC,SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB)
C?        WRITE(6,*) ' GSBBD1 was called '
        IF(I12.EQ.2.AND.NAEL.GE.2) THEN
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
C?        WRITE(6,*) ' GSBBD2A will be called (alpha)'
          CALL GSBBD2A(RHO2(IRHO2AA),
     &         NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &         NGAS,IAOC,JAOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,0,XDUM)
C?        WRITE(6,*) ' GSBBD2A was called '
        END IF
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIB,NIA,C2)
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
        CALL GSBBD2B_AB(RHO2(IRHO2AB),IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                    JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                    IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                    SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                    NOBPTS,IOBPTS,MAXK,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                    NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                    SSCR,CSCR,NACOB,NTEST,SCLFAC,S2_TERM1,
     &                    ISEPAB)
C?      WRITE(6,*) ' GSBBD2B was called '
     &                    
C     GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
C    &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
C    &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
C    &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
C    &                  NOBPTS,IOBPTS,MAXK,
C    &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
C    &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C    &                  CJRES,SIRES,NORB,NTEST)
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
      CALL LUCIAQEXIT('GSDNB')
      RETURN
      END
      SUBROUTINE GSBBD2B_AB(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXK,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  CJRES,SIRES,NORB,NTESTG,SCLFAC,S2_TERM1,
     &                  ISEPAB)
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
      CALL LUCIAQENTER('GSD2B')
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
     &               I4,XI4S,I2,XI2S,IKORD)
*. contributions to Rho2(ij,kl) has been obtained, scatter out
C?              WRITE(6,*) ' Before call to ADTOR2'
C?              WRITE(6,*) ' RHO2B (X) matrix '
C?              call wrtmat(x,ni*nj,nk*nl,ni*nj,nk*nl)
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
                CALL ADTOR2_AB(RHO2,X,2,ISEPAB,
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
      CALL LUCIAQEXIT('GSD2B')
      RETURN
      END
      SUBROUTINE ADTOR2_AB(RHO2,RHO2T,ITYPE,ISEPAB,
     &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
*
* Add contributions to two electron density matrix RHO2
* output density matrix is in the form Rho2(ij,kl),(ij).ge.(kl)
*
* if ISEPAB.EQ.1 Rho2 indices are unrestricted
*
* Jeppe Olsen, Fall of 96
*
*
* Itype = 1 => alpha-alpha or beta-beta loop
*              input is in form Rho2t(ik,jl)
* Itype = 2 => alpha-beta loop
*              input is in form Rho2t(ij,kl)
*               
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION RHO2T(*)
*. Input and output
      DIMENSION RHO2(*)
*
      NTEST = 0000
      IF (ISEPAB.EQ.1.AND.ITYPE.NE.2) THEN
        WRITE(6,*) 'Illegal combination of ISEPAB and ITYPE in ADTOR!',
     &       ISEPAB, ITYPE
        STOP 'ADTOR2'
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to ADTOR2 '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' NI NJ NK NL = ', NI,NJ,NK,NL
        WRITE(6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
        WRITE(6,*) ' ITYPE = ',ITYPE
        IF(NTEST.GE.2000) THEN
          WRITE(6,*) ' Initial two body density matrix '
          CALL PRSYM(RHO2,NORB**2)
        END IF
        WRITE(6,*) ' RHO2T : '
        IF(ITYPE.EQ.1) THEN
          IF(IOFF.EQ.KOFF) THEN
            NROW = NI*(NI+1)/2
          ELSE
            NROW = NI*NK
          END IF
          IF(JOFF.EQ.LOFF) THEN
            NCOL = NJ*(NJ+1)/2
          ELSE
            NCOL = NJ*NL
          END IF
        ELSE IF (ITYPE.EQ.2) THEN
          NROW = NI*NJ
          NCOL = NK*NL
        END IF
        CALL WRTMAT(RHO2T,NROW,NCOL,NROW,NCOL)
      END IF
C?    WRITE(6,*) ' Enforced return in ADTOR2 '
C?    RETURN
      NELMNT = NORB**2*(NORB**2+1)/2
*
      IF(ITYPE.EQ.1) THEN
*
* =======================================
*     Alpha-alpha or beta-beta term
* =======================================
*
*. Four permutations
      DO IPERM = 1, 4
        IF(IPERM.EQ.1) THEN
          NII = NI
          IIOFF = IOFF
          NJJ = NJ
          JJOFF = JOFF
          NKK = NK
          KKOFF = KOFF
          NLL = NL
          LLOFF = LOFF
          SIGN = 1.0D0
          IACTIVE = 1
        ELSE IF(IPERM.EQ.2) THEN
          IF(IOFF.NE.KOFF) THEN
            NII = NK
            IIOFF = KOFF
            NKK = NI
            KKOFF = IOFF
            NJJ = NJ
            JJOFF = JOFF
            NLL = NL
            LLOFF = LOFF
            IACTIVE = 1
          ELSE 
            IACTIVE = 0
          END IF
          SIGN = -1.0D0
        ELSE IF(IPERM.EQ.3) THEN
          IF(JOFF.NE.LOFF) THEN
            NII = NI
            IIOFF = IOFF
            NKK = NK
            KKOFF = KOFF
            NJJ = NL
            JJOFF = LOFF
            NLL = NJ
            LLOFF = JOFF
            SIGN = -1.0D0
            IACTIVE = 1
          ELSE
            IACTIVE = 0
          END IF
        ELSE IF(IPERM.EQ.4) THEN
          IF(IOFF.NE.KOFF.AND.JOFF.NE.LOFF) THEN
            NKK = NI
            KKOFF = IOFF
            NII = NK
            IIOFF = KOFF
            NJJ = NL
            JJOFF = LOFF
            NLL = NJ
            LLOFF = JOFF
            SIGN = 1.0D0
            IACTIVE = 1
          ELSE
            IACTIVE = 0
          END IF
        END IF
*
        IJOFF = (JJOFF-1)*NORB+IIOFF
        KLOFF = (LLOFF-1)*NORB+KKOFF
C       IF(IACTIVE.EQ.1.AND.IJOFF.GE.KLOFF) THEN
        IF(IACTIVE.EQ.1) THEN
          IJOFF = (JJOFF-1)*NORB+IIOFF
          KLOFF = (LLOFF-1)*NORB+LLOFF
            DO II = 1, NII
              DO JJ = 1, NJJ
                DO KK = 1, NKK
                  DO LL = 1, NLL
                    IJ = (JJ+JJOFF-2)*NORB + II+IIOFF - 1
                    KL = (LL+LLOFF-2)*NORB + KK+KKOFF - 1
                    IF(IJ.GE.KL) THEN
                      IJKL = IJ*(IJ-1)/2+KL
                      IF(IPERM.EQ.1) THEN
                        I = II
                        K = KK
                        J = JJ
                        L = LL
                      ELSE IF(IPERM.EQ.2) THEN
                        I = KK
                        K = II
                        J = JJ
                        L = LL
                      ELSE IF(IPERM.EQ.3) THEN
                        I = II
                        K = KK
                        J = LL
                        L = JJ
                      ELSE IF(IPERM.EQ.4) THEN
                        I = KK
                        K = II
                        J = LL
                        L = JJ
                      END IF
                      IF(IOFF.NE.KOFF) THEN
                        IKIND = (K-1)*NI+I
                        NIK = NI*NK
                        SIGNIK = 1.0D0
                      ELSE
                        IKIND = MAX(I,K)*(MAX(I,K)-1)/2+MIN(I,K)
                        NIK = NI*(NI+1)/2
                        IF(I.EQ.MAX(I,K)) THEN
                          SIGNIK = 1.0D0
                        ELSE
                          SIGNIK = -1.0D0
                        END IF
                      END IF
                      IF(JOFF.NE.LOFF) THEN
                        JLIND = (L-1)*NJ+J
                        SIGNJL = 1.0D0
                      ELSE
                        JLIND = MAX(J,L)*(MAX(J,L)-1)/2+MIN(J,L)
                        IF(J.EQ.MAX(J,L)) THEN
                          SIGNJL = 1.0D0
                        ELSE
                          SIGNJL = -1.0D0
                        END IF
                      END IF
                      IKJLT = (JLIND-1)*NIK+IKIND
                      IF(IJKL.GT.NELMNT) THEN
                         WRITE(6,*) ' Problemo 1 : IJKL .gt. NELMNT'
                         WRITE(6,*) ' IJKL, NELMNT',IJKL,NELMNT
                         WRITE(6,*) ' IJ, KL', IJ,KL
                         WRITE(6,*) ' JJ JJOFF ', JJ,JJOFF
                         WRITE(6,*) ' II IIOFF ', II,IIOFF
                         WRITE(6,*) ' IPERM = ', IPERM
                      END IF
                      RHO2(IJKL) = RHO2(IJKL) 
     &                           - SIGN*SIGNJL*SIGNIK*RHO2T(IKJLT)
*. The minus : Rho2t comes as <a+i a+k aj al>, but we want 
* <a+ia+k al aj>
                    END IF
                  END DO
                END DO
              END DO
            END DO
*. End of active/inactive if
        END IF
*. End of loop over permutations
      END DO
      ELSE IF(ITYPE.EQ.2.AND.ISEPAB.EQ.0) THEN
*
* =======================================
*     Alpha-beta term
* =======================================
*
      DO I = 1, NI
       DO J = 1, NJ
         DO K = 1, NK
           DO L = 1, NL
             IJ = (J+JOFF-2)*NORB + I+IOFF - 1
             KL = (L+LOFF-2)*NORB + K+KOFF - 1
             IF(IJ.EQ.KL) THEN
               FACTOR = 2.0D0
             ELSE 
               FACTOR= 1.0D0
             END IF
             IJKL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
             IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI
     &             + (J-1)*NI + I
                      IF(IJKL.GT.NELMNT) THEN
                         WRITE(6,*) ' Problemo 2 : IJKL .gt. NELMNT'
                         WRITE(6,*) ' IJKL, NELMNT',IJKL,NELMNT
                         STOP 'ADTOR2_AB (1)'
                      END IF
             RHO2(IJKL) = RHO2(IJKL) + FACTOR*RHO2T(IJKLT)
            END DO
          END DO
        END DO
      END DO
*
      ELSE IF(ITYPE.EQ.2.AND.ISEPAB.EQ.1) THEN
*
* =======================================
*     Alpha-beta term (full matrix)
* =======================================
*
      NELMNT = NORB**4
      DO I = 1, NI
       DO J = 1, NJ
         DO K = 1, NK
           DO L = 1, NL
             IJ = (J+JOFF-2)*NORB + I+IOFF - 1
             KL = (L+LOFF-2)*NORB + K+KOFF - 1
c I think:
c             IF(IJ.EQ.KL) THEN
c               FACTOR = 2.0D0
c             ELSE 
               FACTOR= 1.0D0
c             END IF
             IJKL = (KL-1)*NORB**2 + IJ
             IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI
     &             + (J-1)*NI + I
                      IF(IJKL.GT.NELMNT) THEN
                         WRITE(6,*) ' Problemo 3 : IJKL .gt. NELMNT'
                         WRITE(6,*) ' IJKL, NELMNT',IJKL,NELMNT
                         STOP 'ADTOR2_AB (2)'
                      END IF
             RHO2(IJKL) = RHO2(IJKL) + FACTOR*RHO2T(IJKLT)
            END DO
          END DO
        END DO
      END DO
*
      ELSE
*
        WRITE(6,*) 'Illegal ITYPE, ISEPAB combination in ADTOR2_AB: ',
     &       ITYPE,ISEPAB
        STOP 'ADTOR2_AB'
*
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Updated two-body density matrix '
         CALL PRSYM(RHO2,NORB**2)
      END IF
*
      RETURN
      END
c $Id$
