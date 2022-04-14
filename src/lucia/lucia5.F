C EXP_BL_MAT
      FUNCTION SUMELM(VEC,NDIM)
*
* Sum elements in Vector
*
      INCLUDE 'implicit.inc'
*
      DIMENSION VEC(NDIM)
*
      SUM = 0.0D0
      DO I = 1, NDIM
        SUM = SUM + VEC(I)
      END DO
      SUMELM = SUM
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Sum of elemnts = ', SUM
      END IF
*
      RETURN
      END
      SUBROUTINE GET_D2_SMBLK(D2,ISM,JSM,KSM,LSM)
*.
*. Obtain  complete symmetry block of D2
*
*. Jeppe Olsen, Dec.4 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Output
      DIMENSION D2(*)
*. Loop over orbitals in symmetry order
      NI = NTOOBS(ISM)
      NJ = NTOOBS(JSM)
      NK = NTOOBS(KSM)
      NL = NTOOBS(LSM)
* 
      IB = IBSO(ISM)
      JB = IBSO(JSM)
      KB = IBSO(KSM)
      LB = IBSO(LSM)
*
      IJKL = 0
      DO L = LB, LB + NL - 1
       DO K = KB, KB + NK - 1
        DO J = JB, JB + NJ - 1
         DO I = IB, IB + NI - 1
*
           IGAS = ITPFSO(I)
           JGAS = ITPFSO(J)
           KGAS = ITPFSO(K)
           LGAS = ITPFSO(L)
*
           IREL = IREOST(I)-IOBPTS(IGAS,ISM)+1
           JREL = IREOST(J)-IOBPTS(JGAS,JSM)+1
           KREL = IREOST(K)-IOBPTS(KGAS,KSM)+1
           LREL = IREOST(L)-IOBPTS(LGAS,LSM)+1
*
           IJKL = IJKL + 1
           D2(IJKL) = GETD2E(IREL,IGAS,ISM,JREL,JGAS,JSM,
     &                       KREL,KGAS,KSM,LREL,LGAS,LSM,1)
         END DO
        END DO
       END DO
      END DO
*
C?    SUM = 0.0D0
C?    NIJKL = NI*NJ*NK*NL 
C?    DO IJKL = 1, NIJKL
C?      SUM = SUM + D2(IJKL)
C?    END DO
C?    WRITE(6,*) ' Sum = ', SUM
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Rho2 block for ISM JSM KSM LSM =',
     &               ISM,JSM,KSM,LSM
        NIJ = NI*NJ
        NKL = NK*NL
        CALL WRTMAT(D2,NIJ,NKL,NIJ,NKL)
      END IF
*
      RETURN
      END 
      SUBROUTINE MERGE_PH_MAT(PMAT,HMAT,PHMAT,NSM,NTOOBS,IHSM)
*
* Merge a particle and a hole matrix to a single matrix
*
* Part of QDOT interface
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*.input
      INTEGER NTOOBS(NSM)
      DIMENSION HMAT(*),PMAT(*)
*. Output 
      DIMENSION PHMAT(*)
*
C     NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN_COL = NDIM_1EL_MAT(IHSM,NTOOBS,NTOOBS,NSM,0)
      ZERO = 0.0D0
      CALL SETVEC(PHMAT,ZERO,LEN_COL)
*
      IOFF_SEP = 1
      IOFF_COL = 1
      DO ISM = 1, NSM
        JSM = MULTD2H(ISM,IHSM)
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)
        NID2 = NI/2
        NJD2 = NJ/2
        DO I = 1, NID2
        DO J = 1, NJD2
*. Hole integrals
          PHMAT(IOFF_COL-1+(J-1)*NI +  I) 
     &  = HMAT(IOFF_SEP-1+(J-1)*NID2+I)
*. Particle integrals
          PHMAT(IOFF_COL-1+(J+NJD2-1)*NI +  I+NID2) 
     &  = PMAT(IOFF_SEP-1+(J-1)*NID2+I)
        END DO
        END DO
*
        IOFF_COL = IOFF_COL + NI*NJ
        IOFF_SEP = IOFF_SEP + NID2*NJD2
      END DO
* 
      RETURN
      END
C      CALL GET_H1AO_QDOT(LABEL,SCR(KLH1AO),IHSM,NBAS,IPERMSM)
      SUBROUTINE GET_H1AO_FROM_DISC_QDOT(LUH,H1AO,IHSM,NBAS,NSM)
*
* Obtain set of property integrals from file written by QDOT 
* program
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NBAS(NSM)
*. Output
      DIMENSION H1AO(*)
      IJ = 0
      IOFF = 1
      DO ISM = 1, NSM
        JSM = MULTD2H(ISM,IHSM)
        NI = NBAS(ISM)
        NJ = NBAS(JSM)
C?      WRITE(6,*) ' ISM JSM NI NJ = ', ISM,JSM,NI,NJ
*. Skip line with symmetry info
        READ(LUH,*)
        DO J = 1, NJ
         DO I = 1, NI
          READ(LUH,*) H1AOXX, II,JJ
C?        WRITE(6,*) ' H1AOXX, II,JJ ', H1AOXX, II,JJ
          H1AO(IOFF-1+(JJ-1)*NI + II) = H1AOXX
         END DO
        END DO
        IOFF = IOFF + NI*NJ
      END DO
*
      RETURN
      END 
C      CALL GET_H1AO_QDOT(LABEL,SCR(KLH1AO),IHSM,NBAS,IPERMSM)
      SUBROUTINE  GET_H1AO_QDOT(LABEL,H1AO,IHSM,NBAS,IPERMSM)
*
* Obtain one-electron integrals from QDOT environment
*
*. Integrals are stored in symmetrypacked form, complete
*. matrices wtih holes first and then particles
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      CHARACTER*8 LABEL
*. Input
      INTEGER NBAS(NSMOB)
*. Output
      DIMENSION H1AO(*)
*. Local scratch
      INTEGER NBASD2(8)
*. Number of basis functions for holes or particles
      DO I = 1, NSMOB
       NBASD2(I) = NBAS(I)/2
      END DO
*
      LUHAO = 73
      IF(LABEL(1:7).EQ.'ZANGMOM') THEN
*. Obtain hole LZ integrals
C       GET_H1AO_FROM_DISC_QDOT(LUH,H1AO,IHSM,NBAS,NSM)
        OPEN(UNIT=LUHAO,FILE='LZ_HOLE',FORM='FORMATTED',
     &       STATUS='OLD')
C?      WRITE(6,*) ' LZ_HOLE opened '
        CALL GET_H1AO_FROM_DISC_QDOT(LUHAO,H1AO,IHSM,NBASD2,NSMOB)
        CLOSE(LUHAO)
*. Particle integrals
C       NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
        LEN = NDIM_1EL_MAT(IHSM,NBASD2,NBASD2,NSMOB,0)
        OPEN(UNIT=LUHAO,FILE='LZ_ELECTRON',FORM='FORMATTED',
     &       STATUS='OLD')
        CALL GET_H1AO_FROM_DISC_QDOT(LUHAO,H1AO(1+LEN),IHSM,NBASD2,
     &       NSMOB)
        CLOSE(LUHAO)
      ELSE
        WRITE(6,'(A,A)') ' Unknown label of integrals : ', LABEL
        STOP             ' Unknown label of integrals '
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,A)') ' Label of integrals = ', LABEL
        WRITE(6,*) ' Hole integrals '
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
        CALL PRHONE(H1AO,NBASD2,IHSM,NSMOB,0)
        WRITE(6,*) ' Particle integrals '
        CALL PRHONE(H1AO(1+LEN),NBASD2,IHSM,NSMOB,0)
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_L3BLKS(LVEC1,LVEC2,LC2)
*
* Get length of the three blocks VEC1, VEC2 and C2 used in sigma ..
*
* Jeppe Olsen, Sept. 2001 from GET_3BLKS
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'csm.inc' 
      INCLUDE 'cstate.inc' 
      INCLUDE 'crun.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'gasstr.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
*. Common block for communicating with sigma
      INCLUDE 'cands.inc'
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_L3')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*
      ISPC = MAX(ICSPC,ISSPC)
      ISM  = ISSM
      CALL IAIBCM(ISPC,WORK(KLCIOIO))
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTP),WORK(KSVST))
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR           
*. alpha and beta strings with an electron removed
      IATPM1 = 3 
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5 
      IBTPM2 = 6
*. Largest number of strings of given symmetry and type
      MAXA = MXNSTR
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
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESCPH(int_mb(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &              NSMST,NSTFSMSPGP,MXPNSMST,
     &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &              NELFSPGP,
     &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &              IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,
     &              MXADKBLK_AS,MX_NSPII)
      IF(ISIMSYM.EQ.1) MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
      IF(IPRCIX.GE.2) THEN
        WRITE(6,*) 'GET_3BL: MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
        WRITE(6,*) 'GET_3BL: MXADKBLK ', MXADKBLK
        WRITE(6,*) ' Space for resolution matrices ',LSCR2
      END IF
*
      IF(ISIMSYM.EQ.0) THEN 
        LBLOCK = MXSOOB
      ELSE
        LBLOCK = MXSOOB_AS
      END IF
*
      LBLOCK = MAX(LBLOCK,LCSBLK)
      LSCR12 = MAX(LBLOCK,2*LSCR2)  
*
      LVEC1 = LBLOCK
      LVEC2 = LBLOCK
      LC2 = LSCR12
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_L3')
*
      RETURN
      END
*
      SUBROUTINE GET_PROP_PERMSM(PROPER,IPERMSM)
*
* Obtain permutational symmetry for operator PROPER
*
* Jeppe Olsen, Feb. 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 PROPER
*
      IF(PROPER(1:6).EQ.'DIPOLE'.OR.
     &   PROPER(1:6).EQ.'DIPLEN'.OR.
     &   PROPER(2:7).EQ.'DIPLEN') THEN
        IPERMSM = 1
      ELSE IF (PROPER(1:6).EQ.'ANGMOM'.OR.
     &         PROPER(2:7).EQ.'ANGMOM') THEN
        IPERMSM = -1
      ELSE IF(PROPER(1:6).EQ.'THETA ' .OR.
     &        PROPER(1:6).EQ.'QUADRU' .OR.
     &        PROPER(1:6).EQ.'SECMOM' .OR.
     &        PROPER(1:3).EQ.'EFG' ) THEN
        IPERMSM = 1
      ELSE
        WRITE(6,'(A,A)') ' Unknown operator ',PROPER
        IPERMSM = 0
      END IF
*
      NTEST = 0
      IF(NTEST.GE.5) THEN
        WRITE(6,'(A,A,I3)') ' Property and perm-sym : ', PROPER,IPERMSM
      END IF
*
      RETURN
      END
      SUBROUTINE EXP_LZ2(KVEC1,KVEC2,RLZEFF,RL2EFF,LUCEFF)
*
* Calculate Expectation value of Lz and L2
*
* Jeppe Olsen, Feb12 2000 - In Dage's living room with the cat sleeping
* Last revision: Feb. 2013: LZEFF and L2EFF returned as arguments
*. Last revision; Feb. 20, 2013; Jeppe Olsen; Opened for use of evaluation using CI
*
* KVEC1, KVEC2, are vectors or subblocks of vectors
* IF ICISTR = 1, then KVEC1 contains input vector
*
* In the present routine the expectation value of LZ is evaluated in 
* two ways :
*
* 1 : Using CI ( assuming space is closed under the action of L_i, i=x,y,z) 
* 2 : Using Density matrices
*
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'multd2h.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cicisp.inc'
*
      CHARACTER*8 LABEL
      REAL*8 INPRDD, INPROD
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' info form EXP_LZ2  '
        WRITE(6,*) ' ================== '
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'EXP_LZ')
*. Pointers to pointers
      KKVEC1 = KVEC1
      KKVEC2 = KVEC2 
*
      CALL MEMMAN(KLZINT,NTOOB**2,'ADDL  ',2,'LZINT ') !done
      CALL MEMMAN(KLSCR,3*NTOOB**2,'ADDL  ',2,'LSCR  ') !done
*
* CI OR density matrix approach:
      I_DO_CI_LI = 1
* Preparations: If CI approach is in use and CSF's are employed, transform 
* to CSF- we do not CSF's for several symmetries pt.
      IF(I_DO_CI_LI.EQ.1) THEN
*
*. Transform to SD's if requested - we do not have CSF expansion
*. for arbitrary symmetry
*
        NOCSF_SAVE = NOCSF
        ICNFBAT_SAVE = ICNFBAT
        LUC2 = LUCEFF
        IF(NOCSF.EQ.0) THEN
         IF(ICNFBAT.EQ.1) THEN
*. In core
           CALL CSDTVCM(KKVEC1,
     &          KCOMVEC1_SD,KCOMVEC2_SD,1,0,ICSM,ICSPC,2)
*. Final vector in now in KCOMVEC1_SD and we need sigma in the SD basis, so
*. use 
           KKVEC1 = KCOMVEC1_SD
           KKVEC2 = KCOMVEC2_SD
         ELSE
*. Not in core  write determinant expansion on scratch unit LUC_SD
          CALL FILEMAN_MINI(LUC_SD,'ASSIGN')
cnw       CALL REWINO(LUC_SD) 
cnw       CALL REWINO(LUCEFF)
C         CSDTVCMN(CSFVEC,DETVEC,SCR,IWAY,ICOPY,ISYM,ISPC,
C    &             IMAXMIN_OR_GAS,ICNFBAT,LUC_DET,LU_CSF,NOCCLS_ACT,
C    &             IOCCLS_ACT,IBLOCK,NBLK_PER_BATCH)
*
C?        WRITE(6,*) ' Bef. CSDTVCMN, LUC_SD, LUCEFF = ', LUC_SD, LUCEFF
C?        WRITE(6,*) ' NCOCCLS == ', NCOCCLS
C?        WRITE(6,*) ' KCLBT, WORK(KCLBT) = ', KCLBT
C?        CALL IWRTMA(WORK(KCLBT),1,1,1,1)
*
          CALL CSDTVCMN(KKVEC1,KKVEC2,dbl_mb(KVEC3),
     &         1,0,ICSM,ICSPC,2,2,LUC_SD,LUCEFF,NCOCCLS,
     &         int_mb(KCIOCCLS_ACT),int_mb(KCIBT),int_mb(KCLBT))
          LUC2 = LUC_SD
C?        WRITE(6,*) ' After CSDTVCMN '
         END IF
         NOCSF = 1
*. Dirty, to get CSF order of blocks
         ICNFBAT = -2
        END IF
*. We now have CI vector in SD/CMB basis on LUC2
      END IF ! CI is used
*
*
* =====================================
*. Obtain LZ integrals in MO basis 
* =====================================
*
      LZ_SYM = MULTD2H(IXYZSYM(1),IXYZSYM(2))
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LENIN = NDIM_1EL_MAT(LZ_SYM,NTOOBS,NTOOBS,NSMOB,0)
      LABEL(1:6) = 'ANGMOM'
C         GET_PROP_PERMSM(PROPER,IPERMSM)
      CALL GET_PROP_PERMSM(LABEL,IPERMSM)
C?    WRITE(6,*) ' IPERMSM from GET_PROP_PERMSM =', IPERMSM
*. (Well I know IPERMSM should be -1, here just testing )
C     GET_PROPINT(H,IHSM,LABEL,SCR,NMO,NBAS,NSM,ILOW,IPERMSM)
      LABEL(1:8) = 'ZANGMOM '
      CALL GET_PROPINT(dbl_mb(KLZINT),LZ_SYM,LABEL,dbl_mb(KLSCR),NTOOBS,
     &             NTOOBS,NSMOB,0,IPERMSM)
*
      IF(I_DO_CI_LI.EQ.1) THEN
* Lz * First vector on LUC
*. Integrals will be stored in symmetryblocked form 
*. without permutational symmetry packing
        IH1FORM_SAVE = IH1FORM
        IH1FORM = 2
        ICSM = IREFSM
        ISSM = MULTD2H(LZ_SYM,ICSM)
*. For later use
        IF(NOCSF.EQ.0) THEN
          LEN = NCM_PER_SYM_GN(ISSM, ICSPC)
        ELSE
          LEN =  XISPSM(ISSM,ICSPC)
        END IF
        IF(NTEST.GE.100) WRITE(6,*) ' ISSM = ', ISSM
C            BVEC(B,IBSM,LUC,LUB,WORK(VEC1),WORK(VEC2))
        CALL BVEC(dbl_mb(KLZINT),LZ_SYM,LUC2,LUHC,
     &       WORK(KKVEC1),WORK(KKVEC2))
*
        IH1FORM = IH1FORM_SAVE
        ISSM = IREFSM
* <C|Lz^2|C>
        IF(ICISTR.EQ.1) THEN
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' LEN = ', LEN
            WRITE(6,*) ' Vec2 as returned from BVEC '
            CALL WRTMAT(WORK(KKVEC2),1,LEN,1,LEN)
          END IF
          EXPV_LZ2 = INPROD(WORK(KKVEC2),WORK(KKVEC2),LEN)
        ELSE
          EXPV_LZ2 = INPRDD(WORK(KKVEC1),WORK(KKVEC2),LUHC,LUHC,1,-1)
        END IF
C?      WRITE(6,*) ' EXPV_LZ2 = ', EXPV_LZ2
        RLZEFF  = SQRT((ABS(EXPV_LZ2)))
*
      ELSE
*
* =================================
*. Evaluate using density matrices 
* =================================
*
C       ABEXP(A,IASM,B,IBSM,AB)
C       CALL ABEXP(WORK(KLZINT),LZ_SYM,WORK(KLZINT),LZ_SYM,EXP_LZ2)
        CALL ABEXP2(dbl_mb(KLZINT),LZ_SYM,dbl_mb(KLZINT),
     &              LZ_SYM,EXPV_LZ2)
*. Well, the LZ integrals are without the imaginary i, so
*. multiply with -1
        EXPV_LZ2 = -EXPV_LZ2
        RLZEFF = SQRT(ABS(EXPV_LZ2))
        IF(NTEST.GE.10) 
     &    WRITE(6,*) ' Expectation value Sqrt<Lz^2> from Densi:', RLZEFF
      END IF
*
*. Perform also calculation for Lx^2 and Ly^2
*
*
*. Lx:
*
*. Lx integrals in MO basis 
      LX_SYM = MULTD2H(IXYZSYM(2),IXYZSYM(3))
      LENIN = NDIM_1EL_MAT(LX_SYM,NTOOBS,NTOOBS,NSMOB,0)
      LABEL(1:6) = 'ANGMOM'
      CALL GET_PROP_PERMSM(LABEL,IPERMSM)
*. (Well I know IPERMSM should be -1, here just testing )
      LABEL(1:8) = 'XANGMOM '
      CALL GET_PROPINT(dbl_mb(KLZINT),LX_SYM,LABEL,dbl_mb(KLSCR),NTOOBS,
     &               NTOOBS,NSMOB,0,IPERMSM)
*. and the expectation value
      IF(I_DO_CI_LI.EQ.1) THEN
* Lx * First vector on LUC
*. Integrals will be stored in symmetryblocked form 
*. without permutational symmetry packing
        IH1FORM_SAVE = IH1FORM
        IH1FORM = 2
        ICSM = IREFSM
        ISSM = MULTD2H(LX_SYM,ICSM)
*. For later use
        IF(NOCSF.EQ.0) THEN
          LEN = NCM_PER_SYM_GN(ISSM, ICSPC)
        ELSE
          LEN =  XISPSM(ISSM,ICSPC)
        END IF
        IF(NTEST.GE.100) WRITE(6,*) ' ISSM = ', ISSM
C            BVEC(B,IBSM,LUC,LUB,WORK(VEC1),WORK(VEC2))
        CALL BVEC(dbl_mb(KLZINT),LX_SYM,LUC2,LUHC,
     &       WORK(KKVEC1),WORK(KKVEC2))
        IH1FORM = IH1FORM_SAVE
        ISSM = IREFSM
* <C|Lx^2|C>
        IF(ICISTR.EQ.1) THEN
          EXPV_LX2 = INPROD(WORK(KKVEC2),WORK(KKVEC2),LEN)
        ELSE
          EXPV_LX2 = INPRDD(WORK(KKVEC1),WORK(KKVEC2),LUHC,LUHC,1,-1)
        END IF
      ELSE
        CALL ABEXP2(dbl_mb(KLZINT),LX_SYM,dbl_mb(KLZINT),LX_SYM,
     &              EXPV_LX2)
        EXPV_LX2 = -EXPV_LX2
C?      WRITE(6,*) ' Expectation value of Lx^2 from Densi:', EXPV_LX2
      END IF
*
*. Ly:
*
*. Ly integrals in MO basis 
        LY_SYM = MULTD2H(IXYZSYM(1),IXYZSYM(3))
        LENIN = NDIM_1EL_MAT(LY_SYM,NTOOBS,NTOOBS,NSMOB,0)
        LABEL(1:6) = 'ANGMOM'
        CALL GET_PROP_PERMSM(LABEL,IPERMSM)
*. (Well I know IPERMSM should be -1, here just testing )
        LABEL(1:8) = 'YANGMOM '
        CALL GET_PROPINT(dbl_mb(KLZINT),LY_SYM,LABEL,dbl_mb(KLSCR),
     &               NTOOBS,NTOOBS,NSMOB,0,IPERMSM)
*. and the expectation value
      IF(I_DO_CI_LI.EQ.1) THEN
* Ly * First vector on LUC
*. Integrals will be stored in symmetryblocked form 
*. without permutational symmetry packing
        IH1FORM_SAVE = IH1FORM
        IH1FORM = 2
        ICSM = IREFSM
        ISSM = MULTD2H(LY_SYM,ICSM)
*. For later use
        IF(NOCSF.EQ.0) THEN
          LEN = NCM_PER_SYM_GN(ISSM, ICSPC)
        ELSE
          LEN =  XISPSM(ISSM,ICSPC)
        END IF
        IF(NTEST.GE.100) WRITE(6,*) ' ISSM = ', ISSM
C            BVEC(B,IBSM,LUC,LUB,WORK(VEC1),WORK(VEC2))
        CALL BVEC(dbl_mb(KLZINT),LY_SYM,LUC2,LUHC,
     &       WORK(KKVEC1),WORK(KKVEC2))
        IH1FORM = IH1FORM_SAVE
        ISSM = IREFSM
* <C|Ly^2|C>
        IF(ICISTR.EQ.1) THEN
          EXPV_LY2 = INPROD(WORK(KKVEC2),WORK(KKVEC2),LEN)
        ELSE
          EXPV_LY2 = INPRDD(WORK(KKVEC1),WORK(KKVEC2),LUHC,LUHC,1,-1)
        END IF
        RLYEFF  = SQRT((ABS(EXPV_LY2)))
      ELSE
        CALL ABEXP2(dbl_mb(KLZINT),LY_SYM,dbl_mb(KLZINT),LY_SYM,
     &              EXPV_LY2)
*. Well, the Ly integrals are without the imaginary i, so
*. multiply with -1
        EXPV_LY2 = -EXPV_LY2
C?      WRITE(6,*) ' Expectation value of Ly^2 from Densi:', EXPV_LY2
      END IF
*
* And the finale
*
      EXPV_L2 =  EXPV_LX2 + EXPV_LY2 + EXPV_LZ2
      RL2EFF = EXPV_L2
*
*
* Clean up time
      IF(NOCSF_SAVE.EQ.0) THEN
          CALL FILEMAN_MINI(LUC_SD,'FREE  ')
          NOCSF = 0
          ICNFBAT = ICNFBAT_SAVE
      END IF
*
      IF(NTEST.GE.10)  THEN
         WRITE(6,*) ' Expectation value of |L_z!:', RLZEFF
         WRITE(6,*) ' Expectation value of L^2 :',  RL2EFF
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'EXP_LZ')
*
      RETURN
      END
      SUBROUTINE  CLS_TO_BASE(CLS_E,EBASC,CLS_C,CBASC,NCLS,NSPC,IBASSPC)
*
* Class info => Base space info
* for energy abd wf correction 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBASSPC(NCLS)
      DIMENSION CLS_C(NCLS),CLS_E(NCLS)
*. Output 
      DIMENSION EBASC(*),CBASC(*)
*
      ZERO = 0.0D0
      CALL SETVEC(EBASC,ZERO,NSPC)
      CALL SETVEC(CBASC,ZERO,NSPC)
*
      DO ICLS = 1, NCLS
        ISPC = IBASSPC(ICLS)
        EBASC(ISPC) = EBASC(ISPC) + CLS_E(ICLS)
        CBASC(ISPC) = CBASC(ISPC) + CLS_C(ICLS)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =============================================='
        WRITE(6,*) ' Contribution  to energy and wf per base space '
        WRITE(6,*) ' =============================================='
        WRITE(6,*)
        WRITE(6,'(A)') '  Class         Energy          wf '
        WRITE(6,'(A)') ' ==========================================='
        DO ISPC = 1, NSPC
          WRITE(6,'(2X,I3,3X,E12.6,2X,E12.6)') 
     &          ISPC,EBASC(ISPC),CBASC(ISPC)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE CKAJJB_PA(CKAJJB,CKAJJBPA,IWAY,NKA,NJ,NJB,
     &                     IREO,NJB_PAS,NJB_ACT,NSMST,SIGN,
     &                     I_ADD_COPY)
*
* Reform C(Ka,j,Jb) between usual and active/passive division of Jb :
*
* C(Ka,j,Jb) <=> C(Ka,Jb_pa,j,Jb_ac)
*
* IWAY : =1 => normal to passive/active form
* IWAY : =2 => passive/active to normal form

* Jeppe Olsen, March 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKAJJB(*) 
      DIMENSION NJB_PAS(*),NJB_ACT(*),IREO(*)
*. Output
      DIMENSION CKAJJBPA(*)
CT    CALL LUCIAQENTER('CKA_PA')
*
      NTEST = 000
*
      IJB = 0
      DO JB_ACT_SM = 1, NSMST
*
        IF(JB_ACT_SM.EQ.1) THEN
         IOFF = 1
         IOFF_PA = 1   
        ELSE
         IOFF = IOFF + NJB_PAS(JB_ACT_SM-1)*NJB_ACT(JB_ACT_SM-1)
         IOFF_PA = IOFF_PA 
     &           + NJB_ACT(JB_ACT_SM-1)*NJ*NJB_PAS(JB_ACT_SM-1)*NKA
        END IF
*
        LJB_ACT = NJB_ACT(JB_ACT_SM)
        LJB_PAS = NJB_PAS(JB_ACT_SM)
        IB_OUT_ADD = LJB_PAS*NKA
* C(Ka,j,Jb) <=> C(Ka,Jb_pa,j,Jb_ac)
        DO JB_ACT = 1, LJB_ACT                
          DO JB_PAS = 1, LJB_PAS
            IJB =  IJB + 1
            IJB_IN = IREO(IJB)
            IB_IN = (IJB_IN-1)*NJ*NKA-NKA
            IB_OUT = (JB_ACT-1)*NJ*LJB_PAS*NKA
     &             + (JB_PAS-1)           *NKA + IOFF_PA - 1
     &             - IB_OUT_ADD
            DO J = 1, NJ
*. Address of C(1,j,Jb)
C             IB_IN = (IJB_IN-1)*NJ*NKA+(J-1)*NKA
              IB_IN = IB_IN + NKA
*. Address of C(1,Jb_pa,j,Jb_ac)
              IB_OUT = IB_OUT + IB_OUT_ADD
C             IB_OUT = (JB_ACT-1)*NJ*LJB_PAS*NKA
C    &               + (J     -1)   *LJB_PAS*NKA
C    &               + (JB_PAS-1)           *NKA + IOFF_PA - 1
C             WRITE(6,*) 'IB_IN, IB_OUT', IB_IN, IB_OUT
C             WRITE(6,*) ' IJB_IN, NJA,NKA ',  IJB_IN, NJA,NKA 
              IF(IWAY.EQ.1) THEN
*. Normal => passive/active form
                DO KA = 1, NKA
                  CKAJJBPA(IB_OUT+KA) = CKAJJB(IB_IN+KA)
                END DO
              ELSE IF(IWAY.EQ.2) THEN
C               IF(I_ADD_COPY.EQ.1) THEN
*. Passive/Active => Normal form
                  DO KA = 1, NKA
                    CKAJJB(IB_IN+KA) =  CKAJJB(IB_IN+KA) +
     &              SIGN*CKAJJBPA(IB_OUT+KA)
                  END DO
C               ELSE IF (I_ADD_COPY.EQ.2) THEN
C                 DO KA = 1, NKA
C                   CKAJJB(IB_IN+KA) =  
C    &              SIGN*CKAJJBPA(IB_OUT+KA)
C                 END DO
C               END IF
*               ^ End of ADD_COPY switch
              END IF 
*             ^ End of switch passive/active <=  => normal
            END DO
*           ^ End of loop over J-orbitals
          END DO
*         ^ End of loop over JB_pa
        END DO
*       ^ End of loop over JB_ac
      END DO
*     ^ End of loop over symmetry of active beta strings
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' =================='
       WRITE(6,*) ' CKAJJB_PA speaking '
       WRITE(6,*) ' =================='
       WRITE(6,*)
       IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' Normal to passive/active form '
       ELSE
         WRITE(6,*) ' passive/active  to normal form '
       END IF
*
       WRITE(6,*) ' C(ka,j,Jb) as C(kaj,Jb) '
       NKAJ = NKA*NJ
       CALL WRTMAT(CKAJJB,NKAJ,NJB,NKAJ,NJB)
*
       WRITE(6,*) ' C(Ka,Jb_pa,j,Jb_ac) as blocks C(KaJB_paj,Jb_ac)'
       DO JB_ACT_SM = 1, NSMST
         WRITE(6,*) ' Symmetry of Jb_act ', JB_ACT_SM
         IF(JB_ACT_SM.EQ.1) THEN
          IOFF_PA = 1   
         ELSE
          IOFF_PA = IOFF_PA 
     &            + NJB_ACT(JB_ACT_SM-1)*NJ*NJB_PAS(JB_ACT_SM-1)*NKA
         END IF
         NROW = NKA*NJB_PAS(JB_ACT_SM)*NJ
         NCOL = NJB_ACT(JB_ACT_SM)
*
C        WRITE(6,*)  ' NJB_ACT(JB_ACT_SM)', NJB_ACT(JB_ACT_SM)
C        WRITE(6,*)  ' NJB_PAS(JB_ACT_SM)', NJB_PAS(JB_ACT_SM)
C        WRITE(6,*)  ' NJ and NKA ', NJ,NKA
C        WRITE(6,*) ' Offset IOFF_PA', IOFF_PA
         CALL WRTMAT(CKAJJBPA(IOFF_PA),NROW,NCOL,NROW,NCOL)
       END DO
      END IF
*
CT    CALL LUCIAQEXIT('CKA_PA')
      RETURN
      END
      SUBROUTINE NXTNUM2_REV(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL,MAXVAL(I). 
*
* In this version rightmost digits are increased first
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen March 1998
*
*. Input
      DIMENSION MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
       NONEW = 1
       GOTO 1001
      END IF
*
      IPLACE = NELMNT + 1
 1000 CONTINUE
        IPLACE = IPLACE - 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.GT.1 ) THEN     
          INUM(IPLACE) = 1               
        ELSE IF ( IPLACE. EQ. 1 ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
      SUBROUTINE NST_SPGRP2(NIGRP,IGRP,ISM_TOT,NSMST,NSTRIN,NDIST)
*
* Number of strings for given combination of groups and 
* symmetry.
*
* Differece to NST_SPGRP : Number of strings per symmetry and group is
* not input paramter. Array NSTFSMGP from /GASSTR/ is used
*
*. Input
*        
*
*   NGRP : Number of active groups 
*   IGRP : The active groups
*   ISM_TOT : Total symmetry of supergroup
*   NSMST   : Number of string symmetries
*
*. Output
*
*  NSTRIN : Number of strings with symmetry ISM_TOT
*  NDIST  : Number of symmetry distributions
*
* Jeppe Olsen, September 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
*. Specific Input
      DIMENSION IGRP(NGRP)
*. General input
      INCLUDE 'gasstr.inc'
C     DIMENSION NSTSGP(NSMST,*)
*. Scratch 
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' NST_SPGP is speaking '
        WRITE(6,*) ' ====================='
*
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
      END IF
*
      IF(NIGRP.EQ.0) THEN
        IF(ISM_TOT.EQ.1) THEN 
          LENGTH = 1
        ELSE
          LENGTH = 0
        END IF
        GOTO 1001
      END IF
*
*. Set up min and max values for symmetries
      CALL MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNSM,MXSM,NDISTX)
*. Loop over symmetry distributions
      IFIRST = 1
      LENGTH = 0 
      NDIST = 0
*. Last group with symmetry differing from total symmetric
      NIGRPL = 1
      DO JGRP = 1, NIGRP
        IF(MXSM(JGRP).GT.1) NIGRPL = JGRP
      END DO
*. Number of strings in groups after NIGRPL
      NSTRL = 1
      DO JGRP = NIGRPL+1,NIGRP
        NSTRL = NSTRL*NSTFSMGP(1,IGRP(JGRP))
      END DO
 1000 CONTINUE
*. Next symmetry distribution
        CALL NEXT_SYM_DISTR(NIGRPL,MNSM,MXSM,ISM,ISM_TOT,IFIRST,NONEW)
        IF(NONEW.EQ.0) THEN
          LDIST = NSTRL
          DO JGRP = 1, NIGRPL
            LDIST = LDIST*NSTFSMGP(ISM(JGRP),IGRP(JGRP))
          END DO
          LENGTH = LENGTH + LDIST
          NDIST = NDIST + 1
      GOTO 1000
        END IF
*
 1001 CONTINUE
      NSTRIN = LENGTH
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings obtained ', LENGTH
        WRITE(6,*) ' Number of symmetry-distributions',NDIST
      END IF
*
      RETURN
      END
      SUBROUTINE ACT_GRPS(ISPGRP,NIGRP,IOP,NOP,
     &                    NACGRP,IACGRP,NPAGRP,IPAGRP)
*
* A  supergroups (ISPGRP) 	
* and a set of NOP creation/annihilation operators working in
* orbital spaces IOP are given.
*
* Divide supergroup into active and  passive parts.
*
* Jeppe Olsen, March 98
*
* Version of March 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
*. Specific input
      INTEGER ISPGRP(NIGRP),IOP(*)
*. Output
      INTEGER IACGRP(*),IPAGRP(*)
*. Local scratch
      INTEGER IACT(MXPNGAS)
*. Output

*
      IZERO = 0  
      CALL ISETVC(IACT,IZERO,NGAS)
*. Active orbital spaces
      DO JOP = 1, NOP
        IIACT = IOP(JOP)
        IACT(IIACT) = IACT(IIACT) + 1
      END DO
*. Active groups in Left and Right supergroups :
*. Output are the actual active groups (not the location in supergroup)
      NACGRP = 0
      NPAGRP = 0
      DO IGRP = 1, NIGRP
        IF(IACT(IGSFGP(ISPGRP(IGRP))).NE.0) THEN
*. active
          NACGRP = NACGRP + 1
          IACGRP(NACGRP) = ISPGRP(IGRP)
        ELSE
*. Passive
          NPAGRP = NPAGRP + 1
          IPAGRP(NPAGRP) = ISPGRP(IGRP)
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)'Division of supergroup into active and passive parts'
        WRITE(6,*)'===================================================='
        WRITE(6,*)
        WRITE(6,*) ' Input supergroup and operators'
        WRITE(6,*)
        CALL IWRTMA(ISPGRP,1,NIGRP,1,NIGRP)
        CALL IWRTMA(IOP,1,NOP,1,NOP)
        WRITE(6,*)
        WRITE(6,*) ' Active  part ', (IACGRP(IAC),IAC=1,NACGRP)
        WRITE(6,*) ' Passive part ', (IPAGRP(IPA),IPA=1,NPAGRP)
      END IF
*
      RETURN 
      END
      SUBROUTINE REPART_NORD_MAT(N,L,NR,IR,NC,IC,IREO)
*
*. Repartion a NORD matrix into a normal two index matrix
*
* Input order C(I1,I2,I3....)
*
* Output order C(IR,IC), where IR is a compound index of the NR indeces in IR
*                              IC is a compound index of the NC indeces in IC
*
* Jeppe Olsen, Sept. 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER L(N),IR(NR),IC(NC)
*. Output
      DIMENSION IREO(*)
*. Local Scratch
      PARAMETER(MXPNORD = 20)
      DIMENSION INUM(MXPNORD)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '***************************'
        WRITE(6,*)  'INFO FROM REPART_NORD_MAT'
        WRITE(6,*) '***************************'
*
        WRITE(6,*) ' Number of index in new row index ', NR
        WRITE(6,*) ' indeces defining row : '
        CALL IWRTMA(IR,1,NR,1,NR)
        WRITE(6,*)
        WRITE(6,*) ' Number of index in new column index ', NC
        WRITE(6,*) ' indeces defining column: '
        CALL IWRTMA(IC,1,NC,1,NC)
*
        WRITE(6,*) ' Dimension of each index '
        CALL IWRTMA(L,1,N,1,N)
      END IF
*
      IF(N.GT.MXPNORD) THEN
        WRITE(6,*) ' Error in REPART_NORD_MAT '
        WRITE(6,*) 
     &  ' Actual order (N) Larger than max order (MXPNORD)'
        WRITE(6,*) ' N, MXPNORD ', N,MXPNORD
        STOP' ERROR in REPART_NORD_MAT'
      END IF
*. Total dimensions of new row and column index
      LR = 1
      DO JR = 1, NR
        LR = LR*L(IR(JR))
      END DO
      LC = 1
      DO JC = 1, NC
        LC = LC*L(IC(JC))
      END DO
*. Loop over elements in input
      NELMNT = 1
      DO I = 1, N
        NELMNT = NELMNT*L(I)
      END DO
      IF(NTEST.GE.100) WRITE(6,*) 'NELMNT = ', NELMNT
      IFIRST = 1
      DO IELMNT = 1, NELMNT
        IF(IFIRST.EQ.1) THEN
          DO IGAS = 1, N
            INUM(IGAS) = 1
          END DO
          IFIRST = 0
        ELSE
C         CALL NXTNUM2(INUM,NORD,1,LEN,NONEW)
          CALL NXTNUM2(INUM,N,1,L,NONEW)
        END IF
*. first index, new ordering
        JR = 1
        MULT = 1
        DO KR = 1, NR
          JR = JR + (INUM(IR(KR))-1)*MULT
          MULT = MULT*L(IR(KR))
        END  DO
*. Second index, new ordering 
        JC = 1
        MULT = 1
        DO KC = 1, NC
          JC = JC + (INUM(IC(KC))-1)*MULT
          MULT = MULT*L(IC(KC))
        END  DO
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Input element, original order '
          CALL IWRTMA(INUM,1,N,1,N)
          WRITE(6,*) ' Corresponding row and column indexes ', 
     &    JR,JC
        END IF
*. row, column order => original order
        IREO(JR+(JC-1)*LR) = IELMNT
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' reordering matrix from REPART.... '
        CALL IWRTMA(IREO,LR,LC,LR,LC)
      END IF
*
      RETURN
      END
      SUBROUTINE MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDIST)
*
* A combination of NIGRP groups are given (IGRP)
*. Find MIN and MAX for symmetry in each group
*
* Jeppe Olsen, September 1997
*              April 1998     From  MINMAX_SM_GP
*
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. Include blocks     
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'strbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'csm.inc'
*. Input
      DIMENSION IGRP(NIGRP)
*.Output
      DIMENSION MNVAL(NIGRP),MXVAL(NIGRP)
*. Local scratch
      DIMENSION LSMGP(MXPOBS,MXPNGAS)
*
      NTEST = 0000
      IF(NTEST.GE.100) WRITE(6,*) ' >> Entering MINMAX_... <<'
*
      DO JGRP = 1, NIGRP
        MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
        MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
      END DO
        
*. Number of strings per sym and group
C     DO JGRP = 1, NIGRP
C       CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(JGRP)-1)*NSMST+1,
C    &               NSMST,LSMGP(1,JGRP))
C     END DO
C     IF(NTEST.GE.1000) THEN
C       WRITE(6,*) ' LSMGP '
C       CALL IWRTMA(LSMGP,NSMST,NIGRP,MXPOBS,NIGRP)
C     END IF
C. Max and min sym in each group
C     DO JGRP = 1, NIGRP
*
C       IMAX = 1
C       DO ISM=1, NSMST
C         IF(LSMGP(ISM,JGRP).GT.0) IMAX = ISM
C       END DO
C       MXVAL(JGRP) = IMAX
*
C       IMIN = NSMST
C       DO ISM = NSMST,1,-1
C        IF(LSMGP(ISM,JGRP).GT.0) IMIN = ISM
C       END DO
C       MNVAL(JGRP) = IMIN
C     END DO
*. Total number of symmetry distributions
      NDIST = 1
      DO JGRP = 1, NIGRP
        NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Group combination : '
        WRITE(6,'(5X,10I3)') (IGRP(JGRP),JGRP=1, NIGRP)
        WRITE(6,*)
        WRITE(6,*) ' Group Minsym Maxsym'
        WRITE(6,*) ' ==================='
        DO JGRP = 1, NIGRP
          WRITE(6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Total number of distributions', NDIST
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving MINMAX_... <<'
*
      RETURN
      END
      SUBROUTINE NEXT_SYM_DISTR(NGAS,MINVAL,MAXVAL,
     &           ISYM,ISYM_TOT,IFIRST,NONEW)
*
* Obtain next distribution of symmetries with given total
* Symmetry. 
*
* Loop over Gaspace2 2-  NGAS spaces are performed, and the symmetry
* of the first space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
* Last modification; May 5, 2013; Jeppe Olsen, switched determined 
*                    space from NGASL to 1
*
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION MINVAL(NGAS),MAXVAL(NGAS)
*. Input and output
      DIMENSION ISYM(NGAS)
*
      INCLUDE 'multd2h.inc'
*
*. Symmetry of space 2- NGAS  spaces
*
      IF(IFIRST.EQ.1) THEN
        DO IGAS = 2, NGAS
          ISYM(IGAS) = MINVAL(IGAS)
        END DO
        NONEW = 0
      END IF
 1001 CONTINUE
      IF(IFIRST.EQ.0) CALL NXTNUM3
     &                (ISYM(2),NGAS-1,MINVAL(2),MAXVAL(2),NONEW)
      IFIRST = 0
*
*. Symmetry of first space
*
      IF(NONEW.EQ.0) THEN
        JSYM = ISYMSTR(ISYM(2),NGAS-1)
        ISYM(1) = MULTD2H(JSYM,ISYM_TOT)
*
        IF(MINVAL(1).GT.ISYM(1).OR.
     &     MAXVAL(1).LT.ISYM(1)    )GOTO 1001
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from NEXT_SYM_DISTR:'
        IF(NONEW.EQ.1) THEN
         WRITE(6,*) ' No new symmetry distributions '
        ELSE
         WRITE(6,*) ' Next symmetry distribution '
         CALL IWRTMA(ISYM,1,NGAS,1,NGAS)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE REO_NORD_MAT(NORD,LEN,IREO,IREO_ARRAY)
*
* A multiorder matrix XIN(I1,I2,I3,... INORD) is considered
* Reorder the indeces according to IREO( OLD order => NEW order)
* Only reorder array is obtained, no reordering of matrix is
* performed
* The dimension of index I is given by LEN(I)
*
* NOTE THE REORDERING ARRAY IREO_ARRAY gives new order => Old order
*
* The Addressing of the multi ordered matrices are
* done with the usual qc ordering :
* Address of (XIN(I1,I2,.....IN) =I1 + (I2-1)*LEN(1)+(I3-1)*LEN(2)*LEN(1) ..
*
* I.e. left most index are the fastest varying
*
* Output is the reorder array IREO_ARRAY
*
* Jeppe Olsen, Sept 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input 
      INTEGER LEN(NORD), IREO(NORD)
*. Output
      DIMENSION IREO_ARRAY(*)
*. Local Scratch
      PARAMETER(MXPNORD=20)
      DIMENSION INUM(MXPNORD),IREONUM(MXPNORD),IREOLEN(MXPNORD)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' REO_NORD_MAT in action '
        WRITE(6,*) ' ======================='   
        WRITE(6,*)  
        WRITE(6,*) ' Input to output order of indeces '
        CALL IWRTMA(IREO,1,NORD,1,NORD)
      END IF
*
      IF(NORD.GE.MXPNORD) THEN
        WRITE(6,*) '  FATAL problem in REO_NORD_MAT'
        WRITE(6,*) 
     &  '  Actual order = ',NORD,'  larger then MAX= ',MXPNORD
        WRITE(6,*) ' Increase MXPNORD in REO_NORD_MAT'
        STOP       ' Increase MXPNORD in REO_NORD_MAT'
      END IF
*. Length of reordered indeces
      DO IELM = 1, NORD
        IREOLEN(IREO(IELM)) = LEN(IELM)
      END DO
*. Number of elements  
      NELMNT = 1
      DO IELM = 1, NORD
        NELMNT = NELMNT*LEN(IELM)
      END DO
*. Ordered loop over elements in input list
      DO IELMNT = 1, NELMNT
        IF(IELMNT.EQ.1) THEN
*. Initialize
          IONE = 1
          CALL ISETVC(INUM,IONE,NORD)
        ELSE
*. Next element
C         NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
          CALL NXTNUM2(INUM,NORD,1,LEN,NONEW)
        END IF
*. To new order
        DO IELM = 1, NORD
          IREONUM(IREO(IELM)) = INUM(IELM)
        END DO
*. Address in reordered matrix
        IMULT = 1
        IADR = 1
        DO IELM = 1, NORD
          IADR = IADR + (IREONUM(IELM)-1)*IMULT
          IMULT = IMULT*IREOLEN(IELM)
        END DO
*.=========================
*. NEW ORDER => OLD ORDER !
*.=========================
        IREO_ARRAY(IADR) = IELMNT
*
        IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Info for element number  ', IELMNT
         WRITE(6,*) ' indeces in input order '
         CALL IWRTMA(INUM,1,NORD,1,NORD)
         WRITE(6,*) ' indeces in output order '
         CALL IWRTMA(IREONUM,1,NORD,1,NORD)
         WRITE(6,*) ' Output address ', IADR
        END IF
*
      END DO
*
      IF(NTEST.GE.10) THEN 
         WRITE(6,*) ' Reorder array from REO_NORD_MAT'             
         CALL IWRTMA(IREO_ARRAY,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END  
*
      FUNCTION IOFF_SYM_DIST(ISYM,NGASL,IOFF,MAXVAL,MINVAL)
*
* A ts block of string is given and the individual 
* symmetrydisrtributions has been obtained ( for example 
* by TS_SYM_PNT)
*
* Obtain offset for symmetrycombination defined by ISYM
*
*. Jeppe Olsen
*. Last modification: May 5, 2013; Jeppe Olsen; implied index is first
*                     instrad of last
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ISYM(*),IOFF(*),MAXVAL(*),MINVAL(*)
* Address in IOFF is 
*     1
*     +  (ISM2-MINVAL(2))
*     +  (ISM3-MINVAL(3))*(MAXVAL(2)-MINVAL(2)+1)
*     +  (ISM4-MINVAL(4))*(MAXVAL(3)-MINVAL(3)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-MINVAL(L))*Prod(i=2,L-1)(MAXVAL(i)-MINVAL(i)+1)           
* 
C     write(6,*) ' Isym and minval '
C     call iwrtma(isym,1,ngasl,1,ngasl)
C     call iwrtma(minval,1,ngasl,1,ngasl)
      I = 1
      IMULT = 1
      DO IGAS = 2, NGASL
        I = I + (ISYM(IGAS)-MINVAL(IGAS)) * IMULT
        IMULT = IMULT*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
C       write(6,*) ' igas i imult ',igas,i,imult
      END DO
      IOFF_SYM_DIST=IOFF(I)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from IOFF_SYM_DIST'
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Address and offset ',I,IOFF_SYM_DIST
        WRITE(6,*) ' Symmetry distribution : ', (ISYM(J),J=1,NGASL)
      END IF
*
      RETURN
      END
      SUBROUTINE TS_SYM_PNT2(IGRP,NIGRP,MAXVAL,MINVAL,ISYM,
     &           IPNT,LPNT)
*
* Construct pointers to start of symmetrydistributions 
* for supergroup of strings with given symmetry ISYM.
* 
* The symmetries are labeled by symmetry of GAS2- GAS L,
* and the symmetry of the first GASpace is then defined to
* give the right overall symmetry
*
* The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
* is thus given as
*     1
*     +  (ISM2-MINVAL(2))
*     +  (ISM3-MINVAL(3))*(MAXVAL(2)-MINVAL(2)+1)
*     +  (ISM4-MINVAL(4))*(MAXVAL(3)-MINVAL(3)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-MINVAL(L))*Prod(i=2,L-1)(MAXVAL(i)-MINVAL(i)+1)           
*
* Where L is the last group of strings with nonvanishing occupation
*
* Jeppe Olsen, September 1997, Shaved in Jan. 2011
*. Last modification: May 5, 2013; Jeppe Olsen; Free group changed to NGASL
*
* Version 2 : Uses IGRP and NIGRP to define supergroup
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
*. Specific Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
*
      INTEGER IGRP_AC(MXPNGAS)
      INTEGER MINVAL_AC(MXPNGAS),MAXVAL_AC(MXPNGAS)
*. Output
      INTEGER MINVAL(*),MAXVAL(*),IPNT(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TS_SYM_PNT2'
        WRITE(6,*) ' ======================='
        WRITE(6,*)
      END IF
*. Info on groups of strings in supergroup
      NGASL = 1
      DO IGAS = 1, NIGRP
        IF(NELFGP(IGRP(IGAS)).GT.0) NGASL = IGAS
*. The following two lines have actually been moved out in ADAST2
        MINVAL(IGAS) = MINMAX_SM_GP(1,IGRP(IGAS))
        MAXVAL(IGAS) = MINMAX_SM_GP(2,IGRP(IGAS))
      END DO
*. Active groups : for consistency, the last active group should be NGASL 
      NGRP_AC = 0
      NSTR_PS = 1
      ISYM_PS = 1
      DO JGRP = 1, NIGRP
        IF(MINVAL(JGRP).NE.MAXVAL(JGRP).OR.JGRP.EQ.1) THEN
          NGRP_AC = NGRP_AC + 1
          IGRP_AC(NGRP_AC) = IGRP(JGRP)
          MINVAL_AC(NGRP_AC) = MINVAL(JGRP)
          MAXVAL_AC(NGRP_AC) = MAXVAL(JGRP)
        ELSE
          ISYM_PS = MULTD2H(ISYM_PS,MINVAL(JGRP))
          NSTR_PS = NSTR_PS*NSTFSMGP(MINVAL(JGRP),IGRP(JGRP))
        END IF
      END DO
*. Total number of symmetry blocks that will be generated
      NBLKS = 1
      DO IGAS = 2, NGRP_AC
       NBLKS = NBLKS*(MAXVAL_AC(IGAS)-MINVAL_AC(IGAS)+1)
      END DO
      IF(NBLKS.GT.LPNT) THEN
        WRITE(6,*) ' Problem in TS_SYM_PNT'
        WRITE(6,*) ' Dimension of IPNT too small'
        WRITE(6,*) ' Actual and required length',NBLKS,LPNT
        WRITE(6,*)
        WRITE(6,*) ' I will Stop and wait for instructions'
        STOP' TS_SYM_PNT too small'
      END IF
*. Loop over symmetry blocks in standard order
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
      JOFF = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 2, NGRP_AC 
            ISMFGS(IGAS) = MINVAL_AC(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in 2-NGAS 
C        CALL NXTNUM3(ISMFGS,NGRP_AC-1,MINVAL_AC,MAXVAL_AC,NONEW)
         IF(NGRP_AC-1.EQ.0) THEN
           NONEW = 1
           GOTO 1001
         END IF
         IPLACE = 1
 1000    CONTINUE
           IPLACE = IPLACE + 1
           IF(ISMFGS(IPLACE).LT.MAXVAL_AC(IPLACE)) THEN
             ISMFGS(IPLACE) = ISMFGS(IPLACE) + 1
             NONEW = 0
             GOTO 1001
           ELSE IF ( IPLACE.LT.NGRP_AC) THEN
             DO JPLACE = 1, IPLACE
               ISMFGS(JPLACE) = MINVAL_AC(JPLACE)
             END DO
           ELSE IF ( IPLACE. EQ. NGRP_AC ) THEN
             NONEW = 1 
             GOTO 1001
           END IF
         GOTO 1000
 1001    CONTINUE

         IF(NONEW.NE.0) GOTO 2001
        END IF
*
        IFIRST = 0
*. Symmetry of spaces 2- NGASL spaces given, symmetry of full space
        IF(NGRP_AC.GT.1) THEN
          ISTSMM1 = 1
        ELSE 
          ISTSMM1 = 1
        END IF
        DO JGAS = 2, NGRP_AC
          ISTSMM1 = MULTD2H(ISTSMM1,ISMFGS(JGAS))
        END DO
        ISTSMM1 = MULTD2H(ISTSMM1,ISYM_PS)
*.  sym of SPACE 1
        ISMGS1 = MULTD2H(ISTSMM1,ISYM)
        ISMFGS(1) = ISMGS1
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Next symmetry distribution: '
        END IF
          
*. Number of strings with this symmetry combination
        NSTRII = NSTR_PS
        DO IGAS = 1, NGRP_AC
          NSTRII = NSTRII*NSTFSMGP(ISMFGS(IGAS),IGRP_AC(IGAS))
        END DO
        JOFF = JOFF + 1
        IPNT(JOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TS_SYM_PNT'
        WRITE(6,*) ' Required total symmetry',ISYM
        WRITE(6,*) ' Number of symmetry blocks ', NBLKS
        WRITE(6,*) 
        WRITE(6,*) ' Offset array  for symmetry blocks'
        CALL IWRTMA(IPNT,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END 
      SUBROUTINE TS_SYM_PNT(ITP,ISYM,ISPGP,MAXVAL,MINVAL,IPNT,LPNT)
*
* Construct pointers to start of symmetrydistributions 
* for supergroup of strings with given symmetry
*
* The start of symmetry block ISYM1 ISYM2 ISYM3 .... ISYMN
* is given as
*     1
*     +  (ISM1-MINVAL(1))
*     +  (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     +  (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)           
*
* Where L is the last group of strings with nonvanishing occupation
*
* Jeppe Olsen, December 1996
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
C     INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
C     INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*. Output
      INTEGER MINVAL(*),MAXVAL(*),IPNT(*)
*
      NTEST = 000
*. Absolute supergroup numbers 
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      IF(NTEST.GE.100) WRITE(6,*) ' Supergroup in action ',ISPGPABS
      
*. Info on groups of strings in supergroup
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
*. Number of strings per symmetry in each gasspace  
        CALL ICOPVE2(int_mb(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
      END DO
C?    WRITE(6,*) ' NNSTSGP'
C?    CALL IWRTMA(NNSTSGP,NSMST,NGAS,MXPNSMST,MXPNGAS)
C     NNSTSGP(MXPNSMST,MXPNGAS)
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' MINVAL and MAXVAL '
        CALL IWRTMA(MINVAL,1,NGAS,1,NGAS)
        CALL IWRTMA(MAXVAL,1,NGAS,1,NGAS)
        WRITE(6,*) ' NGAS = ', NGAS
      END IF
 
*. Total number of strings that will be generated
      NBLKS = 1
      DO IGAS = 1, NGASL-1
       NBLKS = NBLKS*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
      END DO
      IF(NBLKS.GT.LPNT) THEN
        WRITE(6,*) ' Problem in TS_SYM_PNT'
        WRITE(6,*) ' Dimension of IPNT too small'
        WRITE(6,*) ' Actual and required length',NBLKS,LPNT
        WRITE(6,*)
        WRITE(6,*) ' I will Stop and wait for instructions'
        STOP' TS_SYM_PNT too small'
      END IF
*. Loop over symmetry blocks in standard order
      IFIRST = 1
      ISTRBS = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISYM)
        ISMFGS(NGASL) = ISMGSN
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL-1
          IOFF = IOFF + (ISMFGS(IGAS)-MINVAL(IGAS))*MULT
          MULT = MULT * (MAXVAL(IGAS)-MINVAL(IGAS)+1)
        END DO
*
        IPNT(IOFF) = NSTRINT + 1
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOFF, IPNT(IOFF) NSTRII ',
     &                 IOFF, IPNT(IOFF),NSTRII
        END IF
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from TS_SYM_PNT'
        WRITE(6,*) ' Number of symmetry blocks ', NBLKS
        WRITE(6,*) 
        WRITE(6,*) ' Offset array  for symmetry blocks'
        CALL IWRTMA(IPNT,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END 
      FUNCTION ISYMSTR(ISYM,NSTR)
*
* Symmetry of product of NSTR string symmetries
*
* works currently only for D2H and subgroups
*
* Jeppe Olsen, 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER ISYM(*)
*
      IF(NSTR.EQ.0) THEN
        IISYM = 1
      ELSE
        IISYM = ISYM(1)
        DO JSTR = 2, NSTR
           IISYM = MULTD2H(IISYM,ISYM(JSTR))
        END DO
      END IF
*
      ISYMSTR = IISYM
*
      RETURN
      END
      SUBROUTINE GET_DIAG_DOM(AIN,AOUT,NDIM,IREO)
*
* A square matrix AIN is given. Permute columns to obtain
* Largest elements on diagonal
*
* Jeppe Olsen, Feb.98
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION AIN(NDIM,NDIM)
*.Output 
      DIMENSION AOUT(NDIM,NDIM)
      INTEGER IREO(*)
*
*. IREO : Old => New
      IZERO = 0
      CALL ISETVC(IREO,IZERO,NDIM)
      DO ICOL = 1, NDIM
*. Column with largest values at row ICOL
        XVAL = 0.0D0               
        IICOL = 0
        DO JCOL = 1, NDIM
          IF(ABS(AIN(ICOL,JCOL)).GT.XVAL.AND.IREO(JCOL).EQ.0) THEN
            IICOL = JCOL
            XVAL = ABS(AIN(ICOL,JCOL))
          END IF
        END DO
        IF(IICOL.GT.0) THEN
          CALL COPVEC(AIN(1,IICOL),AOUT(1,ICOL),NDIM)
          IREO(IICOL) = ICOL
        ELSE 
          CALL COPVEC(AIN(1,ICOL),AOUT(1,ICOL),NDIM)
          IREO(ICOL) = ICOL
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' OUTPUT FROM GET_DIAG_DOM '
        WRITE(6,*) ' ======================== '
        WRITE(6,*)
        WRITE(6,*) ' IREO array : old => new order '
        CALL IWRTMA(IREO,1,NDIM,1,NDIM)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reordered matrix '
        CALL WRTMAT(AOUT,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE TRACI_CTL
*
* Master routine for transforming CI vectors to new orbital basis
*
* Input vectors on LUC
* output vectors delivered on LUHC
* (LUSC1,LUSC2,LUSC3) used as scratch files
*
* Jeppe Olsen, January 98
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'  
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
*
      INCLUDE 'cands.inc'
*
      NTEST = 0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACIC')
*
      WRITE(6,*)
      WRITE(6,*) ' Wellcome to TRACI_CTL'
      WRITE(6,*) ' ====================='
      WRITE(6,*)
*. The three scratch  blocks
C          GET_3BLKS(KVEC1,KVEC2,KC2)
      CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)

*. MO-MO transformation matrix :
      CALL MEMMAN(KLCMOMO,NTOOB**2,'ADDL  ',2,'CMOMO ')
*. Copy of one-electron integrals
      CALL MEMMAN(KLH1SAVE,NTOOB**2,'ADDL  ',2,'H1SAVE') !done
*. We are going to mess with the one-electron integrals, take a copy
      CALL COPVEC(WORK(KINT1),dbl_mb(KLH1SAVE),NTOOB*NTOOB)
*. Set up block structure of CI space
      IATP = 1
      IBTP = 2
      CALL  Z_BLKFO(ISSPC,ISSM,IATP,IBTP,KLCLBT,KLCLEBT,
     &      KLCI1BT,KLCIBT,KLCBLTP,NBATCH,NBLOCK)
C           Z_BLKFO(ISPC,ISM,IATP,IBTP,KPCLBT,KPCLEBT,
C    &              KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK)
*
      LBLK = -1
      CALL REWINO(LUC)
*. Make sure LUDIA corresponds to original def
      LUDIA = 20
      CALL REWINO(LUDIA)
      DO JROOT = 1, NROOT
*. One-electron density for root JROOT
        CALL REWINO(LUSC1)
        CALL COPVCD(LUC,LUSC1,WORK(KVEC1),0,LBLK)
        CALL COPVCD(LUSC1,LUSC2,WORK(KVEC1),1,LBLK)
        XDUM = 0.0D0
        CALL DENSI2(1,dbl_mb(KRHO1),dbl_mb(KRHO2),
     &        WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC2,EXPS2,
     &        0,XDUM,XDUM,XDUM,XDUM,1)
*. Obtain MO-MO transformation matrix 
        CALL MOROT(IFINMO)
        WRITE(6,*) ' Memcheck after MOROT'
        CALL MEMCHK
*
*. Transform CI vector : Input on LUHC, output on LUDIA (!)
        CALL COPVCD(LUSC1,LUHC,WORK(KVEC1),1,LBLK)
*
        CALL TRACI(WORK(KMOMO),LUHC,LUDIA,ISSPC,ISSM,
     &             WORK(KVEC1),WORK(KVEC2))
C            TRACI(X,LUCIN,LUCOUT,IXSPC,IXSM,VEC1,VEC2)
       WRITE(6,*)
       WRITE(6,*) ' Analysis of rotated state number ', JROOT
       WRITE(6,*) ' ====================================='
       WRITE(6,*)
       CALL GASANA(WORK(KVEC1),NBLOCK,int_mb(KLCIBT),dbl_mb(KLCBLTP),
     &                LUSC1,ICISTR)
      END DO
*     ^ End of loop over roots
      CALL REWINO(LUDIA)
      IF(NTEST.GE.100) THEN
        DO JROOT = 1, NROOT
          CALL WRTVCD(WORK(KVEC1),LUDIA,0,LBLK)
        END DO
      END IF
*
*. clean up time : copy 1-e integrals back in place                  
      CALL COPVEC(dbl_mb(KLH1SAVE),WORK(KINT1),NTOOB*NTOOB)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACIC')
*
      RETURN
      END 
      SUBROUTINE T_ROW_TO_H(T,H,K,TKK)
*
* Set H integrals 
*
*    Column K : H(P,K) = T(P,K)/T(K,K), P.NE.K
*    Other Columns     = 0
* - and return T_{kk} in TKK
* 
*
* Jeppe Olsen, Jan 98
* For rotation of CI vectors
*
c      IMPLICIT REAL*8 (A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input ( in blocked form)
      DIMENSION T(*)
*. Output ( also in blocked form)
      DIMENSION H(*)
*
      KSM = ISMFSO(K)
      KOFF = IBSO(KSM)
      KREL = K - KOFF + 1
      NK = NTOOBS(KSM)

*
      ZERO = 0.0D0
      CALL SETVEC(H,ZERO,NTOOB**2)
*
c     IOFF = IFRMR(WORK(KPGINT1A(1)),1,KSM)
      IOFF = dbl_mb(KPGINT1A(1) + KSM - 1)
      CALL COPVEC(T(IOFF+(KREL-1)*NK),H(IOFF+(KREL-1)*NK),NK)
      TKK = H(IOFF-1+(KREL-1)*NK+KREL)
      IF(TKK .NE. 0.0D0) THEN
        FAC = 1.0D0/TKK                   
        CALL SCALVE(H(IOFF+(KREL-1)*NK),FAC,NK)
C       H(IOFF-1+(K-1)*NK+K) = H(IOFF-1+(K-1)*NK+K) -1.0D0
        H(IOFF-1+(KREL-1)*NK+KREL) = 0.0D0                       
      ELSE
C       TKK = 1.0D0
        TKK = 0.0D0
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output from T_ROW_TO_H '
        WRITE(6,*) ' Row to be transferred ', KREL
        WRITE(6,*) ' TKK = ', TKK
        WRITE(6,*) ' Updated H matrix '
        WRITE(6,*) ' NTOOB, NTOOBS(1)', NTOOB, NTOOBS(1)
        CALL APRBLM2(H,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE EXP_BL_MAT(BLKA,A,IWAY,NBLOCK,LROW,LCOL)
*
* Switch betwen blocked matrix and complete matrix
*
* IWAY = 1 : Blocked to unblocked
*    else  : Unblocked to blocked
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER LROW(*),LCOL(*)
*. Input or output
      DIMENSION A(*),BLKA(*)
*. Dimensions of unblocked matrix
      NROW = IELSUM(LROW,NBLOCK)
      NCOL = IELSUM(LCOL,NBLOCK)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IROFF = 1
          ICOFF = 1
          IBLOFF = 1
        ELSE
          IROFF = IROFF + LROW(IBLOCK-1)
          ICOFF = ICOFF + LCOL(IBLOCK-1)
          IBLOFF = IBLOFF + LROW(IBLOCK-1)*LCOL(IBLOCK-1)
        END IF
        NROWB = LROW(IBLOCK)
        NCOLB = LCOL(IBLOCK)
        DO IC = 1, NCOLB         
         DO IR = 1, NROWB          
           IBL = IBLOFF-1 + (IC-1)*NROWB+IR
           ICM = (IC+ICOFF-1-1)*NROW + IR+IROFF-1
           IF(IWAY.EQ.1) THEN
             A(ICM) = BLKA(IBL)
           ELSE
             BLKA(IBL) = A(ICM)
           END IF
         END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output from EXP_BL_MAT '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' Blocked to complete matrix '
        ELSE
         WRITE(6,*) ' Complete to blocked matrix'
        END IF
        WRITE(6,*) ' Blocked matrix'
C            APRBLM2(A,LROW,LCOL,NBLK,ISYM)
        CALL APRBLM2(BLKA,LROW,LCOL,NBLOCK,0)
        WRITE(6,*)
        WRITE(6,*) ' Complete matrix '
        CALL WRTMAT(A,NROW,NCOL,NROW,NCOL)
      END IF
*
      RETURN
      END
      SUBROUTINE TRACI(X,LUCIN,LUCOUT,IXSPC,IXSM,VEC1,VEC2)
*
* A rotation matrix X is defining expansion from
* old to new orbitals
*        PHI(NEW) = PHI(OLD) * X
*
* change CI coefficients(sym IXSM, space IXSPC ) 
* so they corresponds to PHI(NEW) basis
*
* The input CI vector is on LUCIN and the transformed CI vector
* will be delivered on LUCOUT. 
*
* Transformation as conceived by Per-Aake Malmquist
* (I.J.Q.C. vol XXX, p479 ,1986 (OCTOBER ISSUE ))
*
*  Jeppe Olsen 1988
*
* New LUCIA version of Jan 1998
*
* note The transformation matrix X is supposed to be in complete form
* as a matrix over NTOOB orbitals.
*
c      IMPLICIT REAL*8 (A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
*. Common block for communicating with sigma
      INCLUDE 'cands.inc'
      INCLUDE 'mv7task.inc'
*
      DIMENSION X(*)
*
      CALL LUCIAQENTER('TRACI')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACI ')
*
      NTEST = 0000
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Welcome to TRACI '
        WRITE(6,*) ' ================'
        WRITE(6,*) 
        WRITE(6,*) ' IXSPC,IXSM = ', IXSPC,IXSM
      END IF
*
*. Memory allocation
* for a matrix T
      CALL MEMMAN(KLT,NTOOB**2,'ADDL  ',2,'TMAT  ') !done
cNW     CALL MEMMAN(KLTB,NTOOB**2,'ADDL  ',2,'TMATBL')
*. Scratch in PAMTMT
      LSCR = NTOOB**2 +NTOOB*(NTOOB+1)/2
      CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'KLSCR ') !done
*. Obtain T matrix used for transformation, for each symmetry separately
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN 
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISM-1)**2
        END IF
        IF(NTOOBS(ISM).GT.0) 
     &  CALL PAMTMT(X(IOFF),dbl_mb(KLT-1+IOFF),dbl_mb(KLSCR),
     &             NTOOBS(ISM))
      END DO
      LENT = IOFF + NTOOBS(NSMOB)**2 - 1
*. Save Malmqvist matrix
      CALL COPVEC(dbl_mb(KLT),dbl_mb(KTPAM),LENT)
*. Transform CI-vector
      ICSPC = IXSPC
      ICSM  = ICSM
      ISSPC = IXSPC
      ISSM  = IXSM
*
      I_STAND_OR_BLKDIA = 2
      IF(I_STAND_OR_BLKDIA.EQ.1.OR.NGAS.GT.1) THEN
*. No assumption about the expansion
        CALL TRACID(dbl_mb(KLT),LUCIN,LUCOUT,LUSC1,LUSC2,LUSC3,
     &              VEC1,VEC2)
      ELSE 
*. Special version for block diagonal expansions
        IF(NTEST.GE.10) WRITE(6,*) ' Block TRACI will be called '
        CMV7TASK = 'TRACID'
        CALL MV7(VEC1,VEC2,LUCIN,LUCOUT,XDUM,YDUM)
        CMV7TASK = 'SIGMA '
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Final Biotransformed CI vector from TRACI '
        CALL WRTVCD(VEC1,LUCOUT,1,-1)
      END IF
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACI ')
      CALL LUCIAQEXIT('TRACI')
*
      RETURN
      END
      SUBROUTINE TRACID(T,LUCIN,LUCOUT,LUSC1,LUSC2,LUSC3,
     &           VEC1,VEC2)
*
* Transform CI vector on LUCIN with T matrix after
* Docent Malmquist's recipe. Place result as next vector on LUOUT 
*
* The transformation is done as a sequence of one-electron transformations 
*
* with each orbital transformation being
*
* Sum(k=0,2) ( 1/k! sum(n'.ne.n) S(n'n) E_{n'n} ) Tnn^N_n
* 
* with Sn'n = T(n'n)/Tnn
*
* each transformation is 
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      REAL*8 INPROD
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cands.inc'
      REAL*8 INPRDD
*. Input
      DIMENSION T(*)
*. Scratch blocks ( two of them)
      DIMENSION VEC1(*),VEC2(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACID')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*)
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' Info from TRACID '
        WRITE(6,*) ' ================ '
        WRITE(6,*)
      END IF
*
      LBLK = -1
      IDUM = 1
*. Transfer vector on LUCIN to LUSC1
C           COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      CALL  COPVCD(LUCIN,LUSC1,VEC1,1,LBLK)
*. A bit of info for the sigma routine
      I_RES_AB = 0
*. Do the one-electron update
        I12 = 1
*. With 1-electron integrals in complete block form
        IH1FORM_SAVE = IH1FORM
        IH1FORM = 2
*. Transform each orbital separately
      DO K = 1, NTOOB
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Info for transformation of orbital ', K
        END IF
*
*. Place (T(P,K)/S(K,K)   in one-electron integral list
C                       T_ROW_TO_H(T,H,K)
        CALL T_ROW_TO_H(T,WORK(KINT1),K,TKK)
*. T_{kk}^Nk  
C            T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C)
        CALL T_TO_NK_VEC(TKK,K,ISSM,ISSPC,LUSC1,LUSC2,VEC1)
        CALL COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' output from T_TO_NK'
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
        XXNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
        WRITE(6,*) ' Norm**2 of C(ini) = ', XXNORM
*. For each orbital calculate (1+T+1/2 T^2)|0>
* + T                                      
CERR    CALL TRACI_INCORE_FCI(VEC1,VEC2,LUSC1,LUSC2,T,0,0)
        CALL MV7(VEC1,VEC2,LUSC1,LUSC2,0,0)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction vector'    
          CALL WRTVCD(VEC1,LUSC2,1,LBLK)
        END IF
        ONE = 1.0D0
        CALL VECSMDP(VEC1,VEC2,ONE,ONE,LUSC1,LUSC2,LUSC3,1,LBLK)
        CALL COPVCD(LUSC3,LUSC1,VEC1,1,LBLK)
        XXNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
        WRITE(6,*) ' Norm**2 (1+T)!Prev> = ', XXNORM

        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated vector (1+T)!0> '    
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF                    

*. + 1/2 T^2
        CALL MV7(VEC1,VEC2,LUSC2,LUSC3,0,0)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction vector (1+T+T**2)!0>'    
          CALL WRTVCD(VEC1,LUSC3,1,LBLK)
        END IF                  
        ONE = 1.0D0
        HALF  = 0.5D0
        CALL VECSMDP(VEC1,VEC2,ONE,HALF,LUSC1,LUSC3,LUSC2,1,LBLK)
*. and transfer back to LUSC1
        CALL COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
        XXNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
        WRITE(6,*) ' Norm**2 (1+T+1/2 T^2 )!Prev> = ', XXNORM
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Result af transformation of orbital: '
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
      END DO
*. Clean up tome
      IH1FORM =  IH1FORM_SAVE 
*. And transfer to LUCOUT
      CNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Norm of transformed vector', CNORM
        WRITE(6,*) ' Transformed vector'
        CALL WRTVCD(VEC1,LUSC1,1,LBLK)
      END IF
      CALL REWINO(LUSC1)
C?    WRITE(6,*) ' LUCOUT LUSC1 = ', LUCOUT,LUSC1
      CALL COPVCD(LUSC1,LUCOUT,VEC1,0,LBLK)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACID')
*
      RETURN
      END
      SUBROUTINE TPBLM2(A,AT,LBLR,NBLR,LBLC,NBLC,IFLAG,IORD)
C
C A BLOCKED MATRIX A IS GIVEN .
C
C THE BLOCK STRUCTURE CAN BE OF THE FOLLOWING TYPES
C IORD = 1 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER COLUMNS OF THIS BLOCK
C             LOOP OVER ROWS OF THIS BLOCK
C
C IORD = 2 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER ROWS OF THIS BLOCK
C             LOOP OVER COLUMNS OF THIS BLOCK
C
C     FOR IORD = 2 ARE THE INDIVIDUAL BLOCKS THUS TRANSPOSED
C
C THE COMBINATION OF TWO BLOCKS IABL AND IBBL ARE ALLOWED
C IF IFLAG(IABL,IBBL) = 1
C
C TRANSPOSE THE INDIVIDUAL BLOCKS OF THIS MATRIX TO GIVE AT
C THE ORDER OF THE BLOCKS ARE NOT CHANGED
C
C JEPPE OLSEN , NOVEMBER 1988
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),AT(*)
      DIMENSION LBLR(NBLR),LBLC(NBLC)
      DIMENSION IFLAG(NBLR,NBLC)
C
      IOFF = 1
      DO 200 IBLR = 1, NBLR
        DO 100 IBLC = 1, NBLC
          IF(IFLAG(IBLR,IBLC) .EQ. 1 ) THEN
C?          WRITE(6,*) ' BLOCK INDECES ',IBLR,IBLC
            LR = LBLR(IBLR)
            LC = LBLC(IBLC)
            IF( IORD .EQ. 1 .AND. LR*LC .GT. 0) THEN
              CALL TRPMT3(A(IOFF),LR,LC,AT(IOFF))
            ELSE IF( IORD .EQ. 2.AND. LR*LC .GT.0 ) THEN
              CALL TRPMT3(A(IOFF),LC,LR,AT(IOFF))
            END IF
            IOFF = IOFF + LR * LC
          END IF
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
      SUBROUTINE PRBLM2(A,LBLR,NBLR,LBLC,NBLC,IFLAG,IORD)
C
C A BLOCKED MATRIX A IS GIVEN .
C
C THE BLOCK STRUCTURE CAN BE OF THE FOLLOWING TYPES
C IORD = 1 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER COLUMNS OF THIS BLOCK
C             LOOP OVER ROWS OF THIS BLOCK
C
C IORD = 2 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER ROWS OF THIS BLOCK
C             LOOP OVER COLUMNS OF THIS BLOCK
C
C     FOR IORD = 2 ARE THE INDIVIDUAL BLOCKS THUS TRANSPOSED
C
C THE COMBINATION OF TWO BLOCKS IABL AND IBBL ARE ALLOWED
C IF IFLAG(IABL,IBBL) = 1
C
C PRINT THIS MATRIX !
C
C JEPPE OLSEN , NOVEMBER 1988
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)
      DIMENSION LBLR(NBLR),LBLC(NBLC)
      DIMENSION IFLAG(NBLR,NBLC)
C
      IOFF = 1
      DO 200 IBLR = 1, NBLR
        DO 100 IBLC = 1, NBLC
          IF(IFLAG(IBLR,IBLC) .EQ. 1 ) THEN
            WRITE(6,*) ' BLOCK INDECES ',IBLR,IBLC
            LR = LBLR(IBLR)
            LC = LBLC(IBLC)
            IF( IORD .EQ. 1 .AND. LR*LC .GT. 0) THEN
              CALL WRTMAT(A(IOFF),LR,LC,LR,LC)
            ELSE IF( IORD .EQ. 2.AND. LR*LC.GT.0 ) THEN
              CALL WRTMAT(A(IOFF),LC,LR,LC,LR)
            END IF
            IOFF = IOFF + LR * LC
          END IF
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
      SUBROUTINE PAMTMT(X,T,WORK,NORB)
*
* GENERATE PER AKE'S T MATRIX FROM AN
* ORBITAL ROTATION MATRIX X
*
* T IS OBTAINED AS A STRICTLY LOWER TRIANGULAR
* MATRIX TL AND AN UPPER TRIANGULAR MATRIX TU
*
*         TL = 1 - L
*         TU = U ** -1
*
* WHERE L AND U ARISES FROM A LU DECOMPOSITION OF
* X :
*         X = L * U
* WITH L BEING A LOWER TRIANGULAR MATRIX WITH UNIT ON THE
* DIAGONAL AND U IS AN UPPER TRIANGULAR MATRIX
*
* JEPPE OLSEN OCTOBER 1988
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NORB,NORB),T(NORB,NORB)
      DIMENSION WORK(*)
* DIMENSION OF WORK : NORB ** 2 + NORB*(NORB+1) / 2
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Wellcome to PAMTMT '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' Input matrix '
        CALL WRTMAT(X,NORB,NORB,NORB,NORB)
      END IF
*. Allocate local memory
      KLFREE = 1
C     KLL = KFLREE
      KLL = KLFREE
      KLFREE = KLL + NORB*(NORB+1)/2
      KLU = KLFREE
      KLFREE = KLU + NORB ** 2
*.LU factorize X
      CALL LULU(X,WORK(KLL),WORK(KLU),NORB)
*.Expand U to full matrix
      CALL SETVEC(T,0.0D0,NORB ** 2 )
      DO 10 I = 1,NORB
      DO 10 J = I,NORB
        T(I,J) = WORK(KLU-1+J*(J-1)/2+I)
   10 CONTINUE
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' MATRIX TO BE INVERTED '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*.Invert U
      CALL INVMAT(T,WORK(KLU),NORB,NORB,ISING)
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*.Subtract L
      DO 20 I = 1, NORB
      DO 20 J = 1,I-1
       T(I,J)= - WORK(KLL-1+I*(I-1)/2+J)
   20 CONTINUE
*
      IF( NTEST .GE. 2 ) THEN
        WRITE(6,*) ' INPUT X MATRIX '
        CALL WRTMAT(X,NORB,NORB,NORB,NORB)
        WRITE(6,*) ' T MATRIX '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*
      RETURN
      END
      SUBROUTINE LULU(A,L,U,NDIM)
C
C LU DECOMPOSITION OF MATRIX A
C
C     A = L * U
C
C WHERE L IS A LOWER TRIANGULAR MATRIX WITH A
C UNIT DIAGONAL AND U IS AN UPPER DIAGONAL
C
C L AND U ARE STORED AS ONE DIMENSIONAL ARRAYS
C
C   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
C
C   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
C
C THIS ADRESSING SCHEMES SUPPORTS VECTORIZATION OVER COLUMNS
C FOR L AND  OVER ROWS FOR U .
C
C
C NO PIVOTING IS DONE HERE , SO THE SCHEME GOES :
C
C     LOOP OVER R=1, NDIM
C        LOOP OVER J = R, NDIM
C          U(R,J) = A(R,J) - SUM(K=1,R-1) L(R,K) * U(K,J)
C        END OF LOOP OVER J
C
C        LOOP OVER I = R+1, NDIM
C          L(I,R) = (A(I,R) - SUM(K=1,R-1)L(I,K) * U(K,R) ) /U(R,R)
C        END OF LOOP OVER I
C     END OF LOOP OVER R
C
C JEPPE OLSEN , OCTOBER 1988
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDIM,NDIM)
      DOUBLE PRECISION  L(*),U(*)
      REAL * 8  INPROD
      INTEGER R
      CALL LUCIAQENTER('LULU ') 
C
C
      DO 1000 R = 1, NDIM
C
        DO 100 J = R, NDIM
         U(J*(J-1)/2 + R ) = A(R,J) -
     &   INPROD(L(R*(R-1)/2+1),U(J*(J-1)/2+1),R-1)
  100   CONTINUE
C
        XFACI = 1.0D0/ U(R*(R+1)/2)
        L(R*(R+1)/2 ) = 1.0D0
        DO 200 I = R+1, NDIM
          L(I*(I-1)/2 + R) = (A(I,R) -
     &   INPROD(L(I*(I-1)/2+1),U(R*(R-1)/2+1),R-1) ) * XFACI
  200  CONTINUE
C
 1000 CONTINUE
C
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' L MATRIX '
         CALL PRSYM(L,NDIM)
         WRITE(6,*) ' U MATRIX ( TRANSPOSED ) '
         CALL PRSYM(U,NDIM)
      END IF
*
      CALL LUCIAQEXIT('LULU ')
*
      RETURN
      END
      SUBROUTINE GET_3BLKS(KVEC1,KVEC2,KC2)
*
* Allocate the three blocks VEC1, VEC2, C2 used in sigma, densi etc
*
* Jeppe Olsen, Jan 1998
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'cicisp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'csm.inc' 
      INCLUDE 'cstate.inc' 
      INCLUDE 'crun.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'spinfo.inc'
*
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
*. Common block for communicating with sigma
      INCLUDE 'cands.inc'
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_3B')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ') !done
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ') !done
*
      ISPC = MAX(ICSPC,ISSPC)
      ISM  = ISSM
      CALL IAIBCM(ISPC,dbl_mb(KLCIOIO))
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,dbl_mb(KLCBLTP),int_mb(KSVST))
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR           
*. alpha and beta strings with an electron removed
      IATPM1 = 3 
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5 
      IBTPM2 = 6
*. Largest number of strings of given symmetry and type
      MAXA = MXNSTR
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESCPH(int_mb(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &              NSMST,NSTFSMSPGP,MXPNSMST,
     &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &              NELFSPGP,
     &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &              IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,
     &              MXADKBLK_AS,MX_NSPII)
      IF(ISIMSYM.EQ.1) MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
      IF(IPRCIX.GE.2) THEN
        WRITE(6,*) 'GET_3BL: MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
        WRITE(6,*) 'GET_3BL: MXADKBLK ', MXADKBLK
        WRITE(6,*) ' Space for resolution matrices ',LSCR2
      END IF
*
      IF(ISIMSYM.EQ.0) THEN 
        LBLOCK = MXSOOB
      ELSE
        LBLOCK = MXSOOB_AS
      END IF
      IF(NOCSF.EQ.0) LBLOCK = NSD_FOR_OCCLS_MAX
      IF(IPRCIX.GE.2) WRITE(6,*) ' GET_3BLK: LBLOCK = ', LBLOCK
*
      LBLOCK = MAX(LBLOCK,LCSBLK)
      LSCR12 = MAX(LBLOCK,2*LSCR2)  
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_3B')
*
      CALL MEMMAN(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
      CALL MEMMAN(KVEC2,LBLOCK,'ADDS  ',2,'VEC2  ')
      CALL MEMMAN(KC2,LSCR12,'ADDS  ',2,'C2    ')
*
      RETURN
      END
*
      SUBROUTINE TRIPAK_BLKMG(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK,INFO)
*
* Switch between packed and unpacked form of a blocked matrix
*
* IWAY = 1 => Unpacked to packed
* IWAY = 2 => Packed to unpacked 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input/output
      integer AUTPAK
      DIMENSION APAK(*)
*. Input
      DIMENSION LBLOCK(NBLOCK)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFFU = 1
          IOFFP = 1
        ELSE
          IOFFU = IOFFU + LBLOCK(IBLOCK-1)**2
          IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
        END IF
        L = LBLOCK(IBLOCK)
        IJ = IOFFP
        IJF = IOFFU
        DO I = 1,L
          call ga_get(AUTPAK,IJF,IJF+I-1,1,1,APAK(IJ),1)
          IJ = IJ + I
          IJF = IJF + L
        END DO
      END DO
      call ga_sync()

c    CALL TRIPAK(AUTPAK(IOFFU),APAK(IOFFP),IWAY,L,L)
c    SUBROUTINE TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
c     DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)

c     IF( IWAY .EQ. 1 ) THEN
c       IJ = 0
c       DO I = 1,NDIM
c         DO J = 1, I
c           APAK(IJ+J) = AUTPAK(J,I)
c         END DO
c         IJ = IJ + I
c       END DO

*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from TRIPAK_BLKM'
        WRITE(6,*) 
        WRITE(6,*) ' Unpacked matrix'
        WRITE(6,*) ' ==============='
        call ga_print(autpak)
        WRITE(6,*) ' Packed matrix '
        WRITE(6,*) ' =============='
        CALL APRBLM2(  APAK,LBLOCK,LBLOCK,NBLOCK,1)
      END IF
*
      RETURN
      END
*
      SUBROUTINE TRIPAK_BLKMG2(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK,INFO)
*
* Switch between packed and unpacked form of a blocked matrix
*
* IWAY = 1 => Unpacked to packed
* IWAY = 2 => Packed to unpacked 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input/output
      integer APAK
      DIMENSION AUTPAK(*)
*. Input
      DIMENSION LBLOCK(NBLOCK)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFFU = 1
          IOFFP = 1
        ELSE
          IOFFU = IOFFU + LBLOCK(IBLOCK-1)**2
          IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
        END IF
        L = LBLOCK(IBLOCK)
        IJ = IOFFP
        IJF = IOFFU
        DO I = 1,L
          call ga_put(APAK,IJ,IJ+I-1,1,1,AUTPAK(IJF),1)
          IJ = IJ + I
          IJF = IJF + L
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from TRIPAK_BLKM'
        WRITE(6,*)
        WRITE(6,*) ' Unpacked matrix'
        WRITE(6,*) ' ==============='
        call ga_print(autpack)
        WRITE(6,*) ' Packed matrix '
        WRITE(6,*) ' =============='
        CALL APRBLM2(  APAK,LBLOCK,LBLOCK,NBLOCK,1)
      END IF
*
      RETURN
      END
*
      SUBROUTINE TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK,INFO)
*
* Switch between packed and unpacked form of a blocked matrix
*
* IWAY = 1 => Unpacked to packed
* IWAY = 2 => Packed to unpacked 
*
* Jeppe Olsen, February 1, 1998 (Moensted Kalkgrubber )
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input/output
      DIMENSION AUTPAK(*),APAK(*)
*. Input
      DIMENSION LBLOCK(NBLOCK)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFFU = 1
          IOFFP = 1
        ELSE
          IOFFU = IOFFU + LBLOCK(IBLOCK-1)**2
          IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
        END IF
        L = LBLOCK(IBLOCK)
        CALL TRIPAK(AUTPAK(IOFFU),APAK(IOFFP),IWAY,L,L)                  
C            TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from TRIPAK_BLKM'
        WRITE(6,*)
        WRITE(6,*) ' Unpacked matrix'
        WRITE(6,*) ' ==============='
        CALL APRBLM2(AUTPAK,LBLOCK,LBLOCK,NBLOCK,0)
C            APRBLM2(A,LROW,LCOL,NBLK,ISYM)
        WRITE(6,*) ' Packed matrix '
        WRITE(6,*) ' =============='
        CALL APRBLM2(  APAK,LBLOCK,LBLOCK,NBLOCK,1)
      END IF
*
      RETURN
      END
*
      SUBROUTINE Z_BLKFO(ISPC,ISM,IATP,IBTP,KPCLBT,KPCLEBT,
     &                   KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM,IATP,IBTP.
*
* Output is given in the form of pointers to vectors in WORK
* where the info is stored : 
*
* KPCLBT : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
* KPCI1BT : Length of each block                   
* KPCIBT  : Info on each block
* KPCBLTP : BLock type for each symmetry
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* Jeppe Olsen, Feb. 98
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'cicisp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'crun.inc'
*
      NTEST = 00  
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Output from Z_BLKFO '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Pointers to output arrays
      NTTS = MXNTTS
      CALL MEMMAN(KPCLBT ,MXNTTS,'ADDL  ',1,'CLBT  ') !done
      CALL MEMMAN(KPCLEBT,MXNTTS,'ADDL  ',1,'CLEBT ') !done
      CALL MEMMAN(KPCI1BT,MXNTTS,'ADDL  ',1,'CI1BT ') !done
      CALL MEMMAN(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBT  ') !done
      CALL MEMMAN(KPCBLTP,NSMST,'ADDL  ',2,'CBLTP ') !done
*.    ^ These should be preserved after exit so put mark for flushing here
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'Z_BLKF')
*. Info needed for generation of block info
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ') !done
      CALL IAIBCM(ISPC,dbl_mb(KLCIOIO))
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,dbl_mb(KPCBLTP),int_mb(KSVST))
*. Allowed length of each batch
      IF(ISIMSYM.EQ.0) THEN
        LBLOCK = MXSOOB
      ELSE
        LBLOCK = MXSOOB_AS
      END IF
*
      LBLOCK = MAX(LBLOCK,LCSBLK)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' LBLOCK = ', LBLOCK
      END IF
*
*. Batches  of C vector
      CALL PART_CIV2(IDC,dbl_mb(KPCBLTP),int_mb(KNSTSO(IATP)),
     &              int_mb(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,dbl_mb(KLCIOIO),
     &              ISMOST(1,ISM),
     &              NBATCH,int_mb(KPCLBT),int_mb(KPCLEBT),
     &              int_mb(KPCI1BT),int_mb(KPCIBT),0,ISIMSYM)
*. Number of BLOCKS
      NBLOCK = IFRMR(int_mb(KPCI1BT),1,NBATCH)
     &       + IFRMR(int_mb(KPCLBT),1,NBATCH) - 1
      IF(NTEST.GE.1) THEN
         WRITE(6,*) ' Number of batches', NBATCH
         WRITE(6,*) ' Number of blocks ', NBLOCK
      END IF
*. Length of each block
      CALL EXTRROW(int_mb(KPCIBT),8,8,NBLOCK,int_mb(KPCI1BT))
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'Z_BLKF')
      RETURN
      END
      SUBROUTINE T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C)
*
* Evaluate T**(NK_operator) times vector on file LUIN
* to yield vector on fiel LUOUT
* (NK_operator is number operator for orbital K )             
*
* Note LUCIN and LUCOUT are both rewinded before read/write
* Input 
* =====
*  T : Input constant
*  KORB : Orbital in symmetry order
*
*  ISM,ISPC : Symmetry and space of state on LUIN
*  C : Scratch block
*
*
* Jeppe Olsen, Feb. 98
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'  
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'strinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'csm.inc'

*. Scratch block, must hold a batch of blocks
      DIMENSION C(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T_TO_NK_VEC speaking '
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
        WRITE(6,*) ' T = ', T
      END IF
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'T_TO_N')
*. Set up block and batch structure of vector
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      CALL Z_BLKFO(ISPC,ISM,IATP,IBTP,KLCLBT,KLCLEBT,
     &            KLCI1BT,KLCIBT,KLCBLTP,NBATCH,NBLOCK)
C           Z_BLKFO(ISPC,ISM,IATP,IBTP,KPCLBT,KPCLEBT,
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
      CALL MEMMAN(KLKAOC,MXNSTR,     'ADDL  ',1,'KLKAOC')
      CALL MEMMAN(KLKBOC,MXNSTR,     'ADDL  ',1,'KLKBOC')
*. Orbital K in type ordering
      KKORB = IREOST(KORB)
      CALL T_TO_NK_VECS
     &  (T,KKORB,C,LUCIN,LUCOUT,int_mb(KNSTSO(IATP)),
     &      int_mb(KNSTSO(IBTP)),
     &      NBLOCK,int_mb(KLCIBT),
     &      NAEL,NBEL,int_mb(KLASTR),int_mb(KLBSTR),
     &      dbl_mb(KLCBLTP),NSMST,
     &      ICISTR,NTOOB,int_mb(KLKAOC),int_mb(KLKBOC))
   
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'T_TO_N')
*
      RETURN
      END
      SUBROUTINE T_TO_NK_VECS(T,KORB,C,LUCIN,LUCOUT,NSSOA,NSSOB,
     &                 NBLOCK,IBLOCK,
     &                 NAEL,NBEL,IASTR,IBSTR,IBLTP,NSMST,
     &                 ICISTR,NORB,IKAOCC,IKBOCC)
*
* Multiply Vector in LUCIN with t **NK_op to yield vector on LUCOUT
*
* Both files are initially rewinded
*
*
* Jeppe Olsen, Feb. 1998    
*                                                  

      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)  
*. Scratch
      DIMENSION C(*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IKAOCC(*),IKBOCC(*)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
      DIMENSION IBLTP(*)
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from T_TO_NK_VECS '
        WRITE(6,*) ' ========================'
      END IF
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Input file '
       WRITE(6,*) ' =========== '
       WRITE(6,*)
       CALL WRTVCD(C,LUCIN,1,-1)
      END IF
*
      CALL REWINO(LUCIN)
      CALL REWINO(LUCOUT)
*
      T2 = T**2
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' T and T2 in action ', T, T2
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
C?      WRITE(6,*) ' IATP IBTP IASM IBSM ', IATP,IBTP,IASM,IBSM
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
*. Occupation of orb KORB
        DO JSTR = 1, NASTR1
          KOCC = 0
          DO JAEL = 1, NAEL
            IF(IASTR(JAEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKAOCC(JSTR) = KOCC
        END DO
C?      WRITE(6,*) ' IKAOCC array '
C?      CALL IWRTMA(IKAOCC,1,NASTR1,1,NASTR1)
     
    
*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                           NORB,0,IDUM,IDUM)
C?      WRITE(6,*) ' After GETSTR, NBSTR1=',NBSTR1
*. Occupation of orb KORB
        DO JSTR = 1, NBSTR1
C?        write(6,*) ' JSTR = ', JSTR
          KOCC = 0
          DO JBEL = 1, NBEL
C?          write(6,*) JBEL, IBSTR(JBEL,JSTR)
            IF(IBSTR(JBEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKBOCC(JSTR) = KOCC
        END DO
C?      WRITE(6,*) ' IKBOCC array '
C?      CALL IWRTMA(IKBOCC,1,NBSTR1,1,NBSTR1)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
C?      WRITE(6,*) ' IBLTP ', IBLTP(IASM)
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
C?      WRITE(6,*) ' NIA NIB ', NIA,NIB
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN 
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(LDET,1,-1,LUCIN)
          CALL FRMDSC(C,LDET,-1,LUCIN,IMZERO,IAMPACK)
        END IF
        IF(IMZERO.NE.1) THEN
*
          IDET = 0
          DO  IB = 1,NIB
            IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
              MINIA = IB 
            ELSE
              MINIA = 1     
            END IF
            DO  IA = MINIA,NIA
*
              IDET = IDET + 1
C?            WRITE(6,*) ' IA IB IDET',IA,IB,IDET
              KABOCC = IKAOCC(IA)+IKBOCC(IB)
              IF(KABOCC.EQ.1) THEN
                C(IDET) = T*C(IDET)
              ELSE IF(KABOCC.EQ.2) THEN
                C(IDET) = T2 *C(IDET)
              END IF
            END DO
*           ^ End of loop over alpha strings
          END DO
*         ^ End of loop over beta strings
*
        END IF
*       ^ End of if statement for nonvanishing blocks
*. Save result on LUCOUT
        CALL ITODS(LDET,1,-1,LUCOUT)
        CALL TODSC(C,LDET,-1,LUCOUT)
      END DO
*     ^ End of loop over blocks
*. This is the end, the end of every file my friend, the end
      CALL ITODS(-1,1,-1,LUCOUT)
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Output file '
       WRITE(6,*) ' =========== '
       WRITE(6,*)
       CALL WRTVCD(C,LUCOUT,1,-1)
      END IF
*

      RETURN
      END
*
      SUBROUTINE ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
*
* Advice Sigma routine about best route to take
*
* ITERM : Term  to be studied :  
*         =1 alpha-beta term 
*         ....... ( to be continued )
*
* LADVICE : ADVICE given ( short, an integer !!)
*
* For ITERM = 1 : 
*           LADVICE = 1 : Business as usual, no transpose of matrix
*                         (resolution on alpha strings, direct exc on beta)
*           LADVICE = 2 = Transpose matrices
*                         (resolution on beta strings, direct exc on alpha)
*
* Jeppe Olsen, Tirstrup Airport, Jan 12, 98
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Specific input
      INTEGER IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
*. Local Scratch
       DIMENSION ITP(16),JTP(16),KTP(16),LTP(16)
*
      NTEST = 00
      IF(ITERM.EQ.1) THEN
*.
*. sigma(i,Ka,Ib) = sum(i,kl)<Ib!Eb_kl!Jb>(ij!kl)C(j,Ka,Jb)
*
* Number of ops : Number of sx(kl) N_i*N_j_dimension of C(j,Ka,Jb)
*.No absolute calc of flops is made, only a relative measure
*
* Single excitations connecting the two types
*
C            SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
        CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOCC,JAOCC,IPHGAS)
        CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOCC,JBOCC,IPHGAS)
C?      WRITE(6,*) 'NIJTYP, NKLTYP', NIJTYP,NKLTYP
*. P-h modifications ( I cannot predict these at the moment
        IF(NIJTYP.GE.1.AND.NKLTYP.GE.1) THEN
*
        IF((IPHGAS(ITP(1)).EQ.2.AND.IPHGAS(JTP(1)).EQ.2).OR.
     &     (IPHGAS(KTP(1)).EQ.2.AND.IPHGAS(LTP(1)).EQ.2)     ) THEN
           IPHMODI = 1
         ELSE
           IPHMODI = 0
         END IF
        ELSE
           IPHMODI = 0
        END IF
          
*
        IF(IPHMODI.EQ.1.OR.NIJTYP.NE.1.OR.NKLTYP.NE.1
     &     .OR.IADVICE.EQ.0) THEN
*. Several connections, i.e. the alpha or the beta blocks are identical,
*. or ph modifications
*. just continue
          LADVICE = 1
        ELSE
* =========================================
*.. Index for flops along C(j,Ka,Jb) route
* =========================================
*.Dim of C(j,Ka,Jb) relative to C(Ja,Jb)
*. going from Ja to  Ka reduces occ by one elec, changes dim by n/(N-n+1)
          XNJOB = FLOAT(NOBPT(JTP(1)))
          XNJEL = FLOAT(JAOCC(JTP(1)))
          XCJKAJB = XNJOB*XNJEL/(XNJOB-XNJEL+1)
*. Number of kl excitations per beta string : 
          XNKLSX = FLOAT((NOBPT(KTP(1))-JBOCC(KTP(1)))*JBOCC(LTP(1)))
*. Number of ops (relative to dim of C)
          XNIOB = FLOAT(NOBPT(ITP(1)))
          XFLOPA = XCJKAJB*XNKLSX*XNIOB
* =========================================
*.. Index for flops along C(l,Ja,Kb) route
* =========================================
*.Dim of C(l,Ja,Kb) relative to C(Ja,Jb)
          XNLOB = FLOAT(NOBPT(LTP(1)))
          XNLEL = FLOAT(JBOCC(LTP(1)))
          XCLJAKB = XNLOB*XNLEL/(XNLOB-XNLEL+1)
*. Number of ij excitations per alpha string : 
          XNIJSX = FLOAT((NOBPT(ITP(1))-JAOCC(ITP(1)))*JAOCC(JTP(1)))
*. Number of ops (relative to dim of C)
          XNKOB = FLOAT(NOBPT(KTP(1)))
          XFLOPB = XCLJAKB*XNIJSX*XNKOB
*. Switch to second route if atleast 20 percent less work
          IF(XFLOPB.LE.0.8*XFLOPA) THEN
            LADVICE = 2
          ELSE
            LADVICE = 1
          END IF
*. Well, an additional consideration :
* If the C block involes the smallest allowed number of elecs in hole space,
* and the annihilation is in hole space
* then we do the annihilation in the space with the smallest number of 
* hole electrons.
          LHOLEA =0
          LHOLEB =0
          DO IGAS = 1, NGAS
            IF(IPHGAS(IGAS).EQ.2) THEN
              LHOLEA = LHOLEA + JAOCC(IGAS)
              LHOLEB = LHOLEB + JBOCC(IGAS)
            END IF
          END DO
*
          IF(LHOLEA+LHOLEB.EQ.MNHL.AND.
     &       (IPHGAS(JTP(1)).EQ.2.OR.IPHGAS(LTP(1)).EQ.2))  THEN
*
             IF(IPHGAS(JTP(1)).EQ.2) THEN
              KHOLEA = LHOLEA-1
              KHOLEB = LHOLEB 
             ELSE 
              KHOLEA = LHOLEA
              KHOLEB = LHOLEB - 1
             END IF
*
             IF(KHOLEA.EQ.KHOLEB) THEN
               LLADVICE = LADVICE
             ELSE IF(KHOLEA.LT.KHOLEB) THEN
               LLADVICE= 1
             ELSE
               LLADVICE = 2
             END IF
             IF(NTEST.GE.100.AND.LADVICE.NE.LLADVICE) THEN
               WRITE(6,*) ' Advice changed by hole considetions'
               WRITE(6,*) ' LADVICE, LLADVICE', LADVICE,LLADVICE
             END IF
             LADVICE = LLADVICE  
          END IF
*
*
C         IF(NTEST.GE.100) THEN
          IF(NTEST.GE.100.AND.LADVICE.EQ.2) THEN
            WRITE(6,*) ' ADVICE active '
            WRITE(6,*) ' IAOCC IBOCC JAOCC JBOCC'
            CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
            WRITE(6,*) ' ITP JTP KTP LTP ',ITP(1),JTP(1),KTP(1),LTP(1)
            WRITE(6,'(A,4(2X,E9.3))') 
     &      ' XFLOPA,XFLOPB,XNJOB,XNLOB', XFLOPA,XFLOPB,XNJOB,XNLOB
            WRITE(6,*) ' XCLJAKB*XNKOB, XCJKAJB*XNIOB ', 
     &                   XCLJAKB*XNKOB, XCJKAJB*XNIOB
            WRITE(6,*) ' ADVICE given : ', LADVICE
          END IF
        END IF
*       ^ End if several types/ph modi
      END IF
*     ^ End if ITERM test ( type of excitation)
C     WRITE(6,*) ' MEMCHECK at end of ADVICE'
C     CALL MEMCHK
C     WRITE(6,*) ' MEMCHECK passed '
      RETURN
      END
      SUBROUTINE RSBB2A(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NSCOL,
     &                  NGAS,ISOC,ICOC,
     &                  SB,CB,
     &                  ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGASX,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBSX,
     &                  RIKSX,RKJSX,MXSXST,MXSXBL,IMOC,SCLFAC,NTESTG,
     &                  NSEL2E,ISEL2E,IUSE_PH,IPHGAS,XINT2)
*
* two electron excitations on column strings
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NSCOL : Number of columns in S block
* ISEL1(3) : Number of electrons in RAS1(3) for S block
* ICEL1(3) : Number of electrons in RAS1(3) for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* DXSTST : Sym of sx,!st> => sym of sx !st>
* STSTDX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* SXDXSX : Symmetry of SX1,SX1*SX2 => symmetry of SX2
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals, strings,
*                           single excitations, double excitations
* MAXI   : Largest number of 'spectator strings' treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
*
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nproc
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),DXSTST(NSMDX,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER IPHGAS(NGAS)
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
      INTEGER ISOC(NGAS),ICOC(NGAS)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),XINT(*), XINT2(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
      DIMENSION RIKSX(MXSXBL,4),RKJSX(MXSXBL,4)
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
      INTEGER I4_DIM(4),I4_SM(4),I4_TP(4),I4_REO(4),ISCR(4)
      INTEGER I4_AC(4)
*
      INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
      COMMON/SOMESCR/SCR(MXPTSOB*MXPTSOB*MXPTSOB*MXPTSOB) 
*
      COMMON/XXTEST/ISETVECOPS(10)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'intform.inc'
*
      DIMENSION IACAR(2),ITPAR(2)
      CALL LUCIAQENTER('RS2A') 
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' RSBB2A speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' MXSXST MXSXBL = ', MXSXST,MXSXBL
        WRITE(6,*) ' RSBB2A : IMOC,IUSE_PH ', IMOC, IUSE_PH
        WRITE(6,*) ' ISOC and ICOC : '
        CALL IWRTMA(ISOC,1,NGAS,1,NGAS)
        CALL IWRTMA(ICOC,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' Memcheck at start of RSBB2A '
        CALL MEMCHK 
        WRITE(6,*) ' Memcheck passed '
*
      END IF
      IFRST = 1 
      JFRST = 1
*
*.Types of DX that connects the two strings
*
      IDXSM = STSTDX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0)  GOTO 2001
*. Connecting double excitations
      CALL DXTYP2_GAS(NDXTYP,ITP,JTP,KTP,LTP,NGAS,ISOC,ICOC,IPHGAS)
c     nproc = ga_nnodes()
      DO 2000 IDXTYP = 1, NDXTYP
        ITYP = ITP(IDXTYP)
        JTYP = JTP(IDXTYP)
        KTYP = KTP(IDXTYP)
        LTYP = LTP(IDXTYP)
*. Is this combination of types allowed
         IJKL_ACT = I_DX_ACT(ITYP,KTYP,LTYP,JTYP)
         IF(IJKL_ACT.EQ.0) GOTO 2000
*      
C?      write(6,*) ' test inserted in RSBB2AN'
C?      NPTOT = 0
C?      IF(ITYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(JTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(KTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(LTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(NPTOT.EQ.3) GOTO 2000
*
        ITYP_ORIG = ITYP
        JTYP_ORIG = JTYP
        KTYP_ORIG = KTYP
        LTYP_ORIG = LTYP
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' ITYP_ORIG, JTYP_ORIG, KTYP_ORIG, LTYP_ORIG',
     &                 ITYP_ORIG, JTYP_ORIG, KTYP_ORIG, LTYP_ORIG
        END IF
        NIJKL1 = 0
        IF(ITYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(JTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(KTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(LTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(NIJKL1.EQ.0) CALL LUCIAQENTER('BB2A0')
        IF(NIJKL1.EQ.1) CALL LUCIAQENTER('BB2A1')
        IF(NIJKL1.EQ.2) CALL LUCIAQENTER('BB2A2')
        IF(NIJKL1.EQ.3) CALL LUCIAQENTER('BB2A3')
        IF(NIJKL1.EQ.4) CALL LUCIAQENTER('BB2A4')
*. Optimal ordering of operators 
        I4_AC(1) = 2
        I4_AC(2) = 2
        I4_AC(3) = 1
        I4_AC(4) = 1
        I4_TP(1) = ITYP
        I4_TP(2) = KTYP
        I4_TP(3) = LTYP
        I4_TP(4) = JTYP
        IF(IUSE_PH.EQ.1) THEN
          NOP = 4
          CALL ALG_ROUTERX(ISOC,JSOC,NOP,I4_TP,I4_AC,I4_REO,SIGN4)
        ELSE
          DO IJKL = 1, 4
            I4_REO(IJKL) = IJKL
          END DO
          SIGN4 = 1.0D0
        END IF
*. Type of operators : TP and AC
        DO IJKL = 1, 4
C         ISCR( I4_REO(IJKL) ) = I4_TP(IJKL)
          ISCR(IJKL) = I4_TP( I4_REO(IJKL) )
        END DO
        DO IJKL = 1, 4
          I4_TP(IJKL) = ISCR(IJKL)
        END DO
        ITYP = I4_TP(1)
        KTYP = I4_TP(2)
        LTYP = I4_TP(3)
        JTYP = I4_TP(4)
        DO IJKL = 1, 4
          ISCR(IJKL) = I4_AC( I4_REO(IJKL) )
        END DO
        DO IJKL = 1, 4
          I4_AC(IJKL) = ISCR(IJKL)
        END DO
        IF(NTEST.GE.500) THEN 
          WRITE(6,*) ' I4_AC, IT_TP  defined '
          WRITE(6,*) ' I4_AC, I4_TP '
          CALL IWRTMA(I4_AC,1,4,1,4)
          CALL IWRTMA(I4_TP,1,4,1,4)
        END IF
*
*      ==================================
        IF(I4_AC(1).EQ.I4_AC(2) ) THEN
*      ==================================
*
*. a+ a+ a a or a a a+ a+
*. Largest possible number of orbital pairs
          MI = 0
          MJ = 0
          MK = 0
          ML = 0
          DO IOBSM = 1, NSMST
            MI = MAX(MI,NOBPTS(ITYP,IOBSM))
            MJ = MAX(MJ,NOBPTS(JTYP,IOBSM))
            MK = MAX(MK,NOBPTS(KTYP,IOBSM))
            ML = MAX(ML,NOBPTS(LTYP,IOBSM))
          END DO
          MXPAIR = MAX(MI*MK,MJ*ML)
*. Largest posssible 
*. Symmetry of allowed Double excitation,loop over excitations
          DO 1950 IKOBSM = 1, NSMOB
            JLOBSM = SXDXSX(IKOBSM,IDXSM)
            IF(NTEST.GE.500) WRITE(6,*) ' IKOBSM,JLOBSM', IKOBSM,JLOBSM
            IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined 
            KFRST = 1
*
*. Number of batchs of symmetry pairs IK
*
            LENGTH = 0
            NIKBT = 0
            NBLK = 0 
            NBLKT = 0 
            DO ISM = 1, NSMOB
              KSM = ADSXA(ISM,IKOBSM)
              NI = NOBPTS(ITYP,ISM)
              NK = NOBPTS(KTYP,KSM)
              IF(NTEST.GE.500) write(6,*) ' NI, NK' , NI,NK
*
              IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                NIK = NI*(NI+1)/2
              ELSE IF(ITYP.GT.KTYP.OR.(ITYP.EQ.KTYP.AND.ISM.GT.KSM))THEN
                NIK = NI*NK
              ELSE
                NIK = 0
              END IF
              IF(NIK.NE.0) THEN
                NBLKT = NBLKT + 1
                IF(LENGTH+NIK .GT. MXPAIR) THEN
*. The present batch is complete
                  NIKBT = NIKBT+1      
                  IKBT(1,NIKBT) = NBLKT - NBLK 
                  IKBT(2,NIKBT) = NBLK
                  IKBT(3,NIKBT) = LENGTH
                  LENGTH = 0
                  NBLK = 0
                END IF
                NBLK = NBLK + 1
                LENGTH = LENGTH + NIK
                IKSMBT(1,NBLKT) = ISM
                IKSMBT(2,NBLKT) = KSM
              END IF
            END DO
*. The last batch 
            IF(NBLK.NE.0) THEN
              NIKBT = NIKBT+1     
              IKBT(1,NIKBT) = NBLKT - NBLK + 1
              IKBT(2,NIKBT) = NBLK
              IKBT(3,NIKBT) = LENGTH
            END IF

*. 
            IF(NTEST.GE.2000) THEN 
              WRITE(6,*) ' ITYP, KTYP, IKOBSM,  NIKBT = ',
     &                     ITYP, KTYP, IKOBSM,  NIKBT 
              WRITE(6,*) ' IKBT : Offset, number, length '
              DO JIKBT = 1, NIKBT 
                WRITE(6,'(3i3)') (IKBT(II,JIKBT), II = 1, 3)
              END DO
              WRITE(6,*) ' IKSMBT '
              CALL IWRTMA(IKSMBT,2,NBLKT,2,8)
            END IF
*
*. Number of batchs of symmetry pairs JL
*
            LENGTH = 0
            NJLBT = 0
            NBLK = 0 
            NBLKT = 0 
            DO JSM = 1, NSMOB
              LSM = ADSXA(JSM,JLOBSM)
              NJ = NOBPTS(JTYP,JSM)
              NL = NOBPTS(LTYP,LSM)
*
              IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                NJL = NJ*(NJ+1)/2
              ELSE IF(JTYP.GT.LTYP.OR.(JTYP.EQ.LTYP.AND.JSM.GT.LSM))THEN
                NJL = NJ*NL
              ELSE
                NJL = 0
              END IF
              IF(NJL.NE.0) THEN
                NBLKT = NBLKT + 1
                IF(LENGTH+NJL .GT. MXPAIR) THEN
*. The present batch is complete
                  NJLBT = NJLBT+1      
                  JLBT(1,NJLBT) = NBLKT - NBLK 
                  JLBT(2,NJLBT) = NBLK
                  JLBT(3,NJLBT) = LENGTH
                  LENGTH = 0
                  NBLK = 0
                END IF
                NBLK = NBLK + 1
                LENGTH = LENGTH + NJL
                JLSMBT(1,NBLKT) = JSM
                JLSMBT(2,NBLKT) = LSM
              END IF
            END DO
*. The last batch 
            IF(NBLK.NE.0) THEN
              NJLBT = NJLBT+1     
              JLBT(1,NJLBT) = NBLKT - NBLK + 1
              JLBT(2,NJLBT) = NBLK
              JLBT(3,NJLBT) = LENGTH
            END IF
*. 
            IF(NTEST.GE.2000) THEN 
              WRITE(6,*) ' JTYP, LTYP, JLOBSM,  NJLBT = ',
     &                     JTYP, LTYP, JLOBSM,  NJLBT 
              WRITE(6,*) ' JLBT : Offset, number, length '
              DO JJLBT = 1, NJLBT 
                WRITE(6,'(3i3)') (JLBT(II,JJLBT), II = 1, 3)
              END DO
              WRITE(6,*) ' JLSMBT '
              CALL IWRTMA(JLSMBT,2,NBLKT,2,8)
            END IF
*
*. Loop over batches of IK strings
            DO 1940 IKBTC = 1, NIKBT
              IF(NTEST.GE.1000) WRITE(6,*) ' IKBTC = ', IKBTC
*. Loop over batches of JL strings 
              DO 1930 JLBTC = 1, NJLBT
*. Loop over batches of I strings
*.Obtain two electron integrals xint(ik,jl) = (ij!kl)-(il!kj)
                IXCHNG = 1
* Obtain integrals in ik batch
                NIKT = IKBT(3,IKBTC)
                NJLT = JLBT(3,JLBTC)
                JLOFF = 1
                DO JLPAIR = 1, JLBT(2,JLBTC)
                IKOFF = 1
                DO IKPAIR = 1, IKBT(2,IKBTC)
*
                  ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                  KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                  JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                  LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
*
                  IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                    IKSM = 1
                    NIK =
     &              NOBPTS(ITYP,ISM)*(NOBPTS(ITYP,ISM)+1)/2
                  ELSE
                    IKSM = 0
                    NIK =
     &              NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
                  END IF
*
                  IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                    JLSM = 1
                    NJL =
     &              NOBPTS(JTYP,JSM)*(NOBPTS(JTYP,JSM)+1)/2
                  ELSE
                    JLSM = 0
                    NJL =
     &              NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM)
                  END IF
* ================================================================
*. Required form of integrals : Coulomb - Exchange of just Coulomb
* ================================================================
                  ICOUL = 0
*. Use coulomb - exchange 
                  IXCHNG = 1
*. fetch integrals
                  ONE = 1.0D0
                  IF(IH2FORM.EQ.1) THEN
*. Full conjugation symmetry, do do not worry
                  CALL LUCIAGETINT(SCR,ITYP,ISM,JTYP,JSM,KTYP,KSM,
     &                        LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL,
     &                        ONE,ONE)
                  ELSE IF(IH2FORM.EQ.2) THEN
*. Integrals do not neccessarily have full conjugation symmetry 
                  IF(I4_AC(1).EQ.2) THEN
* a + a+ a a
                    CALL LUCIAGETINT(SCR,ITYP,ISM,JTYP,JSM,KTYP,
     &                          KSM,
     &                          LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL,
     &                          ONE,ONE)
                  ELSE
*. a a a+ a+ : Obtain (jl|ik) and transpose
C?                        WRITE(6,*) ' Memcheck before GETINT '
C?                        CALL MEMCHK
C?                        WRITE(6,*) ' Check passes '
                    CALL LUCIAGETINT(XINT2,JTYP,JSM,ITYP,ISM,LTYP,
     &                          LSM,
     &                          KTYP,KSM,IXCHNG,JLSM,IKSM,ICOUL,
     &                          ONE,ONE)
C?                        WRITE(6,*) ' Memcheck before TRPMT3 '
C?                        CALL MEMCHK
C?                        WRITE(6,*) ' Check passes '
                    CALL TRPMT3(XINT2,NJL,NIK,SCR)
                   END IF
                  END IF
*                 ^ End if Hamiltonian without cc symmetry is used
*                   used
                  DO JL = 1, NJL
                    CALL COPVEC(SCR((JL-1)*NIK+1),
     &                   XINT((JLOFF-1+JL-1)*NIKT+IKOFF),NIK)
                  END DO
                  IKOFF = IKOFF + NIK
                END DO
                JLOFF = JLOFF + NJL
                END DO
*               ^ End if integrals should be fetched
*
                NPART = NROW/MAXI
                IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
c               if (npart.lt.nproc) then
c                  npart = min(nrow,nproc)
c               endif
                IF(NTEST.GE.100)
     &          write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
c               DO 1801 IIPART = ga_nodeid()+1, NPART,nproc
                DO 1801 IIPART = 1, NPART
                  IBOT = 1+(IIPART-1)*MAXI
                  ITOP = MIN(IBOT+MAXI-1,NROW)
                  NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                  KBOT = 1- MAXK
                  KTOP = 0
 1800             CONTINUE
                    KBOT = KBOT + MAXK
                    KTOP = KTOP + MAXK
*
                    IONE = 1
                    JLBOFF = 1
                    NJLT = JLBT(3,JLBTC)
                    DO JLPAIR = 1, JLBT(2,JLBTC)
                      JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                      LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
                      NJ = NOBPTS(JTYP,JSM)
                      NL = NOBPTS(LTYP,LSM)
                      IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                        NJL = NJ*(NJ+1)/2
                        JLSM = 1
                      ELSE
                        NJL = NJ * NL
                        JLSM = 0
                      END IF
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
CT                    CALL LUCIAQENTER('ADADS')
                      II12 = 1
                      K12 = 1
                      IONE = 1
C?       write(6,*) ' Before ADAADAST '
*. Creation / annihilation maps , conjugated of above
                      IF(I4_AC(4).EQ.1) THEN
                        JAC = 2
                      ELSE 
                        JAC = 1
                      END IF 
                      IF(I4_AC(3).EQ.1) THEN
                        LAC = 2
                      ELSE 
                        LAC = 1
                      END IF 
                      CALL ADAADAST_GAS(IONE,JSM,JTYP,NJ,JAC,
     &                                  IONE,LSM,LTYP,NL,LAC,
     &                            ICCTP,ICCSM,IGRP,
     &                            KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                            JFRST,KFRST,II12,K12,SCLFAC)
                      JFRST = 0
                      KFRST = 0
*
CT                    CALL LUCIAQEXIT('ADADS')
                      IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                      J = 0
                      L = 1
*
CT                    CALL LUCIAQENTER('MATCG')
                      DO  IJL = 1, NJL
                        CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                        I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                        JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                        IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                          ZERO = 0.0D0
                          ISETVECOPS(3) = ISETVECOPS(3) + NKBTC*NIBTC
                          CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                        ELSE 
                          CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,
     &                              NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                        END IF
                      END DO
CT                    CALL LUCIAQEXIT ('MATCG')
*
                      JLBOFF = JLBOFF + NJL
                    END DO 
*
*. ( End of loop over jlpair in batch )
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.and now ,to the work
                    LIKB = NIBTC*NKBTC
                    IF(NTEST.GE.3000) THEN
                     WRITE(6,*) ' Integral block '
                     CALL WRTMAT(XINT,NIKT,NJLT,NIKT,NJLT)
                    END IF
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' CSCR matrix '
                      CALL WRTMAT(CSCR,LIKB,NJLT,LIKB,NJLT)
                    END IF
*
C?                  MXACIJO = MXACIJ
                    MXACIJ = MAX(MXACIJ,LIKB*NJLT,LIKB*NIKT)
C?                  IF(MXACIJ.GT.MXACIJO) THEN
C?                    write(6,*) ' New max MXACIJ = ', MXACIJ
C?                    write(6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
C?                    WRITE(6,*) ' ITYP,JTYP,KTYP,LTYP',
C?   &                             ITYP,JTYP,KTYP,LTYP 
C?                    WRITE(6,*)'NIJT, NJLT, NIBTC NKBTC',
C?   &                           NIJT, NJLT,NIBTC,NKBTC
C?                  END IF
*
                    FACTORC = 0.0D0
                    FACTORAB = 1.0D0 
                    CALL MATML7(SSCR,CSCR,XINT,
     &                          LIKB,NIKT,LIKB,NJLT,NIKT,NJLT,
     &                          FACTORC,FACTORAB,2)
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' SSCR matrix '
                      CALL WRTMAT(SSCR,LIKB,NIKT,LIKB,NIKT)
                    END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                    II12 = 2
*
                    IONE = 1
                    IKBOFF = 1
                    DO IKPAIR = 1, IKBT(2,IKBTC)
                      ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                      KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                      NI = NOBPTS(ITYP,ISM)
                      NK = NOBPTS(KTYP,KSM)
                      IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                        NIK = NI*(NI+1)/2
                        IKSM = 1
                      ELSE
                        NIK = NI * NK
                        IKSM = 0
                      END IF
CT                    CALL LUCIAQENTER('ADADS')
                      IF(IFRST.EQ.1) KFRST = 1 
                      ONE = 1.0D0
*
                      IAC = I4_AC(1)
                      KAC = I4_AC(2)
*
                      CALL ADAADAST_GAS(IONE,ISM,ITYP,NI,IAC,
     &                                  IONE,KSM,KTYP,NK,KAC,
     &                                ISCTP,ISCSM,IGRP,
     &                                KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                                IFRST,KFRST,II12,K12,ONE         )
*
                      IFRST = 0
                      KFRST = 0
CT                    CALL LUCIAQEXIT ('ADADS')
*
CT                    CALL LUCIAQENTER('MATCS')
                      I = 0
                      K = 1
                      DO IK = 1, NIK
                        CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                        IKOFF = (K-1)*NI + I
                        ISBOFF = 1+(IKBOFF-1+IK-1)*NIBTC*NKBTC
                        IF(IKSM.EQ.1.AND.I.EQ.k) THEN
* a+ i a+i gives trivially zero
                        ELSE
                          CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                                NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                        END IF
                      END DO
CT                    CALL LUCIAQEXIT ('MATCS')
                      IKBOFF = IKBOFF + NIK
*
                    END DO
*                   ^ End of loop over IKPAIRS in batch
*
                  IF(KEND.EQ.0) GOTO 1800
*.                ^ End of loop over partitionings of resolution strings
 1801           CONTINUE
*               ^ End of loop over partionings of I strings
 1930         CONTINUE
*             ^ End of loop over batches of JL
 1940       CONTINUE
*           ^ End of loop over batches of IK
 1950     CONTINUE
*         ^ End of loop over IKOBSM
*
*
*      ==============================================
        ELSE IF(.NOT.( I4_AC(1).EQ. I4_AC(2)) ) THEN
*      ==============================================
*
*
* Three types of operators :
* a+ a  a+ a  
* a+ a  a  a+
* a  a+ a+ a 
*
* The first end up with 
* -a+ i ak a+l aj X2(ik,jl)
*
* Number two and three end up with
* -a i a k a l aj XC(ik,jl)  ( In coulomb form)
*
          JLSM = 0
          IKSM = 0
*. Symmetry of allowed Double excitation,loop over excitations
          DO 2950 IKOBSM = 1, NSMOB
            JLOBSM = SXDXSX(IKOBSM,IDXSM)
            IF(JLOBSM.EQ.0) GOTO 2950
*. types + symmetries defined => K strings are defined 
            KFRST = 1
            K2FRST = 1
            DO ISM = 1, NSMOB
              KSM = ADSXA(ISM,IKOBSM)
              DO JSM = 1, NSMOB
                LSM = ADSXA(JSM,JLOBSM)
                IF(NTEST.GE.2000) WRITE(6,*) ' ISM KSM LSM JSM',
     &          ISM,KSM,LSM,JSM
                ISCR(I4_REO(1)) = ISM
                ISCR(I4_REO(2)) = KSM
                ISCR(I4_REO(3)) = LSM
                ISCR(I4_REO(4)) = JSM
*
                ISM_ORIG = ISCR(1)             
                KSM_ORIG = ISCR(2)             
                LSM_ORIG = ISCR(3)             
                JSM_ORIG = ISCR(4)             
*
C           DO ISM_ORIG = 1, NSMOB
C             KSM_ORIG = ADSXA(ISM_ORIG,IKOBSM)
C             DO JSM_ORIG = 1, NSMOB
C               LSM_ORIG = ADSXA(JSM_ORIG,JLOBSM)
*
C               ISCR(1) = ISM_ORIG
C               ISCR(2) = KSM_ORIG
C               ISCR(3) = LSM_ORIG
C               ISCR(4) = JSM_ORIG
*
C               ISM = ISCR(I4_REO(1))             
C               KSM = ISCR(I4_REO(2))             
C               LSM = ISCR(I4_REO(3))             
C               JSM = ISCR(I4_REO(4))             
*
                NI = NOBPTS(ITYP,ISM)
                NJ = NOBPTS(JTYP,JSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
                NIK = NI*NK
                NJL = NJ*NL
                IF(NIK.EQ.0.OR.NJL.EQ.0) GOTO 2803
*
                ITPSM_ORIG = (ITYP_ORIG-1)*NSMOB + ISM_ORIG
                JTPSM_ORIG = (JTYP_ORIG-1)*NSMOB + JSM_ORIG
                KTPSM_ORIG = (KTYP_ORIG-1)*NSMOB + KSM_ORIG
                LTPSM_ORIG = (LTYP_ORIG-1)*NSMOB + LSM_ORIG
*
                IF(ITPSM_ORIG.GE.KTPSM_ORIG.AND.
     &             JTPSM_ORIG.GE.LTPSM_ORIG) THEN
*
*.Obtain two electron integrals as xint(ik,jl) = (ij!kl)-(il!kj)
                IKSM = 0
                JLSM = 0
                IF(I4_AC(1).EQ.I4_AC(3)) THEN
* a+ a a+ a
                  ICOUL = 2
                ELSE
* a+ a a a+ or a+ a a a+
                  ICOUL = 1
                END IF
*. Use coulomb - exchange or just coulomb integrals ?
                IF(ITPSM_ORIG.EQ.KTPSM_ORIG
     &            .AND.JTPSM_ORIG.EQ.LTPSM_ORIG)THEN
*. No use of exchange
                  IXCHNG = 0
                  FACX = -0.5D0
                ELSE IF(ITPSM_ORIG.NE.KTPSM_ORIG
     &            .OR.JTPSM_ORIG.NE.LTPSM_ORIG) THEN
*. Exchange used, combines two terms
                  IXCHNG = 1
                  FACX = -0.5D0
                END IF
                IF(ITPSM_ORIG.NE.KTPSM_ORIG
     &            .AND.JTPSM_ORIG.NE.LTPSM_ORIG)THEN
*. Exchange used, combines four terms
                  IXCHNG = 1
                  FACX = -1.0D0
                END IF
           IF( NTEST.GE.1000) WRITE(6,*)
     &   ' ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX',
     &     ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX
*. fetch integrals
* we want the operator in the form a+i ak a+l aj ((ij!lk)-(ik!lj))
                ONE = 1.0D0
                IF(ICOUL.EQ.2) THEN
*. Obtain X2(ik,lj) = (ij!lk)
                  CALL LUCIAGETINT(XINT,ITYP,ISM,JTYP,JSM,LTYP,LSM,
     &                        KTYP,KSM,IXCHNG,IKSM,JLSM,ICOUL,
     &                        ONE,ONE)
                ELSE IF (ICOUL.EQ.1) THEN
                  IF(I_USE_SIMTRH.EQ.0) THEN
                      CALL LUCIAGETINT(XINT,ITYP,ISM,KTYP,KSM,JTYP,
     &                        JSM,
     &                        LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL,
     &                        ONE,ONE)
                  ELSE
                   IF(I4_AC(1).EQ.2) THEN
*. a+i ak al a+j (ik|jl) 
*. obtain integrals (ik!jl) 
                    CALL LUCIAGETINT(XINT,ITYP,ISM,KTYP,KSM,JTYP,
     &                       JSM,
     &                        LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL,
     &                        ONE,ONE)
                   ELSE IF(I4_AC(1).EQ.1) THEN
*. a i a+k a+l a j (ki!lj) 
*. Obtain (ki!lj) and transpose first two and last two indeces 
                    CALL LUCIAGETINT(XINT,KTYP,KSM,ITYP,ISM,LTYP,
     &                      LSM,
     &                      JTYP,JSM,IXCHNG,IKSM,JLSM,ICOUL,
     &                      ONE,ONE)
C                            TRP_H2_BLK(XINT,I12_OR_34,NI,NJ,NK,NL,SCR)
                    CALL TRP_H2_BLK(XINT,46,NK,NI,NL,NJ,XINT2)
                   END IF
                  END IF
                END IF
*               ^ End if integrals should be fetched


*. Loop over batches of I strings
                NPART = NROW/MAXI
                IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
c               if (npart.lt.nproc) then
c                  npart = min(nrow,nproc)
c               endif
                IF(NTEST.GE.100)
     &          write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
C               DO 2801 IIPART = ga_nodeid()+1, NPART,nproc
                DO 2801 IIPART = 1, NPART
                  IBOT = 1+(IIPART-1)*MAXI
                  ITOP = MIN(IBOT+MAXI-1,NROW)
                  NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                  KBOT = 1- MAXK
                  KTOP = 0
 2800             CONTINUE
                    KBOT = KBOT + MAXK
                    KTOP = KTOP + MAXK
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
CT                  CALL LUCIAQENTER('ADADS')
                    II12 = 1
                    K12 = 1
                    IONE = 1
*. Creation / annihilation maps , conjugated of above
                    IF(I4_AC(4).EQ.1) THEN
                      JAC = 2
                    ELSE 
                      JAC = 1
                    END IF 
                    IF(I4_AC(3).EQ.1) THEN
                      LAC = 2
                    ELSE 
                      LAC = 1
                    END IF 
C                   KFRST = 1
                    CALL ADAADAST_GAS(IONE,JSM,JTYP,NJ,JAC,
     &                                IONE,LSM,LTYP,NL,LAC,
     &                          ICCTP,ICCSM,IGRP,
     &                          KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                          JFRST,KFRST,II12,K12,SCLFAC)
                    JFRST = 0
                    KFRST = 0
*
CT                  CALL LUCIAQEXIT('ADADS')
                    IF(NKBTC.EQ.0) GOTO 2801
*. Loop over jl in TS classes and gather
CT                  CALL LUCIAQENTER('MATCG')
                    J = 0
                    L = 1
                    DO  IJL = 1, NJL
                      CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                      I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                      JLOFF = (IJL-1)*NKBTC*NIBTC+1
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,
     &                         NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                    END DO
CT                  CALL LUCIAQEXIT ('MATCG')
*
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
                    IKSM = 0
                    JLSM = 0
*.and now ,to the work
                    LIKB = NIBTC*NKBTC
                    IF(NTEST.GE.3000) THEN
                    WRITE(6,'(A,8(1X,I3))') ' IJKL SM/TP = ',
     &              ISM,ITYP,JSM,JTYP,KSM,KTYP,LSM,LTYP
                     WRITE(6,*) ' Integral block '
                     CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
                    END IF
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' CSCR matrix '
                      CALL WRTMAT(CSCR,LIKB,NJL,LIKB,NJL)
                    END IF
*
C?                  MXACIJO = MXACIJ 
                    MXACIJ = MAX(MXACIJ,LIKB*NJL,LIKB*NIK)
C?                  IF(MXACIJ.GT.MXACIJO) THEN
C?                    write(6,*) ' New max MXACIJ = ', MXACIJ
C?                    write(6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
C?                    WRITE(6,*) ' ITYP,JTYP,KTYP,LTYP',
C?   &                             ITYP,JTYP,KTYP,LTYP 
C?                    WRITE(6,*)'NIJ NJL NIBTC NKBTC',
C?   &                           NIJ,NJL,NIBTC,NKBTC
C?                  END IF
*
                    FACTORC = 0.0D0
                    FACTORAB = FACX  
                    CALL MATML7(SSCR,CSCR,XINT,
     &                          LIKB,NIK,LIKB,NJL,NIK,NJL,
     &                          FACTORC,FACTORAB,2)
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' SSCR matrix '
                      CALL WRTMAT(SSCR,LIKB,NIK,LIKB,NIK)
                    END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                    II12 = 2
*
                    IONE = 1
CT                  CALL LUCIAQENTER('ADADS')
                    IF(IFRST.EQ.1) KFRST = 1 
                    ONE = 1.0D0
*
                    IAC = I4_AC(1)
                    KAC = I4_AC(2)
*
C                   KFRST = 1
                    CALL ADAADAST_GAS(IONE,ISM,ITYP,NI,IAC,
     &                                IONE,KSM,KTYP,NK,KAC,
     &                              ISCTP,ISCSM,IGRP,
     &                              KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                              IFRST,KFRST,II12,K12,ONE          )
*
                    IFRST = 0
                    KFRST = 0
CT                  CALL LUCIAQEXIT ('ADADS')
*
CT                  CALL LUCIAQENTER('MATCS')
                    I = 0
                    K = 1
                    DO IK = 1, NIK
                      CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                      IKOFF = (K-1)*NI + I
                      ISBOFF = 1+(IK-1)*NIBTC*NKBTC
                      CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                     NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                    END DO
C                   write(6,*) ' first element of updated SB', SB(1)
CT                  CALL LUCIAQEXIT ('MATCS')
*
                  IF(KEND.EQ.0) GOTO 2800
*. End of loop over partitionings of resolution strings
 2801           CONTINUE
*               ^ End of loop over batches of I strings
              END IF
*             ^ End of if I. ge. K, J.ge. L
 2803         CONTINUE
              END DO
*             ^ End of loop over KSM
            END DO
*           ^ End of loop over ISM
 2950     CONTINUE
        END IF
*       ^ End of a+ a+ a a/a a a+ a+ versus a+ a a+ a switch

        IF(NIJKL1.EQ.0) CALL LUCIAQEXIT ('BB2A0')
        IF(NIJKL1.EQ.1) CALL LUCIAQEXIT ('BB2A1')
        IF(NIJKL1.EQ.2) CALL LUCIAQEXIT ('BB2A2')
        IF(NIJKL1.EQ.3) CALL LUCIAQEXIT ('BB2A3')
        IF(NIJKL1.EQ.4) CALL LUCIAQEXIT ('BB2A4')
 2000 CONTINUE
*
 2001 CONTINUE
*
C?      WRITE(6,*) ' Memcheck at end of RSBB2A '
C?      CALL MEMCHK 
C?      WRITE(6,*) ' Memcheck passed '
*
      CALL LUCIAQEXIT('RS2A ')
      RETURN
      END
      SUBROUTINE GETINCN2(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,IJKLOF,NSMOB,I2INDX,
     &                  ICOUL,INO12SYM,CFAC,EFAC) 
*
* Obtain integrals 
*
*     ICOUL = 0 : 
*                  XINT(IK,JL) = CFAC*(IJ!KL)         for IXCHNG = 0
*                              = CFAC*(IJ!KL)-EFAC*(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 : 
*                  XINT(IJ,KL) = CFAC*(IJ!KL)         for IXCHNG = 0
*                              = CFAC*(IJ!KL)-EFAC*(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 2 :  XINT(IL,JK) = CFAC*(IJ!KL)         for IXCHNG = 0
*                              = CFAC*(IJ!KL)-EFAC*(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0 
* 
*
* Version for integrals stored in INTLST
*
* If type equals zero, all integrals of given symmetry are fetched 
* ( added aug8, 98)
* Modified July 2010:
* type = -1 => all orbitals of given symmetry
* type =  0 => all inactive orbitals of given symmetry
* type = 1-ngas: all orbitals of a given gas and symmetry
* type = ngas + 1: all secondary orbitals of given symmetry
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLof(NsmOB,NsmOb,NsmOB)
*. Pair of orbital indeces ( symmetry ordered ) => address in symmetry packed 
*. matrix
      Dimension I2INDX(*)
*.Output
      DIMENSION XINT(*)
*. Local scratch      
      DIMENSION IJARR(MXPORB)
*
      IF(ITP.EQ.-1) THEN
        IORB = NTOOBS(ISM)
      ELSE IF(ITP.EQ.0) THEN
        IORB = NINOBS(ISM)
      ELSE IF(ITP.LE.NGAS) THEN
        iOrb=NOBPTS(ITP,ISM)
      ELSE
        IORB = NSCOBS(ISM)
      END IF
*
      IF(JTP.EQ.-1) THEN
        JORB = NTOOBS(JSM)
      ELSE IF(JTP.EQ.0) THEN
        JORB = NINOBS(JSM)
      ELSE IF(JTP.LE.NGAS) THEN
        JOrb=NOBPTS(JTP,JSM)
      ELSE
        JORB = NSCOBS(JSM)
      END IF
*
      IF(KTP.EQ.-1) THEN
        KORB = NTOOBS(KSM)
      ELSE IF(KTP.EQ.0) THEN
        KORB = NINOBS(KSM)
      ELSE IF(KTP.LE.NGAS) THEN
        KOrb=NOBPTS(KTP,KSM)
      ELSE
        KORB = NSCOBS(KSM)
      END IF
*
      IF(LTP.EQ.-1) THEN
        LORB = NTOOBS(LSM)
      ELSE IF(LTP.EQ.0) THEN
        LORB = NINOBS(LSM)
      ELSE IF(LTP.LE.NGAS) THEN
        LOrb=NOBPTS(LTP,LSM)
      ELSE
        LORB = NSCOBS(LSM)
      END IF
*
*. Offsets relative to start of orbitals, symmetry ordered 
*
      IF(ITP.LE.0) THEN
        IOFF = IBSO(ISM)
      ELSE
        IOFF = IBSO(ISM)+NINOBS(ISM)
        DO IITP = 1, ITP -1
          IOFF = IOFF + NOBPTS(IITP,ISM)
        END DO
      END IF
*
      IF(JTP.LE.0) THEN
        JOFF = IBSO(JSM)
      ELSE
        JOFF = IBSO(JSM)+NINOBS(JSM)
        DO JJTP = 1, JTP -1
          JOFF = JOFF + NOBPTS(JJTP,JSM)
        END DO
      END IF
*
      IF(KTP.LE.0) THEN
        KOFF = IBSO(KSM)
      ELSE
        KOFF = IBSO(KSM)+NINOBS(KSM)
        DO KKTP = 1, KTP -1
          KOFF = KOFF + NOBPTS(KKTP,KSM)
        END DO
      END IF
*
*
      IF(LTP.LE.0) THEN
        LOFF = IBSO(LSM)
      ELSE
        LOFF = IBSO(LSM)+NINOBS(LSM)
        DO LLTP = 1, LTP -1
          LOFF = LOFF + NOBPTS(LLTP,LSM)
        END DO
      END IF
*
*     Collect Coulomb terms
*
      ijblk = max(ism,jsm)*(max(ism,jsm)-1)/2 + min(ism,jsm)
      klblk = max(ksm,lsm)*(max(ksm,lsm)-1)/2 + min(ksm,lsm)
*
      IF (INO12SYM.EQ.1) THEN
       IJRELKL = 1
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF(IJBLK.GT.KLBLK) THEN
       IJRELKL = 1
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.EQ.KLBLK) THEN
       IJRELKL = 0
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.LT.KLBLK) THEN
       IJRELKL = -1
       IBLOFF = IJKLOF(MAX(KSM,LSM),MIN(KSM,LSM),MAX(ISM,JSM))
      END IF
*
      itOrb=NTOOBS(iSm)
      jtOrb=NTOOBS(jSm)
      ktOrb=NTOOBS(kSm)
      ltOrb=NTOOBS(lSm)
*
      If(ISM.EQ.JSM) THEN
       IJPAIRS = ITORB*(ITORB+1)/2
      ELSE
       IJPAIRS = ITORB*JTORB
      END IF
*
      IF(KSM.EQ.LSM) THEN
        KLPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KLPAIRS = KTORB*LTORB
      END IF
*
      iInt=0
      Do lJeppe=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lJeppe
        Do jJeppe=jMin,jOff+jOrb-1
*
*
*. Set up array IJ*(IJ-1)/2 
          IF(IJRELKL.EQ.0) THEN 
            DO II = IOFF,IOFF+IORB-1
              IJ = I2INDX((JJEPPE-1)*NTOOB+II)
              IJARR(II) = IJ*(IJ-1)/2
            END DO
          END IF
*
          Do kJeppe=kOff,kOff+kOrb-1
            iMin = iOff
            kl = I2INDX(KJEPPE+(LJEPPE-1)*NTOOB)
            If(IKSM.ne.0) iMin = kJeppe
            IF(ICOUL.EQ.1)  THEN  
*. Address before integral (1,j!k,l)
                IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &               + (KJEPPE-KOFF)*Jorb*Iorb
     &               + (JJEPPE-JOFF)*Iorb
            ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK) 
                IINT = (KJEPPE-KOFF)*JORB*LORB*IORB
     &               + (JJEPPE-JOFF)     *LORB*IORB
     &               + (LJEPPE-LOFF)          *IORB
            END IF
*
            IF(IJRELKL.EQ.1) THEN
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) > (ksm,lsm)
              IJKL0 = IBLOFF-1+(kl-1)*ijPairs
              IJ0 = (JJEPPE-1)*NTOOB         
              Do iJeppe=iMin,iOff+iOrb-1
                  ijkl = ijkl0 + I2INDX(IJEPPE+IJ0)
                  iInt=iInt+1
                  Xint(iInt) = CFAC*Intlst(ijkl)
              End Do
            END IF
*
*. block (ISM JSM !ISM JSM)
            IF(IJRELKL.EQ.0) THEN 
              IJ0 = (JJEPPE-1)*NTOOB         
              KLOFF = KL*(KL-1)/2
              IJKL0 = (KL-1)*IJPAIRS-KLOFF
              Do iJeppe=iMin,iOff+iOrb-1
                ij = I2INDX(IJEPPE+IJ0   )
                If ( ij.ge.kl ) Then
C                 ijkl=ij+(kl-1)*ijPairs-klOff
                  IJKL = IJKL0 + IJ
                Else
                  IJOFF = IJARR(IJEPPE)
                  ijkl=kl+(ij-1)*klPairs-ijOff
                End If
                iInt=iInt+1
                Xint(iInt) = CFAC*Intlst(iblOff-1+ijkl)
              End Do
            END IF
*
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) < (ksm,lsm)
            IF(IJRELKL.EQ.-1) THEN 
              ijkl0 = IBLOFF-1+KL - KLPAIRS
              IJ0 = (JJEPPE-1)*NTOOB         
              Do iJeppe=iMin,iOff+iOrb-1
                IJKL = IJKL0 + I2INDX(IJEPPE + IJ0)*KLPAIRS
                iInt=iInt+1
                Xint(iInt) = CFAC*Intlst(ijkl)
              End Do
            END IF
*
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
*
      IF(ISM.EQ.LSM) THEN
       ILPAIRS = ITORB*(ITORB+1)/2
      ELSE
       ILPAIRS = ITORB*LTORB
      END IF
*
      IF(KSM.EQ.JSM) THEN
        KJPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KJPAIRS = KTORB*JTORB
      END IF
*
        ilblk = max(ism,lsm)*(max(ism,lsm)-1)/2 + min(ism,lsm)
        kjblk = max(ksm,jsm)*(max(ksm,jsm)-1)/2 + min(ksm,jsm)
        IF(ILBLK.GT.KJBLK) THEN
          ILRELKJ = 1
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.EQ.KJBLK) THEN
          ILRELKJ = 0
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.LT.KJBLK) THEN
          ILRELKJ = -1
          IBLOFF = IJKLOF(MAX(KSM,JSM),MIN(KSM,JSM),MAX(ISM,LSM))
        END IF
*
        iInt=0
        Do lJeppe=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=lJeppe
*
          IF(ILRELKJ.EQ.0) THEN
           DO II = IOFF,IOFF+IORB-1
             IL = I2INDX(II+(LJEPPE-1)*NTOOB)
             IJARR(II) = IL*(IL-1)/2
           END DO
          END IF
*
          Do jJeppe=jMin,jOff+jOrb-1
            Do kJeppe=kOff,kOff+kOrb-1
              KJ = I2INDX(KJEPPE+(JJEPPE-1)*NTOOB)
              KJOFF = KJ*(KJ-1)/2
              iMin = iOff
*
              IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                  IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &                  + (KJEPPE-KOFF)*Jorb*Iorb
     &                  + (JJEPPE-JOFF)*Iorb
              ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK) 
                IINT = (KJEPPE-KOFF)*JORB*LORB*IORB
     &               + (JJEPPE-JOFF)     *LORB*IORB
     &               + (LJEPPE-LOFF)          *IORB
              END IF
*
              If(IKSM.ne.0) iMin = kJeppe
*
              IF(ILRELKJ.EQ.1) THEN
                ILKJ0 = IBLOFF-1+( kj-1)*ilpairs
                IL0 = (LJEPPE-1)*NTOOB 
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE + IL0)
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-EFAC*Intlst(ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.0) THEN
                IL0 = (LJEPPE-1)*NTOOB 
                ILKJ0 = (kj-1)*ilPairs-kjOff
                Do iJeppe=iMin,iOff+iOrb-1
                  IL = I2INDX(IJEPPE + IL0 )
                  If ( il.ge.kj ) Then
C                     ilkj=il+(kj-1)*ilPairs-kjOff
                      ILKJ = IL + ILKJ0
                    Else
                      ILOFF = IJARR(IJEPPE)
                      ilkj=kj+(il-1)*kjPairs-ilOff
                    End If
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-EFAC*Intlst(iBLoff-1+ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.-1) THEN
                ILKJ0 = IBLOFF-1+KJ-KJPAIRS
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE+ IL0)*KJPAIRS
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-EFAC*Intlst(ilkj)
                End Do
              END IF
*
            End Do
          End Do
        End Do
      End If
*
      RETURN
      END
      SUBROUTINE ALG_ROUTERX(ILOC,IROC,NOP,IOP_TYP,IOP_AC,IOP_REO,
     &                      SIGN)
*
* Decide route for calculating <ILOC! operator string !ROC>
* with smallest amount of operations/storage
*
* The operator string contains NOP elementary operators
* defined by
*
* IOP_TYP : Orbital type of each operator
* IOP_AC  : Creation/annihilation operator
*
*. Output :
* ==========
*
* IOP_REO : New to old order
* sign    : sign of full rank operator
*
* Method : use IPHGAS , and move creation of holes 
*                       and annihilation of particles to the left
*
* Jeppe Olsen, October 1997
* version  of : dec 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Specific input, I*OC : Contains number of electrons in each space
*  (not used !!!)
       INTEGER ILOC(NGAS),IROC(NGAS)
       INTEGER IOP_TYP(*), IOP_AC(*)
*. Output, new to old order :
       INTEGER IOP_REO(*)
*
*. The loops are assumed to be of the resolution type, so
*. we will insert a resolution of the identity after NOP/2
*  operators(we assume that NOP is even ) :
*
* <ILOC! IOPL !K><K! IOPR ! IROC >
*=<ILOC! IOPL !K><IROC! (IOPR)+ ! K >
*
*
* and the question is now how to order the operators.
* The operation count can in general be mininized by
*. According to above : Should operators stand left or
*                       right of the resolution K
*
      NTEST = 00
*
      IF(IUSE_PH.EQ.0) THEN 
*. No modificatione 
        DO IIIOP = 1, NOP
          IOP_REO(IIIOP) = IIIOP
        END DO
        SIGN = 1.0D0
        GOTO 1001
      END IF
*
*. First time around : not very general
      NMOVE = 0
      IF(NOP.EQ.2) THEN
*
* One-Electron operator 
*
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.1) THEN
* a+1 a 2 : Move around only if both spaces are hole spaces
          NMOVE = 0
          IF(IPHGAS(IOP_TYP(1)).EQ.2) NMOVE = 1
          IF(IPHGAS(IOP_TYP(2)).EQ.2) NMOVE = NMOVE + 1
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) '  IOP_TYP(1), IOP_TYP(2), NMOVE ',
     &                    IOP_TYP(1), IOP_TYP(2), NMOVE
          END IF
        ELSE
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP 
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          STOP ' ALG_ROUTER : Path not implemented '
        END IF
*
        IF(NMOVE.EQ.2) THEN
*. Move the operators around in expansion
          IOP_REO(1) = 2
          IOP_REO(2) = 1
          SIGN = -1.0D0
        ELSE
*. No reorganization adviced
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          SIGN = 1.0D0
        END IF
      ELSE IF(NOP.EQ.4) THEN
*
* Two-electron operator
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.2.AND.
     &     IOP_AC(3).EQ.1.AND.IOP_AC(4).EQ.1) THEN
* a+1 a+2 a3 a4
          IF(IPHGAS(IOP_TYP(1)).EQ.2) THEN
            IMOVE1 = 1
          ELSE
            IMOVE1 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(2)).EQ.2) THEN
            IMOVE2 = 1
          ELSE
            IMOVE2 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(3)).EQ.2) THEN
            IMOVE3 = 1
          ELSE
            IMOVE3 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(4)).EQ.2) THEN
            IMOVE4 = 1
          ELSE
            IMOVE4 = 0
          END IF
        ELSE
*. Not yet implemented 
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP 
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          STOP ' ALG_ROUTER : Path not implemented '
        END IF
*. Number of left and right operators that would like to be moved
        IMOVEL = IMOVE1 + IMOVE2
        IMOVER = IMOVE3 + IMOVE4
        IMOVE = MIN(IMOVEL,IMOVER)
        IF(IMOVE.EQ.2) THEN
*. Well, we are all moving
C         IOP_REO(1) = 3
C         IOP_REO(2) = 4
C         IOP_REO(3) = 1
C         IOP_REO(4) = 2
          IOP_REO(1) = 4
          IOP_REO(2) = 3
          IOP_REO(3) = 2
          IOP_REO(4) = 1
          SIGN = 1.0D0
        ELSE IF(IMOVE.EQ.0) THEN
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
        ELSE IF(IMOVEL.EQ.1.AND.IMOVER.EQ.1) THEN
* a+ a a+ a
* 1 : Position the two operators to be moved as operators 2, 3
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
          IF(IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -SIGN 
          END IF
          IF(IMOVE4.EQ.1) THEN
           IOP_REO(3) = 4
           IOP_REO(4) = 3
           SIGN = - SIGN
          END IF
* 2 : and interchange operators 2 and 3
          ISAVE = IOP_REO(2)                       
          IOP_REO(2) = IOP_REO(3)
          IOP_REO(3) = ISAVE
          SIGN = - SIGN
        ELSE IF (IMOVEL.EQ.2.AND.IMOVER.EQ.1) THEN
* Final operator is a a+ a+ a
          IF(IMOVE4.EQ.1) THEN
            IOP_REO(3) = 4
            IOP_REO(4) = 3
            SIGN = -1.0D0
          ELSE
            IOP_REO(3) = 3
            IOP_REO(4) = 4
            SIGN = 1.0D0
          END IF
          IOP_REO(1) = IOP_REO(3)
          IOP_REO(2) = 1   
          IOP_REO(3) = 2
          SIGN = SIGN
        ELSE IF( IMOVEL.EQ.1 .AND. IMOVER .EQ. 2 ) THEN
* Final operator is a+ a a a+
          IF( IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -1.0D0
          ELSE
            IOP_REO(1) = 1
            IOP_REO(2) = 2
            SIGN = 1.0D0
          END IF
          IOP_REO(4) = IOP_REO(2)
          IOP_REO(2) = 3  
          IOP_REO(3) = 4
          SIGN = SIGN
        END IF
*       ^ End of number of pairs to be moved around
      ELSE 
*. Not one- or two- electron operator
        WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
        WRITE(6,*) ' NOP  = ', NOP 
        WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
        STOP ' ALG_ROUTER : Path not implemented '
      END IF
 1001 CONTINUE
*

      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Information from ALG_ROUTER '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) 
        WRITE(6,*) 
     &  ' Input : anni- or crea- operators + orbital types'
        CALL IWRTMA(IOP_AC,1,NOP,1,NOP)
        CALL IWRTMA(IOP_TYP,1,NOP,1,NOP)
        WRITE(6,*)
        WRITE(6,*) ' Suggested output order '
        CALL IWRTMA(IOP_REO,1,NOP,1,NOP)
        WRITE(6,*)  ' Sign = ', Sign
      END IF
*
      RETURN
      END
      SUBROUTINE ADAADAST_GAS(IOB,IOBSM,IOBTP,NIOB,IAC,
     &                       JOB,JOBSM,JOBTP,NJOB,JAC,
     &                       ISPGP,ISM,ITP,KMIN,KMAX,
     &                       I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12,
     &                       SCLFAC)
*
* Obtain two-operator mappings
* a+/a IORB a+/a JORB !KSTR> = +/-!ISTR>
*
* Whether creation- or annihilation operators are in use depends 
* upon IAC, JAC : 1=> Annihilation, 
*                 2=> Creation
*
* In the form
* I1(KSTR) =  ISTR if a+/a IORB a+/a JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,IGRP.
* (numbering relative to TS start)
*. Only excitations IOB. GE. JOB are included 
* The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
* and IOB (JOB) is the first orbital to be used, and the number of orbitals
* to be checked is NIOB ( NJOB).
*
* Only orbital pairs IOB .gt. JOB are included (if the types are identical)
*
* The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , August of 95   ( adadst)
*               November 1997 : annihilation added
*
*
* ======
*. Input
* ======
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
*./BIGGY
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*./ORBINP/
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
*. Local scratch
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/SSAVE/NELIS(4), NSTRKS(4)
      COMMON/UMMAGUMMA/NSTRIA(4)
      INTEGER KELFGRP(MXPNGAS),KGRP(MXPNGAS)
      COMMON/COMJEP/MXACJ,MXACIJ,MXAADST
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADAADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' IOB,IOBSM,IOBTP,IAC ', IOB,IOBSM,IOBTP,IAC
        WRITE(6,'(A,4I3)') ' JOB,JOBSM,JOBTP,JAC ', JOB,JOBSM,JOBTP,JAC
        WRITE(6,'(A,2I3)') ' I12, K12 ', I12, K12
        WRITE(6,'(A,2I3)') ' IFRST,KFRST', IFRST,KFRST
      END IF
*
*
*. Internal affairs
*
      IF(I12.LE.4.AND.K12.LE.1) THEN
        KLLOC = KLOCSTR(K12)
        KLLZ = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*) 
     &  ' ADST_GAS : Illegal value of I12 or K12 ', I12, K12
        STOP' ADST_GAS : Illegal value of I12 or K12  '
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' KLLOC KLLREO',KLLOC,KLLREO
      END IF

*
*. Supergroup and symmetry of K strings
*
      CALL SYMCOM(2,0,IOBSM,K1SM,ISM)
      CALL SYMCOM(2,0,JOBSM,KSM,K1SM)
      IF(NTEST.GE.100) WRITE(6,*) ' K1SM,KSM : ',  K1SM,KSM
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      IF(IAC.EQ.1) THEN
        IACADJ = 2
        IDELTA =-1
      ELSE IF(IAC.EQ.2) THEN
        IACADJ = 1
        IDELTA = 1
      END IF
      IF(JAC.EQ.1) THEN
        JACADJ = 2
        JDELTA =-1
      ELSE IF(JAC.EQ.2) THEN
        JACADJ = 1
        JDELTA = 1
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' IACADJ, JACADJ', IACADJ,JACADJ
       WRITE(6,*) ' IDELTA, JDELTA', IDELTA, JDELTA
      END IF
*. Occupation of K-strings
      IF(IOBTP.EQ.JOBTP) THEN
        IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA-JDELTA
        JEL = IEL
      ELSE
        IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA
        JEL = NELFSPGP(JOBTP,ISPGPABS)-JDELTA
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' IEL, JEL', IEL,JEL
*. Trivial zero ? (Nice, then mission is complete )
      ITRIVIAL = 0
      IF(IEL.LT.0.OR.JEL.LT.0.OR.
     &   IEL.GT.NOBPT(IOBTP).OR.JEL.GT.NOBPT(JOBTP)) THEN
*. No strings with this number of elecs - be happy : No work 
        NK = 0
        KACT = 0
        KACGRP = 0
        IF(NTEST.GE.100) WRITE(6,*) ' Trivial zero excitations'
        ITRIVIAL = 1
C       GOTO 9999
      ELSE
*. Find group with IEL electrons in IOBTP, JEL in JOBTP
        IIGRP = 0
        DO IGRP = IBGPSTR(IOBTP),IBGPSTR(IOBTP)+NGPSTR(IOBTP)-1
          IF(NELFGP(IGRP).EQ.IEL) IIGRP = IGRP
        END DO
        JJGRP = 0
        DO JGRP = IBGPSTR(JOBTP),IBGPSTR(JOBTP)+NGPSTR(JOBTP)-1
          IF(NELFGP(JGRP).EQ.JEL) JJGRP = JGRP
        END DO
C?      WRITE(6,*) ' ADAADA : IIGRP, JJGRP', IIGRP,JJGRP
*
        IF(IIGRP.EQ.0.OR.JJGRP.EQ.0) THEN
          WRITE(6,*)' ADAADAST : cul de sac, active K groups not found'
          WRITE(6,*)' Active GAS spaces  ' ,IOBTP, JOBTP
          WRITE(6,*)' Number of electrons', IEL, JEL
          STOP      ' ADAADAST : cul de sac, active K groups not found'
        END IF
*
      END IF
*. Groups defining Kstrings
      IF(ITRIVIAL.NE.1) THEN
      CALL ICOPVE(ISPGPFTP(1,ISPGPABS),KGRP,NGAS)
      KGRP(IOBTP) = IIGRP
      KGRP(JOBTP) = JJGRP
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Groups in KGRP '
        CALL IWRTMA(KGRP,1,NGAS,1,NGAS)
      END IF
      END IF
*
* In ADADS1_GAS we need : Occupation of KSTRINGS
*                         lexical => Actual order for I strings
* Generate if required
*
      IF(IFRST.NE.0) THEN
*.. Generate information about I strings
*. Arc weights for ISPGP
        NTEST2 = NTEST
        CALL WEIGHT_SPGP(int_mb(KLLZ),NGAS,
     &                  NELFSPGP(1,ISPGPABS),
     &                  NOBPT,int_mb(KLZSCR),NTEST2)
        NELI = NELFTP(ITP)
        NELIS(I12) = NELI
*. Reorder array for I strings
        CALL GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,
     &                         int_mb(KLLOC),NOCOB,
     &                         1,int_mb(KLLZ),int_mb(KLLREO))
        IF(NTEST.GE.1000) THEN
         write(6,*) ' Info on I strings generated '
         write(6,*) ' NSTRI = ', NSTRI
         WRITE(6,*) ' REORDER array '
         CALL IWRTMA(int_mb(KLLREO),1,NSTRI,1,NSTRI)
       END IF
       NSTRIA(I12) = NSTRI
*
      END IF
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' REORDER array for I STRINGS'
       CALL IWRTMA(int_mb(KLLREO),1,NSTRI,1,NSTRI)
      END IF
*
      IF(ITRIVIAL.EQ.1) GOTO 9999
      NELK = NELIS(I12)
      IF(IAC.EQ.1) THEN
        NELK = NELK + 1
      ELSE
        NELK = NELK - 1
      END IF
      IF(JAC.EQ.1) THEN
        NELK = NELK + 1
      ELSE
        NELK = NELK - 1
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' NELK = ' , NELK
      IF(KFRST.NE.0) THEN
*. Generate occupation of K STRINGS

       CALL GETSTR2_TOTSM_SPGP(KGRP,NGAS,KSM,NELK,NSTRK,
     &                        int_mb(KLLOC),NOCOB,
     &                        0,ISUM,IDUM)
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
       NSTRKS(K12) = NSTRK
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' K strings generated '
         WRITE(6,*) ' Reorder array after generation of K strings'
         CALL IWRTMA(int_mb(KLLREO),1,NSTRI,1,NSTRI)
       END IF
      END IF
*
      NSTRK = NSTRKS(K12)
*
      IIOB = IOBPTS(IOBTP,IOBSM) + IOB - 1
      JJOB = IOBPTS(JOBTP,JOBSM) + JOB - 1
*
      IZERO = 0
      ZERO = 0.0D0
      CALL ISETVC(I1  ,IZERO,LI1*NIOB*NJOB)
COLD  CALL SETVEC(XI1S,ZERO ,LI1*NIOB*NJOB)
*
      CALL ADAADAS1_GAS(NK,I1,XI1S,LI1,
     &          IIOB,NIOB,IAC,JJOB,NJOB,JAC,
     &          int_mb(KLLOC),NELK,NSTRK,int_mb(KLLREO),int_mb(KLLZ),
     &          NOCOB,KMAX,KMIN,IEND,SCLFAC,NSTRIA(I12))
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Reorder array after ADAADAS1'                  
         CALL IWRTMA(int_mb(KLLREO),1,NSTRI,1,NSTRI)
       END IF
 9999 CONTINUE
*
C     WRITE(6,*) ' Memcheck at end of ADAADAS1 '
C     CALL MEMCHK
      RETURN
      END
      SUBROUTINE ADAADAS1_GAS(NK,I1,XI1S,LI1,
     &                IORB,NIORB,IAC,JORB,NJORB,JAC,
     &                KSTR,NKEL,NKSTR,IREO,IZ,
     &                NOCOB,KMAX,KMIN,IEND,SCLFAC,NSTRI)
*
* Obtain I1(KSTR) = +/- a+/a  IORB a+/a JORB !KSTR>
* Only orbital pairs IOB .gt. JOB are included 
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : First I orbital to be added
* NIORB : Number of orbitals to be added : IORB to IORB-1+NIORB
*        are used. They must all be in the same TS group
* JORB : First J orbital to be added 
* LORB : Number of orbitals to be added : JORB to JORB-1+NJORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 =>  a+IORB a+JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
* L.R. Jan 20, 1998
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*. Local scratch, atmost 1000 orbitals in a given TS block)
      DIMENSION ISCR(1000)
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' ADADS1_GAS speaking '
       WRITE(6,*) ' ==================== '
       WRITE(6,'(A,3I4)') ' IORB,NIORB,IAC ', IORB,NIORB,IAC   
       WRITE(6,'(A,3I4)') ' JORB,NJORB,JAC ', JORB,NJORB,JAC   
*
C      WRITE(6,*) ' Kstrings in action (el,string) '
C      WRITE(6,*) ' ==============================='
C      CALL IWRTMA(KSTR,NKEL,NKSTR,NKEL,NKSTR)
       WRITE(6,*) ' 24 elements of reorder array'
       CALL IWRTMA(IREO,1,24,1,24)
*
      END IF
*
      IORBMIN = IORB
      IORBMAX = IORB + NIORB - 1
*
      JORBMIN = JORB
      JORBMAX = JORB + NJORB - 1
*
      NIJ = NIORB*NJORB
*
      KEND = MIN(NKSTR,KMAX)
      IF(KEND.LT.NKSTR) THEN
        IEND = 0
      ELSE
        IEND = 1
      END IF
      NK = KEND-KMIN+1
*
      IF(IAC.EQ.2.AND.JAC.EQ.2) THEN
*
* ==========================
* Creation- creation mapping
* ==========================
*
        DO KKSTR = KMIN,KEND 
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Occupation of string ', KKSTR
           CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
         END IF
*. Loop over electrons after which JORB can be added 
         DO JEL = 0, NKEL
*
           IF(JEL.EQ.0 ) THEN     
             JORB1 = JORBMIN - 1
           ELSE
             JORB1 = MAX(JORBMIN-1,KSTR(JEL,KKSTR))
           END IF
           IF(JEL.EQ.NKEL) THEN
             JORB2 = JORBMAX + 1
           ELSE
             JORB2 = MIN(JORBMAX+1,KSTR(JEL+1,KKSTR))
           END IF
           IF(NTEST.GE.1000)
     &      WRITE(6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
*
           IF(JEL.GT.0.AND.JORB1.GE.JORBMIN.AND.
     &                     JORB1.LE.JORBMAX) THEN
*. vanishing for any IORB
             IJOFF = (JORB1-JORBMIN)*NIORB 
             DO IIORB = 1, NIORB
               IJ = IJOFF + IIORB
               I1(KKSTR-KMIN+1,IJ) = 0   
               XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
             END DO
           END IF
*
           IF(JORB1.LT.JORBMAX.AND.JORB2.GT.JORBMIN) THEN
*. Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
             SIGNJ = (-1) ** JEL * SCLFAC
*. reverse lexical number of the first JEL ELECTRONS
             ILEX0 = 1
             DO JJEL = 1, JEL  
               ILEX0 = ILEX0 + IZ(KSTR(JJEL,KKSTR),JJEL)
             END DO
             DO JJORB = JORB1+1, JORB2-1
* And electron JEL + 1
               ILEX1 = ILEX0 + IZ(JJORB,JEL+1)
*. Add electron IORB
               DO IEL = JEL, NKEL
                 IF(IEL.EQ.JEL) THEN
                   IORB1 = MAX(JJORB,IORBMIN-1)
                 ELSE
                   IORB1 = MAX(IORBMIN-1,KSTR(IEL,KKSTR))
                 END IF
                   IF(IEL.EQ.NKEL) THEN
                   IORB2 = IORBMAX+1
                 ELSE 
                   IORB2 = MIN(IORBMAX+1,KSTR(IEL+1,KKSTR))
                 END IF
                 IF(NTEST.GE.5000)
     &            WRITE(6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
                 IF(IEL.GT.JEL.AND.IORB1.GE.IORBMIN.AND.
     &                             IORB1.LE.IORBMAX) THEN
                   IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
                   I1(KKSTR-KMIN+1,IJ) = 0
                   XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
                 END IF
                 IF(IORB1.LT.IORBMAX.AND.IORB2.GT.IORBMIN) THEN
*. Orbitals IORB1+1 - IORB2 -1 can be added after ELECTRON IEL in KSTR
*. Reverse lexical number of the first IEL+1 electrons
                   ILEX2 = ILEX1
                   DO IIEL = JEL+1,IEL
                     ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+1)
                   END DO
*. add terms for the last electrons
                   DO IIEL = IEL+1,NKEL
                     ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+2)
                   END DO
                   IJOFF = (JJORB-JORBMIN)*NIORB 
                   SIGNIJ =  SIGNJ*(-1.0D0) ** (IEL+1)
                   DO IIORB = IORB1+1, IORB2-1
                     IJ = IJOFF + IIORB - IORBMIN + 1
                     ILEX = ILEX2 + IZ(IIORB,IEL+2)
                     IACT = IREO(ILEX)
                     IF(NTEST.GE.1000) THEN 
                       WRITE(6,*) 'IIORB JJORB', IIORB,JJORB
                       WRITE(6,*) ' ILEX IACT ', ILEX,IACT
                     END IF
                     I1(KKSTR-KMIN+1,IJ) = IACT
                     XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
                   END DO
                 END IF
               END DO
             END DO
           END IF
         END DO
        END DO
      ELSE IF(IAC.EQ.1.AND.JAC.EQ.1) THEN
*
* ===========================================
* annihilation - annihilation mapping (a i a j)
* ===========================================
*
        DO KKSTR = KMIN,KEND 
*. Active range for electrons
         IIELMIN = 0
         IIELMAX = 0 
         JJELMIN = 0
         JJELMAX = 0 
*
         DO KEL = 1, NKEL
          KKORB = KSTR(KEL,KKSTR)
          IF(IIELMIN.EQ.0.AND.KKORB.GE.IORBMIN)IIELMIN = KEL
          IF(JJELMIN.EQ.0.AND.KKORB.GE.JORBMIN)JJELMIN = KEL
          IF(KKORB.LE.IORBMAX) IIELMAX = KEL
          IF(KKORB.LE.JORBMAX) JJELMAX = KEL
         END DO
         IF(IIELMIN.EQ.0) IIELMIN = NKEL  + 1
         IF(JJELMIN.EQ.0) JJELMIN = NKEL  + 1

         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Occupation of string ', KKSTR
           CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
         END IF
*. Loop over first electron to be removed                   
C        DO JEL = 1, NKEL
         DO JEL = JJELMIN,JJELMAX
           JJORB = KSTR(JEL,KKSTR)
*. Loop over second electron to be removed
C          DO IEL = JEL+1, NKEL
           DO IEL = MAX(JEL+1,IIELMIN),IIELMAX
             IIORB = KSTR(IEL,KKSTR)
             IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' IEL JEL IORB JORB ',
     &        IEL,JEL,IORB,JORB
             END IF
             IF(IIORB.GE.IORBMIN.AND.IIORB.LE.IORBMAX.AND.
     &          JJORB.GE.JORBMIN.AND.JJORB.LE.JORBMAX     )THEN
               SIGN = (-1) ** (IEL+JEL-1) * SCLFAC
*. reverse lexical number of the double annihilated string
               ILEX = 1
               DO JJEL = 1, JEL-1
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL)
               END DO
               DO JJEL = JEL+1,IEL-1
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL-1)
               END DO
               DO JJEL = IEL+1, NKEL
                 ILEX = ILEX + IZ(KSTR(JJEL,KKSTR),JJEL-2)
               END DO
               IACT = IREO(ILEX)
               IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                 WRITE(6,*) ' IACT out of bounds, IACT =  ', IACT
                 STOP       ' IACT out of bounds '
               END IF
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' ILEX and IACT ', ILEX, IACT
               END IF
*
               IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
               I1(KKSTR-KMIN+1,IJ) = IACT
               XI1S(KKSTR-KMIN+1,IJ) = SIGN
             END IF 
*            ^ End if orbitals are in correct range
           END DO
*          ^ End of loop over IEL
         END DO
*        ^ End of loop over JEL
        END DO
*       ^ End of loop over Kstrings
*
      ELSE IF(IAC.EQ.2.AND.JAC.EQ.1) THEN
*
* ===================================
* Creation-annihilation map a+ i a j
* ===================================
*
C       DO KKSTR = 1, NKSTR
        DO KKSTR = KMIN,KEND 
*. Indicate where a given orbital i should be added in KKSTR
         IZERO = 0
         CALL ISETVC(ISCR(IORBMIN),IZERO,NIORB)
         IIEL = 1
         DO IIORB = IORBMIN,IORBMAX
 2810      CONTINUE
           IF(IIEL.LE.NKEL) THEN
             IF(IIORB.LT.KSTR(IIEL,KKSTR)) THEN
               ISCR(IIORB)=IIEL
             ELSE IF (IIORB.EQ.KSTR(IIEL,KKSTR)) THEN
               ISCR(IIORB) = - IIEL
               IIEL = IIEL+1
             ELSE IF (IIORB.GT.KSTR(IIEL,KKSTR)) THEN
               IIEL = IIEL + 1
               GOTO 2810
             END IF
           ELSE IF(IIEL.EQ.NKEL+1) THEN
              ISCR(IIORB) = IIEL
           END IF
         END DO
         IF(NTEST.GE.10000) THEN
           WRITE(6,*) ' ISCR from IORBMIN array for KKSTR = ', KKSTR
           WRITE(6,*) ' IORBMIN, NIORB', IORBMIN,NIORB
           CALL IWRTMA(ISCR(IORBMIN),1,NIORB,1,NIORB)
         END IF
         DO JEL = 1, NKEL
           JJORB = KSTR(JEL,KKSTR)
           IF(JJORB.GE.JORBMIN.AND.JJORB.LE.JORBMAX)THEN
             DO IIORB = IORBMIN,IORBMAX
               IEL = ISCR(IIORB)
C?             write(6,*) ' JEL IEL JJORB IIORB',JEL,IEL,JJORB,IIORB
               IACT = 0
               IF(IEL.GT.0.AND.IIORB.GT.JJORB) THEN
*. New string is  a+1 ... a+ jel-1 a+jel+1 ..a+iel-1 a+iiorb a+iel+1 ...
*. Lexical number of new string
                 ILEX = 1
                 DO KEL = 1, JEL-1
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 DO KEL = JEL+1, IEL-1
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL-1)
                 END DO
                 ILEX = ILEX + IZ(IIORB,IEL-1)
                 DO KEL = IEL, NKEL
                  ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 IACT = IREO(ILEX)
                 IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                   WRITE(6,*) ' 1: IACT out of bounds, IACT =  ', IACT
                   WRITE(6,*) ' NSTRI = ', NSTRI
                   WRITE(6,*) 'IIORB,JJORB ',IIORB,JJORB
                   WRITE(6,*) ' Kstring : '
                   CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                   WRITE(6,*) ' ILEX = ', ILEX
                   WRITE(6,*) 'IZ matrix'
                   CALL IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
                   STOP ' IACT out of bounds'
                 END IF
                 SIGN = (-1) ** (IEL+JEL-1) * SCLFAC
               ELSE IF(IEL.GT.0 .AND. IIORB.LT.JJORB) THEN
*. New string is  a+1 ... a+ iel-1 a+ iiorb a+iel+1 ..a+jel-1 a+jel+1 ...
                 ILEX = 1
                 DO KEL = 1, IEL-1
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO
                 ILEX = ILEX + IZ(IIORB,IEL)
                 DO KEL = IEL,JEL-1
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL+1)
                 END DO
                 DO KEL = JEL + 1, NKEL
                   ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                 END DO 
C?               write(6,*) ' ILEX =' , ILEX
                 IACT = IREO(ILEX)
                 IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                   WRITE(6,*) '2 IACT out of bounds, IACT =  ', IACT
                   STOP       ' IACT out of bounds '
                 END IF
C?               write(6,*) ' IACT = ', IACT
                 SIGN = (-1) ** (IEL+JEL  ) * SCLFAC
               ELSE IF (IEL.LT.0. AND. IIORB .EQ. JJORB) THEN
*. Diagonal excitation
                 SIGN = SCLFAC
                 IACT = KKSTR
               END IF
               IF(IACT.NE.0) THEN
                 IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
                 I1(KKSTR-KMIN+1,IJ) = IACT
                 XI1S(KKSTR-KMIN+1,IJ) = SIGN
               END IF
             END DO
*            ^ End of loop over IIORB
           END IF
*          ^ End of  active cases 
         END DO
*        ^ End of loop over electrons to be annihilated
        END DO
*       ^ End of loop over Kstrings 
      ELSE IF(IAC.EQ.1.AND.JAC.EQ.2) THEN
*
* ======================================
* Annihilation-creation  map a i a+ j
* ======================================
*
*. Diagonal excitations ?
        IF(IORBMIN.EQ.JORBMIN) THEN
         IDIAG = 1
        ELSE
         IDIAG = 0
        END IF
C       DO KKSTR = 1, NKSTR
        DO KKSTR = KMIN,KEND 
*. Indicate where a given orbital j should be added in KKSTR
         IZERO = 0
         CALL ISETVC(ISCR(JORBMIN),IZERO,NJORB)
         JJEL = 1
         DO JJORB = JORBMIN,JORBMAX
 0803      CONTINUE
           IF(JJEL.LE.NKEL) THEN
             IF(JJORB.LT.KSTR(JJEL,KKSTR)) THEN
               ISCR(JJORB)=JJEL
             ELSE IF (JJORB.EQ.KSTR(JJEL,KKSTR)) THEN
               ISCR(JJORB) = - JJEL
               JJEL = JJEL+1
             ELSE IF (JJORB.GT.KSTR(JJEL,KKSTR)) THEN
               JJEL = JJEL + 1
               GOTO 0803
             END IF
           ELSE IF(JJEL.EQ.NKEL+1)THEN
              ISCR(JJORB) = JJEL
           END IF
         END DO
         IF(NTEST.GE.10000) THEN
           WRITE(6,*) ' ISCR from JORBMIN array for KKSTR = ', KKSTR
           CALL IWRTMA(ISCR(JORBMIN),1,NJORB,1,NJORB)
         END IF
         DO IEL = 1, NKEL+IDIAG
*. IEL = NKEL + 1 will be used for excitations a+j aj
           IF(IEL.LE.NKEL) IIORB = KSTR(IEL,KKSTR)
           IF((IIORB.GE.IORBMIN.AND.IIORB.LE.IORBMAX).OR.
     &         IEL.EQ.NKEL+1 )THEN
             DO JJORB = JORBMIN,JORBMAX
               JEL = ISCR(JJORB)
C?             WRITE(6,*) ' IEL IIORB JEL JJORB ',
C?   &                      IEL,IIORB,JEL,JJORB
               IACT = 0
               IF(IEL.LE.NKEL) THEN
                 IF(JEL.GT.0.AND.JJORB.GT.IIORB) THEN
*. New string is  a+1 ... a+ iel-1 a+iel+1 ..a+jel-1 a+jjorb a+jel+1 ...
*. Lexical number of new string
                   ILEX = 1
                   DO KEL = 1, IEL-1
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   DO KEL = IEL+1, JEL-1
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL-1)
                   END DO
                   ILEX = ILEX + IZ(JJORB,JEL-1)
                   DO KEL = JEL, NKEL
                    ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   IACT = IREO(ILEX)
                   IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                     WRITE(6,*) '3 IACT out of bounds, IACT =  ', IACT
                     WRITE(6,*) ' ILEX,KKSTR ', ILEX, KKSTR
                     WRITE(6,*) ' occupation of KSTR '
                     CALL IWRTMA(KSTR,1,NKEL,1,NKEL)
                     WRITE(6,*) ' IEL JEL ', IEL,JEL
                     WRITE(6,*) ' IIORB,JJORB',IIORB,JJORB
                     STOP       ' IACT out of bounds '
                   END IF
                   SIGN = (-1) ** (IEL+JEL) * SCLFAC
                 ELSE IF(JEL.GT.0 .AND. JJORB.LT.IIORB) THEN
*. New string is  a+1 ... a+ jel-1 a+ jjorb a+jel+1 ..a+iel-1 a+iel+1 ...
                   ILEX = 1
                   DO KEL = 1, JEL-1
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO
                   ILEX = ILEX + IZ(JJORB,JEL)
                   DO KEL = JEL,IEL-1
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL+1)
                   END DO
                   DO KEL = IEL + 1, NKEL
                     ILEX = ILEX + IZ(KSTR(KEL,KKSTR),KEL)
                   END DO 
                   IACT = IREO(ILEX)
                   SIGN = (-1) ** (IEL+JEL-1) * SCLFAC
                   IF(IACT.LE.0.OR.IACT.GT.NSTRI) THEN
                     WRITE(6,*) '4 IACT out of bounds, IACT =  ', IACT
                     WRITE(6,*) ' NSTRI = ', NSTRI
                     WRITE(6,*) 'IIORB,JJORB ',IIORB,JJORB
                     WRITE(6,*) ' Kstring : '
                     CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                     WRITE(6,*) ' ILEX = ', ILEX
                     WRITE(6,*) 'IZ matrix'
                     CALL IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
                     STOP       ' IACT out of bounds '
                   END IF
                 END IF
                 IF(IACT.NE.0) THEN
                   IJ = (JJORB-JORB)*NIORB + IIORB-IORB+1
                   I1(KKSTR-KMIN+1,IJ) = IACT
                   XI1S(KKSTR-KMIN+1,IJ) = SIGN
                 END IF
               ELSE IF(IEL.EQ.NKEL+1.AND.JEL.GT.0) THEN
*. Diagonal excitations aja+j
                 JJ = (JJORB-JORB)*NJORB + JJORB-JORB+1
                 I1(KKSTR-KMIN+1,JJ) = KKSTR
                 XI1S(KKSTR-KMIN+1,JJ) = SCLFAC
               END IF
             END DO
*            ^ End of loop over JJORB
           END IF
*          ^ End of  active cases 
         END DO
*        ^ End of loop over electrons to be annihilated
        END DO
*       ^ End of loop over Kstrings 
      END IF
*.    ^ End of types of creation mappings
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADADST1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          IJ = 0
          DO  JJORB = JORB,JORB+NJORB-1
            JJORBR = JJORB-JORB+1
            DO  IIORB = IORB, IORB + NIORB - 1
              IJ = IJ + 1
C?            WRITE(6,*) ' IJ = ', IJ
C?            IF(IIORB.GT.JJORB) THEN
                IIORBR = IIORB - IORB + 1
                WRITE(6,'(A,2I4)')
     &          ' Info for orbitals (iorb,jorb) ', IIORB,JJORB
                WRITE(6,*) ' Excited strings and sign '
                CALL IWRTMA(I1(1,IJ),1,NK,1,NK)
                CALL WRTMAT(XI1S(1,IJ),1,NK,1,NK)
C?            END IF
            END DO
          END DO
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
     &                              NORBT,IDOREO,IZ,IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space 
*
*.If  IDOREO .NE. 0 THEN reordering array : lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete 
*                           orbital space, product of strings in
*                           each GAS space 
*
* Compared to GETSTR2_TOTSM_SPGP : Based upon IGRP(NIGRP)
*                                  (Just a few changes in the beginning)
*
* =====
* Input 
* =====
*
* IGRP :  supergroup, here as an array of GAS space 
* NIGRP : Number of active groups 
* ISPGRPSM : Total symmetry of superstrings 
* NEL : Number of electrons 
* IZ  : Reverse lexical ordering array for this supergroup (IF IDOREO.NE.0)
* 
*
* ======
* Output 
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0) 
*
*
* Jeppe Olsen, Written  July 1995
*              Version of Dec 1997
*. Last modification; Jeppe Olsen; May 2013; New ordering of symmetryblocks
*
c      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INTEGER IZ(NORBT,NEL)
      INTEGER IGRP(NIGRP)
*. output
      INTEGER ISTR(NEL,*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
C?    CALL LUCIAQENTER('GETST')
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A)')  ' Strings to be obtained : ' 
        WRITE(6,'(A)') ' **************************'
        WRITE(6,'(A)')
        WRITE(6,'(A,I2)') '   Symmetry : ', ISPGRPSM
        WRITE(6,'(A,16I3)') ' Groups : ', (IGRP(I),I=1,NIGRP)
        WRITE(6,*) ' NEL = ', NEL
        IF(IDOREO.NE.0) THEN 
          WRITE(6,*)
          WRITE(6,*) ' ============= '
          WRITE(6,*) ' The Z array : '
          WRITE(6,*) ' ============= '
          WRITE(6,*)
          WRITE(6,*) ' NORBT,NEL = ',NORBT,NEL
          CALL IWRTMA(IZ,NORBT,NEL,NORBT,NEL)
        END IF
      END IF
*. Absolut number of this supergroup
*. Occupation per gasspace
*. Largest occupied space 
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = IGRP(IGAS)
        NELFGS(IGAS) = NELFGP(IGRP(IGAS))          
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        CALL ICOPVE2(int_mb(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(int_mb(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
* Largest and lowest active symmetries for each GAS space
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF 
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL 
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGASL
         CALL NXTNUM3(ISMFGS(2),NGASL-1,MINVAL(2),MAXVAL(2),NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS(2),NGASL-1,1,NGASL-1,1)
        END IF
        ISTSMM1 = 1
        DO IGAS = 2, NGASL 
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. Required symmetry of GASpace 1
        CALL SYMCOM(2,1,ISTSMM1,ISMGS1,ISPGRPSM)
        ISMFGS(1) = ISMGS1
*. A test that ISFGS(1) is within bounds could be inserted here

*
         DO IGAS = NGASL+1,NGAS
           ISMFGS(IGAS) = 1
         END DO
         IF(NTEST.GE.200) THEN
           WRITE(6,*) ' Next symmetry distribution '
           CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
         END IF
*. Obtain all strings of this symmetry 
CT      CALL LUCIAQENTER('GASSM')
         CALL GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTR(1,ISTRBS),NSTR,NEL,
     &                           NNSTSGP,IISTSGP)
CT       CALL LUCIAQEXIT('GASSM')
*. Reorder Info : Lexical => actual number 
         IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces 
           DO JSTR = ISTRBS, ISTRBS+NSTR-1
             LEX = 1
             DO IEL = 1, NEL 
               LEX = LEX + IZ(ISTR(IEL,JSTR),IEL)
             END DO
             IF(NTEST.GE.100) THEN
               WRITE(6,*) ' string '  
               CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
               WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
             END IF
*
             MAXLEX = MAX(MAXLEX,LEX)
             IREO(LEX) = JSTR
           END DO
         END IF
*
        ISTRBS = ISTRBS + NSTR 
*. ready for next symmetry distribution 
        IF(NGAS-1.NE.0) GOTO 1000
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings : '
        WRITE(6,*)
        CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,NSTR,1,NSTR)
        END IF
      END IF
*
C?    CALL LUCIAQEXIT('GETST')
      RETURN
      END 
      SUBROUTINE RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                   JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                   IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                   SB,CB,ADSXA,STSTSX,MXPNGASX,
     &                   NOBPTS,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                   XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBSX,IUSEAB,
     &                   CJRES,SIRES,SCLFAC,NTESTG,
     &                   NSEL2E,ISEL2E,IUSE_PH,IPHGAS,XINT2)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
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
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT  : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + CJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks 
* October 96 : New routines for accessing annihilation information
*             Cleaned and shaved, only IROUTE = 3 option active
* October   97 : allowing for N-1/N+1 switch
*
* Last change : Aug 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
      REAL*8 INPROD
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*), XINT2(*)
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION H(MXPTSOB*MXPTSOB)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*. Arrays for reorganization 
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C    &          SIGNREO,NADOP,NADDEL,IADDEL,ADSIGN)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'oper.inc'
      CALL LUCIAQENTER('RS2Bn')
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
*
        WRITE(6,*) ' RSBB2BN speaking '
      IF(NTEST.GE.500) THEN
*
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' RSBB2BN speaking '
        WRITE(6,*) ' =============== '
*
        WRITE(6,*) ' Occupation of IA '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of IB '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JA '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JB '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)

*
        WRITE(6,*) ' Memcheck at start of RSBB2BN '
        CALL MEMCHK
        WRITE(6,*) ' Memory check passed '
*
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Groups defining each supergroup
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
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
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)           
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
*
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Should N-1 or N+1 projection be used for alpha strings
          IJ_TYP(1) = ITYP
          IJ_TYP(2) = JTYP
          IJ_AC(1)  = 2
          IJ_AC(2) =  1
          NOP = 2
          IF(IUSE_PH.EQ.1) THEN
            CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,
     &           SIGNIJ)
          ELSE
*. Enforced a+ a
            IJ_REO(1) = 1
            IJ_REO(2) = 2
            SIGNIJ = 1.0D0
          END IF
*. Two choices here :
*  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> ( good old creation mapping)
*  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka>  + delta(i,j)                   
C?        WRITE(6,*) ' RSBB2BN : IOP_REO : ', (IOP_REO(II),II=1,2)
          IF(IJ_REO(1).EQ.1.AND.IJ_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
            IJAC = 2
            SIGNIJ2 = SCLFAC
*
            IJ_DIM(1) = NI
            IJ_DIM(2) = NJ
            IJ_SYM(1) = ISM
            IJ_SYM(2) = JSM
            IJ_TYP(1) = ITYP
            IJ_TYP(2) = JTYP
*
            NOP1   = NI
            IOP1SM = ISM
            IOP1TP = ITYP
            NOP2   = NJ
            IOP2SM = JSM
            IOP2TP = JTYP
          ELSE
*. Terra Nova, annihilation map 
            IJAC = 1
            SIGNIJ2 = -SCLFAC
*
            IJ_DIM(1) = NJ
            IJ_DIM(2) = NI
            IJ_SYM(1) = JSM
            IJ_SYM(2) = ISM
            IJ_TYP(1) = JTYP
            IJ_TYP(2) = ITYP
*
            NOP1   = NJ
            IOP1SM = JSM
            IOP1TP = JTYP
            NOP2   = NI
            IOP2SM = ISM
            IOP2TP = ITYP
          END IF
*
*. Generate creation- or annihilation- mappings for all Ka strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          CALL ADAST_GAS(IJ_SYM(2),IJ_TYP(2),NGAS,JASPGP,JASM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SIGNIJ2,IJAC)
C         CALL ADAST_GAS(JSM,JTYP,JATP,JASM,IAGRP,
C    &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC)
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          CALL ADAST_GAS(IJ_SYM(1),IJ_TYP(1),NGAS,IASPGP,IASM,
     &         I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJAC)
C         CALL ADAST_GAS(ISM,ITYP,NGAS,IASPGP,IASM,
C    &         I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC)
*. Compress list to common nonvanishing elements
          IDOCOMP = 0
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),
     &                        NKASTR,NKAEFF)
          ELSE 
              NKAEFF = NKASTR
          END IF
            
*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK   
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
*
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, IJ_DIM(2)
              WRITE(6,*) ' JJ = ', JJ
              CALL GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
            IF(NTEST.GE.500) THEN
              WRITE(6,*) ' Updated CJRES as C(Kaj,Jb)'
              CALL WRTMAT(CJRES,NKASTR*NJ,NJB,NKASTR*NJ,NJB)
            END IF
*
            MXACJ=MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC*IJ_DIM(1))
* Searching for bug
C           WRITE(6,*) ' JTEST: Dimension of CJRES and SIRES ',
C    &      IJ_DIM(2)*LKABTC*NJB, IJ_DIM(1)*LKABTC*NIB
*
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*. Allowed double excitation ?
              IJKL_ACT = I_DX_ACT(ITYP,KTYP,LTYP,JTYP)
              IF(IJKL_ACT.EQ.0) GOTO 2000
              IF(NTEST.GE.500) THEN
                WRITE(6,*) ' KTYP, LTYP', KTYP, LTYP 
              END IF
*. Should this group of excitations be included 
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              KL_TYP(1) = KTYP
              KL_TYP(2) = LTYP
              KL_AC(1)  = 2
              KL_AC(2) =  1
              NOP = 2
              IF(IUSE_PH.EQ.1) THEN
                CALL ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,
     &               SIGNKL)
              ELSE
*. Enforced a+ a
                KL_REO(1) = 1
                KL_REO(2) = 2
                SIGNKL = 1.0D0
              END IF
*
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = ADSXA(KSM,KLSM)
                IF(NTEST.GE.1500) THEN
                  WRITE(6,*) ' KSM, LSM', KSM, LSM
                END IF
                IF(LSM.EQ.0) GOTO 1930
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*
                IF(KL_REO(1).EQ.1.AND.KL_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
                  KLAC = 2
                  KL_DIM(1) = NK
                  KL_DIM(2) = NL
                  KL_SYM(1) = KSM
                  KL_SYM(2) = LSM
                  KL_TYP(1) = KTYP
                  KL_TYP(2) = LTYP
                ELSE
*. Terra Nova, annihilation map 
                  KLAC = 1
                  KL_DIM(1) = NL
                  KL_DIM(2) = NK
                  KL_SYM(1) = LSM
                  KL_SYM(2) = KSM
                  KL_TYP(1) = LTYP
                  KL_TYP(2) = KTYP
                END IF
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0  
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
*. currently we are using creation mappings for kl
*. (Modify to use ADAST later )
                CALL ADAST_GAS(KL_SYM(2),KL_TYP(2),NGAS,JBSPGP,JBSM,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,SIGNKL,KLAC)
C               CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
C    &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADAST_GAS(KL_SYM(1),KL_TYP(1),NGAS,IBSPGP,IBSM,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,KLAC)
C               CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
C    &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (iop2 iop1 |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                ONE = 1.0D0
                IF(I_USE_SIMTRH .EQ.0 ) THEN
*. Normal integrals with conjugation symmetry
                  CALL LUCIAGETINT(XINT,IJ_TYP(2),IJ_SYM(2),
     &                 IJ_TYP(1),IJ_SYM(1),
     &                 KL_TYP(1),KL_SYM(1),KL_TYP(2),KL_SYM(2),IXCHNG,
     &                 0,0,ICOUL,ONE,ONE)
                ELSE IF (I_USE_SIMTRH.EQ.1) THEN
C?              WRITE(6,*) ' I_USE_SIMTRH = ', I_USE_SIMTRH
*. Integrals does not have conjugation symmetry so be careful...
*. The following is not enough is particle hole symmetry is encountered
*. Obtain ( i j ! k l )
                  CALL LUCIAGETINT(XINT,ITYP,ISM,JTYP,JSM,
     &                             KTYP,KSM,LTYP,LSM,
     &                        IXCHNG,0,0,ICOUL,ONE,ONE)
                  IF(KLAC.EQ.2.AND.IJAC.EQ.2) THEN
*. Transpose to obtain ( j i ! k l )
                    CALL TRP_H2_BLK(XINT,12,NI,NJ,NK,NL,XINT2)
                  ELSE IF(KLAC.EQ.1.AND.IJAC.EQ.2) THEN  
*. Transpose to obtain (j i | l k)
                    CALL TRP_H2_BLK(XINT,46,NI,NJ,NK,NL,XINT2)
                  ELSE IF (KLAC.EQ.1.AND. IJAC .EQ. 1 ) THEN
*. Transpose to obtai (i j | l k)
                    CALL TRP_H2_BLK(XINT,34,NI,NJ,NK,NL,XINT2)
                  END IF
                END IF
*
* S(Ka,i,Ib) = sum(j,k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                IJKL_DIM = IJ_DIM(1)*IJ_DIM(2)*KL_DIM(1)*KL_DIM(2)
                IF(INPROD(XINT,XINT,IJKL_DIM).NE.0.0D0) THEN
                IROUTE = 3
                CALL SKICKJ(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),
     &               KL_DIM(1),KL_DIM(2),
     &               NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE )
                END IF
*
                IF(NTEST.GE.2500) THEN
                  WRITE(6,*) ' Updated Sires as S(Kai,Ib)'
                  CALL WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
                END IF
*
 1930         CONTINUE
*             ^ End of loop over KSM
 2000       CONTINUE
*           ^ End of loop over KLTYP
*
*. Scatter out from s(Ka,Ib,i)
*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' S(Ka,Ib,i) as S(Ka,Ibi)'
              CALL WRTMAT(SIRES,LKABTC,NIB*IJ_DIM(1),LKABTC,IJ_DIM(1))
            END IF
*
            DO II = 1, IJ_DIM(1)
              CALL ADD_SKAIIB(SB,IJ_DIM(1),NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 1940   CONTINUE
*       ^ End of loop over ISM
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL LUCIAQEXIT('RS2Bn')
      RETURN
      END
      SUBROUTINE H1TERMS(SIGMA,C,NR,NK,NI,NJ,MAXK,HIJ,
     &                  KI,XKI,KJ,XKJ,SCLFAC,NCS,NCC,ITRANS)
*
* Update Sigma with one-electron contributions
*
* Sigma(Ir,Ic) = Sigma(Ir,Ic) + Sclfac*H_ij <Ic|a+ ic a jc!Jc> C(Ir,Jc) 
*
* If ITRANS.EQ.1, then Sigma and C are transposed
*
* corresponding to current set of intermediate strings
*
* Jeppe Olsen, Fall of 97    
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION C(*),HIJ(NI,NJ)
*. single operator connections
      DIMENSION KI(MAXK,*),XKI(MAXK,*)
      DIMENSION KJ(MAXK,*),XKJ(MAXK,*)
*. Input and output
      DIMENSION SIGMA(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' H1TERMS speaking '
        WRITE(6,*) ' HIJ matrix '
        CALL WRTMAT(HIJ,NI,NJ,NI,NJ)
        IF(ITRANS.EQ.1) THEN
          WRITE(6,*) ' Initial sigma block (transposed )'
          CALL WRTMAT(SIGMA,NCS,NR,NCS,NR)
        END IF
      END IF
*. Loop over intermediate strings
C?    WRITE(6,*) ' NK = ', NK
      DO KC = 1, NK
*
*. Number of nonvanishing connections from KC
*
       JJ = 0
       DO JORB = 1, NJ
         IF(KJ(KC,JORB).NE.0) JJ = JJ + 1
       END DO
       II = 0
       DO IORB = 1, NI
         IF(KI(KC,IORB).NE.0) II = II + 1
       END DO
C?     WRITE(6,*) ' II, JJ = ', II,JJ
*
       IF(II.NE.0.AND.JJ.NE.0) THEN
         DO IORB = 1, NI
           IC = KI(KC,IORB)
           IF(IC.NE.0) THEN
             SGNI = XKI(KC,IORB)
             DO JORB = 1, NJ
               JC = KJ(KC,JORB)
C?             write(6,*) ' IORB JORB', IORB,JORB
               IF(JC.NE.0) THEN
                 SGNJ = XKJ(KC,JORB)
C?               WRITE(6,*) ' SGNI SGNJ SCLFAC',SGNI,SGNJ,SCLFAC
                 FACTOR = SGNI*SGNJ*HIJ(IORB,JORB)*SCLFAC
C?               WRITE(6,*) ' IC JC FACTOR ', IC,JC,FACTOR
                 ONE = 1.0D0
                 IF(ITRANS.EQ.0) THEN
                   CALL VECSUM(SIGMA(1+(IC-1)*NR),SIGMA(1+(IC-1)*NR),
     &                  C(1+(JC-1)*NR),ONE,FACTOR,NR)
                 ELSE
*. Not pretty, should be done somewhere else or on transposed matrices
                   DO IR = 1, NR
                     SIGMA(IC+(IR-1)*NCS) = SIGMA(IC+(IR-1)*NCS) 
     &              +FACTOR * C(JC+(IR-1)*NCC)
                   END DO
                 END IF
*                ^ End of transpose switch
               END IF
             END DO
           END IF
         END DO
       END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        IF(ITRANS.EQ.1) THEN
          WRITE(6,*) ' Updated sigma block (transposed )'
          CALL WRTMAT(SIGMA,NCS,NR,NCS,NR)
        END IF
      END IF

      RETURN
      END
      SUBROUTINE CHAR_TO_REAL(CHAR_X,REAL_X,L_CHAR_X)
* Character to REAL, using file 8     
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*102 CHAR_X
*
CNW   LU_INTERNAL = 8
CNW   REWIND(LU_INTERNAL)
CNW   WRITE(LU_INTERNAL,'(A)') CHAR_X
CNW   REWIND(LU_INTERNAL)
CNW   READ(LU_INTERNAL,*) REAL_X
      read (CHAR_X,*) REAL_X
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' INFO FROM CHAR_TO_REAL'
        WRITE(6,*) ' CHARACTER INPUT ', CHAR_X
        WRITE(6,*) ' REAL OUTPUT ',REAL_X
      END IF
*
      RETURN
      END
      SUBROUTINE CHAR_TO_INTEGER(CHAR_X,INT_X,L_CHAR_X)
* Character to INTEGER, using file 8     
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*102 CHAR_X
*
C?    WRITE(6,*) ' CHAR_TO_INTEGER, Length of string', L_CHAR_X
C?    WRITE(6,'(A)') ' CHARACTER INPUT ', CHAR_X
*
CNW   LU_INTERNAL = 8
CNW   REWIND(LU_INTERNAL)
CNW   WRITE(LU_INTERNAL,'(A)') CHAR_X
CNW   REWIND(LU_INTERNAL)
CNW   READ(LU_INTERNAL,*) INT_X
      read (CHAR_X,*) INT_X
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' INFO FROM CHAR_TO_INTEGER'
        WRITE(6,'(A)') ' CHARACTER INPUT ', CHAR_X
        WRITE(6,*) ' INTEGER OUTPUT ', INT_X
      END IF
*
      RETURN
      END
      SUBROUTINE DECODE_LINE(LINE,NCHAR,NENTRY,IENTRY,MXENTRY)
*
* A CHAR line is given.
* Find number of separate items, with each item
* being separated by a ,
*
* Jeppe Olsen, Someday in 97 where I really should be doing more
* important things
*
*. Entry
      CHARACTER*(*) LINE
*. Output
      CHARACTER*102 IENTRY(MXENTRY)
*. Local scratch
      CHARACTER*102 CSCR
*
*  
      DO JCHAR = 1, NCHAR
       CSCR(JCHAR:JCHAR) = ' '
      END DO
      JITEM=0
      JEFF = 0
      DO ICHAR = 0, NCHAR
*a pure comment
        IF(ICHAR.EQ.0.AND.LINE(1:1).EQ.'!') THEN
          EXIT
        ELSE IF(ICHAR.EQ.0.OR.LINE(ICHAR:ICHAR).EQ.',') THEN
*Start of new item, make sure there is space and clean up
          JITEM = JITEM + 1
          IF(JITEM .GT.MXENTRY) THEN
            WRITE(6,*) 'DECODE_LINE:MXENTRY too small'
            WRITE(6,*) ' Number of entries larger than MXENTRY'
            WRITE(6,*) ' JITEM, MXENTRY', JITEM, MXENTRY
            STOP'DECODE_LINE:MXENTRY too small'
          END IF
*. Copy previous entry to permanent
          IF(JITEM.NE.1) THEN
            IENTRY(JITEM-1) = CSCR
          END IF
*. and clean
          DO JCHAR = 1, NCHAR
            CSCR(JCHAR:JCHAR) = ' '
          END DO
          JEFF = 0
*. a comment coming up?
        ELSE IF (LINE(ICHAR:ICHAR).EQ.'!') THEN
*. Copy previous entry to permanent
          IF(JITEM.NE.0) THEN
            IENTRY(JITEM) = CSCR
          END IF
*. and basta ...
          EXIT
        ELSE
*. Continuation of previous item
          JEFF = JEFF + 1
          CSCR(JEFF:JEFF) = LINE(ICHAR:ICHAR)
        END IF
*
      END DO
*. Transfer last item to permanant residence
      IF(JEFF.NE.0) IENTRY(JITEM) = CSCR
*
      NENTRY = JITEM
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from Decode line '
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        WRITE(6,*) ' Number of separate entries', NENTRY
        WRITE(6,*)
        DO JENTRY = 1, NENTRY
          WRITE(6,'(A,I3,102A)') 'Entry ',JENTRY,IENTRY(JENTRY)
        END DO
      END IF
*
      RETURN
      END

   
      SUBROUTINE ADAST_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
     &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC,
     &                    IAC)
*
*
* Obtain creation or annihilation mapping
*
* IAC = 2 : Creation map
* a+IORB !KSTR> = +/-!ISTR> 
*
* IAC = 1 : Annihilation map
* a IORB !KSTR> = +/-!ISTR> 
*
* for orbitals of symmetry IOBSM and type IOBTP
* and Istrings defined by the NIGRP groups IGRP and symmetry ISPGPSM
* 
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR> 
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version 
*               October 96   : Improved version
*               September 97 : annihilation mappings added
*                              I groups defined by IGRP
*               May 2013     : New order of symmetrydistributions
*. Last Revision; May 5, 2013; Jeppe Olsen; New order of sym. dist.
*
*
* ======
*. Input
* ======
*
*./BIGGY
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER ISMFGS(MXPNGAS)
      INTEGER MXVLI(MXPNGAS),MNVLI(MXPNGAS)
      INTEGER MXVLK(MXPNGAS),MNVLK(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      INTEGER KGRP(MXPNGAS)
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8)
      DIMENSION IOFFI(LOFFI)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'multd2h.inc'
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension 
* (NKSTR,*), Where NKSTR is the number of K-strings of 
*  correct symmetry . Nk is provided by this routine.
*
       CALL LUCIAQENTER('ADAST')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADAST_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        WRITE(6,*) '  Symmetry of supergroup : ', ISPGPSM
        WRITE(6,*) ' SCLFAC = ', SCLFAC
*
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' Annihilation mapping '
        ELSE IF(IAC.EQ.2) THEN
          WRITE(6,*) ' Creation mapping '
        ELSE 
          WRITE(6,*) ' Unknown IAC parameter in ADAST ',IAC
          STOP       ' Unknown IAC parameter in ADAST '
        END IF
*
      END IF
*. A few preparations
      NORBTS= NOBPTS(IOBTP,IOBSM)
      NORBT= NOBPT(IOBTP)
      IACGAS = IOBTP
*. First orbital of given GASpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1+NINOB
*. First orbital of given GASpace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM) 
       IF(NTEST.GE.1000)
     & WRITE(6,'(A,2I4)') ' IBORBSP, IBORBSPS = ', IBORBSP, IBORBSPS
*
*====================================================
*. K strings : Supergroup, symmetry and distributions
*====================================================
*
      IF(IAC.EQ.1) THEN
       IDELTA = +1
      ELSE
       IDELTA = -1
      END IF
*. Is required mapping contained within current set of maps?
*. a:) Is active GASpace included in IGRP - must be 
      IACGRP = 0
      DO JGRP = 1, NIGRP
       IF(IGSFGP(IGRP(JGRP)).EQ. IACGAS) IACGRP = JGRP
      END DO
*. Note : IACGRP is not the actual active group, it is the address of the
*         active group in IGRP
      IF(IACGRP.EQ.0) THEN
        WRITE(6,*) ' ADAST in problems '
        WRITE(6,*) ' Active GASpace not included in IGRP '
        WRITE(6,*) ' Active GASpace : ', IACGAS
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        STOP       ' ADAST : Active GASpace not included in IGRP '
      END IF
*. b:) active group in K strings
      NIEL = NELFGP(IGRP(IACGRP))
      NKEL = NIEL + IDELTA
      IF(NTEST.GE.1000) WRITE(6,*) ' NIEL and NKEL ',NIEL,NKEL
      IF(NKEL.EQ.-1.OR.NKEL.EQ.NOBPT(IACGAS)+1) THEN
*. No strings with this number of elecs - be happy : No work 
        NKSTR = 0
        KACT = 0
        KACGRP = 0
        GOTO 9999
      ELSE
*. Find group with NKEL electrons in IACGAS
        KACGRP = 0
        DO JGRP = IBGPSTR(IACGAS),IBGPSTR(IACGAS)+NGPSTR(IACGAS)-1
          IF(NELFGP(JGRP).EQ.NKEL) KACGRP = JGRP
        END DO
        IF(NTEST.GE.1000) WRITE(6,*) ' KACGRP = ',KACGRP
*. KACGRP is the Active group itself     
        IF(KACGRP.EQ.0) THEN
          WRITE(6,*)' ADAST : cul de sac, active K group not found'
          WRITE(6,*)' GAS space and number of electrons ',
     &               IACGAS,NKEL
          STOP      ' ADAST : cul de sac, active K group not found'
        END IF
      END IF
*. Okay active K group was found and is nontrivial
C     CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
      KSM = MULTD2H(IOBSM,ISPGPSM)
*. The K supergroup
      CALL ICOPVE(IGRP,KGRP,NIGRP)
      KGRP(IACGRP) = KACGRP
*. Number of strings and symmetry distributions of K strings
      CALL NST_SPGRP(NIGRP,KGRP,KSM,int_mb(KNSTSGP(1)),
     &               NSMST,NKSTR,NKDIST)
C          NST_SPGRP(NGRP,IGRP,ISM_TOT,NSTSGP,NSMST,NSTRIN,NDIST)
      IF(NTEST.GE.1000) WRITE(6,*) 
     & ' KSM, NKSTR : ', KSM, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*. Last active space in K strings and number of strings per group and sym
      NGASL = 1
      DO JGRP = 1, NIGRP
       IF(NELFGP(KGRP(JGRP)).GT.0) NGASL = JGRP
       CALL ICOPVE2(int_mb(KNSTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &              NNSTSGP(1,JGRP))
       CALL ICOPVE2(int_mb(KISTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &              IISTSGP(1,JGRP))
      END DO
C     NGASL = NIGRP
*. MIN/MAX for Kstrings
      CALL MINMAX_FOR_SYM_DIST(NIGRP,KGRP,MNVLK,MXVLK,NKDIST_TOT)
      IF(NTEST.GE.1000) THEN
        write(6,*) 'MNVLK and MXVLK '
        CALL IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
      END IF
*. (NKDIST_TOT is number of distributions, all symmetries )
* ==============
*. I Strings 
* ==============
*. Generate symmetry distributions of I strings with given symmetry
      CALL TS_SYM_PNT2(IGRP,NIGRP,MXVLI,MNVLI,ISPGPSM,
     &                 IOFFI,LOFFI)
*. Offset and dimension for active group in I strings
      CALL ICOPVE2(int_mb(KISTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               IACIST)
C?    WRITE(6,*) ' IACIST for IACGRP,IGRP = ', IACGRP,IGRP(IACGRP)
C?    CALL IWRTMA(IACIST,1,NSMST,1,NSMST)
*
      CALL ICOPVE2(int_mb(KNSTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               NACIST)
*. Number of I strings per group and sym
COLD  DO IGAS = 1, NIGRP
COLD   CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(IGAS)-1)*NSMST+1,NSMST,
COLD &              IISTSGP)
COLD  END DO
*. Last entry in IGRP with a nonvanisking number of strings
      NIGASL = 1
      DO JGRP = 1, NIGRP
        IF(NELFGP(IGRP(JGRP)).GT.0) NIGASL = JGRP
      END DO
C?    WRITE(6,*) ' NIGASL = ', NIGASL
C     NIGASL = NIGRP
       
*. Number of electrons before active space
      NELB = 0
      DO JGRP = 1, IACGRP-1
        NELB = NELB + NELFGP(IGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NELB = ', NELB
*
      ZERO =0.0D0
      IZERO = 0    
COLD  CALL SETVEC(XI1S,ZERO,NORBTS*NKSTR)
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
      DO IGAS = 1, NIGRP
        ISMFGS(IGAS) = 1
      END DO
 1000 CONTINUE
*. Next distribution
C       CALL NEXT_SYM_DISTR(NIGRP,MNVLK,MXVLK,ISMFGS,KSM,KFIRST,NONEW)
        CALL NEXT_SYM_DISTR(NGASL,MNVLK,MXVLK,ISMFGS,KSM,KFIRST,NONEW)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' Symmetry distribution ' 
          call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
        END IF
        IF(NONEW.EQ.1) GOTO 9999
        KFIRST = 0
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IACGRP)
C       CALL  SYMCOM(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        CALL  SYMCOM(3,1,IOBSM,ISMFGS(IACGRP),IACSM)
        ISMFGS(IACGRP) = IACSM
        IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,IOFFI,MXVLI,MNVLI)
        ISMFGS(IACGRP) = ISAVE
C?      WRITE(6,*) ' IBSTRINI ', IBSTRINI
*. Number of strings before active GAS space
        NSTB = 1
C       DO IGAS = 1, IOBTP-1
        DO IGAS = 1, IACGRP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings After active GAS space
        NSTA = 1
C       DO IGAS =  IOBTP +1, NIGRP
        DO IGAS =  IACGRP+1, NIGRP
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group 
        NIAC  = NACIST(IACSM)
        IIAC =  IACIST(IACSM)
C?      WRITE(6,*) ' IIAC, IACSM = ',IIAC, IACSM
*
        NKAC = NNSTSGP(ISMFGS(IACGRP),IACGRP)
        IKAC = IISTSGP(ISMFGS(IACGRP),IACGRP)
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
        IF(NTEST.GE.1000) THEN
        write(6,*) ' nstb nsta niac nkac ',
     &               nstb,nsta,niac,nkac
        END IF
*. Obtain annihilation/creation mapping for all strings of this type
*. Are group mappings in expanded or compact form 
        IF(IAC.EQ.1.AND.ISTAC(KACGRP,2).EQ.0) THEN
          IEC = 2
          LROW_IN = NKEL
        ELSE 
          IEC = 1
C         LROW_IN = NORBTS
          LROW_IN = NORBT
        END IF
        NKACT = NSTFGP(KACGRP)
*
        MXAADST = MAX(MXAADST,NKSTR*NORBTS)
C     COMMON/COMJEP/MXACJ,MXACIJ,MXAADST
        IF(NSTA*NSTB*NIAC*NKAC.NE.0)
     &  CALL ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,   
     &            int_mb(KSTSTM(KACGRP,1)),int_mb(KSTSTM(KACGRP,2)),
     &            IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &            NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC,IAC,
     &            LROW_IN,IEC)
        KSTRBS = KSTRBS + NKSD     
C       IF(NGASL-1.GT.0) GOTO 1000
        GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADAST_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMAT(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
      CALL LUCIAQEXIT('ADAST')
      RETURN
      END
      SUBROUTINE ADAST_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,
     &              IOFFKSP,ICREORB,ICRESTR,
     &              IORBTSF,IORBTF,NORBTS,NSTAK,NSTAKT,NSTAI,
     &              NSTAKTS,ISTAKTS,NELB,NACGSOB,
     *              ISTMAP,SGNMAP,SCLFAC,IAC,LROW_IN,IEC)
*
* Annihilation/Creation mappings from K-strings of given sym in each gasspace
*
* Input 
* NSTB : Number of strings before active gasspace
* NSTA : Number of strings after accive gasspace
* IOFFK : Offset for K group of strings in active gasspace, i.e. start of
*         this symmetry of active K group strings
* IOFFI : Offset for I group of strings in active gasspace, i.e. start of
*         this symmetry of active I group strings
* IOFFISP: Offset for this symmetrydistribution of active I supergroup strings
* IOFFKSP: Offset for this symmetrydistribution of active K supergroup strings
* ICREORB : Orbital part of creation map for active K groupstrings
* ICRESTR : String  part of creation map for active K groupstrings
* IORBTSF   : First active orbital ( first orbital in in active GASspace
*           with required sym)
* IORBTF   : First orbital in active gas space, (can have any sym)
* NORBTS  : Number of orbitals of given symmetry and type            
* NSTAK : Number of K groupstrings with given correct symmetry      
* NSTAKT: Total Number of K groupstrings in active group (all symmetries)
* NSTAKTS: Total Number of K supergroup strings with correct symmetry
* ISTAKTS: Offset for K supergroup strings with hiven symmetrydistribution
* NSTAI : Number of I groupstrings in active gasspace
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION ICREORB(LROW_IN,*), ICRESTR(LROW_IN,*)
*. Output
      DIMENSION ISTMAP(NSTAKTS,*),SGNMAP(NSTAKTS,*)
* 
      IMULTK = NSTAK*NSTB
      IMULTI = NSTAI*NSTB
*
      NTEST = 000
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Reporting from ADAST_GASSM '
        WRITE(6,*) ' ======================== '
        WRITE(6,*) ' ICRESTR '    
        CALL IWRTMA(ICRESTR,LROW_IN,NSTAK,LROW_IN,NSTAK)
        WRITE(6,*) ' ICREORB '    
        CALL IWRTMA(ICREORB,LROW_IN,NSTAK,LROW_IN,NSTAK)
        WRITE(6,*) ' IOFFI = ', IOFFI
        WRITE(6,'(A,2I4)') ' IORBTSF, IORBTF =', IORBTSF, IORBTF
      END IF
      SIGN0 = (-1)**NELB*SCLFAC 
      DO KSTR = IOFFK, NSTAK+IOFFK-1
        DO IORB = IORBTSF, IORBTSF-1+NORBTS
*. Relative to Type-symmetry start
          IORBRTS = IORB-IORBTSF+1
*. Relative to type start
          IORBRT = IORB-IORBTF+1
          IF(NTEST.GE.1000)
     &    write(6,'(A,3I4)') ' IORB, IORBRTS IORBRT',IORB,IORBRTS,IORBRT
*. Change of active group
          I_AM_ACTIVE = 0
          IF(IAC.EQ.2) THEN
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' ICREORB = ', ICREORB(IORBRT,KSTR)
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' ICRESTR = ', ICRESTR(IORBRT,KSTR)
            IF(ICREORB(IORBRT,KSTR) .GT. 0 ) THEN
*. Creation is nonvanishing
              I_AM_ACTIVE = 1
              IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
                SIGN = SIGN0
                ISTR = ICRESTR(IORBRT,KSTR)
              ELSE
                SIGN = -SIGN0
                ISTR = -ICRESTR(IORBRT,KSTR)
              END IF
            END IF
          ELSE IF(IAC.EQ.1) THEN
             IF(IEC.EQ.1) THEN
*. Expanded map
               IF(ICREORB(IORBRT,KSTR) .LT. 0 ) THEN
*. Annihilation is non-vanishing
                 I_AM_ACTIVE = 1
                 IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
                   SIGN = SIGN0
                   ISTR = ICRESTR(IORBRT,KSTR)
                 ELSE
                   SIGN = -SIGN0
                   ISTR = -ICRESTR(IORBRT,KSTR)
                 END IF
               END IF
             ELSE
*. Compressed map 
               DO IROW = 1, LROW_IN
                IF(NTEST.GE.1000) WRITE(6,*)
     &          ' IROW, ICREORB(IROW,KSTR)', IROW, ICREORB(IROW,KSTR)
                 IF(ICREORB(IROW,KSTR) .EQ. -IORB   ) THEN
*. Annihilation is non-vanishing
                   I_AM_ACTIVE = 1
                   IF(ICRESTR(IROW,KSTR) .GT. 0 ) THEN
                     SIGN = SIGN0
                     ISTR = ICRESTR(IROW,KSTR)
                   ELSE
                     SIGN = -SIGN0
                     ISTR = -ICRESTR(IROW,KSTR)
                   END IF
                 END IF
               END DO
             END IF
*            ^ End of expanded/compact switch
           END IF
*          ^ End of Creation/annihilation switch
           
          IF(I_AM_ACTIVE .EQ. 1  ) THEN
*. Excitation is open, corresponding active I string
* Relative to start of given symmetry for this group
            ISTR = ISTR - IOFFI+ 1
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' ISTR, relative = ', ISTR
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
            IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
            IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
            IF(NTEST.GE.1000) WRITE(6,*)
     &      ' IADRK0 IOFFK IOFFKSP ', IADRK0,IOFFK,IOFFKSP
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' IADRI0, IOFFISP ', IADRI0, IOFFISP
*
            NSTAINSTA = NSTAI*NSTA
            NSTAKNSTA = NSTAK*NSTA
*
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' ISTR NSTA NSTB ',ISTR,NSTA,NSTB
            IF(NTEST.GE.1000)
     &      WRITE(6,*) ' NSTAI,NSTAK',NSTAI,NSTAK
            DO IB = 1, NSTB
              DO IA = 1, NSTA
                ISTMAP(IADRK0+IA,IORBRTS) = IADRI0 + IA
                SGNMAP(IADRK0+IA,IORBRTS) = SIGN
              END DO
              IADRI0 = IADRI0 +  NSTAINSTA
              IADRK0 = IADRK0 +  NSTAKNSTA
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output '
        WRITE(6,*) ' ====== '
        NK = NSTB*NSTAK*NSTA
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO IORB = IORBTSF,IORBTSF + NORBTS  - 1 
            IORBR = IORB-IORBTSF+1
            WRITE(6,*) ' Update Info for orbital ', IORB
            WRITE(6,*) ' Mapped strings and sign '
            CALL IWRTMA(ISTMAP(1,IORBR),1,NK,1,NK)
            CALL WRTMAT(SGNMAP(1,IORBR),1,NK,1,NK)
          END DO
        END IF
      END IF

      RETURN
      END
      SUBROUTINE NST_SPGRP(NNGRP,IGRP,ISM_TOT,NSTSGP,NSMST,
     &                     NSTRIN,NDIST)
*
* Number of strings for given combination of groups and 
* symmetry.
*
*. Input
*        
*
*   NNGRP : Number of active groups 
*   IGRP : The active groups
*   ISM_TOT : Total symmetry of supergroup
*   NSTSGP  : Number of strings per symmetry and supergroup
*   NSMST   : Number of string symmetries
*
*. Output
*
*  NSTRIN : Number of strings with symmetry ISM_TOT
*  NDIST  : Number of symmetry distributions
*
* Jeppe Olsen, September 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific Input
      DIMENSION IGRP(NNGRP)
*. General input
      DIMENSION NSTSGP(NSMST,*)
*. Scratch 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'multd2h.inc'
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER ISCR1(8), ISCR2(8)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' NST_SPGP is speaking '
        WRITE(6,*) ' ====================='
*
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NNGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NNGRP)
      END IF
*
      I_NEW_OR_OLD = 1
      NDIST = 1
      IF(I_NEW_OR_OLD.EQ.1) THEN
*. Start with first nontrivial group 
        IFIRST = 0
        DO KGRP = 1, NNGRP
*. A trivial group is identified by one total sym string, 
*. and no strings of the remaining symmetries, so
          IREM = 0
          DO KSM = 2, NSMST
            IREM = IREM + NSTSGP(KSM,IGRP(KGRP))
          END DO
          IF(.NOT.(IREM.EQ.0.AND.NSTSGP(1,IGRP(KGRP)).EQ.1).AND.
     &       IFIRST.EQ.0) THEN
C         IF(MINMAX_SM_GP(2,IGRP(KGRP)).NE.1.AND.IFIRST.EQ.0) THEN
            IFIRST = KGRP
            DO KSM = 1, NSMST
              ISCR1(KSM) = NSTSGP(KSM,IGRP(KGRP))
            END DO
          END IF
        END DO
*
        IF(IFIRST.EQ.0) THEN
*. String type with no electrons so 
          IF(ISM_TOT.EQ.1) THEN
            NSTRIN = 1
            NDIST = 1
          ELSE 
            NSTRIN = 0
            NDIST = 0
          END IF
          GOTO 9999
        END IF
*
        DO I = 1, NSMST
          ISCR2(I) = 0
        END DO
*
        DO KGRP = IFIRST+1, NNGRP
*. Is this group trivial ?
        I_AM_TRIVIAL = 0
        INONSYM= 0
        DO KSM = 2, NSMST
          INONSYM = INONSYM + NSTSGP(KSM,IGRP(KGRP))
        END DO
        IF(INONSYM.EQ.0.AND.NSTSGP(1,IGRP(KGRP)).EQ.1) I_AM_TRIVIAL = 1
        IF(I_AM_TRIVIAL.EQ.0) THEN
C       IF(MINMAX_SM_GP(2,IGRP(KGRP)).NE.1) THEN
*
          NDIST = NDIST 
     &   *(MINMAX_SM_GP(2,IGRP(KGRP))-MINMAX_SM_GP(1,IGRP(KGRP))+1)
*
          DO KSM = MINMAX_SM_GP(1,IGRP(KGRP)),MINMAX_SM_GP(2,IGRP(KGRP))
            DO KSMM1 = 1, NSMST
              KKSM = MULTD2H(KSM,KSMM1)
              ISCR2(KKSM) = ISCR2(KKSM) + 
     &                     ISCR1(KSMM1)*NSTSGP(KSM,IGRP(KGRP))
            END DO
          END DO
          DO KSM = 1, NSMST
            ISCR1(KSM) = ISCR2(KSM)
            ISCR2(KSM) = 0
          END DO
        END IF
        END DO
        NSTRIN = ISCR1(ISM_TOT)
      ELSE IF (I_NEW_OR_OLD.EQ.2) THEN
*
*. Set up min and max values for symmetries
      CALL MINMAX_FOR_SYM_DIST(NNGRP,IGRP,MNSM,MXSM,NDISTX)
*. Loop over symmetry distributions
      IFIRST = 1
      LENGTH = 0 
      NDIST = 0
 1000 CONTINUE
*. Next symmetry distribution
        CALL NEXT_SYM_DISTR(NNGRP,MNSM,MXSM,ISM,ISM_TOT,IFIRST,NONEW)
        IF(NONEW.EQ.0) THEN
          LDIST = 1
          DO JGRP = 1, NNGRP
            LDIST = LDIST*NSTSGP(ISM(JGRP),IGRP(JGRP))
          END DO
          LENGTH = LENGTH + LDIST
          NDIST = NDIST + 1
      GOTO 1000
        END IF
*
      NSTRIN = LENGTH
      END IF
 9999 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings obtained ', NSTRIN
        WRITE(6,*) ' Number of symmetry-distributions',NDIST
      END IF
*
      RETURN
      END
      SUBROUTINE GET_SPGP_INF(ISPGP,ITP,IGRP)
*
* Obtain groups defining supergroup ISPGP,ITP
*
* Jeppe Olsen, November 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
*. Output 
      DIMENSION IGRP(*)
*
      NTEST = 00
*. Absolute group number
C?    WRITE(6,*) ' GET_SPGP_INF : ISPGP, ITP', ISPGP, ITP
      ISPGPABS = ISPGP + IBSPGPFTP(ITP) -1
      DO IGAS = 1, NGAS
        IGRP(IGAS) = ISPGPFTP(IGAS,ISPGPABS)
      END DO
C     CALL ICOPVE(ISPGPFTP(1,ISPGPABS),IGRP,NGAS)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_SPGP_INF : ISPGP ITP ISPGPABS',
     &              ISPGP, ITP, ISPGPABS 
        WRITE(6,*) ' Groups of supergroups'
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE ANNSTR_GAS(STRING,NSTINI,NSTINO,NEL,NORB,IORBOF,
     &                  Z,NEWORD,LSGSTR,ISGSTI,ISGSTO,TI,TTO,NOCOB,
     &                  IEC,LDIM,IPRNT)
*
* A group of strings containing NEL electrons is given
* Set up all possible ways of removing an electron 
*
*========
* Input :
*========
* STRING : Input strings containing NEL electrons
* NSTINI : Number of input  strings
* NSTINO : Number of output strings
* NEL    : Number of electrons in input strings
* NORB   : Number of orbitals
* IORBOF : Number of first orbital
* Z      : Lexical ordering matrix for output strings containing
*          NEL - 1 electrons
* NEWORD : Reordering array for N-1 strings
* LSGSTR : .NE.0 => Include sign arrays ISGSTI,ISGSTO of strings
* ISGSTI : Sign array for NEL   strings
* ISGSTO : Sign array for NEL-1 strings
* IEC    : = 1 Extended map, dimension equals number of orbs
* IEC    : = 2 Compact  map, dimension equals number of elecs
* LDIM   : Row dimension ( see IEC)
*
*=========
* Output :
*=========
*
*TI      : TI(I,ISTRIN) .gt. 0 indicates that orbital I can be added
*          to string ISTRIN .
*TTO     : Resulting NEL + 1 strings
*          if the string have a negative sign
*          then the phase equals - 1
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER STRING,TI,TTO,STRIN2,Z
*.Input
      DIMENSION STRING(NEL,NSTINI),NEWORD(NSTINO),Z(NORB,NEL+1)
      DIMENSION ISGSTI(NSTINI),ISGSTO(NSTINO)
*.Output
      DIMENSION TI(LDIM,NSTINI),TTO(LDIM,NSTINI)
*.Scratch
      DIMENSION STRIN2(500)
*
      NTEST0 =  1
      NTEST = MAX(IPRNT,NTEST0)
      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*)  ' =============== '
        WRITE(6,*)  ' ANNSTR speaking '
        WRITE(6,*)  ' =============== '
        WRITE(6,*)
         WRITE(6,*) ' Number of input electrons ', NEL
      END IF
      LUOUT = 6
*
      DO 1000 ISTRIN = 1,NSTINI
        DO 100 IEL = 1, NEL                  
*. String with electron removed
          DO JEL = 1, IEL-1
           STRIN2(JEL) = STRING(JEL,ISTRIN)
          END DO
          DO JEL = IEL+1, NEL
           STRIN2(JEL-1) = STRING(JEL,ISTRIN)
          END DO
          JSTRIN = ISTRNM(STRIN2,NOCOB,NEL-1,Z,NEWORD,1)
C?        WRITE(6,*) ' anni-string and number '
C?        CALL IWRTMA(STRIN2,1,NEL-1,1,NEL-1)
C?        WRITE(6,*) ' JSTRIN = ', JSTRIN
*
          IORBABS = STRING(IEL,ISTRIN)
          IORB = STRING(IEL,ISTRIN)-IORBOF+1
          IF(IEC.EQ.1) THEN
            IROW = IORB
          ELSE
            IROW = IEL
          END IF
*
          TI(IROW,ISTRIN ) = -IORBABS
C         TI(IROW,ISTRIN ) = -IORB
          TTO(IROW,ISTRIN) = JSTRIN
          IIISGN = (-1)**(IEL-1)
          IF(LSGSTR.NE.0)
     &    IIISGN = IIISGN*ISGSTO(JSTRIN)*ISGSTI(ISTRIN)
          IF(IIISGN .EQ. -1 )
     &    TTO(IROW,ISTRIN) = - TTO(IROW,ISTRIN)
  100   CONTINUE
*
 1000 CONTINUE
*
      IF ( NTEST .GE. 20) THEN
        MAXPR = 60
        NPR = MIN(NSTINI,MAXPR)
        WRITE(LUOUT,*) ' Output from ANNSTR : '
        WRITE(LUOUT,*) '==================='
*
        WRITE(6,*)
        WRITE(LUOUT,*) ' Strings with an electron added or removed'
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' New strings.. ',
     &     (TTO(I,ISTRIN),I = 1,LDIM)
        END DO
        DO ISTRIN = 1, NPR
           WRITE(6,'(2X,A,I4,A,/,(10I5))')
     &     'String..',ISTRIN,' orbitals added or removed ' ,
     &     (TI(I,ISTRIN),I = 1,LDIM)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE EXCCLS2(NCLS,IACTIN,IACTOUT,IEXC,
     &                   IBASSPC_MX,IBASSPC)
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* If IBASSPC_MX .ne. 0, atmost basespaces belonging to this 
*                       space is included
*
* Master routine
*
* Jeppe Olsen, Jan. 1999 - ved siden af ditte, paa MAS efter
*              hendes rygoperation
*
*. Last modification; Nov. 3, 2012; Jeppe Olsen, Aligning with current code
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION IACTIN(NCLS)
*. Output
      DIMENSION IACTOUT(NCLS)
*. From the common blocks
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
*
      CALL EXCCLS2_S(NGAS,int_mb(KIOCCLS),NCLS,IACTIN,IACTOUT,IEXC,
     &               IBASSPC_MX,IBASSPC)
*
      RETURN 
      END
*
      SUBROUTINE EXCCLS2_S(NGAS,ICLS,NCLS,IACTIN,IACTOUT,IEXC,
     &               IBASSPC_MX,IBASSPC)      
*
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* If IBASSPC_MX .ne. 0, atmost basespaces belonging to this 
*                       space is included
* Slave routine
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS(NGAS,NCLS)
      DIMENSION IACTIN(NCLS)
      DIMENSION IBASSPC(*)
*. Output
      DIMENSION IACTOUT(NCLS)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from EXCCLS2_S '
        WRITE(6,*) ' ======================'
      END IF
      IF(NTEST.GE.1000) THEN
        write(6,*) ' ICLS in EXCCLS '
        call iwrtma(icls,ngas,ncls,ngas,ncls)
      END IF
 
      IZERO = 0
      CALL ISETVC(IACTOUT,IZERO,NCLS)
      DO ICLSIN = 1, NCLS
        IF(IACTIN(ICLSIN).NE.0) THEN
        DO ICLSOUT = 1, NCLS
*. Number of anihilations and creations required to connect classes
          NANNI = 0
          NCREA = 0
          DO IGAS = 1, NGAS
*           
            IDIF = ICLS(IGAS,ICLSOUT)-ICLS(IGAS,ICLSIN)
            IF(IDIF .GT. 0 ) THEN                               
              NCREA = NCREA + IDIF  
            ELSE IF (IDIF .LT. 0 ) THEN
              NANNI = NANNI - IDIF
            END IF
          END DO
* 
          IF(NCREA.LE.IEXC) THEN
            IF(IBASSPC_MX.EQ.0.OR.IBASSPC(ICLSOUT).LE.IBASSPC_MX) THEN
              IACTOUT(ICLSOUT) = 1
            END IF
          END IF
*
        END DO
        END IF
      END DO
*
      NACTOUT = 0
      DO ICLSOUT = 1, NCLS
        IF(IACTOUT(ICLSOUT).GT.0) NACTOUT = NACTOUT + 1
      END DO
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Output from EXCCLS '
         WRITE(6,*) ' ==================='
         WRITE(6,*) 
         WRITE(6,*) ' Number of active output classes ',NACTOUT
         WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Active output classes '
         WRITE(6,*) ' ======================'
         DO I = 1, NCLS
           IF(IACTOUT(I).NE.0) WRITE(6,*) I
         END DO
C        CALL IWRTMA(IACTOUT,NCLS,1,NCLS,1)
      END IF
*
      RETURN
      END
*
      SUBROUTINE EXCCLS_S(NGAS,ICLS,NCLS,IACTIN,IACTOUT,IEXC)      
*
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* Slave routine
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS(NGAS,NCLS)
      DIMENSION IACTIN(NCLS)
*. Output
      DIMENSION IACTOUT(NCLS)
*
C     write(6,*) ' ICLS in EXCCLS '
C     call iwrtma(icls,ngas,ncls,ngas,ncls)
 
      IZERO = 0
      CALL ISETVC(IACTOUT,IZERO,NCLS)
      DO ICLSIN = 1, NCLS
        IF(IACTIN(ICLSIN).NE.0) THEN
        DO ICLSOUT = 1, NCLS
*. Number of anihilations and creations required to connect classes
          NANNI = 0
          NCREA = 0
          DO IGAS = 1, NGAS
*           
            IDIF = ICLS(IGAS,ICLSOUT)-ICLS(IGAS,ICLSIN)
            IF(IDIF .GT. 0 ) THEN                               
              NCREA = NCREA + IDIF  
            ELSE IF (IDIF .LT. 0 ) THEN
              NANNI = NANNI - IDIF
            END IF
          END DO
* 
          IF(NCREA.LE.IEXC) THEN
            IACTOUT(ICLSOUT) = 1
          END IF
*
        END DO
        END IF
      END DO
*
      NACTOUT = 0
      DO ICLSOUT = 1, NCLS
        IF(IACTOUT(ICLSOUT).GT.0) NACTOUT = NACTOUT + 1
      END DO
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
         WRITE(6,*) ' Output from EXCCLS '
         WRITE(6,*) ' ==================='
         WRITE(6,*) 
         WRITE(6,*) ' Number of active output classes ',NACTOUT
         WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Active output classes '
         CALL IWRTMA(IACTOUT,NCLS,1,NCLS,1)
      END IF
*
      RETURN
      END
      SUBROUTINE FIND_ACTIVE_BLOCKS_IN_VECTOR(
     &           VEC,BLK_A,NBLOCK,IBLOCK)
*
* Find active (non-vanishing) blocks in vector VEC
*
*. Jeppe Olsen, June 2010
*
C BERT NEW: ONLY LOOK AT BLOCK NODE HAS, LOCAL OPERATION
C ONLY DO LOCAL BLOCK BUT ACCUMULATE RESULT AND REPLICATE 
C ACROSS NODES SO THEY KNOW
C           IF ACTIVE BLK_A = 1
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      REAL*8
     &INPROD
*. Input
      integer VEC,ilo,ihi,idum,ivec,ii,len
      DIMENSION IBLOCK(8,NBLOCK)
CNW   DIMENSION VEC(*), IBLOCK(8,NBLOCK)
*. Output
      DIMENSION BLK_A(NBLOCK)
*
      CALL LUCIAQENTER('ACTBL')
      call ga_distribution(VEC,ga_nodeid(),ilo,ihi,idum,idum)
      call ga_access(VEC, ilo, ihi, 1, 1, ivec, idum)
C?    WRITE(6,*) ' FIND_ACTIVE... ', NBLOCK
      IOFF = 1
      DO IBLK = 1, NBLOCK
        LEN = IBLOCK(8,IBLK)
C?      WRITE(6,*) 'ACT  IBLK, LEN = ', IBLK,LEN,IOFF
        BLK_A(IBLK) = 0.0D0
        if ((ioff.ge.ilo).and.(ioff.le.ihi)) then
         XBLK = INPROD(dbl_mb(ivec+ioff-ilo),dbl_mb(ivec+ioff-ilo),LEN)
cVOG    XBLK = ga_ddot_patch(VEC,'N',IOFF,IOFF+LEN-1,1,1,
cVOG &                       VEC,'N',IOFF,IOFF+LEN-1,1,1)
         IF(XBLK.NE.0.0D0) BLK_A(IBLK) = 1.0D0
        end if
CNW     XBLK = INPROD(VEC(IOFF),VEC(IOFF),LEN)
        IOFF = IOFF + LEN
      END DO
*
      call ga_release(VEC,ilo,ihi,1,1)
      call ga_sync()
      call ga_dgop(1,BLK_A,NBLOCK,'+')
*
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Array of active blocks: '
        CALL WRTMAT(BLK_A,1,NBLOCK,1,NBLOCK)
        do i = 1, nblock
         if (BLK_A(i).eq.1.0d0) write(*,*) 'active blk', i
        end do
      END IF
*
      CALL LUCIAQEXIT('ACTBL')
      RETURN
      END
      SUBROUTINE FIND_ACTIVE_BLOCKS 
     &           (LUIN,LBLK,BLK_A,SEGMNT)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION BLK_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*
      CALL REWINE(LUIN  ,LBLK)
*
      IBLK = 0
      NBLK_A = 0
*. Loop over blocks
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
C         CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN  
           NBLK_A = NBLK_A + 1
           BLK_A(IBLK) = 1.0D0
          ELSE
           BLK_A(IBLK) = 0.0D0
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
      NBLK =  IBLK-1 
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) 
     &  ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active blocks '
        CALL WRTMAT(BLK_A,1,NBLK,1,NBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE FIND_ACTIVE_CLASSES
     &           (LUIN,LBLK,I_BLK_TO_CLS,ICLS_A,NCLS,SEGMNT)
*
*. Find the active (nonvanishing classes ) on LUIN)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION I_BLK_TO_CLS(*)
*. Output
      DIMENSION ICLS_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*
      CALL REWINE(LUIN  ,LBLK)
      IZERO = 0
      CALL ISETVC(ICLS_A,IZERO,NCLS)
*
      IBLK = 0
      NCLS_A = 0
*. Loop over blocks
C?      write(6,*) ' ZAP_BLOCK_VEC :  LBLK = ', LBLK
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
C         CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN  
           IF(ICLS_A(I_BLK_TO_CLS(IBLK)).EQ.0) THEN
             NCLS_A = NCLS_A + 1
           END IF
           ICLS_A(I_BLK_TO_CLS(IBLK)) = 1
*
           IF(I_BLK_TO_CLS(IBLK).EQ.0) THEN
             WRITE(6,*) ' Problem in FIND_ACTIVE_CLASSES : '
             WRITE(6,*) ' IBLK, I_BLK_TO_CLS(IBLK) ', 
     &                    IBLK, I_BLK_TO_CLS(IBLK) 
           END IF
*
C          write(6,*) ' Active block and class ',IBLK,
C    &     I_BLK_TO_CLS(IBLK)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' FIND_A.... Number of active classes ',NCLS_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active classes '
        CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
      END IF
*
      RETURN
      END
      SUBROUTINE REPART_CIV(IBATCH,NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                  MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLKS,
     &                  LBLOCK_VEC)
*
* A CI vector is defined through IBATCH (generated by PART_CIV)
*
*. Divide into batches, length atmost MXLNG, so only blocks
*  that are flagged active by ICLS_A are included
*
*. Output 
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : Inactive blocks are flagged by setting the first element
*          negative ( -1 * the original value )
*
* Input
*
* IBATCH : 
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta sym
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta 
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expandend form                   
*   IBATCH(8,*) : Length of block, packed form 
*
*
*
* Jeppe Olsen, June 1997     
*
* Last modification; Nov. 3, 2012; Jeppe Olsen; I1BATCH
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER IBATCH(8,*)
      INTEGER ICLS_A(*)
*. General input
       INTEGER IBLK_TO_CLS(*)
*. Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER LBLOCK_VEC(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     REPART_CIV      '
        WRITE(6,*) ' =================='
        WRITE(6,*)  
        WRITE(6,*)
C       WRITE(6,*) ' Active classes '
C       CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
C       write(6,*) ' MXLNG = ',MXLNG
      END IF
*
      NBATCH = 0
      LENGTH = 0 
      LENGTHP= 0 
      IIBLOCK = 0
      I1BATCH(1) = 1
*. Loop over blocks in batch
      DO 1000 IBLOCK = 1, NBLKS  
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IF(ICLS_A(IBLK_TO_CLS(IBLOCK)).EQ.0) THEN
*. Block is inactive 
          IBATCH(1,IBLOCK) = - ABS(IBATCH(1,IBLOCK))
          IIBLOCK = IIBLOCK + 1                           
          LBLOCK_VEC(IBLOCK) = -IBATCH(8,IBLOCK)
        ELSE
*. Block belongs to the active !
          IBATCH(1,IBLOCK) =   ABS(IBATCH(1,IBLOCK))
          LBLOCK =  IBATCH(7,IBLOCK) 
          LBLOCKP =  IBATCH(8,IBLOCK) 
          LBLOCK_VEC(IBLOCK) = LBLOCKP
          IF(LENGTH+LBLOCK.LE.MXLNG) THEN
            LENGTH = LENGTH + LBLOCK
            LENGTHP= LENGTHP+ LBLOCKP
            IIBLOCK = IIBLOCK + 1                           
          ELSE IF(LENGTH+LBLOCK.GT.MXLNG) THEN
*. This batch was finished by previous block, goto next batch
            NBATCH = NBATCH + 1
            LEBATCH(NBATCH) = LENGTHP
            LBATCH (NBATCH)  = IIBLOCK
            I1BATCH(NBATCH+1) = IBLOCK 
*. Current block is first block in new batch
            IIBLOCK = 1
            IBATCH(5,IBLOCK) = 1
            IBATCH(6,IBLOCK) = 1
            LENGTHP = LBLOCKP
            LENGTH  = LBLOCK 
          END IF
        END IF
 1000 CONTINUE
*. Final batch
      IF( LENGTH .NE. 0 ) THEN
        NBATCH = NBATCH + 1
        LBATCH(NBATCH) = IIBLOCK
        LEBATCH(NBATCH) = LENGTHP
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' REPART.. : Number of batches ', NBATCH                           
      END IF
      IF(NTEST.GE.100) THEN
        IOFF = 1
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lengths of each block '
          DO IBLOCK = IOFF, IOFF+LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
          IOFF = IOFF + LBATCH(JBATCH)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE FRMDSCN3(VEC,NREC,LBLK,LU,NO_ZEROING,I_AM_ZERO,
     &                    LBLOCK)
*
* Read  VEC as multiple record file, NREC records read
* If NO_ZEROING.EQ.0 then zero blocks are not
* set to zero;  a 1 is instead flagged in the relevant block
* of I_AM_ZERO
*
* IF LBLOCK(IBLK) is less than zero, then no space for this
* block is made
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION LBLOCK(*)
*. Output
      DIMENSION VEC(*),I_AM_ZERO(*)
*
      IOFF = 1
C?    write(6,*) ' FRMDSCN3: Number of records to be read', NREC
      DO IREC = 1, NREC
        CALL IFRMDS(LREC,1,LBLK,LU)
        IF(LBLOCK(IREC).GE.0) THEN
          CALL FRMDSC2(VEC(IOFF),LREC,LBLK,LU,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          I_AM_ZERO(IREC) = IMZERO
          IOFF = IOFF + LREC
        ELSE
          I_AM_ZERO(IREC) = 1
          CALL SKPRCD2(LBLK,-1,LU)
        END IF

      END DO
*
      RETURN
      END
       SUBROUTINE GETOBS_LUCIA(NAOS_ENV,NMOS_ENV)
*
* Obtain info on orbital dimensions from LU91 - LUCIA format
*
* Jeppe Olsen, Feb. 98, AOLABELS added, june 2010
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
*. Output
      INTEGER NMOS_ENV(*),NAOS_ENV(*)
      CHARACTER*4 AO_CENT, AO_TYPE
      COMMON/AOLABELS/AO_CENT(MXPORB),AO_TYPE(MXPORB)
*
      LUH = 91
      CALL REWINO(LUH)
*.
      READ(LUH,*) NSMOB
*. 
      READ(LUH,*) (NMOS_ENV(ISM),ISM=1, NSMOB)
*
      READ(LUH,*) (NAOS_ENV(ISM),ISM=1, NSMOB)
*Skip MO file
      READ(LUH,*) NCMOAO
      IF(NOMOFL.EQ.0) READ(LUH,*) (X, I=1, NCMOAO)
*. And read AO labels
      IF(NOMOFL.EQ.0) THEN
        NAO_TOT = IELSUM(NAOS_ENV,NSMOB)
        READ(LUH,'(20A4)') (AO_CENT(IAO),IAO = 1, NAO_TOT)
        READ(LUH,'(20A4)') (AO_TYPE(IAO),IAO = 1, NAO_TOT)
      END IF
*
      RETURN
      END 
*
* Obtain property integrals with LABEL LABEL from LU91,
* LUCIA format
*
* Jeppe Olsen, Feb.98 

      SUBROUTINE GET_H1AO(LABEL,H1AO,IHSM,NBAS,IPERMSM)
*
* Obtain 1 electron integrals with label LABEL                  
*
* Jeppe Olsen, Feb.98
*              IPERMSM added, Feb 2000
*
*. Note: AO integrals are pr read from DALTON file, even if 
*. environment is LUCIA (quick-fix, June 2012)
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      CHARACTER*8 LABEL
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT_H1A')
*
      IF(ENVIRO(1:6).EQ.'DALTON'.OR.ENVIRO(1:5).EQ.'LUCIA') THEN 
        LSCR = NTOOB**2
        CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'GTH1SC') !done
        CALL GET_H1AO_DALTON(LABEL,H1AO,IHSM,dbl_mb(KLSCR),NBAS,NSMOB,
     &                       IPERMSM)
C            GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS,NSM,IPERMSM)
C!!   ELSE IF (ENVIRO(1:5).EQ.'LUCIA') THEN 
C!!     LU91 = 91
C!!     CALL GET_H1AO_LUCIA(LABEL,H1AO,LU91)
      ELSE IF (ENVIRO(1:4).EQ.'QDOT') THEN
        WRITE(6,*) ' GET_H1AO_QDOT not coded '
        STOP       ' GET_H1AO_QDOT not coded '
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GT_H1A')
*
      RETURN
      END
      SUBROUTINE GET_CMOAO_LUCIA(CMO,NMOS,NAOS,LUH)
*
* Obtain CMOAO expansion matrix from LUCIA formatted file LUH
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
*. Input
      INTEGER NMOS(*),NAOS(*)
*. Output  
      DIMENSION CMO(*)
*
* Structure of file 
* 1 : Number of syms
* 2 : NMO's per sym
* 3 : NAO's per SYM
* 4 : Number of elements in CMOAO
* Note : CMOAO and property integrals written in form 
*     given by ONEEL_MAT_DISC
*
* Jeppe Olsen, Feb. 98
*
C?    WRITE(6,*)  ' GET_CMOAO_LUCIA, LUH = ', LUH
      CALL REWINO(LUH)
*. skip Number of orbital symmetries
      READ(LUH,*) 
*. skip Number of MO's per symmetry
      READ(LUH,*) 
*. skip Number of AO's per symmetry
      READ(LUH,*) 
*. Check that CMO-AO expansion was written to disc
      READ(LUH,*) LCMOAO
      IF(LCMOAO.EQ.0) THEN
        WRITE(6,*) 
     &  ' GET_CMOAO_LUCIA : NO MOAO expansion matrix on FILE 91'
        STOP 'GET_CMOAO_LUCIA : NO MOAO '
      END IF
*. read CMO-AO expansion matrix
      CALL ONEEL_MAT_DISC(CMO,1,NSMOB,NAOS,NMOS,LUH,1)
C          ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MO-AO transformation read in '
        CALL PRHONE(CMO,NMOS,1,NSMOB,0)
C            PRHONE(C,NFUNC,M,NSM,IPACK)
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_CMOAO_NWCHEM(CMO)
*
* Obtain CMOAO expansion matrix from LUCIA formatted file LUH
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "lucinp.inc"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "sym.fh"
#include "rtdb.fh"
#include "nwchem_data.inc"
      logical movecs_gather
      external movecs_gather
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
*. Output  
      DIMENSION CMO(*)
      double precision ctmp(mxporb)
      integer g_moreo
      logical test
*
      g_moreo = ga_create_atom_blocked(luciageom,luciabasis,
     &                                 'LUCIA movecs')
      if(.not.movecs_gather(luciascfmo,1,luciamo(1),reord_nw,ctmp,ctmp,
     &     g_moreo)) call errquit('get_cmoao_nwchem failed',0,DISK_ERR)
c     call ga_clean(g_moreo,1.0d-12)
      icounter=1
      do i = 1, luciamo(1)
         call ga_get(g_moreo, 1, luciaao, i, i, ctmp, 1)
         nir = kirs_nw(reord_nw(i))
         do j = 1, luciaao
            if (nwchem_id(j,nir).gt.0) then
               cmo(icounter) = ctmp(j)*nwchem_sc(j,nir)
               icounter = icounter + 1
             endif
         enddo
      enddo
      if(.not.ga_destroy(g_moreo)) call 
     &        errquit('get_cmoao_nwchem destroy',0,GA_ERR)

      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MO-AO transformation read in '
        CALL PRHONE(CMO,NMOS,1,NSMOB,0)
C            PRHONE(C,NFUNC,M,NSM,IPACK)
      END IF
*
      RETURN
      END 
      subroutine ga_clean(g_input,clean_thresh)
      implicit real*8(a-h,o-z)
      integer g_input
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      call ga_distribution(g_input,ga_nodeid(),ilo,ihi,jlo,jhi)
      if (ihi.gt.0) then
      call ga_access(g_input, ilo, ihi, jlo, jhi, iarray, idum)
      ilen = (ihi-ilo+1)*(jhi-jlo+1)-1
      do i = 0, ilen
         if (dabs(dbl_mb(iarray+i)).le.clean_thresh)
     &      dbl_mb(iarray+i) = 0.0d0
      enddo
      call ga_release(g_input,ilo,ihi,jlo,jhi)
      endif
*
      return
      end
      FUNCTION NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
*
* Length of one-electron matrix with symmetry IHSM
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NRPSM(*),NCPSM(*)
*
      LENGTH = 0
      DO IRSM = 1, NSM
        CALL SYMCOM(2,5,IRSM,ICSM,IHSM)
C       WRITE(6,*) ' IHSM,IRSM,ICSM', IHSM,IRSM,ICSM
        IF(IPACK.EQ.0.OR.(IPACK.EQ.1.AND.IRSM.GT.ICSM)) THEN
          LENGTH = LENGTH + NRPSM(IRSM)*NCPSM(ICSM)
        ELSE IF(IPACK.EQ.1.AND.IRSM.EQ.ICSM) THEN
          LENGTH = LENGTH + NRPSM(IRSM)*(NRPSM(IRSM)+1)/2
        END IF
      END DO
*
      NDIM_1EL_MAT = LENGTH
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' 2 Dim array , sym and length : ', IHSM,LENGTH
      END IF
*
      RETURN
      END
      SUBROUTINE DUMP_1EL_INFO(LUH)
*
*. Dump one-electron information  on file LUH, LUCIA format
*. It is assumed that overlap matrix over AO's are in WORK(KSAO)
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
*
      CHARACTER*4 AO_CENT, AO_TYPE
      COMMON/AOLABELS/AO_CENT(MXPORB),AO_TYPE(MXPORB)
*
      CHARACTER*1 XYZ(3)
      DATA XYZ/'X','Y','Z'/
      CHARACTER*8 LABEL
*
* Structure of file 
* 1 : Number of syms
* 2 : NMO's per sym
* 3 : NAO's per SYM
* 4 : Number of elements in CMOAO
* 5 : CMOAO-expansion matrix (in symmetry packed form) (in KMOAO)
* 6 : Center and type of the AO's
* 7 : Number of property AO lists 
*     Loop over number of properties
*     Label, offset and length of each proprty list
*
*     Property integrals for prop1,prop2 ...
*
* Note : CMOAO and property integrals written in form 
*     given by ONEEL_MAT_DISC
*
* Jeppe Olsen, Feb. 98
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DUMP_1')
C?    WRITE(6,*) ' Entering DUMP_1EL_INFO '
*. A scratch block of length NTOOB ** 2
      LSCR = NTOOB ** 2
      CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'DUMPSC') !done
*. An additional scratch block
      CALL MEMMAN(KLSCR2,3*LSCR,'ADDL  ',2,'DMPSC2') !done
*
      CALL REWINO(LUH)
*. Number of orbital symmetries
      WRITE(LUH,*) NSMOB
*. Number of MO's per symmetry
      WRITE(LUH,*) (NMOS_ENV(ISM),ISM=1,NSMOB)
*. Number of AO's per symmetry
      WRITE(LUH,*) (NAOS_ENV(ISM),ISM=1,NSMOB)
*. Length of CMO-AO expansion
      IF(NOMOFL.EQ.0) THEN
        LENGTH = NDIM_1EL_MAT(1,NAOS_ENV,NMOS_ENV,NSMOB,0)
        WRITE(LUH,*) LENGTH
        DO IJ = 1, LENGTH
          WRITE(LUH,'(E22.15)') WORK(KMOAO-1+IJ)
        END DO
      ELSE
        LENGTH = 0
        WRITE(LUH,*) LENGTH
      END IF
*. Center and type of the various AO's
      NAO_TOT = IELSUM(NAOS_ENV,NSMOB)
      WRITE(LUH,'(20A4)') (AO_CENT(IAO),IAO = 1, NAO_TOT)
      WRITE(LUH,'(20A4)') (AO_TYPE(IAO),IAO = 1, NAO_TOT)
C     WRITE(6,'(20A4)') (AO_CENT(IAO),IAO = 1, NAO_TOT)
C     WRITE(6,'(20A4)') (AO_TYPE(IAO),IAO = 1, NAO_TOT)
*. Total number of properties ( 3 for each rank1, 6 for each rank 2)
      NPROP_COM = 0
      DO IPROP = 1, NPROP
C             GET_PROP_RANK(PROPER,IRANK)
         CALL GET_PROP_RANK(PROPER(IPROP),IRANK)
         IF(IRANK.EQ.0) THEN
           NPROP_COM = NPROP_COM + 1
         ELSE IF (IRANK.EQ.1) THEN
           NPROP_COM =  NPROP_COM + 3
         ELSE IF (IRANK.EQ.2) THEN
           NPROP_COM = NPROP_COM + 6
         END IF
       END DO
       WRITE(LUH,*) NPROP_COM
       IOFF = 1
       DO IPROP = 1, NPROP
         CALL GET_PROP_RANK(PROPER(IPROP),IRANK)
         LABEL = PROPER(IPROP)//'  '
         IF(IRANK.EQ.0) THEN
*. A single  component, total symmetric
           LENGTH = NDIM_1EL_MAT(1,NAOS_ENV,NMOS_ENV,NSMOB,0)
*. Label, offset, length
           WRITE(LUH,'(A,I6,I6)') LABEL,IOFF,LENGTH
           IOFF = IOFF + LENGTH
         ELSE IF(IRANK.EQ.1) THEN
           WRITE(6,'(A)') ' LABEL 1 = ', LABEL
           DO ICOMP = 1, 3
             LABEL = PROPER(IPROP)//'  '
             IF(LABEL(1:6).EQ.'DIPOLE') THEN
               LABEL =XYZ(ICOMP)//'DIPLEN '            
             ELSE
               LABEL =XYZ(ICOMP)//PROPER(IPROP)//' '
             END IF
           WRITE(6,'(A)') ' LABEL 2 = ', LABEL
             CALL SYM_FOR_OP(LABEL,IXYZSYM,IOPSM)
             LENGTH = NDIM_1EL_MAT(IOPSM,NAOS_ENV,NMOS_ENV,NSMOB,0)
             WRITE(LUH,'(A,I6,I6)') LABEL,IOFF,LENGTH
             IOFF = IOFF + LENGTH
           END DO
         ELSE IF (IRANK.EQ.2) THEN
           DO ICOMP = 1, 3
             DO JCOMP = 1, ICOMP
               LABEL = XYZ(JCOMP)//XYZ(ICOMP)//PROPER(IPROP)
               CALL SYM_FOR_OP(LABEL,IXYZSYM,IOPSM)
               LENGTH = NDIM_1EL_MAT(IOPSM,NAOS_ENV,NMOS_ENV,NSMOB,0)
               WRITE(LUH,'(A,I6,I6)') LABEL,IOFF,LENGTH
               IOFF = IOFF + LENGTH
             END DO
           END DO
         END IF
       END DO
*. Fetch integrals and then : Dump them           
       DO IPROP = 1, NPROP
         CALL GET_PROP_RANK(PROPER(IPROP),IRANK)
         CALL GET_PROP_PERMSM(LABEL,IPERMSM)
         LABEL = PROPER(IPROP)//'  '
         IF(IRANK.EQ.0) THEN
CE         CALL GET_H1AO(LABEL,WORK(KLSCR),1,NAOS_ENV)
           CALL GET_PROPINT(dbl_mb(KLSCR),1,LABEL,dbl_mb(KLSCR2),
     &                      NMOS_ENV,NAOS_ENV,NSMOB,0,IPERMSM) 
C               GET_PROPINT(H,IHSM,LABEL,SCR,NMO,NBAS,NSM,ILOW)
           CALL ONEEL_MAT_DISC(dbl_mb(KLSCR),1,NSMOB,NAOS_ENV,
     &          NAOS_ENV,LUH,2)
         ELSE IF(IRANK.EQ.1) THEN
           DO ICOMP = 1, 3
             LABEL = PROPER(IPROP)//'  '
             IF(LABEL(1:6).EQ.'DIPOLE') THEN
               LABEL =XYZ(ICOMP)//'DIPLEN '            
             ELSE
               LABEL =XYZ(ICOMP)//PROPER(IPROP)//' '
             END IF
             CALL SYM_FOR_OP(LABEL,IXYZSYM,IOPSM)
CE           CALL GET_H1AO(LABEL,WORK(KLSCR),IOPSM,NAOS_ENV)
             CALL GET_PROPINT(dbl_mb(KLSCR),IOPSM,LABEL,dbl_mb(KLSCR2),
     &                        NMOS_ENV,NAOS_ENV,NSMOB,0,IPERMSM)
C            GET_PROPINT(H,IHSM,LABEL,SCR,NMO,NBAS,NSM,ILOW)
             CALL ONEEL_MAT_DISC(dbl_mb(KLSCR),IOPSM,NSMOB,
     &            NAOS_ENV,NAOS_ENV,LUH,2)
           END DO
         ELSE IF (IRANK.EQ.2) THEN
           DO ICOMP = 1, 3
             DO JCOMP = 1, ICOMP
               LABEL = XYZ(JCOMP)//XYZ(ICOMP)//PROPER(IPROP)
               CALL SYM_FOR_OP(LABEL,IXYZSYM,IOPSM)
C     GET_PROPINT(H,IHSM,LABEL,SCR,NMO,NBAS,NSM,ILOW)
           CALL GET_PROPINT(dbl_mb(KLSCR),IOPSM,LABEL,dbl_mb(KLSCR2),
     &                      NMOS_ENV,NAOS_ENV,NSMOB,0,IPERMSM)
CE             CALL GET_H1AO(LABEL,WORK(KLSCR),IOPSM,NAOS_ENV)
               CALL ONEEL_MAT_DISC(dbl_mb(KLSCR),IOPSM,NSMOB,
     &              NAOS_ENV,NAOS_ENV,LUH,2)
             END DO
           END DO
         END IF
       END DO
* Overlap matrix over atomic orbitals
      LENGTH = NDIM_1EL_MAT(1,NAOS_ENV,NMOS_ENV,NSMOB,1)
      WRITE(LUH,*) LENGTH
      DO IJ = 1, LENGTH
        WRITE(LUH,'(E22.15)') WORK(KSAO-1+IJ)
      END DO
*
*. Rewind to empty buffer
      CALL REWINO(LUH)
C     WRITE(6,*) ' Enforced stop after end of DUMP_1EL'
C     STOP       ' Enforced stop after end of DUMP_1EL'
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DUMP_1')
*
      RETURN
      END
      SUBROUTINE ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
* Transfer one-electron matrix H between memory and disc file in
* LUCIA format 
*
* IFT = 1 => From disc ( read)
* IFT = 2 => To   disc (write)
*
* Note : File LUH is supposed to be at start of correct integral block
*
* Jeppe Olsen, Feb. 98
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER NRPSM(NSM),NCPSM(*)
      DIMENSION H(*)
*
*. Order of integrals are
*
*     Loop over Symmetry of row index => Symmetry of column  index
*      Loop over columns in symmetry block
*        Loop over rows in symmetry block
*        End of loop over rows in symmetry block
*      End of Loop over columns in symmetry block
*     End of loop over symmetry of row index
*
* Each symmetry block is thus given in complete form
* Note all integrals are in a single record 
*
*. Length of list
C              NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,IPACK)
      LENGTH = NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,0)
*. and read/write
C     WRITE(6,*) ' ONEEL, IFT, LUH = ', IFT, LUH
      IF(IFT.EQ.1) THEN
C       DO IJ = 1, LENGTH
          READ(LUH,*) (H(IJ),IJ=1, LENGTH)
C       END DO
      ELSE IF (IFT.EQ.2) THEN
        DO IJ = 1, LENGTH
          WRITE(LUH,'(E22.15)') H(IJ)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE GET_PROP_RANK(PROPER,IRANK)
*
* Obtain rank for property PROPERTY (gives as CHAR*6 as usual)
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*6 PROPER
*
      IF(PROPER.EQ.'DIPOLE'.OR.
     &   PROPER.EQ.'DIPLEN') THEN
        IRANK = 1
      ELSE IF (PROPER.EQ.'ANGMOM') THEN
        IRANK = 1
      ELSE IF(PROPER.EQ.'THETA ' .OR.
     &        PROPER.EQ.'QUADRU' .OR.
     &        PROPER.EQ.'SECMOM' .OR.
     &        PROPER(1:3).EQ.'EFG' ) THEN
        IRANK  = 2
      ELSE
        WRITE(6,'(A,A)') ' Unknown operator ',PROPER
        IRANK = -1
      END IF
*
      NTEST = 0
      IF(NTEST.GE.5) THEN
        WRITE(6,'(A,A,I3)') ' Property and rank : ', PROPER,IRANK
      END IF
*
      RETURN
      END
*
      SUBROUTINE PROP_NATORB(HDIAG,OCCNUM,NTOOBS,NSMOB)
*
* Analyze property in natural orbital representation
*
* Jeppe Olsen, Feb 98
*
* Input :
* ========
*
* HDIAG : Diagonal of 1-e integrals over nat orbs symmetry order
* OCCNUM: Occupation numbers
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION HDIAG(*),OCCNUM(*)
*
      DIMENSION NTOOBS(*)
*
      SUM = 0.0D0
      WRITE(6,'(A)')
     & ' Orb   Sym    Integral   Occ number   Contribution    Sum '
      WRITE(6,'(A)')
     & ' ========================================================='
      DO ISYM = 1, NSMOB
      WRITE(6,*)
        IF(ISYM.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISYM-1)
        END IF
        DO IORB = 1, NTOOBS(ISYM) 
          IORBEF = IORB+IOFF-1
          CONT = OCCNUM(IORBEF)*HDIAG(IORBEF)
          SUM = SUM + CONT
          WRITE(6,'(1H ,I4, I4,4E13.5)') 
     &    IORB,ISYM,HDIAG(IORBEF),OCCNUM(IORBEF),CONT,SUM
        END DO
      END DO
*
      RETURN 
      END
      SUBROUTINE GET_DIAG_BLOC_MAT(A,ADIAG,NBLOCK,LBLOCK,ISYM)
*
* Obtain diagonal elements from symmetry blocked matrix
*
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION A(*)
      INTEGER LBLOCK(*)
*. Output 
      DIMENSION ADIAG(*)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFF = 1
          LOFF = 1
        ELSE
          IF(ISYM.EQ.1) THEN
            IOFF = IOFF + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
          ELSE
            IOFF = IOFF + LBLOCK(IBLOCK-1)** 2                     
          END IF
          LOFF = LOFF + LBLOCK(IBLOCK-1)
        END IF
*
        L = LBLOCK(IBLOCK)
        CALL COPDIA(A(IOFF),ADIAG(LOFF),L,ISYM)
C            COPDIA(A,VEC,NDIM,IPACK)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        NDIM = IELSUM(LBLOCK,NBLOCK)
        WRITE(6,*) ' output matrix GET_DIAG_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        CALL WRTMAT(ADIAG,1,NDIM,1,NDIM) 
      END IF
*
      RETURN
      END
      SUBROUTINE TRAN_SYM_BLOC_MAT2(AIN,X,NBLOCK,LBLOCK,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LBLOCK(NBLOCK)
*. Output 
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block 
      DIMENSION SCR(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP = 1
         IOFFC = 1
       ELSE
         IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
         IOFFC = IOFFC + LBLOCK(IBLOCK-1)** 2                     
       END IF
       L = LBLOCK(IBLOCK)
       K1 = 1
       K2 = 1 + L **2
*. Unpack block of A
C      TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
         CALL TRIPK3(SCR(K1),AIN(IOFFP),2,L,L,SIGN)
       ELSE
         CALL COPVEC(AIN(IOFFC),SCR(K1),L*L)
       END IF
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL MATML7(SCR(K2),X(IOFFC),SCR(K1),L,L,L,L,L,L,
     &             ZERO,ONE,1)
*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFC),L,L,L,L,L,L,
     &             ZERO,ONE,0)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPK3(SCR(K1),AOUT(IOFFP),1,L,L,SIGN)
       ELSE
         CALL COPVEC(SCR(K1),AOUT(IOFFC),L*L)
       END IF
*
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' Info from TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Transformation matrix '
        CALL APRBLM2(X,LBLOCK,LBLOCK,NBLOCK,0)
        WRITE(6,*) ' Input matrix '
        CALL APRBLM2(AIN,LBLOCK,LBLOCK,NBLOCK,ISYM)
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        CALL APRBLM2(AOUT,LBLOCK,LBLOCK,NBLOCK,ISYM)
      END IF
*
      RETURN
      END
      SUBROUTINE TRAPRP
*
* Calculate one-electron transition properties properties
*
* Jeppe Olsen, June 1997
*              Last modification, Feb. 2000
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cicisp.inc'
*
      REAL*8 INPROD 
      CHARACTER*1 XYZ(3)
      CHARACTER*8 LABEL
      CHARACTER*6 LABEL2
      DATA XYZ/'X','Y','Z'/
*. Common block for communicating with sigma
      INCLUDE 'cands.inc'
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRAPRP')
*
      NCALC = NROOT*NEXCSTATE
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*)
        WRITE(6,*) ' Welcome to TRAPRP'         
        WRITE(6,*) 
        WRITE(6,*) ' Number of transition densities ',NCALC
      END IF
*. A bit of local memory
      LHONE = NTOOB * NTOOB
      WRITE(6,*) ' Max size of one-electron operator',LHONE
      CALL MEMMAN(KLHONE,LHONE,'ADDL  ',2,'KLHONE') 
      CALL MEMMAN(KLCMO,LHONE,'ADDL  ',2,'KLCMO ') 
      CALL MEMMAN(KLRHO1S,LHONE,'ADDL  ',2,'LRHO1S') 
      CALL MEMMAN(KLSCR,4*LHONE,'ADDL  ',2,'KLCMO ') 
      CALL MEMMAN(KLRHOT,NCALC*LHONE,'ADDL  ',2,'RHOT  ')
*. Local memory for density routines
      CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)
*
* 1 : Obtain transition densities
*
*. Right states : current ci vectors on LUC
*. Left states  : states on LUEXC
*
*. Fill /CANDS/
      ICSM = IREFSM
      ISSM = IEXCSYM
*. The space is assumed to be the final space
      ICSPC = NCMBSPC
      ISSPC = NCMBSPC
*. The transition densities
*. stored in KLRHOT as loop over left , Loop over right
      I12 = 1
      ILR = 0
      CALL REWINO(LUEXC)
      DO IL = 1, NEXCSTATE
*. LUEXC => LUHC
        CALL REWINO(LUHC)
        CALL COPVCD(LUEXC,LUHC,WORK(KVEC1),0,-1)
*
        CALL REWINO(LUC)
        DO IR = 1, NROOT
*.LUC => LUSC1
          CALL REWINO(LUSC1)
          CALL COPVCD(LUC,LUSC1,WORK(KVEC1),0,-1)
*
          ILR = ILR + 1
          KLRHOTP = KLRHOT + (ILR-1)*LHONE
          XDUM = 0.0D0
          CALL DENSI2(I12,WORK(KLRHOTP),XRHO2,WORK(KVEC1),WORK(KVEC2),
     &                LUHC,LUSC1,EXPS2,0,XDUM,XDUM,XDUM,XDUM,1)
C              DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR,EXPS2,IDOSRHO1,SRHO1)
*
          WRITE(6,*) ' Transition density beween L and R states ',
     &                 IL, IR
          CALL WRTMAT(WORK(KLRHOTP),NTOOB,NTOOB,NTOOB,NTOOB)
        END DO
      END DO
*
* Transition Properties 
*
      IF(NPROP.GT.0) THEN
*
*. Symmetry of transition densities - and therefore of operators
*
      CALL SYMCOM(3,1,IREFSM,IEXCSYM,ITRASYM) 
      WRITE(6,*) ' Symmetry of transition densities',ITRASYM
*
      IRHO1SM = ITRASYM
      ILR = 0
      DO IL = 1, NEXCSTATE
        DO IR = 1, NROOT
          WRITE(6,*) ' info for IL, IR =', IL,IR
*. Extract symmetry blocks from complete one-electron density 
          ILR = ILR + 1
          KLRHOTP = KLRHOT + (ILR-1)*LHONE
C              REORHO1(RHO1I,RHO1O,IRHO1SM)
          CALL REORHO1(WORK(KLRHOTP),WORK(KLRHO1S),IRHO1SM,1)
*. Number of elements in symmetry blocks of integrals and density
          LRHO1S = 0
          DO ISM = 1, NSMOB
            JSM = MULTD2H(ISM,IRHO1SM)
            LRHO1S = LRHO1S + NTOOBS(ISM)*NTOOBS(JSM)
          END DO
          WRITE(6,*) ' Number of elements in symmetry blocks ',
     &    LRHO1S
*
          DO IPROP =1, NPROP
            WRITE(6,'(A,A)') ' Property to be calculated',
     &                     PROPER(IPROP)
*. one- or two-dimensional tensor ?
            LABEL2 = PROPER(IPROP)
            WRITE(6,'(A,A)') ' LABEL2 =', LABEL2
C           GET_PROP_RANK(PROPER,IRANK)
            CALL GET_PROP_RANK(LABEL2,NRANK)
            CALL GET_PROP_PERMSM(LABEL2,IPERMSM)
C           IF(LABEL2.EQ.'DIPOLE') THEN
C             NRANK = 1
C           ELSE IF(LABEL2.EQ.'THETA ' .OR.
C    &              LABEL2.EQ.'QUADRU' .OR.
C    &              LABEL2(1:3).EQ.'EFG' ) THEN
C             NRANK  = 2
C           ELSE
            IF(NRANK.EQ.-1) THEN
             WRITE(6,'(A,A)') ' Unknown operator ',PROPER(IPROP)
             NRANK = -1
            END IF
            WRITE(6,*) ' Rank of operator ', NRANK
*.  
            IF(NRANK.EQ.1) THEN
              DO ICOMP = 1, 3
                IF(IXYZSYM(ICOMP).EQ.IRHO1SM) THEN
                  WRITE(6,*) ' right symmetry for component',ICOMP
*. Label of integrals on file- currently DALTON FORM !!
                  IF(PROPER(IPROP).EQ.'DIPOLE') THEN
                    LABEL =XYZ(ICOMP)//'DIPLEN '
                  END IF
                  WRITE(6,'(A,A)') ' Label ', LABEL
*. Obtain one-electron integrals
                   CALL GET_PROPINT(WORK(KLHONE),IRHO1SM,LABEL,
     &                  WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,0,IPERMSM)
*. and then : Expectation value
                   EXPVAL = INPROD(WORK(KLRHO1S),WORK(KLHONE),LRHO1S)
                   WRITE(6,'(A,A,E22.14)')
     &             ' Expectation value of ',LABEL,  EXPVAL
                 END IF
               END DO
            ELSE IF(NRANK.EQ.2) THEN
              DO ICOMP = 1,3
                DO JCOMP = 1,ICOMP
                  IF(MULTD2H(IXYZSYM(ICOMP),IXYZSYM(JCOMP))
     &            .EQ.IRHO1SM) THEN
                    WRITE(6,*) ' Right symmetry for components',
     &              ICOMP,JCOMP
*
C                   IF(LABEL2.EQ.'THETA ') THEN
*. Buckinghams traceless quadrupole moment
C                      LABEL=XYZ(JCOMP)//XYZ(ICOMP)//'THETA'
C                   ELSE IF(LABEL2.EQ.'QUADRU') THEN
C                     LABEL=XYZ(JCOMP)//XYZ(ICOMP)//'QUADRU'
C                   ELSE IF (LABEL2(1:3).EQ.'EFG' ) THEN
C
                    LABEL=XYZ(JCOMP)//XYZ(ICOMP)//LABEL2
C                   END IF 
*. Obtain one-electron integrals
                   CALL GET_PROPINT(WORK(KLHONE),IRHO1SM,LABEL,
     &                  WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,0,IPERMSM) 
*. and then : Expectation value
                   EXPVAL = 
     &             INPROD(WORK(KLRHO1S),WORK(KLHONE),LRHO1S)
                   WRITE(6,'(A,A,E22.14)')
     &             ' Expectation value of ',LABEL,  EXPVAL
                  END IF
                END DO
              END DO
            END IF
*
          END DO
* 
        END DO
      END DO
*
      END IF
*     ^ End if Properties should be calculated
 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRAPRP')
      RETURN
      END
*
      SUBROUTINE REORHO1(RHO1I,RHO1O,IRHO1SM,IWAY)
*
* The density matric rho1 is given in complete form
* Extract symmetry blocks with symmetry IRHO1SM
* IWAY = 1: Full matrix to symmetry blocks
* IWAY = 2: Symmetryblocks to full matrix
*
* Jeppe Olsen, June 1997
*              July 2011: Changes from RHO1 restricted to active orbitals
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*. Input or output
      DIMENSION RHO1I(NACOB,NACOB)
*. Output or input
      DIMENSION RHO1O(*)
*
      IF(IWAY.EQ.2) THEN
       ZERO = 0.0D0
       CALL SETVEC(RHO1I,ZERO,NACOB**2)
      END IF
      IMOFF = 1
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IRHO1SM)
        IOFF = IBSO(ISM) + NINOBS(ISM)
        JOFF = IBSO(JSM) + NINOBS(JSM)
*
        NI  = NACOBS(ISM)
        NJ =  NACOBS(JSM)
        DO I = 1, NI
          DO J = 1, NJ
            IP = IREOST(IOFF-1+I)-NINOB
            JP = IREOST(JOFF-1+J)-NINOB
            IF(IWAY.EQ.1) THEN
              RHO1O(IMOFF-1+(J-1)*NI+I) = RHO1I(IP,JP)
            ELSE
              RHO1I(IP,JP) = RHO1O(IMOFF-1+(J-1)*NI+I) 
            END IF
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
        CALL WRTMAT(RHO1I,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*)
        WRITE(6,*) ' Extracted blocks : '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        CALL APRBLM2(RHO1O,NACOBS,NACOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE EXPRHO1(RHO1I,RHO1O,IRHO1SM)
*
* The density matric rho1 is given in symmetry blocked form 
* (symmetry IRHO1SM); expand to full form
*
* Andreas, to reverse the action of Jeppe's routine above
* currently unused (and maybe not useful anyway)
* ----- do not remember why this was needed once ......
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*. Input
      DIMENSION RHO1O(NTOOB,NTOOB)
*. Output
      DIMENSION RHO1I(*)
*
      RHO1O(1:NTOOB,1:NTOOB) = 0D0
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
            RHO1O(IP,JP) = RHO1I(IMOFF-1+(J-1)*NI+I)
          END DO
        END DO
        IMOFF = IMOFF + NI*NJ
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' EXPRHO1 in action '
        WRITE(6,*) ' Symmetry of blocks extracted ',IRHO1SM
        WRITE(6,*) ' Input density '
        CALL PRHONE(RHO1O,NOCOBS,IRHO1SM,NSMOB,0)
        WRITE(6,*) ' Output density '
        CALL WRTMAT(RHO1I,NTOOB,NTOOB,NTOOB,NTOOB)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
      END IF
*
      RETURN
      END
      SUBROUTINE ONE_EL_PROP(I_EXP_OR_TRA,IIRELAX,RELDEN)
*
* Calculate one-electron properties
* One-electron density is assumed in WORK(KRHO1)
*
* Jeppe Olsen, June 1997
*              Updated Feb. 98 ( Natural orbital analysis added )
*              March 1999 : I_EXP_OR_TRA added => distinquish between 
*                           expectation value and transition value
*              April 1999 : Relaxation terms added
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'crun.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
*
      REAL*8 INPROD 
      CHARACTER*1 XYZ(3)
      CHARACTER*8 LABEL
      CHARACTER*6 LABEL2
      DATA XYZ/'X','Y','Z'/
*. Relaxation contribution to density in complete symmetrypacked form
      DIMENSION RELDEN(*)
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ONE_EL')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Welcome to ONE_EL_PROP'
        WRITE(6,*) ' ======================'
        WRITE(6,*)
      END IF
*. A bit of local memory
      LHONE = NTOOB * NTOOB
      IF(NTEST.GE.100)
     & WRITE(6,*) ' Max size of one-electron operator',LHONE
      CALL MEMMAN(KLHONE,LHONE,'ADDL  ',2,'KLHONE') !done
      CALL MEMMAN(KLCMO,LHONE,'ADDL  ',2,'KLCMO ') !nu
      CALL MEMMAN(KLRHO1S,LHONE,'ADDL  ',2,'LRHO1S') !done
      CALL MEMMAN(KLSCR,4*LHONE,'ADDL  ',2,'KLCMO ') !done
*. Natural orbital expansion
      CALL MEMMAN(KLXNAT,LHONE,'ADDL  ',2,'XNAT  ') !done
*. Integrals in natural orbital basis
      CALL MEMMAN(KLHNAT,LHONE,'ADDL  ',2,'HNAT  ') !done
*. Diagonal of integrals in natorb basis
      CALL MEMMAN(KLHDIA,NTOOB,'ADDL  ',2,'HDIA  ') !done
*. Occupation numbers
      CALL MEMMAN(KLOCC,NTOOB,'ADDL  ',2,'OCCNUM') !done
*
*
*. Assumed symmetry for one-electron density- 
*. for making  change to general transitions densities later
*
      IRHO1SM = 1
*. Extract symmetry blocks from complete one-electron density 
      CALL REORHO1(dbl_mb(KRHO1),dbl_mb(KLRHO1S),IRHO1SM,1)
*. Number of elements in symmetry blocks of integrals and density
      LRHO1S = 0
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IRHO1SM)
        LRHO1S = LRHO1S + NTOOBS(ISM)*NTOOBS(JSM)
      END DO
      IF(NTEST.GE.100)
     & WRITE(6,*) ' Number of elements in symmetry blocks ',
     &LRHO1S
*. Natural orbitals
      IF(I_EXP_OR_TRA.EQ.1.AND.IPRPRO.GE.5) THEN
        CALL NATORB(dbl_mb(KRHO1),NSMOB,NTOOBS,NACOBS,NINOBS,
     &              IREOST,dbl_mb(KLXNAT),
     &              dbl_mb(KLHNAT),dbl_mb(KLOCC),NACOB,
     &              dbl_mb(KLSCR),IPRDEN)
      END IF
C     WRITE(6,*) ' memchk 1 '
C     CALL MEMCHK
*
      DO IPROP =1, NPROP
      IF(NTEST.GE.100)
     &   WRITE(6,'(A,A)') ' Property to be calculated',
     &                     PROPER(IPROP)
*. one- or two-dimensional tensor ?
        LABEL2 = PROPER(IPROP)
        IF(NTEST.GE.100)
     & WRITE(6,'(A,A)') ' LABEL2 =', LABEL2
        CALL GET_PROP_RANK(LABEL2,NRANK)
        IF(NRANK.EQ.-1) THEN
         WRITE(6,'(A,A)') ' Unknown operator ',PROPER(IPROP)
         NRANK = 0
        END IF
        IF(NTEST.GE.100)
     & WRITE(6,*) ' Rank of operator ', NRANK
*.  Permutational symmetry of operator
        CALL GET_PROP_PERMSM(LABEL2,IPERMSM)
C       GET_PROP_PERMSM(PROPER,IPERMSM)
   
        IF(NRANK.EQ.1) THEN
          DO ICOMP = 1, 3
            IF(IXYZSYM(ICOMP).EQ.IRHO1SM) THEN
              IF(NTEST.GE.100)
     & WRITE(6,*) ' right symmetry for component',ICOMP
*. Label of integrals on file- currently DALTON FORM !!
              IF(PROPER(IPROP).EQ.'DIPOLE') THEN
                LABEL =XYZ(ICOMP)//'DIPLEN '
              END IF
              IF(NTEST.GE.100)
     & WRITE(6,'(A,A)') ' Label ', LABEL
*. Obtain one-electron integrals
              CALL GET_PROPINT(dbl_mb(KLHONE),IRHO1SM,LABEL,
     &             dbl_mb(KLSCR),NTOOBS,NTOOBS,NSMOB,0,IPERMSM)
*. Testing Hartree-Fock response
              IIITEST = 0
              IF(IIITEST.EQ.1) THEN 
                WRITE(6,*) ' Hartree-Fock Linear Response '
                WRITE(6,*) ' Hartree-Fock Linear Response '
                WRITE(6,*) ' Hartree-Fock Linear Response '
                WRITE(6,*) ' Hartree-Fock Linear Response '
                WRITE(6,*) ' Hartree-Fock Linear Response '
                WRITE(6,*) ' Hartree-Fock Linear Response '
                CALL LIN_RESP(dbl_mb(KLHONE),1,dbl_mb(KLHONE),1)
              END IF
*

*. and then : Expectation value
              EXPVAL = INPROD(dbl_mb(KLRHO1S),dbl_mb(KLHONE),LRHO1S)
              WRITE(6,'(A,A,E22.14)')
     &        ' Expectation value of ',LABEL,  EXPVAL
              IF(IIRELAX.EQ.1) THEN
                RELAVAL = INPROD(RELDEN,dbl_mb(KLHONE),LRHO1S)
                WRITE(6,'(A,A,E22.14)')
     &          ' Expectation + relaxation term of ',LABEL,  
     &          EXPVAL+RELAVAL 
              END IF
*
              IF(I_EXP_OR_TRA.EQ.1.AND.IPRPRO.GE.5) THEN
*. Analysis in terms of natural orbitals
*. Transform  integrals to nat orb basis
               CALL TRAN_SYM_BLOC_MAT2(dbl_mb(KLHONE),dbl_mb(KLXNAT),
     &              NSMOB,NTOOBS,dbl_mb(KLHNAT),dbl_mb(KLSCR),0)
*. Extract diagonal integrals
               CALL GET_DIAG_BLOC_MAT(dbl_mb(KLHNAT),dbl_mb(KLHDIA),
     &              NSMOB,NTOOBS,0)
               CALL PROP_NATORB(dbl_mb(KLHDIA),dbl_mb(KLOCC),NTOOBS,
     &              NSMOB)
              END IF
*
            END IF
          END DO
        ELSE IF(NRANK.EQ.2) THEN
          DO ICOMP = 1,3
            DO JCOMP = 1,ICOMP
              IF(MULTD2H(IXYZSYM(ICOMP),IXYZSYM(JCOMP))
     &        .EQ.IRHO1SM) THEN
                IF(NTEST.GE.100)
     & WRITE(6,*) ' Right symmetry for components',
     &          ICOMP,JCOMP
*
*. Again DALTON format, we assume a efg component for a specific
* nuclei so the label has the form EFGabc 
                LABEL=XYZ(JCOMP)//XYZ(ICOMP)//LABEL2
*. Obtain one-electron integrals
C               WRITE(6,*) ' memchk 2'
C               CALL MEMCHK
                CALL GET_PROPINT(dbl_mb(KLHONE),IRHO1SM,LABEL,
     &               dbl_mb(KLSCR),NTOOBS,NTOOBS,NSMOB,0,IPERMSM)
*. and then : Expectation value
                EXPVAL = INPROD(dbl_mb(KLRHO1S),dbl_mb(KLHONE),LRHO1S)
                WRITE(6,'(A,A,E22.14)')
     &          ' Expectation value of ',LABEL,  EXPVAL
*
                IF(I_EXP_OR_TRA.EQ.1.AND.IPRPRO.GE.5) THEN
*. Analysis in terms of natural orbitals
*. Transform  integrals to nat orb basis
                 CALL TRAN_SYM_BLOC_MAT2(dbl_mb(KLHONE),dbl_mb(KLXNAT),
     &                NSMOB,NTOOBS,dbl_mb(KLHNAT),dbl_mb(KLSCR),0)
C                WRITE(6,*) ' memchk 3'
C                CALL MEMCHK
*. Extract diagonal integrals
                 CALL GET_DIAG_BLOC_MAT(dbl_mb(KLHNAT),dbl_mb(KLHDIA),
     &                NSMOB,NTOOBS,0)
                 CALL PROP_NATORB(dbl_mb(KLHDIA),dbl_mb(KLOCC),
     &                NTOOBS,NSMOB)
                END IF
*
C               WRITE(6,*) ' memchk 4'
C               CALL MEMCHK
              END IF
            END DO
          END DO
        END IF
*
      END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ONE_EL')
      RETURN
      END
*
      SUBROUTINE GET_PROPINT(H,IHSM,LABEL,SCR,NMO,NBAS,NSM,ILOW,
     &                       IPERMSM)
*
*. Obtain Property integrals in MO basis for operator with
*  label LABEL.
*
* If ILOW = 1, only the elements below the diagonal are 
* obtained.
*
* Jeppe Olsen, June 1997
*              September 97 : ILOW added
*              Feb. 2000    : IPERMSM added
*              May 2012     : MOAO matrix changed to MOAO_ACT
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NMO(*),NBAS(*)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'multd2h.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
*. Output
      DIMENSION H(*)
*. Scratch
      DIMENSION SCR(*), NMOSD2(8)
*. Scratch should atleast have length  3 *  NBAST**2
*
      NTEST = 0
      IF(NTEST.GE.10)
     &WRITE(6,*) ' IPERMSM in PROPINT ', IPERMSM
*. Integrals in AO basis, neglect symmetry 
      NBAST = 0
      NMOT = 0
      DO ISM = 1, NSM
        NBAST = NBAST + NBAS(ISM)
        NMOT  = NMOT  + NMO(ISM)
      END DO
C?    WRITE(6,*) ' Total number of basis functions ',NBAST
      LINTMX = NBAST*NBAST
*
      KLH1AO = 1
      KLFREE = KLH1AO + LINTMX
*
      KLC = KLFREE
      KLFREE = KLC + LINTMX
*  
      KLSCR = KLFREE
      KLFREE = KLFREE + LINTMX
      IF(ENVIRO(1:6).EQ.'DALTON') THEN
*. Obtain AO property integrals 
        CALL GET_H1AO(LABEL,SCR(KLH1AO),IHSM,NBAS,IPERMSM)
C*. Obtain MO-AO transformation matrix
C       CALL GET_CMOAO_ENV(SCR(KLC))                
*. Transform from AO to MO basis
C?      WRITE(6,*) ' MOAO_ACT: '
C?      CALL APRBLM2(WORK(KMOAO_ACT),NMO,NMO,NSM,0)
        CALL TRAH1(NBAS,NMO,NSM,SCR(KLH1AO),WORK(KMOAO_ACT),H,IHSM,
     &             SCR(KLSCR))
      ELSE IF (ENVIRO(1:5).EQ.'LUCIA') THEN
C       WRITE(6,*) ' property integrals will be fetched from FORT91'
C            GET_H1AO(LABEL,IHSM,NBAS,IPERMSM)
*. Obtain AO property integrals 
        CALL GET_H1AO(LABEL,SCR(KLH1AO),IHSM,NBAS,IPERMSM)
*. Transform from AO to MO basis
C?      WRITE(6,*) ' MOAO_ACT: '
C?      CALL APRBLM2(WORK(KMOAO_ACT),NMO,NMO,NSM,0)
        CALL TRAH1(NBAS,NMO,NSM,SCR(KLH1AO),WORK(KMOAO_ACT),H,IHSM,
     &             SCR(KLSCR))
      ELSE IF(ENVIRO(1:4).EQ.'QDOT') THEN
C?     WRITE(6,*) ' Calling _QDOT'
       CALL GET_H1AO_QDOT(LABEL,SCR(KLH1AO),IHSM,NBAS,IPERMSM)
       CALL GET_CMOAO_ENV(SCR(KLC))
*. Transform hole integrals
*. Number of mo-s per hole or particle symmetry
       DO KSM = 1, NSM
         NMOSD2(KSM) = NBAS(KSM)/2
       END DO
*. Offset for particle mo's in MOAO transformation
       IB_P = 1
       DO KSM = 1, NSM
         IB_P = IB_P + NMOSD2(KSM)**2
       END DO
*. Transform holes
       CALL TRAH1(NMOSD2,NMOSD2,NSM,SCR(KLH1AO),SCR(KLC),H,IHSM,
     &             SCR(KLSCR))
C      NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
       LH1INT = NDIM_1EL_MAT(IHSM,NMOSD2,NMOSD2,NSM,0)
*. Transform particles
C           TRAH1(NBAS,NORB,NSM,HAO,C,HMO,IHSM,SCR)
       CALL TRAH1(NMOSD2,NMOSD2,NSM,SCR(KLH1AO+LH1INT),SCR(KLC+IB_P-1),
     &            H(1+LH1INT),IHSM,SCR(KLSCR))
C      MERGE_PH_MAT(PMAT,HMAT,PHMAT,NSM,NTOOBS,IHSM)
*. Merge hole and particle integrals
       CALL MERGE_PH_MAT(H(1+LH1INT),H(1),SCR(KLSCR),NSM,NBAS,IHSM)
       LHP1INT = NDIM_1EL_MAT(IHSM,NBAS,NBAS,NSM,0)
       CALL COPVEC(SCR(KLSCR),H,LHP1INT)
      END IF
*
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) 'electron integrals in MO basis, full format '
        CALL PRHONE(H,NMO,IHSM,NSM,0)
      END IF
      IF(ILOW.EQ.1) THEN
*. Complete to lower half form
        IOFF_IN = 1
        IOFF_OUT = 1
        DO ISM = 1, NSM
          JSM = MULTD2H(ISM,IHSM)
          IF(ISM.EQ.JSM) THEN
*. Copy lower half
            LDIM = NMO(ISM)
            NELMNT_IN = LDIM * LDIM        
            NELMNT_OUT = LDIM * (LDIM + 1)/2         
            CALL COPVEC(H(IOFF_IN),SCR(KLSCR),NELMNT_IN)
            CALL TRIPAK(SCR(KLSCR),H(IOFF_OUT),1,LDIM,LDIM)
            IOFF_IN = IOFF_IN + NELMNT_IN
            IOFF_OUT = IOFF_OUT + NELMNT_OUT
          ELSE IF(ISM.LT.JSM) THEN
*. Just skip block in input matrix
            LIDIM = NMO(ISM)
            LJDIM = NMO(JSM)
            IOFF_IN = IOFF_IN + LIDIM*LJDIM
          ELSE IF(ISM.GT.JSM) THEN
*. Copy block to block
            LIDIM = NMO(ISM)
            LJDIM = NMO(JSM)
            NELMNT = LIDIM*LJDIM
C           CALL TRPMT3(H(IOFF_IN),LIDIM,LJDIM,H(IOFF_OUT))
            CALL COPVEC(H(IOFF_IN),H(IOFF_OUT),NELMNT)
            IOFF_IN = IOFF_IN + NELMNT
            IOFF_OUT = IOFF_OUT + NELMNT
          END IF
        END DO
      END IF
*. The one-electron integrals reside in a NMOT X NMOT matrix.
*. Zero trivial integrals
      IF(ILOW.EQ.1) THEN
        NELMNT = IOFF_OUT-1
      ELSE
        NELMNT = 0
        DO ISM = 1, NSM
          JSM = MULTD2H(ISM,IHSM)
          NELMNT = NELMNT + NMO(ISM)*NMO(JSM)
        END DO
        IFREE = NELMNT + 1
      END IF
C?    WRITE(6,*) ' GET_PROP : NELMNT= ', NELMNT
      ZERO = 0.0D0
      NZERO = NMOT*NMOT - NELMNT
      IFREE = NELMNT + 1
C     WRITE(6,*) 'IFREE NZERO ', IFREE,NZERO
      CALL SETVEC(H(IFREE),ZERO,NZERO)
          
      IF(NTEST .GE. 50 ) THEN
        WRITE(6,*) 'Property integrals in MO basis '
        CALL PRHONE(H,NMO,IHSM,NSM,ILOW)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_CMOAO_DALTON(CMO,NBAS,NMO,NSM)
*
* Obtain MO-AO expansion matrix from SIRIUS/DALTON file SIRGEOM
* 
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER NBAS(*), NMO(*)
*. Output
      DIMENSION CMO(*)
*
      
      ITAP30 = 16  
      OPEN(ITAP30,STATUS='OLD',FORM='UNFORMATTED',FILE='SIRIFC')
      REWIND ITAP30
      CALL MOLLAB('TRCCINT ',ITAP30,6)
*. Skip record containing dimensions of orbitals
      READ(ITAP30)
*. And skip record containing eigenvalues etc
      READ(ITAP30)
C     READ (ITAP30) NSYMHF,NORBT,NBAST,NCMOT,(NOCC(I),I=1,NSYMHF),
C    *              (NLAMDA(I),I=1,NSYMHF),(NORB(I),I=1,NSYMHF),
C    *              POTNUC,EMCSCF
C
C
C     READ (ITAP30) (WRK(KEIGVL+I-1),I=1,NORBT),
C    *              (IWRK(KEIGSY+I-1),I=1,NORBT)
*. And then the MO-AO expansion matrix
      NCOEF = 0
      DO ISM = 1, NSM
        NCOEF = NCOEF + NMO(ISM)*NBAS(ISM)
      END DO 
      READ (ITAP30) (CMO(I),I=1,NCOEF)
      CLOSE(ITAP30,STATUS='KEEP')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  MO - AO expansion matrix ' 
        WRITE(6,*) '============================='
        WRITE(6,*)
        CALL APRBLM2(CMO,NBAS,NMO,NSM,0)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS,NSM,IPERMSM)
*
*. Obtain one-electron integrals in ao basis from dalton
*
* Label of integrals LABEL from FILE AORPROPER
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      CHARACTER*8 LABEL
      DIMENSION NBAS(*)  
      INCLUDE 'multd2h.inc'
*. output
      DIMENSION H1AO(*)
*. Scratch
      DIMENSION SCR(*)
*
      LOGICAL FNDLAB
*
      NTEST = 0
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' Fetching one-electron integrals with Label ',
     &  LABEL
        WRITE(6,*) ' IHSM NSM', IHSM,NSM
      END IF
*
*. Number of elements in : Complete lower half array
*                          Symmetry restricted complete matrix
*                          Symmetry restricted lower half matrix
*-- I am not completely sure about the input format of the integrals
      NBAST = 0
      DO ISM = 1, NSM
       NBAST = NBAST + NBAS(ISM)
      END DO
      NINT01 = NBAST*(NBAST+1)/2
C     write(6,*) ' IHSM = ', IHSM
*
      NINT10 = 0
      DO IRSM = 1, NSM
       ICSM = MULTD2H(IHSM,IRSM)
       NINT10 = NINT10 + NBAS(IRSM)*NBAS(ICSM)
      END DO
*
      NINT11 = 0
      DO IRSM = 1, NSM
       ICSM = MULTD2H(IHSM,IRSM)
       IF(IRSM.GT.ICSM) THEN
        NINT11 = NINT11 + NBAS(IRSM)*NBAS(ICSM)
       ELSE IF(IRSM.EQ.ICSM) THEN
        NINT11 = NINT11 + NBAS(IRSM)*(NBAS(IRSM)+1)/2
       END IF
      END DO
*
*. Read in integrals, assumed in complete lower half format
*
         LUPRP = 15
         OPEN (LUPRP,STATUS='OLD',FORM='UNFORMATTED',FILE='AOPROPER')
         REWIND (LUPRP)
         IF (FNDLAB(LABEL,LUPRP)) THEN
           IF(NTEST.GE.2) write(6,*) ' Label obtained'
           READ(LUPRP) (SCR(I),I=1,NINT01)
           IF(NTEST.GE.2) write(6,*) 'integrals readin'
           IF(NTEST.GE.100) call prsym(scr,NBAST)
C           CALL READT(LUPRP,NBAST*(NBAST+1)/2,WRK(KSCR2))
         ELSE
            WRITE(6,*) 'Property label: ',LABEL ,'  not found on file'
            STOP 'Wrong input or integrals not generated'
         ENDIF
        CLOSE(LUPRP,STATUS='KEEP')
*
        IF(NTEST.GE.2)
     &  WRITE(6,*) ' Number of symmetry apdapted integrals',NINT10
*
*. Transfer integrals to symmetry adapted form, complete form
*
         IF(IPERMSM.EQ.1) THEN
          PSIGN_ = 1.0D0
         ELSE
          PSIGN_ = -1.0D0
         END IF
*
         IBINT = 1
*. Loop over symmetry blocks
         DO IRSM = 1, NSM 
           ICSM = MULTD2H(IHSM,IRSM)
           NR = NBAS(IRSM)
           NC = NBAS(ICSM)
*. Offsets
           IBR = 1
           DO ISM = 1, IRSM - 1
             IBR = IBR + NBAS(ISM)
           END DO
           IBC = 1
           DO ISM = 1, ICSM - 1
             IBC = IBC + NBAS(ISM)
           END DO
*. Complete block, stored in usual column wise fashion
           DO ICORB = 1, NC
             DO IRORB = 1, NR
               ICABS = IBC + ICORB -1
               IRABS = IBR + IRORB -1
               ICRMX = MAX(ICABS,IRABS)
               ICRMN = MIN(ICABS,IRABS)
*. I ASSUME that the phase of integrals corresponds to lower half
*  of integrals
               IF(IRABS.GE.ICABS) THEN
                 H1AO(IBINT-1 + (ICORB-1)*NR+IRORB) =
     &           SCR(ICRMX*(ICRMX-1)/2+ICRMN)
               ELSE
                 H1AO(IBINT-1 + (ICORB-1)*NR+IRORB) =
     &           PSIGN_*SCR(ICRMX*(ICRMX-1)/2+ICRMN)
               END IF
*
             END DO
           END DO
           IBINT = IBINT + NR*NC
         END DO
*

*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' One-electron integrals obtained from AOPROPER'
        CALL PRSYM(SCR,NBAST)
*
        WRITE(6,*) ' One-electron integrals in symmetry-packed form'
        CALL PRHONE(H1AO,NBAS,IHSM,NSM,0)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
      END IF
*
      RETURN
      END 
      SUBROUTINE TRAH1(NBAS,NORB,NSM,HAO,C,HMO,IHSM,SCR)
*
*. Transform one-electron integrals from ao's to mo's.
*
*. Symmetry of integrals is IHSM, all integrals blocks assumed complete,
* i.e not packed to lower half 
*
* Jeppe Olsen
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION HAO(*),C(*)
      DIMENSION NORB(*),NBAS(*)
      INCLUDE 'multd2h.inc'

*. Output
      DIMENSION  HMO(*)
*. Scratch
      DIMENSION SCR(*)
*. Loop over integral blocks
      IBHAO = 1
      IBHMO = 1
      DO IRSM = 1, NSM
        ICSM = MULTD2H(IRSM,IHSM)
*. Pointers to offsets in transformation matrices
        IBR = 1
        DO ISM = 1, IRSM-1
          IBR = IBR + NORB(ISM)*NBAS(ISM)
        END DO
        IBC = 1
        DO ISM = 1, ICSM-1
          IBC = IBC + NORB(ISM)*NBAS(ISM)
        END DO
*. 
        LRMO = NORB(IRSM)
        LRAO = NBAS(IRSM)
*
        LCMO = NORB(ICSM)
        LCAO = NBAS(ICSM)
C       write(6,*) ' TRAH1 : IRSM ICSM ',IRSM,ICSM
C       WRITE(6,*) ' LRAO LRMO LCAO LCMO ',LRAO,LRMO,LCAO,LCMO
 
*
C            MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &             NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
        ZERO = 0.0D0
        ONE= 1.0D0
*.C(row)T*Hao
        CALL SETVEC(SCR,ZERO,LRMO*LCAO)
        CALL MATML7(SCR,C(IBR),HAO(IBHAO),
     &       LRMO,LCAO,LRAO,LRMO,LRAO,LCAO,ZERO,ONE,1)
*. (C(row)T*Hao)*C(column)
        CALL SETVEC(HMO(IBHMO),ZERO,LRMO*LCMO)
        CALL MATML7(HMO(IBHMO),SCR,C(IBC),
     &       LRMO,LCMO,LRMO,LCAO,LCAO,LCMO,ZERO,ONE,0)
*
        IBHAO = IBHAO + LRAO*LCAO
        IBHMO = IBHMO + LRMO*LCMO
*.
      END DO
*
      RETURN
      END
      SUBROUTINE PRHONE(H,NFUNC,IHSM,NSM,IPACK)
*
*. Print one-electron integrals with symmetry IHSM
*. If IPACK = 1, then diagonal blocks are assumed packed
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION H(*), NFUNC(*)
      INCLUDE 'multd2h.inc'
*
      IOFF = 1
      DO IRSM = 1, NSM
        ICSM = MULTD2H(IHSM,IRSM)
        NR = NFUNC(IRSM)
        NC = NFUNC(ICSM)
        WRITE(6,'(A,2I3)') 
     &  ' Blocks with row- and column-symmetry',IRSM,ICSM
        WRITE(6,'(A)')
     &  ' =========================================== '
        IF(IRSM.EQ.ICSM.AND.IPACK.EQ.1) THEN
          CALL PRSYM(H(IOFF),NR)
          IOFF = IOFF + NR*(NR+1)/2
        ELSE IF(IRSM.LT.ICSM .AND. IPACK.EQ.1) THEN
*. Upper block, neglected
        ELSE
          CALL WRTMAT(H(IOFF),NR,NC,NR,NC)
          IOFF = IOFF + NR*NC
        END IF
      END DO 
*
      RETURN
      END
      SUBROUTINE RHO1_HH(RHO1,XLR)  
*
* Add terms form hole-hole commutator to one-particle density matrix :
*
*  2*<L!R> to diagonal for Hole orbitals 
*
* Jeppe Olsen, Jan. 1998 (<= Just to show we are geared for the millenium
*                            change)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*Input and output
      DIMENSION RHO1(NACOB,NACOB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' One-body density matrix before hole-commutator'
        WRITE(6,*) ' (See manual)'
        CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
      END IF
*
      DO IORB = 1, NACOB 
       IF(IPHGAS(ITPFTO(IORB+NINOB)).EQ.2) 
     & RHO1(IORB,IORB) = RHO1(IORB,IORB) + 2.0D0*XLR
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' One-body density matrix with hole-commutator term '
        WRITE(6,*) ' (See manual)'
        CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
      END IF
*
      RETURN
      END
      SUBROUTINE DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR,EXPS2,IDOSRHO12,SRHO1,
     &                  RHO2AA,RHO2AB,RHO2BB,IDONATORB)
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
* two-electron spindensities added, Sept. 2004 ( for singularities in IC...)
*
* Restriction of densities to active spaces allowed, Sept. 2005
* Prepared for explicit inactive/secondary orbital spaces, June 2010
*
*. modified, July 8, 2012 (Jeppe) 
*
* Last modification; Oct. 30, 2012; Jeppe Olsen; call to Z_BLKFO changed
*
* CSF and ICISTR = 1, reintroduced in 2011
* ========================================
* IF ICISTR.gt.1, then L and R are two blocks holding a batch
* IF ICISTR .eq. 1, then are the  two vectors holding a vector over
* parameters. Parameters are CSF's if required
*
* IF CSF's are active (NOCSF = 0), then two vectors over SD's
* must be available (KCOMVECX_SD, X = 1, 2)

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
*
* If IDOSRHO12 = 1, spin density is also calculated
*
* if IDOSRHO12 = 2, then the spin-components of the 
* two-body density are also calculated
*
* RHO2AA(i,j,k,l) = <0!a+_ialpha a+_jalpha a_kalpha a_lalpha!0>
*                   i.ge.j, k.ge.l
* RHO2AA(i,j,k,l) = <0!a+_ibeta a+_jbeta a_kbeta a_lbeta!0>
*                   i.ge.j, k.ge.l
* RHO2AB(i,j,k,l) = <0!a+_ialpha a+_jbeta a_kbeta a_lalpha!0>
*                   No restrictions on i,j,k,l
*
* Call tree for densities : 
* =========================
*  DENSI2 - GASDN2 - GSDNBB2 - GSBBD1
*                            - GSBBD2A
*                            - GSBBD2B

      INCLUDE 'wrkspc.inc' 
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      REAL*8 INPRDD
*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
      INCLUDE 'cands.inc'
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
      INCLUDE 'spinfo.inc'
*
      INCLUDE 'csmprd.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'clunit.inc'
*. Scratch for string information
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. Local scratch
      DIMENSION IDACTSPC(MXPNGAS)
* IDTFREORD Should give reordering full TS order to Density TS order 
      
*. Specific input 
CNW   REAL*8 L, INPROD
      REAL*8 INPROD
CNW   DIMENSION L(*),R(*)
      integer L, R
*.Output
      DIMENSION RHO1(*),RHO2(*),SRHO1(*)
*
      DIMENSION RHO2AA(*),RHO2AB(*),RHO2BB(*)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDEN)
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' Information from DENSI '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' NOCSF = ', NOCSF
        WRITE(6,*) ' I12 = ', I12
        WRITE(6,*)
        WRITE(6,*) ' ICSPC, ISSPC = ', ICSPC, ISSPC
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input L and R vectors '
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(L,1,NSVAR,1,NSVAR)
          WRITE(6,*)
          CALL WRTMAT(R,1,NCVAR,1,NCVAR)
        ELSE 
          CALL WRTVCD(L,LUL,1,-1)
          WRITE(6,*)
          CALL WRTVCD(L,LUR,1,-1)
        END IF
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
      CALL LUCIAQENTER('DENSI')
*
*. Some preparations for batching
*
      IF(ISIMSYM.EQ.0) THEN
        LBLOCK = MXSOOB
      ELSE
        LBLOCK = MXSOOB_AS
      END IF
*. Why the below, this is size of 'inner batch'
      IF(NOCSF.EQ.0) THEN
        LBLOCK  = MAX(NSD_FOR_OCCLS_MAX,MXSOOB)
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' TEST, DENSI2: LCSBLK, LBLOCK, MXSOOB  = ',
     &                             LCSBLK, LBLOCK, MXSOOB
      END IF
*
      ISFIRST = 1
      ICFIRST = 1
      ICOMP = 0
      ILTEST = -3006
      IF(ISFIRST.EQ.1) THEN
        CALL Z_BLKFO_FOR_CISPACE(ISSPC,ISSM,LBLOCK,ICOMP,
     &       NTEST,NSBLOCK,NSBATCH,
     &       int_mb(KSIOIO),int_mb(KSBLTP),NSOCCLS_ACT,
     $       int_mb(KSIOCCLS_ACT),
     &       int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSLBLK),int_mb(KSI1BT),
     &       int_mb(KSIBT),
     &       int_mb(KSNOCCLS_BAT),int_mb(KSIBOCCLS_BAT),ILTEST)
      END IF
      IF(ICFIRST.EQ.1) THEN
        CALL Z_BLKFO_FOR_CISPACE(ICSPC,ICSM,LBLOCK,ICOMP,
     &       NTEST,NCBLOCK,NCBATCH,
     &       int_mb(KCIOIO),int_mb(KCBLTP),NCOCCLS_ACT,
     &       int_mb(KCIOCCLS_ACT),
     &       int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &       int_mb(KCIBT),
     &       int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
      END IF
*
*
      IF(NOCSF.EQ.0.AND.ICNFBAT.GE.2) THEN
*. Obtain scratch files for saving combination forms of C and Sigma
C             FILEMAN_MINI(IFILE,ITASK)
         CALL FILEMAN_MINI(LU_LDET,'ASSIGN')
         CALL FILEMAN_MINI(LU_RDET,'ASSIGN')
C?       WRITE(6,*) ' Test: LU_LDET, LU_RDET: ',
C?   &                      LU_LDET, LU_RDET
      END IF
*
      IF(NOCSF.EQ.0) THEN
       IF(ICNFBAT.EQ.1) THEN
*. Incore:
*. Obtain L and R in SD/Combination basis in vector KCOMVEC1_SD, KCOMVEC2_SD
*. A scratch file is used to allow calc to be done with two complete vectors
*L(SD)
         CALL CSDTVCM(L,KCOMVEC1_SD,KCOMVEC2_SD,
     &        1,0,ISSM,ISSPC,2)
         CALL FILEMAN_MINI(LUSCX,'ASSIGN')
*. Save result in LUSCX
cnw      CALL REWINO(LUSCX)
cnw      CALL TODSC(WORK(KCOMVEC1_SD),NSVAR,NSVAR,LUSCX)
         if(dra_write_section(.false.,KCOMVEC1_SD,1,NSVAR,1,1,LUSCX,
     &      1,NSVAR,1,1,itask).ne.0) call errquit('dra err',LUSCX,911)
         if (dra_wait(itask).ne.0) call errquit('dra read err',
     &       itask,911)
*R(SD)
         CALL CSDTVCM(R,KCOMVEC2_SD,KCOMVEC1_SD,
     &        1,0,ICSM,ICSPC,2)
*Retrieve L(SD) from scratch 
cnw      CALL REWINO(LUSCX)
cnw      CALL FRMDSCO(WORK(KCOMVEC1_SD),NSVAR,NSVAR,LUSCX,IMZERO)
         if (dra_read_section(.false.,KCOMVEC2_SD,1,NSVAR,1,1,LUSCX,
     &       1,NSVAR,1,1,itask).ne.0) call errquit('dra err',LUSCX,911)
         if (dra_wait(itask).ne.0) call errquit('dra read err',
     &       itask,911)
         CALL FILEMAN_MINI(LUSCX,'FREE  ')
       ELSE
*. Not in core
C       CSDTVCMN(CSFVEC,DETVEC,SCR,IWAY,ICOPY,ISYM,ISPC,
C    &           IMAXMIN_OR_GAS,ICNFBAT,LU_DET,LU_CSF,NOCCLS_ACT,
C    &           IOCCLS_ACT,IBLOCK,NBLK_PER_BATCH)
*. 
cnw     CALL REWINO(LUR)
cnw     CALL REWINO(LU_RDET)
*. It is assumed that routines defining expansions, for ex.
*. KCIBT have been constructed
        CALL CSDTVCMN(L,R,dbl_mb(KVEC3),
     &       1,0,ICSM,ICSPC,2,2,LU_RDET,LUR,NCOCCLS_ACT,
     &       int_mb(KCIOCCLS_ACT),int_mb(KCIBT),int_mb(KCLBT))
*. 
        CALL REWINO(LUL)
        CALL REWINO(LU_LDET)
        CALL CSDTVCMN(L,R,dbl_mb(KVEC3),
     &       1,0,ISSM,ISSPC,2,2,LU_LDET,LUL,NSOCCLS_ACT,
     &       int_mb(KSIOCCLS_ACT),int_mb(KSIBT),int_mb(KSLBT))
       END IF ! Incore 
      END IF ! CSFs are in use

*
*. Divide orbital space into inactive(hole), active(valence) and secondary(particle)
*  based on space ISSPC. Result is stored in IHPVGAS
*. IF ISSPC and ICSPC differs in division into HPV, and some orbital spaces are
*. excluded one may run into problems
      CALL CC_AC_SPACES(ISSPC,IREFTYP)
      DO IGAS = 1, NGAS
        ITYP = IHPVGAS(IGAS)
        IDACTSPC(IGAS) = 0
        IF(ITYP.EQ.1.AND.IDENS_IN.EQ.1) IDACTSPC(IGAS) = 1
        IF(ITYP.EQ.2.AND.IDENS_SEC.EQ.1) IDACTSPC(IGAS) = 1
        IF(ITYP.EQ.3.AND.IDENS_AC.EQ.1) IDACTSPC(IGAS) = 1
      END DO
      IF(NTEST.GT.5) THEN
       WRITE(6,*) ' Active orbital spaces in DENSI '
       CALL IWRTMA(IDACTSPC,1,NGAS,1,NGAS)
      END IF
*. Number of orbitals for which densities will be constructed 
      NDACTORB = 0
      DO IGAS = 1, NGAS
        IF(IDACTSPC(IGAS).EQ.1) NDACTORB = NDACTORB + NOBPT(IGAS)
      END DO
      IF(NTEST.GT.5) 
     & WRITE(6,*) ' Number of active orbitals for densities',NDACTORB
*. Offsets to restricted set of orbital spaces
      CALL IB_FOR_SEL_ORBSPC(NOBPTS,NOBPS_SEL,IOBPTS_SEL,IDACTSPC,NGAS,
     &     MXPNGAS,NSMST)
C?    WRITE(6,*) ' NOBPS_SEL after call to IB_FOR... '
C?    CALL IWRTMA(NOBPS_SEL,1,NSMOB,1,NSMOB)
C     IB_FOR_SEL_ORBSPC(NOBPTS,NOBPS_SEL,IOBPTS_SEL,I_SEL,
C    &           NGAS,MXPNGAS,NSYM)
*. And reorder arrays for going betweeen  complete set of 
*. orbitals and restricted set of orbitals
      CALL IREO_DACT_TS(NOBPTS,IOBPTS_SEL,IDACTSPC,IREO_SELTF,
     &     IREO_FTSEL,NGAS,MXPNGAS,NSMOB)
C     IREO_DACT_TS(NOBPTS,IOBPTS_SEL,I_SEL,
C    &           IDTFREO,IFTDREO,NGAS,MXPNGAS,NSMOB)
*
      ZERO = 0.0D0
      CALL SETVEC(RHO1,ZERO ,NDACTORB ** 2 )
      IF(I12.EQ.2) 
     &CALL SETVEC(RHO2,ZERO ,NDACTORB ** 2 *(NDACTORB**2+1)/2)
*
      IF(IDOSRHO12.GE.1) THEN
        IDOSRHO1 = 1
      ELSE 
        IDOSRHO1 = 0
      END IF
*
      IF(IDOSRHO12.EQ.2) THEN
        IDOSRHO2 = 1
      ELSE 
        IDOSRHO2 = 0
      END IF
*.
      IF(IDOSRHO1.EQ.1) THEN  
        CALL SETVEC(SRHO1,ZERO,NDACTORB ** 2)
      END IF
*
      IF(IDOSRHO2.EQ.1) THEN 
        LEN_RHO2AA = (NDACTORB*(NDACTORB+1)/2)**2
        LEN_RHO2AB = NDACTORB ** 4
        CALL SETVEC(RHO2AA,ZERO,LEN_RHO2AA)
        CALL SETVEC(RHO2BB,ZERO,LEN_RHO2AA)
        CALL SETVEC(RHO2AB,ZERO,LEN_RHO2AB)
      END IF
*          
C?     WRITE(6,*) ' ISSPC ICSPC in DENSI2 ',ISSPC,ICSPC
*
* Info for this internal space
*. type of alpha and beta strings - as H does not change 
*. the number of electrons, I do not distinguish between spaces for C
*and S
      IF(ICSPC.LE.NCMBSPC) THEN
       IATP = 1
       IBTP = 2
      ELSE
       IATP = IALTP_FOR_GAS(ICSPC)
       IBTP = IBETP_FOR_GAS(ICSPC)
C?     WRITE(6,*) ' DENSI2 : IATP, IBTP = ', IATP, IBTP
      END IF
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. alpha and beta strings with an electron removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-1,IATPM1)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-1,IBTPM1)
*. alpha and beta strings with two electrons removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-2,IATPM2)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-2,IBTPM2)

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
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS') !done
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD') !done
      CALL STSTSM(dbl_mb(KSTSTS),dbl_mb(KSTSTD),NSMST)
*. connection matrices for supergroups
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA') !done
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB') !done
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPB),IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(int_mb(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR          
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(int_mb(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(int_mb(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MAXA = MAX(MAXA,MAXA0)   
      MAXB = MAX(MAXB,MAXB0)   
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRDEN.GT.5 ) WRITE(6,*)
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
      IF(IPRDEN.GT.5) write(6,*) ' DENSI2 : MXSB MXTSOB MXSOOB ',
     &       MXSB,MXTSOB,MXSOOB 
      IF(ISIMSYM.NE.1) THEN
        LSCR1 = MXSOOB
      ELSE
        LSCR1 = MXSOOB_AS
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRDEN.GT.5)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
      IF(ICISTR.EQ.1) THEN
        CALL MEMMAN(KCB,LSCR1,'ADDL  ',2,'KCB   ') !nu
        CALL MEMMAN(KSB,LSCR1,'ADDL  ',2,'KSB   ') !nu
      END IF
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4
      IF(IPRDEN.GT.5)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ') !done
*
      CALL IAIBCM(ISSPC,int_mb(KSIOIO))
      CALL IAIBCM(ISSPC,int_mb(KCIOIO))
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESCPH(int_mb(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSTFSMSPGP,MXPNSMST,
     &     NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &     NELFSPGP,
     &     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &     IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,
     &     MX_NSPII)
      IF(IPRDEN.GT.5) THEN
        WRITE(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF(IPRDEN.GT.5)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
      LSCR12 = MAX(LSCR1,2*LSCR2)
*. It is assumed that the third block already has been allocated, so
      KC2 = KVEC3
      IF(IPRCIX.GT.5)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL MEMMAN(KI1,  LSCR3       ,'ADDL  ',1,'I1    ') !done
      CALL MEMMAN(KI2,  LSCR3       ,'ADDL  ',1,'I2    ') !done
      CALL MEMMAN(KI3,  LSCR3       ,'ADDL  ',1,'I3    ') !done
      CALL MEMMAN(KI4,  LSCR3       ,'ADDL  ',1,'I4    ') !done
      CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ') !done
      CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ') !done
      CALL MEMMAN(KXI3S,LSCR3       ,'ADDL  ',2,'XI3S  ') !done
      CALL MEMMAN(KXI4S,LSCR3       ,'ADDL  ',2,'XI4S  ') !done
*
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,int_mb(KSBLTP),int_mb(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,int_mb(KCBLTP),int_mb(KSVST))
*.0 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
* scratch space containing active one body
      CALL MEMMAN(KRHO1S,NDACTORB ** 2,'ADDL  ',2,'RHO1S ') !done
*. For natural orbitals
      CALL MEMMAN(KRHO1P,NDACTORB*(NDACTORB+1)/2,'ADDL  ',2,'RHO1P ') !d
      CALL MEMMAN(KXNATO,NDACTORB **2,'ADDL  ',2,'XNATO ') !done
*. Natural orbitals in symmetry blocks
      CALL MEMMAN(KRHO1SM,NDACTORB ** 2,'ADDL  ',2,'RHO1S ') !done
      CALL MEMMAN(KXNATSM,NDACTORB ** 2,'ADDL  ',2,'RHO1S ') !nu
      CALL MEMMAN(KOCCSM,NDACTORB ,'ADDL  ',2,'RHO1S ') !done
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
*
      NTTS = MXNTTS
      CALL MEMMAN(KLSCLFCL,NTTS, 'ADDL  ',2,'SCLF_L') !done
      CALL MEMMAN(KLSCLFCR,NTTS, 'ADDL  ',2,'SCLF_R') !done
      CALL MEMMAN(KLSCLNRM,NTTS, 'ADDL  ',2,'NORM_R') !done

      NVARL = NSVAR
      NVARR = NCVAR
      S2_TERM1 = 0.0D0
*
      IF(ICISTR.EQ.1) THEN
       LLUL = 0
       LLUR = 0
      ELSE
       IF(NOCSF.EQ.1) THEN
        LLUL = LUL
        LLUR = LUR
       ELSE
        LLUL = LU_LDET
        LLUR = LU_RDET
       END IF
      END IF
*
      NBATCHL = NSBATCH
      NBATCHR = NCBATCH
     
      ILTEST = -3006
      IF(ICISTR.GE.2) THEN
*. Out of core version
        CALL GASDN2(I12,RHO1,RHO2,L,R,dbl_mb(KC2),
     &       int_mb(KCIOIO),int_mb(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),int_mb(KCBLTP),int_mb(KSBLTP),NACOB,
     &       int_mb(KNSTSO(IATP)),int_mb(KISTSO(IATP)),
     &       int_mb(KNSTSO(IBTP)),int_mb(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,      
     &       MAXK,MAXI,LSCR1,LSCR1,dbl_mb(KCSCR),dbl_mb(KSSCR),
     &       SXSTSM,dbl_mb(KSTSTS),dbl_mb(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &       int_mb(KI1),dbl_mb(KXI1S),int_mb(KI2),dbl_mb(KXI2S),
     &       int_mb(KI3),dbl_mb(KXI3S),int_mb(KI4),dbl_mb(KXI4S),
     &       dbl_mb(KINSCR),
     &       MXPOBS,IPRDEN,dbl_mb(KRHO1S),LLUL,LLUR,
     &       PSSIGN,PSSIGN,dbl_mb(KRHO1P),dbl_mb(KXNATO),
     &       NBATCHL,int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSI1BT),
     &       int_mb(KSIBT),
     &       NBATCHR,int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCI1BT),
     &       int_mb(KCIBT),int_mb(KCONSPA),int_mb(KCONSPB),
     &       dbl_mb(KLSCLFCL),dbl_mb(KLSCLFCR),dbl_mb(KLSCLNRM),
     &       S2_TERM1,IUSE_PH,IPHGAS,
     &       IDOSRHO1,SRHO1,IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &       NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,IOBPTS_SEL,
     &       NINOB,ICISTR,LV,RV,ICISTR,NVARL,NVARR,ILTEST)
*
C    &       NTEST,NSBLOCK,NSBATCH,
C    &       WORK(KSIOIO),WORK(KSBLTP),NSOCCLS_ACT,WORK(KSIOCCLS_ACT),
C    &       WORK(KSLBT),WORK(KSLEBT),WORK(KSLBLK),WORK(KSI1BT),
C    &       WORK(KSIBT),
C    &       WORK(KSNOCCLS_BAT),WORK(KSIBOCCLS_BAT),0,ILTEST)
      ELSE
*, In core version
       IF(NOCSF.EQ.1) THEN
        CALL GASDN2(I12,RHO1,RHO2,dbl_mb(KVEC1P),dbl_mb(KVEC2P),
     &       dbl_mb(KC2),
     &       int_mb(KCIOIO),int_mb(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),int_mb(KCBLTP),int_mb(KSBLTP),NACOB,
     &       int_mb(KNSTSO(IATP)),int_mb(KISTSO(IATP)),
     &       int_mb(KNSTSO(IBTP)),int_mb(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,      
     &       MAXK,MAXI,LSCR1,LSCR1,dbl_mb(KCSCR),dbl_mb(KSSCR),
     &       SXSTSM,dbl_mb(KSTSTS),dbl_mb(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &       int_mb(KI1),dbl_mb(KXI1S),int_mb(KI2),dbl_mb(KXI2S),
     &       int_mb(KI3),dbl_mb(KXI3S),int_mb(KI4),dbl_mb(KXI4S),
     &       dbl_mb(KINSCR),
     &       MXPOBS,IPRDEN,dbl_mb(KRHO1S),-1,-1,
     &       PSSIGN,PSSIGN,dbl_mb(KRHO1P),dbl_mb(KXNATO),
     &       NBATCHL,int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSI1BT),
     &       int_mb(KSIBT),
     &       NBATCHR,int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCI1BT),
     &       int_mb(KCIBT),int_mb(KCONSPA),int_mb(KCONSPB),
     &       dbl_mb(KLSCLFCL),dbl_mb(KLSCLFCR),dbl_mb(KLSCLNRM),
     &       S2_TERM1,IUSE_PH,IPHGAS,
     &       IDOSRHO1,SRHO1,IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &       NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,IOBPTS_SEL,
     &       NINOB,ICISTR,L,R,ICISTR,NVARL,NVARR,ILTEST)
C    &       NTEST,NSBLOCK,NSBATCH,
C    &       WORK(KSIOIO),WORK(KSBLTP),NSOCCLS_ACT,WORK(KSIOCCLS_ACT),
C    &       WORK(KSLBT),WORK(KSLEBT),WORK(KSLBLK),WORK(KSI1BT),
C    &       WORK(KSIBT),
C    &       WORK(KSNOCCLS_BAT),WORK(KSIBOCCLS_BAT),0,ILTEST)
       ELSE
*. CSF's in use
        CALL GASDN2(I12,RHO1,RHO2,dbl_mb(KVEC1P),dbl_mb(KVEC2P),
     &       dbl_mb(KC2),
     &       int_mb(KCIOIO),int_mb(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),int_mb(KCBLTP),int_mb(KSBLTP),NACOB,
     &       int_mb(KNSTSO(IATP)),int_mb(KISTSO(IATP)),
     &       int_mb(KNSTSO(IBTP)),int_mb(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,      
     &       MAXK,MAXI,LSCR1,LSCR1,dbl_mb(KCSCR),dbl_mb(KSSCR),
     &       SXSTSM,dbl_mb(KSTSTS),dbl_mb(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &       int_mb(KI1),dbl_mb(KXI1S),int_mb(KI2),dbl_mb(KXI2S),
     &       int_mb(KI3),dbl_mb(KXI3S),int_mb(KI4),dbl_mb(KXI4S),
     &       dbl_mb(KINSCR),
     &       MXPOBS,IPRDEN,dbl_mb(KRHO1S),-1,-1,
     &       PSSIGN,PSSIGN,dbl_mb(KRHO1P),dbl_mb(KXNATO),
     &       NBATCHL,int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSI1BT),
     &       int_mb(KSIBT),
     &       NBATCHR,int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCI1BT),
     &       int_mb(KCIBT),int_mb(KCONSPA),int_mb(KCONSPB),
     &       dbl_mb(KLSCLFCL),dbl_mb(KLSCLFCR),dbl_mb(KLSCLNRM),
     &       S2_TERM1,IUSE_PH,IPHGAS,
     &       IDOSRHO1,SRHO1,IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &       NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,IOBPTS_SEL,NINOB,
     &       ICISTR,WORK(KCOMVEC1_SD),WORK(KCOMVEC2_SD),
     &       NVARL,NVARR,ILTEST)
C    &       NTEST,NSBLOCK,NSBATCH,
C    &       WORK(KSIOIO),WORK(KSBLTP),NSOCCLS_ACT,WORK(KSIOCCLS_ACT),
C    &       WORK(KSLBT),WORK(KSLEBT),WORK(KSLBLK),WORK(KSI1BT),
C    &       WORK(KSIBT),
C    &       WORK(KSNOCCLS_BAT),WORK(KSIBOCCLS_BAT),0,ILTEST)
       END IF! CSF Switch for incore version
      END IF! ICISTR switch
      call ga_sync()
      call ga_dgop(1,RHO1,NDACTORB ** 2,'+')
      call ga_dgop(1,S2_TERM1,1,'+')
      if (I12.EQ.2) call ga_dgop(2,RHO2,
     &   (NDACTORB**2)*(NDACTORB**2+1)/2,'+')

      IF(NTEST.GE.100) WRITE(6,*) ' Returned from GASDN2' 
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Rho1 after call to GASDN2'
        CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
      END IF
*
*
*. Add terms from hole-hole commutator
      IF(IUSE_PH.EQ.1) THEN 
*. Overlap between left and right vector
       IF(ICSM.NE.ISSM) THEN
        XLR = 0.0D0
       ELSE
        IF(ICISTR.EQ.1) THEN
         XLR = INPROD(L,R,NSVAR)
        ELSE
         XLR = INPRDD(L,R,LLUR,LLUL,1,-1)
        END IF
       END IF
C?     WRITE(6,*) ' XLR = ', XLR
       CALL RHO1_HH(RHO1,XLR)
      END IF
*
      IF(NOCSF.EQ.0.AND.ICISTR.GE.2) THEN
*. Free the temp files
       CALL FILEMAN_MINI(LU_LDET,'FREE  ')
       CALL FILEMAN_MINI(LU_RDET,'FREE  ')
      END IF
*
* Natural Orbitals 
       I_OLD_OR_NEW = 2
       IF(IDONATORB.EQ.1) THEN
         IF(I_OLD_OR_NEW.EQ.1) THEN
           CALL NATORB3(RHO1,NSMOB,NACOBS,NINOBS,NSCOBS,NINOB,NACOB,
     &          IREOST,dbl_mb(KXNATO),dbl_mb(KRHO1SM),dbl_mb(KOCCSM),
     &          dbl_mb(KRHO1P),IPRDEN)
C               NATORB3(RHO1,NSMOB,NACOBS,NINOBS,NSCOBS,
C    &                  NINOB,NACOB,ISTREO,XNAT,RHO1SM,OCCNUM,
C    &                  SCR,IPRDEN)
         ELSE
*. Obtain natural orbitals in blocks over general symmetry
C         NATORB3_GS(RHO1,XNAT,RHO1SM,OCCNUM,SCR,IREO_GS_TO_TS,IPRDEN)
          CALL NATORB3_GS(RHO1,dbl_mb(KXNATO),dbl_mb(KRHO1SM),
     &         dbl_mb(KOCCSM),
     &         dbl_mb(KRHO1P),int_mb(KIREO_GNSYM_TO_TS_ACOB),IPRDEN)
         END IF !old_or_new switch
       END IF ! natural orbitals requested
*
      IF(IPRDEN.GE.20.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' One-electron density matrix '
        WRITE(6,*) ' ============================'
        CALL WRTMAT(RHO1,NDACTORB,NDACTORB,NDACTORB,NDACTORB) 
      END IF
      IF(IPRDEN.GE.100.AND.I12.EQ.2.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) ' Two-electron density '
          CALL PRSYM(RHO2,NDACTORB**2)
      END IF
*
      IF(I12.EQ.2) THEN
* <L!S**2|R>
        EXPS2 = S2_TERM1 + NAEL +
     &          0.5*(NAEL-NBEL)*(0.5*(NAEL-NBEL)-1)
        IF(IPRDEN.GE.5) THEN
COLD      WRITE(6,*) ' Term 1 to S2 ', S2_TERM1
          WRITE(6,*) ' Expectation value of S^2 ', EXPS2
        END IF
      ELSE
        EXPS2 = 0.0D0
      END IF
*
      IF(IDOSRHO1.EQ.1.AND.IPRDEN.GE.5) THEN
        WRITE(6,*) ' One-electron spindensity <0!E(aa) - E(bb)!0> '
        CALL WRTMAT(SRHO1,NDACTORB,NDACTORB,NDACTORB,NDACTORB)
      END IF
*
      IF(IPRDEN.GE.100.AND.IDOSRHO2.EQ.1) THEN
        WRITE(6,*) ' The RHO2AA(ij,kl) spin density '
        NDIM_AA = NDACTORB*(NDACTORB+1)/2
        CALL WRTMAT(RHO2AA,NDIM_AA,NDIM_AA,NDIM_AA,NDIM_AA)
        WRITE(6,*) ' The RHO2BB(ij,kl) spin density '
        CALL WRTMAT(RHO2BB,NDIM_AA,NDIM_AA,NDIM_AA,NDIM_AA)
        WRITE(6,*) ' The RHO2AB(ik,lj) spin density '
        NDIM_AB = NDACTORB*NDACTORB
        CALL WRTMAT(RHO2AB,NDIM_AB,NDIM_AB,NDIM_AB,NDIM_AB)
      END IF
*
      I_CHECK_SRHO2 = 1
      IF(IDOSRHO2.EQ.1.AND.I_CHECK_SRHO2.EQ.1) THEN
*. Obtain standard rho2 from rho2s and check
C             TEST_RHO2S(RHO2,RHO2AA,RHO2AB,RHO2BB,NORB)
         CALL TEST_RHO2S(RHO2,RHO2AA,RHO2AB,RHO2BB,NTOOB)
      END IF

*
*. Eliminate local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'DENSI ')
      CALL LUCIAQEXIT('DENSI')
C     WRITE(6,*) ' Leaving DENSI '
      RETURN
      END
      SUBROUTINE GASDN2(I12,RHO1,RHO2,
     &           CB,SB,C2,ICOCOC,ISOCOC,ICSMOS,ISSMOS,
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
     &           ICONSPA,ICONSPB,SCLFAC_L,SCLFAC_R,RNORM,S2_TERM1,
     &           IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &           IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &           NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,IOBPTS_SEL,
     &           NINOB,ICISTR,LV,RV,NVARL,NVARR,ILTEST)
*
*
* Jeppe Olsen , Winter of 1991
* GAS modificatios, August 1995
*
* Table driven, June 97
*
* Last revision : Jan. 98 (IUSE_PH,IPHGAS added)
*                 Jan. 99 (IDOSRHO1,SRHO1 added)
*                 Sept.04 (IDOSRHO2,RHO2AA,RHO2AB,RHO2BB added)
*                 June 10 (NINOB added)
*                 Jan  11 (LR, LV added, Incore version alive again)
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
*
* ICISTR = 1: 
* L and R vectors are stored in LV, RV. 
* ICISTR > 2:
* The L and R vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
* 
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
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
*. Info on the orbital spaces that are active in density calculations
      INTEGER IDACTSPC(*),IDTFREORD(*),IFTDREORD(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),XI1S(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION RHO1S(*)
      DIMENSION SCLFAC_L(*),SCLFAC_R(*),RNORM(*)
*.
      INTEGER LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      INTEGER RASM(4),RBSM(4),RATP(4),RBTP(4),RSGN(5),RTRP(5)
      LOGICAL ISLOCAL
      REAL * 8 INPROD
*. Vectors holding L and R if ICISTR = 1
CNW   REAL*8
CNW  & LV(*), RV(*)
      integer LV, RV
*.Output
      DIMENSION RHO1(*),RHO2(*)
      DIMENSION RHO2AA(*),RHO2AB(*),RHO2BB(*),SRHO1(*)
      DIMENSION RHO1P(*),XNATO(*)
*
      CALL LUCIAQENTER('GASDN')
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GASDN2 speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' NACOB,MAXK,NGAS,IDC,MXPOBS',
     &             NACOB,MAXK,NGAS,IDC,MXPOBS
        WRITE(6,*) ' LUL, LUR ', LUL,LUR
        WRITE(6,*) ' PSL, PSR = ', PSL, PSR
        WRITE(6,*) ' NVARL, NVARR = ', NVARL, NVARR
        WRITE(6,*) ' ILTEST = ', ILTEST
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial L vector '
        IF(ICISTR.EQ.1) THEN
          CALL WRTRS2(LV,ISSMOS,ISBLTP,ISOCOC,NOCTPA,NOCTPB,
     &                NSSOA,NSSOB,NSMST)
        ELSE
          CALL WRTVCD(CB,LUL,1,-1)
        END IF
        WRITE(6,*) ' Initial R vector '
        IF(ICISTR.EQ.1) THEN
          CALL WRTRS2(RV,ICSMOS,ICBLTP,ICOCOC,NOCTPA,NOCTPB,
     &                NSSOA,NSSOB,NSMST)
        ELSE
          CALL WRTVCD(CB,LUR,1,-1)
        END IF
      END IF
C
C SET SCLFAC_L and SCLFACE_R VALUES FOR BLOCKS.
C GIVEN WE'RE DOING MCSCF THEY ARE THE SAME.
C NOW WE KNOW IF BLOCKS ARE BIG ENOUGH TO USE
C
      NCBLOCK_T = IELSUM(LBATL(1), NBATCHL)
      CALL FIND_ACTIVE_BLOCKS_IN_VECTOR(LV,SCLFAC_L,NCBLOCK_T,
     &                                  IBLOCKL(1,1))
      NBLOCK_C =  I1BATL(NBATCHL)-1+ LBATL(NBATCHL)
      CALL ENFORCE_SIMSYM_ACT_BLK(NBLOCK_C,IBLOCKL,SCLFAC_L)
      DO IIL = 1, NCBLOCK_T
        RNORM(IIL) = SCLFAC_L(IIL)
      END DO 
C 
C LET'S USE LV (LEFT SIDE AS LOCAL). SKIP WHEN BLOCK NOT LOCAL. 
C SIMPLY SET SCLFAC_L TO ZERO FOR NON-LOCAL LV BLOCKS TO SKIP
C
      call ga_distribution(LV,ga_nodeid(),ilo,ihi,idum,idum)
C
      IF(ICISTR.GT.1) THEN
        CALL REWINO(LUL)
      ELSE
        IOFF_L = 1
      END IF
      DO 10001 IBATCHL = 1, NBATCHL
*. Obtain L blocks
        NBLKL = LBATL(IBATCHL)
        IF(LUL.LE.0) THEN
          IF(IBATCHL.EQ.1) THEN
           IOFF_L_BATCH = 1
          ELSE
           IOFF_L_BATCH = IOFF_L_BATCH+LEBATL(IBATCHL-1)
          END IF
        END IF
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
*. Offset for this block in complete vector
          IF(ICISTR.LE.1) THEN
            IOFF_L_BLOCK = IOFF_L_BATCH-1+IBLOCKL(6,IL)
          ELSE
            IOFF_L_BLOCK = 1
          END IF
C
          PLL = 1.0D0
CNW       CALL GSTTBL(LV(IOFF_L_BLOCK),SB(IOFF),IATP,IASM,IBTP,IBSM,
          if (IOFF_L_BLOCK.GE.ilo.AND.
     &        (IOFF_L_BLOCK+IBLOCKL(8,IL)-1).LE.ihi) then
          CALL GSTTBL(LV,RNORM(IL),IOFF_L_BLOCK,SB(IOFF),IATP,IASM,
     &         IBTP,IBSM,ISOCOC,NOCTPA,NOCTPB,NSSOA,NSSOB,PSL,ISOOSC,
     &         IDC,PLL,LUL,C2,NSMST,ISCALE,SCLFAC_L(IL))
          endif
        END DO
*. Loop over batches  of R vector
        IF(ICISTR.GT.1) THEN
          CALL REWINO(LUR)
        END IF
        IOFF_R_BATCH = 0
        DO IR = 1, NCBLOCK_T
          SCLFAC_R(IR) = RNORM(IR)                  
        END DO 
        DO 9001 IBATCHR = 1, NBATCHR
*. Address of start of batch of R
          IF(LUR.LE.0) THEN
            IF(IBATCHR.EQ.1) THEN
             IOFF_R_BATCH = 1
            ELSE
             IOFF_R_BATCH = IOFF_R_BATCH+LEBATR(IBATCHR-1)
            END IF
          END IF
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
*. Offset for this block in complete vector
          IF(ICISTR.LE.1) THEN
            IOFF_R_BLOCK = IOFF_R_BATCH-1+IBLOCKL(6,IR)
          ELSE
            IOFF_R_BLOCK = 1
          END IF
*
          IF(NTEST.GE.2200)
     &    WRITE(6,'(A,4(2X,I4))') ' JATP JBTP JASM JBSM ',
     &                              JATP,JBTP,JASM,JBSM
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
            ISTRFL = 1
            PL = 1.0D0
            CALL PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PL,
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
C?        WRITE(6,*) ' Jeppe sets interact to 1 '
C?        INTERACT = 1
              IF(IABEXC.LE.I12) THEN
                INTERACT = 1
              END IF
            END DO
          END DO
*.          ^ End of checking whether C-block is needed
          ISCALE = 1
          ISLOCAL = IOFF_L_BLOCK.GE.ilo.AND.
     &              IOFF_L_BLOCK+IBLOCKL(8,IL)-1.LE.ihi
          IF(INTERACT.EQ.1.AND.ISLOCAL) THEN
            ISCALE = 1
C?          WRITE(6,*) ' LUR before GSTTBL', LUR
            PLR = 1.0D0
            IF(NTEST.GE.2100) THEN
              WRITE(6,*) ' PSR, PLR before GSTTBL for R ',
     &                     PSR, PLR
            END IF
CNW         CALL GSTTBL(RV(IOFF_R_BLOCK),CB(JOFF),JATP,JASM,JBTP,JBSM,
            CALL GSTTBL(RV,RNORM(IR),IOFF_R_BLOCK,CB(JOFF),JATP,JASM,
     &           JBTP,JBSM,ICOCOC,NOCTPA,NOCTPB,NSSOA,NSSOB,PSR,ICOOSC,
     &           IDC,PLR,LUR,C2,NSMST,ISCALE,SCLFAC_R(IR))
C?        WRITE(6,*) ' IR, SCLFAC_R(IR) = ', IR, SCLFAC_R(IR)
          ELSE
c             WRITE(6,*) ' TTSS for C block skipped  '
c             CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
            IF(LUR.GT.0) THEN
              CALL IFRMDS(LBL,-1,1,LUR)
              CALL SKPRCD2(LBL,-1,LUR)
            END IF
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
            IF(NTEST.GE.100) WRITE(6,'(A,2I4,E8.3)')
     &      '  IIL, IL, SCLFAC_L(IL) = ', 
     &         IIL, IL, SCLFAC_L(IL)
            IF(ICISTR.LE.1) THEN
              IOFF_L_BLOCK = IOFF_L_BATCH-1+IBLOCKL(6,IL)
            ELSE
              IOFF_L_BLOCK = 1
            END IF
            ISLOCAL = IOFF_L_BLOCK.GE.ilo.AND.
     &              IOFF_L_BLOCK+IBLOCKL(8,IL)-1.LE.ihi
          IF(SCLFAC_L(IL).NE.0.0D0.AND.ISLOCAL) THEN
            IATP = IBLOCKL(1,IL)
            IBTP = IBLOCKL(2,IL)
            IASM = IBLOCKL(3,IL)
            IBSM = IBLOCKL(4,IL)
            IOFF = IBLOCKL(5,IL)
*
            NIA = NSSOA(IASM,IATP)
            NIB = NSSOB(IBSM,IBTP)
*. Possible permutations of L blocks
            PL = 1.0D0
            CALL PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PL,
     &           LATP,LBTP,LASM,LBSM,LSGN,LTRP,NLPERM)
            DO 9999 ILPERM = 1, NLPERM
C             write(6,*) ' Loop 9999 ILPERM = ', ILPERM
              IIASM = LASM(ILPERM)
              IIBSM = LBSM(ILPERM)
              IIATP = LATP(ILPERM)
              IIBTP = LBTP(ILPERM)
              IF(NTEST.GE.100) THEN
                WRITE(6,'(A,4(2X,4I4))')
     &          ' L: IIASM, IIBSM, IIATP, IIBTP = ', 
     &               IIASM, IIBSM, IIATP, IIBTP
              END IF
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
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' IIR, IR, SCLFAC_R(IR) = ',
     &                         IIR, IR, SCLFAC_R(IR)
                END IF
              IF(SCLFAC_R(IR).NE.0.0D0) THEN
                JATP = IBLOCKR(1,IR)
                JBTP = IBLOCKR(2,IR)
                JASM = IBLOCKR(3,IR)
                JBSM = IBLOCKR(4,IR)
                JOFF = IBLOCKR(5,IR)
                IF(NTEST.GE.100) THEN
                  WRITE(6,'(A,4(2X,4I4))')
     &            '  R: JASM, JBSM, JATP, JBTP =     ', 
     &                  JASM, JBSM, JATP, JBTP
                END IF
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
     &                  PSR,PL,RATP,RBTP,RASM,RBSM,RSGN,RTRP,
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
                       WRITE(6,*) ' GSDNBB2 will be called for '
                       WRITE(6,'(A,4I5)') 
     &                 '  L: IIASM IIBSM IIATP IIBTP',
     &                       IIASM,IIBSM,IIATP,IIBTP
                       WRITE(6,'(A,4I5)') 
     &                 '  R: JASM JBSM JATP JBTP    ',
     &                       JASM,JBSM,JATP,JBTP
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
     &                    S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &                    IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &                    NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,
     &                    IOBPTS_SEL,NINOB)
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
      IF(NTEST.GE.100) WRITE(6,*) ' Returning from GASDN2 '
      CALL LUCIAQEXIT('GASDN')
      RETURN
      END
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
     &                  SCLFAC,S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &                  IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &                  NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,
     &                  IOBPTS_SEL,NINOB)
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
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*. Input  
      DIMENSION CB(*),SB(*)
      INTEGER IDACTSPC(*),IDTFREORD(*),IFTDREORD(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
      INTEGER IAOC(*), IBOC(*), JAOC(*), JBOC(*)
*. Output
      DIMENSION RHO1(*),RHO2(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION  I1(*),XI1S(*),I2(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION C2(*)
*
      CALL LUCIAQENTER('GSDN2')
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      NTESTO= NTEST
      IF(NTEST.GE.200) THEN
c       WRITE(6,*) ' ================'
c       WRITE(6,*) ' GSDNBB2: R block '
c       WRITE(6,*) ' ================='
c       CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
c       WRITE(6,*) ' ================='
c       WRITE(6,*) ' GSDNBB2: L block '
c       WRITE(6,*) ' ================='
c       CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
*
        WRITE(6,'(A,4(2X,I8))') ' NJA, NJB, NIA, NIB = ',
     &                            NJA, NJB, NIA, NIB
*
        WRITE(6,*)
        WRITE(6,'(A,12(2X,I4))') ' Occupation of alpha strings in L:',
     &  (IAOC(I),I=1,NGAS)
        WRITE(6,'(A,12(2X,I4))') ' Occupation of beta strings  in L:',
     &  (IBOC(I),I=1,NGAS)
        WRITE(6,'(A,12(2X,I4))') ' Occupation of alpha strings in R:',
     &  (JAOC(I),I=1,NGAS)
        WRITE(6,'(A,12(2X,I4))') ' Occupation of beta strings  in R:',
     &  (JBOC(I),I=1,NGAS)
        WRITE(6,*)
*
        WRITE(6,*) ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB
* 
        WRITE(6,*) 'SCLFAC =',SCLFAC
      END IF
C?    WRITE(6,*) ' RHO2 entering GSDNBB2 '
C?    CALL PRSYM(RHO2,NACOB**2)
*
C?    WRITE(6,*) ' IOBPTS_SEL in GSDNBB2'
C?    CALL IWRTMA(IOBPTS_SEL,NGAS,NSMOB,MXPNGAS,NSMOB)
      IACTIVE = 0
*
      IF(IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
*  beta contribution to RHO1
* =============================
*
        IF(NTEST.GE.200) 
     &  WRITE(6,*) ' GSBBD1 will be called (beta)'
        IAB = 2
        CALL GSBBD1(RHO1,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &       NGAS,IBOC,JBOC,
     &       SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB,
     &       NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?      WRITE(6,*) ' GSBBD1 was called '
C?      WRITE(6,*) ' Memory check '
C?      CALL MEMCHK
C?    WRITE(6,*) ' RHO2 in GSDNBB2 after bb(1el)'
C?    CALL PRSYM(RHO2,NACOB**2)
*
* ================================
* beta-beta contribution to RHO2
* ================================
*
        IF(I12.EQ.2.AND.NBEL.GE.2) THEN
        IF(NTEST.GE.200) 
     &  WRITE(6,*) ' GSBBD2A will be called (beta)'
          CALL GSBBD2A(RHO2,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &         NGAS,IBOC,JBOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,IDOSRHO2,RHO2BB,
     &         NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
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
C?    WRITE(6,*) ' RHO2 in GSDNBB2 after bb'
C?    CALL PRSYM(RHO2,NACOB**2)
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
        IF(NTEST.GE.200) 
     &  WRITE(6,*) ' GSBBD1 will be called (alpha)'
        IAB = 1
        CALL GSBBD1(RHO1,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &       NGAS,IAOC,JAOC,SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB,
     &         NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?        WRITE(6,*) ' GSBBD1 was called '
        IF(I12.EQ.2.AND.NAEL.GE.2) THEN
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
        IF(NTEST.GE.200) 
     &  WRITE(6,*) ' GSBBD2A will be called (alpha)'
          CALL GSBBD2A(RHO2,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &         NGAS,IAOC,JAOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,IDOSRHO2,RHO2AA,
     &         NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?        WRITE(6,*) ' GSBBD2A was called '
        END IF
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
C?    WRITE(6,*) ' RHO2 in GSDNBB2 after aa'
C?    CALL PRSYM(RHO2,NACOB**2)
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
        IF(NTEST.GE.200) 
     &  WRITE(6,*) ' GSBBD2B will be called '
        IUSEAB = 0
        CALL GSBBD2BN(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                    JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                    IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                    SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                    NOBPTS,IOBPTS,MAXK,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                    NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                    SSCR,CSCR,NACOB,NTEST,SCLFAC,S2_TERM1,
     &                    IDOSRHO2,RHO2AB,
     &                    NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
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
C?    WRITE(6,*) ' RHO2 in GSDNBB2 after ab'
C?    CALL PRSYM(RHO2,NACOB**2)
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
      IF(NTEST.GE.100) WRITE(6,*) ' Returning from GSDNBB2'
*. A test
C?    WRITE(6,*) ' RHO2 after GSBBD2B '
C?    CALL PRSYM(RHO2,NACOB**2)
*
      CALL LUCIAQEXIT('GSDN2')
      RETURN
      END
      SUBROUTINE GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,SB,CB,
     &                  ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,
     &                  IDOSRHO2,RHO2SS,
     &                  NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
*
* Contributions to two-electron density matrix from column excitations
*
* GAS version, '96 , Jeppe Olsen 
*              Sept. 04, 2-electron spin-densities added 
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
* RHO2SS : Updated alpha-alpha ( or beta-beta) 2e- spindensity (if IDOSRHO2=1)
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
*              Calculating densities only over selected spaces, Sept. 05
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nproc
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST), SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
      INTEGER IDACTSPC(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*.Output
      DIMENSION RHO2(*), RHO2SS(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
C     INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
*
      CALL LUCIAQENTER('GSBDA')
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
     &WRITE(6,*) ' ISCSM,ICCSM IDXSM ', ISCSM,ICCSM,IDXSM
c     nproc = ga_nnodes()
      DO 2000 IDXTP =  1, NDXTP
        ITYP = ITP(IDXTP)
        JTYP = JTP(IDXTP)
        KTYP = KTP(IDXTP)
        LTYP = LTP(IDXTP)
        IF(IDACTSPC(ITYP)+IDACTSPC(JTYP)+IDACTSPC(KTYP)+IDACTSPC(LTYP)
     &     .NE.4) GOTO 2000
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
          IOFF = IOBPTS_SEL(ITYP,ISM)
          JOFF = IOBPTS_SEL(JTYP,JSM)
          KOFF = IOBPTS_SEL(KTYP,KSM)
          LOFF = IOBPTS_SEL(LTYP,LSM)
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
c             if (npart.lt.nproc) then
c                npart = min(nrow,nproc)
c             endif
              IF(NTEST.GE.1000)
     &        write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
c             if (ga_nodeid()+1.gt.npart) goto 1930
c             DO 1801 IIPART = ga_nodeid()+1, NPART,nproc
              DO 1801 IIPART = 1, NPART

                IBOT = 1+(IIPART-1)*MAXI
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
                  CALL LUCIAQENTER('ADADS')
                  II12 = 1
                  K12 = 1
                  IONE = 1
                  CALL ADADST_GAS(IONE,JSM,JTYP,NJ,
     &                            IONE,LSM,LTYP,NL,
     &                        ICCTP,ICCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        JFRST,KFRST,II12,K12,SCLFAC)

                  JFRST = 0
                  KFRST = 0
*
                  CALL LUCIAQEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
                  CALL LUCIAQENTER('MATCG')
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
                  CALL LUCIAQEXIT ('MATCG')
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
CT                CALL LUCIAQENTER('ADADS')
                  II12 = 2
                  K12 = 1
                  IONE = 1
                  IF(IFRST.EQ.1) KFRST = 1
                  ONE = 1.0D0
                  CALL ADADST_GAS(IONE,ISM,ITYP,NI,
     &                            IONE,KSM,KTYP,NK,
     &                        ISCTP,ISCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        IFRST,KFRST,II12,K12,ONE   )
                  IFRST = 0
                  KFRST = 0
*
CT                CALL LUCIAQEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  I = 0
                  K = 1
*
CT                CALL LUCIAQENTER('MATCG')
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
CT                CALL LUCIAQEXIT ('MATCG')
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
              IOFF = IOBPTS_SEL(ITYP,ISM)
              JOFF = IOBPTS_SEL(JTYP,JSM)
              KOFF = IOBPTS_SEL(KTYP,KSM)
              LOFF = IOBPTS_SEL(LTYP,LSM)
              NTESTO = NTEST
C?            IF(IOFF.EQ.3.AND.JOFF.EQ.3.AND.KOFF.EQ.4.AND.LOFF.EQ.4)
C?   &            NTEST = 5000
                  LDUMMY = NKBTC*NIBTC
                  IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
                    WRITE(6,*) ' CSCR matrix '
                    CALL WRTMAT(CSCR,LDUMMY,NJL,LDUMMY,NJL)
                    WRITE(6,*) ' SSCR matrix '
                    CALL WRTMAT(SSCR,LDUMMY,NIK,LDUMMY,NIK)
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
                  IF(NTEST.GE.1000.and.ga_nodeid().eq.1) THEN
                    WRITE(6,*) ' Updated X matrix',IIPART
                    CALL WRTMAT(X,NIK,NJL,NIK,NJL)
                  END IF

*
                IF(KEND.EQ.0) GOTO 1800
*. End of loop over partitionings of resolution strings
 1801         CONTINUE
*. Rho2(ik,jl) has been constructed for ik,jl belonging to
*. Scatter out to density matrix
              IOFF = IOBPTS_SEL(ITYP,ISM)
              JOFF = IOBPTS_SEL(JTYP,JSM)
              KOFF = IOBPTS_SEL(KTYP,KSM)
              LOFF = IOBPTS_SEL(LTYP,LSM)
C?            WRITE(6,*) ' ITYP, ISM, IOFF = ', ITYP, ISM, IOFF
C?            WRITE(6,*) ' JTYP, JSM, JOFF = ', JTYP, JSM, JOFF
C?            WRITE(6,*) ' KTYP, KSM, KOFF = ', KTYP, KSM, KOFF
C?            WRITE(6,*) ' LTYP, LSM, LOFF = ', LTYP, LSM, LOFF
              CALL ADTOR2(RHO2,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,
     &                    NDACTORB)
C                  ADTOR2(RHO2,RHO2T,ITYPE,
C    &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
               IF(IDOSRHO2.EQ.1) THEN
*. and add to spin-density 
              CALL ADTOR2S(RHO2SS,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,
     &                    NDACTORB)
               END IF

 1930       CONTINUE
 1940     CONTINUE
 1950   CONTINUE
 2000 CONTINUE
 2001 CONTINUE
      CALL LUCIAQEXIT('GSBDA')
*
      RETURN
      END
      SUBROUTINE GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXK,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  CJRES,SIRES,NORB,NTESTG,SCLFAC,S2_TERM1,
     &                  IDOSRHO2,RHO2AB,
     &                  NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
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
      INTEGER IOBPTS_SEL(MXPNGAS,*), IDACTSPC(NGAS)
*.Input
      DIMENSION CB(*),SB(*)
*. Output
      DIMENSION RHO2(*),RHO2AB(*)
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
        IF(IDACTSPC(ITYP)+IDACTSPC(JTYP).NE.2) GOTO 2001
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS_SEL(ITYP,ISM)
          JOFF = IOBPTS_SEL(JTYP,JSM)
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
              IF(IDACTSPC(KTYP)+IDACTSPC(LTYP).NE.2) GOTO 2000
*
              DO 1930 KSM = 1, NSMOB
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
C?              WRITE(6,*) ' Loop 1930, KSM LSM ',KSM,LSM
                KOFF = IOBPTS_SEL(KTYP,KSM)
                LOFF = IOBPTS_SEL(LTYP,LSM)
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
C?            WRITE(6,*) ' ITYP, ISM, IOFF = ', ITYP, ISM, IOFF
C?            WRITE(6,*) ' JTYP, JSM, JOFF = ', JTYP, JSM, JOFF
C?            WRITE(6,*) ' KTYP, KSM, KOFF = ', KTYP, KSM, KOFF
C?            WRITE(6,*) ' LTYP, LSM, LOFF = ', LTYP, LSM, LOFF
                CALL ADTOR2(RHO2,X,2,
     &                NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NDACTORB)
                IF(IDOSRHO2.EQ.1) THEN
                  CALL ADTOR2S(RHO2AB,X,2,
     &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NDACTORB)
                END IF 
                IF(NTEST.GE.1000) THEN
                write(6,*) ' updated density matrix '
                call prsym(rho2,NDACTORB*NDACTORB)
                END IF

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
      SUBROUTINE GSBBD1(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &                  IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB,
     &                  NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen 
* Particle-Hole version of Jan. 98
* Active orbital spaces added Sept. 05
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
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nproc
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), ITSOB(*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
      INTEGER IDACTSPC(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*.Output
      DIMENSION RHO1(*), SRHO1(*)
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
*. Add or subtract for spindensity
      IF(IAB.EQ.1) THEN
        XAB = 1.0D0
      ELSE
        XAB = -1.0D0
      END IF
*.Local arrays
      NTEST = 000
      CALL LUCIAQENTER('GSBD1')
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
        WRITE(6,*) ' ISCSM, ICCSM = ', ISCSM, ICCSM
*         
        WRITE(6,*) ' GSBBD1, sclfac ',SCLFAC
      END IF
*
      IFRST = 1
      JFRST = 1
c     nproc = ga_nnodes()
*. Number of partitionings over column strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
c         if (nipart.lt.nproc) then
c            nipart = min(nrow,nproc)
c         endif
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
        IF(IDACTSPC(ITYP)+IDACTSPC(JTYP).NE.2) GOTO 900
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
          IBIORB = IOBPTS_SEL(ITYP,ISM)
          IBJORB = IOBPTS_SEL(JTYP,JSM)
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
c             DO 701 IIPART = ga_nodeid()+1, NIPART,nproc
              DO 701 IIPART = 1, NIPART
                IBOT = (IIPART-1)*MAXI+1
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
                NKI = LKABTC*NIBTC
                IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' CSCR and SSCR '
                 CALL WRTMAT(CSCR,IJ_DIM(2),NKI,IJ_DIM(2),NKI)
                 CALL WRTMAT(SSCR,IJ_DIM(1),NKI,IJ_DIM(1),NKI)
                END IF
*
*. And then the hard  work
                FACTORC = 0.0D0
                FACTORAB = 1.0D0
                CALL MATML7(RHO1S,SSCR,CSCR,IJ_DIM(1),IJ_DIM(2),NKI,
     &               IJ_DIM(1),NKI,IJ_DIM(2),FACTORC,FACTORAB,1)
*
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' Block to one-body density '
                  CALL WRTMAT(RHO1S,IJ_DIM(1),IJ_DIM(2),
     &                              IJ_DIM(1),IJ_DIM(2))
                END IF
*. Scatter out to complete matrix  
                DO JJORB = 1, IJ_DIM(2)
                  JORB = IJ_OFF(2)-1+JJORB
                  DO IIORB = 1, IJ_DIM(1)
                    IORB = IJ_OFF(1)-1+IIORB
                    RHO1((JORB-1)*NDACTORB+IORB) =
     &              RHO1((JORB-1)*NDACTORB+IORB) +
     &              RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)
                    IF(IDOSRHO1.EQ.1) THEN  
                      SRHO1((JORB-1)*NDACTORB+IORB) =
     &                SRHO1((JORB-1)*NDACTORB+IORB) +
     &                XAB*RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)   
                    END IF
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
C!    stop ' enforrced stop in RSBBD1 '
      CALL LUCIAQEXIT('GSBD1')
      RETURN
      END
      FUNCTION IB_H1(ISM,IHSM,NR,NC)
*
*. Offset to symmetryblock H(ISM,*)
*. Jeppe Olsen, Dec. 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*
      INTEGER NR(*),NC(*)
*
      IB = 1
      DO IISM = 1, ISM - 1
        JJSM = MULTD2H(IISM,IHSM)
        IB = IB + NR(IISM)*NC(JJSM)
      END DO
*
      IB_H1 = IB
*
      RETURN
      END
      SUBROUTINE MULT_H1H2(H1,IH1SM,H2,IH2SM,H12,IH12SM)
*. Two set of one-electron integrals H1 and H2 are 
*. given as symmetrypacked complete quadratic matrices. 
*. Obtain product as symmetrypacked complete quadratic matrix
*
* Jeppe Olsen, Dec. 2000
*
      INCLUDE 'implicit.inc'
*.Input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      DIMENSION H1(*),H2(*)
*. Output
      DIMENSION H12(*)
*. 
      INCLUDE 'multd2h.inc'
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Entering MULT_H1H2 '
        WRITE(6,*)' IH1SM, IH2SM = ', IH1SM, IH2SM
      END IF
      IH12SM = MULTD2H(IH1SM,IH2SM)
*. Loop over symmetry blocks of H1H2
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IH12SM)
*. Connecting symmetry A(ISM,KSM)B(KSM,JSM)
        KSM = MULTD2H(ISM,IH1SM)
C?      WRITE(6,*) ' ISM, JSM, KSM =', ISM, JSM, KSM
*. Offsets to A(ISM,KSM) and B(KSM,JSM)
C              IB_H1(ISM,IHSM,NR,NC)
        IB_A = IB_H1(ISM,IH1SM,NTOOBS,NTOOBS)
        IB_B = IB_H1(KSM,IH2SM,NTOOBS,NTOOBS)
        IB_AB = IB_H1(ISM,IH12SM,NTOOBS,NTOOBS)
C?      WRITE(6,*) ' IB_A, IB_B, IB_AB =', IB_A,IB_B,IB_AB
        
*
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)
        NK = NTOOBS(KSM)
*
        FACTORC = 0.0D0
        FACTORAB = 1.0D0
        ZERO = 0.0D0
        CALL SETVEC(H12(IB_AB),ZERO,NI*NJ)
        CALL MATML7(H12(IB_AB),H1(IB_A),H2(IB_B),NI,NJ,NI,NK,NK,NJ,
     &              FACTORC,FACTORAB,0)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input matrices to MULTH1H2 '
        CALL PRHONE(H1,NTOOBS,IH1SM,NSMOB,0)
        CALL PRHONE(H2,NTOOBS,IH2SM,NSMOB,0)
        WRITE(6,*) ' Output matrix from MULTH1H2 '
        CALL PRHONE(H12,NTOOBS,IH12SM,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE ABEXP2(A,IASM,B,IBSM,AB)
*
* Evaluate contribution from active orbitals to expectation value of product of 
* two one-electron operators 
*
* <0!AB!0> = sym(ijkl) A(ij)B(kl) d(ijkl) + sum(ij) rho1(ij) (AB)(ij)
*
* Jeppe Olsen Dec. 4, 2000 in Helsingfors 
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*
      REAL*8 INPROD
*. Input, A and B are required to be in symmetrypacked complete form
      DIMENSION A(*),B(*) 
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ABEXP2 '
        WRITE(6,*) ' ================ '
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ABEXP2')
*. Largest number of orbitals of given sym
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Allocate memory     
      LB_RHO2 = MXSOB**4
      CALL MEMMAN(KLRHO2B,LB_RHO2,'ADDL  ',2,'RHO2B ') !done
      L_VEC = NTOOB**2
      CALL MEMMAN(KLVEC,L_VEC,'ADDL  ',2,'VEC   ') !done
*
      CALL MEMMAN(KLABLK,MXSOB**2,'ADDL  ',2,'ABLK  ') !done
      CALL MEMMAN(KLBBLK,MXSOB**2,'ADDL  ',2,'BBLK  ') !done
*
*. Two-electron contributions
*
      AB2 = 0.0D0
      DO ISM = 1,  NSMOB
       JSM = MULTD2H(ISM,IASM)
       DO IGAS = 1, NGAS
       DO JGAS = 1, NGAS
        NI = NOBPTS(IGAS,ISM)
        NJ = NOBPTS(JGAS,JSM)
        NIJ = NI*NJ
*. Obtain block A(ISM,IGAS,JSM,JGAS)
C       EXTR_SYMGAS_BLK_FROM_ORBMAT(A,ABLK,ISM,IGAS,JSM,JGAS)
        CALL EXTR_SYMGAS_BLK_FROM_ORBMAT
     &  (A,dbl_mb(KLABLK),ISM,IGAS,JSM,JGAS)
        DO KSM = 1, NSMOB
         LSM = MULTD2H(KSM,IBSM)
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          NK = NOBPTS(KGAS,KSM)
          NL = NOBPTS(LGAS,LSM)
          NKL = NK*NL
*. Obtain block B(KSM,KGAS,LSM,LGAS)
          CALL EXTR_SYMGAS_BLK_FROM_ORBMAT
     &    (B,dbl_mb(KLBBLK),KSM,KGAS,LSM,LGAS)
*. Fetch RHO2(ISM,JSM,KSM,LSM)
C         GETD2(RHO2B,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,ISPC)
          CALL GETD2
     &    (dbl_mb(KLRHO2B),ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,1)
* sum(kl) RHO2(ij,kl) B(kl)
          FACTORC = 0.0D0
          FACTORAB = 1.0D0
          ZERO = 0.0D0
          CALL SETVEC(dbl_mb(KLVEC),ZERO,NIJ)
          CALL MATML7(dbl_mb(KLVEC),dbl_mb(KLRHO2B),dbl_mb(KLBBLK),NIJ,
     &                1,NIJ,NKL,NKL,1,FACTORC,FACTORAB,0)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' KLVEC: '
            CALL WRTMAT(dbl_mb(KLVEC),1,NIJ,1,NIJ)
          END IF
*. sum(ij) (A(ij) (sum(kl) RHO2(ij,kl)B(kl) )
          AB2 = AB2 + INPROD(dbl_mb(KLABLK),dbl_mb(KLVEC),NIJ)
          IF(NTEST.GE.100) THEN
            WRITE(6,'(A,4I4,E20.12)') 
     &      ' Two-electron contribution after ISM, JSM, KSM, LSM = ',
     *      ISM, JSM, KSM, LSM, AB2
          END IF
         END DO ! loop over LGAS
         END DO ! loop over KGAS
        END DO ! loop over KSM 
       END DO ! Loop over JGAS
       END DO ! Loop over IGAS
      END DO ! Loop over ISM
      IF(NTEST.GE.10) 
     &WRITE(6,*) ' Two-electron contribution to <0!AB!0> ', AB2
*
*. One-electron part
*
* AB(IJ) RHO1(IJ)
C     MULT_H1H2(H1,IH1SM,H2,IH2SM,H12,IH12SM)
      CALL MULT_H1H2(A,IASM,B,IBSM,dbl_mb(KLVEC),IABSM)
      AB1 = 0.0D0
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IABSM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          NI = NOBPTS(IGAS,ISM) 
          NJ = NOBPTS(JGAS,JSM) 
          NIJ = NI*NJ
          IF(NTEST.GE.100) THEN
            WRITE(6,'(A,4I4)') ' IGAS, JGAS, ISM, JSM = ', 
     &                           IGAS, JGAS, ISM, JSM 
          END IF
*. Obtain block of AB
          CALL EXTR_SYMGAS_BLK_FROM_ORBMAT
     &    (dbl_mb(KLVEC),dbl_mb(KLBBLK),ISM,IGAS,JSM,JGAS)
*. Obtain density RHO1(ISM,JSM)
C         GETD1(RHO1B,ISM,IGAS,JSM,JGAS,ISPIN)
          CALL GETD1(dbl_mb(KLRHO2B),ISM,IGAS,JSM,JGAS,1)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Block of D1 '
            CALL WRTMAT(dbl_mb(KLRHO2B),NI,NJ,NI,NJ)
          END IF
          AB1 = AB1 + INPROD(dbl_mb(KLRHO2B),dbl_mb(KLBBLK),NIJ)
        END DO ! loop over JGAS
        END DO ! loop over IGAS
      END DO ! Loop over ISM
      IF(NTEST.GE.10) 
     &WRITE(6,*) ' One-electron contribution to <0!AB!0> ', AB1
*
      AB = AB1 + AB2
      IF(NTEST.GE.10) 
     &WRITE(6,*) ' <0!AB!0> = ', AB
      CALL MEMCHK2('AFT1EL')
*
      CALL MEMMAN(IDUM,IDUMN,'FLUSM ',IDUM,'ABEXP2')
      RETURN
      END
C        LLBATCHB,LLBATCHE,LBLOCK,LLI1BATCH,ICBLOCK,ISBLOCK,
C LLBATCHB => ICBLBT
C LLBATCHE => ICLEBT 
C LBLOCK => ICLBLK
C LLI1BATCH => ICI1BT
C 
      SUBROUTINE PICO4(VEC1,VEC2,LU1,LU2,LU3,LU4,RNRM,EIG,FINEIG,MAXIT,
     &                 NBATCH,
     &                 ICBLBT,ICLEBT,ICLBLK,ICI1BT,ICBLOCK,
     &                 ISBLBT,ISLEBT,ISLBLK,ISI1BT,ISBLOCK,
     &                 IPRTXX,
     &                 NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                 THRES_ET,THRES_EC,THRES_CC,
     &                 E_CONV,C_CONV,ICLSSEL,
     &                 IBLK_TO_CLS,NCLS,CLS_C,CLS_E,CLS_CT,CLS_ET,
     &                 ICLS_A,ICLS_L,RCLS_L,IBLKS_A,
     &                 CLS_DEL,CLS_DELT,CLS_GAMMA,CLS_GAMMAT,ISKIPEI,
     &                 I2BLK,VEC3,ICLS_A2,MXLNG,IBASSPC,EBASC,CBASC,
     &                 NSPC,IMULTGR,IPAT,LPAT,IREFSPC,CONVER)
*
* Davidson algorithm , requires three blocks in core
*
* Only three vectors in on DISC
*
* Lu4 should only hold a batch of coefficients 
*
*
*
* Block processing version
*
* Jeppe Olsen Winter of 1996
*
* Last modification; Nov. 7, 2012; Jeppe Olsen, Aligned with modern code

*
* Initial version - Only Diagonal preconditioner,
*
* Special version for NROOT = 1, MAXVEC = 2 !!
*
* Input :
* =======
*        LU1 : Initial  vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU2,LU3   : Scatch files
*        MAXIT     : Largest allowed number of iterations
*        NBATCH    : Number of batches of vector
*        ICBLBT  : Number of blocks in each batch
*        ICLEBT  : Number of elements  in each batch
*        ICLBLK    : Length of each block, packed
*        ICI1BT : First block of a given batch
*        ICBLOCK   : Some additional informaition about the blocking 
*        ISBLOCK  :  that this routine does not care about !!!!
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
*   THRES_ET   : Convergence criteria for eigenvalue
*
*   THRES_EC   : Threshold for second order energies for individual terms
*   THRES_CC   : Threshold for first  order wavefunction  for individual terms
*                
*
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL * 8   INPROD, INPRDD, INPRODB
*
      DIMENSION ICBLBT(NBATCH),ICLEBT(NBATCH),
     &          ICI1BT(NBATCH), ICLBLK(*), ICBLOCK(8,*)
      DIMENSION ISBLBT(NBATCH),ISLEBT(NBATCH),
     &          ISI1BT(NBATCH), ISLBLK(*), ISBLOCK(8,*)
*. Cands is needed to save NCBATCH, NSBATCH, not pretty...
      INCLUDE 'cands.inc'
*. I2BLK should atleast have length of number of
*. blocks 
      DIMENSION I2BLK(*)
*
      DIMENSION RNRM(MAXIT,1),EIG(MAXIT,1),FINEIG(1)
      DIMENSION H0(*),IPNTR(1)
      DIMENSION H0SCR(*)
*. Class and block information 
      DIMENSION IBLK_TO_CLS(*) 
      DIMENSION CLS_C(NCLS),CLS_E(NCLS),CLS_CT(NCLS),CLS_ET(NCLS)
      DIMENSION CLS_DEL(*), CLS_DELT(*)
      DIMENSION CLS_GAMMA(*), CLS_GAMMAT(*)
*. Base CI spaces : CI space where a given class is introduced 
      DIMENSION IBASSPC(*),EBASC(*),CBASC(*)
      DIMENSION IPAT(*)
*.Initial VEC3      
      DIMENSION VEC3(*)
      INTEGER ICLS_A(NCLS), ICLS_L(NCLS),IBLKS_A(*),ICLS_A2(NCLS)
      DIMENSION RCLS_L(NCLS)
*
*     H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
      LOGICAL CONVER
*
C?    WRITE(6,*) ' Memchk at start of PICO4'
C?    CALL MEMCHK
      IPICO = 0
      IF(IPICO.NE.0) THEN
C?      WRITE(6,*) ' Perturbative solver '
C       MAXVEC = MIN(MAXVEC,2)
      ELSE IF(IPICO.EQ.0) THEN
C?      WRITE(6,*) ' Variational  solver '
      END IF
*
      WRITE(6,*) ' Number of spaces ', NSPC
      WRITE(6,*) ' Map : Class => Base space '
      CALL IWRTMA(IBASSPC,1,NCLS,1,NCLS)
*
      IF(ICLSSEL.EQ.1) THEN
        WRITE(6,*) ' Class selection will be performed '
        WRITE(6,*) ' Number of classes ', NCLS
        WRITE(6,*) ' Dimension of each class ( Integer )'
        CALL IWRTMA(ICLS_L,1,NCLS,1,NCLS)
        WRITE(6,*) ' Dimension of each class ( Real )'
        CALL WRTMAT(RCLS_L,1,NCLS,1,NCLS)
      END IF
*
      WRITE(6,*) ' IMULTGR = ', IMULTGR
      WRITE(6,*) ' EIGSHF = ', EIGSHF
      IF(IMULTGR.NE.0) THEN
        WRITE(6,*) ' Multispace method in use '
        WRITE(6,*)
        WRITE(6,*) ' Length of pattern ', LPAT
        WRITE(6,*) ' Pattern : '
        CALL IWRTMA(IPAT,1,LPAT,1,LPAT)
        WRITE(6,*) 
        WRITE(6,*) ' Reference space ', IREFSPC
      END IF
*
      IPRT = 20
      IOLSTM = 1
      IF(IPRT.GT.10.AND.IOLSTM.NE.0)
     &WRITE(6,*) ' Inverse iteration modified Davidson '
      IF(IPRT.GT.10.AND.IOLSTM.EQ.0)
     &WRITE(6,*) ' Normal Davidson method '
*
C?    WRITE(6,*) ' LU1 LU2 LU3 LU4 = ', LU1,LU2,LU3,LU4
      IF(IPRT.GE.20) THEN
        WRITE(6,*) ' Convergence threshold for eigenvalues', THRES_ET
        WRITE(6,*)
        WRITE(6,*) ' Elements of trial vectors discarded if '
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        WRITE(6,*) 
     &  '    Estimate of contribution to wavefunction is less than ',
     &  THRES_CC
        WRITE(6,*) 
     &  '    Estimate of contribution to Energy is less than ',
     &  THRES_EC
      END IF
      WRITE(6,*)
      IF(IPRT.GE.10)
     &WRITE(6,*) ' Max number of batches of vector ', NBATCH
*
*. Total number of blocks
      NBLOCKT = 0
      DO IBATCH = 1, NBATCH
        NBLOCKT = NBLOCKT + ICBLBT(IBATCH)
      END DO
      IF(IPRT.GE.10)
     &WRITE(6,*) ' Max number of blocks ', NBLOCKT
      WRITE(6,*)
      TEST = 1.0D-6
      CONVER = .FALSE.
      IROOT = 1
      NROOT = 1
      ZERO = 0.0D0
      MAX_MICRO_IT = 1
*.    ^ Should be moved outside 
*
*. Play around with dynamic allocation of batches ...
      IDYNBATCH = 1
*
* ===================
*.Initial iteration
* ===================
*
      IF(MAXIT.EQ.0) THEN
        WRITE(6,*) ' Max number of iterations is zero'
        WRITE(6,*) ' I will just return from PICO4'
        RETURN
      END IF
*
      ITER = 1
      IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==============================='      
        WRITE(6,*) ' Info from iteration .... ', ITER
        WRITE(6,*) ' ==============================='      
        WRITE(6,*)
      END IF
      CALL LUCIAQENTER('INIIT')
*. Obtain energy of initial vector
      IF(ISKIPEI.EQ.0) THEN
*. active classses of initial vector
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU1,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
*. Mark active blocks and find required number of batches, sigma = c 
          CALL REPART_CIV(ICBLOCK,NBATCHL,ICBLBT,ICLEBT,ICI1BT,
     &         MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ICLBLK)
          CALL REPART_CIV(ISBLOCK,NBATCHL,ISBLBT,ISLEBT,ISI1BT,
     &         MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ISLBLK)
*. NBATCHL is number of batches of Sigma
          NCBATCH = NBATCHL
          NSBATCH = NBATCHL
        ELSE
          NBATCHL = NBATCH
          NCBATCH = NBATCHL
          NSBATCH = NBATCHL
        END IF
*
        WRITE(6,*) ' <0!H!0> to be calculated '
        IF(IPRT.GE.10) WRITE(6,*) ' Number of batches for C ', NBATCHL
        E = 0.0D0
        IRESTRICT = 1
        IOFF = 1
*. We are going to calculate half of sigma and then multiply - correct
        EIGSHFD2 = EIGSHF*0.5D0
        DO IBATCH = 1, NBATCHL
          LBATCHB = ICBLBT(IBATCH)
          LBATCHE = ICLEBT(IBATCH)
          IF(IPRT.GE.10) 
     &    WRITE(6,*) '  <Batch 0!H!0>, Batch : ',IBATCH
          CALL SETVEC(VEC1,ZERO,LBATCHE)
          CALL SBLOCK(LBATCHB,ISBLOCK,IOFF,VEC2,VEC1,LU1,IRESTRICT,0,
     &                0,0,0,0.0D0,EIGSHFD2,'SIGMA ')
C         SBLOCK(NBLOCK,ISBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
C    &                  ECORE,ITASK)

          IF(IPRT.GE.200) THEN
            WRITE(6,*) ' Initial batch of S, number', IBATCH
            CALL WRTBLKN_EP(VEC1,LBATCHB,ICLBLK(IOFF))
          END IF
*. Obtain corresponding C blocks
          CALL GET_BLOCKS_FROM_DISC
     &    (LU1,LBATCHB,IOFF,ICBLOCK,NBLOCKT,VEC2,1)
          E = E + INPROD(VEC1,VEC2,LBATCHE)
          IF(IPRT.GE.200) THEN
            WRITE(6,*) ' Initial batch of C, number', IBATCH
            CALL WRTBLKN_EP(VEC2,LBATCHB,ICLBLK(IOFF))
          END IF
          IOFF = IOFF + LBATCHB
        END DO
        IF(IRESTRICT.EQ.1) E = 2*E
        EIG(1,IROOT) = E                                 
        WRITE(6,*) ' Energy = ', E
      ELSE IF(ISKIPEI.EQ.1) THEN
        E = FINEIG(IROOT) 
        WRITE(6,*) ' Initial energy obtained from previous calc as ',
     &             E
        EAPR = E
        EIG(1,IROOT) = E                                 
      END IF
        
*
      IF(IPRT .GE. 3 ) THEN
        WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
        WRITE(6,'(5F18.13)')
     &  ( EIG(1,IROOT),IROOT=1,NROOT)
      END IF
      ITERX = 1
      CALL LUCIAQEXIT('INIIT')
*
* ======================
*. Loop over iterations
* ======================
*
      DO 1000 ITER = 2, MAXIT
       IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==============================='      
        WRITE(6,*) ' Info from iteration .... ', ITER
        WRITE(6,*) ' ==============================='      
        WRITE(6,*)
       END IF
*. Allow loop over micro-iterations : 
*. In the first micro of a given iteration, the complete Sigma-vector
*. is calculated, in the following micro's, the Sigma vector is 
*. restricted to the space of the C-vectors
      DO 999, MICRO_IT = 1, MAX_MICRO_IT
*
       IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================================='      
        WRITE(6,*) ' Info from micro iteration .... ', MICRO_IT
        WRITE(6,*) ' ====================================='      
        WRITE(6,*)
       END IF
       CALL LUCIAQENTER('PARTA')
*. Largest allowed basespace in this iteration in multispace method
       IF(IMULTGR.GT.0) THEN
         IBASSPC_MX = MAX(1,1-IPAT(MOD(ITER-2,LPAT)+1)+IREFSPC)
         IF(ITER.EQ.MAXIT) THEN
           IBASSPC_MX = IREFSPC
         END IF
         WRITE(6,*) ' Max allowed base space ', IBASSPC_MX
       ELSE
         IBASSPC_MX = 0         
       END IF
*
       ZERO = 0.0D0
       IF(ICLSSEL.EQ.1) THEN
         CALL SETVEC(CLS_CT,ZERO,NCLS)
         CALL SETVEC(CLS_ET,ZERO,NCLS)
         CALL SETVEC(CLS_DELT,ZERO,NCLS)
         CALL SETVEC(CLS_C ,ZERO,NCLS)
         CALL SETVEC(CLS_E ,ZERO,NCLS)
         CALL SETVEC(CLS_DEL,ZERO,NCLS)
         CALL SETVEC(CLS_GAMMA,ZERO,NCLS)
*
         CALL SETVEC(EBASC,ZERO,NSPC)
         CALL SETVEC(CBASC,ZERO,NSPC)
*
       END IF
*. Active classes
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU1,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
*. Mark active blocks and find required number of batches
          CALL REPART_CIV(ICBLOCK,NBATCHL,ICBLBT,ICLEBT,ICI1BT,
     &         MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ICLBLK)
          CALL REPART_CIV(ISBLOCK,NBATCHL,ISBLBT,ISLEBT,ISI1BT,
     &         MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ISLBLK)
          NCBATCH = NBATCHL
        ELSE
          NBATCHL = NBATCH
          NCBATCH = NBATCHL
        END IF
       IF(IPRT.GE.10) WRITE(6,*) ' Number of batches for C ', NBATCHL
       EIGAPR = E
       CALL LUCIAQEXIT('PARTA')
* ===============================
*. Obtain C(T) (H0-E)**-1 C
* ===============================
       CALL LUCIAQENTER('PARTB')
       GAMMA = 0.0D0
       IOFF = 1
       CALL REWINO(LU1)
       DO IBATCH = 1, NBATCHL
         LBATCHB = ICBLBT(IBATCH)
         LBATCHE = ICLEBT(IBATCH)
*. Retrieve Batch of C
         NO_ZEROING = 0
         NO_ZEROING1= 1
         CALL FRMDSCN3(VEC1,LBATCHB,LBLK,LU1,NO_ZEROING1,I2BLK(IOFF),
     &                 ICLBLK(IOFF))
         CALL COPVEC(VEC1,VEC2,LBATCHE)
*. Multiply with (H0-E)** -1 
         FACTOR = -EIGAPR
         ITASK = 1
         CALL DIATERM_GAS(FACTOR,ITASK,VEC2,LBATCHB,ICBLOCK,IOFF,0,
     &         NO_ZEROING1,I2BLK(IOFF))
         IF(NO_ZEROING1.EQ.0) THEN
           GAMMA = GAMMA + INPROD(VEC1,VEC2,LBATCHE)
         ELSE
            GAMMA = GAMMA + INPRODB(VEC1,VEC2,LBATCHB,ICLBLK(IOFF),
     &                      I2BLK(IOFF))
         END IF
         IF(ICLSSEL.EQ.1) THEN
* Contributions to Gamma as <O(Occclass)!(H0-E)**(-1)!0>
           CALL CLASS_PROD3(VEC1,VEC2,IOFF,LBATCHB,ICBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_GAMMA)
         END IF
         IOFF = IOFF + LBATCHB
       END DO
       IF(IPRT.GE.20)
     & WRITE(6,*) ' Gamma  calculated ',GAMMA
       IF(ICLSSEL.EQ.1 .AND. IPRT.GE.100) THEN
         WRITE(6,*) ' Contributions to Gamma from various classes '
         CALL WRTMAT(CLS_GAMMA,1,NCLS,1,NCL)
       END IF
       CALL LUCIAQEXIT('PARTB')

*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = (H0-E)-1 (H*C - EIGAPR*C) : Obtain in batches and save on DISC
*
*
       WRITE(6,*)
       WRITE(6,*) ' ========================'
       WRITE(6,*) ' Delta under construction'
       WRITE(6,*) ' ========================'
       WRITE(6,*)
C      EIGAPR = EIG(ITER-1,1)
       CALL LUCIAQENTER('PARTC')
       XNELMNT = 0.0D0
       XNZERO = 0.0D0
*
       RNORM = 0.0D0
       DELTA = 0.0D0
       CHEDEL =0.0D0
       DELNORM = 0.0D0
*
       DELTAT = 0.0D0
       CHEDELT =0.0D0
       DELNORMT = 0.0D0
       ECC = 0.0D0
*
       CALL REWINO(LU2)
       IOFF = 1
*. Find partitioning of sigma
       IF(MICRO_IT.EQ.1) THEN
*. BLocks obttained by double excitations from classes in C
         NEXC  = 2
       ELSE  
*. Keep sigma-vector equal to C vector 
         NEXC = 0
       END IF
*
       WRITE(6,*) ' Active classes of Sigma '
       CALL EXCCLS2(NCLS,ICLS_A,ICLS_A2,NEXC,IBASSPC_MX,IBASSPC)
*. Partitioning of sigma vector      
       CALL REPART_CIV(ISBLOCK,NBATCHL,ISBLBT,ISLEBT,ISI1BT,
     &      MXLNG,ICLS_A2,IBLK_TO_CLS,NCLS,NBLOCKT,ISLBLK)
       NSBATCH = NBATCHL
    
       IF(IPRT.GE.10) THEN
         WRITE(6,*) ' Number of batches for S ', NBATCHL
         WRITE(6,*)
       END IF
       CALL LUCIAQEXIT('PARTC')
*. Loop over batches of H !0>
       DO IBATCH = 1, NBATCHL
         LBATCHB = ISBLBT(IBATCH)
         LBATCHE = ISLEBT(IBATCH)
         XNELMNT = XNELMNT + FLOAT(LBATCHE)
         IF(IPRT.GE.10) 
     &   WRITE(6,*) '  < Batch delta ! H !0>, Batch : ',IBATCH
* Batch of HC in VEC1
         ZERO = 0.0D0
         CALL SETVEC(VEC1,ZERO,LBATCHE)
         CALL LUCIAQENTER('PARTD')
         CALL SBLOCK(LBATCHB,ISBLOCK,IOFF,VEC2,VEC1,LU1,0,0,
     &                0,0,0,0.0D0,EIGSHF,'SIGMA ')
C     SBLOCK(NBLOCK,ISBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
C    &                  ECORE,ITASK)

         CALL LUCIAQEXIT('PARTD')
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' Batch of H C '
           CALL WRTBLKN_EP(VEC1,LBATCHB,ISLBLK(IOFF)) 
         END IF
*. Retrieve Batch of C
         CALL LUCIAQENTER('PARTE')
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,LBATCHB,IOFF,ISBLOCK,NBLOCKT,VEC2,1)
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' C batch read in '
           CALL WRTBLKN_EP(VEC2,LBATCHB,ISLBLK(IOFF))
         END IF
*. Update energy
         ECC = ECC + INPROD(VEC1,VEC2,LBATCHE)
* Batch of (H-E)C in VEC1
         ONE = 1.0D0
         FACTOR = -EIGAPR
         CALL VECSUM(VEC1,VEC1,VEC2,ONE,FACTOR,LBATCHE)
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' Batch of (H - E ) C '
           CALL WRTBLKN_EP(VEC1,LBATCHB,ISLBLK(IOFF))
         END IF
*. Norm of residual
         RNORM = RNORM + INPROD(VEC1,VEC1,LBATCHE)
*. Batch of (H0-E)-1(H-E)C  in VEC2
         CALL COPVEC(VEC1,VEC2,LBATCHE)
         FACTOR = -EIGAPR
         ITASK = 1
         I12 = 1
         CALL DIATERM_GAS(FACTOR,ITASK,VEC2,LBATCHB,ISBLOCK,IOFF,0,0,0)
         DELNORMT = DELNORMT + INPROD(VEC2,VEC2,LBATCHE)
* C(H-E)(H0-E0)-1(H-E)C
         CHEDELT = CHEDELT + INPROD(VEC1,VEC2,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
*. Contributions divided into occupation classes, complete expansion
*. Wave function correction
           CALL CLASS_PROD3(VEC2,VEC2,IOFF,LBATCHB,ISBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_CT)
*. Energy correction
           CALL CLASS_PROD3(VEC1,VEC2,IOFF,LBATCHB,ISBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_ET)
         END IF
*.[(H0-E0)-1(H-E)C]_{truncated} and 
*.C(H-E) [(H0-E0)-1(H-E)C]_{truncated}
         ZERO = 0.0D0
         CALL SETVEC(VEC3,ZERO,LBATCHE)
         XNZERO = XNZERO + FLOAT(LBATCHE)
         DO I = 1, LBATCHE
           IF(ABS(VEC2(I)*VEC1(I)).GE.THRES_EC.OR.
     &        ABS(VEC2(I)).GE.THRES_CC           ) THEN
             CHEDEL = CHEDEL + VEC1(I)*VEC2(I)
             DELNORM = DELNORM + VEC2(I)*VEC2(I)
             VEC3(I)=VEC2(I)
             XNZERO = XNZERO - 1.0D0
           END IF
         END DO
         CALL LUCIAQEXIT('PARTE')
*
         CALL LUCIAQENTER('PARTF')
         IF(ICLSSEL.EQ.1) THEN
*. Contributions divided into occupation classes, truncated expansion
*. Wave function correction
           CALL CLASS_PROD3(VEC3,VEC3,IOFF,LBATCHB,ISBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_C)
*. Energy correction
           CALL CLASS_PROD3(VEC1,VEC3,IOFF,LBATCHB,ISBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_E)
         END IF
*. retrieve c batch from disc
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,LBATCHB,IOFF,ISBLOCK,NBLOCKT,VEC1,1)
* C(H0-E0)-1(H-E)C
         DELTAT = DELTAT + INPROD(VEC1,VEC2,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
           CALL CLASS_PROD3(VEC1,VEC2,IOFF,LBATCHB,ISBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_DELT)
         END IF
* C[(H0-E0)-1(H-E)C]{truncated}
         DELTA = DELTA + INPROD(VEC1,VEC3,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
           CALL CLASS_PROD3(VEC1,VEC3,IOFF,LBATCHB,ISBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_DEL)
         END IF
*
*. Write packed version to DISC
*. Pack out so zero blocks are given zero entries
         CALL TODSCNP(VEC3,LBATCHB,ISLBLK(IOFF),LBLK,LU2)
         IF(IPRT.GE.200) THEN
           WRITE(6,*) ' Batch of blocks of trial vector '
           CALL WRTBLKN_EP(VEC2,LBATCHB,ISLBLK(IOFF))
         END IF
         IOFF = IOFF + LBATCHB
         CALL LUCIAQEXIT('PARTF')
       END DO
*      /\ End of loop over batches of correction vector
         CALL LUCIAQENTER('PARTG')
         WRITE(6,'(A,E22.15)') 
     &   ' Number of zero elements in delta (before class trunc) ',
     &    XNZERO
         WRITE(6,'(A,E22.15)') 
     &   ' Number of nonzero terms in delta (before class trunc) ',
     &    XNELMNT-XNZERO
       WRITE(6,'(A,E25.16)') 
     & ' ECC(may not be energy...) = ',ECC
       CALL ITODS(-1,1,LBLK,LU2)
       IF(MICRO_IT.EQ.1 ) RNRM(ITER-1,1) = SQRT(RNORM)
*. (End of loop over batches of (H0-E)-1(H-E)C)
*. Predicted energy  
       IF(ICLSSEL.EQ.1) THEN
*. Energy and wave function per base space
         CALL LUCIAQENTER('G1   ')
         CALL CLS_TO_BASE(CLS_E,EBASC,CLS_C,CBASC,NCLS,NSPC,
     &                    IBASSPC) 
         CALL LUCIAQEXIT('G1   ')
*. decide which classes should be truncated
         CALL LUCIAQENTER('PARTT')
C        CALL CLASS_TRUNC(NCLS,ICLS_L,RCLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
C    &                    E_CONV,ICLS_A,N_CLS_TRN,E_CLS_TRN,W_TRN)
*. Test active classes of Sigma (ICLS_A2)
         CALL CLASS_TRUNC(NCLS,ICLS_L,RCLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
     &                    E_CONV,ICLS_A2,N_CLS_TRN,E_CLS_TRN,W_TRN)
         CALL REPART_CIV(ISBLOCK,NBATCHL,ISBLBT,ISLEBT,ISI1BT,
     &        MXLNG,ICLS_A2,IBLK_TO_CLS,NCLS,NBLOCKT,ISLBLK)
         NSBATCH = NBATCHL
*. Perfrom the corresponding partitioning of the Sigma-vector
         CALL LUCIAQEXIT('PARTT')
*. Update corrections for class elimination
         N_ACT_CLS = 0
         CHEDEL = CHEDEL - (-1.0D0) * E_CLS_TRN
         WRITE(6,*) ' Initial Delta, Gamma = ', Delta, GAMMA
         DO JCLS = 1, NCLS
           IF(ICLS_A2(JCLS).EQ.0) THEN
             DELTA = DELTA-CLS_DEL(JCLS) 
             GAMMA = GAMMA-CLS_GAMMA(JCLS) 
           ELSE
             N_ACT_CLS = N_ACT_CLS + 1
           END IF
         END DO
         WRITE(6,*) ' Delta,  Gamma after class trunc', DELTA, GAMMA
           
*. And do the truncation 
*. Truncation of classes => truncation of blocks
         IF(N_CLS_TRN.NE.0) THEN
           CALL LUCIAQENTER('G2   ')
           CALL CLS_TO_BLK(NBLOCKT,IBLK_TO_CLS,ICLS_A,IBLKS_A)
           CALL LUCIAQEXIT('G2   ')
*. from LU2 to LU3 and back to LU2
           CALL LUCIAQENTER('G3   ')
           CALL ZAP_BLOCK_VEC(LU2,LBLK,IBLKS_A,VEC2,LU3)
           CALL LUCIAQEXIT('G3   ')
         END IF
       ELSE
         N_CLS_TRN = 0
         E_CLS_TRN = 0.0D0
         W_CLS_TRN = 0.0D0
       END IF
*. Predicted energy  
       IF(GAMMA.NE.0.0D0) THEN
         E2PREDIT = - CHEDELT + DELTAT**2/GAMMA
         E2PREDI  = - CHEDEL  + DELTA * DELTAT /GAMMA
       ELSE
         E2PREDIT = - CHEDELT
         E2PREDI  = - CHEDEL 
         IF(ICLSSEL.EQ.1) THEN
           CALL ICOPVE(CLS_CT,CLS_C,NCLS)
           CALL ICOPVE(CLS_ET,CLS_E,NCLS)
         END IF
       END IF
*
       WRITE(6,*)
*. 
C?     WRITE(6,*) ' Information for untruncated expansion:'
C?     WRITE(6,*) ' ======================================'
C?     WRITE(6,*) 
C?   & ' CHEDELT DELTAT GAMMA ', CHEDELT,DELTAT,GAMMA
       IF(GAMMA.NE.0.0D0) THEN
         WRITE(6,*) 
     & ' Orthogonalization term to E2 (no trunc.)', DELTAT**2/GAMMA
       END IF
       WRITE(6,'(A,2E25.15)') 
     & ' Predicted energy(no truncation), change and total ', 
     &              E2PREDIT,EIGAPR+E2PREDIT
       IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
C?       WRITE(6,*) 
C?       WRITE(6,*) ' Information for truncated expansion:'
C?       WRITE(6,*) ' ======================================'
C?       WRITE(6,*) 
C?   &   ' CHEDEL DELTA GAMMA ', CHEDEL,DELTA,GAMMA
         WRITE(6,*)
         IF(GAMMA.NE.0.0D0) THEN
           WRITE(6,*) 
     &   ' Orthogonalization term to E2 (trunc.)', DELTA**2/GAMMA
         END IF
         WRITE(6,'(A,2E25.15)') 
     &   ' Predicted energy (truncated), change and total ', 
     &                E2PREDI,EIGAPR+E2PREDI
         WRITE(6,*)
         WRITE(6,*) ' Estimated square-norm of eliminated terms ',
     &   DELNORMT-DELNORM
         WRITE(6,'(A,E25.15)') 
     &   ' Estimated energy contributions of eliminated terms',
     &   E2PREDIT-E2PREDI
C        WRITE(6,'(A,E22.15)') 
C    &   ' Number of zero elements in delta (before class sel) ',
C    &    XNZERO
C        WRITE(6,'(A,E22.15)') 
C    &   ' Number of nonzero terms in delta (before class sel) ',
C    &    XNELMNT-XNZERO
       ELSE
       END IF
       CALL LUCIAQEXIT('PARTG')
*
       IF(ICLSSEL.EQ.1) THEN
         IF(N_ACT_CLS .EQ. 0   ) THEN
           IF(IMULTGR.EQ.0.OR.IBASSPC_MX.EQ.IREFSPC) THEN
*. All classes were eliminated so we are home -and hopefully dry
           WRITE(6,*) ' No active classes  '    
           IF(MICRO_IT.EQ.1) THEN
             WRITE(6,*) ' I will therefore end the diagonalization'
             CONVER = .TRUE.
             GOTO 1001
           ELSE 
             WRITE(6,*) ' I will continue to the next macroiteration'
             GOTO 1000
           END IF
           ELSE
*. No active classes with this IBASSPC_MX, try next 
             EIG(ITER,1) =  EIG(ITER-1,1)
             ITERX = ITER
             WRITE(6,*) ' No active classes  '    
             WRITE(6,*) ' I go to the next iteration '
             GOTO 1000
           END IF
         END IF
       END IF
*
* ============================================
* 1.5 : Inverse Iteration Correction to Delta 
* ============================================
*
*
* Update delta to 
* -(H0-E0)-1(H-E)|0> + delta/gamma * (H0-E0)-1 |0>
*
* ( was +(H0-E0)-1(H-E)|0>)
* 
       CALL REWINO(LU3)
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       IOFF = 1
       IF(IOLSTM.EQ.1.AND.ABS(GAMMA).GT.1.0D-6) THEN
         WRITE(6,*) ' Inverse iteration correction will be added '
         DO IBATCH = 1, NBATCHL
           LBATCHB = ISBLBT(IBATCH)
           LBATCHE = ISLEBT(IBATCH)
*. Retrieve Batch of C
           NO_ZEROING = 0
           CALL FRMDSCN3(VEC1,LBATCHB,LBLK,LU1,NO_ZEROING,I2BLK(IOFF),
     &                   ISLBLK(IOFF))
*. Multiply with (H0-E)** -1 
           FACTOR = -EIGAPR
           ITASK = 1
           CALL DIATERM_GAS(FACTOR,ITASK,VEC1,LBATCHB,ISBLOCK,IOFF,0,
     &                       0,0)
*. Retrieve Batch of Delta
           NO_ZEROING = 0
           CALL FRMDSCN3(VEC2,LBATCHB,LBLK,LU2,NO_ZEROING,I2BLK(IOFF),
     &                   ISLBLK(IOFF))
*. And add
           FAC1 = -1.0D0
           FAC2 = DELTA/GAMMA
           CALL VECSUM(VEC2,VEC2,VEC1,FAC1,FAC2,LBATCHE)
*. Transfer to Disc
           CALL TODSCNP(VEC2,LBATCHB,ISLBLK(IOFF),LBLK,LU3)
           IOFF = IOFF + LBATCHB
         END DO
*        ^ End of loop over batches
         CALL ITODS(-1,1,LBLK,LU3)
*. Well, it is nice to have the correction vector on LU2 so 
         IREW = 1
         CALL COPVCD(LU3,LU2,VEC1,IREW,LBLK)
       END IF
*
       IF(IPRT.GE.1000) THEN
         WRITE(6,*) ' The Complete correction vector Delta '
         CALL WRTVCD_EP(VEC1,LU2,1,-1)
C             WRTVCD(SEGMNT,LU,IREW,LBLK)
       END IF
*
* ===================================
* 2 : Calculate <Delta ! H ! Delta >
* ===================================
*
       WRITE(6,*)
       WRITE(6,*) ' ==================================='
       WRITE(6,*) ' <Delta!H!Delta > under construction'
       WRITE(6,*) ' ==================================='
       WRITE(6,*)
*. Active classes
        CALL LUCIAQENTER('PARTH')
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU2,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
          CALL REPART_CIV(ISBLOCK,NBATCH,ISBLBT,ISLEBT,ISI1BT,
     &         MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ISLBLK)
          CALL REPART_CIV(ICBLOCK,NBATCH,ICBLBT,ICLEBT,ICI1BT,
     &         MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ICLBLK)
          NBATCHL = NBATCH
          NCBATCH = NBATCH
          NSBATCH = NBATCH
        ELSE
          NBATCHL = NBATCH
          NCBATCH = NBATCH
          NSBATCH = NBATCH
        END IF
* Loop over batches of H !delta>
       IOFF = 1
       DELHDEL = 0.0D0
       IRESTRICT = 1
       WRITE(6,*) ' Number of batches for Delta = ', NBATCHL
       IF(IRESTRICT.EQ.1) THEN
         EIGSHFD2 = EIGSHF/2.0D0
       ELSE
         EIGSHFD2 = EIGSHF
       END IF

       DO IBATCH = 1, NBATCHL
         IF(IPRT.GE.10) 
     &   WRITE(6,*) '  <Batch Delta ! H ! Delta >, Batch : ',IBATCH
         LBATCHB = ISBLBT(IBATCH)
         LBATCHE = ISLEBT(IBATCH)
         CALL SETVEC(VEC1,ZERO,LBATCHE)
         CALL SBLOCK(LBATCHB,ISBLOCK,IOFF,VEC2,VEC1,LU2,IRESTRICT,0, 
     &                0,0,0,0.0D0,EIGSHFD2,'SIGMA ')
*. Retrieve Batch of Delta
         CALL GET_BLOCKS_FROM_DISC
     &   (LU2,LBATCHB,IOFF,ISBLOCK,NBLOCKT,VEC2,1)
         IF(IPRT.GE.1000) THEN
           WRITE(6,*) ' Batch of H delta '
           CALL WRTBLKN_EP(VEC1,LBATCHB,ISLBLK(IOFF))
           WRITE(6,*) ' Batch of delta retrieved '
           CALL WRTBLKN_EP(VEC2,LBATCHB,ISLBLK(IOFF))
         END IF
         DELHDEL = DELHDEL + INPROD(VEC1,VEC2,LBATCHE)
         IOFF = IOFF + LBATCHB
       END DO
       CALL LUCIAQEXIT('PARTH')
       IF(IRESTRICT.EQ.1) DELHDEL = 2.0D0*DELHDEL
*
* ===========================================
* 3 : Solve 2 by 2 problem : Nonorthogonal !!
* ===========================================) 
*
*      Norm of delta and overlap between 0 and delta
       S12 = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
       S22 = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
       H11 = EIGAPR
C       H11 = ECC
       WRITE(6,*) ' GAMMA, DELTA, DELTAT = ', GAMMA, DELTA,DELTAT
       IF(IOLSTM.EQ.1.AND.ABS(GAMMA).GT.1.0D-6) THEN
        H12 = -CHEDEL + DELTA*DELTAT/GAMMA
       ELSE 
        H12 = CHEDEL + EIGAPR*DELTA
       END IF
       H22 = DELHDEL
*
       S11 = 1.0D0
*
*.( H11  H12 ) (X1)       (S11    S12 )(X1)
* (          ) (  )   = E (           )( )
* ( H12  H22 ) (X2)       (S12    S22 )(X2)
*
* The eigenvalues
*
        A = S11*S22 -S12 **2
        B = 2.0D0*S12*H12-S11*H22-S22*H11
        C = H11*H22 -H12**2
*
        EA = -B/(2.0D0*A) - SQRT(B**2 - 4.0D0*A*C)/(2.0D0*A)
        EB = -B/(2.0D0*A) + SQRT(B**2 - 4.0D0*A*C)/(2.0D0*A)
*. And the lowest eigenvalue is 
        E = MIN(EA,EB)
*. The corresponding eigenvector
*. Intermediate normalization
        X1 = 1.0D0
COLD    X2 = -(H11-E*S11)/(H12-E*S12)
*. A stable form when H11 approx E
        X2 = -(H12-E*S12)/(H22-E*S22)
*. Normalized
        XNORM2 = S11*X1**2 + S22*X2**2 + 2.0D0*S12*X1*X2
        XNORM = SQRT(XNORM2)
        X1 = X1/XNORM
        X2 = X2/XNORM
*
        IF(IPRT.GE.10) THEN
          WRITE(6,*) 
          WRITE(6,*) ' 2 X 2 Generalized eigenvalue problem, H and S '
          WRITE(6,*) 
          WRITE(6,'(4X,2X,E22.15)') H11
          WRITE(6,'(4X,2(2X,E22.15))') H12,H22
          WRITE(6,*) 
          WRITE(6,'(4X,2X,E22.15)') S11
          WRITE(6,'(4X,2(2X,E22.15))') S12,S22
          WRITE(6,*)
*
          WRITE(6,'(A,E25.15)') 
     &    ' Lowest eigenvalue (with shift)', E
          WRITE(6,*) ' Corresponding eigenvector ', X1,X2
        END IF

*. Save corresponding eigenvector on file LU1 ( first LU3, then COPY)
* VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        IREW = 1
        CALL VECSMDP(VEC1,VEC2,X1,X2,LU1,LU2,LU3,IREW,LBLK)
        CALL COPVCD(LU3,LU1,VEC1,IREW,LBLK)
        IF(IPRT.GE.1000) THEN
          WRITE(6,*)
          WRITE(6,*) ' New approximation to eigenvector '
          WRITE(6,*) ' ================================='
          WRITE(6,*)
          CALL WRTVCD_EP(VEC1,LU1,1,-1)
        END IF
*
        IF(MICRO_IT.EQ.1 ) EIG(ITER,1) = E 
        ITERX = ITER
*. Convergence of complete iteration sequence ?
      IF((IMULTGR.EQ.0.OR.IBASSPC_MX.EQ.IREFSPC).AND.MICRO_IT.EQ.1) THEN
        IF(ABS(EIG(ITER,1) - EIG(ITER-1,1)).LE.THRES_ET)
     &     CONVER = .TRUE.
      END IF
      IF(CONVER) GOTO 1001
  999 CONTINUE
*     ^ End of loop over micro-iterations 
 1000 CONTINUE
* ( End of loop over iterations )
 1001 CONTINUE
      ITER = ITERX
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINO(LU1)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD_EP(VEC1,LU1,0,LBLK)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER-1
          IF(IMULTGR.EQ.0) THEN
            WRITE(6,1340) I,EIG(I,IROOT),RNRM(I,IROOT)
          ELSE 
            IDEL = 1-IPAT(MOD(I-1,LPAT)+1)
            IF(I.EQ.1.OR.IDEL.EQ.0.OR.I.EQ.MAXIT-1) THEN
              WRITE(6,1341) I,EIG(I,IROOT),RNRM(I,IROOT)
 1341         FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5,'  Full resid.')
            ELSE
              WRITE(6,1342) I,EIG(I,IROOT),RNRM(I,IROOT)
 1342         FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5,'  Partial resid.')
            END IF
          END IF
 1330     CONTINUE
          WRITE(6,1340) ITER,EIG(ITER,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)
 1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER-1,IROOT)
 1607   CONTINUE
      END IF
C
*. Clean up : ICBLOCK, ISBLOCK in original form
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      CALL REPART_CIV(ICBLOCK,NBATCHL,ICBLBT,ICLEBT,ICI1BT,
     &     MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ICLBLK)
      CALL REPART_CIV(ISBLOCK,NBATCHL,ICBLBT,ICLEBT,ICI1BT,
     &     MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,ICLBLK)
*
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE RSBB2BN2(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                   JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                   IAGRP,IBGRP,IOCTPA,IOCTPB, 
     &                   NGAS,IAOC,IBOC,JAOC,JBOC,
     &                   SB,CB,ADSXA,STSTSX,MXPNGASX,
     &                   NOBPTS,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                   XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBSX,IUSEAB,
     &                   CJRES,SIRES,SCLFAC,NTESTG,
     &                   NSEL2E,ISEL2E,IUSE_PH,IPHGAS,XINT2,NSTFSMSPGP,
     &                   IASPGP,IBSPGP,JASPGP,JBSPGP)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
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
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT  : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + CJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks 
* October 96 : New routines for accessing annihilation information
*             Cleaned and shaved, only IROUTE = 3 option active
* October   97 : allowing for N-1/N+1 switch
*
* Last change : Aug 2000
*. Some dull optimazation, July 2003
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
      REAL*8 INPROD
      INTEGER NSTFSMSPGP(MXPNSMST,*)
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*), XINT2(*)
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION H(MXPTSOB*MXPTSOB)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
      INTEGER NKA_AS(MXPNSMST), NKB_AS(MXPNSMST)
*. Arrays for reorganization 
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C    &          SIGNREO,NADOP,NADDEL,IADDEL,ADSIGN)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'intform.inc'
*
      COMMON/XXTEST/ISETVECOPS(10)
*
      CALL LUCIAQENTER('RS2B2')
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.500.and. ga_nodeid().eq.0) THEN
*
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' RSBB2BN2 speaking '
        WRITE(6,*) ' ================ '
*
        WRITE(6,*) ' Occupation of IA '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of IB '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JA '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JB '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)

*
        WRITE(6,*) ' Memcheck at start of RSBB2BN '
        CALL MEMCHK
        WRITE(6,*) ' Memory check passed '
*
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Groups defining each supergroup
C     CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
C     CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
C     CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
C     CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600.and.ga_nodeid().eq.0) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)           
*
      ITEST = 0
      IF(ITEST.EQ.1.AND.(NIJTYP.EQ.0.AND.NKLTYP.EQ.0)) THEN
        WRITE(6,*) ' Unneccesary entrance to RSBB2BN2'
        WRITE(6,*) ' IAOC, IBOC,JAOC,JBOC : '
        WRITE(6,*) ' Occupation of IA '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of IB '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JA '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JB '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' IATP, JATP, IBTP, JBTP = ',
     &               IATP, JATP, IBTP, JBTP
        STOP 
      END IF
*
      IAFRST = 1
      IBFRST = 1
      JAFRST = 1
      JBFRST = 1
*
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
c     nproc = ga_nnodes()
      DO 2001 IJTYP = 1, NIJTYP
*
        IJFIRST = 1
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Should N-1 or N+1 projection be used for alpha strings
          IJ_TYP(1) = ITYP
          IJ_TYP(2) = JTYP
          IJ_AC(1)  = 2
          IJ_AC(2) =  1
          NOP = 2
          IF(IUSE_PH.EQ.1) THEN
            CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,
     &           SIGNIJ)
          ELSE
*. Enforced a+ a
            IJ_REO(1) = 1
            IJ_REO(2) = 2
            SIGNIJ = 1.0D0
          END IF
*. Two choices here :
*  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> ( good old creation mapping)
*  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka>  + delta(i,j)                   
C?        WRITE(6,*) ' RSBB2BN : IOP_REO : ', (IOP_REO(II),II=1,2)
          IF(IJ_REO(1).EQ.1.AND.IJ_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
            IJAC = 2
            IOP2AC = 1
            SIGNIJ2 = SCLFAC
*
            IJ_DIM(1) = NI
            IJ_DIM(2) = NJ
            IJ_SYM(1) = ISM
            IJ_SYM(2) = JSM
            IJ_TYP(1) = ITYP
            IJ_TYP(2) = JTYP
*
            NOP1   = NI
            IOP1SM = ISM
            IOP1TP = ITYP
            NOP2   = NJ
            IOP2SM = JSM
            IOP2TP = JTYP
          ELSE
*. Terra Nova, annihilation map 
            IJAC = 1
            IOP2AC = 2
            SIGNIJ2 = -SCLFAC
*
            IJ_DIM(1) = NJ
            IJ_DIM(2) = NI
            IJ_SYM(1) = JSM
            IJ_SYM(2) = ISM
            IJ_TYP(1) = JTYP
            IJ_TYP(2) = ITYP
*
            NOP1   = NJ
            IOP1SM = JSM
            IOP1TP = JTYP
            NOP2   = NI
            IOP2SM = ISM
            IOP2TP = ITYP
          END IF
*
          IF(IJFIRST.EQ.1) THEN
*. Find supergroup type of Kstring
C  NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
            CALL NEWTYP(JATP+IOCTPA-1,IOP2AC,IOP2TP,1,KATP_ABS)
            CALL ICOPVE(NSTFSMSPGP(1,KATP_ABS),NKA_AS,NSMST)
            IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' NKA_AS as delivered '
            CALL IWRTMA(NKA_AS,1,NSMST,1,NSMST)
            END IF
            IJFIRST = 0
          END IF
*
*. Generate creation- or annihilation- mappings for all Ka strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          KASM = MULTD2H(JASM,IJ_SYM(2) )
          NKASTR = NKA_AS(KASM)
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' NKASTR from NKA_AS = ', NKASTR
          CALL ADAST2_GAS(IJ_SYM(2),IJ_TYP(2),NGAS,JASPGP,JASM,
     &         I1,XI1S,NKASTR,IEND,JAFRST,KFRST,KACT,SIGNIJ2,IJAC,
     &         JAACT,1)
          IF(JAACT.EQ.1) JAFRST = 0
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          CALL ADAST2_GAS(IJ_SYM(1),IJ_TYP(1),NGAS,IASPGP,IASM,
     &         I3,XI3S,NKASTR,IEND,IAFRST,KFRST,KACT,ONE,IJAC,
     &         IAACT,2)
          IF(IAACT.EQ.1) IAFRST = 0
*. Compress list to common nonvanishing elements
          IDOCOMP = 0
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),
     &                        NKASTR,NKAEFF)
          ELSE 
              NKAEFF = NKASTR
          END IF
            
*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
c         if (nkabtc.lt.nproc) then
c            nkabtc = min(nkaeff,nproc)
c         endif
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
              KLFIRST = 1
*. Allowed double excitation ?
              IJKL_ACT = I_DX_ACT(ITYP,KTYP,LTYP,JTYP)
              IF(IJKL_ACT.EQ.0) GOTO 2000
              IF(NTEST.GE.500) THEN
                WRITE(6,*) ' KTYP, LTYP', KTYP, LTYP 
              END IF
*. Should this group of excitations be included 
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              KL_TYP(1) = KTYP
              KL_TYP(2) = LTYP
              KL_AC(1)  = 2
              KL_AC(2) =  1
              NOP = 2
              IF(IUSE_PH.EQ.1) THEN
                CALL ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,
     &               SIGNKL)
              ELSE
*. Enforced a+ a
                KL_REO(1) = 1
                KL_REO(2) = 2
                SIGNKL = 1.0D0
              END IF
*
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = ADSXA(KSM,KLSM)
                IF(NTEST.GE.500) THEN
                  WRITE(6,*) ' KSM, LSM', KSM, LSM
                END IF
                IF(LSM.EQ.0) GOTO 1930
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*
                IF(KL_REO(1).EQ.1.AND.KL_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
                  IOP4AC = 1
                  KLAC = 2
                  KL_DIM(1) = NK
                  KL_DIM(2) = NL
                  KL_SYM(1) = KSM
                  KL_SYM(2) = LSM
                  KL_TYP(1) = KTYP
                  KL_TYP(2) = LTYP
                ELSE
*. Terra Nova, annihilation map 
                  IOP4AC = 2
                  KLAC = 1
                  KL_DIM(1) = NL
                  KL_DIM(2) = NK
                  KL_SYM(1) = LSM
                  KL_SYM(2) = KSM
                  KL_TYP(1) = LTYP
                  KL_TYP(2) = KTYP
                END IF
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0  
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
                IF(KLFIRST.EQ.1) THEN
*. Find supergroup type of Kstring
C  NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
                  CALL NEWTYP(JBTP+IOCTPB-1,IOP4AC,KL_TYP(2),1,KBTP_ABS)
                  CALL ICOPVE(NSTFSMSPGP(1,KBTP_ABS),NKB_AS,NSMST)
                  KLFIRST = 0
                END IF
                KBSM = MULTD2H(JBSM,KL_SYM(2) )
                NKBSTR = NKB_AS(KBSM)
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
*. currently we are using creation mappings for kl
*. (Modify to use ADAST later )
                CALL ADAST2_GAS(KL_SYM(2),KL_TYP(2),NGAS,JBSPGP,JBSM,
     &               I2,XI2S,NKBSTR,IEND,JBFRST,KFRST,KACT,SIGNKL,KLAC,
     /               JBACT,3)
                IF(JBACT.EQ.1) JBFRST = 0
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADAST2_GAS(KL_SYM(1),KL_TYP(1),NGAS,IBSPGP,IBSM,
     &               I4,XI4S,NKBSTR,IEND,IBFRST,KFRST,KACT,ONE,KLAC,
     &               IBACT,4)
                IF(IBACT.EQ.1) IBFRST = 0
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (iop2 iop1 |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                ONE = 1.0D0
                IF(IH2FORM.EQ.1) THEN
*. Normal integrals with conjugation symmetry
                  CALL LUCIAGETINT(XINT,IJ_TYP(2),IJ_SYM(2),
     &                 IJ_TYP(1),IJ_SYM(1),
     &                 KL_TYP(1),KL_SYM(1),KL_TYP(2),KL_SYM(2),IXCHNG,
     &                 0,0,ICOUL,ONE,ONE)
                ELSE IF (IH2FORM.EQ.2) THEN
C?              WRITE(6,*) ' I_USE_SIMTRH = ', I_USE_SIMTRH
*. Integrals does not have conjugation symmetry so be careful...
*. The following is not enough is particle hole symmetry is encountered
*. Obtain ( i j ! k l )
                  CALL LUCIAGETINT(XINT,ITYP,ISM,JTYP,JSM,
     &                             KTYP,KSM,LTYP,LSM,
     &                        IXCHNG,0,0,ICOUL,ONE,ONE)
                  IF(KLAC.EQ.2.AND.IJAC.EQ.2) THEN
*. Transpose to obtain ( j i ! k l )
                    CALL TRP_H2_BLK(XINT,12,NI,NJ,NK,NL,XINT2)
                  ELSE IF(KLAC.EQ.1.AND.IJAC.EQ.2) THEN  
*. Transpose to obtain (j i | l k)
                    CALL TRP_H2_BLK(XINT,46,NI,NJ,NK,NL,XINT2)
                  ELSE IF (KLAC.EQ.1.AND. IJAC .EQ. 1 ) THEN
*. Transpose to obtai (i j | l k)
                    CALL TRP_H2_BLK(XINT,34,NI,NJ,NK,NL,XINT2)
                  END IF
                END IF
*
* S(Ka,i,Ib) = sum(j,k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
C               IJKL_DIM = IJ_DIM(1)*IJ_DIM(2)*KL_DIM(1)*KL_DIM(2)
C               IF(INPROD(XINT,XINT,IJKL_DIM).NE.0.0D0) THEN
c         DO 1801 IKABTC = ga_nodeid()+1, NKABTC, nproc
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
CJTEST      WRITE(6,*) ' JTEST: Dimension of CJRES and SIRES ',
C    &      IJ_DIM(2)*LKABTC*NJB, IJ_DIM(1)*LKABTC*NIB
*. Obtain C(ka,J,JB) for Ka in batch
*
            FACS = 1.0D0
            DO JJ = 1, IJ_DIM(2)
              CALL GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
C             GET_CKAJJB(CB,NJ,NJA,CKAJJB,NKA,NJB,J,ISCA,SSCA)

            END DO
            IF(NTEST.GE.2500) THEN
              WRITE(6,*) ' Updated CJRES as C(Kaj,Jb)'
              CALL WRTMAT(CJRES,NKASTR*NJ,NJB,NKASTR*NJ,NJB)
            END IF
*
            ISETVECOPS(2) = ISETVECOPS(2) + NIB*LKABTC*IJ_DIM(1)
            MXACJ=MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC*IJ_DIM(1))
*
            IROUTE = 3
            FACS = 1.0D0
            CALL SKICKJ(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),
     &               KL_DIM(1),KL_DIM(2),
     &               NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE )
C               END IF
*
            IF(NTEST.GE.2500) THEN
               WRITE(6,*) ' Updated Sires as S(Kai,Ib)'
               CALL WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
            END IF
*
*. Scatter out from s(Ka,Ib,i)
*
            DO II = 1, IJ_DIM(1)
              CALL ADD_SKAIIB(SB,IJ_DIM(1),NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
*
 1930     CONTINUE
*         ^ End of loop over KSM
 2000     CONTINUE
*         ^ End of loop over KLTYP
*. Was anything done ?
 1940   CONTINUE
*       ^ End of loop over ISM
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL LUCIAQEXIT('RS2B2')
      RETURN
      END
      SUBROUTINE ADAST2_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
     &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC,
     &                    IAC,IACT,IINDEX)
*
*
* Obtain creation or annihilation mapping
*
* IAC = 2 : Creation map
* a+IORB !KSTR> = +/-!ISTR> 
*
* IAC = 1 : Annihilation map
* a IORB !KSTR> = +/-!ISTR> 
*
* for orbitals of symmetry IOBSM and type IOBTP
* and Istrings defined by the NIGRP groups IGRP and symmetry ISPGPSM
* 
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR> 
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
* If info for the Istrings has been set up, IACT = 1
* IINDEX tells in which IOFFI info on Istrings should go
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version 
*               October 96   : Improved version
*               September 97 : annihilation mappings added
*                              I groups defined by IGRP
*               July 2003,   : Some dull optimization ...
*                              Notice that NKSTR is now input ...
*               May 2013     : Changing induced symmetry index to 
*                              first
*
*. Last modification; May5, 2013; Jeppe Olsen; Changed induced symmetry
*                     index to first
*         
*
*
* ======
*. Input
* ======
*
*./BIGGY
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER ISMFGS(MXPNGAS)
C     INTEGER MXVLI(MXPNGAS),MNVLI(MXPNGAS)
      INTEGER MXVLK(MXPNGAS),MNVLK(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      INTEGER KGRP(MXPNGAS)
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8)
      COMMON/LOC_ADAST2/IOFFI(LOFFI,4),MXVLI(MXPNGAS,4),MNVLI(MXPNGAS,4)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'multd2h.inc'
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension 
* (NKSTR,*), Where NKSTR is the number of K-strings of 
*  correct symmetry . Nk is provided by this routine.
*
C!    CALL LUCIAQENTER('ADAST ')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADAST_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        WRITE(6,*) '  Symmetry of supergroup : ', ISPGPSM
        WRITE(6,*) ' SCLFAC = ', SCLFAC
*
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' Annihilation mapping '
        ELSE IF(IAC.EQ.2) THEN
          WRITE(6,*) ' Creation mapping '
        ELSE 
          WRITE(6,*) ' Unknown IAC parameter in ADAST ',IAC
          STOP       ' Unknown IAC parameter in ADAST '
        END IF
*
      END IF
*. A few preparations
      NORBTS= NOBPTS(IOBTP,IOBSM)
      NORBT= NOBPT(IOBTP)
      IACGAS = IOBTP
*
      IACT = 0
*. First orbital of given GASpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1 + NINOB
*. First orbital of given GASpace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM) 
*
*====================================================
*. K strings : Supergroup, symmetry and distributions
*====================================================
*
      IF(IAC.EQ.1) THEN
       IDELTA = +1
      ELSE
       IDELTA = -1
      END IF
*. Is required mapping contained within current set of maps?
*. a:) Is active GASpace included in IGRP - must be 
      IACGRP = 0
      DO JGRP = 1, NIGRP
       IF(IGSFGP(IGRP(JGRP)).EQ. IACGAS) IACGRP = JGRP
      END DO
*. Note : IACGRP is not the actual active group, it is the address of the
*         active group in IGRP
      IF(IACGRP.EQ.0) THEN
        WRITE(6,*) ' ADAST in problems '
        WRITE(6,*) ' Active GASpace not included in IGRP '
        WRITE(6,*) ' Active GASpace : ', IACGAS
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        STOP       ' ADAST : Active GASpace not included in IGRP '
      END IF
*. b:) active group in K strings
      NIEL = NELFGP(IGRP(IACGRP))
      NKEL = NIEL + IDELTA
      IF(NTEST.GE.1000) WRITE(6,*) ' NIEL and NKEL ',NIEL,NKEL
      IF(NKEL.EQ.-1.OR.NKEL.EQ.NOBPT(IACGAS)+1) THEN
*. No strings with this number of elecs - be happy : No work 
        NKSTR = 0
        KACT = 0
        KACGRP = 0
        GOTO 9999
      ELSE
*. Find group with NKEL electrons in IACGAS
        KACGRP = 0
        DO JGRP = IBGPSTR(IACGAS),IBGPSTR(IACGAS)+NGPSTR(IACGAS)-1
          IF(NELFGP(JGRP).EQ.NKEL) KACGRP = JGRP
        END DO
        IF(NTEST.GE.1000) WRITE(6,*) ' KACGRP = ',KACGRP
*. KACGRP is the Active group itself     
        IF(KACGRP.EQ.0) THEN
          WRITE(6,*)' ADAST : cul de sac, active K group not found'
          WRITE(6,*)' GAS space and number of electrons ',
     &               IACGAS,NKEL
          STOP      ' ADAST : cul de sac, active K group not found'
        END IF
      END IF
*. Okay active K group was found and is nontrivial
      KSM = MULTD2H(IOBSM,ISPGPSM)
*. The K supergroup
      CALL ICOPVE(IGRP,KGRP,NIGRP)
      KGRP(IACGRP) = KACGRP
*. Number of strings and symmetry distributions of K strings
      IF(NTEST.GE.1000) WRITE(6,*) 
     & ' KSM, NKSTR : ', KSM, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*. Last active space in K strings and number of strings per group and sym
      NGASL = 1
      DO JGRP = 1, NIGRP
       IF(NELFGP(KGRP(JGRP)).GT.0) NGASL = JGRP
      END DO
*. MIN/MAX for Kstrings
      DO JGRP = 1, NIGRP
        MNVLK(JGRP) =  MINMAX_SM_GP(1,KGRP(JGRP))
        MXVLK(JGRP) =  MINMAX_SM_GP(2,KGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) THEN
        write(6,*) 'MNVLK and MXVLK '
        CALL IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
      END IF
*. (NKDIST_TOT is number of distributions, all symmetries )
* ==============
*. I Strings 
* ==============
*. Generate symmetry distributions of I strings with given symmetry
      NIGASL = 1
      DO JGRP = 1, NIGRP
        IF(NELFGP(IGRP(JGRP)).GT.0) NIGASL = JGRP
      END DO
*
      IF(IFRST.EQ.1) THEN
        DO IGAS = 1, NIGASL
          MNVLI(IGAS,IINDEX) = MINMAX_SM_GP(1,IGRP(IGAS))
          MXVLI(IGAS,IINDEX) = MINMAX_SM_GP(2,IGRP(IGAS))
        END DO
        CALL TS_SYM_PNT2(IGRP,NIGASL,
     &       MXVLI(1,IINDEX),MNVLI(1,IINDEX),ISPGPSM,
     &       IOFFI(1,IINDEX),LOFFI)
      END IF
      IACT = 1
*. Number of I strings per group and sym
*. Last entry in IGRP with a nonvanishing number of strings
*. Number of electrons before active space
      NELB = 0
      DO JGRP = 1, IACGRP-1
        NELB = NELB + NELFGP(IGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NELB = ', NELB
*
      ZERO =0.0D0
      IZERO = 0    
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
      DO IGAS = 1, NIGRP
        ISMFGS(IGAS) = 1
      END DO
 1000 CONTINUE
*. Next distribution
        CALL NEXT_SYM_DISTR(NGASL,MNVLK,MXVLK,ISMFGS,KSM,KFIRST,NONEW)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' Symmetry distribution ' 
          call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
        END IF
        IF(NONEW.EQ.1) GOTO 9999
        KFIRST = 0
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
C         NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
          NSTRIK = NSTRIK*NSTFSMGP(ISMFGS(IGAS),KGRP(IGAS))
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IACGRP)
        IACSM = MULTD2H(IOBSM,ISMFGS(IACGRP))
        ISMFGS(IACGRP) = IACSM
        IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,IOFFI(1,IINDEX),
     &             MXVLI(1,IINDEX),MNVLI(1,IINDEX))
        ISMFGS(IACGRP) = ISAVE
*. Number of strings before active GAS space
        NSTB = 1
        DO IGAS = 1, IACGRP-1
          NSTB = NSTB*NSTFSMGP(ISMFGS(IGAS),KGRP(IGAS))
        END DO
*. Number of strings After active GAS space
        NSTA = 1
        DO IGAS =  IACGRP+1, NIGRP
          NSTA = NSTA*NSTFSMGP(ISMFGS(IGAS),KGRP(IGAS))
        END DO
*. Number and offset for active group 
        NIAC  = NSTFSMGP(IACSM,IGRP(IACGRP))
        IIAC  = ISTFSMGP(IACSM,IGRP(IACGRP))
C       NIAC  = NACIST(IACSM)
C       IIAC =  IACIST(IACSM)
*
C       NKAC = NNSTSGP(ISMFGS(IACGRP),IACGRP)
C       IKAC = IISTSGP(ISMFGS(IACGRP),IACGRP)
        NKAC = NSTFSMGP(ISMFGS(IACGRP),KGRP(IACGRP))
        IKAC = ISTFSMGP(ISMFGS(IACGRP),KGRP(IACGRP))
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
        IF(NTEST.GE.1000) THEN
        write(6,*) ' nstb nsta niac nkac ',
     &               nstb,nsta,niac,nkac
        END IF
*. Obtain annihilation/creation mapping for all strings of this type
*. Are group mappings in expanded or compact form 
        IF(IAC.EQ.1.AND.ISTAC(KACGRP,2).EQ.0) THEN
          IEC = 2
          LROW_IN = NKEL
        ELSE 
          IEC = 1
          LROW_IN = NORBT
        END IF
        NKACT = NSTFGP(KACGRP)
*
        MXAADST = MAX(MXAADST,NKSTR*NORBTS)
        IF(NSTA*NSTB*NIAC*NKAC.NE.0)
     &  CALL ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,   
     &            int_mb(KSTSTM(KACGRP,1)),int_mb(KSTSTM(KACGRP,2)),
     &            IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &            NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC,IAC,
     &            LROW_IN,IEC)
        KSTRBS = KSTRBS + NKSD     
        GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
*     JTEST = 1
      IF(JTEST.EQ.1) THEN
*. Chasing a bug, Jan 2008
        NNEG = 0
        DO IORB = 1, NORBTS
        DO KSTR = 1, NKSTR
          IF(I1((IORB-1)*NKSTR + KSTR).LT.0) THEN
           WRITE(6,*) ' Problem in ADAST2, negative adress'
           WRITE(6,*) ' IORB, KSTR, I1(KSTR,IORB) = ',
     &                  IORB, KSTR, I1((IORB-1)*NKSTR + KSTR)
          END IF
        END DO
        END DO
*
        IF(NNEG.NE.0) THEN
          CALL MEMCHK2('ADAST')
          STOP ' Problem in ADAST2' 
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADAST_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMAT(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
C!    CALL LUCIAQEXIT('ADAST ')
      RETURN
      END
      SUBROUTINE TEST_RHO2S(RHO2,RHO2AA,RHO2AB,RHO2BB,NORB)
*
* Test two-body spin-density matrices by assembling 
* the standard two-body matrix from these and compare
*
* Jeppe Olsen, sitting in Oak Ridge, Sept. 2004
*
      INCLUDE 'implicit.inc'
*
      DIMENSION RHO2(*),RHO2AA(*),RHO2AB(*), RHO2BB(*)
*
*. Remember the form of the various densities : 
*
* RHO2(ij,kl) = sum_(s,s') <L!a+is a+ks' als' aj s!R> (ij.ge.kl)
* RHO2AA(iklj) = <L!a+iaa+ka aja ala!R> (i.ge.k, j.ge.l)
* RHO2BB(iklj) = <L!a+iba+kb ajb alb!R> (i.ge.k, j.ge.l)
* RHO2AB(iklj) = <L!a+iaa+kb alb aja!R> ( no restrictions)
*
      NTEST = 0000
      WRITE(6,*) ' Welcome to test of 2e-spindensities '
      WRITE(6,*) ' =================================== '
*
      NERROR = 0
      THRES = 1.0D-10
      DO J = 1, NORB
        DO L = 1, J
          DO I = 1, NORB
            IF(J.EQ.L) THEN
              MAXK = I
            ELSE
              MAXK = NORB
            END IF
            DO K = 1, MAXK
*. address in rho2
              IJ = (J-1)*NORB + I
              KL = (L-1)*NORB + K
              IJKL_RHO2 = IJ*(IJ-1)/2 + KL
*. addresses in rho2ab
              IKLJ_RHO2AB = (J-1)*NORB**3 + (L-1)*NORB**2
     &                    + (K-1)*NORB    +  I
              KIJL_RHO2AB = (L-1)*NORB**3 + (J-1)*NORB**2
     &                    + (I-1)*NORB    + K
*. adress in rho2ss
              IF(I.GT.K) THEN
                IK = I*(I-1)/2 + K
                SIGN = -1.0D0
              ELSE
                IK = K*(K-1)/2 + I
                SIGN = 1.0D0
              END IF
C             LJ = L*(L-1)/2 + J
              LJ = J*(J-1)/2 + L
              IKLJ_RHO2SS = (LJ-1)*NORB*(NORB+1)/2 + IK
*
              RHO2_1 = RHO2(IJKL_RHO2)
              RHO2_2 = RHO2AB(IKLJ_RHO2AB)
     &               + RHO2AB(KIJL_RHO2AB)
     &               + SIGN*RHO2AA(IKLJ_RHO2SS)
     &               + SIGN*RHO2BB(IKLJ_RHO2SS)
              IF(ABS(RHO2_1-RHO2_2).GT.THRES) THEN
                NERROR = NERROR + 1
                WRITE(6,*) ' Problem with spinden '
                WRITE(6,*) ' I,J,K,L= ', I,J,K,L
                WRITE(6,*) ' element from RHO2 and from RHO2S ',
     &          RHO2_1,RHO2_2
                WRITE(6,*) ' term from RHO2AB =',
     &          RHO2AB(IKLJ_RHO2AB)+ RHO2AB(KIJL_RHO2AB)
                WRITE(6,*) ' term from RHO2SS ',
     &          SIGN*(RHO2AA(IKLJ_RHO2SS)+RHO2BB(IKLJ_RHO2SS))
                WRITE(6,*) ' IKLJ_RHO2SS = ', IKLJ_RHO2SS
              END IF
            END DO
          END DO
        END DO
      END DO
*
      IF(NERROR.EQ.0) THEN
        WRITE(6,*) ' RHO2 and RHO2S are in agreement '
      ELSE
        WRITE(6,*) ' Number of differences between RHO2 and RHO2S',
     &               NERROR 
      END IF
*
      RETURN
      END
      SUBROUTINE COMHAM(H,NVAR,NBLOCK,LBLOCK,VEC1,VEC2,ECORE)
*
* Construct complete H matrix through a sequence 
* of direct CI iterations 
*
* BLocks assumed to be allocated outside
*
* Jeppe Olsen, April 2003
* Last revision, Oct. 1 2012, Jeppe Olsen, A few cosmetic changes
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
*. Input
      INTEGER LBLOCK(NBLOCK)
*. Output
      DIMENSION H(NVAR,NVAR)
*. Dirty dancing.. no checks
      PARAMETER (MXPVAR = 400)
      DIMENSION SCR(MXPVAR*(MXPVAR+1)/2)
      
*. Scratch through argument list, should be able to hold largest SD TTS block
      DIMENSION VEC1(*), VEC2(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'COMHAM')
*.
      CALL MEMMAN(KLVEC1,NVAR,'ADDL  ',2,'VEC1  ') !done
      CALL MEMMAN(KLVEC2,NVAR,'ADDL  ',2,'VEC2  ') !done
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from COMHAM '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' Number of parameters = ', NVAR
        WRITE(6,*) ' Number of blocks = ', NBLOCK
        WRITE(6,*) ' Length of each block : '
        CALL IWRTMA(LBLOCK,1,NBLOCK,1,NBLOCK)
      END IF
*
       DO I = 1, NVAR
       WRITE(6,*) ' Generating column = ', I
*. Create i'th unit vector
       ZERO = 0.0D0
       CALL SETVEC(dbl_mb(KLVEC1),ZERO,NVAR)
       dbl_mb(KLVEC1-1+I) = 1.0D0
*. Transfer to disc in TTS blocked form 
       IF(ICISTR.GT.1) THEN
         CALL REWINO(LUSC1)
         CALL TODSCN(dbl_mb(KLVEC1),NBLOCK,LBLOCK,-1,LUSC1)
C             TODSCN(VEC,NREC,LREC,LBLK,LU)
         CALL ITODS(-1,1,-1,LUSC1)
       ELSE
         CALL COPVEC(dbl_mb(KLVEC1),VEC1,NVAR)
       END IF
*. H * LUSC1 => LUHC
       CALL MV7(VEC1,VEC2,LUSC1,LUHC,0,0)
C!     STOP ' Jeppe forced me to stop after MV7'
*. Read in He_i and save 
       IF(ICISTR.GT.1) THEN
         CALL REWINO(LUHC)
         CALL FRMDSCN(H(1,I),NBLOCK,-1,LUHC)
C             FRMDSCN(VEC,NREC,LBLK,LU)
       ELSE
          CALL COPVEC(VEC2,H(1,I),NVAR)
       END IF
*. Add core energy to Hamiltonian
       H(I,I) = H(I,I) + ECORE
       IF(NTEST.GE.500) THEN
         WRITE(6,*) ' Output sigma-vector '
         CALL WRTMAT(H(1,I),1,NVAR,1,NVAR)
       END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The complete H matrix '
        CALL WRTMAT_EP(H,NVAR,NVAR,NVAR,NVAR)
      END IF
*
      I_DO_DIAG = 0
      IF(I_DO_DIAG.EQ.1) THEN
*. Diagonalize the Hamiltonian matrix: First pack and then diag
*. Reform to packed matrix
        CALL TRIPAK(H,SCR,1,NVAR,NVAR)
*. Diagonalize
        CALL EIGENL(SCR,H,NVAR,0,1)
*. Pack eigenvalues
        CALL COPDIA(SCR,SCR,NVAR,1)
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMAT_EP(SCR,NVAR,1,NVAR,1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ', IDUM,'COMHAM')
      RETURN
      END 
      FUNCTION IS_ACCOCC_IN_ACCOCC(IACCOCC,IACCOCC_TOT,NGAS,
     &                              NDIM_IACCOCC)
*
* Occupation of an occupation class is given in IACCOCC
* Is this occupation included in accumulated occupation IACCOCC_TOT ?
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IACCOCC(NGAS),IACCOCC_TOT(NDIM_IACCOCC,2)
*
      INCLUDED = 1
      DO IGAS = 1, NGAS
        IF(IACCOCC(IGAS).LT.IACCOCC_TOT(IGAS,1).OR.
     &     IACCOCC(IGAS).GT.IACCOCC_TOT(IGAS,2)    ) INCLUDED = 0
      END DO
*
      IS_ACCOCC_IN_ACCOCC = INCLUDED 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  NDIM_IACCOCC = ',  NDIM_IACCOCC
        WRITE(6,*) ' Occupation to be tested : '
        CALL IWRTMA(IACCOCC,NGAS,1,NGAS,1)
        WRITE(6,*) ' Min and max of accumulated occupation '
        CALL IWRTMA(IACCOCC_TOT,NGAS,2,NDIM_IACCOCC,2)
        IF(INCLUDED.EQ.1) THEN
           WRITE(6,*) ' Occupation class is included '
        ELSE 
           WRITE(6,*) ' Occupation class is not included '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE IAIB_TO_ACCOCC(IATP,IAGRP,IBTP,IBGRP,IACCOCC)
*
* Numbers of IA and IB to accumulated occupation 
*
* Jeppe Olsen, April 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'strbas.inc'
*. Output
      INTEGER IACCOCC(NGAS)
C?    WRITE(6,*) ' IATP, IBTP, IAGRP, IBGRP = ',
C?   &             IATP, IBTP, IAGRP, IBGRP 
      IATP_ABS = IATP + IBSPGPFTP(IAGRP) - 1
      IBTP_ABS = IBTP + IBSPGPFTP(IBGRP) - 1
C?    WRITE(6,*) ' IATP_ABS, IBTP_ABS ', IATP_ABS, IBTP_ABS
*
*. Occupation in each GAS space 
*
      IONE = 1
      CALL IVCSUM(IACCOCC,NELFSPGP(1,IATP_ABS),NELFSPGP(1,IBTP_ABS),
     &            IONE,IONE,NGAS)
C?    WRITE(6,*) ' Occupation before accumulation '
C?    CALL IWRTMA(IACCOCC,1,NGAS,1,NGAS)
*
*. And accumulated occupation 
*
      DO IGAS = 2, NGAS
        IACCOCC(IGAS) = IACCOCC(IGAS) + IACCOCC(IGAS-1)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Accumulated occupation from IAIB_TO_ACCOCC '
        CALL IWRTMA(IACCOCC,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE CMP_INFBLK_OCCSPC(INFBLK,NBLOCK,IOCCSPC,IBLK_IN_OCC,
     &                             NGAS,NVAR_IN_OCC,IFLAG)
*
* A CI space is defined by INFBLK. Check whether the blocks are
* included in accumulated occupation space and report back in
* IBLK_IN_OCC
*
* IF IFLAG = 1, only the number of parameters in space is obtained 
*
* Jeppe Olsen, April 2003
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER INFBLK(8,NBLOCK), IOCCSPC(MXPNGAS,2)
*. Output 
      INTEGER IBLK_IN_OCC(NBLOCK)
*. Local scratch 
      INTEGER IOCCL(MXPNGAS)
*
      NVAR_IN_OCC = 0
      IAGRP = 1
      IBGRP = 2
      DO IBLOCK = 1, NBLOCK
*. Accumulated occupation 
        JATP = INFBLK(1,IBLOCK)
        JBTP = INFBLK(2,IBLOCK)
C?      WRITE(6,*) ' IBLOCK, JATP, JBTP = ', IBLOCK,JATP,JBTP
        CALL IAIB_TO_ACCOCC(JATP,IAGRP,JBTP,IBGRP,IOCCL)
C            IAIB_TO_ACCOCC(IAGRP,IATP,IBGRP,IBTP,IACCOCC)
*. Is accumulated occupation included ? 
        INCLUDED = IS_ACCOCC_IN_ACCOCC(IOCCL,IOCCSPC,NGAS,MXPNGAS)
C                  IS_ACCOCC_IN_ACCOCC(IACCOCC,IACCOCC_TOT,NGAS,
C    &                                              MXPNGAS)
        IF(INCLUDED.EQ.1) NVAR_IN_OCC = NVAR_IN_OCC + INFBLK(IBLOCK,8)
      
        IF(IFLAG.NE.1) IBLK_IN_OCC(IBLOCK) = INCLUDED 
      END DO
*
      NTEST = 00
      IF ( NTEST .GE. 100) THEN
        WRITE(6,*) ' Number of parameters in occspace ', NVAR_IN_OCC
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' IBLK_IN_OCC array '
          CALL IWRTMA(IBLK_IN_OCC,1,NBLOCK,1,NBLOCK)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE SDNUM_FOR_SELBLKS(INFBLK,NBLOCK,ISELBLK,ISELSD,NSELSD)
*
* A CI expansions is defined by INFBLK. Obtain the numbers of 
* the SD's in the blocks selected by nonvanishing numbers in ISELBLK
*
* Jeppe Olsen, April 2003
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER INFBLK(8,NBLOCK), ISELBLK(NBLOCK)
*. Output 
      INTEGER ISELSD(*)
*
      IB=1
      IB_SEL = 1
      DO IBLOCK = 1, NBLOCK
        L = INFBLK(8,IBLOCK)
        IF( ISELBLK(IBLOCK).EQ.1) THEN
          DO I = 1, L
            ISELSD(IB_SEL-1+I) = IB-1+I
          END DO
          IB_SEL = IB_SEL + L
        END IF
        IB = IB + L
      END DO
      NSELSD = IB_SEL - 1
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Addresses of SDs in selected blocks : '
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        CALL IWRTMA(ISELSD,1,NSELSD,1,NSELSD)
      END IF
*
      RETURN
      END
      SUBROUTINE GATMAT(XMATO,XMATI,IGAT,NDIMO,NDIMI)
*
* MATO(I,J) = MATI(IREO(I),IREO(J))
*
      INCLUDE 'implicit.inc'
*. Input matrix
      DIMENSION XMATI(NDIMI,NDIMI)
*. Output
      DIMENSION XMATO(NDIMO,NDIMO)
*. Gather array
      INTEGER IGAT(NDIMO)
*
      DO J = 1, NDIMO
        JREO = IGAT(J)
        DO I = 1, NDIMO
          XMATO(I,J) = XMATI(IGAT(I),JREO)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GATMAT : input and output matrices '
        CALL WRTMAT(XMATI,NDIMI,NDIMI,NDIMI,NDIMI)
        CALL WRTMAT(XMATO,NDIMO,NDIMO,NDIMO,NDIMO)
      END IF
*
      RETURN
      END
      SUBROUTINE ADIAJ_STR(IORB,JORB,ISTR_IN,NEL,LNUM,IZ,IREO,NTOOB,
     &                    ISTR_OUT,INUM,SIGN)
*
* Find string occupation (and address if LNUM.EQ.1) 
* for a+i aj |ISTR_IN> 
*
* It as assumed that the excitation is nonvanishing
*
* Jeppe Olsen, April 99
* 
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISTR_IN(NEL)
      INTEGER IZ(*),IREO(*)
*. Output
      INTEGER ISTR_OUT(NEL)
*
      NTEST = 000
C?    IF(NTEST.GE.100) THEN
C?      WRITE(6,*) 'ADIAJ_STR : Input string '
C?      CALL IWRTMA(ISTR_IN,1,NEL,1,NEL)
C?      WRITE(6,*) ' IORB, JORB ', IORB,JORB
C?    END IF
*
      IEL = 0
      IPLACED = 0
      DO IIEL = 1, NEL
        IF(ISTR_IN(IIEL).GT.IORB .AND. IPLACED.EQ.0 ) THEN
*. Add IORB
          IF(MOD(IIEL-1,2).EQ.0) THEN
            SIGNI =  1.0D0
          ELSE
            SIGNI = -1.0D0
          END IF
*
          IEL = IEL + 1
          ISTR_OUT(IEL) = IORB
          IPLACED = 1
        END IF
*
        IF( ISTR_IN(IIEL).NE.JORB) THEN
          IEL = IEL + 1
          ISTR_OUT(IEL) = ISTR_IN(IIEL)
        ELSE
          IF(MOD(IIEL-1,2).EQ.0) THEN
            SIGNJ = 1.0D0
          ELSE
            SIGNJ = -1.0D0
          END IF
        END IF
      END DO
*. Well, it could be that orbital i should be added as last elec
      IF(IPLACED.EQ.0) THEN
        ISTR_OUT(NEL) = IORB
        IF(MOD(NEL+1-1,2).EQ.0) THEN
          SIGNI = 1.0D0
        ELSE
          SIGNI = -1.0D0
        END IF
      END IF
      SIGN = SIGNI*SIGNJ
      IF(IORB.GT.JORB) SIGN = - SIGN
*
C?    WRITE(6,*) ' Output string '
C?    CALL IWRTMA(ISTR_OUT,1,NEL,1,NEL)
      IF(LNUM.EQ.1) THEN
       INUM = ISTRNM(ISTR_OUT,NTOOB,NEL,IZ,IREO,1)
      END IF
*. And address
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Output from ADIAJ_STR '
        WRITE(6,*) '======================'
        WRITE(6,*) ' Iorb and Jorb (a+(Iorb) a(Jorb)) :',IORB,JORB
        WRITE(6,*) ' Input and output strings '
        CALL IWRTMA(ISTR_IN ,1,NEL,1,NEL)
        CALL IWRTMA(ISTR_OUT,1,NEL,1,NEL)
        WRITE(6,*) ' Sign = ', sign
        IF(LNUM.EQ.1) WRITE(6,*) ' String number = ', INUM
      END IF
*
      RETURN
      END
      SUBROUTINE HCONFDIA_BBM(NAEL,NBEL,IJAGRP,IJBGRP,
     &           IASM,IATP,IAOC,NIA,IBSM,IBTP,IBOC,NIB,
     &           JASM,JATP,JAOC,NJA,JBSM,JBTP,JBOC,NJB,H,CB,SB)
*
* Outer routine for orbital conserving part of Ham times vector
*
* Jeppe Olsen, April 15, 1999 - Snowing in Aarhus
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
      CALL LUCIAQENTER('HCONF')
*. Fetch K matrix
C     CALL GTJK(XDUM1       ,H,NTOOB,XDUM3       ,IREOTS)
      CALL GTJK(int_mb(KLZSCR),H,NTOOB,int_mb(KLZSCR),IREOTS)
*
      CALL HCONFDIA_BBS(NAEL,NBEL,IJAGRP,IJBGRP,
     &           IASM,IATP,IAOC,NIA,IBSM,IBTP,IBOC,NIB,
     &           JASM,JATP,JAOC,NJA,JBSM,JBTP,JBOC,NJB,
     &           NTOOB,H,IOBPTS,NOBPT,IREOTS,int_mb(KLZ(1)),
     &           int_mb(KLZ(2)),int_mb(KLZSCR),ISMFTO,NGAS,
     &           int_mb(KLOCSTR(1)),int_mb(KLOCSTR(2)),
     &           int_mb(KLREO(1)),int_mb(KLREO(2)),CB,SB) 
      CALL LUCIAQEXIT('HCONF')
*
      RETURN
      END
      SUBROUTINE HCONFDIA_BBS(NAEL,NBEL,IJASPGP,IJBSPGP,
     &           IASM,IATP,IAOCC,NIA,IBSM,IBTP,IBOCC,NIB,
     &           JASM,JATP,JAOCC,NJA,JBSM,JBTP,JBOCC,NJB,
     &           NTOOB,RK,IOBPTS,NOBPT,IREOTS,
     &           IAZ,IBZ,IZSCR,ISMFTO,NGAS,
     &           JASTR_OC,JBSTR_OC,IAREO,IBREO,VECIN,VECOUT)
*
* Orbital occupation conserving part of Hamiltonian times vector.
* Part of this hamiltonian that does not conserve spin orbital 
* occupations
*

*
* The part of the Hamiltonian that conserves orbital occupations 
* is 
*
* H = sum_i h_ii + 0.5 sum_{ij}    E_{ii}E_{jj} (ii|jj) 
*                + 0.5 sum(i.ne.j) E_{ij}E_{ji} (ij|ji)
*
*
* and the part that conserves orbital occupatione but not 
* spin orbital occupations is
*
* sum(i.ne.j) a+ia aja a+jb aib (ij!ji) 
*  
* In search for a better preconditioner / H0 for EN pert, 
* Jeppe Olsen, April 99
*
*
* Notice : Present version works with supergroups from lists,
*          can therefore not work with passive/active division 

*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
*.General input
      DIMENSION RK(NTOOB,NTOOB)
      INTEGER IOBPTS(MXPNGAS,*),NOBPT(*),IREOTS(*)
      INTEGER ISMFTO(*)
*. Specific input
      DIMENSION VECIN(NJA,NJB)
      DIMENSION IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
      DIMENSION IAREO(*), IBREO(*)
*. Scratch
      DIMENSION JASTR_OC(NAEL,*),JBSTR_OC(NBEL,*)
      DIMENSION IAZ(*),IBZ(*),IZSCR(*)
*. Local scratch
      INTEGER IASTR_OC(MXPORB),IBSTR_OC(MXPORB)
      INTEGER IXA(MXPORB),IXB(MXPORB),JXA(MXPORB),JXB(MXPORB)
*. Output
      DIMENSION VECOUT(NIA,NIB)
*
      NTEST =  000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from HCONFDIA_BBS'
        WRITE(6,*) ' IASM IATP ', IASM,IATP
        WRITE(6,*) ' JASM JATP ', JASM,JATP
        WRITE(6,*) ' IBSM IBTP ', IBSM,IBTP
        WRITE(6,*) ' JBSM JBTP ', JBSM,JBTP
        WRITE(6,*) ' IJASPGP, IJBSPGP', IJASPGP,IJBSPGP
        WRITE(6,*) ' NGAS = ', NGAS
        WRITE(6,*) ' Input block '
        CALL WRTMAT(VECIN,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' Initial value of sigma block '
        CALL WRTMAT(VECOUT,NIA,NIB,NIA,NIB)
*
        WRITE(6,*) 'IAOCC '
        CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
        WRITE(6,*) 'IBOCC '
        CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
        WRITE(6,*) 'JAOCC '
        CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
        WRITE(6,*) 'JBOCC '
        CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
*
      END IF
*
*
*. Obtain Reordering arrays for Ia,IB strings 
*
*. Arc weights for IA
      NTESTX = 0
      IATP_ABS = IBSPGPFTP(IJASPGP)-1+IATP
      JATP_ABS = IBSPGPFTP(IJASPGP)-1+JATP
      IBTP_ABS = IBSPGPFTP(IJBSPGP)-1+IBTP
      JBTP_ABS = IBSPGPFTP(IJBSPGP)-1+JBTP
      CALL WEIGHT_SPGP(IAZ,NGAS,NELFSPGP(1,IATP_ABS),NOBPT,
     &                 IZSCR,NTESTX)
*. Reorder array for IA strings
      CALL GETSTR_TOTSM_SPGP(IJASPGP,IATP,IASM,NAEL,NIASTR,
     &                       JASTR_OC,NTOOB,1,IAZ,IAREO)
*. Arc weight for IB
      CALL WEIGHT_SPGP(IBZ,NGAS,NELFSPGP(1,IBTP_ABS),NOBPT,
     &                 IZSCR,NTESTX)
*. Reorder array for IA strings
      CALL GETSTR_TOTSM_SPGP(IJBSPGP,IBTP,IBSM,NBEL,NIBSTR,
     &                       JBSTR_OC,NTOOB,1,IBZ,IBREO)
*
* String info for Ja, Jb : Actual string occ
*
*. Arc weight for JA
C     CALL WEIGHT_SPGP(JAZ,NGAS,JAOCC,NOBPT,ZSCR,NTESTX)
*. Occupation for JA strings
      IDUM  = 0
      CALL GETSTR_TOTSM_SPGP(IJASPGP,JATP,JASM,NAEL,NJASTR,
     &                       JASTR_OC,NTOOB,0,IDUM,IDUM)
*. Arc weight for JB
C     CALL WEIGHT_SPGP(JBZ,NGAS,JBOCC,NOBPT,ZSCR,NTEST)
*. Occupation for JB strings
      CALL GETSTR_TOTSM_SPGP(IJBSPGP,JBTP,JBSM,NBEL,NJBSTR,
     &                       JBSTR_OC,NTOOB,0,IDUM,IDUM)
*
      IJSM = MULTD2H(IASM,JASM)
*. Loop over orbital types of i and j
      DO ITP = 1, NGAS
      DO JTP = 1, NGAS
C?    WRITE(6,*) ' ITP JTP = ', ITP,JTP
*. Do a+ia a ja a+jb a ib connect string types
        IF(ITP.EQ.JTP) THEN
          IADEL = 0
          JADEL = 0
        ELSE
          IADEL = 1
          JADEL =-1
        END IF
        IBDEL = - IADEL
        JBDEL = - JADEL
*
        IAMOKAY = 1   
        DO KTP = 1, NGAS
          IF(KTP.NE.ITP.AND.KTP.NE.JTP) THEN
             IF(IAOCC(KTP).NE.JAOCC(KTP) .OR.
     &          IBOCC(KTP).NE.JBOCC(KTP)     ) IAMOKAY = 0
          END IF
        END DO
*
        IF(IAOCC(ITP).NE.JAOCC(ITP)+IADEL) IAMOKAY = 0 
        IF(IAOCC(JTP).NE.JAOCC(JTP)+JADEL) IAMOKAY = 0 
        IF(IBOCC(ITP).NE.JBOCC(ITP)+IBDEL) IAMOKAY = 0 
        IF(IBOCC(JTP).NE.JBOCC(JTP)+JBDEL) IAMOKAY = 0 
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' ITP JTP IAMOKAY : ', ITP,JTP,IAMOKAY
        END IF
        IF(IAMOKAY.EQ.1) THEN
*. Orbital range for I and J, 
          IOFF = IOBPTS(ITP,1)
          NIORB = NOBPT(ITP)
          JOFF = IOBPTS(JTP,1)
          NJORB  = NOBPT(JTP)
*
          IEND = IOFF + NIORB - 1
          JEND = JOFF + NJORB - 1
*
C?        WRITE(6,*) ' NJASTR NJBSTR ', NJASTR,NJBSTR
          DO JBSTR = 1, NJBSTR
            IZERO = 0
            CALL ISETVC(IXB,IZERO,NTOOB)
            CALL ISETVC(JXB,IZERO,NTOOB)
*. Expanded occupation in JB for i- and j- orbitals IXB, JXB
            DO IBEL = 1, NBEL
              IORB = JBSTR_OC(IBEL,JBSTR)
              IF(IOFF.LE.IORB.AND.IORB.LE.IEND) IXB(IORB)= 1
              IF(JOFF.LE.IORB.AND.IORB.LE.JEND) JXB(IORB)= 1
            END DO
*
            DO JASTR = 1, NJASTR
*. I orbitals not occupied in Ja and occupied in JB
             CALL ISETVC(IXA,IZERO,NTOOB)
             DO IEL = 1, NAEL
               IORB = JASTR_OC(IEL,JASTR)
               IF(IOFF.LE.IORB.AND.IORB.LE.IEND) IXA(IORB) = 1
             END DO
             NIACT = 0
             DO IORB = IOFF,IEND
              IF(IXA(IORB).EQ.0.AND.IXB(IORB).EQ.1) THEN
                NIACT = NIACT + 1
                IXA(NIACT) = IORB
              END IF
             END DO
*
*. Loop over J orbitals occupied on JA, Unoccupied in JB
             DO JEL = 1, NAEL
              JORB = JASTR_OC(JEL,JASTR)
              IF(JOFF.LE.JORB.AND.JORB.LE.JEND.AND.
     &           JXB(JORB).EQ.0) THEN
                JOBSM = ISMFTO(JORB)
                IOBSM = MULTD2H(IJSM,JOBSM)
                IF(MOD(JEL,2).EQ.0) THEN
                 SIGNJ = 1.0D0
                ELSE
                 SIGNJ = -1.0D0
                END IF
*. JORB is occupied in JA, unoccupied in JB 
*. Loop over Iorbital, check for sym
                DO IIORB = 1, NIACT
                 IORB = IXA(IIORB)
                 IF(ISMFTO(IORB).EQ.IOBSM.AND.IORB.NE.JORB) THEN
*. We have connection :  Find Ia and Ib strings
C                   ADIAJ_STR(IORB,JORB,ISTR_IN,NEL,LNUM,IZ,IREO,NTOOB,
C    &                    ISTR_OUT,INUM,SIGN)
                  CALL ADIAJ_STR(IORB,JORB,JASTR_OC(1,JASTR),NAEL,
     &                           1,IAZ,IAREO,NTOOB,IASTR_OC,IA,SIGNIA)
                  CALL ADIAJ_STR(JORB,IORB,JBSTR_OC(1,JBSTR),NBEL,
     &                           1,IBZ,IBREO,NTOOB,IBSTR_OC,IB,SIGNIB)
C                 XIJJI = RK(IREOTS(IORB),IREOTS(JORB))
                  XIJJI = RK(IORB,JORB)
                  VECOUT(IA,IB) = VECOUT(IA,IB)             
     &           + SIGNIA*SIGNIB*XIJJI * VECIN(JASTR,JBSTR) 
C?                WRITE(6,*) ' JASTR JBSTR IA IB', JASTR,JBSTR,IA,IB
C?                WRITE(6,*) ' SIGNIA,SIGNIB,XIJJI',SIGNIA,SIGNIB,XIJJI
C?                WRITE(6,*) ' VECOUT,VECIN', VECOUT(IA,IB),
C?   &                         VECIN(JASTR,JBSTR)
                 END IF
*                ^ End if IORB of interest
                END DO
*               ^ End of loop over IORB
              END IF
*             ^ End if JORB of interest
             END DO
*            ^ End of loop over JEL  
            END DO
*           ^ End of loop over JASTR
          END DO
*         ^ End of loop over JBSTR
        END IF
*       ^ End of combination ITP,JTP have connection
      END DO
      END DO
*     ^ End of loop over types of I and J
*
      RETURN
      END 
      SUBROUTINE HCONFINTV(LURHS,LUX,SHIFTG,SHIFT_DIAG,VECIN,VECOUT,
     &                LBLK,LUPROJ,LUPROJ2,LLUDIA)
*
* Solve  (HCONF+Shift)X = RHS
*
* Where HCONF is configuration conserving part of Hamiltonian
*
* If ICISTR.EQ.1 VECIN contains RHS, else RHS is assumed  on LURHS
* Output : solution is on LUX
*
*
* Jeppe Olsen, May 1999           
* 
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'clunit.inc'
      REAL*8  INPRDD
      LOGICAL CONVER
      INCLUDE 'cands.inc'
* SCRATCH files used : (LUSC3,LUSC34,LUSC35,LUSC37 ) <= Old
*                      LUSC36, LUSC37, LUSC38, LUSC39 <= New
*             
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cecore.inc'
      COMMON/H_OCC_CONS/IH_OCC_CONS
      INCLUDE 'cshift.inc'
*
      EXTERNAL H0TVM
      DIMENSION ERROR(100)

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HCNINV')
*
      NTEST = 00
      IH_OCC_CONS = 1
*. Ensure that standard integrals are in WORK(KINT1)
CE    CALL SWAPVE(WORK(KINT1),WORK(KINT1O),NINT1)
*
* 2 : Solve linear set of equations
* ==================================
*
      ZERO = 0.0D0
      IF(LBLK.GT.0 ) THEN
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        STOP       ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
*
      ELSE IF(LBLK.LE.0)   THEN
*
*. Use path allowing segments of vectors
*
* corresponding eigenvalue info
        TEST = 
     &  SQRT(THRES_E) * SQRT(INPRDD(VECIN,VECOUT,LURHS,LURHS,1,-1))
        ILNPRT = NTEST
        CONVER = .FALSE.
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXIT_LOC ',MXIT_LOC
        SHIFT = SHIFTG 
        IPROJ = 0
C?      WRITE(6,*) ' LUX LURHS, LUSC38, LUSC38, LUSC39 LLUDIA',
C?   &               LUX,LURHS, LUSC38, LUSC38, LUSC39,LLUDIA 
C?      WRITE(6,*) ' LUPROJ2, LUPROJ ',  LUPROJ2, LUPROJ A
        WRITE(6,*) ' SHIFTG and SHIFT_DIAG', SHIFT,SHIFT_DIAG
C       SHIFTX = SHIFT_DIAG+ECORE_ORIG-ECORE
        SHIFTX = SHIFT_DIAG
        WRITE(6,*) ' SHIFTX before call to MICGCG ', SHIFTX
        CALL MICGCG(H0TVM,LUX,LURHS,LUSC37,LUSC38,LUSC39,LLUDIA,
     &              VECIN,VECOUT,MXIT_LOC,
     &              CONVER,TEST,SHIFTX,ERROR,NDIM,LUPROJ,LUPROJ2,
     &              VFINAL,ILNPRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Solution to linear set of Equations '
          CALL WRTVCD(VECIN,LUX,1,LBLK)
        END IF
*
      END IF
*
CE    CALL SWAPVE(WORK(KINT1),WORK(KINT1O),NINT1)
      IH_OCC_CONS = 0
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HCNINV')
      RETURN
      END 
      SUBROUTINE IB_FOR_SEL_ORBSPC(NOBPTS,NOBPS_SEL,IOBPTS_SEL,I_SEL,
     &           NGAS,MXPNGAS,NSYM)
*
* Obtain number of orbitals per symmetry and 
* offsets for orbitals in TS ordering when only 
* selected orbital spaces ( as defined by I_SEL) are included
*
*. Jeppe Olsen, Sept 2005
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NOBPTS(MXPNGAS,*),I_SEL(NGAS)
*. Output
      INTEGER IOBPTS_SEL(MXPNGAS,*), NOBPS_SEL(NSYM)
*
      IZERO = 0 
      CALL ISETVC(NOBPS_SEL,IZERO,NSYM)
      IOFF = 1
      DO IGAS = 1, NGAS
        DO ISYM = 1, NSYM
          IOBPTS_SEL(IGAS,ISYM) = IOFF
          IF(I_SEL(IGAS).EQ.1) THEN
             IOFF = IOFF + NOBPTS(IGAS,ISYM)
             NOBPS_SEL(ISYM) = NOBPS_SEL(ISYM) + NOBPTS(IGAS,ISYM)
          END IF
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Offsets for TS ordered selected orbspaces'
        CALL IWRTMA(IOBPTS_SEL,MXPNGAS,NSYM,MXPNGAS,NSYM)
        WRITE(6,*) ' Number of orbitals per sym, selected spaces'
        CALL IWRTMA(NOBPS_SEL,1,NSYM,1,NSYM)
      END IF
*
      RETURN
      END
      SUBROUTINE IREO_DACT_TS(NOBPTS,IOBPTS_SEL,I_SEL,
     &           IDTFREO,IFTDREO,NGAS,MXPNGAS,NSMOB)
*
* A set of active orbital spaces is given in I_SEL. Obtain 
* IDTFREL : Reorder array : Density order => Full TS order 
* IFTDREL : Reorder array : Full ST order => density  order 
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NOBPTS(MXPNGAS,*),I_SEL(NGAS)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*. Output
      INTEGER IDTFREO(*),IFTDREO(*)
*
      IOB_TS = 0
      NDACTORB = 0
      DO IGAS = 1, NGAS
        DO ISYM = 1, NSMOB
          DO IOB = 1, NOBPTS(IGAS,ISYM)
            IOB_TS = IOB_TS + 1
            IF(I_SEL(IGAS).EQ.1) THEN  
              IOB_DACT = IOBPTS_SEL(IGAS,ISYM)-1+ IOB
              IFTDREO(IOB_TS) = IOB_DACT
              IDTFREO(IOB_DACT) = IOB_TS
              NDACTORB = NDACTORB + 1
            ELSE 
              IFTDREO(IOB_TS) = 0
            END IF
          END DO
        END DO
      END DO
      NOB_TOT = IOB_TS
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reorder array : FULL TS => Active '
        CALL IWRTMA(IFTDREO,NOB_TOT,1,NOB_TOT,1)
        WRITE(6,*) ' Reorder array : Active => FULL TS  '
        CALL IWRTMA(IDTFREO,NDACTORB,1,NDACTORB,1)
      END IF
*
      RETURN
      END
      SUBROUTINE EXTR_SYMBLK_ACTMAT(AIN,AOUT,IJSM)
*
* A matrix AIN is given in complete form over active orbitals, 
* symmetry ordered
*
* extract symmetry blocks with symmetry IJSM
*
* Jeppe Olsen, Feb. 2011 from REORHO1 
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*. Input
      DIMENSION AIN(NACOB,NACOB)
*. Output
      DIMENSION AOUT(*)
*
      IMOFF = 1
      DO ISM = 1, NSMOB
       JSM = MULTD2H(ISM,IJSM)
*. Offsets for active orbitals with symmetries ISM, JSM
       IOFF = 1
       DO IISM = 1, ISM -1
        IOFF = IOFF + NACOBS(IISM)
       END DO
       JOFF = 1
       DO JJSM = 1, JSM - 1
        JOFF = JOFF + NACOBS(JJSM)
       END DO
*
        NI  = NACOBS(ISM)
        NJ =  NACOBS(JSM)
        DO I = 1, NI
          DO J = 1, NJ
            AOUT(IMOFF-1+(J-1)*NI+I) = AIN(IOFF-1+I,JOFF-1+J)
          END DO
        END DO
        IMOFF = IMOFF + NI*NJ
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' EXTR_SYMBLK_ACTMAT in action '
        WRITE(6,*) ' Symmetry of blocks extracted ',IJSM
        WRITE(6,*) ' Input matrix '
        CALL WRTMAT(AIN,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*)
        WRITE(6,*) ' extracted blocks : '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
        CALL PRHONE(AOUT,NACOBS,IJSM,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE GETINCN2_A(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,IJKLOF,I2INDX,
     &                  ICOUL,CFAC,EFAC) 
*
* Obtain integrals 
*
*     ICOUL = 0 : 
*                  XINT(IK,JL) = CFAC*(IJ!KL)         for IXCHNG = 0
*                              = CFAC*(IJ!KL)-EFAC*(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 : 
*                  XINT(IJ,KL) = CFAC*(IJ!KL)         for IXCHNG = 0
*                              = CFAC*(IJ!KL)-EFAC*(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 2 :  XINT(IL,JK) = CFAC*(IJ!KL)         for IXCHNG = 0
*                              = CFAC*(IJ!KL)-EFAC*(IL!KJ) for IXCHNG = 1
*
*. Integrals in output block are in unpacked  unless IKSM or JLSM .ne.0
* Storing for ICOUL = 1,2 not working if IKSM or JLSM .ne. 0 
* 
*
* Version for general integral array with dimensions stored in *A arrays
*
*
* Jeppe Olsen, The Lucia growing up campaign, May 2011
*
* Last modification, Oct. 1 2012, Jeppe Olsen, Correcting a bug for ICOUL = 2 (occuring using PH)
*
* type = -1 => all orbitals of given symmetry
* type =  0 => all inactive orbitals of given symmetry
* type = 1-ngas: all orbitals of a given gas and symmetry
* type = ngas + 1: all secondary orbitals of given symmetry
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLof(NsmOB,NsmOb,NsmOB)
*. Pair of orbital indeces ( symmetry ordered ) => address in symmetry packed 
*. matrix ( not active pt)
CNW   Dimension I2INDX(*)
      integer I2INDX
*.Output
      DIMENSION XINT(*)
*. Local scratch      
      DIMENSION IJARR(MXPORB)
*
      NTEST = 000
*. FUSKING
C?    IF(ISM.EQ.2.AND.JSM.EQ.1.AND.KSM.EQ.1.AND.LSM.EQ.2.AND.
C?   &   ITP.EQ.3.AND.JTP.EQ.3.AND.KTP.EQ.2.AND.LTP.EQ.2     ) THEN
C?       NTEST = 10000
C?       WRITE(6,*) ' Jeppe raised NTEST in GETINCN2_A '
C?    END IF
*
      call luciaqenter('gin2a')
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GETINCN2_A'
        WRITE(6,*) ' ======================'
        WRITE(6,'(A,4I3)') ' ITP, JTP, KTP, LTP = ', ITP, JTP, KTP, LTP
        WRITE(6,'(A,4I3)') ' ISM, JSM, KSM, LSM = ', ISM, JSM, KSM, LSM
        WRITE(6,'(A,2I3)') ' IKSM, JLSM = ', IKSM, JLSM
        WRITE(6,'(A,2I3)') ' IXCHNG, ICOUL = ', IXCHNG, ICOUL
        WRITE(6,*) ' CFAC, EFAC = ', CFAC, EFAC
        WRITE(6,*)
        WRITE(6,'(A,3I3)')  ' I12S_A, I34S_A, I1234S_A = ',
     &                        I12S_A, I34S_A, I1234S_A
CNW     WRITE(6,*) ' First integral in input block = ', INTLST(1)
      END IF
*
*Number of orbitals for the four indeces
*
      IF(ITP.EQ.-1) THEN
       NI = NTOOBS_IA(ISM)
      ELSE
       NI = NOBPTS_GN_A(ITP,ISM,1)
      END IF
*
      IF(JTP.EQ.-1) THEN
       NJ = NTOOBS_JA(JSM)
      ELSE
       NJ = NOBPTS_GN_A(JTP,JSM,2)
      END IF
      IF(NTEST.GE.100) WRITE(6,*) 'JSM, JTP, NJ = ', JSM, JTP, NJ
*
      IF(KTP.EQ.-1) THEN
       NK = NTOOBS_KA(KSM)
      ELSE
       NK = NOBPTS_GN_A(KTP,KSM,3)
      END IF
*
      IF(LTP.EQ.-1) THEN
       NL = NTOOBS_LA(LSM)
      ELSE
       NL = NOBPTS_GN_A(LTP,LSM,4)
      END IF
*
*. Offsets relative to start of orbitals, given symmetry
*
      IOFF = 1
      DO IITP = 0, ITP-1
        IOFF = IOFF + NOBPTS_GN_A(IITP,ISM,1)
      END DO
*
      JOFF = 1
      DO JJTP = 0, JTP-1
        JOFF = JOFF + NOBPTS_GN_A(JJTP,JSM,2)
      END DO
*
      KOFF = 1
      DO KKTP = 0, KTP-1
        KOFF = KOFF + NOBPTS_GN_A(KKTP,KSM,3)
      END DO
*
      LOFF = 1
      DO LLTP = 0, LTP-1
        LOFF = LOFF + NOBPTS_GN_A(LLTP,LSM,4)
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,'(A,4I3)') '  NI, NJ, NK, NL = ', NI, NJ, NK, NL
        WRITE(6,'(A,4I3)') 
     & ' IOFF, JOFF, KOFF, LOFF = ', IOFF, JOFF, KOFF, LOFF
      END IF
    
*
* Symmetry of input block
*
      ISM_I = ISM
      JSM_I = JSM
      IF(I12S_A.EQ.1.AND.ISM.LT.JSM) THEN
        ISM_I = JSM
        JSM_I = ISM
      END IF
      IF(I12S_A.EQ.1) THEN
       IJBLK = ISM_I*(ISM_I-1)/2+JSM_I
      ELSE
       IJBLK = (ISM_I-1)*NSMOB + JSM_I
      END IF
*
      KSM_I = KSM
      LSM_I = LSM
      IF(I34S_A.EQ.1.AND.KSM.LT.LSM) THEN
        KSM_I = LSM
        LSM_I = KSM
      END IF
      IF(I34S_A.EQ.1) THEN
       KLBLK = KSM_I*(KSM_I-1)/2 + LSM_I
      ELSE
       KLBLK = (KSM-1)*NSMOB + LSM
      END IF
      I1234P = 0
      IF(I1234S_A.EQ.1.AND.IJBLK.LT.KLBLK) I1234P = 1
      IF(I1234P.EQ.0) THEN
        IBLOFF = IJKLOF(ISM_I,JSM_I,KSM_I)
      ELSE
        IBLOFF = IJKLOF(KSM_I,LSM_I,ISM_I)
      END IF
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' ISM_I, JSM_I, KSM_I, LSM_I ',
     &               ISM_I, JSM_I, KSM_I, LSM_I 
        WRITE(6,*) ' I1234P = ', I1234P 
        WRITE(6,*) ' IBLOFF = ', IBLOFF
      END IF
*
*     Collect Coulomb terms
*
*
      NIS = NTOOBS_IA(ISM)
      NJS = NTOOBS_JA(JSM)
      NKS = NTOOBS_KA(KSM)
      NLS = NTOOBS_LA(LSM)
      IF(NTEST.GE.1000) THEN
        WRITE(6,'(A,4I4)') ' NIS, NJS, NKS, NLS = ',
     &  NIS, NJS, NKS, NLS
      END IF
*
      If(I12S_A.EQ.1.AND.ISM.EQ.JSM) THEN
       NIJS = NIS*(NIS+1)/2
      ELSE
       NIJS = NIS*NJS
      END IF
*
      IF(I34S_A.EQ.1.AND.KSM.EQ.LSM) THEN
        NKLS = NKS*(NKS+1)/2
      ELSE
        NKLS = NKS*NLS
      END IF
*
      IF(IKSM.EQ.1) THEN
       NIK = NI*(NI+1)/2
      ELSE
       NIK = NI*NK
      END IF
*
      IF(JLSM.EQ.1) THEN
       NJL = NJ*(NJ+1)/2
      ELSE
       NJL = NJ*NL
      END IF
      NIJKL = NIK*NJL
*
      ZERO = 0.0d0
      CALL SETVEC(XINT,ZERO,NIJKL)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' Coulomb part '
        WRITE(6,*)
      END IF
*
*. Loop over integrals in output block, order corresponds to 
*. ICOUL = 0, i.e. [IK!JL]
*
      IINT = 0
      DO L = 1, NL
        JMIN=1
        IF (JLSM.NE.0 ) JMIN = L    
        DO J = JMIN, NJ 
          J_ABS = JOFF - 1 + J
*. (Set up array IJ*(IJ-1)/2 )
*
          DO K = 1, NK
            IMIN = 1
            IF(IKSM.NE.0) IMIN = K
*. Index KL in input integral block (IJ!KL)
            K_ABS = KOFF + K - 1
            L_ABS = LOFF + L - 1
            IF(I34S_A.EQ.1.AND.KSM.EQ.LSM) THEN
              KL = MAX(K_ABS,L_ABS)*(MAX(K_ABS,L_ABS)-1)/2
     &           + MIN(K_ABS,L_ABS)
            ELSE IF(KSM.GT.LSM.OR.I34S_A.EQ.0) THEN
              KL = (L_ABS-1)*NKS + K_ABS
            ELSE
              KL = (K_ABS-1)*NLS + L_ABS
            END IF
            IF(NTEST.GE.10000) THEN
              WRITE(6,'(A,4I4)') ' K, L, K_ABS, L_ABS = ',
     &                             K, L, K_ABS, L_ABS
              WRITE(6,*) ' KL = ', KL
            END IF
*
            IF(ICOUL.EQ.1)  THEN  
*. Address in output block before integral (1j!kl)
                IINT = (L-1)*NI*NJ*NK
     &               + (K-1)*NI*NJ    
     &               + (J-1)*NI
            ELSE IF (ICOUL.EQ.2) THEN
*  Address in output block before (1L!JK) 
                IINT = (K-1)*NI*NL*NJ
     &               + (J-1)*NI*NL          
     &               + (L-1)*NI             
            END IF
*. IF ICOUL = 0, the increment of IINT defines this 
*
            IF(NTEST.GE.10000)
     &      WRITE(6,*) ' IJBLK, KLBLK = ', IJBLK, KLBLK
            IF(I1234S_A.EQ.0.OR.IJBLK.GT.KLBLK) THEN
*
*. Contributions from input block (ISM JSM ! KSM LSM ) 
*. with (Ism,jsm) > (ksm,lsm) ( or no restrictions)
*
*. Address -1  of (1 1 !K L) in input block
              IJKL0 = IBLOFF-1+(KL-1)*NIJS
              DO I_ABS = IOFF-1+IMIN, IOFF-1+NI 
*Index IJ in input block
                IF(I12S_A.EQ.0.OR.ISM.GT.JSM) THEN
                 IJ = (J_ABS-1)*NIS + I_ABS
                ELSE IF(ISM.EQ.JSM) THEN
                 IJ = MAX(I_ABS,J_ABS)*(MAX(I_ABS,J_ABS)-1)/2
     &              + MIN(I_ABS,J_ABS)
                ELSE 
*. We enter here if I12S_A.EQ.1.AND.ISM.LT.JSM
                 IJ = (I_ABS-1)*NJS + J_ABS
                END IF
                IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' I_ABS, J_ABS, IJ = ', I_ABS, J_ABS, IJ
                END IF
                IJKL = IJKL0 + IJ
                IINT = IINT + 1
*
                IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' IJKL0, IJKL, IINT = ',
     &                         IJKL0, IJKL, IINT
                  WRITE(6,'(A,4I4)') ' I_ABS, J_ABS, K_ABS, L_ABS = ', 
     &                         I_ABS, J_ABS, K_ABS, L_ABS
                END IF
*
                XINT(IINT) = CFAC*INTLST(IJKL)   
                IF(NTEST.GE.10000)
     &          WRITE(6,*) ' IINT, XINT = ', IINT, XINT(IINT)
              END DO
            END IF
*
*. block (ISM JSM !ISM JSM) with I1234S_A = 1
*
            IF(I1234S_A.EQ.1.AND.IJBLK.EQ.KLBLK) THEN
              KLOFF = KL*(KL-1)/2
              IJKL0 = (KL-1)*NIJS-KLOFF
              DO I_ABS = IOFF-1+IMIN, IOFF + NI - 1
*. IJ in input block
                IF(I12S_A.EQ.0.OR.ISM.GT.JSM) THEN
                  IJ = (J_ABS-1)*NIS + I_ABS
                ELSE IF(ISM.EQ.JSM) THEN
                  IJ = MAX(I_ABS,J_ABS)*(MAX(I_ABS,J_ABS)-1)/2
     &               + MIN(I_ABS,J_ABS)
                ELSE
                  IJ = (I_ABS-1)*NJS + J_ABS
                END IF
                IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' I_ABS, J_ABS = ', I_ABS, J_ABS
                END IF
                IF(IJ.GE.KL) THEN    
C                 IJKL=IJ+(KL-1)*NIJS-KLOFF
                  IJKL = IJKL0 + IJ
                Else
                  IJKL=KL+(IJ-1)*NKLS - IJ*(IJ-1)/2
                END IF
                IINT=IINT+1
                IF(NTEST.GE.10000) THEN
                  WRITE(6,'(A,4I6)') 
     &            ' IINT, IJ, KL, IJKL' ,  IINT, IJ, KL, IJKL
                  WRITE(6,'(A,4I4)') ' I_ABS, J_ABS, K_ABS, L_ABS = ', 
     &                         I_ABS, J_ABS, K_ABS, L_ABS
                END IF
                XINT(IINT) = CFAC*INTLST(IBLOFF-1+IJKL)
                IF(NTEST.GE.10000)
     &          WRITE(6,*) ' IINT, XINT = ', IINT, XINT(IINT)
              END DO
            END IF
*
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) < (ksm,lsm)
            IF(I1234S_A.EQ.1.AND.IJBLK.LT.KLBLK) THEN 
              IJKL0 = IBLOFF-1+KL - NKLS
              DO I_ABS = IOFF-1+IMIN, IOFF + NI - 1
* IJ in input block
                IF(I12S_A.EQ.0.OR.ISM.GT.JSM) THEN
                  IJ = (J_ABS-1)*NIS + I_ABS
                ELSE IF(ISM.EQ.JSM) THEN
                  IJ = MAX(I_ABS,J_ABS)*(MAX(I_ABS,J_ABS)-1)/2
     &               + MIN(I_ABS,J_ABS)
                ELSE 
                  IJ = (I_ABS-1)*NJS + J_ABS
                END IF
                IJKL = IJKL0 + IJ*NKLS
                IINT=IINT+1
                IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' I_ABS, J_ABS, IJ = ', I_ABS, J_ABS, IJ
                  WRITE(6,*) ' IINT, IJKL = ', IINT, IJKL
                END IF
                XINT(IINT) = CFAC*INTLST(IJKL)
                IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' Updated integral = ', XINT(IINT)
                END IF
              END DO
            END IF
*
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Integral block after Coulomb '
        CALL WRTMAT(XINT,1,NIJKL,1,NIJKL)
      END IF
*
*     Collect Exchange terms  [IK!JL] <= (IL!KJ)
*
      IF ( IXCHNG.NE.0 ) THEN
C     WRITE(6,*) ' Warning, Exchange part has been blocked '
C     IF ( IXCHNG.EQ.-127 ) THEN
        IF(NTEST.GE.1000) THEN
          WRITE(6,*)
          WRITE(6,*) ' Exchange part '
          WRITE(6,*)
        END IF
*we will fetch integrals (IL!KJ) so
        NIS = NTOOBS_IA(ISM)
        NJS = NTOOBS_LA(JSM)
        NKS = NTOOBS_KA(KSM)
        NLS = NTOOBS_JA(LSM)
*
*. Offsets relative to start of orbitals, given symmetry
*
      IOFF = 1
      DO IITP = 0, ITP-1
        IOFF = IOFF + NOBPTS_GN_A(IITP,ISM,1)
      END DO
*
      JOFF = 1
      DO JJTP = 0, JTP-1
        JOFF = JOFF + NOBPTS_GN_A(JJTP,JSM,4)
      END DO
*
      KOFF = 1
      DO KKTP = 0, KTP-1
        KOFF = KOFF + NOBPTS_GN_A(KKTP,KSM,3)
      END DO
*
      LOFF = 1
      DO LLTP = 0, LTP-1
        LOFF = LOFF + NOBPTS_GN_A(LLTP,LSM,2)
      END DO
*
        IF(I12S_A.EQ.1.AND.ISM.EQ.LSM) THEN
         NILS = NIS*(NIS+1)/2
        ELSE
         NILS = NIS*NLS         
        END IF
*
        IF(I34S_A.EQ.1.AND.KSM.EQ.JSM) THEN
          NKJS = NKS*(NKS+1)/2
        ELSE
          NKJS = NKS*NJS       
        END IF
*
        ISM_I = ISM
        LSM_I = LSM
        IF(I12S_A.EQ.1.AND.ISM.LT.LSM) THEN
         ISM_I = LSM
         LSM_I = ISM
        END IF
*
        KSM_I = KSM
        JSM_I = JSM
        IF(I34S_A.EQ.1.AND. KSM.LT.JSM) THEN
         KSM_I = JSM
         JSM_I = KSM
        END IF
*
        IF(I12S_A.EQ.1) THEN
          ILBLK = ISM_I*(ISM_I-1)/2 + LSM_I
        ELSE
          ILBLK = (ISM_I-1)*NSMOB + LSM_I
        END IF
        IF(I34S_A.EQ.1) THEN
          KJBLK = KSM_I*(KSM_I-1)/2 + JSM_I
        ELSE
          KJBLK = (KSM_I-1)*NSMOB + JSM_I
        END IF
        IF(NTEST.GE.10000) THEN
          WRITE(6,*) ' ILBLK, KJBLK = ', ILBLK, KJBLK
        END IF
*. Start of block input block (IL!KJ)
        I1234P = 0
        IF(I1234S_A.EQ.0.OR.ILBLK.GE.KJBLK) THEN
          IBLOFF = IJKLOF(ISM_I,LSM_I,KSM_I)
        ELSE 
          IBLOFF = IJKLOF(KSM_I,JSM_I,ISM_I) 
          I1234P = 1
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' ISM_I, JSM_I, KSM_I, LSM_I = ',
     &                 ISM_I, JSM_I, KSM_I, LSM_I  
          WRITE(6,*) ' IBLOFF = ', IBLOFF
        END IF
*
        IINT=0
        DO L = 1, NL     
          L_ABS = LOFF -1 + L
          JMIN = 1  
          IF (JLSM.NE.0) JMIN = L    
*
          DO J = JMIN, NJ              
            J_ABS = JOFF -1 + J
            DO K = 1, NK                   
              K_ABS = KOFF -1 + K
*. Index KJ in (IL!KJ)
              IF(I34S_A.EQ.0.OR.KSM.GT.JSM) THEN
                KJ = (J_ABS-1)*NKS + K_ABS
              ELSE IF (KSM.EQ.JSM) THEN
                KJ = MAX(K_ABS,J_ABS)*(MAX(K_ABS,J_ABS)-1)/2
     &             + MIN(K_ABS,J_ABS)
              ELSE
                KJ = (K_ABS-1)*NJS + J_ABS
              END IF
*  [IJ!KL] <= (IL!KJ)
*
              IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l) in output block
                  IINT = (L-1)*NK*NJ*NI      
     &                 + (K-1)   *NJ*NI    
     &                 + (J-1)      *NI
              ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK)  in output block
                IINT = (K-1)*NI*NL*NJ
     &               + (J-1)*NI*NL          
     &               + (L-1)*NI             
CER             IINT = (K-1)*NJ*NK*NI      
CER  &               + (J-1)   *NK*NI           
CER  &               + (L-1)      *NI
              END IF
*
              IMIN = 1    
              If(IKSM.NE.0) IMIN = K      
*
              IF(I1234S_A.EQ.0.OR.ILBLK.GT.KJBLK) THEN 
                ILKJ0 = IBLOFF-1+(KJ-1)*NILS
                DO I_ABS = IOFF - 1 + IMIN, IOFF - 1 + NI
*. Index IL in (IL!KJ)
                  IF(I12S_A.EQ.0.OR.ISM.GT.LSM) THEN
                    IL = (L_ABS-1)*NIS + I_ABS
                  ELSE IF (ISM.EQ.LSM) THEN
                    IL = MAX(I_ABS,L_ABS)*(MAX(I_ABS,L_ABS)-1)/2
     &                 + MIN(I_ABS,L_ABS)
                  ELSE
                    IL = (I_ABS-1)*NLS + L_ABS
                  END IF
                  ILKJ = ILKJ0 + IL
                  IINT=IINT+1
                  XINT(IINT)=XINT(IINT)-EFAC*INTLST(ILKJ)
CNW               IF(NTEST.GE.10000) THEN 
CNW                 WRITE(6,*) 
CNW  &            ' IINT, ILKJ, INTLST(ILKJ), XINT(IINT) =',
CNW  &              IINT, ILKJ, INTLST(ILKJ), XINT(IINT)
CNW               END IF
                END DO
              ELSE IF(ILBLK.EQ.KJBLK) THEN
                ILKJ0 = (KJ-1)*NILS - KJ*(KJ-1)/2
                DO I_ABS = IOFF-1+IMIN, IOFF-1+NI   
*. Index IL in (IL!KJ)
                  IF(I12S_A.EQ.0.OR.ISM.GT.LSM) THEN
                    IL = (L_ABS-1)*NIS + I_ABS
                  ELSE IF (ISM.EQ.LSM) THEN
                    IL = MAX(I_ABS,L_ABS)*(MAX(I_ABS,L_ABS)-1)/2
     &                 + MIN(I_ABS,L_ABS)
                  ELSE
                    IL = (I_ABS-1)*NLS + L_ABS
                  END IF
                  IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' ISM, LSM, I12S_A =', ISM, LSM, I12S_A
                  WRITE(6,*) ' I_ABS, L_ABS, IL = ', I_ABS, L_ABS, IL
                  END IF
                  IF ( IL.GE.KJ ) THEN
C                     ILKJ = IL + (KJ-1)*NILS - KJ*(KJ-1)/2
                      ILKJ = IL + ILKJ0
                  ELSE
                      ILKJ = KJ + (IL-1)*NKJS - IL*(IL-1)/2
                  END IF
*
                  IINT=IINT+1
                  IF(NTEST.GE.10000) THEN
                   WRITE(6,'(A,5I4)') 
     &                        ' I_ABS, J_ABS, K_ABS, L_ABS, IINT = ',
     &                          I_ABS, J_ABS, K_ABS, L_ABS, IINT
                   WRITE(6,*) ' IL, KJ, ILKJ = ', IL, KJ, ILKJ
                   WRITE(6,*) ' IINT, ILKJ = ', IINT, ILKJ
                  END IF
*
                  XINT(IINT)=XINT(IINT)-EFAC*INTLST(IBLOFF-1+ILKJ)
                  IF(NTEST.GE.10000)
     &            WRITE(6,*) ' updated integral ', XINT(IINT)
                END DO
              ELSE IF(ILBLK.LT.KJBLK) THEN
                ILKJ0 = IBLOFF-1+KJ-NKJS
                DO I_ABS = IOFF - 1 + IMIN, IOFF - 1 + NI
                  IF(I12S_A.EQ.0.OR.ISM.GT.LSM) THEN
                    IL = (L_ABS-1)*NIS + I_ABS
                  ELSE IF (ISM.EQ.LSM) THEN
                    IL = MAX(I_ABS,L_ABS)*(MAX(I_ABS,L_ABS)-1)/2
     &                 + MIN(I_ABS,L_ABS)
                  ELSE
                    IL = (I_ABS-1)*NLS + L_ABS
                  END IF
                  ILKJ = ILKJ0 + IL*NKJS
                  IINT = IINT + 1
                  XINT(IINT)=XINT(IINT)-EFAC*INTLST(ILKJ)
                END DO
              END IF
*
            END DO
          END DO
        END DO
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Integral block after exchange '
          CALL WRTMAT(XINT,1,NIJKL,1,NIJKL)
        END IF
*
      END IF
      call luciaqexit('gin2a')
*
      RETURN
      END
      SUBROUTINE GSBBD2BN(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXK,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  CJRES,SIRES,NORB,NTESTG,SCLFAC,S2_TERM1,
     &                  IDOSRHO2,RHO2AB,
     &                  NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
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
*              Version using ADAST_GAS, January 2011
*
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nproc
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
      INTEGER IOBPTS_SEL(MXPNGAS,*), IDACTSPC(NGAS)
*.Input
      DIMENSION CB(*),SB(*)
*. Output
      DIMENSION RHO2(*),RHO2AB(*)
*.Scratch
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION CJRES(*),SIRES(*)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*
      CALL LUCIAQENTER('GSD2N')
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' GSBBD2BN speaking '
        WRITE(6,*) ' ================== '
        WRITE(6,'(A,4(2X,I4))') ' Left: IASM, IBSM, IATP, IBTP =  ',
     &                                  IASM, IBSM, IATP, IBTP
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
        WRITE(6,'(A,4(2X,I4))') ' Right: JASM, JBSM, JATP, JBTP = ',
     &                                   JASM, JBSM, JATP, JBTP
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
      END IF
*
C?    WRITE(6,*) ' NJAS NJB = ',NJA,NJB
C?    WRITE(6,*) ' IAGRP IBGRP = ', IAGRP,IBGRP
C?    WRITE(6,*) ' MXPNGAS = ', MXPNGAS
C?    WRITE(6,*) ' NSMOB = ', NSMOB
      IROUTE = 3
*
*. Groups defining each supergroup
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
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
c     nproc = ga_nnodes()
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        IF(IDACTSPC(ITYP)+IDACTSPC(JTYP).NE.2) GOTO 2001
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS_SEL(ITYP,ISM)
          JOFF = IOBPTS_SEL(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all Ka strings
*. a+j!ka> = +/-/0 * !Ja>
          IJAC = 2
          KAFRST = 1
          IFRST = 1
          CALL ADAST_GAS(JSM,JTYP,NGAS,JASPGP,JASM,
     &                   I1,XI1S,NKASTR,IEND,IFRST,KAFRST,KACT,
     &                   SCLFAC,IJAC)
C?        WRITE(6,*) ' NKASTR = ', NKASTR
COLD      CALL ADSTN_GAS(JSM,JTYP,JATP,JASM,IAGRP,
COLD &                   I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,
COLD &                   SCLFAC)
*. a+i!ka> = +/-/0 * !Ia>
          ONE    = 1.0D0
          CALL ADAST_GAS(ISM,ITYP,NGAS,IASPGP,IASM,
     &                   I3,XI3S,NKASTR,IEND,IFRST,KAFRST,KACT,
     &                   ONE,IJAC)
COLD      CALL ADSTN_GAS(ISM,ITYP,IATP,IASM,IAGRP,
COLD &                   I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,
COLD &                   ONE   )
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
c         if(nkabtc.lt.nproc) then
c           nkabtc = min(nkaeff,nproc)
c         endif
c         DO 1801 IKABTC = ga_nodeid()+1, NKABTC,nproc
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. explicit zeroin, may be reduced
*. Obtain C(ka,J,JB) for Ka in batch
            ZERO = 0.0D0
            CALL SETVEC(CJRES,ZERO,LKABTC*NJ*NJB)
            DO JJ = 1, NJ
              CALL GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' C(Ka,J,Jb) '
              CALL WRTMAT(CJRES,LKABTC,NJ*NJB,LKABTC,NJ*NJB)
            END IF
*. Obtain S(ka,i,Ib) for Ka in batch
            CALL SETVEC(SIRES,ZERO,LKABTC*NI*NIB)
            DO II = 1, NI
              CALL GET_CKAJJB(SB,NI,NIA,SIRES,LKABTC,NIB,
     &             II,I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' S(Ka,i,Ib) '
              CALL WRTMAT(SIRES,LKABTC,NI*NIB,LKABTC,NI*NIB)
            END IF
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
              IF(IDACTSPC(KTYP)+IDACTSPC(LTYP).NE.2) GOTO 2000
*
              DO 1930 KSM = 1, NSMOB
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
C?              WRITE(6,*) ' Loop 1930, KSM LSM ',KSM,LSM
                KOFF = IOBPTS_SEL(KTYP,KSM)
                LOFF = IOBPTS_SEL(LTYP,LSM)
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
                KLAC = 2
                CALL ADAST_GAS(LSM,LTYP,NGAS,JBSPGP,JBSM,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KAFRST,KACT,ONE,KLAC)
C?        WRITE(6,*) ' NKBSTR = ', NKBSTR
COLD            CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
COLD &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADAST_GAS(KSM,KTYP,NGAS,IBSPGP,IBSM,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KAFRST,KACT,ONE,KLAC)
COLD            CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
COLD &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE)
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
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Before call to ADTOR2'
                  WRITE(6,*) ' RHO2B (X) matrix '
                  call wrtmat(x,ni*nj,nk*nl,ni*nj,nk*nl)
                END IF
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
C?            WRITE(6,*) ' ITYP, ISM, IOFF = ', ITYP, ISM, IOFF
C?            WRITE(6,*) ' JTYP, JSM, JOFF = ', JTYP, JSM, JOFF
C?            WRITE(6,*) ' KTYP, KSM, KOFF = ', KTYP, KSM, KOFF
C?            WRITE(6,*) ' LTYP, LSM, LOFF = ', LTYP, LSM, LOFF
                CALL ADTOR2(RHO2,X,2,
     &                NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NDACTORB)
                IF(IDOSRHO2.EQ.1) THEN
                  CALL ADTOR2S(RHO2AB,X,2,
     &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NDACTORB)
                END IF 
                IF(NTEST.GE.1000) THEN
                write(6,*) ' updated density matrix '
                call prsym(rho2,NDACTORB*NDACTORB)
                END IF

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
      CALL LUCIAQEXIT('GSD2N')
      RETURN
      END
      SUBROUTINE ADVICE_SIGMA3(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE,
     &                         NIA,NIB,NJA,NJB)
*
* Advice Sigma routine about best route to take
*
* ITERM : Term  to be studied :  
*         =1 alpha-beta term 
*         ....... ( to be continued )
*
* LADVICE : ADVICE given ( short, an integer !!)
*
* For ITERM = 1 : 
*           LADVICE = 1 : Business as usual, no transpose of matrix
*                         (resolution on alpha strings, direct exc on beta)
*           LADVICE = 2 = Transpose matrices
*                         (resolution on beta strings, direct exc on alpha)
*
* Jeppe Olsen, Version with info on string dimensions, May 2012
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Specific input
      INTEGER IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
*. Local Scratch
       DIMENSION ITP(16),JTP(16),KTP(16),LTP(16)
*
      NTEST = 00
      IF(ITERM.EQ.1) THEN
*.
*. sigma(i,Ka,Ib) = sum(i,kl)<Ib!Eb_kl!Jb>(ij!kl)C(j,Ka,Jb)
*
* Number of ops : Number of sx(kl) N_i*N_j_dimension of C(j,Ka,Jb)
*.No absolute calc of flops is made, only a relative measure
*
* Single excitations connecting the two types
*
C            SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
        CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOCC,JAOCC,IPHGAS)
        CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOCC,JBOCC,IPHGAS)
C?      WRITE(6,*) 'NIJTYP, NKLTYP', NIJTYP,NKLTYP
*. P-h modifications ( I cannot predict these at the moment
        IF(NIJTYP.GE.1.AND.NKLTYP.GE.1) THEN
*
        IF((IPHGAS(ITP(1)).EQ.2.AND.IPHGAS(JTP(1)).EQ.2).OR.
     &     (IPHGAS(KTP(1)).EQ.2.AND.IPHGAS(LTP(1)).EQ.2)     ) THEN
           IPHMODI = 1
         ELSE
           IPHMODI = 0
         END IF
        ELSE
           IPHMODI = 0
        END IF
          
*
        IF(IPHMODI.EQ.1.OR.NIJTYP.NE.1.OR.NKLTYP.NE.1
     &     .OR.IADVICE.EQ.0) THEN
*. Several connections, i.e. the alpha or the beta blocks are identical,
*. or ph modifications
*. just continue
          LADVICE = 1
        ELSE
* =========================================
*.. Index for flops along C(j,Ka,Jb) route
* =========================================
*.Dim of C(j,Ka,Jb) relative to C(Ja,Jb)
*. going from Ja to  Ka reduces occ by one elec, changes dim by n/(N-n+1)
          XNJOB = FLOAT(NOBPT(JTP(1)))
          XNJEL = FLOAT(JAOCC(JTP(1)))
          XCJKAJB = XNJOB*XNJEL/(XNJOB-XNJEL+1)
*. Number of kl excitations per beta string : 
          XNKLSX = FLOAT((NOBPT(KTP(1))-JBOCC(KTP(1)))*JBOCC(LTP(1)))
*. Number of ops (relative to dim of C)
          XNIOB = FLOAT(NOBPT(ITP(1)))
          XFLOPA = XCJKAJB*XNKLSX*XNIOB
* =========================================
*.. Index for flops along C(l,Ja,Kb) route
* =========================================
*.Dim of C(l,Ja,Kb) relative to C(Ja,Jb)
          XNLOB = FLOAT(NOBPT(LTP(1)))
          XNLEL = FLOAT(JBOCC(LTP(1)))
          XCLJAKB = XNLOB*XNLEL/(XNLOB-XNLEL+1)
*. Number of ij excitations per alpha string : 
          XNIJSX = FLOAT((NOBPT(ITP(1))-JAOCC(ITP(1)))*JAOCC(JTP(1)))
*. Number of ops (relative to dim of C)
          XNKOB = FLOAT(NOBPT(KTP(1)))
          XFLOPB = XCLJAKB*XNIJSX*XNKOB
*. Switch to second route if atleast 20 percent less work
          IF(XFLOPB.LE.0.8*XFLOPA) THEN
            LADVICE = 2
          ELSE
            LADVICE = 1
          END IF
*
* If the flop counts are nearly identical, but one route leads to fewer matrix 
* multiplies, choose this.
         IF(XFLOPB.LE.1.2*XFLOPA.AND.
     &      XNIJSX*FLOAT(NJA).LT.0.9*XNKLSX*FLOAT(NJB))  LADVICE = 2

*. Well, an additional consideration :
* If the C block involes the smallest allowed number of elecs in hole space,
* and the annihilation is in hole space
* then we do the annihilation in the space with the smallest number of 
* hole electrons.
          LHOLEA =0
          LHOLEB =0
          DO IGAS = 1, NGAS
            IF(IPHGAS(IGAS).EQ.2) THEN
              LHOLEA = LHOLEA + JAOCC(IGAS)
              LHOLEB = LHOLEB + JBOCC(IGAS)
            END IF
          END DO
*
          IF(LHOLEA+LHOLEB.EQ.MNHL.AND.
     &       (IPHGAS(JTP(1)).EQ.2.OR.IPHGAS(LTP(1)).EQ.2))  THEN
*
             IF(IPHGAS(JTP(1)).EQ.2) THEN
              KHOLEA = LHOLEA-1
              KHOLEB = LHOLEB 
             ELSE 
              KHOLEA = LHOLEA
              KHOLEB = LHOLEB - 1
             END IF
*
             IF(KHOLEA.EQ.KHOLEB) THEN
               LLADVICE = LADVICE
             ELSE IF(KHOLEA.LT.KHOLEB) THEN
               LLADVICE= 1
             ELSE
               LLADVICE = 2
             END IF
             IF(NTEST.GE.100.AND.LADVICE.NE.LLADVICE) THEN
               WRITE(6,*) ' Advice changed by hole considetions'
               WRITE(6,*) ' LADVICE, LLADVICE', LADVICE,LLADVICE
             END IF
             LADVICE = LLADVICE  
          END IF
*
*
C         IF(NTEST.GE.100) THEN
          IF(NTEST.GE.100.AND.LADVICE.EQ.2) THEN
            WRITE(6,*) ' ADVICE2 active '
            WRITE(6,*) ' IAOCC IBOCC JAOCC JBOCC'
            CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
            WRITE(6,*) ' ITP JTP KTP LTP ',ITP(1),JTP(1),KTP(1),LTP(1)
            WRITE(6,'(A,4(2X,E9.3))') 
     &      ' XFLOPA,XFLOPB,XNJOB,XNLOB', XFLOPA,XFLOPB,XNJOB,XNLOB
            WRITE(6,*) ' XNIJSX*FLOAT(NJA), XNKLSX*FLOAT(NJB) ',
     &                   XNIJSX*FLOAT(NJA), XNKLSX*FLOAT(NJB) 
            WRITE(6,*) ' ADVICE given : ', LADVICE
          END IF
        END IF
*       ^ End if several types/ph modi
      END IF
*     ^ End if ITERM test ( type of excitation)
C     WRITE(6,*) ' MEMCHECK at end of ADVICE'
C     CALL MEMCHK
C     WRITE(6,*) ' MEMCHECK passed '
      RETURN
      END
      SUBROUTINE GET_TPAM(TPAM)
*
*. Obtain ACTIVE PART of Malmqvuist transformation matrix TPAM 
* from global arrays containing complete matrix
*
*. Jeppe Olsen, Amsterdam Airport, May 30
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Output
      DIMENSION TPAM(*)
*
      NTEST = 00
*
C     CALL EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT(
C    &     WORK(KTPAM),TPAM,1)
*
      LENT = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0)
      CALL COPVEC(dbl_mb(KTPAM),TPAM,LENT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' TPAM as copied '
        CALL APRBLM2(TPAM,NACOBS,NACOBS,NSMOB,0)
      END IF
*
      RETURN
      END
c $Id$
