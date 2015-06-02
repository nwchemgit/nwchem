      SUBROUTINE TRA2Q_SIMPLE(CI,CJ,CK,CL,ICCSYM,I12SYM,X2OUT,I2OFF)
*
*
* Trivial 2-electron integral transformation with 
* separate transformation matrices for each index      
*
*. Input integrals in KINT2
*. Output integrals in X2OUT
*  If I12SYM .EQ. 1
*  integrals are packed assuming symmetry between particles 1 and 2
*  If I12SYM .EQ. 0
*  integrals are stored without assuming symmetry between particles 1 and 2
*  
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
*. Specific input
      DIMENSION CI(*),CJ(*),CK(*),CL(*)
      DIMENSION I2OFF(NSMOB,NSMOB,NSMOB)
*. Output
      DIMENSION X2OUT(*)
*
      NTEST  = 000
*
      IF(ICCSYM.EQ.0) THEN
        NOCCSYM = 1
      ELSE
        NOCCSYM = 0
      END IF 

      IF(I12SYM.EQ.0) THEN
        NO12SYM = 1
      ELSE
        NO12SYM = 0
      END IF 

      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA2Q ')
*. Largest symmetry block
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Two symmetry blocks
      LENBL = MXSOB ** 4
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' Size of symmetry block ', LENBL
      CALL  MEMMAN(KLBL1,LENBL,'ADDL ',2,'KLBL1 ')
      CALL  MEMMAN(KLBL2,LENBL,'ADDL ',2,'KLBL2 ')
*. Zero these two blocks ( to avoid the presence of Nan's)
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLBL1),ZERO,LENBL)
      CALL SETVEC(WORK(KLBL2),ZERO,LENBL)
*. Loop over symmetry blocks
      DO ISM = 1, NSMOB
        IF (ICCSYM.EQ.1) THEN
          JSM_MX = ISM
        ELSE
          JSM_MX = NSMOB
        END IF
        DO JSM = 1, JSM_MX   
          IF(I12SYM.EQ.1) THEN
           KSM_MX = ISM
          ELSE
           KSM_MX = NSMOB
          END IF
          DO KSM = 1, KSM_MX
            IF(I12SYM.EQ.1) THEN
              IF(KSM.EQ.ISM) THEN
                LSM_MX = JSM 
              ELSE
                IF (ICCSYM.EQ.1) THEN
                  LSM_MX = KSM
                ELSE
                  LSM_MX = NSMOB
                END IF
              END IF
            ELSE
              IF (ICCSYM.EQ.1) THEN
                LSM_MX = KSM
              ELSE
                LSM_MX = NSMOB
              END IF
            END IF
            DO LSM = 1, LSM_MX
              IF(NTEST.GE.100) 
     &        WRITE(6,*)  'ISM, JSM, KSM, LSM ',ISM,KSM,KSM,LSM
*. Ensure that integrals have correct symmetry
              INTSYM = 1
              IJSM = MULTD2H(ISM,JSM)
              IJKSM = MULTD2H(IJSM,KSM)
              IF(INTSYM.EQ.MULTD2H(IJKSM,LSM)) THEN
*. Fetch 2-e integral block (IJ!KL)
C?            WRITE(6,*)  'ISM, JSM, KSM, LSM ',ISM,KSM,KSM,LSM
              ONE = 1.0D0
              CALL LUCIAGETINT(WORK(KLBL1),-1,ISM,-1,JSM,-1,KSM,-1,LSM,
     &                     0,0,0,1,ONE,ONE)
*. (Type = 0  => Complete symmetryblock)
*. Offsets and dimensions for symmetryblocks in C 
              IOFF = 1
              DO IISM = 1, ISM-1
                IOFF = IOFF + NTOOBS(IISM)**2
              END DO
              NI = NTOOBS(ISM)
*
              JOFF = 1
              DO JJSM = 1, JSM-1
                JOFF = JOFF + NTOOBS(JJSM)**2
              END DO
              NJ = NTOOBS(JSM)
*
              KOFF = 1
              DO KKSM = 1, KSM-1
                KOFF = KOFF + NTOOBS(KKSM)**2
              END DO
              NK = NTOOBS(KSM)
*
              LOFF = 1
              DO LLSM = 1, LSM-1
                LOFF = LOFF + NTOOBS(LLSM)**2
              END DO
              NL = NTOOBS(LSM)
*. Transform 2-electron integral block
C?            WRITE(6,*) ' Before TRA_2EL.. '
              CALL TRA_2EL_BLK_SIMPLE(WORK(KLBL1),
     &        CI(IOFF),NI,CJ(JOFF),NJ,CK(KOFF),NK,CL(LOFF),NL,
     &        WORK(KLBL2))
              IF(NTEST.GE.100) THEN
                WRITE(6,*) ' Transformed 2e- integral block '
                CALL WRTMAT(WORK(KLBL1),NI*NJ,NK*NL,NI*NJ,NK*NL)
              END IF
*. Transfer symmetry block to integral list
             CALL PUTINT(WORK(KLBL1),0,ISM,0,JSM,0,KSM,0,LSM,
     &                    NOCCSYM,NO12SYM,X2OUT,I2OFF)
C?            WRITE(6,*) ' After PUTINT '
            END IF
*           ^ Check if integrals have correct symmetry
            END DO
          END DO
        END DO
      END DO
*     ^ End of loop over symmetries
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRA2Q ')
*
      RETURN
      END
      SUBROUTINE GET_CSIMTRH(T1,CC,CA)
*
* Obtain transformation matrices for Exp(T1) similarity 
* transformation of Hamiltonian
*
* Jeppe Olsen, July 2002
*
*Cc = 1 - T(1)t
*Ca = 1 + T(1)
*
* Matrices are outputted in complete symmetry blocks
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
*. Specific Input 
      DIMENSION T1(*)
*. Output 
      DIMENSION CC(*),CA(*)
*. Length of packed matrices
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN_C = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMST,0)
*. unit matrix in CA
      ZERO = 0.0D0
      CALL SETVEC(CA,ZERO,LEN_C)
      IB = 1
      ONE = 1.0D0
      DO ISM = 1, NSMST
        LEN = NTOOBS(ISM)
        CALL SETDIA(CA(IB),ONE,LEN,0)
        IB = IB + LEN**2
      END DO
C?      WRITE(6,*) ' Assumed unit matrix '
C?      CALL APRBLM2(CA,NTOOBS,NTOOBS,NSMST,0)
*. Transposed T1 in Cc
C          TRP_BLK_MAT(AIN,AOUT,NBLK,LROW,LCOL)
      CALL TRP_BLK_MAT(T1,CC,NSMST,NTOOBS,NTOOBS)
C?      WRITE(6,*) ' Assumed transposed T1 '
C?      CALL APRBLM2(CC,NTOOBS,NTOOBS,NSMST,0)
      ONEM = -1.0D0
*. 1 - T1(T) in Cc
      CALL VECSUM(CC,CA,CC,ONE,ONEM,LEN_C)
*. 1+ T1 in Ca
      CALL VECSUM(CA,CA,T1,ONE,ONE,LEN_C)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformation matrix for creation orbitals'
        CALL APRBLM2(CC,NTOOBS,NTOOBS,NSMST,0)
        WRITE(6,*) ' Transformation matrix for annihilation orbitals'
        CALL APRBLM2(CA,NTOOBS,NTOOBS,NSMST,0)
      END IF
*
      RETURN
      END
      SUBROUTINE TRP_BLK_MAT(AIN,AOUT,NBLK,LROW,LCOL)
*
* Transpose blocked matrix
*
* Jeppe Olsen, August 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LROW(NBLK),LCOL(NBLK)
      DIMENSION AIN(*)
*. Output
      DIMENSION AOUT(*)
*
      IOFF = 1
      DO IBLK = 1, NBLK
        CALL TRPMT3(AIN(IOFF),LROW(IBLK),LCOL(IBLK),AOUT(IOFF))
        IOFF = IOFF + LROW(IBLK)*LCOL(IBLK)
      END DO
*
      RETURN
      END
      SUBROUTINE TRP_H2_BLK(XINT,I12_OR_34,NI,NJ,NK,NL,SCR)
*
* Transform complete block of two-electron integrals 
* SCR must be large enough to hold complete integral block
*
* I12_OR_34 = 12 => transpose first two indeces 
* I12_OR_34 = 34 => transpose last  two indeces 
* I12_OR_34 = 46 => transpose first two and last two  indeces
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION XINT(NI*NJ,NK*NL)
*. Scratch
      DIMENSION SCR(*)
*
C?    WRITE(6,*) ' Memcheck at start of TRP_H2 ... '
C?    CALL MEMCHK
C?    WRITE(6,*) ' Memcheck passed '
*
      IF(I12_OR_34.EQ.12) THEN
        DO IKL = 1, NK*NL
          CALL TRPMT3(XINT(1,IKL),NI,NJ,SCR)
          CALL COPVEC(SCR,XINT(1,IKL),NI*NJ)
        END DO
      ELSE IF (I12_OR_34.EQ.34) THEN
        DO K = 1, NK
        DO L = 1, NL
         KL = (L-1)*NK + K
         LK = (K-1)*NL + L
         CALL COPVEC(XINT(1,KL),SCR(1+(LK-1)*NI*NJ),NI*NJ)
        END DO
        END DO
        CALL COPVEC(SCR,XINT,NI*NJ*NK*NL)
      ELSE IF (I12_OR_34.EQ.46) THEN
        DO K = 1, NK
        DO L = 1, NL
         KL = (L-1)*NK + K
         LK = (K-1)*NL + L
         CALL TRPMT3(XINT(1,KL),NI,NJ,SCR(1+(LK-1)*NI*NJ))
        END DO
        END DO
        CALL COPVEC(SCR,XINT,NI*NJ*NK*NL)
      ELSE 
        WRITE(6,*) ' TRP_H2_BLK :  option not programmed yet '
        STOP       ' TRP_H2_BLK :  option not programmed yet '
      END IF
*  
C?    WRITE(6,*) ' Memcheck at end of TRP_H2 ... '
C?    CALL MEMCHK
C?    WRITE(6,*) ' Memcheck passed '
      RETURN
      END
      SUBROUTINE EXPAND_T1(T,T1EXP,ITSM,IAB)
*
* expand T1 coefficients in complete T array to complete matrix form
*
* Jeppe Olsen, Korshojen, August 2000
* IAB added to argument list July 2002 to allow to access 
* alpha and beta parts 
*
* IAB = 1 => alpha part of one-electron excitations obtained 
* IAB = 2 => beta  part of one-electron excitations obtained
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'csm.inc'
*. Input
      DIMENSION T(*)
*. Output
      DIMENSION T1EXP(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'T1EXP ')
*. Space for creation and annihilation in 2 sets of spinorbital excitations
      CALL MEMMAN(KLT1C,NTOOB,'ADDL  ',2,'T1C   ')
      CALL MEMMAN(KLT1A,NTOOB,'ADDL  ',2,'T1A   ')
      CALL MEMMAN(KLT1MAT,NTOOB**2,'ADDL  ',2,'T1MAT')
*. Expand to complete matrix ( no symmetrypacking) 
C?    WRITE(6,*) ' First element of T = ',T(1)
      CALL EXPAND_T1S(T,WORK(KLT1MAT),ITSM,WORK(KLT1C),WORK(KLT1A),
     &                WORK(KLSOBEX),WORK(KLSOX_TO_OX),
     &                WORK(KLCOBEX_TP),WORK(KLIBSOBEX),
     &                NSPOBEX_TP,IAB)
*. Pack matrix to symmetrypacked form 
C          REORHO1(RHO1I,RHO1O,IRHO1SM)
*. Jeppe notes: The above is not working of there are inactive orbitals:
*  REORHO1 only reforms active orbitals
      CALL REORHO1(WORK(KLT1MAT),T1EXP,ITSM,1)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T1 matrix in symmetrypacked form for IAB =', IAB
        CALL PRHONE(T1EXP,NTOOBS,ITSM,NSMST,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'T1EXP ')
*
      RETURN
      END
      SUBROUTINE EXPAND_T1S(T,T1,ITSM,IC_OCC,IA_OCC,
     &                ISPOBEX_TP,ISOX_TO_OX,IEXC_FOR_OX,
     &                IBSPOBEX_TP,NSPOBEX_TP,IAB)
*
* Slave routine for expanding t1 coefficients to complete matrix
*
* Jeppe Olsen, Korshojen, august 2000
* IAB transferred to argument list, July 2002 to allow for OS
*
* IAB = 1 => alpha excitations are collected 
* IAB = 2 => beta excitations are collected 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cc_exc.inc'
*. General input 
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER ISOX_TO_OX(NSPOBEX_TP)
      INTEGER IEXC_FOR_OX(*)
      INTEGER IBSPOBEX_TP(NSPOBEX_TP)
*. Specific input
      DIMENSION T(*)
*. Output
      DIMENSION T1(NTOOB,NTOOB)
*. Scratch through input
      INTEGER IC_OCC(NTOOB),IA_OCC(NTOOB)
*. Local scratch
      INTEGER IC_GRP(MXPNGAS),IA_GRP(MXPNGAS)
*
      IF(MSCOMB_CC.EQ.1) THEN
*. Combinations, 
        FACTOR = 1.0D0/SQRT(2.0D0)
      ELSE
*. Individual terms
        FACTOR = 1.0D0
      END IF
*
      ZERO = 0.0D0
      CALL SETVEC(T1,ZERO,NTOOB**2)
*
      DO ISOXTP = 1, NSPOBEX_TP
*
      IF(IEXC_FOR_OX(ISOX_TO_OX(ISOXTP)).EQ.1) THEN
*. Is single excitation alpha or beta excitation ?
        NCA_OP = IELSUM(ISPOBEX_TP(1+0*NGAS,ISOXTP),NGAS)
        NCB_OP = IELSUM(ISPOBEX_TP(1+1*NGAS,ISOXTP),NGAS)
        NAA_OP = IELSUM(ISPOBEX_TP(1+2*NGAS,ISOXTP),NGAS)
        NAB_OP = IELSUM(ISPOBEX_TP(1+3*NGAS,ISOXTP),NGAS)
        IELMNT = IBSPOBEX_TP(ISOXTP)
        IF(NCA_OP.EQ.1.AND.IAB.EQ.1.OR.
     &     NCB_OP.EQ.1.AND.IAB.EQ.2     )THEN
          IF(NCA_OP.EQ.1.AND.IAB.EQ.1) THEN
*. Transform from occupations to groups
            CALL OCC_TO_GRP(ISPOBEX_TP(1+0*NGAS,ISOXTP),IC_GRP,1)
            CALL OCC_TO_GRP(ISPOBEX_TP(1+2*NGAS,ISOXTP),IA_GRP,1)
            NC_OP = NCA_OP
            NA_OP = NAA_OP
          ELSE IF(NCB_OP.EQ.1.AND.IAB.EQ.2) THEN
*. Transform from occupations to groups
            CALL OCC_TO_GRP(ISPOBEX_TP(1+1*NGAS,ISOXTP),IC_GRP,1)
            CALL OCC_TO_GRP(ISPOBEX_TP(1+3*NGAS,ISOXTP),IA_GRP,1)
            NC_OP = NCB_OP
            NA_OP = NAB_OP
          END IF
*
          DO ISM_C = 1, NSMST
            ISM_A = MULTD2H(ISM_C,ITSM)
*. Obtain occupations
            IDUM = 0
            CALL GETSTR2_TOTSM_SPGP(IC_GRP,NGAS,ISM_C,NC_OP,NSTR_C,
     &           IC_OCC, NTOOB,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IA_GRP,NGAS,ISM_A,NA_OP,NSTR_A,
     &           IA_OCC, NTOOB,0,IDUM,IDUM)
*. And scatter to complete matrix 
C?           WRITE(6,*) ' ISM_C, ISM_A, NSTR_A, NSTR_C =', 
C?   &                    ISM_C, ISM_A, NSTR_A, NSTR_C 
             DO ISA = 1, NSTR_A
               DO ISC = 1, NSTR_C
                 IORB_A = IA_OCC(ISA)
                 IORB_C = IC_OCC(ISC)
                 T1(IORB_C,IORB_A) = T(IELMNT)*FACTOR
C?               WRITE(6,*) ' IORB_C, IORB_A, IELMNT =',
C?   &                        IORB_C, IORB_A, IELMNT
C?               WRITE(6,*) ' T(IELMNT) = ', T(IELMNT)
                 IELMNT = IELMNT + 1
               END DO
             END DO
          END DO
        END IF
*       ^ End of correct single excitation
      END IF
*     ^ End if single excitation
      END DO
*     ^ End of loop over excitation types
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T1 matrix in complete matrix form '
        CALL WRTMAT(T1,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
      SUBROUTINE TRA_SIMTRH(T,IO_OR_SO_TRA)
*
* Obtain integrals defining similarity transformed Hamiltonian
*
*     Exp(T1) H Exp(-T1)
*
* Jeppe Olsen, Korshojen 53, August 2000
*              Updated July 2002, spin-orbital transformation added
*
*
* IO_OR_SO_TRSA 
* = 1 => orbital transformation assuming that alpha and beta parts of T1
*        are identical
* = 2 => spin-orbital transformation allowing differences in T1(alpha)
*        and T1(beta)
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
*. Specific input
      DIMENSION T(*) 
 
*
      IDUM = 0
      NTEST = 00
      CALL LUCIAQENTER('SIMTR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ET1HE ')
*
      IF(NTEST.GT.0) THEN
      IF(IO_OR_SO_TRA.EQ.1) THEN
        WRITE(6,*) ' Orbital based similarity transformation '
      ELSE 
        WRITE(6,*) ' Spinorbital based similarity transformation '
      END IF
      END IF
*
* 1 : Obtain T1 in matrix-form
*
*. Memory for expanded matrix
      LEN = NTOOB**2
      IF(IO_OR_SO_TRA.EQ.1) THEN
        CALL MEMMAN(KLT1EXP,LEN,'ADDL  ',2,'T1EXP ')
        CALL EXPAND_T1(T,WORK(KLT1EXP),1,1)
      ELSE 
        CALL MEMMAN(KLT1EXPA,LEN,'ADDL  ',2,'T1EXPA')
        CALL MEMMAN(KLT1EXPB,LEN,'ADDL  ',2,'T1EXPB')
        CALL EXPAND_T1(T,WORK(KLT1EXPA),1,1)
        CALL EXPAND_T1(T,WORK(KLT1EXPB),1,2)
      END IF
*
* 2 : Obtain transformation matrices for orbitals in integrals
*
*Cc = 1 - T(1)t
*Ca = 1 + T(1)
*
      IF(IO_OR_SO_TRA.EQ.1) THEN
        CALL MEMMAN(KLCC,LEN,'ADDL  ',2,'Ccrea ')
        CALL MEMMAN(KLCA,LEN,'ADDL  ',2,'Canni ')
        CALL GET_CSIMTRH(WORK(KLT1EXP),WORK(KLCC),WORK(KLCA))
C            GET_CSIMTRH(T1,CC,CA)
      ELSE   
        CALL MEMMAN(KLCCA,LEN,'ADDL  ',2,'CcreaA')
        CALL MEMMAN(KLCAA,LEN,'ADDL  ',2,'CanniA')
        CALL MEMMAN(KLCCB,LEN,'ADDL  ',2,'CcreaB')
        CALL MEMMAN(KLCAB,LEN,'ADDL  ',2,'CanniB')
        CALL GET_CSIMTRH(WORK(KLT1EXPA),WORK(KLCCA),WORK(KLCAA))
        CALL GET_CSIMTRH(WORK(KLT1EXPB),WORK(KLCCB),WORK(KLCAB))
      END IF
*
* Transform one-electron integrals
*
*. Expand one-electron integrals to complete matrix form
      CALL MEMMAN(KLH1EXP,LEN,'ADDL  ',2,'H1EXP ')
      IF(IO_OR_SO_TRA.EQ.1) THEN
        CALL TRIPAK_BLKM(WORK(KLH1EXP),WORK(KINT1O),2,NTOOBS,NSMST)
C             TRA1D_SIMPLE(CCREA,CANNI,HIN,HOUT,ISYM)
        CALL TRA1D_SIMPLE(WORK(KLCC),WORK(KLCA),WORK(KLH1EXP),
     &                    WORK(KINT1_SIMTRH),0)
      ELSE 
        CALL TRIPAK_BLKM(WORK(KLH1EXP),WORK(KINT1O),2,NTOOBS,NSMST)
        CALL TRA1D_SIMPLE(WORK(KLCCA),WORK(KLCAA),WORK(KLH1EXP),
     &                    WORK(KINT1_SIMTRH_A),0)
        CALL TRIPAK_BLKM(WORK(KLH1EXP),WORK(KINT1O),2,NTOOBS,NSMST)
        CALL TRA1D_SIMPLE(WORK(KLCCB),WORK(KLCAB),WORK(KLH1EXP),
     &                    WORK(KINT1_SIMTRH_B),0)
      END IF
*
* Transform two-electron integrals
*
      IF(IO_OR_SO_TRA.EQ.1) THEN
        CALL TRA2D_SIMPLE(WORK(KLCC),WORK(KLCA),WORK(KINT2_SIMTRH))
      ELSE 
*.Alpha-alpha integrals
C            TRA2Q_SIMPLE(CI,CJ,CK,CL,I12SYM,X2OUT)
        CALL TRA2Q_SIMPLE(WORK(KLCCA),WORK(KLCAA),WORK(KLCCA),
     &       WORK(KLCAA),0,1,WORK(KINT2_SIMTRH_AA),WORK(KPINT2_SIMTRH))
*.beta-beta integrals
        CALL TRA2Q_SIMPLE(WORK(KLCCB),WORK(KLCAB),WORK(KLCCB),
     &       WORK(KLCAB),0,1,WORK(KINT2_SIMTRH_BB),WORK(KPINT2_SIMTRH))
*alpa-beta integrals
        CALL TRA2Q_SIMPLE(WORK(KLCCA),WORK(KLCAA),WORK(KLCCB),
     &       WORK(KLCAB),0,0,WORK(KINT2_SIMTRH_AB),
     &       WORK(KPINT2_SIMTRH_AB))
      END IF
*
* hole contributions to core energy and one-electron integrals
* 
      IF(IUSE_PH.EQ.1) THEN
        I_USE_SIMTRH = 1
        IF(IO_OR_SO_TRA.EQ.1) THEN
          CALL FI(WORK(KINT1_SIMTRH),ECORE_HEX,1)
        ELSE
          I_UNRORB = 1
          CALL FI_HS_SM_AB(WORK(KINT1_SIMTRH_A),WORK(KINT1_SIMTRH_B),
     &                  WORK(KFI_AL),WORK(KFI_BE),ECORE_HEX,1) 
          I_UNRORB = 0
        END IF
        I_USE_SIMTRH = 0
      ELSE
         ECORE_HEX = 0.0D0
      END IF
C?    WRITE(6,*) ' ECORE_ORIG, ECORE_HEX =', ECORE_ORIG, ECORE_HEX
      ECORE = ECORE_ORIG + ECORE_HEX
      IF(NTEST.GT.0) 
     &WRITE(6,*) ' TRA_SIMTRH : Updated core energy ',ECORE
     
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ET1HE ')
      CALL LUCIAQEXIT('SIMTR')
      RETURN
      END
      SUBROUTINE TRAINT
*
* Very simple integral transformation for LUCIA
* - for testing the real one
*
* MO => MO transformation defined by matrix C in KMOMO
*
* Jeppe Olsen, One day in November 98, Magistratsvaegen 37D
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cprnt.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      NTEST = 00
      IF(IPRINTEGRAL.GE.100) NTEST = 100
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACTL')
      CALL LUCIAQENTER('TRAIN')
*
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' TRAINT Entered '
       WRITE(6,*) ' ==============='
       IF(ITRA_ROUTE.EQ.1) THEN
         WRITE(6,*) ' Oldfashioned integral trans'
       ELSE
         WRITE(6,*) ' New integrals trans in action'
       END IF
      END IF
*
*
      IF(NTEST.GE.100.AND.ITRA_ROUTE.EQ.1) THEN
        WRITE(6,*) ' Transformation matrix C '
        CALL APRBLM2(dbl_mb(KMOMO),NTOOBS,NTOOBS,NSMOB,0)
C            APRBLM2(A,LROW,LCOL,NBLK,ISYM)
      END IF

*
*. Transform 1-electron integrals : 
*
*     KINT1O => KINT1
      CALL TRA1_SIMPLE(dbl_mb(KKCMO_I),dbl_mb(KKCMO_J))
*
*. Transform 2-electron integrals
*
      IF(ITRA_ROUTE.EQ.1) THEN
*  Integrals in KINT2 => Integrals in KINT2
        CALL TRA2_SIMPLE(dbl_mb(KMOMO),KINT2)
      ELSE
*  Integrals in KINT2 => Integrals in arrays defined by active integral list
*  permutational symmetry of input integrals is defined by I12S, I34S, I1234S
*  If a single array is active, the final integrals are stored in 
*  WORK(KINT2_A(IIE2ARR))
        CALL GEN_TRA_2EI_LIST
      END IF
*
       IDMPINT_LOC = 0
       IF (IDMPINT_LOC.EQ.1 ) THEN
        WRITE(6,*)
     &   ' Integrals written formatted (E22.15) on unit 90'
        LU90 = 90
        REWIND LU90   
*.1 : One-electron integrals
CNW     WRITE(LU90,'(E22.15)')
CNW  &   (WORK(KINT1-1+INT1),INT1=1,NINT1)
        call ga_print(KINT1)
*.2 : Two-electron integrals
CNW     WRITE(LU90,'(E22.15)')
CNW  &   (WORK(KINT2-1+INT2),INT2=1,NINT2)
        call ga_print(KINT2)
*.3. Core energy 
        WRITE(LU90,'(E22.15)')ECORE_ORIG
*.4  Rewind to empty buffer
        REWIND LU90
*.   Symmetry info etc to LU91
        LU91 = 91
        CALL DUMP_1EL_INFO(LU91)
       END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACTL')
      CALL LUCIAQEXIT ('TRAIN')
      RETURN
      END
      SUBROUTINE TRA1D_SIMPLE(CCREA,CANNI,HIN,HOUT,ISYM)
*
* Transform one-electron integrals with dual transformation matrix 
* CCREA, CANNI 
*
* H' = Ccrea(T) H Canni
*
* Input integrals in HIN
* Output integrals also in HOUT
*
* If ISYM = 1, then input and output integrals are lower half packed
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
*. Specific Input
      DIMENSION CCREA(*),CANNI(*)
*
      NTEST = 0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA1_S')
*. Largest symmetry block of orbitals 
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Largest number of orbitals in sym block ', MXSOB
      END IF
*. Space for 1e- integrals and two sym blocks
c weg??
      LENH1 = NTOOB ** 2
      CALL MEMMAN(KLH1,LENH1,'ADDL  ',2,'LH1   ')
      LENSCR = 2 * MXSOB ** 2
      CALL MEMMAN(KLSCR,LENSCR,'ADDL   ',2,'H1SCR ')
C     WRITE(6,*) ' LENH1 LENSCR ', LENH1, LENSCR
*. and do it
     
C?    WRITE(6,*) ' Before TRAN_SYM '
      CALL TRAN_SYM_BLOC_MAT4(HIN,CCREA,CANNI,
     &     NSMOB,NTOOBS,NTOOBS,HOUT,WORK(KLSCR),ISYM) 
C     TRAN_SYM_BLOC_MAT4
C    &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
C     WRITE(6,*) ' After TRAN_SYM '
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformed 1e- integrals '
        WRITE(6,*) ' ========================='
C            WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
        CALL WRTVH1(HOUT,1,NTOOBS,NTOOBS,NSMOB,ISYM)
      END IF
*. Flush memory
C?    WRITE(6,*) ' Returning from TRA1 '
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'TRA1_S')
*
      RETURN
      END
      SUBROUTINE TRA1_SIMPLE(CI,CJ)
*
* Transform one-electron integrals
*
* KINT1O => KINT1
*
* i.e. :
*
* Input integrals in KINT1O 
* Output integrals in KINT1
* MOAO transformation matrices in CI, CJ
* KINT1O and KINT1 may be identical..
*
* Jeppe Olsen
* Modified, July 2011 to allow transformation with different orbital indices
*           Note that input integrals are assumed to be packed
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General Input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'cprnt.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific Input
      DIMENSION CI(*), CJ(*)
*
      NTEST = 000
      NTEST = MAX(IPRINTEGRAL,NTEST)
      IPRNT_TRAINT = 0
      IF(IPRINTEGRAL.GE.1000) IPRNT_TRAINT = 1
      
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA1_S')
*. Largest symmetry block of orbitals 
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Info from TRA1_SIMPLE '
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        WRITE(6,*) ' Largest number of orbitals in sym block ', MXSOB
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Initial set of integrals in KINT1O'
CNW       CALL APRBLM2(WORK(KINT1O),NTOOBS,NTOOBS,NSMOB,1)
          call ga_print(KINT1O)
          WRITE(6,*) ' CI transformation matrix '
          CALL APRBLM2(CI, NTOOBS,NTOOBS,NSMOB,0)
          WRITE(6,*) ' CJ transformation matrix '
          CALL APRBLM2(CJ, NTOOBS,NTOOBS,NSMOB,0)
        END IF
        IF(IH1FORM.EQ.1) THEN
         WRITE(6,*) ' Output symmetryblocks will be packed'
        ELSE
         WRITE(6,*) ' Output symmetryblocks will not be packed'
       END IF
      END IF
*
*. Space for 1e- integrals and two sym blocks
      LENH1 = NTOOB ** 2
      IF(IH1FORM.EQ.1) THEN
        IPACK_OUT = 1
      ELSE
        IPACK_OUT = 0
      END IF
*
      LEN_H = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,IPACK_OUT)
      LEN_HP= NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,1)
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      CALL MEMMAN(KLH1,LEN_H,'ADDL  ',2,'LH1   ')
      LENSCR = 2 * MXSOB ** 2
      CALL MEMMAN(KLSCR,LENSCR,'ADDL   ',2,'H1SCR ')
*. Unpack input integrals if output integrals will be unpacked
      IF(IPACK_OUT.EQ.0) THEN
        call errquit('tra1_simple:should not be here')
        CALL TRIPAK_AO_MAT(WORK(KLH1),WORK(KINT1O),2)
      ELSE
CNW     CALL COPVEC(WORK(KINT1O),WORK(KLH1),LEN_HP)
        if (ga_nodeid().eq.0) 
     &     call ga_get(KINT1O,1,len_hp,1,1,dbl_mb(klh1),1)
      END IF
*. and do it
      if (ga_nodeid().eq.0) 
     &CALL TRAN_SYM_BLOC_MAT4(dbl_mb(KLH1),CI,CJ,NSMOB,NTOOBS,NTOOBS,
     &     KINT1,dbl_mb(KLSCR),IPACK_OUT)
      call ga_sync()
C     TRAN_SYM_BLOC_MAT4(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
C?    WRITE(6,*) ' After TRAN_SYM '
*
      IF(IPRNT_TRAINT.EQ.1.OR.NTEST.GE.100) THEN
       if (ga_nodeid().eq.0) WRITE(6,*) ' Transformed 1e integrals '
CNW    CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,IPACK_OUT)
       call ga_print(KINT1)
      END IF
*
*. Flush memory
      IF(NTEST.GE.1000) WRITE(6,*) ' Returning from TRA1'
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'TRA1_S')
*
      RETURN
      END
      SUBROUTINE TRA2D_SIMPLE(CC,CA,X2OUT)
*
* Trivial 2-electron integral transformation with dual 
* transformation matrices CCREA, CANNI
*
*. Input integrals in KINT2
*. Output integrals in X2OUT
*. Integrals are packed assuming only symmetry between particles 1 and 2
*. symmetry
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
*. Specific input
      DIMENSION CC(*),CA(*)
*. Output
      DIMENSION X2OUT(*)
*
      NTEST  = 00
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA2_S')
*. Largest symmetry block
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Two symmetry blocks
      LENBL = MXSOB ** 4
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' Size of symmetry block ', LENBL
      CALL  MEMMAN(KLBL1,LENBL,'ADDL ',2,'KLBL1 ')
      CALL  MEMMAN(KLBL2,LENBL,'ADDL ',2,'KLBL2 ')
*. Zero these two blocks ( to avoid the presence of Nan's)
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLBL1),ZERO,LENBL)
      CALL SETVEC(WORK(KLBL2),ZERO,LENBL)
C     WRITE(6,*) ' Two integral blocks allocated '
*. Loop over symmetry blocks
      DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB   
          DO KSM = 1, ISM       
            IF(KSM.EQ.ISM) THEN
              LSM_MX = JSM 
            ELSE
              LSM_MX = NSMOB
            END IF
            DO LSM = 1, LSM_MX
              IF(NTEST.GE.100) 
     &        WRITE(6,*)  'ISM, JSM, KSM, LSM ',ISM,KSM,KSM,LSM
*. Ensure that integrals have correct symmetry
              INTSYM = 1
              IJSM = MULTD2H(ISM,JSM)
              IJKSM = MULTD2H(IJSM,KSM)
              IF(INTSYM.EQ.MULTD2H(IJKSM,LSM)) THEN
*. Fetch 2-e integral block (IJ!KL)
              ONE = 1.0D0
C?            WRITE(6,*)  'ISM, JSM, KSM, LSM ',ISM,KSM,KSM,LSM
              CALL LUCIAGETINT(WORK(KLBL1),-1,ISM,-1,JSM,-1,KSM,-1,LSM,
     &                     0,0,0,1,ONE,ONE)
*. (Type = 0  => Complete symmetryblock)
*. Offsets and dimensions for symmetryblocks in C 
              IOFF = 1
              DO IISM = 1, ISM-1
                IOFF = IOFF + NTOOBS(IISM)**2
              END DO
              NI = NTOOBS(ISM)
*
              JOFF = 1
              DO JJSM = 1, JSM-1
                JOFF = JOFF + NTOOBS(JJSM)**2
              END DO
              NJ = NTOOBS(JSM)
*
              KOFF = 1
              DO KKSM = 1, KSM-1
                KOFF = KOFF + NTOOBS(KKSM)**2
              END DO
              NK = NTOOBS(KSM)
*
              LOFF = 1
              DO LLSM = 1, LSM-1
                LOFF = LOFF + NTOOBS(LLSM)**2
              END DO
              NL = NTOOBS(LSM)
*. Transform 2-electron integral block
C?            WRITE(6,*) ' Before TRA_2EL.. '
              CALL TRA_2EL_BLK_SIMPLE(WORK(KLBL1),
     &        CC(IOFF),NI,CA(JOFF),NJ,CC(KOFF),NK,CA(LOFF),NL,
     &        WORK(KLBL2))
              IF(NTEST.GE.100) THEN
                WRITE(6,*) ' Transformed 2e- integral block '
                CALL WRTMAT(WORK(KLBL1),NI*NJ,NK*NL,NI*NJ,NK*NL)
              END IF
*. Transfer symmetry block to integral list
              CALL PUTINT(WORK(KLBL1),0,ISM,0,JSM,0,KSM,0,LSM,1,0,
     &                    X2OUT,WORK(KPINT2_SIMTRH))
C                  PUTINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                  NOCCSYM,NO12SYM,XOUT)
C?            WRITE(6,*) ' After PUTINT '
            END IF
*           ^ Check if integrals have correct symmetry
            END DO
          END DO
        END DO
      END DO
*     ^ End of loop over symmetries
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRA2_S')
*
      RETURN
      END
      SUBROUTINE TRA_2EL_BLK_SIMPLE(XINT,CI,NI,CJ,NJ,CK,NK,CL,NL,SCR)
*
* Transform 2-electron integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CI(NI,NI),CJ(NJ,NJ),CK(NK,NK),CL(NL,NL)
*. Input and output
      DIMENSION XINT(*) 
* Matrix given in complete form XI(NI,NJ,NK,NL)
*. Scratch
      DIMENSION SCR(*)
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' TRA_2EL_BLK_SIMPLE: Input integral block as X(IJ,KL)'
        CALL WRTMAT(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*. Transform first two indeces
      DO K = 1, NK
       DO L = 1, NL
         KLOF = 1 + ((L-1)*NK+K-1)*NI*NJ
         FACTORC = 0.0D0
         FACTORAB = 1.0D0
         CALL MATML7(SCR(KLOF),XINT(KLOF),CJ,NI,NJ,NI,NJ,NJ,NJ,
     &                FACTORC,FACTORAB,0)
         CALL MATML7(XINT(KLOF),CI,SCR(KLOF),NI,NJ,NI,NI,NI,NJ,
     &                FACTORC,FACTORAB,1)
       END DO
      END DO
*. Transpose Block X(IJ,KL)
      NIJ = NI*NJ
      NKL = NK*NL
      CALL TRPMT3(XINT,NIJ,NKL,SCR)
*. Matrix is now SCR(KL,IJ)
*. Transform last two indeces
      DO I = 1, NI
        DO J = 1, NJ
         IJOF = ((J-1)*NI+I-1)*NK*NL+1
         CALL MATML7(XINT(IJOF),SCR(IJOF),CL,NK,NL,NK,NL,NL,NL,
     &               FACTORC,FACTORAB,0)
         CALL MATML7(SCR(IJOF),CK,XINT(IJOF),NK,NL,NK,NK,NK,NL,
     &               FACTORC,FACTORAB,1)
        END DO
      END DO
*. Transpose to obtain XINT(IJ,KL)
      CALL TRPMT3(SCR,NKL,NIJ,XINT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output integral block as X(IJ,KL)'
        CALL WRTMAT(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END
      SUBROUTINE PUTINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  NOCCSYM,NO12SYM,XOUT,I2OFF)
*
* Put integrals in permanent integral list XOUT 
*
* If NOCCSYM = 1, integrals do not possess any
* complex conjugation symmetry
* IF NO12SYM = 1, integrals do not possess symmetry 
* between particle 1 and 2 
*
* Jeppe Olsen, Jan. 1999
*              July 2002 : NO12SYM added for general spinorbital transf
*
c      IMPLICIT REAL*8(A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      DIMENSION XINT(*)
*. Output
      integer XOUT
CNW   DIMENSION XOUT(*)
*
      CALL LUCIAQENTER('PUTIN')
*. Offset and number of integrals
*
      I_AM_KIM = 0
      IF(I_AM_KIM.EQ.1) THEN
        WRITE(6,*) ' Integrals will be dumped in formatted form on LU74'
        WRITE(6,*) ' Integrals will be dumped in formatted form on LU74'
        WRITE(6,*) ' Integrals will be dumped in formatted form on LU74'
        WRITE(6,*) ' Integrals will be dumped in formatted form on LU74'
        WRITE(6,*) ' Integrals will be dumped in formatted form on LU74'
        LU74 = 74
        CALL REWINO(LU74)
*. But first write overlap and one-electron integrals in AO basis 
*. on file 75
        LUSH = 75
        CALL READ_WRITE_SH(LUSH)
      END IF
*
      IF(ITP.EQ.0) THEN
        NI = NTOOBS(ISM)
      ELSE
        NI = NOBPTS(ITP,ISM)
      END IF
*
      IOFF = IBSO(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NOBPTS(IITP,ISM)
      END DO
*
      IF(JTP.EQ.0) THEN
        NJ = NTOOBS(JSM)
      ELSE
        NJ = NOBPTS(JTP,JSM)
      END IF
*
      JOFF = IBSO(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NOBPTS(JJTP,JSM)
      END DO
*
      IF(KTP.EQ.0) THEN
        NK = NTOOBS(KSM)
      ELSE
        NK = NOBPTS(KTP,KSM)
      END IF
*
      KOFF = IBSO(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NOBPTS(KKTP,KSM)
      END DO
*
      IF(LTP.EQ.0) THEN
        NL = NTOOBS(LSM)
      ELSE
        NL = NOBPTS(LTP,LSM)
      END IF
*
      LOFF = IBSO(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NOBPTS(LLTP,LSM)
      END DO
*
      INT_IN = 0
      DO LOB = LOFF,LOFF+NL-1
       DO KOB = KOFF,KOFF+NK-1
        DO JOB = JOFF,JOFF+NJ-1
         DO IOB = IOFF,IOFF+NI-1
C?         WRITE(6,*) ' IOB, JOB, KOB, LOB', IOB,JOB,KOB,LOB

c old:
c           IF(NOCCSYM.EQ.0) THEN
c             INT_OUT = I2EAD(IOB,JOB,KOB,LOB)
c           ELSE 
c             INT_OUT = I2EAD_NOCCSYM(IOB,JOB,KOB,LOB,NO12SYM)
c           END IF
c new:
           INT_OUT = I2ADDR(IOB,JOB,KOB,LOB,I2OFF,NOCCSYM,NO12SYM)

           INT_IN = INT_IN + 1
*
           IF(INT_OUT.EQ.0) THEN
             WRITE(6,*) ' INT_OUT, INT_IN ', INT_OUT,INT_IN
             WRITE(6,*) ' IOB JOB, KOB, LOB = ', 
     &                    IOB,JOB, KOB, LOB
             WRITE(6,*) ' Stop vanishing 2e address '
             STOP       ' Stop vanishing 2e address '
           END IF
*. Dump in formatted form for Kim Vogels VB program
           IF(I_AM_KIM.EQ.1) THEN
             WRITE(LU74,'(4I5,F25.15)') IOB,JOB,KOB,LOB,XINT(INT_IN)
           END IF
CNW        XOUT(INT_OUT) = XINT(INT_IN) 
           call ga_put(XOUT,INT_OUT,INT_OUT,1,1,XINT(INT_IN),1)
C?         WRITE(6,*) ' IOB, JOB, KOB, LOB', IOB,JOB,KOB,LOB
         END DO
        END DO
       END DO
      END DO
*
      CALL LUCIAQEXIT('PUTIN')
      RETURN
      END
      SUBROUTINE DMPINT(LUINT)
*
* Dump integrals in WORK(KINT1),WORK(KINT2) on file LUINT
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*
      INCLUDE 'cintfo.inc'
      INCLUDE 'cecore.inc'
*
      CALL REWINO(LUINT)
*.1 : One-electron integrals
      WRITE(LUINT,'(E22.15)')
     &     (WORK(KINT1-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
      WRITE(LUINT,'(E22.15)')
     &     (WORK(KINT2-1+INT2),INT2=1,NINT2)
*.3 : Core energy
      WRITE(LUINT,'(E22.15)') ECORE_ORIG
*
      RETURN
      END
      SUBROUTINE TRAN_SYM_BLOC_MAT4
     &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrices XL, XR
*  X to yield blocked matrix AOUT
*
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = XL(transposed) A XR
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION AIN(*),LX_ROW(NBLOCK),LX_COL(NBLOCK)
      DIMENSION XL(*),XR(*)
*. Output 
CNW   DIMENSION AOUT(*)
      integer AOUT
*. Scratch : At least twice the length of largest block 
      DIMENSION SCR(*)
*
*. To get rid of annoying and incorrect compiler warnings
      IOFFP_IN = 0
      IOFFC_IN = 0
      IOFFP_OUT = 0 
      IOFFC_OUT = 0
      IOFFX = 0
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP_IN = 1
         IOFFC_IN = 1
         IOFFP_OUT= 1
         IOFFC_OUT= 1
         IOFFX = 1
       ELSE
         IOFFP_IN = IOFFP_IN + LX_ROW(IBLOCK-1)*(LX_ROW(IBLOCK-1)+1)/2
         IOFFC_IN = IOFFC_IN + LX_ROW(IBLOCK-1) ** 2                     
         IOFFP_OUT= IOFFP_OUT+ LX_COL(IBLOCK-1)*(LX_COL(IBLOCK-1)+1)/2
         IOFFC_OUT= IOFFC_OUT+ LX_COL(IBLOCK-1) ** 2                     
         IOFFX = IOFFX + LX_ROW(IBLOCK-1)*LX_COL(IBLOCK-1)
       END IF
       LXR = LX_ROW(IBLOCK)
       LXC = LX_COL(IBLOCK)
       K1 = 1
       K2 = 1 + MAX(LXR,LXC) ** 2
*. Unpack block of A
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK(SCR(K1),AIN(IOFFP_IN),2,LXR,LXR)
C             TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
       ELSE
         CALL COPVEC(AIN(IOFFC_IN),SCR(K1),LXR*LXR)
       END IF
*. XL(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL SETVEC(SCR(K2),ZERO,LXR*LXC)
C?     WRITE(6,*) ' TEST : K2, IOFFX, K1, LXC, LXR ',
C?   &                     K2, IOFFX, K1, LXC, LXR
       CALL MATML7(SCR(K2),XL(IOFFX),SCR(K1),LXC,LXR,LXR,LXC,LXR,LXR,
     &             ZERO,ONE,1)
C?     WRITE(6,*) ' Half transformed matrix '
C?     CALL WRTMAT(SCR(K2),LXC,LXR,LXC,LXR)

*. XL(T) (IBLOCK) A(IBLOCK) XR(IBLOCK)
       CALL SETVEC(SCR(K1),ZERO,LXC*LXC)
       CALL MATML7(SCR(K1),SCR(K2),XR(IOFFX),LXC,LXC,LXC,LXR,LXR,LXC,
     &             ZERO,ONE,0)
C?     WRITE(6,*) ' Transformed matrix '
C?     CALL WRTMAT(SCR(K1),LXC,LXC,LXC,LXC)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
CNW      CALL TRIPAK(SCR(K1),AOUT(IOFFP_OUT),1,LXC,LXC)
         IJ = 0
         IJF = 0
         DO I = 1,LXC 
           DO J = 1, I
             SCR(K2+IJ+J-1) = SCR(K1+IJF+J-1)
           END DO
           IJ = IJ + I
           IJF = IJF + LXC
         END DO
         IBLK = LXC*(LXC+1)/2
         if (IBLK.gt.0)
     &   call ga_put(AOUT,IOFFP_OUT,IOFFP_OUT+IBLK-1,1,1,
     &               SCR(K2),IBLK)
       ELSE
CNW      CALL COPVEC(SCR(K1),AOUT(IOFFC_OUT),LXC*LXC)
         if (LXC.gt.0)
     &   call ga_put(AOUT,IOFFC_OUT,IOFFC_OUT+LXC*LXC-1,1,1,
     &               SCR(K1),LXC*LXC)
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
CNW     CALL APRBLM2(AOUT,LX_COL,LX_COL,NBLOCK,ISYM)      
        call ga_print(AOUT)
      END IF
*
      RETURN
      END
      SUBROUTINE TRAN_SYM_BLOC_MAT4a
     &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrices XL, XR
*  X to yield blocked matrix AOUT
*
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = XL(transposed) A XR
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION AIN(*),LX_ROW(NBLOCK),LX_COL(NBLOCK)
      DIMENSION XL(*),XR(*)
*. Output 
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block 
      DIMENSION SCR(*)
*
*. To get rid of annoying and incorrect compiler warnings
      IOFFP_IN = 0
      IOFFC_IN = 0
      IOFFP_OUT = 0 
      IOFFC_OUT = 0
      IOFFX = 0
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP_IN = 1
         IOFFC_IN = 1
         IOFFP_OUT= 1
         IOFFC_OUT= 1
         IOFFX = 1
       ELSE
         IOFFP_IN = IOFFP_IN + LX_ROW(IBLOCK-1)*(LX_ROW(IBLOCK-1)+1)/2
         IOFFC_IN = IOFFC_IN + LX_ROW(IBLOCK-1) ** 2                     
         IOFFP_OUT= IOFFP_OUT+ LX_COL(IBLOCK-1)*(LX_COL(IBLOCK-1)+1)/2
         IOFFC_OUT= IOFFC_OUT+ LX_COL(IBLOCK-1) ** 2                     
         IOFFX = IOFFX + LX_ROW(IBLOCK-1)*LX_COL(IBLOCK-1)
       END IF
       LXR = LX_ROW(IBLOCK)
       LXC = LX_COL(IBLOCK)
       K1 = 1
       K2 = 1 + MAX(LXR,LXC) ** 2
*. Unpack block of A
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK(SCR(K1),AIN(IOFFP_IN),2,LXR,LXR)
C             TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
       ELSE
         CALL COPVEC(AIN(IOFFC_IN),SCR(K1),LXR*LXR)
       END IF
*. XL(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL SETVEC(SCR(K2),ZERO,LXR*LXC)
C?     WRITE(6,*) ' TEST : K2, IOFFX, K1, LXC, LXR ',
C?   &                     K2, IOFFX, K1, LXC, LXR
       CALL MATML7(SCR(K2),XL(IOFFX),SCR(K1),LXC,LXR,LXR,LXC,LXR,LXR,
     &             ZERO,ONE,1)
C?     WRITE(6,*) ' Half transformed matrix '
C?     CALL WRTMAT(SCR(K2),LXC,LXR,LXC,LXR)

*. XL(T) (IBLOCK) A(IBLOCK) XR(IBLOCK)
       CALL SETVEC(SCR(K1),ZERO,LXC*LXC)
       CALL MATML7(SCR(K1),SCR(K2),XR(IOFFX),LXC,LXC,LXC,LXR,LXR,LXC,
     &             ZERO,ONE,0)
C?     WRITE(6,*) ' Transformed matrix '
C?     CALL WRTMAT(SCR(K1),LXC,LXC,LXC,LXC)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK(SCR(K1),AOUT(IOFFP_OUT),1,LXC,LXC)
       ELSE
         CALL COPVEC(SCR(K1),AOUT(IOFFC_OUT),LXC*LXC)
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL APRBLM2(AOUT,LX_COL,LX_COL,NBLOCK,ISYM)      
      END IF
*
      RETURN
      END
      SUBROUTINE TRAN_SYM_BLOC_MAT4b
     &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrices XL, XR
*  X to yield blocked matrix AOUT
*
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = XL(transposed) A XR
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION LX_ROW(NBLOCK),LX_COL(NBLOCK)
      DIMENSION XL(*),XR(*)
      integer AIN
*. Output 
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block 
      DIMENSION SCR(*)
*
*. To get rid of annoying and incorrect compiler warnings
      IOFFP_IN = 0
      IOFFC_IN = 0
      IOFFP_OUT = 0 
      IOFFC_OUT = 0
      IOFFX = 0
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP_IN = 1
         IOFFC_IN = 1
         IOFFP_OUT= 1
         IOFFC_OUT= 1
         IOFFX = 1
       ELSE
         IOFFP_IN = IOFFP_IN + LX_ROW(IBLOCK-1)*(LX_ROW(IBLOCK-1)+1)/2
         IOFFC_IN = IOFFC_IN + LX_ROW(IBLOCK-1) ** 2                     
         IOFFP_OUT= IOFFP_OUT+ LX_COL(IBLOCK-1)*(LX_COL(IBLOCK-1)+1)/2
         IOFFC_OUT= IOFFC_OUT+ LX_COL(IBLOCK-1) ** 2                     
         IOFFX = IOFFX + LX_ROW(IBLOCK-1)*LX_COL(IBLOCK-1)
       END IF
       LXR = LX_ROW(IBLOCK)
       LXC = LX_COL(IBLOCK)
       K1 = 1
       K2 = 1 + MAX(LXR,LXC) ** 2
*. Unpack block of A
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
CNW      CALL TRIPAK(SCR(K1),AIN(IOFFP_IN),2,LXR,LXR)
         if (LXR.gt.0)
     &   call ga_get(AIN,IOFFP_IN,(IOFFP_IN+LXR*(LXR+1)/2-1),1,1,
     &               SCR(K2),LXR*LXR)
         CALL TRIPAK(SCR(K1),SCR(K2),2,LXR,LXR)
C             TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
       ELSE
         if (LXR.gt.0) call
     &      ga_get(AIN,IOFFC_IN,IOFFC_IN+LXR*LXR-1,1,1,SCR(K1),LXR*LXR)
       END IF
*. XL(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
C?     if (ga_nodeid().eq.0) then
C?     WRITE(6,*) ' Non transformed matrix '
C?     CALL WRTMAT(SCR(K1),LXC,LXC,LXC,LXC)
C?     endif
       CALL SETVEC(SCR(K2),ZERO,LXR*LXC)
C?     WRITE(6,*) ' TEST : K2, IOFFX, K1, LXC, LXR ',
C?   &                     K2, IOFFX, K1, LXC, LXR
       CALL MATML7(SCR(K2),XL(IOFFX),SCR(K1),LXC,LXR,LXR,LXC,LXR,LXR,
     &             ZERO,ONE,1)
C?     if (ga_nodeid().eq.0) then
C?     WRITE(6,*) ' Half transformed matrix '
C?     CALL WRTMAT(SCR(K2),LXC,LXR,LXC,LXR)
C?     endif

*. XL(T) (IBLOCK) A(IBLOCK) XR(IBLOCK)
       CALL SETVEC(SCR(K1),ZERO,LXC*LXC)
       CALL MATML7(SCR(K1),SCR(K2),XR(IOFFX),LXC,LXC,LXC,LXR,LXR,LXC,
     &             ZERO,ONE,0)
C?     if (ga_nodeid().eq.0) then
C?     WRITE(6,*) ' Transformed matrix '
C?     CALL WRTMAT(SCR(K1),LXC,LXC,LXC,LXC)
C?     endif
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK(SCR(K1),AOUT(IOFFP_OUT),1,LXC,LXC)
       ELSE
         CALL COPVEC(SCR(K1),AOUT(IOFFC_OUT),LXC*LXC)
       END IF
*
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL APRBLM2(AOUT,LX_COL,LX_COL,NBLOCK,ISYM)      
      END IF
*
      RETURN
      END
      FUNCTION I2EAD_NOCCSYM(IORB,JORB,KORB,LORB,NO12SYM)
*
* Find adress of integral in integral list without 
* complex conjugation symmetry
*
* If NO12SYM = 1, there is no symmetry between particle 1 and 2 
*
c      IMPLICIT REAL*8           (A-H,O-Z)
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*
*
      IF(NO12SYM.EQ.0) THEN
        I2EAD_NOCCSYM = 
     &  I2EADS_NOCCSYM(IORB,JORB,KORB,LORB,WORK(KPINT2_SIMTRH),0)
      ELSE 
        I2EAD_NOCCSYM = 
     &  I2EADS_NOCCSYM(IORB,JORB,KORB,LORB,WORK(KPINT2_SIMTRH_AB),1)
      END IF
*
      RETURN
      END
      FUNCTION I2EADS_NOCCSYM(IORB,JORB,KORB,LORB,IJKLOF,NO12SYM)
*
* Obtain address of integral (IORB JORB ! KORB LORB) in MOLCAS order 
* IORB JORB KORB LORB corresponds to SYMMETRY ordered indeces 
*
* Routine for integrals without complex conjugation symmetry
* If NO12SYM = 1, there is no symmetry between particle 1 and 2 
* Integrals assumed in core 
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      Dimension IJKLOF(NsmOB,NsmOb,NsmOB)
      Logical  ikSymjl         
*. 
      NTEST = 00
*
      IABS = IORB
      ISM = ISMFTO(IREOST(IORB))
      IOFF = IBSO(ISM)
*
      JABS = JORB
      JSM = ISMFTO(IREOST(JORB))
      JOFF = IBSO(JSM)
*
      KABS = KORB
      KSM = ISMFTO(IREOST(KORB))
      KOFF = IBSO(KSM)
*
      LABS = LORB
      LSM = ISMFTO(IREOST(LORB))
      LOFF = IBSO(LSM)
*
      If( Ntest.ge. 100) THEN
        write(6,*) ' I2EADS_NOCCSYM at your service '
        WRITE(6,*) ' IORB IABS ISM IOFF ',IORB,IABS,ISM,IOFF
        WRITE(6,*) ' JORB JABS JSM JOFF ',JORB,JABS,JSM,JOFF
        WRITE(6,*) ' KORB KABS KSM KOFF ',KORB,KABS,KSM,KOFF
        WRITE(6,*) ' LORB LABS LSM LOFF ',LORB,LABS,LSM,LOFF
      END IF
*
      iSym=iSm
      jSym=jSm
      I = IABS - IOFF + 1
      J = JABS - JOFF + 1
C     ijBlk=jSym+iSym*(iSym-1)/2
      IJBLK = (ISYM-1)*NSMOB + JSYM
      kSym=kSm
      lSym=lSm
      K = KABS - KOFF + 1
      L = LABS -LOFF + 1
C     klBlk=lSym+kSym*(kSym-1)/2
      KLBLK = (KSYM-1)*NSMOB + LSYM
*
      If ( klBlk.gt.ijBlk.AND.NO12SYM.EQ.0 ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(6,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      ikSymjl=.false.         
      IF(NO12SYM.EQ.0.AND.Isym.EQ.Ksym.AND.Jsym.EQ.Lsym) 
     &ikSymjl = .true.
*
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
      IF(NTEST.GE.100) THEN
        print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      END IF
      ijPairs=itOrb*jtOrb
      ij=j + (i-1)*jtOrb
*
      klPairs=ktOrb*ltOrb
      kl=l+(k-1)*ltOrb
      IF(NTEST.GE.100)  THEN
        print *,' ijPairs,klPairs',ijPairs,klPairs
      END IF
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 ) write(6,*) ' ijkl ', ijkl
*
      I2EADS_NOCCSYM = iblOff-1+ijkl
      If( ntest .ge. 100 ) then
        write(6,*) 'i j k l ', i,j,k,l
        write(6,*) ' ibloff ijkl ',ibloff,ijkl
        write(6,*) ' I2EADS_NOCCSYM  = ', I2EADS_NOCCSYM
      END IF
*
      RETURN
      END 
      SUBROUTINE GETINCN2_NOCCSYM(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,IJKLOF,NSMOB,I2INDX,
     &                  ICOUL,INO12SYM) 
*
* Obtain integrals from list without complex conjugation symmetry
*
*     ICOUL = 0 : 
*                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 : 
*                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 2 :  XINT(IL,JK) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0 
* 
*
* Version for integrals stored in INTLST
*
* If type equals zero, all integrals of given type are fetched 
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
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
      IF(ITP.GE.1) THEN
        iOrb=NOBPTS(ITP,ISM)
      ELSE
        IORB = NTOOBS(ISM)
      END IF
*
      IF(JTP.GE.1) THEN
        jOrb=NOBPTS(JTP,JSM)
      ELSE
        JORB = NTOOBS(JSM)
      END IF
*
      IF(KTP.GE.1) THEN
        kOrb=NOBPTS(KTP,KSM)
      ELSE
        KORB = NTOOBS(KSM)
      END IF
*
      IF(LTP.GE.1) THEN
        lOrb=NOBPTS(LTP,LSM)
      ELSE
        LORB = NTOOBS(LSM)
      END IF
*
C?    WRITE(6,*) ' ITP, JTP, KTP, LTP = ', 
C?   &             ITP, JTP, KTP, LTP
C?    WRITE(6,*) ' IKSM, JLSM = ', IKSM, JLSM
C?    WRITE(6,*) ' IORB, JORB, KORB, LORB = ',
C?   &             IORB, JORB, KORB, LORB
*
*. Offsets relative to start of all orbitals, symmetry ordered 
      IOFF = IBSO(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NOBPTS(IITP,ISM)
      END DO
*
      JOFF = IBSO(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NOBPTS(JJTP,JSM)
      END DO
*
      KOFF = IBSO(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NOBPTS(KKTP,KSM)
      END DO
*
      LOFF = IBSO(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NOBPTS(LLTP,LSM)
      END DO
C?    WRITE(6,*) ' IOFF, JOFF, KOFF, LOFF ',
C?   &             IOFF, JOFF, KOFF, LOFF

*
*     Collect Coulomb terms
*
      ijblk = (ism-1)*NSMOB + jsm
      klblk = (ksm-1)*NSMOB + lsm
*
      IJRELKL = 0
      IBLOFF = 0
      IF (INO12SYM.EQ.1) THEN
       IJRELKL = 1
       IBLOFF=IJKLOF(ISM,JSM,KSM)                            
      ELSE IF(IJBLK.GT.KLBLK) THEN
       IJRELKL = 1
       IBLOFF=IJKLOF(ISM,JSM,KSM)                            
      ELSE IF (IJBLK.EQ.KLBLK) THEN
       IJRELKL = 0
       IBLOFF=IJKLOF(ISM,JSM,KSM)                            
      ELSE IF (IJBLK.LT.KLBLK) THEN
       IJRELKL = -1
       IBLOFF = IJKLOF(KSM,LSM,ISM)                                
      END IF
*
      itOrb=NTOOBS(iSm)
      jtOrb=NTOOBS(jSm)
      ktOrb=NTOOBS(kSm)
      ltOrb=NTOOBS(lSm)
*
       IJPAIRS = ITORB*JTORB
       KLPAIRS = KTORB*LTORB
*
      iInt=0
      Do lJeppe=lOff,lOff+lOrb-1
C?      WRITE(6,*) ' ljeppe = ', ljeppe
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lJeppe
C?      WRITE(6,*) ' Range of  J = ', jMin,jOff+jOrb-1
        Do jJeppe=jMin,jOff+jOrb-1
C?      WRITE(6,*) ' jjeppe = ', jjeppe
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
C?      WRITE(6,*) ' kjeppe = ', kjeppe
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
C?      WRITE(6,*) ' ijeppe = ', ijeppe
                  ijkl = ijkl0 + I2INDX(IJEPPE+IJ0)
                  iInt=iInt+1
                  Xint(iInt) = Intlst(ijkl)
C?                WRITE(6,*) ' IINT 1 = ', IINT
C?                WRITE(6,*) ' IJKL, IJKL0, IBLOFF, XINT(iINT) = ',
C?   &                         IJKL, IJKL0, IBLOFF, XINT(iINT)
              End Do
            END IF
*
*. block (ISM JSM !ISM JSM)
            IF(IJRELKL.EQ.0) THEN 
              IJ0 = (JJEPPE-1)*NTOOB         
              KLOFF = KL*(KL-1)/2
              IJKL0 = (KL-1)*IJPAIRS-KLOFF
              Do iJeppe=iMin,iOff+iOrb-1
C?      WRITE(6,*) ' ijeppe = ', ijeppe
                ij = I2INDX(IJEPPE+IJ0   )
                If ( ij.ge.kl ) Then
C                 ijkl=ij+(kl-1)*ijPairs-klOff
                  IJKL = IJKL0 + IJ
                Else
                  IJOFF = IJARR(IJEPPE)
                  ijkl=kl+(ij-1)*klPairs-ijOff
                End If
                iInt=iInt+1
                Xint(iInt) = Intlst(iblOff-1+ijkl)
C?                WRITE(6,*) ' IINT 2 = ', IINT
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
                Xint(iInt) = Intlst(ijkl)
C?                WRITE(6,*) ' IINT 3 = ', IINT
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
       ILPAIRS = ITORB*LTORB
       KJPAIRS = KTORB*JTORB
*
        ilblk = (ism-1)*NSMOB + lsm                                
        kjblk = (ksm-1)*NSMOB + jsm
        ILRELKJ = 0
        IF(ILBLK.GT.KJBLK) THEN
          ILRELKJ = 1
          IBLOFF = IJKLOF(ISM,LSM,KSM)                                
        ELSE IF(ILBLK.EQ.KJBLK) THEN
          ILRELKJ = 0
          IBLOFF = IJKLOF(ISM,LSM,KSM)                                
        ELSE IF(ILBLK.LT.KJBLK) THEN
          ILRELKJ = -1
          IBLOFF = IJKLOF(KSM,JSM,ISM)                               
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
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
C?                WRITE(6,*) ' IINT 4 = ', IINT
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
                  XInt(iInt)=XInt(iInt)-Intlst(iBLoff-1+ilkj)
C?                WRITE(6,*) ' IINT 5 = ', IINT
                End Do
              END IF
*
              IF(ILRELKJ.EQ.-1) THEN
                ILKJ0 = IBLOFF-1+KJ-KJPAIRS
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE+ IL0)*KJPAIRS
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
C?                WRITE(6,*) ' IINT 6 = ', IINT
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
      SUBROUTINE TRA2_SIMPLE(C,X2INTOUT)
*
* Trivial 2-electron integral transformation
*
*. Input integrals in KINT2
*. Output integrals in X2INTOUT
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      DIMENSION C(*)
*. Output
      integer X2INTOUT
CNW   DIMENSION X2INTOUT(*)
*
      NTEST  = 00
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA2_S')
*. Largest symmetry block
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Two symmetry blocks
      LENBL = MXSOB ** 4
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' Size of symmetry block ', LENBL
      CALL  MEMMAN(KLBL1,LENBL,'ADDL ',2,'KLBL1 ')
      CALL  MEMMAN(KLBL2,LENBL,'ADDL ',2,'KLBL2 ')
C     WRITE(6,*) ' Two integral blocks allocated '
*. Loop over symmetry blocks
      DO ISM = 1, NSMOB
        DO JSM = 1, ISM
          DO KSM = 1, ISM   
            IF(KSM.EQ.ISM) THEN
              LSM_MX = JSM 
            ELSE
              LSM_MX = KSM
            END IF
            DO LSM = 1, LSM_MX
*. Ensure that integrals have correct symmetry
              INTSYM = 1
              IJSM = MULTD2H(ISM,JSM)
              IJKSM = MULTD2H(IJSM,KSM)
              IF(INTSYM.EQ.MULTD2H(IJKSM,LSM)) THEN
*. Fetch 2-e integral block (IJ!KL)
              ONE = 1.0D0
C             WRITE(6,*)  'ISM, JSM, KSM, LSM ',ISM,KSM,KSM,LSM
              CALL LUCIAGETINT(dbl_mb(KLBL1),-1,ISM,-1,JSM,-1,KSM,-1,
     &                         LSM,0,0,0,1,ONE,ONE)
*. (Type = 0  => Complete symmetryblock)
*. Offsets and dimensions for symmetryblocks in C 
              IOFF = 1
              DO IISM = 1, ISM-1
                IOFF = IOFF + NTOOBS(IISM)**2
              END DO
              NI = NTOOBS(ISM)
*
              JOFF = 1
              DO JJSM = 1, JSM-1
                JOFF = JOFF + NTOOBS(JJSM)**2
              END DO
              NJ = NTOOBS(JSM)
*
              KOFF = 1
              DO KKSM = 1, KSM-1
                KOFF = KOFF + NTOOBS(KKSM)**2
              END DO
              NK = NTOOBS(KSM)
*
              LOFF = 1
              DO LLSM = 1, LSM-1
                LOFF = LOFF + NTOOBS(LLSM)**2
              END DO
              NL = NTOOBS(LSM)
*. Transform 2-electron integral block
              CALL TRA_2EL_BLK_SIMPLE(dbl_mb(KLBL1),
     &        C(IOFF),NI,C(JOFF),NJ,C(KOFF),NK,C(LOFF),NL,dbl_mb(KLBL2))
*. Transfer symmetry block to integral list
              CALL PUTINT(dbl_mb(KLBL1),0,ISM,0,JSM,0,KSM,0,LSM,
     &                    0,0,X2INTOUT,int_mb(KPINT2))
C     PUTINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                  NOCCSYM,NO12SYM,XOUT)
            END IF
*           ^ Check if integrals have correct symmetry
            END DO
          END DO
        END DO
      END DO
*     ^ End of loop over symmetries
c*. Dump integrals to file LU90 
c      LU90 = 90
c      CALL DMPINT(LU90)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRA2_S')
*
      RETURN
      END
      FUNCTION GTIJKL_SM_AB(I,J,K,L,NA,NB)
*
* Obtain  similarity transformed integral (I J ! K L )
* for integral with NA alpha indeces and NB beta indeces 
* For integrals NA = 2 and NB = 2, it is assumed that I J 
* corresponds to the alpha-indeces
*
* I,J,K L refers to active orbitals in  Type ordering
*                                      ==============
*
* Jeppe Olsen, July 2002
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'oper.inc'
*
      XIJKL = 0.0D0
*
c      IROUTE=1
c      IF (IROUTE==1) THEN
      ISPCAS=0 ! to track unexpected cases
      IF (NA.EQ.4) ISPCAS=1
      IF (NB.EQ.4) ISPCAS=2
      IF (NA.EQ.2) ISPCAS=3

      ! use the standard function
      XIJKL = GTIJKL(I,J,K,L)
*
c      ELSE
c      IF(NA.EQ.4) THEN
c* a a a a integral
c        IADR =  I2EAD_NOCCSYM(IREOTS(I),IREOTS(J),
c     &          IREOTS(K),IREOTS(L),0)
c        XIJKL = WORK(KINT2_SIMTRH_AA-1+IADR)
c      ELSE IF (NB.EQ.4) THEN
c* b b b b integral
c        IADR =  I2EAD_NOCCSYM(IREOTS(I),IREOTS(J),
c     &          IREOTS(K),IREOTS(L),0)
c        XIJKL = WORK(KINT2_SIMTRH_BB-1+IADR)
c      ELSE IF (NA.EQ.2) THEN
c* a a b b integral
c        IADR =  I2EAD_NOCCSYM(IREOTS(I),IREOTS(J),
c     &          IREOTS(K),IREOTS(L),1)
c        XIJKL = WORK(KINT2_SIMTRH_AB-1+IADR)
c      END IF
c      END IF
*
      GTIJKL_SM_AB = XIJKL
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,'(a,6i4,a,g20.10)')
     &       ' 2e int. for I,J,K,L,NA,NB = ', I,J,K,L,NA,NB,' is ',
     &       XIJKL
      END IF
      RETURN
      END
      SUBROUTINE READ_WRITE_SH(LUSH)
*
* Read S and H from DALTON outputfile and write to LUSH
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Title - stored on real .....
      DIMENSION TITLE(24)
      PARAMETER (LBUF=600)
*. Bufffer for reading integrals 
      DIMENSION BUF(LBUF),IBUF(LBUF)
*. Info on number of basis functions per symmetry 
      INTEGER NDEG(8)
*. Max number of orbitals
      PARAMETER (MXPORB = 100)
*. An array for overlap integrals
      DIMENSION S(MXPORB*MXPORB), H(MXPORB*MXPORB)
*
* Obtain info on number of symmetries and basis functions from 
* file  AOONEINT generated by DALTON
      ITAP34 = 66
      OPEN (ITAP34,STATUS='OLD',FORM='UNFORMATTED',FILE='AOONEINT')
      REWIND ITAP34
      READ (ITAP34) TITLE,NST,(NDEG(I),I=1,NST),ENUC
*
*. Report ..
*
*
*. Title 
      WRITE(6,'(//A,2(/12A6)/)')
     *   ' Molecule title from basis set input :',(TITLE(I),I=1,24)
      WRITE(6,*) ' Number of symmetries ', NST
      WRITE(6,*) ' Number of basis functions per symmetry ',
     &           (NDEG(ISM),ISM = 1, NST)
*
* Read in overlap matrix 
*
      LEN = 0
      DO ISM = 1, NST
        LEN = LEN + NDEG(ISM)*(NDEG(ISM)+1)/2 
      END DO
      WRITE(6,*) ' LEN = ', LEN
      ZERO = 0.0D0
      DO I = 1, LEN
        S(I) = ZERO
        H(I) = ZERO
      END DO
*
*. Overlap 
*
      CALL MOLLAB('OVERLAP ',ITAP34,6)
 2100 READ (ITAP34) (BUF(I),I=1,LBUF),(IBUF(I),I=1,LBUF),LENGTH
      DO I = 1,LENGTH
         S(IBUF(I)) = BUF(I)
      END DO
*. Write out overlap matrix 
C          WRITE_H1_WITH_INDEX(H1,NSMOB,NOBPSM,LUOUT)
      CALL WRITE_H1_WITH_INDEX(S,NST,NDEG,LUSH)
*
*. And one-electron Hamiltonian on the same file 
*
      REWIND ITAP34
      CALL MOLLAB('ONEHAMIL',ITAP34,6)
 3100 READ (ITAP34) (BUF(I),I=1,LBUF),(IBUF(I),I=1,LBUF),LENGTH
      DO 3200 I = 1,LENGTH
         H(IBUF(I)) = BUF(I)
 3200 CONTINUE
      IF (LENGTH .GE. 0) GO TO 3100
      CALL WRITE_H1_WITH_INDEX(H,NST,NDEG,LUSH)
      WRITE(LUSH,'(F25.15)')  ENUC
      

*
      RETURN
      END 
      SUBROUTINE WRITE_H1_WITH_INDEX(H1,NSMOB,NOBPSM,LUOUT)
*
* Write one-electron matrix on file LUOUT with symmetry index 
* Integrals are assumed total symmetric and lower half packed 
*
* Jeppe Olsen, Nov. 2003 for Transfer of info to Kims VB program
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION H1(*)
      INTEGER NOBPSM(*)
*
      IOFF = 1
      IJ = 0
      DO IJSM = 1, NSMOB
        DO IOB = 1, NOBPSM(IJSM)
          DO JOB = 1, IOB
            IJ = IJ + 1
            IOB_EFF = IOB + IOFF - 1
            JOB_EFF = JOB + IOFF - 1
            WRITE(LUOUT,'(2I5,F25.15)') IOB_EFF,JOB_EFF,H1(IJ)
          END DO
        END DO
        IOFF = IOFF + NOBPSM(IJSM)
      END DO
*
      RETURN 
      END
      SUBROUTINE TRA2_G_SIMPLE(KLCI,KLCJ,KLCK,KLCL,IP2INTOUT,X2INTOUT)
*
* Simple generel 2-electron integral transformation
*
*.  Jeppe Olsen, Feb. 2011, part of the LUCIA Growing up campaign
*
*. Input integrals in KINT2 and are defined by NTOOBS, I12SM, I34SM, I1234SM
*  with symmetry blocks adressed by KPINT2
*. Output integrals in X2INTOUT  and defined by the *_A indeces in CINTFO
*. and ORBINP and pointers IP2INTOUT
*
* The pointers to the MO-AO transformation matrices are in KLC*, so the
* first transformation matrix is in WORK(KLCI)
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
CNW   DIMENSION X2INTOUT(*)
      integer X2INTOUT
*. Input array
      DIMENSION IP2INTOUT(NSMOB,NSMOB,NSMOB)
*. Local scratch
      INTEGER IPERM_SSSS(7,16)
      integer nxtask, next, itask, icount, nproc
      external nxtask
      NTEST = 00
      NTESTO = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' =========================='
        WRITE(6,*) ' Output form TRA2_G_SIMPLE '
        WRITE(6,*) ' =========================='
        WRITE(6,'(A,3I2)')
     &  ' Input  integrals: I12S, I34S, I1234S ',
     &                      I12S, I34S, I1234S
        WRITE(6,'(A,3I2)') 
     &  ' Output integrals: I12S_A, I34S_A, I1234S_A: ',
     &                      I12S_A, I34S_A, I1234S_A
       END IF
       IF(NTEST.GE.20) THEN
*
        WRITE(6,*) ' CI, CJ, CK, CL = '
        CALL APRBLM2(dbl_mb(KLCI),NTOOBS,NTOOBS,NSMOB,0)
        CALL APRBLM2(dbl_mb(KLCJ),NTOOBS,NTOOBS,NSMOB,0)
        CALL APRBLM2(dbl_mb(KLCK),NTOOBS,NTOOBS,NSMOB,0)
        CALL APRBLM2(dbl_mb(KLCL),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA2_G')
*. Largest symmetry block of input integrals
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Three symmetry blocks of 2ei
      LEN4IBL = MXSOB ** 4
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' Size of Integral symmetry block ', LEN4IBL
      CALL  MEMMAN(KLBL1,LEN4IBL,'ADDL ',2,'KLBL1 ')
      CALL  MEMMAN(KLBL2,LEN4IBL,'ADDL ',2,'KLBL2 ')
      CALL  MEMMAN(KLBL3,LEN4IBL,'ADDL ',2,'KLBL3 ')
C     WRITE(6,*) ' Two integral blocks allocated '
*. four blocks of transformation matrices
      LEN2IBL = MXSOB**2
      CALL MEMMAN(KLINTM1,LEN2IBL,'ADDL  ',2,'INTM1 ')
      CALL MEMMAN(KLINTM2,LEN2IBL,'ADDL  ',2,'INTM2 ')
  
*. Loop over input symmetry blocks
      DO ISM = 1, NSMOB
        IF(I12S.EQ.1) THEN
         JSM_MX = ISM
        ELSE
         JSM_MX = NSMOB
        END IF
        DO JSM = 1, JSM_MX
          IF(I1234S.EQ.1) THEN
            KSM_MX = ISM
          ELSE 
            KSM_MX = NSMOB
          END IF
          DO KSM = 1, KSM_MX
            IF(I34S.EQ.1) THEN
             LSM_MX = KSM
            ELSE
             LSM_MX = NSMOB
            END IF
            IF(I1234S.EQ.1.AND.ISM.EQ.KSM) THEN
             LSM_MX = MIN(LSM_MX,JSM)
            END IF
            DO LSM = 1, LSM_MX
*. Ensure that integrals have correct symmetry
              INTSYM = 1
              IJSM = MULTD2H(ISM,JSM)
              IJKSM = MULTD2H(IJSM,KSM)
              IF(INTSYM.EQ.MULTD2H(IJKSM,LSM)) THEN
* Fusk
C             NTEST = NTESTO
C             IF(ISM.EQ.6.AND.JSM.EQ.2.AND.KSM.EQ.5.AND.LSM.EQ.1) THEN
C               NTEST = 1000
C               WRITE(6,*) ' Jeppe has raised print in Transf '
C             END IF
*
              IF(NTEST.GE.100) THEN
              WRITE(6,'(A,4I2)') ' 1: ISM, JSM, KSM, LSM = ', 
     &                     ISM, JSM, KSM, LSM
CM            WRITE(6,'(A,4I5)') ' IOFF, JOFF, KOFF, LOFF ',
CM   &                     IOFF, JOFF, KOFF, LOFF
              END IF
*. Fetch 2-e integral block (IJ!KL) with full packing of integrals
              IF(NTEST.GE.100) THEN
                WRITE(6,*)
                WRITE(6,'(A)') 
     &          ' ===================================================='
                WRITE(6,'(A,4I2)') 
     &          ' Info on symmetry block (ISM,JSM,KSM,LSM) ',
     &          ISM,JSM,KSM,LSM
                WRITE(6,'(A)') 
     &          ' ===================================================='
                WRITE(6,*)
              END IF
              ONE = 1.0D0
              CALL GET_2EBLK(dbl_mb(KLBL1),ISM,JSM,KSM,LSM,
     &        int_mb(KPINT2))
*. Number of output blocks originating from this input block
              CALL PERM_SYM_SSSSBLK(ISM,JSM,KSM,LSM,I12S,I34S,I1234S,
     &                              I12S_A,I34S_A,I1234S_A,NPERM_SSSS,
     &                              IPERM_SSSS)
*. Loop over the various output blocks arising from this input block
              DO IPERM = 1, NPERM_SSSS
*. The following is put here, instead of in a separate routine as 
*. it may be convenient to do some common work for transformations 
*. containing more than one block.
*
               I12P = IPERM_SSSS(5,IPERM) 
               I34P = IPERM_SSSS(6,IPERM) 
               I1234P = IPERM_SSSS(7,IPERM) 
               IF(NTEST.GE.1000)
     &         WRITE(6,*) ' I12P, I34P, I1234P =', I12P, I34P, I1234P
*
* Obtained permutated indeces, dimensions... for the integrals in
* in the output order
*
C                  DO_2EI_INDEX_PERM(I12P,I34P,I1234P,
C    &                I_IN,J_IN,K_IN,L_IN,
C    &                I_OUT,J_OUT,K_OUT,L_OUT)
               CALL DO_2EI_INDEX_PERM(I12P,I34P,I1234P,
     &         ISM,JSM,KSM,LSM,ISM_O,JSM_O,KSM_O,LSM_O)
C              IF(NTEST.GE.100) 
C    &         WRITE(6,*) ' ISM_O, JSM_O, KSM_O, LSM_O = ',
C    &         ISM_O, JSM_O, KSM_O, LSM_O
*
*. Offsets and dimensions for symmetries ISM_O, JSM_O, KSM_O, LSM_O
*
               IF(NTEST.GE.100) THEN
                WRITE(6,'(A,4I3)') 
     &          ' Info for output symmetries of I,J, K, L = ',
     &           ISM_O, JSM_O, KSM_O, LSM_O
               END IF
               IOFF = 1
               DO IISM = 1, ISM_O-1
                 IOFF = IOFF + NTOOBS(IISM)*NTOOBS_IA(IISM)
               END DO
               NI_O = NTOOBS_IA(ISM_O)
               NI_BAS = NTOOBS(ISM)
               NI_BASO = NTOOBS(ISM_O)
*
               JOFF = 1
               DO JJSM = 1, JSM_O-1
                 JOFF = JOFF + NTOOBS(JJSM)*NTOOBS_JA(JJSM)
               END DO
               NJ_O = NTOOBS_JA(JSM_O)
               NJ_BAS = NTOOBS(JSM)
               NJ_BASO = NTOOBS(JSM_O)
*
               KOFF = 1
               DO KKSM = 1, KSM_O-1
                 KOFF = KOFF + NTOOBS(KKSM)*NTOOBS_KA(KKSM)
               END DO
               NK_O = NTOOBS_KA(KSM_O)
               NK_BAS = NTOOBS(KSM)
               NK_BASO = NTOOBS(KSM_O)
*
               LOFF = 1
               DO LLSM = 1, LSM_O-1
                 LOFF = LOFF + NTOOBS(LLSM)*NTOOBS_LA(LLSM)
               END DO
               NL_O = NTOOBS_LA(LSM_O)
               NL_BAS = NTOOBS(LSM)
               NL_BASO = NTOOBS(LSM_O)

               IF(NTEST.GE.100) 
     &         WRITE(6,*) ' NI_O, NJ_O, NK_O, NL_O = ',
     &         NI_O, NJ_O, NK_O, NL_O
*
C              CALL DO_2EI_INDEX_PERM(I12P,I34P,I1234P,
C    &         IOFF,JOFF,KOFF,LOFF,IOFF_O,JOFF_O,KOFF_O,LOFF_O)
C              CALL DO_2EI_INDEX_PERM(I12P,I34P,I1234P,
C    &         NI,NJ,NK,NL,NI_O,NJ_O,NK_O,NL_O)
C              CALL DO_2EI_INDEX_PERM(I12P,I34P,I1234P,
C    &         NI_BAS,NJ_BAS,NK_BAS,NL_BAS,
C    &         NI_BASO,NJ_BASO,NK_BASO,NL_BASO)
C              CALL DO_2EI_INDEX_PERM(I12P,I34P,I1234P,
C    &         KLCI,KLCJ,KLCK,KLCL,
C    &         KLCI_O,KLCJ_O,KLCK_O,KLCL_O)
*. Permutational symmetries for output  block
               I12S_O = 0
               I34S_O = 0
               I1234S_O = 0 
               IF(I12S_A.EQ.1.AND.ISM_O.EQ.JSM_O) I12S_O = 1
               IF(I34S_A.EQ.1.AND.KSM_O.EQ.LSM_O) I34S_O = 1
               IF(I1234S_A.EQ.1.AND.ISM_O.EQ.KSM_O.AND.JSM_O.EQ.LSM_O)
     &         I1234S_O = 1
*. Symmetry of input integrals with symmetry in current input order
               I12S_I = 0
               I34S_I = 0
               I1234S_I = 0
               IF(I12S.EQ.1.AND.ISM_O.EQ.JSM_O) I12S_I = 1
               IF(I34S.EQ.1.AND.KSM_O.EQ.LSM_O) I34S_I = 1
               IF(I1234S.EQ.1.AND.ISM_O.EQ.KSM_O.AND.JSM_O.EQ.LSM_O)
     &         I1234S_I = 1
*. Symmetry of input integrals in original input list
               I12S_INI = 0
               I34S_INI = 0
               I1234S_INI = 0
               IF(I12S.EQ.1.AND.ISM.EQ.JSM) I12S_INI = 1
               IF(I34S.EQ.1.AND.KSM.EQ.LSM) I34S_INI = 1
               IF(I1234S.EQ.1.AND.ISM.EQ.KSM.AND.JSM.EQ.LSM)
     &         I1234S_INI = 1
*. Obtain input block of integrals, unpacked in index one and 
*. two and in order of the output block IPERM
C              MOD_2E_SSSSBLK(XIN,XOUT,NI_IN,NJ_IN,NK_IN,NL_NL_IN,
C    &                        I12S_IN,I34S_IN,I1234S_IN,
C    &                        I12S_OUT,I34S_OUT,I1234S_OUT,
C    &                        I12P,I34P,I1234P)
*. For the first half transformation, we want the integrals in the 
*. form (ISM_O,JSM_O!KSM_O,LSM_O) without permutational symmetry 
*. between particle one and two
*. Symmetry of the outpacked untransformed integrals
               I1234S_INTM = 0
*. Dimensions of unpacked untransformed integrals
*. Obtain the untransformed integrals (ISM_O JSM_O ! KSM_O LSM_O)
*  without permutational symmetry between 12 and 34 indices, but retained
*  symmetry in 12 and in 34
               IZERO = 0
               CALL MOD_2E_SSSSBLK(dbl_mb(KLBL1),dbl_mb(KLBL2),
     &              NI_BAS,NJ_BAS,NK_BAS,NL_BAS,
     &              I12S_INI,I34S_INI,I1234S_INI,
     &              I12S_I,I34S_I,IZERO,I12P,I34P,I1234P)
*
* ==========================
* First half transformation
* ==========================
*
*. We now have the block (ISM_O,JSM_O!KSM_O,LSM_O) 
*. Transform first two indeces
*
*. Intermediate 1  is untransformed integrals (i' j' ! k' l')
*. Intermediate 2  is half transformed integrals (i j ! k' l')
*. Intermediate 3 is transformed integrals (k l ! i j)
               IF(I34S_I.EQ.0) THEN
                 NKL_I = NK_BASO*NL_BASO
               ELSE
                 NKL_I = NK_BASO*(NK_BASO+1)/2
               END IF
               IF(I12S_I.EQ.0) THEN
                 NIJ_I = NI_BASO*NJ_BASO
               ELSE
                 NIJ_I = NI_BASO*(NI_BASO+1)/2
               END IF
               IF(I12S_O.EQ.0) THEN
                 NIJ_O = NI_O*NJ_O
               ELSE
                 NIJ_O = NI_O*(NI_O+1)/2
               END IF
               IF(I34S_O.EQ.0) THEN
                 NKL_O = NK_O*NL_O
               ELSE
                 NKL_O = NK_O*(NK_O+1)/2
               END IF
*
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Initial integrals as (ij!kl)'
                 CALL WRTMAT(dbl_mb(KLBL2),NIJ_I,NKL_I,NIJ_I,NKL_I)
               END IF
*
* First Half transformation ( i' j'!k' l') => ( i j!k' l')
*                               KLBL2            KLBL3
*
               DO K = 1, NK_BASO
                IF(I34S_I.EQ.1) THEN
                  LMAX = K
                ELSE
                  LMAX = NL_BASO
                END IF
                DO L = 1, LMAX
                 IF(NTEST.GE.10000) WRITE(6,'(A,2I4)') ' K, L = ', K, L
                 IF(I34S_I.EQ.1) THEN
                   KL = K*(K-1)/2 + L
                 ELSE
                   KL = (L-1)*NK_BASO + K
                 END IF
                 KLOF_IN  = KLBL2 + (KL-1)*NIJ_I
                 KLOF_OUT = KLBL3 + (KL-1)*NIJ_O
*
                 FACTORC = 0.0D0
                 FACTORAB = 1.0D0
*
                 KLCI_O = KLCI + IOFF - 1
                 KLCJ_O = KLCJ + JOFF - 1
*. Obtain (i' j' !k' l') without perm symmetry between i' and j'
                 IF(I12S_I.EQ.0) THEN
                   CALL COPVEC(dbl_mb(KLOF_IN),dbl_mb(KLINTM1),NIJ_I)
                 ELSE
C                       TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
                   CALL TRIPAK(dbl_mb(KLINTM1),dbl_mb(KLOF_IN),2,
     &                         NI_BASO,NI_BASO)
                 END IF
*
*(i' j! k' l' ) = Sum(j') (i' j' !k' l') CJ(j' j)
*
                 IF(NTEST.GE.10000) WRITE(6,*) ' Transformation of j '
                 CALL MATML7(dbl_mb(KLINTM2),dbl_mb(KLINTM1),
     &           dbl_mb(KLCJ_O),
     &           NI_BASO,NJ_O,NI_BASO,NJ_BASO,NJ_BASO,NJ_O,
     &           FACTORC,FACTORAB,0)
* ( i j! k' l' ) = sum(i') CI(i',i) (i' j ! k' l')
*
                 IF(NTEST.GE.10000) WRITE(6,*) ' Transformation of i '
                 CALL MATML7(dbl_mb(KLINTM1),dbl_mb(KLCI_O),
     &           dbl_mb(KLINTM2),
     &           NI_O,NJ_O,NI_BASO,NI_O,NI_BASO,NJ_O,
     &           FACTORC,FACTORAB,1)
* Transfer to KLBL3 and pack IJ if required
                 IF(I12S_O.EQ.0) THEN
                   CALL COPVEC(dbl_mb(KLINTM1),dbl_mb(KLOF_OUT),NIJ_O)
                 ELSE
                   CALL TRIPAK(dbl_mb(KLINTM1),dbl_mb(KLOF_OUT),1,
     &                         NI_O,NI_O)
                 END IF
                END DO
               END DO ! End of K L loops
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Half transformed integrals as (ij!kp lp)'
                 CALL WRTMAT(dbl_mb(KLBL3),NIJ_O,NKL_I,NIJ_O,NKL_I)
               END IF
               CALL MEMCHK2('IJTRAN')
*
* =================================================================
* Transpose (i j !k' l') to obtain (k' l' ! i j) and save in KLBL2
* =================================================================
*
C                   TRPMT3(XIN,NROW,NCOL,XOUT)
               CALL TRPMT3(dbl_mb(KLBL3),NIJ_O,NKL_I,dbl_mb(KLBL2))
C?             WRITE(6,*) ' Transposed matrix (kp lp! i j )'
C?             CALL WRTMAT(WORK(KLBL2),NKL_I, NIJ_O,NKL_I, NIJ_O)
               CALL MEMCHK2('TRANSP')
C?             WRITE(6,*) ' I12S_O, I34S_I = ', I12S_O, I34S_I
               DO I = 1, NI_O
                IF(I12S_O.EQ.0) THEN
                  JMAX = NJ_O
                ELSE
                  JMAX = I
                END IF
C?              WRITE(6,*) ' NI_BASO = ', NI_BASO
                DO J = 1, JMAX
                  IF(I12S_O.EQ.1) THEN
                    IJ = I*(I-1)/2 + J
                  ELSE
C                   IJ = (J-1)*NI_BASO + I
                    IJ = (J-1)*NI_O + I
                  END IF
                  KLOF_IN  = KLBL2 + (IJ-1)*NKL_I    
                  KLOF_OUT = KLBL3 + (IJ-1)*NKL_O    
C?                WRITE(6,*) ' I, J, IJ, NKL_I, NKL_O = ',
C?   &            I, J, IJ, NKL_I, NKL_O
*
                  KLCK_O = KLCK + KOFF -1
                  KLCL_O = KLCL + LOFF -1
*. Obtain (k' l' !i j) without perm symmetry between k' and l'
                  IF(I34S_I.EQ.0) THEN
                    CALL COPVEC(dbl_mb(KLOF_IN),dbl_mb(KLINTM1),NKL_I)
                  ELSE
C                        TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
                    CALL TRIPAK(dbl_mb(KLINTM1),dbl_mb(KLOF_IN),2,
     &                          NK_BASO,NK_BASO)
                  END IF
                  IF(NTEST.GE.1000) THEN
                    WRITE(6,*) ' (kp lp ! i j ) for i,j = ', I,J
                    CALL WRTMAT(dbl_mb(KLINTM1),NK_BASO,NL_BASO,
     &                          NK_BASO,NL_BASO)
                  END IF
*
*(k' l! i j) = Sum(l') (k' l' !i j) CL(l' l)
*
                 IF(NTEST.GE.10000) WRITE(6,*) ' Transformation of l '
                  CALL MATML7(dbl_mb(KLINTM2),dbl_mb(KLINTM1),
     &            dbl_mb(KLCL_O),
     &            NK_BASO,NL_O,NK_BASO,NL_BASO,NL_BASO,NL_O,
     &            FACTORC,FACTORAB,0)
                  IF(NTEST.GE.10000) THEN
                    WRITE(6,*) ' (kp l ! i j ) for i,j = ', I,J
                    CALL WRTMAT(dbl_mb(KLINTM2),NK_BASO,NL_O,
     &                          NK_BASO,NL_O)
                  END IF
*
* ( k l! i j ) = sum(k') CK(k',k) (k' l ! i j)
*
                 IF(NTEST.GE.10000) WRITE(6,*) ' Transformation of k '
                  CALL MATML7(dbl_mb(KLINTM1),dbl_mb(KLCK_O),
     &            dbl_mb(KLINTM2),
     &            NK_O,NL_O,NK_BASO,NK_O,NK_BASO,NL_O,
     &            FACTORC,FACTORAB,1)
* Transfer to KLBL3 and pack KL if required
                  IF(I34S_O.EQ.0) THEN
                    CALL COPVEC(dbl_mb(KLINTM1),dbl_mb(KLOF_OUT),NKL_O)
                  ELSE
                    CALL TRIPAK(dbl_mb(KLINTM1),dbl_mb(KLOF_OUT),1,
     &                          NK_O,NK_O)
                  END IF
                END DO
               END DO! End of loop over IJ
*
               CALL MEMCHK2('KLTRAN')
               IF(NTEST.GE.100) THEN
                 WRITE(6,*) 
     &           ' Final transformed integral block as (kl!ij) '
                 CALL WRTMAT(dbl_mb(KLBL3),NKL_O,NIJ_O,NKL_O,NIJ_O)
               END IF
*
*. Save as ( i j ! k l) in packed form
               IOFF_BL = IP2INTOUT(ISM_O,JSM_O,KSM_O)
               IF(NTEST.GE.1000) 
     &         WRITE(6,*) ' IOFF_BL = ', IOFF_BL
C              MOD_2E_SSSSBLK(XIN,XOUT,NI_IN,NJ_IN,NK_IN,NL_NL_IN,
C    &                         I12S_IN,I34S_IN,I1234S_IN,
C    &                         I12S_OUT,I34S_OUT,I1234S_OUT,
C    &                         I12P,I34P,I1234P)
               CALL MOD_2E_SSSSBLK(dbl_mb(KLBL3),
     &              dbl_mb(X2INTOUT+IOFF_BL-1),
     &              NK_O,NL_O,NI_O,NJ_O,I34S_O,I12S_O,0,
     &              I12S_O,I34S_O,I1234S_O,0,0,1)
*. Finito !!!
              END DO ! End of loop over Permutations
            END IF ! End if correct combination of symmetries
            END DO
          END DO
        END DO
      END DO ! End of loop over KSM, LSM, JSM, ISM
      call ga_sync()
*
      CALL MEMCHK('FINOTO')
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRA2_G')
*
      RETURN
      END
C             CALL GET_2EBLK(WORK(KLBL1),ISM,JSM,KSM,LSM,WORK(KPINT2))
      SUBROUTINE GET_2EBLK(XINT,ISM,JSM,KSM,LSM,IPINT2)
*
* Obtain symmetry block of integrals. Integrals are stored in KINT2
* and permutational symmetry is defined by I12S, I34S, I1234S
* Dimension of integrals defined by NTOOBS
*
*. Jeppe Olsen, Feb 2011, For the Lucia growing up campaign
*
C!      INCLUDE 'implicit.inc'
C!    INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IPINT2(NSMOB,NSMOB,NSMOB)
*. Output
      DIMENSION XINT(*)
*. 
      NTEST = 00
*
      IOFF =  IPINT2(ISM,JSM,KSM)
      IF(IOFF.LE.0) THEN
       WRITE(6,*) ' GET_2EBLK trying to fetch block not stored'
       WRITE(6,*) ' ISM, JSM, KSM, LSM = ', ISM,JSM,KSM,LSM
       STOP       ' GET_2EBLK trying to fetch block not stored'
      END IF
*. Number of integrals
      NIJKL = LEN_2EBLK(ISM,NTOOBS,JSM,NTOOBS,KSM,NTOOBS,
     &                  LSM,NTOOBS,I12S,I34S,I1234S)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' KINT2, NIJKL, IOFF = ', 
     &               KINT2, NIJKL, IOFF
        WRITE(6,*) ' KINT_2EMO = ', KINT_2EMO
      END IF
*
      CALL COPVEC(dbl_mb(KINT2-1+IOFF),XINT,NIJKL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Symmetry block of integrals, GET_2EBLK'
        WRITE(6,'(A,4I2)') 
     &  ' ISM, JSM, KSM, LSM = ', ISM, JSM, KSM, LSM
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)
        NK = NTOOBS(KSM)
        NL = NTOOBS(LSM)
        IF(ISM.EQ.JSM) THEN
          NIJ = NI*(NI+1)/2
        ELSE
          NIJ = NI*NJ
        END IF
        IF(KSM.EQ.LSM) THEN
          NKL = NK*(NK+1)/2
        ELSE
          NKL = NK*NL
        END IF
        IF(ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
          CALL PRSYM_GEN(XINT,NIJ,2)
        ELSE
          CALL WRTMAT(XINT,NIJ,NKL,NIJ,NKL)
        END IF
      END IF

      RETURN
      END
*
      FUNCTION LEN_2EBLK(ISM,NI,JSM,NJ,KSM,NK,LSM,NL,I12S,I34S,I1234S)
*
* Length of 2e integral symmetry block
*
*. Input
      INTEGER NI(*),NJ(*),NK(*),NL(*)
*
      IF(ISM.EQ.JSM.AND.I12S.EQ.1) THEN
       NIJ = NI(ISM)*(NI(ISM)+1)/2
      ELSE
       NIJ = NI(ISM)*NJ(JSM)
      END IF
*
      IF(KSM.EQ.LSM.AND.I34S.EQ.1) THEN
       NKL = NK(KSM)*(NK(KSM)+1)/2
      ELSE
       NKL = NK(KSM)*NL(LSM)
      END IF
*
      IF(ISM.EQ.KSM.AND.JSM.EQ.LSM.AND.I1234S.EQ.1) THEN
        NIJKL = NIJ*(NIJ+1)/2
      ELSE
        NIJKL = NIJ*NKL
      END IF
*
      LEN_2EBLK = NIJKL
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from LEN_2EBLK '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' ISM, JSM, KSM, LSM = ', ISM,JSM,KSM,LSM
        WRITE(6,*) ' NIJ, NKL, NIJKL = ', NIJ,NKL, NIJKL
      END IF
*
      RETURN
      END
      SUBROUTINE PRINT_2EBLK_G(XINT,NI,NJ,NK,NL,IJPSM,KLPSM,IJKLPSM)
*
* Print symmetry block of two-electron integralsi XINT
* Dimensions are NI,NJ,NK,NL
*.Permutational symmetries are IJPSM, KLPSM, IJKLPSM
*
*. Jeppe Olsen, Feb. 2011, for the Lucia growing up campaign
*
      INCLUDE 'implicit.inc'
      DIMENSION XINT(*)
      PARAMETER (MXP_ORBBAT = 10)
      INTEGER IJ_BAT(2,MXP_ORBBAT)
      INTEGER KL_BAT(2,MXP_ORBBAT)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to PRINT_2EBLK_G '
        WRITE(6,*) ' ======================== '
        WRITE(6,'(A,4I4)') ' NI, NJ, NK, NL = ', NI, NJ, NK, NL
        WRITE(6,'(A,3I2)') ' IJPSM, KLPSM, IJKLPSM = ',
     &                       IJPSM, KLPSM, IJKLPSM
      END IF
*
      LCOL_BAT = 5
      LROW_BAT = 10
*
      IF(IJPSM.EQ.1) THEN
        NIJ = NI*(NI+1)/2
      ELSE
        NIJ = NI*NJ
      END IF
*
      IF(KLPSM.EQ.1) THEN
        NKL = NK*(NK+1)/2
      ELSE
        NKL = NK*NL
      END IF
*
      NIJ_BAT = NIJ/LROW_BAT
      NKL_BAT = NKL/LCOL_BAT
      IF(NIJ_BAT*LROW_BAT.LT.NIJ) NIJ_BAT = NIJ_BAT + 1
      IF(NKL_BAT*LCOL_BAT.LT.NKL) NKL_BAT = NKL_BAT + 1
      IF(NTEST.GE.100) 
     &WRITE(6,'(A,2I5)') ' NIJ_BAT, NKL_BAT = ', NIJ_BAT, NKL_BAT
      
      IJ_INI = 1
      IJ_IB = 1
      DO IJ_BATL = 1, NIJ_BAT
C     NEXT_BATCH_ORBPAIRS(NI,NJ,LBAT_MX,IJ_PSM,INI,IJPAIRS,LBAT)
C?     WRITE(6,*) ' IJ_INI before call to NEXT ... ', IJ_INI
       CALL NEXT_BATCH_ORBPAIRS(NI,NJ,LROW_BAT,IJPSM,
     &      IJ_INI,IJ_BAT,LIJ)
       IJ_INI = 0
       IF(LIJ.NE.0) THEN
        KL_OFF = 1
        KL_INI = 1
        KL_IB = 1
        DO KL_BATL = 1, NKL_BAT
         CALL NEXT_BATCH_ORBPAIRS(NK,NL,LCOL_BAT,KLPSM,
     &        KL_INI,KL_BAT,LKL)
         KL_INI = 0
         IF(IJKLPSM.EQ.1) THEN
           IJ_START_ABS = MAX(IJ_IB,KL_IB)
           IJ_START = IJ_START-ABS - IJ_IB + 1
         ELSE
           IJ_START = 1
         END IF
         IF(IJ_START.LE.LIJ) THEN
          WRITE(6,'(11X,8(2X,A1,I4,A1,I4,A1,1X))')
     &    ('|',KL_BAT(1,KL),',',KL_BAT(2,KL),')',KL=1, LKL)
          DO IJ = IJ_START, LIJ
           IF(IJKLPSM.EQ.1) THEN
            IJ_ABS = IJ_IB-1+IJ
            KL_MAX_ABS = MIN(KL_IB-1+LKL,IJ_ABS)
            IF(KL_MAX_ABS.GE.KL_IB)
     &      WRITE(6,'(A1,I4,A1,I4,A1,2X,8(E13.7,1X))') 
     &      '(',IJ_BAT(1,IJ),',',IJ_BAT(2,IJ),'|', 
     &       (XINT((KL-1)*NIJ-KL*(KL-1)/2+IJ_ABS),
     &        KL = KL_IB,KL_MAX_ABS)
           ELSE
            WRITE(6,'(A1,I4,A1,I4,A1,2X,8(E13.7,1X))') 
     &      '(',IJ_BAT(1,IJ),',',IJ_BAT(2,IJ),'|', 
     &       (XINT((KL-1)*NIJ+IJ),
     &        KL = KL_IB,KL_IB-1+LKL)
           END IF
          END DO
         END IF !End of KL_BAT was non-empty
         KL_IB = KL_IB + LKL
        END DO ! End of loop over KL_BAT
       END IF !End if IJ_BAT was non-empty
       IJ_IB = IJ_IB + LIJ
      END DO !End of loop over IJ_BAT
*
      RETURN
      END 
      SUBROUTINE NEXT_BATCH_ORBPAIRS(NI,NJ,
     &           LBAT_MX,IJ_PSM,INI,IJPAIRS,LBAT)
*
*. Obtain Next batch of orbital pairs. Orbital pairs are stored row-wise
*
*. Jeppe Olsen, Feb. 2011
*
      INCLUDE 'implicit.inc'
*. Output
      INTEGER IJPAIRS(2,*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' NEXT_BATCH_ORBPAIRS reporting '
        WRITE(6,*)  ' =============================='
        IF(INI.EQ.1) THEN
          WRITE(6,*) ' Initial batch '
        ELSE
          WRITE(6,*) ' Last indeces of latest batch = ', 
     &                 IJPAIRS(1,LBAT_MX),IJPAIRS(2,LBAT_MX)
        END IF
      END IF
*
      IFINITO = 0
*
*. Initial pair I_INI, J_INI
*
      IF(INI.EQ.1) THEN
       I_INI = 1
       J_INI = 1
      ELSE
       I_FIN = IJPAIRS(1,LBAT_MX)
       J_FIN = IJPAIRS(2,LBAT_MX)
       IF(IJ_PSM.EQ.1) THEN
         IF(I_FIN.GT.J_FIN) THEN
           I_INI = I_FIN 
           J_INI = J_FIN + 1
         ELSE IF (I_FIN.LT.NI) THEN
           I_INI = I_FIN + 1
           J_INI = 1
         ELSE
           IFINITO = 1
         END IF
       ELSE 
         IF(I_FIN.LT.NI)  THEN
           I_INI = I_FIN + 1
           J_INI = J_FIN 
         ELSE IF(J_FIN.LT.NJ) THEN
           I_INI = 1
           J_INI = J_FIN + 1
         ELSE
           IFINITO = 1
         END IF
       END IF
      END IF
*
*. And then the batch
*
      LBAT = 0
      IF(IFINITO.EQ.0) THEN
        LBAT = 1
        IJPAIRS(1,LBAT) = I_INI
        IJPAIRS(2,LBAT) = J_INI
        DO IJ = 2, LBAT_MX
         I = IJPAIRS(1,LBAT)
         J = IJPAIRS(2,LBAT)
         IF(IJ_PSM.EQ.1) THEN
           IF(I.GT.J)  THEN
             LBAT = LBAT + 1
             IJPAIRS(1,LBAT) = I 
             IJPAIRS(2,LBAT) = J+1
           ELSE IF(I.LT.NI) THEN
             LBAT = LBAT + 1
             IJPAIRS(1,LBAT) = I+1
             IJPAIRS(2,LBAT) = 1
           ELSE
            IFINITO = 1
            GOTO 1001
           END IF
         ELSE IF(IJ_PSM.EQ.0) THEN
           IF(I.LT.NI)  THEN
             LBAT = LBAT + 1
             IJPAIRS(1,LBAT) = I + 1
             IJPAIRS(2,LBAT) = J
           ELSE IF(J.LT.NJ) THEN
             LBAT = LBAT + 1
             IJPAIRS(1,LBAT) = 1
             IJPAIRS(2,LBAT) = J + 1
           ELSE
             IFINITO = 1
             GOTO 1001
           END IF
         END IF !End of IJ_PSM switch
        END DO
 1001   CONTINUE
      END IF !End of IFINITO check
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Length of next batch: ', LBAT
        IF(LBAT.GT.0) THEN
          WRITE(6,*) ' Indeces of orbital pairs: '
          DO IJ = 1, LBAT
           WRITE(6,'(2I6)') IJPAIRS(1,IJ),IJPAIRS(2,IJ)
          END DO
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE PERM_SYM_SSSSBLK(ISM,JSM,KSM,LSM,I12S,I34S,I1234S,
     &                            I12S_A,I34S_A,I1234S_A,NPERM_SSSS,
     &                              IPERM_SSSS)
*
* Find possible included permutational combinations of symmetry
* block of 2e- integrals
*
* Jeppe Olsen, Lucia growing up campaign, Feb. 2011
*
      INCLUDE 'implicit.inc'
      DIMENSION IPERM_SSSS(7,16)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' PERM_SYM_SSSSBLK, input: '
        WRITE(6,*) ' ======================== '
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' Symmetries  ISM  JSM  KSM  LSM = ',
     &                                   ISM, JSM, KSM, LSM
        WRITE(6,'(A,3I3)') 
     &  ' Input permutational sym  I12S  I34S  I1234S',
     &                             I12S, I34S, I1234S
        WRITE(6,'(A,3I3)') 
     &  ' Output permutational sym  I12S_A  I34S_A  I1234S_A',
     &                              I12S_A, I34S_A, I1234S_A
      END IF
*
      NPERM = 0
      IF(I12S.EQ.1) THEN
       N12S = 2
      ELSE
       N12S = 1
      END IF
      IF(I34S.EQ.1) THEN
       N34S = 2
      ELSE
       N34S = 1
      END IF
      IF(I1234S.EQ.1) THEN
       N1234S = 2
      ELSE
       N1234S = 1
      END IF
*
      DO I12P = 1, N12S
        IF(I12P.EQ.1) THEN
         ISM1 = ISM
         JSM1 = JSM
         KSM1 = KSM
         LSM1 = LSM
        ELSE
         ISM1 = JSM
         JSM1 = ISM
         KSM1 = KSM
         LSM1 = LSM
        END IF
        DO I34P = 1, N34S
          IF(I34P.EQ.1) THEN
            ISM2 = ISM1
            JSM2 = JSM1
            KSM2 = KSM1
            LSM2 = LSM1
          ELSE 
            ISM2 = ISM1
            JSM2 = JSM1
            KSM2 = LSM1
            LSM2 = KSM1
          END IF
          DO I1234P = 1, N1234S
            IF(I1234P.EQ.1) THEN
             ISM3 = ISM2
             JSM3 = JSM2
             KSM3 = KSM2
             LSM3 = LSM2
            ELSE
             ISM3 = KSM2
             JSM3 = LSM2
             KSM3 = ISM2
             LSM3 = JSM2
            END IF
*. Is ISM3, JSM3, KSM3, LSM3 allowed as output block?
            IM_OKAY = 1
            IF(I12S_A.EQ.1.AND.ISM3.LT.JSM3) IM_OKAY = 0
            IF(I34S_A.EQ.1.AND.KSM3.LT.LSM3) IM_OKAY = 0
            IF(I1234S_A.EQ.1.AND.
     &         (ISM3.LT.KSM3.OR.ISM3.EQ.KSM3.AND.JSM3.LT.LSM3))
     &         IM_OKAY = 0
            IF(IM_OKAY.EQ.1) THEN
*. Combination should be included, check that it has not already been included
              IM_NEW = 1
              DO IPERM = 1, NPERM
                IF(IPERM_SSSS(1,IPERM).EQ.ISM3.AND.
     &             IPERM_SSSS(2,IPERM).EQ.JSM3.AND.
     &             IPERM_SSSS(3,IPERM).EQ.KSM3.AND.
     &             IPERM_SSSS(4,IPERM).EQ.LSM3)    IM_NEW = 0
              END DO
              IF(IM_NEW.EQ.1) THEN
*. New permutation, enroll it
                NPERM = NPERM + 1
*. The symmetries
                IPERM_SSSS(1,NPERM) = ISM3
                IPERM_SSSS(2,NPERM) = JSM3
                IPERM_SSSS(3,NPERM) = KSM3
                IPERM_SSSS(4,NPERM) = LSM3
*. The permutations
                IPERM_SSSS(5,NPERM) = I12P - 1
                IPERM_SSSS(6,NPERM) = I34P - 1
                IPERM_SSSS(7,NPERM) = I1234P - 1
              END IF !End if new
            END IF !End if okay
          END DO !I1234P
        END DO !I34P
      END DO !I12P
*
      NPERM_SSSS = NPERM
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output symmetries: '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' ISM JSM KSM LSM 12P 34P 1234P' 
        WRITE(6,*) ' ================================' 
        DO IPERM = 1, NPERM_SSSS
          WRITE(6,'(7I4)') (IPERM_SSSS(I,IPERM),I=1,7)
        END DO
      END IF !Ntest .ge. 100
*
      RETURN
      END 
      SUBROUTINE MOD_2E_SSSSBLK(XIN,XOUT,NI_IN,NJ_IN,NK_IN,NL_IN,
     &           I12S_IN,I34S_IN,I1234S_IN,
     &           I12S_OUT,I34S_OUT,I1234S_OUT,
     &           I12P,I34P,I1234P)
*
* Modify permutational symmetry packing of a 2-electron integral block
*
* I*S_IN: permutational symmetry of input block
* I*S_OUT: Permutational symmetry of output block
* I*P    : Permutations going from input to outputblock
*
* NI_IN, NJ_IN, NK_IN, NL_IN are number of orbitals in input block
*
*. Jeppe Olsen, The Lucia Growing up campaign, Feb. 2011
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION XIN(*)
*. Output
      DIMENSION XOUT(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MOD_2E_SSSSBLK reporting:'
        WRITE(6,*) ' ========================='
        WRITE(6,*)  
        WRITE(6,'(A,4I3)') 
     &  ' Permutational symmetry of input block  (12, 34, 1234) ',
     &    I12S_IN,I34S_IN,I1234S_IN
        WRITE(6,'(A,4I3)') 
     &  ' Permutational symmetry of Output block (12, 34, 1234) ',
     &    I12S_OUT,I34S_OUT,I1234S_OUT
        WRITE(6,'(A,4I3)') 
     &  ' Permutations for going from in to out block (12, 34, 1234)',
     &    I12P,I34P,I1234P
        WRITE(6,'(A,4I4)') ' Dimensions: NI_IN, NJ_IN, NK_IN, NL_IN = ',
     &                                   NI_IN, NJ_IN, NK_IN, NL_IN
      END IF
*
      IF(I12S_IN.EQ.1) THEN
       NIJ_IN = NI_IN*(NI_IN+1)/2
      ELSE
       NIJ_IN = NI_IN*NJ_IN
      END IF
*
      IF(I34S_IN.EQ.1) THEN
        NKL_IN = NK_IN*(NK_IN+1)/2
      ELSE
        NKL_IN = NK_IN*NL_IN
      END IF
*. Dimensions of output blocks
      IF(I12P.EQ.0) THEN
        NI = NI_IN
        NJ = NJ_IN
      ELSE
        NI = NJ_IN
        NJ = NI_IN
      END IF
      IF(I34P.EQ.0) THEN
       NK = NK_IN
       NL = NL_IN
      ELSE
       NK = NL_IN
       NL = NK_IN
      END IF
      IF(I1234P.EQ.1) THEN
       NIX = NK
       NJX = NL
       NK = NI
       NL = NJ
       NI = NIX
       NJ = NJX
      END IF
*
      IF(I12S_OUT.EQ.1) THEN
        NIJ_OUT = NI*(NI+1)/2
      ELSE
        NIJ_OUT = NI*NJ
      END IF
      IF(I34S_OUT.EQ.1) THEN
       NKL_OUT = NK*(NK+1)/2
      ELSE
       NKL_OUT = NK*NL
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,'(A,6I4)') ' NI, NJ, NK, NK, NIJ_OUT, NKL_OUT = ',
     &               NI, NJ, NK, NK, NIJ_OUT, NKL_OUT
      END IF
*
*. Loop over indices in output block 
      DO K = 1, NK
       IF(I34S_OUT.EQ.1) THEN
        MAX_L = K
       ELSE
        MAX_L = NL
       END IF
       DO L = 1, MAX_L
         IF(I34S_OUT.EQ.1) THEN
           KL_OUT = K*(K-1)/2 + L
         ELSE
           KL_OUT = (L-1)*NK + K
         END IF
         DO I = 1, NI
          IF(I12S_OUT.EQ.1) THEN
            MAX_J = I
          ELSE
            MAX_J = NJ
          END IF
          DO J = 1, MAX_J
           IF(I12S_OUT.EQ.1) THEN
             IJ_OUT = I*(I-1)/2 + J
           ELSE
             IJ_OUT = (J-1)*NI + I
           END IF
*. Input orbital indices
* To go from input to output, we first permute 1/2 3/4 ( if requested)
* and then permute 12/34 (if requested). To go from input to output 
* the order should be reversed.
           IF(I1234P.EQ.0) THEN
            I_IN = I
            J_IN = J
            K_IN = K
            L_IN = L
           ELSE
            I_IN = K
            J_IN = L
            K_IN = I
            L_IN = J
           END IF
           IF(I12P.EQ.1) THEN
            I_INX = I_IN
            I_IN = J_IN
            J_IN = I_INX
           END IF
           IF(I34P.EQ.1) THEN
            K_INX = K_IN
            K_IN = L_IN
            L_IN = K_INX
           END IF
COLD
COLD       IF(I12P.EQ.0) THEN
COLD        I_IN = I
COLD        J_IN = J
COLD       ELSE
COLD        I_IN = J
COLD        J_IN = I
COLD       END IF
COLD       IF(I34P.EQ.0) THEN
COLD        K_IN = K
COLD        L_IN = L
COLD       ELSE
COLD        K_IN = L
COLD        L_IN = K
COLD       END IF
COLD       IF(I1234P.EQ.1) THEN
COLD        I_INX = I_IN
COLD        J_INX = J_IN
COLD        I_IN = K_IN
COLD        J_IN = L_IN
COLD        K_IN = I_INX
COLD        L_IN = J_INX
COLD       END IF
           IF(IJ_OUT.GE.KL_OUT.OR.I1234S_OUT.EQ.0) THEN
*
            IF(I1234S_OUT.EQ.1) THEN
              IJKL_OUT = (KL_OUT-1)*NIJ_OUT -KL_OUT*(KL_OUT-1)/2+IJ_OUT
            ELSE
              IJKL_OUT = (KL_OUT-1)*NIJ_OUT + IJ_OUT
            END IF
*. Input indeces
            IF(I12S_IN.EQ.1) THEN
C             IJ_IN = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
              IJ_IN = MAX(I_IN,J_IN)*(MAX(I_IN,J_IN)-1)/2+
     &                MIN(I_IN,J_IN)
            ELSE
              IJ_IN = (J_IN-1)*NI_IN + I_IN
            END IF 
            IF(I34S_IN.EQ.1) THEN
              KL_IN = MAX(K_IN,L_IN)*(MAX(K_IN,L_IN)-1)/2 
     &              + MIN(K_IN,L_IN)
            ELSE
              KL_IN = (L_IN-1)*NK_IN + K_IN
            END IF 
*
            IF(I1234S_IN.EQ.1) THEN
              IJKL_IN = (MIN(IJ_IN,KL_IN)-1)*NIJ_IN
     &                -MIN(IJ_IN,KL_IN)*(MIN(IJ_IN,KL_IN)-1)/2
     &                + MAX(IJ_IN,KL_IN)
            ELSE
              IF(I1234P.EQ.1) THEN
CC              IJKL_IN = (IJ_IN-1)*NKL_IN + KL_IN
C               IJKL_IN = (IJ_IN-1)*NIJ_IN + KL_IN
                IJKL_IN = (KL_IN-1)*NIJ_IN + IJ_IN
              ELSE
                IJKL_IN = (KL_IN-1)*NIJ_IN + IJ_IN
              END IF
            END IF
*
            IF(NTEST.GE.10000) THEN
             WRITE(6,'(A,4I3,2I6,2I8)') 
     &       'I,J,K,L,IJ_IN,KL_IN, IJKL_IN, IJKL_OUT ',
     &        I,J,K,L,IJ_IN,KL_IN,IJKL_IN,IJKL_OUT
            END IF
*
            XOUT(IJKL_OUT) = XIN(IJKL_IN)
            IF(NTEST.GE.10000)
     &      WRITE(6,*) ' XOUT(IJKL_OUT) = ', XOUT(IJKL_OUT)
           END IF ! Check of IJ_OUT versus KL_OUT
          END DO
         END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output symmetry block'
C            PRINT_2EBLK_G(XINT,NI,NJ,NK,NL,IJPSM,KLPSM,IJKLPSM)
        CALL PRINT_2EBLK_G(XOUT,NI,NJ,NK,NL,I12S_OUT,I34S_OUT,
     &                     I1234S_OUT)
      END IF
        
*
      RETURN
      END
      SUBROUTINE MOD_2E_SSSSBLKG(XIN,XOUT,IOFF,NI_IN,NJ_IN,NK_IN,NL_IN,
     &           I12S_IN,I34S_IN,I1234S_IN,
     &           I12S_OUT,I34S_OUT,I1234S_OUT,
     &           I12P,I34P,I1234P)
*
* Modify permutational symmetry packing of a 2-electron integral block
*
* I*S_IN: permutational symmetry of input block
* I*S_OUT: Permutational symmetry of output block
* I*P    : Permutations going from input to outputblock
*
* NI_IN, NJ_IN, NK_IN, NL_IN are number of orbitals in input block
*
*. Jeppe Olsen, The Lucia Growing up campaign, Feb. 2011
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION XIN(*)
*. Output
      integer XOUT
      dimension xbuf(100)
      integer id1(100), id2(100), icount
CNW   DIMENSION XOUT(*)
*
      call luciaqenter('m2ess')
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MOD_2E_SSSSBLK reporting:'
        WRITE(6,*) ' ========================='
        WRITE(6,*)  
        WRITE(6,'(A,4I3)') 
     &  ' Permutational symmetry of input block  (12, 34, 1234) ',
     &    I12S_IN,I34S_IN,I1234S_IN
        WRITE(6,'(A,4I3)') 
     &  ' Permutational symmetry of Output block (12, 34, 1234) ',
     &    I12S_OUT,I34S_OUT,I1234S_OUT
        WRITE(6,'(A,4I3)') 
     &  ' Permutations for going from in to out block (12, 34, 1234)',
     &    I12P,I34P,I1234P
        WRITE(6,'(A,4I4)') ' Dimensions: NI_IN, NJ_IN, NK_IN, NL_IN = ',
     &                                   NI_IN, NJ_IN, NK_IN, NL_IN
      END IF
*
      do icount = 1, 100
         id2(icount) = 1
      enddo
      icount = 0
*
      IF(I12S_IN.EQ.1) THEN
       NIJ_IN = NI_IN*(NI_IN+1)/2
      ELSE
       NIJ_IN = NI_IN*NJ_IN
      END IF
*
      IF(I34S_IN.EQ.1) THEN
        NKL_IN = NK_IN*(NK_IN+1)/2
      ELSE
        NKL_IN = NK_IN*NL_IN
      END IF
*. Dimensions of output blocks
      IF(I12P.EQ.0) THEN
        NI = NI_IN
        NJ = NJ_IN
      ELSE
        NI = NJ_IN
        NJ = NI_IN
      END IF
      IF(I34P.EQ.0) THEN
       NK = NK_IN
       NL = NL_IN
      ELSE
       NK = NL_IN
       NL = NK_IN
      END IF
      IF(I1234P.EQ.1) THEN
       NIX = NK
       NJX = NL
       NK = NI
       NL = NJ
       NI = NIX
       NJ = NJX
      END IF
*
      IF(I12S_OUT.EQ.1) THEN
        NIJ_OUT = NI*(NI+1)/2
      ELSE
        NIJ_OUT = NI*NJ
      END IF
      IF(I34S_OUT.EQ.1) THEN
       NKL_OUT = NK*(NK+1)/2
      ELSE
       NKL_OUT = NK*NL
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,'(A,6I4)') ' NI, NJ, NK, NK, NIJ_OUT, NKL_OUT = ',
     &               NI, NJ, NK, NK, NIJ_OUT, NKL_OUT
      END IF
*
*. Loop over indices in output block 
      DO K = 1, NK
       IF(I34S_OUT.EQ.1) THEN
        MAX_L = K
       ELSE
        MAX_L = NL
       END IF
       DO L = 1, MAX_L
         IF(I34S_OUT.EQ.1) THEN
           KL_OUT = K*(K-1)/2 + L
         ELSE
           KL_OUT = (L-1)*NK + K
         END IF
         DO I = 1, NI
          IF(I12S_OUT.EQ.1) THEN
            MAX_J = I
          ELSE
            MAX_J = NJ
          END IF
          DO J = 1, MAX_J
           IF(I12S_OUT.EQ.1) THEN
             IJ_OUT = I*(I-1)/2 + J
           ELSE
             IJ_OUT = (J-1)*NI + I
           END IF
*. Input orbital indices
* To go from input to output, we first permute 1/2 3/4 ( if requested)
* and then permute 12/34 (if requested). To go from input to output 
* the order should be reversed.
           IF(I1234P.EQ.0) THEN
            I_IN = I
            J_IN = J
            K_IN = K
            L_IN = L
           ELSE
            I_IN = K
            J_IN = L
            K_IN = I
            L_IN = J
           END IF
           IF(I12P.EQ.1) THEN
            I_INX = I_IN
            I_IN = J_IN
            J_IN = I_INX
           END IF
           IF(I34P.EQ.1) THEN
            K_INX = K_IN
            K_IN = L_IN
            L_IN = K_INX
           END IF
COLD
COLD       IF(I12P.EQ.0) THEN
COLD        I_IN = I
COLD        J_IN = J
COLD       ELSE
COLD        I_IN = J
COLD        J_IN = I
COLD       END IF
COLD       IF(I34P.EQ.0) THEN
COLD        K_IN = K
COLD        L_IN = L
COLD       ELSE
COLD        K_IN = L
COLD        L_IN = K
COLD       END IF
COLD       IF(I1234P.EQ.1) THEN
COLD        I_INX = I_IN
COLD        J_INX = J_IN
COLD        I_IN = K_IN
COLD        J_IN = L_IN
COLD        K_IN = I_INX
COLD        L_IN = J_INX
COLD       END IF
           IF(IJ_OUT.GE.KL_OUT.OR.I1234S_OUT.EQ.0) THEN
*
            IF(I1234S_OUT.EQ.1) THEN
              IJKL_OUT = (KL_OUT-1)*NIJ_OUT -KL_OUT*(KL_OUT-1)/2+IJ_OUT
            ELSE
              IJKL_OUT = (KL_OUT-1)*NIJ_OUT + IJ_OUT
            END IF
*. Input indeces
            IF(I12S_IN.EQ.1) THEN
C             IJ_IN = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
              IJ_IN = MAX(I_IN,J_IN)*(MAX(I_IN,J_IN)-1)/2+
     &                MIN(I_IN,J_IN)
            ELSE
              IJ_IN = (J_IN-1)*NI_IN + I_IN
            END IF 
            IF(I34S_IN.EQ.1) THEN
              KL_IN = MAX(K_IN,L_IN)*(MAX(K_IN,L_IN)-1)/2 
     &              + MIN(K_IN,L_IN)
            ELSE
              KL_IN = (L_IN-1)*NK_IN + K_IN
            END IF 
*
            IF(I1234S_IN.EQ.1) THEN
              IJKL_IN = (MIN(IJ_IN,KL_IN)-1)*NIJ_IN
     &                -MIN(IJ_IN,KL_IN)*(MIN(IJ_IN,KL_IN)-1)/2
     &                + MAX(IJ_IN,KL_IN)
            ELSE
              IF(I1234P.EQ.1) THEN
CC              IJKL_IN = (IJ_IN-1)*NKL_IN + KL_IN
C               IJKL_IN = (IJ_IN-1)*NIJ_IN + KL_IN
                IJKL_IN = (KL_IN-1)*NIJ_IN + IJ_IN
              ELSE
                IJKL_IN = (KL_IN-1)*NIJ_IN + IJ_IN
              END IF
            END IF
*
            IF(NTEST.GE.10000) THEN
             WRITE(6,'(A,4I3,2I6,2I8)') 
     &       'I,J,K,L,IJ_IN,KL_IN, IJKL_IN, IJKL_OUT ',
     &        I,J,K,L,IJ_IN,KL_IN,IJKL_IN,IJKL_OUT
            END IF
*
c           icount = icount + 1
c           id1(icount) = IOFF+IJKL_OUT-1
c           xbuf(icount) = xin(ijkl_in)
c           if (icount.eq.100.or.
c    &         (k.eq.nk.and.l.eq.max_l.and.i.eq.ni.and.j.eq.max_j)) then
c              call ga_scatter(XOUT,xbuf,id1,id2,icount)
c              icount = 0
c           endif
            call ga_put(XOUT,IOFF+IJKL_OUT-1,IOFF+IJKL_OUT-1,1,1,
     &                  XIN(IJKL_IN),1)
CNW         XOUT(IJKL_OUT) = XIN(IJKL_IN)
CNW         IF(NTEST.GE.10000)
CNW  &      WRITE(6,*) ' XOUT(IJKL_OUT) = ', XOUT(IJKL_OUT)
           END IF ! Check of IJ_OUT versus KL_OUT
          END DO
         END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output symmetry block'
C            PRINT_2EBLK_G(XINT,NI,NJ,NK,NL,IJPSM,KLPSM,IJKLPSM)
        CALL PRINT_2EBLK_G(XOUT,NI,NJ,NK,NL,I12S_OUT,I34S_OUT,
     &                     I1234S_OUT)
      END IF
      call luciaqexit('m2ess')
*
      RETURN
      END
      SUBROUTINE GEN_TRA_2EI_LIST
*
* Generate transformed two-electron integral list defined by I2ELIST_A,
* IOCOBTP_A, INTSM_A
*
* Jeppe Olsen, April 2011, for the Lucia growing up campaign, 2010-201X
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
*
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cprnt.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GEN2EI')
*
      NTEST = 000
      IPRNT_TRAINT = 0
      IF(IPRINTEGRAL.GE.1000) IPRNT_TRAINT = 1
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' GEN_TRA_2EI_LIST reporting to service '
        WRITE(6,*) ' ======================================'   
        WRITE(6,*)
        WRITE(6,'(A,3I3)') 
     &  ' IE2LIST_A, IOCOBTP_A, INTSM_A (from cintfo)',
     &    IE2LIST_A, IOCOBTP_A, INTSM_A
        WRITE(6,*) ' PNTGRP = ', PNTGRP
        WRITE(6,*) ' WORK(KKCMO_I) '
        CALL PRINT_CMOAO(dbl_mb(KKCMO_I))
      END IF
*
* Set up information for integral list IE2LIST_A in *_A variables 
* and pointers
*
      CALL PREPARE_2EI_LIST 
*. Flag that only this list is active
      CALL FLAG_ACT_INTLIST(IE2LIST_A)
*
*. Construct the various two-electron integral arrays- there are
*. here constructed seperately, but they should be constructed in 
*. a single transformation- will lead to a significant reduction,
*. but that must be another day.
*
*. Scratch space for 4 dimension arrays and MO-AO transformation
*  matrices
      LMOMO = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      IF(NTEST.GE.10000) WRITE(6,*) ' LMOMO = ', LMOMO
      CALL MEMMAN(KLCI,LMOMO,'ADDL  ',2,'CI    ')
      CALL MEMMAN(KLCJ,LMOMO,'ADDL  ',2,'CJ    ')
      CALL MEMMAN(KLCK,LMOMO,'ADDL  ',2,'CK    ')
      CALL MEMMAN(KLCL,LMOMO,'ADDL  ',2,'CL    ')
*
      DO IE2ARR = 1, IE2LIST_N_A
        IIE2ARR = IE2LIST_I_A(IE2ARR)
        IF(NTEST.NE.0) THEN
         WRITE(6,*)
         WRITE(6,'(A,2I3)') 
     &   ' Info from generation of array (IE2ARR, IIE2ARR) ', 
     &    IE2ARR, IIE2ARR
          WRITE(6,'(A)')
     &    ' ========================================================='
          WRITE(6,*)
        END IF
*. Obtain transformation matrices for the four indeces for this 
*. integral list
        I12S_A = I12S_G(IIE2ARR)
        I34S_A = I34S_G(IIE2ARR)
        I1234S_A = I1234S_G(IIE2ARR)
        IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' I12S_A, I34S_A, I1234S_A = ',
     &               I12S_A, I34S_A, I1234S_A 
        WRITE(6,*) ' IOCOBTP_A = ', IOCOBTP_A
        WRITE(6,*) ' KKCMO_I, KKCMO_J, KKCMO_K, KKCMO_L = ',
     &              KKCMO_I, KKCMO_J, KKCMO_K, KKCMO_L
        END IF
*
        IDIM_ONLY = 0
        
        CALL GET_DIM_AND_C_FOR_ORBS(IOCOBTP_A,INT2ARR_G(1,IIE2ARR),
     &       NTOOBS_IA,NTOOBS_JA,NTOOBS_KA,NTOOBS_LA,
     &       NOBPTS_GN_A(0,1,1),NOBPTS_GN_A(0,1,2),
     &       NOBPTS_GN_A(0,1,3),NOBPTS_GN_A(0,1,4),
     &       dbl_mb(KLCI),dbl_mb(KLCJ),dbl_mb(KLCK),dbl_mb(KLCL),
     &       dbl_mb(KKCMO_I),dbl_mb(KKCMO_J),dbl_mb(KKCMO_K),
     &       dbl_mb(KKCMO_L),IDIM_ONLY)
C     GET_DIM_AND_C_FOR_ORBS(IOCOBTP,IOGLIST,
C    &           NIS,NJS,NKS,NLS,
C    &           NITS,NJTS,NKTS,NLTS,CI,CJ,CK,CL,
C    &           CI_IN, CJ_IN, CK_IN, CL_IN,IDIM_ONLY)
*. We know have the transformation matrices, do it..
        NTEST = 000
        IF(NTEST.GE.10000)
     &  WRITE(6,*) ' KINT2_A(IIE2ARR) = ', KINT2_A(IIE2ARR) 
        CALL TRA2_G_SIMPLE(KLCI,KLCJ,KLCK,KLCL,
     &                     int_mb(KPINT2_A(IIE2ARR)),
     &                     KINT2_A(IIE2ARR) )
        call ga_sync()
C            TRA2_G_SIMPLE(KLCI,KLCJ,KLCK,KLCL,IP2INTOUT,X2INTOUT)
*
        IF(IPRNT_TRAINT.EQ.1.OR.NTEST.GE.10000) THEN
          WRITE(6,*) 
          WRITE(6,*) ' ==================================== '
          WRITE(6,*) '    Array of transformed integrals    '
          WRITE(6,*) ' ==================================== '
          WRITE(6,*)
C              PRINT_2EIARR(XINT,I12S,I34S,I1234S,INTSM,NI,NJ,NK,NL,IPNT)
          CALL PRINT_2EIARR(dbl_mb(KINT2_A(IIE2ARR)), I12S_A, I34S_A, 
     &         I1234S_A, INTSM_A,
     &         NTOOBS_IA,NTOOBS_JA,NTOOBS_KA,NTOOBS_LA,
     &         int_mb(KPINT2_A(IIE2ARR)) )
        END IF
*
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' Leaving GEN_TRA_2EI_LIST '
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GEN2EI')
      RETURN
      END
      SUBROUTINE GET_DIM_AND_C_FOR_ORBS(IOCOBTP,IOGLIST,
     &           NIS,NJS,NKS,NLS,
     &           NITS,NJTS,NKTS,NLTS,CI,CJ,CK,CL,
     &           CI_IN, CJ_IN, CK_IN, CL_IN,IDIM_ONLY)
*
* Obtain dimensions and MO expansion coefficients for four set of orbitals
* defined by IOGLIST. The input MO matrices are 
* in GLBBAS
*
* Jeppe Olsen, April 2011, the Lucia growing up campaign
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
*
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER IOGLIST(4)
      DIMENSION CI_IN(*), CJ_IN(*), CK_IN(*), CL_IN(*)

*. Local scratch
      INTEGER NOCOB_L(MXPOBS), ISUBTP(3)
      CHARACTER*1 OGC(2)
*. Output
      INTEGER NIS(NSMOB),NJS(NSMOB),NKS(NSMOB),NLS(NSMOB)
      INTEGER NITS(0:6+MXPR4T,MXPOBS),NJTS(0:6+MXPR4T,MXPOBS)
      INTEGER NKTS(0:6+MXPR4T,MXPOBS),NLTS(0:6+MXPR4T,MXPOBS)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================================'
        WRITE(6,*) ' GET_DIM_AND_C_FOR_ORBS reporting to service '
        WRITE(6,*) ' ============================================'
        WRITE(6,*)
        WRITE(6,*) ' IDIM_ONLY, IOCOBTP = ', IDIM_ONLY, IOCOBTP
        WRITE(6,*) ' IOGLIST: '
        OGC(1:1) = 'O'
        OGC(2:2) = 'G'
          WRITE(6,'(7A2)')
     &    ' (', OGC(IOGLIST(1)),
     &          OGC(IOGLIST(2)),
     &    ' !', OGC(IOGLIST(3)),
     &          OGC(IOGLIST(4)),') '
      END IF
*
* Generate number of occupied orbital per symmetry
*
      ISUBTP(1) = NGAS
      ISUBTP(2) = 0
      ISUBTP(3) = NGAS + 1
      IF(IOCOBTP.EQ.1) THEN
        NOCSUBTP = 1
      ELSE
        NOCSUBTP = 2
      END IF
      NTOSUBTP = 3
*
* And then the four types of orbitals
*
*( To make the following calls fit one line)
      ID = IDIM_ONLY
*. I
      IF(IOGLIST(1).EQ.1) THEN
        CALL CSUB_FROM_C
     &  (CI_IN,CI,NIS,NITS,NOCSUBTP,ISUBTP,ID)
      ELSE
        CALL CSUB_FROM_C
     &  (CI_IN,CI,NIS,NITS,NTOSUBTP,ISUBTP,ID)
      END IF
*. J
      IF(IOGLIST(2).EQ.1) THEN
        CALL CSUB_FROM_C
     &  (CJ_IN,CJ,NJS,NJTS,NOCSUBTP,ISUBTP,ID)
C CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,IONLY_DIM)
      ELSE
        CALL CSUB_FROM_C
     &  (CJ_IN,CJ,NJS,NJTS,NTOSUBTP,ISUBTP,ID)
      END IF
*.  K
      IF(IOGLIST(3).EQ.1) THEN
        CALL CSUB_FROM_C
     &  (CK_IN,CK,NKS,NKTS,NOCSUBTP,ISUBTP,ID)
      ELSE
        CALL CSUB_FROM_C
     &  (CK_IN,CK,NKS,NKTS,NTOSUBTP,ISUBTP,ID)
      END IF
*. L
      IF(IOGLIST(4).EQ.1) THEN
        CALL CSUB_FROM_C
     &  (CL_IN,CL,NLS,NLTS,NOCSUBTP,ISUBTP,ID)
      ELSE
        CALL CSUB_FROM_C
     &  (CL_IN,CL,NLS,NLTS,NTOSUBTP,ISUBTP,ID)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of orbitals per sym for I,J,K,L: '
        CALL IWRTMA(NIS,1,NSMOB,1,NSMOB)
        CALL IWRTMA(NJS,1,NSMOB,1,NSMOB)
        CALL IWRTMA(NKS,1,NSMOB,1,NSMOB)
        CALL IWRTMA(NLS,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' Number of orbitals per type and sym for I,J,K,L'
        CALL IWRTMA(NITS,NGAS+2,NSMOB,7+MXPR4T,MXPOBS)
        CALL IWRTMA(NJTS,NGAS+2,NSMOB,7+MXPR4T,MXPOBS)
        CALL IWRTMA(NKTS,NGAS+2,NSMOB,7+MXPR4T,MXPOBS)
        CALL IWRTMA(NLTS,NGAS+2,NSMOB,7+MXPR4T,MXPOBS)
*
        IF(IDIM_ONLY.EQ.0) THEN
        WRITE(6,*) ' CI, CJ, CK, CL: '
          CALL APRBLM2(CI,NTOOBS,NIS,NSMOB,0)
          CALL APRBLM2(CJ,NTOOBS,NJS,NSMOB,0)
          CALL APRBLM2(CK,NTOOBS,NKS,NSMOB,0)
          CALL APRBLM2(CL,NTOOBS,NLS,NSMOB,0)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE DO_2EI_INDEX_PERM(I12P,I34P,I1234P,
     &           I_IN,J_IN,K_IN,L_IN,
     &           I_OUT,J_OUT,K_OUT,L_OUT)
*
* Perform Permutations defined by I12P, I34P, I1234P
* on some indeces,*_IN of two-electron integrals to
* produce output indeces *_OUT
*
*. Jeppe Olsen, April 2011
*
      INCLUDE 'implicit.inc'
*
      IF(I12P.EQ.0) THEN
       I_OUT = I_IN
       J_OUT = J_IN
      ELSE
       I_OUT = J_IN
       J_OUT = I_IN
      END IF
*
      IF(I34P.EQ.0) THEN
        K_OUT = K_IN
        L_OUT = L_IN
      ELSE
        K_OUT = L_IN
        L_OUT = K_IN
      END IF
*
      IF(I1234P.EQ.1) THEN
        K_SAVE = K_OUT
        L_SAVE = L_OUT
        K_OUT = I_OUT
        L_OUT = J_OUT
        I_OUT = K_SAVE
        J_OUT = L_SAVE
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from DO_2EI_INDEX_PERM '
        WRITE(6,'(A,4I3)') 
     &  ' Required permutations: I12P,I34P,I1234P =',
     &  I12P,I34P, I1234P
        WRITE(6,'(A,4I5)') 
     &  ' Input indices I, J, K, L ', I_IN,J_IN,K_IN,L_IN
        WRITE(6,'(A,4I5)') 
     &  ' Output indices : ', I_OUT, J_OUT, K_OUT, L_OUT
      END IF
*
      RETURN
      END
      SUBROUTINE PRINT_2EIARR(XINT,I12S,I34S,I1234S,INTSM,NI,NJ,NK,NL,
     &                        IPNT)
*
* Print two-electron integral array XINT
*
*. Jeppe Olsen, May 2011, part of the lucia growing up campaign
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
*. Input: intgrals, pointer to start of blocks, and dimensions
      DIMENSION XINT(*)
      INTEGER IPNT(NSMOB,NSMOB,NSMOB)
      INTEGER NI(*), NJ(*), NK(*), NL(*)
*
      WRITE(6,*) 
      WRITE(6,*) ' ==================================================='
      WRITE(6,*)  ' Integral array: '
      WRITE(6,'(A,I2)')  '   Spatial symmetry = ', INTSM
      WRITE(6,'(A,3I2)') '   Perm. symmetries (12, 34, 1234) = ',
     &             I12S,I34S,I1234S
      WRITE(6,*) ' ==================================================='
      WRITE(6,*)
*
      NTEST = 00
*. Loop over input symmetry blocks
      DO ISM = 1, NSMOB
        IF(I12S.EQ.1) THEN
         JSM_MX = ISM
        ELSE
         JSM_MX = NSMOB
        END IF
        DO JSM = 1, JSM_MX
          IF(I1234S.EQ.1) THEN
            KSM_MX = ISM
          ELSE 
            KSM_MX = NSMOB
          END IF
          DO KSM = 1, KSM_MX
            IF(I34S.EQ.1) THEN
             LSM_MX = KSM
            ELSE
             LSM_MX = NSMOB
            END IF
            IF(I1234S.EQ.1.AND.ISM.EQ.KSM) THEN
             LSM_MX = MIN(LSM_MX,JSM)
            END IF
            DO LSM = 1, LSM_MX
*. Ensure that integrals have correct symmetry
              IF(NTEST.GE.10000) THEN
              WRITE(6,'(A,4I2)') ' 1: ISM, JSM, KSM, LSM = ', 
     &                     ISM, JSM, KSM, LSM
CM            WRITE(6,'(A,4I5)') ' IOFF, JOFF, KOFF, LOFF ',
CM   &                     IOFF, JOFF, KOFF, LOFF
              END IF
              IJSM = MULTD2H(ISM,JSM)
              IJKSM = MULTD2H(IJSM,KSM)
              IF(INTSM.EQ.MULTD2H(IJKSM,LSM)) THEN
                WRITE(6,*)
                WRITE(6,'(A,4I2)')  
     &          ' Block with symmetries ', ISM, JSM, KSM, LSM
                WRITE(6,*) 
     &          ' ================================'
                WRITE(6,*)
*. Correct symmetry combination
                NIS = NI(ISM)
                NJS = NJ(JSM)
                NKS = NK(KSM)
                NLS = NL(LSM)
                IPNTS = IPNT(ISM,JSM,KSM)
                IF(NTEST.GE.100) WRITE(6,'(A,4I3,I7)')
     &          ' NIS, NJS, NKS, NLS, IPNTS = ', 
     &            NIS, NJS, NKS, NLS, IPNTS
*
                I12S_B = 0
                IF(I12S.EQ.1.AND.ISM.EQ.JSM) I12S_B = 1
                I34S_B = 0
                IF(I34S.EQ.1.AND.KSM.EQ.LSM) I34S_B = 1
                I1234S_B = 0
                IF(I1234S.EQ.1.AND.ISM.EQ.KSM.AND.JSM.EQ.LSM)
     &          I1234S_B = 1
*
C     PRINT_2EBLK_G(XINT,NI,NJ,NK,NL,IJPSM,KLPSM,IJKLPSM)
                CALL PRINT_2EBLK_G(XINT(IPNTS),NIS,NJS,NKS,NLS,
     &                             I12S_B,I34S_B,I1234S_B)
              END IF! End of block has correct symmetry
            END DO
          END DO
        END DO
      END DO !End of loops over symmetries
*
      RETURN
      END
      SUBROUTINE GET_INTARR_F4TP(INTARR,IGAS,JGAS,KGAS,LGAS)
*
* Four types of orbital spaces are given. Determine first
* array where this space is included
*
* To check whether a given integral array is active, it is
* checked whether the corresponding pointer to list of integrals
* is positive.
*
* Jeppe Olsen, May 2011, for the LUCIA growing up campaign
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
*. Local scratch
      INTEGER IJKL_GAS(4),IJKL_TP(4)
*
      NTEST = 00
*. The form of the various integral arrays are defined in 
      IJKL_GAS(1) = IGAS
      IJKL_GAS(2) = JGAS
      IJKL_GAS(3) = KGAS
      IJKL_GAS(4) = LGAS
*
      IF(NTEST.GE.10) THEN 
        WRITE(6,*) ' Info from GET_INTARR_F4TP '
        WRITE(6,'(A,4I4)') ' IGAS, JGAS, KGAS, LGAS = ',
     &                       IGAS, JGAS, KGAS, LGAS
      END IF
*
      IARR_FOUND = 0
      DO IARR = 1, NE2ARR
*. Is this array active
        IF(KINT2_A(IARR).GT.0) THEN
           IOCOBTP = IOCOBTP_G(IARR)
           DO INDEX = 1, 4
             IG = IJKL_GAS(INDEX)
             IF((IOCOBTP.EQ.2.AND.IG.EQ.0).OR.
     &          (1.LE.IG.AND.IG.LE.NGAS)     ) THEN
*. Occupied index
               IJKL_TP(INDEX) = 1
             ELSE
*. General index
               IJKL_TP(INDEX) = 2
             END IF
           END DO
*
           IF(NTEST.GE.100) THEN
            WRITE(6,*) ' IJKL_TP = '
            CALL IWRTMA(IJKL_TP,1,4,1,4)
            WRITE(6,*) ' IARR, INT2ARR_G(*,IARR) ', IARR
            CALL IWRTMA(INT2ARR_G(1,IARR),1,4,1,4)
           END IF
*
           INCLUDED = 1
*. Check that a general index in input is not matched by occupied in array
           DO INDEX  = 1, 4
            IF(IJKL_TP(INDEX).EQ.2.AND.INT2ARR_G(INDEX,IARR).EQ.1)
     &      INCLUDED = 0
           END DO
           IF(INCLUDED.EQ.1) IARR_FOUND = IARR
           IF(IARR_FOUND.NE.0) GOTO 101
        END IF
      END DO
  101 CONTINUE
      INTARR = IARR_FOUND
*
      IF(INTARR.EQ.0) THEN
        WRITE(6,*) ' Problem in GET_INTARR_F4TP '
        WRITE(6,*) ' Integral array not found '
        WRITE(6,'(A,4I3)') ' IGAS, JGAS, KGAS, LGAS = ',
     &               IGAS,JGAS,KGAS,LGAS
cVOG    WRITE(6,'(A,10(2X,I9))') ' KINT2_A(*): ',
        WRITE(6,'(A,10(2X,I17))') ' KINT2_A(*): ', 
     &  (KINT2_A(IARR),IARR=1,NE2ARR)
        STOP ' Problem in GET_ANTARR_F4TP'
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from GET_INTARR_F4TP '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' Input: IGAS, JGAS, KGAS, LGAS = ',
     &               IGAS,JGAS,KGAS,LGAS
        WRITE(6,'(A,I2)') ' First active array containing block ',
     &  INTARR
      END IF
*
      RETURN
      END
       
        
      

      



    

     
      
      
    

      
             

          
          
      

c $Id$
