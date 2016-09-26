
* single- and double- amplitudes 
*
* =================
* Singles (E_{IJ}) 
* =================
*
* form 1 : Only symmetry correct terms included
*          organization :
*          Loop over symmetry of I
*          Find symmetry of J
*           Loop over type of I (a+)
*            Loop over type of J ( a)
*             IF (a+i aj is allowed ) :
*              Loop over J
*               Loop over I
*               End of loop over I
*              End of loop over J
*             End if ( allowed types )
*            End of loop over type of J
*           End of loop over type of I
*          End of loop over symmetry of I
*
* Form 2 : in complete matrix form IJ = (J-1)*NORB + I,
*          where I and J are the absolute orbital numbers
*          and the orbitals are type ordered
*
*
* ===================
* Doubles E(IJ) E(KL)
* ===================
*
* Form 1 : Nonredundant form without singlet/triplet splitting
*
* Loop over symmetry of I
*  Loop over symmetry of J
*   Loop over symmetry of K
*     => Symmetry of L
*    Loop over type of I
*     Loop over type of J
*      Loop over type of K
*       Loop over type of L
*        If (type is allowed ) then
*         Loop over L
*          Loop over K
*           Loop over J
*            Loop over I
*            End of loop over I
*           End of Loop over J
*          End of Loop over K
*         End of Loop over L
*        End if (types are allowed )
*       End of loop over types of L
*      End of loop over types of K
*     End of loop over type over J
*    End of loop over types of I
*
* Form 2 : Nonredundant form with singlet/triplet splitting
*
*   A singlet-singlet coupled
*     Loop structure as above
*   B triplet-triplet coupled
*     Loop structure as above
*
* Form 3 : redundant form in symmetric matrix (IJ,KL),IJ.GT.KL
*          and the orbitals are type ordered
*          ( this form will probably disappear later but ... )

*
      SUBROUTINE BIO_TO_STANDARD(V_BIO,V_STANDARD,IWAY)
*
* Transform s vector between standard format and biorthgonal 
*
* Jeppe Olsen, Sept. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Specific input/output
      DIMENSION V_BIO(*),V_STANDARD(*)
*. For the 2s He case ( STANDARD => BIO)
      V_BIO(1) = 0.5*V_STANDARD(1)
      V_BIO(2) = 0.25*V_STANDARD(2)
*
      RETURN
      END
C     CALL GET_DIAG_BLMAT(WORK(KFI),WORK(KFDIA),NSMOB,NTOOB,1)
      SUBROUTINE GET_DIAG_BLMAT(A,DIAG,NBLK,LBLK,ISYM)
*
* Extract diagonal from blocked matrix
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION A(*)  
      INTEGER LBLK(*)
*. Output
      DIMENSION DIAG(*)
*
      IOFF1 = 1
      IOFF2 = 1
C?    WRITE(6,*) '  GET_DIA..  NBLK = ', NBLK
      DO IBLK = 1, NBLK
        L = LBLK(IBLK)
C?      WRITE(6,*) ' IBLK and L ', IBLK,L
        CALL COPDIA(A(IOFF2),DIAG(IOFF1),L,ISYM)
*
        IOFF1 = IOFF1 + L
        IF(ISYM.EQ.1) THEN
          IOFF2 = IOFF2 + L*(L+1)/2
        ELSE
          IOFF2 = IOFF2 + L ** 2
        END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Diagonal elements of blocked matrix '
        LEN = IOFF1 - 1
        CALL WRTMAT(DIAG,1,LEN,1,LEN)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_H0_DIAG(DIAG,I_SYM,F,ISCALE)              
*
* Obtain SCF H0 matrix  :
*
*           <HF![F,E]!HF>
*           <HF![F,EE]!HF>
* If Iscale .ne. 0, the diagonal is scaled as
*
* diag(a,i,b,j) => 1/2*(2-delta(ai,bj)) diag(a,i,b,j)
*
* Jeppe Olsen, Sept. 98
*
* Blocks are packed column wise
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
*. Input : Diagonal elements of F
      DIMENSION F(*)
*. Output
      DIMENSION DIAG(*)         
*. To get rid of compiler warnings
      IJBL_GT_KLBL = 0
      I_MIN = 0
*
* One-electron part
*
      IJ = 0
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,I_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*. Offsets
            I_OFF = IOBPTS(IGAS,ISM)
            J_OFF = IOBPTS(JGAS,JSM)
*. Numbers
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            DO J = J_OFF,J_OFF+NJ-1
            DO I = I_OFF,I_OFF+NI-1
              IJ = IJ + 1             
              DIAG(IJ) = F(I) - F(J)
C?            write(6,*) ' I J IJ DIAG(IJ) ', I,J,IJ,DIAG(IJ)
            END DO
            END DO
          END IF
        END DO
        END DO
      END DO
*
* Two-electron part
*
      IJKL = IJ
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,I_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
C?        WRITE(6,*) ' IGAS,JGAS,KGAS,LGAS,IJKL_ACT',
C?   &                 IGAS,JGAS,KGAS,LGAS,IJKL_ACT
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1  
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1  
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0  
            ELSE
              IJBL_GT_KLBL = -1 
            END IF
          END IF
C?        WRITE(6,*) ' IJBL_GT_KLBL' ,  IJBL_GT_KLBL
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            DO L = L_OFF,L_OFF+NL-1
            DO K = K_OFF,K_OFF+NK-1
C
            IF(IJBL_GT_KLBL .EQ. 0 ) THEN
             J_MIN = L
            ELSE
             J_MIN = J_OFF
            END IF
            DO J = J_MIN,J_OFF+NJ-1
*
            IF( IJBL_GT_KLBL .EQ. 1 ) THEN
              I_MIN = I_OFF
            ELSE IF ( IJBL_GT_KLBL .EQ. 0 ) THEN
              IF(J.GT.L) THEN
                I_MIN = I_OFF 
              ELSE
                I_MIN = K
              END IF
            END IF
            DO I = I_MIN,I_OFF+NI-1
C
*
              IJKL = IJKL + 1   
C?            WRITE(6,*) 'I,J,K,L,IJKL',
C?   &                    I,J,K,L,IJKL        
*. I and K corresponds to creation ops, J and L to annihilation ops 
              DIAG(IJKL) = F(I) + F(K) - F(J) - F(L)
*. Scale
              IF(ISCALE.EQ.1) THEN
                IF(I.EQ.K.AND.J.EQ.L) DIAG(IJKL) = 0.5D0*DIAG(IJKL)
              END IF
            END DO
            END DO
            END DO
            END DO
*           ^ End of loop over orbitals over given TS
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(ISCALE.EQ.0) THEN
          WRITE(6,*) ' Unscaled <HF![F,t]!HF> matrix '
        ELSE
          WRITE(6,*) '   scaled <HF![F,t]!HF> matrix '
        END IF
        CALL WRT_CC_VEC(DIAG,LU)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_CC_VEC2(CC,LU,CCTYPE)
*
* Print vector of CC amplitudes, type is governed by CCTYPE
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input 
c      INCLUDE 'mxpdim.inc' 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'ctcc.inc'
C     COMMON/CTCC/KLSOBEX,NSPOBEX_TP,KLLSOBEX,KLIBSOBEX,LEN_T_VEC,
C    &             MX_ST_TSOSO,MX_ST_TSOSO_BLK

      CHARACTER*6 CCTYPE
      DIMENSION CC(*) 
*
      IF(CCTYPE(1:2).EQ.'CC') THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' Single excitation part '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        CALL WRT_SX(CC,1)
*
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' Double excitation part '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        CALL WRT_DX1(CC(1+NSXE),1)
      ELSE IF (CCTYPE(1:6).EQ.'GEN_CC') THEN
        CALL WRTBLKN(CC,NSPOBEX_TP,WORK(KLLSOBEX))
      ELSE
        STOP 'UNKNOWN TYPE IN CC_WRT_VEC2'
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_CC_VEC(CC,LU)
*
* Print vector of CC amplitudes 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Generel input 
      INCLUDE 'mxpdim.inc' 
      INCLUDE 'crun.inc'
      DIMENSION CC(*) 
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Single excitation part '
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_SX(CC,1)
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Double excitation part '
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_DX1(CC(1+NSXE),1)
*
      RETURN
      END 
      SUBROUTINE CC_VEC_FNC(CC_AMP,CC_VEC,E_CC,E_CC_A,
     &                      VEC1,VEC2,IBIO,CCTYPE,
     &                      CC_VEC2) 
*
* Calculate energy and CC vector function for 
* a set of CC amplitudes defined by CC_AMP
*
* Jeppe Olsen, Initiated summer of 98
*              Last modified : June 2001
*
*
* The CC vector function reads
*
*     <\my! Exp(-T) H Exp (T) !0>
*
* A note on VEC1 VEC2 (VEC3/C2) (June 2001)
* These arrays are only needed for the older versions 
* of the code, and may be dummy indeces for the newer versions
*
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*.Input
      DIMENSION CC_AMP(*)
      CHARACTER*6 CCTYPE
*. Output
      DIMENSION CC_VEC(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*), CC_VEC2(*)
*. 
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
      COMMON/CMXCJ/MXCJ,MAXK1_MX
*. Local scratch 
      INTEGER ISX(2*MXPNGAS*MXPNGAS)
*
      REAL*8 INPRDD, INPROD
* 
      STOP 'FATAL: CALL TO OLD CC_VEC_FNC'

      RETURN
      END

      SUBROUTINE CC_VEC_FNC2(CC_AMP,CC_VEC,E_CC,E_CC_A,
     &                       XAMPNRM,XVECNRM,XLAMPNRM,
     &                       VEC1,VEC2,IBIO,CCTYPE,
     &                       CC_VEC2,
     &                       LU_AMP,LU_VECF,LULAMP,
     &                       LUINT1,LUINT2,LUINT3
     &     )
*
* Calculate energy and CC vector function for 
* a set of CC amplitudes defined by CC_AMP
*
* Jeppe Olsen, Initiated summer of 98
*              Last modified : June 2001
*
*
* The CC vector function reads
*
*     <\my! Exp(-T) H Exp (T) !0>
*
* A note on VEC1 VEC2 (VEC3/C2) (June 2001)
* These arrays are only needed for the older versions 
* of the code, and may be dummy indeces for the newer versions
*
* LUINT1, LUINT2: files for passing intermediates 
*  old vector function: pass   e^T|ref> on LUINT1=LU_EXRF,
*                       pass  He^T|ref> on LUINT2=LUHEXRF
*
*  ECC                  pass   e^L|ref> on LUINT3=LU_LXRF
*
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*.Input
      DIMENSION CC_AMP(*)
      CHARACTER*6 CCTYPE
*. Output
      DIMENSION CC_VEC(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*), CC_VEC2(*)
*. 
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
      COMMON/CMXCJ/MXCJ,MAXK1_MX
*. Local scratch 
      INTEGER ISX(2*MXPNGAS*MXPNGAS)
*
      REAL*8 INPRDD, INPROD
* 
      CALL ATIM(CPU0,WALL0)

      LBLK = -1
      NTEST = 0
      NTEST = MAX(NTEST,IPRCC)
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' CC_VEC_FNC in action '
        WRITE(6,*) ' ====================='
        WRITE(6,*) '  CCFORM : ',CCFORM(1:6)
      END IF
*
* Load T amplitudes to CC_AMP
      CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LU_AMP)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' T-amplitudes at start of CC_VEC_FNC'
        CALL WRT_CC_VEC2(CC_AMP,6,CCTYPE)
      END IF
       
      XAMPNRM = SQRT(INPROD(CC_AMP,CC_AMP,N_CC_AMP))
*
      IF(I_DO_NEWCCV.EQ.1.OR.I_DO_NEWCCV.EQ.2) THEN
*----------------------------------------------------------------------*
*. Use new CC vector function 
*----------------------------------------------------------------------*
*.     
C?      WRITE(6,*) ' ISIMTRH set to zero '
C?      ISIMTRH = 0
C?      WRITE(6,*) ' I_USE_SIMTRH set to zero '
C?      I_USE_SIMTRH = 0
*..
*..
        IF(ISIMTRH.EQ.1.OR.I_BCC.EQ.1) THEN
*. Inactivate single excitations   
          CALL GET_SPOBTP_FOR_EXC_LEVEL(1,WORK(KLCOBEX_TP),NSPOBEX_TP,
     &         NSXTP,ISX,WORK(KLSOX_TO_OX))
          IZERO = 0

          CALL ISCASET(WORK(KLSPOBEX_AC),IZERO,ISX,NSXTP)
        END IF

        IF(ISIMTRH.EQ.1) THEN
*. Transform integrals to Exp(-T1)HExp(T1) basis
          ICC_EXC = 0
          I_USE_SIMTRH = 0
          I_UNRORB = 0
*. The inactive Fock matrix will be calculated in TRA_SIM, 
*. Ensure that the initial PH division is used 
          CALL ISWPVE(IPHGAS,IPHGAS1,NGAS)
          IF(IREFTYP.NE.2) THEN
*. Orbital integral transformation 
            IO_OR_SO_TRA = 1
          ELSE 
*. Spinorbital integral transformation 
            IO_OR_SO_TRA = 2
          END IF
C?        WRITE(6,*) ' Before call to TRA_SIMTRH : '
C?        WRITE(6,*) ' IREFTYP,  IO_OR_SO_TRA = ', IREFTYP,IO_OR_SO_TRA
*
          CALL TRA_SIMTRH(CC_AMP,IO_OR_SO_TRA)
          CALL ISWPVE(IPHGAS,IPHGAS1,NGAS)
*. Tell GETINT to use similarity transformed integrals
          I_USE_SIMTRH = 1
        END IF
*.      ^ End if similarity transformed Hamiltonian is used 
        IF (IREFTYP.EQ.2) I_UNRORB = 1

        IF(I_DO_NEWCCV.EQ.1) THEN
          WRITE(6,*) ' New CC_VECTOR function in use '
          CALL EXP_MT_H_EXP_T(CC_AMP,CC_VEC,CC_VEC2,E_CC_A)
          E_CC = CC_VEC(N_CC_AMP+1)
        ELSE IF(I_DO_NEWCCV.EQ.2) THEN
          WRITE(6,*) ' Very new CC_VECTOR function in use '
          CALL EXP_MT_H_EXP_TC(CC_AMP,CC_VEC)
          E_CC = CC_VEC(N_CC_AMP+1)+ECORE
        END IF

*. Clean up, make single excitations active again
        IF(ISIMTRH.EQ.1.OR.I_BCC.EQ.1) THEN
          IONE = 1
          CALL ISCASET(WORK(KLSPOBEX_AC),IONE,ISX,NSXTP)
          ECORE = ECORE_INI
          I_USE_SIMTRH = 0
        END IF
      ELSE
*----------------------------------------------------------------------*
*. Use old CC vector function
*.  All kinds of gradients and vector functions of prototype
*.  implementations happen to be here:
*----------------------------------------------------------------------*
*
* Get some scratch files
        LUVS1 = IOPEN_NUS('CCVFSCR1') ! -> lusc1
        LUVS2 = IOPEN_NUS('CCVFSCR2') ! -> lusc2
        LUVS3 = IOPEN_NUS('CCVFSCR3') ! -> lusc3
        LUVS4 = IOPEN_NUS('CCVFSCR4') ! -> lusc34
        LUVS5 = IOPEN_NUS('CCVFSCR5') ! -> lusc35
        LUVS6 = IOPEN_NUS('CCVFSCR6') 
*
* ==========================
* .1 : Exp(T) !0> (on LUHC)
* ==========================
*
        LU_EXRF = LUINT1
        LUHEXRF = LUINT2
        LU_LXRF = LUINT3
        IF (LUINT2.LT.0) LUHEXRF=-LUINT2
* Zero order state is assumed on LUC, results will be stored on LUHC
* CC amplitudes are NO MORE assumed located in KCC1
       MX_TERM = 10
       ICC_EXC = 1
       THRES_E = 0.0D0   
*
       ICSPC = IETSPC
       ISSPC = IETSPC
*
       IF(ISIMTRH.EQ.1) THEN
*. Inactivate single excitations   
         CALL GET_SPOBTP_FOR_EXC_LEVEL(1,WORK(KLCOBEX_TP),NSPOBEX_TP,
     &        NSXTP,ISX,WORK(KLSOX_TO_OX))
         IZERO = 0
         CALL ISCASET(WORK(KLSPOBEX_AC),IZERO,ISX,NSXTP)
       END IF
*
       IOPTYP=0
       IF (CCFORM(1:4).EQ.'UCC2') THEN
         IOPTYP=-1
         MX_TERM=1
       ELSE IF (CCFORM(1:3).EQ.'UCC') THEN
         IOPTYP=-1
         MX_TERM=-150
       END IF

       CALL EXPT_REF2(LUC,LU_EXRF,LUVS1,LUVS2,LUVS3,THRES_E,MX_TERM,
     &             CC_AMP,CC_VEC,VEC1,VEC2,N_CC_AMP,CCTYPE,IOPTYP)
       ! CC_VEC used as scratch if IOPTYP.NE.0

       IF(NTEST.GE.10) THEN
         X  = INPRDD(VEC1,VEC2,LU_EXRF,LU_EXRF,1,LBLK)
         WRITE(6,*) ' Square norm of Exp(t) |ref>  = ', X
       END IF
CM     WRITE(6,*) ' Memcheck after first EXPT_REF '
CM     CALL MEMCHK
CM     WRITE(6,*) ' Memcheck passsed '
*
* ==========================
*. 2 H Exp(T) !0> (On LUVS1)
* ==========================
*
       ICC_EXC = 0
*
       IF(ISIMTRH.EQ.1) THEN
c*. Ensure that the initial PH division is used 
c         CALL ISWPVE(IPHGAS,IPHGAS1,NGAS)
*. Transform integrals to Exp(-T1)HExp(T1) basis
         I_USE_SIMTRH = 0
         I_UNRORB = 0
         IO_OR_SO_TRA = 1
         IF(IREFTYP.NE.2) THEN
*. Orbital integral transformation 
           IO_OR_SO_TRA = 1
         ELSE 
*. Spinorbital integral transformation 
           IO_OR_SO_TRA = 2
         END IF

         CALL TRA_SIMTRH(CC_AMP,IO_OR_SO_TRA)

         IF (IUSE_PH.EQ.1.AND.IREFTYP.EQ.2) THEN
           CALL COPVEC(WORK(KFI_AL),WORK(KINT1_SIMTRH_A),NTOOB**2)
           CALL COPVEC(WORK(KFI_BE),WORK(KINT1_SIMTRH_B),NTOOB**2)
         END IF
c         CALL ISWPVE(IPHGAS,IPHGAS1,NGAS)
*. Tell GETINT to use similarity transformed integrals
         I_USE_SIMTRH = 1
         IF (IREFTYP.EQ.2) I_UNRORB = 1
       END IF
C?    WRITE(6,*) ' I12 after TRA_SIMTRH = ', I12
*
       IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Input file to MV7 '
         CALL WRTVCD(VEC1,LU_EXRF,1,LBLK)
       END IF
       MAXK1_MX = 0
       IF(CCTYPE(1:6).EQ.'GEN_CC'.AND.CCFORM(1:3).EQ.'TCC'.AND.
     &      LUINT2.GT.0) THEN
         ISSPC = ITSPC 
       END IF
C?    WRITE(6,*) ' I12 before MV7 = ', I12
       CALL MV7(VEC1,VEC2,LU_EXRF,LUHEXRF,0,0)
*. clean up time
       IF(ISIMTRH.EQ.1) THEN 
         I_USE_SIMTRH = 0
       END IF
C?    WRITE(6,*) ' I12 after MV7 = ', I12
CM     WRITE(6,*) ' Memcheck after MV7 '
CM     CALL MEMCHK
CM     WRITE(6,*) ' Memcheck passsed '
*
       IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Output file from MV7 '
         CALL WRTVCD(VEC1,LUHEXRF,1,LBLK)
       END IF
       LBLK = -1
       IF(NTEST.GE.10) THEN
         X  = INPRDD(VEC1,VEC2,LUHEXRF,LUHEXRF,1,LBLK)
         WRITE(6,*) ' Square norm of sigma = ', X
       END IF
*
* =================================
* 3 Exp(-T) H Exp(T) !0> (On LUVS2)
* =================================
*
       IF(CCFORM(1:3).EQ.'TCC'.OR.CCFORM(1:3).EQ.'ECC') THEN
         ONEM = -1.0D0
         ICC_EXC = 1
         CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
         IF (CCFORM(1:3).EQ.'TCC') THEN
           ICSPC = ITSPC
         ELSE
           ICSPC = IETSPC
         END IF
         CALL EXPT_REF2(LUHEXRF,LUVS2,
     &               LUVS3,LUVS4,LUVS5,THRES_E,MX_TERM,
     &               CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
         CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
         IF(NTEST.GE.10) THEN
           X  = INPRDD(VEC1,VEC2,LUVS2,LUVS2,1,LBLK)
           WRITE(6,*) ' Square norm of exp(-T)H exp(T) |HF>  = ', X
         END IF
         
         IF (CCFORM(1:3).EQ.'TCC') THEN
*. vector was delivered in  ITSPC, Obtain |0> in ITSPC on LUVS4
           LBLK = -1
           CALL EXPCIV(IREFSM,IETSPC,LUC,ITSPC,LUVS4,
     &          LBLK,LUVS5,1,0,IDC,NTEST)
           ISSPC = ITSPC
           ICSPC = ITSPC
*. Coupled cluster energy
           E_CC = INPRDD(VEC1,VEC2,LUVS4,LUVS2,1,LBLK)
*
           I_SUBTRACT_E = 1
           IF(I_SUBTRACT_E.EQ.1) THEN
*  exp(-T)H exp(T) |HF> =>  exp(-T)H exp(T) |HF> - E |HF>
             ONE = 1.0D0
             EMINUS = (-1.0D0)*E_CC
             CALL VECSMD(VEC1,VEC2,ONE,EMINUS,LUVS2,LUVS4,LUVS5,1,LBLK)
C               VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
             CALL COPVCD(LUVS5,LUVS2,VEC1,1,LBLK)
           END IF
*          ^ End if -E|Ref> should be subtracted
         END IF ! TCC
       END IF ! CCFORM TCC or ECC
*
       IF(ISIMTRH.EQ.1) THEN 
*. Reactivate all excitations
         IONE = 1
         CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TP)
       END IF
*
* ========================================================
* 4   <0!T+mu Exp(-T) H Exp(T) !0>  for TCC
*   2*<0!Exp(T)+ T+mu  (H-E) Exp(T) !0>/<0!Exp(T)+Exp(T)!0> for VCC
*  
*   and energy
* ========================================================
*
*
       IF(CCFORM(1:3).EQ.'TCC') THEN
*. TCC : Densities with !0> and  Exp(-T) H Exp(T) !0>
         CALL DEN_GCC(VEC1,VEC2,LUVS4,LUVS2,CC_VEC)
       ELSE IF(CCFORM(1:3).EQ.'ECC') THEN
*. ECC : generate modified projection manifold
         
*  <ref|e^{L^+} :
         MXTERM = 100
         IF (CCFORM(1:4).EQ.'ECC2') MX_TERM=1
         IF (CCFORM(1:4).EQ.'ECC3') MX_TERM=2
         IF (CCFORM(1:4).EQ.'ECC4') MX_TERM=3
         IF (CCFORM(1:4).EQ.'ECC5') MX_TERM=4
         IF (CCFORM(1:4).EQ.'ECC6') MX_TERM=5
         IF (CCFORM(1:4).EQ.'ECC7') MX_TERM=6
         IF (CCFORM(1:4).EQ.'ECC8') MX_TERM=7
         IF (CCFORM(1:4).EQ.'ECC9') MX_TERM=8

         CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LULAMP)
         XLAMPNRM=SQRT(INPROD(CC_AMP,CC_AMP,N_CC_AMP))
         CALL EXPT_REF2(LUC,LU_LXRF,LUVS1,LUVS3,LUVS4,THRES_E,MX_TERM,
     &        CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)

*  and the density
         CALL DEN_GCC(VEC1,VEC2,LU_LXRF,LUVS2,CC_VEC)

         IF (CCFORM(4:4).NE.' '.AND.CCFORM(4:4).NE.'-') THEN

* to be solved better: we now need one more term in the expansion
*     we rely on the fact, that the last 1/N!L^N|ref> is on LUVS1
c           FAC = 1D0/DBLE(MXTERM + 1)
c           CALL SIG_GCC(VEC1,VEC2,LUVS1,LUVS3,CC_AMP)
c           CALL VECSMD(VEC1,VEC2,FAC,1D0,LUVS3,LU_LXRF,LUVS1,1,LBLK)
c           CALL COPVCD(LUVS1,LU_LXRF,VEC1,1,LBLK)
* does not work for some reason, so brute force:
           IF (CCFORM(1:4).EQ.'ECC2') MX_TERM=2
           IF (CCFORM(1:4).EQ.'ECC3') MX_TERM=3
           IF (CCFORM(1:4).EQ.'ECC4') MX_TERM=4
           IF (CCFORM(1:4).EQ.'ECC5') MX_TERM=5
           IF (CCFORM(1:4).EQ.'ECC6') MX_TERM=6
           IF (CCFORM(1:4).EQ.'ECC7') MX_TERM=7
           IF (CCFORM(1:4).EQ.'ECC8') MX_TERM=8
           IF (CCFORM(1:4).EQ.'ECC9') MX_TERM=9
           
c           CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LULAMP)
           CALL EXPT_REF2(LUC,LU_LXRF,LUVS1,LUVS3,LUVS4,THRES_E,MX_TERM,
     &          CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
           

         END IF

*.  E(ECC) = <ref|exp(L^+) exp(-T)Hexp(T)|ref>:
         E_CC = INPRDD(VEC1,VEC2,LU_LXRF,LUVS2,1,LBLK)

       ELSE IF(CCFORM(1:3).EQ.'VCC') THEN

         XNORM = INPRDD(VEC1,VEC2,LU_EXRF,LU_EXRF,1,LBLK)
         E_CC  = INPRDD(VEC1,VEC2,LU_EXRF,LUHEXRF,1,LBLK)/XNORM
         FAC1 = 1.0D0/XNORM
         FAC2 = (-E_CC)/XNORM

         CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUHEXRF,LUVS2,LUVS3,1,LBLK)

*. VCC : Densities with <0!Exp(T)+ and (H-E) Exp(T) !0> in CC_VEC
         CALL DEN_GCC(VEC1,VEC2,LU_EXRF,LUVS3,CC_VEC)
       ELSE IF(CCFORM(1:4).EQ.'UCC2') THEN
* UCC2
*     we still miss T^2 |0> ....
         ! to be sure
         CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LU_AMP)
         CALL SIG_GCC_U(VEC1,VEC2,LUC,LUVS3,LUVS1,LUVS2,CC_AMP,CC_VEC)
         CALL SIG_GCC_U(VEC1,VEC2,LUVS3,LUVS4,LUVS1,LUVS2,CC_AMP,CC_VEC)

*     ... H |0> ...
         CALL MV7(VEC1,VEC2,LUC,LUVS1,0,0)

*     ... and T H |0> 
         CALL SIG_GCC_U(VEC1,VEC2,LUVS1,LUVS2,LUVS6,LUVS5,CC_AMP,CC_VEC)

c T|0> on LUVS3
c TT|0> on LUVS4
c TH|0> on LUVS2
c  H|0> on LUVS1
         
* E = <0|(1+T^+) H (1+T)|0> + 1/2<0|H T^2|0> + 1/2<0|T^2 H|0>
         E_CC = INPRDD(VEC1,VEC2,LU_EXRF,LUHEXRF,1,LBLK) +
     &          INPRDD(VEC1,VEC2,LUVS1,LUVS4,1,LBLK)

* dE/dt =   <0|(1+T^+) H (tau-tau^+) |0> 
*         - <0|(tau-tau^+) H (1+T)|0> 
*         + 1/2(  <0| H T (tau-tau^+)|0>
*               + <0|(tau-tau^+) T H |0>
*               + <0| H (tau-tau^+)T |0>
*               + <0| T (tau-tau^+)H |0>   )

         CALL DEN_GCC_S(VEC1,VEC2,LUHEXRF,LUC,CC_VEC,CC_AMP,-1)
         CALL VEC_TO_DISC(CC_VEC,N_CC_AMP,1,LBLK,LUVS4)

         CALL DEN_GCC_S(VEC1,VEC2,LUVS2,LUC,CC_VEC,CC_AMP,1)
         CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUVS4)
         CALL VECSUM(CC_VEC,CC_VEC,CC_AMP,-1D0,2D0,N_CC_AMP)
         CALL VEC_TO_DISC(CC_VEC,N_CC_AMP,1,LBLK,LUVS4)

         CALL DEN_GCC_S(VEC1,VEC2,LUVS1,LUVS3,CC_VEC,CC_AMP,1)
         CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,LBLK,LUVS4)
         CALL VECSUM(CC_VEC,CC_VEC,CC_AMP,-1D0,1D0,N_CC_AMP)

       ELSE IF(CCFORM(1:3).EQ.'UCC') THEN
         E_CC = INPRDD(VEC1,VEC2,LU_EXRF,LUHEXRF,1,LBLK)
         
         CALL RELUNIT(LUVS1,'delete')
         CALL RELUNIT(LUVS2,'delete')
         CALL RELUNIT(LUVS3,'delete')
         CALL RELUNIT(LUVS4,'delete')
         CALL RELUNIT(LUVS5,'delete')
         CALL RELUNIT(LUVS6,'delete')

         ! gradient using trick 17 (Wilson formula)
         CALL UCC_GRAD(CC_VEC,CC_AMP,
     &                 VEC1,VEC2,N_CC_AMP,
     &                 LU_AMP,LU_EXRF,LUHEXRF)

       END IF ! CCFORM

       IF (CCFORM(1:3).NE.'UCC'.OR.CCFORM(1:4).EQ.'UCC2') THEN
         CALL RELUNIT(LUVS1,'delete')
         CALL RELUNIT(LUVS2,'delete')
         CALL RELUNIT(LUVS3,'delete')
         CALL RELUNIT(LUVS4,'delete')
         CALL RELUNIT(LUVS5,'delete')
         CALL RELUNIT(LUVS6,'delete')
       END IF

*
* ================================
* 7    Modifications for CC3
* ================================
*
       IF(I_DO_CC3.EQ.1) THEN 
C             CC3_VECFNC(CCVEC,CCDIA,CCAMP,VEC1,VEC2)
         CALL CC3_VECFNC(CC_VEC,WORK(KDIA),CC_AMP,VEC1,VEC2) 
       END IF
       X = INPROD(CC_VEC,CC_VEC,N_CC_AMP)
       IF(NTEST.GE.5) WRITE(6,*) ' Norm of CC vector function =',X

      END IF
*     ^ End of switch between new and old CC vector function
      IF(I_DO_MASK_CC.EQ.1) THEN
*. set unwanted elements of CCvectorfunction to zero
C?      WRITE(6,*) ' CC vector function before masking '
C?      CALL WRT_CC_VEC2(CC_VEC,6,CCTYPE)
        VALUE = 0.0
        CALL MASK_CCVEC(WORK(KLSOBEX),NSPOBEX_TP,CC_VEC,1,
     &       MASK_SD,MSK_AEL,MSK_BEL,VALUE,MX_ST_TSOSO_BLK_MX)
      END IF
*     ^ End if masking was required

      CALL VEC_TO_DISC(CC_VEC,N_CC_AMP,1,LBLK,LU_VECF)
      XVECNRM = SQRT(INPROD(CC_VEC,CC_VEC,N_CC_AMP))

*
*. Mission over, print 
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' CC-vector function '
        WRITE(6,*) ' ==================='
        WRITE(6,*)  
        CALL WRT_CC_VEC2(CC_VEC,6,CCTYPE)
      END IF
C?    WRITE(6,*) ' I12 at END of CC_VEC = ', I12
*. 
C     WRITE(6,*) ' Enforced stop after CC_VEC_FNC'
C     STOP ' Enforced stop after CC_VEC_FNC'

      CALL ATIM(CPU,WALL)
      CALL PRTIM(6,'time in vector-function',cpu-cpu0,wall-wall0)

      RETURN
      END

      SUBROUTINE GET_SX_BLK(HBLK,H,IGAS,ISM,JGAS,JSM)
*
* Fetch block of one-electron excitations from H
*
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION H(*)
*. Output
      DIMENSION HBLK(*)
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
*
      IJ_ACT = I_SX_ACT(IGAS,JGAS)
      IF(IJ_ACT.EQ.0) THEN
*. Just another empty block 
       ZERO = 0.0D0
       CALL SETVEC(HBLK,ZERO,NI*NJ) 
      ELSE
*. Block assumed total symmetric 
       ISX_SYM = 1
       CALL I_OFF_SX(IOFF,IGAS,ISM,JGAS,JSM,ISX_SYM)
       CALL COPVEC(H(IOFF),HBLK,NI*NJ) 
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,'(A,A, 4I3)') 
     &   ' Block of single excitations with sym and type ',
     &   ' (ism,itp,jsm,jtp)', ISM,IGAS,JSM,JGAS
         CALL WRTMAT(HBLK,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END 
      SUBROUTINE I_OFF_SX(IOFF,IIGAS,IISM,JJGAS,JJSM,ISX_SYM)
*
* Offset for single excitation block
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
* 
      IJ_OFF = 1
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,ISX_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
C                 I_SX_CCACT(IGAS,JGAS)
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            IF(ISM.EQ.IISM.AND.IGAS.EQ.IIGAS.AND.
     &         JSM.EQ.JJSM.AND.JGAS.EQ.JJGAS     ) THEN 
               IOFF = IJ_OFF
            END IF
            IJ_OFF = IJ_OFF + NI*NJ
          END IF
        END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Offset for sx block (gas and sym for i and j )',
     &             IISM,IIGAS,JJSM,JJGAS ,' is ', IOFF
      END IF
*
      RETURN
      END
      SUBROUTINE GET_DX_BLK(IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM, 
     &                      C,CBLK,IEXP,IXCHNG,IKLJ,IKSM,JLSM,SCR,
     &                      IJ_TRNSP)
*
* Fetch block of double excitation coefficients
*
* IXCNG : Obtain coulom - exchange
* IKLJ  : = 1 => Obtain in dirac form <IK!LJ> 
*         = 0 => Obtain in Coulomb form (IJ!KL) 
*
* IJ_TRNSP = 1 : Block of interest is C(ji,kl) in form CBLK(IJ,KL)
*
* If the integrals are exported in dirac form 
* there is the additional possibility of supplying integrals in
* symmetric forms 
* IKSM : use i.ge.k
* IKSM : use j.ge.l
*
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION C(*)
*.output 
      DIMENSION CBLK(*)
*. Scratch
      DIMENSION SCR(*) 
*. General information 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
*. Double excitations assumed symmetrix
      IDX_SYM = 1
*
*. Obtain C_{ijkl} in usual (IJ!KL) form
*
      CALL FETCH_DX_BLK(IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,C,CBLK,
     &                  IDX_SYM,IEXP,IJ_TRNSP)
      IF(IXCHNG.EQ.1) THEN
*. Obtain C_{iljk} 
        CALL FETCH_DX_BLK(IGAS,ISM,LGAS,LSM,KGAS,KSM,JGAS,JSM,C,SCR ,
     &                    IDX_SYM,IEXP,IJ_TRNSP)
*. C(I,J,K,L) = C_{ijkl} - C_{ilkj}
        CALL EXCHN_TO_COUL(CBLK,SCR,NI,NJ,NK,NL,1)
C            EXCHN_TO_COUL(C_COUL,C_EXCHN,NI,NJ,NK,NL,ISC)
      END IF
*. Into Dirac form ?
      IF(IKLJ.EQ.1) THEN
        CALL REO_DXBLK_MUL_TO_DIR(NI,NJ,NK,NL,IKSM,JLSM,CBLK,SCR,1)
      END IF
*
      NTEST = 00 
      IF(ISM+JSM+KSM+LSM.EQ.16) THEN
        NTEST = 00
      ELSE
        NTEST = 0
      END IF
      IF(NTEST.GE.100) THEN 
*
        WRITE(6,*)
        WRITE(6,*) ' ============================='
        WRITE(6,*) ' Output block from GET_DX_BLK '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
*
        IF(IKLJ.EQ.0) THEN
          WRITE(6,*) ' ISM JSM JSM LSM ', ISM,JSM,KSM,LSM
          WRITE(6,*) ' Coulomb form C(IJ,KL) '
          WRITE(6,*)
          CALL WRTMAT(CBLK,NI*NJ,NK*NL,NI*NJ,NK*NL)
        ELSE IF (IKLJ .EQ. 1 ) THEN
          WRITE(6,*) ' Exchange form C(IK,JL) '
          WRITE(6,*)
          IF(IKSM.EQ.1) THEN
            NIK = NI*(NI+1)/2 
          ELSE
            NIK = NI*NK
          END IF
          IF(JLSM.EQ.1) THEN
            NJL = NJ*(NJ+1)/2
          ELSE
            NJL = NJ*NL
          END IF
          CALL WRTMAT(CBLK,NIK,NJL,NIK,NJL)
        END IF
*
      END IF
*
      RETURN
      END 
      SUBROUTINE REO_DXBLK_MUL_TO_DIR(NI,NJ,NK,NL,IKSM,JLSM,DIN,DOUT,
     &                                ICOPY)
*
* Reorganize block of 4-electron terms from Mullikan to Dirac form 
*
* If Icopy .eq. 1, then output block is copied over inputblock
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input D(I,J,K,L)
      DIMENSION DIN(NI,NJ,NK,NL)
*. output D(IK,JL)
      DIMENSION DOUT(*)
*
      IF(IKSM.EQ.0) THEN
        NIK = NI*NK
      ELSE
        NIK = NI*(NI+1)/2
      END IF
*
      IF(JLSM.EQ.0) THEN
        NJL = NJ*NL
      ELSE
        NJL = NJ*(NJ+1)/2
      END IF
*
      DO L = 1, NL
       IF(JLSM.EQ.1) THEN
         JMIN = L
       ELSE
         JMIN = 1
       END IF
       DO J = JMIN, NJ
        IF(JLSM.EQ.0) THEN
         JL = (L-1)*NJ + J
        ELSE
C        JL = J*(J-1)/2+L
         JL = (L-1)*NJ + J - L*(L-1)/2
        END IF
        DO K = 1, NK
         IF(IKSM.EQ.1) THEN
          IMIN = K
         ELSE
          IMIN = 1
         END IF
         DO I = IMIN,NI
          IF(IKSM.EQ.0) THEN
           IK = (K-1)*NI+I
          ELSE
C          IK = I*(I-1)/2+K
           IK = (K-1)*NI + I - K*(K-1)/2
          END IF
*
          DOUT((JL-1)*NIK+IK) = DIN(I,J,K,L)
*
         END DO
        END DO
       END DO
      END DO
*
      IF(ICOPY.EQ.1) THEN
       CALL COPVEC(DOUT,DIN,NIK*NJL)
      END IF
*
C     WRITE(6,*) ' NI NJ NK NL NIK NJL', NI,NJ,NK,NL,NIK,NJL
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Double coefficients in DIRAC format D(IK,JL)'
        WRITE(6,*)
        CALL WRTMAT(DOUT,NIK,NJL,NIK,NJL)
      END IF
*
      RETURN
      END 
      SUBROUTINE EXCHN_TO_COUL(C_COUL,C_EXCHN,NI,NJ,NK,NL,ISC)
*
* Form exchange to coulomb format for block of doubles coefs
*
* ISC = 1 C_COUL(I,J,K,L) =  C_COUL(I,J,K,L) - C_EXCHN(I,L,K,J)
* ISC = 2 C_COUL(I,J,K,L) =  C_EXCHN(I,L,K,J)
*
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION C_COUL(NI,NJ,NK,NL)
*. Output
      DIMENSION C_EXCHN(NI,NL,NK,NJ)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The two input matrices to EXCHN.... '
        CALL WRTMAT(C_COUL ,NI*NJ,NK*NL,NI*NJ,NK*NL)
        CALL WRTMAT(C_EXCHN,NI*NL,NK*NJ,NI*NL,NK*NJ)
      END IF
*
      DO I = 1, NI
       DO J = 1, NJ
        DO K = 1, NK
         DO L = 1, NL
           IF(ISC.EQ.2) THEN
             C_COUL(I,J,K,L) = C_EXCHN(I,L,K,J)
           ELSE
             C_COUL(I,J,K,L) =  C_COUL(I,J,K,L) - C_EXCHN(I,L,K,J)
           END IF
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from COUL_TO_EXCHN '
        WRITE(6,*)
        CALL WRTMAT(C_COUL,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END 
      SUBROUTINE FETCH_DX_BLK(IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
     &                        C,CBLK,IDX_SYM,IEXP,IJ_TRNSP)
*
* Fetch block with given type and sym from list of double excitations
*
* If block is packed and IEXP = 1, then the block is expanded   
*
* IF IJ_TRNSP = 1, then we are after the doubles block C(ji,kl) stored 
* as Cblk(ij,kl). Introduced to comply with call from RSBB2BN 
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input
      DIMENSION C(*)
*. Output
      DIMENSION CBLK(*)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Local scratch (dirty, but simple)
      DIMENSION SCR(200,200)
*
C?    WRITE(6,*) ' First 5 elements of C in FETCH_DX.. '
C?    CALL WRTMAT(C,1,5,1,5)
      IF(IJ_TRNSP.EQ.1) THEN
        L = IGAS
        IGAS = JGAS
        JGAS = L
        L = ISM
        ISM = JSM
        JSM = L
      END IF
*. Is combination of types allowed
      IJKL_ACT = I_DX_ACT(IGAS,KGAS,LGAS,JGAS)
*. Complete or packed block
      IF(IGAS.EQ.KGAS.AND.ISM.EQ.KSM .AND. 
     &   JGAS.EQ.LGAS.AND.JSM.EQ.LSM       ) THEN
        ISYM = 1
      ELSE
        ISYM = 0
      END IF
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
      IF(ISYM.EQ.0) THEN
        LEN = NI*NJ*NK*NL
      ELSE
        LEN = (NI*NJ+1)*NI*NJ/2
      END IF
      IF(IJKL_ACT.EQ.0) THEN
*. Zero block
        ZERO = 0.0D0
        CALL SETVEC(CBLK,ZERO,LEN)          
      ELSE
*. Obtain offset for block
        CALL I_OFF_DX(IOFF,ITRNSP,
     &         IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,IDX_SYM) 
        IF(IEXP.EQ.0.OR.ISYM.EQ.0) THEN
          IF(ITRNSP.EQ.0) THEN
            CALL COPVEC(C(IOFF),CBLK(1),LEN)         
          ELSE IF(ITRNSP.EQ.1 ) THEN
            CALL TRPMAT(C(IOFF),NK*NL,NI*NJ,CBLK(1) )
          END IF
        ELSE 
          SIGN = 1.0D0
          CALL TRIPK2(CBLK(1),C(IOFF),2,NI*NJ,NI*NJ,SIGN)
        END IF
      END IF
*
      IF(IJ_TRNSP.EQ.1) THEN
* C(JI,KL) obtained, transpose to C(IJ,KL) 
       DO KL = 1, NK*NL
         KLOFF = 1 + (KL-1)*NI*NJ
         CALL TRPMAT(CBLK(KLOFF),NI,NJ,SCR)
         CALL COPVEC(SCR,CBLK(KLOFF),NI*NJ)
        END DO
*. Clean up     
        L = IGAS
        IGAS = JGAS
        JGAS = L
        L = ISM
        ISM = JSM
        JSM = L
      END IF
*
      NTEST = 00
      IF(ISM+JSM+KSM+LSM.EQ.16) THEN
        NTEST = 00
      ELSE
        NTEST = 0  
      END IF
*
      IF(NTEST.GE.100) THEN  
        WRITE(6,*) ' Block of double-coefficients '
        IF(IJ_TRNSP.EQ.0) THEN
          WRITE(6,*) ' Type and symmetry of indeces (I,J,K,L) ' ,
     &                 IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
        ELSE IF(IJ_TRNSP.EQ.1) THEN
          WRITE(6,*) ' Type and symmetry of indeces (J,I,K,L) ' ,
     &                 JGAS,JSM,IGAS,ISM,KGAS,KSM,LGAS,LSM
         WRITE(6,*) ' (Transposed in indeces I and J ) '
        END IF
        CALL WRTMAT(CBLK,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END
      SUBROUTINE I_OFF_DX(IOFF,ITRNSP,
     &           IXGAS,IXSM,JXGAS,JXSM,KXGAS,KXSM,LXGAS,LXSM,IDX_SYM)
*
* Obtain offset for  block of double excitations. 
*
* Coefficients are stored in blocks C(ai,bj), (ij).ge.(kl) 
*
* If ijkl corresponds to (ij).lt.(kl), the ITRNSP flag
* is set, and the offset to the corresponding klij block is returned
* Blocks with (kl).gt. (ij) are referred to the corresponding 
* klij block, and the transpose flag is set.
*
*
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
*
      IIGAS = IXGAS
      IISM  = IXSM
*
      JJGAS = JXGAS
      JJSM  = JXSM
*
      KKGAS = KXGAS
      KKSM  = KXSM
*
      LLGAS = LXGAS
      LLSM  = LXSM
*. C(ij,kl) or C(kl,ij) ( Transposed)
      ITRNSP = 0
      IJSM_INDEX = (IISM-1)*NSMOB + JJSM
      KLSM_INDEX = (KKSM-1)*NSMOB + LLSM
      IJGAS_INDEX = (IIGAS-1)*NGAS+ JJGAS
      KLGAS_INDEX = (KKGAS-1)*NGAS+ LLGAS
*. Transpose if   
      IF(IJSM_INDEX.LT.KLSM_INDEX .OR.
     &   IJSM_INDEX.EQ.KLSM_INDEX.AND.IJGAS_INDEX.LT.KLGAS_INDEX) THEN
        ITRNSP = 1
*
        IIGAS = KXGAS
        IISM  = KXSM
*
        JJGAS = LXGAS
        JJSM  = LXSM
*
        KKGAS = IXGAS
        KKSM  = IXSM
*
        LLGAS = JXGAS
        LLSM  = JXSM
*
      END IF
*
      IJBL_GT_KLBL = 0  
      IJKL_OFF = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
C?     WRITE(6,*) ' ISM JSM KSM LSM', ISM,JSM,KSM,LSM
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1  
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1  
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0  
            ELSE
              IJBL_GT_KLBL = -1 
            END IF
          END IF
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            IF(IIGAS.EQ.IGAS.AND.IISM.EQ.ISM.AND.
     &         JJGAS.EQ.JGAS.AND.JJSM.EQ.JSM.AND.
     &         KKGAS.EQ.KGAS.AND.KKSM.EQ.KSM.AND.
     &         LLGAS.EQ.LGAS.AND.LLSM.EQ.LSM ) THEN
               IOFF = IJKL_OFF
            END IF
*
            IF(IJBL_GT_KLBL.EQ.1) THEN
              IJKL_OFF = IJKL_OFF + NI*NJ*NK*NL
            ELSE IF (IJBL_GT_KLBL .EQ. 0) THEN
              IJKL_OFF = IJKL_OFF + NI*NJ*(NI*NJ+1)/2
            END IF
C?          WRITE(6,*) ' IGAS,JGAS,KGAS,LGAS,IJKL_OFF',
C?   &                   IGAS,JGAS,KGAS,LGAS,IJKL_OFF
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Offset for block of double excitations' 
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        WRITE(6,*) ' Type and Sym for orbitals I,J,K,L ', 
     &               IIGAS,IISM,JJGAS,JJSM,KKGAS,KKSM,LLGAS,LLSM  
        WRITE(6,*) ' IDX_SYM and IOFF = ', IDX_SYM,IOFF
      END IF
*
      RETURN
      END
      SUBROUTINE RENRM_DX(C,IWAY,IDX_SYM)
*
* Switch normalization of double excitations corresponding
* to switching between restricted and unrestricted summation 
*
* Restricted   summation : T2 = sum(ij.ge.kl) CR_{ijkl} e_{ijkl}
* Unrestricted summation : T2 = sum(ij    kl) CU_{ijkl} e_{ijkl}
*
* Relation between restricted and unrestricted summation
*
* CU(ijkl) = (1+delta((ij),kl))/2 CR(ijkl)
*
* IWAY = 1 restricted to unrestricted
* IWAY = 2 unrestricted to restricted

* Jeppe Olsen, Summer of 98
*
*. Input and output
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc' 
*. Specific input/output
      DIMENSION C(*)
*
      IJBL_GT_KLBL = 0  
      IF(IWAY.EQ.1) THEN
        FACTOR = 0.5D0
      ELSE
        FACTOR = 2.0D0
      END IF
      FACTORI = 1.0D0/FACTOR
*
      IJKL_OFF = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1  
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1  
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0  
            ELSE
              IJBL_GT_KLBL = -1 
            END IF
          END IF
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
            NK = NOBPTS(KGAS,KSM)
            NL = NOBPTS(LGAS,LSM)
*
            WRITE(6,*) ' REN.. ACTIVE IGAS,JGAS,KGAS,LGAS,IJKL_OFF',
     &      IGAS,JGAS,KGAS,LGAS,IJKL_OFF
            IF(IJBL_GT_KLBL.EQ.1) THEN
              CALL SCALVE(C(IJKL_OFF),FACTOR,NI*NJ*NK*NL)
              IJKL_OFF = IJKL_OFF + NI*NJ*NK*NL
            ELSE IF (IJBL_GT_KLBL .EQ. 0) THEN
              CALL SCALVE(C(IJKL_OFF),FACTOR,(NI*NJ+1)*NI*NJ/2)
              CALL SCLDIA(C(IJKL_OFF),FACTORI,NI*NJ,1)
C                  SCLDIA(A,FACTOR,NDIM,IPACK)
              IJKL_OFF = IJKL_OFF + NI*NJ*(NI*NJ+1)/2
            END IF
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ================================='
       WRITE(6,*) ' Renormalized double coefficients '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       WRITE(6,*) ' Symmetry = ', IDX_SYM
       CALL WRT_DX1(C,IDX_SYM)
      END IF
*
      RETURN
      END
      FUNCTION I_DX_ACT(IGAS,KGAS,LGAS,JGAS)
*
* Is double excitation a+ igas a+ kgas a lgas a jgas active
*
* Hiding cc restrictions etc
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'cc_exc.inc'
*
      IF(ICC_EXC.EQ.1) THEN
*. Check for coupled cluster restrictions
        I_DX_ACT = I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
      ELSE
*. Normal CI, extensions for perturbation etc can be inseted here 
        I_DX_ACT = 1
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(I_DX_ACT.EQ.1) THEN 
          WRITE(6,*) 
     &    ' allowed excitation a+(igas) a (jgas) a+(kgas) a (lgas) for'
     &    ,IGAS,JGAS,KGAS,LGAS
        ELSE IF(I_DX_ACT.EQ.0) THEN 
          WRITE(6,*) 
     &    ' excluded excitation a+(igas) a (jgas) a+(kgas) a (lgas) for'
     &    ,IGAS,JGAS,KGAS,LGAS
        END IF
      END IF
*
      RETURN
      END
      FUNCTION I_SX_ACT(IGAS,JGAS)
*
* Is single excitation a+ igas a jgas active
*
* Hiding cc restrictions etc
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'cc_exc.inc'
*
      IF(ICC_EXC.EQ.1) THEN
*. Check for coupled cluster restrictions
        I_SX_ACT = I_SX_CCACT(IGAS,JGAS)
      ELSE
*. Normal CI, extensions for perturbation etc can be inseted here 
        I_SX_ACT = 1
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(I_SX_ACT.EQ.1) THEN 
          WRITE(6,*) 
     &    ' allowed excitation a+ (igas) a (jgas) for igas,jgas=',
     &    IGAS,JGAS
        ELSE IF(I_SX_ACT.EQ.0) THEN 
          WRITE(6,*) 
     &    ' exluded excitation a+ (igas) a (jgas) for igas,jgas=',
     &    IGAS,JGAS
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_DX1(C,IDX_SYM)
*
*
* Print list of double excitations in compressed form
* - without singlet-singlet, triplet-triplet separation
* Jeppe Olsen, summer of 98
*

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
*. Input 
      DIMENSION C(*)
*
      IJBL_GT_KLBL = 0  
*
C     WRITE(6,*)
C     WRITE(6,*) ' ================================================'
C     WRITE(6,*) ' List of double excitations without S/T splitting'
C     WRITE(6,*) ' ================================================'
C     WRITE(6,*) 
      WRITE(6,*) 
     & '  (Blocks for E(IJ) E(KL) written as matrices C(IJ,KL) )'
      IJKL_OFF = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1  
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1  
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0  
            ELSE
              IJBL_GT_KLBL = -1 
            END IF
          END IF
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            IF(NI*NJ*NK*NL.GT.0) THEN
              WRITE(6,'(A,8I3)') 
     &        ' Orbital indeces I,J,K,L (type and sym)',
     &          IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM 
              IF(IJBL_GT_KLBL.EQ.1) THEN
                CALL WRTMAT(C(IJKL_OFF),NI*NJ,NK*NL,NI*NJ,NL*NL)
                IJKL_OFF = IJKL_OFF + NI*NJ*NK*NL
              ELSE IF (IJBL_GT_KLBL .EQ. 0) THEN
                CALL PRSM2(C(IJKL_OFF),NI*NJ)
                IJKL_OFF = IJKL_OFF + NI*NJ*(NI*NJ+1)/2
              END IF
            END IF
*           ^ End if nonvanishing block
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      RETURN
      END
      SUBROUTINE REF_DX_EXP_COMP(CEXP,CCOM,IDX_SYM,IWAY,IBIO)
*
* Reform double excitations between expanded form and
* compressed form - without singlet-triplet separation
*
* IWAY = 1 : Expanded to compressed form
* IWAY = 2 : Compressed to expanded form
*
* Jeppe Olsen, summer of 98
*
*. Modified to column wise packing
*
*. Assumes IBIO = 1
*
* Diagonal elements are divided with two.
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'lucinp.inc'
*. Input and output
      DIMENSION CEXP(*),CCOM(*)
*
      FACTOR1 = 1.0D0/3.0D0
      FACTOR2 = 1.0D0/6.0D0
      HALF = 0.5D0
*
C?    WRITE(6,*) ' IBIO = ', IBIO
*
      IJKL_COM = 0
      IJBL_GT_KLBL = 0  
      I_MIN = 0
*. To eliminate warning 
      ILKJ_EXP = 0
*
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
C?     WRITE(6,*) ' ISM,JSM,KSM,LSM,IJSM_GT_KLSM',
C?   &              ISM,JSM,KSM,LSM,IJSM_GT_KLSM
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
C?        WRITE(6,*) ' IGAS,JGAS,KGAS,LGAS,IJKL_ACT',
C?   &                 IGAS,JGAS,KGAS,LGAS,IJKL_ACT
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1  
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1  
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0  
            ELSE
              IJBL_GT_KLBL = -1 
            END IF
          END IF
C?        WRITE(6,*) ' IJBL_GT_KLBL' ,  IJBL_GT_KLBL
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            DO L = L_OFF,L_OFF+NL-1
            DO K = K_OFF,K_OFF+NK-1
C
            IF(IJBL_GT_KLBL .EQ. 0 ) THEN
             J_MIN = L
            ELSE
             J_MIN = J_OFF
            END IF
            DO J = J_MIN,J_OFF+NJ-1
*
            IF( IJBL_GT_KLBL .EQ. 1 ) THEN
              I_MIN = I_OFF
            ELSE IF ( IJBL_GT_KLBL .EQ. 0 ) THEN
              IF(J.GT.L) THEN
                I_MIN = I_OFF 
              ELSE
                I_MIN = K
              END IF
            END IF
            DO I = I_MIN,I_OFF+NI-1
*
              IJ = (J-1)*NTOOB+I
              KL = (L-1)*NTOOB+K
*
              IF(IBIO.EQ.1) THEN
               IL = (L-1)*NTOOB + I
               KJ = (J-1)*NTOOB + K
               ILKJ_EXP = MAX(IL,KJ)*(MAX(IL,KJ)-1)/2 + MIN(IL,KJ)
              END IF
              IJKL_EXP = IJ*(IJ-1)/2+KL
              IJKL_EXP = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
              IJKL_COM = IJKL_COM + 1   
C?            WRITE(6,*) 'I,J,K,L,IJ,KL,IJKL_EXP,IJKL_COM',
C?   &                    I,J,K,L,IJ,KL,IJKL_EXP,IJKL_COM        
              IF(IWAY.EQ.1) THEN
*. expanded to compressed
                IF(IBIO.EQ.0) THEN
                  CCOM(IJKL_COM) = CEXP(IJKL_EXP)
                  IF(IJ.EQ.KL) CCOM(IJKL_COM) = HALF*CCOM(IJKL_COM)
                ELSE
                  CCOM(IJKL_COM) = FACTOR1*CEXP(IJKL_EXP) +
     &                             FACTOR2*CEXP(ILKJ_EXP)
                  IF(IJ.EQ.KL) CCOM(IJKL_COM) = HALF*CCOM(IJKL_COM)
                END IF
              ELSE
                CEXP(IJKL_EXP) = CCOM(IJKL_COM)
              END IF
*
            END DO
            END DO
            END DO
            END DO
*           ^ End of loop over orbitals over given TS
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      RETURN
      END
      SUBROUTINE WRT_SX(C,ISX_SYM)
*
* Write list of single excitations given in  compact form
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
*. Input 
      DIMENSION C(*)
* 
      WRITE(6,*)
      WRITE(6,*) ' List of single excitations in compressed form '
      WRITE(6,*) ' ============================================= '
      WRITE(6,*)
      IJ_OFF = 1
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,ISX_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
C                 I_SX_CCACT(IGAS,JGAS)
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            IF(NI*NJ.GT.0) THEN
              WRITE(6,'(A,A, 4I3)') 
     &        ' Block of single excitations with sym and type ',
     &        ' (ism,itp,jsm,jtp)', ISM,IGAS,JSM,JGAS
              CALL WRTMAT(C(IJ_OFF),NI,NJ,NI,NJ)
              IJ_OFF = IJ_OFF + NI*NJ
            END IF
*           ^ End if nonvanishing block
          END IF
*         ^ End if active block
        END DO
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE REF_SX(CIN,COUT,INFRM,IOUTFRM,ISX_SYM,IBIO)
*
* Reform single excitations
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
*. Input and output
      DIMENSION CIN(*),COUT(*)
*
      IF(INFRM.EQ.1.AND.IOUTFRM.EQ.2) THEN
*. restricted to unrestricted
        IWAY = 1
      ELSE 
        IWAY = 2
      END IF
*
      IF(IBIO.EQ.0) THEN
        FACTOR = 1.0D0
      ELSE
        FACTOR = 0.5D0
      END IF
*
      IJ_COMP = 0
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,ISX_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
C                 I_SX_CCACT(IGAS,JGAS)
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*. Offsets
            I_OFF = IOBPTS(IGAS,ISM)
            J_OFF = IOBPTS(JGAS,JSM)
*. Numbers
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            DO J = J_OFF,J_OFF+NJ-1
            DO I = I_OFF,I_OFF+NI-1
              IJ_COMP = IJ_COMP + 1
              IF(IWAY.EQ.1) THEN
               COUT((J-1)*NTOOB + I) = CIN(IJ_COMP)*FACTOR
              ELSE
               COUT(IJ_COMP) = CIN((J-1)*NTOOB + I)*FACTOR
              END IF
            END DO
            END DO
          END IF
        END DO
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE SD_TO_EXC(IREFA,IREFB,IEXCA,IEXCB,NAEL,NBEL,
     &                     NAEXC,IAEXC,NBEXC,IBEXC)
*
* Transfer between slater determinant representation and
* excitation representation
*
* Jeppe Olsen, summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IREFA(NAEL),IEXCA(NAEL)
      INTEGER IREFB(NBEL),IEXCB(NBEL)
*.Output
      INTEGER IAEXC(*),IBEXC(*)
*. For alpha string
*. COding interruptus
      RETURN
      END
      SUBROUTINE ST_TO_EXC(IREF,IEXC,NEL,NEXC,IEXC_OP)
*
* Find excitation that Ex | IREF > = | IEXC >
*
*. Jeppe Olsen, Summer of 98
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IREF(NEL),IEXC(NEL)
*. Output 
      INTEGER IEXC_OP(*)
*. Scratch
      INCLUDE 'mxpdim.inc'
      INTEGER IORB(MXPORB),IEL(MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reference string '
        CALL IWRTMA(IREF,1,NEL,1,NEL)
        WRITE(6,*) ' Excited string '
        CALL IWRTMA(IEXC,1,NEL,1,NEL)
      END IF
*
*. Strings are assumed in ascending order
*. Find orbital occuring only in one string. 
      JEEL = 1
      JDIF = 1
      DO JREL = 1, NEL
       IF(IREF(JREL).GT.IEXC(JEEL)) THEN
*. IEXC(JEEL) is only in IEXC
        IORB(JDIF) = IEXC(JEEL)
        IEL (JDIF) = JEEL
        JEEL = JEEL + 1
        JDIF = JDIF + 1
       ELSE IF(IREF(JREL).LT.IEXC(JEEL) ) THEN
*. IREF(JREL) is only in IREF
        IORB(JDIF) = -IREF(JREL)
        IEL(JDIF) = JREL
        JDIF = JDIF + 1
       END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of differences detected ', JDIF
        WRITE(6,*) ' IORB and IEL arrays '
        CALL IWRTMA(IORB,1,JDIF,1,JDIF)
        CALL IWRTMA(IEL ,1,JDIF,1,JDIF)
      END IF
*. We will write the excitation operator in descending order
* !IEXC> = +/-a+(IEXC(1)) ... a+(IEXC(JDIF)) 
*             a (IEXC(NDIF+1)) ... a (IEXC(2*JDIF)) ! IREF>
*  with IEXC(I).GT.IEXC(I+1),I=1, NEXC-1,..., NEXC+1, ...., 2*NEXC
*
       ICREA = 0
       IANNI = 0
       DO JEL = 1, 2*NEXC
         IF(IORB(JEL).GT.0) THEN
           ICREA = ICREA + 1
           IEXC_OP(NEXC+1 -ICREA) = IORB(JEL)
         ELSE IF(IORB(JEL).LT.0) THEN
           IANNI = IANNI + 1
           IEXC_OP(2*NEXC+1-IANNI) = ABS(IORB(JEL))
         END IF
       END DO
*. Then, only the sign is missing
*. 1 : sign for moving orbitals to be annihilated outside
*      (and ordering them in ascending order outside )
       ISIGN = 1
       IANNI = 0
       DO JEL = 1, 2*NEXC
         IF(IORB(JEL).LT.0) THEN
           ISIGN = ISIGN*(-1)**(IEL(JEL)-1-IANNI)
           IANNI = IANNI + 1
         END IF
       END DO 
*. 2 : Sign for adding electrons to annihilated list,
*      done from list of creation operators in 
*.1 : Sign for changing moving all differing electrons in 
       NEXC = NDIF/2
*. We have now tabulated all differences
*
* Coding iterruptus
      RETURN
      END
      SUBROUTINE CC_AC_SPACES(ISPC,IREFTYP)
*
* Divide orbital spaces  into 
*  Hole spaces     : Only annihilation allowed
*  particle spaces : Only creation allowed
*  valence  spaces : Both annihilation and creation allowed
*
* Result is stored in IHPVGAS
*
* Division based upon occupation in CI space ISPC
* Used for Coupled Cluster Calculations
*
* Find type of referene state 
*
* IREFTYP = 1 : Closed Shell Hartree-Fock state
* IREFTYP = 2 : Highspin open shell state 
* IREFTYP = 3 : More generel reference space (CAS, RAS, GAS ..)
*
* Jeppe Olsen, Summer of 98 ( not much of an summer !)
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc' 
      INCLUDE 'strinp.inc'
*
      NEL_REF = NELEC(1) + NELEC(2)
C     WRITE(6,*) ' NELEC(1), NELEC(2) ', NELEC(1),NELEC(2)
*
*. To get rid of annoying and incorrect compiler warnings
      NEL_MAX = 0
*
      NHOLE = 0
      NVAL  = 0
      I_NEW_OR_OLD = 1
      IF(I_NEW_OR_OLD.EQ.2) THEN
*
* Old route : Assumes all hole spaces occurs first 
*
      DO IGAS = 1, NGAS
*.. Occupation in CI space 1 ;
*. hole space : doubly occupied =>2
*. part space : un     occupied 
*. val  space : Various occupations  
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS) 
       END IF
*
       IF(IGSOCCX(IGAS,1,ISPC) .EQ. NEL_MAX  .AND.
     &    IGSOCCX(IGAS,2,ISPC) .EQ. NEL_MAX       ) THEN 
*. hole space
          IHPVGAS(IGAS) = 1
          NHOLE = NHOLE + NGSOBT(IGAS)
       ELSE IF(IGAS.GT.1.AND.IGSOCCX(IGAS-1,1,ISPC) .EQ. NEL_REF) THEN
*. Particle space
          IHPVGAS(IGAS) = 2
       ELSE 
*. Valence space
          IHPVGAS(IGAS) = 3
          NVAL = NVAL + NGSOBT(IGAS)
       END IF
*
      END DO
*
      ELSE 
*
* A more recent code without assuming that the hole spaces 
* occur first ( for the QD project) 
       DO IGAS = 1, NGAS
*. Minimum number of electrons in this space
         IF(IGAS.EQ.1) THEN
           NEL_MIN = IGSOCCX(1,1,ISPC)
         ELSE 
           NEL_MIN = IGSOCCX(IGAS,1,ISPC)-IGSOCCX(IGAS-1,2,ISPC)
           NEL_MIN = MAX(0,NEL_MIN)
         END IF
*. Largest number of electrons in this space 
         IF(IGAS.EQ.1) THEN
           NEL_MAX = IGSOCCX(1,2,ISPC)
         ELSE 
           NEL_MAX = IGSOCCX(IGAS,2,ISPC)-IGSOCCX(IGAS-1,1,ISPC)
           NEL_MAX = MIN(NEL_MAX,2*NGSOBT(IGAS))
         END IF
C?       WRITE(6,*) ' IGAS, NEL_MAX, NEL_MIN = ', 
C?   &                IGAS, NEL_MAX, NEL_MIN
         IF(NEL_MAX.EQ.0) THEN
*. particle space
          IHPVGAS(IGAS) = 2
         ELSE IF (NEL_MIN.EQ.2*NGSOBT(IGAS)) THEN
*. Hole space 
          IHPVGAS(IGAS) = 1
          NHOLE = NHOLE + NGSOBT(IGAS)
         ELSE 
*. Valence space 
          IHPVGAS(IGAS) = 3
          NVAL = NVAL + NGSOBT(IGAS)
         END IF
       END DO
      END IF
           
           
      NEL_AL = NELEC(1)
      NEL_BE = NELEC(2)
      IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE) THEN
*. Closed shell Hartree-Fock
        IREFTYP = 1
      ELSE IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE+NVAL.OR.
     &        NEL_BE.EQ.NHOLE.AND.NEL_AL.EQ.NHOLE+NVAL) THEN
*. Highspin openshell 
        IREFTYP = 2
      ELSE 
*. More general, not analyzed in detail p.t. 
        IREFTYP = 3
      END IF 

*. Test open-shell code (with closed-shell reference:
      i_test_oscode_with_csref = 0
      if(i_test_oscode_with_csref.eq.1) then
        if (ireftyp.ne.1) stop 'ireftyp.ne.1'
        ireftyp = 2
        do ii = 1, 40
          write(6,*) 'TEST TEST TEST : open shell test'
        end do
      end if

*
*. IHPVGAS_AB : Hole, particle, valence for alpha and beta orbitals
*
*. Differs from IHPVGAS for high spin open shell case
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) THEN
          IHPVGAS_AB(IGAS,1) = 1
          IHPVGAS_AB(IGAS,2) = 1
        ELSE IF (IHPVGAS(IGAS).EQ.2) THEN 
          IHPVGAS_AB(IGAS,1) = 2
          IHPVGAS_AB(IGAS,2) = 2
        ELSE IF (IHPVGAS(IGAS).EQ.3) THEN 
          IF(IREFTYP.EQ.2.AND.NEL_AL.GT.NEL_BE) THEN
            IHPVGAS_AB(IGAS,1) = 1
            IHPVGAS_AB(IGAS,2) = 2
          ELSE IF(IREFTYP.EQ.2.AND.NEL_AL.LT.NEL_BE) THEN
            IHPVGAS_AB(IGAS,1) = 2
            IHPVGAS_AB(IGAS,2) = 1
          ELSE 
            IHPVGAS_AB(IGAS,1) = 3
            IHPVGAS_AB(IGAS,2) = 3
          END IF
        END IF
      END DO
*
c      IREFTYP = 2
c      DO I = 1, 20
c       WRITE(6,*) ' IREFTYP enforced to 2 '
c      END DO
          
        
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' CC division of orbitals '
        WRITE(6,*) ' ======================= ' 
        WRITE(6,*)
        WRITE(6,*) ' Hole =>1, Part=>2, Val=>3 '
        WRITE(6,*)
        CALL IWRTMA(IHPVGAS,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' CC division for alpha-spinorbitals'
        CALL IWRTMA(IHPVGAS_AB(1,1),1,NGAS,1,NGAS)
*
        WRITE(6,*) ' CC division for betaa-spinorbitals'
        CALL IWRTMA(IHPVGAS_AB(1,2),1,NGAS,1,NGAS)
*
        IF(IREFTYP.EQ.1) THEN
          WRITE(6,*) ' Reference state is closed shell single SD'
        ELSE IF( IREFTYP.EQ.2) THEN
          WRITE(6,*)
     &    ' Reference state is high spin open shell single SD'
        ELSE 
          WRITE(6,*)
     &    ' Reference state is general multireference state '
        END IF
        
      END IF
*
      RETURN
      END
      FUNCTION I_SX_CCACT(IGAS,JGAS)
*
*  Is excitation a+Igas a Jgas active
*
      IMPLICIT REAL*8(A-H,O-Z)
C     COMMON/CC_EXC/ICC_EXC
      INCLUDE 'cc_exc.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*
      IACT = 0
      IF(ICC_EXC .EQ. 0 ) THEN
        IACT = 1
      ELSE IF (ICC_EXC .EQ. 1 ) THEN
        IF(IHPVGAS(IGAS) .GE.2 .AND. IHPVGAS(JGAS) .NE. 2 
     &     .AND..NOT.(IHPVGAS(IGAS).EQ.3.AND.IHPVGAS(JGAS).EQ.3) ) THEN
          IF(I_IAD(IGAS).EQ.2.AND.I_IAD(JGAS).EQ.2) THEN
            IACT = 1
          ELSE
            IACT = 0
          END IF
        ELSE
          IACT = 0
        END IF
      END IF
*
      I_SX_CCACT = IACT
*
      RETURN
      END
      FUNCTION I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
*
*  Is excitation a+Igas a+Kgas a Lgas a Jgas active
*
      IMPLICIT REAL*8(A-H,O-Z)
C     COMMON/CC_EXC/ICC_EXC
      INCLUDE 'cc_exc.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*
      IF(ICC_EXC .EQ. 0 ) THEN
        IACT = 1
      ELSE IF (ICC_EXC .EQ. 1 ) THEN
*. Not allowed to excite into inactive, annihilate from 
*  secondary, furthermore not all four indeces can be active
        IACT = 0
        IF(I_IAD(IGAS).EQ.2.AND.I_IAD(JGAS).EQ.2.AND.
     &     I_IAD(KGAS).EQ.2.AND.I_IAD(LGAS).EQ.2     ) THEN
          IACT = 1
          IF(IHPVGAS(IGAS).EQ.1.OR.IHPVGAS(KGAS).EQ.1.OR.
     &       IHPVGAS(JGAS).EQ.2.OR.IHPVGAS(LGAS).EQ.2) THEN 
             IACT = 0
          END IF
*. 4 indeces in valence space not allowed
          IF(IHPVGAS(IGAS).EQ.3.AND.IHPVGAS(JGAS).EQ.3.AND.
     &       IHPVGAS(KGAS).EQ.3.AND.IHPVGAS(LGAS).EQ.3) THEN
            IACT = 0
          END IF
        END IF
      END IF
*
      I_DX_CCACT = IACT
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' IGAS(c) JGAS(a) KGAS(c) LGAS(a) and IACT',
     &   IGAS,JGAS,KGAS,LGAS,IACT
      END IF
*
      RETURN
      END
      SUBROUTINE FIND_N_CC_AMP(IEXSYM,NSXA,NSXB,NDXAA,NDXBB,NDXAB,
     &                    NSXE,NDXEE)
* 
* Number of coupled cluster amplitudes
*
* Jeppe Olsen, Summer of 98 
*
* Two sets of amplitudes : (Spin-) Restricted and unrestricted.
*. Currently only unrestricted set is implemented
*
* Restricted
* T = sum(ai) C(ai)E(ai) + 1/2 sum(ai>=bj) C(aibj) E(ai) E(bj)

* Unrestricted
* T = sum(ai) Ca(ai) Ea(ai) + sum(ai) Cb(ai)Eb(ai)
*   + sum(a>b,i>j) Caa(abij) Ea(ai)Ea(bj)
*   + sum(a>b,i>j) Cbb(abij) Eb(ai)Eb(bj)
*   + sum(ab,ij )  Cab(abij) Ea(ai)Eb(bj)
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
C     COMMON/CC_EXC/ICC_EXC
      INCLUDE 'cc_exc.inc'
*. ICC_EXC should be set to one before call
*
*. Single excitations
*
      NSX = 0
      DO IGAS = 1, NGAS
       DO JGAS = 1, NGAS
*. Is this excitation allowed
        IACT = I_SX_ACT(IGAS,JGAS)
        IF(IACT.EQ.1) THEN
          DO ISYM = 1, NSMST
            JSYM = MULTD2H(ISYM,IEXSYM)
            NSX = NSX + NGSOB(ISYM,IGAS)*NGSOB(JSYM,JGAS)
          END DO
        END IF
       END DO
      END DO
      NSXA = NSX
      NSXB = NSX
      NSXE = NSX
*
* Double excitations
*
      NDXAA = 0
      NDXAB = 0
      NDXEE = 0
      DO IGAS = 1, NGAS 
       DO JGAS = 1, NGAS
        DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
           IACT = I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
           IF(IACT.EQ.1) THEN
            DO ISYM = 1, NSMST
             DO JSYM = 1, NSMST
              DO KSYM = 1, NSMST
*
               IJSYM = MULTD2H(ISYM,JSYM)
               IJKSYM = MULTD2H(IJSYM,KSYM)
               LSYM = MULTD2H(IJKSYM,IEXSYM)
*
               I_INDEX = (ISYM-1)*NGAS + IGAS
               J_INDEX = (JSYM-1)*NGAS + JGAS
               K_INDEX = (KSYM-1)*NGAS + KGAS
               L_INDEX = (LSYM-1)*NGAS + LGAS
*
               NI = NGSOB(ISYM,IGAS)
               NJ = NGSOB(JSYM,JGAS)
               NK = NGSOB(KSYM,KGAS)
               NL = NGSOB(LSYM,LGAS)
*. Alpha-alpha and beta-beta excitations
               NIK = 0
               NJL = 0
               IF(I_INDEX.GT.K_INDEX) THEN
                 NIK =  NI*NK
               ELSE IF(I_INDEX.EQ.K_INDEX) THEN
                 NIK =  NK*(NK-1)/2
               END IF
               IF(J_INDEX.GT.L_INDEX) THEN
                 NJL = NJ*NL
               ELSE IF(J_INDEX.EQ.L_INDEX) THEN
                 NJL = NJ*(NJ-1)/2
               END IF
               NDXAA = NDXAA + NIK*NJL
*. Alpha-beta excitations
               NDXAB = NDXAB + NI*NJ*NK*NL
*. Restricted EE excitations
               IF(I_INDEX.GT.K_INDEX.OR.
     &            I_INDEX.EQ.K_INDEX.AND.J_INDEX.GT.L_INDEX) THEN
*. Allowed sym/types, no restrictions 
                   NDXEE = NDXEE + NI*NJ*NK*NL
               ELSE IF(I_INDEX.EQ.K_INDEX.AND.J_INDEX.EQ.L_INDEX) THEN
*. Allowed sym/types, restriced sum 
                   NDXEE = NDXEE + (NI*NJ+1)*NI*NJ/2
               END IF
              END DO
             END DO
            END DO
          END IF
         END DO
        END DO
       END DO
      END DO
      NDXBB = NDXAA
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' NSXA,NSXB = ', NSXA,NSXB
       WRITE(6,*) ' NDXAA,NDXBB,NDXAB = ', NDXAA,NDXBB,NDXAB
       WRITE(6,*) ' NDXEE = ', NDXEE
      END IF
*
      RETURN
      END 
*
      SUBROUTINE INI_CC_AMP(LUCCAMP,CC,IFORM)
* 
* Initialize Coupled Cluster amplitudes
*
*. Method of initialization depends upon IFORM
*
* IFORM = 1 => Set to zero
* IFORM = 2 => Read in from LUCCAMP
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
*(    ^contains number of single and double excitations and N_CC_AMP)
      INCLUDE 'clunit.inc'
*. Amplitudes to be defined
      DIMENSION CC(*)
*
      STOP 'do not use this routine no more never ever ...'
      IF(IFORM.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(CC,ZERO,N_CC_AMP)
        WRITE(6,*) ' Initial set of amplitudes set to zero '
*
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
*
      ELSE IF(IFORM.EQ.2) THEN
        ZERO = 0.0D0
        CALL SETVEC(CC,ZERO,N_CC_AMP)
        WRITE(6,*) ' Reading in CC amplitudes from ', LUCCAMP
        I_FORMATTED = 0
        IF(I_FORMATTED.EQ.1) THEN
*. Formatted 
          CALL REWINO(LUCCAMP)
          READ(LUCCAMP,*) N_CC_AMPP
          N_CC_AMP_READ = MIN(N_CC_AMPP,N_CC_AMP)
          DO I = 1, N_CC_AMP_READ
            READ(LUCCAMP,*) CC(I)
          END DO
        ELSE
*. Unformatted
          CALL REWINO(LUCCAMP)
          READ(LUCCAMP) N_CC_AMPP
          N_CC_AMP_READ = MIN(N_CC_AMPP,N_CC_AMP)
          LBLK = -1
          CALL VEC_FROM_DISC(CC,N_CC_AMP_READ,1,LBLK,LUCCAMP)
c          READ(LUCCAMP) (CC(I),I=1,N_CC_AMP_READ)
        END IF
      ELSE
        WRITE(6,*) ' Unknown parameter in INI_CC_AMP ',IFORM
        STOP       ' Unknown parameter in INI_CC_AMP '             
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Initial set of amplitudes '
        WRITE(6,*) ' ========================= '
        WRITE(6,*)
C            WRT_CC_VEC(CC,LU)
        CALL WRT_CC_VEC(CC,6)
      END IF
*
      RETURN
      END
      SUBROUTINE LUCIA_CC(ISM,ISPC,IPRNT,ECC,II_RESTRT_CC,I_TRANS_WF)
*
* Coupled Cluster calculations with LUCIA 
*
* Jeppe Olsen, March 1998
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPRDD
      CHARACTER*6 CCTYPE
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc' 
      INCLUDE 'cands.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
C     COMMON/CC_EXC/ICC_EXC
      INCLUDE 'cc_exc.inc'
*
      stop 'call to obsolete routine'
*
      RETURN
      END
C       EXPT_REF(LUC,LUHC,LUSC1,LUSC2,THRES_E,MX_TERM,
C    &             WORK(KVEC1),WORK(KVEC2))
      SUBROUTINE EXPT_REF(LUC,LUHC,LUSC1,LUSC2,LUSC3,THRES_C,MX_TERM,
     &                    VEC1,VEC2,CCTYPE)
*
* Obtain Exp (T) !ref> by Taylor expansion of exponential
*
* Jeppe Olsen, March 1998
*
* Extended to include general CC, summer of 99
*
* T is in work(KCC1) 
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPRDD
      CHARACTER*6 CCTYPE
*
      INCLUDE 'glbbas.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'crun.inc'
*
      DIMENSION VEC1(*),VEC2(*)
      COMMON/CINT_CC/INT_CC
* 
      LBLK = -1
*
C?    MX_TERM = 10
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF' 
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRCC)
*
      IF(NTEST.GE.5) THEN
       WRITE(6,*)
       WRITE(6,*) '==================='
       WRITE(6,*) 'EXPT_REF in action '
       WRITE(6,*) '==================='
       WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' LUC,LUHC,LUSC1,LUSC2',LUC,LUHC,LUSC1,LUSC2
       WRITE(6,*) ' Initial vector on LUC '
       CALL WRTVCD(VEC1,LUC,1,LBLK)
      END IF
* Tell integral fetcher to fetch cc amplitudes, not integrals
      INT_CC = 1
*. Loop over orders of expansion
      N = 0
      XFACN = 1.0D0
*
      CALL COPVCD(LUC,LUSC1,VEC1,1,LBLK)
      CALL COPVCD(LUC,LUSC3,VEC1,1,LBLK)
*
 1000 CONTINUE
       N = N+1
       IF(NTEST.GE.5) THEN
         WRITE(6,*) ' Info for N = ', N
       END IF
*. T^N  times vector on LUSC1
C?     WRITE(6,*) ' Input vector to MV7 '
C?     CALL WRTVCD(VEC1,LUSC1,1,LBLK)
       IF(CCTYPE(1:2).EQ.'CC') THEN
         CALL MV7(VEC1,VEC2,LUSC1,LUHC,0,0)
       ELSE IF(CCTYPE(1:6).EQ.'GEN_CC') THEN
         XNORM = INPROD(WORK(KCC1),WORK(KCC1),N_CC_AMP)
         WRITE(6,'(A,E22.15)') ' square-norm of T = ', XNORM
         CALL SIG_GCC(VEC1,VEC2,LUSC1,LUHC,WORK(KCC1))
C  SIG_GCC(C,HC,LUC,LUHC,T)
* General CC sigma
       END IF

       CALL COPVCD(LUHC,LUSC1,VEC1,1,LBLK)
       IF(NTEST.GE.500) THEN
         WRITE(6,*) ' T**(N) |0> '
         WRITE(6,*) ' ==========='
         CALL WRTVCD(VEC1,LUSC1,1,LBLK)
       END IF
*. Norm of this correction term
       LBLK = -1
       XNORM2 = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
       XFACN = XFACN/FLOAT(N)
       XNORM = SQRT(XNORM2)/XFACN
       IF(NTEST.GE.5) THEN
         WRITE(6,*) ' Norm of correction ', XNORM
       END IF
*. Update output file with 1/N! T^N !ref>
       ONE = 1.0D0
       CALL VECSMD(VEC1,VEC2,ONE,XFACN,LUSC3,LUSC1,LUSC2,1,LBLK)
       CALL COPVCD(LUSC2,LUSC3,VEC1,1,LBLK)
*. Take another turn ?
      IF(XNORM.GT. THRES_C.AND. N .LT. MX_TERM) GOTO 1000
*
*. Result on LUHC
      IF (XNORM.GT.THRES_C) THEN
        WRITE(6,'(x,a,i5,a)')
     $        'Fatal: No convergence in EXPT_REF (max. iter.:',
     $        MX_TERM, ' )'
        STOP 'No convergence in EXPT_REF!'
      END IF
      CALL COPVCD(LUSC3,LUHC,VEC1,1,LBLK)
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Convergence obtained in ', N, ' iterations'
        WRITE(6,*) ' Norm of last correction ', XNORM
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============'
        WRITE(6,*) ' Exp T |ref> '
        WRITE(6,*) ' ============'
        WRITE(6,*)
        CALL WRTVCD(VEC1,LUHC,1,LBLK)
      END IF
*
      RETURN
      END 
*------------------------------------------------------------------------*
*     clone of EXPT_REF with T-operator as explicit argument on TAMP
*------------------------------------------------------------------------*
      SUBROUTINE EXPT_REF2(LUC,LUHC,LUSC1,LUSC2,LUSC3,THRES_C,MXTERM,
     &                    TAMP,TSCR,VEC1,VEC2,N_CC_AMP,CCTYPE,IOPTYP)
*
* Obtain Exp (T) !ref> by Taylor expansion of exponential
*
* Jeppe Olsen, March 1998 
*
* Extended to include general CC, summer of 99
*
* IOPTYP defines symmetry of operator: 
*
*    +1 Hermitian
*    -1 unitary
*     0 general
*
* MXTERM: If a negative number is provided, exp(T)|R> is expected to be
*         and infinite series that should be approximated to THRSH_C
*         within -MXTERM iterations. If not, the program will cancel with
*         an error message
*         If the value is positive, it is assumed that only a finite expansion
*         up to T^N, N=MXTERM is requested. If the contributions drop below
*         THRSH_C at an earlier T^N, the routine will return before N==MXTERM
*
* TSCR is only needed in the first two cases.
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'

      REAL*8 INPRDD, INPROD
      CHARACTER*6 CCTYPE
*
      INCLUDE 'glbbas.inc'
      INCLUDE 'cprnt.inc'
*
      DIMENSION VEC1(*),VEC2(*),TAMP(*),TSCR(*)
      COMMON/CINT_CC/INT_CC
*
      LOGICAL STRICT
* 
      LBLK = -1
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRCC)
*
      IF (MXTERM.LT.0) THEN
        MXTERM_L = -MXTERM
        STRICT = .TRUE.
      ELSE
        MXTERM_L = MXTERM
        STRICT = .FALSE.
      END IF
*
      IF (IOPTYP.EQ.1) THEN
        SFAC = 1d0
      ELSE IF(IOPTYP.EQ.-1) THEN
        SFAC = -1d0
      ELSE IF (IOPTYP.NE.0) THEN
        WRITE(6,*) 'Indigestible input in EXPT_REF2!!!'
        STOP 'IOPTYP in EXPT_REF2'
      END IF
*
      IF(NTEST.GE.5) THEN
       WRITE(6,*)
       WRITE(6,*) '==================='
       WRITE(6,*) 'EXPT_REF in action '
       WRITE(6,*) '==================='
       WRITE(6,*) ' ioptyp = ',ioptyp
       WRITE(6,*) ' MXTERM = ',MXTERM
       IF (MXTERM.LT.0) THEN
         WRITE(6,*) ' infinite series expansion with max. ',
     &        -MXTERM,' terms'
       ELSE
         WRITE(6,*) ' series expansion truncated after ',MXTERM,' terms'
       END IF
       WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' LUC,LUHC,LUSC1,LUSC2',LUC,LUHC,LUSC1,LUSC2
       WRITE(6,*) ' Initial vector on LUC '
       CALL UNIT_INFO(LUC)
       CALL UNIT_INFO(LUHC)
       CALL UNIT_INFO(LUSC1)
       CALL UNIT_INFO(LUSC2)
       IF (NTEST.GE.1000) THEN
         CALL WRTVCD(VEC1,LUC,1,LBLK)
       ELSE
         CALL WRTVSD(VEC1,LUC,1,LBLK)
       END IF
      END IF
* Tell integral fetcher to fetch cc amplitudes, not integrals
      INT_CC = 1
*. Loop over orders of expansion
      N = 0
      XFACN = 1.0D0
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) 'TAMP:'
        CALL WRT_CC_VEC2(TAMP,6,CCTYPE)
      END IF

      IF (IOPTYP.NE.0) THEN
        CALL CONJ_CCAMP(TAMP,1,TSCR)
        CALL SCALVE(TSCR,SFAC,N_CC_AMP)
        IF(NTEST.GE.500) THEN
          WRITE(6,*) 'TAMP+:'
          CALL WRT_CC_VEC2(TSCR,6,CCTYPE)
        END IF        
      END IF
*
      CALL COPVCD(LUC,LUSC1,VEC1,1,LBLK)
      CALL COPVCD(LUC,LUHC,VEC1,1,LBLK)
*
 1000 CONTINUE
       N = N+1
       IF(NTEST.GE.5) THEN
         WRITE(6,*) ' Info for N = ', N
       END IF
*. T^N  times vector on LUSC1
C?     WRITE(6,*) ' Input vector to MV7 '
C?     CALL WRTVCD(VEC1,LUSC1,1,LBLK)
       IF(CCTYPE(1:2).EQ.'CC') THEN
         STOP 'DO WE STILL NEED THAT ROUTE?'
         CALL MV7(VEC1,VEC2,LUSC1,LUSC2,0,0)
       ELSE IF(CCTYPE(1:6).EQ.'GEN_CC') THEN
         CALL SIG_GCC(VEC1,VEC2,LUSC1,LUSC2,TAMP)
         IF(NTEST.GE.500.AND.IOPTYP.NE.0) THEN
           WRITE(6,*) ' 1/(N-1)! T (T +/- T^+)**(N-1) |0> '
           WRITE(6,*) ' =================================='
           CALL WRTVCD(VEC1,LUSC2,1,LBLK)
         END IF
C  SIG_GCC(C,HC,LUC,LUHC,T)
* General CC sigma
       END IF
* NEW: scale result directly with 1/N which should result in
*      higher numerical stability for large exponents
       FAC = 1D0/DBLE(N)
       IF(IOPTYP.NE.0) THEN
         CALL SCLVCD(LUSC2,LUSC3,FAC,VEC1,1,LBLK)
         CALL CONJ_T
         CALL SIG_GCC(VEC1,VEC2,LUSC1,LUSC2,TSCR)
         CALL CONJ_T
         IF(NTEST.GE.500) THEN
           WRITE(6,*) ' 1/(N-1)! T^+ (T +/- T^+)**(N-1) |0> '
           WRITE(6,*) ' =================================='
           IF (NTEST.GE.5000) THEN
             CALL WRTVCD(VEC1,LUSC2,1,LBLK)
           ELSE
             CALL WRTVSD(VEC1,LUSC2,1,LBLK)
           END IF
         END IF
c                                      in1   in2   res         
         CALL VECSMD(VEC1,VEC2,FAC,1d0,LUSC2,LUSC3,LUSC1,1,LBLK)
       ELSE
         CALL SCLVCD(LUSC2,LUSC1,FAC,VEC1,1,LBLK)
       END IF
       IF(NTEST.GE.500) THEN
         WRITE(6,*) ' 1/N! T**(N) |0> '
         WRITE(6,*) ' ================'
         IF (NTEST.GE.5000) THEN
           CALL WRTVCD(VEC1,LUSC1,1,LBLK)
         ELSE
           CALL WRTVSD(VEC1,LUSC1,1,LBLK)
         END IF
       END IF
*. Norm of this correction term
       LBLK = -1
       XNORM2 = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
       XMXNRM = FDMNXD(LUSC1,2,VEC1,1,LBLK)
       XFACN=1D0
       XNORM = SQRT(XNORM2)
       IF(NTEST.GE.5) THEN
         WRITE(6,*) ' Norm of correction ', XNORM, XMXNRM
       END IF
*. Update output file with 1/N! T^N !ref>
       ONE = 1.0D0
       CALL VECSMD(VEC1,VEC2,ONE,ONE,LUSC1,LUHC,LUSC2,1,LBLK)
       CALL COPVCD(LUSC2,LUHC,VEC1,1,LBLK)
*. give up?
      IF (XNORM.GT.1d+100) THEN
        WRITE(6,*) 'Wavefunction blows up! Take a step back :-)'
        WRITE(6,*) ' Norm of last 1/N! T^N !ref>: ',XNORM,' for N=',N
        XNORM=SQRT(INPROD(TAMP,TAMP,N_CC_AMP))
        WRITE(6,*) ' Norm of T was: ', XNORM
        STOP 'WOOMM!'
      END IF
*. Take another turn ?
      IF(XNORM.GT. THRES_C.AND. N .LT. MXTERM_L) GOTO 1000
*
*. Result on LUHC
      IF (XNORM.GT.THRES_C.AND.STRICT) THEN
        WRITE(6,'(x,a,i5,a)')
     $        'Fatal: No convergence in EXPT_REF (max. iter.:',
     $        MXTERM_L, ' )'
        STOP 'No convergence in EXPT_REF!'
      END IF
C      CALL COPVCD(LUSC3,LUHC,VEC1,1,LBLK)
      IF(NTEST.GE.2.AND.STRICT) THEN
        WRITE(6,*) ' Convergence obtained in ', N, ' iterations'
        WRITE(6,*) ' Norm of last correction ', XNORM
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============'
        WRITE(6,*) ' Exp T |ref> '
        WRITE(6,*) ' ============'
        WRITE(6,*)
         IF (NTEST.GE.1000) THEN
           CALL WRTVCD(VEC1,LUHC,1,LBLK)
         ELSE
           CALL WRTVSD(VEC1,LUHC,1,LBLK)
         END IF
      END IF
*
      RETURN
      END 
*------------------------------------------------------------------------*
      SUBROUTINE VEC_FROM_DISC(VEC,LENGTH,IREW,LBLK,LU)
*
* Read vector VEC of length LENGTH from discfile LU in
* standard LUCIA form, as for example written by VEC_TO_DISC
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VEC(LENGTH)
*
      IF(IREW.EQ.1) THEN  
        CALL REWINO(LU)
      END IF
      CALL IFRMDS(LENGTH2,1,LBLK,LU)
      IF (LENGTH.NE.LENGTH2) THEN
        WRITE(6,*)' Length of vector on file does not match with input!'
        WRITE(6,*) ' File number = ', LU
        WRITE(6,*) ' I take the shorter one of ',LENGTH, LENGTH2
        LENGTH2 = MIN(LENGTH2,LENGTH)
      END IF
      CALL FRMDSC(VEC,LENGTH2,LBLK,LU,IAMZERO,IAMPACKED)
*. Skip the -1 at end 
      CALL IFRMDS(IONEM,1,LBLK,LU)
*
      RETURN
      END
      SUBROUTINE VEC_FROM_DISC_E(VEC,LENGTH,IREW,LBLK,LU,IERR)
*
* Read vector VEC of length LENGTH from discfile LU in
* standard LUCIA form, as for example written by VEC_TO_DISC
*
* version with error code (andreas, 2004)
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VEC(LENGTH)
*
      IF(IREW.EQ.1) THEN  
        CALL REWINO(LU)
      END IF
      CALL IFRMDSE(LENGTH2,1,LBLK,LU,IERR)
      IF (IERR.NE.0) RETURN
      IF (LENGTH.NE.LENGTH2) THEN
        WRITE(6,*) ' Length of vector on file do not match with input!'
        WRITE(6,*) ' I take the shorter one of ',LENGTH, LENGTH2
        LENGTH2 = MIN(LENGTH2,LENGTH)
      END IF
      CALL FRMDSCE(VEC,LENGTH2,LBLK,LU,IAMZERO,IAMPACKED)
      IF (IERR.NE.0) RETURN
*. Skip the -1 at end 
      CALL IFRMDSE(IONEM,1,LBLK,LU,IERR)
*
      RETURN
      END
      SUBROUTINE VEC_TO_DISC(VEC,LENGTH,IREW,LBLK,LU)
*
* Write vector VEC of length LENGTH to discfile LU in
* standard LUCIA form
*
* 1 : Length as written by ITODSC
* 2 : The vector as written by TODSC
* 3 : End of vector, -1 written by ITODSC
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VEC(LENGTH)
*
      IF(IREW.EQ.1) THEN  
        CALL REWINO(LU)
      END IF
      CALL ITODS(LENGTH,1,LBLK,LU)
      CALL TODSC(VEC,LENGTH,LBLK,LU)
      IONEM = -1 
      CALL ITODS(IONEM,1,LBLK,LU)
*
      RETURN
      END
      subroutine optim_cc_new(irestrt_cc,ecc,lgrad,cc_amp,ccvec1,ccvec2,
     &     vec1,vec2,
     &     ittact,iooexcc,iooexc,nooexcl,cctype,
     &     luccamp,lu_omg,lu_lamp,lu_1den,lu_2den,
     &     lucmo,lu1into,lu2into,E_FINAL,ERROR_NORM_FINAL,CONV_F)
*
* Optimize CC wave function
*
* New version, 2004
*
* Input :
* ========
*
* LUCCAMP : Initial set of CC amplitiudes
*
* Output
* ========
*
* LUCCAMP : Final set of amplitudes 
*
* prelim: LUDIA will contain the diagonal
*
* Scratch
* ========
*
* VEC1, VEC2 : Two vectors for CI calculations (old CCVEC)
* CCVEC1, CCVEC2 : Two complete CC vectors
*
c      implicit real*8(a-h,o-z)
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
      real*8 inprod, inprdd
      character*6 cctype
      LOGICAL CONV_F
      include 'crun.inc'
      include 'glbbas.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cecore.inc' 
      include 'clunit.inc'
      include 'ctcc.inc'
      include 'cprnt.inc'
      include 'cc_exc.inc'
      include 'cands.inc'
      include 'cstate.inc'
      include 'cintfo.inc'
      include 'cgas.inc'
      include 'oper.inc'
      include 'opti.inc'
      include 'symrdc.inc'
cVOG  common/cntdl_max/ntdl_max_act
*. scratch 
      dimension vec1(*),vec2(*), ccvec1(*), ccvec2(*)
*. input and output
      logical lgrad
      dimension cc_amp(*)
      dimension ittact(ngas*ngas,2)
      dimension iooexcc(*)
      dimension nooexcl(2)     
*. 5 little local arrays for informations
      logical lrestart
      dimension eccinf(maxit), ecainf(maxit),
     &          vfninf(maxit), ampinf(maxit), dtninf(maxit)

*
*
      call atim(cpu0,wall0)
*
      idummy = 0
      call memman(idummy,idummy,'MARK  ',idummy,'CC_OPT')
*
      lblk = -1
      ntdl_max_act = 0
*
      ntest = 50 
      ntest = max(ntest,iprcc)
      if(ntest.ge.5) then
        write(6,'(a,a)') ' optim_cc: form of CC-expansion ',cctype
        write(6,*) 'luccamp,lu_lamp: ',luccamp,lu_lamp
      end if
*. Optimization of coupled cluster wavefunction defines symmetry as
      icsm = irefsm
      issm = irefsm
      itex_sm = 1
*      
*. Coupled cluster flag --- should be redundant now, ikke?
      icc_exc = 1
*
* ===========================================
*             Construct diagonal 
* ===========================================
*
* Obtain Inactive Fock matrix 
*
      call copvec(work(kint1o),work(kfi),nint1)
      call iswpve(iphgas,iphgas1,ngas)
*
      i_use_simtrh = 0
      i_unrorb_merk = i_unrorb
      i_unrorb = 0
      call fi(work(kfi),ecc,1)
      i_unrorb = i_unrorb_merk
      call copvec(work(kfi),work(kfio),nint1)
      call iswpve(iphgas,iphgas1,ngas)
*
      ivcc = 0
      if (ccform(1:3).eq.'VCC') ivcc = 1
      iucc = 0
      if (ccform(1:3).eq.'UCC') iucc = 1
      iecc = 0
      if (ccform(1:3).eq.'ECC') iecc = 1

      if(i_do_newccv.ge.1) then
        if(ccform(1:3).eq.'VCC'.or.ccform(1:3).eq.'UCC') then
          write(6,*)'ERROR: VCC cannot go along with i_do_newccv.eq.1'
          stop ' inconsistency in optim_cc2'
        end if
*. new cc vector function uses particle hole form so 
*. use inactive fock matrix and modified core-energy 
        ecore_h = ecore_ab
        ecore = ecore_orig + ecore_h
        write(6,*) ' Updated core energy ',ecore 
        call copvec(work(kfi),work(kint1),nint1)
      else
        ecore_h = 0.0d0
        ecore_orig = ecore_ini
      end if
*. and the diagonal 

      ! preliminary position:
      nspin = 1
      if (ireftyp.eq.2) nspin = 2
      
      norbhss = 0
      if (i_bcc.eq.1.or.i_oocc.eq.1.or.i_obcc.eq.1) then
        norbhss = nspin
        nooexc_tot = nooexcl(1)+nooexcl(2)
        call memman(idum,idum,'MARK  ',idum,'OHSS  ')
        call memman(korbhss,nooexc_tot,'ADDL  ',2,'ORBHSS')
        lu_odia = iopen_nus('OCC_ODIAG')
      end if

      imod = 2 ! use alpha/beta fock-matrix
      call gencc_f_diag_m(imod,work(klsobex),nspobex_tp,ccvec1,1,
     &      work(korbhss),nooexcl,iooexcc,norbhss,
     &      work(kvec1p),work(kvec2p),mx_st_tsoso_mx,
     &      mx_st_tsoso_blk_mx)

      if (i_bcc.eq.1.or.i_oocc.eq.1.or.i_obcc.eq.1) then
        call vec_to_disc(work(korbhss),
     &       nooexcl(1)+nooexcl(2),1,lblk,lu_odia)
        call memman(idum,idum,'FLUSM ',idum,'OHSS  ')
      end if
*
      if(ntest.ge.100) then
        write(6,*) ' optim_cc : Diagonal '
        call wrtmat(ccvec1,1,n_cc_amp,1,n_cc_amp)
      end if
      if(i_do_mask_cc.eq.1) then
*. mask the excitations not arising from mask det
        value = 1.0d+15
        call mask_ccvec(work(klsobex),nspobex_tp,ccvec1,1,
     &       mask_sd,msk_ael,msk_bel,value,mx_st_tsoso_blk_mx)
      end if
*. save diagonal on disc
c prelim: shift the diagonal, if neg. eigenvalues appear:
      xmin = 1d13
      do ii = 1, n_cc_amp
        xmin = min(ccvec1(ii),xmin)
      end do
      if (xmin.lt.0d0) then
        print *,' min. element of diagonal: ',xmin
        print *,' shift diagonal by ',-xmin+0.01d0
        do ii = 1, n_cc_amp
          ccvec1(ii) = ccvec1(ii)-xmin+0.01d0
        end do
      end if

c end prelim
      call vec_to_disc(ccvec1,n_cc_amp,1,lblk,ludia)
      
      if (iuse_ph.eq.1.and.ireftyp.eq.2.and.i_do_newccv.eq.0) then
        do ii = 1, 10
          print *,' IPHGAS GEFUMMEL!!!'
        end do
        do ii = 1, ngas
          if (iphgas(ii).eq.1) then
            ihpvgas_ab(ii,1) = 2
            ihpvgas_ab(ii,2) = 2
          end if
          if (iphgas(ii).eq.2) then
            ihpvgas_ab(ii,1) = 1
            ihpvgas_ab(ii,2) = 1
          end if
        end do
      END IF
*
* set some variables on common /opti/
      ivar = 0
      ilin = 0
*. (Maxit is obtained from CRUN)
      maxmacit = maxit
      micifac  = 20
      maxmicit = maxmacit*micifac

      DTNINF(1)=0.0d0 ! step to first vector was zero

      xngrad = 1000
      xnomg  = 1000
      itask = 0
      imacit = 0
      imicit = 0
      imicit_tot = 0
      itaskl = 0
      imacitl = 0
      imicitl = 0
      imicit_totl = 0
      energy = 0d0
      iprint = 100

      lulsig = iopen_nus('CCOPTLSIG')
      lursig = iopen_nus('CCOPTRSIG')
      if (lu_omg.le.0) then
        luomg = iopen_nus('CCOPTRES')
      else
        luomg = lu_omg
      end if
      lutrvec = iopen_nus('CCOPTVEC')
      luint1 = iopen_nus('CC_INTM1')
      luint2 = iopen_nus('CC_INTM2')
      luint3 = iopen_nus('CC_INTM3')

      i_lr_sim = 0

      i_packvec = 0
      if (iecc.eq.1) i_lr_sim = 1
      if (iecc.eq.1) i_packvec = 1

      if ((i_obcc.eq.1.or.i_oocc.eq.1).and.ivcc.eq.0) i_lr_sim = 1
c TEST
      if ((i_obcc.eq.1.or.i_oocc.eq.1).and.ivcc.eq.0) i_packvec = 2

      if (i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1) then
        lukappa = iopen_nus('OCC_KAPPA')
        lutrmat = iopen_nus('OCC_TRMAT')
        lu_ogrd = iopen_nus('OCC_OGRAD')

        lutrv_l = iopen_nus('CCTRVECL')
        lutrv_o = iopen_nus('CCTRVECO')

        lusig_l = iopen_nus('CCSIG_L')
        lusig_o = iopen_nus('CCSIG_O')

      end if

      if (i_lr_sim.eq.0) then
        luomg_c = luomg
        lutrvec_c = lutrvec
        luamp_c = luccamp
        ludia_c = ludia
      else
        ! combined amplitudes for optimizer
        lurhs = iopen_nus('CCOPTRHS')
        luomg_c = iopen_nus('CCOPTRESC')
        lutrvec_c = iopen_nus('CCOPTVECC')
        luamp_c = iopen_nus('CCOPTAMPC')
        ludia_c = iopen_nus('CCDIAC')
      end if

      lrestart = .not.symred.and.(irestrt_cc.ne.0.or.lgrad)
      call init_vec(luccamp,n_cc_amp,ccvec1,i_bcc,lrestart)
      
      if (i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1) then

        call init_vec(lukappa,nooexc_tot,ccvec1,0,lrestart)

        idum = 0
        call memman(idum,idum,'MARK  ',idum,'LUMATB')
        call memman(klumat,ntoob*ntoob,'ADDL  ',2,'UMAT')
        work(klumat:klumat-1+ntoob*ntoob) = 0d0
        ioff = 0
        do ism = 1, nsmob
          ndim = ntoobs(ism)
          idx = 1
          do ii = 1, ndim
            work(klumat-1+ioff+idx) = 1d0
            idx = idx + ndim+1
          end do
          ioff = ioff + ndim*ndim
        end do
        ndimtot = ioff
c        call memchk2('abuse ')
        call vec_to_disc(work(klumat),ndimtot,1,lblk,lutrmat)
        if (nspin.eq.2)
     &       call vec_to_disc(work(klumat),ndimtot,0,lblk,lutrmat)

        call memman(idum,idum,'FLUSM ',idum,'LUMATB')

        ccvec1(1:n_cc_amp) = 0d0
        xnrm = sqrt(inprdd(ccvec1,ccvec1,lukappa,lukappa,1,lblk))
        if (nspin.eq.2) then
          call copvec(work(kint1),work(kint1b),nint1)
        end if
        ! nspin.eq.2: use routine to init the other 2int integrals
        if (xnrm.gt.1d-10.or.nspin.eq.2) then
          call kap2u(3,ludum,lukappa,lutrmat,iooexcc,nooexcl,nspin)
          call tra_kappa(lukappa,lutrmat,iooexcc,nooexcl,nspin,
     &                   1,lu1into,lu2into)
        end if
        
      end if

      if (i_lr_sim.eq.1) then
        call init_vec(lu_lamp,n_cc_amp,ccvec1,0,lrestart)
        if (i_packvec.eq.1.and.i_oocc.eq.0.and.i_oocc.eq.0) then
          ! pack it
          imode = 11
          ludum = -1
          call cmbamp(imode,luccamp,lu_lamp,ludum,luamp_c,
     &         ccvec1,n_cc_amp,n_cc_amp,0)
          ! double the diagonal
          imode = 11
          call cmbamp(imode,ludia,ludia,ludum,ludia_c,
     &         ccvec1,n_cc_amp,n_cc_amp,0)
        else if (i_packvec.eq.1.and.(i_oocc.eq.1.or.i_obcc.eq.1)) then
          ! pack it
          imode = 11
          ludum = -1
          call cmbamp(imode,luccamp,lu_lamp,lukappa,luamp_c,
     &         ccvec1,n_cc_amp,n_cc_amp,nooexcl)
          ! double the diagonal
          imode = 11
          call cmbamp(imode,ludia,ludia,lu_odia,ludia_c,
     &         ccvec1,n_cc_amp,n_cc_amp,nooexcl)
        end if
      end if

      call atim(cpu00i,wall00i)

      itransf = 0

      if (i_lr_sim.eq.0) then
        if (i_bcc.eq.1.or.(i_oocc.eq.1.and.ivcc.eq.1)) then
          write(6,'(">>>",A)')
     &     '  Iter.               energy     norm(Omg)  grad(kappa)'
          write(6,'(">>>",A)')
     &     '-------------------------------------------------------'
        else
          write(6,'(">>>",A)')
     &     '  Iter.               energy     norm(Omg)     norm(T) '
          write(6,'(">>>",A)')
     &     '-------------------------------------------------------'
        end if
      else if (i_obcc.eq.1.or.i_oocc.eq.1) then
          write(6,'(">>>",A)')
     &     '  Iter.               energy     norm(Omg)     grad(T) '//
     &                                                    ' grad(kappa)'
          write(6,'(">>>",A)')
     &     '-------------------------------------------------------'//
     &                                                    '------------'
      else
          write(6,'(">>>",A)')
     &     '  Iter.        Lambda energy          proj. CC energy  '//
     &             '   norm(Omg)     grad(T)    norm(T)     norm(L)'
          write(6,'(">>>",A)')
     &     '-------------------------------------------------------'//
     &             '-----------------------------------------------'
      end if

*----------------------------------------------------------------------*
*     Optimization loop:
*----------------------------------------------------------------------*
      DO WHILE (ITASK.LT.8)
        call atim(cpu0i,wall0i)

        if(ntest.ge.2) then
           write(6,*)
           write(6,*) ' ==========================================='
           write(6,*) '  Information from iteration ', imacit
           write(6,*) ' ==========================================='
           write(6,*)
        end if
          
*----------------------------------------------------------------------*
*     optimization kernel
*----------------------------------------------------------------------*
        if ((i_lr_sim.eq.0.or.i_packvec.eq.1)
     &       .and.i_bcc.ne.1.and..not.(i_oocc.eq.1.and.ivcc.eq.1)) then
          if (i_lr_sim.eq.0) then
            luamp_c = luccamp
            luomg_c = luomg
            ludia_c = ludia
          end if
          call optcont(imacit,imicit,imicit_tot,iprint,
     &                   itask,iconv,
     &                   luamp_c,lutrvec,
     &                   ener,
     &                   ccvec1,ccvec2,n_cc_amp,
     &                   luomg_c,lursig,ludia_c,
     &                   0,ludum)
          if (i_oocc.eq.1.and.imacit.gt.1) itransf = 1
        else
          ikapmod = 1
          if (iorder.eq.2) ikapmod = 1 !0
          if (ikapmod.eq.0) then
            ! do not keep track of total kappa, just get new
            ! transformation matrix as exp(dkappa)exp(kappa)
            ! does not allow for DIIS
            if (i_obcc.eq.1.or.i_oocc.eq.1) iorbopt = -3
            if ((i_obcc.eq.1.or.i_oocc.eq.1)
     &                 .and.i_packvec.eq.2) iorbopt = -2
          else if (ikapmod.eq.1) then
            ! CEP gradient: get new kappa as ln(exp(dkappa)exp(kappa0)) 
            ! should be stable now with new ln(U)-routine
            if (i_obcc.eq.1.or.i_oocc.eq.1) iorbopt = 3
            if ((i_obcc.eq.1.or.i_oocc.eq.1)
     &                 .and.i_packvec.eq.2) iorbopt = 2
          else if (ikapmod.eq.2) then
            ! real gradient at kappa0: get new kappa as kappa+dkappa
            if (i_obcc.eq.1.or.i_oocc.eq.1) iorbopt = -3
            if ((i_obcc.eq.1.or.i_oocc.eq.1)
     &                 .and.i_packvec.eq.2) iorbopt = -2
          end if

          if (i_bcc.eq.1) iorbopt = 1
          if (i_oocc.eq.1.and.ivcc.eq.1) iorbopt = 1

          call optcont_orbopt
     &                  (imacit,imicit,imicit_tot,iprint,
     &                   itask,iconv,iorbopt,ikapmod,itransf,
     &                   luccamp,lu_lamp,lukappa,lutrmat,
     &                   lutrvec,lutrv_l,lutrv_o,
     &                   ener,
     &                   ccvec1,ccvec2,n_cc_amp,nooexcl,iooexcc,nspin,
     &                   luomg,lulsig,lu_ogrd,lursig,lusig_l,lusig_o,
     &                   ludia,lu_odia,0,ludum)
        end if

        if (i_packvec.eq.1.and.i_oocc.eq.0.and.i_obcc.eq.0) then
          ! unpack the amplitudes for the workers
          imode = 01
          ludum = -1
          call cmbamp(imode,luccamp,lu_lamp,ludum,luamp_c,
     &         ccvec1,n_cc_amp,n_cc_amp,0)
        else if (i_packvec.eq.1.and.(i_oocc.eq.1.or.i_obcc.eq.1)) then
          imode = 01
          call cmbamp(imode,luccamp,lu_lamp,lukappa,luamp_c,
     &         ccvec1,n_cc_amp,n_cc_amp,nooexcl)
        end if

*----------------------------------------------------------------------*
*     a: orbital transformation, if required
*----------------------------------------------------------------------*
        if ((i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1)
     &       .and.itransf.eq.1) then
          call tra_kappa(lukappa,lutrmat,iooexcc,nooexcl,nspin,
     &                   1,lu1into,lu2into)
        end if

*----------------------------------------------------------------------*
*     b: get vector-function (CC-residual)
*----------------------------------------------------------------------*
        if (iand(itask,1).eq.1.or.iand(itask,2).eq.2) then
          call cc_vec_fnc2(ccvec1,ccvec2,ecc1,eccl,
     &                     xampnrm,xomgnrm,xlampnrm,
     &                     vec1,vec2,ibio,cctype,
     &                     cc_amp,
     &                     luccamp,luomg,lu_lamp,
     &                     luint1,luint2,luint3)

*----------------------------------------------------------------------*
*     c: get left-hand Jacobi-contraction (gives gradient wrt T)
*----------------------------------------------------------------------*
          if (i_lr_sim.eq.1) then
            if (iecc.eq.0)
     &        call zero_ord_rhs(ccvec1,vec1,vec2,luccamp,lurhs,luint1)

            lr_switch = 1
            iadd_rhs = 1
            if (iecc.eq.1) iadd_rhs = 0
            eccl = ecc1
            itex_sm = 1
            call jac_t_vec2(lr_switch,iadd_rhs,0,1,1,
     &           ccvec1,ccvec2,vec1,vec2,
     &           n_cc_amp,n_cc_amp,
     &           eccl,xlampnrm,xlresnrm,
     &           luccamp,luomg,lu_lamp,lulsig,lurhs,
     &           luint1,luint2,luint3)
          end if

*----------------------------------------------------------------------*
*     d: orbital gradients
*----------------------------------------------------------------------*
          if (i_obcc.eq.1.or.i_bcc.eq.1) then
            ! Brueckner says: take Omega_1 as orbital gradient
            ! NOTE: the T1-part on luomg is set to zero
            !       therefore call is after jac_t_vec
            call omg2ogrd(luomg,lu_ogrd,ccvec1,
     &             iooexcc,nooexcl,n_cc_amp,nspin,
     &             xkresnrm,xomgnrm)
          end if

          if (i_obcc.eq.1.or.i_oocc.eq.1) then
            ! get the orbital gradient
            ! for the hybrid method, we have to modify the
            ! singles part on lulsig
            call occ_orbgrad(i_obcc,ivcc,
     &             ccvec1,ccvec2,vec1,vec2,
     &             ittact,iooexcc,iooexc,
     &             nooexcl,n_cc_amp,nspin,
     &             xkapnrm,xkresnrm,xlresnrm,
     &             luccamp,lu_lamp,lulsig,
     &             lukappa,lu_ogrd,lu_odia,
     &             lu_1den,lu_2den,
     &             lu1into,lu2into,
     &             luint1,luint3,luc)

            imod = 2            ! use alpha/beta fock-matrix
            call gencc_f_diag_m(imod,work(klsobex),
     &           nspobex_tp,ccvec1,1,
     &           xdummy,nooexcl,iooexcc,0,
     &           vec1,vec2,mx_st_tsoso_mx,
     &           mx_st_tsoso_blk_mx)

            xmin = 1d13
            do ii = 1, n_cc_amp
              xmin = min(ccvec1(ii),xmin)
            end do
            if (xmin.lt.0d0) then
              print *,' min. element of diagonal: ',xmin
              print *,' shift diagonal by ',-xmin+0.01d0
              do ii = 1, n_cc_amp
                ccvec1(ii) = ccvec1(ii)-xmin+0.01d0
              end do
            end if
            
            call vec_to_disc(ccvec1,n_cc_amp,1,lblk,ludia)

          end if

          if (i_lr_sim.eq.1) then
            ! pack everything for optimizer
            if (i_packvec.eq.1.and.i_oocc.eq.0.and.i_obcc.eq.0) then
              imode = 11
              ludum = -1
              if (iecc.eq.0) then
                call cmbamp(imode,luomg,lulsig,ludum,luomg_c,
     &               ccvec1,n_cc_amp,n_cc_amp,0)
              else
                call cmbamp(imode,luomg,lulsig,ludum,luomg_c,
     &               ccvec1,n_cc_amp,n_cc_amp,0)
              end if
            else if (i_packvec.eq.1.and.(i_oocc.eq.1.or.i_obcc.eq.1))
     &             then
              imode = 11
              call cmbamp(imode,luomg,lulsig,lu_ogrd,luomg_c,
     &             ccvec1,n_cc_amp,n_cc_amp,nooexcl)
            end if
          end if


*----------------------------------------------------------------------*
*     some output
*----------------------------------------------------------------------*
          if(i_do_newccv.ge.1) then
            ecc = ecc1
          else
            ecc = ecc1 + ecore
            eccl = eccl + ecore
          end if

          eccinf(imacit)=ecc
          ecainf(imacit)=eccl
          vfninf(imacit)=xomgnrm
          ampinf(imacit)=xampnrm
*
          if (i_lr_sim.eq.0) then
            if (i_bcc.eq.1.or.(i_oocc.eq.1.and.ivcc.eq.1)) then
              ener = ecc
              write(6,'(">>>",i5,f25.12,2(2x,e10.4))') 
     &             imacit,ecc,xomgnrm,xkresnrm
            else
              ener = ecc
              write(6,'(">>>",i5,f25.12,2(2x,e10.4))') 
     &             imacit,ecc,xomgnrm,xampnrm
            end if
          else if (i_oocc.eq.1.or.i_obcc.eq.1) then
            ener = eccl
            write(6,'(">>>",i5,f25.12,3(2x,e10.4))') 
     &           imacit,eccl,xomgnrm,xlresnrm,xkresnrm
          else
            ener = eccl
            write(6,'(">>>",i5,2f25.12,4(2x,e10.4))') 
     &           imacit,eccl,ecc,xomgnrm,xlresnrm,xampnrm,xlampnrm
          end if
*
        end if

*----------------------------------------------------------------------*
*     e: Hessian-vector contractions for 2nd-order methods follow
*----------------------------------------------------------------------*
        if (iand(itask,4).eq.4) then
          if (i_bcc.eq.1) stop 'no Newton meets Brueckner'
          if (i_oocc.eq.1.or.i_obcc.eq.1) then
            call occ_scnd_num(i_obcc,
     &           ccvec1,ccvec2,vec1,vec2,
     &           ittact,iooexcc,iooexc,
     &           nooexcl,n_cc_amp,nspin,
     &           lu1into,lu2into,luc,
     &           luccamp,lu_lamp,lukappa,
     &           lutrvec,lutrv_l,lutrv_o,
     &           lursig,lusig_l,lusig_o,
     &           1,1,1)
          else
            lr_switch = 2
            iadd_rhs = 0
            itex_sm = 1

            call jac_t_vec2(lr_switch,iadd_rhs,0,1,1,
     &         ccvec1,ccvec2,vec1,vec2,
     &         n_cc_amp,n_cc_amp,
     &         eccl,xlampnrm,xlresnrm,
     &         luccamp,luomg,lutrvec,lursig,lurhs,
     &         luint1,luint2,luint3)

            numtest=0
            if (numtest.eq.1) then
              call jac_t_vec_num(ccvec1,ccvec2,vec1,vec2,
     &             luccamp,lutrvec,lursig,lursig)
            end if
          end if

        end if

        call atim(cpui,walli)
        
        call prtim(6,'time for current iteration',
     &       cpui-cpu0i,walli-wall0i)

      END DO
*----------------------------------------------------------------------*
*     ^ End of loop over iterations
*----------------------------------------------------------------------*

      call atim(cpu,wall)
      call prtim(6,'time for CC optimization',
     &     cpu-cpu0,wall-wall0)
      call prtim(6,'average time for iterations',
     &     (cpu-cpu00i)/dble(imacit),(wall-wall00i)/dble(imacit))
*
* for some CC-variants we may now easily calculate the density
      if (idensi.ne.0.and.(ccform(1:3).eq.'ECC'.or.
     &                     ccform(1:3).eq.'VCC'.or.
     &                     ccform(1:3).eq.'UCC')) then

        if (ccform(1:3).eq.'ECC') then
          luleft  = luint3
          luright = luint1
        else 
          luleft  = luint3
          luright = luint1
          call copvcd(luright,luleft,vec1,1,lblk)
        end if

        write(6,*) ' calculating densities now ... '

        call atim(cpu0,wall0)
        call densi2(idensi,work(krho1),work(krho2),
     &         vec1,vec2,luleft,luright,exps2,
     &         0,work(ksrho1))

        lrho1 = ntoob**2
        lrho2 = ntoob**2*(ntoob**2+1)/2
        
        if (ccform(1:3).eq.'VCC') then
          ovl = inprdd(vec1,vec1,luright,luright,1,lblk)
          call scalve(work(krho1),1d0/ovl,lrho1)
          if (idensi.eq.2)
     &         call scalve(work(krho2),1d0/ovl,lrho2)
        end if

        if (ccform(1:3).eq.'ECC') then
          call sym_blmat(work(krho1),1,ntoob)
          if (idensi.eq.2)
     &         call sym_2dens(work(krho2),ntoob,0)
        end if

        call vec_to_disc(work(krho1),lrho1,1,lblk,lu_1den)

        if (idensi.eq.2) then
          call vec_to_disc(work(krho2),lrho2,1,lblk,lu_2den)
        end if

        call atim(cpu,wall)

        call prtim(6,'time for densities',
     &       cpu-cpu0,wall-wall0)

      end if

      if ((i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1).and.lgrad) then
        ! make and save new CMO-coefficients
        call mknewcmo(lukappa,lutrmat,-lucmo,iooexcc,nooexcl,nspin)
        call relunit(lukappa,'keep')
      else if (i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1) then
        call relunit(lukappa,'keep')
      end if

      call relunit(lursig,'delete')
      call relunit(lulsig,'delete')
      if (lu_omg.le.0) 
     &    call relunit(luomg,'delete')
      call relunit(lutrvec,'delete')
      call relunit(luint1,'delete')
      call relunit(luint2,'delete')
      call relunit(luint3,'delete')

      if (i_lr_sim.ne.0) then
        call relunit(lurhs,'delete')
        call relunit(luomg_c,'delete')
        call relunit(lutrvec_c,'delete')
        call relunit(luamp_c,'delete')
        call relunit(ludia_c,'delete')
      end if

      if (i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1) then
        call relunit(lutrmat,'delete')
        call relunit(lu_ogrd,'delete')
        call relunit(lu_odia,'delete')
      end if
*
      if(iconv.eq.0) then
        write(6,'(a,i4,a)') 
     &  ' Convergence not obtained in ', IMACIT, ' iterations'
      else
        write(6,'(a,i4,a)') 
     &  ' Convergence obtained in ', IMACIT, ' iterations'
        if (ccform(1:3).eq.'TCC') then
          if (i_oocc.eq.1) then
        write(6,'(/,a,/,a,/,a,f25.12,a,/,a,/,a,/)')
     &   ' *****************************************************',
     &   ' *                                                   *',
     &   ' *   Final OCC energy:',ECCL,                 ' E_h  *',
     &   ' *                                                   *',
     &   ' *****************************************************'
          else if (i_obcc.eq.1.or.i_bcc.eq.1) then
        write(6,'(/,a,/,a,/,a,f25.12,a,/,a,/,a,/)')
     & ' *************************************************************',
     & ' *                                                           *',
     & ' *   Final Brueckner CC energy:',ECC,                  ' E_h *',
     & ' *                                                           *',
     & ' *************************************************************'
          else
        write(6,'(/,a,/,a,/,a,f25.12,a,/,a,/,a,/)')
     &   ' *****************************************************',
     &   ' *                                                   *',
     &   ' *   Final CC energy:',ECC,                  ' E_h   *',
     &   ' *                                                   *',
     &   ' *****************************************************'
          end if
        else if (ccform(1:3).eq.'VCC') then
          if (i_oocc.eq.1) then
        write(6,'(/,a,/,a,/,a,f25.12,a,/,a,/,a,/)')
     &' **************************************************************',
     &' *                                                            *',
     &' *  Final variational OCC energy:',ECC,         '          E_h*',
     &' *                                                            *',
     &' **************************************************************'
          else
        write(6,'(/,a,/,a,/,a,f25.12,a,/,a,/,a,/)')
     &' **************************************************************',
     &' *                                                            *',
     &' *  Final variational CC energy:',ECC,                  ' E_h *',
     &' *                                                            *',
     &' **************************************************************'
          end if
        else if (ccform(1:3).eq.'UCC') then
          if (ccform(1:4).eq.'UCC2') then
            write(6,'(/,a,/,a,/,a,f25.12,a,/,a,/,a,/)')
     &' **************************************************************',
     &' *                                                            *',
     &' *  Final UCC2 (CEPA-0) energy:',ECC,                 ' E_h   *',
     &' *                                                            *',
     &' **************************************************************'
          else
            write(6,'(/,a,/,a,/,a,f25.12,a,/,a,/,a,/)')
     &' **************************************************************',
     &' *                                                            *',
     &' *  Final variational UCC energy:',ECC,                 ' E_h *',
     &' *                                                            *',
     &' **************************************************************'
          end if
        else if (ccform(1:3).eq.'ECC') then
        write(6,'(/,a,/,a,/,3a,f25.12,a,/,a,/,a,/)')
     &' **************************************************************',
     &' *                                                            *',
     &' *  Final ',CCFORM(1:4),' energy:',ECC,       ' E_h           *',
     &' *                                                            *',
     &' **************************************************************'
        else
          write(6,*) 'ccform: ',ccform(1:3)
          write(6,*) 'I am quite surprised! Please edit my source code!'
          stop 'programmer''s error in OPTIM_CC2'
        end if
      end if
*
      write(6,'(2(/x,a))')
     & ' Convergence Information:',
     & ' ========================'
      if (ccform(1:3).eq.'TCC'.and.i_do_newccv.eq.1) then
        write(6,'(3(/x,a),/,2(/x,a))')
     & '    E_CC(T)     = <R| H e^T|R>',
     & '    L_CC(T,T)   = <R|e^{T^+}P e^{-T} H e^T|R>'//
     &                                          '/<R|e^{T^+} P e^T|R>|',
     & '    Omega       = <R|tau^+ e^{-T} H e^T|R> ',
     & '  Iter         E_CC(T)        L_CC(T,T)     |Omega|  '//
     &                                          '|delta T|      |T|',
     & ' --------------------------------------------------------'//
     &                                              '-----------------'
        do iter = 1, imacit
          write(6,'(2x,i4,2x,f16.8,x,f16.8,3x,e8.2,3x,e8.2,3x,e8.2)')
     &       iter,eccinf(iter),ecainf(iter),
     &       vfninf(iter),dtninf(iter),ampinf(iter)
        end do
        write(6,'(x,a,//)')
     &' --------------------------------------------------------'//
     &                                              '-----------------'
      else
        write(6,'(/,2(/x,a))')
     & '  Iter         E_CC         |Omega|  |delta T|      |T|',
     & ' --------------------------------------------------------'
        do iter = 1, imacit
          write(6,'(2x,i4,2x,f16.8,3x,e8.2,3x,e8.2,3x,e8.2)')
     &       iter,eccinf(iter),
     &       vfninf(iter),dtninf(iter),ampinf(iter)
        end do
        write(6,'(x,a,//)')
     & ' --------------------------------------------------------'
      end if
*. For transfer to main program a
      E_FINAL = ECCINF(IMACIT)
      ERROR_NORM_FINAL=VFNINF(IMACIT)
      IF(ICONV.EQ.0) THEN
        CONV_F = .FALSE.
      ELSE
        CONV_F = .TRUE.
      END IF
*
      WRITE(6,*) ' Test: E_FINAL, ERROR_NORM_FINAL, CONV_F =',
     &                   E_FINAL, ERROR_NORM_FINAL, CONV_F
*
      if(ntest.ge.100) then
        write(6,*) ' Final set of CC-amplitudes '
        call vec_from_disc(ccvec1,n_cc_amp,1,lblk,luccamp)
        call wrt_cc_vec2(ccvec1,lu,cctype)
      end if
*
      write(6,*) ' Largest T(D,L) block in use = ', ntdl_max_act
      call memman(idummy,idummy,'FLUSM ',idummy,'CC_OPT')

      return
      end


      subroutine init_vec(lu,namp,vec,izero_t1,lrestart)

      implicit none

      logical, intent(in) ::
     &     lrestart
      integer, intent(in) ::
     &     lu, namp, izero_t1
      real(8), intent(inout) ::
     &     vec(namp)

      integer ::
     &     namp_read, lblk

      lblk = -1
      
      if (.not.lrestart) then
        write(6,*) 'no restart indicated ...'
        write(6,*) 'I have to start from scratch (sorry) ...'
        goto 100
      end if

      rewind(lu)
      read(lu,err=100,end=100) namp_read

      ! reuse the file as is
      if (namp_read.eq.namp.and.izero_t1.ne.1) return

      namp_read = min(namp,namp_read)
      ! read as much as there is and set rest to zero
      vec(1:namp) = 0d0
      call vec_from_disc(vec,namp_read,1,lblk,lu)
      
      if (izero_t1.eq.1)
     &     call zero_t1(vec)

      goto 200

 100  continue

      ! init with zero's
      vec(1:namp) = 0d0

 200  continue

      call vec_to_disc(vec,namp,1,lblk,lu)

      return
      end 

c $Id$
