      SUBROUTINE LUCIA_HF(IHFSM,IHFSPC,MAXIT_HF,
     &                    E_HF,E1_FINAL,CONVER_HF)
*
* Master routine for Hartree-Fock optimization using LUCIA
*
* Written for playing around with the HF optimization, Fall of 03
*
* New version allowing use of supersymmetry, May 23, 2012
* Last revision, Sept 4 2012, Jeppe Olsen, cecore.inc added
*
* Final orbitals are saved in KMOAOIN+KMOAOUT+KMOAO_CUR and integrals are 
* transformed to this basis
*
* Jeppe Olsen
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc' 
      INCLUDE 'cecore.inc'
*
      LOGICAL CONVER_HF
      CONVER_HF = .FALSE.
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HF_SEC')
*
      WRITE(6,*) ' ================================== '
      WRITE(6,*) '                                    '
      WRITE(6,*) ' Hartree-Fock optimization section  '
      WRITE(6,*) '                                    '
      WRITE(6,*) '            Version of May  2012    '
      WRITE(6,*) '                                    '
      WRITE(6,*) ' ================================== '
      WRITE(6,*)
      WRITE(6,'(A,I3)') 
     &' Largest allowed number of iterations ', MAXIT_HF
* 
*. Number of occupied orbitals per symmetry 
*
      CALL GET_NOC
*. Obtain Info on Hamiltonian in AO basis 
*. S,H,D,F in AO basis 
      LEN1E = NTOOB **2
      LENC =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0) 
      CALL MEMMAN(KSAO,LENC,'ADDL  ',2,'S_AO  ')
      CALL MEMMAN(KHAO,LENC,'ADDL  ',2,'H_AO  ')
      CALL MEMMAN(KDAO,LENC,'ADDL  ',2,'D_AO  ')
      CALL MEMMAN(KFAO,LENC,'ADDL  ',2,'D_AO  ')
*. Current MO-AO expansion  matrix
      CALL MEMMAN(KCCUR,LENC,'ADDL  ',2,'CCUR  ')
*. Extra copy of MO-AO matrix
      CALL MEMMAN(KLCMOAO2,LENC,'ADDL  ',2,'MOAO2 ')
*. And for saving old(er) MO-AO expansion
      CALL MEMMAN(KCOLD,LENC,'ADDL  ',2,'COLD  ')
*. Saving all AO densities in expanded form
      CALL MEMMAN(KDAO_COLLECT,LEN1E*MAXIT_HF,'ADDL  ',2,'DAO_CO')
*. Save all AO Fock matrices in expanded form 
      CALL MEMMAN(KFAO_COLLECT,LEN1E*MAXIT_HF,'ADDL  ',2,'FAO_CO')
*. And space for the energies in the various iterations 
      CALL MEMMAN(KEITER,(MAXIT_HF+1),'ADDL  ',2,'EITER ')

*. Integral treatment 
*. Obtain AO integrals HAO, SAO
C     GETHSAO(HAO,SAO,IGET_HAO,IGET_SAO)
      CALL GET_HSAO(WORK(KHAO),WORK(KSAO),1,1)
      I2E_AOINT_CORE = 1
      IF(I2E_AOINT_CORE .EQ.1 ) THEN
        WRITE(6,*) ' All AO-integrals are stored in core '
        CALL GET_H2AO
      END IF
*
*. Initial guess to MO-coefficients
*
* INI_HF_MO = 2 => Read coefficients in from fil LUMOIN
* INI_HF_MO = 1 => Diagonalize one-electron Hamiltonian
      CALL GET_INI_GUESS(WORK(KCCUR),INI_HF_MO)
      IF(MAXIT_HF.EQ.0) GOTO 3006
*. And do the optimization 
*. Roothaan-Hall or EOPD
      IF(IHFSOLVE.EQ.1.OR.IHFSOLVE.EQ.2) THEN
        CALL OPTIM_SCF(MAXIT_HF,E_HF,CONVER_HF,E1_FINAL,NIT_HF,
     &                WORK(KEITER))
      ELSE IF(IHFSOLVE.EQ.3) THEN
C            OPTIM_SCF_USING_ONE_STEP(MAXIT_HF)       
        CALL OPTIM_SCF_USING_ONE_STEP(MAXIT_HF,
     &       E_HF,CONVER_HF,E1_FINAL,NIT_HF,WORK(KEITER))
      ELSE 
*. Second order method
        CALL OPTIM_SCF_USING_SECORDER(MAXIT_HF,
     &       E_HF,CONVER_HF,E1_FINAL,NIT_HF,WORK(KEITER))
      END IF
*
* We now have the final orbitals in KCCUR.  If supersymmetry,
* is active, the ordering inside a given symmetry is important.
* In  this case, the orbitals in KCCUR are ordered according to the *DVS* arrays
* and is therefore not in standard order. The question is now
* what order, the orbitals should be returned in.
* The orbitals are ordered depending on 
* whether the HF is the only calculation or not:
*   Only one calculation: Standard order 
*   two or more calulations: Order according to the GAS_SP
*
      IF(I_USE_SUPSYM.EQ.1) THEN
*
*
       WRITE(6,*) 
     & ' Jeppe Testing, bring MO''s in supersymmetry and shell order'
       I_DO_ATEST = 0
       IF(I_DO_ATEST.EQ.1) THEN
*
*. Reorder to standard (symmetry-order)
         CALL REO_CMOAO(WORK(KCCUR),WORK(KLCMOAO2),
     &        WORK(KMO_STA_TO_ACT_REO),0,2)
*. Reorder to supersymmetry order
         CALL REFORM_CMO_STA_GEN(WORK(KLCMOAO2),WORK(KMOAOIN),0,IDUM,1)
*. Reorder to shell order
         CALL REFORM_CMO_SUP_SHL(WORK(KMOAOIN),WORK(KLCMOAO2),1)
*. And print
         CALL PRINT_CSHELL(WORK(KLCMOAO2))
*
* Some tests of reforms: 1:  Actual => Standard => Actual
C REFORM_CMO(C_IN,IFORM_IN, C_OUT, IFORM_OUT)
         WRITE(6,*) ' REFORM: actual => standard => actual '
         CALL REFORM_CMO(WORK(KCCUR),2,WORK(KLCMOAO2),2)
       END IF
  
       WRITE(6,*)
       WRITE(6,*) ' Print shell form using general routines '
       WRITE(6,*) ' ======================================== '
       WRITE(6,*)
       CALL PRINT_CMO_AS_SHELLS(WORK(KCCUR),2)
*
*. return to standard order and save in KCCUR
*
       CALL REO_CMOAO(WORK(KCCUR),WORK(KMOAOIN),
     &      WORK(KMO_STA_TO_ACT_REO),1,2)
*. If there gas calculation in action, sort to the 
*. order specified by NGAS_IRREP_SUPSYM
       IF(I_DO_GAS.EQ.0) THEN
         WRITE(6,*) 
     &   ' Final HF orbitals have standard super-symmetry order'
         CMO_ORD = 'STA'
       ELSE 
         CALL ORDER_GAS_SUPSYM_ORBITALS
         CALL REO_CMOAO(WORK(KCCUR),WORK(KMOAOIN),
     &        WORK(KMO_STA_TO_ACT_REO),1,1)
         WRITE(6,*) 
     &  ' Final HF orbitals ordered according to GASspecs'
        CMO_ORD = 'OCC'
       END IF
*
      END IF
*
 3006 CONTINUE
*
      WRITE(6,*) ' Final list of MO coefficients '
      WRITE(6,*) ' =============================='
C     CALL APRBLM2(WORK(KCCUR),NTOOBS,NTOOBS,NSMOB,0)
      CALL PRINT_CMOAO(WORK(KCCUR))
*. And save in KMOAOIN, KMOAOUT, KMOAO_ACT ...
      CALL COPVEC(WORK(KCCUR),WORK(KMOAOIN),LENC)
      CALL COPVEC(WORK(KCCUR),WORK(KMOAOUT),LENC)
      CALL COPVEC(WORK(KCCUR),WORK(KMOAO_ACT),LENC)
*
* It could very easily that the below should be moved outside
*
*  Set MO to current HF i.e.
*      a) Save MOs in MOAOIN
*      b) Transform integrals to MOAOIN
*      c) Set MOMO to 1
*      d) Construct new inactive Fock matrix
*
* In case of supersymmetry, the order of the MO's depends 
* on whether the HF is the only calculation or not:
*   Only one calculation: Standard order 
*   two or more calulations: Order according to the GAS_SP

* NCMBSPC
*. Transform integrals to converged (at least final) basis '
*. Flag type of integral list to be obtained: Pt complete list of integrals
      IE2LIST_A = IE2LIST_FULL
      IOCOBTP_A = 1
      INTSM_A = 1
      IH1FORM = 1
      IH2FORM = 1
*
      KKCMO_I = KMOAOIN
      KKCMO_J = KMOAOIN
      KKCMO_K = KMOAOIN
      KKCMO_L = KMOAOIN
      CALL TRAINT
      WRITE(6,*) ' Integral transformation completed '
*. And overwrite two-electron integrals
*. Move 2e- integrals to KINT_2EMO
      IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
      NINT2_F = NINT2_G(IE2ARR_F)
      KINT2_F = KINT2_A(IE2ARR_F)
      CALL COPVEC(WORK(KINT2_F),WORK(KINT_2EMO),NINT2_F)
C?      WRITE(6,*) ' NINT2_F = ', NINT2_F
C?      WRITE(6,*) ' Integrals transformed to KINT_2EMO'
C?    CALL WRTMAT(WORK(KINT_2EMO),1,NINT2_F,1,NINT2_F)
*. one-electron integrals to KINT1O
      CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
*. And to KH
      CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
*. The integrals corresponds now to the new initial orbitals, reset MOMO
*. matrix to one
      ZERO = 0.0D0
      CALL SETVEC(WORK(KMOMO),ZERO,LENC)
*
      ONE = 1.0D0
C          SETDIA_BLM(B,VAL,NBLK,LBLK,IPCK)
      CALL SETDIA_BLM(WORK(KMOMO),ONE,NSMOB,NTOOBS,0)
*
*
*. Construct inactive Fock matrix
*
      IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
      KINT2_FSAVE = KINT2_A(IE2ARR_F)
      KINT2_A(IE2ARR_F) = KINT_2EMO
C          FI_FROM_INIINT(FI,CINI,H,EINAC,IHOLETP)
      CALL FI_FROM_INIINT(WORK(KFI),WORK(KMOMO),WORK(KH),
     &                    ECORE_HEX,3)
      ECORE = ECORE_ORIG + ECORE_HEX
      WRITE(6,*) ' TEST: ECORE, ECORE_ORIG, ECORE_HEX = ',
     &                   ECORE, ECORE_ORIG, ECORE_HEX
      CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' ECORE_ORIG, ECORE_HEX, ECORE(2) ',
     &               ECORE_ORIG, ECORE_HEX, ECORE
      END IF
*. and   redirect integral fetcher back to actual integrals
      KINT2 = KINT_2EMO
      KINT2_A(IE2ARR_F) = KINT2_FSAVE
*
*. Print summary
      CALL PRINT_SUMMARY_HF(WORK(KEITER),CONVER_HF,NIT_HF)
*. Transform integrals to current basis
      RETURN
      END 
        
      SUBROUTINE GET_HSAO(HAO,SAO,IGET_HAO,IGET_SAO)
*
* Obtain list of AO integrals in lower half packed form from environment 
* IGET_HAO = 1 => get 1-electron integrals
* IGET_SAO = 1 => get overlap integrals '
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      DIMENSION HAO(*), SAO(*)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
        WRITE(6,*) ' Info from GET_HSAO '
        WRITE(6,*) ' ================== '
        WRITE(6,*) 
      END IF
*
*. Dalton environment
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GETHSA')
*. We only have information on the MO basis, CMO**(-1) will be generated
*. to transform to the AO basis
C                 NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN_H = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLH1,LEN_H,'ADDL  ',2,'H1AO  ')
      CALL MEMMAN(KLH2,LEN_H,'ADDL  ',2,'H2AO  ')
      CALL MEMMAN(KLH3,2*LEN_H,'ADDL  ',2,'H3AO  ')
      CALL MEMMAN(KLH4,LEN_H,'ADDL  ',2,'H4AO  ')
*
      IF(IGET_HAO.EQ.1) THEN
*. C ** -1 in KLH2
        CALL COPVEC(dbl_mb(KMOAOIN),dbl_mb(KLH1),LEN_H)
C            INV_BLKMT(A,AINV,SCR,NBLK,LBLK,IPROBLEM)
        CALL INV_BLKMT(dbl_mb(KLH1),dbl_mb(KLH2),dbl_mb(KLH3),
     &       NSMOB,NTOOBS,IPROBLEM)
        IF(IPROBLEM.NE.0) THEN
          WRITE(6,*) ' GET_HSAO: Warning: problem in block-inverse '
        END IF
*. Read integrals in MO basis in HAO
        REWIND LU2INT
        READ(LU2INT,'(E22.15)') (dbl_mb(KLH1-1+INT1),INT1=1,NINT1)
*. Obtain HAO = C** -1(T) HMO C** -1, in packed form 
C            TRAN_SYM_BLOC_MAT3(AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
        CALL TRAN_SYM_BLOC_MAT3(dbl_mb(KLH1),dbl_mb(KLH2),
     &       NSMOB,NTOOBS,NTOOBS,HAO,dbl_mb(KLH3),1)
      END IF
*
      IF(IGET_SAO.EQ.1) THEN
*. Obtain SAO as ( C CT) ** -1
C                 NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
        LEN_H = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
        CALL MEMMAN(KLH1,LEN_H,'ADDL  ',2,'H1AO  ')
        CALL MEMMAN(KLH2,LEN_H,'ADDL  ',2,'H2AO  ')
        CALL MEMMAN(KLH3,2*LEN_H,'ADDL  ',2,'H3AO  ')
        IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input C matrix '
          CALL APRBLM2(dbl_mb(KMOAOIN),NTOOBS,NTOOBS,NSMOB,0)
        END IF
* C CT in H1
C              MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
        CALL MULT_BLOC_MAT(dbl_mb(KLH1),dbl_mb(KMOAOIN),dbl_mb(KMOAOIN),
     &       NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,2)
* ( C CT) ** -1 in H2
C              INV_BLKMT(A,AINV,SCR,NBLK,LBLK,IPROBLEM)
        CALL INV_BLKMT(dbl_mb(KLH1),dbl_mb(KLH2),dbl_mb(KLH3),
     &       NSMOB,NTOOBS,IPROBLEM)
        IF(IPROBLEM.NE.0) THEN
          WRITE(6,*) ' Warning: GET_HSAO, problem in block-inverse '
        END IF
*.pack to lower half
C               TRIPAK_AO_MAT(AUTPAK,APAK,IWAY)
        CALL TRIPAK_AO_MAT(dbl_mb(KLH2),SAO,1)
      END IF
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GETHSA')
*
      IF(NTEST.GE.100) THEN
        IF(IGET_HAO.EQ.1) THEN
          WRITE(6,*)
          WRITE(6,*) ' H(ao): '
          CALL APRBLM2(HAO,NTOOBS,NTOOBS,NSMOB,1)
        END IF
        IF(IGET_SAO.EQ.1) THEN
          WRITE(6,*)
          WRITE(6,*) ' S(ao): '
          CALL APRBLM2(SAO,NTOOBS,NTOOBS,NSMOB,1)
        END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE GETHAO_DALTON(HAO)
*
* Obtain one-electron integrals from Dalton environment 
*
* Jeppe Olsen, Sept 2003
*
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Local scratch
      PARAMETER ( LBUF = 600 )
      DIMENSION BUF(LBUF), IBUF(LBUF)
*
      DIMENSION HAO(*)
*. Number of integrals in symmetrypacked triangular matrix
C           NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      ZERO = 0.0D0
      CALL SETVEC(HAO,ZERO,LEN)
C
C
C     Read information on file AONEINT from HERMIT.
C
      ITAP34 = 66
      OPEN (ITAP34,STATUS='OLD',FORM='UNFORMATTED',FILE='AOONEINT')
*
      CALL MOLLAB('ONEHAMIL',ITAP34,6)
 2100 READ (ITAP34) (BUF(I),I=1,LBUF),(IBUF(I),I=1,LBUF),LENGTH
      DO 2200 I = 1,LENGTH
         HAO(IBUF(I)) = BUF(I)
 2200 CONTINUE
      IF (LENGTH .GE. 0) GO TO 2100
*
      CLOSE(ITAP34,STATUS='KEEP')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' 1-e integrals in AO basis '
C        WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
        CALL WRTVH1(HAO,1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      END IF
*
      RETURN
      END
      SUBROUTINE GETSAO_DALTON(SAO)
*
* Obtain overlap integrals from Dalton environment 
*
* Jeppe Olsen, Sept 2003
*
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Local scratch
      PARAMETER ( LBUF = 600 )
      DIMENSION BUF(LBUF), IBUF(LBUF)
*
      DIMENSION SAO(*)
*. Number of integrals in symmetrypacked triangular matrix
C           NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      ZERO = 0.0D0
      CALL SETVEC(SAO,ZERO,LEN)
C
C
C     Read information on file AONEINT from HERMIT.
C
      ITAP34 = 66
      OPEN (ITAP34,STATUS='OLD',FORM='UNFORMATTED',FILE='AOONEINT')
*
      REWIND(ITAP34)
      CALL MOLLAB('OVERLAP ',ITAP34,6)
 2100 READ (ITAP34) (BUF(I),I=1,LBUF),(IBUF(I),I=1,LBUF),LENGTH
      DO 2200 I = 1,LENGTH
         SAO(IBUF(I)) = BUF(I)
 2200 CONTINUE
      IF (LENGTH .GE. 0) GO TO 2100
*
      CLOSE(ITAP34,STATUS='KEEP')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap integrals in AO basis '
C        WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
        CALL WRTVH1(SAO,1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_H2AO
*
* Obtain AO 2-electron integrals in core  and save in KINT_2EMO
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cintfo.inc'
*
      IF(ENVIRO(1:6).EQ.'DALTON'.OR.ENVIRO(1:5).EQ.'LUCIA') THEN
         CALL GET_H2AO_INNER
      ELSE 
         WRITE(6,'(A,A)') ' Unknown environment parameter ', ENVIRO
      END IF
*
      RETURN
      END
      SUBROUTINE GET_H2AO_INNER
*
* 2 e Ints in AO basis
*
*. are obtained by back transforming MO integrals !!
*
* It is assumed that integrals in MO basis have been read in
*
* Updated, May 2012
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
*. Output : in WORK(KINT_2EMO)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' Info from GET_H2AO_INNER '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
       END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'H2AO_D')
*
*. Obtain C**(-1) as C(T) * S
*
C  NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
*. Obtain CMOAO
      LENE = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLC,LENE,'ADDL  ',2,'CMOAO ')
      CALL GET_CMOAO_ENV(WORK(KLC))
*. Expand S to  full form
      CALL MEMMAN(KLSE,LENE,'ADDL  ',2,'CMOAO ')
      CALL TRIPAK_BLKM(WORK(KLSE),WORK(KSAO),2,NAOS_ENV,NSMOB)
C          TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
*. and Multiply C(T) and S
      CALL MEMMAN(KLCINV,LENE,'ADDL  ',2,'CINV ')
      CALL MULT_BLOC_MAT(WORK(KLCINV),WORK(KLC),WORK(KLSE),
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,1)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' C(inv) matrix '
        CALL WRTVH1(WORK(KLCINV),1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      END IF
*
* Perform integral transformation to AO basis ...
* Tell program to work with full two-electron integral list 
      IE2LIST_A = IE2LIST_FULL
      IOCOBTP_A = 1
      INTSM_A = 1
      CALL PREPARE_2EI_LIST
      I12S_A = 1
      I34S_A = 1
      I1234S_A = 1
      CALL FLAG_ACT_INTLIST(IE2LIST_FULL)
*
      KKCMO_I = KLCINV
      KKCMO_J = KLCINV
      KKCMO_K = KLCINV
      KKCMO_L = KLCINV
      IH1FORM = 1
      IH2FORM = 1
      CALL TRAINT
      WRITE(6,*) ' Integral transformation to AO basis completed '
*. And overwrite two-electron integrals
*. Move 2e- integrals to KINT_2EMO
      IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
      NINT2_F = NINT2_G(IE2ARR_F)
      KINT2_F = KINT2_A(IE2ARR_F)
      CALL COPVEC(WORK(KINT2_F),WORK(KINT_2EMO),NINT2_F)
C?    WRITE(6,*) ' NINT2_F = ', NINT2_F
C?    WRITE(6,*) ' Integrals transformed to KINT_2EMO'
C?    CALL WRTMAT(WORK(KINT_2EMO),1,NINT2_F,1,NINT2_F)
*. one-electron integrals to KINT1O
      CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
*. And to KH
        CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
*. The integrals corresponds now to the new initial orbitals, reset MOMO
*. matrix to one
C            SETDIA_BLM(B,VAL,NBLK,LBLK,IPCK)
        ONE = 1.0D0
      CALL SETDIA_BLM(WORK(KMOMO),ONE,NSMOB,NTOOBS,0)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'H2AO_D')
*
      RETURN
      END
      SUBROUTINE GET_INI_GUESS(CINI,INI_MO)
*
*. Obtain initial MO-AO expansion matrix 
*
* INI_MO = 1 => Read coefficients in from fil LUMOIN
* INI_MO = 2 => Diagonalize one-electron Hamiltonian
*
* Jeppe Olsen, Sept. 2003
*              Modified to use general symmetry, May 24, 2012
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      IDUM = 0
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GET_INI_GUESS '
        WRITE(6,*) ' ======================='
        WRITE(6,*) 
        WRITE(6,*) ' INI_MO = ', INI_MO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'INI_MO')
 
      LENC   = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0)
      IF(INI_MO.EQ.2) THEN
         WRITE(6,*) ' Initial MO-AO transformation readin '
*. Is already stored in KMOAOIN, so just copy
C        CALL GET_CMOAO_ENV(CINI)
*. Reorder - in case of supersymmetry
         CALL REO_CMOAO(WORK(KMOAOIN),CINI,
     &        WORK(KMO_STA_TO_ACT_REO),0,1)
C            REO_CMOAO(CIN,COUT,IREO,ICOPY,IWAY)
COLD     CALL COPVEC(WORK(KMOAOIN),CINI,LENC)
      ELSE IF (INI_MO.EQ.1) THEN
*. Obtain Initial MO-coefficients by diagonalizing H matrix 
*. in blocks defined by general symmetry
*. Local copies of S and H
C        LEN_BLMAT(NBLK,LROW,LCOL,IPACK)
         LEN_HS = LEN_BLMAT(NGENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,1)
         LEN_HSU = LEN_BLMAT(NGENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,0)
         WRITE(6,*) ' LEN_HS = ', LEN_HS  
         WRITE(6,*) ' NGENSMOB= ', NGENSMOB
         WRITE(6,*) ' NBAS_GENSMOB: '
         CALL IWRTMA3(NBAS_GENSMOB,1,NGENSMOB,1,NGENSMOB)
         CALL MEMMAN(KLH,LEN_HS,'ADDL  ',2,'HLOC  ')
         CALL MEMMAN(KLS,LEN_HS,'ADDL  ',2,'SLOC  ')
         CALL MEMMAN(KLC,LENC,'ADDL  ',2,'CLOC  ')
*. Reform from standard to general symmetry 
C             REFORM_MAT_STA_GEN(ASTA,AGEN,IPACK,IWAY)
         CALL REFORM_MAT_STA_GEN(WORK(KHAO),WORK(KLH),1,1)
         CALL REFORM_MAT_STA_GEN(WORK(KSAO),WORK(KLS),1,1)
*.space for eigenvalues
         NAOT = IELSUM(NBAS_GENSMOB,NGENSMOB)
         CALL MEMMAN(KLEPS,NAOT,'ADDL  ',2,'EPS   ')
*
         CALL SCF_DIA(WORK(KLH),WORK(KLS),WORK(KLC),WORK(KLEPS))
*. The C-coefficient are given in the general symmetry basis,
*  transfer to standard
C     REFORM_CMO_STA_GEN(CMO_STA,CMO_GEN,
C    &           IDO_REORDER,IREO,IWAY)
         CALL REFORM_CMO_STA_GEN(CINI,WORK(KLC),0,IDUM,2)
*. And reorder internally in each symmetry
C           REO_CMOAO(CIN,COUT,IREO,ICOPY,IWAY)
         CALL REO_CMOAO(CINI,WORK(KLC),WORK(KMO_STA_TO_ACT_REO),1,1)
C. 

      END IF
* 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial MO-AO expansion matrix (gen sym)'
        CALL APRBLM2(WORK(KLC),NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,0)
        WRITE(6,*) ' Initial MO-AO expansion matrix (standard)'
        CALL APRBLM2(CINI,NTOOBS,NTOOBS,NSMOB)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'INI_MO')
      IF(NTEST.GE.100)  WRITE(6,*) ' Leaving GET_INI_GUESS '
*
      RETURN
      END
      SUBROUTINE SCF_DIA_CTL(F,S,C,E,D,COLD,WOCOCMIN,ALPHA_INI,NOCPSM)
*. Solve the SCF equations with the restriction that 
*. the projection of each occupied orbital in the space of 
*. old occupied orbitals must atleast have the weight WOCOCMIN
*. Occupied orbitals stored in COLD.
*. If the change is larger, a constant times the density (SDS)
*. is added, and the constant is choosen, so the 
*. above criteria eventually is fulfilled
*. If ALPHA_INI .ne. 0 this value is used as initial shift
* 
*
* Jeppe Olsen, Homeward bound from the Gwatt meeeting 
* Oct. 2. 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION F(*),S(*),COLD(*),D(*)
*. Number of occupied orbitals per symmetry
      DIMENSION NOCPSM(*)
*. Output
      DIMENSION C(*),E(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'SCF_CT')
*. A bit of local scratch
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      CALL MEMMAN(KLSCR,LENC,'ADDL  ',2,'SCR   ')
      CALL MEMMAN(KLCOSCN,LENC,'ADDL  ',2,'COSCN ')
      CALL MEMMAN(KLW,NTOOB,'ADDL  ',2,'WCO_CN')
      CALL MEMMAN(KLSDS,LENC,'ADDL  ',2,'SDS   ')
      CALL MEMMAN(KLS_E,LENC,'ADDL  ',2,'SDS   ')
*
*. Loop over iterations of shifts
      MXSHFTIT = 5
      ALPHA_TOT = 0.0D0
      I_HAVE_DONE_SDS = 0
      ALPHA = 0.0D0
      DO ISHFTIT = 1, MXSHFTIT
        IF(ISHFTIT.EQ.1) THEN 
          ALPHA = ALPHA_INI
          ALPHA_DEL = ALPHA_INI
        ELSE 
          ALPHA = ALPHA + ALPHA_DEL 
        END IF
*
        IF(ALPHA_DEL.NE.0.0D0) THEN
          IF(I_HAVE_DONE_SDS.EQ.0) THEN
*. Construct SDS : Unpack S, CALC SDS, and pack
*. S(unpack) in KLS_E
            CALL TRIPAK_BLKM(WORK(KLS_E),S,2,NAOS_ENV,NSMOB)
*. DS in KLSCR
            CALL MULT_BLOC_MAT(WORK(KLSCR),D,WORK(KLS_E),
     &      NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &      1)
*. SDS in KLSDS
            CALL MULT_BLOC_MAT(WORK(KLSDS),WORK(KLS_E),WORK(KLSCR),
     &      NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &      1)
            CALL TRIPAK_BLKM(WORK(KLSDS),WORK(KLSCR),1,NAOS_ENV,NSMOB)
            CALL COPVEC(WORK(KLSCR),WORK(KLSDS),LENP)
            I_HAVE_DONE_SDS = 1
          END IF
*. F => F + ALPHA_DEL * D
          ONE = 1.0D0
          CALL VECSUM(F,F,WORK(KLSDS),ONE,ALPHA_DEL,LENP)
        END IF
*. Diagonalize
        CALL SCF_DIA(F,S,C,E)
*. Find overlap between new and old vectors 
        CALL OVERLAP_COLD_CNEW(COLD,C,S,WORK(KLCOSCN),WORK(KLSCR),
     &       NAOS_ENV,NOCPSM,NSMOB)
*. Find weight of old occupied MO's in new MO's
        CALL WEIGHT_NEW_OCCOLD(WORK(KLCOSCN),NOCPSM,WORK(KLW),WMIN)
*
        IF(ISHFTIT.LT.MXSHFTIT.AND.WMIN.LT.WOCOCMIN) THEN
*. Obtain HOMO-LUMO gap ( is actually obtained for shifted F, 
*. but that is fine
C              GET_HOMO_LUMO(EPSIL,NAOS,NOCC,NSMOB,EHOMO,ELUMO)
          CALL GET_HOMO_LUMO(E,NAOS_ENV,NOCPSM,NSMOB,
     &         EHOMO,ELUMO,EHL_GAP)
*. We have the weight of the old occupied orbitals in each new 
* in WORK(KLW) and we have the orbital energies in E.
*. If the minimal weight is too small just add HOMO-LUMO-GAP
C         ALPHA_DEL = MIN(-EHL_GAP,-0.5D0)
*. The F matrix will be shifted to reduce the steps.
*. It is assumed for now that E(HOMO) .LT. E(LUMO). that is
*. the MO's have been chosen according to the aufbau principle
C         ALPHA_DEL = -SQRT((WOCOCMIN-WMIN))/(EHL_GAP)
C          ALPHA_DEL = -0.8D0*((1-SQRT(WMIN))/(1-SQRT(WOCOCMIN))-1)*
C    &                        EHL_GAP
C         ALPHA_DEL=-1.5D0*(SQRT(1.0D0-WMIN)/SQRT(1.0D0-WOCOCMIN)-1)*
C    &                     EHL_GAP
C     SELECT_RH_SHIFT(W,WMIN,E,NOBPSM,NOCPSM,NSMOB,SHIFT)
          CALL SELECT_RH_SHIFT(WORK(KLW),WOCOCMIN,E,NAOS_ENV,NOCPSM,
     &                         NSMOB,SHIFT)
          ALPHA_DEL = SHIFT
          WRITE(6,*) ' ALPHA_DEL = ', ALPHA_DEL
      
        ELSE
*. We are finished so 
         GOTO 1001
        END IF
      END DO
1001  CONTINUE
*. Clean up : remove SDS term from shifted F
      CALL VECSUM(F,F,WORK(KLSDS),ONE,-ALPHA,LENP)
          
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'SCF_CT')
      RETURN
      END 
*
      SUBROUTINE SCF_DIA(F,S,C,E)   
*
* Solve SCF equations for given matrices 
*
*     F C = S C E
*
* F and S : input
* C and E : output 
*
*. Obtained from Dage Sundholms version in QDOT program 
*
* Jeppe Olsen, Sept. 2003
*. Modified to general symmetry division, May 24, 2012
*. General input
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION F(*),S(*)
*. output
      DIMENSION C(*),E(*)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' =================='
        WRITE(6,*) ' SCF_DIA reporting '
        WRITE(6,*) ' =================='
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' INPUT to SCF_DIA '
        WRITE(6,*) ' Blocks of overlapmatrix  '
        CALL APRBLM2(S,NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,1)      
        WRITE(6,*) ' Blocks of fock matrix '
        CALL APRBLM2(F,NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,1)      
      END IF
*
      IDUMMY = 0
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'SCFDIA') 
*. A bit of scratch 
      LENGTH = NDIM_1EL_MAT(1,NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,0)
      CALL MEMMAN(KLHIN,LENGTH ,'ADDL  ',2,'KLHIN ')
      CALL MEMMAN(KLSIN,LENGTH ,'ADDL  ',2,'KLSIN ')
      CALL MEMMAN(KLPVEC ,LENGTH ,'ADDL  ',2,'PVEC  ')
      CALL MEMMAN(KLVOUT,LENGTH ,'ADDL  ',2,'VOUT  ')
*. Unpack F and S
      ONE = 1.0D0
      CALL TRIPAK_BLKM(WORK(KLSIN),S,2,NBAS_GENSMOB,NGENSMOB)
      CALL TRIPAK_BLKM(WORK(KLHIN),F,2,NBAS_GENSMOB,NGENSMOB)
C          TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Blocks of overlapmatrix  '
        CALL APRBLM2(WORK(KLSIN),NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,0)      
        WRITE(6,*) ' Blocks of fock matrix '
        CALL APRBLM2(WORK(KLHIN),NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,0)      
      END IF
*. Generalized diagonalization (sorted)
      ISORT=1
      CALL GENDIA_BLMAT(WORK(KLHIN),WORK(KLSIN),C,E,
     &                  WORK(KLPVEC),NBAS_GENSMOB,NGENSMOB,ISORT)

*. Report to home
      IF(NTEST.GE.10) THEN
        CALL WRITE_ORBENERGIES(E,NBAS_GENSMOB,NGENSMOB)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Blocks of eigenvectors '
        CALL APRBLM2(C,NBAS_GENSMOB,NBAS_GENSMOB,NGENSMOB,0)
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'SCFDIA') 
*
      RETURN
      END
      SUBROUTINE DIAG_BLK_SYMMAT(A,NBLK,LBLK,X,EIGENLV,SCR,ISYM)
*
* Diagonalize blocked symmetric matrix. 
* 
*
* Last modification: July 8, 2012, JO
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION A(*)
      INTEGER LBLK(*)
*. Output X : eigenvectors, eigenv : eigenvalues
      DIMENSION X(*), EIGENLV(*)
*. Scratch ( dim of largest block)
      DIMENSION SCR(*)
*
      IOFFA = 0
      IOFFX = 0
      IOFFE = 0
      DO IBLK = 1, NBLK
*. Offsets to current block
        IF(IBLK.EQ.1) THEN
         IOFFA = 1
         IOFFX = 1
         IOFFE = 1
        ELSE
         IF(ISYM.EQ.1) THEN
           IOFFA = IOFFA + LBLK(IBLK-1)*(LBLK(IBLK-1)+1)/2
         ELSE
           IOFFA = IOFFA + LBLK(IBLK-1) ** 2 
         END IF
         IOFFX = IOFFX + LBLK(IBLK-1) ** 2
         IOFFE = IOFFE + LBLK(IBLK-1)
        END IF
        LEN = LBLK(IBLK)
        IF(ISYM.EQ.1) THEN
          CALL COPVEC(A(IOFFA),SCR,LEN*(LEN+1)/2)
        ELSE
          CALL TRIPAK(A(IOFFA),SCR,1,LEN,LEN)      
        END IF
*. block is now in symmetric packed form in SCR, diagonalize
*DS     WRITE(6,*) ' Input matrix to EIGENL '
*DS     CALL PRSYM(SCR,LEN)
*DS  no sort = 0,0)  CALL EIGENL(SCR,X(IOFFX),LEN,0,0) 
COLD    write(6,*) 'no sort in DIAG_BLK_SYMMAT'
COLD    CALL EIGENL(SCR,X(IOFFX),LEN,0,0)
*. Sorting reinstated, JO, July 2012
        CALL EIGENL(SCR,X(IOFFX),LEN,0,1) 
*. Eigenvalues are diagonal elements of scr, copy to eigen
        CALL COPDIA(SCR,EIGENLV(IOFFE),LEN,1)
C            COPDIA(A,VEC,NDIM,IPACK)
      END DO
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info about diag. of blocked matrix '
       WRITE(6,*) ' Eigenvectors - in blocked form '
       CALL APRBLM2(X,LBLK,LBLK,NBLK,0   )
C           APRBLM2(AOUT,LBLOCK,LBLOCK,NBLOCK,ISYM)
       WRITE(6,*) ' Eigenvalues ' 
       LEN_TOT = IELSUM(LBLK,NBLK)
       CALL WRTMAT(EIGENLV,1,LEN_TOT,1,LEN_TOT)
      END IF
*
      RETURN
      END
      SUBROUTINE TRAN_SYM_BLOC_MAT_2(AIN,X,NBLOCK,LBLOCK,AOUT,SCR,ISYM)
*
* Transform a blocked symmetric matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
* Daughter of TRAN_SYM_BLOC_MAT : ISYM added to allow use of complete blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LBLOCK(NBLOCK)
*. Output 
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block 
      DIMENSION SCR(*)
*
C?    WRITE(6,*) ' TRAN.... ISYM = ', ISYM
      IOFFP = 0
      IOFFC = 0
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
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN 
C?       WRITE(6,*) ' IBLOCK IOFFP L ',IBLOCK,IOFFP,L
         CALL TRIPAK(SCR(K1),AIN(IOFFP),2,L,L)
       ELSE
         CALL COPVEC(AIN(IOFFC),SCR(K1),L*L) 
       END IF
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
C?     WRITE(6,*) ' TRAN ... IBLOCK,L ', IBLOCK,L
C?     WRITE(6,*) ' TRAN ... : Input to MATML7 '
C?     CALL WRTMAT(X(IOFFC),L,L,L,L)
C?     CALL WRTMAT(SCR(K1) ,L,L,L,L)
       CALL MATML7(SCR(K2),X(IOFFC),SCR(K1),L,L,L,L,L,L,
     &             ZERO,ONE,1)
C?     WRITE(6,*) ' Half Transformed block '
C?     CALL WRTMAT(SCR(K2),l,l,l,l)
*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFC),L,L,L,L,L,L,
     &             ZERO,ONE,0)
C?     WRITE(6,*) ' Transformed block '
C?     CALL WRTMAT(SCR(K1),l,l,l,l)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK(SCR(K1),AOUT(IOFFP),1,L,L)
       ELSE 
         CALL COPVEC(SCR(K1),AOUT(IOFFC),L*L)
       END IF
*
      END DO
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT_2 '
        WRITE(6,*) ' ==================================='
        CALL APRBLM2(AOUT,LBLOCK,LBLOCK,NBLOCK,ISYM)
      END IF
*
      RETURN
      END
      SUBROUTINE SQRT_BLMAT(A,NBLK,LBLK,ITASK,ASQRT,AMSQRT,SCR,ISYM)
* Obtain square root -and inverse square root- of blocked matrix 
*
* input
* =====
*
* A : Input blocked matrix
* NBLK : Number of blocks in A
* LBLK : Length of each block
* ITASK = 1 => Just square root of matrix
*       = 2 => Calculate also inverse square root
* ISYM  = 1 => Matrix is packed as lower half
*       = 0 => Matris is not packed as lower half
*
* Output
* ======
*
* ASQRT : square root in blocked form
* AMSQRT : Inverse square root in blocked form ( IF ITASK = 2 )
*
* Scratch
* ========
*
* SCR : Should atleast be of length : 6 times largest block dim
*         
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION A(*)
      INTEGER LBLK(*)
*. Output
      DIMENSION ASQRT(*),AMSQRT(*)
*. Scratch
      DIMENSION SCR(*)
*
      DO IBLK = 1, NBLK
        IF(IBLK.EQ.1) THEN
          IOFF = 1
        ELSE
          IF(ISYM.EQ.1) THEN
            IOFF = IOFF + LBLK(IBLK-1)*(LBLK(IBLK-1)+1)/2
          ELSE
            IOFF = IOFF + LBLK(IBLK-1)**2
          END IF
        END IF
        LEN = LBLK(IBLK)
C?      WRITE(6,*) ' IBLK IOFF LBLK ', IBLK,IOFF,LEN        
        IF(ISYM.EQ.1) THEN
          K1 = 1     
          K2 = K1 + LEN*LEN
          K3 = K2 + LEN*LEN      
          K4 = K3 + LEN*LEN      
          SIGN = 1.0D0
          CALL TRIPAK(SCR(K1),A(IOFF),2,LEN,LEN)
          CALL SQRTMT(SCR(K1),LEN,ITASK,SCR(K2),SCR(K3),SCR(K4))
          CALL TRIPAK(SCR(K2),ASQRT(IOFF),1,LEN,LEN)
          CALL TRIPAK(SCR(K3),AMSQRT(IOFF),1,LEN,LEN)
        ELSE
          CALL SQRTMT(A(IOFF),LEN,ITASK,ASQRT(IOFF),AMSQRT(IOFF),SCR)
C             SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
        END IF
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SQRT_BLMT '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,*) ' Square root matrix : '
        CALL APRBLM2(ASQRT,LBLK,LBLK,NBLK,ISYM)
C       CALL APRBLM2(WORK(KFOCK),NTOOBS,NTOOBS,NSMOB,ISM)
        IF(ITASK.EQ.2) THEN
          WRITE(6,*) ' Inverse square root matrix '
          CALL APRBLM2(AMSQRT,LBLK,LBLK,NBLK,ISYM)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE SORT_BLOC_MAT(E,C,NBLK,LBLK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(*),E(*),NBLK(LBLK)

      IOFF = 0
      IEOFF = 0
      DO I=1,LBLK
        IF(I.EQ.1) THEN 
         IOFF=1
         IEOFF=1
        ELSE 
         IOFF=IOFF+NBLK(I-1)*NBLK(I-1)
         IEOFF=IEOFF+NBLK(I-1)
       END IF 
       CALL SORTERA(E(IEOFF),C(IOFF),NBLK(I))
      END DO

      RETURN
      END
      SUBROUTINE GENDIA_BLMAT(HIN,SIN,C,E,PVEC,NBLK,LBLK,ISORT)
*
      INCLUDE 'implicit.inc'
*
      DIMENSION HIN(*),SIN(*),C(*),E(*),PVEC(*),NBLK(LBLK)

*
      IOFF = 0
      IEOFF = 0
      DO I=1,LBLK
        IF(I.EQ.1) THEN 
          IOFF=1
          IEOFF=1
        ELSE 
          IOFF=IOFF+NBLK(I-1)*NBLK(I-1)
          IEOFF=IEOFF+NBLK(I-1)
        END IF 

* Solve first
        NOSORT=0
        CALL GENDIA(HIN(IOFF),SIN(IOFF),C(IOFF),E(IEOFF),PVEC(IOFF),
     &              NBLK(I),NOSORT)

        END DO 

* Sort later
      DO I=1,LBLK
        IF(I.EQ.1) THEN 
          IOFF=1
          IEOFF=1
        ELSE 
          IOFF=IOFF+NBLK(I-1)*NBLK(I-1)
          IEOFF=IEOFF+NBLK(I-1)
        END IF 

        IF(ISORT.EQ.1) THEN 
          CALL SORTERA(E(IEOFF),C(IOFF),NBLK(I))
        END IF
      END DO

      RETURN
      END
      SUBROUTINE SORTERA(E,C,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(NDIM),C(NDIM,NDIM)

      DO I=1,NDIM
        DO J=I+1,NDIM
          IF(E(I)-E(J).GT.0.D0) THEN
            X=E(I)
            E(I)=E(J)
            E(J)=X
            DO K=1,NDIM
              X=C(K,I)
              C(K,I)=C(K,J)
              C(K,J)=X
            END DO
          END IF
        END DO
      END DO

      RETURN
      END
      SUBROUTINE GENDIA(HIN,SIN,VOUT,EIGENLV,PVEC,NDIM,ISORT)
*
*   Diagonalization of a symmetric eigenvalue problem
*   Eigenvalues and eigenvectors are obtained
*
*   HIN    = Matrix to be diagonalized, Destroyed
*   SIN    = Metric, Destroyed
*   VOUT   = Eigenvector
*   EIGENLV = Eigenvalues
*   PVEC   = Transformation matrix between the two bases
*   NDIM   = Size of the matrices
*   ISORT  = Sorting parameter (=0 no sorting and =1 sorting)
 
*   Subroutines called:
*                       TRANSH
*                       PACKHM
*                       EIGENL
*                       MATML7
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION HIN(NDIM,NDIM),SIN(NDIM,NDIM),EIGENLV(NDIM)
      DIMENSION VOUT(NDIM,NDIM),PVEC(NDIM,NDIM)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to TRANSH, HIN and SIN '
        CALL WRTMAT(HIN,NDIM,NDIM,NDIM,NDIM)
        CALL WRTMAT(SIN,NDIM,NDIM,NDIM,NDIM)
      END IF
 
* Transform to the diagonal basis (SIN = Unit matrix after transform)
* The transformation matrix is saved in PVEC
* VOUT is used as Trash area in TRANSH
 
      CALL TRANSH_DAGE(NDIM,HIN,SIN,PVEC,VOUT)
 
      IF(NTEST.GT.20) THEN
 
        WRITE(6,*)
        WRITE(6,*) 'TRANSFORMED H-MATRIX'
        WRITE(6,*)
        CALL WRTMAT(HIN,NDIM,NDIM,NDIM,NDIM)
C        DO 10 I=1,NDIM
C          DO 10 J=1,NDIM
C            WRITE(6,1000) HIN(J,I),J,I
C10      CONTINUE
      END IF
 
      IF(NTEST.GT.30) THEN
 
        WRITE(6,*)
        WRITE(6,*) 'TRANSFORMATION MATRIX'
        WRITE(6,*)
        CALL WRTMAT(PVEC,NDIM,NDIM,NDIM,NDIM)
C        DO 20 I=1,NDIM
C          DO 20 J=1,NDIM
C            WRITE(6,1000) PVEC(J,I),J,I
C20      CONTINUE
      END IF
 
* Pack the matrix for diagonalization using subroutine EIGENL
* Pack HIN into the SIN matrix
 
      CALL PACKHM(HIN,SIN,NDIM)
 
* Both eigenvalues (in every I*(I+1)/2 element of SIN) and
* Eigenvectors (unsorted if MFKR = 0 and sorted if MFKR = 1) in hin
 
      MV=0
      MFKR=ISORT
*   
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to EIGENL '
        WRITE(6,*) ' Dimension = ', NDIM
        CALL PRSYM(SIN,NDIM)
      END IF
      CALL EIGENL(SIN,HIN,NDIM,MV,MFKR)
 
      IF(NTEST.GT.40) THEN
 
        WRITE(6,*)
        WRITE(6,*) 'TRANSFORMED EIGENLVECTORS'
        WRITE(6,*)
        CALL WRTMAT(HIN,NDIM,NDIM,NDIM,NDIM)
C        DO 30 I=1,NDIM
C          DO 30 J=1,NDIM
C            WRITE(6,1000) HIN(J,I),J,I
C30      CONTINUE
      END IF
 
* Fetch the eigenvalues from SIN to EIGENLV
 
      DO 40 I=1,NDIM
        INDEX=I*(I+1)/2
        IN=INT(INDEX/NDIM)+1
        JN=INDEX-(IN-1)*NDIM
        IF(JN.EQ.0) THEN
          IN=IN-1
          JN=NDIM
        END IF
        EIGENLV(I)=SIN(JN,IN)
40    CONTINUE
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(VOUT,PVEC,HIN,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,
     &            FACTORC,FACTORAB,0)
 
C     CALL DGEMUL(PVEC,NDIM,'N',HIN,NDIM,'N',VOUT,NDIM,
C    &            NDIM,NDIM,NDIM)
 
      IF(NTEST.GT.1) THEN
 
        WRITE(6,*)
        WRITE(6,*) 'EIGENLVALUES FROM GENDIA'
        WRITE(6,*)
        DO 100 I=1,NDIM
          WRITE(6,1000) EIGENLV(I),I
100     CONTINUE
      END IF
 
      IF(NTEST.GT.10) THEN
 
        WRITE(6,*)
        WRITE(6,*) 'EIGENLVECTORS FROM GENDIA'
        WRITE(6,*)
        CALL WRTMAT(VOUT,NDIM,NDIM,NDIM,NDIM)
C        DO 200 I=1,NDIM
C          DO 200 J=1,NDIM
C            WRITE(6,1000) VOUT(J,I),J,I
C200     CONTINUE
      END IF
*
 1000 FORMAT(1X,D20.12,3X,I6,3X,I6)
 
      RETURN
      END
      SUBROUTINE PACKHM(HIN,HPACK,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION HIN(NDIM,NDIM),HPACK(NDIM*(NDIM+1)/2)
*
*  Set the symmetric HIN-matrix into packed form HPACK
*
      DO 100 I=1,NDIM
        DO 100 J=I,NDIM
          HPACK(J*(J-1)/2+I)=HIN(J,I)
100   CONTINUE
 
      RETURN
      END
      SUBROUTINE TRANSH_DAGE(N,H,S,P,WORK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
C Symmetric matrices are assumed
C Transform H to H' obtain the transformation matrix P
C (HC=ESC) => (H'C'=EC') and (C=PC')
 
C N     : Dimension of the problem
C H     : Hamilton matrix, in H out H' (full matrix)
C S     : Overlap matrix,  in S out I  (full matrix)
C P     : Transformation matrix in trash out P (full matrix)
 
      DIMENSION S(N,N),H(N,N),P(N,N),WORK(N)
 
C Neglect matrix elements less than DEPS
 
      DEPS=0.5D-14
      ONE=1.D0
 
C Set P to unit matrix
 
      CALL SETVEC(P,0.0D0, N ** 2 )
      CALL SETDIA(P,ONE,N,0)
 
C First part of the transformation of the H and P matrices
 
      DO 20 K=1,N-1
        DO 20 J=N,K+1,-1
 
          D=S(K,J)/S(K,K)
          IF(ABS(D).GT.DEPS) THEN
 
            DO 30 I=K+1,J
30          S(I,J)=S(I,J)-D*S(K,I)
            DO 31 I=K+1,J
31          H(I,J)=H(I,J)-D*H(K,I)
 
            DO 40 I=1,K
40          H(I,J)=H(I,J)-D*H(I,K)
 
            DO 50 I=J,N
50          H(J,I)=H(J,I)-D*H(K,I)
 
            DO 60 I=1,K
60          P(I,J)=P(I,J)-D*P(I,K)
 
          END IF
20    CONTINUE
 
C Second part of the transformation obtaining the final H and P matrices
C but just the upper triangle.
 
      DO 70 I=1,N
        E=SQRT(S(I,I))
 
        DO 80 J=1,N
80      H(I,J)=H(I,J)/E
 
        DO 90 J=1,I
90      H(J,I)=H(J,I)/E
 
        DO 100 J=1,I
        P(J,I)=P(J,I)/E
100     CONTINUE
70    CONTINUE
 
C To be sure, copy the upper triangle to the lower triangle
C set the S matrix to be unit matrix
C (Just in case)
 
      DO 200 I=1,N-1
        DO 200 J=I+1,N
200     H(J,I)=H(I,J)
 
      CALL SETVEC(S,0.0D0, N ** 2 )
      CALL SETDIA(S,1.0D0,N,0)
 
      RETURN
      END
      SUBROUTINE GET_FOCK(FOCK,H,P,ENERGY,ENUC,I12) 
*
* Obtain FOCK matrix for given density P
*
* FOCK is delivered in symmetry packed, lower half form
*
* Total energy is also obtained as 1/2 Tr (H + F )
*
* I12 = 1 => only 1e term
* I12 = 2 => only 2e term
* I12 = 3 => 1 and 2 e terms
*. General input
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      DIMENSION H(*)
*. Specific input ( P is delivered as packed complete(not lower half) matrix)
      DIMENSION P(*)
*. Output
      DIMENSION FOCK(*)
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Information from GET_FOCK'
        WRITE(6,*) ' ========================='
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' I12 = ', I12
      END IF
*
      IDUMMY = 0
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ', IDUMMY,'GET_FO')
      CALL LUCIAQENTER('GT_FO')
* Largest number of orbitals belonging to given sym/sys
C IMNMX(IVEC,NDIM,MINMAX)
      MXSSOB = IMNMX(NTOOBS,NSMOB,2)
C?    WRITE(6,*) ' Largest number of orbital of given sym', MXSSOB
      LEN = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
*. Local memory : two and four index block over orbitals 
      CALL MEMMAN(KLM2,MXSSOB**2,'ADDL  ',2,'MAT2I  ')
      CALL MEMMAN(KLM2P,MXSSOB**2,'ADDL  ',2,'MAT2I  ')
      CALL MEMMAN(KLM4,MXSSOB**4,'ADDL  ',2,'MAT4I  ')
      CALL MEMMAN(KLM4B,MXSSOB**4,'ADDL  ',2,'MAT4IB ')
C?    WRITE(6,*) ' KLM2 KLM4 = ',KLM2,KLM4
      ZERO = 0.0D0
C?    WRITE(6,*) ' Len = ', LEN
*.loop over symmetries
      IOFFP = 0 
      IOFFC = 0
      DO ISYM = 1, NSMOB
        IF(ISYM.EQ.1) THEN
          IOFFP = 1
          IOFFC = 1
         ELSE 
          IOFFP = IOFFP + NTOOBS(ISYM-1)*(NTOOBS(ISYM-1)+1)/2
          IOFFC = IOFFC + NTOOBS(ISYM-1) ** 2
         END IF
         LOB = NTOOBS(ISYM)
*
*. One-electron terms : Just copy H-block   
*
        IF(I12.NE.2) THEN
           CALL COPVEC(H(IOFFP),FOCK(IOFFP),LOB*(LOB+1)/2 )
        ELSE
           ZERO = 0.0D0
           CALL SETVEC(FOCK(IOFFP),ZERO,LOB*(LOB+1)/2 )
        END IF
*. Two-electron terms : 2*Coulomb - Exchange 
        IF(I12.NE.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLM2),ZERO,LOB**2)
        JOFFC = 0
        JOFFP = 0
        DO JSYM = 1, NSMOB
          IF(JSYM.EQ.1) THEN
            JOFFP = 1
            JOFFC = 1
          ELSE 
            JOFFP = JOFFP + NTOOBS(JSYM-1)*(NTOOBS(JSYM-1)+1)/2
            JOFFC = JOFFC + NTOOBS(JSYM-1) ** 2
          END IF
          JLEN = NTOOBS(JSYM)
*
*. Obtain 2*(ij!kl) - (il!kj)
*. =========================
*.  As the codes 
*. only can deliver (ij!kl)-(il!kj) and (ij!kl) this is done in 
*. two steps - very unefficient
*
*. (ij!kl) - (il!kj) 
          CALL GETINT_AO(WORK(KLM4),ISYM,ISYM,JSYM,JSYM,1)
*. (ij!kl)
          CALL GETINT_AO(WORK(KLM4B),ISYM,ISYM,JSYM,JSYM,0)
*
          NIJKL = JLEN*JLEN*LOB*LOB
          ONE = 1.0D0
          CALL VECSUM(WORK(KLM4),WORK(KLM4),WORK(KLM4B),
     &                ONE,ONE,NIJKL)
*
          CALL MATVCC(WORK(KLM4),P(JOFFC),WORK(KLM2P),
     &                LOB*LOB,JLEN*JLEN,0)
C              MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
          ONE = 1.0D0
          CALL VECSUM(WORK(KLM2),WORK(KLM2),WORK(KLM2P),
     &         ONE,ONE,LOB*LOB)
        END DO
*       ^ End of loop over symmetries JSYM
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Two electron terms in unpacked form'
          CALL WRTMAT(WORK(KLM2),LOB,LOB,LOB,LOB)
        END IF
*. The 2-electron terms came out in complete matrix form, 
*  pack and add to fock matrix
        XDUM = 1.0D0
        CALL TRIPAK(WORK(KLM2),WORK(KLM2P),1,LOB,LOB)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Two-electron terms in packed form'
          CALL PRSYM(WORK(KLM2P),LOB)
C         CALL PRSYM(FOCK(IOFFP),LOB)
        END IF
        HALF = 0.5D0
        CALL VECSUM(FOCK(IOFFP),FOCK(IOFFP),WORK(KLM2P),ONE,HALF,
     &              LOB*(LOB+1)/2 )
        END IF
*       ^ End if two-electron terms should be added
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Fock matrix  in packed form'
          CALL PRSYM(FOCK(IOFFP),LOB)
        END IF
      END DO
*       ^ End of loop over symmetries ISYM
*
* E = 1/2 Tr(H+F)P + ENUC
*
      E1 = 0.0D0
      EF = 0.0D0
      DO ISYM = 1, NSMOB
        LEN = NTOOBS(ISYM)
        IF(ISYM.EQ.1) THEN
          IOFFP = 1
          IOFFC = 1
        ELSE 
          IOFFP = IOFFP 
     &          + NTOOBS(ISYM-1)*(NTOOBS(ISYM-1)+1 )/2
          IOFFC = IOFFC 
     &          + NTOOBS(ISYM-1)**2
        END IF
*. One-electron contribution 
C TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
        IF(I12.NE.2) THEN
          CALL TRIPAK(WORK(KLM2),H(IOFFP),2,LEN,LEN)      
          E1 = E1 + INPROD(WORK(KLM2),P(IOFFC),LEN*LEN)
        END IF
*. Fock energy
        CALL TRIPAK(WORK(KLM2),FOCK(IOFFP),2,LEN,LEN)      
        EF = EF + INPROD(WORK(KLM2),P(IOFFC),LEN*LEN)
C?      WRITE(6,*) ' GET_FOCK : Block of F and P : '
C?      CALL WRTMAT(WORK(KLM2),LEN,LEN,LEN,LEN)
C?      CALL WRTMAT(P(IOFFC),LEN,LEN,LEN,LEN)
      END DO
      E = 0.5D0*(E1+EF)
      ENERGY = ENUC + E
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Fock matrix in AO basis '
        CALL APRBLM2(FOCK,NAOS_ENV,NAOS_ENV,NSMOB,1)      
      END IF
*
      IF(NTEST.GE.1) THEN
       WRITE(6,*) ' One- and two-electron contributions to energy ', 
     & E1, E - E1
C      WRITE(6,*) ' Two-electron contribution to energy ', E2
C      WRITE(6,*) ' Total energy ', E
       WRITE(6,*) ' Total energy ', ENERGY
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ', IDUMMY,'GET_FO')
      CALL LUCIAQEXIT('GT_FO')
      RETURN
      END
*
      SUBROUTINE GET_P_FROM_C(C,P,NOC,NAO,NSMOB)
*
* Obtain density matrix P in AO basis . P is obtained 
* in complete symmetry packed form 
*
* P(mu,nu) = 2*sum(i_occ)C(mu,i)C(nu,i)
*
* Jeppe Olsen, Sept. 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION C(*)
      DIMENSION NOC(NSMOB),NAO(NSMOB)
*. Output
      DIMENSION P(*)
*
      IOFFC = 0
      DO ISYM = 1, NSMOB
        IF(ISYM.EQ.1) THEN
          IOFFC = 1
        ELSE 
          IOFFC = IOFFC + NAO(ISYM-1)**2
        END IF
        LOB = NAO(ISYM) 
        ZERO = 0.0D0
C?      WRITE(6,*) ' Input C block : '
C?      CALL WRTMAT(C(IOFFC),LOB,LOB,LOB,LOB)
        CALL SETVEC(P(IOFFC),ZERO,LOB**2)
C?      WRITE(6,*) ' ISYM, NOC(ISYM) = ', ISYM, NOC(ISYM)
        DO I = 1, NOC(ISYM)
          DO NU = 1, LOB
            CNUI = C(IOFFC-1+(I-1)*LOB + NU)
            FACTOR = 2.0D0*CNUI
            DO MU = 1, LOB
              P(IOFFC-1+(NU-1)*LOB+MU) = P(IOFFC-1+(NU-1)*LOB+MU)
     &      + FACTOR*C(IOFFC-1+(I-1)*LOB + MU)
            END DO
          END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Density matrix in AO basis '
        WRITE(6,*) ' ========================== '
C APRBLM2(A,LROW,LCOL,NBLK,ISYM)
        CALL APRBLM2(P,NAO,NAO,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE GETINT_AO(XINT,ISYM,JSYM,KSYM,LSYM,IXCHNG)
*
* Obtain block of 2-e integrals with given symmetry 
* of the four blocks. The integrals are delivered 
* as an unpacked matrix XINT(I,J,K,L)
*
* IXCHNG = 0 : XINT(I,J,K,L) = (IJ!KL)
* IXCHNG = 1 : XINT(I,J,K,L) = (IJ!KL) - (IL!KJ)
*
* Jeppe Olsen, Sept. 2003
*              May 2012, aligned with new fetch of integrals
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      IKSM = 0
      JLSM = 0
      ICOUL = 1
      CFACX = 1.0D0
      IF(IXCHNG.EQ.1) THEN
        EFACX = 1.0D0
      ELSE
        EFACX = 0.0D0
      END IF
*
C     CALL GETINCN2(XINT,0,ISYM,0,JSYM,0,KSYM,0,LSYM,
C    &              IXCHNG,IKSM,JLSM,WORK(KINT2AO),
C    &              WORK(KPINT2),NSMOB,WORK(KINH1),ICOUL,0)
      CALL LUCIAGETINT(XINT,-1,ISYM,-1,JSYM,-1,KSYM,-1,LSYM,
     &              IXCHNG,IKSM,JLSM,ICOUL,CFACX,EFACX)
C
C          GETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                  IXCHNG,IKSM,JLSM,ICOUL,CFACX,EFACX)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GETINT_AO '
        WRITE(6,*) ' ISYM, JSYM, KSYM, LSYM = ',
     &               ISYM, JSYM, KSYM, LSYM
        IF(IXCHNG.EQ.0) THEN
          WRITE(6,*) ' Output is (IJ!KL)'
        ELSE IF (IXCHNG.EQ.1) THEN
          WRITE(6,*) ' Output is (IJ!KL)-(IL!KJ)'
        END IF
        NI = NAOS_ENV(ISYM)
        NJ = NAOS_ENV(JSYM)
        NK = NAOS_ENV(KSYM)
        NL = NAOS_ENV(LSYM)
        CALL WRTMAT(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_NOC
*
* Find number of doubly occupied orbitals per symmetry 
* for Hartree-Fock.
*  
* Routine name has been kept, but output and procedure has been changed
*
* Output is now stored in NHFD_GNSYM and NHFD_STASYM
*
* Jeppe Olsen, Sept. 2003, heavily modified May 23, 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc' 
      INCLUDE 'lucinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================='
        WRITE(6,*) ' Info from GET_NOC'
        WRITE(6,*) ' ================='
        WRITE(6,*)
      END IF
*
      IF(I_USE_SUPSYM.EQ.0) THEN
*. Just use input in NHFD_IRREP_SUPSYM, NHFS_IRREP_SUPSYM
        CALL ICOPVE(NHFD_IRREP_SUPSYM,NHFD_SUPSYM,NSMOB)
        CALL ICOPVE(NHFS_IRREP_SUPSYM,NHFS_SUPSYM,NSMOB)
*
        CALL ICOPVE(NHFD_SUPSYM,NHFD_GNSYM,NSMOB)
        CALL ICOPVE(NHFS_SUPSYM,NHFS_GNSYM,NSMOB)
*
        CALL SET_HF_DIST_STASYM
*. No reordering, so
C            ISTVC2(IVEC,IBASE,IFACT,NDIM)
        CALL ISTVC2(WORK(KMO_STA_TO_ACT_REO),0,1,NTOOB)       
        CALL ICOPVE(NHFD_SUPSYM,NHFD_GNSYM,NSMOB)
        CALL ICOPVE(NHFS_SUPSYM,NHFS_GNSYM,NSMOB)
*
        CALL ICOPVE(NHFD_SUPSYM,NHFD_STASYM,NSMOB)
        CALL ICOPVE(NHFS_SUPSYM,NHFS_STASYM,NSMOB)
      ELSE
*. expand from super symmetry irreps to components
        CALL N_SUPSYM_IRREP_TO_SUPSYM(NHFD_IRREP_SUPSYM,NHFD_SUPSYM)
        CALL N_SUPSYM_IRREP_TO_SUPSYM(NHFS_IRREP_SUPSYM,NHFS_SUPSYM)
*. and expand modify from supersymmetry to standard symmetry
        CALL N_SUPSYM_TO_STASYM(NHFD_SUPSYM,NHFD_STASYM)
        CALL N_SUPSYM_TO_STASYM(NHFS_SUPSYM,NHFS_STASYM)
*
        CALL ICOPVE(NHFD_SUPSYM,NHFD_GNSYM,NSMOB)
        CALL ICOPVE(NHFS_SUPSYM,NHFS_GNSYM,NSMOB)
*
        CALL SET_HF_DIST_SUPSYM
C       ORDER_SUPSYM_ORBITALS(NSPC,ISPC,MO_SUPSYM,IREO,
C    &           ISUPSYM_FOR_BAS)
        CALL ORDER_SUPSYM_ORBITALS(3,NHF_DSV_GNSYM,
     &       WORK(KMO_SUPSYM),
     &       WORK(KMO_STA_TO_ACT_REO),
     &       WORK(KISUPSYM_FOR_BAS)                )
      END IF
*
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Number of occupied orbitals per general symmetry '
        WRITE(6,*) ' ================================================='
        WRITE(6,*)
        CALL IWRTMA3(NHFD_GNSYM,1,NGENSMOB,1,NGENSMOB)
        WRITE(6,*) ' Number of occupied orbitals per standard symmetry'
        WRITE(6,*) ' ================================================='
        WRITE(6,*)
        CALL IWRTMA3(NHFD_STASYM,1,NSMOB,1,NSMOB)
      END IF
*
*
      RETURN
      END
      SUBROUTINE OPTIM_SCF(MAXIT_HF,
     &                    E_HF,CONVER_HF,E1_FINAL,NIT_HF,E_ITER)       
*
* Optimize SCF wavefunction
*
* Current version assumes closed shell states
*
* Modified to allow general symmetries, May 24, 2012
* Note: A major part of the work is done using the standard
* form of matrices. It is only for the diagonalization 
* that generalized symmetry is invoked.
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      LOGICAL LOGDIIS, CONVER_HF
*
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'gasstr.inc'
*. Energy in the various iterations
      DIMENSION E_ITER(*)
*
      NTEST= 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Info from OPTIM_SCF '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
      END IF
*
      LOGDIIS=.TRUE.
*. Max number of Fock matrices and densities
      MAX_NMAT = 5
      WRITE(6,*) ' Largest allowed number of saved F-matrices',
     &           MAX_NMAT
      IDUMMY = 0
      CALL LUCIAQENTER('OPTHF')
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ', IDUMMY,'OPTIM ')
*
      CALL MEMMAN(KDSNTO,NSMOB,'ADDL  ',1,'DSNTO  ')
      CALL MEMMAN(KDSENE,MAXIT_HF,'ADDL  ',2,'KDSENE ')
*
      CALL MEMMAN(KEPSIL,NTOOB,'ADDL  ',2,'EPSIL ')
      CALL MEMMAN(KEPSILT,NTOOB,'ADDL  ',2,'EPSILT')
*. And two one-electron matrices
      LENC = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,1)
*
      CALL MEMMAN(KLHSCR1,LENC,'ADDL  ',2,'HSCR1 ')
      CALL MEMMAN(KLHSCR2,LENC,'ADDL  ',2,'HSCR2 ')
*
      CALL MEMMAN(KLCBAR ,LENC,'ADDL  ',2,'CBAR  ')
      CALL MEMMAN(KLOCCBAR ,LENC,'ADDL  ',2,'OCCBAR  ')
*. Old density matrix 
        CALL MEMMAN(KLDAO_OLD,LENC,'ADDL  ',2,'DAO_OL')
*. Matrix for Improved virtual orbitals 
      CALL MEMMAN(KLFIVO,LENC,'ADDL  ',2,'F_IVO ')
*
      THRES_E = 1.0D-10
*
      I_ANALYZE_DENSI = 0
*. Allow subiterations with projected densities - outdated !!
      I_ALLOW_SUBIT = 0
      IF(IHFSOLVE.EQ.1) THEN
*. Standard Roothaan-Hall
        I_DO_SUBOPT = 0
      ELSE
*. eopd optimization
        I_DO_SUBOPT = 1
        IF(I_DO_SUBOPT.EQ.0) THEN
          WRITE(6,*) ' No EOPD optimization '
        ELSE 
          WRITE(6,*) ' EOPD optimization '
        END IF
      END IF
*. Max number of subiterations per Fock evaluation
      MAX_SUBIT = 0
*. Min Ratio between total and projected density for 
*. allowing subspace iterations 
      XPT = 0.50
      IF(I_DO_SUBOPT.EQ.1) THEN
*
*. Prepare for analysis of density matrices
*
*. Allocate space for some extra matrices
*. SAO in expanded form 
        CALL MEMMAN(KLSAO_E,LENC,'ADDL  ',2,'SAO_E ')
*. HAO in expanded form
        CALL MEMMAN(KLHAO_E,LENC,'ADDL  ',2,'HAO_E ')
*. FAO in symmetry packed form
        CALL MEMMAN(KLFAO_E,LENC,'ADDL  ',2,'SAO_E ')
*. And extra test matrices
        CALL MEMMAN(KLCTEST,LENC,'ADDL  ',2,'CTEST ')
        CALL MEMMAN(KLCTEST2,LENC,'ADDL  ',2,'CTEST2')
*. A few more scratch matrices 
        LEN = MAX(2*LENC,MAXIT_HF**2)
        WRITE(6,*) ' LEN = ', LEN
        CALL MEMMAN(KLSCR1,LEN,'ADDL  ',2,'SCR1  ')
        CALL MEMMAN(KLSCR2,LEN,'ADDL  ',2,'SCR2  ')
        CALL MEMMAN(KLSCR3,LEN,'ADDL  ',2,'SCR3  ')
        CALL MEMMAN(KLSCR4,LEN,'ADDL  ',2,'SCR4  ')
*. Inverse overlap of densities
        CALL MEMMAN(KLSINV,MAXIT_HF**2,'ADDL  ',2,'SINV  ')
*. Projected density
        CALL MEMMAN(KLPROJD,LENC,'ADDL  ',2,'PROJD ')
*. And two vectors of length MAXIT_HF
        CALL MEMMAN(KLPROJC,MAXIT_HF,'ADDL  ',2,'PROJC ')
        CALL MEMMAN(KLPROJC_OLD,MAXIT_HF,'ADDL  ',2,'PROJC ')
        CALL MEMMAN(KLSCRVEC,MAXIT_HF,'ADDL  ',2,'SCRVEC')
        CALL MEMMAN(KLSCRVEC2,MAXIT_HF,'ADDL  ',2,'SCRVEC')
*. Change in density matrices 
        CALL MEMMAN(KLDDEL,LENC,'ADDL  ',2,'DDEL  ')
*. And an extra density to play around with
        CALL MEMMAN(KLDAO_X,LENC,'ADDL  ',2,'DAO_X ')
*. Overlap matrix between old and new set of orbitals
        CALL MEMMAN(KLCOSCN,LENC,'ADDL  ',2,'COSCN ')
*. Weight of old orbitals in new 
        CALL MEMMAN(KLW,NTOOB,'ADDL  ',2,'W_ON  ')
*. Obtain SAO and HAO matrices in expanded form '
        CALL TRIPAK_BLKM(WORK(KLSAO_E),WORK(KSAO),2,NTOOBS,NSMOB)
        CALL TRIPAK_BLKM(WORK(KLHAO_E),WORK(KHAO),2,NTOOBS,NSMOB)
        WRITE(6,*) ' Unpacked AO overlap matrix '
        CALL APRBLM2(WORK(KLSAO_E),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. And still a few
        CALL MEMMAN(KLHX,LENC,'ADDL  ',2,'HX    ')
        CALL MEMMAN(KLSX,LENC,'ADDL  ',2,'SX    ')
        CALL MEMMAN(KLFX,LENC,'ADDL  ',2,'FX    ')
        CALL MEMMAN(KLCX,LENC,'ADDL  ',2,'CX    ')
* 
*
*. On input an initial matrix is assumed residing in work(kccur)
*
      NFOCK = 0
      WRITE(6,*) ' Largest allowed number of iterations', MAXIT_HF
      DO ITER = 1, MAXIT_HF + 1
        WRITE(6,*)
        WRITE(6,*) ' =========================================='
        WRITE(6,*) ' Information from Fock iteration ', ITER
        WRITE(6,*) ' =========================================='
        WRITE(6,*)
        NIT_HF = ITER
*. Copy current C expansion to old
        CALL COPVEC(WORK(KCCUR),WORK(KCOLD),LENC)
*. Construct density matrix for given mo-ao expansion
        IF(ITER.GT.1) CALL COPVEC(WORK(KDAO),WORK(KLDAO_OLD),LENC)
        CALL GET_P_FROM_C(WORK(KCCUR),WORK(KDAO),NHFD_STASYM,
     &       NTOOBS,NSMOB)
*. For comparison with Lea : Test various energy approximations for current 
*. density 
C            GET_APPROX_HFEM(E_APPROX,IAPPROX,DENSI,NDENSI)
        I_DO_LEA_TESTS = 0
        IF(I_DO_LEA_TESTS.EQ.1) THEN
          CALL GET_APPROX_HFEM(EA0, 0, WORK(KDAO),NFOCK) 
          WRITE(6,*) ' Energy, approximation 0 = ', EA0
          CALL GET_APPROX_HFEM(EA1, 1, WORK(KDAO),NFOCK) 
          WRITE(6,*) ' Energy, approximation 1 = ', EA1
          CALL GET_APPROX_HFEM(EA2, 2, WORK(KDAO),NFOCK) 
          WRITE(6,*) ' Energy, approximation 2 = ', EA2
          CALL GET_APPROX_HFEM(EA3, 3, WORK(KDAO),NFOCK) 
          WRITE(6,*) ' Energy, approximation 3 = ', EA3
        END IF
*. Obtain Fock matrix for current density 
        WRITE(6,*) 
     &  ' New Fock matrix will be generated '
        CALL GET_FOCK(WORK(KFAO),WORK(KHAO),WORK(KDAO),ENERGY,
     &                ECORE_ORIG,3)
        WRITE(6,'(A,I3, E25.15)') ' Iteration, current energy = ', 
     &                          ITER, ENERGY
        E_HF = ENERGY
*
*. Test for convergence 
*
        E_ITER(ITER)=ENERGY
        IF(ITER.GT.1) THEN
          DELTAE=DABS(E_ITER(ITER)-E_ITER(ITER-1))
        END IF
*
        IF(ITER.GT.1.AND.DELTAE.LT.THRES_E)THEN 
          WRITE(6,*) ' Energy converged, last energy change ',
     &    DELTAE
          CONVER_HF = .TRUE.
          GOTO 1111
        END IF
        IF(ITER.EQ.MAXIT_HF + 1) GOTO 1111
*
*. Save current density 
        NFOCK = NFOCK + 1
*. Eliminate, if required, the older matrices 
        IF(NFOCK.GT.MAX_NMAT.AND.I_DO_SUBOPT.EQ.1) THEN
*. eliminate oldest matrix - written first
          DO IMAT = 2, MAX_NMAT
*
            ID_FROM = KDAO_COLLECT+(IMAT-1)*LENC
            ID_TO   = KDAO_COLLECT+(IMAT-1-1)*LENC
            CALL COPVEC(WORK(ID_FROM),WORK(ID_TO),LENC)
*
            IF_FROM = KFAO_COLLECT + (IMAT-1)*LENP
            IF_TO   = KFAO_COLLECT + (IMAT-1-1)*LENP
            CALL COPVEC(WORK(IF_FROM),WORK(IF_TO),LENP)
          END DO
*. Number of Fock matrices ( when new is added ..)
          NFOCK = MAX_NMAT
          CALL MEMCHK2('CHCKA ') 
        END IF
*. and new matrices to our collection 
        IF(I_DO_SUBOPT.EQ.1) THEN
          CALL COPVEC(WORK(KDAO),WORK(KDAO_COLLECT+(NFOCK-1)*LENC),
     &               LENC)
          CALL COPVEC(WORK(KFAO),WORK(KFAO_COLLECT+(NFOCK-1)*LENP),
     &               LENP)
          ZERO = 0.0D0
          CALL SETVEC(WORK(KLPROJC_OLD),ZERO,NFOCK)
        END IF
*
        IF(NFOCK.GT.1.AND.I_DO_SUBOPT.EQ.1) THEN
*. Find optimal combination of available densities 
          IPROJ_DOPT = 1
          CALL OPTIM_EHFAPR(WORK(KDAO),NFOCK,NFOCK,WORK(KLPROJC),
     &    IPROJ_DOPT)
*. When IPROJ_DOPT = 1, then the returned is the optimal density (Dtilde)
*. projected in the space, and WORK(KLPROJC) contains the coefs of the 
*. optimized density
*. Obtain orbitals corresponding to DENSITY DOPT, and the 
*. corresponding orbitals 
C              DIAG_SDS(D,S,C,XOCCNUM)
          CALL DIAG_SDS(WORK(KDAO),WORK(KSAO),WORK(KLCBAR),
     &                 WORK(KLOCCBAR))
*. Analyze overlap between orbitals coming out from OPTIM_EHFAPR
*. and the input orbitals this routine
          WRITE(6,*) ' Info on overlap of Cbar and Cold orbitals'
          CALL OVERLAP_COLD_CNEW(WORK(KCCUR),WORK(KLCBAR),
     &         WORK(KSAO),WORK(KLCOSCN),WORK(KLSCR1),
     &         NBAS_GENSMOB,NOCOBS,NSMOB)
*. Find weight of old occupied MO's in new MO's
          CALL WEIGHT_NEW_OCCOLD(WORK(KLCOSCN),NOCOBS,WORK(KLW),WMIN)
*. And make the orbitals from DBAR the new ones
          CALL COPVEC(WORK(KLCBAR),WORK(KCCUR),LENC)
*
          IF(IPROJ_DOPT.EQ.1) THEN
*. Find Fock matrix corresponding to specified sum of densities 
            CALL GET_F_FOR_SUM_DENSI(WORK(KFAO),WORK(KLPROJC),NFOCK)
          ELSE
             WRITE(6,*) ' You asked for NEW Fock matrix for D_TILDE '
             CALL GET_FOCK(WORK(KFAO),WORK(KHAO),WORK(KDAO),ENERGY,
     &                 ECORE_ORIG,3)
          END IF
        END IF
*       ^ End of diagonalization should be done with control
*
*. Diagonalize to obtain new improved MO's
*
C?      WRITE(6,*) ' Info from diagonalization of F '
*. Will we shift to ensure that the overlap between new and 
* old orbitals are larger than some threshold
        I_DO_SHIFT_IN_SCFDIA = 0
        IF( I_DO_SHIFT_IN_SCFDIA.EQ.0) THEN
*. Reform from standard to general symmetry 
C             REFORM_MAT_STA_GEN(ASTA,AGEN,IPACK,IWAY)
         CALL REFORM_MAT_STA_GEN(WORK(KFAO),WORK(KLFX),1,1)
         CALL REFORM_MAT_STA_GEN(WORK(KSAO),WORK(KLSX),1,1)
          CALL SCF_DIA(WORK(KLFX),WORK(KLSX),
     &                 WORK(KLCX),WORK(KEPSIL)) 
*. The C-coefficient are given in the general symmetry basis,
*  transfer to standard
C     REFORM_CMO_STA_GEN(CMO_STA,CMO_GEN,
C    &           IDO_REORDER,IREO,IWAY)
         CALL REFORM_CMO_STA_GEN(WORK(KCCUR),WORK(KLCX),0,IDUM,2)
         CALL REO_CMOAO(WORK(KCCUR),WORK(KLCX),
     &        WORK(KMO_STA_TO_ACT_REO),1,1)
        ELSE 
          WOCOCMIN = 0.90
          ALPHA_INI = 0.0D0
C              SCF_DIA_CTL(F,S,C,E,D,COLD,WOCOCMIN,ALPHA_INI,NOCPSM)
          CALL SCF_DIA_CTL(WORK(KFAO),WORK(KSAO),WORK(KCCUR),
     &         WORK(KEPSIL),WORK(KDAO),WORK(KCOLD),
     &           WOCOCMIN,ALPHA_INI,NOCOBS)
*. Scan the solution of the shifted Fock equations for various shifts
           I_SCAN = 1
           IF(I_SCAN.EQ.1) THEN
C                 SCF_DIA_SHIFT_SCAN (F,S,C,E,D,COLD,NOCPSM)
*
             I_DO_IVO = 1
             IF(I_DO_IVO.EQ.1) THEN
*. construct Fock matrix for improved virtual orbitals
               NEL = NELFTP(1) + NELFTP(2)
               CALL GET_FIVO(WORK(KFAO),WORK(KHAO),WORK(KLFIVO),
     &               NOCOBS,NBAS_GENSMOB,NSMOB,NEL,WORK(KCCUR),
     &               WORK(KLSAO_E),
     /               WORK(KLSCR1),WORK(KLSCR2),WORK(KLSCR3) )
               CALL SCF_DIA_SHIFT_SCAN(WORK(KLFIVO),WORK(KSAO),
     &               WORK(KLCTEST),WORK(KCCUR),
     &               WORK(KEPSILT),WORK(KDAO),WORK(KCOLD),NOCOBS,
     &               NFOCK)
             ELSE 
               CALL SCF_DIA_SHIFT_SCAN(WORK(KFAO),WORK(KSAO),
     &               WORK(KLCTEST),WORK(KCCUR),
     &               WORK(KEPSILT),WORK(KDAO),WORK(KCOLD),NOCOBS,
     &               NFOCK)
             END IF
           END IF
        END IF
*
      END DO
1111  CONTINUE 
*
C      WRITE(6,*) ' Number of Fock matrix evaluations ', NFOCK
      I_ANA_STABILITY = 0
      IF(I_ANA_STABILITY.EQ.1) THEN
C     ANA_STABILITY_FOR_HF(C)
*. Note: Jeppe. at the moment the program is doing MO-MO 
* transformation, gives problem in DALTON environment
        XDUM = 0
        NEW_MO_OBTAINED = -1
        CALL ANA_STABILITY_FOR_HF(WORK(KCCUR),1,0,INEG_TOTSYM,XDUM,
     &                            NEW_MO_OBTAINED)
C            ANA_STABILITY_FOR_HF(C,I_WILL_DO_TRANS,I_FOLLOW_MODE,
C    &                                INEG_TOTSYM,CNEW,NEW_MO_OBTAINED)
      END IF
*
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) ' Final energy in au  ', ENERGY

C     IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Final list of MO-AO expansion coefficients '
        CALL APRBLM2(WORK(KCCUR),NTOOBS,NTOOBS,NSMOB,0)
*
        WRITE(6,*) ' Final list of MO eigenvalues over gen sym'
        CALL WRITE_ORBENERGIES(WORK(KEPSIL),NBAS_GENSMOB,NGENSMOB)
C       CALL WRTMAT(WORK(KEPSIL),1,NTOOB,1,NTOOB)
C     END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ', IDUMMY,'OPTIM ')
      CALL LUCIAQEXIT('OPTHF')

      RETURN
      END
      SUBROUTINE WRITE_ORBENERGIES(EPSIL,NAOS,NSMOB)
*
* Print list of orbital energies 
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION EPSIL(*)
      INTEGER NAOS(NSMOB)
*
      WRITE(6,*)
      WRITE(6,*) ' ================='
      WRITE(6,*) ' Orbital energies '
      WRITE(6,*) ' ================='
      WRITE(6,*)
*
      IOFF = 1
      DO ISYM = 1, NSMOB
        IF(ISYM.GT.1) THEN
          IOFF = IOFF + NAOS(ISYM-1)
        END IF
        WRITE(6,*) ' Symmetry ', ISYM
        LEN = NAOS(ISYM)
        CALL WRTMAT(EPSIL(IOFF),1,LEN,1,LEN)
      END DO
*
      RETURN
      END
      SUBROUTINE OVERLAP_COLD_CNEW(COLD,CNEW,S,COSCN,SCR,
     &          NAOS,NOC,NSMOB)
*
* Obtain overlap matrix between old and new mo-coefficients
*
* COSCN = COLD(T) S CNEW
*
*. Jeppe Olsen, Sept. 2003
*
      INCLUDE 'implicit.inc'
*. General input
      INTEGER NAOS(NSMOB),NOC(NSMOB)
*. Specific input
      DIMENSION COLD(*),CNEW(*),S(*)
*. Output
      DIMENSION COSCN(*)
*. Scratch : Should be able to contain complete total symmetric matrix
      DIMENSION SCR(*)
*. Expand S matrix to complete form
C TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
      CALL TRIPAK_BLKM(COSCN,S,2,NAOS,NSMOB)
*. Obtain  COLD (T) S in SCR
C     MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
      CALL MULT_BLOC_MAT(SCR,COLD,COSCN,NSMOB,
     &     NAOS,NAOS,NAOS,NAOS,NAOS,NAOS,1)
*.Obtain COLD(T) S CNEW 
      CALL MULT_BLOC_MAT(COSCN,SCR,CNEW,NSMOB,
     &     NAOS,NAOS,NAOS,NAOS,NAOS,NAOS,0)
*. Obtain the part of the new orbitals, that are 
*. spanned by the old occupied orbitals
C      IOFF = 1
C      IOFFM = 1
C      DO ISYM = 1, NSMOB
C        DO IORB = 1, NAOS(ISYM)
C          X = 0.0D0
C          DO JORB = 1, NOC(ISYM)
C            X = X + COSCN(IOFFM-1+(
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' all-occupied part of COLD(T) S CNEW '
        IOFF = 1
        DO ISYM = 1, NSMOB
          WRITE(6,*) ' Symmetry ', ISYM
          LOC = NOC(ISYM)
          LOB = NAOS(ISYM)
          CALL WRTMAT(COSCN(IOFF),LOB,LOC,LOB,LOB)
          IOFF = IOFF + LOB*LOB
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE GET_HOMO_LUMO(EPSIL,NAOS,NOCC,NSMOB,EHOMO,ELUMO,
     &                         EHL_GAP)
*
* Obtain energy of HOMO and LUMO from a list of orbital energies
* THe occupied are assumed to be the lowest in each symmetry
*
* EHLGAP is the smallest difference between occupied and 
* unoccupied orbitals for orbitals with a given symmetry
*
*. Jeppe Olsen, Oct. 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION EPSIL(*)
      INTEGER NAOS(NSMOB),NOCC(NSMOB)
*. Initial guesses
      EHOMO = EPSIL(1)
      ELUMO = EPSIL(1+NOCC(1))
*
      IOFF = 1
      EHL_GAP = 100000000.0
      DO ISM = 1, NSMOB
        LOC = NOCC(ISM) 
        LOB = NAOS(ISM)
        DO IOB = 1, LOC
          EHOMO = MAX(EHOMO,EPSIL(IOFF-1+IOB))
        END DO
        DO IVIRT = LOC+1,LOB
          ELUMO = MIN(ELUMO,EPSIL(IOFF-1+IVIRT))
        END DO
C?      WRITE(6,*) ' ,EPSIL(IOFF-1+LOC+1)-EPSIL(IOFF-1+LOC)',
C?   &                EPSIL(IOFF-1+LOC+1)-EPSIL(IOFF-1+LOC)
        IF(LOC.NE.0) 
     &  EHL_GAP = MIN(EHL_GAP,EPSIL(IOFF-1+LOC+1)-EPSIL(IOFF-1+LOC))
        IOFF = IOFF + LOB
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' EHOMO = ', EHOMO
        WRITE(6,*) ' ELUMO = ', ELUMO
        WRITE(6,*) ' EHL_GAP = ', EHL_GAP
      END IF
*
      RETURN
      END
      SUBROUTINE WEIGHT_NEW_OCCOLD(COSCN,NOCPSM,W_NEW_OCCOLD,WOCC_MIN)
*
*. Obtain weight of old occupied orbitals in new orbitals
*. The information about the overlap between old and new orbitals 
*. is given in COSCN
*. The smallest weight of the new occupied orbitals is returned 
* in WOCC_MIN
*
* Jeppe Olsen, Oct. 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION COSCN(*), NOCPSM(*)
*. Output
      DIMENSION W_NEW_OCCOLD(*)
*
      IOFF = 1
      IORB_ABS = 0
      DO ISYM = 1, NSMOB
        NI_ORB = NAOS_ENV(ISYM)
        NI_OCC = NOCPSM(ISYM)
        DO I = 1, NI_ORB
          X = 0.0D0
          DO IOCC = 1, NI_OCC
            X = X + COSCN(IOFF-1+(I-1)*NI_ORB+IOCC)**2
          END DO
          IORB_ABS = IORB_ABS + 1
          W_NEW_OCCOLD(IORB_ABS) = X
        END DO
        IOFF = IOFF + NI_ORB**2
      END DO
*
      WOCC_MIN = 1.0D0
      IOFF = 1
      DO ISYM = 1, NSMOB
        DO IOCC = 1, NOCPSM(ISYM)
          WOCC_MIN = MIN(WOCC_MIN,W_NEW_OCCOLD(IOFF-1+IOCC))
        END DO
        IOFF = IOFF + NAOS_ENV(ISYM)
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The weight of old occupied orbitals in new MO''s '
        CALL PRINT_ORBVEC(W_NEW_OCCOLD,NSMOB,NAOS_ENV)
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Smallest overlap between new occupied orbitals '
        WRITE(6,*) ' and space of previous occupied orbitals is ',
     &             WOCC_MIN
      END IF
*
      RETURN
      END 
      SUBROUTINE PRINT_ORBVEC(XVEC,NSMOB,NOBPSM)
*
* Print vector over orbitals, separating the various symmetries
*
*. Jeppe Olsen, Oct. 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION XVEC(*)
      INTEGER NOBPSM(NSMOB)
*
      IOFF = 1
      DO ISM = 1, NSMOB
        WRITE(6,*) ' Symmetry .... ', ISM
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        LOB = NOBPSM(ISM)
        CALL WRTMAT(XVEC(IOFF),1,LOB,1,LOB)
        IOFF = IOFF + LOB
      END DO
*
      RETURN
      END 
      SUBROUTINE PROJECT_DENSI_ON_DENSI(DENSIIN,DENSI,LDENSI,NDENSI,
     &           SAO,SINV,PROJ_DENSI,PROJ_COEF,NSMOB,NOBPSM,
     &           SCR1,SCR2,SCRVEC,I_DO_SINV,DORIG_NORM, DPROJ_NORM)
*
* An input density DENSIIN is given. Project this density on 
* a set of NDENSI known densities in  DENSI.
* Obtain the projected density, PROJ_DENSI, and 
* the expansion coefficients of the projected density
*
* If I_DO_SINV .ne. 0 then the inverted overlap matrix of the 
* input densities are generated
*
*. Densities and AO matrices are assumed in complete 
*. (i.e. not lower half) symmetry packed form 
* Jeppe Olsen, Gwatt, Schweiz, Sept. 2003
*
      INCLUDE 'implicit.inc'
      REAL*8  INNER_PRODUCT_MAT, INPROD
*
*. Input
*
*. Density to be projected 
      DIMENSION DENSIIN(LDENSI)
*. Densities on which we will project 
      DIMENSION DENSI(LDENSI,NDENSI)
*. Overlap matrix in AO basis in complete packed form
      DIMENSION SAO(*)
*. Inverse overlap matrix of densities on which we will project 
*. (is constructed if I_DO_SINV .eq. 1)
      DIMENSION SINV(NDENSI,NDENSI)
*. Basis set specification
      INTEGER NOBPSM(NSMOB)
*
*. Output
*
*. The projected density
      DIMENSION PROJ_DENSI(LDENSI)
*. The expansion coefficients of the projected density
      DIMENSION PROJ_COEF(NDENSI)
*
*. Scratch : 2 matrices of length max(LDENSI,NDENSI**2)
*            and one vector of length NDENSI
*
      DIMENSION SCR1(*),SCR2(*)
    
      DIMENSION SCRVEC(NDENSI)
*
      NTEST = 00
*. Construct inverse overlap matrix
      IF(I_DO_SINV.EQ.1) THEN
*. First the overlap matrix 
C     GET_OVERLAPMAT_AO_DENS(DENSI,LDENSI,NDENSI,SAO,
C    &           SOVERLAP,SCR1,SCR2,NSMOB,NOBPSM)
        CALL GET_OVERLAPMAT_AO_DENS(DENSI,LDENSI,NDENSI,SAO,
     &       SINV,SCR1,SCR2,NSMOB,NOBPSM)
C            INVERT_BY_DIAG2(A,B,SCR,VEC,NDIM)
        CALL INVERT_BY_DIAG2(SINV,SCR1,SCR2,SCRVEC,NDENSI)
      END IF
*
*. Obtain overlap <DENSIIN!DENSI(I)>
      DO I = 1, NDENSI
       SCRVEC(I) = INNER_PRODUCT_MAT(DENSIIN,DENSI(1,I),SAO,SCR1,SCR2,
     &             NSMOB,NOBPSM)
C     INNER_PRODUCT_MAT(A,B,S,SCR1,SCR2,NBLK,LBLK)
      END DO
*. Obtain expansion coefficients as SINV(I,J) * <DENSI(J)!DENSIIN>
      CALL MATVCB(SINV,SCRVEC,PROJ_COEF,NDENSI,NDENSI,0)
*. norm of projected matrix 
      DPROJ_NORM = INPROD(SCRVEC,PROJ_COEF,NDENSI)
      DPROJ_NORM = SQRT(DPROJ_NORM)
C?    IF(NTEST.GE.5) WRITE(6,*) ' Norm of projected part of density ', 
C?   &           SQRT(DPROJ_NORM)
*. norm of Original matrix 
      DORIG_NORM = INNER_PRODUCT_MAT(DENSIIN,DENSIIN,SAO,SCR1,SCR2,
     &             NSMOB,NOBPSM)
      DORIG_NORM = SQRT(DORIG_NORM)
      IF(NTEST.GE.5) 
     &   WRITE(6,*) ' Norm of input and projected density ', 
     &   DORIG_NORM, DPROJ_NORM
*
C     INNER_PRODUCT_MAT(A,B,S,SCR1,SCR2,NBLK,LBLK)
*  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
*. And the projected density 
C  MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(DENSI,PROJ_COEF,PROJ_DENSI,LDENSI,NDENSI,0)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info on projected density : '
        WRITE(6,*) ' Matrix to be projected '
        CALL APRBLM2(DENSIIN,NOBPSM,NOBPSM,NSMOB,0)      
        WRITE(6,*) ' Expansion coefficients of projected density'
        CALL WRTMAT(PROJ_COEF,1,NDENSI,1,NDENSI)
        WRITE(6,*) ' Projected density '
        CALL APRBLM2(PROJ_DENSI,NOBPSM,NOBPSM,NSMOB,0)      
      END IF
*
      RETURN
      END
*
      SUBROUTINE GET_OVERLAPMAT_AO_DENS(DENSI,LDENSI,NDENSI,SAO,
     &           SOVERLAP,SCR1,SCR2,NSMOB,NOBPSM)
*
* Obtain overlap matrix of NDENSI AO densities DENSI
*
* Overlap of two densities is defineed as Tr (D1 S) (D2 S)
*
* Densities and overlap matrix are supposed to be in 
* complete symmetry packed form
*
* Jeppe Olsen, Gwatt, Schweiz, Sept. 3003
*
      INCLUDE 'implicit.inc'
      REAL * 8 INNER_PRODUCT_MAT
*. Input
      INTEGER NOBPSM(NSMOB)
      DIMENSION DENSI(LDENSI,NDENSI), SAO(*)
*. Output
      DIMENSION SOVERLAP(NDENSI,NDENSI)
*. Scratch through argument list , each must be of length atleast LDENSI
      DIMENSION SCR1(*),SCR2(*)
*
C?    WRITE(6,*) ' LDENSI = ', LDENSI
      DO I = 1, NDENSI
       DO J = 1, I
*
C?       WRITE(6,*) ' I,J = ', I,J
C?       WRITE(6,*) ' DENSI(I) '
C?       CALL APRBLM2(DENSI(1,I),NOBPSM,NOBPSM,NSMOB,0)      
C?       WRITE(6,*) ' DENSI(J) '
C?       CALL APRBLM2(DENSI(1,J),NOBPSM,NOBPSM,NSMOB,0)      
C?       WRITE(6,*) ' Overlap matrix '
C?       CALL APRBLM2(SAO,NOBPSM,NOBPSM,NSMOB,0)      
*
         SOVERLAP(I,J) = INNER_PRODUCT_MAT(DENSI(1,I),DENSI(1,J),SAO,
     &            SCR1, SCR2,NSMOB,NOBPSM)
         SOVERLAP(J,I) = SOVERLAP(I,J)
       END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap matrix of densities '
        CALL WRTMAT(SOVERLAP,NDENSI,NDENSI,NDENSI,NDENSI)
      END IF
*
      RETURN
      END
* 
      REAL*8 FUNCTION INNER_PRODUCT_MAT(A,B,S,SCR1,SCR2,NBLK,LBLK)
*
* Obtain general inner product of two blocked matrices, 
* matrices are in complete block form 
*
* <A|B> = Tr (ASBS)
*
* This form corresponds to consider matrices as vectors, and 
* using a standard (nonorthogonal) vector norm.
*
* Jeppe Olsen, Sept. 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*),B(*),S(*)
      INTEGER LBLK(NBLK)
*. Scratch : Must each hold complete matrix
      DIMENSION SCR1(*),SCR2(*)
C       MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)

*. AS in SCR1
      CALL MULT_BLOC_MAT(SCR1,A,S,NBLK,LBLK,LBLK,LBLK,LBLK,
     &                   LBLK,LBLK,0)
*. ASB in SCR2
      CALL MULT_BLOC_MAT(SCR2,SCR1,B,NBLK,LBLK,LBLK,LBLK,LBLK,
     &                   LBLK,LBLK,0)
*. ASBS in SCR1
      CALL MULT_BLOC_MAT(SCR1,SCR2,S,NBLK,LBLK,LBLK,LBLK,LBLK,
     &                   LBLK,LBLK,0)
*. And then the trace
      TRACE = TRACE_BLK_MAT(SCR1,NBLK,LBLK,0)
*
      INNER_PRODUCT_MAT = TRACE 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inner product of two matrices = ' ,
     &               INNER_PRODUCT_MAT
      END IF
*
      RETURN
      END
      FUNCTION TRACE_BLK_MAT(A,NBLK,LBLK,IPAK)
*
* Obtain trace of a blocked matrix 
*
* Jeppe Olsen, Sept. 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*) 
      INTEGER LBLK(NBLK)
*
      TRACE = 0.0D0
      IOFF = 1
      DO IBLK = 1, NBLK
        TRACE = TRACE + TRACE_MAT(A(IOFF),LBLK(IBLK),IPAK)
        IF(IPAK.EQ.1) THEN
          IOFF  = IOFF + LBLK(IBLK)*(LBLK(IBLK)+1)/2
        ELSE 
          IOFF  = IOFF + LBLK(IBLK)*LBLK(IBLK)
        END IF
      END DO
*
      TRACE_BLK_MAT = TRACE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Trace of blocked matrix = ', TRACE
      END IF
*
      RETURN
      END
      FUNCTION TRACE_MAT(A,NDIM,IPAK)
*
* Obtain trace of matrix 
*
* IPAK = 1 => matrix is packed in lower triangular form
*
*. Jeppe Olsen, Sept. 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*)
*
      TRACE = 0.0D0
*
      IF(IPAK.EQ.1) THEN
        DO I = 1, NDIM
          II = I*(I+1)/2
          TRACE = TRACE + A(II)
        END DO
      ELSE 
        DO I = 1, NDIM
          II = (I-1)*NDIM + I
          TRACE = TRACE + A(II)
        END DO
      END IF
*
      TRACE_MAT = TRACE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Trace of matrix = ', TRACE
      END IF
*
      RETURN
      END
      SUBROUTINE OBTAIN_DBAR_COEF(COEF,DBAR_COEF,IREFDENS,NDENSI,INORM)
*
*. Old way :
* ===========
* IREFDENS .EQ.0 : DBAR_COEF(I) = COEF(I)
*
* IREFDENS .NE.0 : DBAR_COEF(I) = COEF(I) for I .lt. IREFDENS
*                  DBAR_COEF(IREFDENS) = 1
*                  DBAR_COEF(I) = COEF(I-1) for I .ge. IREFDENS
*
* If INORM = 1, then the coefficients are normalized in the 
*               standard DIIS sense : sum_i DBAR_COEF(I) = 1
*
* New way : 
* =========
* 
* D_BAR = D_IREFDENS + Sum(I.ne.IREFDENS) COEF_I(D_I-D_IREFDENS)
*
* The output, DBAR_COEF, are thus in general the expansion coefficients of 
* the individual densities

*
*. Jeppe Olsen, Sept. 2003
*               New form added Sept. 2004
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION COEF(*)
*. Output
      DIMENSION DBAR_COEF(*)
*
      INEW_OR_OLD = 1
      IF(INEW_OR_OLD.EQ.2) THEN
        IF(IREFDENS.EQ.0) THEN
          CALL COPVEC(COEF,DBAR_COEF,NDENSI)
        ELSE 
          DO IDENS = 1, IREFDENS-1
            DBAR_COEF(IDENS) = COEF(IDENS)
          END DO
          DBAR_COEF(IREFDENS) = 1
          DO IDENS = IREFDENS+1,NDENSI
            DBAR_COEF(IDENS) = COEF(IDENS-1)
          END DO
        END IF
*
        IF(INORM.EQ.1) THEN
           XNORM = ELSUM(DBAR_COEF,NDENSI)
           FACTOR = 1.0D0/XNORM
           CALL SCALVE(DBAR_COEF,FACTOR,NDENSI)
        END IF
      ELSE 
*. New way 
        CSUM = ELSUM(COEF,NDENSI-1)
        DO IDENS = 1, IREFDENS-1
          DBAR_COEF(IDENS) = COEF(IDENS)
        END DO
        DBAR_COEF(IREFDENS) = 1.0D0 - CSUM
        DO IDENS = IREFDENS+1,NDENSI
          DBAR_COEF(IDENS) = COEF(IDENS-1)
        END DO
      END IF
*.    ^ End of switch between new and old ways
*
      NTEST = 00 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input coefficients '
        CALL WRTMAT(COEF,1,NDENSI-1,1,NDENSI-1)
        WRITE(6,*) ' Coefficients for DBAR : '
        CALL WRTMAT(DBAR_COEF,1,NDENSI,1,NDENSI)
      END IF
*
      RETURN
      END
      SUBROUTINE OBTAIN_D_TILDE(D_TILDE,IREFDENS,DCOEF,NDENSI,
     &                          LDENSI,DENSI,
     &                          INORMA,NEL,S,NSMOB,NAOS,SCR1,SCR2,SCR3,
     &                          SCR4)
*
* Obtain D_TILDE densi by purifying, and normalizing 
* a D_BAR density ... 
*
* D_BAR density is a sum of densities 
*     DBAR = D(IREFDENS) + sum(i) DCOEF(i) DENSI(I)  (Old form)
*     DBAR = D(IREFDENS) 
*          + sum(i.ne.IRFEDENS) DCOEF(I)(DENSI(I)-DENSI(IREFDENS))
*. and calling the density purifier as P, D_TILDE is obtained as 
*
* INORMA = 0 D_TILDE = P(D_BAR)
* 
* INORMA = 1 D_TILDE = P(D_BAR) * NEL/TR(P_D_BAR S )
*
* The densities are assumed in symmetry packed complete form, 
* and LDENSI refers to the length of this array
*
*
* The density is supposed to be a orbital density 
*. Jeppe Olsen, Sept. 2003
*
      INCLUDE 'implicit.inc'
*. General input
      INTEGER NAOS(NSMOB)
      DIMENSION S(*)
*. Specific input
      DIMENSION DENSI(LDENSI,NDENSI)
      DIMENSION DCOEF(*)
*. Output
      DIMENSION D_TILDE(*)
*. Local scratch
      DIMENSION SCR1(*), SCR2(*), SCR3(*), SCR4(*)
*
      NTEST = 00
*
*. Obtain coefficients defining D_BAR
C     OBTAIN_DBAR_COEF(COEF,DBAR_COEF,IREFDENS,NDENSI,INORM)
      CALL OBTAIN_DBAR_COEF(DCOEF,SCR2,IREFDENS,NDENSI,1)
*. Obtain D_BAR in SCR1
C  MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(DENSI,SCR2,SCR1,LDENSI,NDENSI,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The D_BAR density '
        CALL APRBLM2(SCR1,NAOS,NAOS,NSMOB,0)      
      END IF
*. And purify 
*. First expand D_BAR and S to complete form 
* D_BAR expanded in SCR3
* S     expanded in SCR2
      CALL TRIPAK_BLKM(SCR2,S,2,NAOS,NSMOB)
      CALL COPVEC(SCR1,SCR3,LDENSI)
C     CALL TRIPAK_BLKM(SCR3,SCR1,2,NAOS,NSMOB)
      MAXIT_PUR = 10
      CALL MCWEENY_PUR(SCR3,D_TILDE,SCR2,SCR1,SCR4,MAXIT_PUR,NSMOB,NAOS)
C     MCWEENY_PUR(DIN,DOUT,S,SCR1,SCR2,MAXIT,NSMOB,NOBPSM)
*. And normalize - if required
      IF(INORMA.EQ.1) THEN
*. Obtain D_TILDE S in SCR1
        CALL MULT_BLOC_MAT(SCR1,D_TILDE,SCR2,NSMOB,NAOS,NAOS,NAOS,NAOS,
     &                     NAOS,NAOS,0)
C                TRACE_BLK_MAT(A,NBLK,LBLK,IPAK)
        TRACE =  TRACE_BLK_MAT(SCR1,NSMOB,NAOS,0)
        COEF = DBLE(NEL)/TRACE
        CALL SCALVE(D_TILDE,COEF,LDENSI)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================== ' 
        WRITE(6,*) ' The D_tilde density '
        WRITE(6,*) ' ==================== ' 
        WRITE(6,*)
        CALL APRBLM2(D_TILDE,NAOS,NAOS,NSMOB,0)      
      END IF
*
      RETURN
      END 
      SUBROUTINE MCWEENY_PUR(DIN,DOUT,S,SCR1,SCR2,MAXIT,NSMOB,NOBPSM)
*
* Purify input density DIN by repeated use of McWeeny purifucation
* On input, DIN and S are  supposed to be in symmetryblocked complete form
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. General Input
      INTEGER NOBPSM(NSMOB)
*. Input
      DIMENSION DIN(*),S(*)
*. Output
      DIMENSION DOUT(*)
*. Scratch
      DIMENSION SCR1(*),SCR2(*)
*
      NTEST = 00
      LEND = NDIM_1EL_MAT(1,NOBPSM,NOBPSM,NSMOB,0)
C?    WRITE(6,*) 'LEND = ', LEND
*. For mo-density, it is the density/2 that is idempotent, so scale
      I_MO_OR_MOS = 1
      IF(I_MO_OR_MOS.EQ.1) THEN
         FACTOR = 0.5D0
      ELSE 
         FACTOR = 1.0D0
      END IF
      CALL COPVEC(DIN,DOUT,LEND)
      CALL SCALVE(DOUT,FACTOR,LEND)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial (scaled) D '
        CALL  APRBLM2(DOUT,NOBPSM,NOBPSM,NSMOB,0)
      END IF


      ICONVER = 0
      TEST = 1.0D-12
      DO IT = 1, MAXIT
*. Obtain 3DSD - 2 DSDSD with D in DOUT
*. DS in SCR1
C     MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
C?      WRITE(6,*) ' input D, S to MATMULT '
C?        CALL  APRBLM2(DOUT,NOBPSM,NOBPSM,NSMOB,0)
C?        CALL  APRBLM2(S   ,NOBPSM,NOBPSM,NSMOB,0)
  

        CALL MULT_BLOC_MAT(SCR1,DOUT,S,NSMOB,NOBPSM,NOBPSM,NOBPSM,
     &                     NOBPSM,NOBPSM,NOBPSM,0)
C?      WRITE(6,*) ' NEW DS '
C?        CALL  APRBLM2(SCR1,NOBPSM,NOBPSM,NSMOB,0)
*.DSDS in SCR2
        CALL MULT_BLOC_MAT(SCR2,SCR1,SCR1,NSMOB,NOBPSM,NOBPSM,NOBPSM,
     &                     NOBPSM,NOBPSM,NOBPSM,0)
C?      WRITE(6,*) ' NEW DSDS '
C?        CALL  APRBLM2(SCR2,NOBPSM,NOBPSM,NSMOB,0)
*
*.3DS-2DSDS in SCR1
        CALL VECSUM(SCR1,SCR1,SCR2,3.0D0,-2.0D0,LEND)
*.3DSD-2DSDSD in SCR2
        CALL MULT_BLOC_MAT(SCR2,SCR1,DOUT,NSMOB,NOBPSM,NOBPSM,NOBPSM,
     &                     NOBPSM,NOBPSM,NOBPSM,0)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Purified density in iteration ', IT
          CALL  APRBLM2(SCR2,NOBPSM,NOBPSM,NSMOB,0)
        END IF
*. Difference between 3DSD and 2DSDSD 
        CALL VECSUM(SCR1,SCR2,DOUT,1.0D0,-1.0D0,LEND)
        DIFF2 = INPROD(SCR1,SCR1,LEND)
        DIFF = SQRT(DIFF2)
*
        CALL COPVEC(SCR2,DOUT,LEND)
C?      WRITE(6,*) ' DIFF = ', DIFF
        IF(DIFF.LT.TEST) ICONVER = 1
        IF(ICONVER.EQ.1) GOTO 1001
      END DO
 1001 CONTINUE
*
*. Scale to original convention
      FACTOR2 = 1.0D0/FACTOR
      CALL SCALVE(DOUT,FACTOR2,LEND)
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Output from McWeeny purification '
        IF(ICONVER.EQ.1) THEN
           WRITE(6,*) ' Convergence was obtained in ',IT,
     &                ' iterations '
        ELSE
           WRITE(6,*) ' Convergence was not obtained in ',MAXIT,
     &                ' iterations '
        END IF
        WRITE(6,*) ' Norm of last correction : ', DIFF
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Final purified matrix '
C              APRBLM2(S,NAOS_ENV,NAOS_ENV,NSMOB,1)
        CALL  APRBLM2(DOUT,NOBPSM,NOBPSM,NSMOB,0)
      END IF
*
      IF(NTEST.GE.1 .AND. ICONVER.EQ.0) THEN
        WRITE(6,*) ' Warning : McWeeny purification did not converge'
      END IF
*
      RETURN
      END
   
*.

C     SUBROUTINE DIIS_HF(NVECOLD,IREFDENSG,NVEC,B,
*
* Perform standard DIIS extrapolation 
*
* Jeppe Olsen, Sept. 2003
*
*. An density,IREFDENSG is choosen 
*. ( if IREFDENSG = 0, then the last density is choosen)
*. and the new density is expanded as 
* Cbar = 
* (C(IREFDENS) + sum(k.ne.irefdens) C(k) D(K))/(1 + sum(k.ne.irefdens) C(k))
*. 
*. The coefficients are obtained by minimizing the norm of the gradients 
*. 
*
* is choosen as the 
C      INCLUDE 'implicit.inc'
      SUBROUTINE GET_APPROX_HFE(E_APPROX,IAPPROX,
     &           DENSIIN,DENSI,NDENSI,
     &           HAO,FAO,SAO,SINV,I_DO_SINV,
     &           SCR1,SCR2,SCR3,SCR4,SCRVEC1,SCRVEC2,
     &           LENC,LENP,NSMOB,NOBPSM,ENUC)
*
* Obtain Approximate Hartree-Fock energy for a given density
* DENSIIN. Various approximations are possible :
*
* IAPPROX = 0 : E_APPROX = Tr FAO()*D   (FAO = WORK(KFAO))
* IAPPROX = 1 : E_APPROX = Tr(H D_proj) + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 2 : E_APPROX = Tr(H D ) +  Tr D_ort G(D_proj) 
*                        + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 3 : E_APPROX is exact energy of Current density
*               E_APPROX = E_APPROX(IAPPROX=2) 
*                        + 0.5D0*D_proj G(D_proj)
* Where G(D_proj) is the two-electron part of F(D_proj)
* where D_proj is the projection of the density DENSIIN on 
* a set of NDENSI densities, DENSI, and D_ort is the part 
* of DENSIIN that is orthogonal to the densities in DENSI
* G(D_proj) is the two-electron part of F(D_proj)
*
* Jeppe Olsen, Oct. 1 2003, sitting in GWATT working 
*              on Jette and Jeppes 28 year anniversary
*
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc' 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*
*. Input
*  ===== 
*
      DIMENSION DENSIIN(LENC),DENSI(LENC,NDENSI)
*. S and H in expanded form
      DIMENSION SAO(LENC), HAO(LENC)
*. F matrices in compact form ... 
      DIMENSION FAO(LENP,NDENSI)
*. S inverted in space of densities - or space for this matrix
      DIMENSION SINV(NDENSI,NDENSI)
*
*  ======
*. Scratch through argument list 
*  ======
*. Two matrices of length MAX(LENC,NDENSI**2)
      DIMENSION SCR1(*), SCR2(*), SCR3(*), SCR4(*)
*. Two vectors of length NDENSI
      DIMENSION SCRVEC1(NDENSI),SCRVEC2(NDENSI)
*
      E1 = 666666.0D0
*
      IF(IAPPROX.EQ.0) THEN
* E_APPROX = Tr (F(NDENSI) * D
*. Expand last Fock matrix to complete form in SCR3
        CALL TRIPAK_BLKM(SCR3,WORK(KFAO),2,NOBPSM,NSMOB)
*. Multiply D and F in SCR2
        CALL MULT_BLOC_MAT(SCR2,SCR3,DENSIIN,NSMOB,NOBPSM,NOBPSM,
     &       NOBPSM,NOBPSM,NOBPSM,NOBPSM,0)
*. and Trace fD
C                 TRACE_BLK_MAT(A,NBLK,LBLK,IPAK)
        E = TRACE_BLK_MAT(SCR2,NSMOB,NOBPSM,0)
        E_APPROX = E
      END IF
*
      IF(IAPPROX.GT.0) THEN
*. Get the projection of DENSIIN on the input densities 
*. Projected density in SCR1, and expansion coefficients in 
*. SCRVEC1
      CALL PROJECT_DENSI_ON_DENSI(DENSIIN,DENSI,LENC,NDENSI,
     &     SAO,SINV,SCR1,SCRVEC1,NSMOB,NOBPSM,SCR2,SCR3,
     &     SCRVEC2,I_DO_SINV,DORIG_NORM,DPROJ_NORM)
*
      END IF
*
      IF(IAPPROX.EQ.1.OR.IAPPROX.EQ.2.OR.IAPPROX.EQ.3) THEN
*
* E_APPROX = Tr(H D_proj) + 0.5* Tr D_proj G(D_proj)
*          = 0.5* Tr (h + F(D_proj) * D_proj
*. Obtain F(D_proj) in SCR4 in packed form 
C  MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
        CALL MATVCC(FAO,SCRVEC1,SCR4,LENP,NDENSI,0)
*. Expand F(D_proj) in expanded form in SCR3
        CALL TRIPAK_BLKM(SCR3,SCR4,2,NOBPSM,NSMOB)
*. In SCR3 we now have sum(i) scrvec1(i) h + G(D_PROJ), but we 
*. wanted h + G(D_proj) so :
        SUM = 0.0D0
        DO I = 1, NDENSI
          SUM = SUM + SCRVEC1(I)
        END DO
        COEF = 1.0D0 - SUM
        ONE = 1.0D0
        CALL VECSUM(SCR3,SCR3,HAO,ONE,COEF,LENC)
*
C?      WRITE(6,*) ' F(Dproj) obtained as linear comb '
C?      CALL APRBLM2(SCR3,NOBPSM,NOBPSM,NSMOB,0)      
*. Save SCR4 for future work
        CALL COPVEC(SCR3,SCR4,LENC)
*. (h + f(D_proj)) in SCR3
        ONE = 1.0D0
        CALL VECSUM(SCR3,SCR3,HAO,ONE,ONE,LENC)
C?      WRITE(6,*) ' H + F(DPROJ) in ...HFE '
C?      CALL WRT_AOMATRIX(SCR3,0)
*. (h + f(D_proj))D_proj in SCR2
        CALL MULT_BLOC_MAT(SCR2,SCR3,SCR1,NSMOB,NOBPSM,NOBPSM,
     &       NOBPSM,NOBPSM,NOBPSM,NOBPSM,0)
*. and Trace ( h + f(D_proj) D_proj
C                 TRACE_BLK_MAT(A,NBLK,LBLK,IPAK)
        E1 = 0.5D0*TRACE_BLK_MAT(SCR2,NSMOB,NOBPSM,0)
        E_APPROX = E1
      END IF
      IF (IAPPROX.EQ.2.OR.IAPPROX.EQ.3) THEN
* IAPPROX = 2 : E_APPROX2 = Tr(H D ) +  Tr D_ort G(D_proj) 
*                         + 0.5* Tr D_proj G(D_proj)
*                         = E_APPROX1 + Tr D_ort F(D_proj)
*. D_ort in SCR2
        ONE = 1.0D0
        ONEM = -1.0D0
        CALL VECSUM(SCR2,DENSIIN,SCR1,ONE,ONEM,LENC)
C?      WRITE(6,*) ' D_ort : '
C?      CALL APRBLM2(SCR2,NOBPSM,NOBPSM,NSMOB,0)      
*. D_ort F(D_proj) in SCR3
        CALL MULT_BLOC_MAT(SCR3,SCR2,SCR4,NSMOB,NOBPSM,NOBPSM,
     &       NOBPSM,NOBPSM,NOBPSM,NOBPSM,0)
        EDEL = TRACE_BLK_MAT(SCR3,NSMOB,NOBPSM,0)
        E_APPROX = E1 + EDEL
      END IF
      IF(IAPPROX.EQ.3) THEN
*. obtain the second order contribution 
      CALL GET_FOCK(SCR4,SCR3,SCR2,DOGDO,0,2)
C          GET_FOCK(FOCK,H,P,ENERGY,ENUC,I12) 
C?    WRITE(6,*) ' DOGDO term = ', DOGDO
      END IF

*
* TEST ZONE 
*
*. For test : Obtain also D_ORT G(DORT) term
CT    CALL GET_FOCK(SCR4,SCR3,SCR2,DOGDO,0,2)
CT    GET_FOCK(FOCK,H,P,ENERGY,ENUC,I12) 
CT    WRITE(6,*) ' DOGDO term = ', DOGDO
*. Tr H D(proj)
*. Pack H - required by GET_FOCK
CT     CALL TRIPAK_BLKM(HAO,SCR4,1,NOBPSM,NSMOB)
CT     CALL GET_FOCK(SCR3,SCR4,SCR1,E1PROJ,0,1)
*. Tr H D(ort)
CT     CALL GET_FOCK(SCR3,SCR4,SCR2,E1ORT ,0,1)
CT    GET_FOCK(FOCK,H,P,ENERGY,ENUC,I12) 
CT         TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
*. Dproj G(Dproj)
CT     CALL GET_FOCK(SCR3,SCR4,SCR1,DPGDP,0,2)
*
*. And add the nuclear energy
      E_APPROX = E_APPROX + ENUC
*
      IF(IAPPROX.EQ.3) E_APPROX = E_APPROX + DOGDO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Approximate HF energy evaluator : '
        IF(IAPPROX.EQ.0) 
     &   WRITE(6,*) ' Tr F(NDENSI) D = ', E_APPROX
        IF(IAPPROX.GE.1)
     &   WRITE(6,*) ' 0.5* Tr (h + F(D_proj) * D_proj = ',E1
        IF(IAPPROX.GE.2) WRITE(6,*) '  Tr D_ort F(D_proj) = ', EDEL
        IF(IAPPROX.EQ.3) WRITE(6,*) ' 1/2 D_ort G(D_ort)  = ',
     &        DOGDO
        WRITE(6,'(A,E25.15)') '  Final energy approx      = ', 
     &  E_APPROX
*. From the test zone
CT      WRITE(6,*) ' E1PROJ = ', E1PROJ
CT      WRITE(6,*) ' E1ORT  = ', E1ORT
CT      WRITE(6,*) ' DPGDP  = ',DPGDP
      END IF
*
      RETURN
      END
      SUBROUTINE GET_APPROX_HFEM(E_APPROX,IAPPROX,DENSI,NDENSI)
*
* Master routine for calculating approximate HF energy for a 
* given density by projecting density DENSI on a set of 
* NDENSI densities 
* The collected densities and Fock matrices are stored in 
* KDAO_COLLECT and KFAO_COLLECT
*
* IAPPROX = 1 : E_APPROX = Tr(H D_proj) + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 2 : E_APPROX = Tr(H D ) +  Tr D_ort G(D_proj) 
*                        + 0.5* Tr D_proj G(D_proj)
*
* IAPPROX = 3 : Exact energy of density 
* Jeppe Olsen, Oct. 5, 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cecore.inc'
*. Input
      DIMENSION DENSI(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'GET_EM')
*. Local scratch
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
*. Inverse overlap of densities
      CALL MEMMAN(KLSINV,NDENSI*NDENSI,'ADDL  ',2,'SINV  ')
*. SAO in expanded form 
      CALL MEMMAN(KLSAO_E,LENC,'ADDL  ',2,'SAO_E ')
*. HAO in expanded form
      CALL MEMMAN(KLHAO_E,LENC,'ADDL  ',2,'HAO_E ')
*. A few more scratch matrices 
      LEN = MAX(LENC,NDENSI**2)
      CALL MEMMAN(KLSCR1,LEN,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LEN,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LEN,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSCR4,LEN,'ADDL  ',2,'SCR4  ')
*. And two vectors of length NDENSI
      CALL MEMMAN(KLSCRVEC, NDENSI,'ADDL  ',2,'SCRVEC')
      CALL MEMMAN(KLSCRVEC2,NDENSI,'ADDL  ',2,'SCRVEC')
*. Expand H and S to complete symmetrypacked form
      CALL TRIPAK_BLKM(WORK(KLSAO_E),WORK(KSAO),2,NAOS_ENV,NSMOB)
      CALL TRIPAK_BLKM(WORK(KLHAO_E),WORK(KHAO),2,NAOS_ENV,NSMOB)
*. and call the routine doing the job
      CALL GET_APPROX_HFE(E_APPROX,IAPPROX,DENSI,WORK(KDAO_COLLECT),
     &     NDENSI,WORK(KLHAO_E),WORK(KFAO_COLLECT),WORK(KLSAO_E),
     &     WORK(KLSINV),
     &     1,WORK(KLSCR1),WORK(KLSCR2),WORK(KLSCR3),WORK(KLSCR4),
     &     WORK(KLSCRVEC),WORK(KLSCRVEC2),LENC,LENP,NSMOB,NAOS_ENV,
     &     ECORE_ORIG)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'GET_EM')
*
      RETURN
      END 
      SUBROUTINE GET_EAPR_FOR_DTILDE(E_APPROX,IAPPROX,
     &           D_COEF,IREFDN,NDENSI)
*. Obtain approximate energy for DTILDE density 
*. (i.e. purified density) for a given linear combination 
*. of input densities. The expansion of densities is 
*. defined by D_COEF and IREFDN, the latter being 
*. the reference density
*
*. Jeppe Olsen, Oct. 5, 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'gasstr.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E_DTIL')
*. Local scratch
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENGTH = MAX(LENC,NDENSI**2)
*
      CALL MEMMAN(KLSCR1,LENGTH,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LENGTH,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LENC,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSCR4,LENC,'ADDL  ',2,'SCR4  ')
      CALL MEMMAN(KLSCR5,LENC,'ADDL  ',2,'SCR5  ')
      CALL MEMMAN(KLDTILDE,LENC,'ADDL  ',2,'DTILDE')
*
      NEL = NELFTP(1) + NELFTP(2)
C?    WRITE(6,*) ' Number of electrons = ', NEL
      CALL OBTAIN_D_TILDE(WORK(KLDTILDE),IREFDN,D_COEF,NDENSI,LENC,
     &     WORK(KDAO_COLLECT),1,NEL,WORK(KSAO),NSMOB,NAOS_ENV,
     &     WORK(KLSCR1),WORK(KLSCR2),WORK(KLSCR3),WORK(KLSCR4))
*. And the  energy
C     GET_APPROX_HFEM(E_APPROX,IAPPROX,DENSI,NDENSI)
      CALL GET_APPROX_HFEM(E_APPROX,IAPPROX,WORK(KLDTILDE),NDENSI)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E_DTIL')
*
      RETURN
      END 
      SUBROUTINE GET_E1E2_FOR_EHF_APPROX_FD(E0,E1,E2,REFC,NDENSI,IREFDN,
     &           ISCALE_STEP,SCALE_STEP)
*
*. Obtain gradient and Hessian for approximate HF energy function
*. at point defined by REFC 
*
* If ISCALE_STEP = 1, then the finite difference evalaluation 
*                     currently used is done with 
*                     steps DELTA*SCALE_STEP(I) for variable I
*
* Jeppe Olsen, Oct. 2003
*              ISCALE_STEP,SCALE_STEP added Oct. 2004
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'gasstr.inc'
*. Input
      DIMENSION REFC(NDENSI-1), SCALE_STEP(*)
*. Output
      DIMENSION E1(NDENSI-1), E2(NDENSI-1,NDENSI-1)
*
      IDUM = 0
*. Method use to approximate energy 2 : Standard approx E(DBAR)
*                                   3 : Exact E(DBAR)
      IM = 2
      IF(IM.EQ.3) THEN
        WRITE(6,*) ' Note : exact Gradient and Hessian constructed '
      END IF
      IF(ISCALE_STEP.EQ.0) THEN
       WRITE(6,*) ' Steps are unscaled '
      ELSE 
       WRITE(6,*) ' Steps are scaled '
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_E ')
*. Local scratch for storing expansion coefficients
      CALL MEMMAN(KLDCOEF,NDENSI,'ADDL  ',2,'DCOEF ')
*.  Energy at point of expansion
      CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
      CALL GET_EAPR_FOR_DTILDE(E0,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*. step for finite difference 
      DELTA = 0.003D0
*. Gradient and diagonal Hessian elements 
      ID_EFF = 0
      DO ID = 1, NDENSI
        IF(ID.NE.IREFDN) THEN
          ID_EFF = ID_EFF + 1
          IF(ISCALE_STEP.EQ.0) THEN
            DELTAI = DELTA
          ELSE 
            DELTAI = DELTA*SCALE_STEP(ID_EFF)
          END IF
* E(+Delta)
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) + DELTAI
C     GET_EAPR_FOR_DTILDE(EAPPROX,IAPPROX,
C    &           D_COEF,IREFDN,NDENSI)
          CALL GET_EAPR_FOR_DTILDE(EP1,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*. E(-Delta)
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) - DELTAI
          CALL GET_EAPR_FOR_DTILDE(EM1,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*. E(+2 Delta)
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) + 2.0D0*DELTAI
          CALL GET_EAPR_FOR_DTILDE(EP2,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*. E(-2 Delta)
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) -2.0D0*DELTAI
          CALL GET_EAPR_FOR_DTILDE(EM2,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*. And we can obtain gradient and diagonal elements
           E1(ID_EFF) 
     &  = (8.0D0*EP1-8.0D0*EM1-EP2+EM2)/(12.0D0*DELTAI)
C?         WRITE(6,*) ' ID, ID_EFF, E1(ID_EFF) = ',
C?   &                  ID, ID_EFF, E1(ID_EFF)
           E2(ID_EFF,ID_EFF)  
     &  =  (16.0D0*(EP1+EM1-2.0D0*E0)-EP2-EM2+2.0D0*E0)/
     &     (12.0D0*DELTAI**2)
        END IF
      END DO
*. And then the non-diagonal Hessian elements
      ID_EFF = 0
      DO ID = 1, NDENSI
       IF(ID.NE.IREFDN) THEN
        ID_EFF = ID_EFF + 1
        IF(ISCALE_STEP.EQ.0) THEN
          DELTAI = DELTA
        ELSE 
          DELTAI = DELTA*SCALE_STEP(ID_EFF)
        END IF
        JD_EFF = 0
        DO JD = 1, ID-1
         IF(JD.NE.IREFDN) THEN
*EP1P1
          JD_EFF = JD_EFF + 1
          IF(ISCALE_STEP.EQ.0) THEN
            DELTAJ = DELTA
          ELSE 
            DELTAJ = DELTA*SCALE_STEP(JD_EFF)
          END IF
C?        WRITE(6,*) ' ID_EFF, JD_EFF = ', ID_EFF, JD_EFF
*
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) + DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) + DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EP1P1,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*EM1M1
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) -DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) -DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EM1M1,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*EP1M1
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) + DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) - DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EP1M1,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*EM1P1
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) - DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) + DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EM1P1,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*EP2P2
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) + 2.0D0*DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) + 2.0D0*DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EP2P2,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*EM2M2
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) - 2.0D0*DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) - 2.0D0*DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EM2M2,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*EP2M2
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) + 2.0D0*DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) - 2.0D0*DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EP2M2,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*EM2P2
          CALL COPVEC(REFC,WORK(KLDCOEF),NDENSI-1)
          WORK(KLDCOEF-1+ID_EFF) = WORK(KLDCOEF-1+ID_EFF) - 2.0D0*DELTAI
          WORK(KLDCOEF-1+JD_EFF) = WORK(KLDCOEF-1+JD_EFF) + 2.0D0*DELTAJ
          CALL GET_EAPR_FOR_DTILDE(EM2P2,IM,WORK(KLDCOEF),IREFDN,NDENSI)
*
          G1 = EP1P1-EP1M1-EM1P1+EM1M1 
          G2 = EP2P2-EP2M2-EM2P2+EM2M2 
C?        WRITE(6,*) ' EP1P1 EP1M1 EM1P1 EM1M1 = ',
C?   &                 EP1P1,EP1M1,EM1P1,EM1M1
C?        WRITE(6,*) ' EP2P2 EP2M2 EM2P2 EM2M2 = ',
C?   &                 EP2P2,EP2M2,EM2P2,EM2M2
C?        WRITE(6,*) ' G1, G2 = ', G1,G2
C?        WRITE(6,*) ' G1/4*DELTA**2 = ',  G1/(4*DELTAI*DELTAJ)
C?        WRITE(6,*) ' G2/16*DELTA**2 = ',  G2/(16*DELTAI*DELTAJ)
          E2(ID_EFF,JD_EFF) = (16.0D0*G1-G2)/(48*DELTAI*DELTAJ)
          E2(JD_EFF,ID_EFF) = E2(ID_EFF,JD_EFF) 
         END IF
        END DO
       END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Gradient obtained by finite difference '
        WRITE(6,*) ' ======================================='
        CALL WRTMAT(E1,1,NDENSI-1,1,NDENSI-1)
        WRITE(6,*)
        WRITE(6,*) ' Hessian obtained by finite difference '
        WRITE(6,*) ' ======================================'
        CALL WRTMAT(E2,NDENSI-1,NDENSI-1,NDENSI-1,NDENSI-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_E ')
*
      RETURN
      END
      SUBROUTINE OPTIM_EHFAPR(DOPT,NDENSI,IREFDNS,PROJ_COEF,
     &                        IPROJ_DOPT)
*
* Optimize approximate HF energy function in space 
* of purified linear combination of known densities
*
* Jeppe Olsen, Oct. 2003  
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
COLD  REAL*8 INPROD
      REAL*8 INNER_PRODUCT_MAT
*. Output
      DIMENSION DOPT(*), PROJ_COEF(NDENSI)
*. 
      NTEST = 10  
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'OP_HFA')
*
* Set up gradient and Hessian 
*. Space for gradient and Hessian
      CALL MEMMAN(KLE1,NDENSI-1     ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLE2,(NDENSI-1)**2,'ADDL ',2,'E2_HFA')
*. For diagonalization 
      LENP = (NDENSI-1)*NDENSI/2
      CALL MEMMAN(KLE2_P,LENP,'ADDL  ',2,'E2_P  ')
      LENC = (NDENSI-1)**2
      CALL MEMMAN(KLU,LENC,'ADDL  ',2,'EIGV  ')
      CALL MEMMAN(KLE1TR,NDENSI-1    ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLDELTATR,NDENSI   ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLDELTA  ,NDENSI   ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLREFC   ,NDENSI-1 ,'ADDL ',2,'REFC  ')
      CALL MEMMAN(KLSCALE  ,NDENSI-1 ,'ADDL ',2,'SCALE ')
*. Set up metric giving norm of change of D_bar compared to 
* reference density : 
* ||D_bar - D_irefdens!! = sum_ij C_I S_IJ C_J
      LEND = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENDBAR = (NDENSI-1)**2
      CALL MEMMAN(KLSDBAR,LENDBAR,'ADDL  ',2,'SDBAR ')
      CALL DBAR_CHANGE_METRIC(WORK(KDAO_COLLECT),NDENSI,LEND,IREFDNS,
     &                        WORK(KLSDBAR)) 
C          DBAR_CHANGE_METRIC(DENSI,NDENSI,LDENS,IREFDENS,SDBAR)
*. Obtain diaoginal elements in WORK(KLSCALE)
C  COPDIA(A,VEC,NDIM,IPACK)
      CALL COPDIA(WORK(KLSDBAR),WORK(KLSCALE),NDENSI-1,0)
*. Obtain scaling pt used for finite difference - used to 
*. make sure all changes corresponds to the same change 
*. of density matrix 
      DO I = 1, NDENSI-1
       WORK(KLSCALE-1+I) = 1.0D0/SQRT(WORK(KLSCALE-1+I))
      END DO
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' 1/SQRT(||D_I-D_0||) array '
        CALL WRTMAT(WORK(KLSCALE),1,NDENSI-1,1,NDENSI-1)
      END IF
*
C     MAX_NR_IT = 16
*. Number of allowed NR iterations
      MAX_NR_IT = 6
*. Initial allowed total max step 
      NEL = NELFTP(1) + NELFTP(2)
      XMAX_TOTSTEP = 0.05D0*FLOAT(NEL)
*. Initial allowed max step in NR
      XMAX_NRSTEP = 0.03D0*FLOAT(NEL)
*. It semms like starting out with a small step and increase is better so
      XMAX_TOTSTEP = 0.5
      XMAX_NRSTEP = 0.3
*. The above continue is a branching point to 
*. Which we will step if the XMAX_TOTSTEP should be changed
      MAX_REDEF_XMAX = 1
      N_REDEF_XMAX = 0
*. Start by setting reference coefs to zero
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLREFC),ZERO,NDENSI-1)
  999 CONTINUE
      DO IT = 1, MAX_NR_IT
*.. Obtain gradient and Hessian
*
*. Use scaled parameters for finite difference currrently 
*. used to evaluate E1 and E2
       I_SCALE_STEP = 1
       CALL GET_E1E2_FOR_EHF_APPROX_FD(E0,WORK(KLE1),WORK(KLE2),
     &      WORK(KLREFC),NDENSI,IREFDNS,I_SCALE_STEP,WORK(KLSCALE))
*. Solve Trust region controlled NR equations 
       CALL SOLVE_NR_WITH_GEN_TRCTL(WORK(KLE1),WORK(KLE2),
     &      WORK(KLSDBAR),WORK(KLDELTA),NDENSI-1,XMAX_NRSTEP)
*. Update coefficients
       ONE = 1.0D0
       CALL VECSUM(WORK(KLREFC),WORK(KLREFC),WORK(KLDELTA),
     &             ONE,ONE,NDENSI-1)
       IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Updated coefficients '
         CALL WRTMAT(WORK(KLREFC),1,NDENSI-1,1,NDENSI-1)
       END IF
*. Test length of total step - Use SDBAR metric for change
       SSS = V1T_MAT_V2(WORK(KLREFC),WORK(KLREFC),WORK(KLSDBAR),
     &                  NDENSI-1)
       XNORMT = SQRT(SSS)
       I_RESCALED = 0
       IF(XNORMT.GT.XMAX_TOTSTEP) THEN
*. Scale total step and exit
         SCALE = XMAX_TOTSTEP/XNORMT
         CALL SCALVE(WORK(KLREFC),SCALE,NDENSI-1)
         I_RESCALED = 1
         WRITE(6,*) ' Exit as max step exceeded allow max '
         GOTO 1001
       END IF
*
      END DO
 1001 CONTINUE
*     ^ End of loop over iterations
*
* Analysis of result and perhaps projection of Dopt on 
* current space of densities 
*
      LEN_MAT = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LEN_MATX = LEN_MAT
      LENGTH = MAX(LEN_MAT,NDENSI**2)
      CALL MEMMAN(KLSCR1,LENGTH,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LENGTH,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LEN_MAT,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSCR4,LEN_MAT,'ADDL  ',2,'SCR4  ')
      CALL MEMMAN(KLD_TILDE,LEN_MAT,'ADDL  ',2,'D_TILD')
      CALL MEMMAN(KLX,LEN_MAT,'ADDL  ',2,'D_TILD')
      CALL MEMMAN(KLSAOE,LEN_MAT,'ADDL  ',2,'SAOE  ')
      CALL MEMMAN(KLCDBAR,NDENSI,'ADDL  ',2,'CDBAR ')
      CALL MEMMAN(KLDBAR,LEN_MAT,'ADDL  ',2,'DBAR  ')
      CALL MEMMAN(KLDCOEF,NDENSI,'ADDL  ',2,'DCOEF ')
*. Coefficients of Dbar
      CALL OBTAIN_DBAR_COEF(WORK(KLREFC),WORK(KLCDBAR),IREFDNS,
     &       NDENSI,1)
*. Unpack S 
        CALL TRIPAK_BLKM(WORK(KLSAOE),WORK(KSAO),2,NAOS_ENV,NSMOB)
*
      I_CHECK_NORMS = 1
      IF(I_CHECK_NORMS.EQ.1) THEN
*. Obtain DBAR - DREF and check norm ( something fishy is happening I think)
        CALL OBTAIN_DBAR_COEF(WORK(KLREFC),WORK(KLDCOEF),IREFDNS,
     &       NDENSI,1)
        WORK(KLDCOEF-1+IREFDNS) = WORK(KLDCOEF-1+IREFDNS) - 1.0D0
        CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLDCOEF),WORK(KLSCR1),
     &              LEN_MAT,NDENSI,0)
*. and obtain D_bar ( for later use)
        ONE = 1.0D0
        KLDREF = KDAO_COLLECT + (IREFDNS-1)*LEN_MAT
        CALL VECSUM(WORK(KLDBAR),WORK(KLSCR1),WORK(KLDREF),
     &              1.0D0,1.0D0,LEN_MAT)
*. Unpack S 
        CALL TRIPAK_BLKM(WORK(KLSAOE),WORK(KSAO),2,NAOS_ENV,NSMOB)
*. and the square-norm of Dbar-Dref
         DELNORM1 = INNER_PRODUCT_MAT(WORK(KLSCR1),WORK(KLSCR1),
     &              WORK(KLSAOE),WORK(KLSCR2),WORK(KLSCR3),
     &              NSMOB,NAOS_ENV)
         WRITE(6,*) 'Square norm of Dbar - Dref = ', DELNORM1
      END IF
*
      NEL = NELFTP(1) + NELFTP(2)
      CALL OBTAIN_D_TILDE(WORK(KLD_TILDE),IREFDNS,WORK(KLREFC),
     &     NDENSI,LEN_MATX,WORK(KDAO_COLLECT),1,NEL,WORK(KSAO),
     &     NSMOB,NAOS_ENV,
     &     WORK(KLSCR1),WORK(KLSCR2),WORK(KLSCR3),WORK(KLSCR4))
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' Optimized D-tilde density '
        CALL APRBLM2(WORK(KLD_TILDE),NAOS_ENV,NAOS_ENV,NSMOB,0)
      END IF
*
      IF(I_CHECK_NORMS.EQ.1) THEN
*. Obtain DTILDE - DREF and check norm ( something fishy is happening I think)
        ONE = 1.0D0
        ONEM = -1.0D0
        KLDREF = KDAO_COLLECT + (IREFDNS-1)*LEN_MAT 
        CALL VECSUM(WORK(KLSCR1),WORK(KLD_TILDE),WORK(KLDREF),
     &              ONE,ONEM,LEN_MAT)
*. and the square-norm
         DELNORM2 = INNER_PRODUCT_MAT(WORK(KLSCR1),WORK(KLSCR1),
     &              WORK(KLSAOE),WORK(KLSCR2),WORK(KLSCR3),
     &              NSMOB,NAOS_ENV)
         WRITE(6,*) 'Square norm of Dtilde - Dref = ', DELNORM2
*. D_tilde - D_bar  in KLSCR1
         CALL VECSUM(WORK(KLSCR1),WORK(KLD_TILDE),WORK(KLDBAR),
     &               ONE,ONEM,LEN_MAT)
*
         DELNORM3 = INNER_PRODUCT_MAT(WORK(KLSCR1),WORK(KLSCR1),
     &              WORK(KLSAOE),WORK(KLSCR2),WORK(KLSCR3),
     &              NSMOB,NAOS_ENV)
         WRITE(6,*) 'Square norm of Dtilde - Dbar = ', DELNORM3
      END IF
      IF(IPROJ_DOPT.EQ.1) THEN
        I_SKIP_PROJ = 0
        IF(I_SKIP_PROJ.EQ.1) THEN
          WRITE(6,*) 
     &    ' Dbar taken as optimal density(projection of Dtilde skipped)'
          WRITE(6,*) 
     &    ' Dbar taken as optimal density(projection of Dtilde skipped)'
          WRITE(6,*) 
     &    ' Dbar taken as optimal density(projection of Dtilde skipped)'
          CALL OBTAIN_DBAR_COEF(WORK(KLREFC),WORK(KLDCOEF),IREFDNS,
     &         NDENSI,1)
C            MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
          CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLDCOEF),DOPT,
     &              LEN_MAT,NDENSI,0)
          CALL COPVEC(WORK(KLDCOEF),PROJ_COEF,NDENSI)
        ELSE
*. Project the obtained density on the actual densities
          CALL MEMMAN(KLSINV,NDENSI**2,'ADDL   ',2,'SINV  ')
          CALL MEMMAN(KLPROJ_COEF,NDENSI,'ADDL   ',2,'PROJ_C')
          CALL MEMMAN(KLSCRVEC,NDENSI,'ADDL   ',2,'SCRVEC')
*On return DOPT contains Dtilde projected on densities 
          CALL PROJECT_DENSI_ON_DENSI(WORK(KLD_TILDE),
     &          WORK(KDAO_COLLECT),LEN_MAT,NDENSI,
     &          WORK(KLSAOE),WORK(KLSINV),DOPT,WORK(KLPROJ_COEF),
     &          NSMOB,NAOS_ENV,WORK(KLSCR1),WORK(KLSCR2),WORK(KLSCRVEC),
     &          1,DORIG_NORM,DPROJ_NORM)
          CALL COPVEC(WORK(KLPROJ_COEF),PROJ_COEF,NDENSI)
          IF(I_CHECK_NORMS.EQ.1) THEN
*. Obtain Dtilde(proj) - DREF and check norm 
*. ( something fishy is happening I think)

            ONE = 1.0D0
            ONEM = -1.0D0
            KLDREF = KDAO_COLLECT + (IREFDNS-1)*LEN_MAT 
            CALL VECSUM(WORK(KLSCR1),DOPT,WORK(KLDREF),
     &                  ONE,ONEM,LEN_MAT)
*. and the square-norm
             DELNORM4 = INNER_PRODUCT_MAT(WORK(KLSCR1),WORK(KLSCR1),
     &                  WORK(KLSAOE),WORK(KLSCR2),WORK(KLSCR3),
     &                  NSMOB,NAOS_ENV)
             WRITE(6,*) 'Square norm of Dtilde(proj)-Dref = ', DELNORM4
*. Norm of Dtilde(proj) - Dtilde 
             CALL VECSUM(WORK(KLSCR1),DOPT,WORK(KLD_TILDE),
     &                   ONE,ONEM,LEN_MAT)
             DELNORM5 = INNER_PRODUCT_MAT(WORK(KLSCR1),WORK(KLSCR1),
     &                  WORK(KLSAOE),WORK(KLSCR2),WORK(KLSCR3),
     &                  NSMOB,NAOS_ENV)
             WRITE(6,*) 'Square norm of Dtilde(proj)-Dtilde = ',DELNORM5
*. Norm of Dtilde(proj) - Dbar 
             CALL VECSUM(WORK(KLSCR1),DOPT,WORK(KLDBAR),
     &                   ONE,ONEM,LEN_MAT)
             DELNORM6 = INNER_PRODUCT_MAT(WORK(KLSCR1),WORK(KLSCR1),
     &                  WORK(KLSAOE),WORK(KLSCR2),WORK(KLSCR3),
     &                  NSMOB,NAOS_ENV)
             WRITE(6,*) 'Square norm of Dtilde(proj)-Dbar = ', DELNORM6
* 
          END IF
*         ^ End if norms should be checked 
        END IF
*       ^ End of switch if projection should be skipped - even if requested 
      ELSE
        CALL COPVEC(WORK(KLD_TILDE),DOPT,LEN_MAT)
      END IF
*     ^ End of density should be projected 
      IF(I_CHECK_NORMS.EQ.1) THEN
*. Various norms were calculated in the above.
*. Use this to see if the total step  was too small - or too large
* ||Dtilde-D_bar||/||Dbar-D_ref|| gives info about the ratio
* between the total step and purified step. 
         RATIO = SQRT(DELNORM3/DELNORM1)
         WRITE(6,*) ' ||Dtilde-D_bar||/||D_bar-D_ref|| = ',RATIO
         I_WILL_REDO = 0
         IF((I_RESCALED.EQ.1) .AND. (RATIO.LT.0.4) ) THEN
           XMAX_TOTSTEP = 1.5*XMAX_TOTSTEP
           I_WILL_REDO = 1
           N_REDEF_XMAX = N_REDEF_XMAX + 1
         END IF
         IF( RATIO.GT.0.7) THEN
*. The ratio was too large. It is my experience that 
*. this usually is helped by reducing the step, so I will skip this
C          XMAX_TOTSTEP = 0.5*XMAX_TOTSTEP
C          I_WILL_REDO = 1
C          XMAX_NRSTEP = 0.5*XMAX_NRSTEP
C          N_REDEF_XMAX = N_REDEF_XMAX + 1
         END IF
         IF(I_WILL_REDO.EQ.1.AND.N_REDEF_XMAX.LE.MAX_REDEF_XMAX)THEN  
          WRITE(6,*) ' Iterations will be redone with '
          WRITE(6,*) '  XMAX_TOTSTEP = ', XMAX_TOTSTEP
          WRITE(6,*) '  XMAX_NRSTEP = ', XMAX_NRSTEP
          GOTO 999
         END IF
      END IF
*     ^ End if norms were checked
*
*. Find approximate energy of optimized tilde density
        CALL GET_APPROX_HFEM(E_APPROX,2,WORK(KLD_TILDE),NDENSI)
        WRITE(6,*) ' Expected energy of D_tilde ', E_APPROX
        CALL GET_APPROX_HFEM(E_APPROX,2,DOPT,NDENSI)
        WRITE(6,*) ' Expected energy of D_bar ', E_APPROX
*. Find exact energy of optimized D_tilde -- constructs F(D_ort)
        CALL GET_APPROX_HFEM(E_EX,3,WORK(KLD_TILDE),NDENSI)
        WRITE(6,*) ' Exact energy of D_tilde( F constructed) ', E_EX
*. Find exact energy of optimized D_bar e -- constructs F(D_ort)
        CALL GET_APPROX_HFEM(E_EX,3,DOPT,NDENSI)
        WRITE(6,*) ' Exact energy of D_bar( F constructed) ', E_EX
*
      IITEST = 0
      IF(IITEST.EQ.1) THEN
        WRITE(6,*) ' Energy at point of expansion ', E0
*. Find approximate energy of optimized tilde density
        CALL GET_APPROX_HFEM(E_APPROX,2,WORK(KLD_TILDE),NDENSI)
        WRITE(6,*) ' Expected energy of D_tilde ', E_APPROX
*. Find exact energy of optimized D_tilde -- constructs F(D_ort)
        CALL GET_APPROX_HFEM(E_APPROX,3,WORK(KLD_TILDE),NDENSI)
        WRITE(6,*) ' Exact energy of D_tilde ', E_APPROX
*. Obtain second order approximation to energy of optimized D_tilde
C       VALUE_SECORD_TAYLOR(X,E0,E1,E2,NDIM,SCR)
        E2APR = VALUE_SECORD_TAYLOR(WORK(KLDELTA),E0,WORK(KLE1),
     &          WORK(KLE2),NDENSI-1,WORK(KLDELTATR))
        WRITE(6,*) ' Taylor estimate of E(D_TILDE) ', E2APR
        IF(IPROJ_DOPT.EQ.1) THEN
*. And approximate energy for projected D-tilde
          CALL GET_APPROX_HFEM(E_APPROX,2,DOPT,NDENSI)
          WRITE(6,*) ' Expected energy of D_bar ', E_APPROX
        END IF
      END IF
*
      NTEST = 10
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Optimized density '
        CALL APRBLM2(DOPT,NAOS_ENV,NAOS_ENV,NSMOB,0)
C?      CALL APRBLM2(SCR2,NOBPSM,NOBPSM,NSMOB,0)      
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Coefficients of Dbar as obtained from optim'
        CALL WRTMAT(WORK(KLCDBAR),1,NDENSI,1,NDENSI)
      END IF
      IF(IPROJ_DOPT.EQ.1.AND.NTEST.GE.10) THEN
          WRITE(6,*) ' Coefficients of projected D_TILDE aka D_BAR'
          CALL WRTMAT(PROJ_COEF,1,NDENSI,1,NDENSI)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'OP_HFA')
*
      RETURN
      END
      FUNCTION VALUE_SECORD_TAYLOR(X,E0,E1,E2,NDIM,SCR)
*
* Evaluate second-order taylor expansion for given expansion coefficients X
*
*. Jeppe Olsen, Oct. 2003
*
      INCLUDE 'implicit.inc'
      REAL*8   INPROD
*. input
      DIMENSION E1(NDIM),E2(NDIM,NDIM)
*. Scratch
      DIMENSION SCR(NDIM)
*
      VALUE = E0
*. First order term X(T) * E1
      E1TERM =  INPROD(X,E1,NDIM)
      VALUE = VALUE + E1TERM
*. Second order term
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(E2,X,SCR,NDIM,NDIM,0)
      E2TERM = 0.5D0*INPROD(SCR,X,NDIM)
      VALUE = VALUE + E2TERM
*
      VALUE_SECORD_TAYLOR = VALUE
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Value of second order expansion = ', VALUE
        WRITE(6,*)
        WRITE(6,*) ' Zero-order contribution    = ', E0
        WRITE(6,*) ' First-order contributions  = ', E1TERM
        WRITE(6,*) ' Second-order contributions = ', E2TERM
      END IF
*
      RETURN
      END
      SUBROUTINE GET_F_FOR_SUM_DENSI(FOUT,COEF,NDENSI)
*
* Obtain Fock matrix for a linear combination of densities 
*
* F = h + G(sum_i COEF_i D_i)
*
* The Fock matrices are stored in KFAO_COLLECT
*
*. Jeppe Olsen, October 12, 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      DIMENSION COEF(NDENSI)
*. Output : in packed form as standard Fock matrices
      DIMENSION FOUT(*)
*
      IDUM = 0
C     CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'F_SUMD')
* 
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
*. Obtain sum_i COEF_i F(D_i) in FOUT
C  MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(WORK(KFAO_COLLECT),COEF,FOUT,LENP,NDENSI,0)
*. We now have sum(i) COEF_i h + sum_i COEF_i G(D_i), but we 
*. wanted h + sum_i COEF_i G(D_i) so :
      SUM = 0.0D0
      DO I = 1, NDENSI
        SUM = SUM + COEF(I)
      END DO
      FACTOR = 1.0D0 - SUM
      ONE = 1.0D0
      CALL VECSUM(FOUT,FOUT,WORK(KHAO),ONE,FACTOR,LENP)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Fock matrix generated for density weights '
        CALL WRTMAT(COEF,1,NDENSI,1,NDENSI)
        WRITE(6,*) ' Obtained Fock matrix '
        CALL APRBLM2(FOUT,NAOS_ENV,NAOS_ENV,NSMOB,1) 
      END IF
*
      RETURN
      END
      FUNCTION EHFAPR_FOR_X(X,DEN_REF,NDENSI,IAPPROX)
*
* Obtain approximate HF energy for given antisymmetrix matrix X
* And reference density DEN_REF.
*
* The density is obtained as 
* D = Exp (-XS) DEN_REF Exp (SX)
*
* and the approximate energy is governed by IAPPROX : 
*
* IAPPROX = 1 : E_APPROX = Tr(H D_proj) + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 2 : E_APPROX = Tr(H D ) +  Tr D_ort G(D_proj) 
*                        + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 3 : E_APPROX is exact energy of Current density
*               E_APPROX = E_APPROX(IAPPROX=2) 
*                        + 0.5D0*D_proj G(D_proj)
* Where G(D_proj) is the two-electron part of F(D_proj)
* where D_proj is the projection of the density DENSIIN on 
* a set of NDENSI densities, DENSI, and D_ort is the part 
* of DENSIIN that is orthogonal to the densities in DENSI
* G(D_proj) is the two-electron part of F(D_proj)
*
*. Jeppe Olsen, Oct 14, 2003 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
       DIMENSION X(*), DEN_REF(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EAPR_X')
*. Obtain Actual density
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLDEN,LENC,'ADDL  ',2,'DEN_NW')
      CALL NEW_DENSI_FROM_BCH(DEN_REF,WORK(KLDEN),X,0)
C     NEW_DENSI_FROM_BCH(DIN,DOUT,X,IPURIFY)
*. and approximate energy for density
C          GET_APPROX_HFEM(E_APPROX,IAPPROX,DENSI,NDENSI)
      CALL GET_APPROX_HFEM(E_APPROX,IAPPROX,WORK(KLDEN),NDENSI)
*
       EHFAPR_FOR_X = E_APPROX
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EAPR_X')
      RETURN
      END 
*
      SUBROUTINE NEW_DENSI_FROM_BCH(DIN,DOUT,X,IPURIFY)
*
* Obtain new density from generalized BCH expansion 
*
* DOUT = Exp(-XS ) DIN Exp( SX)
*
* In input X is assumed to be in lower half packed form - without diagonal
*
* If IPURIFY = 1, then the obtained density is purified 
*
*. Jeppe Olsen, October 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      DIMENSION DIN(*), X(*)
*. Output
      DIMENSION DOUT(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DN_BCH')
* 
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to NEW_DENSI_FROM_BCH '
        WRITE(6,*) ' Reference density '
        CALL APRBLM2(DIN,NAOS_ENV,NAOS_ENV,NSMOB,0)
      END IF
*
*. Expand X to complete antisymmetric form 
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLXEXP,LENC,'ADDL  ',2,'XEXP  ')
      CALL REFORM_ANTISYM_BLMAT(X,WORK(KLXEXP),NAOS_ENV,NSMOB,
     &                          0,1)
*. Project X so only nonredundant part remains - use DIN in projection
      CALL MEMMAN(KLXPROJ,LENC,'ADDL  ',2,'XPROJ')
      CALL PROJ_NONRED(WORK(KLXEXP),WORK(KLXPROJ),DIN,0)
C          PROJ_NONRED(X,XPROJ,D)
*. Expand S to complete form 
      CALL MEMMAN(KLSAO_E,LENC,'ADDL  ',2,'SAO_E ')
      CALL TRIPAK_BLKM(WORK(KLSAO_E),WORK(KSAO),2,NAOS_ENV,NSMOB)
*. Obtain exp(-XS) DIN exp( SX) 
      MAXORD = 20
C     GENBCH(A,S,X,OUT,NBLK,LBLK,MAXORD)
      CALL GENBCH(DIN,WORK(KLSAO_E),WORK(KLXPROJ),DOUT,NSMOB,NAOS_ENV,
     &            MAXORD)
*
      IF(IPURIFY.EQ.1) THEN
C     MCWEENY_PUR(DIN,DOUT,S,SCR1,SCR2,MAXIT,NSMOB,NOBPSM)
        CALL MEMMAN(KLSCR1,LENC,'ADDL  ',2,'SCR1  ')
        CALL MEMMAN(KLSCR2,LENC,'ADDL  ',2,'SCR2  ')
        CALL MEMMAN(KLDOUT,LENC,'ADDL  ',2,'LDOUT ')
        MAXIT = 10
*. Note : Scaling of MO density to idempotent form is done in MCWEENY_PUR
        CALL MCWEENY_PUR(DOUT,WORK(KLDOUT),WORK(KLSAO_E),WORK(KLSCR1),
     &                   WORK(KLSCR2),MAXIT,NSMOB,NAOS_ENV)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from NEW_DENSI_FROM_BCH '
        WRITE(6,*) ' Input density : '
        CALL APRBLM2(DIN,NAOS_ENV,NAOS_ENV,NSMOB,0) 
        WRITE(6,*) ' Output density '
        CALL APRBLM2(DOUT,NAOS_ENV,NAOS_ENV,NSMOB,0) 
        WRITE(6,*) ' X-matrix (without diagonal) '
        IOFFP = 1
        DO IBLK = 1, NSMOB
            WRITE(6,*) ' Block ', IBLK
            LEN = NAOS_ENV(IBLK)-1
            CALL PRSYM(X(IOFFP),LEN)
            IOFFP = IOFFP + LEN*(LEN-1)/2
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DN_BCH')
      RETURN
      END

      SUBROUTINE REFORM_ANTISYM_BLMAT(APAK,ACOM,LBLK,NBLK,
     &           IDIAG_IS_IN_APAK,IWAY)
* Reform antisymmetric blocked matrix between full and lower half
* packed form.
*
* IWAY = 1 : Packed to complete form
* IWAY = 2 : complete to packed form
*
*. If IDIAG_IS_IN_APAK = 0, APAK does not contain the (vanishing) diagonal
*
*. Jeppe Olsen, October 2003
*
      INCLUDE 'implicit.inc'
*. Input 
      INTEGER LBLK(NBLK)
*. Input and output
      DIMENSION APAK(*),ACOM(*)
*
      IOFFP = 1
      IOFFC = 1
      DO IBLK = 1, NBLK
        LEN = LBLK(IBLK)
C?      WRITE(6,*) ' IBLK, LEN, IOFFP = ', IBLK,LEN,IOFFP
        CALL REFORM_ANTISYM_MAT(APAK(IOFFP),ACOM(IOFFC),LEN,
     &       IDIAG_IS_IN_APAK,IWAY)
        IOFFC = IOFFC + LEN**2
        IF(IDIAG_IS_IN_APAK.EQ.1) THEN
          IOFFP = IOFFP + LEN*(LEN+1)/2
        ELSE 
          IOFFP = IOFFP + LEN*(LEN-1)/2
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Antisymmetric packed matrix => complete matrix'
        ELSE 
          WRITE(6,*) ' Complete matrix => antisymmetric packed matrix'
        END IF
        IF(IDIAG_IS_IN_APAK.EQ.0) THEN
          WRITE(6,*)
     &    ' Packed antisymmetric matrix does not contain diagonal '
        ELSE 
          WRITE(6,*)
     &    ' Packed antisymmetric matrix does     contain diagonal '
        END IF
*
        WRITE(6,*) ' Complete matrix : '
        CALL APRBLM2(ACOM,LBLK,LBLK,NBLK,0) 
        WRITE(6,*) ' Packed matrix : '
        IF(IDIAG_IS_IN_APAK.EQ.0) THEN
*. Packed matrix has no diagonal - I have no routine to print this so 
          IOFFP = 1
          DO IBLK = 1, NBLK
            WRITE(6,*) ' Block ', IBLK
            LEN = LBLK(IBLK)-1
            CALL PRSYM(APAK(IOFFP),LEN)
            IOFFP = IOFFP + LEN*(LEN-1)/2
          END DO
        ELSE 
          CALL APRBLM2(APAK,LBLK,LBLK,NBLK,1) 
        END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE REFORM_ANTISYM_MAT(APAK,ACOM,NDIM,
     &           IDIAG_IS_IN_APAK,IWAY)
*
* Reform antisymmetric matrix between full and lower half form 
* Lower half form is packed in standard rowwise form
*
* IWAY = 1 : Packed to complete form
* IWAY = 2 : complete to packed form
*
*. If IDIAG_IS_IN_APAK = 0, APAK does not contain the (vanishing) diagonal
*
*. Jeppe Olsen, October 2003
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION APAK(*), ACOM(NDIM,NDIM)
*
      IJ = 0
      DO I = 1, NDIM
        IF(IDIAG_IS_IN_APAK.EQ.1) THEN
          MAXJ = I 
        ELSE 
          MAXJ = I-1
        END IF
        IF(IWAY.EQ.1) THEN
          DO J = 1, MAXJ
            IJ = IJ + 1
            ACOM(I,J) = APAK(IJ)
            ACOM(J,I) =-APAK(IJ)
          END DO
        ELSE 
          DO J = 1, MAXJ
            IJ = IJ + 1
            APAK(IJ) = ACOM(I,J) 
          END DO
        END IF
      END DO
*
      IF(IDIAG_IS_IN_APAK.EQ.0.AND.IWAY.EQ.1) THEN
*. Set diagonal in complete matrix to zero
C  SETDIA(MATRIX,VALUE,NDIM,IPACK)
        ZERO = 0.0D0
        CALL SETDIA(ACOM,ZERO,NDIM,0)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Antisymmetric packed matrix => complete matrix'
        ELSE 
          WRITE(6,*) ' Complete matrix => antisymmetric packed matrix'
        END IF
        IF(IDIAG_IS_IN_APAK.EQ.0) THEN
          WRITE(6,*)
     &    ' Packed antisymmetric matrix does not contain diagonal '
        ELSE 
          WRITE(6,*)
     &    ' Packed antisymmetric matrix does     contain diagonal '
        END IF
*
        WRITE(6,*) ' Complete matrix : '
        CALL WRTMAT(ACOM,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Packed matrix '
        IF(IDIAG_IS_IN_APAK.EQ.0) THEN
          LENGTH = NDIM -1
        ELSE
          LENGTH = NDIM
        END IF
        CALL PRSYM(APAK,LENGTH)
      END IF
*
      RETURN
      END
      SUBROUTINE S_COMMU_BLKMT(A,B,S,COMMU,SCR,NBLK,LBLK)
*
* Obtain S-commutator of a a blocked matrix containing 
* square submatrices 
*
* The S-matrix must be in complete blocked form
*
* Jeppe Olsen, Oct. 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*),B(*),S(*)
      INTEGER LBLK(NBLK)
*. Output
      DIMENSION COMMU(*)
*. Scratch : Must be able to hold largest submatrix
      DIMENSION SCR(*)
*
      IOFF = 1
      DO IBLK = 1, NBLK
         LENGTH = LBLK(IBLK)
         CALL S_COMMU(A(IOFF),B(IOFF),S(IOFF),COMMU(IOFF),SCR,LENGTH)
         IOFF = IOFF + LENGTH**2
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' S-commutator of blocked matrices '
        CALL APRBLM2(COMMU,LBLK,LBLK,NBLK,0)
      END IF
*
      RETURN
      END
      SUBROUTINE GENBCH(A,S,X,OUT,NBLK,LBLK,MAXORD)
*
* General BCH expansion of a matrix :
*
* OUT = exp(-XS) A exp(SX)
*
*. The exponentials are expanded through  order MAXORD
*
* Scratch is located internally
*
* Jeppe Olsen, October 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*
*. Input
      DIMENSION A(*),S(*),X(*)
      INTEGER LBLK(NBLK)
*. Output
      DIMENSION OUT(*)
*
      NTEST = 00
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GENBCH')
* Local scratch : Two matrices, all blocks and one block -
*                 for simplicity/laziness two complete matrices allocated
C?    WRITE(6,*) ' GENBCH : NBLK = ', NBLK
C?    WRITE(6,*) ' LBLK : '
C?    CALL IWRTMA(LBLK,1,NBLK,1,NBLK)
      LENBLM = NDIM_1EL_MAT(1,LBLK,LBLK,NBLK,0)
C?    WRITE(6,*) ' LENBLM = ', LENBLM
*
      CALL MEMMAN(KLSCR1,LENBLM,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LENBLM,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LENBLM,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSCR4,LENBLM,'ADDL  ',2,'SCR4  ')
*. Norm of A - used for convergence 
      ANORM = SQRT(INPROD(A,A,LENBLM))
*. Required threshold for convergence - currently very strict
      THRES = 1.0D-14*ANORM
*. zero-order term
      CALL COPVEC(A,OUT,LENBLM)
      CALL COPVEC(A,WORK(KLSCR4),LENBLM)
      XFAC = 1.0D0
      DO IORD = 1, MAXORD
        XFAC = XFAC*DBLE(IORD)
        FACTOR = 1.0D0/XFAC
*. Obtain next commutator :  SCR4 S X - X S SCR4,  store in SCR2
C            S_COMMU_BLKMT(A,B,S,COMMU,SCR,NBLK,LBLK)
        CALL S_COMMU_BLKMT(WORK(KLSCR4),X,S,WORK(KLSCR2),WORK(KLSCR3),
     &                     NBLK,LBLK)
*. Norm of change
        XNORM = SQRT(INPROD(WORK(KLSCR2),WORK(KLSCR2),LENBLM))*FACTOR
*. and add to result 
        ONE = 1.0D0
        CALL VECSUM(OUT,OUT,WORK(KLSCR2),ONE,FACTOR,LENBLM)
        IF(XNORM.LE.THRES) GOTO 1001
*. and prepare for next it
        CALL COPVEC(WORK(KLSCR2),WORK(KLSCR4),LENBLM)
      END DO
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' General BCH expansion, exp(-XS) A exp( SX) '
        WRITE(6,*) ' Input X and A matrices '
        CALL APRBLM2(X,LBLK,LBLK,NBLK,0)
        CALL APRBLM2(A,LBLK,LBLK,NBLK,0)
        WRITE(6,*)  ' Output matrix : '
        CALL APRBLM2(OUT,LBLK,LBLK,NBLK,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GENBCH')
*
      RETURN
      END
      SUBROUTINE S_COMMU(A,B,S,COMMU,SCR,NDIM)
*
* Obtain S-commutator of two matrices
*
* [A,B]_S = ASB - BSA
*
* Jeppe Olsen, October 2003
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM),S(NDIM,NDIM)
*. Output
      DIMENSION COMMU(NDIM,NDIM)
*. Scratch
      DIMENSION SCR(NDIM,NDIM)
*. AS in SCR
C      MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(SCR,A,S,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,
     &            FACTORC,FACTORAB,0)
*. ASB in COMMU
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(COMMU,SCR,B,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,
     &            FACTORC,FACTORAB,0)
*. BS in SCR
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(SCR,B,S,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,
     &            FACTORC,FACTORAB,0)
*. ASB - BSA in COMMU
      FACTORC = 1.0D0
      FACTORAB = -1.0D0
      CALL MATML7(COMMU,SCR,A,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,
     &            FACTORC,FACTORAB,0)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' S-commutator ASB - BSA will be evaluated for '
        WRITE(6,*) ' A matrix : '
        CALL WRTMAT(A    ,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' B matrix : '
        CALL WRTMAT(B    ,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' S matrix : '
        CALL WRTMAT(S    ,NDIM,NDIM,NDIM,NDIM)
    
        WRITE(6,*) ' S-commutator matrix : '
        CALL WRTMAT(COMMU,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE E1E2_FOR_APPROX_EHF_EXPSX(E0,E1,E2,DEN_REF,X_REF,LENX,
     &                                     NDENSI,IAPPROX)
*
* Obtain gradient and Hessian for approximate HF energy function
* using the Exp SX formalism for the density. 
*
* The energyfunction assumes that D0 density used for the expansion is
* in WORK(KDAO) and the corresponding Fock matrix is on WORK(KFAO)
*
* IAPPROX defines energy approximation 
*
* IAPPROX = 1 : E_APPROX = Tr(H D_proj) + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 2 : E_APPROX = Tr(H D ) +  Tr D_ort G(D_proj) 
*                        + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 3 : E_APPROX is exact energy of current density
*
*. Jeppe Olsen, October 2003, updated with analytical method, may 2005
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
*. Input : Reference density and reference X used for density generation
*          (Densi = Exp-(X_REF + X)S DEN_REF Exp S(X_REF+X)
*
* LENX is the number of redundant parameters in the AO parameterization 
* of the density
      DIMENSION DEN_REF(*),X_REF(*)
*. Output : gradient and Hessian 
      DIMENSION E1(LENX),E2(LENX,LENX)
*
      I_NUM_OR_ANA = 2 
* 1 => Numerical 
* 2 => Analytical
* 3 => Numerical and and analytical 
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E12_ESX')
*
      IF(I_NUM_OR_ANA.EQ.1.OR.I_NUM_OR_ANA.EQ.3) THEN
*. Obtain gradient and Hessian by analytical method
*. Memory for local copy of X coefficients
        CALL MEMMAN(KLX,LENX,'ADDL  ',2,'X_COPY')
*.  Energy at point of expansion
        E0 =  EHFAPR_FOR_X(X_REF,DEN_REF,NDENSI,IAPPROX)
C             EHFAPR_FOR_X(X,DEN_REF,NDENSI,IAPPROX)
*. step for finite difference 
        DELTA = 0.005D0
*. Gradient and diagonal Hessian elements 
        DO I = 1, LENX
* E(+Delta)
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) + DELTA
            EP1 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*. E(-Delta)
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) - DELTA
            EM1 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*. E(+2 Delta)
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) + 2.0D0*DELTA
            EP2 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*. E(-2 Delta)
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) -2.0D0*DELTA
            EM2 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*. And we can obtain gradient and diagonal elements
*
C?         WRITE(6,*) ' E0, EP1, EP2, EM1, EM2 = ',
C?   &                  E0, EP1, EP2, EM1, EM2
             E1(I) 
     &    = (8.0D0*EP1-8.0D0*EM1-EP2+EM2)/(12.0D0*DELTA)
C?           WRITE(6,*) ' ID, ID_EFF, E1(ID_EFF) = ',
C?   &                    ID, ID_EFF, E1(ID_EFF)
             E2(I,I)  
     &    =  (16.0D0*(EP1+EM1-2.0D0*E0)-EP2-EM2+2.0D0*E0)/
     &       (12.0D0*DELTA**2)
      END DO
*. And the non-diagonal Hessian elements
        DO I = 1, LENX
          DO J = 1, I-1
*EP1P1
*
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) + DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) + DELTA
            EP1P1 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*EM1M1
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) -DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) -DELTA
            EM1M1 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*EP1M1
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) + DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) - DELTA
            EP1M1 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*EM1P1
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) - DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) + DELTA
            EM1P1 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*EP2P2
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) + 2.0D0*DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) + 2.0D0*DELTA
            EP2P2 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*EM2M2
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) - 2.0D0*DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) - 2.0D0*DELTA
            EM2M2 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*EP2M2
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) + 2.0D0*DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) - 2.0D0*DELTA
            EP2M2 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*EM2P2
            CALL COPVEC(X_REF,WORK(KLX),LENX)
            WORK(KLX-1+I) = WORK(KLX-1+I) - 2.0D0*DELTA
            WORK(KLX-1+J) = WORK(KLX-1+J) + 2.0D0*DELTA
            EM2P2 =  EHFAPR_FOR_X(WORK(KLX),DEN_REF,NDENSI,IAPPROX)
*
            G1 = EP1P1-EP1M1-EM1P1+EM1M1 
            G2 = EP2P2-EP2M2-EM2P2+EM2M2 
C?        WRITE(6,*) ' EP1P1 EP1M1 EM1P1 EM1M1 = ',
C?   &                 EP1P1,EP1M1,EM1P1,EM1M1
C?        WRITE(6,*) ' EP2P2 EP2M2 EM2P2 EM2M2 = ',
C?   &                 EP2P2,EP2M2,EM2P2,EM2M2
C?        WRITE(6,*) ' G1, G2 = ', G1,G2
C?        WRITE(6,*) ' G1/4*DELTA**2 = ',  G1/(4*DELTA**2)
C?        WRITE(6,*) ' G2/16*DELTA**2 = ',  G2/(16*DELTA**2)
            E2(I,J) = (16.0D0*G1-G2)/(48*DELTA**2)
            E2(J,I) = E2(I,J) 
          END DO
        END DO
      END IF
*     ^ End if gradient and Hessian should be calculated using finite difference
      IF(I_NUM_OR_ANA.EQ.2.OR.I_NUM_OR_ANA.EQ.3) THEN
*. Calculate gradient and Hessian using analytical methods
*. It is *. assumed that X=0
        LENX = LEN_AO_MAT(-1)
        CALL MEMMAN(KLE1AN,LENX,'ADDL  ',2,'E1_AN ')
        CALL MEMMAN(KLE2AN,LENX**2,'ADDL  ',2,'E2_AN ')
        CALL GET_E1_HF_ANALYTICAL(WORK(KLE1AN),WORK(KFAO),DEN_REF,
     &       WORK(KDAO),IAPPROX,NDENSI)
C       GET_E1_HF_ANALYTICAL(E1,FNOT,DSEED,DNOT,IAPPROX,
C    &                                  NDENSI)
        CALL GET_E2_ONESTEP(WORK(KLE2AN),LENX,WORK(KFAO),DEN_REF,
     &       WORK(KDAO),NDENSI,IAPPROX)
C       GET_E2_ONESTEP(E2.LENX,FNOT,DSEED,DNOT,NDENSI,IAPPROX)
      END IF
      IF(I_NUM_OR_ANA.EQ.2) THEN
*. Transfer to permanent arrays
        CALL COPVEC(WORK(KLE1AN),E1,LENX)
        CALL COPVEC(WORK(KLE2AN),E2,LENX**2)
      END IF
*
      IF(I_NUM_OR_ANA.EQ.3) THEN
*. Compare numerical and analytical gradient and Hessian 
*. Check projection af E1(NUM) and E1(ANA)
        WRITE(6,*) ' Check of E1(num) and P(T) E1(num) '
        CALL CHECK_PROJ_VEC(E1,-1,DEN_REF,1)
        WRITE(6,*) ' Check of E1(ana) and P(T) E1(ana) '
        CALL CHECK_PROJ_VEC(WORK(KLE1AN),-1,DEN_REF,1)
C       CHECK_PROJ_VEC(X,IMATFORM,D,ITRNSP)
*. Compare numerical and analytical density 
        WRITE(6,*) ' Comparison of numerical and analytical gradient'
        THRES = 1.0D-7
        CALL CMP2VC(E1,WORK(KLE1AN),LENX,THRES)
        WRITE(6,*) ' Comparison of numerical and analytical Hessian'
        CALL CMP2VC(E2,WORK(KLE2AN),LENX**2,THRES)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Gradient obtained by finite difference '
        WRITE(6,*) ' ======================================='
        CALL WRTMAT(E1,1,LENX,1,LENX)              
        WRITE(6,*)
        WRITE(6,*) ' Hessian obtained by finite difference '
        WRITE(6,*) ' ======================================'
        CALL WRTMAT(E2,LENX,LENX,LENX,LENX)                   
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E12_ESX')
      RETURN
      END 
      SUBROUTINE OPTIM_EHFAPR_EXPSX(DOPT,NDENSI,DEN_REF)
*
* Optimize approximate HF energy using Exp XS formalism.
*
* The energy function being optimized is of the form
*
* E = Tr H D(X) + 0.5*Tr(D_proj) G(D_proj) + Tr(D_ort G(D_proj)
*
* where D_proj is the projection of D(X) on a set of known 
* densities.
*
* The resulting optimized density is returned in DOPT
*
* Jeppe Olsen, Oct. 2003  
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
      REAL*8 INPROD
*. Input : Reference density used as initial point for optimization
      DIMENSION DEN_REF(*)
*. Output
      DIMENSION DOPT(*)
*. 
      NTEST = 10  
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'OP_EXS')
*. Number of parameters  in lower half matrix without diagonal
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      LENX = LENP - NTOOB
C?    WRITE(6,*) ' Number of parameters in antisymmetrix matrix', LENX
*. Space for gradient and Hessian
      CALL MEMMAN(KLE1,LENX     ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLE2,LENX**2,  'ADDL ',2,'E2_HFA')
*. For diagonalization 
      LEN_E2P = LENX*(LENX+1)/2
      CALL MEMMAN(KLE2_P,LEN_E2P,'ADDL  ',2,'E2_P  ')
      CALL MEMMAN(KLU,LENX**2,'ADDL  ',2,'EIGV  ')
      CALL MEMMAN(KLE1TR,LENX    ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLDELTATR,LENX ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLDELTA  ,LENX ,'ADDL ',2,'E1_HFA')
      CALL MEMMAN(KLREFX   ,LENX ,'ADDL ',2,'REFX  ')
*. and a for holding a complete, symmetrypacked matrix 
      CALL MEMMAN(KLMATC   ,LENC ,'ADDL ',2,'MATC  ')
*. The point of expansion will change and is in general stored in DOPT, so 
*. start by copying initial density to DOPT
      CALL COPVEC(DEN_REF,DOPT,LENC)
*
      MAX_NR_IT = 3
      DO IT = 1, MAX_NR_IT
*. Check DSD condition on input density matrix 
        CALL CHECK_IDEMP(DOPT)
C            CHECK_IDEMP(D)
*. Set reference coefs to zero - we are changing point of expansion 
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLREFX),ZERO,LENX)
*. Choose energy form for one-step method
       IAPPROX = 2
       IF(IAPPROX.EQ.3) THEN
         WRITE(6,*) ' Exact Hessian will be calculated '
       ELSE IF (IAPPROX.EQ.0) THEN
         WRITE(6,*) ' Hessian is obtained for Tr F(LAST) D '
       END IF
       CALL E1E2_FOR_APPROX_EHF_EXPSX(E0,WORK(KLE1),WORK(KLE2),DOPT,
     &      WORK(KLREFX),LENX,NDENSI,IAPPROX)
*
*. Diagonalize Hessian 
       CALL DIAG_SYM_MAT(WORK(KLE2),WORK(KLU),WORK(KLE2_P),LENX,0)
*. On output eigenvalues are in WORK(KLE2_P) and eigenvectors are in KLU
*. Transform gradient to diagonal basis 
       CALL MATVCC(WORK(KLU),WORK(KLE1),WORK(KLE1TR),
     &             LENX,LENX,1)         
*
       IF(NTEST.GE.100) THEN
         WRITE(6,'(A,E25.12)') 
     &   ' Approximate energy at current point ', E0
         WRITE(6,*) ' eigenvalues and gradient '
         DO I = 1, LENX
           WRITE(6,*)  WORK(KLE2_P-1+I), WORK(KLE1TR-1+I)
         END DO
       END IF
*
       THRES = 1.0D-3
*. Use new routine for soloving potentially damped NR equations
       TOLER = 1.1D0
       STEP_MAX = 0.6D0
*. In the Hessian there are a number of singularities,
*. redefine these as a large number 
       XLARGE = 1.0D10
       DO I = 1, LENX
         IF(ABS(WORK(KLE2_P-1+I)).LT.THRES) 
     &          WORK(KLE2_P-1+I) = XLARGE
       END DO
C     SOLVE_SHFT_NR_IN_DIAG_BASIS(
C    &           E1,E2,NDIM,STEP_MAX,TOLERANCE,X,ALPHA,DELTA_F_PRED)
       CALL SOLVE_SHFT_NR_IN_DIAG_BASIS(
     &      WORK(KLE1TR),WORK(KLE2_P),LENX,STEP_MAX,
     &      TOLER,WORK(KLDELTATR),ALPHA_OUT,DELTA_F_PRED)
       XSTEP = SQRT(INPROD(WORK(KLDELTATR),WORK(KLDELTATR),LENX))
       WRITE(6,*) ' Length of step ', XSTEP
*. Transform to orginal basis 
       CALL MATVCC(WORK(KLU),WORK(KLDELTATR),WORK(KLDELTA),
     &             LENX,LENX,0)         
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Step in original basis '
         CALL WRTMAT(WORK(KLDELTA),1,LENX,1,LENX)          
       END IF
*. Update coefficients
       ONE = 1.0D0
       CALL VECSUM(WORK(KLREFX),WORK(KLREFX),WORK(KLDELTA),
     &             ONE,ONE,LENX)
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Updated coefficients '
         CALL WRTMAT(WORK(KLREFX),1,LENX,1,LENX)
       END IF
*. Generate new density 
C          NEW_DENSI_FROM_BCH(DIN,DOUT,X,IPURIFY)
      CALL NEW_DENSI_FROM_BCH(DOPT,WORK(KLMATC),WORK(KLREFX),1)
      CALL COPVEC(WORK(KLMATC),DOPT,LENC)
*
      END DO
*     ^ End of loop over Newton iterations
*
*. Approximate energy of optimized density
C     GET_APPROX_HFEM(E_APPROX,IAPPROX,DENSI,NDENSI)
      CALL GET_APPROX_HFEM(E_APPROX,IAPPROX,DOPT,NDENSI)
      WRITE(6,*) ' Approximate energy of new density ', E_APPROX
*. Find exact energy of optimized D_tilde -- constructs F(D_ort)
C?      CALL GET_APPROX_HFEM(E_EXACT,3,DOPT,NDENSI)
C?      WRITE(6,*) ' Exact energy of D_tilde(F constructed) ', E_EXACT
*
      NTEST = 10
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Optimized density '
        CALL APRBLM2(DOPT,NAOS_ENV,NAOS_ENV,NSMOB,0)
      END IF
** TEST ZONE 
*. :. Obtain Fock matrix for current density 
C?    CALL GET_FOCK(WORK(KLMATC),WORK(KHAO),DOPT,ENERGY,ECORE_ORIG,3)
C?    WRITE(6,*) ' Fock matrix obtained directly from DOP '
C?    CALL WRT_AOMATRIX(WORK(KLMATC),1)
*. Obtain energy for current energy
C?      CALL MEMMAN(KLESCR2  ,LENC ,'ADDL ',2,'ESCR2 ')
C?      CALL MEMMAN(KLESCR3  ,LENC ,'ADDL ',2,'ESCR3 ')
C?      ENERGY = EHF_FROM_HFD(WORK(KHAO),WORK(KLMATC),DOPT,
C?   &           WORK(KLESCR2),WORK(KLESCR3)) + ECORE_ORIG
C?      WRITE(6,*) ' Energy obtained in TEST zone ', ENERGY
*. END OF TEST ZONE
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'OP_EXS')
*
      RETURN
      END
      SUBROUTINE PROJ_NONRED(X,XPROJ,D,ITRNSP)
*
* Project redundant part of matrix X out :
*
* ITRNSP = 0 : 
* X = P X Q(T) + Q X P(T) = DS X (1 - SD) + (1 - DS) X SD
*                         = DSX + XSD - 2 DSXSD
* ITRNSP = 1 : 
* X = P(T) X Q) + Q(T) X P = SD X (1 - DS) + (1 - SD) X DS
*                         = SDX + XDS - 2 SDXDS
*
* where P = DS, where D is the density matrix scaled to 
* idempotent form.
*
* Input D is assumed to be standard orbital matrix in complete form, so 
* scaling is performed here 
*
*. Jeppe Olsen, October 2003
*               ITRNSP added, May 2005
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
*. Input
      DIMENSION X(*), D(*)
*. Output
      DIMENSION XPROJ(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'PROJNR')
*. Space for scratch matrices 
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLSCR1,LENC,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LENC,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LENC,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSAO_E,LENC,'ADDL  ',2,'SAO_E ')
*. Expand S to complete form 
      CALL TRIPAK_BLKM(WORK(KLSAO_E),WORK(KSAO),2,NAOS_ENV,NSMOB)
      IF(ITRNSP.EQ.0) THEN
*. Obtain DS in SCR1
        CALL MULT_BLOC_MAT(WORK(KLSCR1),D,WORK(KLSAO_E),
     &       NSMOB,NAOS_ENV,NAOS_ENV,
     &       NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,0)
      ELSE
*. Obtain SD in SCR1
        CALL MULT_BLOC_MAT(WORK(KLSCR1),WORK(KLSAO_E),D,
     &       NSMOB,NAOS_ENV,NAOS_ENV,
     &       NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,0)
      END IF
*. Remember the factor 0.5 that should be multiplied with D
      HALF = 0.5D0
      CALL SCALVE(WORK(KLSCR1),HALF,LENC)
*. DSX/SDX in SCR2
      CALL MULT_BLOC_MAT(WORK(KLSCR2),WORK(KLSCR1),X,
     &     NSMOB,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,0)
*. XSD/XDS in SCR3
      CALL MULT_BLOC_MAT(WORK(KLSCR3),X,WORK(KLSCR1),
     &     NSMOB,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,2)
*. DSX + XSD /SDX + XDS in XPROJ
      ONE = 1.0D0
      CALL VECSUM(XPROJ,WORK(KLSCR2),WORK(KLSCR3),ONE,ONE,LENC)
C?    WRITE(6,*) ' DSX + XSD/SDX+XDS in PROJ_NR'
C?    CALL APRBLM2(XPROJ,NAOS_ENV,NAOS_ENV,NSMOB,0)
*. DSXSD/SDXDS in SCR2
      CALL MULT_BLOC_MAT(WORK(KLSCR2),WORK(KLSCR1),WORK(KLSCR3),
     &     NSMOB,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,0)
       TWOM = -2.0D0
       CALL VECSUM(XPROJ,XPROJ,WORK(KLSCR2),ONE,TWOM,LENC)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(ITRNSP.EQ.0) THEN
          WRITE(6,*) ' Matrix and nonredundant part of matrix '
        ELSE 
          WRITE(6,*) ' Matrix and P(script)(T) times matrix '
        END IF
        CALL APRBLM2(X    ,NAOS_ENV,NAOS_ENV,NSMOB,0)      
        CALL APRBLM2(XPROJ,NAOS_ENV,NAOS_ENV,NSMOB,0)      
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'PROJNR')
*
      RETURN
      END
      SUBROUTINE OPTIM_SCF_USING_ONE_STEP(MAXIT_HF,
     &  E_HF,CONVER_HF,E1_FINAL,NIT_HF,E_ITER)       
C       CALL OPTIM_SCF_USING_ONE_STEP(MAXIT_HF,
C    &       E_HF,CONVER_HF,E1_FINAL,NIT_HF,WORK(KEITER))
*
* Optimize SCF wavefunction using one-step method
*
* Current version assumes closed shell states
*
*. On input an initial MO-AO is assumed residing in work(kccur)
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INNER_PRODUCT_MAT
      LOGICAL CONVER_HF
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cecore.inc'
*
      DIMENSION E_ITER(*)
*
      CALL LUCIAQENTER('OPTHF')
*
      IDUMMY = 0
      CONVER_HF = .FALSE.
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ', IDUMMY,'OPTIM ')
*
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
*
      CALL MEMMAN(KLDNEW,LENC,'ADDL  ',2,'DNEW  ')
*
      LEN = MAX(LENC,MAXIT_HF**2)
      CALL MEMMAN(KLSCR1,LEN,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LEN,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LEN,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSCRV1,MAXIT_HF,'ADDL  ',2,'SCRV1 ')
      CALL MEMMAN(KLC_PROJ,MAXIT_HF,'ADDL  ',2,'C_PROJ')
      CALL MEMMAN(KLSDD,LEN,'ADDL  ',2,'SDD   ')
      CALL MEMMAN(KLSDDINV,LEN,'ADDL  ',2,'SDDINV')
      CALL MEMMAN(KLSAO_E,LEN,'ADDL  ',2,'SAO_E  ')
      CALL MEMMAN(KLD_PROJ,LEN,'ADDL  ',2,'D_PROJ')
      CALL MEMMAN(KLD_EXP,MAXIT_HF**2,'ADDL  ',2,'D_EXP ')
      CALL MEMMAN(KLDORTN,MAXIT_HF,'ADDL  ',2,'D_ORTN')
      CALL MEMMAN(KLDCHAN,MAXIT_HF,'ADDL  ',2,'D_CHAN')
*
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLD_EXP),ZERO,MAXIT_HF**2)
*
      NTEST= 000
      THRES_E = 1.0D-10
*. Obtain input density from input orbitals
        CALL GET_P_FROM_C(WORK(KCCUR),WORK(KDAO),
     &       NHFD_STASYM,NAOS_ENV,NSMOB)
*. S matrix in expanded form
        CALL TRIPAK_BLKM(WORK(KLSAO_E),WORK(KSAO),2,NAOS_ENV,NSMOB)
      NFOCK = 0
      WRITE(6,*) ' Largest allowed number of iterations', MAXIT_HF
      DO ITER = 1, MAXIT_HF
        WRITE(6,*)
        WRITE(6,*) ' =========================================='
        WRITE(6,*) ' Information from iteration ', ITER
        WRITE(6,*) ' =========================================='
        WRITE(6,*)
        NIT_HF = ITER
*. For numerical stability, the Fock matrices are 
*. constructed as incremental matrices
*. Obtain projection of current density on previous densities,
*. and the part of the density that cannot be projected
C     PROJECT_DENSI_ON_DENSI(DENSIIN,DENSI,LDENSI,NDENSI,
C    &           SAO,SINV,PROJ_DENSI,PROJ_COEF,NSMOB,NOBPSM,
C    &           SCR1,SCR2,SCRVEC,I_DO_SINV,DORIG_NORM, DPROJ_NORM)
        CALL PROJECT_DENSI_ON_DENSI(WORK(KDAO),WORK(KDAO_COLLECT),LENC,
     &       NFOCK,WORK(KLSAO_E),WORK(KLSDDINV),WORK(KLD_PROJ),
     &       WORK(KLC_PROJ),NSMOB,NAOS_ENV,WORK(KLSCR1),WORK(KLSCR2),
     &       WORK(KLSCRV1),1,DORIG_NORM,D_PROJ_NORM)
*. Save expansion coefficients for current density in 
* D_EXP dimensioned as D_EXP(MAXIT_HF,MAXIT_HF)
        CALL COPVEC(WORK(KLC_PROJ),WORK(KLD_EXP+NFOCK*MAXIT_HF),NFOCK)
        WORK(KLD_EXP+NFOCK*MAXIT_HF+NFOCK) = 1.0D0
*. Part of density that could not be projected 
        ONE = 1.0D0
        ONEM = -1.0D0
        CALL VECSUM(WORK(KDAO_COLLECT+NFOCK*LENC),WORK(KDAO),
     &              WORK(KLD_PROJ),ONE,ONEM,LENC)
*. Norm of part of new density that could not be projected on 
*. Previous densities 
      KLDAO_EFF = KDAO_COLLECT+NFOCK*LENC
      DORTN = INNER_PRODUCT_MAT(WORK(KLDAO_EFF),WORK(KLDAO_EFF),
     &       WORK(KLSAO_E),
     &       WORK(KLSCR1),WORK(KLSCR2),NSMOB,NAOS_ENV)
      IF(NFOCK.GT.0) WORK(KLDORTN+NFOCK-1) = SQRT(DORTN)
*. Norm of part of 
*. Construct Fock matrix for current density-difference
        CALL GET_FOCK(WORK(KFAO_COLLECT+NFOCK*LENP),WORK(KHAO),
     &       WORK(KDAO_COLLECT+NFOCK*LENC),ENERGYX,
     &       ECORE_ORIG,3)
*. Obtain Fock matrix for current total density
C     GET_F_FOR_SUM_DENSI(FOUT,COEF,NDENSI)
        CALL GET_F_FOR_SUM_DENSI(WORK(KFAO),
     &       WORK(KLD_EXP+NFOCK*MAXIT_HF),NFOCK+1)
*
C?      WRITE(6,*) ' Current density and Fock matrix '
C?      CALL WRT_AOMATRIX(WORK(KDAO),0)
C?      CALL WRT_AOMATRIX(WORK(KFAO),1)
*. Obtain energy for current energy
        ENERGY = EHF_FROM_HFD(WORK(KHAO),WORK(KFAO),WORK(KDAO),
     &           WORK(KLSCR1),WORK(KLSCR2)) + ECORE_ORIG
        E_ITER(ITER)=ENERGY 
        E_HF = ENERGY
*. Well a new fock matrix has been constructed and saved, 
*. it is now time to register it
        NFOCK = NFOCK + 1
        WRITE(6,'(A,I3, E25.15)') ' Iteration, current energy = ', 
     &                          ITER, ENERGY
*. study overlap of densities
C       GET_OVERLAPMAT_AO_DENS(DENSI,LDENSI,NDENSI,SAO,
C    &           SOVERLAP,SCR1,SCR2,NSMOB,NOBPSM)
        CALL GET_OVERLAPMAT_AO_DENS(WORK(KDAO_COLLECT),LENC,
     &       NFOCK,WORK(KLSAO_E),WORK(KLSDD),WORK(KLSCR1),
     &       WORK(KLSCR2),NSMOB,NAOS_ENV)
        WRITE(6,*) ' The <D_i!D_j> matrix '
        CALL WRTMAT(WORK(KLSDD),NFOCK,NFOCK,NFOCK,NFOCK)
        CALL DIAG_SYM_MAT(WORK(KLSDD),WORK(KLSCR1),WORK(KLSCR2),
     &  NFOCK,0)
        WRITE(6,*) ' Eigenvalues of <D_i!D_j> '
        CALL WRTMAT(WORK(KLSCR2),1,NFOCK,1,NFOCK)

*. Test for convergence 
        IF(ITER.GT.1) THEN
          DELTAE=DABS(E_ITER(ITER)-E_ITER(ITER-1))
          WRITE(6,*) ' Change of energy ', DELTAE
        END IF
*
        IF(ITER.GT.1.AND.DELTAE.LT.THRES_E)THEN 
          WRITE(6,*) ' Energy converged, last energy change ',
     &    DELTAE
          CONVER_HF = .TRUE.
          GOTO 1111
        END IF
*. Obtain improved density by optimizing energy function 
*. Defined by current and previous densities and Fock matrices
        CALL OPTIM_EHFAPR_EXPSX(WORK(KLDNEW),NFOCK,WORK(KDAO))
*. Find norm of change of density
        ONE  =  1.0D0
        ONEM = -1.0D0
        CALL VECSUM(WORK(KLSCR3),WORK(KLDNEW),WORK(KDAO),
     &       ONE,ONEM,LENC)
        DCHAN = INNER_PRODUCT_MAT(WORK(KLSCR3),WORK(KLSCR3),
     &       WORK(KLSAO_E),
     &       WORK(KLSCR1),WORK(KLSCR2),NSMOB,NAOS_ENV)
        WORK(KLDCHAN-1+NFOCK) = SQRT(DCHAN)
*. And enroll new density 
        CALL COPVEC(WORK(KLDNEW),WORK(KDAO),LENC)
      END DO
1111  CONTINUE 
*. Obtain orbitals corresponding to final density (KLDNEW) and sort into 
*. occupied and virtual orbitals, save in KCCUR
C HER ER JEG
C     DIAG_SDS(D,S,C,XOCCNUM)
C     CALL DIAG_SDS(WORK(KLDNEW),WORK(KSAO),WORK(KCCUR),WORK(KLSCR1))
      CALL DIAG_SDS_GEN(WORK(KLDNEW),WORK(KSAO),WORK(KCCUR),
     &     WORK(KLSCR1))
*
      WRITE(6,*) ' Number of Fock matrix evaluations ', NFOCK
*
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) ' Final energy in au  ', ENERGY
      WRITE(6,*) ' Overlap between final and earlier densities : '
*
      WRITE(6,*)
      WRITE(6,*) ' =============================================='
      WRITE(6,*)
      CALL DOVLAP(WORK(KLSCRV1),WORK(KDAO),WORK(KDAO_COLLECT),
     &            NFOCK,1,WORK(KLD_EXP),MAXIT_HF)
C     DOVLAP(VOVLP,DTARGET,DOTHER,NDENSI,ITRANS,CTRANS,
C    &                  NDIMC)
*. Find norm of differences <DFinal - D_i!DFinal - D_i> =
*                           2<Dfinal!DFinal> - 2<D_i!DFinal>
* (assuming all densities have the same trace )
C     INNER_PRODUCT_MAT(A,B,S,SCR1,SCR2,NBLK,LBLK)
      DFDF = INNER_PRODUCT_MAT(WORK(KDAO),WORK(KDAO),WORK(KLSAO_E),
     &       WORK(KLSCR1),WORK(KLSCR2),NSMOB,NAOS_ENV)
      WRITE(6,*) ' DFDF = ', DFDF
      DO I = 1, NFOCK-1
        WORK(KLSCRV1-1+I) = 
     &  SQRT(MAX(0.0D0,2.0D0*DFDF-2.0D0*WORK(KLSCRV1-1+I)))
      END DO
*
      WRITE(6,*) ' Norms of Dfinal - D_i '
      WRITE(6,*) ' ======================'
      CALL WRTMAT_EP(WORK(KLSCRV1),1,NFOCK-1,1,NFOCK-1)
*
      WRITE(6,*) ' Norm of change of density'
      CALL WRTMAT_EP(WORK(KLDCHAN),1,NFOCK-1,1,NFOCK-1)
*
      WRITE(6,*) ' Norm of orthogonal part of change of density'
      CALL WRTMAT_EP(WORK(KLDORTN),1,NFOCK-1,1,NFOCK-1)
*
      WRITE(6,*) ' Summary of iteration sequence '
      WRITE(6,*)  
     &' i  E(i+1) - E(i)  E(i) - Econv |D(i+1) - D(i)| |D(i+1)-D(i)| P',
     &  ' |D(i) - Dconv|'
      WRITE(6,*) 
     & '==============================================================',
     & '================'
      DO I = 1, NFOCK-1
          WRITE(6,'(I3,5(2X,E13.6))') I,
     &    WORK(KEITER+I)-WORK(KEITER+I-1),
     &    WORK(KEITER-1+I)-ENERGY,
     &    WORK(KLDCHAN-1+I),
     &    WORK(KLDCHAN-1+I)-WORK(KLDORTN-1+I),
     &    WORK(KLSCRV1-1+I)
      END DO
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ', IDUMMY,'OPTIM ')
      CALL LUCIAQEXIT('OPTHF')

      RETURN
      END
      FUNCTION EHF_FROM_HFD(HAO,FAO,DAO,SCR1,SCR2)
*
* Obtain energy as 1/2 Tr (H + F) D
*
*. Jeppe Olsen, October 2003
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
* H and F in compact form, D in expanded form
      DIMENSION HAO(*), FAO(*), DAO(*)
*. Two scratch matrices, length of full matrices
      DIMENSION SCR1(*), SCR2(*)
*
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
*. 0.5*( H + F)
      HALF = 0.5D0
      CALL VECSUM(SCR1,HAO,FAO,HALF,HALF,LENP)
C?    WRITE(6,*) ' 0.5*(HAO + FAO) in ...HFD'
C?    CALL WRT_AOMATRIX(SCR1,1)
* Expand H + F to complete form in SCR2
      CALL TRIPAK_BLKM(SCR2,SCR1,2,NAOS_ENV,NSMOB)
*. Product 0.5(H+F)D in SCR1
      CALL MULT_BLOC_MAT(SCR1,SCR2,DAO,NSMOB,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,0)
*. and Trace fD
C                 TRACE_BLK_MAT(A,NBLK,LBLK,IPAK)
      E = TRACE_BLK_MAT(SCR1,NSMOB,NAOS_ENV,0)
      EHF_FROM_HFD = E
*
      RETURN
      END 
      SUBROUTINE SOLVE_SHFT_NR_IN_DIAG_BASIS(
     &           E1,E2,NDIM,STEP_MAX,TOLERANCE,X,ALPHA,DELTA_F_PRED)
*
* A gradient and Hessian is given in a diagonal basis
*
* Solve shifted Newton-Raphson equations :
*
* (E2+ alpha*1) X = -E1
*
* so the norm of X is atmost STEP (with a relative tolerance TOLERANCE)
*
* The shifted Hessian, E2 + alpha 1 is required to be positive definite
*
* The step is requrired to be in the range 
* STEP_MAX/TOLERANCE to STEP_MAX*TOLERANCE
* so the closer TOLERANCE (>1)  is to 1.0D0, the smaller is the 
* allowed tolerance
*
* The predicted change of the function value is returned in DELTA_F_PRED
*
*. Jeppe Olsen, Oct 2003
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
*. Input
      DIMENSION E1(NDIM),E2(NDIM)
*. Output 
      DIMENSION X(NDIM)
*
      NTEST = 00
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Output from optimizer with step-control '
      END IF
*
*. Smallest ACTIVE eigenvalue
*. ACTIVE ?? Well, there can be eigenvalue where the 
*. corresponding gradient vanishes, f.ex. due to symmetry 
*. reasons.
*. 
*. Norm of gradient 
      E1_NORM = SQRT(INPROD(E1,E1,NDIM))
*. Threshold for defining gradient to be vanishing
      ZERO_EFF = E1_NORM*1.0D-7
*. It is assumed that the eigenvalues are sorted 
      EIG_MIN = 1.0D0
      DO I = 1, NDIM
        IF(ABS(E1(I)).GT.ZERO_EFF) THEN
          EIG_MIN = E2(I)
          GOTO 1111
        END IF
      END DO
 1111 CONTINUE
 
C     EIG_MIN = XMNMX(E2,NDIM,1)
*. If Hessian is positive definite, try first unshifted solution
      IF(EIG_MIN.GT.0.0D0) THEN
        ZERO = 0.0D0
        CALL DIAVC2(X,E1,E2,ZERO,NDIM)
        ONEM = -1.0D0
        CALL SCALVE(X,ONEM,NDIM)
        XNORM = SQRT(INPROD(X,X,NDIM))
        ALPHA_0 = 0.0D0
        ALPHA = 0.0D0
        IF(XNORM.GT.STEP_MAX*TOLERANCE) THEN
          I_MUST_SHIFT = 1
        ELSE 
          I_MUST_SHIFT = 0
          IF(NTEST.GE.5) THEN
           WRITE(6,*) ' Unshifted step, norm of gradient and step ',
     &     E1_NORM, XNORM
          END IF
        END IF
      ELSE
*. There are negative eigenvalues, so shift is neccesary
        I_MUST_SHIFT = 1
        ALPHA_0 = - EIG_MIN
      END IF
*
      IF(I_MUST_SHIFT.EQ.1) THEN
*. Choose initial shift assuming that all eigenvalues equals the 
*. smallest - this will lead to an underestimate of the shift 
       ALPHA = E1_NORM/STEP_MAX - EIG_MIN - ALPHA_0
       MAXIT = 6
       DO IT = 1, MAXIT
         IF(IT.GT.1) THEN
*.  We model the step function as xk/(eig_min + alpha_0 + alpha). 
*.  We have a step, so use this to determine xk, and then 
*. determine a new step
           XK = XNORM*(EIG_MIN + ALPHA_0 + ALPHA )
           ALPHA = XK/STEP_MAX - EIG_MIN - ALPHA_0
         END IF 
*. New step 
         ALPHA_TOT = ALPHA_0 + ALPHA
         CALL DIAVC2(X,E1,E2,ALPHA_TOT,NDIM)
         ONEM = -1.0D0
         CALL SCALVE(X,ONEM,NDIM)
         XNORM = SQRT(INPROD(X,X,NDIM))
         IF(XNORM.LT.STEP_MAX*TOLERANCE .AND. 
     &      XNORM.GT.STEP_MAX/TOLERANCE      ) THEN
           ICONVER = 1
         ELSE 
           ICONVER = 0
         END IF
         IF(NTEST.GE.5) THEN
           WRITE(6,'(A, 2(2X,E8.3))') 
     &     ' Shift and norm of step', ALPHA,XNORM
         END IF
C?       WRITE(6,*) ' IT, XK, ALPHA, XNORM = ', 
C?   &                IT, XK, ALPHA, XNORM 
         IF(ICONVER.EQ.1) GOTO 1001
       END DO
*      ^ end of loop over iterations
      END IF
*     ^ End if shift was neccessary
 1001 CONTINUE
*
*. Predicted change of function value
*
      E1TERM = INPROD(X,E1,NDIM)
      E2TERM = 0.0D0
      DO I = 1, NDIM
        E2TERM = E2TERM + X(I)*E2(I)*X(I)
      END DO
      E2TERM = 0.5D0*E2TERM
      DELTA_F_PRED= E1TERM + E2TERM
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SOLVE_NR ... '
        WRITE(6,*) ' Eigenvalue, gradient, step '
        DO I = 1, NDIM
         WRITE(6,'(I4,3E18.10)') I, E2(I),E1(I),X(I)
        END DO
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,'(A,E18.10)')
     &  ' Predicted change of function value ', DELTA_F_PRED 
      END IF
*
      RETURN
      END
      FUNCTION XMNMX(VEC,NDIM,MINMAX)
*
* Find smallest (MINMAX = 1), largest (MINMAX = 2), largest absolute (MINMAX = 3)
* value of elements in VEC. 
*
*. Jeppe Olsen, Oct. 2003
*
* Last revision, Sept. 3 2012, Jeppe Olsen, MINMAX = 3 added
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VEC(*)
*
      IF(MINMAX.LE.2) THEN
        VALUE = VEC(1)
      ELSE
        VALUE = ABS(VEC(1))
      END IF
*
      IF(MINMAX.EQ.1) THEN
        DO I = 2, NDIM
          VALUE = MIN(VALUE,VEC(I))
        END DO
      ELSE IF(MINMAX.EQ.2) THEN
        DO I = 2, NDIM
          VALUE = MAX(VALUE,VEC(I))
        END DO
      ELSE IF(MINMAX.EQ.3) THEN
        DO I = 2, NDIM
          VALUE = MAX(VALUE,ABS(VEC(I)))
        END DO
      ELSE
        WRITE(6,*) ' Unknown parameter MINMAX in XMNMX = ', MINMAX
        STOP ' Unknown parameter MINMAX in XMNMX ' 
      END IF
*
      XMNMX = VALUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(MINMAX.EQ.1) THEN
         WRITE(6,*) ' Smallest element in vector ', VALUE
        ELSE 
         WRITE(6,*) ' Largest element in vector ', VALUE
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE SELECT_RH_SHIFT(W,WMIN,E,NOBPSM,NOCPSM,NSMOB,SHIFT)
*
* Select shift in RH iteration to ensure 
* a given minimal overlap  between old and new occ orbitals
*
* A shift can contribute to stabilization in two ways
*  1 : Complete interchange of occupied orbitals, by 
*      shifting the orbital energies of the previous occ orbs.
*      down.
*  2 : Damping the change without shifting order of orbitals.
*
* It is not obvious how to distinguish between the two cases, and 
* currently I am only taking care of situation 1
*
* It could be possible to use the weight of a given orb to
* distinguish between the 2 cases
*
* Jeppe Olsen, Sept. 2004
*
      INCLUDE 'implicit.inc'
*. Input : Orbital energies and weight of old occ in new orbitals
      DIMENSION E(*), W(*)
*. Number of orbitals per symmetry and number of occ per symmetry
      DIMENSION NOBPSM(NSMOB),NOCPSM(NSMOB)
*
      SHIFT = 0.0D0
      IOFF = -137
      DO ISMOB = 1, NSMOB
        IF(ISMOB.EQ.1) THEN
          IOFF = 1
        ELSE 
          IOFF = IOFF + NOBPSM(ISMOB-1)
        END IF
        SHIFTL = 0.0D0
        DO IOC = 1, NOCPSM(ISMOB)
          IF(W(IOFF-1+IOC).LT.WMIN.AND.SHIFTL.EQ.0.0D0) THEN
*. shift so this orbital becomes above lowest unoccupied of this 
*. symmetry
            EPSIL = E(IOFF-1+IOC)
            ELUMO = E(IOFF-1+NOCPSM(ISMOB)+1)
            SHIFTL = -(ELUMO-EPSIL) - 0.001D0
          END IF
        END DO
        SHIFT = MIN(SHIFTL,SHIFT)
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Suggested shift ', SHIFT 
      END IF
*
      RETURN
      END 
      SUBROUTINE DBAR_CHANGE_METRIC(DENSI,NDENSI,LDENS,IREFDENS,SDBAR)
*
* Obtain metric for changes of DBAR
*
* ||Dbar-D_irefdens!!^2 = sum_ij C_I S_IJ C_J
*
* Jeppe Olsen, Sept. 2004
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      REAL*8 INNER_PRODUCT_MAT
*
*. Input in symmetry-blocked form, complete unpacked blocks
      DIMENSION DENSI(LDENS,NDENSI)
*. Output
      DIMENSION SDBAR(NDENSI-1,NDENSI-1)
*
C?    WRITE(6,*) ' DBAR_CHANGE.... : NDENSI, IREFDENS = ',
C?   &             NDENSI,IREFDENS
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DBAR_C')
*. Space for expanded AO metric and two scratch matrices
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLSAO_E,LENC,'ADDL  ',2,'SAO_E ')
      CALL MEMMAN(KLSCR1 ,LENC,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2 ,LENC,'ADDL  ',2,'SCR2  ')
*. Obtain S in expanded form 
      CALL TRIPAK_BLKM(WORK(KLSAO_E),WORK(KSAO),2,NAOS_ENV,NSMOB)
*
      ONE = 1.0D0
      ONEM = -1.0D0
*. Modify all densities(different from REFDENS) to D(I)-D_(REFDENS)
      DO I = 1, NDENSI
        IF(I.NE.IREFDENS) THEN
          CALL VECSUM(DENSI(1,I),DENSI(1,I),DENSI(1,IREFDENS),
     &                ONE,ONEM,LDENS)
        END IF
      END DO
*
      I = 0
      DO II = 1, NDENSI
        IF(II.NE.IREFDENS) THEN
          I = I + 1
          J = 0
          DO JJ = 1, II
            IF(JJ.NE.IREFDENS) THEN
              J = J + 1
*. And product
C?            WRITE(6,*) ' II, JJ, I,J = ', II,JJ,I,J
              SDBAR(I,J) = INNER_PRODUCT_MAT(DENSI(1,II),DENSI(1,JJ),
     &                     WORK(KLSAO_E),WORK(KLSCR1),WORK(KLSCR2),
     &                     NSMOB,NAOS_ENV)
              SDBAR(J,I) = SDBAR(I,J)
            END IF
          END DO 
        END IF
      END DO
*. Clean up by get good old densities back
      DO I = 1, NDENSI
        IF(I.NE.IREFDENS) THEN
          CALL VECSUM(DENSI(1,I),DENSI(1,I),DENSI(1,IREFDENS),
     &                ONE,ONE,LDENS)
        END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Metric for changes of Dbar '
        CALL WRTMAT(SDBAR,NDENSI-1,NDENSI-1,NDENSI-1,NDENSI-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DBAR_C')
      RETURN
      END 
      SUBROUTINE SOLVE_NR_WITH_GEN_TRCTL(E1,E2,S,STEP,NVAR,XMAXNORM)
*
*. Solve Newton-Raphson equations with general trust region control, 
* i.e. the constraint x^T S x  .le. XMAXNORM**2
*
*. Jeppe Olsen, Sept. 2004
*
*. Small matrix version, assuming diagonalizations are no problem ..
*
*. Presently we are assuming minimization, so the Hessian is shifted 
*. to become positive definite
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*. Input
      DIMENSION E1(NVAR),E2(NVAR,NVAR),S(NVAR,NVAR)
*. Output
      DIMENSION STEP(NVAR)
*
      NTEST = 10
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'NR_TRC')
*
*. Memory allocations
*. Local copies of S and E2, as these are destroyed during 
*. diagonalization 
      CALL MEMMAN(KLE2,NVAR**2,'ADDL  ',2,'L_E2  ')
      CALL MEMMAN(KLS ,NVAR**2,'ADDL  ',2,'L_S   ')
*. eigenvectors 
      CALL MEMMAN(KLX ,NVAR**2,'ADDL  ',2,'L_X   ')
*. Eigenvalues 
      CALL MEMMAN(KLEIGVAL ,NVAR,'ADDL  ',2,'EIGVAL')
*. And a scratch matrix
      CALL MEMMAN(KLSCRM ,NVAR**2,'ADDL  ',2,'SCRMAT')
*. Gradient in eigenvector basis 
      CALL MEMMAN(KLE1T,NVAR,'ADDL  ',2,'E1T   ')
*. Step in eigenvector basis 
      CALL MEMMAN(KLSTEPT,NVAR,'ADDL  ',2,'STEPT ')
*
* ===================================================
*. Solve general eigenvalueproblem  E2 X = S X Eigval
* ===================================================
*
C     GENDIA_BLMAT(HIN,SIN,C,E,PVEC,NBLK,LBLK,ISORT)
      CALL COPVEC(E2,WORK(KLE2),NVAR**2)
      CALL COPVEC(S ,WORK(KLS ),NVAR**2)
      CALL GENDIA_BLMAT(WORK(KLE2),WORK(KLS),WORK(KLX),WORK(KLEIGVAL),
     &                  WORK(KLSCRM),NVAR,1,1)
*. Transform gradient to eigenvector basis 
      CALL MATVCC(WORK(KLX),E1,WORK(KLE1T),NVAR,NVAR,1)
      E1NORM = SQRT(INPROD(WORK(KLE1T),WORK(KLE1T),NVAR))
*. Negative ACTIVE eigenvalues, use that the eigenvalues are sorted, so 
*. the lowest comes first 
*. Define an effective zero. THis is not without complications
*. As convergence is approached the gradient becomes zero and 
* the distinction between truly vanishing gradients and 
* converged may become blurred.
*
      ZERO_EFF = MAX(E1NORM*1.0D-12,1.0D-10)
      EIGVMIN = 1.0D0
      E1MIN = 666666.0D0
      WRITE(6,*) ' Lowest eigenvalue = ', WORK(KLEIGVAL)
      DO I = 1, NVAR
        IF(ABS(WORK(KLE1T-1+I)).GT.ZERO_EFF) THEN
          EIGVMIN = WORK(KLEIGVAL-1+I)
*. gradient of smallest ACTIVE eigenvalue
          E1MIN = WORK(KLE1T-1+I)
          WRITE(6,*) 'Lowest eigenvalue with nonvanishing gradient',
     &    EIGVMIN
C?        WRITE(6,*) ' ZERO_EFF = ', ZERO_EFF
          GOTO 1111
        END IF
      END DO
 1111 CONTINUE
C     EIGVMIN = WORK(KLEIGVAL)
      THRES = 1.0D-10
      IF(EIGVMIN.LT.-THRES) THEN
        INEG = 1
        ISHIFT = 1
      ELSE
        INEG = 0
        ISHIFT = 0
      END IF
*
      IF(INEG.EQ.0) THEN
        SHIFT = 0.0D0
*. Solve unshifted equations in eigen basis 
        CALL NR_STEP_IN_EIGBASIS(WORK(KLE1T),WORK(KLEIGVAL),SHIFT,
     &       WORK(KLSTEPT),NVAR)
*. Norm of step  - easy as we are in basis with unit metric
        STEP_NORM = SQRT(INPROD(WORK(KLSTEPT),WORK(KLSTEPT),NVAR))
        IF(STEP_NORM.LT.XMAXNORM) THEN
          ISHIFT = 0
          SHIFT = 0.0D0
        ELSE 
          ISHIFT = 1
        END IF
      END IF
*
      IF(ISHIFT.EQ.1) THEN
*. Iterate to obtain shift that gives allowed steplength
*. we can put lower and upper bounds on the shift as we know 
*. the lowest eigenvalue. The lower bound is obtained by assuming
*. that only the lowest eigenvalue contribute. The upper bound is 
*. obtained by assuming all eigenvalues are identical to the lowest
*
        SHIFT_MAX = -ABS(E1MIN)/XMAXNORM + EIGVMIN
        E1NRM = SQRT(INPROD(WORK(KLE1T),WORK(KLE1T),NVAR))
        SHIFT_MIN = -E1NRM/XMAXNORM + EIGVMIN
        WRITE(6,*) ' E1NRM, XMAXNORM = ', E1NRM, XMAXNORM
*
        IF(NTEST.GE.5) THEN
          WRITE(6,'(A,2E15.8)') ' Lower and upper bound on shift ', 
     &    SHIFT_MIN,SHIFT_MAX
          WRITE(6,*) '(shifts are subtracted) '
        END IF
*. Iterate over shifts and obtain converged value
        MAXIT = 10
        DO IT = 1, MAXIT
*. Take simple average
          SHIFT = (SHIFT_MIN + SHIFT_MAX)/2.0D0
*. And solve NR equations with this shift
          CALL NR_STEP_IN_EIGBASIS(WORK(KLE1T),WORK(KLEIGVAL),SHIFT,
     &         WORK(KLSTEPT),NVAR)
          STEP_NORM = SQRT(INPROD(WORK(KLSTEPT),WORK(KLSTEPT),NVAR))
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' SHIFT and STEP_NORM ', SHIFT, STEP_NORM
          END IF
          IF(STEP_NORM.LT.XMAXNORM*0.9) THEN
*. Step was to small, so current shift is to negative 
            SHIFT_MIN = SHIFT
          ELSE IF (STEP_NORM.GT.XMAXNORM/0.9) THEN
*. Stqp is too large, so current shift is not sufficient negative
            SHIFT_MAX = SHIFT
          ELSE 
*. Step is within bound, i.e. converged so 
            GOTO 1001
          END IF
        END DO
 1001   CONTINUE
*.      ^ End of loop over iterations for determing shift
      END IF
*     ^ End if shift was required
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Final shift and stepnorm ',
     &  SHIFT,STEP_NORM
      END IF
*. The step was obtained in eigenvector basis, transform to 
*. original basis 
      CALL MATVCC(WORK(KLX),WORK(KLSTEPT),STEP,NVAR,NVAR,0)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Solution to TR controlled NR equations '
        CALL WRTMAT(STEP,1,NVAR,1,NVAR)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'NR_TRC')
      RETURN
      END 

      SUBROUTINE NR_STEP_IN_EIGBASIS(E1,EIGVAL,SHIFT,STEP,NVAR)
*
*. Solve shifted NR problem in eigenvalue basis 
*
* Step(i) = -E1(i)/ (eigval(i) - shift )
*
* Modes corresponding to singularities ( eigval(i)=0) are inactive
*
* Jeppe Olsen, Sept. 2004
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION E1(NVAR),EIGVAL(NVAR)
*. Output
      DIMENSION STEP(NVAR)
*. Threshold for singularities 
      THRES = 1.0D-10
      DO I = 1, NVAR
        IF(ABS(EIGVAL(I)).GT.THRES) THEN
          STEP(I) = -E1(I)/(EIGVAL(I)-SHIFT)
        ELSE
          STEP(I) = 0.0D0
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Solution of NR equation in diagonal basis : '
        WRITE(6,*) ' Step, Eigenvalue and gradient '
        WRITE(6,*) ' ==============================='
        DO I = 1, NVAR
          WRITE(6,*) STEP(I),EIGVAL(I),E1(I)
        END DO
      END IF
*
      RETURN
      END
      FUNCTION V1T_MAT_V2(V1,V2,XMAT,NVAR)
*
* Obtain inner product V1T*XMAT*V2 = Sum(IJ) V1(I) XMAT(I,J) V2(J)
*
*. Jeppe Olsen, Sept. 2004
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION V1(NVAR),V2(NVAR),XMAT(NVAR,NVAR)
*
      SUM = 0.0D0
      DO J = 1, NVAR
        SUMI = 0.0D0
        DO I = 1, NVAR
          SUMI = SUMI + V1(I)*XMAT(I,J)
        END DO
        SUM = SUM + SUMI*V2(J)
      END DO
*
      V1T_MAT_V2 = SUM
*
      RETURN
      END 
      SUBROUTINE SCF_DIA_SHIFT_SCAN
     &(F,S,C,CBEST,E,D,COLD,NOCPSM,NDENSI)
*
* Scan the SCF diagonalization as function of shift parameter
* a number of full energy calculations are performed..
* 
*
* Jeppe Olsen, Sept. 2004
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION F(*),S(*),COLD(*),D(*)
*. Number of occupied orbitals per symmetry
      DIMENSION NOCPSM(*)
*. Output
      DIMENSION C(*),E(*),CBEST(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'SCF_SC')
*. A bit of local scratch
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      CALL MEMMAN(KLSCR,LENC,'ADDL  ',2,'SCR   ')
      CALL MEMMAN(KLCOSCN,LENC,'ADDL  ',2,'COSCN ')
      CALL MEMMAN(KLW,NTOOB,'ADDL  ',2,'WCO_CN')
      CALL MEMMAN(KLSDS,LENC,'ADDL  ',2,'SDS   ')
      CALL MEMMAN(KLS_E,LENC,'ADDL  ',2,'SDS   ')
      CALL MEMMAN(KLP,LENC,'ADDL  ',2,'LP   ')
*
*. Obtain matrix SDS
*
*. Construct SDS : Unpack S, CALC SDS, and pack
*. S(unpack) in KLS_E
      CALL TRIPAK_BLKM(WORK(KLS_E),S,2,NAOS_ENV,NSMOB)
*. DS in KLSCR
      CALL MULT_BLOC_MAT(WORK(KLSCR),D,WORK(KLS_E),
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &     1)
*. SDS in KLSDS
      CALL MULT_BLOC_MAT(WORK(KLSDS),WORK(KLS_E),WORK(KLSCR),
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &     1)
      CALL TRIPAK_BLKM(WORK(KLSDS),WORK(KLSCR),1,NAOS_ENV,NSMOB)
      CALL COPVEC(WORK(KLSCR),WORK(KLSDS),LENP)
           I_HAVE_DONE_SDS = 1
*. Shifts are obtained as 0, STEP*1.25**(N-1)
      STEP = -0.0001
*. Loop over positive and negative shifts '
      WRITE(6,*) '# Shiftinfo, exact and aproximate energy : '
      INI = 1
      DO INP = 1,2
      MXSHFTIT = 46
      IF(INP.EQ.1) THEN
        FACTOR = STEP*(1.25**MXSHFTIT)
      ELSE 
        FACTOR = -STEP/1.25
      END IF
*. Loop over iterations of shifts
      ALPHA = -3006
      DO ISHFTIT = 1, MXSHFTIT
*. Determine new shift
        IF(INI.EQ.1) THEN
          ALPHA_OLD  = 0.0D0
        ELSE 
          ALPHA_OLD = ALPHA
        END IF
        IF(INP.EQ.1) THEN
          FACTOR = FACTOR/1.25D0
        ELSE
          FACTOR = FACTOR*1.25D0
        END IF
        ALPHA = FACTOR
        ALPHA_DEL = ALPHA - ALPHA_OLD
*. F => F + ALPHA_DEL * SDS
        ONE = 1.0D0
        CALL VECSUM(F,F,WORK(KLSDS),ONE,ALPHA_DEL,LENP)
*. Diagonalize
        CALL SCF_DIA(F,S,C,E)
*. Find exact energy for these coefficients
C            GET_P_FROM_C(C,P,NOC,NAO,NSMOB)
        CALL GET_P_FROM_C(C,WORK(KLP),NOCPSM,NAOS_ENV,NSMOB)
        CALL GET_APPROX_HFEM(E_EX,3,WORK(KLP),0)
        CALL GET_APPROX_HFEM(E_APR,2,WORK(KLP),NDENSI)
        IF(INI.EQ.1) THEN
          E_EX_MIN = E_EX
          CALL COPVEC(C,CBEST,LENC)
          SHIFT_MIN = ALPHA
          INI = 0
        ELSE 
          IF(E_EX.LT.E_EX_MIN) THEN
            E_EX_MIN = E_EX
            CALL COPVEC(C,CBEST,LENC)
            SHIFT_MIN = ALPHA
          END IF
        END IF
        WRITE(6,'(A,3E21.10)') 'Shiftinfo ', ALPHA,E_EX,E_APR
      END DO
*     ^ End of loop over shifts
      END DO
*.    ^ End of loop over positive/negative shifts
      WRITE(6,'(A,2E22.12)') 
     & ' Shiftinfo: Min shift and energy ',SHIFT_MIN,E_EX_MIN
*. Clean up : remove SDS term from shifted F
      CALL VECSUM(F,F,WORK(KLSDS),ONE,-ALPHA,LENP)
          
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'SCF_SC')
      RETURN
      END 
      SUBROUTINE DIAG_SDS(D,S,C,XOCCNUM)
*
* Diagonalize Density matrix in AO basis (SDS !!),
* and sort resulting orbitals according to 
* occupation numbers, so the occupied orbitals
* occur first in each symmetry class
*
*. Jeppe Olsen, Oct. 2004
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input : D unpacked, S packed 
      DIMENSION D(*),S(*)
*. Output
      DIMENSION  C(*),XOCCNUM(*)
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIAG_S')
*. Scratch
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      CALL MEMMAN(KLS_E,LENC,'ADDL  ',2,'S_E   ')
      CALL MEMMAN(KLSCR,LENC,'ADDL  ',2,'SCR   ')
      CALL MEMMAN(KLSDS,LENC,'ADDL  ',2,'SCR   ')
*. Construct SDS : Unpack S, CALC SDS, and pack
*. S(unpack) in KLS_E
      CALL TRIPAK_BLKM(WORK(KLS_E),S,2,NAOS_ENV,NSMOB)
*. DS in KLSCR
      CALL MULT_BLOC_MAT(WORK(KLSCR),D,WORK(KLS_E),
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &      1)
*. SDS in KLSDS
      CALL MULT_BLOC_MAT(WORK(KLSDS),WORK(KLS_E),WORK(KLSCR),
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &     1)
CER   CALL TRIPAK_BLKM(WORK(KLSDS),WORK(KLSCR),1,NAOS_ENV,NSMOB)
CER   CALL COPVEC(WORK(KLSCR),WORK(KLSDS),LENP)
*. Multiply with -1 to get largest occupation numbers first
      ONEM = -1.0D0
      CALL SCALVE(WORK(KLSDS),ONEM,LENC)
*. and diagonalize
      ISORT=1
      CALL GENDIA_BLMAT(WORK(KLSDS),WORK(KLS_E),C,XOCCNUM,
     &                  WORK(KLSCR),NAOS_ENV,NSMOB,ISORT)
*. And multiply eigenvalues to get occupation numbers 
*. ( to compensate for the -1 we multiplied by earlier )
      CALL SCALVE(XOCCNUM,ONEM,NTOOB)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Obtained occupation numbers form DIAG_SDS... '
        CALL WRITE_ORBVECTOR(XOCCNUM,NAOS_ENV,NSMOB)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DIAG_S')
*
      RETURN 
      END 
      SUBROUTINE WRITE_ORBVECTOR(ORBVECTOR,NAOS,NSMOB)
*
* Print vector with elements running over orbitals
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION ORBVECTOR(*)
      INTEGER NAOS(NSMOB)
*
*
      IOFF = 1
      DO ISYM = 1, NSMOB
        IF(ISYM.GT.1) THEN
          IOFF = IOFF + NAOS(ISYM-1)
        END IF
        WRITE(6,*) ' Symmetry ', ISYM
        LEN = NAOS(ISYM)
        CALL WRTMAT(ORBVECTOR(IOFF),1,LEN,1,LEN)
      END DO
*
      RETURN
      END
      SUBROUTINE SELECT_DENS_IN_SUBSPC(DENSI,NDENSI,LDENSI,IREFDNS,
     &           NALLOWED,IALLOWED)
*
* A set of densities are given in DENSI. Find which of 
* these may be used in subspace minimization.
* 
      RETURN
      END 
      SUBROUTINE DEL_OV_ORBMAT(A,NSMOB,NOBPSM,NOCPSM,ITASK,IPACK)
*
* Delete OCC/VIRT blocks of an orbital matrix A in MO basis
*
* ITASK = 1 : Delete evrything but VIR-VIR block
* ITASK = 2 : Delete evrything but OCC-OCC block
* ITASK = 3 : Delete evrything but OCC-OCC and VIR-VIR blocks
*
*. Jeppe Olsen, Oct. 2004
*
      INCLUDE 'implicit.inc'
*. Input 
      INTEGER NOBPSM(*),NOCPSM(*)
*. Matrix  in complete blocked form (IPACK = 0) or lower half packed 
*. form 
      DIMENSION A(*)
*
      IOFF_SM = 2810
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF_SM = 1
        ELSE 
          IF(IPACK.EQ.0) THEN 
            IOFF_SM = IOFF_SM + NOBPSM(ISM-1)**2
          ELSE 
            IOFF_SM = IOFF_SM + NOBPSM(ISM-1)*(NOBPSM(ISM-1)+1)/2
          END IF
        END IF
*
        NOC = NOCPSM(ISM)
        NOB = NOBPSM(ISM)
*. Delete OC-VI and VI-OC blocks
        DO I = 1, NOC
        DO J = NOC+1,NOB
          IF(IPACK.EQ.0) THEN
            A(IOFF_SM-1+(J-1)*NOB+I) = 0.0D0
            A(IOFF_SM-1+(I-1)*NOB+J) = 0.0D0
          ELSE 
            A(IOFF_SM-1+J*(J-1)/2+I) = 0.0D0
          END IF
        END DO
        END DO
*.
        IF(ITASK.EQ.1) THEN
*. Delete also OC-OC block
          DO I = 1, NOC
            IF(IPACK.EQ.0) THEN
              JMAX = NOC
            ELSE 
              JMAX = I
            END IF
            DO J = 1, JMAX
              IF(IPACK.EQ.0) THEN
                A(IOFF_SM-1+(J-1)*NOB + I) = 0.0D0
              ELSE 
                A(IOFF_SM-1+I*(I-1)/2 + J) = 0.0D0
              END IF
            END DO
         END DO
        END IF
*
        IF(ITASK.EQ.2) THEN
*. Delete virtual-virtual block
          DO I = NOC+1,NOB
            IF(IPACK.EQ.0) THEN
              JMAX = NOB
            ELSE 
              JMAX = I
            END IF
            DO J = NOC+1, JMAX
              IF(IPACK.EQ.0) THEN
                A(IOFF_SM-1+(J-1)*NOB + I) = 0.0D0
              ELSE 
                A(IOFF_SM-1+I*(I-1)/2 + J) = 0.0D0
              END IF
            END DO
          END DO
        END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output matrix from DEL_OV_ORBMAT '
        CALL APRBLM2(A,NOBPSM,NOBPSM,NSMOB,IPACK)
      END IF
*
      RETURN
      END 
C              CALL GET_FIVO(WORK(KFAO),WORK(KHAO),WORK(KLFIVO),
C    &               NOCOBS,NAOS_ENV,NEL,WORK(KCCUR),WORK(KSAO),
C    /               WORK(KLSCR1),WORK(KLSCR2))
      SUBROUTINE GET_FIVO(F,H,FIVO,NOCPSM,NOBPSM,NSMOB,NEL,C,SAO_E,
     &           SCR1,SCR2,SCR3)
*
*. Modify Virt-Virt blocks of Fock matrix so it corresponds 
*. to an (NEL-1) average field
*
*. Jeppe Olsen, Oct 2004
*
      INCLUDE 'implicit.inc'
*. Input : matrices in packed form, in AO basis
      DIMENSION F(*), H(*),SAO_E(*),C(*)
*. output 
      DIMENSION FIVO(*)
*. Scratch
      DIMENSION SCR1(*),SCR2(*),SCR3(*)
*
      LENC = NDIM_1EL_MAT(1,NOBPSM,NOBPSM,NSMOB,0)
C?    WRITE(6,*) ' LENC = ', LENC
*. obtain effective e-e field in FIVO
      ONE = 1.0D0
      ONEM = -1.0D0
      CALL VECSUM(FIVO,F,H,ONE,ONEM,LENC)
*. Transform to MO basis
C     TRAN_SYM_BLOC_MAT4
C    &(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      CALL TRAN_SYM_BLOC_MAT4
     &     (FIVO,C,C,NSMOB,NOBPSM,NOBPSM,SCR1,SCR2,1)
      CALL COPVEC(SCR1,FIVO,LENC)
*. Eliminate OC-OC,OC-VI,VI-OC blocks from FIVO
C     DEL_OV_ORBMAT(A,NSMOB,NOBPSM,NOCPSM,ITASK,IPACK)
      CALL DEL_OV_ORBMAT(FIVO,NSMOB,NOBPSM,NOCPSM,1,1)
*. construct (NEL-1)/NEL -1)*F_EE
      FACTOR = FLOAT(NEL-1)/FLOAT(NEL) - 1
      WRITE(6,*) ' FACTOR = ', FACTOR
      CALL SCALVE(FIVO,FACTOR,LENC)
C?    WRITE(6,*) ' Scaled F(e-e) '
C?    CALL APRBLM2(FIVO,NOBPSM,NOBPSM,NSMOB,1)
*. Transform back with C**(-1) = C^T S
*. Transform with C^T
      CALL TRP_BLK_MAT(C,SCR2,NSMOB,NOBPSM,NOBPSM)
C  TRP_BLK_MAT(AIN,AOUT,NBLK,LROW,LCOL)
      CALL TRAN_SYM_BLOC_MAT4
     &     (FIVO,SCR2,SCR2,NSMOB,NOBPSM,NOBPSM,SCR1,SCR3,1)
C?    WRITE(6,*) ' After transformation with C^T '
C?    CALL APRBLM2(SCR1,NOBPSM,NOBPSM,NSMOB,1)
*. And with S
      WRITE(6,*) ' Input SAO_E matrix' 
      CALL APRBLM2(SAO_E,NOBPSM,NOBPSM,NSMOB,0)
      WRITE(6,*) ' Input SCR1 matrix '
      CALL APRBLM2(SCR1 ,NOBPSM,NOBPSM,NSMOB,0)
      WRITE(6,*) ' Input SCR2 matrix '
      CALL APRBLM2(SCR2 ,NOBPSM,NOBPSM,NSMOB,0)
      CALL TRAN_SYM_BLOC_MAT4
     &     (SCR1,SAO_E,SAO_E,NSMOB,NOBPSM,NOBPSM,FIVO,SCR2,1)
C?    WRITE(6,*) ' After transformation with S  '
C?    CALL APRBLM2(FIVO,NOBPSM,NOBPSM,NSMOB,1)
*. And add F to FIVO
      CALL VECSUM(FIVO,FIVO,F,ONE,ONE,LENC)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' F-matrix for improved virtual orbitals '
       CALL APRBLM2(FIVO,NOBPSM,NOBPSM,NSMOB,1)
      END IF
*
      RETURN
      END 
      SUBROUTINE ANA_STABILITY_FOR_HF(C,I_WILL_DO_TRANS,I_FOLLOW_MODE,
     &                                INEG_TOTSYM,CNEW,NEW_MO_OBTAINED)
*
*. Perform stability analysis by obtaining and diagonalizing SCF Hessian 
*. for all symmetries
*
*. Integrals are assumed in place
*
*
* If I_FOLLOW_MODE, then the code will do a linesearch along the 
* lowest negative eigenvector of total symmetrix type, and the 
* code will return the corresponding new mo coefs in CNEW.
*
* Absence or presence of negative total symmetric modes are flagged by
* INEG_TOTSYM
*
*. Jeppe Olsen, Oct. 2004
*
* Last revision, Aug. 30 2012, Jeppe Olsen
*
* 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
*. Input (not active..)
       DIMENSION C(*)
*
      NTEST = 10
      IF(NTEST.GE.1) THEN
        WRITE(6,*) 
        WRITE(6,*) ' Information about stability analysis of Hessian:'
        WRITE(6,*) ' ================================================'
        WRITE(6,*)
      END IF
*
      IDUM = 0
      NEW_MO_OBTAINED = 0
*. Tell integral-fetchers to get normal integrals
      IH1FORM = 1
      IH2FORM = 1
      I_USE_SIMTRH = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'ANA_S')
*
* Transform integrals with C. At the moment MO-MO  transformation 
* is used, so it is assumed  the C is MO-MO transf
*
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      IF(I_WILL_DO_TRANS.EQ.1) THEN
*. Is not correctly working pt
*. KMOMO is used by the transformation program
        LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
        CALL COPVEC(C,WORK(KMOMO),LENC)
        CALL TRAINT
*. Obtain inactive Fock matrix - if no transformation it is assumed in place
        CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
        CALL FI(WORK(KINT1),ECORE_HEX,1)
      END IF
*. Nonredundant excitation between different types
*. Nonredundant type-type excitations
      CALL MEMMAN(KLTTACT,(NGAS+2)**2,'ADDL  ',1,'TTACT ')
      CALL NONRED_TT_EXC(WORK(KLTTACT),0,0)
*. Space for non-redundant operators - max
*. Nonredundant orbital excitations
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
      CALL MEMMAN(KLOOEXCC,2*NTOOB*NTOOB,'ADDL  ',1,'OOEXCC')
*. Space for saving lowest eigenvector 
      CALL MEMMAN(KLEIGV,NTOOB*NTOOB,'ADDL  ',2,'EIGV  ')
*
      EIGVMIN = 0.0D0
      INEG_TOTSYM = 0
      DO IJSM = 1, NSMOB
*. Generate orbital excitations of this symmetry
        CALL NONRED_OO_EXC(NOOEXC,WORK(KLOOEXC),WORK(KLOOEXCC),
     &                     IJSM,WORK(KLTTACT),2)
*. Obtain Hessian for this symmetry
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HESGEN')
        CALL MEMMAN(KLE2,NOOEXC*(NOOEXC+1)/2,'ADDL  ',2,'E2    ')
        CALL E2_FUSK(WORK(KLE2),NOOEXC,WORK(KLOOEXCC),0)
*. Diagonalize
        CALL MEMMAN(KLX,NOOEXC*NOOEXC,'ADDL  ',2,'X     ')
        CALL MEMMAN(KLSCR,NOOEXC*(NOOEXC+1)/2,'ADDL  ',2,'EIGVL ')
        CALL DIAG_SYM_MAT(WORK(KLE2),WORK(KLX),WORK(KLSCR),NOOEXC,1)
C             DIAG_SYM_MAT(A,X,SCR,NDIM,ISYM)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Eigenvalues for symmetry = ', IJSM
          CALL WRTMAT(WORK(KLSCR),1,NOOEXC,1,NOOEXC)
        END IF
        NNEG = 0
        DO I = 1, NOOEXC
          IF(WORK(KLSCR-1+I).LT.0.0D0) NNEG = NNEG + 1
        END DO
        IF(NNEG.EQ.0) THEN
          WRITE(6,*) ' No negative eigenvalues for symmetry ', IJSM
        ELSE
          WRITE(6,*) ' Symmetry and number of negative eigenvalues ',
     &    IJSM,  NNEG
          WRITE(6,*) ' The negative eigenvalues '
          CALL WRTMAT(WORK(KLSCR),1,NNEG,1,NNEG)
           WRITE(6,*) ' Major components of lowest mode '
           THRES = 0.1D0
           CALL WRT_IOOEXCOP_WTHRES(WORK(KLX),WORK(KLOOEXCC),
     &               NOOEXC,THRES)
C               WRT_IOOEXCOP_WTHRES(XOOEXC,IOOEXC,NOOEXC,THRESREL)
        END IF
*. Save eigvector for lowest total symmetric negative eigenvalue, 
*. as usual I assume that the eigenvalues are returned ordered
        IF(IJSM.EQ.1.AND.WORK(KLSCR).LT.0.0D0) THEN
          EIGVMIN = WORK(KLSCR)
          CALL COPVEC(WORK(KLX),WORK(KLEIGV),NOOEXC)
          INEG_TOTSYM = 1
        END IF
        CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HESGEN')
      END DO
*
      IF(I_FOLLOW_MODE.EQ.1.AND.INEG_TOTSYM.EQ.1) THEN
*
*. A negative eigenvalue was obtained in the space of 
* totally symmetric excitations, minimize energy 
* along this path
*
*. Regenerate total symmetric excitations
        CALL NONRED_OO_EXC(NOOEXC,WORK(KLOOEXC),WORK(KLOOEXCC),
     &                     1,WORK(KLTTACT),2)
         IF(NTEST.GE.10) THEN
           WRITE(6,*) ' Major components of mode to be followed '
C               WRT_IOOEXCOP_WTHRES(XOOEXC,IOOEXC,NOOEXC,THRESREL)
           THRES = 0.1D0
           CALL WRT_IOOEXCOP_WTHRES(WORK(KLEIGV),WORK(KLOOEXCC),
     &               NOOEXC,THRES)
         END IF
*. Put 1-electron integrals back in KINT1
        CALL COPVEC(WORK(KH),WORK(KINT1),LENP)
        CALL MIN_EHF_ALONG_MODE(WORK(KLEIGV),WORK(KLOOEXCC),NOOEXC,
     &                          EMIN,ALPHAMIN)
        WRITE(6,*) ' ALPHAMIN returned from MIN... ', ALPHAMIN
*. Obtain new MO-coefs associated with step

        CALL SCALVE(WORK(KLEIGV),ALPHAMIN,NOOEXC)
C     GET_NEWMO_FROM_KAPPA(CNEW,COLD,KAPPA,IOOEXC,NOOEXC,ICOLD_IS_UNI)
        ICOLD_IS_UNI = 0
        CALL GET_NEWMO_FROM_KAPPA(CNEW,WORK(KCCUR),WORK(KLEIGV),
     &       WORK(KLOOEXCC),NOOEXC,ICOLD_IS_UNI)
*. and flag that new MO's have been obtained
        NEW_MO_OBTAINED = 1
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'ANA_S')
*
*. 
      RETURN
      END
      SUBROUTINE OPTIM_SCF_USING_SECORDER(MAXIT_HF,E_HF,CONVER_HF,
     &           E1_FINAL,NIT_HF,E_ITER)       
*
* Optimize SCF wavefunction using second-order method
* 
* Current code is EXTREMELY INEFFICIENT using
* complete integral transformations and 
* explicit construction of all matrices including Hessian 
* in MO basis
*
* Current version assumes closed shell states
*
*. On input an initial MO-AO is assumed residing in work(kccur)
*
*. Jeppe Olsen, 2004
*. Last revision, Aug. 29 2012, Jeppe Olsen - updated to modern times..
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      REAL*8 INPROD
      DIMENSION E_ITER(*)
      LOGICAL CONVER_HF
*
      CALL LUCIAQENTER('OPTHF')
*
      IDUMMY = 0
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ', IDUMMY,'OPTIM ')
*
      CONVER_HF = .FALSE.
*
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
*
      CALL MEMMAN(KLTTACT,(NGAS+2)*(NGAS+2),'ADDL  ',1,'TTACT ')
      CALL NONRED_TT_EXC(WORK(KLTTACT),0,0)
*. Space for non-redundant operators
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
*. Space for Nonredundant orbital excitations
      CALL MEMMAN(KLOOEXCC,2*NTOOB*NTOOB,'ADDL  ',1,'OOEXCC')
*. And the nonredundant operators of symmetry 1
      CALL NONRED_OO_EXC(NOOEXC,WORK(KLOOEXC),WORK(KLOOEXCC),
     &                     1,WORK(KLTTACT),2)
*. Space for gradient and Hessian and expanded Hessian
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1    ')
      LEN_E2 = NOOEXC*(NOOEXC+1)/2
      CALL MEMMAN(KLE2,LEN_E2,'ADDL  ',2,'E2    ')
      CALL MEMMAN(KLE2E,NOOEXC**2,'ADDL  ',2,'E2_E  ')
*. Space for Metric ( is actually just unit matrix ..)
      CALL MEMMAN(KLS2E,NOOEXC**2,'ADDL  ',2,'S2_E  ')
*. Space for Kappa parameters
      CALL MEMMAN(KLKAPPA,NOOEXC,'ADDL  ',2,'KAPPA ')
*. Space for a new set of MO's
      CALL MEMMAN(KLCNEW,LENC,'ADDL  ',2,'CNEW  ')
*. Space for saving number of inactive orbitals per symmetry
      CALL MEMMAN(KLNINOBSS,NSMOB,'ADDL  ',2,'NINOSS')
*. Tell the integrals get normal integrals
      IH1FORM = 1
      IH2FORM = 1
      I_USE_SIMTRH = 0
*. We will set the inactive orbitals to the doubly occupied, save the
*. original number of inactive orbitals
      CALL ICOPVE(NINOBS,WORK(KLNINOBSS),NSMOB)
      CALL ICOPVE(NHFD_STASYM,NINOBS,NSMOB)
*. Set metric to unit 
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLS2E),ZERO,NOOEXC**2)
      ONE = 1.0D0
      CALL SETDIA(WORK(KLS2E),ONE,NOOEXC,0)

*. Branch point to which I can return if stability analysis
*. reveals that stability exists
*. Max number of stability analysis 
      MAX_MACRO = 5
      IMACRO = 0
00001 CONTINUE
      IMACRO = IMACRO + 1
*
      THRES_E = 1.0D-10
      EOLD = -0803
      EHF  = -0803
*. Obtain input density from input orbitals
      WRITE(6,*) ' Largest allowed number of iterations', MAXIT_HF
      DO ITER = 1, MAXIT_HF + 1
        WRITE(6,*)
        WRITE(6,*) ' =========================================='
        WRITE(6,*) ' Information from iteration ', ITER
        WRITE(6,*) ' =========================================='
        WRITE(6,*)
        NIT_HF = ITER
*
        E_HF = EHF
        IF(ITER.EQ.MAXIT_HF+1) WRITE(6,*) ' (just transf and energy) '
        WRITE(6,*)
        IF(ITER.GT.1) EOLD = EHF
*. Integral transformation with MO-coefs in KMOMO
        CALL COPVEC(WORK(KCCUR),WORK(KMOMO),LENC)
        KKCMO_I = KMOMO
        KKCMO_J = KMOMO
        KKCMO_K = KMOMO
        KKCMO_L = KMOMO
        CALL TRAINT
        CALL COPVEC(WORK(KINT1),WORK(KH),LENP)
COLD    CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
*. Obtain density matrix in AO basis
        CALL GET_P_FROM_C(WORK(KCCUR),WORK(KDAO),NHFD_STASYM,
     &       NTOOBS,NSMOB)
*. Obtain inactive Fock matrix
        WRITE(6,*)
     &  ' New Fock matrix will be generated '
COLD    CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
*. Obtain inactive Fock matrix
        CALL FISM(WORK(KINT1),EHF)
*. Add core-energy 
        EHF = EHF + ECORE_ORIG
        E_ITER(NIT_HF) = EHF
        E_HF = EHF
        IF(ITER.EQ.MAXIT_HF+1) GOTO 1111
        IF(ITER.GT.1) THEN
          DELTA_E = ABS(EHF-E_ITER(ITER-1))
          WRITE(6,*) ' Change of energy = ', EHF-EOLD
          IF(DELTA_E.LT.THRES_E) THEN 
            WRITE(6,*) ' Energy converged '
            CONVER_HF = .TRUE.
            GOTO 1111
          END IF
        END IF
*. Obtain Gradient
        CALL GET_E1HF_FROM_FI(WORK(KLE1),NOOEXC,WORK(KLOOEXCC))
        E1NORM = SQRT(INPROD(WORK(KLE1),WORK(KLE1),NOOEXC))
        WRITE(6,'(A,I5,2E23.12)') ' ITER, Energy and grad-norm', 
     &  ITER,EHF,E1NORM
*. Obtain Hessian 
*
        IONLYF = 0
        IF(IONLYF.EQ.1) THEN
          WRITE(6,*) ' Hessian only including Fterms '
        END IF
        CALL E2_FUSK(WORK(KLE2),NOOEXC,WORK(KLOOEXCC),IONLYF)
*. Expand  Hessian to complete form 
        CALL TRIPAK_BLKM(WORK(KLE2E),WORK(KLE2),2,NOOEXC,1)
*. Solve Trust region NR equations
C     SOLVE_NR_WITH_GEN_TRCTL(E1,E2,S,STEP,NVAR,XMAXNORM)
         XMAXNORM = 0.60D0
         CALL SOLVE_NR_WITH_GEN_TRCTL(WORK(KLE1),WORK(KLE2E),
     &        WORK(KLS2E),WORK(KLKAPPA),NOOEXC,XMAXNORM)
*. Obtain new MO  coeffients in WORK(KCCUR), As we pt are 
*  transforming from MO's to MO's the old set oc MO's are 
*  simple the unit matrix 
C     GET_NEWMO_FROM_KAPPA(CNEW,COLD,KAPPA,IOOEXC,NOOEXC,
C    &                                ICOLD_IS_UNI)
         CALL GET_NEWMO_FROM_KAPPA(WORK(KCCUR),WORK(KMOMO),
     &        WORK(KLKAPPA),WORK(KLOOEXCC),NOOEXC,0)
      END DO
 1111 CONTINUE
*
      WRITE(6,*)
      WRITE(6,*) ' Final energy in au  ', EHF
      WRITE(6,*) ' Final set of MOs: '
      CALL PRINT_CMOAO(WORK(KCCUR))
*
*. Check for negative modes, and do a mode-following along 
*. lowest negative mode
*
* IT WOULD BE BETTER TO TAKE THIS OUTSIDE THE OPTIMIZER 
* TO THE LEVEL ABOVE
      I_ANA_STA = 1
      I_FOLLLOW_MODE = 0
      IF(I_ANA_STA.EQ.1) THEN
        I_WILL_DO_TRANS = 0
        I_FOLLOW_MODE = 1
        CALL ANA_STABILITY_FOR_HF(WORK(KCCUR),I_WILL_DO_TRANS,
     &       I_FOLLOW_MODE,INEG_TOTSYM,WORK(KLCNEW),NEW_MO_OBTAINED)
        IF(NEW_MO_OBTAINED.EQ.1) THEN
*. We should go back and redo SCF with new mo's
            CALL COPVEC(WORK(KLCNEW),WORK(KCCUR),LENC)
            WRITE(6,*) ' CNEW copied to CCUR '
            WRITE(6,*) ' Updated MO coefficients after ANA_STABILITY'
            CALL PRINT_CMOAO(WORK(KLCNEW))
            IF(IMACRO.LT.MAX_MACRO) GOTO 1
        END IF
C     ANA_STABILITY_FOR_HF(C,I_WILL_DO_TRANS,I_FOLLOW_MODE,
C    &                                INEG_TOTSYM,CNEW,NEW_MO_OBTAINED)
      END IF
*. Restore number of inactive orbitals
      CALL ICOPVE(WORK(KLNINOBSS),NINOBS,NSMOB)
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ', IDUMMY,'OPTIM ')
      CALL LUCIAQEXIT('OPTHF')

      RETURN
      END
      SUBROUTINE GET_NEWMO_FROM_KAPPA(CNEW,COLD,KAPPA,IOOEXC,NOOEXC,
     &                                ICOLD_IS_UNI)
*
* Obtain New mo's corresponding to a given set of Kappa-parameters
*
*. Jeppe Olsen, Oct. 2004
*
* If ICOLD_IS_UNI = 1, then the old MO coefs, are assumed to 
* be a unit matrix, and do not need to be explicitly given
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      DIMENSION IOOEXC(2,NOOEXC),COLD(*)
      REAL*8 KAPPA(*)
*. Output
      DIMENSION CNEW(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTMOKA')
*. Obtain Exponential of -Kappa
      LEN_E =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)     
      CALL MEMMAN(KLEXPMK,LEN_E,'ADDL  ',2,'EXPMK ')
C?    WRITE(6,*) ' NOOEXC, 1 = ', NOOEXC
      CALL GET_EXP_MKAPPA(WORK(KLEXPMK),KAPPA,IOOEXC,NOOEXC)
*. Cnew = Cold * Exp -Kappa
C      MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)

      IF(ICOLD_IS_UNI.EQ.0) THEN
        CALL MULT_BLOC_MAT(CNEW,COLD,WORK(KLEXPMK),NSMOB,NTOOBS,NTOOBS,
     &                     NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
      ELSE
        CALL COPVEC(WORK(KLEXPMK),CNEW,LEN_E)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' New set of MO''s '
        CALL APRBLM2(CNEW,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTMOKA')
*
      RETURN
      END
      SUBROUTINE GET_EXP_MKAPPA(EXPMK,KAPPAP,IOOEXC,NOOEXC)
*
* Obtain Exp(-Kappa) where Kappa is given in packed form
*
*. Jeppe Olsen, Oct. 2004
*. Last revision: Aug. 30 2012, Jeppe Olsen, Some print added 
*
c      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific Input
      DIMENSION IOOEXC(2,NOOEXC)
      REAL*8 KAPPAP(NOOEXC)
*. Output
      DIMENSION EXPMK(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_EXP_MKAPPA '
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Input: Kappa in packed form '
        CALL WRT_EXCVEC(KAPPAP,IOOEXC,NOOEXC)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTEXKA')
*. Obtain Kappa matrix in expanded form
      LEN_E =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)     
      CALL MEMMAN(KLKAPPAE,LEN_E,'ADDL  ',2,'KAPPAE')
      CALL REF_AS_KAPPA(KAPPAP,dbl_mb(KLKAPPAE),1,1,IOOEXC,NOOEXC)
*
*. And obtain the exponential of -Kappa - done in symmetryblocks
      NOBBLK_MAX = IMNMX(NTOOBS,NSMOB,2)
      LSCR = 4*NOBBLK_MAX ** 2 + 3*NOBBLK_MAX
      CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'KAPSCR')
      IOFF = 1
      DO ISM = 1, NSMOB
       IF(ISM.GT.1) IOFF = IOFF + NTOOBS(ISM-1) ** 2
       NOB = NTOOBS(ISM)
       CALL EXPMA(EXPMK(IOFF),dbl_mb(KLKAPPAE-1+IOFF),NOB,
     &            dbl_mb(KLSCR),0)
C           EXPMA(EMA,A,NDIM,SCR,ISUB)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Exp -Kappa '
        CALL APRBLM2(EXPMK,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTEXKA')
*
      RETURN
      END
      SUBROUTINE REF_AS_KAPPA(KAPPAP,KAPPAE,ISM,IWAY,IOOEX,NOOEX)
*
* Switch between packed and expanded form of an antisymmetrix kappa matrix 
* of symmetry ISM. 
* KAPPA_E is expanded matrix in ST order
* KAPPA_P is packed matrix in TS order 
*
*
* IWAY = 1 Packed => expanded form
* IWAY = 2 Expanded to packed form 
*
*. Jeppe Olsen, Oct. 2004
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input and output
      INTEGER IOOEX(2,NOOEX)
      REAL*8  KAPPAP(*), KAPPAE(*)
*. local scratch
      INTEGER IOFFSM(MXPOBS)
*. Set up array giving offset for block of given symmetry
*  GET_OFFSET_FOR_SYMBLK(NROW,NCOL,IRC,ISYM,NSMST,IB)
      CALL GET_OFFSET_FOR_SYMBLK(NTOOBS,NTOOBS,1,ISM,NSMOB,IOFFSM)
*. Number of elements in expanded matrix 
      LENGTH_E = IOFFSM(NSMOB) + NTOOBS(NSMOB)**2 - 1
*
      ZERO = 0.0D0
      CALL SETVEC(KAPPAE,ZERO,LENGTH_E)
      DO IEX = 1, NOOEX
*. Orb numbers in TS- order 
        ICREA_TS = IOOEX(1,IEX)
        IANNI_TS = IOOEX(2,IEX)
*. and in ST-order
        ICREA_ST = IREOTS(ICREA_TS)
        IANNI_ST = IREOTS(IANNI_TS)
*
        ICREA_SM = ISMFSO(ICREA_ST)
        IANNI_SM = ISMFSO(IANNI_ST)
*
        ICREA_OFF = IBSO(ICREA_SM)
        IANNI_OFF = IBSO(IANNI_SM)
*
        IAC_E = IOFFSM(IANNI_SM) -1  
     &       + (ICREA_ST-ICREA_OFF)*NTOOBS(IANNI_SM)
     &       + IANNI_ST - IANNI_OFF + 1
        ICA_E = IOFFSM(ICREA_SM) - 1
     &        + (IANNI_ST-IANNI_OFF)*NTOOBS(ICREA_SM)
     &        + ICREA_ST - ICREA_OFF + 1
C?      WRITE(6,*) ' IAC_E, ICA_E = ', IAC_E, ICA_E
        IF(IWAY.EQ.1) THEN
          KAPPAE(ICA_E) = KAPPAP(IEX)
          KAPPAE(IAC_E) = -KAPPAP(IEX)
        ELSE
          KAPPAP(IEX) = KAPPAE(ICA_E)
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from  REF_AS_KAPPA'
        IF(IWAY.EQ.1) THEN 
          WRITE(6,*) ' Packed => expanded form '
        ELSE IF (IWAY.EQ.2) THEN
          WRITE(6,*) ' Expanded => packed form '
        END IF
        WRITE(6,*) ' Kappa-matrix in expanded form '
        CALL APRBLM2(KAPPAE,NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*) ' Kappa-matrix in packed form '
        CALL WRT_IOOEXCOP(KAPPAP,IOOEX,NOOEX)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_IOOEXCOP(XOOEXC,IOOEXC,NOOEXC)
*
* Print orbital-excitation vector
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
      DIMENSION XOOEXC(NOOEXC),IOOEXC(2,NOOEXC)
*
C?    WRITE(6,*) 'NOOEXC in WRT... ', NOOEXC
*
      WRITE(6,*) ' Crea   Anni           Value '
      WRITE(6,*) ' ============================'
      DO IEX = 1, NOOEXC
        WRITE(6,'(2I6,E22.12)') IOOEXC(1,IEX),IOOEXC(2,IEX),XOOEXC(IEX)
      END DO
*
      RETURN
      END
      SUBROUTINE WRT_IOOEXCOP_WTHRES(XOOEXC,IOOEXC,NOOEXC,THRESREL)
*
* Print orbital-excitation vector
*
* Only elements that are larger than THRESREL* largest value is printed
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
      DIMENSION XOOEXC(NOOEXC),IOOEXC(2,NOOEXC)
*. Largest value
      XMAX = XMNMX(XOOEXC,NOOEXC,3)
*
C?    WRITE(6,*) 'NOOEXC in WRT... ', NOOEXC
*
      WRITE(6,*) ' Crea   Anni           Value '
      WRITE(6,*) ' ============================'
      DO IEX = 1, NOOEXC
        IF(ABS(XOOEXC(IEX)).GE.THRESREL*XMAX) 
     &  WRITE(6,'(2I6,E22.12)') IOOEXC(1,IEX),IOOEXC(2,IEX),XOOEXC(IEX)
      END DO
*
      RETURN
      END
      SUBROUTINE EXPMA(EMA,A,NDIM,SCR,ISUB)
*
* Obtain exponential of antisymmetric matrix -A, exp(-A)
*
* IF ISUB .NE. 0 exp(-A) - I is constructed.
*
      INCLUDE 'implicit.inc'
      DIMENSION EMA(NDIM,NDIM),A(NDIM,NDIM),SCR(*)
*
* SCR should at least be of length 4*NDIM**2 + 3*NDIM
*
* Jeppe Olsen , June 1989 for LUCAS 
*               Oct. 2004, adopted for LUCIA
*
*
      NTEST = 00
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Output from EXPMA '
        WRITE(6,*) ' ================= '
        WRITE(6,*)
        WRITE(6,*) ' Input matrix '
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
** 1 : Local memory
*
      KLFREE = 1
* A ** 2
      KLA2 =  KLFREE
      KLFREE = KLFREE + NDIM ** 2
* Eigenvectors of A ** 2
      KLA2VC = KLFREE
      KLFREE = KLFREE + NDIM ** 2
* Eigenvalues of A ** 2
      KLA2VL = KLFREE
      KLFREE = KLFREE + NDIM
* Extra matrix
      KLMAT1 = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*. And an extra, extra matrix
      KLSCR = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*. And two vectors 
      KLAR1 = KLFREE
      KLFREE = KLFREE + NDIM
*
      KLAR2 = KLFREE
      KLFREE = KLFREE + NDIM
C?    WRITE(6,*) ' Memory used in EXPMA = ', KLFREE
*
** Obtain A ** 2 and diagonalize
*
      CALL MATML4(SCR(KLA2),A,A,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0)
      CALL TRIPAK(SCR(KLA2),SCR(KLSCR),1,NDIM,NDIM)
      CALL EIGENL(SCR(KLSCR),SCR(KLA2VC),NDIM,0,1)
      CALL COPDIA(SCR(KLSCR),SCR(KLA2VL),NDIM,1)
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Eigenvalues of A squared '
        CALL WRTMAT(SCR(KLA2VL),NDIM,1,NDIM,1)
      END IF
*
** 3 Obtain arrays sum(n) e ** n /(2n)! and e **n /(2n+1)!
*
* Max number of terms required for iterative method
      NTERM = -1
      XMAX = FNDMNX(SCR(KLA2VL),NDIM,2)
      IF(XMAX.LT.1.0D-1) THEN
        NTERM = 0
        TEST = 1.0D-15
        ELMNT = XMAX
        X2NP1 = 1.0D0
  230   CONTINUE
          NTERM = NTERM + 1
          X2NP1 = X2NP1 + 2.0D0
          ELMNT = ELMNT*XMAX/(X2NP1*(X2NP1-1))
          IF( XMAX .EQ. 0.0D0 ) GOTO 231
C       IF(ELMNT/XMAX.GT.TESt) GOTO 230
        IF(ELMNT     .GT.TESt) GOTO 230
  231   CONTINUE
      END IF
      IF(NTEST.GE.10)
     &WRITE(6,*) ' XMAX NTERM ', XMAX,NTERM
*
      IF(XMAX.GE.1.0D-1) THEN
* Explicit use of sine and cosine formulaes
        DO 2810  I = 1 , NDIM
          IF(SCR(KLA2VL-1+I).GE.0.0D0) THEN
* All eigenvalues should be nonpositive,
* set small positive eigenvalues to zero
             SCR(KLAR1-1+I) = 1.0D0
             IF(ISUB.EQ.0) THEN
               SCR(KLAR2-1+I) = 1.0D0
             ELSE
               SCR(KLAR2-1+I) = 0.0D0
             END IF
           ELSE
             DELTA = SQRT(-SCR(KLA2VL-1+I))
             SCR(KLAR1-1+I) = SIN(DELTA)/DELTA
             IF( ISUB.EQ.0) THEN
               SCR(KLAR2-1+I) = COS(DELTA)
             ELSE
               SCR(KLAR2-1+I) = COS(DELTA) - 1.0D0
             END IF
           END IF
 2810   CONTINUE
      ELSE
* Iterative formulaes ( actually only needed for stability of array 2, if 
*                       ISUB = 1 )
*
        CALL SETVEC(SCR(KLAR1),1.0D0,NDIM)
        CALL SETVEC(SCR(KLSCR),1.0D0,NDIM)
        X2NP1 = 1.0D0
        DO 300 N = 1, NTERM
          X2NP1 = X2NP1 + 2.0D0
          FACTOR = 1.0D0/(X2NP1*(X2NP1-1))
          CALL VVTOV(SCR(KLA2VL),SCR(KLSCR),SCR(KLSCR),NDIM)
          CALL SCALVE(SCR(KLSCR),FACTOR,NDIM)
          CALL VECSUM(SCR(KLAR1),SCR(KLAR1),SCR(KLSCR),
     &                1.0D0,1.0D0,NDIM)
  300   CONTINUE
*
        IF(ISUB.EQ.0) THEN
          CALL SETVEC(SCR(KLAR2),1.0D0,NDIM)
        ELSE
          CALL SETVEC(SCR(KLAR2),0.0D0,NDIM)
        END IF
        CALL  SETVEC(SCR(KLSCR),1.0D0,NDIM)
        X2N = 0.0D0
        DO 330 N = 1, NTERM
          X2N = X2N + 2.0D0
          FACTOR = 1.0D0/(X2N*(X2N-1))
          CALL VVTOV(SCR(KLA2VL),SCR(KLSCR),SCR(KLSCR),NDIM)
          CALL SCALVE(SCR(KLSCR),FACTOR,NDIM)
          CALL VECSUM(SCR(KLAR2),SCR(KLAR2),SCR(KLSCR),
     &                1.0D0,1.0D0,NDIM)
  330   CONTINUE
      END IF
*
** 4  Obtain Exp(-A)
*
* A * U * Dia1 * U(t)
      KLMAT2 = KLA2
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(SCR(KLMAT1),SCR(KLA2VC),SCR(KLAR1),
     &            NDIM,SCR(KLSCR) )
      CALL MATML4(SCR(KLMAT2),A,SCR(KLMAT1),
     &            NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0)
* U * Dia2 * U(T)
      CALL XDIAXT(SCR(KLMAT1),SCR(KLA2VC),SCR(KLAR2),
     &            NDIM,SCR(KLSCR) )
*
      CALL VECSUM(EMA,SCR(KLMAT1),SCR(KLMAT2),
     &            1.0D0,-1.0D0,NDIM ** 2 )
*
      IF( NTEST .GE. 5 ) THEN
        IF(ISUB .EQ. 0 ) THEN
          WRITE(6,*) '    Exp(-A)'
        ELSE
          WRITE(6,*) ' Exp(-A) - I '
        END IF
        WRITE(6,*) '  ============='
        CALL WRTMAT_EP(EMA,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_E1HF_FROM_FI(E1,NOOEXC,IOOEXC)
*
* Obtain HF gradient in MO basis from FI matrix also in MO basis
* FI is assumed residing in WORK(KINT1)
*
* E1(A,I) = -4*FI(I,A)
*
*. Jeppe Olsen, Oct. 2004 ( Yes, I admit...)
*. Last Modification, Sept. 24 2012, Jeppe Olsen, Modification of print
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc' 
*. Specific input
      DIMENSION IOOEXC(2,NOOEXC)
*. Output
      DIMENSION E1(NOOEXC)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GET_E1HF_FROM_FI '
      END IF
*
C?    WRITE(6,*) ' NOOEXC in GET_E1HF.. ', NOOEXC
      DO IEX = 1, NOOEXC
        II = IOOEXC(1,IEX)
        IA = IOOEXC(2,IEX)
        E1(IEX) = -4.0D0*GETH1_B(IA,II)
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' IEX, II, IA = ', IEX, II, IA 
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orbital gradient of HF '
        WRITE(6,*) ' ======================='
        CALL WRT_IOOEXCOP(E1,IOOEXC,NOOEXC)
      END IF
*
      RETURN
      END 
      FUNCTION EHF_FROM_KAPPAM(KAPPA)
*. Obtain energy from a set of kappa-parameters.
*. the orbital excitations are defined by transfer common block
*. CLOOXC
*
*. Outer routine for EHF_FROM_KAPPA, so reference to 
*. orbital excitations are hidden
*
*. Jeppe Olsen
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      COMMON/CLOOXC/KLIOOXC,NOOXC
*. input
      REAL*8 KAPPA(*)
*
      EHF_FROM_KAPPAM = EHF_FROM_KAPPA(KAPPA,WORK(KLIOOXC),NOOXC)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Energy obtained in  EHF_FROM_KAPPAM = ',
     &                EHF_FROM_KAPPAM
      END IF
*
      RETURN
      END
      FUNCTION EHF_FROM_KAPPA(KAPPA,IOOEXC,NOOEXC)
*
*. Obtain energy from a set of kappa-parameters
*. Fock matrix is calculated 
*
*. Jeppe Olsen, Oct. 2004
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'glbbas.inc'
*. Specific input
      REAL*8 KAPPA(NOOEXC)
      INTEGER IOOEXC(NOOEXC)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EHF_KA')
*
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLCNEW,LENC,'ADDL  ',2,'CNEW  ')
      CALL MEMMAN(KLP   ,LENC,'ADDL  ',2,'P     ')
      CALL MEMMAN(KLF   ,LENC,'ADDL  ',2,'F     ')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input Kappa parameters to EHF_FROM_KAPPA '
        CALL WRT_IOOEXCOP(KAPPA,IOOEXC,NOOEXC)
C            WRT_IOOEXCOP(XOOEXC,IOOEXC,NOOEXC)
      END IF
*. Obtain MO-matrix for current Kappa-parameters, relative to 
*. current orbitals
      ICOLD_IS_UNI = 1
      CALL GET_NEWMO_FROM_KAPPA(WORK(KLCNEW),WORK(KCCUR),KAPPA,
     &     IOOEXC,NOOEXC,ICOLD_IS_UNI)
C          GET_NEWMO_FROM_KAPPA(CNEW,COLD,KAPPA,IOOEXC,NOOEXC,
C    &                          ICOLD_IS_UNI)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Updated MO coefficients '
        CALL APRBLM_F7(WORK(KLCNEW),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. density for current MO-coefficients
      CALL GET_P_FROM_C(WORK(KLCNEW),WORK(KLP),
     &     NINOBS,NAOS_ENV,NSMOB)
      CALL GET_FOCK(WORK(KLF),WORK(KINT1),WORK(KLP),EHF,
     &                  ECORE_ORIG,3)
      EHF_FROM_KAPPA = EHF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EHF_KA')
      RETURN
      END
      SUBROUTINE LINES_SEARCH_BY_BISECTION(FUNC,REF,DIR,NVAR,XINI,
     &           XFINAL,FFINAL,IKNOW,F0,FXINI)
*
* Minimize function FUNC(REF + X*DIR) to determine X
* using only function-values
*
* IKNOW = 2: F(X=0) is given in F0, F(X=XINI) is given in FXINI
* IKNOW = 1: F(X=0) is given in F0
* IKNOW = 0: I know nothing..
*
*. Simple linesearch only using bi-section
*
*. The game to play in bisectioning is to work 
*. with three points XMIN, XMID, XMAX with XMIN < XMID < XMAX
*. and f(xmin) > f(xmid) < f(xmax). 
*. one then knows that there is a minimum between xmin and xmax,
*. and it is then a question of reducing the size of this interval
*
*. Jeppe Olsen, Oct. 2004
*
      INCLUDE 'implicit.inc'
      DIMENSION REF(NVAR),DIR(NVAR)
      EXTERNAL FUNC
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Entering LINES_SEARCH_BY_BISECT'
      END IF
*
*. Functionvalue for X=0
      IF(IKNOW.GE.1) THEN
       FX = F0
      ELSE 
       FX = FUNC(REF)
      END IF
      FMIN = FX
*
      X1 = 0.0D0
      F1 = FX
*
      I_HAVE_BRACKETING = 0
      I_AM_CONVERGED = 0
*
      MAXIT = 20
      XOLD = 0.0D0
*. Initial alpha
      IF(XINI.EQ.0.0D0) THEN
        X = 1.0D0
      ELSE
        X = XINI
      END IF
*
      ONE = 1.0D0
      DO ITER = 1, MAXIT
*. A x has been obtained, obtain function value for this alpha
        XDEL = X-XOLD
        CALL VECSUM(REF,REF,DIR,ONE,XDEL,NVAR)
        FOLD = FX
        IF(IKNOW.EQ.2.AND.X.EQ.XINI) THEN
          FX = FXINI
        ELSE
          FX = FUNC(REF)
        END IF
        XOLD = X
*. So we now have function value FX for step XOLD
*. And analyze to obtain new X
        IF(ITER.EQ.1) THEN
*. We have 2 x-values and function values
          X2 = X
          F2 = FX
          IF(FX.LE.F1) THEN
*. We went downhills, increase alpha
            X = 1.5D0*X
          ELSE
*. We went uphills - the initial value was too large
            X = X/2
          END IF
        ELSE IF (ITER.GE.2) THEN
         IF(I_HAVE_BRACKETING.EQ.0) THEN
*. We do not yet have bracketing, but have info from two previous points
*. Enroll third point
           X3 = X
           F3 = FX
*. Sort the three x's to ensure increasing order
           CALL ORD_3DATA(X1,X2,X3,F1,F2,F3,
     &          X1O,X2O,X3O,F1O,F2O,F3O)
           X1 = X1O
           X2 = X2O
           X3 = X3O
           F1 = F1O
           F2 = F2O
           F3 = F3O
C       ORD_3DATA(X1IN,X2IN,X3IN,F1IN,F2IN,F3IN,
C    &                       X1UT,X2UT,X3UT,F1UT,F2UT,F3UT)
*. Analyze these three data to see if they constitute bisection
           CALL ANA_3DATA(F1O,F2O,F3O,ITYP)
C              ANA_3DATA(F1,F2,F3,ITYP)
* ITYP = 1 : constitute an bisection : F1 > F2 < F3, 
*     i.e a minimum exist between X1 and X3
* ITYP = 2 : F1 < F2 : a minimum exist before X2, we can actually 
*            say before X1 as function starts as decreasing
* ITYP = 3 : F1 > F2 > F3 : a minumim exists after X2
           IF(ITYP.EQ.1) THEN
*. Bracketing has been found
             I_HAVE_BRACKETING = 1
             XMIN = X1
             XMID = X2
             XMAX = X3
             FMIN = F1
             FMID = F2
             FMAX = F3
             IF(NTEST.GE.100) THEN
              WRITE(6,'(A,6(2X,E13.5))') 
     &        ' Input to QUADPOL: X/F MIN MID MAX ',
     &        X1,X2,X3,F1,F2,F3
             END IF
*. Obtain new X by quadratic interpolation 
C                 FIND_QUADPOL_FOR_3PT(A,B,C,X1,X2,X3,F1,F2,F3)
             CALL FIND_QUADPOL_FOR_3PT(A,B,C,X1,X2,X3,F1,F2,F3)
C                 EXTREMUM_QUADPOL(A,B,C,XEXTR,IMNMX)
             CALL EXTREMUM_QUADPOL(A,B,C,X,IMNMX)
*. There may be a risk that X is too close to one of the points, 
*. move it 
             IF(X.EQ.XMIN) X = XMIN + 0.01*(XMID-XMIN)
             IF(X.EQ.XMAX) X = XMAX - 0.01*(XMAX-XMID)
             IF(X.EQ.XMID) X = XMID + 0.01*(XMID-XMIN)
           ELSE IF(ITYP.EQ.2) THEN
*. We stepped to far, minimum, take a point between X1 and X2 
CE           X = X1/2
             X = (X1+X2)/2
*. The points to keep in the following are X1 and X2, so 
*. no reorg.
           ELSE IF(ITYP.EQ.3) THEN
*. We should continue our search
             X = 1.5*X3
             XMIN = X2
*. The points to keep in the following are X2 and X3 so
             X1 = X2
             F1 = F2
             X2 = X3
             F2 = F3
           END IF
         ELSE
*. We are already doing bracketing, so combine new info with old 
*. to obtain improved bracketing
C     REFINE_BISECTION(XMIN,XMID,XMAX,FMIN,FMID,FMAX,X,F)
          CALL REFINE_BISECTION(XMIN,XMID,XMAX,FMIN,FMID,FMAX,X,FX)
          X1 = XMIN
          F1 = FMIN
          X2 = XMID
          F2 = FMID
          X3 = XMAX
          F3 = FMAX
*
          
*. and quadratic interpolation 
C              FIND_QUADPOL_FOR_3PT(A,B,C,X1,X2,X3,F1,F2,F3)
          CALL FIND_QUADPOL_FOR_3PT(A,B,C,X1,X2,X3,F1,F2,F3)
C              EXTREMUM_QUADPOL(A,B,C,XEXTR,IMNMX)
          CALL EXTREMUM_QUADPOL(A,B,C,X,IMNMX)
*. There may be a risk that X is too close to one of the points, 
*. move it 
          IF(X.EQ.XMIN) X = XMIN + 0.01*(XMID-XMIN)
          IF(X.EQ.XMAX) X = XMAX - 0.01*(XMAX-XMID)
          IF(X.EQ.XMID) X = XMID + 0.01*(XMID-XMIN)
         END IF
*.       ^ End if bisectioning was active
       END IF
*.     ^ End of switch depending on iteration numbers
*. Compare X and XOLD to see whether it is worth taking 
*. another round
       RATIO = ABS(X-XOLD)/X
       IF (NTEST.GT.10)
     & WRITE(6,*) ' X, XOLD, RATIO = ', X,XOLD,RATIO
       THRES = 0.0001
*. A rough estimate is sufficient
       THRES = 0.1D0
       IF(RATIO.LT.THRES) THEN 
         I_AM_CONVERGED = 1
         GOTO 1111
       END IF
      END DO
*     ^ End of loop over iterations
 1111 CONTINUE
*
      XFINAL = XOLD
      FFINAL = FX
*. Clean up, by setting reference to original reference
      XM = -XOLD
      CALL VECSUM(REF,REF,DIR,ONE,XM,NVAR)
*
      IF(NTEST.GE.100) THEN
        IF(I_AM_CONVERGED.EQ.1) THEN 
           WRITE(6,*) ' Bisection converged '
        ELSE 
           WRITE(6,*) ' Bisection not converged '
        END IF
        WRITE(6,*) ' Final value of linesearch parameter ', X
        WRITE(6,*) ' Function value at final point ', FFINAL
      END IF
      RETURN
      END 
      SUBROUTINE ORD_3DATA(X1IN,X2IN,X3IN,F1IN,F2IN,F3IN,
     &                       X1UT,X2UT,X3UT,F1UT,F2UT,F3UT)
*
* Three values of X and corresponding function values are given.
* Order so they come in increasing order of X
*
*. Jeppe Olsen ( I admit)
*
      INCLUDE 'implicit.inc'
*
      DIMENSION F(3),X(3)
*
      X(1) = X1IN
      F(1) = F1IN
      X(2) = X2IN
      F(2) = F2IN
      X(3) = X3IN
      F(3) = F3IN
*
      DO I = 1, 3
        DO J = I+1,3
          IF(X(I).GT.X(J)) THEN
           XI = X(I)
           FI = F(I)
           X(I) = X(J)
           F(I) = F(J)
           X(J) = XI
           F(J) = FI
          END IF
        END DO
      END DO
*
      X1UT = X(1)
      F1UT = F(1)
      X2UT = X(2)
      F2UT = F(2)
      X3UT = X(3)
      F3UT = F(3)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Sorting x, f(x) '
         WRITE(6,*) ' IN : X1, X2, X3, F(X1),F(X2),F(X3) = '
         WRITE(6,*)   X1IN,X2IN,X3IN,F1IN,F2IN,F3IN 
         WRITE(6,*) ' OUT : X1, X2, X3, F(X1),F(X2),F(X3) = '
         WRITE(6,*)   X1UT,X2UT,X3UT,F1UT,F2UT,F3UT
      END IF
*
      RETURN
      END
      SUBROUTINE ANA_3DATA(F1,F2,F3,ITYP)
*
*. Analyze three function values F1,F2,F3 corresponding 
*  to F(X1), F(X2), F(X3) with X1 < X2 < X3 to see if they
*
* 1 : constitute an bisection : F1 > F2 < F3, 
*     i.e a minimum exist between X1 and X3
* 2 : F1 < F2 : a minimum exist before X2
*
* 3 : F1 > F2 > F3 : a minimum exists after X2
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
*
      IF(F1.GT.F2.AND.F2.LT.F3) THEN
        ITYP = 1
      ELSE IF( F1.LT.F2) THEN
        ITYP = 2
      ELSE 
        ITYP = 3
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,3(1X,E22.15))') 
     &  ' Analyzing the three function values:',  F1,F2,F3
        IF(ITYP.EQ.1) WRITE(6,*) ' Defines bracketing '
        IF(ITYP.EQ.2) WRITE(6,*) ' Minimum before second value'
        IF(ITYP.EQ.3) WRITE(6,*) ' Minimum after second value'
      END IF
*
      RETURN
      END
      SUBROUTINE FIND_QUADPOL_FOR_3PT(A,B,C,X1,X2,X3,F1,F2,F3)
*
* Find a quadratic poly, that passes through the points FI = F(XI)
*
*. Jeppe Olsen, Oct. 2004
*
      INCLUDE 'implicit.inc'
*
* We start by writing the poly as 
* a(X-X1)(X-X2) + B'(X-X1) + c' so
      CP = F1
      BP = (F2-CP)/(X2-X1)
      AP = (F3-BP*(X3-X1)-CP)/((X3-X1)*(X3-X2))
*. And then rewrite to A* X² + B*X + C
      A = AP
      B = BP - A*(X1+X2)
      C = CP - BP*X1+ AP*X1*X2
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Quadratic poly fitting X1,X2,X3, F1,F2,F3 '
        WRITE(6,*)  X1,X2,X3,F1,F2,F3
        WRITE(6,*) 'Coefficients of poly, A,B,C ', A,B,C
      END IF
*
      RETURN
      END
      SUBROUTINE EXTREMUM_QUADPOL(A,B,C,XEXTR,IMNMX)
*
* A quadratic poly is given by A*x^2 + B*x + C
* Find extremum  and report on whether it is a min or max
*
* Jeppe Olsen
*
      INCLUDE 'implicit.inc'
*
      XEXTR = -B/(2.0D0*A)
      IF(A.GT.0.0D0) THEN
        IMNMX = 1
      ELSE
        IMNMX = 2
      END IF
*
      RETURN
      END
      SUBROUTINE REFINE_BISECTION(XMIN,XMID,XMAX,FMIN,FMID,FMAX,X,F)
*
* a bisection is given by the points XMIN, XMID, XMAX 
* and the function values FMIN,FMID,FMAX, FMIN > FMID < FMAX.
* The function has been evaluated at a new point, X with function 
* value F. Use this to refince bisection
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' REFINING bi-section '
        WRITE(6,*) ' New point, X, F = ', X,F
        WRITE(6,*) ' Old XMIN, XMID, XMAX, FMIN, FMID,FMAX'
        WRITE(6,*)   XMIN,XMID,XMAX,FMIN,FMID,FMAX
      END IF
*
      IF(X.LT.XMID) THEN
* XMIN X XMID XMAX
        IF(F.GT.FMID) THEN
*. The minimum is between X and XMAX, so 
          XMIN = X
          FMIN = F
        ELSE
*. The minimum is between XMIN and XMID
          XMAX = XMID
          FMAX = FMID
          XMID = X
          FMID = F
        END IF
      ELSE 
* XMIN XMID X MAX
        IF(F.GT.FMID) THEN
*. The minimum is between XMIN and X
          XMAX = X
          FMAX = F
        ELSE
*. The minimum is between XMID and XMAX
          XMIN = XMID
          FMIN = FMID
          XMID = X
          FMID = F
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' New XMIN, XMID, XMAX, FMIN, FMID,FMAX'
        WRITE(6,*)   XMIN,XMID,XMAX,FMIN,FMID,FMAX
      END IF
*
      RETURN
      END
      SUBROUTINE  TEST_BISECTION
*
* testing bisection routine
*
      INCLUDE 'implicit.inc'
*
       PARAMETER(MAXNVAR = 100)
       DIMENSION REF(MAXNVAR),DIR(MAXNVAR)
*
       EXTERNAL FUNKY
*. First test, a one-dimensional function 
      REF(1)  = 0.0D0
      DIR(1) = 1.0D0
      XINI = 0.1D0
      NVAR = 1
      RDUM = 0.0D0
      CALL LINES_SEARCH_BY_BISECTION(FUNKY,REF,DIR,NVAR,XINI,X,FX,
     %                               0, RDUM,RDUM)
C     LINES_SEARCH_BY_BISECTION(FUNC,REF,DIR,NVAR,XINI,XFINAL,FFINAL)
      STOP ' Enforced stop in TEST_BISECTION '
* 
      RETURN
      END
      FUNCTION FUNKY(X)
*
* F(X) = (X-1)(X-2)(X-3)(X-4)
*
      INCLUDE 'implicit.inc'
*
      FUNKY = (X-1.0D0)*(X-2.0D0)*(X-3.0D0)*(X-4.0D0)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'X and F(X) ', X, FUNKY
      END IF
*
      RETURN
      END
C       MIN_EHF_ALONG_MODE(WORK(KLEIGV),WORK(KLOOEXCC),NOOEXC,
C    &                          EMIN,ALPHAMIN)
      SUBROUTINE MIN_EHF_ALONG_MODE(XMODE,IOOXC,NOOXCX,EMIN,ALPHAMIN)
*
* Minimize HF energy along a mode given by XMODE
*
*. Jeppe Olsen, Oct. 2004
*
*. Last revision, Aug. 30 2012, Jeppe Olsen, updated to modern times
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
*
      EXTERNAL EHF_FROM_KAPPAM
*. Input
      DIMENSION XMODE(NOOXCX),IOOXC(2*NOOXCX)
*
      COMMON/CLOOXC/KLIOOXC,NOOXC
*
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from MIN_EHF_ALONG_MODE'
        WRITE(6,*) ' ============================'
        WRITE(6,*) 
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input mode '
        CALL WRTMAT(XMODE,1,NOOXCX,1,NOOXCX)
      END IF
*
      NOOXC = NOOXCX
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MNE_MO')
*
      CALL MEMMAN(KLXNOT,NOOXC,'ADDL  ',2,'LXNOT ')
      CALL MEMMAN(KLIOOXC,2*NOOXC,'ADDL  ',1,'LIOOXC')
*. As the bisection routine is not allowed to know anything but 
*. direction, info about orbital excitations must be transferred 
*. through a common block
      CALL ICOPVE(IOOXC,WORK(KLIOOXC),2*NOOXC)
*. We are expanding around Kappa = 0, so
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLXNOT),ZERO,NOOXC)
*. Start by obtaining the first three derivatives
C                         (EFUNC,XNOT,XMODE,NVAR,E1,E2,E3)
      CALL E123_ALONG_MODE(EHF_FROM_KAPPAM,WORK(KLXNOT),
     &       XMODE,NOOXC,E1,E2,E3)
*. We will use a rather small initial steplength to 
*. reduce change of overlooking minimum
      XINI = 0.1D0
      RDUM = 0.0D0
      CALL LINES_SEARCH_BY_BISECTION(EHF_FROM_KAPPAM,
     &     WORK(KLXNOT),XMODE,NOOXC,XINI,XFINAL,EMIN,
     &     0,RDUM,RDUM) 
      ALPHAMIN = XFINAL
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Results from search along mode '
        WRITE(6,*) ' Final step and energy ', ALPHAMIN,EMIN
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'MNE_MO')
      RETURN
      END 
      SUBROUTINE GET_E1_HF_ANALYTICAL(E1,FNOT,DSEED,DNOT,IAPPROX,NDENSI)
*
* Obtain Hartree-Fock gradient by analytical means at
* point of expansion 
* ( No more fooling around with finite difference jeppe)
*
* IAPPROX defines energy approximation 
* IAPPROX = 0 : E_APPROX = Tr FNOT D
* IAPPROX = 1 : E_APPROX = Tr(H D_proj) + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 2 : E_APPROX = Tr(H D ) +  Tr D_ort G(D_proj) 
*                        + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 3 : E_APPROX is exact energy of current density
*
* Jeppe Olsen, May 2005 - for the onestep project with Stinne
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*
*. Input
      DIMENSION FNOT(*), DSEED(*),DNOT(*) 
*. Output
      DIMENSION E1(*)
*
      NTEST = 00
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E1OSAN')
*
      LENC = LEN_AO_MAT(0)
      LENX = LEN_AO_MAT(-1)
*. Expand S to complete form 
*. 
      IF(IAPPROX.EQ.0.OR.IAPPROX.EQ.2) THEN
*.  Calculate gradient of E_epsilon = Tr FNOT D : 
*   -4 P(script)^T[(FNOT D_seed S)]^A_{\mu\vu} 
        CALL GET_E1_OF_EEPSILON(E1,FNOT,DSEED)
      END IF
*
      IF(IAPPROX.EQ.2) THEN
*. Add contribution  from E_extra to gradient 
        CALL MEMMAN(KLE1EXTRA,LENX,'ADDL  ',2,'E1_EXT')
C            GET_E1_OF_EEXTRA(E1,FNOT,DSEED,DNOT,NDENSI)
        CALL GET_E1_OF_EEXTRA(WORK(KLE1EXTRA),FNOT,DSEED,DNOT,NDENSI)
        ONE = 1.0D0
        CALL VECSUM(E1,E1,WORK(KLE1EXTRA),ONE,ONE,LENX)
      END IF
      IF(IAPPROX.EQ.3) THEN
*. Obtain Fock matrix corresponding to DSEED
        CALL MEMMAN(KLFOCK,LENP,'ADDL  ',2,'FOCK  ')
        XDUM = 0.0D0
        CALL GET_FOCK(WORK(KLFOCK),WORK(KHAO),DSEED,ENERGY,
     &                XDUM,3)
*.  Calculate gradient of E_epsilon = Tr F D : 
*   -4 P(script)^T[(F D_seed S)]^A_{\mu\vu} 
        CALL GET_E1_OF_EEPSILON(E1,WORK(KLFOCK),DSEED)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Gradient '
        CALL WRT_AOMATRIX(E1,-1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E1OSAN')
*
      RETURN
      END 
      SUBROUTINE GET_E1_OF_EEPSILON(E1,F,D)
*
* Obtain gradient of epsilon energy =  
*   -4 P(script)^T[(F D S)]^A_{\mu\vu}
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
*
*. Input F is in packed form whereas D is in complete form
*
*. Input
      DIMENSION F(*), D(*)
*. Output
      DIMENSION E1(*)
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' F matrix as delivered '
        CALL WRT_AOMATRIX(F,1)
      END IF
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E1_EEP')
* a : Expand F and S to complete form 
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLF_E,LENC,'ADDL  ',2,'FNOT_E')
      CALL TRIPAK_BLKM(WORK(KLF_E),F,2,NAOS_ENV,NSMOB)
      CALL MEMMAN(KLSAO_E,LENC,'ADDL  ',2,'SAO_E ')
      CALL TRIPAK_BLKM(WORK(KLSAO_E),WORK(KSAO),2,NAOS_ENV,NSMOB)
*. b : Obtain FNOT DSEED S
      CALL MEMMAN(KLSMAT1,LENC,'ADDL  ',2,'SCRMT1')
      CALL MEMMAN(KLSMAT2,LENC,'ADDL  ',2,'SCRMT2')
*.  b1 : F D in SMAT1
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' F_E, D, S_E '
        CALL WRT_AOMATRIX(WORK(KLF_E),0)
        CALL WRT_AOMATRIX(D,0)
        CALL WRT_AOMATRIX(WORK(KLSAO_E),0)
      END IF
      CALL MULT_BLOC_MAT(WORK(KLSMAT1),WORK(KLF_E),D,
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,0)
*. b2 : F D S in SMAT2
      CALL MULT_BLOC_MAT(WORK(KLSMAT2),WORK(KLSMAT1),WORK(KLSAO_E),
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,0)
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' F D S matrix '
         CALL WRT_AOMATRIX(WORK(KLSMAT2),0)
       END IF
*. c : obtain antisymmetric part of F D S in SMAT2
      CALL GET_SYMPART_OF_BLOCKED_MATRIX(WORK(KLSMAT2),NSMOB,NAOS_ENV,2)
*. d : Project with P(script) (T) and save in SMAT1
      CALL PROJ_NONRED(WORK(KLSMAT2),WORK(KLSMAT1),D,1)
C     PROJ_NONRED(X,XPROJ,D,ITRNSP)
*. and scale with -4
      FACTOR = -4.0D0
      CALL SCALVE(WORK(KLSMAT1),FACTOR,LENC)
*. And obtain the lower half
C     CALL REFORM_ANTISYM_BLMAT(X,WORK(KLXEXP),NAOS_ENV,NSMOB,0,1)
      CALL REFORM_ANTISYM_BLMAT(E1,WORK(KLSMAT1),NAOS_ENV,NSMOB,0,2)
C     REFORM_ANTISYM_BLMAT(APAK,ACOM,LBLK,NBLK,IDIAG_IS_IN_APAK,IWAY)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Gradient of E_epsilon : '
        WRITE(6,*) ' ========================'
        CALL PR_SYM_BLMAT_NODIAG(E1,NSMOB,NAOS_ENV)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E1_EEP')
*
      RETURN
      END
      SUBROUTINE PR_SYM_AOMAT_NODIAG(A)
*
* Print an blockes lower half AO matrix without an diagonal
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      CALL PR_SYM_BLMAT_NODIAG(A,NSMOB,NAOS_ENV)
*
      RETURN
      END
      SUBROUTINE PR_SYM_BLMAT_NODIAG(A,NBLOCK,LBLOCK)
*
*. Print a blocked lower half matrix without diagonal
*
*. Jeppe Olsen, May 2005 
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*),LBLOCK(NBLOCK)
*
      IOFFP = 1
      DO IBLOCK = 1, NBLOCK
        WRITE(6,*) ' Block ', IBLOCK
        LEN = LBLOCK(IBLOCK)-1
        CALL PRSYM(A(IOFFP),LEN)
        IOFFP = IOFFP + LEN*(LEN-1)/2
      END DO
*
      RETURN
      END
      SUBROUTINE GET_SYMPART_OF_AO_MATRIX(A,ISA)
*
* Obtain symmetric (ISA = 1) or antisymmetric (ISA = 2) part of 
* AO matrix A
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      CALL GET_SYMPART_OF_BLOCKED_MATRIX(A,NSMOB,NAOS_ENV,ISA)
*
      RETURN
      END

      SUBROUTINE GET_SYMPART_OF_BLOCKED_MATRIX(A,NBLOCK, LBLOCK,ISA)
*
* Obtain symmetric(ISA=1) or antisymmetric part(ISA=2) of blocked matrix
*. 
* Symmetric part     = 0.5 * (A + A(T))
* Antisymmetric part = 0.5 * (A - A(T))
*
      INCLUDE 'implicit.inc'
*. Matrix is in complete form
*. Input 
      INTEGER LBLOCK(NBLOCK)
*. Input and output
      DIMENSION A(*)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFF = 1
        ELSE 
          IOFF = IOFF + LBLOCK(IBLOCK-1)**2
        END IF 
        LENGTH = LBLOCK(IBLOCK)
        CALL GET_SYMPART_OF_MATRIX(A(IOFF),LENGTH,ISA)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(ISA.EQ.1) THEN
          WRITE(6,*) ' Symmetric part of matrix  : '
        ELSE IF (ISA.EQ.2) THEN
          WRITE(6,*) ' Antisymmetrix part of matrix '
        END IF
        CALL APRBLM2(A,LBLOCK,LBLOCK,NBLOCK,0)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_SYMPART_OF_MATRIX(A,LENGTH,ISA)
*
* Obtain symmetric(ISA = 1) or antisymmetric (ISA = 2) part of matrix
*
      INCLUDE 'implicit.inc'
      DIMENSION A(LENGTH,LENGTH)
*
      IF(ISA.EQ.1) THEN
        DO I = 1, LENGTH
         DO J = 1, I-1
           A(I,J) = 0.5D0*(A(I,J) + A(J,I))
           A(J,I) = A(I,J)
         END DO
        END DO
      ELSE IF (ISA.EQ.2) THEN
        DO I = 1, LENGTH
         DO J = 1, I-1
           A(I,J) = 0.5D0*(A(I,J) - A(J,I))
           A(J,I) = -A(I,J)
         END DO
         A(I,I) = 0.0D0
        END DO
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       IF(ISA.EQ.1) THEN 
         WRITE(6,*) ' Symmetric part of matrix '
       ELSE IF (ISA.EQ.2) THEN
         WRITE(6,*) ' Antisymmetric part of matrix '
       END IF
       CALL WRTMAT(A,LENGTH,LENGTH,LENGTH,LENGTH)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_E1_OF_EEXTRA(E1,FNOT,DSEED,DNOT,NDENSI)
*
* Obtain contribution from E_Extra to gradient 
*
* = -factor { \script P(T) [S DSEED (G(Dproj(0) + S(Dcc-0.5Dc)S)]^A_{\mu\vu} 
*
* where 
* Dc = sum_{ij} T_{ij}^{-1} Tr (D_i G(Dproj(0)) D_j
* DCc = sum_{ij} T_{ij}^{-1} Tr ((Dseed - DNOT - 0.5 Dproj(0)) G(D_i) ) D_j
* Dproj = sum_{ij} T_{ij}^{-1} D_i < D_j!Dseed-Dnot>
*
* and the factor depends on the normalization  of the density 
* and the G-term. If the density is the standard MO density
* with eigenvalues 2 and 0, the factor is 4. 
* In case the density is the idempotent density with eigenvalues 
* 0 and 1, the factor is 8.
*
* So : If we are using our standard book/paper normalization  factor = 8
*      If we are using the normalization used in this program factor = 4
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION FNOT(*), DSEED(*), DNOT(*)
*
C?    WRITE(6,*) ' Entering GET_E1_OF_EEXTRA '
      NTEST = 000
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E1_EEX')
C     LEN_AO_MAT(IPAK)
      LENC = LEN_AO_MAT(0)
*. Scratch matrices 
      CALL MEMMAN(KLSCR1,LENC,'ADDL  ',2,'SCR1  ')
      LENSCR = MAX(LENC,NDENSI**2)
      CALL MEMMAN(KLSCR2,LENSCR,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LENSCR,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSCR4,LENSCR,'ADDL  ',2,'SCR4  ')
*. Space for Dc and Dcc
      CALL MEMMAN(KLDC  ,LENC,'ADDL  ',2,'DC     ')
      CALL MEMMAN(KLDCC ,LENC,'ADDL  ',2,'DCC    ')
*
      CALL MEMMAN(KLSINV,NDENSI**2,'ADDL  ',2,'SINV  ')
      CALL MEMMAN(KLCPROJ,NDENSI,'ADDL  ',2,'CPROJ ')
      CALL MEMMAN(KLSVEC,NDENSI,'ADDL  ',2,'SVEC  ')
*. Expand S to complete form
      CALL MEMMAN(KLSAO_E,LENC,'ADDL  ',2,'SAO_E ')
C     TRIPAK_AO_MAT(AUTPAK,APAK,IWAY)
      CALL TRIPAK_AO_MAT(WORK(KLSAO_E),WORK(KSAO),2)
*
* ========================
* Dproj(0) and G(Dproj(0))
* ========================
*
*. Dseed - Dnot in SCR1
      ONE = 1.0D0
      ONEM = -1.0D0
      CALL VECSUM(WORK(KLSCR1),DSEED,DNOT,ONE,ONEM,LENC)
*. Projected density in SCR4, expansion coefficient in CPROJ
      CALL PROJECT_DENSI_ON_DENSI(WORK(KLSCR1),WORK(KDAO_COLLECT),LENC,
     &     NDENSI,WORK(KLSAO_E),WORK(KLSINV),WORK(KLSCR4),WORK(KLCPROJ),
     &     NSMOB,NTOOBS,WORK(KLSCR2),WORK(KLSCR3),WORK(KLSVEC),1,
     &     DORIG_NORM,DPROJ_NORM)
C     PROJECT_DENSI_ON_DENSI(DENSIIN,DENSI,LDENSI,NDENSI,
C    &           SAO,SINV,PROJ_DENSI,PROJ_COEF,NSMOB,NOBPSM,
C    &           SCR1,SCR2,SCRVEC,I_DO_SINV,DORIG_NORM, DPROJ_NORM)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Dseed - Dnot and Dproj '
        CALL WRT_AOMATRIX(WORK(KLSCR1),0)
        CALL WRT_AOMATRIX(WORK(KLSCR4),0)
      END IF
*. Obtain G(D_proj)
      CALL MEMMAN(KLGPROJ,LENC,'ADDL  ',2,'GPROJ ')
      CALL GET_G_DSUM(WORK(KLGPROJ),NDENSI,WORK(KLCPROJ))
*
* =======
*    Dc
* =======
*. Obtain Tr (D_i G(DPROJ)) in CPROJ
      DO I = 1, NDENSI
        KLDI = KDAO_COLLECT + (I-1)*LENC
        WORK(KLCPROJ-1+I) = 
     &  TRACE_PROD_AOMATRICES(WORK(KLDI),WORK(KLGPROJ))
C       TRACE_PROD_AOMATRICES(A,B)
      END DO
*. T^{-1} CPROJ in SVEC
      CALL MATVCB(WORK(KLSINV),WORK(KLCPROJ),WORK(KLSVEC),
     &            NDENSI,NDENSI,0)
*. Dc = sum_i svec(i) densi(i) in DC
      CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLSVEC),WORK(KLDC),
     &            LENC,NDENSI,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The Dc matrix '
        CALL WRT_AOMATRIX(WORK(KLDC),0)
      END IF
*
* =======
*   Dcc
* =======
*
* Dseed - Dnot - 0.5Dproj in SCR2
* (Dseed - Dnot is in SCR1, Dproj in SCR4)
      ONE = 1.0D0
      HALFM = -0.5D0
      CALL VECSUM(WORK(KLSCR2),WORK(KLSCR1),WORK(KLSCR4),ONE,HALFM,LENC)
*. Obtain Tr(Dseed - Dnot - 0.5Dproj)G(D_i)
      DO IDENSI = 1, NDENSI
*. Obtain G(D_i) in SCR1
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLSVEC),ZERO,NDENSI)
        WORK(KLSVEC-1+IDENSI) = 1.0D0
        CALL GET_G_DSUM(WORK(KLSCR1),NDENSI,WORK(KLSVEC))
C            GET_G_DSUM(G,NDENSI,CVEC)
*. Tr((Dseed - Dnot - 0.5Dproj)G(D_i) in CPROJ(i)
C     TRACE_PROD_AOMATRICES(A,B)
        WORK(KLCPROJ-1+IDENSI) = 
     &  TRACE_PROD_AOMATRICES(WORK(KLSCR2),WORK(KLSCR1))
      END DO
C
*. T^{-1}  CPROJ in SVEC
      CALL MATVCB(WORK(KLSINV),WORK(KLCPROJ),WORK(KLSVEC),
     &            NDENSI,NDENSI,0)
C
*. and Dcc = sum(i) svec(i) densi(i)
      CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLSVEC),WORK(KLDCC),
     &            LENC,NDENSI,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The Dcc matrix '
        CALL WRT_AOMATRIX(WORK(KLDCC),0)
      END IF
*
* =========================================
*. G(D(proj) + S (Dcc - 0.5 Dc ) S in SCR1
* =========================================
*
      CALL VECSUM(WORK(KLDCC),WORK(KLDCC),WORK(KLDC),ONE,HALFM,LENC)
      CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLSAO_E),WORK(KLDCC),0)
      CALL MULT_AO_MATRICES(WORK(KLSCR3),WORK(KLSCR2),WORK(KLSAO_E),0)
      CALL VECSUM(WORK(KLSCR1),WORK(KLGPROJ),WORK(KLSCR3),ONE,ONE,LENC)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' G(D(proj) + S (Dcc - 0.5 Dc ) S '
        CALL WRT_AOMATRIX(WORK(KLSCR1),0)
C            WRT_AOMATRIX(A,IPAK)
      END IF
* ====================================================
*. S Dseed ( G(D(proj) + S (Dcc - 0.5 Dc ) S) in SCR1
* ====================================================
      CALL MULT_AO_MATRICES(WORK(KLSCR2),DSEED,WORK(KLSCR1),0)
      CALL MULT_AO_MATRICES(WORK(KLSCR1),WORK(KLSAO_E),WORK(KLSCR2),0)
* ===============================================================================
*. And antisymmetrize :  [ S Dseed ( G(D(proj) + S (Dcc - 0.5 Dc ) S)]^A in SCR1
* ===============================================================================
C           GET_SYMPART_OF_AO_MATRIX(A,ISA)
      CALL  GET_SYMPART_OF_AO_MATRIX(WORK(KLSCR1),2)
* ====================================================
*. Project  with  P(script)^T and save in KLSCR2 
* ====================================================
C          PROJ_NONRED(X,XPROJ,D,ITRNSP)
      CALL PROJ_NONRED(WORK(KLSCR1),WORK(KLSCR2),DSEED,1)
*. Multiply with factor
      FACTOR = 4.0D0
      CALL SCALVE(WORK(KLSCR2),FACTOR,LENC)
*. And extract lower half - diagonal
C     REFORM_ANTISYM_AOMAT(APAK,ACOM,IDIAG_IS_IN_APAK,IWAY)
      CALL REFORM_ANTISYM_AOMAT(E1,WORK(KLSCR2),0,2)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Contribution to gradient from Eextra '
C            PR_SYM_AOMAT_NODIAG(A)
        CALL PR_SYM_AOMAT_NODIAG(E1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E1_EEX')
*
      RETURN
      END
      SUBROUTINE GET_G_DSUM(G,NDENSI,CVEC)
*
* A density is given as Dproj = sum_{i=1,ndensi} CVEC(i) Densi(i)
* Obtain the corresponding two-electron interaction 
* G(Dproj) = G( sum_{i=1,ndensi} CVEC(i) Densi(i)) from 
* the fock matrices calculated for the various densities
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
*. Input
      DIMENSION CVEC(NDENSI)
*. Output
      DIMENSION G(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_G ')
*
*. Scratch matrices 
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      CALL MEMMAN(KLSCR1,LENC,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LENC,'ADDL  ',2,'SCR2  ')
* Obtain sum_i CVEC(i) F(D_i) in SCR1
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      CALL MATVCC(WORK(KFAO_COLLECT),CVEC,WORK(KLSCR1),LENP,NDENSI,0)
*. Expand in SCR2
      CALL TRIPAK_BLKM(WORK(KLSCR2),WORK(KLSCR1),2,NTOOBS,NSMOB)
*. In SCR2 we now have sum(i) cvec(i) h + G(D_PROJ), but we 
*. wanted G(D_proj) so :
      SUM = 0.0D0
      DO I = 1, NDENSI
        SUM = SUM + CVEC(I)
      END DO
      COEF = - SUM
      ONE = 1.0D0
* h in expanded form 
      CALL MEMMAN(KLHAO_E,LENC,'ADDL  ',2,'HAO_E ')
      CALL TRIPAK_BLKM(WORK(KLHAO_E),WORK(KHAO),2,NAOS_ENV,NSMOB)
* G(Dproj) in G
      CALL VECSUM(G,WORK(KLSCR2),WORK(KLHAO_E),ONE,COEF,LENC)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_G ')
      RETURN
      END 
      FUNCTION TRACE_PROD_AOMATRICES(A,B)
*
*. Find trace of product of two complete AO matrices
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      DIMENSION A(*), B(*)
     
      TRACE = TRACE_PROD_BLMAT(A,B,NAOS_ENV,NSMOB)
      TRACE_PROD_AOMATRICES = TRACE
*
      RETURN
      END
      FUNCTION TRACE_PROD_BLMAT(A,B,LBLOCK,NBLOCK)
*
* Trace of product of two blocked complete matrices 
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*),B(*)
      INTEGER LBLOCK(NBLOCK)
*
      SUM = 0.0D0
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFF = 1
        ELSE 
          IOFF = IOFF + LBLOCK(IBLOCK-1)**2
        END IF
        LENGTH = LBLOCK(IBLOCK)
        SUM  = SUM + TRACE_PROD_MAT(A(IOFF),B(IOFF),LENGTH)
      END DO
*
      TRACE_PROD_BLMAT = SUM
*
      RETURN
      END 
      FUNCTION TRACE_PROD_MAT(A,B,NDIM)
*
* Trace of product of two square matrices 
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM)
*
* Trace = sum(i,j) A(i,j)*B(j,i)
*
      TRACE = 0.0D0
      DO I = 1, NDIM
       DO J = 1, NDIM
        TRACE = TRACE + A(I,J)*B(J,I)
       END DO
      END DO
*
      TRACE_PROD_MAT = TRACE
*
      RETURN
      END
      SUBROUTINE MULT_AO_MATRICES(C,A,B,ITRNSP)
*
* Multiply two complete and symmetry-blocked matrices 
*
*  Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION A(*), B(*)
*. Output
      DIMENSION C(*)
*
      CALL MULT_BLOC_MAT(C,A,B,
     &     NSMOB,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,NAOS_ENV,
     &     NAOS_ENV,ITRNSP)
*
      RETURN
      END
      SUBROUTINE REFORM_ANTISYM_AOMAT(APAK,ACOM,IDIAG_IS_IN_APAK,IWAY)
*
* Reform antisymmetric blocked AO matrix between full and lower 
* half form
*
* IWAY = 1 : Packed to complete form
* IWAY = 2 : complete to packed form
*
*. If IDIAG_IS_IN_APAK = 0, APAK does not contain the (vanishing) diagonal
*
* Jeppe Olsen, May 2005
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input and output
      DIMENSION APAK(*), ACOM(*)
*
      CALL REFORM_ANTISYM_BLMAT(APAK,ACOM,NAOS_ENV,NSMOB,
     &           IDIAG_IS_IN_APAK,IWAY)
*
      RETURN
      END
      FUNCTION LEN_AO_MAT(IPAK)
*
* Length of an packed (IPAK = 1) or unpacked (IPAK = 0) AO matrix with 
* symmetry 1
* IPAK = -1 => packed without diagonal
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      IABS_IPAK = ABS(IPAK)
      LEN_AO_MAT = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,IABS_IPAK)
      IF(IPAK.EQ.-1) LEN_AO_MAT = LEN_AO_MAT - NTOOB
*
      RETURN
      END
      SUBROUTINE TRIPAK_AO_MAT(AUTPAK,APAK,IWAY)
*
* Reformat AO matrix
*
* IWAY = 1 => Unpacked to packed
* IWAY = 2 => Packed to unpacked
*
*. Jeppe Olsen, May 2005
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input and output
      DIMENSION AUTPAK(*), APAK(*)
*
C          TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
      CALL TRIPAK_BLKM(AUTPAK,APAK,IWAY,NAOS_ENV,NSMOB)
*
      RETURN
      END 
      SUBROUTINE WRT_AOMATRIX(A,IPAK)
*
* Print AO matrix of symmetry 1
*
* IPAK = 0 => Not packed 
* IPAK = 1 => Packed 
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION A(*)
*
      IF(IPAK.GE.0) THEN
        CALL APRBLM2(A,NAOS_ENV,NAOS_ENV,NSMOB,IPAK)      
      ELSE IF (IPAK.EQ.-1) THEN 
       CALL PR_SYM_BLMAT_NODIAG(A,NSMOB,NAOS_ENV) 
      END IF
*
      RETURN
      END
      SUBROUTINE INPROD_DENSITIES_MAT(V,DENSI,LDENSI,NDENSI,A)
*
* Obtain inner products V_i = < A Densi_i> = TR (S A S Densi_i), i=1, ndensi
* Densities and A are assumed delivered as complete, symmetrypacked matrices
*
*. Jeppe Olsen, still May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION DENSI(LDENSI,NDENSI),A(*)
*. Output
      DIMENSION V(NDENSI)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'INPDNM')
      LENC = LEN_AO_MAT(0)
*. Space for unpacked S and scratch matrix 
      CALL MEMMAN(KLSCR1,LENC,'ADDL  ',2,'SCR1  ')
      CALL MEMMAN(KLSCR2,LENC,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LENC,'ADDL  ',2,'SCR3  ')
*. Obtain unpacked S matrix in SCR1
      CALL TRIPAK_AO_MAT(WORK(KLSCR1),WORK(KSAO),2)
*. Calculate S A in SCR2
C MULT_AO_MATRICES(C,A,B,ITRNSP)
      CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLSCR1),A,0)
*. Calculate S A S in SCR3
      CALL MULT_AO_MATRICES(WORK(KLSCR3),WORK(KLSCR2),WORK(KLSCR1),0)
*. And then the inner product
      DO IDENSI = 1, NDENSI
C  TRACE_PROD_AOMATRICES(A,B)
       V(IDENSI) = TRACE_PROD_AOMATRICES(DENSI(1,IDENSI),WORK(KLSCR3))
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' <D_i | A > : '
        WRITE(6,*) ' ============ '
        CALL WRTMAT(V,1,NDENSI,1,NDENSI)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'INPDNM')
      RETURN
      END
      SUBROUTINE E2_TV_ONESTEP(E2TV,X,FNOT,DSEED,DNOT,IAPPROX,NDENSI)
*
* Obtain E2 times vector for onestep approach
*
*
* IAPPROX defines energy approximation 
* IAPPROX = 0 : E_APPROX = Tr FNOT D
* IAPPROX = 1 : E_APPROX = Tr(H D_proj) + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 2 : E_APPROX = Tr(H D ) +  Tr D_ort G(D_proj) 
*                        + 0.5* Tr D_proj G(D_proj)
* IAPPROX = 3 : E_APPROX is exact energy of current density
*
* Jeppe Olsen, May 2005 - for the onestep project with Stinne
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
*. Input : X is in packed form
      DIMENSION FNOT(*), DSEED(*),DNOT(*),X(*)
*. Output
      DIMENSION E2TV(*)
*
*
* E2 X = 
*   4 P(T) [(F + G(Dproj) + S (Dcc-0.5Dcc) S) P(X) S Dseed S]^A_{\mu\nu}
* - 4 P(T) [[(F + G(Dproj) + S (Dcc-0.5Dcc) S) Dseed S]^S P(X) S ]^A_{\mu\nu}
* + 4 P(T) [S Dseed (G(Dproj(1)) + S (Dcc(1) -0.5 Dc(1) ) S)]^A_{\mu\nu}
*
* where 
* Dc = sum_{ij} T_{ij}^{-1} Tr (D_i G(Dproj(0)) D_j
* Dcc = sum_{ij} T_{ij}^{-1} Tr ((Dseed - DNOT - 0.5 Dproj) G(D_i) ) D_j
* Dproj = sum_{ij} T_{ij}^{-1} D_i < D_j!Dseed-Dnot>
*
* Dproj(1) = sum_{ij} T_{ij}^{-1} D_i <D_j[Dseed,P(X)]_S>
* Dc(1) = sum_{ij} T_{ij}^{-1} Tr (D_i G(Dproj(1)) D_j
* Dcc(1) = sum_{ij} T_{ij}^{-1} Tr ([Dseed,P(X)]_S - 0.5Dproj(1))G(D_i)
*
* If IAPPROX = 0, only the terms containing F should be calculated
*
      NTEST = 00
      IF(NTEST.GE.10) WRITE(6,*) '  E2_TV_ONESTEP entered '
*
      ONE = 1.0D0
      HALFM = -0.5D0
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E2_TV ')
*
      LENC = LEN_AO_MAT(0)
      LENX = LEN_AO_MAT(-1)
*. Scratch matrices 
      CALL MEMMAN(KLSCR1,LENC,'ADDL  ',2,'SCR1  ')
      LENSCR = MAX(LENC,NDENSI**2)
      CALL MEMMAN(KLSCR2,LENSCR,'ADDL  ',2,'SCR2  ')
      CALL MEMMAN(KLSCR3,LENSCR,'ADDL  ',2,'SCR3  ')
      CALL MEMMAN(KLSCR4,LENSCR,'ADDL  ',2,'SCR4  ')
      CALL MEMMAN(KLSCR5,LENSCR,'ADDL  ',2,'SCR4  ')
*. Space for Dc and Dcc
      CALL MEMMAN(KLDC  ,LENC,'ADDL  ',2,'DC     ')
      CALL MEMMAN(KLDCC ,LENC,'ADDL  ',2,'DCC    ')
*
      CALL MEMMAN(KLSINV,NDENSI**2,'ADDL  ',2,'SINV  ')
      CALL MEMMAN(KLCPROJ,NDENSI,'ADDL  ',2,'CPROJ ')
      CALL MEMMAN(KLSVEC,NDENSI,'ADDL  ',2,'SVEC  ')
*.
      CALL MEMMAN(KLS_E,LENC,'ADDL  ',2,'SAO_E ')
      CALL MEMMAN(KLX_E,LENC,'ADDL  ',2,'X_E ')
      CALL MEMMAN(KLF_E,LENC,'ADDL  ',2,'F_E ')
*. Expand S to complete form
C     TRIPAK_AO_MAT(AUTPAK,APAK,IWAY)
      CALL TRIPAK_AO_MAT(WORK(KLS_E),WORK(KSAO),2)
*. Expand F to complete form 
      CALL TRIPAK_AO_MAT(WORK(KLF_E),FNOT,2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Expanded Fock matrix '
        CALL WRT_AOMATRIX(WORK(KLF_E),0)
      END IF
*. Expand X to complete form
C          REFORM_ANTISYM_AOMAT(APAK,ACOM,IDIAG_IS_IN_APAK,IWAY)
      CALL REFORM_ANTISYM_AOMAT(X,WORK(KLX_E),0,1)
*. And project with P
C     PROJ_NONRED(X,XPROJ,D,ITRNSP)
      CALL PROJ_NONRED(WORK(KLX_E),WORK(KLSCR1),DSEED,0)
      CALL COPVEC(WORK(KLSCR1),WORK(KLX_E),LENC)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' P(X) in expanded form '
        CALL WRT_AOMATRIX(WORK(KLX_E),0)
      END IF
*
      IF(IAPPROX.NE.0.AND.IAPPROX.NE.3) THEN
*
* ========================
* Dproj(0) and G(Dproj(0))
* ========================
*
*. Dseed - Dnot in SCR1
       ONE = 1.0D0
       ONEM = -1.0D0
       CALL VECSUM(WORK(KLSCR1),DSEED,DNOT,ONE,ONEM,LENC)
*. Projected density in SCR4, expansion coefficient in CPROJ
       CALL PROJECT_DENSI_ON_DENSI(WORK(KLSCR1),WORK(KDAO_COLLECT),LENC,
     &      NDENSI,WORK(KLS_E),WORK(KLSINV),WORK(KLSCR4),
     &      WORK(KLCPROJ),
     &      NSMOB,NTOOBS,WORK(KLSCR2),WORK(KLSCR3),WORK(KLSVEC),1,
     &      DORIG_NORM,DPROJ_NORM)
C      PROJECT_DENSI_ON_DENSI(DENSIIN,DENSI,LDENSI,NDENSI,
C    &            SAO,SINV,PROJ_DENSI,PROJ_COEF,NSMOB,NOBPSM,
C    &            SCR1,SCR2,SCRVEC,I_DO_SINV,DORIG_NORM, DPROJ_NORM)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Dseed - Dnot and Dproj '
         CALL WRT_AOMATRIX(WORK(KLSCR1),0)
         CALL WRT_AOMATRIX(WORK(KLSCR4),0)
       END IF
*. Obtain G(D_proj)
       CALL MEMMAN(KLGPROJ,LENC,'ADDL  ',2,'GPROJ ')
       CALL GET_G_DSUM(WORK(KLGPROJ),NDENSI,WORK(KLCPROJ))
*
* =======
*    Dc
* =======
*. Obtain Tr (D_i G(DPROJ)) in CPROJ
       DO I = 1, NDENSI
         KLDI = KDAO_COLLECT + (I-1)*LENC
         WORK(KLCPROJ-1+I) = 
     &   TRACE_PROD_AOMATRICES(WORK(KLDI),WORK(KLGPROJ))
C        TRACE_PROD_AOMATRICES(A,B)
       END DO
*. T^{-1} CPROJ in SVEC
       CALL MATVCB(WORK(KLSINV),WORK(KLCPROJ),WORK(KLSVEC),
     &             NDENSI,NDENSI,0)
*. Dc = sum_i svec(i) densi(i) in DC
       CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLSVEC),WORK(KLDC),
     &             LENC,NDENSI,0)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' The Dc matrix '
         CALL WRT_AOMATRIX(WORK(KLDC),0)
       END IF
*
* =======
*   Dcc
* =======
*
* Dseed - Dnot - 0.5Dproj in SCR2
* (Dseed - Dnot is in SCR1, Dproj in SCR4)
       ONE = 1.0D0
       HALFM = -0.5D0
       CALL VECSUM(WORK(KLSCR2),WORK(KLSCR1),WORK(KLSCR4),ONE,HALFM,
     &             LENC)
*. Obtain Tr(Dseed - Dnot - 0.5Dproj)G(D_i)
       DO IDENSI = 1, NDENSI
*. Obtain G(D_i) in SCR1
         ZERO = 0.0D0
         CALL SETVEC(WORK(KLSVEC),ZERO,NDENSI)
         WORK(KLSVEC-1+IDENSI) = 1.0D0
         CALL GET_G_DSUM(WORK(KLSCR1),NDENSI,WORK(KLSVEC))
C             GET_G_DSUM(G,NDENSI,CVEC)
*. Tr((Dseed - Dnot - 0.5Dproj)G(D_i) in CPROJ(i)
C      TRACE_PROD_AOMATRICES(A,B)
         WORK(KLCPROJ-1+IDENSI) = 
     &   TRACE_PROD_AOMATRICES(WORK(KLSCR2),WORK(KLSCR1))
       END DO
C
*. T^{-1}  CPROJ in SVEC
       CALL MATVCB(WORK(KLSINV),WORK(KLCPROJ),WORK(KLSVEC),
     &             NDENSI,NDENSI,0)
C
*. and Dcc = sum(i) svec(i) densi(i)
       CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLSVEC),WORK(KLDCC),
     &             LENC,NDENSI,0)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' The Dcc matrix '
         CALL WRT_AOMATRIX(WORK(KLDCC),0)
       END IF

*
* =============================================
*. F + G(D(proj) + S (Dcc - 0.5 Dc ) S in SCR1
* =============================================
*

       CALL VECSUM(WORK(KLDCC),WORK(KLDCC),WORK(KLDC),ONE,HALFM,LENC)
       CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLS_E),WORK(KLDCC),0)
       CALL MULT_AO_MATRICES(WORK(KLSCR3),WORK(KLSCR2),WORK(KLS_E),0)
       CALL VECSUM(WORK(KLSCR1),WORK(KLGPROJ),WORK(KLSCR3),ONE,ONE,LENC)
       CALL VECSUM(WORK(KLSCR1),WORK(KLSCR1),WORK(KLF_E),ONE,ONE,LENC)
      ELSE IF (IAPPROX.EQ.0) THEN
*. Just transfer F to SCR1
       CALL COPVEC(WORK(KLF_E),WORK(KLSCR1),LENC)
      END IF
*     ^ End if IAPPROX ne, 0, 3
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' F + G(D(proj) + S (Dcc - 0.5 Dc ) S '
         CALL WRT_AOMATRIX(WORK(KLSCR1),0)
C             WRT_AOMATRIX(A,IPAK)
       END IF
*. 
*
*. =================================================================
*. 4 [(F + G(D(proj) + S (Dcc - 0.5 Dc ) S) X S Dseed S]^A in SCR4
*. =================================================================
*
      CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLSCR1),WORK(KLX_E),0)
      CALL MULT_AO_MATRICES(WORK(KLSCR3),WORK(KLSCR2),WORK(KLS_E),0)
      CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLSCR3),DSEED,0)
      CALL MULT_AO_MATRICES(WORK(KLSCR4),WORK(KLSCR2),WORK(KLS_E),0) 
      CALL GET_SYMPART_OF_AO_MATRIX(WORK(KLSCR4),2)
      FACTOR = 4.0D0
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '(F+G(D(proj)+S(Dcc-0.5 Dc)S)XSDseed S]^A'
        CALL WRT_AOMATRIX(WORK(KLSCR4),0)
      END IF
      CALL SCALVE(WORK(KLSCR4),FACTOR,LENC)
*
*. ===================================================================
*.-4 [[(F + G(D(proj) + S (Dcc - 0.5 Dc ) S) Dseed S]^S X S]^A in SCR3
*. ===================================================================
      CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLSCR1),DSEED,0)
      CALL MULT_AO_MATRICES(WORK(KLSCR3),WORK(KLSCR2),WORK(KLS_E),0)
      CALL GET_SYMPART_OF_AO_MATRIX(WORK(KLSCR3),1)
      CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLSCR3),WORK(KLX_E),0)
      CALL MULT_AO_MATRICES(WORK(KLSCR3),WORK(KLSCR2),WORK(KLS_E),0)
      CALL GET_SYMPART_OF_AO_MATRIX(WORK(KLSCR3),2)
      FACTOR = -4.0D0
      CALL SCALVE(WORK(KLSCR3),FACTOR,LENC)
*. and add term 1 and term 2 in SCR5
      CALL VECSUM(WORK(KLSCR5),WORK(KLSCR4),WORK(KLSCR3),ONE,ONE,LENC)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Term 1 + 2 : '
        CALL WRT_AOMATRIX(WORK(KLSCR5),0)
      END IF
*
      IF(IAPPROX.EQ.2) THEN
*
* =========================
*  Dproj(1) and G(Dproj(1))
* =========================
*
* Dproj(1) = sum_{ij} T_{ij}^{-1} D_i <D_j[Dseed,X]_S>
*
*. [Dseed,X]_S = 2 [Dseed S X]^S in  SCR1
       CALL MULT_AO_MATRICES(WORK(KLSCR2),DSEED,WORK(KLS_E),0)
       CALL MULT_AO_MATRICES(WORK(KLSCR1),WORK(KLSCR2),WORK(KLX_E),0)
       CALL GET_SYMPART_OF_AO_MATRIX(WORK(KLSCR1),1)
       FACTOR = 2.0D0
       CALL SCALVE(WORK(KLSCR1),FACTOR,LENC)
*. Projected density in SCR4, expansion coefficient in CPROJ
       CALL PROJECT_DENSI_ON_DENSI(WORK(KLSCR1),WORK(KDAO_COLLECT),LENC,
     &      NDENSI,WORK(KLS_E),WORK(KLSINV),WORK(KLSCR4),
     &      WORK(KLCPROJ),
     &      NSMOB,NTOOBS,WORK(KLSCR2),WORK(KLSCR3),WORK(KLSVEC),1,
     &      DORIG_NORM,DPROJ_NORM)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' D(1) and Dproj(1) '
         CALL WRT_AOMATRIX(WORK(KLSCR1),0)
         CALL WRT_AOMATRIX(WORK(KLSCR4),0)
       END IF
*. Obtain G(D_proj(1))
       CALL GET_G_DSUM(WORK(KLGPROJ),NDENSI,WORK(KLCPROJ))
*
* ===========
*    Dc(1) 
* ===========
*
* Dc(1) = sum_{ij} T_{ij}^{-1} Tr (D_i G(Dproj(1)) D_j
*. Obtain Tr (D_i G(DPROJ(1))) in CPROJ
       DO I = 1, NDENSI
         KLDI = KDAO_COLLECT + (I-1)*LENC
         WORK(KLCPROJ-1+I) = 
     &   TRACE_PROD_AOMATRICES(WORK(KLDI),WORK(KLGPROJ))
C        TRACE_PROD_AOMATRICES(A,B)
       END DO
*. T^{-1} CPROJ in SVEC
       CALL MATVCB(WORK(KLSINV),WORK(KLCPROJ),WORK(KLSVEC),
     &             NDENSI,NDENSI,0)
*. Dc(1) = sum_i svec(i) densi(i) in DC
       CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLSVEC),WORK(KLDC),
     &             LENC,NDENSI,0)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' The Dc(1) matrix '
         CALL WRT_AOMATRIX(WORK(KLDC),0)
       END IF
*.
*. =======
*.  Dcc(1)
*. =======
*.
* Dcc(1) = sum_{ij} T_{ij}^{-1} Tr ([Dseed,X]_S - 0.5Dproj(1))G(D_i)
* [Dseed,X]_S - 0.5Dproj(1) in SCR2
* ([Dseed,X]_S is in SCR1, Dproj(1) in SCR4)
       ONE = 1.0D0
       HALFM = -0.5D0
       CALL VECSUM(WORK(KLSCR2),WORK(KLSCR1),WORK(KLSCR4),ONE,HALFM,
     &             LENC)
*. Obtain Tr([Dseed,X]_S - 0.5Dproj(1))G(D_i)
       DO IDENSI = 1, NDENSI
*. Obtain G(D_i) in SCR1
         ZERO = 0.0D0
         CALL SETVEC(WORK(KLSVEC),ZERO,NDENSI)
         WORK(KLSVEC-1+IDENSI) = 1.0D0
         CALL GET_G_DSUM(WORK(KLSCR1),NDENSI,WORK(KLSVEC))
C             GET_G_DSUM(G,NDENSI,CVEC)
*. Tr(([Dseed,X]_S - 0.5Dproj(1))G(D_i) in CPROJ(i)
C      TRACE_PROD_AOMATRICES(A,B)
         WORK(KLCPROJ-1+IDENSI) = 
     &   TRACE_PROD_AOMATRICES(WORK(KLSCR2),WORK(KLSCR1))
       END DO
*. T^{-1}  CPROJ in SVEC
       CALL MATVCB(WORK(KLSINV),WORK(KLCPROJ),WORK(KLSVEC),
     &             NDENSI,NDENSI,0)
*. and Dcc(1) = sum(i) svec(i) densi(i)
       CALL MATVCC(WORK(KDAO_COLLECT),WORK(KLSVEC),WORK(KLDCC),
     &             LENC,NDENSI,0)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' The Dcc(1) matrix '
         CALL WRT_AOMATRIX(WORK(KLDCC),0)
       END IF
*
* ==================================================
*. G(D(proj(1)) + S (Dcc(1) - 0.5 Dc(1) ) S in SCR1
* ==================================================
*
       CALL VECSUM(WORK(KLDCC),WORK(KLDCC),WORK(KLDC),ONE,HALFM,LENC)
       CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLS_E),WORK(KLDCC),0)
       CALL MULT_AO_MATRICES(WORK(KLSCR3),WORK(KLSCR2),WORK(KLS_E),0)
       CALL VECSUM(WORK(KLSCR1),WORK(KLGPROJ),WORK(KLSCR3),ONE,ONE,LENC)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' G(D(proj(1)) + S (Dcc(1) - 0.5 Dc(1) ) S '
         CALL WRT_AOMATRIX(WORK(KLSCR1),0)
       END IF
*
* ===============================================================
*. 4 [ S Dseed G(D(proj(1)) + S (Dcc(1) - 0.5 Dc(1) ) S]^A in SCR1
* ===============================================================
*
       CALL MULT_AO_MATRICES(WORK(KLSCR2),DSEED,WORK(KLSCR1),0)
       CALL MULT_AO_MATRICES(WORK(KLSCR1),WORK(KLS_E),WORK(KLSCR2),0)
       CALL GET_SYMPART_OF_AO_MATRIX(WORK(KLSCR1),2)
       FACTOR = 4.0D0
       CALL SCALVE(WORK(KLSCR1),FACTOR,LENC)
*
* ======================================
* Add terms 1,2,3 and project with P(T) and save in KLSCR1
* ======================================
*
       CALL VECSUM(WORK(KLSCR5),WORK(KLSCR1),WORK(KLSCR5),ONE,ONE,LENC)
      END IF
*     ^ End if IAPPROX = 2
       CALL PROJ_NONRED(WORK(KLSCR5),WORK(KLSCR1),DSEED,1)
C          PROJ_NONRED(X,XPROJ,D,ITRNSP)
* ===================================
*. And extract lower half - diagonal
* ===================================
C     REFORM_ANTISYM_AOMAT(APAK,ACOM,IDIAG_IS_IN_APAK,IWAY)
      CALL REFORM_ANTISYM_AOMAT(E2TV,WORK(KLSCR1),0,2)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' E2 times vector '
        CALL PR_SYM_AOMAT_NODIAG(E2TV)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E2_TV ')
*
      WRITE(6,*) '  E2_TV_ONESTEP will be exited '
      RETURN
      END
      SUBROUTINE GET_E2_ONESTEP(E2,LENX,FNOT,DSEED,DNOT,NDENSI,IAPPROX)
*
* Obtain Hessian for onestep approach
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'wrkspc.inc'
*. Input
      DIMENSION FNOT(*), DNOT(*), DSEED(*)
*. Output
      DIMENSION E2(LENX,LENX)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E2_ONE')
*
      CALL MEMMAN(KLX,LENX,'ADDL  ',2,'X     ')
*
      DO I = 1, LENX
        ZERO = 0.0D0
        CALL SETVEC(WORK(KLX),ZERO,LENX)
        WORK(KLX-1+I) = 1.0D0
C            E2_TV_ONESTEP(E2TV,X,FNOT,DSEED,DNOT,IAPPROX,NDENSI)
        CALL E2_TV_ONESTEP(E2(1,I),WORK(KLX),FNOT,DSEED,DNOT,
     &                     IAPPROX,NDENSI)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Obtained Hessian'
        CALL WRTMAT(E2,LENX,LENX,LENX,LENX)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E2_ONE')
      RETURN
      END 
      SUBROUTINE CHECK_IDEMP(D)
*
*. A density matrix D is given as a complete matrix 
*. Check idempotency condition 
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*
      REAL*8 INPROD
*. Input
      DIMENSION D(*)
*
*. Factor for putting D into form required for idempotency 
* ( = 1/2 for orbital density, 1 for spin-orbital density)
*
      FACTOR_IDEM = 0.5D0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CHK_D ')
*.
      LENC = LEN_AO_MAT(0)
      CALL MEMMAN(KLS_E,LENC,'ADDL  ',2,'S_E    ')
      CALL MEMMAN(KLSCR1,LENC,'ADDL  ',2,'LSCR1  ')
      CALL MEMMAN(KLSCR2,LENC,'ADDL  ',2,'LSCR2  ')
*. Expand S to complete form 
      CALL TRIPAK_AO_MAT(WORK(KLS_E),WORK(KSAO),2)
*. D S in SCR1
      CALL MULT_AO_MATRICES(WORK(KLSCR1),D,WORK(KLS_E),0)
*. DS D in SCR2
      CALL MULT_AO_MATRICES(WORK(KLSCR2),WORK(KLSCR1),D,0)
*. Scaled DSD - D in SCR1
      FAC_DSD = FACTOR_IDEM**2
      FAC_D =   -FACTOR_IDEM
      CALL VECSUM(WORK(KLSCR1),WORK(KLSCR2),D,
     &            FAC_DSD,FAC_D,LENC)
      XNORM = INPROD(WORK(KLSCR1),WORK(KLSCR1),LENC)
*
      WRITE(6,*) ' Norm of DSD - D = ',SQRT(XNORM)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CHK_D ')
      RETURN
      END 
      SUBROUTINE PROJ_PCKMAT(X,XPROJ,D,IMATFORM,ITRNSP)
*
* Project matrix using P(script) (ITRNSP=0), P(script)(t) (ITRNSP=1)
*
* IMATFORM determines form of input and output matrix :
* IMATFORM = 1  Lower packed symmetric matrix 
* IMATFORM = -1 Lower packed antisymmetrix matrix without diagonal
* IMATFORM = 0 : Complete matrix
*
*. Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*. Input
      DIMENSION X(*),D(*)
*. Output
      DIMENSION XPROJ(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PROJPM')
*. Space for unpacked unprojected and projected matrix 
      LENC = LEN_AO_MAT(0)
      CALL MEMMAN(KLM     ,LENC,'ADDL  ',2,'M     ')
      CALL MEMMAN(KLM_PROJ,LENC,'ADDL  ',2,'M_PROJ ')
*. Unpack intial matrix 
      IF(IMATFORM.EQ.-1) THEN
C          REFORM_ANTISYM_AOMAT(APAK,ACOM,IDIAG_IS_IN_APAK,IWAY)
        CALL REFORM_ANTISYM_AOMAT(X,WORK(KLM),0,1)
      ELSE IF(IMATFORM.EQ.1) THEN
        CALL TRIPAK_AO_MAT(WORK(KLM),X,2)
      ELSE IF(IMATFORM.EQ.0) THEN
        CALL COPVEC(X,WORK(KLM),LENC)
      END IF
*. And project unpacked matrix 
      CALL PROJ_NONRED(WORK(KLM),WORK(KLM_PROJ),D,ITRNSP)
*. And pack again 
      IF(IMATFORM.EQ.-1) THEN
        CALL REFORM_ANTISYM_AOMAT(XPROJ,WORK(KLM_PROJ),0,2)
      ELSE IF(IMATFORM.EQ.1) THEN
        CALL TRIPAK_AO_MAT(WORK(KLM_PROJ),XPROJ,1)
      ELSE IF (IMATFORM.EQ.0) THEN
        CALL COPVEC(WORK(KLM_PROJ),XPROJ,LENC)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PROJPM')
      RETURN
      END 
      SUBROUTINE CHECK_PROJ_VEC(X,IMATFORM,D,ITRNSP)
*
* Check whether a given matrix VEC is invariant to
* Projection with P(script) (ITRNSP = 0) or P(script)(T) ( ITRNSP=1)
* i.e. find the norm of 
*    P(script) matrix - matrix  (ITRNSP = 0)
*    P(script)(T) matrix - matrix (ITRNSP = 1)
*
* IMATFORM = 0 : Full matrix form 
* IMATFORM = 1 : Packed form, symmetric matrix with diagonal
* IMATFORM =-1 : packed form, antisymmetric matrix without diagonal
*
*     Jeppe Olsen, May 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      REAL*8 INPROD
*. Input
      DIMENSION X(*), D(*)
*. 
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CHK_PV')
*
      LENC = LEN_AO_MAT(0)
      CALL MEMMAN(KLXPROJ,LENC,'ADDL  ',2,'XPROJ ')
C     PROJ_PCKMAT(X,XPROJ,D,IMATFORM,ITRNSP)
      CALL PROJ_PCKMAT(X,WORK(KLXPROJ),D,IMATFORM,ITRNSP)
*
      LENGTH = LEN_AO_MAT(IMATFORM)
      CALL VECSUM(WORK(KLXPROJ),WORK(KLXPROJ),X,1.0D0,-1.0D0,LENGTH)
*
      DIFF_NORM = INPROD(WORK(KLXPROJ),WORK(KLXPROJ),LENGTH)
      X_NORM    = INPROD(X,X,LENGTH)
*
      IF(ITRNSP.EQ.0) THEN
        WRITE(6,*) ' Norm of X and P X - X ', X_NORM, DIFF_NORM
      ELSE IF (ITRNSP.EQ.1) THEN
        WRITE(6,*) ' Norm of X and P(T) X - X ', X_NORM, DIFF_NORM
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CHK_PV')
      RETURN
      END
      SUBROUTINE DOVLAP(VOVLP,DTARGET,DOTHER,NDENSI,ITRANS,CTRANS,
     &                  NDIMC)
*
* A density matrix Dtarget is given. Obtain overlap 
* between this density and a set of densities in DOTHER.
* If ITRANS = 1, then the overlap with transformed densities 
* (transformation defined by CTRANS ) are obtained
*
*. Jeppe Olsen, June 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
*. Input 
      DIMENSION DTARGET(*), DOTHER(*),CTRANS(NDIMC,NDIMC)
*. Output
      DIMENSION VOVLP(NDENSI)
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'DOVLAP ')
*
      LENC = LEN_AO_MAT(0)
      CALL MEMMAN(KLOVLAP,NDENSI,'ADDL  ',2,'OVLAP ')
      CALL INPROD_DENSITIES_MAT(WORK(KLOVLAP),DOTHER,LENC,NDENSI,
     &                          DTARGET)
C   DENSITIES_MAT(V,DENSI,LDENSI,NDENSI,A)
      IF(ITRANS.EQ.1) THEN
        CALL MATVCB(CTRANS,WORK(KLOVLAP),VOVLP,NDIMC,NDENSI,1)
      ELSE 
        CALL COPVEC(WORK(KLOVLAP),VOVLP,NDENSI)
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100)  THEN
        WRITE(6,*) ' Overlap between target and other densities '
        CALL WRTMAT(VOVLP,1,NDENSI,1,NDENSI)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'DOVLAP ')
      RETURN
      END
      SUBROUTINE ANA_ORBS_IN_IT(F,D,NMAT)
*
* A set of densities and fock-matrices are given.
* Obtain canonical occupied orbitals and analyze these
*
*. Jeppe Olsen, June 2005
*
*. Canonical occupied orbitals in it I are obtained by 
*. diagonaling F(I) in the occupied orbitals from D(I).
*. In RH iterations this corresponds to having as matrix I
*  both the Fockmatrix, and the density obtained by diagonalizing this.
*. In RH iterations this corresponds to having as matrix I
*  both the Fockmatrix, and the density obtained by diagonalizing this.
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
*
*. The densitymatrices are in complete symmetry-blocked form 
*. and the Fock matrices are in lower-half symmetry packed form 
      DIMENSION F(*),D(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ANA_OI')
*
* Part I : obtain the canonical occupied orbitals in the 
*          various iterations
*
      LENC = LEN_AO_MAT(0)
      LENP = LEN_AO_MAT(1)
*. Total number of occupied orbitals 
      NOCOBT = IELSUM(NOCOBS,NSMOB)
* The symmetry of the occupied orbitals may vary in each 
* iteration, so we only know the total number of occupied 
* orbitals
      CALL MEMMAN(KLCOCC,LENC*NMAT,'ADDL  ',2,'COCC  ')
      CALL MEMMAN(KLSCR ,LENC,'ADDL  ',2,'SCR   ')
      CALL MEMMAN(KLNOCCFSM,NMAT*NSMOB,'ADDL  ',2,'NOCCFS')
*
      DO IMAT = 1, NMAT
         IB_D = (IMAT-1)*LENC + 1
         IB_F = (IMAT-1)*LENP + 1
         IB_OCC = KLCOCC + (IMAT-1)*LENC 
         IB_NOCC = KLNOCCFSM + (IMAT-1)*NSMOB 
C             GET_OCC_FOR_D(F,D,S,OCC,NOCC,IWAY)
         CALL GET_OCC_FOR_D(F(IB_F),D(IB_D),WORK(KSAO),
     &        WORK(IB_OCC),WORK(IB_NOCC),2)
*. We now have the various occupied orbitals in WORK(IB_OCC)
      END DO
* 1 : Analyze changes of orbitals in each iteration 
      CALL MEMMAN(KLCOSCN,LENC*NMAT,'ADDL  ',2,'COSCN ')
      CALL MEMMAN(KLWNO,LENC*NMAT,'ADDL  ',2,'WNO   ')
      DO IMAT = 1, NMAT-1
        WRITE(6,*) ' Comparison of occupied orbital sets ', IMAT,IMAT+1
        IB_OLD = KLCOCC + (IMAT-1)*LENC
        IB_NEW = KLCOCC + (IMAT+1-1)*LENC
        IB_COSCN = KLCOSCN + (IMAT-1)*LENC
        IB_NOLD = KLNOCCFSM + (IMAT-1)*NSMOB
        IB_NNEW = KLNOCCFSM + (IMAT+1-1)*NSMOB
        CALL OVERLAP_COCOLD_COCNEWB(WORK(IB_OLD),WORK(IB_NEW),
     &       WORK(KSAO),WORK(IB_COSCN),WORK(KLSCR),NAOS_ENV,
     &       WORK(IB_NOLD),WORK(IB_NNEW),NSMOB)
C       OVERLAP_COCOLD_COCNEWB(COLD,CNEW,S,COSCN,SCR,
C    &            NAOS,NOCOLD,NOCNEW,NSMOB)
* The weight are collected in SCR, should be put into array ..
         CALL WEIGHT_NEWOCC_OCCOLD(WORK(IB_COSCN),WORK(IB_NOLD),
     &        WORK(IB_NEW),WORK(KLSCR),WOCC_MIN)
C       WEIGHT_NEWOCC_OCCOLD(COSCN,NOCC_O,NOCC_N,W_NEW_OCCOLD,
C    &             WOCC_MIN)
       END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ANA_OI')
      RETURN
      END 
      SUBROUTINE GET_OCC_FOR_D(F,D,S,OCC,NOCC,IWAY)
*
* Obtain occupied orbitals corresponding to a given density 
* D.
* IWAY = 1 : Just diagonalize D
* IWAY = 2 : Diagonalize D to obtain occupied orbitals,
*            Obtain F in the space of occ orbitals and diagonalize
*
* Jeppe Olsen, June 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
* Input : F, S in symmetry-packed lower-half form, D in symmetry-packed 
*         complete form 
      DIMENSION F(*),D(*),S(*)
*. Output
      DIMENSION OCC(*)
*
      NTEST = 100
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT_OCD')
*. Obtain occupied and unoccupied orbitals by solving SDS C = S C
      CALL MEMMAN(KLXOCC,NOCC,'ADDL  ',2,'XOCC  ')
      LENC = LEN_AO_MAT(0)
      CALL MEMMAN(KLC,LENC,'ADDL  ',2,'LC    ')
      CALL DIAG_SDS(D,S,WORK(KLC),WORK(KLXOCC))
C       DIAG_SDS(D,S,C,XOCCNUM)
*. Extract occupied orbitals
      CALL MEMMAN(KLNOCCFSM,NSMOB,'ADDL  ',1,'NOCCFS')
      CALL MEMMAN(KLCOCC,LENC,'ADDL  ',2,'COCC  ')
      CALL GET_OCC_ORB_FROM_LIST(WORK(KLC),WORK(KLCOCC),
     &     WORK(KLXOCC),WORK(KLNOCCFSM))
*
      IF(IWAY.EQ.2) THEN 
*. Obtain in KLC Fock matrix in space of occupied orbitals
C     TRAN_SYM_BLOC_MAT3
C    (AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
        CALL MEMMAN(KLSCR ,LENC,'ADDL  ',2,'LSCR  ')
        CALL MEMMAN(KLEVEC ,LENC,'ADDL  ',2,'EVEC  ')
        CALL MEMMAN(KLEVAL ,LENC,'ADDL  ',2,'EVAL  ')
        CALL MEMMAN(KLF_E ,LENC,'ADDL  ',2,'F_E   ')
*. Expand Fock matrix to complete form 
C     CALL TRIPAK_AO_MAT(WORK(KLS_E),WORK(KSAO),2)
        CALL TRIPAK_AO_MAT(WORK(KLF_E),F,2)
        CALL TRAN_SYM_BLOC_MAT3(WORK(KLF_E),WORK(KLCOCC),NSMOB,NAOS_ENV,
     &       WORK(KLNOCCFSM),WORK(KLC),WORK(KLSCR),0)
*. Diagonalize fock-matrix in occupied space
C     DIAG_BLK_SYMMAT(A,NBLK,LBLK,X,EIGENLV,SCR,ISYM)
         CALL DIAG_BLK_SYMMAT(WORK(KLC),NSMOB,WORK(KLNOCCFSM),
     &        WORK(KLEVEC),WORK(KLEVAL),WORK(KLSCR),0)
*. And transform to the occupied canonical basis 
C     MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
         CALL MULT_BLOC_MAT(WORK(KLC),WORK(KLCOCC),WORK(KLEVEC),
     &       NSMOB,NAOS_ENV,WORK(KLNOCCFSM),NAOS_ENV,WORK(KLNOCCFSM),
     &       WORK(KLNOCCFSM),0)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Generated occupied orbitals '
        CALL APRBLM2(WORK(KLC),NAOS_ENV,WORK(KLNOCCFSM),NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GT_OCD')
      RETURN
      END
      SUBROUTINE GET_OCC_ORB_FROM_LIST(CTOT,COCC,XOCC,NOCC_FSM)
*
* Obtain list of occupied orbitals from list of occupation 
* numbers and list of orbitals. 
*
* Occupied orbitals are identified as those 
* with occ numbers close to 1 or 2 
*
*. Jeppe Olsen, June 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION CTOT(*),XOCC(*)
*. Output
      DIMENSION COCC(*),NOCC_FSM(*)
*
      IB_IN = 1
      IB_OUT = 1
      THRES = 1.0D-6
      NOCC_ACT = 0
      IORB_TOT = 0
*
      DO ISM = 1, NSMOB
        NMO = NAOS_ENV(ISM)
        NOCC_FSM(ISM) = 0
        DO IMO = 1, NMO
          IORB_TOT = IORB_TOT + 1
          IF(ABS(XOCC(IORB_TOT)-1.0D0).LT.THRES.OR.
     &       ABS(XOCC(IORB_TOT)-2.0D0).LT.THRES   ) THEN
             CALL ICOPVE(CTOT(IB_IN),COCC(IB_OUT),NMO)
             NOCC_FSM(ISM) = NOCC_FSM(ISM) + 1
             IB_OUT = IB_OUT + NMO
             NOCC_ACT = NOCC_ACT + 1
          END IF
          IB_IN = IB_IN + NMO
        END DO
      END DO
*
      NTEST = 100
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Number of occupied orbitals ', NOCC_ACT
       WRITE(6,*) ' Number of occpied orbitals per symmetry'
       CALL IWRTMA(NOCC_FSM,1,NSMOB,1,NSMOB)
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' The set of occupied orbitals '
       CALL APRBLM2(COCC,NAOS_ENV,NOCC_FSM,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE OVERLAP_COCOLD_COCNEWB(COLD,CNEW,S,COSCN,SCR,
     &          NAOS,NOCOLD,NOCNEW,NSMOB)
*
* Obtain overlap matrix between old and new mo-coefficients
* for OCCUPIED orbitals 
*
* COSCN = COLD(T) S CNEW
*
*. Jeppe Olsen, June 2005
*.              Extension to OVERLAP_COLD_CNEW : Allows 
*               different number of occ orbs in new and old.
*               Includes only occupies orbitals
*
      INCLUDE 'implicit.inc'
*. General input
      INTEGER NAOS(NSMOB)
*. Specific input
      DIMENSION COLD(*),CNEW(*),S(*)
      INTEGER NOCOLD(NSMOB),NOCNEW(NSMOB)
*. Output
      DIMENSION COSCN(*)
*. Scratch : Should be able to contain complete total symmetric matrix
      DIMENSION SCR(*)
*. Expand S matrix to complete form
C TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
      CALL TRIPAK_BLKM(COSCN,S,2,NAOS,NSMOB)
*. Obtain  COLD (T) S in SCR
C     MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
      CALL MULT_BLOC_MAT(SCR,COLD,COSCN,NSMOB,
     &     NOCOLD,NAOS,NAOS,NOCOLD,NAOS,NAOS,1)
*.Obtain COLD(T) S CNEW 
      CALL MULT_BLOC_MAT(COSCN,SCR,CNEW,NSMOB,
     &     NOCOLD,NOCNEW,NOCOLD,NAOS,NAOS,NOCNEW,0)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occ-Occ part of COLD(T) S CNEW '
        CALL APRBLM2(COSCN,NOCOLD,NOCNEW,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE WEIGHT_NEWOCC_OCCOLD(COSCN,NOCC_O,NOCC_N,W_NEW_OCCOLD,
     &           WOCC_MIN)
*
*. Obtain for each new occupied orbital weight of this in the space of 
*. of all old orbitals.
*. The information about the overlap between old and new orbitals 
*. is given in COSCN
*. The smallest weight of a new occupied orbital is returned 
* in WOCC_MIN
*
* Jeppe Olsen, June 2005
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION COSCN(*), NOCC_O(NSMOB),NOCC_N(NSMOB)
*. Output
      DIMENSION W_NEW_OCCOLD(*)
*
      IOFF = 1
      IORB_ABS = 0
      DO ISYM = 1, NSMOB
        DO I = 1, NOCC_N(ISYM)
          X = 0.0D0
          DO IOCC = 1, NOCC_O(ISYM)
            X = X + COSCN(IOFF-1+(I-1)*NOCC_O(ISYM)+IOCC)**2
          END DO
          IORB_ABS = IORB_ABS + 1
          W_NEW_OCCOLD(IORB_ABS) = X
        END DO
        IOFF = IOFF + NOCC_O(ISYM)*NOCC_N(ISYM)
      END DO
*
      WOCC_MIN = 1.0D0
      DO IOCC = 1, IORB_ABS
        WOCC_MIN = MIN(WOCC_MIN,W_NEW_OCCOLD(IOCC))
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Part of each new occ orb that is in space of old occ orbs'
        CALL WRTMAT(W_NEW_OCCOLD,NSMOB,1,IORB_ABS,1,IOB_ABS)
      END IF
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Smallest overlap between new occupied orbitals '
        WRITE(6,*) ' and space of previous occupied orbitals is ',
     &             WOCC_MIN
      END IF
*
      RETURN
      END 
      SUBROUTINE REF_GN_KAPPA(KAPPAP,KAPPAE,IAS,ISM,IWAY,IOOEX,NOOEX)
*
* Switch between packed and expanded form of an kappa matrix 
* of symmetry ISM. 
* IAS = 1 => Antisymmetric Kappa
* IAS = 2 => Symmetric Kappa
* KAPPA_E is expanded matrix in ST order
* KAPPA_P is packed matrix in TS order 
*
*
* IWAY = 1 Packed => expanded form
* IWAY = 2 Expanded to packed form 
*
*. Jeppe Olsen, July 2011, generalization of code Oct. 2004
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input and output
      INTEGER IOOEX(2,NOOEX)
      REAL*8  KAPPAP(*), KAPPAE(*)
*. local scratch
      INTEGER IOFFSM(MXPOBS)
*
*. Check input parameters
*
      ISTOP = 0
      IF(1.GT.IAS.OR.IAS.GT.2) THEN
       WRITE(6,*)  ' IAS out of range, value = ', IAS
       ISTOP = 1
      END IF
      IF(1.GT.IWAY.OR.IWAY.GT.2) THEN
       WRITE(6,*)  ' IWAY out of range, value = ', IWAY
       ISTOP = 1
      END IF
      IF(ISTOP.EQ.1) THEN
       WRITE(6,*) ' Input to REF_GN_KAPPA out of range '
       STOP       ' Input to REF_GN_KAPPA out of range '
      END IF
*
*. Set up array giving offset for block of given symmetry
*  GET_OFFSET_FOR_SYMBLK(NROW,NCOL,IRC,ISYM,NSMST,IB)
      CALL GET_OFFSET_FOR_SYMBLK(NTOOBS,NTOOBS,1,ISM,NSMOB,IOFFSM)
*. Number of elements in expanded matrix 
      LENGTH_E = IOFFSM(NSMOB) + NTOOBS(NSMOB)**2 - 1
*
      ZERO = 0.0D0
      CALL SETVEC(KAPPAE,ZERO,LENGTH_E)
*. Sign for distinguishing betweeen antisymmetric and symmetric Kappa
      IF(IAS.EQ.1) THEN
        SIGN = -1.0D0
      ELSE 
        SIGN = 1.0D0
      END IF
*
      DO IEX = 1, NOOEX
*. Orb numbers in TS- order 
        ICREA_TS = IOOEX(1,IEX)
        IANNI_TS = IOOEX(2,IEX)
*. and in ST-order
        ICREA_ST = IREOTS(ICREA_TS)
        IANNI_ST = IREOTS(IANNI_TS)
*
        ICREA_SM = ISMFSO(ICREA_ST)
        IANNI_SM = ISMFSO(IANNI_ST)
*
        ICREA_OFF = IBSO(ICREA_SM)
        IANNI_OFF = IBSO(IANNI_SM)
*
        IAC_E = IOFFSM(IANNI_SM) -1  
     &       + (ICREA_ST-ICREA_OFF)*NTOOBS(IANNI_SM)
     &       + IANNI_ST - IANNI_OFF + 1
        ICA_E = IOFFSM(ICREA_SM) - 1
     &        + (IANNI_ST-IANNI_OFF)*NTOOBS(ICREA_SM)
     &        + ICREA_ST - ICREA_OFF + 1
C?      WRITE(6,*) ' IAC_E, ICA_E = ', IAC_E, ICA_E
        IF(IWAY.EQ.1) THEN
          KAPPAE(ICA_E) = KAPPAP(IEX)
          KAPPAE(IAC_E) = SIGN*KAPPAP(IEX)
        ELSE
          KAPPAP(IEX) = KAPPAE(ICA_E)
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from  REF_GN_KAPPA'
        IF(IWAY.EQ.1) THEN 
          WRITE(6,*) ' Packed => expanded form '
        ELSE IF (IWAY.EQ.2) THEN
          WRITE(6,*) ' Expanded => packed form '
        END IF
*
        IF(IAS.EQ.1) THEN
         WRITE(6,*) ' Antisymmetric kappa '
        ELSE 
         WRITE(6,*) ' Symmetric kappa '
        END IF
*
        WRITE(6,*) ' Kappa-matrix in expanded form '
        CALL APRBLM2(KAPPAE,NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*) ' Kappa-matrix in packed form '
        CALL WRT_IOOEXCOP(KAPPAP,IOOEX,NOOEX)
      END IF
*
      RETURN
      END
      SUBROUTINE PRINT_SUMMARY_HF(E_ITER,CONVER_HF,NIT_HF)
*
* Initial summerizer of SCF convergence
*
*. Jeppe Olsen, May 2012 - when I should be preparing for the 
*  Palermo meeting
*
      INCLUDE 'implicit.inc'
      LOGICAL CONVER_HF
*
      DIMENSION E_ITER(*)
*
      WRITE(6,*) ' Summary of HF convergence: '
      WRITE(6,*) ' ==========================='
      WRITE(6,*)
      IF(CONVER_HF) THEN
        WRITE(6,'(A, I3, A )')
     &  ' Convergence was obtained in ', NIT_HF, ' iterations '
      ELSE 
        WRITE(6,'(A, I3, A )')
     &  ' Convergence was NOT obtained in ', NIT_HF, ' iterations '
      END IF
*
      WRITE(6,*) ' Iter          E_HF           Delta_E  '
      WRITE(6,*) ' ======================================'
      DO ITER = 1, NIT_HF
        IF(ITER.EQ.1) THEN
          WRITE(6,'(2X,I3,2X,F22.12)') 
     &    ITER, E_ITER(1)
        ELSE
          WRITE(6,'(2X,I3,2X,F22.12,2X,E10.3)') 
     &    ITER, E_ITER(ITER), E_ITER(ITER)-E_ITER(ITER-1)
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE LUCIA_HF_OLD(IHFSM,IHFSPC,MAXIT_HF,
     &                    E_HF,E1_FINAL,CONVER_HF)
*
* Master routine for Hartree-Fock optimization using LUCIA
*
* Written for playing around with the HF optimization, Fall of 03
*
* Jeppe Olsen
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'crun.inc'
*
      LOGICAL CONVER_HF
      CONVER_HF = .FALSE.
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HF_SEC')
*
      WRITE(6,*) ' ================================== '
      WRITE(6,*) '                                    '
      WRITE(6,*) ' Hartree-Fock optimization section  '
      WRITE(6,*) '                                    '
      WRITE(6,*) '            Version of May  2012    '
      WRITE(6,*) '                                    '
      WRITE(6,*) ' ================================== '
* 
C?    WRITE(6,*) ' I will skip to test-section for bisection '
C?    CALL TEST_BISECTION
*
*
*. Check that ci space is a single det space 
*
      IF(XISPSM(IHFSM,IHFSPC).NE.1.0D0) THEN
        WRITE(6,*) ' Hartree-Fock section entered ' 
        WRITE(6,*) ' for multi-determinant wavefunction '
        WRITE(6,*) ' Please use MCSCF keyword instead of HF '
        WRITE(6,*) ' I will stop for now '
        STOP ' Hartree-Fock section with MC wavefunction '
      END IF 
*
*. Type of wf : closed shell or high spin open shell 
*
* 
*   IREFTP = 1 => CLOSED Shell HARTRE-FOCK, 
*   IREFTP = 2 => High spin open shell single det state
*   IREFTP = 3 => Cas state or more general multireference state
      WRITE(6,*) ' IHFSPC = ', IHFSPC
      CALL CC_AC_SPACES(IHFSPC,IREFTYP)
      IF(IREFTYP.EQ.1) THEN
        WRITE(6,*) ' Closed shell single reference state '
      ELSE IF (IREFTYP.EQ.2) THEN
        WRITE(6,*) ' High spin open shell reference state '
      END IF
*. Check that the required number of iterations is less 
* than MXPHFIT
COLD  IF(MAXIT_HF.GT.MXPHFIT) THEN
COLD    WRITE(6,*) ' Too many HF iterations required '
COLD    WRITE(6,*) ' Actual and max = ', MAXIT_HF, MXPHIF
COLD    WRITE(6,*) ' Reduce number of HF iterations '
COLD    WRITE(6,*) ' or split up in several runs '
COLD    WRITE(6,*) ' or increase MXPHFIT in mxpdim.inc'
COLD    STOP '  Too many HF iterations required '
COLD  END IF
*. Number of occupied orbitals per symmetry 
      CALL GET_NOC
*. Obtain Info on Hamiltonian in AO basis 
*. S,H,D,F in AO basis 
      LEN1E = NTOOB **2
      CALL MEMMAN(KSAO,LEN1E,'ADDL  ',2,'S_AO  ')
      CALL MEMMAN(KHAO,LEN1E,'ADDL  ',2,'H_AO  ')
      CALL MEMMAN(KDAO,LEN1E,'ADDL  ',2,'D_AO  ')
      CALL MEMMAN(KFAO,LEN1E,'ADDL  ',2,'D_AO  ')
*. Current MO-AO expansion  matrix
      CALL MEMMAN(KCCUR,LEN1E,'ADDL  ',2,'CCUR  ')
*. And for saving old(er) MO-AO expansion
      CALL MEMMAN(KCOLD,LEN1E,'ADDL  ',2,'COLD  ')
*. Saving all AO densities in expanded form
      CALL MEMMAN(KDAO_COLLECT,LEN1E*MAXIT_HF,'ADDL  ',2,'DAO_CO')
*. Save all AO Fock matrices in expanded form 
      CALL MEMMAN(KFAO_COLLECT,LEN1E*MAXIT_HF,'ADDL  ',2,'FAO_CO')
*. And space for the energies in the various iterations 
      CALL MEMMAN(KEITER,(MAXIT_HF+1),'ADDL  ',2,'EITER ')

*. Integral treatment 
      I2E_AOINT_CORE = 1
      IF(I2E_AOINT_CORE .EQ.1 ) THEN
        WRITE(6,*) ' All AO-integrals are stored in core '
        WRITE(6,*) ' Number of 2e integrals = ', NINT2
        CALL MEMMAN(KINT2AO,NINT2,'ADDL  ',2,'INT2AO')
      END IF
*. Obtain AO integrals HAO, SAO
C     GETHSAO(HAO,SAO,IGET_HAO,IGET_SAO)
      CALL GET_HSAO(WORK(KHAO),WORK(KSAO),1,1)
*. And if requested 2-electron integrals in AO basis
      IF( I2E_AOINT_CORE.EQ.1) THEN
        CALL GET_H2AO(WORK(KINT2AO))
      END IF
*
*. Initial guess to MO-coefficients
*
* INI_HF_MO = 2 => Read coefficients in from fil LUMOIN
* INI_HF_MO = 1 => Diagonalize one-electron Hamiltonian
      CALL GET_INI_GUESS(WORK(KCCUR),INI_HF_MO)
*. And do the optimization 
*. Roothaan-Hall or EOPD
      IF(IHFSOLVE.EQ.1.OR.IHFSOLVE.EQ.2) THEN
        CALL OPTIM_SCF(MAXIT_HF,E_HF,CONVER_HF,E1_FINAL,NIT_HF,
     &                WORK(KEITER))
      ELSE IF(IHFSOLVE.EQ.3) THEN
C            OPTIM_SCF_USING_ONE_STEP(MAXIT_HF)       
        CALL OPTIM_SCF_USING_ONE_STEP(MAXIT_HF,
     &       E_HF,CONVER_HF,E1_FINAL,NIT_HF,WORK(KEITER))
      ELSE 
*. Second order method
        CALL OPTIM_SCF_USING_SECORDER(MAXIT_HF,
     &       E_HF,CONVER_HF,E1_FINAL,NIT_HF,WORK(KEITER))
      END IF
*. Print summary
      CALL PRINT_SUMMARY_HF(WORK(KEITER),CONVER_HF,NIT_HF)
*
      RETURN
      END 
      SUBROUTINE DIAG_SDS_GEN(D,S,C,XOCCNUM)
*
* Diagonalize Density matrix in AO basis (SDS !!),
* and sort resulting orbitals according to 
* occupation numbers, so the occupied orbitals
* occur first in each symmetry class
*
* The procedure is carried out in the generalized basis to ensure
* supersymmetry
*
*. Jeppe Olsen, May 24, 2012, from Oct. 2004
*
c      INCLUDE 'implicit.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input : D unpacked, S packed 
      DIMENSION D(*),S(*)
*. Output
      DIMENSION  C(*),XOCCNUM(*)
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIAG_S')
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from DIAG_SDS_GEN '
       WRITE(6,*) ' ======================='
       WRITE(6,*) ' First elements of S and D ', D(1), S(1)
      END IF
*. Scratch
      LENC = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,0)
      LENP = NDIM_1EL_MAT(1,NAOS_ENV,NAOS_ENV,NSMOB,1)
      CALL MEMMAN(KLS_E,LENC,'ADDL  ',2,'S_E   ')
      CALL MEMMAN(KLD_E,LENC,'ADDL  ',2,'D_E   ')
      CALL MEMMAN(KLSCR,LENC,'ADDL  ',2,'SCR   ')
      CALL MEMMAN(KLSDS,LENC,'ADDL  ',2,'SDS   ')
*. Obtain( in KLD_E) D packed in generalized symmetry
      WRITE(6,*) ' Info on reform of D '
C     REFORM_MAT_STA_SUP(ASTA,ASUP,IPACK,IWAY)
      CALL REFORM_MAT_STA_SUP(D,WORK(KLD_E),0,1)
*. Obtain in KLSCR, S in packed generalized form
      WRITE(6,*) ' Info on reform of S '
      CALL REFORM_MAT_STA_SUP(S,WORK(KLSCR),1,1)
*. Obtain (in KLS_E) S in unpacked generalized form
      CALL TRIPAK_BLKM(WORK(KLS_E),WORK(KLSCR),2,NBAS_GENSMOB,
     &     NGENSMOB)
*. DS in KLSCR
      CALL MULT_BLOC_MAT(WORK(KLSCR),WORK(KLD_E),WORK(KLS_E),
     &     NGENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,
     &     NBAS_GENSMOB,NBAS_GENSMOB,1)
*. SDS in KLSDS
      CALL MULT_BLOC_MAT(WORK(KLSDS),WORK(KLS_E),WORK(KLSCR),
     &     NGENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,
     &     NBAS_GENSMOB,NBAS_GENSMOB,1)
*. Multiply with -1 to get largest occupation numbers first
C LEN_BLMAT(NBLK,LROW,LCOL,IPACK)
      LEN_GEN = LEN_BLMAT(NGENSMOB,NBAS_GENSMOB,NBAS_GENSMOB,0)
      ONEM = -1.0D0
      CALL SCALVE(WORK(KLSDS),ONEM,LEN_GEN)
*. and diagonalize
      ISORT=1
      CALL GENDIA_BLMAT(WORK(KLSDS),WORK(KLS_E),C,XOCCNUM,
     &                  WORK(KLSCR),NBAS_GENSMOB,NGENSMOB,ISORT)
*. And multiply eigenvalues to get occupation numbers 
*. ( to compensate for the -1 we multiplied by earlier )
      CALL SCALVE(XOCCNUM,ONEM,NTOOB)
*. The C-coefficients are returned in generalized symmetry, reform to standard
C       REFORM_CMO_STA_GEN(CMO_STA,CMO_GEN,IDO_REORDER,IREO,IWAY)
        CALL REFORM_CMO_STA_GEN(WORK(KLSCR),C,0,IDUM,2)
        CALL COPVEC(WORK(KLSCR),C,LENC)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Obtained occupation numbers form DIAG_SDS... '
        CALL WRITE_ORBVECTOR(XOCCNUM,NAOS_ENV,NSMOB)
        WRITE(6,*) ' Obtained MO-coefficients '
        CALL APRBLM2(C,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DIAG_S')
*
      RETURN 
      END 
     

               

        
        


c $Id$
