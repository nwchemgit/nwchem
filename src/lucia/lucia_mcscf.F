      SUBROUTINE LUCIA_MCSCF(IREFSPC_MCSCF,MAXMAC,MAXMIC,
     &                       EFINAL,CONVER,VNFINAL)
*
* Master routine for MCSCF optimization.
*
* Jeppe Olsen
* Last revision; Oct. 10, 2012; Jeppe Olsen; Cleaning + CI in linesearch
* Last revision; Feb. 15, 2013; Jeppe Olsen; Allowing no restart of CI + selection of root
*
* Sept. 2011: Option to calculate Fock matrices from 
*             transformed integrals removed - assumed complete
*             list of transformed integrals
* Oct. 2011: reorganization of code 
*
* Initial MO-INI transformation matrix is assumed set outside and is in MOMO
* Initial MO-AO transformation matrix is in MOAOIN
*
*. Output matrix is in
*   1) MOAOUT   - as it is the output matrix
*   2) MOAO_ACT - as it is the active matrix
*   3) MOAOIN   - as the integrals are in this basis ...
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cstate.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Some indirect transfer - for communicating with EMCSCF_FROM_KAPPA
      COMMON/EXCTRNS/KLOOEXCC,KINT1_INI,KINT2_INI, IREFSPC_MCSCFL,
     &               IPRDIALL,IIUSEH0PL,MPORENP_EL,
     &               ERROR_NORM_FINALL,CONV_FL,
     &               I_DO_CI_IN_INNER_ACT
C           GASCI(IREFSM,IREFSPC_MCSCF,IPRDIAL,IIUSEH0P,
C    &           MPORENP_E,EREF,ERROR_NORM_FINAL,CONV_F)  
* A bit of local scratch
      INTEGER ISCR(2), ISCR_NTS((7+MXPR4T)*MXPOBS),only_dens
*
      REAL*8
     &INPROD
      LOGICAL DISCH, CONV_INNER
*
      LOGICAL CONV_F,CONV_FL,CONVER
      EXTERNAL EMCSCF_FROM_KAPPA
*. A bit of local scratch
C     INTEGER I2ELIST_INUSE(MXP2EIARR),IOCOBTP_INUSE(MXP2EIARR)
*
* Removing (incorrect) compiler warnings
      KINT2_FSAVE = 0
      IE2ARR_F = -1

      IDUMMY = 0
      CALL MEMMAN(IDUMMY, IDUMMY, 'MARK ', IDUMMY,'MCSCF ') 
      CALL LUCIAQENTER('MCSCF')
*
*. Local parameters defining optimization
*
*. reset kappa to zero in each inner or outer iteration
*
* IRESET_KAPPA_IN_OR_OUT = 1 => Reset kappa in each inner iteration
* IRESET_KAPPA_IN_OR_OUT = 2 => Reset kappa in each outer iteration
*
*. Use gradient or Brillouin vector (differs only when gradient is 
*  evaluated for Kappa ne. 0, ie. IRESET_KAPPA = 2
*
* I_USE_BR_OR_E1 = 1 => Use Brilloin vector
* I_USE_BR_OR_E2 = 2 => Use E1
      IRESET_KAPPA_IN_OR_OUT = 2
      I_USE_BR_OR_E1 = 2 
*. Largest allowed number of vectors in update
      NMAX_VEC_UPDATE = 50
*
*. Super symmetry options
*
      IF(I_USE_SUPSYM.EQ.1) THEN
*. Average over orbital excitations of shell excitations
        I_AVERAGE_ORBEXC = 1
*. Restrict orbital excitations in case of super-symmetry
        INCLUDE_ONLY_TOTSYM_SUPSYM = 1
      ELSE
        I_AVERAGE_ORBEXC = 0
        INCLUDE_ONLY_TOTSYM_SUPSYM = 0
      END IF
*
CNW   WRITE(6,*) ' TESTY: I_USE_SUPSYM, I_AVERAGE_ORBEXC =',
CNW  &                    I_USE_SUPSYM, I_AVERAGE_ORBEXC
CNW   WRITE(6,*) ' INCLUDE_ONLY_TOTSYM_SUPSYM = ',
CNW  &             INCLUDE_ONLY_TOTSYM_SUPSYM
*
*. Eliminate restart of CI (put a zero..)
      I_RESTART_OF_CI = 1

*
      if (ga_nodeid().eq.0) then
      WRITE(6,*) ' **************************************'
      WRITE(6,*) ' *                                    *'
      WRITE(6,*) ' * MCSCF optimization control entered *'
      WRITE(6,*) ' *                                    *'
      WRITE(6,*) ' * Version 1.3, Jeppe Olsen, Oct. 12  *'
      WRITE(6,*) ' **************************************'
      WRITE(6,*)
      WRITE(6,*) ' Occupation space: ', IREFSPC_MCSCF
      WRITE(6,*) ' Allowed number of outer iterations ', MAXMAC
      WRITE(6,*) ' Allowed number of inner iterations ', MAXMIC
*
      IF(I_USE_SUPSYM.EQ.1) THEN
        IF(INCLUDE_ONLY_TOTSYM_SUPSYM.EQ.1) THEN
          WRITE(6,*) 
     &   ' Excitations only between orbs with the same supersymmetry'
        ELSE
          WRITE(6,'(2X,A)') 
     &   'Excitations only between orbs with the same standard symmetry'
        END IF
*
        IF(I_AVERAGE_ORBEXC.EQ.1) THEN
          WRITE(6,*)
     &   'Average over orbital components of shell-excitations '
        ELSE
          WRITE(6,*)
     &   'No average over orbital components of shell-excitations '
        END IF
*
      END IF
*
      WRITE(6,*)
      WRITE(6,*) ' MCSCF optimization method in action:'
      IF(IMCSCF_MET.EQ.1) THEN
        WRITE(6,*) '    One-step method NEWTON'
      ELSE  IF (IMCSCF_MET.EQ.2) THEN
        WRITE(6,*) '    Two-step method NEWTON'
      ELSE  IF (IMCSCF_MET.EQ.3) THEN
        WRITE(6,*) '    One-step method Update'
      ELSE  IF (IMCSCF_MET.EQ.4) THEN
        WRITE(6,*) '    Two-step method Update'
      END IF
      IF(I_RESTART_OF_CI .EQ. 0) THEN
        WRITE(6,*) '    No restart of CI '
      END IF
*
      IF(IOOE2_APR.EQ.1) THEN
        WRITE(6,*) '    Orbital-Orbital Hessian constructed'
      ELSE IF(IOOE2_APR.EQ.2) THEN
        WRITE(6,*) 
     &  '    Diagonal blocks of Orbital-Orbital Hessian constructed'
      ELSE IF(IOOE2_APR.EQ.3) THEN
        WRITE(6,*) 
     &  '    Approx. diagonal of Orbital-Orbital Hessian constructed'
      END IF
*
*. Linesearch
*
      IF(IMCSCF_MET.LE.2) THEN
       IF(I_DO_LINSEA_MCSCF.EQ.1) THEN 
         WRITE(6,*) 
     &   '    Line search for Orbital optimization '
       ELSE IF(I_DO_LINSEA_MCSCF.EQ.0) THEN
         WRITE(6,*) 
     &   '    Line search when energy increases '
       ELSE IF(I_DO_LINSEA_MCSCF.EQ.-1) THEN
         WRITE(6,*) 
     &   '    Line search never carried out '
       END IF
      ELSE
*. Update method linesearch always used
        WRITE(6,*) 
     &  '    Line search for Orbital optimization '
      END IF
      IF(IMCSCF_MET.EQ.3.OR.IMCSCF_MET.EQ.4) THEN
        WRITE(6,'(A,I4)') 
     &  '     Max number of vectors in update space ', NMAX_VEC_UPDATE
      END IF
*
      IF(IRESET_KAPPA_IN_OR_OUT .EQ.1 ) THEN
        WRITE(6,*) 
     &  '       Kappa is reset to zero in each inner iteration '
      ELSE IF( IRESET_KAPPA_IN_OR_OUT .EQ.2 ) THEN
        WRITE(6,*) 
     &  '    Kappa is reset to zero in each outer iteration '
      END IF
*
      IF(I_USE_BR_OR_E1.EQ.1) THEN
        WRITE(6,*) '    Brillouin vector in use'
      ELSE IF(I_USE_BR_OR_E1 .EQ.2) THEN
        WRITE(6,*) '    Gradient vector in use'
      END IF
      endif
*
      NFRZ_ORB_ACT = 0
      IF(NFRZ_ORB.NE.0.AND.IREFSPC_MCSCF.GE.IFRZFST) THEN
        if (ga_nodeid().eq.0) then
        WRITE(6,*) ' Orbitals frozen in MCSCF optimization: '
        CALL IWRTMA3(IFRZ_ORB,1,NFRZ_ORB,1,NFRZ_ORB)
        endif
        NFRZ_ORB_ACT = NFRZ_ORB
      END IF
      
      I_MAY_DO_CI_IN_INNER_ITS = 1
      XKAPPA_THRES = 1.0D0
      MIN_OUT_IT_WITH_CI = 4
      I_MAY_DO_CI_IN_INNER_ITS = 0
      IF(IMCSCF_MET.EQ.4) THEN
        I_MAY_DO_CI_IN_INNER_ITS = 1
        XKAPPA_THRES = 1.0D0
        if (ga_nodeid().eq.0) then
        WRITE(6,'(A)') 
     &  '     CI - optimization in inner iterations starts when: '
        WRITE(6,'(A)')
     &  '       Hessian approximation is not shifted'
        WRITE(6,'(A,E8.2)') 
     &  '       Initial step is below ',  XKAPPA_THRES
        WRITE(6,'(A,I3)') 
     &  '     Outer iteration is atleast number ', MIN_OUT_IT_WITH_CI
        endif
      END IF
*
*. Initial allowed step length 
      STEP_MAX = 0.75D0
C     WRITE(6,*) ' Jeppe has reduced step to ', STEP_MAX
      TOLER = 1.1D0
      NTEST = 0000
      IPRNT= MAX(NTEST,IPRMCSCF)
*
      last_ci = 0
      I_DO_NEWTON = 0
      I_DO_UPDATE = 0
      I_UPDATE_MET = 0
      IF(IMCSCF_MET.LE.2) THEN
        I_DO_NEWTON = 1
      ELSE IF (IMCSCF_MET.EQ.3.OR.IMCSCF_MET.EQ.4) THEN
        I_DO_UPDATE = 1
*. use BFGS update
        I_UPDATE_MET = 2
*. Update vectors will be kept in core
        DISCH = .FALSE.
        LUHFIL = -2810
      END IF
COLD  WRITE(6,*) ' I_DO_NEWTON, I_DO_UPDATE = ', 
COLD &             I_DO_NEWTON, I_DO_UPDATE
      I_DO_CI_IN_INNER_ACT = 0
*
*. Memory for information on convergence of iterative procedure
      NITEM = 4
      LEN_SUMMARY = NITEM*(MAXMAC+1)
      CALL MEMMAN(KL_SUMMARY,LEN_SUMMARY,'ADDL  ',2,'SUMMRY')
*. Memory for the initial set of MO integrals
CNW   CALL MEMMAN(KINT1_INI,NINT1,'ADDL  ',2,'INT1_IN')
      CALL MEMMAN(KINT2_INI,NINT2,'ADDL  ',2,'INT2_IN')
      if(.not.ga_create(mt_dbl,NINT1,1,'INT1_IN',0,1,KINT1_INI))
     &   call errquit('cannot create lucia_mc: kint1_ini',0,GA_ERR)
*. And for two extra MO matrices 
      LEN_CMO =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLMO1,LEN_CMO,'ADDL  ',2,'MO1   ')
      CALL MEMMAN(KLMO2,LEN_CMO,'ADDL  ',2,'MO2   ')
      CALL MEMMAN(KLMO3,LEN_CMO,'ADDL  ',2,'MO3   ')
      CALL MEMMAN(KLMO4,LEN_CMO,'ADDL  ',2,'MO4   ')
*. And for storing MO coefficients from outer iteration
      CALL MEMMAN(KLMO_OUTER,LEN_CMO,'ADDL  ',2,'MOOUTE')
*. And initial set of MO's
      CALL MEMMAN(KLCMOAO_INI,LEN_CMO,'ADDL  ',2,'MOINI ')
*. Normal integrals accessed
      IH1FORM = 1
      I_RES_AB = 0
      IH2FORM = 1
*. CI not CC
      ICC_EXC = 0
* 
*. Non-redundant orbital excitations
*
*. Nonredundant type-type excitations
      CALL MEMMAN(KLTTACT,(NGAS+2)**2,'ADDL  ',1,'TTACT ')
      CALL NONRED_TT_EXC(int_mb(KLTTACT),IREFSPC_MCSCF,0)
*. Nonredundant orbital excitations
*.. Number : 
      KLOOEXC = 1
      KLOOEXCC= 1
*
      IF(I_USE_SUPSYM.EQ.1.AND.INCLUDE_ONLY_TOTSYM_SUPSYM.EQ.1) THEN
        I_RESTRICT_SUPSYM = 1
      ELSE 
        I_RESTRICT_SUPSYM = 0
      END IF
      CALL NONRED_OO_EXC2(NOOEXC,int_mb(KLOOEXC),int_mb(KLOOEXCC),
     &     1,int_mb(KLTTACT),I_RESTRICT_SUPSYM,int_mb(KMO_SUPSYM),
     &     N_INTER_EXC,N_INTRA_EXC,1)
*
      IF(NOOEXC.EQ.0) THEN
        WRITE(6,*) ' STOP: zero orbital excitations in MCSCF '
        STOP       ' STOP: zero orbital excitations in MCSCF '
      END IF
*.. And excitations
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
      CALL MEMMAN(KLOOEXCC,2*NOOEXC,'ADDL  ',1,'OOEXCC')
*. Allow these parameters to be known outside
      KIOOEXC = KLOOEXC
      KIOOEXCC = KLOOEXCC
      CALL NONRED_OO_EXC2(NOOEXC,int_mb(KLOOEXC),int_mb(KLOOEXCC),
     &     1,int_mb(KLTTACT),I_RESTRICT_SUPSYM,int_mb(KMO_SUPSYM),
     &     N_INTER_EXC,N_INTRA_EXC,2)
*. Memory for gradient 
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1_MC ')
*. And Brilluoin matrix in complete form
      CALL MEMMAN(KLBR,LEN_CMO,'ADDL  ',2,'BR_MAT')
*. And an extra gradient
      CALL MEMMAN(KLE1B,NOOEXC,'ADDL  ',2,'E1B   ')
*. and a scratch vector for excitation
      CALL MEMMAN(KLEXCSCR,NOOEXC,'ADDL  ',2,'EX_SCR')
*. Memory for gradient and orbital-Hessian - if  required
      IF(IOOE2_APR.EQ.1) THEN
        LE2 = NOOEXC*(NOOEXC+1)/2
        CALL MEMMAN(KLE2,LE2,'ADDL  ',2,'E2_MC ')
*. For eigenvectors of orbhessian
        LE2F = NOOEXC**2
        CALL MEMMAN(KLE2F,LE2F,'ADDL  ',2,'E2_MCF')
*. and eigenvalues, scratch, kappa
        CALL MEMMAN(KLE2VL,NOOEXC,'ADDL  ',2,'EIGVAL')
      ELSE
        KLE2 = -1
        KLE2F = -1
        KLE2VL = -1
      END IF
      KLIBENV = -2810
      KCLKSCR = -2810
* Shell-excitations, if supersymmetry is in use 
      IF(I_USE_SUPSYM.EQ.1) THEN
*. Obtain nonredundant shell excitations - output is in pointers defined in NONRED_SS
C            NONRED_SS_EXC(NOOEX,IOOEXC,NSSEX)
        CALL NONRED_SS_EXC(NOOEXC,int_mb(KIOOEXCC), NSSEX)
      END IF
*
      
      IF(I_DO_UPDATE.EQ.1) THEN
*. Space for update procedure
*. Array defining envelope and a scratch vector
* and matrix
        CALL MEMMAN(KLIBENV,NOOEXC,'ADDL  ',2,'IBENV')
        CALL MEMMAN(KLCLKSCR,NOOEXC,'ADDL  ',2,'CLKSCR')
*. rank 2 matrices
        CALL MEMMAN(KLRANK2,4*NMAX_VEC_UPDATE,'ADDL  ',2,'RNK2MT')
* two vectors defining each rank two-space
        LENGTH_V = 2*NMAX_VEC_UPDATE*NOOEXC
        CALL MEMMAN(KLUPDVEC,LENGTH_V,'ADDL  ',2,'RNK2VC')
*. Vectors for saving previous kappa and gradient
        CALL MEMMAN(KLE1PREV,NOOEXC,'ADDL  ',2,'E1PREV')
        CALL MEMMAN(KLKPPREV,NOOEXC,'ADDL  ',2,'KPPREV')
C KLRANK2, KLUPDVEC, KLCLKSCR,KLE1PREV,KLKPPREV
      END IF
*. 
*. and scratch, kappa
      CALL MEMMAN(KLE2SC,NOOEXC,'ADDL  ',2,'EIGSCR')
      CALL MEMMAN(KLKAPPA,NOOEXC,'ADDL  ',2,'KAPPA ')
      CALL MEMMAN(KLSTEP, NOOEXC,'ADDL  ',2,'STEP  ')
*. Save the initial set of MO integrals 
CNW   CALL COPVEC(WORK(KINT1O),WORK(KINT1_INI),NINT1)
      CALL COPVEC(dbl_mb(KINT_2EMO),dbl_mb(KINT2_INI),NINT2)
      call ga_copy_patch('N',KINT1O,1,nint1,1,1,
     &                   KINT1_INI,1,nint1,1,1)
*. Print will be reduced for densities
      IPRDEN_SAVE = IPRDEN
      IPRDEN = 0
      IRESTR_SAVE = IRESTR
*
      IIUSEH0P = 0
      MPORENP_E = 0
      IPRDIAL = IPRMCSCF
      IPRDIAL = 0
*. Transfer to common block for communicating with EMCSCF
C     COMMON/EXCTRNS/KLOOEXCC,KINT1_INI,KINT2_INI, IREFSPC_MCSCFL,
C    &               IPRDIALL,IIUSEH0PL,MPORENP_EL,
C    &               ERROR_NORM_FINALL,CONV_FL,
C    &               I_DO_CI_IN_INNER_ACT
      IREFSPC_MCSCFL = IREFSPC_MCSCF
      IPRDIALL = IPRDIAL
      IIUSEH0PL = IIUSEH0P
      MPORENP_EL = MPORENP_E

*
      CONVER = .FALSE.
      CONV_F = .FALSE.
*. The various types of integral lists- should probably be made in
* def of lists
      IE2LIST_0F = 1
      IE2LIST_1F = 2
      IE2LIST_2F = 3
      IE2LIST_4F = 5
*. For integral transformation: location of MO coefs
      KKCMO_I = KMOMO
      KKCMO_J = KMOMO
      KKCMO_K = KMOMO
      KKCMO_L = KMOMO
*
      IF(I_DO_UPDATE.EQ.1) THEN
*. Define envelope for used orbital Hessian - pt complete
* is constructed so
        IONE = 1
        CALL ISETVC(dbl_mb(KLIBENV),IONE,NOOEXC)
      END IF
*
*. Loop over outer iterations
*
* In summery
* 1: Norm of orbgradient
* 2: Norm of orbstep
* 3: Norm of CI after iterative procedure
* 4: Energy
*
*. Convergence is pt  energy change le THRES_E
*
      ZERO = 0.0D0
      NMAT_UPD = 0
*. Line search is not meaning full very close to convergence
      THRES_FOR_ENTER_LINSEA = 1.0D-8

      N_INNER_TOT = 0
      DO IOUT = 1, MAXMAC
*
        call ga_sync()
        IF(IPRNT.GE.1.and.ga_nodeid().eq.0) THEN
          WRITE(6,*)
          WRITE(6,*) ' ----------------------------------'
          WRITE(6,*) ' Output from outer iteration', IOUT
          WRITE(6,*) ' ----------------------------------'
          WRITE(6,*)
        END IF
        NOUTIT = IOUT
*
* Save initial integrals
*
CNW     CALL COPVEC(WORK(KINT1_INI),WORK(KINT1O),NINT1)
*
* Bert: Not for first iteration, as these should be the same
*
        if (IOUT.GT.1) then
        CALL COPVEC(dbl_mb(KINT2_INI),dbl_mb(KINT_2EMO),NINT2)
        call ga_copy_patch('N',KINT1_INI,1,nint1,1,1,
     &                   KINT1O,1,nint1,1,1)
        endif
*
*. Transform integrals to current set of MO's
*
        IF(IPRNT.GE.10) WRITE(6,*) ' Integral transformation:' 
*. Where integrals should be read from
        KINT2 = KINT2_INI
*. Flag type of integral list to be obtained
*. Flag for integrals with two  free index: energy + gradient+orb-Hessian
*. Check problem: raise!!
        IE2LIST_AL = IE2LIST_2F
        IE2LIST_AL = IE2LIST_4F
        IOCOBTP_AL = 1
        INTSM_AL = 1
*. Perform integral transformation and construct inactive Fock matrix
C       DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
        CALL DO_ORBTRA(1,1,0,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
        call ga_sync()
CNW     CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
        call ga_copy(KFI,KINT1)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' MCSCF: ECORE_ORIG, ECORE_HEX, ECORE(2) ',
     &                 ECORE_ORIG, ECORE_HEX, ECORE
        END IF
*. Prepare for reading integrals
        IE2ARRAY_A = IE2LIST_I(IE2LIST_IB(IE2LIST_A))
*
*. Perform CI - and calculate densities
*
        IF(IPRNT.GE.10) WRITE(6,*) ' CI: '
        IF(IOUT.NE.1) IRESTR = 1
        MAXIT_SAVE = MAXIT
C       MAXIT = MAXMIC
        IF(I_RESTART_OF_CI.EQ.0) THEN
*. This is CI without restart. Should special settings be used?
*. Root selection is not used in  iterative procedure
          IRESTR = 0
          IROOT_SEL_SAVE = IROOT_SEL
          IROOT_SEL = 0
          WRITE(6,*) ' INI_SROOT, INI_NROOT = ',
     &                 INI_SROOT, INI_NROOT
          IF(INI_SROOT.NE.INI_NROOT) THEN
            NROOT_SAVE = NROOT
            MXCIV_SAVE = MXCIV
            NROOT = INI_NROOT
            MXCIV = MAX(2*NROOT,MXCIV_SAVE)
            WRITE(6,*) ' INI_*ROOT option in action '
          END IF ! special setting for initial CI
        END IF
        IPRDIA_SAVE = IPRDIA
C       IPRDIA = 1
        call ga_sync()
        CALL GASCI(IREFSM,IREFSPC_MCSCF,IPRDIAL,IIUSEH0P,
     &             MPORENP_E,EREF,ERROR_NORM_FINAL,CONV_F)  
        IROOT_SEL = IROOT_SEL_SAVE
        IF(I_RESTART_OF_CI.EQ.0.AND.INI_SROOT.NE.NROOT) THEN
*. Reset parameters
              NROOT = NROOT_SAVE
              MXCIV = MXCIV_SAVE
              IROOT_SEL = IROOT_SEL_SAVE
        END IF
C       WRITE(6,*) ' TESTY, NROOT(b) = ', NROOT
*
        MAXIT = MAXIT_SAVE
        if (ga_nodeid().eq.0) 
     &  WRITE(6,*) ' Energy and residual from CI :', 
     &  EREF,ERROR_NORM_FINAL
        dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+3) = ERROR_NORM_FINAL
        dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+4) = EREF
        EOLD = EREF
        ENEW = EREF
*. (Sic)
*
        IF(IOUT.GT.1) THEN
*. Check for convergence
          DELTA_E = dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+4)-
     &              dbl_mb(KL_SUMMARY-1+(IOUT-1-1)*NITEM+4)
          IF(IPRNT.GE.2.and.ga_nodeid().eq.0) WRITE(6,'(A,E9.3)') 
     &    '  Change of energy between outer iterations = ', DELTA_E
          IF(ABS(DELTA_E).LE.THRES_E) CONVER = .TRUE.
        END IF
        IF(CONVER) THEN
          NOUTIT = NOUTIT-1
          IF(IPRNT.GE.1.and.ga_nodeid().eq.0) THEN
            WRITE(6,*) ' MCSCF calculation has converged'
          END IF
          GOTO 1001
        END IF
*. A test
C       CALL EN_FROM_DENS(ENERGY,2,0)
        CALL EN_FROM_DENS(ENERGY2,2,0)
        if (ga_nodeid().eq.0)
     &  WRITE(6,*) ' Energy from density matrices ', ENERGY2
*. The active Fock matrix
C       DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
        CALL DO_ORBTRA(0,0,1,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
*
*.======================================
*. Exact or approximate orbital Hessian 
*.======================================
*
*
*. Fock matrix in KF
        CALL FOCK_MAT_STANDARD(KF,2,KFI,KFA)
        IOOSM = 1
C            ORBHES(OOHES,IOOEXC,NOOEXC,IOOSM,ITTACT)
        IF(IOOE2_APR.EQ.1) THEN
          CALL ORBHES(dbl_mb(KLE2),int_mb(KLOOEXC),NOOEXC,IOOSM,
     &         int_mb(KLTTACT))
          IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
           WRITE(6,*) ' The orbital Hessian '
           CALL PRSYM(dbl_mb(KLE2),NOOEXC)
          END IF
        END IF
*
*. Diagonalize to determine lowest eigenvalue
*
*. Outpack to complete form
        CALL TRIPAK(dbl_mb(KLE2F),dbl_mb(KLE2),2,NOOEXC,NOOEXC)
C            TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
*. Lowest eigenvalue
C            DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
        CALL DIAG_SYMMAT_EISPACK(dbl_mb(KLE2F),dbl_mb(KLE2VL),
     &       dbl_mb(KLE2SC),NOOEXC,IRETURN)
        IF(IRETURN.NE.0) THEN
           WRITE(6,*) 
     &     ' Problem with diagonalizing E2, IRETURN =  ', IRETURN
        END IF
        IF(IPRNT.GE.1000.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) ' Eigenvalues: '
          CALL WRTMAT(dbl_mb(KLE2VL),1,NOOEXC,1,NOOEXC)
        END IF
*. Lowest eigenvalue
        E2VL_MN = XMNMX(dbl_mb(KLE2VL),NOOEXC,1)
        IF(IPRNT.GE.2.and.ga_nodeid().eq.0)  
     &  WRITE(6,*) ' Lowest eigenvalue of E2(orb) = ', E2VL_MN
*
*. Cholesky factorization orbital Hessian if required
*
        I_SHIFT_E2 = 0
        IF(I_DO_UPDATE.EQ.1) THEN
*. Cholesky factorization requires positive matrices.
*. add a constant to diagonal if needed
          XMINDIAG = 1.0D-4
          IF(E2VL_MN.LE.XMINDIAG) THEN
           ADD = XMINDIAG - E2VL_MN 
C               ADDDIA(A,FACTOR,NDIM,IPACK)
           CALL ADDDIA(dbl_mb(KLE2),ADD,NOOEXC,1)
           I_SHIFT_E2 = 1
          END IF
C CLSKHE(AL,X,B,NDIM,IB,IALOFF,ITASK,INDEF)
C         WRITE(6,*) ' NOOEXC before CLSKHE = ', NOOEXC 
          CALL CLSKHE(dbl_mb(KLE2),XDUM,XDUM,NOOEXC,dbl_mb(KLIBENV),
     &         dbl_mb(KLCLKSCR),1,INDEF)
          IF(INDEF.NE.0) THEN
            WRITE(6,*) ' Indefinite matrix in CKSLHE '
            STOP ' Indefinite matrix in CKSLHE '
          END IF
        END IF! Cholesky decomposition required
*
*
*. Finite difference check
*
        I_DO_FDCHECK = 0
        IF(I_DO_FDCHECK.EQ.1) THEN
*. First: Analytic gradient from Fock matrix - As kappa = 0, Brillouin vector
* = gradient
          CALL E1_FROM_F(dbl_mb(KLE1),KF,1,int_mb(KLOOEXC),
     &                   int_mb(KLOOEXCC),
     &                   NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
*
          CALL MEMMAN(KLE1FD,NOOEXC,'ADDL  ',2,'E1_FD ')
          LE2 = NOOEXC*NOOEXC
          CALL MEMMAN(KLE2FD,LE2,   'ADDL  ',2,'E2_FD ')
          CALL SETVEC(dbl_mb(KLE2VL),ZERO,NOOEXC)
          CALL GENERIC_GRA_HES_FD(E0,dbl_mb(KLE1FD),dbl_mb(KLE2FD),
     &         dbl_mb(KLE2VL),NOOEXC,EMCSCF_FROM_KAPPA)
C              GENERIC_GRA_HES_FD(E0,E1,E2,X,NX,EFUNC)
*. Compare gradients
          ZERO = 0.0D0
          CALL CMP2VC(dbl_mb(KLE1FD),dbl_mb(KLE1),NOOEXC,ZERO)
*. transform Finite difference Hessian to packed form
          CALL TRIPAK(dbl_mb(KLE2FD),dbl_mb(KLE2F),1,NOOEXC,NOOEXC)
          LEN = NOOEXC*(NOOEXC+1)/2
          CALL CMP2VC(dbl_mb(KLE2),dbl_mb(KLE2F),LEN,ZERO)
              STOP ' Enforced stop after FD check'
        END IF
*       ^ End of finite difference check
*. Initialize sum of steps for outer iteration
        dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+2) = 0.0D0
*. Loop over Inner iterations, where orbitals are optimized
*. Initialize Kappa as zero
        IF(IRESET_KAPPA_IN_OR_OUT.EQ.2) THEN
          CALL SETVEC(dbl_mb(KLKAPPA),ZERO,NOOEXC)
        END IF
*. Save MO's from start of each outer iteration
        CALL COPVEC(dbl_mb(KMOMO),dbl_mb(KMOREF),LEN_CMO)
*. Convergence Threshold for inner iterations
*. At the moment just chosen as the total convergence threshold
        THRES_E_INNER = THRES_E
        CONV_INNER = .FALSE.
        I_DID_CI_IN_INNER = 0
*
        DO IINNER = 1, MAXMIC
          N_INNER_TOT = N_INNER_TOT + 1
*
          IF(IPRNT.GE.5.and.ga_nodeid().eq.0) THEN
            WRITE(6,*)
            WRITE(6,*) ' Info from inner iteration = ', IINNER
            WRITE(6,*) ' ===================================='
            WRITE(6,*)
          END IF
*
          IF(IRESET_KAPPA_IN_OR_OUT.EQ.1) THEN
            CALL SETVEC(dbl_mb(KLKAPPA),ZERO,NOOEXC)
          END IF
          E_INNER_OLD = EREF
          EOLD = ENEW
*
          IF(IINNER.NE.1) THEN
*
*. gradient integral transformation and Fock matrices
*
*. Flag type of integral list to be obtained:
*. Flag for integrals with one free index: energy + gradient
           IE2LIST_AL = IE2LIST_1F
           IE2LIST_AL = IE2LIST_4F
           IOCOBTP_AL = 1
           INTSM_AL = 1
*. Integral transformation, inactive and active Fock-matrices
           CALL DO_ORBTRA(1,1,1,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
CNW        CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
CNW        CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
           call ga_sync()
           call ga_copy(KINT1,KH)
           call ga_copy(KFI,KINT1)
*
           IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' ECORE_ORIG, ECORE_HEX, ECORE(2) ',
     &                    ECORE_ORIG, ECORE_HEX, ECORE
           END IF
*. Fock matrix in KF
          CALL FOCK_MAT_STANDARD(KF,2,KFI,KFA)
          END IF ! IINNER .ne.1
*
*. Construct orbital gradient
*
          IF(IPRNT.GE.10.and.ga_nodeid().eq.0) 
     &    WRITE(6,*) ' Construction of E1: '
          XKAPPA_NORM = SQRT(INPROD(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),
     &                              NOOEXC))
          IF(I_USE_BR_OR_E1.EQ.1.OR.XKAPPA_NORM.EQ.0.0D0) THEN
*. Brillouin vector from Fock matrix is used
           CALL E1_FROM_F(dbl_mb(KLE1),KF,1,int_mb(KLOOEXC),
     &                   int_mb(KLOOEXCC),
     &                   NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
          ELSE
*. Calculate gradient at non-vanishing Kappa
*. Complete Brillouin matrix
C              GET_BRT_FROM_F(BRT,F)
          CALL GET_BRT_FROM_F(dbl_mb(KLBR),KF)
C              E1_MCSCF_FOR_GENERAL_KAPPA(E1,F,KAPPA)
          CALL E1_MCSCF_FOR_GENERAL_KAPPA(dbl_mb(KLE1),dbl_mb(KLBR),
     &         dbl_mb(KLKAPPA))
          END IF
          IF(I_AVERAGE_ORBEXC.EQ.1) THEN
*. Average over orbital excitations belonging to a given shell excitation
C                SHELL_AVERAGE_ORBEXC(VECIN,NSSEX,NOOFSS,IBOOFSS,
C    &                                IOOFSS,VECUT,NOOEX,ICOPY)
            CALL SHELL_AVERAGE_ORBEXC(dbl_mb(KLE1),NSSEX,
     &           int_mb(KNOOFSS),int_mb(KIBOOFSS),int_mb(KIOOFSS),
     &           dbl_mb(KLEXCSCR),NOOEXC,1) 
          END IF


          IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
            WRITE(6,*) ' E1, Gradient: '
            CALL WRTMAT(dbl_mb(KLE1),1,NOOEXC,1,NOOEXC)
          END IF
*
          E1NRM = SQRT(INPROD(dbl_mb(KLE1),dbl_mb(KLE1),NOOEXC))
          IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &    WRITE(6,*) ' Norm of orbital gradient ', E1NRM
          dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+1) = E1NRM
*
* ==========================
* Two step Newton procedure
* ==========================
*
          IF(I_DO_NEWTON.EQ.1) THEN
*
*. Transform gradient to diagonal basis
*
*. (save original gradient)
            CALL COPVEC(dbl_mb(KLE1),dbl_mb(KLE1B),NOOEXC)
            CALL MATVCC(dbl_mb(KLE2F),dbl_mb(KLE1),dbl_mb(KLE2SC),
     &           NOOEXC,NOOEXC,1)
            CALL COPVEC(dbl_mb(KLE2SC),dbl_mb(KLE1),NOOEXC)
*
*. Solve shifted NR equations with step control
*
*           SOLVE_SHFT_NR_IN_DIAG_BASIS(
*    &            E1,E2,NDIM,STEP_MAX,TOLERANCE,X,ALPHA)A
            CALL SOLVE_SHFT_NR_IN_DIAG_BASIS(dbl_mb(KLE1),
     &           dbl_mb(KLE2VL),NOOEXC,STEP_MAX,TOLER,
     &           dbl_mb(KLSTEP),ALPHA,DELTA_E_PRED)
*
*
            XNORM_STEP = SQRT(INPROD(dbl_mb(KLSTEP),dbl_mb(KLSTEP),
     &                               NOOEXC))
*. Is step close to max
            I_CLOSE_TO_MAX = 0 
            IF(0.8D0.LE.XNORM_STEP/STEP_MAX) I_CLOSE_TO_MAX  = 1
*
            dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+2) = 
     &      dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+2) + XNORM_STEP
            IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &      WRITE(6,'(A,2(2X,E12.5))')
     &      ' Norm of step and predicted energy change = ',
     &       XNORM_STEP, DELTA_E_PRED
*. transform step to original basis
            CALL MATVCC(dbl_mb(KLE2F),dbl_mb(KLSTEP),dbl_mb(KLE2SC),
     &           NOOEXC,NOOEXC,0)
            CALL COPVEC(dbl_mb(KLE2SC),dbl_mb(KLSTEP),NOOEXC)
*
            IF(I_AVERAGE_ORBEXC.EQ.1) THEN
*. Average kappa- a bit unclean, as energy prediction is not strictly valid anymore
CM           WRITE(6,*) ' MEMCHECK before second SHELL_AVERAGE'
CM           CALL MEMCHK2('BEAVE2')
             CALL SHELL_AVERAGE_ORBEXC(dbl_mb(KLSTEP),NSSEX,
     &            int_mb(KNOOFSS),int_mb(KIBOOFSS),int_mb(KIOOFSS), 
     &            dbl_mb(KLEXCSCR),NOOEXC,1) 
CM           CALL MEMCHK2('AFAVE2')
            END IF
            IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
              WRITE(6,*) ' Step in original basis:'
              CALL WRTMAT(dbl_mb(KLSTEP),1,NOOEXC,1,NOOEXC)
            END IF
*. Is direction down-hills
            E1STEP = INPROD(dbl_mb(KLSTEP),dbl_mb(KLE1B),NOOEXC)
            IF(IPRNT.GE.2) WRITE(6,'(A,E12.5)')
     &      ' < E1!Step> = ', E1STEP
            IF(E1STEP.GT.0.0D0.AND.IPRNT.GE.2.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' Warning: step is in uphill direction '
            END IF
*. Energy for rotated orbitals
*
            ONE = 1.0D0
            CALL VECSUM(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),dbl_mb(KLSTEP),
     &      ONE,ONE,NOOEXC)
            XNORM2 = SQRT(INPROD(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),
     &                           NOOEXC))
            if(ga_nodeid().eq.0)
     &      WRITE(6,*) ' Norm of updated kappa step =', XNORM2
            call ga_sync()
            ENERGY1 = EMCSCF_FROM_KAPPA(dbl_mb(KLKAPPA))
            ENEW = ENERGY1
            if(ga_nodeid().eq.0)
     &      WRITE(6,*) ' Energy for rotated orbitals', ENERGY1
*. Compare old and new energy to decide with to do
            DELTA_E_ACT = ENEW-EOLD
            E_RATIO = DELTA_E_ACT/DELTA_E_PRED  
            IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &      WRITE(6,'(A,3(2X,E12.5))') 
     &      ' Actual and predicted energy change, ratio ', 
     &      DELTA_E_ACT, DELTA_E_PRED,E_RATIO
*
            IF(E_RATIO.LT.0.0D0) THEN
             IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &          WRITE(6,*) ' Trustradius reduced '
             RED_FACTOR = 2.0D0
             STEP_MAX = STEP_MAX/RED_FACTOR
             IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &          WRITE(6,*) ' New trust-radius ', STEP_MAX
            END IF
            IF(IOUT.GT.1.AND.E_RATIO.GT.0.8D0.AND.I_CLOSE_TO_MAX.EQ.1) 
     &      THEN
             IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &          WRITE(6,*) ' Trustradius increased '
             XINC_FACTOR = 1.5D0
             STEP_MAX = STEP_MAX*XINC_FACTOR
             IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &         WRITE(6,*) ' New trust-radius ', STEP_MAX
            END IF
C?          WRITE(6,*) ' ABS(E_RATIO-1.0D0) = ', ABS(E_RATIO-1.0D0)
*
            IF((ABS(DELTA_E_ACT).GT.THRES_FOR_ENTER_LINSEA.AND.
     &         ABS(E_RATIO-1.0D0).GT.0.1D0).AND.
     &         (I_DO_LINSEA_MCSCF.EQ.1.OR.
     &         I_DO_LINSEA_MCSCF.EQ.2.AND.EOLD.GT.ENEW)) THEN
*
*. line-search for orbital optimization
*
C                 LINES_SEARCH_BY_BISECTION(FUNC,REF,DIR,NVAR,XINI,
C    &            XFINAL,FFINAL,IKNOW,F0,FXINI)
*. Step was added to Kappa when calculating energy, get Kappa back
              ONE = 1.0D0
              ONEM = -1.0D0
              CALL VECSUM(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),
     &        dbl_mb(KLSTEP),
     &        ONE,ONEM,NOOEXC)
              CALL LINES_SEARCH_BY_BISECTION(EMCSCF_FROM_KAPPA,
     &             dbl_mb(KLKAPPA),dbl_mb(KLSTEP),NOOEXC,ONE,
     &             XFINAL,FFINAL,2, EOLD, ENEW)
              ENEW = FFINAL
              IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &          WRITE(6,*) ' Line search value of X = ',
     &        XFINAL
              XKAPPA_NORM2 = 
     &        SQRT(INPROD(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),NOOEXC))
              CALL VECSUM(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),
     &             dbl_mb(KLSTEP),
     &             ONE, XFINAL,NOOEXC)
            END IF! line search should be entered
            IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
              WRITE(6,*) ' Updated total Kappa '
              CALL WRTMAT(dbl_mb(KLKAPPA),1,NOOEXC,1,NOOEXC)
            END IF
          END IF! Newton method
          call ga_sync()
          CALL MEMCHK2('AF_NEW')
          IF(I_DO_UPDATE.EQ.1) THEN
*
* ====================
*  Update procedure
* ====================
*
*. Update Hessian
            IF(IINNER.EQ.1) THEN
*. Just save current info
              CALL COPVEC(dbl_mb(KLE1),dbl_mb(KLE1PREV),NOOEXC)
              CALL COPVEC(dbl_mb(KLKAPPA),dbl_mb(KLKPPREV),NOOEXC)
              NMAT_UPD = 0
            ELSE
C             HESUPV (E2,AMAT,AVEC,
C    &                 X,E1,VEC2,
C    &                 VEC3,NVAR,IUPDAT,IINV,VEC1,NMAT,
C    &                 LUHFIL,DISCH,IHSAPR,IBARR,E2,VEC4)
C            HESUPV (HDIAG,A,AVEC,X,G,XPREV,GPREV,NVAR,
C    &                   IUPDAT,IINV,SCR,NMAT,LUHFIL,DISCH,
C    &                   IHSAPR,IB,E2,VEC4)

*. Update on inverse
              IINV = 1
*. Initial approximation is a cholesky factorized matrix
              IHSAPR = 3
              CALL HESUPV(dbl_mb(KLE2),dbl_mb(KLRANK2),dbl_mb(KLUPDVEC),
     &             dbl_mb(KLKAPPA),dbl_mb(KLE1),dbl_mb(KLKPPREV),
     &             dbl_mb(KLE1PREV),NOOEXC,I_UPDATE_MET,IINV,
     &             dbl_mb(KLCLKSCR),NMAT_UPD,LUHFIL,DISCH,IHSAPR,
     &             dbl_mb(KLIBENV),dbl_mb(KLE2),dbl_mb(KLEXCSCR))
*. Forget the first(when starting out with exact Hessian)
              NMAT_UPD = NMAT_UPD + 1
COLD          IF(IOUT.GE.2) THEN
COLD            NMAT_UPD = 0
COLD            WRITE(6,*) ' Fusk:  update removed '
COLD            WRITE(6,*) ' Fusk:  update removed '
COLD            WRITE(6,*) ' Fusk:  update removed '
COLD            WRITE(6,*) ' Fusk:  update removed '
COLD            WRITE(6,*) ' Fusk:  update removed '
COLD          END IF
            END IF! IINNER = 1
*
*. New search direction = step
*==============================
*
*. Inverse initial Hessian approximation times gradient
            IF(IHSAPR.EQ.1) THEN
*. Just inverse diagonal (in E2) times gradient
              CALL VVTOV(dbl_mb(KLE2),dbl_mb(KLE1),dbl_mb(KLSTEP),
     &                   NOOEXC)
            ELSE
              CALL COPVEC(dbl_mb(KLE1),dbl_mb(KLCLKSCR),NOOEXC)
C                  CLSKHE(AL,X,B,NDIM,IB,IALOFF,ITASK,INDEF)
              CALL CLSKHE(dbl_mb(KLE2),dbl_mb(KLSTEP),dbl_mb(KLCLKSCR),
     &             NOOEXC,dbl_mb(KLIBENV),dbl_mb(KLEXCSCR),2,INDEF)
            END IF
            IF(NTEST.GE.10000.and.ga_nodeid().eq.0) THEN
              WRITE(6,*) ' Contribution from H(ini) to (-1) step:'
              CALL WRTMAT(dbl_mb(KLSTEP),1,NOOEXC,1,NOOEXC)
            END IF
*. And the rank-two updates
            IF(NMAT_UPD.NE.0) THEN
C                SLRMTV(NMAT,NVAR,A,AVEC,NRANK,VECIN,VECOUT,IZERO,
C    &                  DISCH,LUHFIL)
              IZERO = 0
              CALL SLRMTV(NMAT_UPD,NOOEXC,dbl_mb(KLRANK2),
     &                    dbl_mb(KLUPDVEC),2,dbl_mb(KLE1),
     &                    dbl_mb(KLSTEP),IZERO,DISCH,LUHFIL)
            END IF
*. And the proverbial minus 1
            ONEM = -1.0D0
            CALL SCALVE(dbl_mb(KLSTEP),ONEM,NOOEXC)
*. Check norm and reduce to STEP_MAX if required
            STEP_NORM = SQRT(INPROD(dbl_mb(KLSTEP),dbl_mb(KLSTEP),
     &                              NOOEXC))
            IISCALE = 0
            IF(STEP_NORM.GT.STEP_MAX) THEN
              FACTOR = STEP_MAX/STEP_NORM
              IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &        WRITE(6,'(A,E8.2)') ' Step reduced by factor = ', FACTOR
              CALL SCALVE(dbl_mb(KLSTEP),FACTOR,NOOEXC)
              IISCALE = 1
            END IF
*
            IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
              WRITE(6,*) ' Step:'
              CALL WRTMAT(dbl_mb(KLSTEP),1,NOOEXC,1,NOOEXC)
            END IF
*. Is direction down-hills
            E1STEP = INPROD(dbl_mb(KLSTEP),dbl_mb(KLE1),NOOEXC)
            IIREVERSE = 0
            IF(IPRNT.GE.2.and.ga_nodeid().eq.0) WRITE(6,'(A,E12.5)')
     &      '  < E1!Step> = ', E1STEP
            IF(E1STEP.GT.0.0D0) THEN
             IF(IPRNT.GE.2.and.ga_nodeid().eq.0) THEN
             WRITE(6,*) ' Warning: step is in uphill direction '
             WRITE(6,*) ' Sign of step is changed '
             ENDIF
             ONEM = -1.0D0
             CALL SCALVE(dbl_mb(KLSTEP),ONEM,NOOEXC)
             IIREVERSE = 1
            END IF
            XNORM_STEP = SQRT(INPROD(dbl_mb(KLSTEP),dbl_mb(KLSTEP),
     &                               NOOEXC))
            IF(IPRNT.GE.2.and.ga_nodeid().eq.0) WRITE(6,'(A,E12.5)')
     &      '  Norm of step  = ', XNORM_STEP
*. Are conditions for performing CI in inner its satisfied
            I_DO_CI_IN_INNER_ACT = 0
            IF(I_MAY_DO_CI_IN_INNER_ITS.EQ.1.AND.I_SHIFT_E2.EQ.0.AND.
     &         IOUT.GE.MIN_OUT_IT_WITH_CI.AND.
     &         XNORM_STEP.LT.XKAPPA_THRES.AND.IISCALE.EQ.0.AND.
     &         IIREVERSE.EQ.0) THEN
                 I_DO_CI_IN_INNER_ACT = 1
            END IF
*. Perform CI if required for this step
            IF(I_DO_CI_IN_INNER_ACT.EQ.1) THEN
            IF(IPRNT.GE.10.and.ga_nodeid().eq.0) 
     &         WRITE(6,*) ' CI in inner it '
            I_DID_CI_IN_INNER = 1
            I_DO_CI_IN_INNER_ITS = 1
            if(ga_nodeid().eq.0) WRITE(6,*) ' CI in inner it '
*. Integrals are in place
*
*. Perform CI - and calculate densities
*
            IF(IPRNT.GE.10.and.ga_nodeid().eq.0) WRITE(6,*) ' CI: '
            IRESTR = 1
            MAXIT_SAVE = MAXIT
            MAXIT = 5
C           WRITE(6,*) ' Number of CI-iterations reduced to 1 '
CNW         CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
            call ga_copy(KFI,KINT1)
*
            IF(I_RESTART_OF_CI.EQ.0) THEN
*. This is CI without restart. Should special settings be used?
*. Root selection is not used in  iterative procedure
              IRESTR = 0
              IROOT_SEL_SAVE = IROOT_SEL
              IROOT_SEL = 0
              IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &                     WRITE(6,*) ' INI_SROOT, INI_NROOT = ',
     &                     INI_SROOT, INI_NROOT
              IF(INI_SROOT.NE.INI_NROOT) THEN
                NROOT_SAVE = NROOT
                MXCIV_SAVE = MXCIV
                NROOT = INI_NROOT
                MXCIV = MAX(2*NROOT,MXCIV_SAVE)
                IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &            WRITE(6,*) ' INI_*ROOT option in action '
              END IF ! special setting for initial CI
            END IF
*
            CALL GASCI(IREFSM,IREFSPC_MCSCF,IPRDIAL,IIUSEH0P,
     &           MPORENP_E,EREF,ERROR_NORM_FINAL,CONV_F)  
            IROOT_SEL = IROOT_SEL_SAVE
            IF(IRESTR.EQ.0.AND.INI_SROOT.NE.NROOT) THEN
*. Reset parameters
              NROOT = NROOT_SAVE
              MXCIV = MXCIV_SAVE
              IROOT_SEL = IROOT_SEL_SAVE
            END IF
*
            MAXIT = MAXIT_SAVE
            if(ga_nodeid().eq.0)
     &      WRITE(6,*) ' Energy and residual from CI :', 
     &      EREF,ERROR_NORM_FINAL
            ENEW  = EREF
          END IF! CI in inner iterations
*
*. Determine step length along direction
*. ======================================
*
*. Energy for rotated orbitals
*
            ONE = 1.0D0
            CALL VECSUM(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),dbl_mb(KLSTEP),
     &      ONE,ONE,NOOEXC)
            XNORM2 = SQRT(INPROD(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),
     &                           NOOEXC))
            IF(IPRNT.GE.2.and.ga_nodeid().eq.0) WRITE(6,'(A,E12.5)') 
     &      '  Norm of total kappa = ', XNORM2
            ENERGY1 = EMCSCF_FROM_KAPPA(dbl_mb(KLKAPPA))
            ENEW = ENERGY1
            IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &         WRITE(6,*) ' Energy for rotated orbitals', ENERGY1
*. Compare old and new energy to decide with to do
            DELTA_E_ACT = ENEW-EOLD
            IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &      WRITE(6,'(A,3(2X,E9.3))') 
     &      '  Actual energy change without linesearch ', DELTA_E_ACT
*
            IF((ABS(DELTA_E_ACT).GT.THRES_FOR_ENTER_LINSEA).AND.
     &         (I_DO_LINSEA_MCSCF.EQ.1.OR.
     &         I_DO_LINSEA_MCSCF.EQ.2.AND.EOLD.GT.ENEW)) THEN
*
*. line-search for orbital optimization
*
*. Step was added to Kappa when calculating energy, get Kappa back
              ONE = 1.0D0
              ONEM = -1.0D0
              CALL VECSUM(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),
     &                    dbl_mb(KLSTEP),ONE,ONEM,NOOEXC)
              CALL LINES_SEARCH_BY_BISECTION(EMCSCF_FROM_KAPPA,
     &             dbl_mb(KLKAPPA),dbl_mb(KLSTEP),NOOEXC,ONE,XFINAL,
     &             FFINAL,2, EOLD, ENEW)
              ENEW = FFINAL
              IF(IPRNT.GE.2.and.ga_nodeid().eq.0) WRITE(6,'(A,E9.3)') 
     &        '  Step-scaling parameter from lineseach = ', XFINAL
              XKAPPA_NORM2 = 
     &        SQRT(INPROD(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),NOOEXC))
              CALL VECSUM(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),
     &                    dbl_mb(KLSTEP),ONE, XFINAL,NOOEXC)
              DELTA_E_ACT = ENEW-EOLD
              IF(IPRNT.GE.2.and.ga_nodeid().eq.0) 
     &        WRITE(6,'(A,3(2X,E9.3))') 
     &        '  Actual energy change with  linesearch ', DELTA_E_ACT
            END IF! line search should be entered
*    
            IF(ABS(DELTA_E_ACT).LT.THRES_E_INNER) THEN
             if(ga_nodeid().eq.0) 
     &         WRITE(6,*) ' Inner iterations converged '
             CONV_INNER = .TRUE.
            END IF
*
            IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
               WRITE(6,*) ' Updated total Kappa '
               CALL WRTMAT(dbl_mb(KLKAPPA),1,NOOEXC,1,NOOEXC)
            END IF
            XNORM_IT = INPROD(dbl_mb(KLKAPPA),dbl_mb(KLKAPPA),NOOEXC)
            dbl_mb(KL_SUMMARY-1+(IOUT-1)*NITEM+2) = XNORM_IT
          END IF ! Update method
*
*=======================================
*. The new and improved MO-coefficients
*=======================================
*
*. Obtain exp(-kappa)
          CALL MEMCHK2('BE_NWM')
C              GET_EXP_MKAPPA(EXPMK,KAPPAP,IOOEXC,NOOEXC)
          CALL GET_EXP_MKAPPA(dbl_mb(KLMO1),dbl_mb(KLKAPPA),
     &                        int_mb(KLOOEXCC),NOOEXC)
          CALL MEMCHK2('AF_EMK')
*. 
          CALL MULT_BLOC_MAT(dbl_mb(KLMO2),dbl_mb(KMOREF),dbl_mb(KLMO1),
     &         NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
          CALL COPVEC(dbl_mb(KLMO2),dbl_mb(KMOMO),LEN_CMO)
          CALL MEMCHK2('AF_ML1')
*. And the new MO-AO coefficients
C?        WRITE(6,*) '  KMOAO_ACT = ', KMOAO_ACT
          CALL MULT_BLOC_MAT(dbl_mb(KMOAO_ACT),dbl_mb(KMOAOIN),
     &       dbl_mb(KMOMO),
     &       NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
          CALL MEMCHK2('AF_ML2')
          IF(IPRNT.GE.100.and.ga_nodeid().eq.0) THEN
            WRITE(6,*) ' Updated MO-coefficients'
            CALL APRBLM2(dbl_mb(KMOMO),NTOOBS,NTOOBS,NSMOB,0)
          END IF
          IF(IRESET_KAPPA_IN_OR_OUT.EQ.1) THEN
            CALL COPVEC(dbl_mb(KMOMO),dbl_mb(KMOREF),LEN_CMO)
          END  IF
          CALL MEMCHK2('AF_NWM')
*
*
*  ===========================================================
*. CI in inner its- should probably be moved (but not removed)
*  ===========================================================
*
          IF(I_MAY_DO_CI_IN_INNER_ITS.EQ.1.AND.I_SHIFT_E2.EQ.0.AND.
     &      XNORM2.LT.XKAPPA_THRES.AND.IOUT.GE.MIN_OUT_IT_WITH_CI) THEN
            IF(IPRNT.GE.10.and.ga_nodeid().eq.0) 
     &        WRITE(6,*) ' CI in inner it '
            I_DID_CI_IN_INNER = 1
C           WRITE(6,*) ' CI in inner it '
            I_DO_CI_IN_INNER_ITS = 10
*
*. Transform integrals to current set of MO's
*
            IF(IPRNT.GE.10.and.ga_nodeid().eq.0) 
     &         WRITE(6,*) ' Integral transformation:' 
            KINT2 = KINT2_INI
            IE2LIST_AL = IE2LIST_4F
C                DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,
C                          INTSM_IN)
            CALL DO_ORBTRA(1,1,0,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
CNW         CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
CNW         CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
CBERT Could probably just switch pointers
            call ga_copy(KINT1,KH)
            call ga_copy(KFI,KINT1)
*. The diagonal will fetch J and K integrals using GTIJKL_GN,* 
*. prepare for this routine
            IE2ARRAY_A = IE2LIST_I(IE2LIST_IB(IE2LIST_A))
*
*. Perform CI - and calculate densities
*
            IF(IPRNT.GE.10.and.ga_nodeid().eq.0) WRITE(6,*) ' CI: '
            IRESTR = 1
            MAXIT_SAVE = MAXIT
            MAXIT = 5
C           WRITE(6,*) ' Number of CI-iterations reduced to 1 '
            CALL GASCI(IREFSM,IREFSPC_MCSCF,IPRDIAL,IIUSEH0P,
     &           MPORENP_E,EREF,ERROR_NORM_FINAL,CONV_F)  
            MAXIT = MAXIT_SAVE
            if(ga_nodeid().eq.0)
     &      WRITE(6,*) ' Energy and residual from CI :', 
     &      EREF,ERROR_NORM_FINAL
            ENEW  = EREF
          END IF! CI in inner iterations
*
*. Obtain and block diagonalize FI+FA
*
          I_DIAG_FIFA = 0
          IF(I_DIAG_FIFA.EQ.1) THEN
*. Obtain FI +FA
            CALL DO_ORBTRA(0,1,1,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
CBERT INCORRECT DIMENSIONALITY OF KLMO1 VS FKI/KFA!!!!
            stop 'vecsum not done'
            CALL VECSUM(dbl_mb(KLMO1),KFI,KFA,ONE,ONE,
     &                  NINT1)
*. Diagonalize FI+FA and save in KLMO2
            CALL DIAG_GASBLKS(dbl_mb(KLMO1),dbl_mb(KLMO2),
     &           IDUM,IDUM,IDUM,dbl_mb(KLMO3),dbl_mb(KLMO4),2)
*. And new MO-coefficients
            CALL MULT_BLOC_MAT(dbl_mb(KLMO3),dbl_mb(KMOMO),
     &           dbl_mb(KLMO2),
     &           NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
            CALL COPVEC(dbl_mb(KLMO3),dbl_mb(KMOMO),LEN_CMO)
          END IF !FIFA should be diagonalized
*
          IF(CONV_INNER.AND.I_DO_CI_IN_INNER_ITS.EQ.1) THEN
            CONVER = .TRUE.
            GOTO 1001
          END IF
          IF(CONV_INNER) GOTO 901
        END DO !End of loop over inner iterations
 901    CONTINUE
        CALL MEMCHK2('EN_OUT')
      END DO
*     ^ End of loop over outer iterations
 1001 CONTINUE
      if(ga_nodeid().eq.0) then
      IF(CONVER) THEN
        WRITE(6,*) 
     &  ' Convergence of MCSCF was obtained in ', NOUTIT,' iterations'
      ELSE
        WRITE(6,*) 
     &  ' Convergence of MCSCF was not obtained in ', NOUTIT, 
     &  'iterations'
      END IF
      WRITE(6,'(A,I4)') 
     &'  Total number of inner iterations ', N_INNER_TOT
      endif
*
*
*. Finalize: Transform integrals to final MO's, obtain
*  norm of CI- and orbital gradient
*
*
*. Expansion of final orbitals in AO basis, pt in KLMO2
      MAXIT_SAVE = MAXIT
cVOG  IF(NCISPC.EQ.2.AND.NOMOFL.EQ.0.AND.IFINMO.NE.0) THEN
      IF(NOMOFL.EQ.0.AND.IFINMO.NE.0) THEN
       CALL MOROT_GS(IFINMO)
       call copvec(dbl_mb(KMOAOUT),dbl_mb(KMOAO_ACT),len_cmo)
       call copvec(dbl_mb(KMOAOUT),dbl_mb(KLMO2),len_cmo)
       IFINMO = 0
      ELSE
      CALL MULT_BLOC_MAT(dbl_mb(KLMO2),dbl_mb(KMOAOIN),dbl_mb(KMOMO),
     &       NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
      CALL COPVEC(dbl_mb(KLMO2),dbl_mb(KMOAO_ACT),LEN_CMO)
      CALL COPVEC(dbl_mb(KLMO2),dbl_mb(KMOAOUT),LEN_CMO)
*. And 0 CI iterations with new integrals
      MAXIT = 1
      END IF
c     CALL MULT_BLOC_MAT(dbl_mb(KLMO2),dbl_mb(KMOAOIN),dbl_mb(KMOMO),
c    &       NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
c     CALL COPVEC(dbl_mb(KLMO2),dbl_mb(KMOAO_ACT),LEN_CMO)
c     CALL COPVEC(dbl_mb(KLMO2),dbl_mb(KMOAOUT),LEN_CMO)
      if(ga_nodeid().eq.0) WRITE(6,*) 
     &' Final MO-AO transformation stored in MOAOIN, MOAO_ACT, MOAOUT'
*. Integral transformation
      KINT2 = KINT2_INI
*. Flag for integrals with one free index: energy + gradient
      IE2LIST_AL = IE2LIST_1F
      IE2LIST_AL = IE2LIST_4F
      IOCOBTP_AL = 1
      INTSM_AL = 1
*. Integral transf and FI 
      CALL DO_ORBTRA(1,1,0,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
CNW   CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
CNW   CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
CBERT depending on what happens above we could just switch pointers
      call ga_copy(KINT1,KH)
      call ga_copy(KFI,KINT1)
      IRESTR = 1
*. and normal density print
      IPRDEN = IPRDEN_SAVE 
      last_ci = 1
      CALL GASCI(IREFSM,IREFSPC_MCSCF,IPRDIAL,IIUSEH0P,
     &            MPORENP_E,EREF,ERROR_NORM_FINAL_CI,CONV_F)
      EFINAL = EREF
      MAXIT = MAXIT_SAVE
*. Current orbital gradient
!*. VOG 4 Density
!      only_dens = 1
!      LRHO2 = NACOB**2*(NACOB**2+1)/2
!      CALL MEMMAN(KRHO1BU,NACOB**2,'ADDL  ',2,'RHO1BU  ')
!      CALL MEMMAN(KRHO2BU,LRHO2   ,'ADDL  ',2,'RHO2BU  ')
!      call reduce_gas_spaces(trgas,dbl_mb(KRHO1),dbl_mb(krho1bu),
!     &                       dbl_mb(KRHO2),dbl_mb(krho2bu),only_dens)
      CALL DO_ORBTRA(0,0,1,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
!      call restore_gas_spaces(dbl_mb(KRHO1),dbl_mb(krho1bu),
!     &                        dbl_mb(KRHO2),dbl_mb(krho2bu),only_dens)
      CALL FOCK_MAT_STANDARD(KF,2,KINT1,KFA)
      IF(IPRNT.GE.100.and.ga_nodeid().eq.0) WRITE(6,*) ' F constructed '
      CALL E1_FROM_F(dbl_mb(KLE1),KF,1,int_mb(KLOOEXC),
     &               int_mb(KLOOEXCC),
     &               NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
*. VOG 4 Destiny End
      IF(I_AVERAGE_ORBEXC.EQ.1) THEN
*. Average over orbital excitations belonging to a given shell excitation
        CALL SHELL_AVERAGE_ORBEXC(dbl_mb(KLE1),NSSEX,int_mb(KNOOFSS),
     &       int_mb(KIBOOFSS),int_mb(KIOOFSS),dbl_mb(KLEXCSCR),NOOEXC,1)
        END IF
      E1NRM_ORB = SQRT(INPROD(dbl_mb(KLE1),dbl_mb(KLE1),NOOEXC))
      VNFINAL = E1NRM_ORB + ERROR_NORM_FINAL_CI
*
      IF(IPRORB.GE.2.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) 
     &  ' Final MOs in initial basis (not natural or canonical)'
        CALL APRBLM2(dbl_mb(KMOMO),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      IF(IPRORB.GE.1.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) 
     &  ' Final MOs in AO basis (not natural or canonical)'
        CALL PRINT_CMOAO(dbl_mb(KLMO2))
      END IF
*
*. Projection of final occupied orbitals on initial set of occupied orbitals
*
*. Obtain initial and final occupied orbitals
      ISCR(1) = 0
      ISCR(2) = NGAS
      CALL MEMMAN(KLCOCC_INI,LEN_CMO,'ADDL  ',2,'COCC_IN')
      CALL MEMMAN(KLCOCC_FIN,LEN_CMO,'ADDL  ',2,'COCC_FI')
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,
      CALL CSUB_FROM_C(dbl_mb(KMOAOIN),dbl_mb(KLCOCC_INI),NOCOBS,
     &                 ISCR_NTS,2,ISCR,0)
      CALL CSUB_FROM_C(dbl_mb(KLMO2),dbl_mb(KLCOCC_FIN),NOCOBS,ISCR_NTS,
     &                 2,ISCR,0)
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,IONLY_DIM)
      if(ga_nodeid().eq.0) WRITE(6,*) 
     &' Projecting final (MO2) on initial (MO1) occupied orbitals'
      CALL PROJ_ORBSPC_ON_ORBSPC(dbl_mb(KLCOCC_INI),dbl_mb(KLCOCC_FIN),
     &     NOCOBS,NOCOBS)
      if(ga_nodeid().eq.0) WRITE(6,*)
C     PROJ_ORBSPC_ON_ORBSPC(CMOAO1,CMOAO2,NMO1PSM,NMO2PSM)
*
*. Projection of final active orbitals on initial set of active orbitals
*
*. Obtain initial and final active orbitals
      ISCR(1) = NGAS
C     CALL MEMMAN(KLCOCC_INI,LEN_CMO,'ADDL  ',2,'COCC_IN')
C     CALL MEMMAN(KLCOCC_FIN,LEN_CMO,'ADDL  ',2,'COCC_FI')
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,
      CALL CSUB_FROM_C(dbl_mb(KMOAOIN),dbl_mb(KLCOCC_INI),NACOBS,
     &                 ISCR_NTS,1,ISCR,0)
      CALL CSUB_FROM_C(dbl_mb(KLMO2),dbl_mb(KLCOCC_FIN),NACOBS,ISCR_NTS,
     &                 1,ISCR,0)
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,IONLY_DIM)
      if(ga_nodeid().eq.0) WRITE(6,*) 
     &' Projecting final (MO2) on initial (MO1) active orbitals'
      CALL PROJ_ORBSPC_ON_ORBSPC(dbl_mb(KLCOCC_INI),dbl_mb(KLCOCC_FIN),
     &     NACOBS,NACOBS)
      if(ga_nodeid().eq.0) then
      WRITE(6,*)
C     PROJ_ORBSPC_ON_ORBSPC(CMOAO1,CMOAO2,NMO1PSM,NMO2PSM)
*. Print summary
      CALL PRINT_MCSCF_CONV_SUMMARY(dbl_mb(KL_SUMMARY),NOUTIT)
      WRITE(6,'(A,F20.12)') ' Final energy = ', EFINAL
      WRITE(6,'(A,F20.12)') ' Final norm of orbital gradient = ', 
     &                        E1NRM_ORB
      endif
*
C?    WRITE(6,*) ' E1NRM_ORB, ERROR_NORM_FINAL_CI = ',
C?   &             E1NRM_ORB, ERROR_NORM_FINAL_CI
C?    WRITE(6,*) ' Final energy = ', EFINAL

cVOG  Open file where backup information will be stored for the second
cVOG  LUCIA call for AkBk
      IF(ga_nodeid().eq.0) THEN
       iunit = 668
       open(iunit,file='ecore_bu',status='unknown')
       write(iunit,*) ECORE
       write(iunit,*) ECORE_ORIG
       write(iunit,*) ECORE_H
       write(iunit,*) ECORE_HEX
       write(iunit,*) ECORE_AB
       write(iunit,*) ECORE_INI
       write(iunit,*) ECORE_EXT
       write(iunit,*) ECORE_INA
       write(iunit,*) ECORE_FI
       close(iunit)
      end if
cVOG
      CALL MEMMAN(IDUMMY, IDUMMY, 'FLUSM', IDUMMY,'MCSCF ') 
      if(.not.ga_destroy(KINT1_INI))
     &   call errquit('cannot destroy: kint1_ini',0,GA_ERR)
      CALL LUCIAQEXIT('MCSCF')
      RETURN
      END
      SUBROUTINE SETDIA_BLM(B,VAL,NBLK,LBLK,IPCK)
*
* Set a blocked matrix to a diagonal matrix with diagonal values VAL
* (and off diagonal  elements = 0)
*
*. Jeppe Olsen, April 2010
*
      INCLUDE 'wrkspc.inc'
*. input
      INTEGER LBLK(NBLK)
*. Output
      DIMENSION B(*)
*
      NTEST = 00
*
      IOFF = -1
      DO IBLK = 1, NBLK
       LEN = LBLK(IBLK)
C?     WRITE(6,*) ' IBLK, LEN =', IBLK, LEN
       IF(IBLK.EQ.1) THEN
         IOFF = 1
       ELSE
         IF(IPCK.EQ.1) THEN
           IOFF = IOFF + LBLK(IBLK-1)*(LBLK(IBLK-1)+1)/2
         ELSE
           IOFF = IOFF + LBLK(IBLK-1)**2
         END IF
       END IF
       IF(IPCK.EQ.0) THEN
         LENB = LEN**2
       ELSE
         LENB = LEN*(LEN+1)/2
       END IF
       ZERO = 0.0D0
C?     WRITE(6,*) ' IOFF = ', IOFF
       CALL SETVEC(B(IOFF),ZERO,LENB)
       CALL SETDIA(B(IOFF),VAL,LEN,IPCK)
      END DO
*
      IF(NTEST.GE.100) THEN
C APRBLM2(A,LROW,LCOL,NBLK,ISYM)
        WRITE(6,*) ' Output matrix from SETDIA_BLM'
        CALL APRBLM2(B,LBLK,LBLK,NBLK,IPCK)
      END IF
*
      RETURN
      END
      FUNCTION EMCSCF_FROM_KAPPA(XKAPPA)
*
* Obtain MCSCF energy for orbital rotations defined by kappa
*
*. Notice: Integrals and inactive Fock matrix is new basis on return.
*
*. Jeppe Olsen, April 2010
*
*. Last revision; Oct 2012: Jeppe Olsen; clean up + CI
*
* NOTE: The reference orbitals- defining the orbitals together with Kappa
*       are now required to reside in KMOREF (and not in KMOMO) - Nov. 2011
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      LOGICAL CONV_FL
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
*
*. Some indirect transfer
      COMMON/EXCTRNS/KLOOEXCC,KINT1_INI,KINT2_INI, IREFSPC_MCSCFL,
     &               IPRDIALL,IIUSEH0PL,MPORENP_EL,
     &               ERROR_NORM_FINALL,CONV_FL,
     &               I_DO_CI_IN_INNER_ACT
*. Orbital rotations in compact form
      DIMENSION  XKAPPA(*)
*
      NTEST = 000
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',ADDL,'GTEMC ')
*. Space for matrices for MO transformations
      LEN_CMO =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLMO1,LEN_CMO,'ADDL  ',2,'LMO1  ')
      CALL MEMMAN(KLMO2,LEN_CMO,'ADDL  ',2,'LMO2  ')
      CALL MEMMAN(KLMO3,LEN_CMO,'ADDL  ',2,'LMO3  ')
*. Save current set of MO's
      CALL COPVEC(dbl_mb(KMOMO),dbl_mb(KLMO3),LEN_CMO)
*. Exp(-Kappa)
C GET_EXP_MKAPPA(EXPMK,KAPPAP,IOOEXC,NOOEXC)
      CALL GET_EXP_MKAPPA(dbl_mb(KLMO1),XKAPPA,
     &                    int_mb(KLOOEXCC),NOOEXC)
*. And new MO-coefficients: MO2 = MO Exp(-Kappa)
      CALL MULT_BLOC_MAT(dbl_mb(KLMO2),dbl_mb(KMOREF),dbl_mb(KLMO1),
     &     NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
      CALL COPVEC(dbl_mb(KLMO2),dbl_mb(KMOMO),LEN_CMO)
*
*. Integral transformation and FI for current MO expansion
*
*. Flag for integrals with one free index: energy + gradient
      IE2LIST_0F = 1
      IE2LIST_4F = 5
      IE2LIST_AL = IE2LIST_0F
      IE2LIST_AL = IE2LIST_4F
      IOCOBTP_AL = 1
      INTSM_AL = 1
      CALL DO_ORBTRA(1,1,0,IE2LIST_AL,IOCOBTP_AL,INTSM_AL)
CNW   CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
CNW   CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
      call ga_copy(KINT1,KH)
      call ga_copy(KFI,KINT1)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ECORE_ORIG, ECORE_HEX, ECORE(2) ',
     &               ECORE_ORIG, ECORE_HEX, ECORE
       END IF
*. CI if requested - calculations of densities may be eliminated
      IF(I_DO_CI_IN_INNER_ACT.EQ.1) THEN
        CALL GASCI(IREFSM,IREFSPC_MCSCFL,IPRDIALL,IIUSEH0PL,
     &             MPORENP_EL,EREFL,ERROR_NORM_FINALL,CONV_FL)  
        WRITE(6,*) ' Energy from CI = ', EREFL
      END IF
*
*. Energy for these MO-coefficients and densities
*
      CALL EN_FROM_DENS(ENERGY,2,0)
*. Clean up time: restore MO coefficients - but not integrals and FI
      CALL COPVEC(dbl_mb(KLMO3),dbl_mb(KMOMO),LEN_CMO)
*. And the conclusion
      EMCSCF_FROM_KAPPA = ENERGY
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',ADDL,'GTEMC ')
      RETURN
      END
      SUBROUTINE GENERIC_GRA_HES_FD(E0,E1,E2,X,NX,EFUNC)
*
* Obtain gradient and Hessian for for a general function 
* depending on NX parameters in X
*
* The function values are obtained through an external function EFUNC
*
*. Jeppe Olsen, April 2010
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input: current set of parameters
      DIMENSION X(NX)
*. Output
      DIMENSION E1(NX),E2(NX,NX)
      EXTERNAL EFUNC
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GEN_FD')
      CALL MEMMAN(KLX,NX,'ADDL  ',2,'XLOCAL')
*
      NTEST = 100
*.
*.  Energy at point of expansion
      E0 =  EFUNC(X)
      WRITE(6,*) ' Energy at reference', E0
*. step for finite difference 
      DELTA = 0.0010D0
*. Gradient and diagonal Hessian elements 
      DO I = 1, NX
* E(+Delta)
        CALL COPVEC(X,dbl_mb(KLX),NX)
        dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) + DELTA
        EP1 =  EFUNC(dbl_mb(KLX))
*. E(-Delta)
        CALL COPVEC(X,dbl_mb(KLX),NX)
        dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) - DELTA
        EM1 =  EFUNC(dbl_mb(KLX))
*. E(+2 Delta)
        CALL COPVEC(X,dbl_mb(KLX),NX)
        dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) + 2.0D0*DELTA
        EP2 =  EFUNC(dbl_mb(KLX))
*. E(-2 Delta)
        CALL COPVEC(X,dbl_mb(KLX),NX)
        dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) -2.0D0*DELTA
        EM2 =  EFUNC(dbl_mb(KLX))
*. And we can obtain gradient and diagonal elements
C?      WRITE(6,*) ' E0, EP1, EP2, EM1, EM2 = ',
C?   &               E0, EP1, EP2, EM1, EM2
*
        E1(I) 
     &  = (8.0D0*EP1-8.0D0*EM1-EP2+EM2)/(12.0D0*DELTA)
*
C?       WRITE(6,*) ' ID, ID_EFF, E1(ID_EFF) = ',
C?   &                ID, ID_EFF, E1(ID_EFF)
        E2(I,I)  
     &  =  (16.0D0*(EP1+EM1-2.0D0*E0)-EP2-EM2+2.0D0*E0)/
     &     (12.0D0*DELTA**2)
      END DO
*
*. And the non-diagonal Hessian elements
*
      DO I = 1, NX
        DO J = 1, I-1
*EP1P1
*
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) + DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) + DELTA
          EP1P1 =  EFUNC(dbl_mb(KLX))
*EM1M1
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) -DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) -DELTA
          EM1M1 =  EFUNC(dbl_mb(KLX))
*EP1M1
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) + DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) - DELTA
          EP1M1 =  EFUNC(dbl_mb(KLX))
*EM1P1
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) - DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) + DELTA
          EM1P1 =  EFUNC(dbl_mb(KLX))
*EP2P2
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) + 2.0D0*DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) + 2.0D0*DELTA
          EP2P2 =  EFUNC(dbl_mb(KLX))
*EM2M2
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) - 2.0D0*DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) - 2.0D0*DELTA
          EM2M2 =  EFUNC(dbl_mb(KLX))
*EP2M2
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) + 2.0D0*DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) - 2.0D0*DELTA
          EP2M2 =  EFUNC(dbl_mb(KLX))
*EM2P2
          CALL COPVEC(X,dbl_mb(KLX),NX)
          dbl_mb(KLX-1+I) = dbl_mb(KLX-1+I) - 2.0D0*DELTA
          dbl_mb(KLX-1+J) = dbl_mb(KLX-1+J) + 2.0D0*DELTA
          EM2P2 =  EFUNC(dbl_mb(KLX))
*
          G1 = EP1P1-EP1M1-EM1P1+EM1M1 
          G2 = EP2P2-EP2M2-EM2P2+EM2M2 
*
          E2(I,J) = (16.0D0*G1-G2)/(48*DELTA**2)
          E2(J,I) = E2(I,J) 
*
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Gradient obtained by finite difference '
        WRITE(6,*) ' ======================================='
        CALL WRTMAT(E1,1,NX,1,NX)              
        WRITE(6,*)
        WRITE(6,*) ' Hessian obtained by finite difference '
        WRITE(6,*) ' ======================================'
        CALL WRTMAT(E2,NX,NX,NX,NX)                   
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GEN_FD')
      RETURN
      END 
      SUBROUTINE PRINT_MCSCF_CONV_SUMMARY(SUMMARY,NIT)
*
* Print summary of MCSCF calculations
*
*. Jeppe Olsen, April 2010
*
*. Last modification; Oct. 2012; Jeppe Olsen, allowing for +999 its
      INCLUDE 'implicit.inc'
*. Current number of items per iteration
      PARAMETER(NITEM = 4)
*
      DIMENSION SUMMARY(NITEM,NIT)
*
      WRITE(6,*)
      WRITE(6,*) ' Summary of MCSCF convergence: '
      WRITE(6,*) ' =============================='
      WRITE(6,*) 
     & ' Iter Orb-gradient  Orb-step  CI-gradient     Energy '
      WRITE(6,*) 
     & ' =========================================================='
      DO IT = 1, NIT
        WRITE(6,'(I4,4X,E8.3,4X,E8.3,3X,E8.3,F20.12)')
     &  IT, (SUMMARY(ITEM,IT),ITEM =1, NITEM)
      END DO
      WRITE(6,*)
      WRITE(6,*)
*
      RETURN
      END
      SUBROUTINE CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,
     &                       IONLY_DIM)
*
* Obtain MO-INI transformation for a subset of molecular orbitals
* IF IONLY_DIM .ne. 0, then only the dimensions of the subset
* is generated.
*
* The types of the subsets are defined by ISUBTP
*
*. Jeppe Olsen, October 2010, Modified May 2011
*
*. General input
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      DIMENSION C(*)
      INTEGER ISUBTP(NSUBTP)
*
* ISUBTP = 0       => Inactive
* ISUBTP = NGAS    => Active
* ISUBTP = NGAS + 1=> Secondary
*
*. Output
      DIMENSION CSUB(*) 
      INTEGER LENSUBS(*),LENSUBTS(0:6+MXPR4T,MXPOBS)
*
      NTEST = 00
*
      IZERO = 0
      CALL ISETVC(LENSUBTS,IZERO,(7+MXPR4T)*MXPOBS)
*
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
         IOBOFF = 1
         ICOFF = 1
         ICSOFF = 1
        ELSE
         IOBOFF = IOBOFF + NTOOBS(ISM-1)
         ICOFF = ICOFF + NTOOBS(ISM-1)**2
         ICSOFF = ICSOFF + LENSUBS(ISM-1)*NTOOBS(ISM-1)
        END IF
        LENS = NTOOBS(ISM)
        IF(NTEST.GE.1000) WRITE(6,*) ' ISM, LENS = ', ISM, LENS
        LENSUBS(ISM) = 0
*
*. Dimensions of the various orbitalsubspaces 
*
*. Loop over inactive/active/secondary
        DO IAS = 1, 3
          IF(IAS.EQ.1) THEN
*. Looking for inactive
           ITARGET = 0
           ISTART = 0
           ISTOP = 0
          ELSE IF(IAS.EQ.2) THEN
*. Looking for active
           ISTART = 1
           ISTOP = NGAS
           ITARGET = NGAS
          ELSE 
*. Looking for secondary
           ITARGET = NGAS + 1
           ISTART = NGAS + 1
           ISTOP = NGAS + 1
          END IF
          IACTIVE = 0
          DO JSUBTP = 1, NSUBTP
            IF(ISUBTP(JSUBTP).EQ.ITARGET) IACTIVE = 1
          END DO 
C?        WRITE(6,*) ' ITARGET, IACTIVE = ', ITARGET,IACTIVE
          IF(IACTIVE.EQ.1) THEN
            DO IGAS = ISTART, ISTOP
              LENSUBTS(IGAS,ISM) = NOBPTS_GN(IGAS,ISM)
            END DO
          END IF
        END DO ! End of loop over IAS
*
COLD    DO IGAS = 0, NGAS + 1
        DO IORB = 1, LENS
         ITP = ITPFSO(IOBOFF-1+IORB)
*. Modify GAS-types to NGAS
         IF(0.LT.ITP.AND.ITP.LE.NGAS) ITP = NGAS 
         I_AM_RIGHT_TYPE = 0
         DO JSUBTP = 1, NSUBTP
          IF(ISUBTP(JSUBTP).EQ.ITP) I_AM_RIGHT_TYPE = 1
         END DO
         IF(I_AM_RIGHT_TYPE.EQ.1) THEN
           LENSUBS(ISM) = LENSUBS(ISM) + 1
           ICOFF2 = ICOFF + (IORB-1)*LENS
           ICSOFF2 = ICSOFF + (LENSUBS(ISM)-1)*LENS
C?         WRITE(6,*) ' Orbital included, IORB = ', IORB
C?         WRITE(6,*) ' ICOFF, ICSOFF = ', ICOFF,ICSOFF
C?         WRITE(6,*) ' ICOFF2, ICSOFF2 = ',ICOFF2, ICSOFF2 
           IF(IONLY_DIM.EQ.0)
     &     CALL COPVEC(C(ICOFF2),CSUB(ICSOFF2),LENS)
         END IF
        END DO
COLD    END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from CSUB_FROM_C '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,'(A,3I3)') ' Requested types:',
     &  (ISUBTP(I),I=1,NSUBTP)
        WRITE(6,*) ' Number of MOs per symmetry in CSUB '
        CALL IWRTMA(LENSUBS,NSMOB,1,NSMOB,1)
        WRITE(6,*) ' Number of MOs per type and sym '
        CALL IWRTMA(LENSUBTS,NGAS+2,NSMOB,7+MXPR4T,NSMOB)
*
        IF(IONLY_DIM.EQ.0) THEN
          WRITE(6,*) ' Resulting CSUB'
C              APRBLM2(A,LROW,LCOL,NBLK,ISYM)
          CALL APRBLM2(CSUB,NTOOBS,LENSUBS,NSMOB,0)
        END IF 
      END IF
*
      RETURN
      END 
      SUBROUTINE PREPARE_2EI_LIST
*
* Prepare for using two-electron integral list IE2LIST_A, IOCOBTP_A,INTSM_A
* (from cintfo.inc)  
*
* i.e.  set up relevant arrays - which are assumed to have been allocated
*
* Jeppe Olsen, April 2011, for the LUCIA growing up campaign
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
*. Local scratch
      INTEGER NOCOBS_L(MXPOBS),NOBS_L(MXPOBS,4), ISUBTP(2)
      INTEGER NOCOBTS_L(MXPOBS*(7+MXPR4T))

*
      INCLUDE 'cintfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csmprd.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'intform.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
      NTEST = 00
*
      IE2LIST_N_A =IE2LIST_N(IE2LIST_A)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' PREPA... IE2LIST_A, IE2LIST_N_A = ',
     &                        IE2LIST_A, IE2LIST_N_A
      END IF
*
      IB = IE2LIST_IB(IE2LIST_A)
      DO IARR = 1, IE2LIST_N_A
       IE2LIST_I_A(IARR) = IE2LIST_I(IB-1+IARR)
      END DO
*
*. Complex conjugation symmetry of the one-electron integrals in a form used by
*. integral fetch routines
*
      IF(IE1_CCSM_G(IE2LIST_A).EQ.1) THEN
*. Permutational symmetry of integrals
       IH1FORM = 1
      ELSE
*. No permutational symmetry of integrals
       IH1FORM = 2
      END IF
*. Some routines also need to know about the complex conjugation symmetry
* of two-electron integrals so:
      IF(IE2_CCSM_G(IE2LIST_A).EQ.1) THEN
*. Permutational symmetry of integrals
       IH2FORM = 1
      ELSE
*. No permutational symmetry of integrals
       IH2FORM = 2
      END IF
C?    WRITE(6,*) ' PREPARE.. IH1FORM, IH2FORM = ', IH1FORM, IH2FORM
*
*. Number of occupied per symmetry
*
      IF(IOCOBTP_A.EQ.1) THEN
        NSUBTP = 1
        ISUBTP(1) = NGAS
      ELSE
        NSUBTP = 2
        ISUBTP(1) = 0
        ISUBTP(2) = NGAS
      END IF
      CALL CSUB_FROM_C(XDUM,XDUM,NOCOBS_L,NOCOBTS_L,NSUBTP,ISUBTP,1)
*
      DO IARR = 1, IE2LIST_N_A
        IIARR = IE2LIST_I_A(IARR)
        DO INDEX = 1, 4
          IF(INT2ARR_G(INDEX,IIARR).EQ.1) THEN
            CALL ICOPVE(NOCOBS_L,NOBS_L(1,INDEX),NSMOB)
          ELSE
            CALL ICOPVE(NTOOBS,NOBS_L(1,INDEX),NSMOB)
          END IF
        END DO
        I12S_L = I12S_G(IIARR)
        I34S_L = I34S_G(IIARR)
        I1234S_L = I1234S_G(IIARR)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 
     &   ' Before call to PNT4DM: IIARR, KPINT2_A, KPLSM2_A = ',
     &    IIARR,KPINT2_A(IIARR),KPLSM2_A(IIARR)
        END IF
        CALL PNT4DM(NSMOB,NSMSX,MXPOBS,
     &       NOBS_L(1,1),NOBS_L(1,2),NOBS_L(1,3),
     &       NOBS_L(1,4),INTSM_A,ADSXA,SXDXSX,I12S_L,I34S_L,I1234S_L,
     &       int_mb(KPINT2_A(IIARR)), int_mb(KPLSM2_A(IIARR)),
     &       ADASX,NINT4D)
        IE2ARR_L_A(IIARR) = NINT4D
      END DO
*
      RETURN
      END
      SUBROUTINE FLAG_ACT_INTLIST(IACT_LIST)
*
* Flag that only the integral arrays of integral list IACT_LIST is active
* This is realized by setting the pointers to all other integral arrays
* to their negative value
*
* Jeppe Olsen, May 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      NTEST = 000
*
      DO IE2LIST = 1, NE2LIST
        N_IA = IE2LIST_N(IE2LIST)
        IB_IA = IE2LIST_IB(IE2LIST)
        DO II_AR = IB_IA, IB_IA-1+N_IA
          I_AR = IE2LIST_I(II_AR)
          IF(IE2LIST.EQ.IACT_LIST) THEN
           KINT2_A(I_AR) = IABS(KINT2_A(I_AR))
          ELSE
           KINT2_A(I_AR) = -IABS(KINT2_A(I_AR))
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from FLAG_ACT_INTLIST'
        WRITE(6,*) ' Integral list to be flagged positive ', IACT_LIST
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Modified integral pointers '
        CALL IWRTMA(KINT2_A,1,NE2ARR,1,NE2ARR)
      END IF
*
      RETURN
      END
      SUBROUTINE E1_MCSCF_FOR_GENERAL_KAPPA(E1,BRT,KAPPA)
*
* A kappa matrix in  packed form and a Brilloun vector, BRT, in expanded from are given. 
* Obtain gradient and save in E1
*
*. Jeppe Olsen October 2011 - Reusing old LUCAS routines
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc' 
      INCLUDE 'crun.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*. Input
      REAL*8 KAPPA(*), BRT(*)
*. Output
      REAL*8 E1(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'E1_GEN')
*. Matrix for Brilloin vector for  general/occupied index
      N_GG = IINPROD(NTOOBS,NTOOBS,NSMOB)
      N_GO = IINPROD(NTOOBS,NOCOBS,NSMOB)
      CALL MEMMAN(KLE1,N_GO,'ADDL  ',2,'E1_EXP')
*. Obtain Kappa in expanded form
      CALL MEMMAN(KLKAPPAE,N_GG,'ADDL  ',2,'KAP_E  ')
      CALL REF_AS_KAPPA(KAPPA,dbl_mb(KLKAPPAE),1,1,int_mb(KIOOEXCC),
     &                  NOOEXC)
*. If we want to do FUSK (explicit calculation of gradient)
      I_DO_FUSK_E1 = 0
      IF(I_DO_FUSK_E1.EQ.1) THEN
*. Scratch space for fusk calculation, and storage of result
        CALL MEMMAN(KLE1_FUSK_E,N_GG,'ADDL  ',2,'E1EXPF')
        CALL MEMMAN(KLE1_FUSK,NOOEXC,'ADDL  ',2,'E1FSK')
        MXSOB_BLK = IMNMX(NTOOBS,NSMOB,2)
        CALL MEMMAN(KLE1_FUSK_SCR,3*MXSOB_BLK**2,'ADDL  ',2,'E1FSCR')
      END IF
*
      IOFF_GO = 1
      IOFF_GG = 1
*
      DO ISM  = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF_GO = 1
          IOFF_GG = 1
        ELSE
          IOFF_GO = IOFF_GO + NTOOBS(ISM-1)*NOCOBS(ISM-1)
          IOFF_GG = IOFF_GG + NTOOBS(ISM-1)*NTOOBS(ISM-1)
        END IF
        NOCC = NOCOBS(ISM)
        NORB = NTOOBS(ISM)
*. Fusk. set elements Kappa(7,1),(4,1) to 1, rest to zero
C?      IF(ISM.EQ.1) THEN
C?        ZERO = 0.0D0
C?        CALL SETVEC(WORK(KLKAPPAE),ZERO,NORB**2)
C?        WORK(KLKAPPAE-1+(1-1)*NORB+7) = 1.0D0
C?        WORK(KLKAPPAE-1+(1-1)*NORB+4) = 1.0D0
C?        WORK(KLKAPPAE-1+(7-1)*NORB+1) = -1.0D0
C?        WORK(KLKAPPAE-1+(4-1)*NORB+1) = -1.0D0
C?      END IF
*. End of fusk
        CALL LINGRA_FOR_SYM(dbl_mb(KLKAPPAE-1+IOFF_GG),
     &       BRT(IOFF_GG),NOCC,NORB,dbl_mb(KLE1-1+IOFF_GO))
        IF(I_DO_FUSK_E1.EQ.1) THEN
*. Fusk(explicit calculation of gradient)
C             COMMUP(C,A,B,NDIM,SCR,ISKIP0)
         ONEM = -1.0D0
         CALL SCALVE(dbl_mb(KLKAPPAE-1+IOFF_GG),ONEM,NORB**2)
         CALL COMMUP(dbl_mb(KLE1_FUSK_E-1+IOFF_GO),BRT(IOFF_GG),
     &        dbl_mb(KLKAPPAE-1+IOFF_GG),NORB,dbl_mb(KLE1_FUSK_SCR),0)
         CALL SCALVE(dbl_mb(KLKAPPAE-1+IOFF_GG),ONEM,NORB**2)
        END IF
        
       END DO
*. We have now the gradient in the matrix WORK(KLE1) in expanded form
*. as a matrix with general/occupied index. Obtain the corresponding 
*. gradient in standard packed form
C     EXC_VEC_FROM_GO_MAT(EXC_VEC, GOMAT,IJSM,
C    &           NOOEXC,IOOEXCC,IOOEXC,
C    &           NSMOB,NOCOBS,NTOOBS,NTOOB,IBSO,IREOST)
       CALL EXC_VEC_FROM_GO_MAT(E1,dbl_mb(KLE1),1,NOOEXC,
     &      int_mb(KIOOEXCC),int_mb(KIOOEXC),NSMOB,NOCOBS,NTOOBS,
     &      NTOOB,IBSO,IREOST)
       IF(I_DO_FUSK_E1.EQ.1) THEN
*. Obtain also the Fusk gradient
       CALL EXC_VEC_FROM_GO_MAT(dbl_mb(KLE1_FUSK),dbl_mb(KLE1_FUSK_E),1,
     &      NOOEXC,int_mb(KIOOEXCC),int_mb(KIOOEXC),NSMOB,NOCOBS,NTOOBS,
     &      NTOOB,IBSO,IREOST)
       END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Gradient for general Kappa '
        CALL WRT_EXCVEC(E1,int_mb(KIOOEXCC),NOOEXC)
        IF(I_DO_FUSK_E1.EQ.1) THEN
          WRITE(6,*) ' Gradient(FUSK) for general Kappa '
          CALL WRT_EXCVEC(dbl_mb(KLE1_FUSK),int_mb(KIOOEXCC),NOOEXC)
        END IF
*
      END IF
*
      IF(I_DO_FUSK_E1.EQ.1) THEN
*. Compare the two gradients 
        WRITE(6,*)  ' Comparison of Analytical and Fusk gradients'
*. Use a threshold relevant for larger gradients
        THRES = 1.0D-15
        CALL CMP2VC(E1,dbl_mb(KLE1_FUSK),NOOEXC,THRES)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'E1_GEN')
*
      RETURN
      END
      SUBROUTINE LINGRA_FOR_SYM(KAPPA,BRT,NOCC,NORB,E1)
*
* For a given symmetry, calculate orbital gradient at a general point KAPPA,
* from Brillouin vector BRT, given in full form
*
* Simplified version of LUCAS code
*
* The Gradient is calculated as E1(R,S) where R is a general 
* orbital index and S is an occupied orbital
*
*=========
* Input :
*=========
*
* KAPPA : actual kappa parameters
* BRT   : Brillouin vector in current basis, expanded form: BRT(NORB,NOCC)
* NOCC : Number of occupied orbitals of this symmetry
* NORB : Number of orbitals of this symmetry
*
*=========
* Output :
*=========
*
* E1 : Gradient E1(NORB,NOCC)
*
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
C     REAL*8 INPROD
      REAL*8 KAPPA(NORB,NOCC)
      DIMENSION BRT(NORB,NORB),E1(NORB,NOCC)
** test arrays
      dimension xjep1(700),xjep2(700),xjep3(700),xjep4(700)
COLD  dimension xjep5(700),xjep6(700),xjep7(700),xjep8(700)
*
      IDUM = 0
      ONE = 1.0D0
      ONEM= -1.0D0
      ZERO = 0.0D0
      TWO = 2.0D0
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LINGRA')
      NTEST = 000
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Output from LINGRA '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' Initial Brillouin matrix'
        CALL WRTMAT(BRT,NORB,NORB,NORB,NORB)
      END IF
* Dimension of nonsingular Kappa  matrix
      NRDIM  = MIN(NORB,2 * NOCC)
*
** 1 : some memory allocation
*
* Kappa in subspace
      CALL MEMMAN(KLKPR,NRDIM**2,'ADDL  ',2,'KAPRED')
* Vectors defining subspace
      CALL MEMMAN(KLV,NRDIM*NORB,'ADDL  ',2,'VVEC  ')
*. Brilloins vector in current basis (tilde basis)
COLD  CALL MEMMAN(KLBRT,NOCC*NORB,'ADDL  ',2,'BRTVEC')
*. The C-matrix
      CALL MEMMAN(KLCFORA,NRDIM**2,'ADDL  ',2,'CFORA  ')
*. Scratch for generation of C-matrix
      LEN_CFORA_SCR = 4*NRDIM**2 + 4*NRDIM
      CALL MEMMAN(KLCFORA_SCR,LEN_CFORA_SCR,'ADDL  ',2,'S_CFORA')
*. A matrix
      CALL MEMMAN(KLMAT1,NORB**2,'ADDL  ',2,'MAT1  ')
      CALL MEMMAN(KLMAT2,NORB**2,'ADDL  ',2,'MAT2  ')
*
* Brillouins vector: Is now input.... and is therefore commented out
*
* <0![E(RS)-E(SR),H]!0> = 2(F(RS) - F(SR))
* ========================================
*
C     DO IR = 1, NORB
C      DO IS = 1, NOCC
C       WORK(KLBRT-1+(IS-1)*NORB+IR) = -TWO*F(IS,IR)
C       IF(IR.LE.NOCC) WORK(KLBRT-1+(IS-1)*NORB+IR) =
C    &  WORK(KLBRT-1+(IS-1)*NORB+IR) + TWO*F(IR,IS)
C      END DO
C     END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Brillouin vector in original basis '
        CALL WRTMAT(BRT,NORB,NOCC,NORB,NOCC)
      END IF
* The Brillouin vector is the first contribution to gradient so
      CALL COPVEC(BRT,E1,NOCC*NORB)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Gradient including the first term '
        CALL WRTMAT(E1,NORB,NOCC,NORB,NOCC)
      END IF
*
** 2 : Obtain subspace for Kappa and kappa in subspace
*      Kappa = V Kappa(red) V(T)
*
COLD  CALL MEMMAN(KLSCR1,NRDIM**2,'ADDL  ',2,'SCRRED')
      CALL REDKAP(KAPPA,NOCC,NORB,NREDVC,dbl_mb(KLV),dbl_mb(KLKPR))
*. Update NRDIM
      NRDIM = NREDVC
C          REDKAP(KAPPA,NOCC,NORB,NREDVC,REDVEC,REDKP,SCR)
      IF( NTEST .GE. 10 ) THEN
         WRITE(6,*) ' Reduced Kappa matrix '
         WRITE(6,*) ' ==================== '
         CALL WRTMAT(dbl_mb(KLKPR),NRDIM,NRDIM,NRDIM,NRDIM)
      END IF
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' The V basis '
         WRITE(6,*) ' =========== '
         CALL WRTMAT(dbl_mb(KLV),NORB,NRDIM,NORB,NRDIM)
      END IF 
*
*
* ==============================================
* Second term to gradient: (1-VV(T)) B C(Kappa)
* ==============================================
*
* C(Kappa) = Sum(n) 1/(n+1)! Kappa^n
C          CFORA(C,A,NDIM,SCR)
      CALL CFORA(dbl_mb(KLCFORA),dbl_mb(KLKPR),NRDIM,
     &           dbl_mb(KLCFORA_SCR))
*. Obtain C in occupied subspace and save in MAT1
C COPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
      CALL COPMT2(dbl_mb(KLCFORA),dbl_mb(KLMAT1),NRDIM,NRDIM,
     &            NOCC,NOCC,1)
*. B C(Kappa) in MAT2
      FACTORAB = 1.0D0
      FACTORC  = 0.0D0
      CALL MATML7(dbl_mb(KLMAT2),BRT,dbl_mb(KLMAT1),
     &     NORB,NOCC,NORB,NOCC,NOCC,NOCC,FACTORC,FACTORAB,0)
      CALL VECSUM(E1,E1,dbl_mb(KLMAT2),ONE,ONE,NOCC*NORB)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' BC(Kappa) matrix '
        CALL WRTMAT(dbl_mb(KLMAT2),NORB,NOCC,NORB,NOCC)
      END IF
*. V V(T) B C in MAT2
      CALL MATML7(dbl_mb(KLMAT1),dbl_mb(KLV),dbl_mb(KLMAT2),
     &            NRDIM,NOCC,NORB,NRDIM,NORB,NOCC,
     &            FACTORC,FACTORAB,1)
      CALL MATML7(dbl_mb(KLMAT2),dbl_mb(KLV),dbl_mb(KLMAT1),
     &            NORB,NOCC,NORB,NRDIM,NRDIM,NOCC,
     &            FACTORC,FACTORAB,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' V V(T) BC(Kappa) matrix '
        CALL WRTMAT(dbl_mb(KLMAT2),NORB,NOCC,NORB,NOCC)
      END IF
      CALL VECSUM(E1,E1,dbl_mb(KLMAT2),ONE,ONEM,NOCC*NORB)
*
*. B C(Kappa) in MAT2
*
C?    FACTORAB = 1.0D0
C?    FACTORC  = 0.0D0
C?    CALL MATML7(WORK(KLMAT2),BRT,WORK(KLCFORA),
C?   &     NORB,NRDIM,NORB,NRDIM,NRDIM,NRDIM,FACTORC,FACTORAB,0)
C?    CALL VECSUM(E1,E1,WORK(KLMAT2),ONE,ONE,NOCC*NORB)
C?    IF(NTEST.GE.1000) THEN
C?      WRITE(6,*) ' BC(Kappa) matrix '
C?      CALL WRTMAT(WORK(KLMAT2),NORB,NRDIM,NORB,NRDIM)
C?    END IF
*. V V(T) B C in MAT2
C?    CALL MATML7(WORK(KLMAT1),WORK(KLV),WORK(KLMAT2),
C?   &            NRDIM,NRDIM,NORB,NRDIM,NORB,NRDIM,
C?   &            FACTORC,FACTORAB,1)
C?    CALL MATML7(WORK(KLMAT2),WORK(KLV),WORK(KLMAT1),
C?   &            NORB,NRDIM,NORB,NRDIM,NRDIM,NRDIM,
C?   &            FACTORC,FACTORAB,0)
C?    IF(NTEST.GE.1000) THEN
C?      WRITE(6,*) ' V V(T) BC(Kappa) matrix '
C?      CALL WRTMAT(WORK(KLMAT2),NORB,NRDIM,NORB,NRDIM)
C?    END IF
C?    CALL VECSUM(E1,E1,WORK(KLMAT2),ONE,ONEM,NOCC*NORB)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Gradient including the first two terms'
        CALL WRTMAT(E1,NORB,NOCC,NORB,NOCC)
      END IF
* 
* ==============================================
* Third term to gradient: X* (G BB) X(T)
* ==============================================
*
* U diagonalized kappa(red): kappa(red) = i U epsil_i U^\dagger
* X = VU gives Kappa in reduced diagonal basis: Kappa = i X epsil X^\dagger
* g(r,s) = sum_i=1^\infty i^n/(n+1)! (epsil_r-epsil_s)^n
* BRB is Brilloin vector in bar basis - i.e. the X-basis diagonalizing Kappa 
* G BRB (R,S) =  G(R,S) BRB(R,S)
* Eigenvalues, epsilon:
      CALL MEMMAN(KLER,NRDIM,'ADDL  ',2,'EIGVLR')
      CALL MEMMAN(KLEI,NRDIM,'ADDL  ',2,'EIGVLI')
* eigenvectors, UR, UI
      CALL MEMMAN(KLUR,NRDIM**2,'ADDL  ',2,'EIGVCR')
      CALL MEMMAN(KLUI,NRDIM**2,'ADDL  ',2,'EIGVCI')
* XR, XI
      CALL MEMMAN(KLXR,NRDIM*NORB,'ADDL  ',2,'XVECR ')
      CALL MEMMAN(KLXI,NRDIM*NORB,'ADDL  ',2,'XVECI ')
* GR, GI
      CALL MEMMAN(KLGR,NRDIM**2,'ADDL  ',2,'GMATR ')
      CALL MEMMAN(KLGI,NRDIM**2,'ADDL  ',2,'GMATI ')
      CALL MEMMAN(KLBRBI,NRDIM**2,'ADDL  ',2,'BRBI  ')
* A couple of matrices
      CALL MEMMAN(KLMAT3,NORB**2,'ADDL  ',2,'MAT3  ')
      CALL MEMMAN(KLMAT4,NORB**2,'ADDL  ',2,'MAT4  ')
      CALL MEMMAN(KLMAT5,NRDIM**2,'ADDL  ',2,'MAT5  ')
*. The Brillouin vector in the  bar basis (x-basis)
      CALL MEMMAN(KLBRBR,NRDIM**2,'ADDL  ',2,'BRXR  ')
      CALL MEMMAN(KLBRBI,NRDIM**2,'ADDL  ',2,'BRXI  ')
* G BRB
      CALL MEMMAN(KLGBRBR,NRDIM**2,'ADDL  ',2,'GBRBR ')
      CALL MEMMAN(KLGBRBI,NRDIM**2,'ADDL  ',2,'GBRBI ')

*
** Diagonalize Kappa(red): Kappa(red) = i U esilon(imag) U^{\dagger}
*
*
      CALL MEMMAN(KLZ,2*NRDIM**2,'ADDL  ',2,'Z_DIA ')
      CALL MEMMAN(KLW,2*NRDIM,'ADDL  ',2,'W_DIA ')
      CALL MEMMAN(KLSCRDIA,2*NRDIM,'ADDL  ',2,'SCRDIA')
*
C     EIGGMT(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,Z,W,SCR)
      CALL EIGGMTN(dbl_mb(KLKPR),NRDIM,dbl_mb(KLER),dbl_mb(KLEI),
     &           dbl_mb(KLUR),dbl_mb(KLUI),dbl_mb(KLZ),dbl_mb(KLW),
     &           dbl_mb(KLSCRDIA))
*. Check that X(dagger)KappaR X = i \epsilon
      I_DO_TEST = 0
      IF(I_DO_TEST.EQ.1) THEN
        WRITE(6,*) ' Checking U(dagger)KappaR U = i \epsilon '
*. Regenerate KappaR
        CALL REDKAP(KAPPA,NOCC,NORB,NREDVC,dbl_mb(KLV),dbl_mb(KLKPR))
        NRDIM = NREDVC
        WRITE(6,*) ' KappaR before call '
        CALL WRTMAT(dbl_mb(KLKPR),NRDIM,NRDIM,NRDIM,NRDIM)
        CALL SETVEC(XJEP1,0.0D0,NRDIM**2)
C       CTRN_MAT(XAXR,XAXI,XR,XI,AR,AI,NXR,NXC,NAR,NAC,SCR)
        CALL CTRN_MAT(XJEP2,XJEP3,dbl_mb(KLUR),dbl_mb(KLUI),
     &                dbl_mb(KLKPR),XJEP1,NRDIM,NRDIM,NRDIM,NRDIM,
     &                XJEP4)
         WRITE(6,*) ' U(Dagger) KappaR U '
         CALL WRTCMAT(XJEP2,XJEP3,NRDIM,NRDIM)
         WRITE(6,*) ' U after call '
         CALL WRTCMAT(dbl_mb(KLUR),dbl_mb(KLUI),NRDIM,NRDIM)
         WRITE(6,*) ' V after call '
         CALL WRTMAT(dbl_mb(KLV),NRDIM,NRDIM,NRDIM,NRDIM)
       END IF
* End of test zone
*
*. Matrices XR = V UR, XI = V UI - actual MO to reduced diagonal
*
      CALL MATML4(dbl_mb(KLXR),dbl_mb(KLV),dbl_mb(KLUR),
     &            NORB,NRDIM,NORB,NRDIM,
     &            NRDIM,NRDIM,0)
      CALL MATML4(dbl_mb(KLXI),dbl_mb(KLV),dbl_mb(KLUI),
     &            NORB,NRDIM,NORB,NRDIM,
     &            NRDIM,NRDIM,0)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The X matrix '
        WRITE(6,*) ' ============='
        CALL WRTCMAT(dbl_mb(KLXR),dbl_mb(KLXI),NORB,NRDIM)
      END IF
*
* test zone: test that X(dagger) Kappa X = i\epsilon
*
      IF(I_DO_TEST.EQ.1) THEN
        WRITE(6,*) ' Test that X(dagger)KappaR X = i \epsilon '
C       CTRN_MAT(XAXR,XAXI,XR,XI,AR,AI,NXR,NXC,NAR,NAC,SCR)
        CALL SETVEC(XJEP1,0.0D0,NRDIM**2)
        CALL CTRN_MAT(XJEP2,XJEP3,dbl_mb(KLXR),dbl_mb(KLXI),
     &                KAPPA,XJEP1,NRDIM,NRDIM,NRDIM,NRDIM,
     &                XJEP4)
         WRITE(6,*) ' X(Dagger) KappaR X '
         CALL WRTCMAT(XJEP2,XJEP3,NRDIM,NRDIM)
      END IF
*. End of test zone
*
*. Obtain Brilluoin matrix in X basis: B(bar) = X(T) B(Tilde) X*
*. B(bar)(rs) = sum_(r's') X(r'r) B(tilde)(r's') X*(s's)
* X(T) B: real part in MAT3, imag part in MAT4
      FACTORAB = 1.0D0
      FACTORC = 0.0D0
*. X(T)R B in MAT3
      CALL MATML7(dbl_mb(KLMAT3),dbl_mb(KLXR),BRT,NRDIM,NORB,
     &            NORB,NRDIM,NORB,NORB,FACTORC,FACTORAB,1)
*. X(T)I B in MAT4
      CALL MATML7(dbl_mb(KLMAT4),dbl_mb(KLXI),BRT,NRDIM,NORB,
     &            NORB,NRDIM,NORB,NORB,FACTORC,FACTORAB,1)
*. X(T) B X* = (X(T)(R) B X(R) + X(T)(I) B X(I)) 
*            +i(X(T)(I) B X(R) - X(T)(R) B X(I))
* (X(T)R B) XR
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(dbl_mb(KLBRBR),dbl_mb(KLMAT3),dbl_mb(KLXR),
     &            NRDIM,NRDIM,
     &            NRDIM,NORB,NORB,NRDIM,FACTORC,FACTORAB,0)
      FACTORC =  1.0D0
      FACTORAB = 1.0D0
* (X(T)I B) XI
      CALL MATML7(dbl_mb(KLBRBR),dbl_mb(KLMAT4),dbl_mb(KLXI),
     &            NRDIM,NRDIM,
     &            NRDIM,NORB,NORB,NRDIM,FACTORC,FACTORAB,0)
*
      FACTORC =  0.0D0
      FACTORAB = 1.0D0
*. (X(T)I B) XR 
      CALL MATML7(dbl_mb(KLBRBI),dbl_mb(KLMAT4),dbl_mb(KLXR),
     &            NRDIM,NRDIM,
     &            NRDIM,NORB,NORB,NRDIM,FACTORC,FACTORAB,0)
      FACTORC =  1.0D0
      FACTORAB =-1.0D0
*. (X(T)R B) XI
      CALL MATML7(dbl_mb(KLBRBI),dbl_mb(KLMAT3),dbl_mb(KLXI),
     &            NRDIM,NRDIM,
     &            NRDIM,NORB,NORB,NRDIM,FACTORC,FACTORAB,0)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Brillioun matrix in the X-basis '
        WRITE(6,*) ' ==================================== '
        CALL WRTCMAT(dbl_mb(KLBRBR),dbl_mb(KLBRBI),NRDIM,NRDIM)
      END IF
*
** Obtain G(i,j)
*
C          GPQ(E,NDIM,GR,GI)
      CALL GPQ(dbl_mb(KLEI),NRDIM,dbl_mb(KLGR),dbl_mb(KLGI))
*
** Obtain GB(i,j) = G(I,J)*BRB(I,J)
C          CVVTOV(VRIN1,VIIN1,VRIN2,VIIN2,VROUT,VIOUT,NDIM)
      CALL CVVTOV(dbl_mb(KLBRBR),dbl_mb(KLBRBI),dbl_mb(KLGR),
     &            dbl_mb(KLGI),
     &            dbl_mb(KLGBRBR),dbl_mb(KLGBRBI),NRDIM**2)
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' The element product of G and BRB '
        WRITE(6,*) '=========================== ===== '
        CALL WRTCMAT(dbl_mb(KLGBRBR),dbl_mb(KLGBRBI),NRDIM,NRDIM)
      END IF
*
*. Back transform GB to obtain contribution to gradient in actual basis
*
* X* GB X(T)
      CALL TRPMT3(dbl_mb(KLXR),NORB,NRDIM,dbl_mb(KLMAT1))
      CALL COPVEC(dbl_mb(KLMAT1),dbl_mb(KLXR),NORB*NRDIM)
      CALL TRPMT3(dbl_mb(KLXI),NORB,NRDIM,dbl_mb(KLMAT1))
      CALL COPVEC(dbl_mb(KLMAT1),dbl_mb(KLXI),NORB*NRDIM)
C     CTRN_MAT(XAXR,XAXI,XR,XI,AR,AI,NXR,NXC,NAR,NAC,SCR)
      CALL CTRN_MAT(dbl_mb(KLMAT1),dbl_mb(KLMAT2),dbl_mb(KLXR),
     &              dbl_mb(KLXI),
     &              dbl_mb(KLGBRBR),dbl_mb(KLGBRBI),NRDIM,NORB,
     &              NRDIM,NRDIM,dbl_mb(KLMAT3))
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GB matrix transformed to original basis '
        CALL WRTCMAT(dbl_mb(KLMAT1),dbl_mb(KLMAT2),NORB,NORB)
      END IF
*. testzone: Back transform B(Bar)
      IF(I_DO_TEST.EQ.1) THEN
        WRITE(6,*) ' B(Bar) will be back-transformed to B(tilde)'
        CALL CTRN_MAT(XJEP1,XJEP2,dbl_mb(KLXR),dbl_mb(KLXI),
     &                dbl_mb(KLBRBR),dbl_mb(KLBRBI),NRDIM,NORB,
     &                NRDIM,NRDIM,WORK(KLMAT3))
        WRITE(6,*) ' B(tilde) from back-transformation'
        CALL WRTCMAT(XJEP1,XJEP2,NORB,NORB)
      END IF
*. End of test-zone
       
*. And update E1
      CALL VECSUM(E1,E1,dbl_mb(KLMAT1),ONE,ONE,NORB*NOCC)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Linear transformed gradient as E1(R,S) '
       CALL WRTMAT(E1,NORB,NOCC,NORB,NOCC)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'LINGRA')
*
      RETURN
      END
      SUBROUTINE CMATML(CR,CI,AR,AI,BR,BI,
     &                  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL,
     &                  ITP,ICC,SCR)
 
*
* Multiply complex matrix (AR+i*AI) with complex matrix BR+i*BI
* to give CR + i*CI
*
* IF ITP=1 then A is transposed before multiplication
* IF ITP=2 then B is transposed before multiplication
*
* IF ICC=1 then A is complex conjugated before multiplication
* IF ICC=2 then B is complex conjugated before multiplication
*
* Jeppe Olsen , June 1989
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      DIMENSION CR(NCROW,NCCOL),CI(NCROW,NCCOL)
      DIMENSION AR(NAROW,NACOL),AI(NAROW,NACOL)
      DIMENSION BR(NBROW,NBCOL),BI(NBROW,NBCOL)
      DIMENSION SCR(NCROW,NCCOL)
*
      IF(ICC .LT. 0 .OR. ICC .GT. 2 ) THEN
        WRITE(6,*) ' Wrong input to CMATML, ICC = ',ICC
        STOP ' CMATML, ICC out of range '
      END IF
      IF(ITP .LT. 0 .OR. ICC .GT. 2 ) THEN
        WRITE(6,*) ' Wrong input to CMATML, ITP = ',ITP
        STOP ' CMATML, ITP out of range '
      END IF
* CR = AR*BR - SIGN1*SIGN2*AI*BI
      CALL MATML4(CR,AR,BR,NCROW,NCCOL,NAROW,NACOL,
     &            NBROW,NBCOL,ITP)
      CALL MATML4(SCR,AI,BI,NCROW,NCCOL,NAROW,NACOL,
     &            NBROW,NBCOL,ITP)
       IF(ICC.EQ.0) THEN
         SIGN = -1.0D0
       ELSE
         SIGN = 1.0D0
       END IF
       CALL VECSUM(CR,CR,SCR,1.0D0,SIGN,NCROW*NCCOL)
* CI = SIGN2*AR*BI + SIGN1*AI*BR
      CALL MATML4(CI,AR,BI,NCROW,NCCOL,NAROW,NACOL,
     &            NBROW,NBCOL,ITP)
      CALL MATML4(SCR,AI,BR,NCROW,NCCOL,NAROW,NACOL,
     &            NBROW,NBCOL,ITP)
      IF(ICC .EQ. 0 ) THEN
        SIGN1 =  1.0D0
        SIGN2 =  1.0D0
      ELSE IF(ICC .EQ. 1 ) THEN
        SIGN1 =  1.0D0
        SIGN2 = -1.0D0
      ELSE IF(ICC .EQ. 2 ) THEN
        SIGN1 = -1.0D0
        SIGN2 =  1.0D0
      END IF
*
      CALL VECSUM(CI,CI,SCR,SIGN1,SIGN2,NCROW*NCCOL)
*
      NTEST = 0
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Output from CMATML '
        WRITE(6,*) '===================='
        write(6,*) ' real and imaginary part of matrix A '
        CALL WRTMAT(AR,NAROW,NACOL,NAROW,NACOL)
        CALL WRTMAT(AI,NAROW,NACOL,NAROW,NACOL)
        write(6,*) ' real and imaginary part of matrix B '
        CALL WRTMAT(BR,NBROW,NBCOL,NBROW,NBCOL)
        CALL WRTMAT(BI,NBROW,NBCOL,NBROW,NBCOL)
        WRITE(6,*) ' Real part of product matrix : '
        CALL WRTMAT(CR,NCROW,NCCOL,NCROW,NCCOL)
        WRITE(6,*) ' Imaginary  part of product matrix : '
        CALL WRTMAT(CI,NCROW,NCCOL,NCROW,NCCOL)
      END IF
*
      RETURN
      END
      SUBROUTINE LINSE3(ALPHA,DIRECT,NVAR,MAXIT,ITNUM,X,E1,F,IFLAG,
     &                  IQUACI)
*
*
* Adapted for LUCIA, Aug. 2011
*
* Linesearch , version of July 92     
* uses derivatives as well as function values
*
* Possibility of dividing direction in groups have been removed
*
* Jeppe Olsen
*
*=======
* Input
*=======
*   DIRECT     : Search direction
*   X          : initial point
* IFLAG .NE. 0 : E1 and E0 contains already initial gradient anf value
*   ALPHA      : initial step size
*   MAXIT      : max. number of energy and gradient calculations
*
*   ISEPOP     : Differs from zero if individual optimization of 
*                different groups of parameters
*   NGRP       : Number of groups
*   IGRPO      : Offset for each group
*   IGRPN      : Number of elements in each group
*========
* Output
*========
*  E0         : Function value at final point
*  E1         : gradient at final point
*  X          : final parameter set
*  ALPHA      : final step size
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DIRECT(NVAR),X(NVAR),E1(NVAR)
      LOGICAL INRANG
      REAL * 8  INPROD
*
**.    Initialize
*
      NTEST = 10
      IF(NTEST.GE.2)WRITE( 6,*)  ' INFO from linesearch ( LINSE3 ) '
      ITNUM = 0
      IF ( MAXIT .LE. 0 ) MAXIT = 3
      MAXIT2 = 12
      MINIT = 1
* ?
      ALPINI = 1.0D0
      ALPHAN = 1.0D0
      FPMAX  = 1.0D0
*
      INRANG = .FALSE.
      ALPMN1 = 1.0D-4
      ALPMN2 = 1.0D-11
      ALPMIN = 0.0D0
      ALPLOW = 0.0D0
      ALPMAX = 0.0D0
      DIV = 3.0D0
      RHO = 0.25D0
*. Convergence criterion  for function value
      DELMIN = 1.0D-12
*. Numbert of sweeps over groups
*. Gradient and function value in initial point
      IF( IFLAG .EQ. 0 ) THEN
        CALL GRADIE(X,E1,F,1)
        IF(IQUACI .NE. 0 ) THEN
          CALL NEWCI
          CALL GRADIE(X,E1,F,1)
        END IF
      END IF
*
          KGRPN = NVAR
      ALPHAP = 0.0D0
      ALPHA = 1.0D0
      FINI=F
      FLOW = F
      FMIN = F
      FP = INPROD(DIRECT,E1,NVAR)
C?    write(6,*) ' FP ', FP
      IF (FP.GT. 0.0D0) THEN
         WRITE(6,*) ' ***LINSEA , WARNING :',
     &     ' *** DIRECTION REVERSED '
         CALL SCALVE(DIRECT,-1.0D0,NVAR)
         FP = - FP
      END IF
      FPLOW = FP
      FPINI = FP
      FPMIN = FP
      ALPMIN = 0.0D0
      ALPLOW = 0.0D0
      ALPMAX = 0.0D0
*
      VECNRM =SQRT ( INPROD(DIRECT,DIRECT,NVAR) )
 
      IF(ALPHA.EQ.0.0D0) ALPHA = 1.0D0
      IF(ALPHA .LE. ALPMN1 ) ALPHA = 0.05
      IF (NTEST .GE. 30 ) THEN
         write(6,*) ' info from linsea '
         write(6,*) ' ================='
         WRITE(6,*) ' INITIAL GRADIENT'
         CALL WRTMAT(E1,1,NVAR,1,NVAR)
         WRITE(6,*) ' INITIAL VECTOR '
         CALL WRTMAT(X ,1,NVAR,1,NVAR)
         WRITE(6,*) ' INITIAL DIRECTION '
         CALL WRTMAT(DIRECT ,1,NVAR,1,NVAR)
      END IF
*
**.  Loop over iterations
*
      ITNUM = 0
100   CONTINUE
       ITNUM = ITNUM + 1
*   ========================================================
*.. Energy and gradient corresponding to new iteration point
*   ========================================================
       ALPHAD = ALPHA - ALPHAP
C?     IF(NTEST.GE.2)
C?   & write(6,*) ' Iteration ALPHA ALPHAD ',ITNUM,ALPHA,ALPHAD
*.     X(*) + ALPHA*VECTOR(*)
       CALL VECSUM(X,X,DIRECT,1.0D0,ALPHAD,NVAR)
*.     New gradient and function value
       CALL GRADIE(X,E1,F,1)
       IF(IQUACI .NE. 0 ) THEN
         CALL NEWCI
         CALL GRADIE(X,E1,F,1)
       END IF
       IF( NTEST .GE. 1 ) THEN
         WRITE(6,'(1X,A,I4,2E15.7)') '  ** LINSE3: It alpha F    ',
     &            ITNUM,ALPHA,F
C        WRITE(6,*) 'PREDICTED AND ACTUAL CHANGE : ',
C    &                (FPINI*ALPHA),(F-FINI)
       END IF
*=========================
*.   Check for convergence
*=========================
       FP = INPROD(DIRECT,E1,NVAR)
C?     write(6,*) ' FP ' , FP
       IF(( ABS(FP).LT.RHO*ABS(FPINI).AND.
     &    (FINI-F).GT.(-RHO*FPINI*ALPHA  )) .OR.
     &    ABS(F-FINI).LT.DELMIN) THEN
*. Convergence
         INRANG = .TRUE.
       ELSE
*. no convergence so
*===============
*.     New alpha
*===============
       IF(F.LT.FLOW) THEN
*. New point is the lowest obtained
         IF(ALPHA.GT.ALPLOW) THEN
*(This combination implies that FPLOW.LT.0
           IF(FPLOW.LT.FP) THEN
*=================================
*. F<FLOW,FP>FPLOW,ALPHA>ALPLOW quadratic interpolation can be used
*=================================
             ALPHAN = (ALPHA-ALPLOW)*FP/(FPLOW-FP) + ALPHA
             IF(FP.GT.0.0D0) THEN
*. right bound is ALPHA
                ALPMAX = ALPHA
                FMAX   = F
                FPMAX  = FP
             ELSE IF(FP.LT.0.0D0) THEN
*. Minimum is to the left of alpha,check new alpha
                IF(ALPHAN.GT.ALPMAX.AND.ALPMAX.NE.0) THEN
*. use quadratic interpolation using current point
*. and end point instead
                  ALPHAN = (ALPHA-ALPMAX)*FP/(FPMAX-FP) + ALPHA
                ELSE IF (ALPMAX.EQ.0.0D0) THEN
                  IF(ALPHAN.GT.3*ALPHA) ALPHAN = 3*ALPHA
                END IF
                ALPMIN = ALPHA
                FMIN   = F
                FPMIN  = FP
             END IF
           ELSE IF(FPLOW.GE.FP) THEN
*=================================
*. F<FLOW,FP<FPLOW,ALPHA>ALPLOW  minimum is to the left of alpha
*=================================
             IF(ALPMAX.NE.0) THEN
*. bisection ( quad interpol could be used )
                ALPHAN = (ALPHA+ALPMAX)/2
             ELSE
                ALPHAN = 2* ALPHA
             END IF
             ALPMIN = ALPHA
             FMIN   = F
             FPMIN  = FP
           END IF
         ELSE IF ( ALPHA.LT.ALPLOW) THEN
*(Implies that FPLOW.GT.0
           IF(FP.LT.0.0D0) THEN
*=================================
*. F<FLOW,FP<0.0D0,ALPHA<ALPLOW quadratic interpolation can be used
*=================================
*. minimum is bounded by alpha and alplow
              ALPHAN = (ALPHA-ALPLOW)*FP/(FPLOW-FP) + ALPHA
              ALPMIN = ALPHA
              FMIN   = F
              FPMIN  = FP
           ELSE IF( FP.GT.0.0D0) THEN
*=================================
*. F<FLOW,FP>0.0D0,ALPHA<ALPLOW quadratic interpolation can be used
*=================================
*. Minimum as bounded by alpmin and alpha
             ALPHAN = (ALPHA-ALPMIN)*FP/(FPMIN-FP) + ALPHA
             ALPMAX = ALPHA
             FMAX = F
             FPMAX = FP
           END IF
         END IF
         ALPLOW = ALPHA
         FPLOW = FP
         FLOW =  F
       ELSE IF( F.GT.FLOW) THEN
         IF(ALPHA.GT.ALPLOW) THEN
           IF(FP.GT.0.0D0) THEN
*=================================
*. F>FLOW,FP>0 ,ALPHA>ALPLOW  quad interpolation
*================================
             ALPHAN = (ALPHA-ALPLOW)*FP/(FPLOW-FP) + ALPHA
           ELSE
*=================================
*. F>FLOW,FP0 ,ALPHA>ALPLOW  bisection
*=================================
             ALPHAN = (ALPHA+ALPLOW)/2
           END IF
           ALPMAX = ALPHA
           FMAX = F
           FPMAX = FP
         ELSE IF(ALPHA.LT.ALPLOW) THEN
*. a mess , a minimum has been overlooked,rectify this
           IF(FP.GT.0.0D0) THEN
*=================================
*. F>FLOW,FP>0 ,ALPHA<ALPLOW  quad interpolation using 0 and alpha
*=================================
             ALPHAN = (ALPHA-ALPINI)*FP/(FPINI-FP) + ALPHA
           ELSE
*=================================
*. F>FLOW,FP<0 ,ALPHA<ALPLOW  bisection
*=================================
             ALPHAN = ALPHA/2
           END IF
           ALPMIN = ALPINI
           FMIN =  FINI
           FPMIN =  FPINI
           ALPMAX = ALPHA
           FMAX =   F
           FPMAX = FP
         END IF
       END IF
       ALPHAP = ALPHA
       ALPHA  = ALPHAN
      IF(ITNUM.LT.MAXIT)  GOTO 100
      IF(F.GT.FINI .AND. ITNUM.LT.MAXIT2) GOTO 100
      IF(ITNUM .LT. MINIT ) GOTO 100
       END IF
*
*
      RETURN
      END
      SUBROUTINE EIGGMTN_LUCAS(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,
     &                   Z,W,SCR)
*
* Outer routine for calculating eigenvectors and eigenvalues
* of a general real matrix
*
* Version employing EISPACK path RG
*
* Current implementation is rather wastefull with respect to
* memory but at allows one to work with real arithmetic
* outside this routine
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8 INPROD
      DIMENSION AMAT(NDIM,NDIM),SCR(*)
      DIMENSION ARVAL(NDIM),AIVAL(NDIM)
      DIMENSION ARVEC(NDIM,NDIM),AIVEC(NDIM,NDIM)
      DIMENSION Z(NDIM,NDIM),W(NDIM)
*
* Diagonalize
*
      NSCR = 2*NDIM
      CALL RG(NDIM,NDIM,AMAT,ARVAL,AIVAL,1,Z,SCR(1),SCR(1+NDIM),IERR)
      IF( IERR.NE.0) THEN
        WRITE(6,*) ' Problem in EIGGMTN, no convergence '
        WRITE(6,*) ' I have to stop '
        STOP ' No convergence in EIGGMTN '
      END IF
*
* Extract real and imaginary parts according to Eispack manual p.89
*
      DO 150 K = 1, NDIM
*
        IF(AIVAL(K).NE.0.0D0) GOTO 110
        CALL COPVEC(Z(1,K),ARVEC(1,K),NDIM)
        CALL SETVEC(AIVEC(1,K),0.0D0,NDIM)
        GOTO 150
*
  110   CONTINUE
        IF(AIVAL(K).LT.0.0D0) GOTO 130
        CALL COPVEC(Z(1,K),ARVEC(1,K),NDIM)
        CALL COPVEC(Z(1,K+1),AIVEC(1,K),NDIM)
        GOTO 150
*
  130   CONTINUE
        CALL COPVEC(ARVEC(1,K-1),ARVEC(1,K),NDIM)
        CALL VECSUM(AIVEC(1,K),AIVEC(1,K),AIVEC(1,K-1),
     &              0.0D0,-1.0D0,NDIM)
*
  150 CONTINUE
 
 
*
* explicit orthogonalization of eigenvectors with
* (degenerate eigenvalues are not orthogonalized by DGEEV)
*
      DO 200 IVEC = 1, NDIM
         RNORM = INPROD(ARVEC(1,IVEC),ARVEC(1,IVEC),NDIM)
     &         + INPROD(AIVEC(1,IVEC),AIVEC(1,IVEC),NDIM)
         FACTOR = 1.0d0/SQRT(RNORM)
         CALL SCALVE(ARVEC(1,IVEC),FACTOR,NDIM)
         CALL SCALVE(AIVEC(1,IVEC),FACTOR,NDIM)
         DO 190 JVEC = IVEC+1,NDIM
* orthogonalize jvec to ivec
           OVLAPR = INPROD(ARVEC(1,IVEC),ARVEC(1,JVEC),NDIM)
     &            + INPROD(AIVEC(1,JVEC),AIVEC(1,IVEC),NDIM)
           OVLAPI = INPROD(ARVEC(1,IVEC),AIVEC(1,JVEC),NDIM)
     &            - INPROD(AIVEC(1,IVEC),ARVEC(1,JVEC),NDIM)
           CALL VECSUM(ARVEC(1,JVEC),ARVEC(1,JVEC),ARVEC(1,IVEC),
     &                 1.0D0,-OVLAPR,NDIM )
           CALL VECSUM(AIVEC(1,JVEC),AIVEC(1,JVEC),AIVEC(1,IVEC),
     &                 1.0D0,-OVLAPR,NDIM )
           CALL VECSUM(ARVEC(1,JVEC),ARVEC(1,JVEC),AIVEC(1,IVEC),
     &                 1.0D0,OVLAPI,NDIM )
           CALL VECSUM(AIVEC(1,JVEC),AIVEC(1,JVEC),ARVEC(1,IVEC),
     &                 1.0D0,-OVLAPI,NDIM )
  190    CONTINUE
  200 CONTINUE
 
*
* Normalize eigenvectors
*
      DO 300 L = 1, NDIM
        XNORM = INPROD(ARVEC(1,L),ARVEC(1,L),NDIM)
     &        + INPROD(AIVEC(1,L),AIVEC(1,L),NDIM)
        FACTOR = 1.0D0/SQRT(XNORM)
        CALL SCALVE(ARVEC(1,L),FACTOR,NDIM)
        CALL SCALVE(AIVEC(1,L),FACTOR,NDIM)
  300 CONTINUE
      NTEST = 0
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Output from EIGGMT '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Real and imaginary parts of eigenvalues '
        CALL WRTMAT(ARVAL,1,NDIM,1,NDIM)
        CALL WRTMAT(AIVAL,1,NDIM,1,NDIM)
        WRITE(6,*) ' real part of eigenvectors '
        CALL WRTMAT(ARVEC,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' imaginary part of eigenvectors '
        CALL WRTMAT(AIVEC,NDIM,NDIM,NDIM,NDIM)
      END IF
*
* Test : check orthonormality
C     kl1 = 1
C     kl2 = 1 + ndim ** 2
C     kl3 = kl2 + ndim ** 2
C     call cmatml(scr(kl1),scr(kl2),arvec,aivec,arvec,aivec,
C    &            ndim,ndim,ndim,ndim,ndim,ndim,1,1,scr(kl3))
C
C      write(6,*) ' real and imaginary parts of u* u '
C      call wrtmat(scr(kl1),ndim,ndim,ndim,ndim)
C      call wrtmat(scr(kl2),ndim,ndim,ndim,ndim)
      RETURN
      END
      SUBROUTINE COMMAT(A,B,ACOMB,SCR,NDIM)
*
* ACOMB = (A,B)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),B(*),ACOMB(*),SCR(*)
*
      NTEST = 0
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' COMMAT speaking '
        WRITE(6,*) ' ================'
      END IF
*
      CALL MATML4(ACOMB,A,B,
     &            NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' AB = '
        CALL WRTMAT(ACOMB,NDIM,NDIM,NDIM,NDIM)
      END IF
      
*
      CALL MATML4(SCR,B,A,
     &            NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' BA = '
        CALL WRTMAT(SCR,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      CALL VECSUM(ACOMB,ACOMB,SCR,1.0D0,-1.0D0,NDIM ** 2 )
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' AB-BA = '
        CALL WRTMAT(ACOMB,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE COMMUP(C,A,B,NDIM,SCR,ISKIP0)
*
*  SIMPLE TEST ROUTINE !!!!!
*
*  C = Sum(N = 0,Infinity) (B Super) ** N A  /(N+1)!
*
*       With B Super A = (B,A)
*
* If iskip0 differs from zero, the first term (A) is skipped
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*),B(*),C(*),SCR(*)
      REAL*8 INPROD
*
      NTEST = 000
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Information from COMMUP '
        WRITE(6,*) ' ========================'
      END IF
      IF(NTEST.GE.100 ) THEN
        WRITE(6,*) ' commutator matrix '
        CALL WRTMAT(B,NDIM,NDIM,NDIM,NDIM)
        WRITE(6,*) ' Matrix to be commuted '
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
* Length of SCR should at least be 3*NDIM ** 2
*
      KLFREE = 1
*
      KLM1 = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      KLM2 = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      KLM3 = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      IF(ISKIP0 .EQ. 0 ) THEN
        CALL COPVEC(A,C,NDIM ** 2 )
      ELSE
        CALL SETVEC(C,0.0D0,NDIM ** 2 )
      END IF
      CALL COPVEC(A,SCR(KLM1),NDIM ** 2 )
      XNP1 = 1.0D0
      MAXN = 20
*
      DO 100 N = 1, MAXN
        XNP1 = XNP1 + 1.0D0
C       COMMAT(A,B,ACOMB,SCR,NDIM)
        CALL COMMAT(B,SCR(KLM1),SCR(KLM2),SCR(KLM3),NDIM)
        CALL COPVEC(SCR(KLM2),SCR(KLM1),NDIM ** 2 )
        CALL SCALVE(SCR(KLM1),1.0D0/XNP1,NDIM ** 2 )
        XNORM = SQRT(INPROD(SCR(KLM1),SCR(KLM1),NDIM**2))
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' C before VECSUM '
          CALL WRTMAT(C,NDIM,NDIM,NDIM,NDIM)
        END IF
        CALL VECSUM(C,C,SCR(KLM1),1.0D0,1.0D0, NDIM ** 2 )
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Output matrix after N = ', N
          CALL WRTMAT(C,NDIM,NDIM,NDIM,NDIM)
        END IF
  100 CONTINUE
*
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Norm of last added correction = ', XNORM
      END IF
      IF(NTEST.GE.100 ) THEN
        WRITE(6,*) ' Output matrix '
        CALL WRTMAT(C,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE EXPFSK(B,EXPMB,NDIM,SCR)
*
* Straight forward calculation of EXP(-B) for a general matrix b
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(*),EXPMB(*),SCR(*)
*
      KLFREE = 1
*
      KLM1  = KLFREE
      KLFREE = KLM1 + NDIM ** 2
*
      KLM2  = KLFREE
      KLFREE = KLM2 + NDIM ** 2
*
C          ADDDIA(A,FACTOR,NDIM,IPACK)
      CALL SETVEC(EXPMB,0.0D0,NDIM ** 2 )
      CALL ADDDIA(EXPMB,1.0D0,NDIM,0)
      CALL SETVEC(SCR(KLM1),0.0D0, NDIM ** 2)
      CALL ADDDIA(SCR(KLM1),1.0D0,NDIM,0)
*
      XFAC = 1.0D0
      XN = 0.0D0
      DO 100 N = 1, 20
        XN = XN + 1.0D0
        XFAC = XFAC * XN
          CALL MATML4(SCR(KLM2),SCR(KLM1),B,
     &         NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0)
          CALL COPVEC(SCR(KLM2),SCR(KLM1),NDIM ** 2 )
          CALL SCALVE(SCR(KLM1),-1.0D0/XN,NDIM ** 2 )
          CALL VECSUM(EXPMB,EXPMB,SCR(KLM1),1.0D0,1.0D0,NDIM ** 2 )
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output form EXPMB '
        WRITE(6,*) ' =========================== '
        CALL WRTMAT(B,NDIM,NDIM,NDIM,NDIM)
        CALL WRTMAT(EXPMB,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE CVVTOV(VRIN1,VIIN1,VRIN2,VIIN2,VROUT,VIOUT,NDIM)
*
* Vout(i) = Vin1(i)*Vin2(i)
*
* Vout,Vin1,Vin2 Complex vectors stored as real and imaginary parts
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VRIN1(NDIM),VIIN1(NDIM),VRIN2(NDIM),VIIN2(NDIM)
      DIMENSION VROUT(NDIM),VIOUT(NDIM)
*
      DO 100 I = 1, NDIM
        VROUT(I) = VRIN1(I)*VRIN2(I)-VIIN1(I)*VIIN2(I)
        VIOUT(I) = VRIN1(I)*VIIN2(I)+VIIN1(I)*VRIN2(I)
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Real and imaginary vector from CVVTOV'
        WRITE(6,*) ' ====================================='
        CALL WRTMAT(VROUT,1,NDIM,1,NDIM)
        CALL WRTMAT(VIOUT,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE PROJVC(VECIN,X,NDIM,NVEC,VECOUT,IPQ)
*
* Project a vector :
*
* IPQ = 1 :VECOUT = SUM(I) X(I) * (X(I)(T) * VECIN)
* IPQ = 2 :VECOUT = VECIN - SUM(I) X(I) * (X(I)(T) * VECIN)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECIN(NDIM),X(NDIM,NVEC),VECOUT(NDIM)
      REAL*8 INPROD
*
      CALL SETVEC(VECOUT,0.0D0,NDIM)
      DO 100 I = 1, NVEC
        OVLAP = INPROD(X(1,I),VECIN,NDIM)
        CALL VECSUM(VECOUT,VECOUT,X(1,I),1.0D0,OVLAP,NDIM)
  100 CONTINUE
*
      IF(IPQ.EQ.2) THEN
        CALL VECSUM(VECOUT,VECOUT,VECIN,-1.0D0,1.0D0,NDIM)
      END IF
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from PROJVC'
        WRITE(6,*) ' ===================================='
        CALL WRTMAT(VECIN,1,NDIM,1,NDIM)
        CALL WRTMAT(VECOUT,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE CFORA(C,A,NDIM,SCR)
*
* Obtain the matrix sum
*
* C(A) = Sum(n) (A) ** n /(N+1)!
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION C(NDIM,NDIM),A(NDIM,NDIM),SCR(*)
*
* SCR should at least be of length 4*NDIM**2 + 4*NDIM
*
* Jeppe Olsen , June 1989
*
* Analytical evaluation added,assuming that A is antisymmetric
* Jeppe Olsen, April 1990
*
      NTEST = 0
      IF( NTEST .GE. 5 ) THEN
        WRITE(6,*) ' Output from CFORA '
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' Input matrix : '
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
*. yet a matrix
      KLMAT2 = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      KLAR1 = KLFREE
      KLFREE = KLFREE + NDIM
*
      KLAR2 = KLFREE
      KLFREE = KLFREE + NDIM
*
      KLAR3 = KLFREE
      KLFREE = KLFREE + NDIM
*
** Obtain A ** 2 and diagonalize
*
      CALL MATML4(SCR(KLA2),A,A,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0)
      CALL TRIPAK(SCR(KLA2),SCR(KLMAT1),1,NDIM,NDIM)
      CALL EIGENL(SCR(KLMAT1),SCR(KLA2VC),NDIM,0,1)
      CALL COPDIA(SCR(KLMAT1),SCR(KLA2VL),NDIM,1)
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Eigenvalues of A squared '
        CALL WRTMAT(SCR(KLA2VL),NDIM,1,NDIM,1)
        WRITE(6,*) ' eigenvectors '
        CALL WRTMAT(SCR(KLA2VC),NDIM,NDIM,NDIM,NDIM)
      END IF
*
** 3 Obtain arrays sum(n) e ** n /(2n+1)! and e **n /(2n+2)!
*
* Max number of terms required
c     REAL FUNCTION FNDMNX*8(VECTOR,NDIM,MINMAX)
      XMAX = FNDMNX(SCR(KLA2VL),NDIM,2)
      XMIN = FNDMNX(SCR(KLA2VL),NDIM,1)
      XMAX = MAX(ABS(XMAX),ABS(XMIN))
      NTERM = 0
*. Threshold for switching between analytical and expansion evaluation
C     THRES = 1.0D-1
      THRES = 1.0D+2
      IF(ABS(XMAX).LT.THRES) THEN
        NTERM = 0
        TEST = 1.0D-15
        ELMNT = XMAX
        X2NP1 = 1.0D0
  230   CONTINUE
          NTERM = NTERM + 1
          X2NP1 = X2NP1 + 2.0D0
          ELMNT = ELMNT*XMAX/(X2NP1*(X2NP1-1))
        IF(XMAX .EQ. 0.0D0 ) GOTO 231
        IF(ELMNT/XMAX.GT.TEST) GOTO 230
  231   CONTINUE
        IF(NTEST.GE.5) write(6,*) ' XMAX NTERM ', XMAX,NTERM
      END IF
*
*. First term:   
*
      IF(ABS(XMAX).GE.THRES) THEN
        IF(NTEST.GE.5) 
     &  WRITE(6,*) ' analytical expansion used in CFORA '
*. Use analytical formulaes
        DO 311 N = 1,NDIM
          EPSIL = SCR(KLA2VL-1+N)
          IF(EPSIL.GE.(-1.0D-12)) THEN
            SCR(KLAR1-1+N) = 0.0D0
          ELSE
            SCR(KLAR1-1+N) = SIN(SQRT(-EPSIL))/SQRT(-EPSIL) - 1
          END IF
  311   CONTINUE
      ELSE
*. Series expansion
        IF(NTEST.GE.5)
     &  WRITE(6,*) ' series expansion used in CFORA '
        CALL SETVEC(SCR(KLAR1),0.0D0,NDIM)
        CALL SETVEC(SCR(KLAR3),1.0D0,NDIM)
        X2NP1 = 1.0D0
        DO 300 N = 1, NTERM
          X2NP1 = X2NP1 + 2.0D0
          FACTOR = 1.0D0/(X2NP1*(X2NP1-1))
          CALL VVTOV(SCR(KLA2VL),SCR(KLAR3),SCR(KLAR3),NDIM)
          CALL SCALVE(SCR(KLAR3),FACTOR,NDIM)
          CALL VECSUM(SCR(KLAR1),SCR(KLAR1),SCR(KLAR3),
     &                1.0D0,1.0D0,NDIM)
  300   CONTINUE
      END IF
*
*. Second term
*
      IF(ABS(XMAX).GE.THRES) THEN
*. Use analytical formulaes
        DO 333 N = 1, NDIM
          EPSIL = SCR(KLA2VL-1+N)
          IF(EPSIL.GE.(-1.0D-12)) THEN
            SCR(KLAR2-1+N)= 0.5D0
          ELSE
            SCR(KLAR2-1+N)= (COS(SQRT(-EPSIL))-1)/EPSIL
          END IF
  333   CONTINUE
      ELSE
*. Use expansion
        CALL SETVEC(SCR(KLAR2),0.5D0,NDIM)
        CALL SETVEC(SCR(KLAR3),0.5D0,NDIM)
        X2NP2 = 2.0D0
        DO 330 N = 1, NTERM
          X2NP2 = X2NP2 + 2.0D0
          FACTOR = 1.0D0/(X2NP2*(X2NP2-1))
          CALL VVTOV(SCR(KLA2VL),SCR(KLAR3),SCR(KLAR3),NDIM)
          CALL SCALVE(SCR(KLAR3),FACTOR,NDIM)
          CALL VECSUM(SCR(KLAR2),SCR(KLAR2),SCR(KLAR3),
     &                1.0D0,1.0D0,NDIM)
  330   CONTINUE
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The sin(epsil)/epsil -1 array '
        CALL WRTMAT(SCR(KLAR1),1,NDIM,1,NDIM)
        WRITE(6,*) ' The -(cos(epsil)-1)/epsil**2  array '
        CALL WRTMAT(SCR(KLAR2),1,NDIM,1,NDIM)
      END IF
*
** 4  C(A)
*
* A * U * Dia2 * U(t)
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(SCR(KLMAT1),SCR(KLA2VC),SCR(KLAR2),
     &            NDIM,SCR(KLMAT2))
      CALL MATML4(SCR(KLMAT2),A,SCR(KLMAT1),
     &            NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,0)
      IF(NTEST.GE.100) THEN
      write(6,*) 'a u dia2 u(t)'
      CALL WRTMAT(SCR(KLMAT2),NDIM,NDIM,NDIM,NDIM)
      END IF
      CALL COPVEC(SCR(KLMAT2),C,NDIM**2)
* U * Dia1 * U(T)
      CALL XDIAXT(SCR(KLMAT1),SCR(KLA2VC),SCR(KLAR1),
     &            NDIM,SCR(KLMAT2) )
*
      IF(NTEST.GE.100) THEN
      write(6,*) ' u dia1 u(t)'
      CALL WRTMAT(SCR(KLMAT1),NDIM,NDIM,NDIM,NDIM)
      END IF
*
CT    CALL VECSUM(C,SCR(KLMAT1),C,
CT   &            +1.0D0,-1.0D0,NDIM ** 2 )
      CALL VECSUM(C,SCR(KLMAT1),C,
     &            +1.0D0,+1.0D0,NDIM ** 2 )
*
      IF( NTEST .GE. 5 ) THEN
        WRITE(6,*) '  C(A)'
        WRITE(6,*) ' ======'
        CALL WRTMAT(C,NDIM,NDIM,NDIM,NDIM)
      END IF
*
* TEST : straight forward way of obtaining C
*
C     CALL COPVEC(A,SCR(KLMAT1),NDIM ** 2)
C     CALL SCALVE(SCR(KLMAT1),-0.5D0,NDIM ** 2 )
C     CALL COPVEC(SCR(KLMAT1),SCR(KLMAT2),NDIM ** 2 )
C     DO 1236 N = 2, 2*NTERM
C       XNP1 =  DBLE(N+1)
C       CALL MATML4(SCR(KLFREE),SCR(KLMAT2),A,NDIM,NDIM,
C    &              NDIM,NDIM,NDIM,NDIM,0)
C       CALL COPVEC(SCR(KLFREE),SCR(KLMAT2),NDIM ** 2 )
C       FACTOR = -1.0D0/XNP1
C       CALL SCALVE(SCR(KLMAT2),FACTOR,NDIM ** 2 )
C       CALL VECSUM(SCR(KLMAT1),SCR(KLMAT1),SCR(KLMAT2),
C    &              1.0D0,1.0D0,NDIM ** 2 )
C1236 CONTINUE
C     IF( NTEST .GE. 5 ) THEN
C       WRITE(6,*) ' C(A) obtained with  matrix multiplications'
C       CALL WRTMAT(SCR(KLMAT1),NDIM,NDIM,NDIM,NDIM)
C     END IF
      RETURN
      END
      SUBROUTINE LRMTV2(VECIN,ARED,VRED,NFULL,NRED,VECOUT,SCR,ITRNSP)
*
* Calculate product of low rank matrix and vector
*
* Vecout = Vred * Ared    * Vred (T) * Vecin  ( Itrnsp = 0 )
* Vecout = Vred * Ared(T) * Vred (T) * Vecin  ( Itrnsp = 1 )
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECIN(NFULL),VECOUT(NFULL)
      DIMENSION VRED(NFULL,NRED)
      DIMENSION ARED(NRED,NRED)
      DIMENSION SCR(*)
      REAL * 8 INPROD
* 1 : Vred(T) * Vecin
      DO 100 I = 1, NRED
        SCR(I) = INPROD(VRED(1,I),VECIN(1),NFULL)
  100 CONTINUE
* 2 : Ared*Vred(T)*Vecin
C    MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(ARED,SCR(1),SCR(1+NRED),NRED,NRED,ITRNSP)
* 3 : Vred*Ared*Vred(T)*Vecin
      CALL SETVEC(VECOUT,0.0D0,NFULL)
      DO 200 J = 1, NRED
        FACTOR = SCR(NRED+J)
        CALL VECSUM(VECOUT,VECOUT,VRED(1,J),1.0D0,
     &              FACTOR,NFULL)
  200 CONTINUE
*
      NTEST = 0
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Input and output vector from LRMTVC '
        CALL WRTMAT(VECIN,1,NFULL,1,NFULL)
        CALL WRTMAT(VECOUT,1,NFULL,1,NFULL)
      END IF
*
      RETURN
      END
      SUBROUTINE REDKAP(KAPPA,NOCC,NORB,NREDVC,REDVEC,REDKP)
*
* For a given symmetry, 
* a kappa matrix is given in form of KAPPA in expanded form
* (the lower half of kappa)
*
*========
* Purpose
*=========
*
* Exploit the low rank nature of Kappa to obtain a 2*NOCC
* Dimensional space so
*
* Kappa = V Kappa(red) V(T)
* Where
*       V : NPOINT X 2*NOCC  matrix giving orthonormal basis for
*           reduced space
*       Kappa(red): Kappa in subspace
*
*=========
* Input
*=========
*
* KAPPA : Independent parameters of kappa in expanded form,
*         Lower triangular matrix stored in NORB * NOCC matrix
* NOCC  : Number of occupied orbitals
* NORB : Number of orbitals
*
*=========
* output
*=========
*
* NREDVC : Dimension of subspace
* REDVEC : Expansion of basis vectors of subspace , V
* REDKP  : Reduced kappa matrix, NREDVC X NREDVC matrix ,kappa (red)
*
* Jeppe Olsen, Oct. 2011
*              Simplified version of REDKS2 from LUCAS
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL * 8 INPROD
      REAL * 8  KAPPA(NORB,NOCC),REDVEC(NORB,*)
      DIMENSION REDKP(*)
*. Length of scratch: NRDIM**2 : Not in use anymore
C     DIMENSION SCR(*)
*
      DIMENSION XJEP1(2000),XJEP2(2000)
*
*
      NTEST = 000
*
      ZERO = 0.0D0
      ONE = 1.0D0
*. Expected number of linear independent vectors 
      NRDIM_I = MIN(NORB,2*NOCC)
      IF(NTEST.GE.100) write(6,*) 
     &' Initial NRDIM in REDKAP ',NRDIM_I
*
* Threshold for linear dependency ( squared norm )
      XDEP = 1.0D-18
      NZERO = 0
*
* Scratch space
*
      KL1 = 1
      KLEND = KL1 + NRDIM_I**2
*
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' Input kappa matrix '
        CALL WRTMAT(KAPPA,NORB,NOCC,NORB,NOCC)
      END IF
*
*=============================================================
* 1 :        Basis vectors of reduced space                  *
*=============================================================
*
*
* The first new nocc vectors equals the old- the occupied orbitals 
*
      NREDVC = 0
      DO IVEC = 1, NOCC
        CALL SETVEC(REDVEC(1,IVEC),0.0D0,NORB)
        REDVEC(IVEC,IVEC) = ONE
        NREDVC = NREDVC + 1
      END DO
*
* Following NOCC vectors span Kappa ( except occ - occ ) block
*
      DO IVEC = 1, NRDIM_I - NOCC
        IVECEF = IVEC + NOCC
        CALL COPVEC(KAPPA(1,IVEC),REDVEC(1,IVECEF),NORB)
        CALL SETVEC(REDVEC(1,IVECEF),ZERO,NOCC)
        XINI = INPROD(REDVEC(1,IVECEF),REDVEC(1,IVECEF),NORB)
* Orthogonalize to  previous vectors
        DO JVEC = 1, IVEC-1
          JVECEF = JVEC+NOCC
          OVLAP = INPROD(REDVEC(1,IVECEF),REDVEC(1,JVECEF),NORB)
          FACTOR = - OVLAP
          CALL VECSUM(REDVEC(1,IVECEF),REDVEC(1,IVECEF),
     &               REDVEC(1,JVECEF),1.0D0,FACTOR,NORB)
        END DO
* Normalize or kill
        XNORM = INPROD(REDVEC(1,IVECEF),REDVEC(1,IVECEF),NORB)
        IF(XNORM .GT. XDEP*XINI) THEN
          FACTOR = 1.0D0/SQRT(XNORM)
          CALL SCALVE(REDVEC(1,IVECEF),FACTOR,NORB)
          NREDVC = NREDVC + 1
        ELSE
          CALL SETVEC(REDVEC(1,IVECEF),0.0D0,NORB)
          NZERO = NZERO + 1
        END IF
      END DO ! End of loop over IVEC
*
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information form REDKAP'
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' Number of zero vectors ', NZERO
        WRITE(6,*) ' Dimension of reduced space ', NREDVC
        WRITE(6,*) ' Input kappa matrix '
        CALL WRTMAT(KAPPA,NORB,NOCC,NORB,NOCC)
        WRITE(6,*) ' Basis for reduced kappa matrix '
        CALL WRTMAT(REDVEC,NORB,NREDVC,NORB,NREDVC)
      END IF
*
*
*=============================================================
* 2 :          Kappa reduced space                           *
*=============================================================
*
* Kappa matrix in reduced basis
*
      CALL SETVEC(REDKP,ZERO,NREDVC** 2 )
* Occupied - Occupied block
      DO JOCC = 1,NOCC
        DO IOCC =JOCC+1,NOCC
          REDKP((JOCC-1)*NREDVC+IOCC) = KAPPA(IOCC,JOCC)
          REDKP((IOCC-1)*NREDVC+JOCC) = -KAPPA(IOCC,JOCC)
        END DO
      END DO
* virtual - occupied block and occupied-virtual blocks
      DO JOCC = 1, NOCC
        DO IVIR = 1,NREDVC-NOCC
          REDKP((JOCC-1)*NREDVC+NOCC+IVIR) = 
     &    INPROD(KAPPA(1,JOCC),REDVEC(1,NOCC+IVIR),NORB)
          REDKP((IVIR+NOCC-1)*NREDVC+JOCC) = 
     &   -REDKP((JOCC-1)*NREDVC+NOCC+IVIR)
        END DO
      END DO
*
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Reduced Kappa matrix '
        WRITE(6,*) ' ========================='
        CALL WRTMAT(REDKP,NREDVC,NREDVC,NREDVC,NREDVC)
      END IF
*
*=============================================================
* 3 :          Preliminary tests ( to be removed )           *
*=============================================================
*
*
      I_DO_TESTS = 0
      IF(I_DO_TESTS.EQ.1) THEN
        WRITE(6,*) ' Tests of reduced kappa: '
*  reproduce Kappa vector from reduced space
* Kappa(red) * V(T)
        CALL MATML4(XJEP1,REDKP,REDVEC,NREDVC,NORB,NREDVC,NREDVC,
     &              NORB,NREDVC,2)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Kappa(ref) * V(T) '
          CALL WRTMAT(XJEP1,NREDVC,NORB,NREDVC,NORB)
        END IF
* V *  Kappa(red) * V(T)
        CALL MATML4(XJEP2,REDVEC,XJEP1,NORB,NORB,NORB,NREDVC,
     &              NREDVC,NORB,0)
        write(6,*) ' Reexpanded Kappa matrix '
        call wrtmat(Xjep2,NORB,NORB,NORB,NORB)
*
* V(T) V (should be 1)
*
        CALL MATML4(XJEP2,REDVEC,REDVEC,
     &  NREDVC,NREDVC,NORB,NREDVC,NORB,NREDVC,1)
        WRITE(6,*) ' V(T) V (should be 1 ) '
        CALL WRTMAT(XJEP2,NREDVC,NREDVC,NREDVC,NREDVC)
      END IF
      RETURN
      END
      SUBROUTINE GPQ(E,NDIM,GR,GI)
*
* GENERATE THE COMPLEX ARRAY
*
* G(P,Q) = SUM(N=1) i**N ( E(P) - E(Q))**N/(N+1)!
*        = i * (1-COS(E(P)-E(Q)))/(E(P)-E(Q))
*        + SIN(E(P)-E(Q))/(E(P)-E(Q)) -1
*
* Jeppe Olsen,Summer of '89 ,
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER P,Q
      DIMENSION E(NDIM),GR(NDIM,NDIM),GI(NDIM,NDIM)
*
      NMAX = 0
      TEST1 = 1.0D-15
*
      ZERO = 0.0D0
      CALL SETVEC(GR,ZERO,(NDIM)**2)
      CALL SETVEC(GI,ZERO,(NDIM)**2)
*
      DO 100 P = 1, NDIM
      DO 100 Q = 1, NDIM
       EPQ = E(P) - E(Q)
       IF(ABS(EPQ).GE.1.0D-3) THEN
* Use analytical formulaes
         GR(P,Q) = SIN(EPQ)/EPQ - 1.0D0
         GI(P,Q) = (1.0D0-COS(EPQ))/EPQ
       ELSE
* use series expansion
        EPQN = 1.0D0
        XNP1F = 1.0D0
        NPL1 = 1
*       LOOP OVER 4 * N
   99   CONTINUE
*
          EPQN = EPQN * EPQ
          NPL1 = NPL1 + 1
          XNP1F = XNP1F * DBLE(NPL1)
          GI(P,Q) = GI(P,Q) + EPQN/XNP1F
*
          EPQN = EPQN * EPQ
          NPL1 = NPL1 + 1
          XNP1F = XNP1F * DBLE(NPL1)
          GR(P,Q) = GR(P,Q) - EPQN/XNP1F
*
          EPQN = EPQN * EPQ
          NPL1 = NPL1+1
          XNP1F = XNP1F * DBLE(NPL1)
          GI(P,Q) = GI(P,Q) - EPQN/XNP1F
*
          EPQN = EPQN * EPQ
          NPL1= NPL1 + 1
          XNP1F = XNP1F * DBLE(NPL1)
          GR(P,Q) = GR(P,Q) + EPQN/XNP1F
*
          IF (NPL1 .GT.NMAX) NMAX = NPL1
C         IF(ABS(EPQN/XNP1F) .GT. TEST1 ) GOTO 99
        IF(EPQ .EQ. 0.0D0 ) GOTO 100
        IF(ABS(EPQN/(XNP1F*EPQ)).GT.TEST1 ) GOTO 99
       END IF
  100 CONTINUE
*
C?    WRITE(6,*) ' MAXIMAL N NEEDED TO GENERATE G ..', (NMAX-1)
*
      NTEST = 00
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Output from routine GPQ '
        WRITE(6,*) ' ========================'
        WRITE(6,*) ' Input vector of eigenvalues '
        CALL WRTMAT(E,1,NDIM,1,NDIM)
        WRITE(6,*) ' Real and imaginary parts of G obtained in GPQ '
        CALL WRTMAT(GR,NDIM,NDIM,NDIM,NDIM)
        CALL WRTMAT(GI,NDIM,NDIM,NDIM,NDIM)
      END IF
 
      RETURN
      END
      SUBROUTINE QUASI( MAXIT,X,E1,E2,SCR,NVAR,IBARR,IDOE2,IQUACI,
     &             ISEPOP,NGRP,IGRPO,IGRPN,                   
     &                  ICONV)
C
C OUTER ROUTINE FOR QUASI NEWTON METHODS
C
C
C.. INPUT : X : SPACE FOR SOLUTION VECTOR
C           E1 : SPACE FOR GRADIENT
C           E2 : SPACE FOR APPROXIMATION TO INVERSE HESSIAN
C              : CAN BE DIAGONAL OR GENERAL . OBTAINED ON HESINV.
C           SCR : SCRATCH SPACE, MINIMUM LENGTH :
C                 6*NVAR +4 * MAXIT
C           NVAR : TOTAL NUMBER OF VARIABLES
C           IB(I) :  First column in row I of Hessian approximation
C           IDOE2 : RECALCULATE HESSIAN APPROXIMATION
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL DISCH
      DIMENSION SCR(*),E1(*),E2(*),X(*),IBARR(*)
C
C
C *** ALLOCATION OF SCRATCH MEMORY FOR QUASIN
C
C..1 THREE SCRATCH VECTORS
      KVEC1 = 1
      KVEC2 = KVEC1 + NVAR
      KVEC3 = KVEC2 + NVAR
      KVEC4 = KVEC3 + NVAR
C..2  MAXIT TWO BY TWO MATRICES
      KAMAT = KVEC4 + NVAR
C..3  2 VECTORS FOR ( DISCH IN USE )
      KAVEC = KAMAT + 4*MAXIT
C..4
      KFREE = KAVEC + 2*NVAR
C?    WRITE(6,*) ' SCRATCH SPACE USED IN QUASIN ', KFREE - 1
C
C** TRANSFER CONTROL TO QUASI NEWTON PART
C
      CALL SETVEC(X,0.0D0,NVAR)
      DISCH = .TRUE.
      NRESET = 21
*. Convergence threshold for norm of gradient 
      CNVNRM = 1.0D-8
      CALL QUASIN(DISCH ,X,E0,E1,SCR(KVEC1),SCR(KVEC2),
     &            SCR(KVEC3),MAXIT,NRESET,CNVNRM,37,
     &            SCR(KAVEC),SCR(KAMAT), NVAR ,IBARR,E2,IDOE2,
     &            SCR(KVEC4),IQUACI,
     &            ISEPOP,NGRP,IGRPO,IGRPN,ICONV)                   
C
      RETURN
      END
      SUBROUTINE QUASIN( DISCH,X,E0,E1,VEC1,VEC2,VEC3,MAXIT,NRESET,
     &                        CNVNRM,LUHFIL,AVEC,AMAT,NVAR,IBARR,
     &                        E2 ,IDOE2,VEC4,IQUACI,
     &                        ISEPOP,NGRP,IGRPO,IGRPN,
     &                        ICONV)
C     Master routine for optimization with quasi - newton methods
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL * 8  INPROD
      LOGICAL CONVRG,DISCH
      DIMENSION X(NVAR),E1(NVAR),VEC1(NVAR),VEC2(NVAR),VEC3(NVAR)
      DIMENSION AMAT(4*MAXIT)
      DIMENSION AVEC(*),VEC4(*)
      DIMENSION IBARR(*)
C
      DIMENSION E2(*)
C E2 IS SPACE FOR HESSIAN APPROXIMATION
C
C IF DISCH IS TRUE VECTORS DEFINING HESSIAN APPROXIMATION ARE WRITTEN
C ON DISC ON FILE  LUHFIL
C
      IF ( DISCH ) REWIND LUHFIL
*
*     INPUT: X   INITIAL GUESS
*     Convergence defined as obtained when
*     norm of gradient .le. CNVNRM
*
      IUPDAT=2
      IPRT = 0
*
      WRITE (*,*)
      IF(IPRT .GT. 0 ) THEN
      IF(IUPDAT.EQ.1) WRITE(6,*) ' Quasi - Newton method (DFP) in use'
      IF(IUPDAT.EQ.2) WRITE(6,*) ' Quasi - Newton method (BFGS) in use'
      IF(IUPDAT.EQ.3) WRITE(6,*) ' Initial approximation not updated '
      END IF
*
**   Initialize
*
      IINV=1
      ALPHA = 1.0D0
      CONVRG = .FALSE.
      E0  =0.0
      IRESET = 0
      ICONV = 0
*
** Max number of iterations in linesearch
*
* dfp needs fairly accurate linesearch
      IF ( IUPDAT .EQ. 1 ) MAXITL = 3
* bfgs is not so sensitive, experience shows 3 anyhow
      IF ( IUPDAT .EQ. 2 )  MAXITL = 3
* scaled steepest descent depends on linesearch so
      IF( IUPDAT.EQ.3) MAXITL = 4
*
**   Loop over iterations
*
      ITNUM = 0
  100 CONTINUE
         ITNUM = ITNUM +1
         IRESET = IRESET + 1
         IF(ITNUM.EQ.1 ) THEN
*           Initial hessian approximation
            IHSAPR= 3
            IF ( IDOE2 .NE. 0 ) CALL HESINV(E2,IBARR)
            NMAT = 0
            IFLAG = 1
*           GRADIENT AND FUNCTION VALUE AT INITIAL POINT
            CALL GRADIE(X,E1,E0,1)
C!          CALL GFUSK(E1,X)
C                GFUSK (E1FD,STEP)
            CALL COPVEC(X,VEC2,NVAR)
            CALL COPVEC(E1,VEC3,NVAR)
         ELSE
**          Perform update of hessian approximation
            IF(IRESET.EQ.NRESET.AND.NRESET.NE.0) THEN
*              Hessian reset
                WRITE(*,*) ' HESSIAN RESET '
                IRESET = 0
                NMAT = 0
C SPECIFIC PROGRAMMING
                CALL COPVEC(X,VEC3,NVAR)
                CALL NEWMOS(VEC3,1)
                CALL SETVEC(X,0.0D0,NVAR)
                CALL HESINV(E2,IBARR)
C END OF SPECIFIC PROGRAMMING
                CALL COPVEC(X,VEC2,NVAR)
                CALL COPVEC(E1,VEC3,NVAR)
             ELSE
                IF (IUPDAT.NE.3) THEN
                  CALL HESUPV (E2,AMAT,AVEC,
     &                 X,E1,VEC2,
     &                 VEC3,NVAR,IUPDAT,IINV,VEC1,NMAT,
     &                 LUHFIL,DISCH,IHSAPR,IBARR,E2,VEC4)
                  NMAT = NMAT + 1
                END IF
             END IF
         IFLAG = 1
         END IF
*
**      New search direction : -H*G
*
         IF(IHSAPR .EQ. 3 ) THEN
           CALL COPVEC(E1,AVEC,NVAR)
           CALL CLSKHE(E2,VEC1,AVEC,NVAR,IBARR,VEC4,2,INDEF)
         ELSE
           CALL VVTOV(E2,E1,VEC1,NVAR)
         END IF
C
         IZERO = 0
         IF(NMAT.NE.0)
     &   CALL SLRMTV(NMAT,NVAR,AMAT,AVEC,2,E1,
     &               VEC1,IZERO,DISCH,LUHFIL)
         CALL SCALVE(VEC1,-1.0D0,NVAR)
*
**      Line search along search direction
*
         E0PREV = E0
*. Set alpha to one always ! ?
         ALPHA = 1.0D0
*. Scale step to have length atmost one in iteration one
         IF(ITNUM.EQ.1) THEN
           STPNRM = SQRT(INPROD(VEC1,VEC1,NVAR))
           IF(STPNRM.GT.1.0D0) THEN
             CALL SCALVE(VEC1,1.0D0/STPNRM,NVAR)
             WRITE(6,*) ' Initial step scaled to norm 1 !!! '
           END IF
         END IF
*
           CALL LINSE3(ALPHA,VEC1,NVAR,MAXIT,ITNUML,X,E1,E0,IFLAG,
     &                 IQUACI,ISEPOP,NGRP,IGRPO,IGRPN)
*
**      Test for convergence
*
         GNRM = SQRT(INPROD(E1,E1,NVAR))
C        IF(GNRM.LE.CNVNRM.OR. ABS(E0-E0PREV) .LT. 1.0D-11 ) 
         IF(GNRM.LE.CNVNRM) 
     &   THEN
            CONVRG = .TRUE.
            ICONV = 1
         END IF
         WRITE(6,'(A,1X,E16.10,2(6X,E8.2))')
     &   '  >>> E0 E1NRM DELTA',E0,GNRM,E0-E0PREV
 
         IF(.NOT.CONVRG.AND.ITNUM.LT.MAXIT ) GOTO 100
*.       End of loop over iterations.
      IF(CONVRG) THEN
         WRITE(6,*) ' >>> Convergence obtained with quasin <<< '
         WRITE(6,'(A,I4)') ' Number of iterations used ',ITNUM
      END IF
*
      IF( IPRT .GT. 1 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' FINAL PARAMETER SET '
        CALL WRTMAT(X,1,NVAR,1,NVAR)
        WRITE(6,*)
C       WRITE(6,*) ' FINAL GRADIENT '
C       CALL WRTMAT(E1,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
      SUBROUTINE HESUPV (HDIAG,A,AVEC,X,G,XPREV,GPREV,NVAR,
     &                   IUPDAT,IINV,SCR,NMAT,LUHFIL,DISCH,
     &                   IHSAPR,IB,E2,VEC4)
 
*     Routine for variable metric update of hessian.
*     VECTOR BASED VERSION
*
*     On input :  current approximation to hessian CONSISTS OF
*                    1) hessian diagonal hdiag
*                    2) nmat rank-2 matrices defined by a and avec
*
*                x   current parameter set
*                g   current gradient
*              xprev : previous parameter set
*              gprev : previous gradient
*              nvar  : number of variables
*              iupdat: update method to be used
*                     = 1 :  DFP formulae
*                     = 2 : BFGS formulae
*              iinv  :  = 1 : update on inverse hessian (in use)
*                       = 1 : update on hessian ( to be published)
*               (scr : array of length at least nvar)
*     On output
*              new approximation to hessian STORED IN HDIAG AND
*              NMAT + 1 RANK-2 MATRICES
*             xprev : present parameter set( for use in next iteration)
*             gprev : present gradient  ( for use in next iteration)
*
* in current version are the hessian stored as an diagonal plus
* a sum of rank-two matrices.
* according to mister fletcher himself it is poissible to write
*
*   h(k+1) = h(k) + (delta, h(k)*gamma ) * a * (delta, h(k)*gamma) (t)
*            where delta and gamma are vectors defined by
*           delta = x - xprev
*           gamma = g - gprev
*
* and a is a  2 x 2 matrix with elements
*           a(1,1) = (1 +phi * gamhgam/delgem)/delgam
*           a(1,2) = -phi /delgam
*           a(2,1) = a(1,2)
*           a(2,2) = (phi - 1 ) /gamhgam
* with
*
* delgam = delta(t)*gamma
* gamhgam = gamma(t)*h(k)*gamma
* and phi depends on the method
*      for dfp :  phi = 0
*      for bfgs : phi = 1
*** ADDED ( APRIL 15 '87 ) : POSSIBILITY FOR SAVING VECTORS IN DISC
*   IF  DISCH IS TRUE VECTORS ARE READ IN FROM FILE LUHFIL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL * 8  INPROD
      LOGICAL DISCH
      DIMENSION HDIAG(NVAR),X(NVAR),G(NVAR),XPREV(NVAR),GPREV(NVAR),
     &          SCR(NVAR),AVEC(*),A(*)
      DIMENSION E2(*), VEC4(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from HESUPV: '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' IUPDAT,IINV = ', IUPDAT, IINV
        WRITE(6,*) ' NVAR, NMAT = ', NVAR, NMAT
        WRITE(6,*) ' DISCH ', DISCH
      END IF
*
      IF(IINV.EQ.1) THEN
*       inverse update
*       common preparations for  bfgs and dfp update
*       delta in xprev and gamma in gprev
        CALL VECSUM(XPREV,XPREV,X,-1.0D0,1.0D0,NVAR)
        CALL VECSUM(GPREV,GPREV,G,-1.0D0,1.0D0,NVAR)
*
*       h(k)*gamma in scr
*
        IF(IHSAPR .EQ. 3 ) THEN
          IADDR_FREE = 2*NMAT*NVAR + 1
          CALL COPVEC(GPREV,AVEC(IADDR_FREE),NVAR)
          CALL CLSKHE(E2,SCR,AVEC(IADDR_FREE),NVAR,IB,VEC4,2)
        ELSE
          CALL VVTOV(HDIAG,GPREV,SCR,NVAR)
        END IF
*
        IZERO = 0
        CALL SLRMTV(NMAT,NVAR,A,AVEC,2,GPREV,SCR,IZERO,DISCH,LUHFIL)
        DELGAM = INPROD(GPREV,XPREV,NVAR)
*       update only if delta(t)*gamma is positive ; only
*       under this assumption is positive definiteness conserved.
        ISKIP = 0
        IF(DELGAM .GT. 0.0D0 .OR. ISKIP .EQ. 0 ) THEN
          GAMHGA = INPROD(GPREV,SCR,NVAR)
          IF ( DELGAM .LT. 0.0D0) WRITE(6,*)
     &    '  WARNING QUASIN : POSITIVE DEFINITE H NOT ENSURED '
          IF ( IUPDAT .EQ. 1 ) PHI = 0
          IF ( IUPDAT .EQ. 2 ) PHI = 1.0D0
*
          IF(DISCH) THEN
            WRITE(LUHFIL) (XPREV(II),II=1,NVAR)
            WRITE(LUHFIL) (SCR(II),II=1,NVAR)
          ELSE
            IADR_DELTA = 2*NMAT*NVAR+1
            IADR_GAMMA = 2*NMAT*NVAR+NVAR+1
            CALL COPVEC(XPREV,AVEC(IADR_DELTA),NVAR)
            CALL COPVEC(SCR,  AVEC(IADR_GAMMA),NVAR)
          END IF
*
          IA = 4*NMAT
          A(IA+1) = (1 + PHI*GAMHGA/DELGAM)/DELGAM
          A(IA+2) = -PHI/DELGAM
          A(IA+3) = A(IA+2)
          A(IA+4) = (PHI-1)/GAMHGA
        ELSE
          WRITE(6,*) ' NO UPDATE PERFORMED SINCE ,DELGAM : ',DELGA
        END IF! DELGAM > 0
      END IF! IINV = 1
*
*   Prepare for next iteration, copy x to xprev, g to gprev
*
      CALL COPVEC(X,XPREV,NVAR)
      CALL COPVEC(G,GPREV,NVAR)
*
      RETURN
      END
      SUBROUTINE APRMAT(VECA,VECB,AB,NINBL,NBLK,NDIM)
C
C OBTAINED BLOCKED PROJECTED MATRIX
C AB(I,J) = Inprod(VECA(*,I),VECB(*,J))
C ONLY VECTORS BELONGING TO THE SAME BLOCK ARE MULTIPLIED
C AND RESULTING MATRIX IS ALSO BLOCKED
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NINBL(NBLK)
      DIMENSION AB(*)
      DIMENSION VECA(NDIM,*),VECB(NDIM,*)
      REAL * 8 INPROD
C
      IBLKBS = 1
      IABBS  = 1
      DO 300 IBLK = 1, NBLK
        LBLK = NINBL(IBLK)
        IF(IBLK.NE.1) IBLKBS = IBLKBS+NINBL(IBLK-1)
        IF(IBLK.NE.1) IABBS  = IABBS +NINBL(IBLK-1)**2
        DO 200 I = 1,LBLK
        DO 200 J = 1,LBLK
          AB(IABBS-1+(J-1)*LBLK+I) =
     &    INPROD(VECA(1,IBLKBS-1+I),VECB(1,IBLKBS-1+J),NDIM)
  200   CONTINUE
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE CTRN_MAT(XAXR,XAXI,XR,XI,AR,AI,NXR,NXC,NAR,NAC,SCR)
*
* Transform complex matrix A with complex matrix X
* XAX = (XR + i XI)^\dagger (AR + i AI)(XR + i XI) giving
*
* XAX_R = XR^T AR XR + XI^T AR XI + XI AI XR - XR AI XI
* XIX_I = XR^T AI XR + XR^R AI XR - XI^T AR XR +  XI AI XI
*
* Jeppe Olsen, implementing  gradient at general point in LUCIA, oct 2011
*
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION XR(NXR,NXC),XI(NXR,NXC)
      DIMENSION AR(NAR,NAC),AI(NAR,NAC)
*. Output
      DIMENSION XAXR(NXC,NXC), XAXI(NXC,NXC)
*. Scratch:
      DIMENSION SCR(NAR,NXC)
*. check that dimensions are consistent:
      IF(NXR.NE.NAR.OR.NAC.NE.NXR) THEN
       WRITE(6,*) ' Inconsistent dimensions in CTRN_MAT'
       WRITE(6,'(A,4I6)') ' NXR, NXC, NAR, NAC = ',
     &                      NXR, NXC, NAR, NAC
       STOP 'Inconsistent dimensions in CTRN_MAT'
      END IF
*
* AR XR
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(SCR,AR,XR,NAR,NXC,NAR,NAC,NXR,NXC,
     &     FACTORC, FACTORAB,0)
* XR(T) AR XR => XAXR
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(XAXR,XR,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
* -XI(T) AR XR => XAXI
      FACTORC = 0.0D0
      FACTORAB = -1.0D0
      CALL MATML7(XAXI,XI,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
* 
* AR XI
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(SCR,AR,XI,NAR,NXC,NAR,NAC,NXR,NXC,
     &     FACTORC, FACTORAB,0)
* XR(T) AR XI => XAXI
      FACTORC = 1.0D0
      FACTORAB = 1.0D0
      CALL MATML7(XAXI,XR,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
* XI(T) AR XI => XAXR
      FACTORC = 1.0D0
      FACTORAB = 1.0D0
      CALL MATML7(XAXR,XI,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
* 
* AI XR
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(SCR,AI,XR,NAR,NXC,NAR,NAC,NXR,NXC,
     &     FACTORC, FACTORAB,0)
* XR(T) AI XR => XAXI
      FACTORC = 1.0D0
      FACTORAB = 1.0D0
      CALL MATML7(XAXI,XR,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
* XI(T) AI XR => XAXR
      FACTORC = 1.0D0
      FACTORAB = 1.0D0
      CALL MATML7(XAXR,XI,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
* 
* AI XI
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(SCR,AI,XI,NAR,NXC,NAR,NAC,NXR,NXC,
     &     FACTORC, FACTORAB,0)
*-XR(T) AI XI => XAXR
      FACTORC = 1.0D0
      FACTORAB = -1.0D0
      CALL MATML7(XAXR,XR,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
* XI(T) AI XI => XAXI
      FACTORC = 1.0D0
      FACTORAB = 1.0D0
      CALL MATML7(XAXI,XI,SCR,NXC,NXC,NXR,NXC,NAR,NXC,
     &     FACTORC,FACTORAB,1)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from CTRN_MAT '
       WRITE(6,*) ' ====================='
       WRITE(6,*) ' Real part of transformed matrix '
       CALL WRTMAT(XAXR,NXC,NXC)
       WRITE(6,*) ' Imaginary part of transformed matrix '
       CALL WRTMAT(XAXI,NXC,NXC)
      END IF
*
      RETURN
      END
      SUBROUTINE WRTCMAT(CMATR, CMATI, NROW, NCOL)
*
* Write real and imaginary parts of general complex matrix CMAT
*
* Jeppe Olsen, Oct. 2011
*
      INCLUDE 'implicit.inc'
*
      DIMENSION CMATR(NROW,NCOL), CMATI(NROW,NCOL)
*
      WRITE(6,*) ' Real part of matrix: '
      CALL WRTMAT(CMATR,NROW,NCOL,NROW,NCOL)
      WRITE(6,*)
      WRITE(6,*) ' Imaginary part of matrix: '
      CALL WRTMAT(CMATI,NROW,NCOL,NROW,NCOL)
*
      RETURN
      END
      SUBROUTINE EXC_VEC_FROM_GO_MAT(EXC_VEC, GOMAT,IJSM,
     &           NOOEXC,IOOEXCC,IOOEXC,
     &           NSMOB,NOCOBS,NTOOBS,NTOOB,IBSO,IREOST)
*
* A matrix GOMAT(R,S) of symmetry IJSM with R:  general, S: occupied index
* is given
* Fetch the elements specified by the orbital excitation array
* IOOEXC, and save in EXC_VEC
*
*. Jeppe Olsen, October 2011- Never to old to reform matrices
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. General input
      INTEGER IOOEXC(NTOOB,NTOOB), IBSO(NSMOB), IREOST(*)
      INTEGER IOOEXCC(2,NOOEXC)
      INTEGER NOCOBS(NSMOB),NTOOBS(NSMOB)
*. Specific input
      DIMENSION GOMAT(*)
*. Output
      DIMENSION EXC_VEC(*)
*. Loop over orbitals in GOMAT- symmetry order
      IJOFF = 1
      DO ISM = 1, NSMOB
       IF(ISM.EQ.1) THEN
         IJOFF = 1
       ELSE
         JSM_PREV = MULTD2H(ISM-1,IJSM)
         IJOFF = IJOFF + NTOOBS(ISM-1)*NOCOBS(JSM_PREV)
       END IF
*
       JSM = MULTD2H(ISM,IJSM)
       NI = NTOOBS(ISM)
       NJ = NOCOBS(JSM)
*
       IB = IBSO(ISM)
       JB = IBSO(JSM)
*
       DO IORB_S = 1, NI
       DO JORB_S = 1, NJ
*. Type-order - as used in IOOEXC
        IORB_T = IREOST(IB - 1 + IORB_S)
        JORB_T = IREOST(JB - 1 + JORB_S)
        IF(IOOEXC(IORB_T,JORB_T).GT.0) THEN
          IJEXC = IOOEXC(IORB_T,JORB_T)
          EXC_VEC(IJEXC) = GOMAT(IJOFF-1+(JORB_S-1)*NI+IORB_S)
        END IF
       END DO
       END DO
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' EXC_VEC from EXC_VEC_FROM_GO_MAT '
        CALL WRT_EXCVEC(EXC_VEC,IOOEXCC,NOOEXC)
      END IF
*
      RETURN
      END
      SUBROUTINE NEWCI 
        STOP ' Entered Dummy NEWCI routine '
      RETURN
      END 
      SUBROUTINE GRADIE
        STOP ' Entered Dummy GRADIE routine '
      RETURN
      END 
      SUBROUTINE NEWMOS
        STOP ' Entered Dummy NEWMOS routine '
      RETURN
      END 
      SUBROUTINE HESINV
        STOP ' Entered Dummy HESINV routine '
      RETURN
      END 
      SUBROUTINE GET_BRT_FROM_F(BRT,F)
*
* Obtain Brillouin Vector in complete symmetry-blocked
* matrix form from Fock matrix
*
* B(R,S) = 2F(R,S) - 2F(S,R)
*
* Jeppe Olsen, November 2011
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
CNW   DIMENSION F(*)
      integer F
*. Output
      DIMENSION BRT(*)
*
      call luciaqenter('gbrtf')
      IOFF = 1
      call ga_distribution(F,ga_nodeid(),ilo,ihi,idum,idum)
      if (ihi.gt.0) then
      call ga_access(F, ilo, ihi, 1, 1, ikf, idum)
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
         IOFF = 1
        ELSE
         IOFF = IOFF + NTOOBS(ISM-1)**2
        END IF
        NOB = NTOOBS(ISM)
        DO IR = 1, NOB
        DO IS = 1, NOB 
           index1 = IOFF-1+(IR-1)*NOB+IS
           index2 = IOFF-1+(IS-1)*NOB+IR
           XVAL = 0.0d0
           if (index1.ge.ilo.and.index1.le.ihi)
     &     XVAL = XVAL + dbl_mb(ikf+index1-ilo+1)
           if (index2.ge.ilo.and.index1.le.ihi)
     &     XVAL = XVAL - dbl_mb(ikf+index2-ilo+1)
           BRT(IOFF-1+(IS-1)*NOB+IR) = 2.0D0*(XVAL)
        ENDDO
        ENDDO
CNW     DO IR = 1, NOB
CNW     DO IS = 1, NOB
CNW       BRT(IOFF-1+(IS-1)*NOB + IR) = 2.0D0*
CNW  &   (F(IOFF-1+(IS-1)*NOB + IR)-F(IOFF-1+(IR-1)*NOB + IS))
CNW     END DO
CNW     END DO
      END DO
      IOFF = IOFF + NOB*NOB - 1
      call ga_release(F,ilo,ihi,1,1)
      endif
      call ga_sync()
      call ga_dgop(1, BRT, IOFF, '+')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Brillouin matrix in complete form'
       WRITE(6,*) ' ================================='
       CALL APRBLM2(BRT,NTOOBS,NTOOBS,NSMOB,0)
      END IF
      call luciaqexit('gbrtf')
*
      RETURN
      END
      SUBROUTINE PROJ_ORBSPC_ON_ORBSPC(CMOAO1,CMOAO2,NMO1PSM,NMO2PSM)
*
* Project orbitals CMOAO1 on Orbitals CMOAO2 and find norm 
* of resulting orbitals
*
* 1: Obtain X = CMOAO1 CMOAO1(T) SAO CMOOA2
* 2: Obtain X(T)SX and print norms (X(T)SX)
*
*. Jeppe Olsen, November 2011
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input:
*. Dimension  of orbital sets
      INTEGER NMO1PSM(NSMOB),NMO2PSM(NSMOB)
*. And the orbital expansions
      DIMENSION CMOAO1(*), CMOAO2(*)
*
      NTEST = 10
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Info from PROJ_ORBSPC_ON_ORBSPC '
        WRITE(6,*) ' =============================== '
      END IF
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' MO1: ' 
        CALL APRBLM2(CMOAO1,NTOOBS,NMO1PSM,NSMOB,0)
        WRITE(6,*) ' MO2: ' 
        CALL APRBLM2(CMOAO2,NTOOBS,NMO2PSM,NSMOB,0)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PROBOB')
*
      LEN_CMO = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLMAT1,LEN_CMO,'ADDL  ',2,'MAT1  ')
      CALL MEMMAN(KLMAT2,LEN_CMO,'ADDL  ',2,'MAT2  ')
      CALL MEMMAN(KLMAT3,LEN_CMO,'ADDL  ',2,'MAT3  ')
      CALL MEMMAN(KLSAOE,LEN_CMO,'ADDL  ',2,'SAOE  ')
*. Obtain SAO in expanded form
      XDUM = 2810.1979
      CALL GET_HSAO(XDUM,dbl_mb(KSAO),0,1)
C          GETHSAO(HAO,SAO,IGET_HAO,IGET_SAO)
*. Obtain SAO in expanded (unpacked form)
      CALL TRIPAK_AO_MAT(dbl_mb(KLSAOE),dbl_mb(KSAO),2)
*S C2 IN MAT1
C     MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &              LAROW,LACOL,LBROW,LBCOL,ITRNSP)
      CALL MULT_BLOC_MAT(dbl_mb(KLMAT1),dbl_mb(KLSAOE),CMOAO2,NSMOB,
     &     NTOOBS,NMO2PSM,NTOOBS,NTOOBS,NTOOBS,NMO2PSM,0)
* C1(T) (S C2) in MAT2
      CALL MULT_BLOC_MAT(dbl_mb(KLMAT2),CMOAO1,dbl_mb(KLMAT1),NSMOB,
     &     NMO1PSM,NMO2PSM,NTOOBS,NMO1PSM,NTOOBS,NMO2PSM,1)
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' C1(T) S C2: '
        CALL APRBLM2(dbl_mb(KLMAT2),NMO1PSM,NMO2PSM,NSMOB,0)
      END IF
* C1 (C1(T) S C2) in MAT1
      CALL MULT_BLOC_MAT(dbl_mb(KLMAT1),CMOAO1,dbl_mb(KLMAT2),NSMOB,
     &     NTOOBS,NMO2PSM,NTOOBS,NMO1PSM,NMO1PSM,NMO2PSM,0)
*
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' The matrix C1 C1(T) S C2 '
        CALL APRBLM2(dbl_mb(KLMAT1),NTOOBS,NMO2PSM,NSMOB,0)
      END IF
*. X = C1 C1(T) S C2 is in MAT1, obtain S X in MAT2
      CALL MULT_BLOC_MAT(dbl_mb(KLMAT2),dbl_mb(KLSAOE),dbl_mb(KLMAT1),
     &     NSMOB,NTOOBS,NMO2PSM,NTOOBS,NTOOBS,NTOOBS,NMO2PSM,0)
* X(T) S X in MAT3
      CALL MULT_BLOC_MAT(dbl_mb(KLMAT3),dbl_mb(KLMAT1),dbl_mb(KLMAT2),
     &     NSMOB,NMO2PSM,NMO2PSM,NTOOBS,NMO2PSM,NTOOBS,NMO2PSM,1)
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' The matrix X(T) S X '
        CALL APRBLM2(dbl_mb(KLMAT3),NMO2PSM,NMO2PSM,NSMOB,0)
      END IF
*
C GET_DIAG_BLMAT(A,DIAG,NBLK,LBLK,ISYM)
      CALL GET_DIAG_BLMAT(dbl_mb(KLMAT3),dbl_mb(KLMAT2),NSMOB,NMO2PSM,0)
*. And print
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Part of CMO2 that is spanned by CMO1'
        CALL PRINT_SCALAR_PER_ORB(dbl_mb(KLMAT2),NMO2PSM)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PROBOB')
*
      RETURN
      END
      SUBROUTINE PRINT_SCALAR_PER_ORB(SCALAR,NOBPSM)
*
* Print a scalar for each orbital from orbital set with 
* dimension NOBPSM
*
* Jeppe Olsen
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
*. General Input
      INTEGER NOBPSM(NSMOB)
*. Specific input
       DIMENSION SCALAR(*)
*
      IOFF = 1
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF  = 1
        ELSE
          IOFF = IOFF + NOBPSM(ISM-1)
        END IF
        IF(NOBPSM(ISM).NE.0) THEN 
          WRITE(6,'(A,I5)') ' Symmetry = ', ISM
          WRITE(6,*)        ' ==================='
          NOB = NOBPSM(ISM)
          CALL WRTMAT(SCALAR(IOFF),1,NOB,1,NOB)
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE PRINT_SCALAR_PER_ORB2(SCALAR,NOBPSMX,NOBSMX)
*
* Print a scalar for each orbital from orbital set with 
* dimension NOBPSMX and NOSMX symmetries
*
* Jeppe Olsen
*
*. Last modification: July 8, 2012
      INCLUDE 'implicit.inc'
*. General Input
      INTEGER NOBPSMX(NOBSMX)
*. Specific input
       DIMENSION SCALAR(*)
*
      IOFF = 1
      DO ISM = 1, NOBSMX
        IF(ISM.EQ.1) THEN
          IOFF  = 1
        ELSE
          IOFF = IOFF + NOBPSMX(ISM-1)
        END IF
        IF(NOBPSMX(ISM).NE.0) THEN 
          WRITE(6,'(A,I5)') ' Symmetry = ', ISM
          WRITE(6,*)        ' ==================='
          NOB = NOBPSMX(ISM)
          CALL WRTMAT(SCALAR(IOFF),1,NOB,1,NOB)
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE E123_ALONG_MODE(EFUNC,XNOT,XMODE,NVAR,E1,E2,E3)
*
* A function,EFUNC,  depending on a set of parameter X is
* defined. Obtain first three directional derivatives
* along direction XMODE for initial values of X given 
* by X
*
*. Jeppe Olsen, Aug. 30 2012
*. Last revision, Aug. 30. 2012, Written
*
* For mode walking out from stationary point
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION XNOT(NVAR),XMODE(NVAR)
*
      EXTERNAL EFUNC
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from E123_ALONG_MODE'
       WRITE(6,*) ' ========================='
       WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Input mode '
       CALL WRTMAT(XMODE,1,NVAR,1,NVAR)
      END IF
*. Steplength
      DELTA = 1.0D-2
*. Energy at point of expansion
      E0 = EFUNC(XNOT)
*.E(XNOT+XMODE)
      CALL VECSUM(XNOT,XNOT,XMODE,1.0D0,DELTA,NVAR)
      EP1 = EFUNC(XNOT)
*.E(XNOT+2*XMODE)
      CALL VECSUM(XNOT,XNOT,XMODE,1.0D0,DELTA,NVAR)
      EP2 = EFUNC(XNOT)
*.E(XNOT-XMODE)
      CALL VECSUM(XNOT,XNOT,XMODE,1.0D0,-3.0D0*DELTA,NVAR)
      EM1 = EFUNC(XNOT)
*.E(XNOT-2*XMODE)
      CALL VECSUM(XNOT,XNOT,XMODE,1.0D0,-1.0D0*DELTA,NVAR)
      EM2 = EFUNC(XNOT)
*. and restore
      CALL VECSUM(XNOT,XNOT,XMODE,1.0D0,2.0D0*DELTA,NVAR)
*. Elementary finite difference equations
      E1 = (8.0D0*EP1-8.0D0*EM1-EP2+EM2)/(12.0D0*DELTA)
      E2 = (16.0D0*(EP1+EM1-2.0D0*E0)-EP2-EM2+2.0D0*E0)/
     &     (12.0D0*DELTA**2)
      E3 = (EP2-EM2-2.0D0*(EP1-EM1))/2*DELTA**3
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from E123_ALONG_MODE '
        WRITE(6,*) ' Finite difference to the first 3 derivatives ',
     &             E1,E2,E3
      END IF
*
      RETURN
      END
c $Id$
