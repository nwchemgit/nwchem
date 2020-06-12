      SUBROUTINE LUCIA_NORT(I_DO_NONORT_MCSCF,
     &           JCMBSPC,E_FINAL,CONV_F,ERROR_NORM_FINAL,INI_NORT,
     &           IVBGNSP,IVBGNSP_PREV)
*
* Perform Nonorthogonal CI calculation 
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'vb.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'fragmol.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
      LOGICAL CONV_F, CONV_NORTCI
*
      IPRVB = 10
      NTEST = 1000
      NTEST = MAX(IPRVB, NTEST)
*
      WRITE(6,*) ' *************************************** '
      WRITE(6,*) ' *                                     * '
      WRITE(6,*) ' *  Nonorthogonal section entered      * '
      WRITE(6,*) ' *                                     * '
      WRITE(6,*) ' *  Jeppe Olsen                        * ' 
      WRITE(6,*) ' *                                     * '
      WRITE(6,*) ' *  Version of June 2013 ( 0.96)       * '
      WRITE(6,*) ' *************************************** '
*
      WRITE(6,*) ' TEST: INI_NORT, IRESTR = ', INI_NORT, IRESTR
      WRITE(6,*) ' TEST: IVBGNSP, IVBGNSP_PREV = ', 
     &                   IVBGNSP, IVBGNSP_PREV
*
      IF(IVBGNSP.NE.0) THEN
*. Copy general space to reference VB space
        NORBVBSPC = NOBPT(NORTCIX_SCVB_SPACE)
        CALL ICOPVE(VB_GNSPC_MIN(1,IVBGNSP),VB_REFSPC_MIN(1),NORBVBSPC)
        CALL ICOPVE(VB_GNSPC_MAX(1,IVBGNSP),VB_REFSPC_MAX(1),NORBVBSPC)
      END IF
*
      IF(NTEST.GE.0) THEN
       WRITE(6,*) ' Information on nonorthogonal calculation: '
       WRITE(6,*) ' ==========================================' 
       WRITE(6,*)
*
       IF(NORT_MET.EQ.1) THEN
        WRITE(6,'(5X,A)') 
     &  ' Non-orthogonal wave function will be expanded in CI space'
       ELSE IF( NORT_MET.EQ.2) THEN
        WRITE(6,'(5X,A)') 
     &  ' Non-orthogonal wave function will be expanded configurations'
       ELSE
         WRITE(6,*) ' Currently unknown NORT_MET = ', NORT_MET
         STOP ' Currently unknown NORT_MET '
       END IF
*
       WRITE(6,'(5X,A,I3)') 
     &  ' Orbital space for non-orthogonal calculation:', 
     &    NORTCIX_SCVB_SPACE
       WRITE(6,'(5X,A,I3)') 
     & ' Allowed excitation level from Spin-coupled valence space',
     &   NORTCI_SCVB_EXCIT
       WRITE(6,'(5X,A,I3)') 
     &  ' Spanning CI-space:', JCMBSPC
*
       WRITE(6,'(5X,A)') 
     &' Min and max accumulated occupation in valence ref CI space: '
       NORBVBSPC = NOBPT(NORTCIX_SCVB_SPACE)
       DO IORB = 1, NORBVBSPC
        WRITE(6,'(10X,2I3)') VB_REFSPC_MIN(IORB),VB_REFSPC_MAX(IORB)
       END DO
      END IF !NTEST is large enough for printing
*
      IF(NTEST.GE.10) THEN
       IF(INI_NORT.EQ.1) THEN
C       WRITE(6,*) ' INI_MO_TP, INI_MO_ORT = ', INI_MO_TP, INI_MO_ORT
        WRITE(6,*)
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Initial set of orbitals '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
*
        IF(INI_MO_TP.EQ.1) THEN
          WRITE(6,'(4X,A)') ' Atomic orbitals will be used '
        ELSE IF (INI_MO_TP.EQ.2) THEN
          WRITE(6,'(4X,A)') 
     &    ' Input MOs in VB space rotated  to give diagonal block'
        ELSE IF (INI_MO_TP.EQ.3) THEN
          WRITE(6,'(4X,A)') 
     &    ' Initial MO orbitals from SIRIFC will be used'
        ELSE IF (INI_MO_TP.EQ.4) THEN
          WRITE(6,'(4X,A)') 
     &    ' Constructed from fragment orbitals'
        END IF
        WRITE(6,'(4X,A)') 
     &  ' Orbitals in inactive and secondary space will be ort.'
        WRITE(6,'(4X,A)') ' Orbitals in GAS orbital spaces(.ne. VB ): '
        IF(INI_MO_ORT.EQ.0) THEN
          WRITE(6,'(6X,A)') ' No orthogonalization  '
        ELSE IF (INI_MO_ORT.EQ.1) THEN
          WRITE(6,'(6X,A)') ' Orthogonalized'
        END IF
        WRITE(6,'(4X,A)') ' Orbitals in VB orbital space: '
        IF(INI_ORT_VBGAS.EQ.0) THEN
          WRITE(6,'(6X,A)') ' No orthogonalization  '
        ELSE IF (INI_ORT_VBGAS.EQ.1) THEN
          WRITE(6,'(6X,A)') ' Orthogonalized'
        END IF
*
        IF(INI_MO_TP.EQ.4) THEN
         WRITE(6,*) ' Distribution of orbitals from fragments:'
         DO IFRAG = 1, NFRAG_MOL
          NSMOB_L = NSMOB_FRAG(IFRAG)
          WRITE(6,'(A,I3)') ' For fragment ', IFRAG
          WRITE(6,*)        ' ===================='
          WRITE(6,*) ' Number of orbitals per GAS (row) and sym (col) '
          CALL IWRTMA
     &    (N_GS_SM_BAS_FRAG(0,1,IFRAG),NGAS+2,NSMOB_L,MXPNGAS+1,MXPOBS)
         END DO
        END IF ! End if INI_MO_TP.eq.4
       ELSE
        WRITE(6,*) ' Start from orbitals in place '
       END IF
*
        IF(IRESTR.EQ.0) THEN
         WRITE(6,*)
         WRITE(6,*) ' ======================= '
         WRITE(6,*) ' Initial configuration: '
         WRITE(6,*) ' ======================= '
         WRITE(6,*)
         IF(I_HAVE_INI_CONF.EQ.0) THEN
           WRITE(6,'(5X,A)') ' None given '
         ELSE
           WRITE(6,'(5X,A)') ' In compressed form '
           CALL IWRTMA(INI_CONF,1,NOB_INI_CONF,1,NOB_INI_CONF)
         END IF
        ELSE
          WRITE(6,*) ' Restarted calculation '
        END IF
*
      END IF ! NTEST is large enough for testoutput
*
* Some general info on configuration expansions
*
*. First orbital and number of electrons in VB orbital space
      IB_VBOBSPC= NINOB + 1
      DO IOBSPC = 1, NORTCIX_SCVB_SPACE-1
        IB_VBOBSPC = IB_VBOBSPC + NOBPT(IOBSPC)
      END DO
      NORBVBSPC = NOBPT(NORTCIX_SCVB_SPACE)
      IF(NTEST.GE.10)
     &WRITE(6,*) ' Dimension and offset for orbitals in VB-space',
     &             NORBVBSPC,IB_VBOBSPC
*. Number of electrons
      NELEC = VB_REFSPC_MIN(NORBVBSPC)
      IF(NTEST.GE.10) WRITE(6,*) ' Test: NELEC = ', NELEC
*. Save for communication with configuration routines
      IB_ORB_CONF = IB_VBOBSPC
      N_ORB_CONF = NORBVBSPC
      N_EL_CONF = NELEC
*
*. Check number of electrons in initial configuration
*
      IF(I_HAVE_INI_CONF.EQ.1) THEN
       NEL_INI = NEL_IN_COMPACT_CONF(INI_CONF,NOB_INI_CONF)
       IF(NEL_INI.NE.NELEC) THEN
         WRITE(6,*) 
     &   ' Incorrect number of electrons in initial configuration'
         WRITE(6,*) ' Actual and required number of electrons ',
     &               NEL_INI, NELEC
         STOP
     &   ' Incorrect number of electrons in initial configuration'
       END IF
      END IF
*
* =========================================================
* information about prototype configurations in  CI space
* =========================================================
*
*
*. Max. and min. number of open orbitals- based only number of orbitals
* and electrons
*. And the prototype information
*
* ======================================
* The various min-max occupation spaces
* ======================================
*
* Space 1: The reference space for |0>
       ICSPC_CNF = 1
       CALL ICOPVE(VB_REFSPC_MIN,IOCC_MIN_GN(1,ICSPC_CNF),NORBVBSPC)
       CALL ICOPVE(VB_REFSPC_MAX,IOCC_MAX_GN(1,ICSPC_CNF),NORBVBSPC)
*. Space 2: Space where Hamiltonian vector will be calculated, currently
*      also reference space   
       ISSPC_CNF = 2
       CALL ICOPVE(VB_REFSPC_MIN,IOCC_MIN_GN(1,ISSPC_CNF),NORBVBSPC)
       CALL ICOPVE(VB_REFSPC_MAX,IOCC_MAX_GN(1,ISSPC_CNF),NORBVBSPC)
*. Space 3: Intermediate space where |0> is expanded in biothonormal basis,
*. Must interact with final space (2) through a given level of excit
       IMSPC_CNF = 3
*. For atmost two-body operators
       IF(NORT_M.EQ.1) THEN
         NEXCIT = 2
         NEXCIT = NELEC
*. I have been having some errors with orb gradient when reordering
*. orbitals, so I have increased this in the aboce
         WRITE(6,*) ' IMPORTANT: NEXCIT raised to NELEC for test'
       ELSE
         NEXCIT = 2
       END IF
       CALL MINMAX_EXCIT(
     &      IOCC_MIN_GN(1,ISSPC_CNF),IOCC_MAX_GN(1,ISSPC_CNF),NEXCIT,
     &      IOCC_MIN_GN(1,IMSPC_CNF),IOCC_MAX_GN(1,IMSPC_CNF),
     &      NORBVBSPC)
      NVBCISPC = 3
      NVBCNSPC =  NVBCISPC
      IB_INTM_SPC = NVBCISPC + 1

      IF(NORT_MET.EQ.2) THEN
*. The bioorthogonal C vector will be obtained as a 
*. sequence of one-orbital transformations. Generate spaces for these
         N_INTM_SPC = N_ORB_CONF
*.Is there enough space for pointers
         IF(NVBCISPC+N_INTM_SPC.GE.MXPICI) THEN
           WRITE(6,*) ' Too many intermediate MAXMIN spaces required'
           WRITE(6,*) ' Needed number of spaces ',  N_ORB_CONF
           WRITE(6,*) ' Present number of spaces ',  MXPICI - NVBCISPC
           WRITE(6,*) ' Increase MXPICI and recompile '
           STOP       ' Too many intermediate MAXMIN spaces required'
         END IF
*. Generate the various MAXMIN spaces and their dimensions
C        MINMAX_FOR_ORBTRA(MIN_IN,MAX_IN,MIN_OUT,MAX_OUT,
C    &   MIN_INTM,MAX_INTM,MIN_INTMS,MAX_INTMS,ISYM,IDODIM)
         IDODIM = 1
         WRITE(6,*) ' ICSPC_CNF, IOCC_MIN_GN(1,ICSPC_CNF) = ',
     &                ICSPC_CNF, IOCC_MIN_GN(1,ICSPC_CNF) 
         WRITE(6,*) 
     &  ' Configuration information for orbital transformation'
         WRITE(6,*) 
     &  ' ===================================================='
         CALL MINMAX_FOR_ORBTRA(
     &                   IOCC_MIN_GN(1,ICSPC_CNF),
     &                   IOCC_MAX_GN(1,ICSPC_CNF),
     &                   IOCC_MIN_GN(1,IMSPC_CNF),
     &                   IOCC_MAX_GN(1,IMSPC_CNF),
     &                   IOCC_MIN_GN(1,IB_INTM_SPC),
     &                   IOCC_MAX_GN(1,IB_INTM_SPC),
     &                   ISYM,IDODIM,NCONF_GN(IB_INTM_SPC),
     &                   NCSF_GN(IB_INTM_SPC),
     &                   NSD_GN(IB_INTM_SPC))
*. In and out spaces for the orbital transformation
        DO IORB = 1, N_ORB_CONF
         IF(IORB.EQ.1) THEN
            IORBTRA_SPC_IN(IORB) = ICSPC_CNF
         ELSE
            IORBTRA_SPC_IN(IORB) = IORBTRA_SPC_OUT(IORB-1)
         END IF
         IF(IORB.LT.N_ORB_CONF) THEN
           IORBTRA_SPC_OUT(IORB) = IB_INTM_SPC - 1 + IORB
         ELSE
           IORBTRA_SPC_OUT(IORB) = IMSPC_CNF
         END IF
        END DO
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' In and out spaces for the orbital trans '
          WRITE(6,*) ' ======================================= '
          WRITE(6,*) 
          WRITE(6,*) ' Orbital Inspace Outspace '
          WRITE(6,*) ' ========================='
          DO IORB = 1, N_ORB_CONF
            WRITE(6,'(3(I3,4X))') 
     &      IORB, IORBTRA_SPC_IN(IORB), IORBTRA_SPC_OUT(IORB)
          END DO
        END IF
*
        NVBCNSPC = NVBCNSPC+N_INTM_SPC
*
*. Largest number of CSFs of given sym in a CI space
*
      END IF! NORT_MET = 2
*
* ==================================================
* Generate configurations for the active CN spaces
* ==================================================
*
      NCSF_MNMX_MAX = 0
      DO ISPC = 1, NVBCNSPC
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 
          WRITE(6,*) ' ========================================'
          WRITE(6,*) ' Information about MINMAX space= ', ISPC
          WRITE(6,*) ' ========================================'
          WRITE(6,*)
        END IF
*
        CALL GEN_CONF_FOR_MINMAX_SPC(
     &      IOCC_MIN_GN(1,ISPC),IOCC_MAX_GN(1,ISPC),
     &       NORBVBSPC, IREFSM,IB_VBOBSPC,ISPC)
*. Configurations are returned in WORK(KICONF_OCC_GN(IREFSM,ISPC))
*. Number of SD's ..
C            NPARA_FOR_MINMAX_SPC(NCONF_OP,NCSF,NSD,NCMB)
        CALL NPARA_FOR_MINMAX_SPC(NCONF_PER_OPEN_GN(1,IREFSM,ISPC),
     &       NCSF,NSD,NCMB,NCNF)
        NSD_PER_SYM_GN(IREFSM,ISPC) = NSD
        NCSF_PER_SYM_GN(IREFSM,ISPC) = NCSF
        NCONF_PER_SYM_GN(IREFSM,ISPC) = NCNF
*
        NCSF_MNMX_MAX  = MAX(NCSF_MNMX_MAX,NCSF)
*
        IF(NORT_MET.EQ.1) THEN
*
* =======================================================
*. Generate mapping of SD's from configuration order to 
*. standard string order
* =======================================================
*
*. Obtain information about reexpansion in CI space
* Reorder array for determinants, index and sign
          CALL MEMMAN(KSDREO_I_GN(IREFSM,ISPC),NSD,'ADDL  ',1,'SDREOI')
*. Offsets for determinants with a given numbner of open orbitals
*. The code  below is a but confusing, I am not sure of its use..
          IZERO = 0
          CALL ISETVC(IB_SD_OPEN_GN(1,ISPC),IZERO,MAXOP+1)
          IB = 1
          DO IOPEN = MINOP, MAXOP
            IB_SD_OPEN_GN(IOPEN+1,ISPC) = IB
            IF(MOD(IOPEN-MS2,2).EQ.0) THEN
              IB = IB +
     &        NCONF_PER_OPEN_GN(IOPEN+1,IREFSM,ISPC)*NPCMCNF(IOPEN+1)
            END IF
          END DO

*. Reorder array for determinants, index and sign
          CALL MEMMAN(KSDREO_I_GN(IREFSM,ISPC),NSD,'ADDL  ',1,'SDREOI')
*. And then the reordering
C     CNFORD2(ISM,ICTSDT,ICONF_OCC,NCONF_PER_OP,
C    &           IDFTP,ICONF_ORBSPC)
          CALL CNFORD2(IREFSM,WORK(KSDREO_I_GN(IREFSM,ISPC)),
     &                 WORK(KICONF_OCC_GN(IREFSM,ISPC)),
     &                 NCONF_PER_OPEN_GN(1,IREFSM,ISPC),
     &                 WORK(KDFTP),NORTCIX_SCVB_SPACE,
     &                 JCMBSPC) 
        ENDIF ! End if NORTCI = 1
      END DO ! End of loop over CI spaces
*
      WRITE(6,*) ' Largest number of CSF''s in a space ',
     &             NCSF_MNMX_MAX
*
* =============================================================
* Generate atom orbitals and integrals over these orbitals
* ==============================================================
*
* At the moment: It is assumed that integrals have been 
* delivered in an orthogonal basis defined by C(MOAO) in WORK(KMOAOIN). 
* Obtain matrix for transforming from MO's to AO's
* and backtransform integrals....
*
* IN MOAOIN we actually have the actual expansion of the set of non-orthoginal
* orbitals that we will use. Save this, and read in original copy of C(MOAO)
      LENC = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0) 
      CALL COPVEC(WORK(KMOAOIN),WORK(KMOAO_ACT),LENC)
      CALL GET_CMOAO_ENV(WORK(KMOAOIN))
*
*. Allocate space for  H in AO basis
      LEN1E = NTOOB **2
      IF(NTEST.GE.1000)
     &WRITE(6,*) ' NTOOB, LEN1E = ', NTOOB, LEN1E
      CALL MEMMAN(KLHAO,LEN1E,'ADDL  ',2,'H_AO  ')
*. Allocate space for inverse of C(MOAO)
      CALL MEMMAN(KLCAOMO,LEN1E,'ADDL  ',2,'CAOMO ')
*. Obtain AO integrals SAO
      XDUM = 2810.1979
      CALL GET_HSAO(XDUM,WORK(KSAO),0,1)
C          GETHSAO(HAO,SAO,IGET_HAO,IGET_SAO)
*. Obtain SAO in expanded (unpacked form)
C?    WRITE(6,*) ' LEN1E = ', LEN1E
      CALL MEMMAN(KLSAOE,LEN1E,'ADDL  ',2,'S_AO_E')
C TRIPAK_AO_MAT(AUTPAK,APAK,IWAY)
*
      CALL TRIPAK_AO_MAT(WORK(KLSAOE),WORK(KSAO),2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' SAOE: '
        CALL APRBLM2(WORK(KLSAOE),NTOOBS,NTOOBS,NSMOB,0)
        WRITE(6,*) ' MOAOIN: '
        CALL APRBLM2(WORK(KMOAOIN),NTOOBS,NTOOBS,NSMOB,0)
      END IF

*. CMOAO(T) * SAO - it is assumed that CMOAO is in KMOAOIN
      CALL MULT_BLOC_MAT(WORK(KLCAOMO),WORK(KMOAOIN),WORK(KLSAOE),
     &     NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,1)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' C(AOMO) matrix: '
        CALL WRTVH1(WORK(KLCAOMO),1,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. And clean up
      CALL COPVEC(WORK(KMOAO_ACT),WORK(KMOAOIN),LENC)
*. 
*
*.The two-electron integrals in the AO basis - only done in initial NORT
*
      IF(INI_NORT.EQ.1) THEN
*
       IF(NOMOFL.EQ.1) THEN
        WRITE(6,*) 
     &  ' Lucia is trying to make a MO=>AO transformation of integrals'
        WRITE(6,*)
     &  ' But there is no AO=> MO transformation present'
        STOP ' NORTCI: NO AO => MO transformation matrix present'
       END IF
*
       IF(NTEST.GE.10) WRITE(6,*) ' Integral transformation:'
*. Input integrals in place for integral transformation
       KINT2 = KINT_2EMO
       CALL COPVEC(WORK(KH),WORK(KINT1O), NINT1)
*. Flag type of integral list to be obtained: Pt complete list of integrals
       IE2LIST_A = IE2LIST_FULL
       IOCOBTP_A = 1
       INTSM_A = 1
       KKCMO_I = KLCAOMO
       KKCMO_J = KLCAOMO
       KKCMO_K = KLCAOMO
       KKCMO_L = KLCAOMO
       IH1FORM = 1
       IH2FORM = 1
       CALL TRAINT
*. Move integrals in AO basis to KINT_2EMO (sorry for the name..)
       IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
       NINT2_F = NINT2_G(IE2ARR_F)
       KINT2_F = KINT2_A(IE2ARR_F)
       CALL COPVEC(WORK(KINT2_F),WORK(KINT_2EMO),NINT2_F)
C?     WRITE(6,*) ' NINT2_F = ', NINT2_F
C?     WRITE(6,*) ' Integrals transformed to KINT_2EMO'
C?     CALL WRTMAT(WORK(KINT_2EMO),1,NINT2_F,1,NINT2_F)
*. one-electron AO integrals to KINT1O
       CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
*
* End of generation of integrals over atomic orbitals: We have now in KINT_2EMO the 
* two-electron integrals over AO's and in KINT1O, the one-electron integrals in the AP basis.
      ELSE
       WRITE(6,*) ' AO integrals assumed in place '
      END IF
*
* ======================================
*. Obtain initial set of  orbitals
* ======================================
*
* Two steps : 1) Obtain a set of (nonorthogonal) initial orbitals
*             2) Perform (partial) orthonormalization to obtain 
*                Final initial orbitals
*
*. 1: Generate/Read in the initial orbitals
* Generate set of (nonorthogonal) initial orbitals
*
      IF(INI_NORT.EQ.1) THEN
       CALL GET_INIMO(WORK(KMOAOUT))
C           GET_INIMO(CMO_INI)
      ELSE
       WRITE(6,*) ' Starting from MOAOUT orbitals '
      END IF
*
*. Obtain, if required, supersymmetry of MO's
*
      IF(I_USE_SUPSYM.EQ.1) THEN
*. Supersymmetry of orbital in MOAOUT
         WRITE(6,*) ' Supersymmetry of orbitals in MOAOUT: '
         CALL SUPSYM_FROM_CMOAO(WORK(KMOAOUT),WORK(KISUPSYM_FOR_BAS),
     &                         WORK(KMO_ACT_SUPSYM))
*. Obtain reorder array going from correct order to actual order
         CALL REO_2SUPSYM_ORDERS(WORK(KMO_OCC_SUPSYM),
     &        WORK(KMO_ACT_SUPSYM),WORK(KIREO_INI_OCC))
*. Reorder to obtain the occ order of supersymmetry
         CALL REO_CMOAO(WORK(KMOAOUT),WORK(KMOAO_ACT),
     &        WORK(KIREO_INI_OCC),1,2)
*. Check that we now have correct supersymmetry (Jeppe has been messing up...)
         CALL SUPSYM_FROM_CMOAO(WORK(KMOAOUT),WORK(KISUPSYM_FOR_BAS),
     &                         WORK(KMO_ACT_SUPSYM))
         CALL ICOPVE(WORK(KMO_ACT_SUPSYM), WORK(KMO_SUPSYM), NTOOB)
         IDENT = IS_I1_EQ_I2(WORK(KMO_OCC_SUPSYM),
     &                       WORK(KMO_SUPSYM),NTOOB)
         IF(IDENT.EQ.0) THEN
           WRITE(6,*) ' Error: Reordered orbitals are not in occ order'
           WRITE(6,*) ' Obtained symmetry of reordered orbitals '
           CALL IWRTMA3(WORK(KMO_SUPSYM),1,NTOOB,1,NTOOB)
           WRITE(6,*) ' Required order '
           CALL IWRTMA3(WORK(KMO_OCC_SUPSYM),1,NTOOB,1,NTOOB)
           STOP ' Error: Jeppe is STILL messing supersymmetry up!!! '
         END IF



      END IF

      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Expansion of initial MOs in AOs '
        WRITE(6,*) ' ================================'
        CALL APRBLM2(WORK(KMOAOUT),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. Calculate metric over MO's in KLCMOAO2)..
      CALL GET_SMO(WORK(KMOAOUT),WORK(KLSAOE),0)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Metric in final initial orbitals ... '
        WRITE(6,*) ' ===================================='
        CALL APRBLM2(WORK(KLSAOE),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. Obtain CBIO: expansion of orbitals in MO's, CBIO2: expansion
*  of orbitals in AO's orbitals
      CALL GET_CBIO(WORK(KMOAOUT),WORK(KCBIO),WORK(KCBIO2))
*
* =======================================================================
* Bioorthogonal integral transformation with indices corresponding to 
* annihilation indices being in bioorthonormal basis
* =======================================================================
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Bioorthogonal integral transformation '
      END IF
      IE2LIST_A = IE2LIST_FULL_BIO
      IOCOBTP_A = 1
      INTSM_A = 1
      CALL PREPARE_2EI_LIST
*. Two forms 1: Operator acts on bioorthonormal expansion
*               creation operators are in bio, annihilation are in orig,
*               integral indices converse
*            2: Operator acts on origonal expansion
*. 
      I_STRINGS_BIO_OR_ORIG = 1
      IF(I_STRINGS_BIO_OR_ORIG.EQ.1) THEN
        KKCMO_I = KMOAOUT
        KKCMO_J = KCBIO2
        KKCMO_K = KMOAOUT
        KKCMO_L = KCBIO2
      ELSE
        KKCMO_I = KCBIO2
        KKCMO_J = KMOAOUT
        KKCMO_K = KCBIO2
        KKCMO_L = KMOAOUT
      END IF
C     DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
      CALL DO_ORBTRA(1,1,0,IE2LIST_FULL_BIO,IOCOBTP_A,INTSM_A)
      CALL FLAG_ACT_INTLIST(IE2LIST_FULL_BIO)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' one-electron integrals in biobase'
        WRITE(6,*) ' ================================='
        CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. Transfer the inactive Fock-matrix to feeder matrix for integral fetchers
      NINT1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1_F)
*
      IPRNTL = IPRDIA
      CALL NORTCALC(IREFSM,JCMBSPC,ICSPC_CNF,I_DO_NONORT_MCSCF,IPRNTL,
     &            E_FINAL,ERROR_NORM_FINAL,CONV_F,IVBGNSP_PREV)

      RETURN
      END
      SUBROUTINE NORTCALC(ISM,ISPC_GAS,ISPC_CNF,I_DO_NONORT_MCSCF,
     &                 IPRNT,
     &                 E_FINAL,VN_NORTCI,CONV_NORTCI,IVBGNSP_PREV)
*
* CI optimization in GAS space number ISPC for symmetry ISM              
*
* Information about the number of SD, CSF's is assumed to have
* been determined outside this routine
*
*
* Jeppe Olsen, June 2011
*  
*. Last modifications; Jeppe 2013; Analytic orbital Hessian and more
*
      INCLUDE 'wrkspc.inc'
      LOGICAL CONVER_NORTCI, CONVER_NORTMC
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
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'comjep.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'vb.inc'
*. Common block for communicating with sigma
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3, 
     &       LUCBIO_SAVE, LUHCBIO_SAVE,LUC_SAVE
*. Common block for transferring info to finite difference routines.
      COMMON/EVB_TRANS/KLIOOEXC_A, KLKAPPA_A,
     &                 KLIOOEXC_S,KLKAPPA_S,
     &                 KL_C,KL_VEC2,KL_VEC3,
     &                 KLOOEXC
*
      INCLUDE 'cecore.inc'
      COMMON/CMXCJ/MXCJ,MAXK1_MX,LSCMAX_MX
*
      COMMON/H_OCC_CONS/IH_OCC_CONS
*
      REAL*8 INPRDD, INPROD
*
      EXTERNAL SIGMA_NORTCI, PRECOND_NORTCI, E_VB_FROM_KAPPA_WRAP
*

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'NORTCI')
      NTEST = 10
      NTEST = MAX(NTEST,IPRNT)
*
      IF(NORT_MET.GE.3) THEN
         WRITE(6,*) ' Currently unknown NORT_MET = ', NORT_MET
         STOP ' Currently unknown NORT_MET '
       END IF
*
      IF(NTEST.GT.1) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ======================================='
        WRITE(6,*) ' Control has been transferred to NORTCI'
        WRITE(6,*) ' ======================================='
        WRITE(6,*) 
        IF(NORT_MET.EQ.1) THEN
          WRITE(6,'(5X,A)')
     &   ' Non-orthogonal wave function will be expanded in CI space'
       ELSE IF (NORT_MET.EQ.2) THEN
          WRITE(6,'(5X,A)')
     &   ' Initial suite of non-orthogonal configuration codes '
         ELSE
           WRITE(6,*) ' Currently unknown NORT_MET = ', NORT_MET
           STOP ' Currently unknown NORT_MET '
        END IF
*
        IF(I_DO_NONORT_MCSCF.EQ.1) THEN
          WRITE(6,*) ' I will also do MCSCF.... '
          WRITE(6,*) ' ========================='
        END IF
*
        WRITE(6,'(5X,A,I3)')
     &   ' Configuration space', ISPC_CNF
        WRITE(6,'(5X,A,I3)')
     &   ' Spanning CI-space:', ISPC_GAS
        WRITE(6,'(5X,A,I3)')
     & ' Orbital space containing non-orthogonal orbitals ',
     &   NORTCIX_SCVB_SPACE
        WRITE(6,'(5X,A,I3)')
     & ' Allowed excitation level from Spin-coupled valence space',
     &  NORTCI_SCVB_EXCIT
        WRITE(6,*) ' Orbital Min. occ Max. occ '
        WRITE(6,*) ' =========================='
        DO IORB = 1, NOBPT(NORTCIX_SCVB_SPACE)
          WRITE(6,'(3X,I4,2I3)')
     &    IORB, IOCC_MIN_GN(IORB,ISPC_CNF), IOCC_MAX_GN(IORB,ISPC_CNF)
        END DO
      END IF
*
*. Prepare for integral handling for complete array
*. (needed for codes where integrals are accessed individually)
      IE2ARRAY_A = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL_BIO))
      I12S_A = I12S_G(IE2ARRAY_A)
      I34S_A = I34S_G(IE2ARRAY_A)
      I1234S_A = I1234S_G(IE2ARRAY_A)
      IOCOBTP_A = IOCOBTP_G(IE2ARRAY_A)
      KINT2_LA = KINT2_A(IE2ARRAY_A)
      KPINT2_LA = KPINT2_A(IE2ARRAY_A)
*
*. Number of dets, csf's and configs in CI expansion
*
      NDET = NSD_PER_SYM_GN(ISM,ISPC_CNF)
      NCSF = NCSF_PER_SYM_GN(ISM,ISPC_CNF)
      NCONF = NCONF_PER_SYM_GN(ISM,ISPC_CNF)

      IF(IPRNT.GT.1) THEN
       WRITE(6,'(A,I9)') 
     & ' Number of determinants/combinations  ',NDET
       WRITE(6,'(A,I9)') 
     & ' Number of CSFs  ',NCSF
       WRITE(6,'(A,I9)') 
     & ' Number of Confs  ',NCONF
      END IF
*.Transfer to CANDS
      ICSM = ISM
      ISSM = ISM
*. Complete operator 
      I12 = 2
*
      ICSPC_CN = ICSPC_CNF
      ISSPC_CN = ISSPC_CNF
      IMSPC_CN = IMSPC_CNF
*
      IF(NORT_MET.EQ.1) THEN
*
*.Initial version with standard CI behind the scene
*
*. allocate memory for this
         ICSPC = ISPC_GAS
         ISSPC = ISPC_GAS
         IMSPC = ISPC_GAS
*
         WRITE(6,*) ' NORTCI: ICSPC_CNF,ISSPC_CNF,IMSPC_CNF',
     &                        ICSPC_CNF,ISSPC_CNF,IMSPC_CNF
         WRITE(6,*) ' NORTCI: ICSPC_CN,ISSPC_CN,IMSPC_CN',
     &                        ICSPC_CN,ISSPC_CN,IMSPC_CN
*
         CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)
         KVEC1P = KVEC1
         KVEC2P = KVEC2
      END IF ! if NORT_MET = 1
      IF(NORT_MET.EQ.2) THEN
* We will use a number of different spaces, vectors should be 
* able to store max space
      END IF
*
*
* Set up complete H and S for test
*
      I_DO_COMHAM = 0
      IF(I_DO_COMHAM .EQ. 1) THEN
        CALL COMHAM_HS_GEN(SIGMA_NORTCI,NCSF)
C            COMHAM_HS_GEN(MSTV,NDIM)
        STOP ' Enforced stop after COMHAM_HS_GEN '
      END IF
*
*. CI diagonal - if required
* 
*. Not yet implemented
*
      WRITE(6,*) ' Diagonal in Non-orthogonal CI not yet implemented'
*
      I_DO_PRECOND = 0
      IPREC_FORM = 0
      I_ER_CONV = 2
      THRES_R = SQRT(THRES_E)
      SHIFT = 0.0D0
*
      MAXITL = MAXIT
      MAXVECL = MXCIV
*
*. Allocate space for iterative solver: 
*. Four scratch vectors
C     CALL MEMMAN(KL_VEC1,NCSF,'ADDL  ',2,'EXTVC1')
C     CALL MEMMAN(KL_VEC2,NCSF,'ADDL  ',2,'EXTVC2')
C     CALL MEMMAN(KL_VEC3,NCSF,'ADDL  ',2,'EXTVC3')
*. Increased for CONF approach
      CALL MEMMAN(KL_VEC1,NCSF_MNMX_MAX,'ADDL  ',2,'EXTVC1')
      CALL MEMMAN(KL_VEC2,NCSF_MNMX_MAX,'ADDL  ',2,'EXTVC2')
      CALL MEMMAN(KL_VEC3,NCSF_MNMX_MAX,'ADDL  ',2,'EXTVC3')
*. Space for subsspace matrices
      CALL MEMMAN(KL_RNRM,MAXITL*NROOT,'ADDL  ',2,'RNRM  ')
      CALL MEMMAN(KL_EIG ,MAXITL*NROOT,'ADDL  ',2,'EIG   ')
      CALL MEMMAN(KL_FINEIG,NROOT,'ADDL  ',2,'FINEIG')
      CALL MEMMAN(KL_APROJ,MAXVECL**2,'ADDL  ',2,'APROJ ')
      CALL MEMMAN(KL_SPROJ,MAXVECL**2,'ADDL  ',2,'SPROJ ')
      CALL MEMMAN(KL_AVEC ,MAXVECL**2,'ADDL  ',2,'AVEC  ')
      LLWORK = 5*MAXVECL**2 + 2*MAXVECL
      CALL MEMMAN(KL_WORK ,LLWORK   ,'ADDL  ',2,'WORK  ')
      CALL MEMMAN(KL_AVEC ,MAXVECL**2,'ADDL  ',2,'AVECP ')
      CALL MEMMAN(KL_AVECP,MAXVECL**2,'ADDL  ',2,'AVECP ')
*. And a matrix over active orbitals
      CALL MEMMAN(KL_MACT,NACOB**2,'ADDL  ',2,'M_ACT ')
*
*. Initial approximation to CI-vector
*
      IF(IRESTR.EQ.0) THEN
C       INI_CSFEXP(CINI)
        CALL INI_CSFEXP(WORK(KL_VEC1))
*. and transfer initial guess to DISC
        CALL VEC_TO_DISC(WORK(KL_VEC1),NCSF,1,-1,LUSC54)
       ELSE
        WRITE(6,*) ' Restart from previous CI vector '
*. Expand from previous to current cnf-space
        NCSF_PREV = NVB_CSF(IVBGNSP_PREV)
        WRITE(6,*) 'IVBGNSP_PREV, NCSF_PREV = ',
     &              IVBGNSP_PREV, NCSF_PREV 
        CALL VEC_FROM_DISC(WORK(KL_VEC1),NCSF_PREV,1,-1,LUSC54)
C       EXP_CNFSPC(CIVECIN,CIVECUT,ICONF_OCC,NCONF_FOR_OPEN,
C    &           MINOCC_IN,MAXOC_IN,NOBCNF)
        IF(IVBGNSP_PREV.EQ.0) THEN
        CALL EXP_CNFSPC(WORK(KL_VEC1), WORK(KL_VEC2), 
     &       WORK(KICONF_OCC_GN(ICSM,ISPC_CNF)),
     &       NCONF_PER_OPEN_GN(1,ICSM,ISPC_CNF),
     &       VB_REFSPCO_MIN, VB_REFSPCO_MAX,
     &       NACOB)
        ELSE
        CALL EXP_CNFSPC(WORK(KL_VEC1), WORK(KL_VEC2), 
     &       WORK(KICONF_OCC_GN(ICSM,ISPC_CNF)),
     &       NCONF_PER_OPEN_GN(1,ICSM,ISPC_CNF),
     &       VB_GNSPC_MIN(1,IVBGNSP_PREV),VB_GNSPC_MAX(1,IVBGNSP_PREV),
     &       NACOB)
        END IF ! Test of VBGNSP_PREV
        CALL VEC_TO_DISC(WORK(KL_VEC2),NCSF,1,-1,LUSC54)
       END IF
*. And diagonalize
      NTESTL = 10000
      SHIFT = 0.0D0
      CALL MINGENEIG(SIGMA_NORTCI,PRECOND_NORTCI,
     &     IPREC_FORM,THRES_E,THRES_R,I_ER_CONV,
     &     WORK(KL_VEC1),WORK(KL_VEC2),WORK(KL_VEC3),
     &     LUSC54, LUSC37,
     &     WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &     NCSF,LUSC38,LUSC39,LUSC40,LUSC53,LUSC51,LUSC52,
     &     NROOT,MAXVECL,NROOT,WORK(KL_APROJ),
     &     WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &     NTESTL,SHIFT,WORK(KL_AVECP),I_DO_PRECOND,
     &     CONV_NORTCI,E_NORTCI,VN_NORTCI)
      E_FINAL = E_NORTCI
*
      IF(I_DO_NONORT_MCSCF.EQ.0) CONV_F = CONV_NORTCI
*
      WRITE(6,*) ' Final energy in non-orthogonal CI ', E_NORTCI
      WRITE(6,*) ' Final residual norm in non-orthogonal CI',
     &           VN_NORTCI
      IF(NTEST.GE.10000) THEN
       WRITE(6,*) ' Final approximation to eigenvector from MINGENEIG'
       CALL WRTMAT(WORK(KL_VEC1),1,NCSF,1,NCSF)
      END IF
*
* Analyze the CI- coefficients of the resulting wave function
*
C    ANACSF(CIVEC,ICONF_OCC,NCONF_FOR_OPEN,IPROCS,THRES,
C    &           MAXTRM,IOUT)
      MAXTRM = 1000
      THRES = 0.03
      IOUT = 6
*. The analyzer assumes full set of active electrons, adjust for this
      NACTEL = NACTEL - 2*(IB_ORB_CONF-NINOB-1)
      CALL ANACSF(WORK(KL_VEC1),WORK(KICONF_OCC_GN(ICSM,ISPC_CNF)), 
     &     NCONF_PER_OPEN_GN(1,ICSM,ISPC_CNF),
     &     WORK(KCFTP),THRES, MAXTRM,IOUT)
      NACTEL = NACTEL + 2*(IB_ORB_CONF-NINOB-1)
*. Density etc not implemented for NORT_MET = 2, so
      IF(NORT_MET.EQ.2) GOTO 9999
*
* And construct density matrix
*
      XDUM = 0.0D0
      CALL VB_DENSI(WORK(KRHO1),XDUM,1,WORK(KL_VEC1),WORK(KL_VEC2),
     &              WORK(KL_VEC3))
*. Obtain natural orbitals and natural occupation numbers
*. 1: Metric over active orbitals
C     SACT(SACT,C)
      CALL GET_SACT(WORK(KL_MACT),WORK(KMOAOUT))
*. 2: and diagonalize using metric of active orbitals
C     NONORT_NATORB(SACT,RHO1)
      CALL NONORT_NATORB(WORK(KL_MACT),WORK(KRHO1))
* 
      IF(I_DO_NONORT_MCSCF.EQ.1) THEN
        IREFSPC_MCSCF = ISPC_GAS
*
        IF(NORT_MET.NE.1) THEN
          WRITE(6,*) ' MCSCF works only for NORT_MET = 1'
          STOP       ' MCSCF works only for NORT_MET = 1'
        END IF
*. Allowed number of micro and macro's

        WRITE(6,*) ' MCSCF part entered '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' Allowed number of macroiterations ', MAXIT_MAC
        WRITE(6,*) ' Allowed number of microiterations ', MAXIT_MIC
*
*  ====================
*  Generate excitations
*  ====================
* 
* Two types:
* Interspace excitations: only antisymmtric conformal operators
* Active-Active exciations: both symmetric and antisymmetric operators
*
* For historical reasons, there is a flag for eliminating the 
* interspace excitations
*
*. Number of excitations
* ======================
*
*
        INCLUDE_ONLY_TOTSYM_SUPSYM = 1
        IF(I_USE_SUPSYM.EQ.1.AND.INCLUDE_ONLY_TOTSYM_SUPSYM.EQ.1) THEN
          I_RESTRICT_SUPSYM = 1
        ELSE
          I_RESTRICT_SUPSYM = 0
        END IF
        I_DO_INTER_EXC = 1
*
*. Number of internal excitations in active space
C       ORB_EXCIT_INT_SPACE(IORBSPC,ITOTSYM,NOOEXC,IOOEXC,NUMONLY)
        CALL ORB_EXCIT_INT_SPACE
     &  (NORTCIX_SCVB_SPACE,1,NOOEXC_AA,IDUM,1,1,
     &  I_RESTRICT_SUPSYM,WORK(KMO_SUPSYM))
*.Number of interspace excitations
        IF(I_DO_INTER_EXC.EQ.1) THEN
*. Number of interspace excitations
*. Nonredundant type-type excitations
          CALL MEMMAN(KLTTACT,(NGAS+2)**2,'ADDL  ',1,'TTACT ')
          CALL NONRED_TT_EXC(WORK(KLTTACT),IREFSPC_MCSCF,0)
*. Nonredundant interspace orbital excitations
          KLOOEXC = 1
          KLOOEXCC= 1
*. Number of interspace excitations
          CALL NONRED_OO_EXC2(NOOEXC_IS,WORK(KLOOEXC),WORK(KLOOEXCC),
     &         1,WORK(KLTTACT),I_RESTRICT_SUPSYM,WORK(KMO_SUPSYM),
     &         N_INTER_EXC,N_INTRA_EXC,1)
        END IF
*. Number of symmetric rotations
        NOOEXC_S = NOOEXC_AA
*. Number of antisymmetric rotations
        NOOEXC_A = NOOEXC_IS + NOOEXC_AA
*. The total number of excitations 
        NOOEXC =  NOOEXC_S + NOOEXC_A 
*
*. Allocate space
* ======================
*
*. Separate arrays are set up for all and symmetric excitations(??)
*. 
*. For all excitations
        CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
        CALL MEMMAN(KLOOEXCC,2*NOOEXC,'ADDL  ',1,'OOEXCC')
*. For the symmetric excitations
        CALL MEMMAN(KLOOEXCC_S,2*NOOEXC_S,'ADDL  ',1,'OOEXCS')
*. Allow these parameters to be known outside
        KIOOEXC = KLOOEXC
        KIOOEXCC = KLOOEXCC
*
*. And the excitations: The active- active are added twice..
* ======================
*
        IF(I_DO_INTER_EXC.EQ. 1) THEN
*. The interspace excitations 
          CALL NONRED_OO_EXC2(NOOEXC_IS,WORK(KLOOEXC),WORK(KLOOEXCC),
     &     1,WORK(KLTTACT),I_RESTRICT_SUPSYM,WORK(KMO_SUPSYM),
     &     N_INTER_EXC,N_INTRA_EXC,2)
        END IF
*. The internal excitations 
        CALL ORB_EXCIT_INT_SPACE(NORTCIX_SCVB_SPACE,1,NOOEXC_S,
     &       WORK(KLOOEXCC),0,NOOEXC_IS+1,
     &       I_RESTRICT_SUPSYM,WORK(KMO_SUPSYM))
        CALL ORB_EXCIT_INT_SPACE(NORTCIX_SCVB_SPACE,1,NOOEXC_S,
     &       WORK(KLOOEXCC),0,NOOEXC_IS+NOOEXC_S+1,
     &       I_RESTRICT_SUPSYM,WORK(KMO_SUPSYM))
*. Save also the internal excitations in KLOOEXCC_S
        CALL ICOPVE3(WORK(KLOOEXCC),NOOEXC_IS*2+1,
     &       WORK(KLOOEXCC_S),1,2*NOOEXC_S)
C ICOPVE3(IIN,IOFFIN,IOUT,IOFFOUT,NDIM)
        WRITE(6,*) ' NOOEXC after ORB_EXCIT.. ', NOOEXC
C PRINT_ORBEXC_LIST(IOOEXC,NOOEXC_A,NOOEXC_S) 
        WRITE(6,*) ' The list of orbital excitations'
        CALL  PRINT_ORBEXC_LIST(WORK(KLOOEXCC),NOOEXC_A,NOOEXC_S)
        WRITE(6,*) ' The list of symmetric excitations'
        CALL  PRINT_ORBEXC_LIST(WORK(KLOOEXCC_S),0,NOOEXC_S)
*
* Allocate space for gradient, kappa, Hessian, etc
* ================================================
*
        WRITE(6,*) ' NOOEXC before MEMMAN ', NOOEXC
        CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1_MC ')
        CALL MEMMAN(KLKAPPA,NOOEXC,'ADDL  ',2,'LKAPPA')
        CALL MEMMAN(KLE2SC,NOOEXC,'ADDL  ',2,'E2SC  ')
*. Memory for orbital-Hessian - if  required
        LE2 = NOOEXC*(NOOEXC+1)/2
        CALL MEMMAN(KLE2,LE2,'ADDL  ',2,'E2P_MC')
*. For eigenvectors of orbhessian
        LE2F = NOOEXC**2
        CALL MEMMAN(KLE2F,LE2F,'ADDL  ',2,'E2_MC ')
*. and eigenvalues, scratch, kappa
        CALL MEMMAN(KLE2VL,NOOEXC,'ADDL  ',2,'EIGVAL')
*. Space for two one-bodydensity matrices
        CALL MEMMAN(KLS,NTOOB**2,'ADDL  ',2,'SMO   ')
*. KMOAOIN will be used for storing MO's that should be transformed
        I_STRINGS_BIO_OR_ORIG = 1
        IF(I_STRINGS_BIO_OR_ORIG.EQ.1) THEN
          KKCMO_I = KMOAOIN
          KKCMO_J = KCBIO2
          KKCMO_K = KMOAOIN
          KKCMO_L = KCBIO2
        ELSE
          KKCMO_I = KCBIO2
          KKCMO_J = KMOAOIN
          KKCMO_K = KCBIO2
          KKCMO_L = KMOAOIN   
        END IF
*
        LEN1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
        LEN1_A = NDIM_1EL_MAT(1,NACOBS,NACOBS,NSMOB,0)
        CALL COPVEC(WORK(KMOAOUT),WORK(KMOAOIN),LEN1_F)
*. For summary
        NITEM = 6
        CALL MEMMAN(KL_SUMMARY,NITEM*MAXIT_MAC,'ADDL  ',2,'SUMMAR')
*
*. Finished with the preparations
*
        CONVER_NORTMC = .FALSE.
        XKAP_THRES = 1.0D-6
        STEP_MAX = 0.75D0
*
        DO IMAC = 1, MAXIT_MAC
          NMAC = IMAC
          WRITE(6,*) ' Output from Macroiteration', IMAC
          WRITE(6,*) ' =================================='
          INIMIC = 1
          DO IMIC = 1, MAXIT_MIC
            WRITE(6,*) ' Output from Microiteration', IMIC
            WRITE(6,*) ' =================================='
*
*. The current expansion of the AOs is in KMOAOIN. Obtain the
* bioorbitals
            CALL GET_CBIO(WORK(KMOAOIN),WORK(KCBIO),WORK(KCBIO2))
            IF(NTEST.GE.100) THEN
               WRITE(6,*) 
     &         ' Current set of orbitals'
               CALL APRBLM2(WORK(KMOAOIN),NTOOBS,NTOOBS,NSMOB,0)
            END IF
            IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' Current set of bioorbitals '
               CALL APRBLM2(WORK(KCBIO),NTOOBS,NTOOBS,NSMOB,0)
            END IF
*
            IF(NTEST.GE.100) THEN
*             Calculate and print metric 
              CALL GET_SMO(WORK(KMOAOIN),WORK(KLS),0)
              WRITE(6,*) ' Metric in Current MO basis '
              CALL APRBLM2(WORK(KLS),NTOOBS,NTOOBS,NSMOB,0)
            END IF
*
* =====================================================
* Integral transformation to current set of orbitals
* =====================================================
*
            IF(NTEST.GE.10) THEN
              WRITE(6,*) 
     &        ' Bioorthogonal integral transformation '
            END IF
            IOCOBTP_A = 1
            INTSM_A = 1
            CALL FLAG_ACT_INTLIST(IE2LIST_FULL_BIO)
            CALL DO_ORBTRA(1,1,0,IE2LIST_FULL_BIO,
     &                     IOCOBTP_A,INTSM_A)
*

*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' one-electron integrals in biobase'
              WRITE(6,*) ' ================================='
              CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,0)
            END IF
            NINT1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
            CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1_F)
*
* ==============================
*. Perform CI in current basis in first Micro of each macro
* ==============================
*
            IF(IMIC.EQ.1) THEN
             IF(NTEST.GE.1000) WRITE(6,*) ' CI optimization '
             CALL MINGENEIG(SIGMA_NORTCI,PRECOND_NORTCI,
     &            IPREC_FORM,THRES_E,THRES_R,I_ER_CONV,
     &            WORK(KL_VEC1),WORK(KL_VEC2),WORK(KL_VEC3),
     &            LUSC54, LUSC37,
     &            WORK(KL_RNRM),WORK(KL_EIG),WORK(KL_FINEIG),MAXITL,
     &            NCSF,LUSC38,LUSC39,LUSC40,LUSC53,LUSC51,LUSC52,
     &            NROOT,MAXVECL,NROOT,WORK(KL_APROJ),
     &            WORK(KL_AVEC),WORK(KL_SPROJ),WORK(KL_WORK),
     &            NTESTL,SHIFT,WORK(KL_AVECP),I_DO_PRECOND,
     &            CONV_NORTCI,E_NORTCI,VN_NORTCI)
             WORK(KL_SUMMARY-1+(IMAC-1)*NITEM + 1) = E_NORTCI
             IF(IMAC.EQ.1) THEN
               DELTA = 0.0D0
             ELSE 
               DELTA = E_NORTCI - WORK(KL_SUMMARY-1+(IMAC-2)*NITEM + 1)
             END IF
             WORK(KL_SUMMARY-1+(IMAC-1)*NITEM + 2) = DELTA
*. Preliminary final energy
             E_FINAL = E_NORTCI
*
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) 
     &          ' Final approximation to eigenvector from MINGENEIG'
               CALL WRTMAT(WORK(KL_VEC1),1,NCSF,1,NCSF)
             END IF
*
* And construct the one- and two-body density matrices
*
             CALL VB_DENSI(WORK(KRHO1),WORK(KRHO2),2,WORK(KL_VEC1),
     &            WORK(KL_VEC2),WORK(KL_VEC3))
*. Construct Active Fock-matrix
             CALL DO_ORBTRA(1,1,1,IE2LIST_FULL_BIO,
     &                     IOCOBTP_A,INTSM_A)
            END IF !micit = 1
*
* =====================================
*. Construct Gradient at current point
* =====================================
*
*. The Fock-matrices for biorthogonal expansion
C FOCK_MAT_NORT(F1,F2,I12,FI,FA)
            CALL FOCK_MAT_NORT(WORK(KF),WORK(KF2),2,WORK(KFI),WORK(KFA))
*. And the interspace gradient
C     E1_FROM_F_NORT(E1,F1,F2,IOPSM,IOOEXC,IOOEXCC,
C    &           NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
            CALL E1_FROM_F_NORT(WORK(KLE1),WORK(KF),WORK(KF2),1,
     &           WORK(KLOOEXC),WORK(KLOOEXCC),NOOEXC_A,NTOOB,
     &           NTOOBS,NSMOB,IBSO,IREOST)
*. And add the active-active gradient
* The interspace excitations
C           VB_GRAD_ORBVBSPC(NOOEXC,IOOEXC,E1,C,VEC1_CSF,VEC2_CSF)
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' Active-active gradient will be calculated '
            CALL VB_GRAD_ORBVBSPC(NOOEXC_S,WORK(KLOOEXCC_S),
     &      WORK(KLE1+NOOEXC_IS),
     &      WORK(KL_VEC1),WORK(KL_VEC2),WORK(KL_VEC3))
            IF(NTEST.GE.0) WRITE(6,*) ' Gradient calculated '
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Gradient vector '
              CALL WRTMAT(WORK(KLE1),1,NOOEXC,1,NOOEXC)
            END IF
            E1NORM = INPROD(WORK(KLE1),WORK(KLE1),NOOEXC)
            WORK(KL_SUMMARY-1+(IMAC-1)*NITEM + 3) = E1NORM
C           STOP ' Enforced stop '
*
            I_DO_DIFTEST = 0
            IF(I_DO_DIFTEST.EQ.1) THEN
*. Finite difference test of gradient at Kappa = 0
              KLIOOEXC_A = KIOOEXCC
              KLIOOEXC_S = KLOOEXCC_S
              KLKAPPA_A = KLKAPPA
              KLKAPPA_S = KLKAPPA+NOOEXC_A
              KL_C = KL_VEC1
*
C     COMMON/EVB_TRANS/KLIOOEXC_A, KLKAPPA_A,
C    &                 KLIOOEXC_S,KLKAPPA_S,
C    &                 KL_C,KL_VEC2,KL_VEC3
              ZERO = 0.0D0
              CALL SETVEC(WORK(KLKAPPA),ZERO,NOOEXC)
              CALL MEMMAN(KLE1_EXTRA,NOOEXC,'ADDL  ',2,'E1_EXT')
              CALL GENERIC_GRAD_FROM_F(WORK(KLE1_EXTRA),NOOEXC,
     &            E_VB_FROM_KAPPA_WRAP, WORK(KLKAPPA))
*
*. Clean up: recalculate integrals corresponding to MO's in MOAOIN
*
              IF(I_STRINGS_BIO_OR_ORIG.EQ.1) THEN
                KKCMO_I = KMOAOIN
                KKCMO_J = KCBIO2
                KKCMO_K = KMOAOIN
                KKCMO_L = KCBIO2
              ELSE
                KKCMO_I = KCBIO2
                KKCMO_J = KMOAOIN
                KKCMO_K = KCBIO2
                KKCMO_L = KMOAOIN   
              END IF
              IF(NTEST.GE.10) THEN
                WRITE(6,*) 
     &          ' Bioorthogonal integral transformation '
              END IF
              CALL TRAINT
              CALL FLAG_ACT_INTLIST(IE2LIST_FULL_BIO)
              STOP ' Enforced stop after FD calc of gradient'
            END IF ! Finite Difference test
*
            IF(IMIC.EQ.1) THEN
*
*. Obtain new orbital Hessian
*
*. Complete orbital Hessian
              IE2FORM = 1
*. Prepare for finite difference calc of Hessian
              KLIOOEXC_A = KIOOEXCC
              KLIOOEXC_S = KLOOEXCC_S
              KLKAPPA_A = KLKAPPA
              KLKAPPA_S = KLKAPPA+NOOEXC_A
              KL_C = KL_VEC1
*
*. IE2FORM Is not active at the moment
              IE2FORM = 1
              CALL ORBHES_VB(WORK(KLE2),IE2FORM)
*
*. Diagonalize to determine lowest eigenvalue
*
*. Outpack to complete form
              CALL TRIPAK(WORK(KLE2F),WORK(KLE2),2,NOOEXC,NOOEXC)
C                  TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
*. Lowest eigenvalue
C             DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
              CALL DIAG_SYMMAT_EISPACK(WORK(KLE2F),WORK(KLE2VL),
     &             WORK(KLE2SC),NOOEXC,IRETURN)
              IF(IRETURN.NE.0) THEN
                WRITE(6,*)
     &          ' Problem with diagonalizing E2, IRETURN =  ', IRETURN
              END IF
              IF(IPRNT.GE.1000) THEN
                WRITE(6,*) ' Eigenvalues: '
                CALL WRTMAT(WORK(KLE2VL),1,NOOEXC,1,NOOEXC)
              END IF
              IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' Eigenvectors of Hessian '
               CALL WRTMAT(WORK(KLE2F),NOOEXC,NOOEXC,
     &                     NOOEXC,NOOEXC)
              END IF
*. Lowest eigenvalue
C                       XMNMX(VEC,NDIM,MINMAX)
              E2VL_MN = XMNMX(WORK(KLE2VL),NOOEXC,1)
              IF(NTEST.GE.2)
     &        WRITE(6,*) ' Lowest eigenvalue of E2(orb) = ', E2VL_MN
            END IF !imic = 1
*. Transform gradient to diagonal basis
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Gradient in original basis before MATVCC'
              CALL WRTMAT(WORK(KLE1),1,NOOEXC,1,NOOEXC)
            END IF
            CALL MATVCC(WORK(KLE2F),WORK(KLE1),WORK(KLE2SC),
     &           NOOEXC,NOOEXC,1)
            CALL COPVEC(WORK(KLE2SC),WORK(KLE1),NOOEXC)
*. Solve shifted NR equations with step control
  666       CONTINUE
            TOLER = 1.1D0
*           SOLVE_SHFT_NR_IN_DIAG_BASIS(
*    &              E1,E2,NDIM,STEP_MAX,TOLERANCE,X,ALPHA)A
            CALL SOLVE_SHFT_NR_IN_DIAG_BASIS(WORK(KLE1),WORK(KLE2VL),
     &           NOOEXC,STEP_MAX,TOLER,WORK(KLKAPPA),ALPHA,DELTA_E_PRED)
            XNORM_STEP = 
     &      SQRT(INPROD(WORK(KLKAPPA),WORK(KLKAPPA),NOOEXC))
            WORK(KL_SUMMARY-1+(IMAC-1)*NITEM + 4) = XNORM_STEP
            IF(NTEST.GE.2) WRITE(6,*) ' Norm of step = ', XNORM_STEP
*. transform step to original basis
            CALL MATVCC(WORK(KLE2F),WORK(KLKAPPA),WORK(KLE2SC),
     &           NOOEXC,NOOEXC,0)
            CALL COPVEC(WORK(KLE2SC),WORK(KLKAPPA),NOOEXC)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Kappa in original basis '
              CALL WRTMAT(WORK(KLKAPPA),1,NOOEXC,1,NOOEXC)
            END IF
*. Energy for new step:
            ENEW =  E_VB_FROM_KAPPA_WRAP(WORK(KLKAPPA))
            WRITE(6,*) ' Energy at new iteration point ', ENEW
*. Preliminary E_FINAL ..
            E_FINAL = ENEW
            THRESD = 1.0D-7
            IF(ENEW.GT.E_NORTCI+THRES) THEN
*. Step was to large, Decrease steplength and recalculate step
               STEP_MAX = STEP_MAX/2.0D0
               GOTO 666
            END IF

*
*. Obtain New MO coefficients in MOAOUT: MOAOIN* Exp(-Kappa_A S) Exp(-Kappa_S S)
*
C     NEWMO_FROM_KAPPA_NORT(
C    &           NOOEXC_A,IOOEXC_A,KAPPA_A,
C    &           NOOEXC_S,IOOEXC_S,KAPPA_S,CMOAO_IN,CMOAO_OUT)
            CALL NEWMO_FROM_KAPPA_NORT(
     &           NOOEXC_A,WORK(KIOOEXCC),WORK(KLKAPPA),
     &           NOOEXC_S,WORK(KLIOOEXC_S),WORK(KLKAPPA+NOOEXC_A),
     &           WORK(KMOAOIN),WORK(KMOAOUT))
*. And copy to KMOAOIN for next round
            CALL COPVEC(WORK(KMOAOUT),WORK(KMOAOIN),LEN1_F)
*
            IF(IPRNT.GE.100) THEN
              WRITE(6,*) ' Updated MOAO-coefficients'
              CALL APRBLM2(WORK(KMOAOUT),NTOOBS,NTOOBS,NSMOB,0)
            END IF
         IF(XNORM_STEP.LT.XKAP_THRES) THEN
           CONVER_NORTMC = .TRUE.
           GOTO 1001
         END IF
*
         END DO ! End of loop over microiterations
       END DO ! End of loop over macroiterations
 1001 CONTINUE
*
      IF(CONVER_NORTMC) THEN
        WRITE(6,*) ' Convergence was obtained in ', NMAC , ' iterations'
      ELSE 
        WRITE(6,*) ' Convergence was not obtained '
      END IF
*
      WRITE(6,*) ' Final energy = ', E_FINAL
*
      IF(IPRNT.GE.2) THEN
        WRITE(6,*) ' Optimized MOAO-coefficients:'
        WRITE(6,*) ' ============================'
        CALL PRINT_CMOAO(WORK(KMOAOUT))
      END IF
*
* And construct the final density matrices
*
      CALL VB_DENSI(WORK(KRHO1),WORK(KRHO2),1,WORK(KL_VEC1),
     &     WORK(KL_VEC2),WORK(KL_VEC3))
*. Obtain natural orbitals and natural occupation numbers
*. 1: Metric over active orbitals
C     SACT(SACT,C)
      CALL GET_SACT(WORK(KL_MACT),WORK(KMOAOUT))
*. 2: and diagonalize using metric of active orbitals
C     NONORT_NATORB(SACT,RHO1)
      CALL NONORT_NATORB(WORK(KL_MACT),WORK(KRHO1))
*
* Analyze the CI- coefficients of the resulting wave function
*
C    ANACSF(CIVEC,ICONF_OCC,NCONF_FOR_OPEN,IPROCS,THRES,
C    &           MAXTRM,IOUT)
      MAXTRM = 1000
      THRES = 0.03
      IOUT = 6
      NACTEL = NACTEL - 2*(IB_ORB_CONF-NINOB-1)
      CALL ANACSF(WORK(KL_VEC1),WORK(KICONF_OCC_GN(ICSM,ISPC_CNF)), 
     &     NCONF_PER_OPEN_GN(1,ICSM,ISPC_CNF),
     &     WORK(KCFTP),THRES, MAXTRM,IOUT)
      NACTEL = NACTEL + 2*(IB_ORB_CONF-NINOB-1)
*
      WRITE(6,*) ' ======================='
      WRITE(6,*) ' Summary of convergence '
      WRITE(6,*) ' ======================='
      WRITE(6,*)
      WRITE(6,*) 
     &' Iter            Energy   Delta E   E1-norm    Step '
      WRITE(6,*) 
     &' ===================================================='
      DO IMAC = 1, NMAC
        II = KL_SUMMARY + (IMAC-1)*NITEM-1
        WRITE(6,'(2X, I3, 1X, F18.10,1X, E10.3, E10.3, E10.3)')
     &  IMAC, WORK(II+1),WORK(II+2),WORK(II+3), WORK(II+4)
      END DO
*
      END IF ! I do MCSCF
*
 9999 CONTINUE
*
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUM,'NORTCI')
      RETURN
      END
      SUBROUTINE SIGMA_NORTCI(C,HC,SC,IDOHC,IDOSC)
*
* Routine for sigma-generation, nonorthogonal CI, using biortogonal 
* approach. Integrals in biobasis assumed in place
*
* Initial version, Jeppe Olsen June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'vb.inc'
*. Two local scratch files
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3, 
     &       LUCBIO_SAVE, LUHCBIO_SAVE,LUC_SAVE
*. Input:  C in CSF basis, configuration space ICSPC_CN
*. Output: Sigma CSF basis, configuration space  ISSPC_CN
*. Output files
      DIMENSION HC(*), SC(*)
*
      DIMENSION C(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'NOCISI')
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from SIGMA_NORTCI'
        WRITE(6,*) ' ========================'
        WRITE(6,*)
        IF(IDOHC.EQ.1) WRITE(6,*) ' HC will be calculated '
        IF(IDOSC.EQ.1) WRITE(6,*) ' SC will be calculated'
        WRITE(6,*)
        WRITE(6,*) ' CI and MINMAX space for C ', ICSPC,ICSPC_CN
        WRITE(6,*) ' CI and MINMAX space for S ', ISSPC,ISSPC_CN
        WRITE(6,*) ' CI and MINMAX space, Intermediate ', IMSPC,IMSPC_CN
*
        IF(NORT_MET.EQ.1) THEN
          WRITE(6,*) 'Approach based on reexpansion in GASpace '
        ELSE IF (NORT_MET.EQ.2) THEN 
          WRITE(6,*) ' Approach using initial configuration routines'
        END IF
*
      END IF
*. This routine does all the CSF-SD transformation explicitly, 
*. so fool inner routines, especially MV7, to believe that we 
*. are working with SD's
      NOCSF_SAVE = NOCSF
      NOCSF = 1
*
      NCSF_CSPC_CNF = NCSF_PER_SYM_GN(ICSM,ICSPC_CN)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input vector in CSF basis '
        CALL WRTMAT(C,1,NCSF_CSPC_CNF,1,NCSF_CSPC_CNF)
      END IF
*
      IF(NORT_MET.EQ.1) THEN
*
* Initial version, with standard FCI behind the screen (and you are pt 
* behind the screen...)
*  The route
*
* 1) Obtain Input state in biortogonal basis in space IMSPC_CN - 
*    in Slater determinants
* 2) Obtain Hamiltonian times input state in bioorthogonal basis
* 
* In the initial version step 1 consists of 
* 1.1) Transform C from CSF to SD in CI ICSPC_CN
* 1.2) Expand C in SD from from SPACE ICSPC_CN to space ICSPC
* 1.3) Calculate biortogonal C-vector
*      in CI space IMSPC 
* 1.4) Contract bioorthogonal C-vector from space IMSPC to IMSPC_CN
* whereas step 2 consists of 
* 2.1) Expand Bioorthogonal C-vector from space IMSPC_CN to IMSPC
* 2.2) calculate biorthogonal sigma-vector in space ISSPC
* 2.2) Contract bioothogonal sigma-vector to space ISSPC_CN
* 2.3) Transform sigma-vector from SD to CSF-basis
*
* The SC-vector is obtained from step 1.4
*
      ICSPC_ORIG = ICSPC
      ISSPC_ORIG = ISSPC
      IMSPC_ORIG = IMSPC
      ICSM_ORIG = ICSM
      ISSM_ORIG = ISSM
*
*. The CI space are actually assumed to be identical 

* =========
*. Step 1
* =========
*
* 1.1) Transform C from CSF to SD in CI ICSPC_CN
*
* Allocate space for output CI vector in SD basis in config basis
*
      NSD_CSPC_CNF = NSD_PER_SYM_GN(ICSM,ICSPC_CN)
*
      NSD =  NSD_CSPC_CNF
      NCSF = NCSF_CSPC_CNF
*
      CALL MEMMAN(KLC_SD,NSD_CSPC_CNF,'ADDL  ',2,'C_SD  ')
      CALL MEMMAN(KLVCI ,NSD_CSPC_CNF,'ADDL  ',2,'VCI   ')
*. Expand Input C vector from CSD to SD form
      CALL COPVEC(C,WORK(KLVCI),NCSF_CSPC_CNF)
C          CSDTVCM(CSFVEC,DETVEC,IWAY,ICOPY,ISYM,ICSPC,IMAXMIN_OR_GAS)
      XDUM = 0.0D0
      CALL CSDTVCM(WORK(KLVCI),WORK(KLC_SD),XDUM,1,0,ICSM,ICSPC_CN,1)
*
* 1.2) Expand C in SD from SPACE ICSPC_CN to space ICSPC
*
*.Obtain number and length of blocks of expansion
      CALL MEMMAN(KLBLK,MXNTTS,'ADDL  ',1,'LBLKCI')
C     LBLOCK_FOR_CIXP(LBLOCK,NBLOCK,ICISPC,ISYM)
      CALL LBLOCK_FOR_CIXP(WORK(KLBLK),NBLOCK_C,ICSPC,ICSM)
C     SCA_VEC_TO_BLKV_DISC(VEC,ISCA,NELMNT,LUOUT,NBLOCK,LBLOCK,VECBLK,IREW)  
      IF(NTEST.GE.10) WRITE(6,*) ' I will CALL SCA_VEC ... '
      CALL SCA_VEC_TO_BLKV_DISC(WORK(KLC_SD),
     &     WORK(KSDREO_I_GN(ICSM,ICSPC_CN)),
     &     NSD_CSPC_CNF,LUSCR1,NBLOCK_C,WORK(KLBLK),WORK(KVEC1P),1)
      IF(NTEST.GE.10) WRITE(6,*) ' Home from SCAVEC .... '
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input C-vector in GAS space form '
        CALL WRTVCD(WORK(KVEC1P),LUSCR1,1,-1)
      END IF
*
* 1.3) Calculate biortogonal C-vector in CI space IMSPC 
*
*. Expand CI vector to space IMSPC
      CALL EXPCIV(ICSM,ICSPC,LUSCR1,IMSPC,LUSCR2,-1,LUSCR3,1,1,IDC,
     &     NTEST)
C          EXPCIV(ISM,ISPCIN,LUIN,ISPCUT,LUUT,LBLK,LUSCR,NROOT,ICOPY,IDC,NTESTG)
      IF(NTEST.GE.10) WRITE(6,*) ' Back from EXPCIV I'
*. And then do the transformation defined by KCBIO
*. Save one-electron integrals
      IF(IH1FORM.EQ.1) THEN
        IPACK_H1 = 1
      ELSE
        IPACK_H1 = 0
      END IF
C NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN_H1 =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,IPACK_H1)
      CALL MEMMAN(KLH1SAVE,NTOOB**2,'ADDL  ',2,'H1SAVE')
      CALL COPVEC(WORK(KINT1),WORK(KLH1SAVE),LEN_H1)
*. Ecore is now adays included in MV7( called in TRACI) --hide it
      ECORE_SAVE = ECORE
      ECORE = 0.0D0
*
      IF(LUC_SAVE.NE.0) THEN
C?     WRITE(6,*) ' C in orig base will be saved in unit ', LUC_SAVE
       CALL COPVCD(LUSCR1,LUC_SAVE,WORK(KVEC1P),1,-1)
      END IF
      CALL REWINO(LUSCR2)
*. biotransform and save result in LUSCR2
      IF(NTEST.GE.1000) WRITE(6,*) ' Traci will be called '
      CALL TRACI(WORK(KCBIO),LUSCR1,LUSCR2,IMSPC,ICSM,
     &           WORK(KVEC1P),WORK(KVEC2P))
      IF(NTEST.GE.1000) WRITE(6,*) ' Home from Traci '
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' C in biobase, SD expansion '
        CALL WRTVCD(WORK(KVEC1P),LUSCR2,1,-1)
      END IF
C TRACI(X,LUCIN,LUCOUT,IXSPC,IXSM,VEC1,VEC2)
      ECORE = ECORE_SAVE
      CALL COPVEC(WORK(KLH1SAVE),WORK(KINT1),LEN_H1)
      IF(NTEST.GE.1000) WRITE(6,*) ' Back from TRACI '
      CALL LBLOCK_FOR_CIXP(WORK(KLBLK),NBLOCK_M,IMSPC,ICSM)
*
      IF(LUCBIO_SAVE.NE.0) THEN
C?     WRITE(6,*) ' C in biobase will be saved in unit ', LUCBIO_SAVE
       CALL COPVCD(LUSCR2,LUCBIO_SAVE,WORK(KVEC1P),1,-1)
      END IF
      IF(IDOSC.EQ.1) THEN
* Obtain the metric vector = <i!0> in space ICSPC_CNF
        CALL GAT_VEC_FROM_BLKV_DISC(WORK(KLC_SD),
     &       WORK(KSDREO_I_GN(ICSM,ICSPC_CN)),
     &       NSD_CSPC_CNF,LUSCR2,NBLOCK_M,WORK(KLBLK),WORK(KVEC1P),1)
      IF(NTEST.GE.10) WRITE(6,*) ' Back from GAT_VEC '
*. And transform to CSF basis
C       CSDTVCM(CSFVEC,DETVEC,IWAY,ICOPY,ISYM,ISPC,IMAXMIN_OR_GAS)
        XDUM = 0.0D0
        CALL CSDTVCM(WORK(KLVCI),WORK(KLC_SD),XDUM,2,0,ICSM,ICSPC_CN,1)
        CALL COPVEC(WORK(KLVCI),SC,NCSF)
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Back from CSDTVCM ' 
        END IF
      END IF
      IF(IDOHC.EQ.1) THEN
*. Obtain the transformed vector for the determinants of space IMSPC_CN
        CALL LBLOCK_FOR_CIXP(WORK(KLBLK),NBLOCK_M,IMSPC,ICSM)
        NSD_MSPC_CNF = NSD_PER_SYM_GN(ICSM,IMSPC_CN)
        CALL MEMMAN(KLCM_SD,NSD_MSPC_CNF,'ADDL  ',2,'CM_SD ')
        CALL GAT_VEC_FROM_BLKV_DISC(WORK(KLCM_SD),
     &       WORK(KSDREO_I_GN(ICSM,IMSPC_CN)),
     &       NSD_MSPC_CNF,LUSCR2,NBLOCK_M,WORK(KLBLK),WORK(KVEC1P),1)
C       GAT_VEC_FRM_BLKV_DISC(VEC,ISCA,NELMNT,LUIN,NBLOCK,LBLOCK,VECBLK,IREW)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Biotransformed C in Intermediate CN space '
          CALL WRTMAT(WORK(KLCM_SD),1,NSD_MSPC_CNF,1,NSD_MSPC_CNF)
          WRITE(6,*) ' SIGMA_NORTCI speaking, end of step 1: '
        END IF
      END IF
*
* ======
* Step 2
* ======
*
      IF(IDOHC.EQ.1) THEN
* 2.1) Expand Bioorthogonal C-vector from space IMSPC_CN to IMSPC
        CALL SCA_VEC_TO_BLKV_DISC(WORK(KLCM_SD),
     &       WORK(KSDREO_I_GN(ICSM,IMSPC_CN)),
     &       NSD_MSPC_CNF,LUSCR1,NBLOCK_M,WORK(KLBLK),WORK(KVEC1P),1)
* 2.2) calculate biorthogonal sigma-vector in space ISSPC
*
       ICSPC = IMSPC_ORIG
       ISSPC = ISSPC_ORIG
       I12 = 2
       XDUM = 3006.1956D0
       CALL MV7(WORK(KVEC1P),WORK(KVEC2P),LUSCR1,LUSCR2,XDUM,XDUM)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' HC in Biobase '
         CALL WRTVCD(WORK(KVEC1P),LUSCR2,1,-1)
       END IF
       IF(LUHCBIO_SAVE.NE.0) THEN
C?      WRITE(6,*) ' HC in biobase will be saved in unit ', LUHCBIO_SAVE
        CALL COPVCD(LUSCR2,LUHCBIO_SAVE,WORK(KVEC1P),1,-1)
       END IF
* 2.2) Contract biothogonal sigma-vector to space ISSPC_CN
        CALL LBLOCK_FOR_CIXP(WORK(KLBLK),NBLOCK_S,ISSPC,ISSM)
        NSD_SSPC_CNF = NSD_PER_SYM_GN(ISSM,ISSPC_CN)
        CALL MEMMAN(KLSS_SD,NSD_SSPC_CNF,'ADDL  ',2,'SS_SD ')
        CALL GAT_VEC_FROM_BLKV_DISC(WORK(KLSS_SD),
     &       WORK(KSDREO_I_GN(ISSM,ISSPC_CN)),
     &       NSD_SSPC_CNF,LUSCR2,NBLOCK_S,WORK(KLBLK),WORK(KVEC1P),1)
* 2.3) Transform sigma-vector from SD to CSF-basis
        XDUM = 0.0D0
        CALL CSDTVCM(WORK(KLVCI),WORK(KLSS_SD),XDUM,2,0,ISSM,ISSPC_CN,1)
C       CSDTVCM(CSFVEC,DETVEC,IWAY,ICOPY,ISYM,ISPC,IMAXMIN_OR_GAS)
        CALL COPVEC(WORK(KLVCI),HC,NCSF)
      END IF
*  
      ELSE IF (NORT_MET.EQ.2) THEN
*
* 1:. Perform Bioorthogonal transformation of C from space ISPC_CN to 
*  space IMSPC_CN. It is required that the spaces for for the individual steps
*  have been defined in IORBTRA_SPC_IN, IORBTRA_SPC_OUT
*. Pt evrything is in CORE
*
*
        ICISTR = 1
        IF(NTEST.GE.10) 
     &  WRITE(6,*) 
     &  ' TRACI_CONF will be called to perform orbital transformation'
        LUC = 0
        LUS = 0
*. TRACI_CONF modifies input vector, so
        NCSF_C = NCSF_PER_SYM_GN(ICSM,ICSPC_CN)
        NCSF_S = NCSF_PER_SYM_GN(ISSM,ISSPC_CN)
        NCSF_CS = MAX(NCSF_C,NCSF_S)
        CALL MEMMAN(KLC_CSF,NCSF_MNMX_MAX,'ADDL  ',2,'C_CSF ')
        WRITE(6,*) ' TEST: ICSM, ICSPC_CN, NCSF_C = ', 
     &                     ICSM, ICSPC_CN, NCSF_C
*.
*. Obtained transformed vector in C_CSF
        CALL COPVEC(C,SC,NCSF_C)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Input C vector (CSF basis) '
          CALL WRTMAT(C, NCSCF_C, 1, NCSCF_C, 1)
          WRITE(6,*)
        END IF
        CALL TRACI_CONF(SC,WORK(KLC_CSF),LUC,LUS)
C            TRACI_CONF(C,S,LUC,LUHC)
        IF(NTEST.GE.1000) THEN
          NCSF_M =  NCSF_PER_SYM_GN(ISSM,IMSPC_CNF)
          WRITE(6,*) ' IMSPC_CNF, NCSF_M = ', IMSPC_CNF, NCSF_M
          CALL WRTMAT(WORK(KLC_CSF),1,NCSF_M,1,NCSF_M)
        END IF
*. Extract metric in initial space
        IF(NTEST.GE.10)  WRITE(6,*) 
     &  ' REF_CNFVEC will be called to get metric times initial vector'
C            REF_CNFVEC(VECIN,ISPCIN,VECOUT,ISPCOUT,ISYM)
        CALL REF_CNFVEC(WORK(KLC_CSF),IMSPC_CN,SC,ICSPC_CN,ICSM)
        IF(NTEST.GE.10)  WRITE(6,*) ' Returned from REF_CNFVEC'
*. And then do the Sigma from M space to S space
        IF(IDOHC.EQ.1) THEN
          ICSPC_CN_SAVE = ICSPC_CN
          ISSPC_CN_SAVE = ISSPC_CN
          ICSPC_CN = IMSPC_CN
          ISSPC_CN = ISSPC_CN
          LUC = 0
          LUHC = 0
          CALL SIGMA_CONF(WORK(KLC_CSF),HC,LUC,LUHC)
          IF(NTEST.GE.1000) WRITE(6,*) ' Home from SIGMA_CONF'
*. And restore
          ICSPC_CN = ICSPC_CN_SAVE 
          ISSPC_CN = ISSPC_CN_SAVE
        END IF !DOHC = 1
      END IF! switch between different algorithms 
*
C?    WRITE(6,*) ' TEST, NTEST = ', NTEST
      IF(NTEST.GE.100) THEN
        NCSF_SSPC_CNF = NCSF_PER_SYM_GN(ISSM,ISSPC_CN)
*
        WRITE(6,*) ' Final vectors from SIGMA_NORTCI '
        WRITE(6,*) ' ================================'
        IF(IDOSC.EQ.1) THEN
          WRITE(6,*) ' Metric times C vector: '
          CALL WRTMAT(SC,1,NCSF_SSPC_CNF,1,NCSF_SSPC_CNF)
        END IF
*
        IF(IDOHC.EQ.1) THEN
          WRITE(6,*) ' Hamiltonian times C vector: '
          CALL WRTMAT(HC,1,NCSF_SSPC_CNF,1,NCSF_SSPC_CNF)
       END IF
      END IF
*. And clean up 
      NOCSF = NOCSF_SAVE 
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'NOCISI')
      WRITE(6,*) ' Returning from SIGMA_NORTCI'
COLD  STOP ' Enforced stop at end of SIGMA_NORTCI' 
      RETURN
      END 
      SUBROUTINE SCA_VEC_TO_BLKV_DISC(VEC,ISCA,NELMNT,LUOUT,NBLOCK,
     &          LBLOCK,VECBLK,IREW)
* A vector is given in compact as elements and scatter array.
* Write this vector to disc, FILE LUOUT, in blocked form as defined by LBLOCK.
* Vecblk shoul be able to hold largest block
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      DIMENSION VEC(NELMNT), VECBLK(*)
      INTEGER LBLOCK(NBLOCK), ISCA(NELMNT)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from SCA_VEC_TO_BLKV_DISC '
        WRITE(6,*) ' LUOUT, IREW = ', LUOUT, IREW
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ISCA: '
        CALL IWRTMA(ISCA,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' LBLOCK: '
        CALL IWRTMA(LBLOCK,1,NBLOCK,1,NBLOCK)
        WRITE(6,*) ' Input vector: '
        CALL WRTMAT(VEC,1,NELMNT,1,NELMNT)
      END IF
*
      IF(IREW.EQ.1) THEN
       CALL REWINO(LUOUT)
      END IF
*
      IB_BL = 1
      DO IBLOCK = 1, NBLOCK
        LEN_BL = LBLOCK(IBLOCK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IBLOCK, LEN_BLK = ', IBLOCK, LEN_BL
        END IF
*
        ZERO = 0.0D0
        CALL SETVEC(VECBLK,ZERO,LEN_BL)
*
*. Find and copy elements in input vector that are in block IBLOCK
        DO IELMNT = 1, NELMNT
          JSCA     = ISCA(IELMNT)
          JSCA_ABS = IABS(JSCA)
          IF(NTEST.GE.10000) WRITE(6,*) ' IELMNT, JSCA, JSCA_ABS =', 
     &                 IELMNT, JSCA, JSCA_ABS
          IF(IB_BL.LE.JSCA_ABS.AND.JSCA_ABS.LE.IB_BL+LEN_BL-1) THEN
* Element is in block
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' Element in block, IELMNT, JSCA', IELMNT,JSCA
              WRITE(6,*) ' Output address = ', JSCA_ABS-IB_BL + 1
            END IF
            IF(JSCA.GT.0) THEN
              VECBLK(JSCA_ABS-IB_BL + 1) = VEC(IELMNT)
            ELSE
              VECBLK(JSCA_ABS-IB_BL + 1) =-VEC(IELMNT)
            END IF
          END IF
        END DO
*. Write block  to disc
        CALL ITODS(LEN_BL,1,-1,LUOUT)
        CALL TODSCP(VECBLK,LEN_BL,-1,LUOUT)
        IB_BL = IB_BL + LEN_BL
      END DO! End of loop over blocks
*. Write end of file
      CALL ITODS(-1,1,-1,LUOUT)
*
      RETURN
      END
      SUBROUTINE GAT_VEC_FROM_BLKV_DISC(VEC,ISCA,NELMNT,LUIN,NBLOCK,
     &          LBLOCK,VECBLK,IREW)
* A vector is given in disc, file LUIN,  with block-structure defined by LBLOCK
* Obtain elements given by scatter vector ISCA
* Vecblk shoul be able to hold largest block
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      DIMENSION VEC(NELMNT), VECBLK(*)
      INTEGER LBLOCK(NBLOCK), ISCA(NELMNT)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from GAT_VEC_TO_BLKV_DISC '
        WRITE(6,*) ' LUIN, IREW = ', LUIN, IREW
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ISCA: '
        CALL IWRTMA(ISCA,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' LBLOCK: '
        CALL IWRTMA(LBLOCK,1,NBLOCK,1,NBLOCK)
        WRITE(6,*) ' Initial vector on disc '
        CALL WRTVCD(VECBLK,LUIN,1,-1)
      END IF
*
      IF(IREW.EQ.1) THEN
       CALL REWINO(LUIN)
      END IF
*
      IB_BL = 1
      DO IBLOCK = 1, NBLOCK
        LEN_BL = LBLOCK(IBLOCK)
*. Obtain block 
        CALL IFRMDS(LBL,1,-1,LUIN)
        IF(LBL.NE.LEN_BL) THEN
          WRITE(6,*) 
     &    ' Difference between expected and actual block sizes',
     &   LEN_BL, LBL
         STOP
     &    ' Difference between expected and actual block sizes'
        END IF
        NO_ZEROING = 0
        CALL FRMDSC2(VECBLK,LBL,-1,LUIN,IMZERO,IAMPACK,
     &         NO_ZEROING)
*. Find and copy elements from input vector that are in block IBLOCK
        DO IELMNT = 1, NELMNT
          JSCA     = ISCA(IELMNT)
          JSCA_ABS = IABS(JSCA)
          IF(IB_BL.LE.JSCA_ABS.AND.JSCA_ABS.LE.IB_BL+LEN_BL-1) THEN
* Element is in block
            IF(JSCA.GT.0) THEN
              VEC(IELMNT) = VECBLK(JSCA_ABS-IB_BL + 1) 
            ELSE
              VEC(IELMNT) =-VECBLK(JSCA_ABS-IB_BL + 1) 
            END IF
          END IF
        END DO
        IB_BL = IB_BL + LEN_BL
      END DO! End of loop over blocks
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Vector gathered from DISC '
       CALL WRTMAT(VEC,1,NELMNT,1,NELMN)
      END IF
*
      RETURN
      END
      SUBROUTINE LBLOCK_FOR_CIXP(LBLOCK,NBLOCK,ICISPC,ISYM)
*
* Obtain number of blocks and lengths of blocks for CI expansion 
* in space ICISPC and symmetry ISYM
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'csm.inc'
*
      NTEST = 000
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LBLCIX')
*
* Number of occupation classes
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      NAEL = NELFTP(IATP)
      NBEL = NELFTP(IBTP)
      NEL = NAEL + NBEL
*
      IWAY = 1
      CALL OCCLS(1,NOCCLS,IOCCLS,NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),0,0,NOBPT)
*. and the occupation classes
      CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL')
      CALL MEMMAN(KLBASSPC,NOCCLS,'ADDL  ',1,'BASSPC')
      IWAY = 2
      CALL OCCLS(2,NOCCLS,WORK(KLOCCLS),NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),1,WORK(KLBASSPC),NOBPT)
*. Allocate space for largest encountered number of TTSS blocks
      NTTS = MXNTTS
C     WRITE(6,*) ' GASCI : NTTS = ', NTTS
*.
      CALL MEMMAN(KLCLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLCLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLCI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLCIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
      CALL MEMMAN(KLC2B  ,  NTTS,'ADDL  ',1,'C2BT  ')
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Matrix giving allowed combination of alpha- and beta-strings
      CALL IAIBCM(ICISPC,WORK(KLCIOIO))
*. option KSVST not active so
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,WORK(KLCBLTP),WORK(KSVST))
*. Blocks of  CI vector, using a single batch for complete  expansion
      ICOMP = 1
      ISIMSYM = 1
      CALL PART_CIV2(IDC,WORK(KLCBLTP),WORK(KNSTSO(IATP)),
     &              WORK(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLCIOIO),
     &              ISMOST(1,ISYM),
     &              NBATCH,WORK(KLCLBT),WORK(KLCLEBT),
     &              WORK(KLCI1BT),WORK(KLCIBT),ICOMP,ISIMSYM)
*. Number of BLOCKS
        NBLOCK = IFRMR(WORK(KLCI1BT),1,NBATCH)
     &         + IFRMR(WORK(KLCLBT),1,NBATCH) - 1
        IF(NTEST.GE.1000) WRITE(6,*) ' Number of blocks ', NBLOCK
*. And the lengths of the various blocks
*
      CALL EXTRROW(WORK(KLCIBT),8,8,NBLOCK,LBLOCK)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info in CI space ', ICISPC, ' with sym ', ISYM
        WRITE(6,*) ' =============================================='
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks: ', NBLOCK
        WRITE(6,*) ' Length of each block: '
        CALL IWRTMA(LBLOCK,1,NBLOCK,1,NBLOCK)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'LBLCIX')
*
      RETURN
      END
      SUBROUTINE MINMAX_EXCIT(IOCC_MIN_IN,IOCC_MAX_IN,NEXCIT,
     &                        IOCC_MIN_OUT,IOCC_MAX_OUT,NORB)
*
* A CI space is defined by accumulated occations IOCC_MIN_IN, IOCC_MAX_IN.
* Apply NEXCIT excitations to this space to obtain IOCC_MIN_OUT,IOCC_MAX_OUT
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC_MIN_IN(NORB),IOCC_MAX_IN(NORB)
*. Output
      INTEGER IOCC_MIN_OUT(NORB),IOCC_MAX_OUT(NORB)
*
      NTEST = 00
*
      NELEC = IOCC_MIN_IN(NORB)
      DO IORB = 1, NORB
       IOCC_MIN_OUT(IORB) = 
     & MAX(0,IOCC_MIN_IN(IORB)-NEXCIT,NELEC-2*(NORB-IORB))
       IOCC_MAX_OUT(IORB) = MIN(2*IORB,NELEC,IOCC_MAX_IN(IORB)+NEXCIT)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from MINMAX_EXCIT '
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' allowed excitation level = ', NEXCIT
        WRITE(6,*) ' Input occupation constraints '
        CALL WRT_MINMAX_OCC(IOCC_MIN_IN,IOCC_MAX_IN,NORB)
        WRITE(6,*) ' Output occupation constraints '
        CALL WRT_MINMAX_OCC(IOCC_MIN_OUT,IOCC_MAX_OUT,NORB)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_MINMAX_OCC(IOCC_MIN,IOCC_MAX,NORB)
*
* Write min and max accumulated occupation arrays
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC_MIN(NORB),IOCC_MAX(NORB)
*
      WRITE(6,*) ' Min and Max accumulated occupations: '
      WRITE(6,*)
      WRITE(6,*) ' Orbital Min. occ Max. occ '
      WRITE(6,*) ' =========================='
      DO IORB = 1, NORB
        WRITE(6,'(3X,I4,2(4X,I4))')
     &  IORB, IOCC_MIN(IORB), IOCC_MAX(IORB)
      END DO
*
      RETURN
      END
      SUBROUTINE PRECOND_NORTCI
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*
      WRITE(6,*) ' Dummy PRECOND_NORTCI entered'
      STOP       ' Dummy PRECOND_NORTCI entered'
*
      
      END
      SUBROUTINE GET_CBIO(C,CBIOMO,CBIOAO)
* A MO-AO transformation matrix, C,  to a (non-orthogonal) basis is given.
* Obtain the corresponding bioorthogonal transformation matrix CBIOMO
* (Bio = > MO's in C) and CBIOAO (Bio => AO'S)
*
* Jeppe Olsen, July 2011 for the nonorthogonal CI work
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
*. Input
      DIMENSION C(*)
*. Output
      DIMENSION CBIOMO(*), CBIOAO(*)
*
      NTEST = 00
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_CB')
*
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN_M = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLSAOE,LEN_M,'ADDL  ',2,'SAO_E ')
      CALL MEMMAN(KLMSCR,LEN_M,'ADDL  ',2,'MSCR  ')
      CALL MEMMAN(KLSCR,2*LEN_M,'ADDL  ',2,'SCR   ')
*. Expand SAO
      CALL TRIPAK_AO_MAT(WORK(KLSAOE),WORK(KSAO),2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Expanded SAO '
        CALL APRBLM2(WORK(KLSAOE),NTOOBS,NTOOBS,NSMOB,0) 
      END IF
*. Obtain Metric in MO-basis, SMO,  in CBIOAO
C  TRAN_SYM_BLOC_MAT4(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      CALL TRAN_SYM_BLOC_MAT4(WORK(KLSAOE),C,C,NSMOB,NTOOBS,NTOOBS, 
     &CBIOAO,WORK(KLSCR),0)
C CBIOMO = SMO ** -1
      IPROBLEM = 0
C     INV_BLKMT(A,AINV,SCR,NBLK,LBLK,IPROBLEM)
      CALL INV_BLKMT(CBIOAO,CBIOMO,WORK(KLSCR),NSMOB,NTOOBS,
     &               IPROBLEM)
      IF(IPROBLEM.NE.0) THEN
        WRITE(6,*) ' Problem inverting matrix C(T) S(AO) '
        STOP       ' Problem inverting matrix C(T) S(AO) '
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' CBIOMO = SMO ** -1 '
        CALL APRBLM2(CBIOMO,NTOOBS,NTOOBS,NSMOB,0) 
      END IF
* CBIOAO = C * CBIOMO
      CALL MULT_BLOC_MAT(CBIOAO,C,CBIOMO,
     &     NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
*
* Check: Calculate C(T) S CBIO 
C  TRAN_SYM_BLOC_MAT4(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      I_DO_CHECK = 0
      IF(I_DO_CHECK.EQ.1) THEN
       CALL TRAN_SYM_BLOC_MAT4(WORK(KLSAOE),C,CBIOAO,NSMOB,NTOOBS,
     & NTOOBS, 
     & WORK(KLMSCR),WORK(KLSCR),0)
       WRITE(6,*) ' C(T) S CBIO '
       CALL APRBLM2(WORK(KLMSCR),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*

      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' Bioorthogonal MOAO expansion matrix '
       WRITE(6,*) ' =================================== '
       WRITE(6,*)
       CALL APRBLM2(CBIOAO,NTOOBS,NTOOBS,NSMOB,0) 
       WRITE(6,*)
       WRITE(6,*) ' Bioorthogonal MOMO expansion matrix '
       WRITE(6,*) ' =================================== '
       WRITE(6,*)
       CALL APRBLM2(CBIOMO,NTOOBS,NTOOBS,NSMOB,0) 
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_CB')
      RETURN
      END
      SUBROUTINE INV_BLKMT(A,AINV,SCR,NBLK,LBLK,IPROBLEM)
*
* Invert blocked matrix  A to give AINV
* Problems with inversion is flagged by IPROBLEM.NE. 0
* IPROBLEM gives last block with problems
*
* SCR should at least be twice the size of the largest block
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(*)
      INTEGER LBLK(NBLK)
*. Output
      DIMENSION AINV(*)
*Scratch
      DIMENSION  SCR(*)
*
      NTEST = 000
*
      IPROBLEM = 0
      DO IBLK = 1, NBLK
        IF(IBLK.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + LBLK(IBLK-1)**2
        END IF
        LEN_BLK = LBLK(IBLK)
        CALL COPVEC(A(IOFF),SCR,LEN_BLK**2)
        IF(NTEST.GE.1000) THEN 
          WRITE(6,*) ' Matrix to be inverted '
          CALL WRTMAT(SCR,LEN_BLK,LEN_BLK,LEN_BLK,LEN_BLK)
        END IF
C       INVMAT(A,B,MATDIM,NDIM,ISING)
        CALL INVMAT(SCR,SCR(1+LEN_BLK**2),LEN_BLK,LEN_BLK,ISING)
        IF(ISING.NE.0) IPROBLEM = IBLK
        CALL COPVEC(SCR,AINV(IOFF),LEN_BLK**2)
      END DO
*
      IF(IPROBLEM.NE.0) THEN
        WRITE(6,*) 
     &  ' Problem in INV_BLKMAT, number of last singular block =', 
     &   IPROBLEM
        WRITE(6,*) ' Complete input block matrix '
        CALL APRBLM2(A,LBLK,LBLK,NBLK,0)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverted block matrix:'
        CALL APRBLM2(AINV,LBLK,LBLK,NBLK,0)
C            APRBLM2(A,LROW,LCOL,NBLK,ISYM)
      END IF
*
      RETURN
      END
      SUBROUTINE COMHAM_HS_GEN(MSTV,NDIM)
*
* Set up Complete Hamiltonian matrices using external
* routine MSTV
*
* Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'clunit.inc'
*
      PARAMETER(MXLDIM = 200)
      DIMENSION H(MXLDIM*MXLDIM), S(MXLDIM*MXLDIM)
      DIMENSION VEC1(MXLDIM),VEC2(MXLDIM),VEC3(MXLDIM)
      DIMENSION SCR(5*MXLDIM**2+2*MXLDIM), EIGVEC(MXLDIM**2)
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3, 
     &       LUCBIO_SAVE, LUHCBIO_SAVE
*
      EXTERNAL MSTV
*
      LUSCR1 = LUSC34
      LUSCR2 = LUSC35
      LUSCR3 = LUSC36
      LUCBIOSAVE = 0
      LUHCBIOSAVE = 0
*
      NTEST = 1000
      IF(NDIM.GT.MXLDIM) THEN
         WRITE(6,*) 
     &  ' COMHAM_HS_GEN called with larger dimension than allowed '
         WRITE(6,*) ' ALlowed (MXLDIM) and actual (NDIM) dimensions ',
     &   MXLDIM, NDIM
        WRITE(6,*) 'LUCIA suggests that you increase MXLDIM '
        STOP
     &  ' COMHAM_HS_GEN called with larger dimension than allowed '
      END IF
*. Restrict
      NDIML = NDIM
*
      ZERO = 0.0D0
      ONE = 1.0D0
      DO IVEC = 1, NDIML
       CALL SETVEC(VEC1,ZERO,NDIM)
       VEC1(IVEC) = ONE
       CALL MSTV(VEC1,VEC2,VEC3,1,1)
*
       IOFF = (IVEC-1)*NDIML+1
       CALL COPVEC(VEC2,H(IOFF),NDIML)
       CALL COPVEC(VEC3,S(IOFF),NDIML)
      END DO
*
      IF(NTEST.GE.1000) THEN
      WRITE(6,*) ' matrices from COMHAM_HS_GEN'
      CALL WRTMAT(H,NDIML,NDIML,NDIML,NDIML)
      CALL WRTMAT(S,NDIML,NDIML,NDIML,NDIML)
      END IF
*
      I_DO_DIAG = 1
      IF(I_DO_DIAG.EQ.1) THEN
C     GENEIG_WITH_SING_CHECK(A,S,EIGVEC,EIGVAL,NVAR,NSING,
C    &                                  WORK,IASPACK)
        CALL GENEIG_WITH_SING_CHECK(H,S,EIGVEC,VEC1,NDIML,
     &       NSING,SCR,0)
      END IF
*
      RETURN
      END  
      SUBROUTINE EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT
     &           (A,AGAS,IGAS,JGAS,I_EX_OR_CP)
*
* A symmetryblocked (not lower half packed) matrix A over orbitals is given
* Extract blocks referring to GASpaCE IGAS, JGAS
*
* I_EX_OR_CP = 1 => Extract from A to IGAS
* I_EX_OR_CP = 1 => Copy from IGAS to A
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Specific input or output
      DIMENSION A(*), AGAS(*)
*. Scratch- for output
      DIMENSION IDIM(MXPNGAS), JDIM(MXPNGAS)
*
      DO ISYM = 1, NSMOB
       IF(ISYM.EQ.1) THEN
        IOFF_IN = 1
        IOFF_OUT = 1
       ELSE
        IOFF_IN = IOFF_IN + NTOOBS(ISYM-1)**2
        IOFF_OUT = 
     &  IOFF_OUT + NOBPTS_GN(IGAS,ISYM-1)*NOBPTS_GN(JGAS,ISYM-1)
       END IF
*
       IIOFF = 1
       DO IIGAS = 0, IGAS -1
         IIOFF = IIOFF + NOBPTS_GN(IIGAS,ISYM)
       END DO
*
       IJOFF = 1
       DO IIGAS = 0, JGAS -1
         IJOFF = IJOFF + NOBPTS_GN(IIGAS,ISYM)
       END DO
*
       NI = NOBPTS_GN(IGAS,ISYM)
       NJ = NOBPTS_GN(JGAS,ISYM)
       NIS = NTOOBS(ISYM)
       NJS = NTOOBS(ISYM)
       DO J = 1, NJ
         DO I = 1, NI
           IJ_OUT = IOFF_OUT -1 + (J-1)*NI + I 
           IJ_IN  = IOFF_IN -1 
     &            + (IJOFF+J-1-1)*NIS + IIOFF+I-1
           IF(I_EX_OR_CP.EQ.1) THEN
             AGAS(IJ_OUT) = A(IJ_IN)
           ELSE
             A(IJ_IN) = AGAS(IJ_OUT)
           END IF
         END DO
       END DO
      END DO ! End of loop over symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Submatrix with IGAS, JGAS = ',
     &   IGAS, JGAS
         CALL EXTRROW(NOBPTS_GN,IGAS+1,7+MXPR4T,NSMOB,IDIM)
         CALL EXTRROW(NOBPTS_GN,JGAS+1,7+MXPR4T,NSMOB,JDIM)
C             EXTRROW(INMAT,IROW,NROW,NCOL,IOUTVEC)
C APRBLM2(A,LROW,LCOL,NBLK,ISYM)
         CALL APRBLM2(AGAS,IDIM,JDIM,NSMOB,0)
         WRITE(6,*) ' Full matrix '
         CALL APRBLM2(A,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE PREPARE_CMOAO_INI
     &(INI_MO_TP, CMOAO_OUT,CMOAO_IN,IVBGAS)
*
* Obtain initial orbitals for Lucia calculation
*
* INI_MO_TP = 1 => CMOAO_OUT = 1
*       = 2 => Transform MO's so diagonal block IVBGAS is a unit matrix
*       = 3 => CMOAO_OUT = CMOAO_IN
*       = 4 => from fragment MO's
* Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
*. Input
      DIMENSION CMOAO_IN(*)
*. Output
      DIMENSION CMOAO_OUT(*)
*. Local scratch
      DIMENSION IDIMV(MXPOBS), IDIMI(MXPOBS)
*
      IDUM = 0
      NTEST = 10
*
      IF(NTEST.GE.10) 
     &WRITE(6,*) ' PREPARE..., INI_MO_TP = ', INI_MO_TP
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PREPMO')
*
      IF(INI_MO_TP.EQ.1) THEN
*
*. CMOAO_OUT = 1
*
        ONE = 1.0D0
        CALL SETDIA_BLM(CMOAO_OUT,ONE,NSMOB,NTOOBS,0)      
      ELSE IF ( INI_MO_TP.EQ.2) THEN
*
* Rotate orbitals in GASpace IVBGAS, so the diagonal IVBGAS block 
* become diagonal- could require pivoting
*
        LEN1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
C                NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
        CALL MEMMAN(KLMO1,LEN1_F,'ADDL  ',2,'MO1   ')
        CALL MEMMAN(KLMO2,LEN1_F,'ADDL  ',2,'MO2   ')
        CALL MEMMAN(KLSCR,2*LEN1_F,'ADDL  ',2,'SCR   ')
*
        CALL COPVEC(CMOAO_IN,CMOAO_OUT,LEN1_F)
*. Extract block (IVBGAS,IVBGAS) of CMO
C     EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT
C    &           (A,AGAS,IGAS,JGAS,I_EX_OR_CP)
        CALL EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT
     &           (CMOAO_IN,WORK(KLMO1),IVBGAS,IVBGAS,1)
*. Number of orbitals per sym in this space
         CALL EXTRROW(NOBPTS_GN,IVBGAS+1,7+MXPR4T,NSMOB,IDIMV)
*. Invert block and save in KLMO2
C             INV_BLKMT(A,AINV,SCR,NBLK,LBLK,IPROBLEM)
         CALL INV_BLKMT(WORK(KLMO1),WORK(KLMO2),WORK(KLSCR),NSMOB,
     &                   IDIMV,IPROBLEM)
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Inverted diagonal GAS block'
           CALL APRBLM2(WORK(KLMO2),IDIMV,IDIMV,NSMOB,0)
         END IF
*. Multiply inverted block on ini MO's in space IVBGAS
         DO IGAS = 0, NGAS +1
*. Extract block (IGAS,IVBGAS) in KLMO1
          CALL EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT
     &         (CMOAO_IN,WORK(KLMO1),IGAS,IVBGAS,1)
*. Dimensions of block IGAS
          CALL EXTRROW(NOBPTS_GN,IGAS+1,7+MXPR4T,NSMOB,IDIMI)
*         CMOAO_IN(IGAS,IVBGAS)*CMOAO_IN(IGAS,IGAS)**(-1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' C(IGAS,IVGAS) block '
          CALL APRBLM2(WORK(KLMO1),IDIMI,IDIMV,NSMOB,0)
          WRITE(6,*) ' C(IVGAS,IVGAS)** (-1) block'
          CALL APRBLM2(WORK(KLMO2),IDIMV,IDIMV,NSMOB,0)
        END IF
C       MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
          CALL MULT_BLOC_MAT(WORK(KLSCR),WORK(KLMO1),WORK(KLMO2),
     &    NSMOB,IDIMI,IDIMV, IDIMI,IDIMV,IDIMV,IDIMV,0)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) 
     &   ' C(IGAS,IVGAS)*C**(-1)(IVGAS,IVGAS) for IGAS = ', IGAS
           CALL APRBLM2(WORK(KLSCR),IDIMI,IDIMV,NSMOB,0)
          END IF

*. And copy to CMOAO_OUT
          CALL EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT
     &         (CMOAO_OUT,WORK(KLSCR),IGAS,IVBGAS,2)
         END DO
*
      ELSE IF(INI_MO_TP.EQ.3.OR.INI_MO_TP.EQ.5) THEN
*
* CMOAO_OUT = CMOAO_IN
*
        LEN1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
        CALL COPVEC(CMOAO_IN,CMOAO_OUT,LEN1_F)
      ELSE IF(INI_MO_TP.EQ.4) THEN
* obtain MO's from Fragment AO's
        CALL GET_CMO_FROM_FRAGMENTS(CMOAO_OUT)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from  PREPARE CMOAO_INI_NORTCI'
        WRITE(6,*) ' ====================================='
        WRITE(6,*)
        WRITE(6,*) ' INI_MO_TP = ', INI_MO_TP
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output set of MOs '
        CALL APRBLM2(CMOAO_OUT,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PREPMO')
*
      RETURN
      END
      SUBROUTINE GET_CMO_FROM_FRAGMENTS(CMO)
* Obtain MOAO coefficients CMO from fragments as specified by
* N_GS_SM_BAS_FRAG
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'fragmol.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
*. Molecule to fragment symmetry
      INTEGER LSYMEXP(8) 
*. Output
      DIMENSION CMO(*)
*
      NTEST = 10
*
*
* 1: Check information in fragments with total number of orbitals
*    and basis functions 
*
*
*. Total number of orbitals per symmetry
      NERROR = 0
      DO ISYM = 1, NSMOB
       NNTOOBS = 0
*. Loop over equivalent groups of atoms
       DO IEQV = 1, NEQVGRP_FRAG
         IFRAG = IEQVGRP_FRAG(1,IEQV)
         IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IEQV, IFRAG = ', IEQV, IFRAG
            WRITE(6,*) ' LEQVGRP_FRAG(IEQV) = ', LEQVGRP_FRAG(IEQV)
         END IF
         IF(LEQVGRP_FRAG(IEQV).EQ.1) THEN
*. No expansion of symmetries
           DO JSYM = 1, NSMOB
             LSYMEXP(JSYM) = JSYM
           END DO
         ELSE IF(LEQVGRP_FRAG(IEQV).EQ.2) THEN
           IF(NSMOB.EQ.4) THEN
*. Assumed expansion from Cs to C2V
             LSYMEXP(1) = 1
             LSYMEXP(2) = 2
             LSYMEXP(3) = 1
             LSYMEXP(4) = 2
           ELSE IF(NSMOB.EQ.8) THEN
*. Assumed expansion from C2V to D2H
             LSYMEXP(1) = 1
             LSYMEXP(2) = 2
             LSYMEXP(3) = 3
             LSYMEXP(4) = 4
             LSYMEXP(5) = 1
             LSYMEXP(6) = 2
             LSYMEXP(7) = 3
             LSYMEXP(8) = 4
           ELSE
              WRITE(6,*) ' Combination not programmed(2) '
              WRITE(6,*) ' IEQV, LEQVGRP_FRAG, NSMOB = ', 
     &                     IEQV, LEQVGRP_FRAG(IEQV), NSMOB
              STOP       ' Combination not programmed '
           END IF
         ELSE IF(LEQVGRP_FRAG(IEQV).EQ.4) THEN
           IF(NSMOB.EQ.8) THEN
*. Assumed expansion from CS to D2H
             LSYMEXP(1) = 1
             LSYMEXP(2) = 2
             LSYMEXP(3) = 3
             LSYMEXP(4) = 4
             LSYMEXP(5) = 1
             LSYMEXP(6) = 2
             LSYMEXP(7) = 3
             LSYMEXP(8) = 4
           ELSE
              WRITE(6,*) ' Combination not programmed(3) '
              WRITE(6,*) ' LEQVGRP_FRAG, NSMOB = ', LEQVGRP_FRAG, NSMOB
              STOP       ' Combination not programmed '
           END IF
         END IF ! Switch between dimension of equivalence class
         NNTOOBS = NNTOOBS + NBAS_FRAG(LSYMEXP(ISYM),IFRAG)
         IF(NTEST.GE.1000) WRITE(6,*) ' ISYM, LSYM, IFRAG, NBAS = ',
     &     ISYM,LSYMEXP(ISYM),IFRAG,NBAS_FRAG(LSYMEXP(ISYM),IFRAG)
       END DO ! Loop over equivalent groups of atoms
*
       IF(NNTOOBS.NE.NTOOBS(ISYM)) THEN
        WRITE(6,*) 
     &  ' Number of basisfunctions from fragments is not correct '
        WRITE(6,'(A,3I3)') ' ISYM, NTOOBS, Sum of fragments: ',
     &  ISYM, NTOOBS(ISYM),NNTOOBS
        NERROR = NERROR + 1
       END IF
      END DO
*. Check internal consistency for each fragment
      DO IIFRAG = 1, NFRAG_MOL
       IFRAG = IFRAG_MOL(IIFRAG)
       NSMOB_L = NSMOB_FRAG(IFRAG)
       DO ISYM = 1, NSMOB_L
        NNTOOBS_FRAG = 0
        DO IGAS = 0, NGAS + 1
         NNTOOBS_FRAG  = 
     &   NNTOOBS_FRAG + N_GS_SM_BAS_FRAG(IGAS,ISYM,IFRAG)
        END DO
        IF(NNTOOBS_FRAG.NE.NBAS_FRAG(ISYM,IFRAG)) THEN
          WRITE(6,*) 
     &    ' Inconsistency between N_GS_SM_BAS_FRAG and NBAS_FRAG'
          WRITE(6,'(A,4I3)') 
     &    ' IFRAG, ISYM, Sum over gaspaces and Required ',
     &    IFRAG, ISYM, NNTOOBS_FRAG, NBAS_FRAG(ISYM,IFRAG)
          NERRROR = NERROR + 1
        END IF
       END DO
      END DO
*. Check consistency for each GASpace and symmetry
      WRITE(6,*) ' Warning: some consistency checks skipped '
      WRITE(6,*) ' Warning: some consistency checks skipped '
      WRITE(6,*) ' Warning: some consistency checks skipped '
      WRITE(6,*) ' Warning: some consistency checks skipped '
      WRITE(6,*) ' Warning: some consistency checks skipped '
      WRITE(6,*) ' Warning: some consistency checks skipped '
      WRITE(6,*) ' Warning: some consistency checks skipped '
CTEMP DO IGAS = 0, NGAS + 1
CTEMP  DO ISYM = 1, NSMOB
CTEMP   NNTOOBS_GS_SM = 0
CTEMP   DO IIFRAG =1, NFRAG_MOL
CTEMP     IFRAG = IFRAG_MOL(IIFRAG)
CTEMP     NNTOOBS_GS_SM = 
CTEMP&    NNTOOBS_GS_SM + N_GS_SM_BAS_FRAG(IGAS,ISYM,IFRAG)
CTEMP   END DO
*
CTEMP   IF(NNTOOBS_GS_SM.NE.NOBPTS_GN(IGAS,ISYM)) THEN
CTEMP    WRITE(6,*) 
CTEMP&   ' Inconsistency in number of orbitals of given SYM and GAS'
CTEMP    WRITE(6,'(A,4I4)') ' ISYM, IGAS, Sum over fragments, Total ',
CTEMP&   ISYM, IGAS, NNTOOBS_GS_SM, NOBPTS_GN(IGAS,ISYM) 
CTEMP    NERROR = NERROR + 1
CTEMP   END IF
CTEMP  END DO
CTEMP END DO
*
      IF(NERROR.NE.0) THEN
       WRITE(6,*) 
     & ' Inconsistency between info on fragments and molecule '
C!     STOP
C!   & ' Inconsistency between info on fragments and molecule '
      END IF
*
* 2: And then set up the CMO matrix from fragment info
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ================================== '
        WRITE(6,*) ' CMO(FRAGMENTS) => CMO(MOLECULE)(1) '
        WRITE(6,*) ' ================================== '
      END IF
*
      DO ISYM = 1, NSMOB
       IF(NTEST.GE.1000) WRITE(6,*) ' ISYM = ', ISYM
       IF(ISYM.EQ.1) THEN
        IB_CMOL = 1
       ELSE
        IB_CMOL = IB_CMOL + NTOOBS(ISYM-1)**2
       END IF
       NOB_SM = NTOOBS(ISYM)
       ZERO = 0.0D0
       CALL SETVEC(CMO(IB_CMOL),ZERO,NOB_SM**2)
       IOFF_ORB = 1
       IOFF_BAS = 1
       JMO = 0
       DO IGAS = 0, NGAS + 1
         IB_BAS_MOL = 1
         IF(NTEST.GE.1000) WRITE(6,*) ' IGAS = ', IGAS
*. Loop over equivalent set of fragment orbitals
         DO IEQV = 1, NEQVGRP_FRAG
* First fragment of class 
           IFRAG = IEQVGRP_FRAG(1,IEQV)
           IF(NTEST.GE.1000) 
     &     WRITE(6,*) ' IEQV, IFRAG = ', IEQV, IFRAG
*
           XL = DBLE(LEQVGRP_FRAG(IEQV))
           SCALE = 1.0D0/SQRT(XL)
*. Symmetry in fragment
     
           IF(LEQVGRP_FRAG(IEQV).EQ.1) THEN
            ISYML = ISYM
           ELSE IF (LEQVGRP_FRAG(IEQV).EQ.2) THEN
            IF(NSMOB.EQ.8) THEN
              ISYML = ISYM
              IF(ISYM.GT.4) ISYML = ISYM-4
            ELSE
              WRITE(6,*) ' Symmetry reduction not programmed(1) '
              WRITE(6,*) ' ISYM, NSMOB, LEQVGRP_FRAG = ',
     &                     ISYM, NSMOB, LEQVGRP_FRAG(IEQV)
              STOP ' Symmetry reduction not programmed '
            END IF
           END IF
           IF(NTEST.GE.1000) 
     &     WRITE(6,*) ' ISYM,ISYML = ', ISYM,ISYML
*. Address of symmetryblock in C for fragment
           IB_C_FRAG = 1
           DO JSYM = 1, ISYML-1
             IB_C_FRAG = IB_C_FRAG + NBAS_FRAG(JSYM,IFRAG)**2
           END DO
*. Start and number of orbitals in input fragment
           IB_OB_FRAG = 1
           DO JGAS = 0, IGAS - 1
             IB_OB_FRAG = IB_OB_FRAG 
     &     + N_GS_SM_BAS_FRAG(JGAS,ISYML,IFRAG)
           END DO
           IF(NTEST.GE.1000) WRITE(6,*) ' IB_OB_FRAG = ',
     &     IB_OB_FRAG
           N_OB_GS_SM_FRAG = N_GS_SM_BAS_FRAG(IGAS,ISYML,IFRAG)
           N_OB_SM_FRAG = NBAS_FRAG(ISYML,IFRAG)
           DO JJMO = 1, N_OB_GS_SM_FRAG
           JMO = JMO + 1
           IF(NTEST.GE.1000) WRITE(6,*) ' Info for Orbital ', JMO
           DO IIMO = 1, N_OB_SM_FRAG
            IF(NTEST.GE.1000) WRITE(6,*) ' JJMO, IIMO = ', 
     &      JJMO, IIMO
            IF(NTEST.GE.1000) WRITE(6,*) ' IB_BAS_MOL = ',
     &      IB_BAS_MOL
            IADR_OUT = IB_CMOL-1+(JMO-1)*NOB_SM +IB_BAS_MOL-1 + IIMO 
            IADR_IN = IB_C_FRAG-1
     &            + (JJMO+IB_OB_FRAG-1-1)*N_OB_SM_FRAG
     &            + IIMO 
            IF(NTEST.GE.1000) WRITE(6,*) ' IADR_IN, IADR_OUT ',
     &      IADR_IN, IADR_OUT
            CMO(IADR_OUT) =  WORK(KCMOAO_FRAG(IFRAG)-1+IADR_IN)*SCALE
           END DO !loop over IIMO
           END DO !loop over JJMO
*. Start of basis functions for given sym and fragment in molecule
           IB_BAS_MOL = IB_BAS_MOL + N_OB_SM_FRAG
           IF(NTEST.GE.1000) 
     &     WRITE(6,*) ' IB_BAS_MOL, N_OB_SM_FRAG', 
     &                  IB_BAS_MOL, N_OB_SM_FRAG
         END DO ! End of loop over fragments
       END DO ! End of loop over GAspaces
      END DO ! End of loop over Symmetries
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' CMO matrix from fragments(not orthogonalized) '
       WRITE(6,*) ' =============================================='
       WRITE(6,*)
       CALL APRBLM_F7(CMO,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
COLD  STOP ' Jeppe enforced me to stop after CMO '
*
      RETURN
      END
C     ORT_ORB(WORK(KLCMOAO1),CMOAO_OUT,INTER_ORT, 
C    &     INTERGAS_ORT,INTRAGAS_OUT,IORT_VB)
      SUBROUTINE ORT_ORB(CMOAO_IN, CMOAO_OUT, 
     &           INTER_ORT,INTERGAS_ORT,
     &           INTRAGAS_ORT,IORT_VB)
*
* Two parts
* 1: Inter Gas orthogonaliztion
* 2: Intra Gas orthonormalization:
*
*. The inter gas orthogonalization: CMOAO_IN => CMOAO_OUT
* ==================================
* INTER_ORT = 1 => All GA Spaces are orthogonalized to inactive and 
*                    secondary space
*
* INTERGAS_ORT = 1 => Gaspaces are orthogonalized to each other
*
*. The Intra gas orthonormalization: CMOAO_OUT => CMOAO_OUT
* ====================================
* INTRAGAS_ORT = 0 => no Intra gas orthogonalization
*              = 1 => Intra gas orthogonalization using symmetric orthog
*              = 2 => Intra gas orthogonalization using orthog by diag
* IORT_VB   = 0 => No orthogonalization of space VB space
*           = 1 => orthog using method specified by INTRAGAS_ORT
*
* Note: If INTRAGAS_ORT = 1, then the VB orb space is left untouched, 
*       irrespectively of IORT_VB
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'crun.inc'
*. Input
      DIMENSION CMOAO_IN(*)
*. Output 
      DIMENSION CMOAO_OUT(*)
*. Local scratch
      INTEGER IDIM(MXPOBS)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ORTOBV')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from ORT_ORB '
        WRITE(6,*) ' ====================='
        WRITE(6,*)
        WRITE(6,'(A,2I4)') 
     &  ' INTER_ORT, INTERGAS_ORT = ',
     &    INTER_ORT, INTERGAS_ORT
        WRITE(6,'(A,2I4)')
     &  ' INTRAGAS_ORT, IORT_VB ',
     &    INTRAGAS_ORT, IORT_VB 
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input CMO coefficients '
        CALL APRBLM2(CMOAO_IN,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      IDUM = 0
*. Obtain metric over molecular orbitals
      LEN_1F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
C              NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      CALL MEMMAN(KLSMO,LEN_1F,'ADDL  ',2,'SMO   ')
      CALL MEMMAN(KLCMOAO2,LEN_1F,'ADDL  ',2,'MOAO2 ')
*. Obtain metric in MO basis in SMO
      IPACK_OUT = 0
      CALL GET_SMO(CMOAO_IN,WORK(KLSMO),IPACK_OUT)
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Overlap matrix for initial orbitals '
       CALL APRBLM2(WORK(KLSMO),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
* ==============================
* The intergas orthogonalization
* ==============================
*
*. Resulting MOAO transformation will be saved in CMOAO_OUT
      CALL COPVEC(CMOAO_IN,CMOAO_OUT,LEN_1F)
      CALL COPVEC(CMOAO_IN,WORK(KLCMOAO2),LEN_1F)
C?    WRITE(6,*) ' INTRAGAS_ORT after COPVEC(1)', INTRAGAS_ORT
*
      IF(INTER_ORT.EQ.1) THEN
*
*. Orthogonalize GAS spaces for inactive
*
        IF(NINOB.NE.0) THEN
          DO IGAS = 1, NGAS+1
*. Orthogonalize GAS IGAS to inactive
C                ORT_GAS_TO_GAS(IGAS,JGAS,SIN,CIN,COUT)
            CALL ORT_GAS_TO_GAS(0,IGAS,WORK(KLSMO),WORK(KLCMOAO2),
     &           CMOAO_OUT)
*. Test..
CT          CALL ORT_GAS_TO_GAS(IGAS,0,WORK(KLSMO),WORK(KLCMOAO2),
CT   &           CMOAO_OUT)
            CALL COPVEC(CMOAO_OUT,WORK(KLCMOAO2),LEN_1F)
            CALL COPVEC(CMOAO_OUT,WORK(KLCMOAO2),LEN_1F)
*. Update metric
            CALL GET_SMO(CMOAO_OUT,WORK(KLSMO),IPACK_OUT)
          END DO
        END IF
*
CM      IF(NSCOB.NE.0) THEN
*
*. Orthogonalize Secondary space to GASpaces
*
CM        DO IGAS = 1, NGAS
*. Orthogonalize Secondary to GAS IGAS 
CM          CALL ORT_GAS_TO_GAS(IGAS,NGAS+1,WORK(KLSMO),WORK(KLCMOAO2),
CM   &           CMOAO_OUT)
CM          CALL COPVEC(CMOAO_OUT,WORK(KLCMOAO2),LEN_1F)
*. Update metric
CM          CALL GET_SMO(CMOAO_OUT,WORK(KLSMO),IPACK_OUT)
CM        END DO
CM      END IF
*
C?    WRITE(6,*) ' INTRAGAS_ORT after INTER(1)', INTRAGAS_ORT
        IF(INTERGAS_ORT.EQ.1) then
* Orthogonalize JGAS to IGAS with JGAS > IGAS
          IF(NTEST.GE.10000) THEN
            WRITE(6,*) ' SMO before GAS GAS orthog'
            CALL APRBLM2(WORK(KLSMO),NTOOBS,NTOOBS,NSMOB,0)
            WRITE(6,*) ' MOAO2 before GAS GAS orthog '
            CALL APRBLM2(WORK(KLCMOAO2),NTOOBS,NTOOBS,NSMOB,0)
          END IF
          DO JGAS = 2, NGAS
            DO IGAS = 1, JGAS -1
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) 
     &          ' InterGAS orthogonalization for IGAS, JGAS ',
     &          IGAS, JGAS
              END IF
              CALL ORT_GAS_TO_GAS(IGAS,JGAS,WORK(KLSMO),WORK(KLCMOAO2),
     &             CMOAO_OUT)
C                 ORT_GAS_TO_GAS(IGAS,JGAS,SIN,CIN,COUT)
              CALL COPVEC(CMOAO_OUT,WORK(KLCMOAO2),LEN_1F)
*. Update metric
              CALL GET_SMO(CMOAO_OUT,WORK(KLSMO),IPACK_OUT)
            END DO
          END DO
        END IF ! End if intergas orthogonalization was called
*
        IF(NSCOB.NE.0) THEN
*
*. Orthogonalize Secondary space to GASpaces
*
          DO IGAS = 1, NGAS
*. Orthogonalize Secondary to GAS IGAS 
            CALL ORT_GAS_TO_GAS(IGAS,NGAS+1,WORK(KLSMO),WORK(KLCMOAO2),
     &           CMOAO_OUT)
            CALL COPVEC(CMOAO_OUT,WORK(KLCMOAO2),LEN_1F)
*. Update metric
            CALL GET_SMO(CMOAO_OUT,WORK(KLSMO),IPACK_OUT)
          END DO
        END IF
      END IF ! End if interspace orthogonalization was called
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' MOAO transformation matrix after INTERORT'
        CALL APRBLM2(CMOAO_OUT,NTOOBS,NTOOBS,NSMB,0)
      END IF

*
*
* ==============================
* The intragas orthogonalization
* ==============================
*
*
      IF(NTEST.GE.1000)
     &WRITE(6,*) ' INTRAGAS_ORT after INTERORT', INTRAGAS_ORT
      IF(INTRAGAS_ORT .NE.0) THEN
*. Space for Metric in MO basis MO-MO transformation, blocks of S and C,
*. and scratch
       CALL MEMMAN(KLSMO,LEN_1F,'ADDL  ',2,'SAOE  ')
       CALL MEMMAN(KLCMOMO,LEN_1F,'ADDL  ',2,'CMOMO ')
       CALL MEMMAN(KLSBLK,MXTOB**2,'ADDL  ',2,'SBLK  ')
       CALL MEMMAN(KLCBLK,MXTOB**2,'ADDL  ',2,'CBLK  ')
       LSCR = 2*LEN_1F + 6*MXTOB**2
       CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'SCRORT')
*. Initialize MOMO- transformation matrix to 1
       ZERO = 0.0D0
       CALL SETVEC(WORK(KLCMOMO),ZERO,LEN_1F)
       ONE = 1.0D0
       CALL SETDIA_BLM(WORK(KLCMOMO),ONE,NSMOB,NTOOBS,0)
C           SETDIA_BLM(B,VAL,NBLK,LBLK,IPCK)
*. Obtain metric in MO basis
       IPACK_OUT = 0
       CALL GET_SMO(CMOAO_OUT,WORK(KLSMO),IPACK_OUT)
*. Loop over gas-spaces
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Information from GAS-GAS orthog '
         WRITE(6,*) ' Metric in MO basis after INTERGAS part'
         CALL APRBLM2(WORK(KLSMO),NTOOBS,NTOOBS,NSMOB,0)
       END IF
       DO IGAS = 0, NGAS+1
*. Number of orbitals per sym of this GASpace
        CALL EXTRROW(NOBPTS_GN,IGAS+1,7+MXPR4T,NSMOB,IDIM)
        NTOB = IELSUM(IDIM,NSMOB)
        IF(NTOB.NE.0) THEN
         IF(NTEST.GE.1000)
     &   WRITE(6,*) ' Orthonormalization of GAS space = ', IGAS
*. Extract block (IGAS,IGAS) of S-matrix and save in KLSBLK
C             EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT(A,AGAS,IGAS,JGAS,I_EX_OR_CP)
         CALL EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT
     &        (WORK(KLSMO),WORK(KLSBLK),IGAS,IGAS,1)
*. And obtain transformation matrix giving othogonal basis  
*. Orthogonalization method defined differently in ORTHGNORM..
         IF(IGAS.NE.NORTCIX_SCVB_SPACE) THEN
           IORTMET_L = INTRAGAS_ORT
         ELSE
           IF(IORT_VB.EQ.0) THEN
             IORTMET_L = 0
           ELSE
             IORTMET_L = INTRAGAS_ORT
           END IF
         END IF
*
C?       WRITE(6,*) ' IORTMET_L = ', IORTMET_L
         IF(IORTMET_L.NE.0) THEN
          CALL ORTHNORM_BLKMT(WORK(KLSBLK),WORK(KLCBLK),NSMOB,IDIM,
     &         WORK(KLSCR),IORTMET_L)
C              ORTHNORM_BLKMT(S,C,NBLK,LBLK,SCR,IORTMET)
*. Copy transformation matrix to complete matrix
C              EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT(A,AGAS,IGAS,JGAS,I_EX_OR_CP)
          CALL EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT
     &         (WORK(KLCMOMO),WORK(KLCBLK),IGAS,IGAS,2)
         END IF ! IORTMET_L .ne. 0
         END IF ! transformation should be done
       END DO ! loop over GASpaces
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) 
     &  ' Intra-gas  MO-MO transformation matrix '
         CALL APRBLM2(WORK(KLCMOMO),NTOOBS,NTOOBS,NSMOB,0)
       END IF
* CMOAO_OUT = "CMOAO_IN " * CMOMO 
C           MULT_BLOC_MAT
C           (C,A,B,NBLOCK,LCROW,LCCOL,LAROW,LACOL,LBROW,LBCOL,ITRNSP)
       CALL COPVEC(CMOAO_OUT,WORK(KLCMOAO2),LEN_1F)
       CALL MULT_BLOC_MAT(CMOAO_OUT,WORK(KLCMOAO2),WORK(KLCMOMO),
     &      NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
      END IF ! End if intragas orthogonalization was required.
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MO-AO transformation matrix '
        CALL APRBLM2(CMOAO_OUT,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ORTOBV')
*
      RETURN
      END
      SUBROUTINE GET_SMO(CMO,SMO,IPACK_OUT)
*
*. Obtain Metric, SMO, over a set of orbitals, CMO.
*. Metric is given in packed form if IPACK_OUT = 1
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      DIMENSION CMO(*)
*. Output
      DIMENSION SMO(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from GET_SMO'
        WRITE(6,*) ' ================='
        WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input CMO basis '
        CALL APRBLM2(CMO,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GETSMO')
*
C             NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
      LEN_M = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLSAOE,LEN_M,'ADDL  ',2,'SAO_E ')
      CALL MEMMAN(KLSCR,2*LEN_M,'ADDL  ',2,'SCR   ')
*. Expand SAO
      CALL TRIPAK_AO_MAT(dbl_mb(KLSAOE),dbl_mb(KSAO),2)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Expanded SAO '
        CALL APRBLM2(dbl_mb(KLSAOE),NTOOBS,NTOOBS,NSMOB,0) 
      END IF
*. Obtain Metric in MO-basis, SMO  = CMO(T) SAO CMO
C          TRAN_SYM_BLOC_MAT4(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      CALL TRAN_SYM_BLOC_MAT4(dbl_mb(KLSAOE),CMO,CMO,
     &     NSMOB,NTOOBS,NTOOBS,SMO,dbl_mb(KLSCR),0)
*
      IF(IPACK_OUT.EQ.1) THEN
*. Pack output matrix
        CALL COPVEC(SMO,dbl_mb(KLSAOE),LEN_M)
        CALL TRIPAK_AO_MAT(dbl_mb(KLSAOE),SMO,1)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Metric in MO basis '
        CALL APRBLM2(SMO,NTOOBS,NTOOBS,NSMOB,IPACK_OUT)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GETSMO')
*
      RETURN
      END
      SUBROUTINE ORTHNORM_BLKMT(S,C,NBLK,LBLK,SCR,IORTMET)
*
* Obtain transformation matrix that orthonormalizes basis 
* defined by blocked metric S
*
* IMET = 1:  Symmetric orthonormalization
* IMET = 2:  Orthonormalize by diagonalization 
* 
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input: S is given in packed form
      INTEGER LBLK(NBLK)
      DIMENSION S(*)
*. Output
      DIMENSION C(*)
*. Scratch: Should atleast be: 2* Dimension of matrix + 6 times largest block
      DIMENSION SCR(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ORTHNORM_BLKMT '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' Number of elements per block '
        CALL IWRTMA(LBLK,1,NBLK,1,NBLK)
        IF(IORTMET.EQ.1) THEN
          WRITE(6,*) ' Symmetric orthogonalization '
        ELSE IF (IORTMET.EQ.2) THEN
          WRITE(6,*) ' Orthonormalization by diagonalization of metric'
        END IF
        WRITE(6,*) ' (IORTMET = ', IORTMET
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input metric: '
        CALL APRBLM2(S,LBLK,LBLK,NBLK,0)
      END IF
*
      LEN_MAT = LEN_BLMAT(NBLK,LBLK,LBLK,0)
      IF(IORTMET.EQ.1) THEN
*. Obtain S ** (-1/2)
        KLSQRT = 1
        KLSCR = KLSQRT + LEN_MAT
        CALL SQRT_BLMAT(S,NBLK,LBLK,2,SCR(1),C,SCR(KLSCR),0)
C            SQRT_BLMAT(A,NBLK,LBLK,ITASK,ASQRT,AMSQRT,SCR,ISYM)
      ELSE
        CALL GET_ON_BASIS_BY_DIAG_BLKMT(S,NBLK,LBLK,C,SCR,1)
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) 
     &  ' ORTHNORM_BLKMT: Matrix defining orthonormal basis '
        WRITE(6,*) 
     &  ' ==================================================='
        WRITE(6,*)
        CALL APRBLM2(C,LBLK,LBLK,NBLK,0)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_ON_BASIS_BY_DIAG_BLKMT(S,NBLK,LBLK,C,SCR,IPACK)
*
* A blocked metric S is given (lower half packed if IPACK = 1)
* Obtain block form of transformation matrix giving the orthonormal basis
* that is obtained by diagonalization
* S = U(T) Sigma  U, C = U  Sigma**(-1/2)
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LBLK(NBLK)
      DIMENSION S(*)
*. Output
      DIMENSION C(*)
*. Scratch:  Should at least be of length L**2 + 2L, where L is dimension
*            of largest block
      DIMENSION SCR(*)
*
      NTEST = 0
*
      NSING = 0
*. To get rid of compiler warninf
      IOFF = 0
      DO IBLK = 1, NBLK
        IF(IBLK.EQ.1) THEN
          IOFF = 1
          IOFFS = 1
        ELSE
          IOFF = IOFF + LBLK(IBLK-1)**2
          IF(IPACK.EQ.0) THEN 
            IOFFS = IOFF
          ELSE
            IOFFS = IOFFS + LBLK(IBLK-1)*(LBLK(IBLK-1)-1)/2
          END IF
        END IF
*
        KLS = 1
        KLVEC1 = KLS + LBLK(IBLK)**2
        KLVEC2=  KLVEC1 + LBLK(IBLK)
*. Obtain unpacked, but blocked, matrix in SCR(KLS)
        IF(IPACK.EQ.0) THEN
          LL = LBLK(IBLK)**2
          CALL COPVEC(S(IOFFS),SCR(KLS),LL)
        ELSE
          CALL TRIPAK_BLKM(SCR(KLS),S,2,LBLK,NBLK)
        END IF
*. And obtain orthonormal basis
        THRES_SINGU = 1.0D-14
C            GET_ON_BASIS2(S,NVEC,NSING,X,SCRVEC1,SCRVEC2,THRES_SINGU)
        CALL GET_ON_BASIS2(SCR(KLS),LBLK(IBLK),NSING_BLK,C(IOFF),
     &                     SCR(KLVEC1),SCR(KLVEC2), THRES_SINGU)
        NSING = NSING + NSING_BLK 
        IF(NSING_BLK.NE.0) THEN
          WRITE(6,*) ' Singularities in metric block ', IBLK
          WRITE(6,*) ' Number of singularities ',       NSING_BLK
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orthonormalization matrix from diagonalization'
        CALL APRBLM2(C,LBLK,LBLK,NBLK,0)
      END IF
*
      IF(NSING.NE.0) THEN 
        WRITE(6,*) ' Singularities in metric '
        WRITE(6,*) ' Number of singularities in metric ', NSING
        STOP       ' Singularities in metric '
      END IF
*
      RETURN
      END
      FUNCTION LEN_BLMAT(NBLK,LROW,LCOL,IPACK)
*
* Determine number of elements in packed matrix with NBLK blocks
* with dimensions LROW, LCOL.
* IPACK = 1 => matrix is packed
*
* Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LROW(NBLK),LCOL(NBLK)
*
      LEN = 0
      IF(IPACK.EQ.0) THEN
        DO IBLK = 1, NBLK
          LEN = LEN + LROW(IBLK)*LCOL(IBLK)
        END DO
      ELSE
        DO IBLK = 1, NBLK
          LEN = LEN + LROW(IBLK)*(LROW(IBLK)+1)/2
        END DO
      END IF
*
      LEN_BLMAT = LEN
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Dimension of block matrix ', LEN
      END IF
*
      RETURN
      END
      SUBROUTINE ORT_GAS_TO_GAS(IGAS,JGAS,SIN,CIN,COUT)
*
* Orthogonalize Orbitals in space JGAS to orbitals in space IGAS, i.e.
* modify orbitals in space JGAS so they are orthogonal to 
* orbitals in space IGAS
*
* Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc' 
      INCLUDE 'wrkspc-static.inc'
*. Specific Input: Expansion of input MO's in AO's 
      DIMENSION CIN(*),SIN(*)
*. Output: Expansion of output MO's in AO's
      DIMENSION COUT(*)
*. Local scratch
      INTEGER IDIM(MXPOBS),JDIM(MXPOBS)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ORT_GAS_TO_GAS'
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' IGAS, JGAS = ', IGAS, JGAS
      END IF
      IF(NTEST.GE.10000) THEN 
        WRITE(6,*) ' CIN entering ORT_GAS_TO_GAS '
        CALL APRBLM2(CIN,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'ORTGAS')
*. A bit of scratch
      LSCR = 2 * MXTOB **2
C?    WRITE(6,*) ' Test: MXTOB = ',MXTOB
      CALL MEMMAN(KLSCR,LSCR, 'ADDL  ', 2, 'SCRORT')
*
      CALL MEMMAN(KLSII, MXTOB**2, 'ADDL  ',2,'SJJ   ')
      CALL MEMMAN(KLSIJ, MXTOB**2, 'ADDL  ',2,'SJI   ')
      CALL MEMMAN(KLC, MXTOB**2, 'ADDL  ',2,'SJI   ')
*
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
      LSCR = MXTOB*MXSOB
      CALL MEMMAN(KLCI, LSCR, 'ADDL  ',2,'CIMOAO')
      CALL MEMMAN(KLCJ, LSCR, 'ADDL  ',2,'CJMOAO')
      CALL MEMMAN(KLCJT, LSCR, 'ADDL  ',2,'CJMOAO')
*
*. Dimensions of IGAS, JGAS (over symmetries)
*
      CALL EXTRROW(NOBPTS_GN,IGAS+1,7+MXPR4T,NSMOB,IDIM)
      CALL EXTRROW(NOBPTS_GN,JGAS+1,7+MXPR4T,NSMOB,JDIM)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Number of orbitals per sym in IGAS = ', IGAS
        CALL IWRTMA(IDIM,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' Number of orbitals per sym in JGAS = ', JGAS
        CALL IWRTMA(JDIM,1,NSMOB,1,NSMOB)
      END IF
*
*. Extract S(IGAS,IGAS),S(IGAS,JGAS)
*
C     EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT(A,AGAS,IGAS,JGAS,I_EX_OR_CP)
      CALL EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT(SIN,WORK(KLSII),IGAS,IGAS,1)
      CALL EXTR_OR_CP_GAS_BLKS_FROM_ORBMAT(SIN,WORK(KLSIJ),IGAS,JGAS,1)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' S(IGAS,IGAS) for IGAS = ', IGAS
        CALL APRBLM2(WORK(KLSII),IDIM,IDIM,NSMOB,0)
        WRITE(6,*) ' S(IGAS,JGAS) for IGAS, JGAS = ', IGAS, JGAS
        CALL APRBLM2(WORK(KLSIJ),IDIM,JDIM,NSMOB,0)
      END IF
*
*. Obtain coefficient matrix of I-vectors to obtain orthogonality
*
C     ORT_SPCY_TO_SPCX_BLK(NX,NY,NBLK,SXX,SXY,C,SCR)
      CALL ORT_SPCY_TO_SPCX_BLK(IDIM,JDIM,NSMOB,
     &     WORK(KLSII), WORK(KLSIJ),WORK(KLC),WORK(KLSCR))
*
*. Obtain MO-orbitals of space I and J
*
      CALL EX_OR_CP_MO_FOR_GAS(CIN,WORK(KLCI),IGAS,1)
      CALL EX_OR_CP_MO_FOR_GAS(CIN,WORK(KLCJ),JGAS,1)
*
*. Update MO- coefficients for JGAS
*
C           MULT_BLOC_MAT
C           (C,A,B,NBLOCK,LCROW,LCCOL,LAROW,LACOL,LBROW,LBCOL,ITRNSP)
* C(JGAS) =  C(IGAS)*C
      CALL MULT_BLOC_MAT(WORK(KLCJT),WORK(KLCI),WORK(KLC),NSMOB,
     &     NTOOBS,JDIM,NTOOBS,IDIM,IDIM,JDIM,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Correction to Y_j = sum_k X_k C(k,j) '
        CALL APRBLM2(WORK(KLCJT),NTOOBS,JDIM,NSMOB,0)
      END IF
      LEN = LEN_BLMAT(NSMOB,JDIM,NTOOBS,0)
C     LEN_BLMAT(NBLK,LROW,LCOL,IPACK)
      ONE = 1.0D0
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input block CIN for JGAS = ', JGAS
        CALL APRBLM2(WORK(KLCJ),NTOOBS,JDIM,NSMOB,0)
      END IF
      CALL VECSUM(WORK(KLCJ),WORK(KLCJ),WORK(KLCJT),ONE,ONE,LEN)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Updated matrix C(JGAS) '
        CALL APRBLM2(WORK(KLCJ),NTOOBS,JDIM,NSMOB,0)
      END IF
*
*. And transfer to COUT
*
      LEN_1F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL COPVEC(CIN,COUT,LEN_1F)
      CALL EX_OR_CP_MO_FOR_GAS(COUT,WORK(KLCJ),JGAS,2)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' MO expansion after orthogonalization of GAS ', JGAS , ' TO ',
     &    IGAS
        CALL APRBLM2(COUT,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'ORTGAS')
*
      RETURN
      END
      SUBROUTINE ORT_SPCY_TO_SPCX_BLK(NX,NY,NBLK,SXX,SXY,C,SCR)
* A space X with metric SXX 
* and a space Y with overlap SXY with X
* is given.  The space and metrics are divided into NBLK blocks.
*
* Obtain the matrix C so (Y_ai + sum_k C_ki X_ak) is 
* orthogonal to space X
*
* C(IBLK) = -SXX(IBLK)(-1) SXY(IBLK)
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NX(NBLK),NY(NBLK)
      DIMENSION SXX(*),SXY(*)
*. Output
      DIMENSION C(*)
*. Scratch: Should be length 2*NXM*NXM where NX is dim of largest block
      DIMENSION SCR(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from ORT_SPCY_TO_SPCX_BLK '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input matrix SXX '
        CALL APRBLM2(SXX,NX,NX,NBLK,0)
        WRITE(6,*) ' Input matrix SXY '
        CALL APRBLM2(SXY,NX,NY,NBLK,0)
      END IF
*
      NXM = IMNMX(NX,NBLK,2)
      KLS = 1
      KLSCR = KLS + NXM*NXM
*
      DO IBLK = 1, NBLK
        IF(IBLK.EQ.1) THEN
         IOFFXX = 1
         IOFFXY = 1
        ELSE
         IOFFXX = IOFFXX + NX(IBLK-1)**2
         IOFFXY = IOFFXY + NX(IBLK-1)*NY(IBLK-1)
        END IF
        NNX = NX(IBLK)
        NNY = NY(IBLK)
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' IBLK, NNX, NNY = ', IBLK, NNX, NNY
*. Obtain SXX(IBLK)  (-1)
        CALL COPVEC(SXX(IOFFXX),SCR(KLS),NNX**2)
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' BLOCK SXX: '
         CALL WRTMAT(SCR(KLS),NNX,NNX,NNX,NNX)
        END IF
        ISING = 0
        CALL INVMAT(SCR(KLS),SCR(KLSCR),NNX,NNX,ISING)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Inverted SXX block'
          CALL WRTMAT(SCR(KLS),NNX,NNX,NNX,NNX)
        END IF
        IF(ISING.GT.0) THEN
         WRITE(6,*) ' Problem inverting  SXX '
        END IF
*. And multiply
C         MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &           NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
        FACTORC = 0.0D0
        FACTORAB = -1.0D0
        CALL MATML7(C(IOFFXY),SCR(KLS),SXY(IOFFXY),
     &       NNX,NNY,NNX,NNX,NNX,NNY,FACTORC, FACTORAB,0)
      END DO! End of loop over blocks
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' C matrix for space-space orthogonalization '
       CALL APRBLM2(C,NX,NY,NBLK,0)
      END IF
*
      RETURN
      END
C     CALL EX_OR_CP_MO_FOR_GAS(CMO,WORK(KLCI),IGAS,I_EX_OR_CP)
      SUBROUTINE EX_OR_CP_MO_FOR_GAS(CMO_TOT, CMO_GAS, IGAS,
     &            I_EX_OR_CP)
*
* Extract from or copy to CMO_TOT orbitals belonging to GASpace IGAS  
* to/from  CMO_GAS
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      DIMENSION CMO_TOT(*)
*. Output
      DIMENSION CMO_GAS(*)
*. Local scratch
      DIMENSION IDIM(MXPOBS)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from EX_OR_CP_MO_FOR_GAS'
        WRITE(6,*) ' ============================='
        WRITE(6,*) ' I_EX_OR_CP = ', I_EX_OR_CP
        WRITE(6,*) ' IGAS = ',IGAS
      END IF
      IF(NTEST.GE.10000) THEN 
        WRITE(6,*) ' Input complete matrix '
        CALL APRBLM2(CMO_TOT,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      DO ISM = 1, NSMOB
*. Start of symmetry block
       IF(ISM .EQ. 1 ) THEN
         IOFF_TOT = 1
         IOFF_GAS = 1
       ELSE
         IOFF_TOT = IOFF_TOT + NTOOBS(ISM-1)**2
         IOFF_GAS = IOFF_GAS + NTOOBS(ISM-1)*NOBPTS_GN(IGAS,ISM-1)
       END IF
*. First orbital in GASpace IGAS in sym ISM - relative to start of sym
       IOFF_REL = 1
       DO JGAS = 0, IGAS-1
         IOFF_REL = IOFF_REL + NOBPTS_GN(JGAS,ISM)
       END DO
       NOB_GAS = NOBPTS_GN(IGAS,ISM)
       NOB_SM  = NTOOBS(ISM)
       IF(I_EX_OR_CP.EQ.1) THEN
         CALL COPVEC(CMO_TOT(IOFF_TOT-1+(IOFF_REL-1)*NOB_SM+1),
     &               CMO_GAS(IOFF_GAS),NOB_GAS*NOB_SM)
       ELSE
         CALL COPVEC(CMO_GAS(IOFF_GAS),
     &               CMO_TOT(IOFF_TOT-1+(IOFF_REL-1)*NOB_SM+1),
     &               NOB_GAS*NOB_SM)
       END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        CALL EXTRROW(NOBPTS_GN,IGAS+1,7+MXPR4T,NSMOB,IDIM)
        IF(I_EX_OR_CP.EQ.1) THEN
         WRITE(6,*) ' Extracted MO coefficients for IGAS = ', IGAS
         CALL APRBLM2(CMO_GAS,NTOOBS,IDIM,NSMOB,0)
        ELSE
         WRITE(6,*) ' Updated matrix of MO coefficients ' 
         CALL APRBLM2(CMO_TOT,NTOOBS,NTOOBS,NSMOB,0)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE INI_CSFEXP(CINI)
*
* Obtain initial CI expansion in terms of the CSF expansion CINI
* Configuration space is specified as ICSPC_CN
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc' 
      INCLUDE 'spinfo.inc'
*. Output
      DIMENSION CINI(*)
*
*. If an initial configuration has been specified use thus
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'INICSF')
*
      NTEST = 1000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' INI_CSFEXP reporting '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' I_HAVE_INI_CONF = ', I_HAVE_INI_CONF
      END IF
*
      NCSF_TOT = NCSF_PER_SYM_GN(ICSM, ICSPC_CN)
*. Initialize by zero
      ZERO = 0.0D0
      CALL SETVEC(CINI,ZERO,NCSF_TOT)
*
      IF(I_HAVE_INI_CONF.EQ.1) THEN
        WRITE(6,*) ' Initial configuration used as initial guess '
*
*. Find address of configuration
C              ILEX_FOR_CONF_G(ICONF,NOCC_ORB,ICONF_SPC,IDOREO)
        ILEX = ILEX_FOR_CONF_G(INI_CONF,NOB_INI_CONF,ICSPC_CN,1)
        IF(NTEST.GE.1000) WRITE(6,*) ' Address of config = ', ILEX
*. Number of CSF's for this configuration
        IOPEN = 2*NOB_INI_CONF-N_EL_CONF
*. Address in CSFVEC of first CSF with this number of open orbitals
        IB_OPEN = IB_OPEN_CSF(IOPEN+1,ICSM,ICSPC_CN)
        IF(NTEST.GE.1000) WRITE(6,*) ' IB_OPEN = ', IB_OPEN
*. Address of first configuration with this number of open orbitals
        IB_CONF = IB_CONF_REO_GN(IOPEN+1,ICSM,ICSPC_CN)
        IF(NTEST.GE.1000) WRITE(6,*) ' IB_CONF = ', IB_CONF
*. Address of first CSF belonging to this configuration
        IADDR = IB_OPEN + (ILEX-IB_CONF)*NPCSCNF(IOPEN+1)
        IF(NTEST.GE.1000) WRITE(6,*) ' IADDR = ', IADDR
   
*. Equal contribution to all CSF's of config
        NCSF_CONF = NPCSCNF(IOPEN+1)
        XNCSF_CONF = DBLE(NCSF_CONF)
C?      WRITE(6,*) ' IOPEN, NPCSCNF(IOPEN+1) = ',
C?   &               IOPEN, NPCSCNF(IOPEN+1)
        FACTOR = 1.0D0/SQRT(XNCSF_CONF)
C?      WRITE(6,*) ' XNCSF_CONF, FACTOR = ', 
C?   &              XNCSF_CONF, FACTOR
        CALL SETVEC(CINI(IADDR),FACTOR,NCSF_CONF)
      ELSE
*. Set configuration one to 1
       IF(NCSF_TOT.GE.7) THEN
         CINI(7) = 1.0D0
         WRITE(6,*) ' Initial guess set to CSF 7 !!!! '
       ELSE
         CINI(1) = 1.0D0 
         WRITE(6,*) ' Initial guess set to CSF 1  '
       END IF
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial CI vector '
        WRITE(6,*) ' ================='
        CALL WRTMAT(CINI,1,NCSF_TOT,1,NCSF_TOT)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'INICSF')
*
      RETURN
      END
      FUNCTION NEL_IN_COMPACT_CONF(ICONF,NOCOB)
*
*. Number of electrons in configuration, compact form
*
* Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INTEGER ICONF(NOCOB)
*
      NEL = 0
      DO IORB = 1, NOCOB
       IF(ICONF(IORB).GT.0) THEN
         NEL = NEL + 1
       ELSE
         NEL = NEL + 2
       END IF
      END DO
*
      NEL_IN_COMPACT_CONF = NEL
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from NEL_IN_COMPACT_CONF '
        WRITE(6,*) ' Configuration: '
        CALL IWRTMA(ICONF,1,NOCOB,1,NOCOB)
        WRITE(6,*) ' Number of electrons = ', NEL
      END IF
*
      RETURN
      END
      SUBROUTINE SIGMA_CONF(C,HC,LUC,LUHC)
*
* Configuration driven Sigma routine
* Jeppe Olsen, July 2011
*
*. The input and output CI spaces in action are defined by the 
* ICPSC_CN, ISSPC_CN parameters in cands
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'spinfo.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGCNF')
*
      NTEST = 1000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from SIGMA_CONF '
        WRITE(6,*) ' ===================== '
        WRITE(6,'(A,2I3)') ' Config space and sym for C ', ICSPC_CN,ICSM
        WRITE(6,'(A,2I3)') ' Config space and sym for S ', ISSPC_CN,ISSM
      END IF
*. 
      NCONF_C = NCONF_PER_SYM_GN(ICSM,ICSPC_CN)
      NSD_C   = NSD_PER_SYM_GN(ICSM,ICSPC_CN)
      NCSF_C = NCSF_PER_SYM_GN(ICSM,ICSPC_CN)
*
      NCONF_S = NCONF_PER_SYM_GN(ISSM,ISSPC_CN)
      NSD_S   = NSD_PER_SYM_GN(ISSM,ISSPC_CN)
      NCSF_S = NCSF_PER_SYM_GN(ISSM,ISSPC_CN)
*
      NCONF_MAX = MAX(NCONF_C,NCONF_S)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,3I8)') ' Number of confs, SDs and CSFs for C ',
     &  NCONF_C, NSD_C, NCSF_C
        WRITE(6,'(A,3I8)') ' Number of confs, SDs and CSFs for S ',
     &  NCONF_S, NSD_S, NCSF_S
      END IF
*
*. Number of batches for configuration expansions (each batch atmost dim LCSBLK)
* ================================================================================
*
*. Allowed length of batch:
* ==========================
*. IF LCSBLK has not been specified, a default batch size is used
*
      LCSBLK_L = LCSBLK
      IF(LCSBLK_L.LE.0) THEN
        WRITE(6,*) ' SIGMA_CONF will define length of batch '
        LCSBLK_DEFAULT = 2000000
*. Compare with dimension of largest single configuration
        LCONF_MAX = IMNMX(NPCSCNF,MAXOP+1,2)
        IF(LCONF_MAX.GT.LCSBLK_DEFAULT) LCSBLK_DEFAULT = LCONF_MAX
        LCSBLK_L = LCSBLK_DEFAULT
      END IF
*. If ICISTR = 1, vectors are stored in one batch, so
      IF(ICISTR.EQ.1) LCSBLK_L = MAX(NCSF_S,NCSF_S)
      
      IF(NTEST.GE.1000) WRITE(6,*) ' Allowed size of batch ', LCSBLK_L
*. Batches of C
*. ==============
*. One could here either use CSF's or SD's. As memory maybe the defining parameter,
* I opt for CSF's and will then expand/contract each configuration when needed.
*. Length of each configuration
      CALL MEMMAN(KLLCNFEXP,NCONF_MAX,'ADDL  ',1,'LCNFEX')
*. For C
C     CONF_EXP_LEN_LIST(ILEN,NCONF_PER_OPEN,NELMNT_PER_OPEN,MAXOP)
      CALL CONF_EXP_LEN_LIST(WORK(KLLCNFEXP),
     &     NCONF_PER_OPEN_GN(1,ICSM,ICSPC_CN),NPCSCNF,MAXOP)
C     PART_VEC(LBLK,NBLK,MAXSTR,LBAT,NBAT,IONLY_NBAT)
      CALL PART_VEC(WORK(KLLCNFEXP),NCONF_C,LCSBLK_L,IDUM,NBAT_C,1)
      CALL MEMMAN(KLLBAT_C,NBAT_C,'ADDL  ',1,'LBAT_C')
      CALL PART_VEC(WORK(KLLCNFEXP),NCONF_C,LCSBLK_L,WORK(KLLBAT_C),
     &     NBAT_C,0)
*. And for Sigma
C     CONF_EXP_LEN_LIST(ILEN,NCONF_PER_OPEN,NELMNT_PER_OPEN,MAXOP)
      CALL CONF_EXP_LEN_LIST(WORK(KLLCNFEXP),
     &     NCONF_PER_OPEN_GN(1,ISSM,ISSPC_CN),NPCSCNF,MAXOP)
C     PART_VEC(LBLK,NBLK,MAXSTR,LBAT,NBAT,IONLY_NBAT)
      CALL PART_VEC(WORK(KLLCNFEXP),NCONF_S,LCSBLK_L,IDUM,NBAT_S,1)
      CALL MEMMAN(KLLBAT_S,NBAT_S,'ADDL  ',1,'LBAT_S')
      CALL PART_VEC(WORK(KLLCNFEXP),NCONF_S,LCSBLK_L,WORK(KLLBAT_S),
     &     NBAT_S,0)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Number of batches for C and S ', NBAT_C, NBAT_S
      END IF
*. Largest number of configurations in a given batch
      MAX_CONF_BATCH_C = IMNMX(WORK(KLLBAT_C),NBAT_C,2)
      MAX_CONF_BATCH_S = IMNMX(WORK(KLLBAT_S),NBAT_S,2)
      MAX_CONF_BATCH = MAX(MAX_CONF_BATCH_C,MAX_CONF_BATCH_S)
*
      IF(NTEST.GE.100)
     &WRITE(6,*) ' Largest number of configs in batch ', MAX_CONF_BATCH
      CALL MEMMAN(KLLBLK_BAT_C,MAX_CONF_BATCH ,'ADDL  ',2,'LBLBTC')
      CALL MEMMAN(KLLBLK_BAT_S,MAX_CONF_BATCH ,'ADDL  ',2,'LBLBTS')
*. Two vectors for holding expansion in SD of given config
      LEN_SD_CONF_MAX = IMNMX(NPDTCNF,MAXOP+1,2)
      CALL MEMMAN(KLCONF_SD_C,LEN_SD_CONF_MAX,'ADDL  ',2,'CN_SDC')
      CALL MEMMAN(KLCONF_SD_S,LEN_SD_CONF_MAX,'ADDL  ',2,'CN_SDS')
*. Scratch space in routine for evuluating H for configurations (allowing combs)
*. Scratch: Length: INTEGER: (NDET_C + NDET_S)*N_EL_CONF + NDET_C + 6*NORB
      L_CNHCN = LEN_SD_CONF_MAX*(1+2*N_EL_CONF) + 6*N_ORB_CONF
      CALL MEMMAN(KL_CNHCN, L_CNHCN,'ADDL  ',1,'LCNHCN')
*. Space for two integers arrays for signs
      CALL MEMMAN(KLISIGNC,LEN_SD_CONF_MAX,'ADDL  ',1,'ISIGNC')
      CALL MEMMAN(KLISIGNS,LEN_SD_CONF_MAX,'ADDL  ',1,'ISIGNS')
*
C?    WRITE(6,*) ' KDFTP, KL_CNHCN = ', KDFTP, KL_CNHCN
C?    WRITE(6,*) ' KLLBLK_BAT_C, KLLBLK_BAT_S = ',
C?   &            KLLBLK_BAT_C, KLLBLK_BAT_S
C?    WRITE(6,*) ' KLCONF_SD_C, KLCONF_SD_S = ',
C?   &             KLCONF_SD_C, KLCONF_SD_S
C?    WRITE(6,*) ' KLLBAT_C, KLLBAT_S = ',
C?   &             KLLBAT_C, KLLBAT_S
*
      IADOB = IB_ORB_CONF - 1
      CALL SIGMA_CONF_SLAVE(C,HC,LUC,LUHC,ICISTR,
     &     NCONF_PER_OPEN_GN(1,ICSM,ICSPC_CN),
     &     NCONF_PER_OPEN_GN(1,ISSM,ISSPC_CN),
     &     NBAT_C,WORK(KLLBAT_C),
     &     NBAT_S,WORK(KLLBAT_S),
     &     WORK(KLLBLK_BAT_C),WORK(KLLBLK_BAT_S),
     &     WORK(KLCONF_SD_C),WORK(KLCONF_SD_S),
     &     IADOB,WORK(KDFTP),WORK(KL_CNHCN),
     &     WORK(KLISIGNC),WORK(KLISIGNS))
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SIGCNF')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Final sigma-vector from SIGMA_CONF'
        CALL WRTMAT(HC,1,NCSF_S,1,NCSF_S)
      END IF
       
      IF(NTEST.GE.1000) WRITE(6,*) ' SIGMA_CONF finished '
      RETURN
      END 
      SUBROUTINE SIGMA_CONF_SLAVE(C,S,LUC,LUS,ICISTR,
     &           NCONF_PER_OPEN_C,NCONF_PER_OPEN_S,
     &           NBAT_C,LBAT_C,NBAT_S,LBAT_S,
     &           LBLK_BAT_C,LBLK_BAT_S,
     &           CONF_SD_C,CONF_SD_S,IADOB,IPRODT,
     &           ISCR_CNHCN,ISIGN_C,ISIGN_S)
*
* Inner (aka slave) routine for direct CI in configuration based methodsø
*
*. Jeppe Olsen,July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cecore.inc'
*. Input
*. C-vector or space for batch of C-vector
      DIMENSION C(*)
*. Info on the two configuration expansions
       INTEGER NCONF_PER_OPEN_C(*), NCONF_PER_OPEN_S(*)
*. Number of blocks in the batches of C and S
      INTEGER LBAT_C(*), LBAT_S(*)
*. Scratch for Info on batches of C and S: Length of each block (configuration in batch)
      INTEGER LBLK_BAT_C(*),LBLK_BAT_S(*)
*. Space for SD expansion of single configurations
      DIMENSION CONF_SD_C(*), CONF_SD_S(*)
*. Space for signs for phase change for dets of a configurations
      INTEGER ISIGN_C(*),ISIGN_S(*)
*. CSF info: proto type dets 
      INTEGER IPRODT(*)

*. Output
      DIMENSION S(*)
*. Scratch transferred through to CNHCN
      INTEGER ISCR_CNHCN(*)
*. Local scratch
      INTEGER IOCC_C(MXPORB),IOCC_S(MXPORB)
*
      NTEST = 000
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGCNI')
*. Initialization of some parameters for controlling loop over configurations
      IOPEN_S = 0
      INUM_OPS = 0
      IOPEN_C = 0
      INUM_OPC = 0
      IB_CSF_C = 1
      IB_CSF_S = 1
*
      CALL MEMCHK2('INISIG')
*
*. Loop over batches of S
      INI_S = 1
      IF(NTEST.GE.1000) WRITE(6,*) ' NBAT_C, NBAT_S = ', 
     &NBAT_C, NBAT_S
      DO IBAT_S = 1, NBAT_S
       IF(NTEST.GE.1000) 
     & WRITE(6,'(A,I3)') ' >>> Start of sigma batch ', IBAT_S
*
       IF(IBAT_S.EQ.1) THEN
         IB_CONF_S = 1
         IB_CSF_S = 1
       ELSE
         IB_CONF_S = IB_CONF_S + LBAT_S(IBAT_S-1)
       END IF
C?     WRITE(6,*) ' LBAT_S(1) = ', LBAT_S(1)
       N_CONF_S = LBAT_S(IBAT_S) 
*. Number of CSF's per config in S-batch
C           GET_LBLK_CONF_BATCH(ICNF_INI,NCNF,LBLK_BAT,ISYM,ISPC,
C    &      NSD_BAT_TOT,NCSF_BAT_TOT) 
       CALL GET_LBLK_CONF_BATCH(IB_CONF_S,N_CONF_S,LBLK_BAT_S,ISSM,
     &      ISSPC_CN,NSD_BAT_TOT_S,NCSF_BAT_TOT_S)
       CALL MEMCHK2('AFGTL1')
       IF(NTEST.GE.100) THEN
         WRITE(6,'(A,2I9)') 
     &   ' Number of CSFs and SDs in S-batch ', NCSF_BAT_TOT_S,
     &     NSD_BAT_TOT_S
       END IF
*. Initialize sigma batch
       ZERO = 0.0D0
C?     WRITE(6,*) ' IB_CSF_S, NCSF_BAT_TOT_S = ',
C?   &              IB_CSF_S, NCSF_BAT_TOT_S
       CALL SETVEC(S(IB_CSF_S),ZERO,NCSF_BAT_TOT_S)
*. Loop over batches of C
C      IF(ICISTR.NE.1) REWIND LUHC
       INI_C = 1
*. First time in this batch
       ISBAT_FIRST_TIME =1
       DO IBAT_C = 1, NBAT_C
        IF(NTEST.GE.1000) 
     &  WRITE(6,'(A,I3)') ' >>> Start of C batch ', IBAT_C
        CALL MEMCHK2('STCBAT')
        IF(IBAT_C.EQ.1) THEN
          IB_CONF_C = 1
          IB_CSF_C = 1
        ELSE
          IB_CONF_C = IB_CONF_C + LBAT_C(IBAT_C-1)
        END IF
        N_CONF_C = LBAT_C(IBAT_C)
*. Number of configs per config in S-batch
        CALL GET_LBLK_CONF_BATCH(IB_CONF_C,N_CONF_C,LBLK_BAT_C,ICSM,
     &      ICSPC_CN,NSD_BAT_TOT_C,NCSF_BAT_TOT_C)
        IF(NTEST.GE.100) THEN
          WRITE(6,'(A,2I9)') 
     &   ' Number of CSFs and SDs in C-batch ', NCSF_BAT_TOT_C,
     &     NSD_BAT_TOT_C
        END IF
      CALL MEMCHK2('AFGTLB')
*. Read, if required, next batch of C- Each configuration stored in a record by itself
        IF(ICISTR.NE.1) THEN
          CALL FRMDSCN(C,N_CONF_C,-1,LUC)
C              FRMDSCN(VEC,NREC,LBLK,LU)
        END IF
*. And then to the configurations of the C and sigma
*. First time in this batch
        IF(ISBAT_FIRST_TIME.EQ.1) THEN
* Save pointers to start of configuration
          IOPEN_S_SAVE = IOPEN_S
          INUM_OPS_SAVE = INUM_OPS
          IB_CSF_S_SAVE = IB_CSF_S
          INI_S_SAVE = INI_S
        ELSE
          IOPEN_S = IOPEN_S_SAVE
          INUM_OPS = INUM_OPS
          IB_CSF_S = IB_CSF_S_SAVE
          INI_S = INI_S_SAVE
        END IF
        ISBAT_FIRST_TIME  = 0
*
        ICBAT_FIRST_TIME =1
        DO ICONF_S = IB_CONF_S, IB_CONF_S + N_CONF_S -1
*. Obtain occupation in IOCC_S and iopen for this sigma-configuration
C             NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
         IF(NTEST.GE.1000) WRITE(6,*) ' Requesting next S-conf: '
         CALL NEXT_CONF_IN_CONFSPC(IOCC_S,IOPEN_S,INUM_OPS,INI_S,
     &        ISSM,ISSPC_CN,NEW_S)
         INI_S = 0
         IOCOB_S = (IOPEN_S + N_EL_CONF)/2
*. Signs for going between configuration and interaction order of dets
C     SIGN_CONF_SD(ICONF,NOB_CONF,IOP,ISGN,IPDET_LIST,ISCR)
         CALL SIGN_CONF_SD(IOCC_S,IOCOB_S,IOPEN_S,ISIGN_S,IPRODT,
     &                      ISCR_CNHCN)
      CALL MEMCHK2('AFSIGN')
*
         IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Sigma configuration number ', ICONF_S
          CALL IWRTMA(IOCC_S,1,IOCOB_S,1,IOCOB_S)
         END IF
         NCSF_S = NPCSCNF(IOPEN_S+1)
         NSD_S = NPDTCNF(IOPEN_S+1)
*
         ZERO = 0.0D0
*. The contribution to a given sigma conf from all C-conf in C-batch
* will be stored in CONF_SD_S(1)
         CALL SETVEC(CONF_SD_S,ZERO,NSD_S)
*
         IF( ICBAT_FIRST_TIME .EQ. 1) THEN
           IOPEN_C_SAVE = IOPEN_C
           INUM_OPC_SAVE = INUM_OPC
           IB_CSF_C_SAVE = IB_CSF_C
         ELSE
           IOPEN_C = IOPEN_C_SAVE
           INUM_OPC = INUM_OPC_SAVE
           IB_CSF_C = IB_CSF_C_SAVE
         END IF
         ICBAT_FIRST_TIME = 0
         DO ICONF_C = IB_CONF_C, IB_CONF_C + N_CONF_C -1
*. Obtain occupation in IOCC_C and iopen for this C-configuration
C             NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
          IF(NTEST.GE.1000) WRITE(6,*) ' Requesting next C-conf: '
          CALL NEXT_CONF_IN_CONFSPC(IOCC_C,IOPEN_C,INUM_OPC,INI_C,
     &         ICSM,ICSPC_CN,NEW_C)
*. Signs for going between configuration and interaction order of dets
C     SIGN_CONF_SD(ICONF,NOB_CONF,IOP,ISGN,IPDET_LIST,ISCR)
          IOCOB_C = (IOPEN_C + N_EL_CONF)/2
          CALL SIGN_CONF_SD(IOCC_C,IOCOB_C,IOPEN_C,ISIGN_C,IPRODT,
     &                      ISCR_CNHCN)
          INI_C = 0
          IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' C configuration number ', ICONF_C
           IOCOB_C = (IOPEN_C + N_EL_CONF)/2
           CALL IWRTMA(IOCC_C,1,IOCOB_C,1,IOCOB_C)
          END IF
*. Expand coefficients for configuration from CSF to SD basis
          NCSF_C = NPCSCNF(IOPEN_C+1)
          NSD_C = NPDTCNF(IOPEN_C+1)
C              CSDTVC_CONF(C_SD,C_CSF,NOPEN,ISIGN,IAC,IWAY)
          CALL CSDTVC_CONF(CONF_SD_C,C(IB_CSF_C),IOPEN_C,ISIGN_C,2,1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' C(ICONF_C)  in SD'
            CALL WRTMAT(CONF_SD_C,1,NSD_C,1,NSD_C)
          END IF
          IF(NTEST.GE.1000)  THEN
            WRITE(6,'(A,2I6)') 
     &      ' Info on sigma for ICONF_C, ICONF_S = ',
     &      ICONF_C, ICONF_S
          END IF
*. Core energy is pt added in DIHDJ2, so the code below is outcommented
C!        IF(ICONF_C.EQ.ICONF_S) THEN
*. Add core energy
C!          ONE = 1.0D0
C!          CALL VECSUM(CONF_SD_S,CONF_SD_S,CONF_SD_C,
C!   &            ONE,ECORE,NSD_C)
C!        END IF
*. Update: S(I) = S(I) + Sum(J) <I!H!J> C(J)
C         CNHCN_LUCIA(ICNL,IOPL,ICNR,IOPR,C,SIGMA,
C    &                IADOB,IPRODT,I12OP,IORBTRA,IORB,IAB,ISCR)
          I12OP = 2
          I_DO_ORBTRA = 0
          IORB = 0
C     CNHCN_LUCIA(ICNL,IOPL,ICNR,IOPR,C,CNHCNM,SIGMA,
C    &           IADOB,IPRODT,I12OP,I_DO_ORBTRA,IORBTRA,
C    &           ECORE,ISCR)
          CALL CNHCN_LUCIA(IOCC_S,IOPEN_S,IOCC_C,IOPEN_C,
     &                     CONF_SD_C,XDUM,CONF_SD_S,IADOB,
     &                     IPRODT,I12OP,I_DO_ORBTRA,IORB,
     &                     ECORE,2,0,RJ,RK, ISCR_CNHCN)
C     CNHCN_LUCIA(ICNL,IOPL,ICNR,IOPR,C,CNHCNM,SIGMA,
C    &           IADOB,IPRODT,I12OP,I_DO_ORBTRA,IORBTRA,
C    &           ECORE,IHORS,ISYM,RJ,RK,ISCR)
*. Update address of C in action
          IB_CSF_C = IB_CSF_C + NCSF_C
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Updated Sigma(ICONF_S)  in SD'
            CALL WRTMAT(CONF_SD_S,1,NSD_S,1,NSD_S)
          END IF
         END DO ! over configs in batch of C
*. And transform sigma part to CSF and update sigma vector
         CALL CSDTVC_CONF(CONF_SD_S,S(IB_CSF_S),IOPEN_S,ISIGN_S,1,2)
         IB_CSF_S = IB_CSF_S + NCSF_S
         IF(NTEST.GE.1000) WRITE(6,*) ' End of conf for S-batch '
        END DO ! over configs in batch of S
        IF(NTEST.GE.1000) WRITE(6,*) ' End of C-batch '
       END DO ! Over batches of C
       IF(NTEST.GE.1000) WRITE(6,*) ' End of S-batch '
      END DO ! over batches of Sigma
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SIGCNI')
      IF(NTEST.GE.1000) WRITE(6,*) 'SIGMA_CONF_SLAVE finished'
      RETURN
      END
      SUBROUTINE PART_VEC(LBLK,NBLK,MAXSTR,LBAT,NBAT,IONLY_NBAT)
*
* A vector consists of NBLK BLocks with lengths  given by LBLK(IBLK).
* Partition the vector into batches of blocks, so each batch has atmost 
* length MAXSTR. 
* IF IONLY_NBAT = 1, then only the number of batches is calculated
*
* Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER LBLK(NBLK)
*.Output
      DIMENSION LBAT(*)
*
      NBAT = 1
      NBLK_B = 0
      LBLK_B = 0
      DO IBLK = 1, NBLK
       IF(LBLK(IBLK)+LBLK_B.LE.MAXSTR) THEN
*. Can be included in current batch
         NBLK_B = NBLK_B + 1
         LBLK_B = LBLK_B + LBLK(IBLK)
       ELSE
         IF(IONLY_NBAT.EQ.0) LBAT(NBAT) = NBLK_B
*. Start new batch
         NBAT = NBAT + 1
         LBLK_B = LBLK(IBLK)
         NBLK_B = 1
       END IF
      END DO
*. Save last Batch
      IF(IONLY_NBAT.EQ.0) LBAT(NBAT) = NBLK_B
      
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from PART_VEC '
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' Largest allowed batchsize ', MAXSTR
        WRITE(6,*) ' Number of batches ', NBAT  
        IF(IONLY_NBAT.EQ.0) THEN
          WRITE(6,*) ' Number of blocks in each batch '
          CALL IWRTMA(LBAT,1,NBAT,1,NBAT)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE MATVCC2(A,VIN,VOUT,NROW,NCOL,ITRNS,FACIN)
*
* ITRNS = 0 : VOUT(I) = FACIN*VOUT(I) + A(I,J)*VIN(J)
* ITRNS = 1 : VOUT(I) = FACIN*VOUT(I) + A(J,I)*VIN(J)
*
* NROW, NCOL are rows and column of input matrix (not transposed)
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(NROW,NCOL)
      DIMENSION VIN(*)
*. Output
      DIMENSION VOUT(*)
*
      IF(ITRNS.EQ.0) THEN
*
        IF(FACIN.EQ.0.0D0) THEN
          ZERO = 0.0D0
          CALL SETVEC(VOUT,ZERO,NROW)
        ELSE
          CALL SCALVE(VOUT,FACIN,NROW)
        END IF
*
        DO J = 1, NCOL
         VINJ = VIN(J)
         DO I = 1, NROW
           VOUT(I) = VOUT(I) + A(I,J)*VINJ
         END DO
        END DO
*
      ELSE IF( ITRNS.EQ.1) THEN
*
        DO I = 1, NCOL
          IF(FACIN.EQ.0.0D0) THEN
            X = 0.0D0
          ELSE
            X = FACIN*VOUT(I)
          END IF
*
          DO J = 1, NROW
            X = X + A(J,I)*VIN(J)
          END DO
          VOUT(I) = X
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(ITRNS.EQ.0) THEN
          WRITE(6,*) ' Vectorout = matrix * vectorin (MATVCC) '
          WRITE(6,*) ' Input and output vectors '
          CALL WRTMAT(VIN,1,NCOL,1,NCOL)
          CALL WRTMAT(VOUT,1,NROW,1,NROW)
          WRITE(6,*) ' Matrix '
          CALL WRTMAT(A,NROW,NCOL,NROW,NCOL)
        ELSE
          WRITE(6,*) ' Vectorout = matrix(T) * vectorin (MATVCC) '
          WRITE(6,*) ' Input and output vectors '
          CALL WRTMAT(VIN,1,NROW,1,NROW)
          CALL WRTMAT(VOUT,1,NCOL,1,NCOL)
          WRITE(6,*) ' Matrix (untransposed)'
          CALL WRTMAT(A,NROW,NCOL,NROW,NCOL)
        END IF
      END IF

*
      RETURN
      END
      SUBROUTINE ISIGN_TIMES_REAL(ISIGN,VEC,NDIM)
*
* VEC(I) = ISIGN(I)*VEC(I)
*
* X X 
*
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER ISIGN(*)
      DIMENSION VEC(*)
*
      DO I = 1, NDIM
        IF(ISIGN(I).EQ.-1) VEC(I) = -VEC(I)
      END DO
*. (No NTEST here, as it could identify programmer....)
      RETURN
      END
      SUBROUTINE MINMAX_FOR_ORBTRA(MIN_IN,MAX_IN,MIN_OUT,MAX_OUT,
     &           MIN_INTM,MAX_INTM,MIN_INTMS,MAX_INTMS,ISYM,IDODIM,
     &           NCONF_INTM,NCSF_INTM,NSD_INTM)
*
* Obtain intermediate MINMAX spaces for transforming between
* initial (MIN/MAX_IN) and final (MIN/MAX_OUT) spaces.
*
* Two intermediate spaces are produced
*
* _INTM:   Just overall occupations are considered
* _INTMS:  Also occupations in each orbital symmetry is 
*           considered
* (INTMS arrays not activated yet...)
* 
* IF IDODIM.EQ.1, the number of configs, CSF's and SD's 
* is calculated for the various spaces and SYM ISYM.
*
* IP_SPC is the first space in MIN
* 
*
* Jeppe Olsen, July 16 2011 (55 years birthday- still programming)
*
* No distinction is made here of the two operators used to
* transform a given orbital. The IORB array should be 
* used as final space for both operators for this orbital
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'spinfo.inc'
*. Input
      INTEGER MIN_IN(MXPORB), MAX_IN(MXPORB) 
      INTEGER MIN_OUT(MXPORB), MAX_OUT(MXPORB)
*. Output
      INTEGER MIN_INTM(MXPORB,N_ORB_CONF),
     &        MAX_INTM(MXPORB,N_ORB_CONF)
      INTEGER MIN_INTMS(MXPORB,N_ORB_CONF),
     &        MAX_INTMS(MXPORB,N_ORB_CONF)
*
      INTEGER NCONF_INTM(N_ORB_CONF)
      INTEGER NCSF_INTM(N_ORB_CONF)
      INTEGER NSD_INTM(N_ORB_CONF)
*. Local scratch
      INTEGER NOCPSM_IN(MXPOBS,2),NOCPSM_INTM(MXPOBS,2),
     &        NOCPSM_OUT(MXPOBS,2), NREM(MXPOBS)
*
* The occupations of the intermediate codes is based on the 
* following considerations:
*  In each step of the transformation one orbital is transformed
*  from initial to final basis. In step IORB, electrons in
*  orbitals 1 - IORB may this be added, but never removed
*  Note also that the transformation is symmetry conserving.
*  So restrictions does not only hold for complete
*  electron occupations, but also for occupations in each 
*  orbital symmetry
*
*. Note that MINMAX_INTM(*,*,IORB) refers to the occupations
*  after orbital IORB has been transformed
*
      NTEST = 1000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Info from MINMAX_FOR_ORBTRA '
        WRITE(6,*) ' ============================ '
        WRITE(6,*)
        WRITE(6,*) ' MINMAX for IN: '
        CALL WRT_MINMAX_OCC(MIN_IN,MAX_IN,N_ORB_CONF)
        WRITE(6,*) ' MINMAX for OUT: '
        CALL WRT_MINMAX_OCC(MIN_OUT,MAX_OUT,N_ORB_CONF)
      END IF
*
      IZERO = 0
*
*. Number of electrons per symmetry in IN  and OUT
*
C     MINMAX_PER_SYM(MIN_OCC,MAX_OCC,MIN_PER_SYM,MAX_PER_SYM)
      CALL MINMAX_PER_SYM(MIN_IN,MAX_IN,
     &     NOCPSM_IN(1,1),NOCPSM_IN(1,2))
      CALL MINMAX_PER_SYM(MIN_OUT,MAX_OUT,
     &     NOCPSM_OUT(1,1),NOCPSM_OUT(1,2))
*
*
* For convenience, during debugging
      INUM = -55
      DO IORB = 1, N_ORB_CONF-1
        CALL ISETVC(MIN_INTM(1,IORB),INUM,N_ORB_CONF)
        CALL ISETVC(MAX_INTM(1,IORB),INUM,N_ORB_CONF)
      END DO
*. Loop over orbitals to be transformed
      DO ITORB = 1, N_ORB_CONF
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Orbital to be transformed ', ITORB
        END IF
        IF(ITORB.EQ.N_ORB_CONF) THEN
*. Just copy final list
          N = N_ORB_CONF
          CALL ICOPVE(MIN_OUT(1),MIN_INTM(1,N),N)
          CALL ICOPVE(MAX_OUT(1),MAX_INTM(1,N),N)
        ELSE IF(ITORB.EQ.1) THEN
*. The number of electrons in orbital 1 cannot be increased
          MAX_INTM(1,ITORB) = MAX_IN(1)
          MIN_INTM(1,ITORB) = 0
*. The accumulated occupation for the remaining orbitals may be decreased by
*. the number of electrons in orbital 1
          MAX_AC = MAX_IN(1)
          MIN_AC = MIN_IN(1)
          DO IORB = ITORB+1, N_ORB_CONF
            MIN_INTM(IORB,ITORB) = 
     &      MAX(0,MIN_IN(IORB)-MAX_AC)
          END DO
*. In the untransformed orbitals: Never less in IORB-N_ORB_CONF than in INI
           DO IORB = ITORB+1, N_ORB_CONF
             MAX_INTM(IORB,ITORB) = MAX_IN(IORB)
           END DO
        ELSE
*. Max in ITORB
           MAX_AC  = 
     &     MIN((MAX_INTM(ITORB,ITORB-1) - MIN_INTM(ITORB-1,ITORB-1)),2)
*. Orbital IORB .le. ITORB:  Never more in these orbitals than in the end
*. Occupations once created are never annihilated
           DO IORB = 1, ITORB-1 
             MAX_INTM(IORB,ITORB) = MAX_OUT(IORB)
             MIN_INTM(IORB,ITORB) = MAX(0,MIN_INTM(IORB,ITORB-1)-MAX_AC)
           END DO
*. Orbital ITORB: Accumulated in 1 - ITORB can never be more than 
*  in the initial space
           MAX_INTM(ITORB,ITORB) = MAX_IN(ITORB)
           MIN_INTM(ITORB,ITORB) = MAX(0,MIN_INTM(ITORB,ITORB-1)-MAX_AC)

*. Orbital IORB .gt. ITORB: Never less in orbitals IORB- N_ORB_CONF than in INI
           DO IORB = ITORB+1, N_ORB_CONF
             MAX_INTM(IORB,ITORB) = MAX_IN(IORB)
             MIN_INTM(IORB,ITORB) = MAX(0,MIN_INTM(IORB,ITORB-1)-MAX_AC)
           END DO
         END IF ! switch between orbitals
        END DO ! loop over orbitals to be transformed
*
*. Ensure that the MINMAX arrays are consistent with atmost
*. two electrons in each orb
*
      IZEROSPC = 0
      DO ITORB = 1, N_ORB_CONF
        CALL CHECK_MINMAX(MIN_INTM(1,ITORB),MAX_INTM(1,ITORB),
     &       N_ORB_CONF,IZEROSPC)
C     CHECK_MINMAX(MIN_OCC,MAX_OCC,NORB,IZEROSPC)
        IF(IZEROSPC.EQ.1) THEN
          WRITE(6,*) ' Vanishing space detected by CHECK_MINMAX'
          STOP       ' Vanishing space detected by CHECK_MINMAX'
        END IF
      END DO
*
*
* Test: Set evrything to Max space
*
      IFUSK = 0
      IF(IFUSK .EQ.1) THEN
       DO I = 1, 100
         WRITE(6,*) ' MINMAX spaces set to largest possible space'
       END DO
       NELECT = MAX_IN(N_ORB_CONF)
       DO ITORB = 1, N_ORB_CONF
        DO IORB = 1, N_ORB_CONF
         MIN_INTM(IORB,ITORB) = 0
         MAX_INTM(IORB,ITORB) = NELECT
        END DO
        CALL CHECK_MINMAX(MIN_INTM(1,ITORB),MAX_INTM(1,ITORB),
     &       N_ORB_CONF,IZEROSPC)
       END DO
      END IF ! FUSK
*
      IF(IDODIM.EQ.1) THEN
       DO IORB = 0, N_ORB_CONF
C             GET_DIM_MINMAX_SPACE(MIN_OCC,MAX_OCC,NORB,ISYM,
         CALL GET_DIM_MINMAX_SPACE(MIN(1,IORB),MAX(1,IORB),
     &   IREO_MNMX_OB_NO,N_ORB_CONF,ISYM,NCONFL,NCSFL,NSDL)
         NCONF_INTM(IORB) = NCONFL
         NCSF_INTM(IORB) = NCSFL
         NSD_INTM(IORB) = NSDL
       END DO
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ========================================'
       WRITE(6,*) ' MINMAX arrays for orbital transformation'
       WRITE(6,*) ' ========================================'
       WRITE(6,*) 
       DO IORB = 1, N_ORB_CONF
         WRITE(6,'(A,I4)') ' After transforming orbital ', IORB
         WRITE(6,*)        ' =================================='
         WRITE(6,*) 
         CALL WRT_MINMAX_OCC(
     &   MIN_INTM(1,IORB),MAX_INTM(1,IORB),N_ORB_CONF)
         IF(IDODIM.EQ.1) WRITE(6,'(A,3I9)')
     &   ' Number Confs, CSFs, SDs ',
     &   NCONF_INTM(IORB),NCSF_INTM(IORB),NSD_INTM(IORB)
       END DO
      END IF
*
      RETURN
      END
      SUBROUTINE CHECK_MINMAX(MIN_OCC,MAX_OCC,NORB,IZEROSPC)
* 
* Accumulated occupations for configuration space is 
* given in the form or a min max space. Ensure that the space
* is physically reasonable:
* 1: each orbital may contain atmost two electrons
*
* The spaces are corrected to produce the same space as input
* Therefore: Min_occ may be increased and max_occ may be decreased
* 
* A vanisning space is flagged by IZEROSPC = 1
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER MIN_OCC(NORB),MAX_OCC(NORB)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from CHECK_MINMAX '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' MINMAX to be examined '
        CALL WRT_MINMAX_OCC(MIN_OCC,MAX_OCC,NORB)
      END IF
*
      IZEROSPC = 0
*. Check that MAX is larger to or equal to MIN
      DO IORB = 1, NORB
        IF(MIN_OCC(IORB).GT. MAX_OCC(IORB)) IZEROSPC = 1
      END DO
*. Ensure that lower bounds are non-negative
      DO IORB = 1, NORB
        IF(MIN_OCC(IORB).LT.0) MIN_OCC(IORB) = 0
      END DO
*. Upper bound negative => vanishing space
      DO IORB = 1, NORB
        IF(MAX_OCC(IORB).LT.0) IZEROSPC = 1
      END DO
*. Upper bound .le. number of electrons
      NELEC = MAX_OCC(NORB)
      DO IORB = 1, NORB
        IF(MAX_OCC(IORB).GT.NELEC) MAX_OCC(IORB) = NELEC
      END DO
*. Ensure non-decreasing upper and lowe bounds
      DO IORB = NORB, 2, -1
        IF(MAX_OCC(IORB-1).GT.MAX_OCC(IORB)) 
     &  MAX_OCC(IORB-1) = MAX_OCC(IORB)
      END DO
      DO IORB = 2, NORB
       IF(MIN_OCC(IORB-1).GT.MIN_OCC(IORB))
     &    MIN_OCC(IORB) = MIN_OCC(IORB-1)
      END DO
*. Atmost two electrons may be added in each orbital
      DO IORB =1, NORB
        IF(MAX_OCC(IORB).GT.2*IORB) MAX_OCC(IORB) = 2*IORB
        IF(MIN_OCC(IORB).GT.2*IORB) IZEROSPC = 1
      END DO
*. Atleast two electrons may be added in each of the remaining orbitals
      DO IORB = NORB,1,-1
        MAXLEFT = (NORB-IORB)*2
        IF(MIN_OCC(IORB).LE.NELEC-MAXLEFT) MIN_OCC(IORB) = NELEC-MAXLEFT
      END DO
*
      IF(NTEST.GE.100.OR.IZEROSPC.EQ.1) THEN
        IF(IZEROSPC.EQ.1) THEN
          WRITE(6,*) ' CHECK_MINMAX was presented for a vanishing space'
          WRITE(6,*) ' Space, perhaps partly cleaned up'
        ELSE
          WRITE(6,*) ' MINMAX space after shaving by CHECK_MINMAX'
        END IF
        CALL WRT_MINMAX_OCC(MIN_OCC,MAX_OCC,NORB)
      END IF
*
      RETURN
      END
      FUNCTION IINPROD(IA,IB,NDIM)
*
* Inner product of two integer arrays IA, IB
*
* Jeppe Olsen, July 16, 2011
*
      INCLUDE 'implicit.inc'
      INTEGER IA(NDIM), IB(NDIM)
*
      IPROD = 0
      DO I = 1, NDIM
        IPROD = IPROD + IA(I)*IB(I)
      END DO
*
      IINPROD = IPROD
*
      RETURN
      END
      SUBROUTINE TRACI_CONF(C,S,LUC,LUHC)
*
*. Perform orbital transformation in the configuration approach
*. Initial version some 40 hours before take of to WATOC 2011
*
* Note: Routine uses C as scratch so this is modified during calc.
*
*. The MO-MO transformation matrix is stored in KCBIO
*. The spaces defining the in and out spaces are defined by 
*  
* ICPSC_CN, ISSPC_CN parameters in cands
*
*
*. Last modification; Jeppe Olsen; June 18, 2013; Allowing inactive orbitals
*.                    and several symmetries(sic)
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'vb.inc'
      INCLUDE 'lucinp.inc'
*. Input and scratch
      DIMENSION C(*)
*. Output
      DIMENSION S(*)
*. Local scratch
      DIMENSION FUSK(1000)
*
      IDUM = 0
      CALL LUCIAQENTER('TRACNF')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACNF')
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from TRACI_CONF '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' LUC, LUHC = ', LUC, LUHC
        WRITE(6,*) ' ICISTR = ', ICISTR
        WRITE(6,*) ' ICSM, ISSM = ',ICSM, ISSM
      END IF
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' Initial vector to be transformed'
        NCSF_C = NCSF_PER_SYM_GN(ICSM,ICSPC_CN)
        CALL WRTMAT(C,1,NCSF_C,1,NCSF_C)
      END IF
*
* 1:  Obtain the matrix T defining the steps of the orbital transformation 
*     using the approach of PAM
* T 
      CALL MEMMAN(KLT,NTOOB**2,'ADDL  ',2,'TMAT  ')
      CALL MEMMAN(KLTB,NTOOB**2,'ADDL  ',2,'TMATBL')
*. Scratch in PAMTMT
      LSCR = NTOOB**2 +NTOOB*(NTOOB+1)/2
      CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'KLSCR ')
*. Each symmetry separate
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISM-1)**2
        END IF
        IF(NTOOBS(ISM).GT.0)
     &  CALL PAMTMT(WORK(KCBIO-1+IOFF),WORK(KLT-1+IOFF),
     &       WORK(KLSCR),NTOOBS(ISM))
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The T-matrix for the orbital trans '
        CALL APRBLM2(WORK(KLT),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. LUCIA will use space for one-electron integrals for orbital transformation.
*. save a copy of original KINT1
      LEN_1F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLINT1_ORIG,LEN_1F,'ADDL  ',2,'INT1_O')
      CALL COPVEC(WORK(KINT1),WORK(KLINT1_ORIG),LEN_1F)
*. Default block size
      LCSBLK_L = LCSBLK
      IF(LCSBLK_L.LE.0) THEN
        WRITE(6,*) ' SIGMA_CONF will define length of batch '
        LCSBLK_DEFAULT = 2000000
*. Compare with dimension of largest single configuration
        LCONF_MAX = IMNMX(NPCSCNF,MAXOP+1,2)
        IF(LCONF_MAX.GT.LCSBLK_DEFAULT) LCSBLK_DEFAULT = LCONF_MAX
        LCSBLK_L = LCSBLK_DEFAULT
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' LCSBLK_L = ', LCSBLK_L
      END IF
*
      ICSPC_CN_SAVE = ICSPC_CN
      ISSPC_CN_SAVE = ISSPC_CN
*
*. Now do the transformation for each orbital
*
      
      DO IORB = 1, N_ORB_CONF
*. We are looping over orbitals in the configurations, i.e. 
*. in type-order 
         IIORB = IB_ORB_CONF -1 + IORB
         IIORB_SO = IREOTS(IIORB)
         IF(NTEST.GE.1000) THEN
           WRITE(6,'(A,I2,I3) ') 
     &     ' >>>> Info for orb. transformation for orbital',
     &     IORB
         END IF
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' IORB, IIORB,IIORB_SO = ',
     &                  IORB, IIORB,IIORB_SO
         END IF
* For each orbital I we will calculate 
*( 1+ \hat T(I) + 1/2\hat T(I)^2)) TII^\hat N_I C(I-1), 
* where C(I-1) is result of all previous transformations.
*. We will collect the contributions for each orb in KLCSFVC  
*. At start we have the transformed operator so far in C
*
* Prepare for transforming orbital IORB
*
*. Place (T(P,I)/S(I,I)   in one-electron integral list
C            T_ROW_TO_H(T,H,K)
        CALL T_ROW_TO_H(WORK(KLT),WORK(KINT1),IIORB_SO,TII)
*. T_{II}^Ni C in ICSPC_CN, save in C
C           T_TO_NK_T_VEC_CONF(T,K,VEC,ISPC,ISYM)
        CALL T_TO_NK_T_VEC_CONF(TII,IORB,C,ICSPC_CN,ICSM)
*
        ICSPC_CN = IORBTRA_SPC_IN(IORB)
        ISSPC_CN = IORBTRA_SPC_OUT(IORB)
        NCSF_C = NCSF_PER_SYM_GN(ICSM,ICSPC_CN)
        NCSF_S = NCSF_PER_SYM_GN(ISSM,ISSPC_CN)
        NCSF_CS = MAX(NCSF_S,NCSF_C)
*. A scratch CSF vector
        CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACNI')
        CALL MEMMAN(KLCSFVC,NCSF_MNMX_MAX,'ADDL  ',2,'CSFVC ')
*. Loop over the two operators needed for each orbitaltransf
        DO IPOT = 1, 2
         CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACNI')
* The input (C) and output spaces
         IF(IPOT.EQ.1) THEN
           ICSPC_CN = IORBTRA_SPC_IN(IORB)
           ISSPC_CN = IORBTRA_SPC_OUT(IORB)
         ELSE
           ICSPC_CN = IORBTRA_SPC_OUT(IORB)
           ISSPC_CN = IORBTRA_SPC_OUT(IORB)
         END IF
         IF(NTEST.GE.100)
     &   WRITE(6,*) ' ICSPC_CN, ISSPC_CN = ', ICSPC_CN, ISSPC_CN
         IF(IPOT.EQ.1) THEN
*. Expand TII^\hat N_I C(I-1) in CSFVC
C          REF_CNFVEC(VECIN,ISPCIN,VECOUT,ISPCOUT,ISYM)
           CALL REF_CNFVEC(C,ICSPC_CN,WORK(KLCSFVC),ISSPC_CN,ICSM)
         END IF
*
         NCONF_C = NCONF_PER_SYM_GN(ICSM,ICSPC_CN)
         NSD_C   = NSD_PER_SYM_GN(ICSM,ICSPC_CN)
         NCSF_C = NCSF_PER_SYM_GN(ICSM,ICSPC_CN)
*
         NCONF_S = NCONF_PER_SYM_GN(ISSM,ISSPC_CN)
         NSD_S   = NSD_PER_SYM_GN(ISSM,ISSPC_CN)
         NCSF_S = NCSF_PER_SYM_GN(ISSM,ISSPC_CN)
*
         NCONF_MAX = MAX(NCONF_C,NCONF_S)
*
         IF(NTEST.GE.1000) THEN
           WRITE(6,'(A,3I8)') 
     &     ' Number of confs, SDs and CSFs for C ',
     &       NCONF_C, NSD_C, NCSF_C
           WRITE(6,'(A,3I8)') 
     &     ' Number of confs, SDs and CSFs for S ',
     &      NCONF_S, NSD_S, NCSF_S
         END IF
*
*
*. If ICISTR = 1, vectors are stored in one batch, so
         IF(ICISTR.EQ.1) LCSBLK_L = MAX(NCSF_C,NCSF_S)
         IF(NTEST.GE.100) WRITE(6,*) ' Size of batch ', LCSBLK_L
*. Batches of C
*. ==============
*. One could here either use CSF's or SD's. As memory maybe the defining parameter,
* I opt for CSF's and will then expand/contract each configuration when needed.
*. Length of each configuration
         CALL MEMMAN(KLLCNFEXP,NCONF_MAX,'ADDL  ',1,'LCNFEX')
*. For C
C        CONF_EXP_LEN_LIST(ILEN,NCONF_PER_OPEN,NELMNT_PER_OPEN,MAXOP)
         CALL CONF_EXP_LEN_LIST(WORK(KLLCNFEXP),
     &        NCONF_PER_OPEN_GN(1,ICSM,ICSPC_CN),NPCSCNF,MAXOP)
C        PART_VEC(LBLK,NBLK,MAXSTR,LBAT,NBAT,IONLY_NBAT)
         CALL PART_VEC(WORK(KLLCNFEXP),NCONF_C,LCSBLK_L,IDUM,NBAT_C,1)
         CALL MEMMAN(KLLBAT_C,NBAT_C,'ADDL  ',1,'LBAT_C')
         CALL PART_VEC(WORK(KLLCNFEXP),NCONF_C,LCSBLK_L,WORK(KLLBAT_C),
     &        NBAT_C,0)
*. And for Sigma
C        CONF_EXP_LEN_LIST(ILEN,NCONF_PER_OPEN,NELMNT_PER_OPEN,MAXOP)
         CALL CONF_EXP_LEN_LIST(WORK(KLLCNFEXP),
     &        NCONF_PER_OPEN_GN(1,ISSM,ISSPC_CN),NPCSCNF,MAXOP)
C        PART_VEC(LBLK,NBLK,MAXSTR,LBAT,NBAT,IONLY_NBAT)
         CALL PART_VEC(WORK(KLLCNFEXP),NCONF_S,LCSBLK_L,IDUM,NBAT_S,1)
         CALL MEMMAN(KLLBAT_S,NBAT_S,'ADDL  ',1,'LBAT_S')
         CALL PART_VEC(WORK(KLLCNFEXP),NCONF_S,LCSBLK_L,WORK(KLLBAT_S),
     &      NBAT_S,0)
*
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Number of batches for C and S ', NBAT_C, NBAT_S
         END IF
*. Largest number of configurations in a given batch
         MAX_CONF_BATCH_C = IMNMX(WORK(KLLBAT_C),NBAT_C,2)
         MAX_CONF_BATCH_S = IMNMX(WORK(KLLBAT_S),NBAT_S,2)
         MAX_CONF_BATCH = MAX(MAX_CONF_BATCH_C,MAX_CONF_BATCH_S)
*
         IF(NTEST.GE.1000)
     &   WRITE(6,*) ' Largest number of configs in batch ',
     &   MAX_CONF_BATCH
         CALL MEMMAN(KLLBLK_BAT_C,MAX_CONF_BATCH ,'ADDL  ',2,'LBLBTC')
         CALL MEMMAN(KLLBLK_BAT_S,MAX_CONF_BATCH ,'ADDL  ',2,'LBLBTS')
*. Two vectors for holding expansion in SD of given config
         LEN_SD_CONF_MAX = IMNMX(NPDTCNF,MAXOP+1,2)
         CALL MEMMAN(KLCONF_SD_C,LEN_SD_CONF_MAX,'ADDL  ',2,'CN_SDC')
         CALL MEMMAN(KLCONF_SD_S,LEN_SD_CONF_MAX,'ADDL  ',2,'CN_SDS')
*. Scratch space in routine for evuluating H for configurations (allowing combs)
*. Scratch: Length: INTEGER: (NDET_C + NDET_S)*N_EL_CONF + NDET_C + 6*NORB
         L_CNHCN = LEN_SD_CONF_MAX*(1+2*N_EL_CONF) + 6*N_ORB_CONF
         CALL MEMMAN(KL_CNHCN, L_CNHCN,'ADDL  ',1,'LCNHCN')
*. Space for two integers arrays for signs
         CALL MEMMAN(KLISIGNC,LEN_SD_CONF_MAX,'ADDL  ',1,'ISIGNC')
         CALL MEMMAN(KLISIGNS,LEN_SD_CONF_MAX,'ADDL  ',1,'ISIGNS')
*
       ZERO = 0.0D0
       CALL SETVEC(S,ZERO,NCSF_S)
*
        IADOB = IB_ORB_CONF - 1
        CALL MEMCHK2('BETRAC')
        IF(IPOT.EQ.1) THEN
          XXNORM = INPROD(C,C,NCSF_C)
          WRITE(6,*) ' Norm**2 C(ini) = ', XXNORM
        END IF
        CALL TRACI_CONF_SLAVE(C,S,LUC,LUHC,ICISTR,
     &       NCONF_PER_OPEN_GN(1,ICSM,ICSPC_CN),
     &       NCONF_PER_OPEN_GN(1,ISSM,ISSPC_CN),
     &       NBAT_C,WORK(KLLBAT_C),
     &       NBAT_S,WORK(KLLBAT_S),
     &       WORK(KLLBLK_BAT_C),WORK(KLLBLK_BAT_S),
     &       WORK(KLCONF_SD_C),WORK(KLCONF_SD_S),
     &       IADOB,WORK(KDFTP),WORK(KL_CNHCN),
     &       WORK(KLISIGNC),WORK(KLISIGNS),IORB)
        CALL MEMCHK2('AFTRAC')
*
*. And copy output to input for next round..
*
        ONE = 1.0D0
        IF(IPOT.EQ.1) THEN
*. Collecting (1 + T ) !C(K-1)> in KLCSFVC
          FACTOR = 1.0D0
          CALL VECSUM(WORK(KLCSFVC),WORK(KLCSFVC),S,ONE,FACTOR,NCSF_S)
          XXNORM = INPROD(WORK(KLCSFVC),WORK(KLCSFVC),NCSF_S)
          WRITE(6,*) ' Norm**2 (1+T)!Prev> = ', XXNORM
CD        IF(NTEST.GE.1000) THEN
CD         WRITE(6,*) ' Fusk Updated Sigma vector (1+T)!Prev> '
CD         WRITE(6,*) ' Fusk Updated Sigma vector (1+T)!Prev> '
CD         CALL CSDTVC_CONFSPACE(NCONF_S,WORK(KLCSFVC),
CD   &          FUSK,ISSM,ISSPC_CN,1)
CD        END IF
*. And prepare for next op
          CALL COPVEC(S,C,NCSF_S)
        ELSE
*. Collecting (1 + T + 1/2T^2) !C(K-1)> in KLCSFVC
          FACTOR = 0.5D0
          ONE = 1.0D0
          CALL VECSUM(WORK(KLCSFVC),WORK(KLCSFVC),S,ONE,FACTOR,NCSF_S)
          CALL COPVEC(WORK(KLCSFVC),C,NCSF_S)
          CALL COPVEC(WORK(KLCSFVC),S,NCSF_S)
          XXNORM = INPROD(WORK(KLCSFVC),WORK(KLCSFVC),NCSF_S)
          WRITE(6,*) ' Norm**2 (1+T+1/2 T^2)!Prev> = ', XXNORM
*
CD        WRITE(6,*) ' Fusk, (1 + T + 1/2T^2) !C(K-1)> SD basis '
CD        WRITE(6,*) ' Fusk, (1 + T + 1/2T^2) !C(K-1)> SD basis '
CD        CALL CSDTVC_CONFSPACE(NCONF_S,S,FUSK,ISSM,ISSPC_CN,1)
        END IF
        CALL MEMCHK2('AFTSUM')
*
*
        CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACNI')
       END DO !End of loop over the two powers of the operator
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Updated TRACI vector after a orbtrans'
         CALL WRTMAT(S,1,NCSF_S,1,NCSF_S)
       END IF
       CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACNI')
      END DO ! Loop over orbitals to be transformed
*. and restors defs
      ICSPC_CN = ICSPC_CN_SAVE 
      ISSPC_CN = ISSPC_CN_SAVE 
      CALL COPVEC(WORK(KLINT1_ORIG),WORK(KINT1),LEN_1F)
*
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' Final PAM transformed CI vector '
        WRITE(6,*) ' ================================'
        NCSF_S = NCSF_PER_SYM_GN(ISSM,ISSPC_CN)
        CALL WRTMAT(S,1,NCSF_S,1,NCSF_S)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACNF')
      CALL LUCIAQEXIT('TRACNF')
      RETURN
      END 
      SUBROUTINE TRACI_CONF_SLAVE(C,S,LUC,LUS,ICISTR,
     &           NCONF_PER_OPEN_C,NCONF_PER_OPEN_S,
     &           NBAT_C,LBAT_C,NBAT_S,LBAT_S,
     &           LBLK_BAT_C,LBLK_BAT_S,
     &           CONF_SD_C,CONF_SD_S,IADOB,IPRODT,
     &           ISCR_CNHCN,ISIGN_C,ISIGN_S,IORB)
*
* Inner (aka slave) routine for orbital transformation in configuration based methods
*
* Transform Orbital IORB
*
*. Jeppe Olsen,July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cecore.inc'
*. Input
*. C-vector or space for batch of C-vector
      DIMENSION C(*)
*. Info on the two configuration expansions
       INTEGER NCONF_PER_OPEN_C(*), NCONF_PER_OPEN_S(*)
*. Number of blocks in the batches of C and S
      INTEGER LBAT_C(*), LBAT_S(*)
*. Scratch for Info on batches of C and S: Length of each block (configuration in batch)
      INTEGER LBLK_BAT_C(*),LBLK_BAT_S(*)
*. Space for SD expansion of single configurations
      DIMENSION CONF_SD_C(*), CONF_SD_S(*)
*. Space for signs for phase change for dets of a configurations
      INTEGER ISIGN_C(*),ISIGN_S(*)
*. CSF info: proto type dets 
      INTEGER IPRODT(*)

*. Output
      DIMENSION S(*)
*. Scratch transferred through to CNHCN
      INTEGER ISCR_CNHCN
*. Local scratch
      INTEGER IOCC_C(MXPORB),IOCC_S(MXPORB)
*
*. TEMP SCRATCH
      DIMENSION SFUSK(2000), SFUSK2(2000)
*
      NTEST = 0010
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from TRACI_CONF_SLAVE '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' ICISTR = ', ICISTR
        WRITE(6,'(A,I4)') ' IORB = ', IORB
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACNI')
*. Initialization of some parameters for controlling loop over configurations
      IOPEN_S = 0
      INUM_OPS = 0
      IOPEN_C = 0
      INUM_OPC = 0
      IB_CSF_C = 1
      IB_CSF_S = 1
*
      CALL MEMCHK2('INISIG')
*
*. Loop over batches of S
      INI_S = 1
C?    WRITE(6,*) ' NBAT_S = ', NBAT_S
      DO IBAT_S = 1, NBAT_S
       IF(NTEST.GE.1000) 
     & WRITE(6,'(A,I3)') ' >>> Start of sigma batch ', IBAT_S
*
       IF(IBAT_S.EQ.1) THEN
         IB_CONF_S = 1
         IB_CSF_S = 1
       ELSE
         IB_CONF_S = IB_CONF_S + LBAT_S(IBAT_S-1)
       END IF
C?     WRITE(6,*) ' LBAT_S(1) = ', LBAT_S(1)
       N_CONF_S = LBAT_S(IBAT_S) 
*. Number of CSF's per config in S-batch
C           GET_LBLK_CONF_BATCH(ICNF_INI,NCNF,LBLK_BAT,ISYM,ISPC,
C    &      NSD_BAT_TOT,NCSF_BAT_TOT) 
       CALL GET_LBLK_CONF_BATCH(IB_CONF_S,N_CONF_S,LBLK_BAT_S,ISSM,
     &      ISSPC_CN,NSD_BAT_TOT_S,NCSF_BAT_TOT_S)
       CALL MEMCHK2('AFGTL1')
       IF(NTEST.GE.100) THEN
         WRITE(6,'(A,2I9)') 
     &   ' Number of CSFs and SDs in S-batch ', NCSF_BAT_TOT_S,
     &     NSD_BAT_TOT_S
       END IF
*. Initialize sigma batch
       ZERO = 0.0D0
C?     WRITE(6,*) ' IB_CSF_S, NCSF_BAT_TOT_S = ',
C?   &              IB_CSF_S, NCSF_BAT_TOT_S
       CALL SETVEC(S(IB_CSF_S),ZERO,NCSF_BAT_TOT_S)
*. Loop over batches of C
       IF(ICISTR.NE.1) REWIND LUS
       INI_C = 1
*. First time in this batch
       ISBAT_FIRST_TIME =1
       DO IBAT_C = 1, NBAT_C
        IF(NTEST.GE.1000) 
     &  WRITE(6,'(A,I3)') ' >>> Start of C batch ', IBAT_C
        CALL MEMCHK2('STCBAT')
        IF(IBAT_C.EQ.1) THEN
          IB_CONF_C = 1
          IB_CSF_C = 1
        ELSE
          IB_CONF_C = IB_CONF_C + LBAT_C(IBAT_C-1)
        END IF
        N_CONF_C = LBAT_C(IBAT_C)
*. Number of configs per config in S-batch
        CALL GET_LBLK_CONF_BATCH(IB_CONF_C,N_CONF_C,LBLK_BAT_C,ICSM,
     &      ICSPC_CN,NSD_BAT_TOT_C,NCSF_BAT_TOT_C)
        IF(NTEST.GE.100) THEN
          WRITE(6,'(A,2I9)') 
     &   ' Number of CSFs and SDs in C-batch ', NCSF_BAT_TOT_C,
     &     NSD_BAT_TOT_C
        END IF
      CALL MEMCHK2('AFGTLB')
*. Read, if required, next batch of C- Each configuration stored in a record by itself
        IF(ICISTR.NE.1) THEN
          CALL FRMDSCN(C,N_CONF_C,-1,LUC)
C              FRMDSCN(VEC,NREC,LBLK,LU)
        END IF
*. And then to the configurations of the C and sigma
*. First time in this batch
        IF(ISBAT_FIRST_TIME.EQ.1) THEN
* Save pointers to start of configuration
          IOPEN_S_SAVE = IOPEN_S
          INUM_OPS_SAVE = INUM_OPS
          IB_CSF_S_SAVE = IB_CSF_S
        ELSE
          IOPEN_S = IOPEN_S_SAVE
          INUM_OPS = INUM_OPS
          IB_CSF_S = IB_CSF_S_SAVE
        END IF
        ISBAT_FIRST_TIME  = 0
*
        ICBAT_FIRST_TIME =1
        DO ICONF_S = IB_CONF_S, IB_CONF_S + N_CONF_S -1
*. Obtain occupation in IOCC_S and iopen for this sigma-configuration
C             NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
         CALL NEXT_CONF_IN_CONFSPC(IOCC_S,IOPEN_S,INUM_OPS,INI_S,
     &        ISSM,ISSPC_CN,NEW_S)
         INI_S = 0
         IOCOB_S = (IOPEN_S + N_EL_CONF)/2
*. Signs for going between configuration and interaction order of dets
C     SIGN_CONF_SD(ICONF,NOB_CONF,IOP,ISGN,IPDET_LIST,ISCR)
         CALL SIGN_CONF_SD(IOCC_S,IOCOB_S,IOPEN_S,ISIGN_S,IPRODT,
     &                      ISCR_CNHCN)
      CALL MEMCHK2('AFSIGN')
*
         NCSF_S = NPCSCNF(IOPEN_S+1)
         NSD_S = NPDTCNF(IOPEN_S+1)
C?       IF(NSD_S.EQ.6) THEN
C?         WRITE(6,*) ' Fusk NTEST increased '
C?         WRITE(6,*) ' Fusk NTEST increased '
C?         WRITE(6,*) ' Fusk NTEST increased '
C?         NTEST = 10000
C?       END IF
*
         IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Sigma configuration number ', ICONF_S
          CALL IWRTMA(IOCC_S,1,IOCOB_S,1,IOCOB_S)
         END IF
*
         ZERO = 0.0D0
*. The contribution to a given sigma conf from all C-conf in C-batch
* will be stored in CONF_SD_S(1)
         CALL SETVEC(CONF_SD_S,ZERO,NSD_S)
*
         IF( ICBAT_FIRST_TIME .EQ. 1) THEN
           IOPEN_C_SAVE = IOPEN_C
           INUM_OPC_SAVE = INUM_OPC
           IB_CSF_C_SAVE = IB_CSF_C
         ELSE
           IOPEN_C = IOPEN_C_SAVE
           INUM_OPC = INUM_OPC_SAVE
           IB_CSF_C = IB_CSF_C_SAVE
         END IF
         ICBAT_FIRST_TIME = 0
         DO ICONF_C = IB_CONF_C, IB_CONF_C + N_CONF_C -1
           IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' ICONF_C, ICONF_S: ', ICONF_C, ICONF_S
           END IF
*. Obtain occupation in IOCC_C and iopen for this C-configuration
C             NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
          CALL NEXT_CONF_IN_CONFSPC(IOCC_C,IOPEN_C,INUM_OPC,INI_C,
     &         ICSM,ICSPC_CN,NEW_C)
*. Signs for going between configuration and interaction order of dets
C     SIGN_CONF_SD(ICONF,NOB_CONF,IOP,ISGN,IPDET_LIST,ISCR)
          IOCOB_C = (IOPEN_C + N_EL_CONF)/2
          CALL SIGN_CONF_SD(IOCC_C,IOCOB_C,IOPEN_C,ISIGN_C,IPRODT,
     &                      ISCR_CNHCN)
          INI_C = 0
          IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' C configuration number ', ICONF_C
           IOCOB_C = (IOPEN_C + N_EL_CONF)/2
           CALL IWRTMA(IOCC_C,1,IOCOB_C,1,IOCOB_C)
          END IF
*. Expand coefficients for configuration from CSF to SD basis
          NCSF_C = NPCSCNF(IOPEN_C+1)
          NSD_C = NPDTCNF(IOPEN_C+1)
C              CSDTVC_CONF(C_SD,C_CSF,NOPEN,ISIGN,IAC,IWAY)
          CALL CSDTVC_CONF(CONF_SD_C,C(IB_CSF_C),IOPEN_C,ISIGN_C,2,1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' C(ICONF_C)  in SD'
            CALL WRTMAT(CONF_SD_C,1,NSD_C,1,NSD_C)
          END IF
          IF(NTEST.GE.1000)  THEN
            WRITE(6,'(A,2I6)') 
     &      ' Info on sigma for ICONF_C, ICONF_S = ',
     &      ICONF_C, ICONF_S
          END IF
*. Update: S(I) = S(I) + Sum(J) sum p <I!a+_(P,IAB)a_(IORB,IAB)!J> C(J)
          I12OP = 1
          I_DO_ORBTRA = 1
*. As want to add S(I), we set a local core-energy to one
          ECORE_L = 0.0D0
C     CNHCN_LUCIA(ICNL,IOPL,ICNR,IOPR,CNHCNM,SIGMA,
C    &           IADOB,IPRODT,I12OP,IORBTRA,ECORE,ISCR)
          CALL CNHCN_LUCIA(IOCC_S,IOPEN_S,IOCC_C,IOPEN_C,
     &         CONF_SD_C,XDUM,CONF_SD_S,IADOB,
     &         IPRODT,I12OP, I_DO_ORBTRA, IORB, ECORE_L,2,
     &         0,RJ,RK,ISCR_CNHCN)
*. Update address of C in action
          IB_CSF_C = IB_CSF_C + NCSF_C
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Updated Sigma(ICONF_S)  in SD'
            CALL WRTMAT(CONF_SD_S,1,NSD_S,1,NSD_S)
          END IF
         END DO ! over configs in batch of C
*. And transform sigma part to CSF and update sigma vector
         CALL CSDTVC_CONF(CONF_SD_S,S(IB_CSF_S),IOPEN_S,ISIGN_S,1,2)
         IB_CSF_S = IB_CSF_S + NCSF_S
        END DO ! over configs in batch of S
       END DO ! Over batches of C
      END DO ! over batches of Sigma
*
*. Test transformation back to CSF 
*
C?    WRITE(6,*) ' FUSK: back transf to SD basis at end of TRACI..'
*
C?    INI_S = 1
C?    IB_CSF = 1
C?    IB_SD = 1
C?    N_CONF_S = LBAT_S(1)
C?    DO ICONF_S = 1, N_CONF_S
C?      CALL NEXT_CONF_IN_CONFSPC(IOCC_S,IOPEN_S,INUM_OPS,INI_S,
C?   &       ISSM,ISSPC_CN,NEW_S)
C?      INI_S = 0
C?      IOCOB_S = (IOPEN_S + N_EL_CONF)/2
*. Signs for going between configuration and interaction order of dets
C?      CALL SIGN_CONF_SD(IOCC_S,IOCOB_S,IOPEN_S,ISIGN_S,IPRODT,
C?   &                     ISCR_CNHCN)
C?      NCSF_S = NPCSCNF(IOPEN_S+1)
C?      NSD_S = NPDTCNF(IOPEN_S+1)
C?      CALL CSDTVC_CONF(SFUSK(IB_SD),S(IB_CSF),IOPEN_S,ISIGN_S,2,1)
C?      IB_CSF = IB_CSF + NCSF_S
C?      IB_SD = IB_SD + NSD_S
C?    END DO
C?    WRITE(6,*) ' Resulting vector transformed to SD''s '
C?    NSD_TOT =  IB_SD - 1
C?    CALL WRTMAT(SFUSK,1,NSD_TOT,1,NSD_TOT)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACNI')
      IF(NTEST.GE.10)  WRITE(6,*) ' Returning from TRACI_CONF_SLAVE '
      RETURN
      END
      SUBROUTINE REF_CNFVEC(VECIN,ISPCIN,VECOUT,ISPCOUT,ISYM)
* 
* A vector VECIN is given in configuration space IVECIN.
* Obtain corresponding vector in configuratin space ISPCOUT
* Terms that are in ISPCIN, but not in ISPCOUT are eliminated
* Terms that are in ISPCOUT, but not in ISPCIN are set to 0
*
*
*. Jeppe Olsen, July 17, 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'vb.inc'
*. Input
      DIMENSION VECIN(*)
*. Output
      DIMENSION VECOUT(*)
*. Local scratch
      DIMENSION IOCC(MXPORB), IOCC2(MXPORB), IOCC3(MXPORB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from REF_CNFVEC '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,'(A,2I5)') ' In- and Out-spaces: ', ISPCIN,ISPCOUT
        WRITE(6,*) ' ISYM = ', ISYM
      END IF
*
      NCSF_IN = NCSF_PER_SYM_GN(ISYM,ISPCIN)
      NCSF_OUT = NCSF_PER_SYM_GN(ISYM,ISPCOUT)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCSF_IN, NCSF_OUT = ', NCSF_IN, NCSF_OUT
      END IF
      ZERO = 0.0D0
      CALL SETVEC(VECOUT,ZERO,NCSF_OUT)
*
      INI = 1
      IB_IN = 1
      NEW = 1
      DO IOPEN = 0, MAXOP
        NCSF_PT = NPCSCNF(IOPEN+1)
        NCNF_OPEN_IN = NCONF_PER_OPEN_GN(IOPEN+1,ISYM,ISPCIN)
        NOCOBL = (IOPEN+N_EL_CONF)/2
*. First configuration in out space with given number of open orbs and sym
        IF_OPEN_OUT = IB_CONF_REO_GN(IOPEN+1,ISYM,ISPCOUT)
*. Offset in CSF vector to first elements with given sym and number of orbs
        IB_OPEN_OUT = IB_OPEN_CSF(IOPEN+1,ISYM,ISPCOUT)
        DO ICNF = 1, NCNF_OPEN_IN
*. Obtain occupation of configuration
C              NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
          CALL NEXT_CONF_IN_CONFSPC(IOCC,IOPENX,INUM_OP_IN,INI,ISYM,
     &         ISPCIN,NEW)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*)  ' Next config from NEXT_CONF.... '
            NOCOBL = (N_EL_CONF + IOPEN)/2
            CALL IWRTMA(IOCC,1,NOCOBL,1,NOCOBL)
          END IF
          INI = 0
*.Is IOCC in output space?
* Reform from compact to occ number form
C  REFORM_CONF_OCC2(ICONF_EXP,ICONF_PACK,NORBL,NOCOBL,IWAY)
          CALL REFORM_CONF_OCC2(IOCC2,IOCC,N_ORB_CONF,NOCOBL,2)
*.occ number to accumulated
C         REFORM_CONF_ACCOCC(IACOCC,IOCC,IWAY,NORB)
          CALL REFORM_CONF_ACCOCC(IOCC3,IOCC2,2,N_ORB_CONF)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Next configuration in accumulated form '
            CALL IWRTMA(IOCC3,1,N_ORB_CONF,1,N_ORB_CONF)
          END IF
*. Check to see if configuration is within bounds
          IN_OUT = IS_IACC_CONF_IN_MINMAX_SPC(IOCC3,
     &             IOCC_MIN_GN(1,ISPCOUT),IOCC_MAX_GN(1,ISPCOUT),
     &             N_ORB_CONF)
          IF(IN_OUT.EQ.1) THEN
*. Find number of this configuration
C                  ILEX_FOR_CONF_G(ICONF,NOCC_ORB,ICONF_SPC,IDOREO)
            ILEX = ILEX_FOR_CONF_G(IOCC,NOCOBL,ISPCOUT,1)
            IB_OUT = IB_OPEN_OUT  + (ILEX-IF_OPEN_OUT)*NCSF_PT
            CALL COPVEC(VECIN(IB_IN),VECOUT(IB_OUT),NCSF_PT)
          END IF ! conf was in out space
            IB_IN = IB_IN + NCSF_PT
        END DO ! End of loop over input configs with a given number of open orbs
      END DO ! End of loop over number of open orbitals
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input and output vectors '
        CALL WRTMAT(VECIN,1,NCSF_IN,1,NCSF_IN)
        WRITE(6,*)
        CALL WRTMAT(VECOUT,1,NCSF_OUT,1,NCSF_OUT)
      END IF
*
      RETURN
      END
      FUNCTION IS_IACC_CONF_IN_MINMAX_SPC(IOCC,MIN_OCC,MAX_OCC,NORB)
*
* An accumulated configuration IOCC is given. Check if this configuration
* in in space defined by MIN_OCC, MAX_OCC.
* Returns 1/0 as answer
*
*. Jeppe Olsen, July 2011
*
      INTEGER MIN_OCC(NORB),MAX_OCC(NORB)
      INTEGER IOCC(NORB)
*
      INBOUND = 1
      DO IORB = 1, NORB
        IF(MIN_OCC(IORB).GT.IOCC(IORB).OR.
     &     IOCC(IORB).GT.MAX_OCC(IORB)) INBOUND = 0
      END DO
*
      IS_IACC_CONF_IN_MINMAX_SPC = INBOUND
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Configuration: '
        CALL IWRTMA(IOCC,1,NORB,1,NORB)
        IF(INBOUND.EQ.1) THEN
          WRITE(6,*) ' Configuration is in space '
        ELSE
          WRITE(6,*) ' Configuration is not in space '
        END IF
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Min Max space tested: '
        CALL WRT_MINMAX_OCC(MIN_OCC,MAX_OCC,NORB)
      END IF
*
      RETURN
      END
      SUBROUTINE T_TO_NK_T_VEC_CONF(T,K,VEC,ISPC,ISYM)
*
* A vector VEC is given in CI space ISPC.
* Multiply with T^(\hat N_k), where \hat N_k is the 
* number operator for orbital K
*
*. Jeppe Olsen, July 17, 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'vb.inc'
*. Input and output
      DIMENSION VEC(*)
*. Local scratch
      INTEGER IOCC(MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from T_TO_NK_T_VEC_CONF '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        WRITE(6,'(A,I5)') ' Confspaces: ', ISPC
        WRITE(6,'(A,I3,2X,E13.7)') ' K and T ', K, T
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input  vector to T_TO_NK_T_VEC_CONF '
        NCSF = NCSF_PER_SYM_GN(ISYM,ISPC)
        CALL WRTMAT(VEC,1,NCSF,1,NCSF)
      END IF
*
      TT = T*T
*
      INI = 1
      IB = 1
      NEW = 1
      DO IOPEN = 0, MAXOP
        NCSF_PT = NPCSCNF(IOPEN+1)
        NCNF_FOR_IOPEN = NCONF_PER_OPEN_GN(IOPEN+1,ISYM,ISPC)
        NOCOBL = (IOPEN+N_EL_CONF)/2
        DO ICNF = 1, NCNF_FOR_IOPEN
*. Obtain occupation of next configuration
C              NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
          CALL NEXT_CONF_IN_CONFSPC(IOCC,IOPENX,INUM_OP,INI,ISYM,
     &       ISPC,NEW)
          INI = 0
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*)  ' Next config from NEXT_CONF.... '
            CALL IWRTMA(IOCC,1,NOCOBL,1,NOCOBL)
          END IF
*
*. Number of electrons in K
          NKOCC = 0
          DO IORB = 1, NOCOBL
            IF(IOCC(IORB).EQ.K) THEN
*. Singly occupied
              NKOCC = 1
            ELSE IF(IOCC(IORB).EQ.-K) THEN
*. Doubly occupied
              NKOCC = 2
            END IF
          END DO
*
          IF(NKOCC.EQ.1) THEN
            CALL SCALVE(VEC(IB),T,NCSF_PT)
          ELSE IF (NKOCC.EQ.2) THEN 
            CALL SCALVE(VEC(IB),TT,NCSF_PT)
          END IF
*
          IB = IB + NCSF_PT
        END DO ! End of loop over input configs with a given number of open orbs
      END DO ! End of loop over number of open orbitals
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output vector from T_TO_NK_T_VEC_CONF '
        NCSF = NCSF_PER_SYM_GN(ISYM,ISPC)
        CALL WRTMAT(VEC,1,NCSF,1,NCSF)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_EXPMKS(EXPMKS,KAPPA_S, KAPPA_A,S,NOBPS,NSMOB)
*
* A symmetric and an antisymmetric kappa-matrix, KAPPA_S, KAPPA_A, 
* respectively, are given for a orbital space, in complete form  
* Obtain Exp (-Kappa_A S) Exp(-Kappa_S S)
*
* By varying the choice of NOBPS, the code can be used both for 
* a complete and for a subspace matrix.
*
* Jeppe Olsen, July 19 in Santiago de COmpostela, 24 hours before talk
* (I decided on the plane to make a MCSCF program for the VB code ...)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
*. Input
      REAL*8 KAPPA_S(*), KAPPA_A(*), S(*)
      INTEGER NOBPS(NSMOB)
*. Output
      DIMENSION EXPMKS(*)
*.
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2, 'GTEMKS')
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output from GET_EXPMKS '
        WRITE(6,*) ' ======================='
        WRITE(6,*) 
        WRITE(6,*) ' Input matrix KAPPA_A '
        CALL APRBLM2(KAPPA_A,NOBPS,NOBPS,NSMOB,0)
        WRITE(6,*) ' Input matrix KAPPA_S '
        CALL APRBLM2(KAPPA_S,NOBPS,NOBPS,NSMOB,0)
      END IF
*
* Exp (-Kappa_x S ) =  S^(-1/2) Exp(-S^(1/2) Kappa_x S^1/2) S(-1/2)
*. Scratch: Should atleast be: 2* Dimension of matrix + 6 times largest block
*
*. Obtain S^1/2, S^-1/2
*
        LEN_1 =  NDIM_1EL_MAT(1,NOBPS,NOBPS,NSMOB,0)
        CALL MEMMAN(KLSQRT,LEN_1,'ADDL  ',2,'SQRT  ')
        CALL MEMMAN(KLSQRTI,LEN_1,'ADDL  ',2,'SQRTI ')
        CALL MEMMAN(KLMAT,LEN_1,'ADDL  ',2,'MAT   ')
        CALL MEMMAN(KLMAT2,LEN_1,'ADDL  ',2,'MAT2  ')
        CALL MEMMAN(KLMAT3,LEN_1,'ADDL  ',2,'MAT3  ')
        NOB_MAX = IMNMX(NOBPS,NSMOB,2)
        LSCR = 6*NOB_MAX**2
        CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'LSQRT ')
        CALL COPVEC(S,WORK(KLMAT),LEN_1)
C            SQRT_BLMAT(A,NBLK,LBLK,ITASK,ASQRT,AMSQRT,SCR,ISYM)
        CALL SQRT_BLMAT(WORK(KLMAT),NSMOB,NOBPS,2,
     &       WORK(KLSQRT),WORK(KLSQRTI),WORK(KLSCR),0)
*
* ==========================================
* Exp( S^1/2 Kappa A S^1/2) in WORK(KLMAT2)
* ==========================================
*
C  TRAN_SYM_BLOC_MAT4(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
         CALL TRAN_SYM_BLOC_MAT4(KAPPA_A,WORK(KLSQRT),WORK(KLSQRT),
     &        NSMOB,NOBPS,NOBPS,WORK(KLMAT),WORK(KLSCR),0)
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' The matrix  S^1/2 Kappa A S^1/2 '
           CALL APRBLM2(WORK(KLMAT),NOBPS,NOBPS,NSMOB,0)
         END IF
* Exp(S^1/2) Kappa A S^1/2)
         LSCR_EXP = 4*NOB_MAX**2 + 3*NOB_MAX
C?       WRITE(6,*) ' LSCR_EXP, NOB_MAX = ', LSCR_EXP, NOB_MAX
         CALL MEMMAN(KLSCR_EXP,LSCR_EXP,'ADDL  ',2,'SCR_EX')
         DO ISYM = 1, NSMOB
           IF(ISYM .EQ.1) THEN
             IOFF = 1
           ELSE
             IOFF = IOFF + NOBPS(ISYM-1)**2
           END IF
* Exp(S^1/2) Kappa A S^1/2) in KLMAT2
C                EXPMA(EMA,A,NDIM,SCR,ISUB)
           CALL EXPMA(WORK(KLMAT2+IOFF-1),WORK(KLMAT+IOFF-1),
     &          NOBPS(ISYM),WORK(KLSCR_EXP),0)
C?         WRITE(6,*) ' After EXPMA '
         END DO
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' The matrix Exp( S^1/2 Kappa A S^1/2) '
           CALL APRBLM2(WORK(KLMAT2),NOBPS,NOBPS,NSMOB,0)
         END IF
*
* ===========================================
* Exp( S^1/2 Kappa S S^1/2) in WORK(KLMAT3)
* ===========================================
*
C  TRAN_SYM_BLOC_MAT4(AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
*. S^1/2 Kappa S S^1/2 in KLMAT
         CALL TRAN_SYM_BLOC_MAT4(KAPPA_S,WORK(KLSQRT),WORK(KLSQRT),
     &        NSMOB,NOBPS,NOBPS,WORK(KLMAT),WORK(KLSCR),0)
C?         WRITE(6,*) ' After TRAN_SYM_BLOC_MAT(2) '
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' The matrix  S^1/2 Kappa S S^1/2 '
           CALL APRBLM2(WORK(KLMAT),NOBPS,NOBPS,NSMOB,0)
         END IF
*
* Exp( S^1/2 Kappa S S^1/2) in KLMAT3
         DO ISYM = 1, NSMOB
           IF(ISYM .EQ.1) THEN
             IOFF = 1
           ELSE
             IOFF = IOFF + NOBPS(ISYM-1)**2
           END IF
C               EXP_MAS(EMA,A,NDIM,SCR)
           CALL EXP_MAS(WORK(KLMAT3+IOFF-1),WORK(KLMAT+IOFF-1),
     &          NOBPS(ISYM),WORK(KLSCR_EXP))
C?         WRITE(6,*) ' After EXP_MAS' 
         END DO
*
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' The matrix Exp( S^1/2 Kappa S S^1/2) '
           CALL APRBLM2(WORK(KLMAT3),NOBPS,NOBPS,NSMOB,0)
         END IF
* Exp( S^1/2) Kappa A S^1/2) Exp( S^1/2) Kappa S S^1/2) in KLMAT
C      SUBROUTINE MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
           CALL MULT_BLOC_MAT(WORK(KLMAT),WORK(KLMAT2),WORK(KLMAT3),
     &          NSMOB,NOBPS,NOBPS,NOBPS,NOBPS,NOBPS,NOBPS,0)
*. Premultipy with S^-1/2 and save on KLMAT3
           CALL MULT_BLOC_MAT(WORK(KLMAT3),WORK(KLSQRTI),WORK(KLMAT),
     &          NSMOB,NOBPS,NOBPS,NOBPS,NOBPS,NOBPS,NOBPS,0)
*. Postmultiply with S^1/2 and save in EXPMKS
           CALL MULT_BLOC_MAT(EXPMKS,WORK(KLMAT3),WORK(KLSQRT),
     &          NSMOB,NOBPS,NOBPS,NOBPS,NOBPS,NOBPS,NOBPS,0)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Matrix Exp(-K_A S) Exp(-K_S S) '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL APRBLM2(EXPMKS,NOBPS,NOBPS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2, 'GTEMKS')
*
      RETURN
      END
      SUBROUTINE EXP_MAS(EMA,A,NDIM,SCR)
*
* Expontial of minus a symmetric matrix A
* The matrix is given in complete form
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION A(NDIM,NDIM)
*. Output
      DIMENSION EMA(NDIM,NDIM)
*. Scratch: Length should be 2*NDIM**2 + NDIM*(NDIM+1)/2+ NDIM
      DIMENSION SCR(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from EXP_MAS '
       WRITE(6,*) ' ================= '
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Symmetrix matrix to be exponentialized '
        CALL WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
* Diagonalize matrix A
*
      KLX = 1
      KLSCR = KLX + NDIM**2
      KLMAT2 = KLSCR + NDIM*(NDIM+1)/2
      KLFREE = KLMAT2 + NDIM*NDIM
*
*. Obtain eigenvalues and eigenvectors of A
C          DIAG_SYM_MAT(A,X,SCR,NDIM,ISYM)
      CALL DIAG_SYM_MAT(A,SCR(KLX),SCR(KLSCR),NDIM,0)
*. Eigenvalues have been returned in SCR(KLSCR) and the eigenvectors V
*  in SCR(KLX)
*. The exponential of the eigenvalues -and remember the - from Exp(-A)
      DO I = 1, NDIM
       SCR(KLSCR-1+I) = EXP(-SCR(KLSCR-1+I))
      END DO
* V Exp(eigenvalues)
      DO J = 1, NDIM
       EPSILJ = SCR(KLSCR-1+J)
       CALL COPVEC(SCR(KLX + (J-1)*NDIM),
     &             SCR(KLMAT2+(J-1)*NDIM),NDIM)
       CALL SCALVE(SCR(KLMAT2+(J-1)*NDIM),EPSILJ,NDIM)
      END DO
* V Exp(eigenvalues) V+
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(EMA,SCR(KLMAT2),SCR(KLX),
     &            NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,
     &            FACTORC,FACTORAB,2)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Exponential of symmetrix matrix '
       CALL WRTMAT(EMA,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE ORB_EXCIT_INT_SPACE(IORBSPC,ITOTSYM,
     &           NOOEXC,IOOEXC,NUMONLY,IOFF_EXC,
     &           I_RESTRICT_SUPSYM,MO_SUPSYM)
*
* Number of orbital excitations of symmetry ITOTSYM in orbitals space 
* IORBSPC.
* NUMONLY = 1 => Only number is calculated
*         = 0 => Also the excitations are set up, starting at IOFF_EXC
*
* Jeppe Olsen, July 19, 2011, the IOFF parameter added June 2012
* Last modification; Jeppe Olsen; June 3 2013; Supersymmetry added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
*.Input
      INTEGER MO_SUPSYM(*)
*. Output
      INTEGER IOOEXC(2,*)
*
      NTEST = 10
      NOOEXC = 0
*. First orbital of space  IORBSPC
      IOFF = NINOB + 1
      DO IGAS = 0, IORBSPC-1
        IOFF = IOFF + NOBPT(IGAS)
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' Offset for orbital excitations ',
     & IOFF
      NORB = NOBPT(IORBSPC)
      DO IORB = IOFF, IOFF + NORB - 1
        DO JORB = IOFF, IORB - 1
          ISM = ISMFTO(IORB)
          JSM = ISMFTO(JORB)
          IF(NTEST.GE.100) WRITE(6,*) ' IORB, JORB, ISM, JSM = ',
     &    IORB, JORB, ISM, JSM
          IF(MULTD2H(ISM,JSM).EQ.ITOTSYM) THEN
            IMOKAY2 = 1
            IF(I_RESTRICT_SUPSYM.EQ.1) THEN
*. Check that supersymmetries are identical
              IF(MO_SUPSYM(IREOTS(IORB)).NE.MO_SUPSYM(IREOTS(JORB)))THEN
               IMOKAY2 = 0
               IF(NTEST.GE.10) THEN
                 WRITE(6,*) 
     &           ' Excitation eliminated by supersym: IORB, JORB = ',
     &             IORB, JORB
               END IF
              END IF
            END IF! Supersymmetry restrictions are active
            IF(IMOKAY2.EQ.1) THEN
              NOOEXC = NOOEXC + 1
              IF(NUMONLY.EQ.0) THEN
                IOOEXC(1,NOOEXC+IOFF_EXC-1) = IORB
                IOOEXC(2,NOOEXC+IOFF_EXC-1) = JORB
              END IF
            END IF
          END IF ! Symmetry was right
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Number of active- active orbital excitations ', 
     & NOOEXC
      END IF
      IF(NTEST.GE.100) THEN
       IF(NUMONLY.EQ.0) THEN
         WRITE(6,*)  ' And the orbital excitations '
         CALL IWRTMA(IOOEXC(1,IOFF_EXC),2,NOOEXC,2,NOOEXC)
       END IF
      END IF
*
      RETURN
      END
      SUBROUTINE E1_VB_FROM_ACTMAT(E1,IOOEXC,NOOEXC,E,RHOA,RHOB)
*
* Obtain VB gradient in active space from densities 
* RHOB = <c!a+iaj!c(bio)>/<0!0>, RHOA = <c!a+aj!hc(bio)>/<0!0>
*. (note that the densities are in the original basis)
*
*. Jeppe Olsen, July19, 2011
*     
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER IOOEXC(2,NOOEXC)
      DIMENSION RHOA(NACOB,NACOB),RHOB(NACOB,NACOB)
*. Output
      DIMENSION E1(2*NOOEXC)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from E1_VB_FROM_ACTMAT '
        WRITE(6,*) ' ============================'
        WRITE(6,*) 
        WRITE(6,*) ' Energy = ', E
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' <0!E(ij)!H0>/<0!0> '
        CALL WRTMAT(RHOA,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*)
        WRITE(6,*) ' <0!E(ij)!0>/<0!0> '
        CALL WRTMAT(RHOB,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*)
      END IF
*. The antisymmetric part of the gradient
      DO JOO = 1, NOOEXC
        IORB = IOOEXC(1,JOO)-NINOB
        JORB = IOOEXC(2,JOO)-NINOB
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' JOO, IORB, JORB = ', IORB, JORB
*. Antisymmetric part
        E1(JOO) = 2.0D0*(RHOA(IORB,JORB)-RHOA(JORB,IORB))
*. Symmetric part
        E1(JOO+NOOEXC) = 
     &  -2.0D0*     (RHOA(IORB,JORB)+RHOA(JORB,IORB)
     &           -E*(RHOB(IORB,JORB)+RHOB(JORB,IORB)))
      END DO
*
      IF(NTEST.GE.100) THEN
*
       WRITE(6,*) ' Active-active gradient for nonorthogonal MCSCF '
       WRITE(6,*) ' ==============================================='
       WRITE(6,*)   
       CALL WRTMAT(E1,1,2*NOOEXC,1,2*NOOEXC)
      END IF
*
      RETURN
      END
      SUBROUTINE DO_ORBTRA(IDOTRA,IDOFI,IDOFA,
     &           IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
*
* Perform orbital transformations on integrals and Inactive/active Fock
* matrix
* 
* IDOTRA = 1 => Transformed one- and two-electron integrals
* IDOFI  = 1 => Inactive Fock-matrix
* IDOFA  = 1 => Active Fock-matrix
*
* Jeppe Olsen, July 2011 - In a hotel room in Santiago de Compostella
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
      IDUM = 0
      CALL LUCIAQENTER('ORBTR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ORBTRA')
*
      NTEST = 0000
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' Info from DO_ORBTRA '
       WRITE(6,*) ' ====================='
       WRITE(6,*) 
       WRITE(6,*) ' Tasks: IDOTRA, IDOFI, IDOFA = ', 
     &                     IDOTRA, IDOFI, IDOFA
       IF(IDOTRA.EQ.1) THEN
       WRITE(6,*)  ' IE2LIST_IN, IOCOBTP_IN, INTSM_IN = ',
     &               IE2LIST_IN, IOCOBTP_IN, INTSM_IN
       END IF
      END IF! NTEST .ge. 100
* 
      IE2LIST_A = IE2LIST_IN
      IOCOBTP_A = IOCOBTP_IN
      INTSM_A = INTSM_IN
*
      call ga_sync()
      IF(IDOTRA.EQ.1) THEN
*. Perform one- and two-electron transformations.
* The pointers to the mo-ao transformation matrices KKCMO_X, X=I,J,K,L
* must have been set up outside.
        CALL PREPARE_2EI_LIST
        CALL TRAINT
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' one-electron transformed integrals'
          WRITE(6,*) ' ================================='
          IPACK_H1 = IE1_CCSM_G(IE2LIST_IN)
CNW       CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,IPACK_H1)
          call ga_print(KINT1)
        END IF
*
      END IF! Integral transformation should be performed
*
      call ga_sync()
      IF(IDOFI.EQ.1) THEN
*
*.      Calculate inactive Fock matrix in basis defined by  KKCMI, KKCMJ
*       ================================================================
*
*. Use AO integrals in KINT_2EMO
*
        IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
        KINT2_FSAVE = KINT2_A(IE2ARR_F)
        KINT2_A(IE2ARR_F) = KINT_2EMO
*. The permutational symmetry of the inactive Fock-matrix is inherited from
*. the complex conjugation symmetry of the one-electron integrals
        IPACK_F = IE1_CCSM_G(IE2LIST_IN)
*
        CALL FI_FROM_INIINT_G(KFI,dbl_mb(KKCMO_I),dbl_mb(KKCMO_J),
     &                        KINT1,ECORE_HEX,3,IPACK_F)
        ECORE = ECORE_ORIG + ECORE_HEX
        IF(NTEST.GE.100)
     &  WRITE(6,*) ' Updated core energy =  ', ECORE
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Inactive Fock-matrix '
          call ga_print(kfi)
CNW       CALL APRBLM2(WORK(KFI),NTOOBS,NTOOBS,NSMOB,IPACK_F)
        END IF
*. And clean up
        KINT2_A(IE2ARR_F) = KINT2_FSAVE
      END IF !  FI  should be calculated
*
      call ga_sync()
      IF(IDOFA.EQ.1) THEN
*
*.      Calculate active Fock matrix in basis defined by  KKCMI, KKCMJ
*        =============================================================
*
*
*. Use AO integrals in KINT_2EMO
*
        IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
        KINT2_FSAVE = KINT2_A(IE2ARR_F)
        KINT2_A(IE2ARR_F) = KINT_2EMO
*. The permutational symmetry of the inactive Fock-matrix is inherited from
*. the complex conjugation symmetry of the one-electron integrals
        IPACK_F = IE1_CCSM_G(IE2LIST_IN)
*
* A bit dirty: I will use IPACK_F to decide whether it is an 
* normal or bio-calculation- will probably give my trouble later..
        IF(IPACK_F.EQ.0) THEN
         IBIO_CALC = 1
        ELSE
         IBIO_CALC = 0
        END IF
*
C            FA_FROM_INIINT(FA,CINI,CINIB,D,IPACK)
        IF(IBIO_CALC.EQ.1) THEN
*. transform RHO1 to bio-actual MO basis
*
*. Obtain first in symmetry block form
          LEN_R = NDIM_1EL_MAT(1,NACOBS,NACOBS,NSMOB,0)
          CALL MEMMAN(KLRHO1,NACOB**2,'ADDL  ',2,'RHO1L ')
          CALL MEMMAN(KLRHO1B,NACOB**2,'ADDL  ',2,'RHO1S ')
          CALL MEMMAN(KLCBIOA,LEN_R,'ADDL  ',2,'CBIOAC')
C              REORHO1(RHO1I,RHO1O,IRHO1SM)
          CALL REORHO1(dbl_mb(KRHO1),dbl_mb(KLRHO1),1,1)
*.  Obtain CBIO over active orbitals only
C            EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT(A,AGAS,I_EX_OR_CP)
          CALL EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT
     &         (dbl_mb(KCBIO),dbl_mb(KLCBIOA),1)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' CBIO in active orbitals '
            CALL APRBLM2(dbl_mb(KLCBIOA),NACOBS,NACOBS,NSMOB,0)
          END IF
          CALL TR_BIOMAT(dbl_mb(KLRHO1),dbl_mb(KLRHO1B),dbl_mb(KLCBIOA),
     &                   NACOBS,1,2,1,1)
*. Transfer back to full matrix over active orbitals
          CALL REORHO1(dbl_mb(KLRHO1),dbl_mb(KLRHO1B),1,2)
        ELSE
          KLRHO1 = KRHO1
        END IF
*
        CALL FA_FROM_INIINT(KFA,dbl_mb(KKCMO_I),dbl_mb(KKCMO_J),
     &                      dbl_mb(KLRHO1),IPACK_F)
*. And clean up
        KINT2_A(IE2ARR_F) = KINT2_FSAVE
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Active Fock-matrix '
          call ga_print(kfa)
CNW       CALL APRBLM2(WORK(KFA),NTOOBS,NTOOBS,NSMOB,IPACK_F)
        END IF
      END IF ! Active Fock matrix should be calculated
*
      call ga_sync()
      IF(NTEST.GE.100) WRITE(6,*) ' Leaving DO_ORBTRA '
*
      CALL LUCIAQEXIT('ORBTR')
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ORBTRA')
*
      RETURN
      END
      SUBROUTINE GET_INIMO(CMO_INI)
*
* Obtain initial set of Molecular orbitals in CMO_INI as specified by 
* parameters INI_MO_TP,INI_MO_ORT, INI_ORT_VBGAS in crun
*
*   Two steps : 1) Obtain a set of (nonorthogonal) initial orbitals
*                  according to INI_MO_TP
*               2) Perform (partial) orthonormalization to obtain 
*                  Final initial orbitals according to INI_MO_ORT,
*                  and INI_ORT_VBGAS,IGAS_SEL
*
* The INI_MO_TP parameter defines the raw (nonorthogonal) initial orbitals:
*
* INI_MO_TP = 1 => Unit matrix
* INI_MO_TP = 2 => Rotate orbitals from environment so 
*                    Diagonal block in GAS IGAS_SEL is diagonal
* INI_MO_TP = 3 => Use orbitals read in from environment
* INI_MO_TP = 4 => Read in fragment orbitals
* INI_MO_TP = 5 => Read in from LUCINF_O
*
*
* INI_MO_ORT = 0 => No orthonormalization
*            = 1 => symmetric orthogonalization
*            = 2 => orthonormalization by biagonalization  
*
* INI_ORT_VBGAS = 0 => No orthonormalization of VB gas space
*               = 1 => Orthonornormalization of VB gas space according to
*                      INI_MO_ORT
*
* Jeppe Olsen, Restructuring some code in a Hotel room in Santiago De
*              Compostella, July 2011
*              June 2012, INI_MO_TP = 5 added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'  
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'fragmol.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
*.
      CHARACTER*6 CSAVE
*. Output
      DIMENSION CMO_INI(*)
*
      IDUM = 0 
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' Info from GET_INIMO: '
        WRITE(6,*) ' ====================='
        WRITE(6,*)
      END IF
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'INIMO ')
*
*. 0: Obtain some input information if required
*
      IF(INI_MO_TP.EQ.4) THEN
*
*.      Set up fragment information
*
        IF(NFRAG_TP.EQ.0) THEN
          WRITE(6,*) 
     &    ' Input orbitals from fragment MOs requested(INI_MO_TP=4)'
          WRITE(6,*) 
     &    ' But no fragment information provided (keyword: MOFRAG)'
          WRITE(6,*)   ' Specify keyword MOFRAG '
          STOP         ' Specify keyword MOFRAG '
        ELSE
         CALL MOINF_FRAG 
        END IF
      END IF ! Iform = 4
*
      LEN_1F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLCMOAO1,LEN_1F,'ADDL  ',2,'CMOAO1')
      CALL MEMMAN(KLCMOAO2,LEN_1F,'ADDL  ',2,'CMOAO2')
*
      IF(INI_MO_TP.EQ.2.OR.INI_MO_TP.EQ.3) THEN
*
*.      Obtain MOAO transformation matrix from environment 
*
        CALL GET_CMOAO_ENV(WORK(KLCMOAO1))
      END IF
*
      IF(INI_MO_TP.EQ.5) THEN
*
*. Read in from LUCINF_O which is a fort.91 output file, but
*. perhaps from another geometry.
*
*. a bit of dirty dancing: let LUCINF_O be the standard fort.91 
*. for a few microseconds.
*. Obtain a free unit-number
       LU91_SAVE = LU91
       CALL FILEMAN_MINI(LU91,'ASSIGN')
       OPEN(LU91,STATUS='OLD',FORM='FORMATTED',FILE='LUCINF_O')
*. Fool also environment to think it is LUCIA
       CSAVE = ENVIRO
       ENVIRO(1:6) = 'LUCIA '
*. Obtain CMO as usual from environment - with changed LU91
*. 
       CALL GET_CMOAO_ENV(WORK(KLCMOAO1))
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' MOAO for INI_MO_TP = 5 '
         CALL APRBLM2(WORK(KLCMOAO1),NTOOBS,NTOOBS,NSMOB,0)
       END IF
*. And restore order
       CLOSE(LU91,STATUS='KEEP')
       CALL FILEMAN_MINI(LU91,'FREE  ')
       LU91 = LU91_SAVE
       ENVIRO = CSAVE
      END IF
*
*
*. 1: Generate/Read in the 'initial initial' orbitals and store in CMOAO2
*
*. The split of work between current routine and PREPARE_CMOAO_INI is 
*. strange, but works..
C     PREPARE CMOAO_INI(INI_MO_TP_L, CMOAO_OUT,CMOAO_IN,IVBGAS)
      CALL PREPARE_CMOAO_INI
     &     (INI_MO_TP,WORK(KLCMOAO2),WORK(KLCMOAO1),
     &     NORTCIX_SCVB_SPACE)
      CALL COPVEC(WORK(KLCMOAO2),WORK(KLCMOAO1),LEN_1F)
*
*. 2. Orthonormalize parts of the orbital spaces
*
*.
*. Orthogonalize Active to inactive and secondary to active- always done
      INTER_ORT = 1
*. Between GA spaces
      IF(INI_MO_ORT.EQ.0) THEN
        INTERGAS_ORT = 0
        INI_ORT_VBGASL = 0
      ELSE
        INTERGAS_ORT = 1
        INI_ORT_VBGASL = INI_ORT_VBGAS
      END IF
*. Intragas orthogonalization
      INTRAGAS_ORT = INI_MO_ORT
*. Orthogonalization in VB space- defined by parameter INI_MO_ORT
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,4I4)')
     &  ' INTER_ORT, INTERGAS_ORT, INTRAGAS_ORT, INI_ORT_VBGASL',
     &    INTER_ORT, INTERGAS_ORT, INTRAGAS_ORT, INI_ORT_VBGASL
      END IF
C     ORT_ORB(CMOAO_IN, CMOAO_OUT, INTER_ORT,INTERGAS_ORT,
C    &        INTRAGAS_ORT,IORT_VBSPC)
      CALL ORT_ORB(WORK(KLCMOAO1),CMO_INI,INTER_ORT, 
     &     INTERGAS_ORT,INTRAGAS_ORT,INI_ORT_VBGASL)
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Expansion of final initial MOs in AOs '
        WRITE(6,*) ' ======================================'
        CALL APRBLM_F7(CMO_INI,NTOOBS,NTOOBS,NSMOB,0)
C       CALL PRINT_CMOAO(CMO_INI)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'INIMO ')
*
      RETURN
      END
      SUBROUTINE BLK_CHECK_UNI_MAT
     &           (UNI,NBLK,LBLK,XMAX_DIFF_DIAG,XMAX_DIFF_OFFD)
*
* A full blocked matrix UNI is given. Find largest deviation of 
* matrix from unit matrix as 
*    The largest deviation of diagonal element from one
*    The largest deviation of block-diagonal element from zero
*
* Jeppe Olsen, May 2012
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION UNI(*)
      INTEGER LBLK(NBLK)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ==========================='
        WRITE(6,*) ' Info from BLK_CHECK_UNI_MAT'
        WRITE(6,*) ' ==========================='
      END IF
*
      IB = 1
      XMAX_DIFF_DIAG = 0.0D0
      XMAX_DIFF_OFFD = 0.0D0
*
      DO IBLK = 1, NBLK
        L = LBLK(IBLK)
C       CHECK_UNIT_MAT(UNI,NDIM,XMAX_DIFF_DIAG,XMAX_DIFF_OFFD)
        CALL CHECK_UNIT_MAT(UNI(IB),L,XDIAG_LOC, XOFFD_LOC,0)
        XMAX_DIFF_DIAG = MAX(XMAX_DIFF_DIAG,XDIAG_LOC)
        XMAX_DIFF_OFFD = MAX(XMAX_DIFF_OFFD,XOFFD_LOC)
        IB = IB + L**2
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' Deviations of block matrix from unit matrix: '
        WRITE(6,*) 
     &  '   Largest deviation of diagonal element from 1:',
     &  XMAX_DIFF_DIAG
        WRITE(6,*) 
     &  ' Largest deviation of off-diagonal element from 1:',
     &    XMAX_DIFF_OFFD
      END IF
*
      RETURN
      END
      SUBROUTINE CHECK_UNIT_MAT(UNI,NDIM,XMAX_DIFF_DIAG,XMAX_DIFF_OFFD,
     &           ISYM)
*
* A matrix UNI is given. Check difference between UNI and UNIT matrix
* and report in:
*     XMAX_DIFF_DIAG: Max absolute difference between between diagonal 
*                     element and 1
*     XMAX_DIFF_OFFD: Max absolute difference between off diagonal and zero
*
*. Jeppe Olsen, July 2011 (Thought I had written this routine before...)
*  Last modification; Feb 27, 2013; Jeppe Olsen; ISYM added
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION UNI(*)
*. Diagonal element
      XMAX_DIFF_DIAG = 0.0D0
      DO I = 1, NDIM
         IF(ISYM.EQ.0) THEN
           II = (I-1)*NDIM + I
         ELSE
          II = I*(I-1)/2 + I
        END IF
        XMAX_DIFF_DIAG = MAX(XMAX_DIFF_DIAG,ABS(UNI(II)-1.0D0))
      END DO
*. Off diagonal elements
      XMAX_DIFF_OFFD = 0.0D0
      DO I = 1, NDIM
       DO J = 1, I-1
        IF(ISYM.EQ.0) THEN
          JI = (I-1)*NDIM + J
          IJ = (J-1)*NDIM + I
          XMAX_DIFF_OFFD = MAX(XMAX_DIFF_OFFD,ABS(UNI(IJ)),ABS(UNI(JI)))
        ELSE
         IJ = I*(I-1)/2 + J
         XMAX_DIFF_OFFD = MAX(XMAX_DIFF_OFFD,ABS(UNI(IJ)))
        END IF
       END DO
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Comparison of matrix with unit matrix: '
       WRITE(6,*) '   Largest deviation of diagonal elements ', 
     & XMAX_DIFF_DIAG
       WRITE(6,*) '   Largest deviation of of-diagonal elements ', 
     & XMAX_DIFF_OFFD
      END IF
*
      RETURN
      END
      SUBROUTINE TR_BIOMAT(XIN,XOUT,CBIO,NORB_PSM, 
     &            INB_IN,INB_OUT,JNB_IN,JNB_OUT)
*
* An orbital matrix XIN(I,J) is given  in symmetry blocked form
* with NORB_PSM orbitals per symmetry
* INB_IN = 1 => I is in normal basis
*        = 2 => I is in bioorthogonal basis
* JNB_IN = 1 => J is in normal basis
*        = 2 => J is in bioorthogonal basis
* Obtain the matrix in the representation XOUT(I,J) defined by 
* 
* INB_OUT = 1 => I is in normal basis
*         = 2 => I is in bioorthogonal basis
* JNB_OUT = 1 => J is in normal basis
*         = 2 => J is in bioorthogonal basis
* The matrix CBIO giving the transformation from the normal to the 
* bioorthogonal basis is in the same basis.
*
* Note: The use of locally defined NORB_PSM, allows the restriction 
*       of the matrice to for example the active orbitals.
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION XIN(*)
      DIMENSION NORB_PSM(NSMOB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info form TR_BIOMAT '
        WRITE(6,*) ' =================== '
        WRITE(6,'(A,4I2)') ' INB_IN,INB_OUT,JNB_IN,JNB_OUT = ',
     &                       INB_IN,INB_OUT,JNB_IN,JNB_OUT
*
        WRITE(6,*) ' NORB_PSM = '
        CALL IWRTMA3(NORB_PSM,1,NSMOB,1,NSMOB)
*
        WRITE(6,*) ' The Input Cbio matrix '
        CALL APRBLM2(CBIO,NORB_PSM,NORB_PSM,NSMOB,0)
      END IF
*
*. Check that input parameters are in range
*
      INB_IN_OK = 1
      JNB_IN_OK = 1
      INB_OUT_OK = 1
      JNB_OUT_OK = 1
      IF(1.GT.INB_IN.OR.INB_IN.GT.2) INB_IN_OK = 0
      IF(1.GT.JNB_IN.OR.JNB_IN.GT.2) JNB_IN_OK = 0
      IF(1.GT.INB_OUT.OR.INB_OUT.GT.2) INB_OUT_OK = 0
      IF(1.GT.JNB_OUT.OR.JNB_OUT.GT.2) JNB_OUT_OK = 0
*
      IF(INB_IN_OK.EQ.0.OR.JNB_IN_OK.EQ.0.OR.
     &   INB_OUT_OK.EQ.0.OR.JNB_OUT_OK.EQ.0) THEN
       WRITE(6,*) ' Error in input to TR_BIOMAT'
       WRITE(6,*) ' Input parameter out or range (1,2)'
       WRITE(6,'(A,4(2X,I2))') 
     & ' INB_IN,JNB_IN,INB_OUT,JNB_OUT = ',
     &   INB_IN,JNB_IN,INB_OUT,JNB_OUT
       STOP ' Error in input to TR_BIOMAT' 
      END IF     
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TR_BIO')
*
      LEN_1F = NDIM_1EL_MAT(1,NORB_PSM,NORB_PSM,NSMOB,0)
*. Local copy of CBIO
      CALL MEMMAN(KLCBIO,LEN_1F,'ADDL  ',2,'CBIOL ')
      CALL COPVEC(CBIO,WORK(KLCBIO),LEN_1F)
*
      NOBS_MX = IMNMX(NORB_PSM,NSMOB,2)
      LSCR = 2*NOBS_MX**2
      KLCBIOINV = 0
*. Obtain transformation from BIO to normal basis if required
      IF(INB_IN.EQ.2.AND.INB_OUT.EQ.1.OR.
     &   JNB_IN.EQ.2.AND.JNB_OUT.EQ.1) THEN
         CALL MEMMAN(KLCBIOINV,LEN_1F,'ADDL  ',2,'CBIINV')
         CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'CBIOSC')
*
C             INV_BLKMT(A,AINV,SCR,NBLK,LBLK,IPROBLEM)
         CALL INV_BLKMT(WORK(KLCBIO),WORK(KLCBIOINV),WORK(KLSCR),
     &        NSMOB,NORB_PSM,IPROBLEM)
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Inverted CBIO '
           CALL APRBLM2(WORK(KLCBIOINV),NORB_PSM,NORB_PSM,NSMOB,0)
         END IF
         IF(IPROBLEM.NE.0) THEN
           WRITE(6,*) ' Problem inverting CBIO(MO,MO) '
         END IF
      END IF
*
*. Local pointers to pointers to transformations matrices for I and J
*
      IF(INB_IN.EQ.INB_OUT) THEN
       KKLI = 0
      ELSE 
       IF(INB_IN.EQ.1.AND.INB_OUT.EQ.2) THEN
* Normal => BIO
        KKLI = KLCBIO
       ELSE
        KKLI = KLCBIOINV
       END IF
      END IF
*
      IF(JNB_IN.EQ.JNB_OUT) THEN
       KKLJ = 0
      ELSE 
       IF(JNB_IN.EQ.1.AND.JNB_OUT.EQ.2) THEN
* Normal => BIO
        KKLJ = KLCBIO
       ELSE
        KKLJ = KLCBIOINV
       END IF
      END IF
*. And do the transformation as requested
      IF(INB_IN.EQ.INB_OUT.AND.JNB_IN.EQ.JNB_OUT) THEN
*       No transformation, just copy
        CALL COPVEC(XIN,XOUT,LEN_1F)
      ELSE IF( INB_IN.NE.INB_OUT.AND.JNB_IN.EQ.JNB_OUT) THEN
*.      Transformation of first index I
C            MULT_BLOC_MAT(C,A,B,
C            NBLOCK,LCROW,LCCOL,LAROW,LACOL,LBROW,LBCOL,ITRNSP)
        CALL MULT_BLOC_MAT(XOUT,WORK(KKLI),XIN,
     &       NSMOB,NORB_PSM,NORB_PSM,NORB_PSM,NORB_PSM,NORB_PSM,
     &       NORB_PSM,1)
      ELSE IF(INB_IN.EQ.INB_OUT.AND.JNB_IN.NE.JNB_OUT) THEN
*       Transformation of second index, J
        CALL MULT_BLOC_MAT(XOUT,XIN,WORK(KKLJ), 
     &       NSMOB,NORB_PSM,NORB_PSM,NORB_PSM,NORB_PSM,NORB_PSM,
     &       NORB_PSM,0)
      ELSE
*. Transform both I and J indeces
C           TRAN_SYM_BLOC_MAT4
C           (AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
         CALL MEMMAN(KLSCRTRA,LSCR,'ADDL  ',2,'SCRTRA')
         CALL TRAN_SYM_BLOC_MAT4(XIN,WORK(KKLI),WORK(KKLJ),
     &        NSMOB,NORB_PSM,NORB_PSM,XOUT,WORK(KLSCRTRA),0)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from TR_BIOMAT'
        WRITE(6,*) ' ====================== '
        WRITE(6,*) 
        WRITE(6,'(A,4(2X,I3))') ' INB_IN, JNB_IN, INB_OUT, JNB_OUT =',
     &                           INB_IN, JNB_IN, INB_OUT, JNB_OUT 
        WRITE(6,*) ' Input matrix: '
        CALL APRBLM2(XIN,NORB_PSM,NORB_PSM,NSMOB,0)
        WRITE(6,*) ' Output matrix: '
        CALL APRBLM2(XOUT,NORB_PSM,NORB_PSM,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TR_BIO')
      RETURN
      END
      SUBROUTINE EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT
     &           (A,AGAS,I_EX_OR_CP)
*
* A symmetryblocked (not lower half packed) matrix A over orbitals is given
* Extract all blocks referring to the GASpaces (i.e. 1-ngas)
*
* Matrix is assumed total symmetric wrt pointgroup
*
* I_EX_OR_CP = 1 => Extract from A to AGAS
* I_EX_OR_CP = 1 => Copy from AGAS to A
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Specific input or output
      DIMENSION A(*), AGAS(*)
*
      DO ISYM = 1, NSMOB
       IF(ISYM.EQ.1) THEN
        IOFF_IN = 1
        IOFF_OUT = 1
       ELSE
        IOFF_IN = IOFF_IN + NTOOBS(ISYM-1)**2
        IOFF_OUT = IOFF_OUT + NACOBS(ISYM-1)**2
       END IF
*
       IOFF = NINOBS(ISYM)+1
       NIA= NACOBS(ISYM)
       NIT= NTOOBS(ISYM)
*
       DO J = 1, NIA
         DO I = 1, NIA
           IJ_OUT = IOFF_OUT -1 + (J-1)*NIA + I 
           IJ_IN  = IOFF_IN -1 
     &            + (IOFF+J-1-1)*NIT + IOFF+I-1
           IF(I_EX_OR_CP.EQ.1) THEN
             AGAS(IJ_OUT) = A(IJ_IN)
           ELSE
             A(IJ_IN) = AGAS(IJ_OUT)
           END IF
         END DO
       END DO
      END DO ! End of loop over symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Submatrix Over active orbitals' 
         CALL APRBLM2(AGAS,NACOBS,NACOBS,NSMOB,0)
         WRITE(6,*) ' Full matrix '
         CALL APRBLM2(A,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE VB_GRAD_ORBVBSPC(NOOEXC_AA,IOOEXC_AA,E1,C,
     &           VEC1_CSF,VEC2_CSF)
*
* Obtain gradient over orbitals in active space
*
* E1(A)(IJ) = 2 ( <0!(E(ij) - E(ji))H!0>  
* E1(S)(IJ) =-2 ( <0!(E(ij) + E(ji))(H-E)!0>
*
* The number of active-active excitations is NOOEXC_AA
* and the corresponding excitations are IOOEXC_AA 
*
* So to obtain gradient 
* 1: construct bioorthogonal expansion of S = H!0> and !0>
* 2: Set up density matrices <0!E(ij)!s> <0!E(ij)!0> 
*    where i is in biobase and j in normal
* 3: Transform density matrices to standard basis
* To accomplish 1, the sigma routine is called with the current set of 
* CI coefficients
*
* The current CI coefficients in the CSF basis are in C, where 
* VEC1_CSF, VEC2_CSF, must be able to hold these expansions
*
* This is an initial version, for initial calculations and checks
*
* Jeppe Olsen, July 2011, for the initial NORTMCSCF program
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3, 
     &       LUCBIO_SAVE, LUHCBIO_SAVE,LUC_SAVE
      REAL*8 INPRDD
*. Input
      DIMENSION C(*)
      INTEGER IOOEXC_AA(2,NOOEXC_AA)
*. Scratch 
      DIMENSION VEC1_CSF(*), VEC2_CSF(*)
*. Output
      DIMENSION E1(2*NOOEXC_AA)
*
      NTEST = 000
*. CSFs are handled explicitly, so
      NOCSF = 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==========================='
        WRITE(6,*) ' Input from VB_GRAD_ORBVBSPC'
        WRITE(6,*) ' ==========================='
        WRITE(6,*)
        WRITE(6,*) ' NOOEXC_AA = ', NOOEXC_AA
        WRITE(6,*) ' The active-active excitations '
        CALL PRINT_ORBEXC_LIST(IOOEXC_AA,0,NOOEXC_AA)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'VBGRAD')
*
      LUSCR1 = LUSC34
      LUSCR2 = LUSC35
      LUSCR3 = LUSC36
      LUCBIO_SAVE = 110
      LUHCBIO_SAVE = 111
      LUC_SAVE = 112
*
* A bit of scratch
*
      LEN_1A = NDIM_1EL_MAT(1,NACOBS,NACOBS,NSMOB,0)
      CALL MEMMAN(KLRHOA,NACOB**2,'ADDL  ',2,'RHOA  ')
      CALL MEMMAN(KLRHOB,NACOB**2,'ADDL  ',2,'RHOB  ')
      CALL MEMMAN(KLSCR ,NACOB**2,'ADDL  ',2,'SCR   ')
      CALL MEMMAN(KLCBIOA,LEN_1A,'ADDL  ',2,'CBIOAC')
*. Preparation: Obtain CBIO over active orbitals only
C          EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT(A,AGAS,I_EX_OR_CP)
      CALL EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT
     &     (WORK(KCBIO),WORK(KLCBIOA),1)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' CBIO in active orbitals '
        CALL APRBLM2(WORK(KLCBIOA),NACOBS,NACOBS,NSMOB,0)
      END IF
*
*. Sigma with the current C 
*
C          SIGMA_NORTCI(C,HC,SC,IDOHC,IDOSC)
      CALL SIGMA_NORTCI(C,VEC1_CSF,VEC2_CSF,1,1)
      IF(NTEST.GE.1000) WRITE(6,*) ' Back from SIGMA_NORTCI'
* calculate energy from vectors on file
      CHC = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE,LUHCBIO_SAVE,1,-1)
      CC  = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE, LUCBIO_SAVE,1,-1)
      EVB = CHC/CC
      IF(NTEST.GE.10) WRITE(6,*) ' Energy is ', EVB
*
*. Set up density <0! a+i(bio) aj!0(bio)> in RHOB
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' C in SD expansion '
        CALL WRTVCD(WORK(KVEC1P),LUC_SAVE,1,-1)
        WRITE(6,*) ' C(bio) in SD expansion '
        CALL WRTVCD(WORK(KVEC1P),LUCBIO_SAVE,1,-1)
        WRITE(6,*) ' HC(bio) in SD expansion '
        CALL WRTVCD(WORK(KVEC1P),LUHCBIO_SAVE,1,-1)
      END IF
      XDUM = 0.0D0
      CALL DENSI2(1 ,WORK(KLRHOB),XDUM,
     &WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE,LUCBIO_SAVE,EXPS2,
     &0,XDUM,XDUM,XDUM,XDUM,0)
*. Scale with 1/<0!0>
      FACTOR = 1.0D0/CC
      CALL SCALVE(WORK(KLRHOB),FACTOR,NACOB**2)
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Density matrix <0! a+i(bio) aj!bio 0>/<0!0> '
       CALL WRTMAT(WORK(KLRHOB),NACOB,NACOB,NACOB,NACOB)
      END IF
*. Obtain density as blocked matrix over symmetry blocks of active orbitals
C          REORHO1(RHO1I,RHO1O,IRHO1SM)
      CALL REORHO1(WORK(KLRHOB),WORK(KLSCR),1,1)
      CALL COPVEC(WORK(KLSCR),WORK(KLRHOB),LEN_1A)
*. Transform the densities from bio, normal to the normal, normal basis
C     TR_BIOMAT(XIN,XOUT,CBIO,NORB_PSM, 
C    &            INB_IN,INB_OUT,JNB_IN,JNB_OUT)
      CALL TR_BIOMAT(WORK(KLRHOB),WORK(KLSCR),WORK(KLCBIOA),
     &     NACOBS,2,1,1,1)
*. Transfer back to full matrix over active orbitals
      CALL REORHO1(WORK(KLRHOB),WORK(KLSCR),1,2)
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Density matrix <0! a+i aj!bio 0> '
       CALL WRTMAT(WORK(KLRHOB),NACOB,NACOB,NACOB,NACOB)
      END IF
*
*. Set up density <0! a+i(bio) aj!H0(bio)>   in RHOA
*
      CALL DENSI2(1 ,WORK(KLRHOA),XDUM,
     &     WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE,LUHCBIO_SAVE,EXPS2,
     &     0,XDUM,XDUM,XDUM,XDUM,0)
*. Scale with 1/<0!0>
      FACTOR = 1.0D0/CC
      CALL SCALVE(WORK(KLRHOA),FACTOR,NACOB**2)
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Density matrix <0! a+i(bio) aj!bio H0>/<0!0> '
       CALL WRTMAT(WORK(KLRHOA),NACOB,NACOB,NACOB,NACOB)
      END IF
*. Obtain density as blocked matrix over symmetry blocks of active orbitals
      CALL REORHO1(WORK(KLRHOA),WORK(KLSCR),1,1)
      CALL COPVEC(WORK(KLSCR),WORK(KLRHOA),LEN_1A)
*. Transform the densities from bio, normal to the normal, normal basis
      CALL TR_BIOMAT(WORK(KLRHOA),WORK(KLSCR),WORK(KLCBIOA),
     &     NACOBS,2,1,1,1)
*. Transfer back to full matrix over active orbitals
      CALL REORHO1(WORK(KLRHOA),WORK(KLSCR),1,2)
*. and construct the gradient
C     E1_VB_FROM_ACTMAT(E1,IOOEXC_S,NOOEXC_AA,E,RHOA,RHOB)
       CALL E1_VB_FROM_ACTMAT(E1,IOOEXC_AA,
     &      NOOEXC_AA,EVB, WORK(KLRHOA),WORK(KLRHOB))
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'VBGRAD')
      RETURN
      END
      SUBROUTINE GET_VB_VF_VBSPC_FROM_KAPPA(E1,
     &           KAPPA_A,NOOEXC_A,IOOEXC_A,
     &           KAPPA_S,NOOEXC_S,IOOEXC_S,CCI,
     &           VEC1_CSF,VEC2_CSF)
*
* Obtain gradient-like Vector function E1 in VB orbital space from 
* given Kappa and S
*
* Using method with expansion in complete VI space
*
*. It is assumed that the current MO-AO coefficients are in KMOAOIN.
* Integrals etc are overwritten, so the exit from this routine is 
* not clean.
*
*. Jeppe Olsen, July 24 2011 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'spinfo.inc'
*. Specific input
      INTEGER IOOEXC_A(2,NOOEXC_A), IOOEXC_S(2,NOOEXC_S)
      REAL*8 KAPPA_A(*), KAPPA_S(*)
*. Scratch
      DIMENSION VEC1_CSF(*),VEC2_CSF(*)
*. Output
      DIMENSION E1(NOOEXC_S+NOOEXC_A)
*
      NTEST = 100
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Input from GET_VB_VF_VBSPC_FROM_KAPPA '   
        WRITE(6,*) ' ======================================'
        WRITE(6,*) 
        WRITE(6,*) ' Input Kappa_A, Kappa_S: '
        CALL WRTMAT(KAPPA_A,1,NOOEXC_A,1,NOOEXC_A)
        WRITE(6,*)
        CALL WRTMAT(KAPPA_S,1,NOOEXC_S,1,NOOEXC_S)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTVBVF')
*
*. Obtain New MO coefficients in MOAOUT: MOAOIN* Exp(-Kappa_A S) Exp(-Kappa_S S)
*
C     NEWMO_FROM_KAPPA_NORT(
C    &           NOOEXC_A,IOOEXC_A,KAPPA_A,
C    &           NOOEXC_S,IOOEXC_S,KAPPA_S,CMOAO_IN,CMOAO_OUT)
      CALL NEWMO_FROM_KAPPA_NORT(
     &     NOOEXC_A,IOOEXC_A,KAPPA_A,NOOEXC_S,IOOEXC_S,KAPPA_S,
     &     WORK(KMOAOIN),WORK(KMOAOUT))
*
* Obtain the set of biorthonormal orbitals
*
      CALL GET_CBIO(WORK(KMOAOIN),WORK(KCBIO),WORK(KCBIO2))
*
* Biorthonormal integral transformaion
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Bioorthogonal integral transformation '
      END IF
*
      IE2LIST_A = IE2LIST_FULL_BIO
      IOCOBTP_A = 1
      INTSM_A = 1
      CALL PREPARE_2EI_LIST
*
      KKCMO_I = KMOAOUT
      KKCMO_J = KCBIO2
      KKCMO_K = KMOAOUT
      KKCMO_L = KCBIO2
*
C     DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
      CALL DO_ORBTRA(1,1,1,IE2LIST_FULL_BIO,IOCOBTP_A,INTSM_A)
      NINT1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1_F)
      CALL FLAG_ACT_INTLIST(IE2LIST_FULL_BIO)
*. The antisymmetric part of gradient
      CALL FOCK_MAT_NORT(WORK(KF),WORK(KF2),2,WORK(KFI),WORK(KFA))
*. And the interspace gradient
C     E1_FROM_F_NORT(E1,F1,F2,IOPSM,IOOEXC,IOOEXCC,
C    &           NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
      CALL E1_FROM_F_NORT(E1,WORK(KF),WORK(KF2),1,
     &     WORK(KLOOEXC),WORK(KLOOEXCC),NOOEXC_A,NTOOB,
     &     NTOOBS,NSMOB,IBSO,IREOST)
*. And add the active-active gradient
* The interspace excitations
C           VB_GRAD_ORBVBSPC(NOOEXC,IOOEXC,E1,C,VEC1_CSF,VEC2_CSF)
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' Active-active gradient will be calculated '
            CALL VB_GRAD_ORBVBSPC(NOOEXC_S,WORK(KLOOEXCC_S),
     &      WORK(KLE1+NOOEXC_IS),
     &      WORK(KL_VEC1),WORK(KL_VEC2),WORK(KL_VEC3))


C     VB_GRAD_ORBVBSPC(NOOEXC,IOOEXC,E1,C,
C    &           VEC1_CSF,VEC2_CSF)
*. Assuming just optimization in the VB space,
      CALL VB_GRAD_ORBVBSPC(NOOEXC_S,IOOEXC_S,E1,CCI,
     &     VEC1_CSF,VEC2_CSF)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) 
     & ' Orbital vector function from GET_VB_VF_VBSPC_FROM_KAPPA'
       WRITE(6,*) 
     & ' ======================================================='
       WRITE(6,*)
       WRITE(6,*) ' Part referring to antisymmetric operators: '
       CALL WRT_IOOEXCOP(E1,IOOEXC_S,NOOEXC_S)
       WRITE(6,*) ' Part referring to symmetric operators: '
       CALL WRT_IOOEXCOP(E1(1+NOOEXC_S),IOOEXC_S,NOOEXC_S)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTVBVF')
      RETURN
      END
      SUBROUTINE NEWMO_FROM_KAPPA_NORT(
     &           NOOEXC_A,IOOEXC_A,KAPPA_A,
     &           NOOEXC_S,IOOEXC_S,KAPPA_S,CMOAO_IN,CMOAO_OUT)
*
* Obtain New MO coefficients from symmetric and anti-symmetric
* kappa for VB calculation:
*
* CMOAO_OUT = CMOAO_IN * Exp(-Kappa_A S) Exp(-Kappa_S S)
*
* Jeppe Olsen, July 24, 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER IOOEXC_A(2,NOOEXC_A),IOOEXC_S(2,NOOEXC_S)
*. Antisymmetric and symmetric part of Kappa in packed form
      REAL*8
     &KAPPA_A(*),KAPPA_S(*)
      DIMENSION CMOAO_IN(*)
*. Output
      DIMENSION CMOAO_OUT(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from NEWMO_FROM_KAPPA_NORT'
        WRITE(6,*) ' ================================ '
        WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input KAPPA_A, KAPPA_S: '
        CALL WRTMAT(KAPPA_A,1,NOOEXC_A,1,NOOEXC_A)
        WRITE(6,*)
        CALL WRTMAT(KAPPA_S,1,NOOEXC_S,1,NOOEXC_S)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'NWMONO')
*
* Obtain Kappa_A and Kappa_S in full form
*
      NDIM_1F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLKAPPA_AE,NDIM_1F,'ADDL  ',2,'KAPPAE')
      CALL MEMMAN(KLKAPPA_SE,NDIM_1F,'ADDL  ',2,'KAPPSE')
C          REF_GN_KAPPA(KAPPAP,KAPPAE,IAS,ISM,IWAY,IOOEX,NOOEX)
      CALL REF_GN_KAPPA(KAPPA_A,WORK(KLKAPPA_AE),1,1,1,
     &     IOOEXC_A,NOOEXC_A)
      CALL REF_GN_KAPPA(KAPPA_S,WORK(KLKAPPA_SE),2,1,1,
     &     IOOEXC_S,NOOEXC_S)
*, Obtain metric in MO basis
      CALL MEMMAN(KLS,NDIM_1F,'ADDL  ',2,'SMOMO ')
      CALL GET_SMO(CMOAO_IN,WORK(KLS),0)
*
*. Obtain Exp (-Kappa_A S) Exp(-Kappa_S S)
*
      CALL MEMMAN(KLEXPMKS,NDIM_1F,'ADDL  ',2,'SMOMO ')
C     GET_EXPMKS(EXPMKS,KAPPA_S, KAPPA_A,S,NOBPS,NSMOB)
      CALL GET_EXPMKS(WORK(KLEXPMKS),WORK(KLKAPPA_SE),
     &     WORK(KLKAPPA_AE),WORK(KLS),
     &     NTOOBS,NSMOB)
*
* CMOAO_OUT = CMOAO_IN (-Kappa_A S) Exp(-Kappa_S S)
*
C  MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,LAROW,LACOL,LBROW,LBCOL,ITRNSP)
      CALL MULT_BLOC_MAT(CMOAO_OUT,CMOAO_IN,WORK(KLEXPMKS),NSMOB,
     &     NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) 
       WRITE(6,*) 
       WRITE(6,*) ' CMOAO_OUT: '
       CALL APRBLM2(CMOAO_OUT,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'NWMONO')
*
      RETURN
      END
      SUBROUTINE GENERIC_JAC_FROM_VF(JAC,NDIM,E1FUNC,E1,X,IDOSYM,
     &           ISTART,ISTOP)
*
*. Obtain Jacobian around X using external gradient function E1FUNC
*
* The Jacobian is assumed full, but is only calculated for
* the Columns ISTART to ISTOP
*
*. Jeppe Olsen, July 2011
*. Last modification; Jeppe Olsen; June 2013; ISTART, ISTOP added
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION X(NDIM)
*. Output
      REAL*8 JAC(NDIM,NDIM)
*. External
      EXTERNAL E1FUNC
*
* IORDER = 2: - second order formulae:
*
*  J Delta X = E1(X+Delta X) - E1(X-Delta)
*
* IORDER = 4: Fourth order formulae:
*
*  J Delta X = (-1/12) (       (E1(X+2Delta X) - E1(X-2Delta X))
*                       -8.0D0*(E1(X+Delta X)  - E1(X-Delta X) ))
*
*
*
      IORDER = 2
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from GENERIC_JAC_FROM_VF '
        WRITE(6,*) ' =============================='
        WRITE(6,*)
        WRITE(6,*) ' Order of method in use ', IORDER
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial set of parameters '
        CALL WRTMAT(X,1,NDIM,1,NDIM)
      END IF
*Evaluate vector function at point of expansion for check
C?    CALL E1FUNC(X,E1)
C?    WRITE(6,*) ' Vector function at initial point'
C?    CALL WRTMAT(E1,1,NDIM,1,NDIM)
C?    STOP ' After initial test '
*
      
*
*. Shift and constants for finite difference
      IF(IORDER.EQ.2) THEN
        DELTA = 1.0D-4
        FAC1 = 0.5D0/DELTA
      ELSE IF (IORDER.EQ.4) THEN
        DELTA = 1.0D-2
        FAC1 = 1.0D0/(12.0D0*DELTA)
        FAC2 = 8.0D0/(12.0D0*DELTA)
      END IF
*
      IF(IORDER.EQ.2) THEN
        DO J = ISTART, ISTOP
* E1(X+Delta X)
         X(J) = X(J) + DELTA
         CALL E1FUNC(X,E1)
         CALL COPVEC(E1,JAC(1,J),NDIM)
         CALL SCALVE(JAC(1,J),FAC1,NDIM)
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' E1(X+Delta)*FAC1: '
           CALL WRTMAT(JAC(1,J),1,NDIM,1,NDIM)
         END IF
* E1(X-Delta X)
         X(J) = X(J) - DELTA - DELTA
         CALL E1FUNC(X,E1)
         ONE = 1.0D0
         FAC1M = -FAC1
         CALL VECSUM(JAC(1,J),JAC(1,J),E1,ONE,FAC1M,NDIM)
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' (E1(X+Delta)-E1(X-Delta))*FAC1: '
           CALL WRTMAT(JAC(1,J),1,NDIM,1,NDIM)
         END IF
*. Clean up
         X(J) = X(J) + DELTA 
        END DO
      ELSE IF (IORDER.EQ.4) THEN
        DO J = 1, NDIM
* E1(X+2Delta X)
         X(J) = X(J) + 2.0D0*DELTA
         CALL E1FUNC(X,E1)
         CALL COPVEC(E1,JAC(1,J),NDIM)
         CALL SCALVE(JAC(1,J),-FAC1,NDIM)
* E1(X-2Delta X)
         X(J) = X(J) - 2.0D0*DELTA - 2.0D0*DELTA
         CALL E1FUNC(X,E1)
         ONE = 1.0D0
         CALL VECSUM(JAC(1,J),JAC(1,J),E1,ONE,FAC1,NDIM)
* E1(X+ Delta X)
         X(J) = X(J) + 2.0D0*DELTA + DELTA
         CALL E1FUNC(X,E1)
         ONE = 1.0D0
         CALL VECSUM(JAC(1,J),JAC(1,J),E1,ONE,FAC2,NDIM)
* E1(X- Delta X)
         X(J) = X(J) - DELTA - DELTA
         CALL E1FUNC(X,E1)
         ONE = 1.0D0
         CALL VECSUM(JAC(1,J),JAC(1,J),E1,ONE,-FAC2,NDIM)
*. Clean up
         X(J) = X(J) + DELTA 
        END DO
      END IF !Switch between the two procedures
*
      IF(IDOSYM.EQ.1) THEN
*. Symmetrize Jacobian
       DO I = 1, NDIM
        DO J = 1, I
         JAC(I,J) = 0.5D0*(JAC(I,J) + JAC(J,I))
         JAC(J,I) = JAC(I,J)
        END DO
       END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GENERIC_JAC_FROM_VF '
        WRITE(6,*) ' ================================='
        WRITE(6,*)
        CALL WRTMAT(JAC,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE GENERIC_GRAD_FROM_F(GRAD,NDIM,EFUNC,X)
*
*. Obtain gradient around  X using external function EFUNC
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION X(NDIM)
*. Output
      REAL*8 GRAD(NDIM)
*. External
      EXTERNAL EFUNC
*. The Gradient is obtained from finite difference using
*  Gradient Delta X = (-1/12) (       (E(X+2Delta X) - E(X-2Delta X))
*                       -8.0D0*(E(X+Delta X)  - E(X-Delta X) ))
*
*. Shift for finite difference
      DELTA = 1.0D-3
*
      DO J = 1, NDIM
* E(X+2Delta X)
       X(J) = X(J) + 2.0D0*DELTA
       EP2D = EFUNC(X)
* E1(X-2Delta X)
       X(J) = X(J) - 2.0D0*DELTA - 2.0D0*DELTA
       EM2D =  EFUNC(X)
* E(X+ Delta X)
       X(J) = X(J) + 2.0D0*DELTA + DELTA
       EP1D =  EFUNC(X)
* E1(X- Delta X)
       X(J) = X(J) - DELTA - DELTA
       EM1D = EFUNC(X)
*. And the synthesis
       GRAD(J) = -1.0D0/(12.0D0*DELTA)*(EP2D-EM2D) 
     &         +8.0D0/(12.0D0*DELTA)*(EP1D-EM1D)
*. Clean up
       X(J) = X(J) + DELTA 
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GENERIC_GRAD_FROM_F '
        WRITE(6,*) ' ================================='
        WRITE(6,*)
        CALL WRTMAT(GRAD,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
      FUNCTION E_VB_FROM_KAPPA_WRAP(KAPPA)
* 
* Wrapper routine for calculating Valence bond energy 
* from Kappa
* It is required that common /EVB_TRANS/ has been defined
*
*. Jeppe Olsen, July 25, 2011, on the train to Fjerritslev- cannot get the
*               code out of my head..
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      COMMON/EVB_TRANS/KLIOOEXC_A, KLKAPPA_A,
     &                 KLIOOEXC_S,KLKAPPA_S,
     &                 KL_C,KL_VEC2,KL_VEC3
*. Input
      REAL*8 KAPPA(*)
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from E_VB_FROM_KAPPA_WRAP'
        WRITE(6,*) ' ================================'
        WRITE(6,*) 
        WRITE(6,*) ' Kappa_A, Kappa_S '
        WRITE(6,*)
        CALL WRTMAT(KAPPA(1),NOOEXC_A,1,NOOEXC_A)
        WRITE(6,*)
        CALL WRTMAT(KAPPA(1+NOOEXC_A),1,NOOEXC_S,1,NOOEXC_S)
      END IF
*
      E_VB_FROM_KAPPA_WRAP = 
     &E_VB_FROM_KAPPA(KAPPA,NOOEXC_A,WORK(KLIOOEXC_A),
     &                KAPPA(1+NOOEXC_A),NOOEXC_S,WORK(KLIOOEXC_S),
     &                WORK(KL_C),WORK(KL_VEC2),WORK(KL_VEC3))
C     E_VB_FROM_KAPPA(
C    &           KAPPA_A,NOOEXC_A,IOOEXC_A,
C    &           KAPPA_S,NOOEXC_S,IOOEXC_S,CCI,
C    &           VEC1_CSF,VEC2_CSF)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Energy from E_VB_FROM_KAPPA_WRAP '
        WRITE(6,'(A,E15.8)') ' E = ', E_VB_FROM_KAPPA_WRAP
      END IF
*
      RETURN
      END
      FUNCTION E_VB_FROM_KAPPA(
     &         KAPPA_A,NOOEXC_A,IOOEXC_A,
     &         KAPPA_S,NOOEXC_S,IOOEXC_S,CCI,
     &         VEC1_CSF,VEC2_CSF)
*
* Obtain Valence bond energy from Kappa_A, Kappa_S
* Using method with expansion in complete VI space
*
*. It is assumed that the current MO-AO coefficients are in KMOAOIN.
* Integrals etc are overwritten, so the exit from this routine is 
* not clean.
*
*. Jeppe Olsen, July 24 2011 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'spinfo.inc'
*. Common block for communicating with sigma
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3, 
     &       LUCBIO_SAVE, LUHCBIO_SAVE,LUC_SAVE
*. Specific input
      INTEGER IOOEXC_A(2,NOOEXC_A), IOOEXC_S(2,NOOEXC_S)
      REAL*8 KAPPA_A(*), KAPPA_S(*)
      DIMENSION CCI(*)
      REAL*8 INPRDD
*. Scratch
      DIMENSION VEC1_CSF(*),VEC2_CSF(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' info from E_VB_FROM_KAPPA '
        WRITE(6,*) ' =========================='
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input Kappa_A and Kappa_S '
        CALL WRTMAT(KAPPA_A,1,NOOEXC_A,1,NOOEXC_A)
        WRITE(6,*)
        CALL WRTMAT(KAPPA_S,1,NOOEXC_S,1,NOOEXC_S)
       END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EVBFKA')
*
*. Obtain New MO coefficients in MOAOUT: MOAOIN* Exp(-Kappa_A S) Exp(-Kappa_S S)
*
C     NEWMO_FROM_KAPPA_NORT(
C    &           NOOEXC_A,IOOEXC_A,KAPPA_A,
C    &           NOOEXC_S,IOOEXC_S,KAPPA_S,CMOAO_IN,CMOAO_OUT)
C?    WRITE(6,*) ' NOOEXC_A, NOOEXC_S before call to NEWMO' ,
C?   &             NOOEXC_A, NOOEXC_S
      CALL NEWMO_FROM_KAPPA_NORT(
     &     NOOEXC_A,IOOEXC_A,KAPPA_A,NOOEXC_S,IOOEXC_S,KAPPA_S,
     &     WORK(KMOAOIN),WORK(KMOAOUT))
*
* Obtain the set of biorthonormal orbitals
*
      CALL GET_CBIO(WORK(KMOAOUT),WORK(KCBIO),WORK(KCBIO2))
*
* Biorthonormal integral transformaion
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Bioorthogonal integral transformation '
      END IF
*
      IE2LIST_A = IE2LIST_FULL_BIO
      IOCOBTP_A = 1
      INTSM_A = 1
      CALL PREPARE_2EI_LIST
*
      KKCMO_I = KMOAOUT
      KKCMO_J = KCBIO2
      KKCMO_K = KMOAOUT
      KKCMO_L = KCBIO2
*
C          DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
      CALL DO_ORBTRA(1,1,0,IE2LIST_FULL_BIO,IOCOBTP_A,INTSM_A)
      CALL FLAG_ACT_INTLIST(IE2LIST_FULL_BIO)
      NINT1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1_F)
*
*. Sigma with the current C 
*
C          SIGMA_NORTCI(C,HC,SC,IDOHC,IDOSC)
      CALL SIGMA_NORTCI(CCI,VEC1_CSF,VEC2_CSF,1,1)
      IF(NTEST.GE.1000) WRITE(6,*) ' Back from SIGMA_NORTCI'
* calculate energy from vectors on file
      CHC = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE,LUHCBIO_SAVE,1,-1)
      CC =  INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE,LUCBIO_SAVE,1,-1)
      EVB = CHC/CC
*
      E_VB_FROM_KAPPA = EVB
*
      WRITE(6,'(A,3(2X,E14.8))') ' Energy: CHC, CC, CHC/CC ', 
     &                        CHC,CC,EVB
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EVBFKA')
      RETURN
      END
      SUBROUTINE ORBHES_VB(E2,IFORM)
*
*Obtain complete or part of Orbital Hessian for VB approach
*
* IFORM = 1 => Complete orbital Hessian
*
*. Jeppe Olsen, July 26, 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'wrkspc-static.inc'
*. Output: Complete orbital Hessian in lower packed form 
      DIMENSION E2(*)
      EXTERNAL VB_BR_FOR_KAPPA_WRAP
*
* Method for calculating orbital Hesssian
      I_ORBHES_MET = 2
* IORBHES_MET = 1 => Finite difference based on energy
* IORBHES_MET = 2 => Finite difference bases on Vector function
* IORBHES_MET = 2 => Analytic calc of antisym, FD calc of symmetric part
      NTEST = 000
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Info from ORBHES_FD'
       WRITE(6,*) ' ================== '
       WRITE(6,*) 
       IF(I_ORBHES_MET.EQ.1) THEN
         WRITE(6,*)
     &   ' Orbital Hessian obtained from energy finite difference'
       ELSE IF( I_ORBHES_MET.EQ.2) THEN
         WRITE(6,*)
     &   ' Orbital Hessian obtained from gradient finite difference'
       END IF
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'OBE2VB')
*  
      KKCMO_I = KMOAOUT
      KKCMO_J = KCBIO2
      KKCMO_K = KMOAOUT
      KKCMO_L = KCBIO2
*
      NOOEXC_TOT = NOOEXC_A + NOOEXC_S
*
      IF(I_ORBHES_MET .EQ. 1) THEN
        CALL ORBHES_VB_FD(E2)
      ELSE
* A local copy of complete Hessian, BR-vector and kappa
       CALL MEMMAN(KLE2,NOOEXC_TOT**2,'ADDL  ',2,'E2FULL')
       CALL MEMMAN(KLBR,NOOEXC_TOT,'ADDL  ',2,'BRVEC ')
       CALL MEMMAN(KLKAP,NOOEXC_TOT,'ADDL  ',2,'KLKAP ')
*
* We will evaluate Hessian at current expansion point, so
       ZERO = 0.0D0
       CALL SETVEC(WORK(KLKAP),ZERO,NOOEXC_TOT)
*. FUSK
       IREADJ = 0
       IF(IREADJ.EQ.1) THEN
*. Jacobian is read in rather than constructed '
        WRITE(6,*) ' WARNING: JACO READ IN FROM LU95 '
        LU95 = 95
        CALL REWINO(LU95)
        NELMNT = NOOEXC_TOT*(NOOEXC_TOT+1)/2
        READ(LU95,*) (E2(IJ), IJ = 1, NELMNT)
       ELSE
*
       CALL GENERIC_JAC_FROM_VF(WORK(KLE2),NOOEXC_TOT,
     &      VB_BR_FOR_KAPPA_WRAP, WORK(KLBR),WORK(KLKAP),1,
     &      1, NOOEXC_TOT)
C      GENERIC_JAC_FROM_VF(JAC,NDIM,E1FUNC,E1,X,IDOSYM)
C      TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
       CALL TRIPAK(WORK(KLE2),E2,1,NOOEXC_TOT,NOOEXC_TOT)
      END IF ! JACO read in
      END IF ! I_ORBHES_MET = 1
*
      IDUMPJ = 1
      IF(IDUMPJ.EQ.1) THEN
        WRITE(6,*) ' Jacobian is dumped to file 95 '
        LU95 = 95
        CALL REWINO(LU95)
        NELMNT = NOOEXC_TOT*(NOOEXC_TOT+1)/2
        WRITE(LU95,*) (E2(IJ), IJ = 1, NELMNT)
      END IF
       
*
      IF(I_ORBHES_MET.LE.2) THEN
*. Restore order- and integrals
        IE2LIST_A = IE2LIST_FULL_BIO
        IOCOBTP_A = 1
        INTSM_A = 1
        CALL PREPARE_2EI_LIST
        CALL GET_CBIO(WORK(KMOAOIN),WORK(KCBIO),WORK(KCBIO2))
*
        KKCMO_I = KMOAOIN
        KKCMO_J = KCBIO2
        KKCMO_K = KMOAOIN
        KKCMO_L = KCBIO2
*
C          DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
        CALL DO_ORBTRA(1,1,0,IE2LIST_FULL_BIO,IOCOBTP_A,INTSM_A)
        NINT1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
        CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1_F)
        CALL FLAG_ACT_INTLIST(IE2LIST_FULL_BIO)
      END IF
*
      IF(NTEST.GE.1000) THEN
        NOOEXC_TOT = NOOEXC_A + NOOEXC_S
        WRITE(6,*) ' Orbital Hessian '
        CALL APRBLM2(E2,NOOEXC_TOT,NOOEXC_TOT,1,1)
      END IF
* 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'OBE2VB')
      RETURN
      END
      SUBROUTINE ORBHES_VB_FD(E2)
*
*. Obtain Orbital Hessian for VB by energy Finite difference
*
*. Jeppe Olsen, July 28, 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'crun.inc'
* EVB_TRANS must have been set outside
      COMMON/EVB_TRANS/KLIOOEXC_A, KLKAPPA_A,
     &                 KLIOOEXC_S,KLKAPPA_S,
     &                 KL_C,KL_VEC2,KL_VEC3
      EXTERNAL E_VB_FROM_KAPPA_WRAP
*
*. Output: Hessian in lower packed form
*
      DIMENSION E2(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'OBE2FD')
*. Copy of Hessian in complete form
      NOOEXC_T = NOOEXC_A + NOOEXC_S
      CALL MEMMAN(KLE2F,NOOEXC_T**2,'ADDL  ',2,'E2F   ')
*
      CALL MEMMAN(KLE1,NOOEXC_T,'ADDL  ',2,'KLE1')
C          GENERIC_GRA_HES_FD(E0,E1,E2,X,NX,EFUNC)
      KLKAPPA = KLKAPPA_A
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLKAPPA),ZERO,NOOEXC_T)
      CALL GENERIC_GRA_HES_FD(E0,WORK(KLE1),WORK(KLE2F),WORK(KLKAPPA),
     &     NOOEXC_T,E_VB_FROM_KAPPA_WRAP)
*. Pack to lower half
            CALL TRIPAK(WORK(KLE2F),E2,1,NOOEXC_T,NOOEXC_T)
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ORBHES_VB_FD '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        WRITE(6,'(A,E15.8)') ' Current energy = ', E0
        WRITE(6,'(A)') ' Gradient: '
        CALL WRTMAT(WORK(KLE1),1,NOOEXC_T,1,NOOEXC_T)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,'(A)') ' Hessian: '
        CALL PRSYM(E2,NOOEXC_T)
C?      CALL WRTMAT(E2,NOOEXC_T,NOOEXC_T,NOOEXC_T,NOOEXC_T)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'OBE2FD')
      RETURN
      END
      SUBROUTINE GET_CMOINI_GEN(CINIAO_UT,CINIUT_INIIN,CINIAO_IN)
*
* Obtain starting set of orbitals.
* May be obtained from fragment orbitals
*
*. Output:
*     Expansion of starting orbitals in AO: CINIAO_UT
*     Expansion of starting orbitals in initial orbitals: CINIUT_INIIN
*  Input:
*     Expansion of initial initial orbitals: CINIAO_IN

*
*. Jeppe Olsen, April 2012, extended June 2012
*               March 2013, added a bit for supersymmetry
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'fragmol.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Input
       DIMENSION CINIAO_IN(*)
*. Output
       DIMENSION CINIAO_UT(*), CINIUT_INIIN(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUN,'MOING')
      NTEST = 10
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Wellcome to GET_CMOINI_GEN'
       WRITE(6,*) ' =========================='
      END IF
*
      IF(NTEST.GE.2) THEN
C       WRITE(6,*) ' INI_MO_TP, INI_MO_ORT = ', INI_MO_TP, INI_MO_ORT
        WRITE(6,*)
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Initial set of orbitals '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
*
        IF(INI_MO_TP.EQ.1) THEN
          WRITE(6,'(4X,A)') ' Atomic orbitals will be used '
        ELSE IF (INI_MO_TP.EQ.2) THEN
          WRITE(6,'(4X,A)')
     &    ' Input MOs in VB space rotated  to give diagonal block'
        ELSE IF (INI_MO_TP.EQ.3) THEN
          WRITE(6,'(4X,A)')
     &    ' Initial MO orbitals from SIRIFC/91 will be used'
        ELSE IF (INI_MO_TP.EQ.4) THEN
          WRITE(6,'(4X,A)')
     &    ' Constructed from fragment orbitals'
        ELSE IF (INI_MO_TP.EQ.5) THEN
          WRITE(6,'(4X,A)')
     &    ' Initial MO orbitals from LUCINF_O will be used'
        END IF
*
        IF(INI_MO_TP.NE.3) THEN
         WRITE(6,'(4X,A)')
     &   ' Orbitals in inactive and secondary space will be ort.'
         WRITE(6,'(4X,A)') ' Orbitals in GAS orbital spaces(.ne. VB ): '
         IF(INI_MO_ORT.EQ.0) THEN
           WRITE(6,'(6X,A)') ' No orthogonalization  '
         ELSE IF (INI_MO_ORT.EQ.1) THEN
           WRITE(6,'(6X,A)') ' Orthogonalized'
         END IF
         WRITE(6,'(4X,A)') ' Orbitals in VB orbital space: '
         IF(INI_ORT_VBGAS.EQ.0) THEN
           WRITE(6,'(6X,A)') ' No orthogonalization  '
         ELSE IF (INI_ORT_VBGAS.EQ.1) THEN
           WRITE(6,'(6X,A)') ' Orthogonalized'
         END IF
        END IF
*
*. In general, the output form of the orbitals are unknown
*
        CMO_ORD = 'UNK'
*
        IF(INI_MO_TP.EQ.4) THEN
         WRITE(6,*) ' Distribution of orbitals from fragments:'
         DO IFRAG = 1, NFRAG_MOL
          NSMOB_L = NSMOB_FRAG(IFRAG)
          WRITE(6,'(A,I3)') ' For fragment ', IFRAG
          WRITE(6,*)        ' ===================='
          WRITE(6,*) ' Number of orbitals per GAS (row) and sym (col) '
          CALL IWRTMA
     &    (N_GS_SM_BAS_FRAG(0,1,IFRAG),NGAS+2,NSMOB_L,MXPNGAS+1,MXPOBS)
         END DO
        END IF ! End if INI_MO_TP.eq.4
      END IF !NTEST test
*
* Two steps : 0) Orthogonalize to frozen orbitals
*             1) Obtain a set of (nonorthogonal) initial orbitals
*             2) Perform (partial) orthonormalization to obtain
*                Final initial orbitals
*
* Generate set of (nonorthogonal) initial orbitals
*
       CALL GET_INIMO(CINIAO_UT)
C           GET_INIMO(CMO_INI)
*
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Expansion of initial MOs in AOs '
        WRITE(6,*) ' ================================'
        CALL APRBLM_F7(CINIAO_UT,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*. MO_TP = 3 => we are done...
      IF(INI_MO_TP.EQ.3) GOTO 9999
*
*. Orthogonalize to frozen orbitals
*. Jeppe, I am not sure if this is working in connection with supersymmetry reordering...
*. (What are the numbers defining the localized orbitals?)
      IF(NFRZ_ORB.NE.0) THEN
        CALL ORT_CMO_TO_FROZEN_ORBITALS(CINIAO_UT)
        IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Orbitals orthogonalized to frozen '
         CALL APRBLM_F7(CINIAO_UT,NTOOBS,NTOOBS,NSMOB,0)
        END IF
      END IF
*
      CMO_ORD = 'UNK'
*
* New initial orbitals in terms of initial initial orbitals(KMOAOIN)
*
* CINIUT_INIIN = CINIAO_UT* CINIAO_IN**-1
*
*. Invert CINIAO_IN 
      LMOMO = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLCMOS,2*LMOMO,'ADDL ',2,'CMOS  ')
      CALL MEMMAN(KLCMOI,  LMOMO,'ADDL ',2,'CMOI  ')
      IPROBLEM = 0
      CALL INV_BLKMT(CINIAO_IN,WORK(KLCMOI),WORK(KLCMOS),NSMOB,
     &               NTOOBS,IPROBLEM)
C          INV_BLKMT(A,AINV,SCR,NBLK,LBLK,IPROBLEM)
      IF(IPROBLEM.NE.0) THEN
        WRITE(6,*) ' Problem inverting CMOAOUT '
        STOP       ' Problem inverting CMOAOUT '
      END IF
*. And multiply
C  MULT_H1H2(H1,IH1SM,H2,IH2SM,H12,IH12SM)
      CALL MULT_H1H2(WORK(KLCMOI),1,CINIAO_UT,1,CINIUT_INIIN,IUTSM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Expansion of initial MOs in Initial initial MOs '
        WRITE(6,*) ' ====================================='
        CALL APRBLM_F7(CINIUT_INIIN,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
* Check of orthogonality of reexpansion of initial orbitals
*
C    MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
C    &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
      CALL MULT_BLOC_MAT(WORK(KLCMOS),CINIUT_INIIN,CINIUT_INIIN,NSMOB,
     &     NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,1)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CINIUT_INIIN*CINIUT_INIIN(T) '
        WRITE(6,*) ' ============================='
        CALL APRBLM2(WORK(KLCMOS),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
 9999 CONTINUE
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUN,'MOING')
*
      RETURN
      END 
      SUBROUTINE EXTR_SYMGAS_BLK_FROM_ORBMAT
     &           (A,ABLK,ISM,IGAS,JSM,JGAS)
*
* A symmetryblocked (not lower half packed) matrix A over orbitals is given
* Extract block referring to GASpaCE IGAS, JGAS and symmetry ISM,JSM
*
* I_EX_OR_CP = 1 => Extract from A to IGAS
* I_EX_OR_CP = 1 => Copy from IGAS to A
*
*. Jeppe Olsen, May 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
*. Specific input and output
      DIMENSION A(*), ABLK(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' EXTR_SYMGAS_BLK_FROM_ORBMAT '
        WRITE(6,*) ' =========================== '
      END IF
        
*. Symmetry of matrix
      IJSM = MULTD2H(ISM,JSM)
*. Offsets to symmetry block in full matrix matrix
      IOFF_IN = 1
      DO IISM = 1, ISM-1
        JJSM = MULTD2H(IISM,IJSM)
        IOFF_IN = IOFF_IN + NTOOBS(IISM)*NTOOBS(JJSM)
      END DO
*. Offset to start of orbitals in given gas
      IOFF = 1
      DO IIGAS = 0, IGAS -1
        IOFF = IOFF + NOBPTS_GN(IIGAS,ISM)
      END DO
*
      JOFF = 1
      DO JJGAS = 0, JGAS -1
        JOFF = JOFF + NOBPTS_GN(JJGAS,JSM)
      END DO
*
      NI = NOBPTS_GN(IGAS,ISM)
      NJ = NOBPTS_GN(JGAS,JSM)
      NIS = NTOOBS(ISM)
      NJS = NTOOBS(JSM)
      DO J = 1, NJ
        DO I = 1, NI
          IJ_OUT = (J-1)*NI + I 
          IJ_IN  = IOFF_IN -1 + (JOFF+J-1-1)*NIS + IOFF+I-1
          ABLK(IJ_OUT) = A(IJ_IN)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Submatrix with ISM, JSM, IGAS, JGAS = ',
     &   ISM, JSM, IGAS, JGAS
         CALL WRTMAT(ABLK,NI,NJ,NI,NJ)
      END IF
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Full matrix '
         CALL APRBLM2(A,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE VB_DENSI(RHO1,RHO2,IR12,C,VEC1_CSF,VEC2_CSF)
*
*
* Obtain one-body density matrix over active space for VB function
*
*
* E(IJ) =  <0!(E(ij))!0> /<0!0>
*
* and if IR12 = 2 also the two-body density matrix in mixed basis
* E(IJ,KL) = <0!\tilde a+i \sigma \tilde a+k sigma' a l sigma' a j sigma!0>
*
* Note that whereas the one-eletron density is transformed to the 
* actual MO-basis, the two-body density is kept in the mixed basis
*
* So to obtain gradient 
* 1: construct bioorthogonal expansion of  !0>
* 2: Set up density matrices <0!E(ij)!0> 
*    where i is in biobase and j in normal
* 3: Transform density matrices to standard basis
*
* The current CI coefficients in the CSF basis are in C, where 
* VEC1_CSF, VEC2_CSF, must be able to hold these expansions
*
* This is an initial version, for initial calculations and checks
*
* Jeppe Olsen, May 2012, for the initial NORTMCSCF program
*
* Sitting in Palermo, preparing for a talk ...
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3, 
     &       LUCBIO_SAVE, LUHCBIO_SAVE,LUC_SAVE
      REAL*8 INPRDD
*. Input
      DIMENSION C(*)
*. Scratch 
      DIMENSION VEC1_CSF(*), VEC2_CSF(*)
*. Output
      DIMENSION RHO1(*), RHO2(*)
*
      NTEST = 10
*. CSFs are handled explicitly, so
      NOCSF = 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========'
        WRITE(6,*) ' VB_DENSI'
        WRITE(6,*) ' ========'
        WRITE(6,*)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',2,'VBDENS')
*
      LUSCR1 = LUSC34
      LUSCR2 = LUSC35
      LUSCR3 = LUSC36
      LUCBIO_SAVE = 110
      LUC_SAVE = 112
*
* A bit of scratch
*
      LEN_1A = NDIM_1EL_MAT(1,NACOBS,NACOBS,NSMOB,0)
      CALL MEMMAN(KLRHOB,NACOB**2,'ADDL  ',2,'RHOB  ')
      CALL MEMMAN(KLSCR ,NACOB**2,'ADDL  ',2,'SCR   ')
      CALL MEMMAN(KLCBIOA,LEN_1A,'ADDL  ',2,'CBIOAC')
*. Preparation: Obtain CBIO over active orbitals only
C          EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT(A,AGAS,I_EX_OR_CP)
      CALL EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT
     &     (WORK(KCBIO),WORK(KLCBIOA),1)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' CBIO in active orbitals '
        CALL APRBLM2(WORK(KLCBIOA),NACOBS,NACOBS,NSMOB,0)
      END IF
*
*. Biotransform C 
*
C          SIGMA_NORTCI(C,HC,SC,IDOHC,IDOSC)
      CALL SIGMA_NORTCI(C,VEC1_CSF,VEC2_CSF,0,1)
      IF(NTEST.GE.1000) WRITE(6,*) ' Back from SIGMA_NORTCI'
* calculate Overlap from vectors on file - for check
      CC  = INPRDD(WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE, LUCBIO_SAVE,1,-1)
      IF(NTEST.GE.100) WRITE(6,*) ' <0!0> =', CC
*
*. Set up density <0! a+i(bio) aj!0(bio)> in RHOB
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' C in SD expansion '
        CALL WRTVCD(WORK(KVEC1P),LUC_SAVE,1,-1)
        WRITE(6,*) ' C(bio) in SD expansion '
        CALL WRTVCD(WORK(KVEC1P),LUCBIO_SAVE,1,-1)
      END IF
      XDUM = 0.0D0
      CALL DENSI2(IR12 ,WORK(KLRHOB),RHO2,
     &WORK(KVEC1P),WORK(KVEC2P),LUC_SAVE,LUCBIO_SAVE,EXPS2,
     &0,XDUM,XDUM,XDUM,XDUM,0)
*
      FACTOR = 1.0D0/CC
C?    WRITE(6,*) ' CC = ', CC 
      CALL SCALVE(WORK(KLRHOB),FACTOR,NACOB**2)
      IF(IR12.EQ.2) THEN
        LRHO2 = NACOB**2*(NACOB**2+1)/2
        CALL SCALVE(RHO2,FACTOR,LRHO2)
      END IF
        
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Density matrix <0! a+i(bio) aj!bio 0>/<0!0> '
       CALL WRTMAT(WORK(KLRHOB),NACOB,NACOB,NACOB,NACOB)
      END IF
*. Obtain density as blocked matrix over symmetry blocks of active orbitals
C          REORHO1(RHO1I,RHO1O,IRHO1SM)
      CALL REORHO1(WORK(KLRHOB),WORK(KLSCR),1,1)
      CALL COPVEC(WORK(KLSCR),WORK(KLRHOB),LEN_1A)
*. Transform the densities from bio, normal to the normal, normal basis
C     TR_BIOMAT(XIN,XOUT,CBIO,NORB_PSM, 
C    &            INB_IN,INB_OUT,JNB_IN,JNB_OUT)
      CALL TR_BIOMAT(WORK(KLRHOB),WORK(KLSCR),WORK(KLCBIOA),
     &     NACOBS,2,1,1,1)
*. Transfer back to full matrix over active orbitals
      CALL REORHO1(RHO1,WORK(KLSCR),1,2)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Density matrix <0! E(ij) !> '
       CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',2,'VBDENS')
      RETURN
      END
      SUBROUTINE GET_SACT(SACT,C)
*
*. Obtain the overlap matrix of the active orbitals for a MO-AO expansion 
* given by the MO-AO expansion matric C
*
*. Jeppe Olsen, May 31 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
*. Specific input
      DIMENSION C(*)
*. Specific output: in symmetry-packed lower half form
      DIMENSION SACT(*)
*. Obtain expansion of active orbitals only
    
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from SACT '
        WRITE(6,*) ' ============== '
      END IF
*
*. It is assumed that SAO resides in WORK(KSAO)
      IDUM = 0 
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'VB_SAC')
*
*. Two ways of which IWAY = 1 has a bug...
*
      LEN_CACT = LEN_BLMAT(NSMOB,NACOBS,NTOOBS,0)
      LEN_C = LEN_BLMAT(NSMOB,NTOOBS,NTOOBS,0)
*
      IWAY = 2
      IF(IWAY.EQ.1) THEN
      CALL MEMMAN(KLCACT,LEN_CACT,'ADDL  ',2,'C_AC  ')
      CALL MEMMAN(KLSCR ,2*LEN_C,'ADDL  ',2,'SCR   ')
*. Obtain C over active orbitals only
C          EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT(A,AGAS,I_EX_OR_CP)
        CALL EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT
     &       (C,WORK(KLCACT),1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' C over active orbitals '
          CALL APRBLM2(WORK(KLCACT),NACOBS,NACOBS,NSMOB,0)
        END IF
* 
        CALL TRAN_SYM_BLOC_MAT4(WORK(KSAO),WORK(KLCACT),WORK(KLCACT),
     &       NSMOB,NTOOBS,NACOBS,SACT,WORK(KLSCR),1)
C            TRAN_SYM_BLOC_MAT4
C    &  (AIN,XL,XR,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      ELSE
*. Obtain full SMO and extract active blocks
       CALL MEMMAN(KLS1,LEN_C,'ADDL  ',2,'S_FULL')
       CALL MEMMAN(KLS2,LEN_C,'ADDL  ',2,'S2FULL')
       CALL GET_SMO(WORK(KMOAOUT),WORK(KLS1),0)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Full S matrix '
         CALL APRBLM2(WORK(KLS1),NTOOBS,NTOOBS,NSMOB,0)
       END IF
*. Extract active blocks
       CALL EXTR_OR_CP_ACT_BLKS_FROM_ORBMAT
     &       (WORK(KLS1),WORK(KLS2),1)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' S matrix over activt orbitals'
         CALL APRBLM2(WORK(KLS2),NACOBS,NACOBS,NSMOB,0)
       END IF
*. And pack these
C  TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
       CALL TRIPAK_BLKM(WORK(KLS2),SACT,1,NACOBS,NSMOB)
      END IF !switch between routes
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Overlap matrix over active orbitals '
       CALL APRBLM2(SACT,NACOBS,NACOBS,NSMOB,1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'VB_SAC')
*
      RETURN
      END
      SUBROUTINE NONORT_NATORB(SACT,RHO1)

* Obtain natural orbitals for a density matrix in a 
* nonorthogonal basis
*
*. Jeppe Olsen, May 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Specific input: SACT in symmetry-blocked lower half packed form 
*. and RHO1 over all active orbitals in standard type-symmetry order
*
      DIMENSION SACT(*),RHO1(NACOB,NACOB)
*
      NTEST = 10
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from NONORT_NATORB '
        WRITE(6,*) ' ========================'
        WRITE(6,*)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'NORNAT')
*
*. Some scratch space
*  ==================
*. Density matrix in symmetry-packed complete form
      LEN_CACT = LEN_BLMAT(NSMOB,NACOBS,NACOBS,0)
      CALL MEMMAN(KLRH_SYM,LEN_CACT,'ADDL  ',2,'RH_SYM')
*. Unpacked overlap matrix
      CALL MEMMAN(KLSUNP,LEN_CACT,'ADDL  ',2,'S_UNP ')
*. Expansion coefficient of natural orbitals
      CALL MEMMAN(KLCNAT,LEN_CACT,'ADDL  ',2,'C_NAT ')
*. Matrix for going to orthonormal basis
      CALL MEMMAN(KLP,LEN_CACT,'ADDL  ',2,'P_TRA ')
*. Natural occupation numbers
      CALL MEMMAN(KLOCC,NACOB,'ADDL  ',2,'P_TRA ')

*. Obtain density in blocks of symmetry
*. Loop over active orbitals in output order: symmetry type
      IOBOFF = 0
      IMTOFF = 0
      IADD_ST = 0
      IADD_TS = NINOB
      DO ISMOB = 1, NSMOB
        IF(ISMOB.EQ.1) THEN
          IOBOFF     = 1
          IMTOFF     = 1
          IADD_ST    = NINOBS(1)
        ELSE
          IOBOFF     = IOBOFF + NACOBS(ISMOB-1)
          IMTOFF     = IMTOFF + NACOBS(ISMOB-1)**2
          IADD_ST    = IADD_ST + NINOBS(ISMOB) + NSCOBS(ISMOB-1)
        END IF
        LOB = NACOBS(ISMOB)
C?      WRITE(6,*) ' ISMOB, LOB, = ', ISMOB, LOB
C?      WRITE(6,*) ' IADD_TS = ', IADD_TS
*
*. Extract symmetry block of density matrix
*
*. Loop over active orbitals of symmetry ISMOB in ST order
        DO IOB = IOBOFF,IOBOFF + LOB-1
           IOB_ABS = IOB + IADD_ST
C          IOB_TS = ISTREO(IOB_ABS) - IADD_TS
           IOB_TS = IREOST(IOB_ABS) - IADD_TS
           IOB_REL = IOB  - IOBOFF + 1
           DO JOB = IOBOFF,IOBOFF + LOB-1
               JOB_ABS = JOB + IADD_ST
               JOB_TS = IREOST(JOB_ABS) - IADD_TS
               JOB_REL = JOB  - IOBOFF + 1
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' JOB, JOB_ABS, JOB_TS, IREOST() = ',
     &                        JOB, JOB_ABS, JOB_TS, IREOST(JOB_ABS)
                 WRITE(6,*) ' IOB_TS, JOB_TS = ', IOB_TS, JOB_TS
                 WRITE(6,'(A,6I3)')
     &           ' IOB_TS, JOB_TS, IOB, JOB, IOB_REL, JOB_REL  = ',
     &             IOB_TS, JOB_TS, IOB, JOB, IOB_REL, JOB_REL
               END IF
               WORK(KLRH_SYM-1+IMTOFF-1+(JOB_REL-1)*LOB+IOB_REL)
     &       = RHO1(IOB_TS,JOB_TS)
           END DO !Job
        END DO ! Iob
      END DO! Loop over symmetries of orbitals
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' One-body density matrix in symmetry-blocks ' 
        CALL APRBLM2(WORK(KLRH_SYM),NACOBS,NACOBS,NSMOB,0)
      END IF
*. Unpack overlapmatrix
C TRIPAK_BLKM(AUTPAK,APAK,IWAY,LBLOCK,NBLOCK)
      CALL TRIPAK_BLKM(WORK(KLSUNP),SACT,2,NACOBS,NSMOB)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Overlap matrix in unpacked form '
        CALL APRBLM2(WORK(KLSUNP),NACOBS,NACOBS,NSMOB,0)
      END IF
*. Multiply density with -1 to get highest occupation numbers first
      ONEM = -1.0D0
      CALL SCALVE(WORK(KLRH_SYM),ONEM,LEN_CACT)
*. Diagonalize
C     GENDIA_BLMAT(HIN,SIN,C,E,PVEC,NBLK,LBLK,ISORT)
      CALL GENDIA_BLMAT(WORK(KLRH_SYM),WORK(KLSUNP),WORK(KLCNAT),
     &     WORK(KLOCC),WORK(KLP),NACOBS,NSMOB,1)
*. Multiply occupation numbers with -1 to counteract previous multiply
      CALL SCALVE(WORK(KLOCC),ONEM,NACOB)
*
      WRITE(6,*) ' Natural occupation numbers: '
      WRITE(6,*) ' =========================== '
      WRITE(6,*)
*
      DO ISYM = 1, NSMOB
       IF(ISYM.EQ.1) THEN
         IOFF_I = 1
         IOFF_IJ = 1
       ELSE
         IOFF_I = IOFF_I + NACOBS(ISYM-1)
         IOFF_IJ = IOFF_IJ + NACOBS(ISYM-1)**2
       END IF
       WRITE(6,*)
       WRITE(6,*)
     & ' Natural occupation numbers for symmetry = ', ISYM
       WRITE(6,*)
     & ' ==================================================='
       L = NACOBS(ISYM)
       CALL WRTMAT(WORK(KLOCC-1+IOFF_I),1,L,1,L)
       WRITE(6,*)
     & ' Expansion of natural orbitals for symmetry = ', ISYM
       WRITE(6,*)
     & ' ==================================================='
       CALL WRTMAT(WORK(KLCNAT-1+IOFF_IJ),L,L,L,L)
      END DO! Loop over symmetries
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'NORNAT')
      RETURN
      END
      SUBROUTINE VB_BR_FOR_KAPPA_WRAP(KAPPA,BR)
*
* Outer routine for obtaining generalized Brillouin vector
* at a given point
*
*. Jeppe Olsen, May 31, 2012 in Palermo, (18 hours to talk)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'crun.inc'
*
      COMMON/EVB_TRANS/KLIOOEXC_A, KLKAPPA_A,
     &                 KLIOOEXC_S,KLKAPPA_S,
     &                 KL_C,KL_VEC2,KL_VEC3,
     &                 KLOOEXC
*
*
*. Input
      REAL*8 KAPPA(*)
*. And output
      DIMENSION BR(*)
      
      NTEST = 01
      IF(NTEST.GE.1) WRITE(6,*) ' Entering VB_BR_FOR_KAPPA_WRAP'
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from VB_BR_FOR_KAPPA_WRAP'
        WRITE(6,*) ' =============================='
        WRITE(6,*)
        WRITE(6,*) ' Kappa_A, Kappa_S '
        WRITE(6,*)
        WRITE(6,*) ' NOOEXC_A, NOOEXC_S = ', 
     &               NOOEXC_A, NOOEXC_S
        CALL WRTMAT(KAPPA(1),NOOEXC_A,1,NOOEXC_A)
        WRITE(6,*)
        CALL WRTMAT(KAPPA(1+NOOEXC_A),1,NOOEXC_S,1,NOOEXC_S)
      END IF
*. And call the routine that does the job
      CALL VB_BR_FROM_KAPPA(BR,
     &     NOOEXC_A,WORK(KLIOOEXC_A),KAPPA(1),
     &     NOOEXC_S,WORK(KLIOOEXC_S),KAPPA(1+NOOEXC_A),
     &     WORK(KLOOEXC),
     &     WORK(KL_C),WORK(KL_VEC2),WORK(KL_VEC3))
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Brillouin vector from VB_BR_FOR_KAPPA_WRAP'
        WRITE(6,*) ' ========================================= '
        WRITE(6,*)
        N = NOOEXC_A + NOOEXC_S
        CALL WRTMAT(BR,1,N,1,N)
      END IF
*
      RETURN
      END 
      SUBROUTINE VB_BR_FROM_KAPPA(BR,
     &           NOOEXC_A,IOOEXC_A, KAPPA_A,
     &           NOOEXC_S,IOOEXC_S, KAPPA_S,
     &           IOOEXC,
     &           C,VEC2,VEC3)
*
* Obtain VB Brillouin vector for a given set of Kappa parameters
*
*. Jeppe Olsen, May 31, Palermo  - Finished June 3, Zurich
*
*.It is assumed that the current MO-AO coefficients are in KMOAOIN.
* Integrals etc are overwritten, so the exit from this routine is
* not clean.
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'spinfo.inc'
*. Explicit input
      REAL*8 KAPPA_A(NOOEXC_A),KAPPA_S(NOOEXC_S)
      INTEGER IOOEXC_A(2,NOOEXC_A), IOOEXC_S(2,NOOEXC_S), IOOEXC(*)
*. Coefficients
      DIMENSION C(*)
*. Output
      DIMENSION BR(*)
*. Scratch vectors
      DIMENSION VEC2(*),VEC3(*)
*
*. Common block for communicating with sigma
      COMMON/SCRFILES_MATVEC/LUSCR1,LUSCR2,LUSCR3, 
     &       LUCBIO_SAVE, LUHCBIO_SAVE,LUC_SAVE
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from VB_BR_FROM_KAPPA '
        WRITE(6,*) ' ==========================='
        WRITE(6,*) ' NOOEXC_S, NOOEXC_A = ',
     &               NOOEXC_S, NOOEXC_A
      END IF
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Antisymmetric and symmetric part of Kappa '
       CALL WRTMAT(KAPPA_A,1,NOOEXC_A,1,NOOEXC_A)
       CALL WRTMAT(KAPPA_S,1,NOOEXC_S,1,NOOEXC_S)
*
       WRITE(6,*) ' IOOEXC: '
       CALL IWRTMA3(IOOEXC,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'VBBRKA')
*
*. Obtain New MO coefficients in MOAOUT: MOAOIN* Exp(-Kappa_A S) Exp(-Kappa_S S)
*
      CALL NEWMO_FROM_KAPPA_NORT(
     &     NOOEXC_A,IOOEXC_A,KAPPA_A,NOOEXC_S,IOOEXC_S,KAPPA_S,
     &     WORK(KMOAOIN),WORK(KMOAOUT))
*
* Obtain the set of biorthonormal orbitals
*
      CALL GET_CBIO(WORK(KMOAOUT),WORK(KCBIO),WORK(KCBIO2))
*
* Biorthonormal integral transformaion
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Bioorthogonal integral transformation '
      END IF
*
C     IE2LIST_A = IE2LIST_FULL_BIO
C     IOCOBTP_A = 1
C     INTSM_A = 1
      IE2LIST_A = IE2LIST_1G_BIO
C     IOCOBTP_A = 2
      IOCOBTP_A = 1
      INTSM_A = 1
      CALL PREPARE_2EI_LIST
*
      KKCMO_I = KMOAOUT
      KKCMO_J = KCBIO2
      KKCMO_K = KMOAOUT
      KKCMO_L = KCBIO2
*
C          DO_ORBTRA(IDOTRA,IDOFI,IDOFA,IE2LIST_IN,IOCOBTP_IN,INTSM_IN)
C     CALL DO_ORBTRA(1,1,1,IE2LIST_FULL_BIO,IOCOBTP_A,INTSM_A)
C     CALL FLAG_ACT_INTLIST(IE2LIST_FULL_BIO)
      CALL DO_ORBTRA(1,1,1,IE2LIST_1G_BIO,IOCOBTP_A,INTSM_A)
      CALL FLAG_ACT_INTLIST(IE2LIST_1G_BIO)

      NINT1_F = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1_F)
*
* And construct the one- and two-body density matrices
*
      CALL VB_DENSI(WORK(KRHO1),WORK(KRHO2),2,C,VEC2,VEC3)
*. Construct Active Fock-matrix
      CALL DO_ORBTRA(1,1,1,IE2LIST_FULL_BIO,
     &     IOCOBTP_A,INTSM_A)
*
      CALL FOCK_MAT_NORT(WORK(KF),WORK(KF2),2,WORK(KFI),WORK(KFA))
*. And the interspace gradient
C     E1_FROM_F_NORT(E1,F1,F2,IOPSM,IOOEXC,IOOEXCC,
C    &           NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
            CALL E1_FROM_F_NORT(BR,WORK(KF),WORK(KF2),1,
     &           IOOEXC,IOOEXC_A,NOOEXC_A,NTOOB,
     &           NTOOBS,NSMOB,IBSO,IREOST)
*. And add the active-active gradient
* The interspace excitations
C           VB_GRAD_ORBVBSPC(NOOEXC,IOOEXC,E1,C,VEC1_CSF,VEC2_CSF)
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' Active-active gradient will be calculated '
            CALL VB_GRAD_ORBVBSPC(NOOEXC_S,IOOEXC_S,
     &      BR(1+NOOEXC_A-NOOEXC_S),C,VEC2,VEC3)
                  
* And calculate gradient
C     VB_GRAD_ORBVBSPC(NOOEXCA,IOOEXC,E1,C,
C    &           VEC1_CSF,VEC2_CSF)
COLD  CALL VB_GRAD_ORBVBSPC(NOOEXC_A,IOOEXC_A,BR,C,VEC2,VEC3)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The Brilloin vector as delivered by VEC_BR_FRO..'
        WRITE(6,*) ' ================================================='
        CALL WRTMAT(BR,NOOEXC_A+NOOEXC_S,1,NOOEXC_A+NOOEXC_S,1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'VBBRKA')
      RETURN
      END
      SUBROUTINE CSDTVC_CONFSPACE(NCONF,VCSF,VSD,ISYM,ISPC,IWAY)
*
* Transform a CI vector between CSF and SD form for configuration 
* expansion using on-flight generation of info
*
*. Jeppe Olsen, Kristiansand, June 11, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'spinfo.inc'
*
      PARAMETER (LSCR = 1000)
*. Input / output
      DIMENSION VCSF(*), VSD(*)
*
*. Local scratch - is not general pt....
*
      DIMENSION IOCC(LSCR), ISIGN(LSCR), ISCR(LSCR)
*
      NTEST = 100
*
      WRITE(6,*) ' CSDTVC_CONFSPACE, Preliminary version '
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CSDTVC_CONFSPACE '
        WRITE(6,*) ' ============================ '
        WRITE(6,*)
        WRITE(6,*) ' Space and sym: ', ISPC, ISYM
        WRITE(6,*) ' IWAY = ', IWAY
      END IF
*
      INI = 1
      IB_CSF = 1
      IB_SD = 1
*
      DO ICONF = 1, NCONF
C            NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
        CALL NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,ISYM,ISPC,NEW)
        INI = 0
        IOCOB = (IOPEN + N_EL_CONF)/2
*. Signs for going between configuration and interaction order of dets
C            SIGN_CONF_SD(ICONF,NOB_CONF,IOP,ISGN,IPDET_LIST,ISCR)
        CALL SIGN_CONF_SD(IOCC,IOCOB,IOPEN,ISIGN,WORK(KDFTP),ISCR)
        NCSF = NPCSCNF(IOPEN+1)
        NSD  = NPDTCNF(IOPEN+1)
C             CSDTVC_CONF(C_SD,C_CSF,NOPEN,ISIGN,IAC,IWAY)
        CALL CSDTVC_CONF(VSD(IB_SD),VCSF(IB_CSF),IOPEN,ISIGN,2,IWAY)
        IB_CSF = IB_CSF + NCSF
        IB_SD = IB_SD + NSD
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CSDTVC_CONFSPACE:'
        NCSFT = IB_CSF-1
        NSDT = IB_SD - 1
        WRITE(6,*) ' CSF expansion: '
        CALL WRTMAT(VCSF,1,NCSFT,1,NCSFT)
        WRITE(6,*) ' SD expansion '
        CALL WRTMAT(VSD,1,NSDT,1,NSDT)
      END IF
*
      RETURN
      END
        
c $Id$
