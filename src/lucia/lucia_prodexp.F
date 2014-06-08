* Some more nomenclature 
* first of all, I will use the word atom for orbitalsubspace, thus 
* instead of talking about orbital subspace K, I will say atom K.
*
* So the input is a set of configurations and spins for each atom.
* An atomic state (AS) for atom K specifies : Configuration 
*                                            Spin-multiplicity 
*                                            2*MS value
* An atomic term (AT) for atom K specifies :  configuration 
*                                            Spin-multiplicity 
* Thus a given atomic term typically contains more than 
* one atomic state, differing in their 2*MS value 
* An atomic configuration (AC) specifies : configuration ( surprised ?)
* An atomic configuration may comprise several several (active)
* atomic  terms
*
* A PAS (product of atomic states)  is a product of atomic states
* A PAT (product of atomic terms) is a product of atomic terms.
* A PAT therefore typically contains many PAS'es, differing in 
* the MS2-values of the various atoms.
* A PAC (product of atomic configurations) is a product of 
* atomic configurations
* 
* The CI-expansion will thus be 
* Loop over the product atomic configurations 
*   Loop over the product atomic terms this product atomic configuration 
*     Loop over product atomic states for this product atomic term
*     End of loop over PAS
*   End of loop over PAT
* End of loop over PAC

      SUBROUTINE LUCIA_PRODEXP
*
* Control routine for Wavefunctions built as products of wave functions for each subspace
*
*. Jeppe Olsen, October 05
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'prdwvf.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Local scratch 
      INTEGER IPACSPIN(MXPNGAS)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PRODEX')
      NTEST = 1000
      IF(NTEST.GE.2) THEN 
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*) '***********************************************'
        WRITE(6,*) '* LUCIA : Expansion of wavefunction as        *'
        WRITE(6,*) '*         products of subspace wave functions *'
        WRITE(6,*) '*                                             *'
        WRITE(6,*) '*         Version of Oct 05                   *'
        WRITE(6,*) '***********************************************'
        WRITE(6,*)
        WRITE(6,*)
      END IF
*
      IF(NTEST.GE.1) THEN
       WRITE(6,*) ' Included configurations for the various subspaces :'
       WRITE(6,*) ' ==================================================='
       WRITE(6,*)
       DO IGAS = 1, NGAS
         WRITE(6,*) '   Orbital subspace ', IGAS 
         WRITE(6,*) '   *************************'
         WRITE(6,*) '   Number of configurations included ', 
     &   NWF_PER_SUBSPC(IGAS)
         DO ICONF = 1,  NWF_PER_SUBSPC(IGAS)
            WRITE(6,*) '     Configuration : ', ICONF
            WRITE(6,'(A,20I3)') '        Occupation : ',  
     &      (ISUBSPCWF_OCC(IORB,ICONF,IGAS),IORB=1,NOBPT(IGAS))
            WRITE(6,*) '       Multiplicity ', 
     &      MULT_FOR_SUBSPCWF(ICONF,IGAS)
         END DO
       END DO
       WRITE(6,*) ' Allowed excitation level between subspaces',
     &            INTRA_EXC_PRWF
      END IF
*. Pt I will assume only a singe wf per atom and highspin
*. Check number of wf 
      NPAT_INPUT = 1
      DO IAT = 1, NGAS
        NPAT_INPUT = NPAT_INPUT*NWF_PER_SUBSPC(IAT)
      END DO
      IF(NPAT_INPUT.NE.1) THEN
        WRITE(6,*) ' Input corresponds to more than one PAT '
        WRITE(6,*) ' I am only programmed for one PAT '
        STOP       ' Input corresponds to more than one PAT '
      END IF
*. Check also that input corresponds to high-spin case
      IHIGH_SPIN = 1
      DO IAT = 1, NGAS
        NOPEN = 0
        DO IOB = 1, NOBPT(IAT)
          IF(ISUBSPCWF_OCC(IOB,1,IAT).EQ.1) NOPEN = NOPEN + 1
        END DO
*
        IMAX_MULT = NOPEN + 1

        IF(MULT_FOR_SUBSPCWF(1,IAT).NE.IMAX_MULT) THEN
          WRITE(6,*) ' Not high-spin for atom ', IAT
          WRITE(6,*) ' Actual mult and high-spin mult ',
     &    MULT_FOR_SUBSPCWF(1,IAT), IMAX_MULT
          IHIGH_SPIN = 0
        END IF
      END DO
      IF(IHIGH_SPIN.EQ.0) THEN
        WRITE(6,*) ' Reference PAT is not high-spin '
        STOP       ' Reference PAT is not high-spin '
      END IF
*. We are now sure that there is just one reference PAT and 
*. it is high-spin, construct/reconstruct reference PAC and PAT
      CALL MEMMAN(KPACREF,NTOOB,'ADDL  ',1,'PACREF')
      CALL MEMMAN(KLSCR1,NGAS,'ADDL  ',1,'ISCR1 ')
      IONE = 1
      CALL ISETVC(WORK(KLSCR1),IONE,NGAS)
      CALL GET_PAC_FROM_INI_AC(WORK(KPACREF),NGAS,NOBPT,WORK(KLSCR1))
*. Number of open shells in reference
C     GET_PAC_FROM_INI_AC(IPAC,NAT,NOBAT,IACFAT)
*
* Generate all PAC's that may be obtained by applying excitations
* to the reference PAC
C     EXC_PAC(IFLAG,IOC_PAC,NAT,NOBPAT,NOBAC,NEXC,NEXC_PAC,
C    &                   IOC_EXPAC,INTRA,ICHKSM,IREFSM)
*. At the moment : no check of symmetry, no intraatomic excitations
      INTRA = 0
      ICHKSM = 0
      IADREF = 1
*. Number of PACS
      CALL EXC_PAC(1,WORK(KPACREF),NGAS,NOBPT,NTOOB,INTRA_EXC_PRWF,
     &             NEXC_PAC,IDUMMY,INTRA,ICHKSM,0,IADREF)
*. And the actual PACS
      LEN_PACS = NTOOB*NEXC_PAC
      CALL MEMMAN(KPACS,LEN_PACS,'ADDL  ',1,'PACS  ')
      CALL EXC_PAC(0,WORK(KPACREF),NGAS,NOBPT,NTOOB,INTRA_EXC_PRWF,
     &             NEXC_PAC,WORK(KPACS),INTRA,ICHKSM,0,IADREF)
*. Give the various PACs high spin giving one PAT for each PAC
      LEN_SPINAT = NEXC_PAC*NGAS
      CALL MEMMAN(KPACSPIN,LEN_SPINAT,'ADDL  ',1,'PACSPN')
      CALL Z_HIGH_SPIN_PACS(NEXC_PAC,WORK(KPACS),NTOOB,NOBPT,NGAS,
     &     WORK(KPACSPIN))
      NPAT = NEXC_PAC
C     Z_HIGH_SPIN_PACS(NPAC,IPACOCC,NACOB,NOBAT,NAT,IPATSPIN)
*
* =============================================================
*. Dimension of the PAS expansions for various values of MS2
* =============================================================
*
      NOPEN_REF = NOPEN_FOR_PAC(WORK(KPACREF),NTOOB)
C     NOPEN_FOR_PAC(IPACOCC,NACOB)
      WRITE(6,*) ' Number of singly occupied orbitals in ref PAC',
     &           NOPEN_REF
C     DIM_FOR_PASEXPANSION(IPATSPIN,NPAT,
C    &           NAT,MS2TOT,NPAS,XPAS)
      MS2_MIN = MOD(NOPEN_REF,2)
      CALL MEMMAN(KNPAS_PER_PAT,NEXC_PAC,'ADDL  ',2,'PASPPA')
      
      MXLENPAT = 0
      DO  MS2TOT = MS2_MIN, NOPEN_REF,2
        CALL DIM_FOR_PASEXPANSION(WORK(KPACSPIN),NEXC_PAC,NGAS,MS2TOT,
     &               NPAS,XPAS,MXLENPAT_L,WORK(KNPAS_PER_PAT))
        MXLENPAT = MAX(MXLENPAT,MXLENPAT_L)
        WRITE(6,*) ' Dimension for MS2 = ', MS2TOT, ' is ', NPAS, XPAS
        WRITE(6,*) ' Largest dimension of a single PAT ', MXLENPAT_L
      END DO
*
*
* ========================================
* Find the various atomic prototype states 
* ========================================
*
C     FIND_PROTO_AS(NPAT,IPATSPIN,IPATOCC,NAT,IFLAG,
C    &                       NOBAC,NOBPAT,NAS_PROTO,IAS_PROTO)
*. Number 
       CALL FIND_PROTO_AS(NPAT,WORK(KPACSPIN),WORK(KPACS),NGAS,1,
     &                    NTOOB,NOBPT,NAS_PROTO, IDUMMY)
*. And the actual prototype atomic states
       CALL MEMMAN(KPROTO_AS,3*NAS_PROTO,'ADDL  ',2,'AS_PRO')
       CALL FIND_PROTO_AS(NPAT,WORK(KPACSPIN),WORK(KPACS),NGAS,2,
     &                    NTOOB,NOBPT,NAS_PROTO, WORK(KPROTO_AS))
*. 
* ===================================================
* Prototype information about the various prototypes 
* ===================================================
*. 
*. The info concerning SD's and coupling and transformation
*  will be constructed for each prototype atomic state. 
*. Define pointers to the arrays of pointers ...
*. ( So I am defining pointers to pointers, I nearly feel 
*    that I am a computer scientist..)
      CALL MEMMAN(KNSD_PROTO_AS,NAS_PROTO,'ADDL  ',2,'NSDPRO')
      CALL MEMMAN(KNCSF_PROTO_AS,NAS_PROTO,'ADDL  ',2,'NCSFPR')
      CALL MEMMAN(KSD_PROTO_AS,NAS_PROTO,'ADDL  ',2,'SD_PRO')
      CALL MEMMAN(KCSF_PROTO_AS,NAS_PROTO,'ADDL  ',2,'CSFPRO')
      CALL MEMMAN(KCSDCSF_PROTO_AS,NAS_PROTO,'ADDL  ',2,'CSDCSF')
C     GEN_SPINFO_FOR_ALL_PROTO_AS(NPROTO_AS,IPROTO_AS,
C    &           NSD_PROTO,NCSF_PROTO,KSD_PROTO,KCSF_PROTO,
C    &           KCSDCSF_PROTO)
      CALL GEN_SPINFO_FOR_ALL_PROTO_AS(NAS_PROTO,WORK(KPROTO_AS),
     &     WORK(KNSD_PROTO_AS),WORK(KNCSF_PROTO_AS),
     &     WORK(KSD_PROTO_AS),WORK(KCSF_PROTO_AS),WORK(KCSDCSF_PROTO_AS)
     &     )
*
*. Generate the various PAS'ses for each PAT -testing 
*
      LEN_PASMS2 = MXLENPAT*NGAS
      CALL MEMMAN(KIPAS,LEN_PASMS2,'ADDL  ',2,'PASOCC')
      CALL MEMMAN(KIPASSCR,LEN_PASMS2,'ADDL  ',2,'PASCR')
      DO MS2TOT = MS2_MIN, NOPEN_REF,2
*. Dimensions for this MS2TOT
        CALL DIM_FOR_PASEXPANSION(WORK(KPACSPIN),NEXC_PAC,NGAS,MS2TOT,
     &               NPAS,XPAS,MXLENPAT,WORK(KNPAS_PER_PAT))
       DO IPAS = 1, NEXC_PAC
*. 
         DO IAT = 1, NGAS
           IPACSPIN(IAT) = IFRMR(WORK,KPACSPIN,(IPAS-1)*NGAS+IAT)
         END DO
         NPAS_FOR_PAC = IFRMR(WORK,KNPAS_PER_PAT,IPAS)
C             GENPAS_FOR_PAT(IMULT,NAT,MS2TOT,NPAS,IPAS,ISCR)
         CALL GENPAS_FOR_PAT(IPACSPIN,NGAS,MS2TOT,NPAS_FOR_PAC,
     &        WORK(KIPAS),WORK(KIPASSCR))
       END DO
      END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PRODEX')
      RETURN
      END
      SUBROUTINE ATCNF_FOR_ATTRM(IFLAG,NATCNF,NATTRM,IATTRM,NORB,
     &           IATCNF)
*
* Atomic terms => Atomic configurations
*
* ( There may be several atomic terms with the same configuration, and 
*   this routine thus identifies the common configurations)
*
* Jeppe Olsen, Warwick Oct 5, 2005
* 
* ======
*  Input
* ======
*
* IFLAG : =1 : Number of atomic configurations 
*         =2 : and the atomic configurations 
* NATTRM : Number of atomic terms
* IATTRM : The configuration part of the various terms
* NORB   : Number of orbitals
*
*
* ========
*  Output
* ========
* NATCNF : Number of (different) atomic configurations
* IATCNF : The various atomic configurations.
*
      INCLUDE 'implicit.inc'
*. Input : configuration part of the various terms
      DIMENSION IATTRM(NORB,NATTRM)
*. Output ( if IFLAG = 2)
      DIMENSION IATCNF(NORB,*)
*
      NATCNF = 0
      DO JATTRM1 = 1, NATTRM
*. Is configuration part of this term different from previous confs
       INEW = 1
       DO JATTRM2 = 1, JATTRM1 -1
         CALL COMPARE_TWO_INTARRAYS(IATTRM(1,JATTRM1),IATTRM(1,JATTRM2),
     &                              NORB,IDENT)
         IF(IDENT.EQ.1) INEW = 0
       END DO
       IF(INEW.EQ.1) THEN
         NATCNF = NATCNF + 1
         IF(IFLAG.EQ.2) THEN
           CALL ICOPVE(IATTRM(1,JATTRM1),IATCNF(1,NATCNF),NORB)
         END IF
       END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Atomic terms => atomic configurations '
        WRITE(6,*) ' Number of atomic configurations ', NATCNF
        IF(IFLAG.EQ.2) THEN
          WRITE(6,*) ' The atomic configurations (as columns)'
          CALL IWRTMA(IATCNF,NORB,NATCNF,NORB,NATCNF)
        END IF
      END IF
*
      RETURN 
      END
      
   
      SUBROUTINE COMPARE_TWO_INTARRAYS(IA,IB,NAB,IDENT)
*
* Compare two arrays of integers, IA, IB and return IDENT =1
* if they are identical, IDENT = 0 if they differ
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IA(NAB),IB(NAB)
*
      IDENT = 1
      DO I = 1, NAB
        IF(IA(I).NE.IB(I)) IDENT = 0
      END DO
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Integer arrays IA, IB '
        CALL IWRTMA(IA,1,NAB,1,NAB)
        CALL IWRTMA(IB,1,NAB,1,NAB)
        IF(IDENT.EQ.1) THEN
          WRITE(6,*) ' Arrays are identical '
        ELSE 
          WRITE(6,*) ' Arrays differ '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE EXC_PAC(IFLAG,IOC_PAC,NAT,NOBPAT,NOBAC,NEXC,NEXC_PAC,
     &                   IOC_EXPAC,INTRA,ICHKSM,IREFSM,IADREF)
*. A product of atomic configurations is given in IOC_PAT.
*. Find all other atomic configurations that my be 
*. obtained by 1-NEXC-fold excitations.
*. Only excitations between different atoms are included
*. if ICHKSM = 1, then only configurations with symmetry = IREFSM
*. are included
*. if INTRA = 1, intraatomic excitations are allowed
*
* IFLAG = 1 : Only the number of excited configurations is returned
*
* ======
*. Input
* ======
* IFLAG = 1 : Calculate only number of excited configurations
* IFLAG.ne.1 : Calculate number of excited configurations and 
*              return occupations of these configurations
* NAT       : Number of atoms (number of orbital subspaces)
* NOBPAT    : Number of orbitals per atom
* NOBAC     : Total number of orbitals 
* NEXC      : Allowed number of excitations 
* INTRA = 1 : Allow intraatomic excitations
*      ne 1 : no intraatomic excitations
* ICHKSM    : Check symmetry of the various configurations 
* IREFSM    : Required symmetry of configurations (if ICHKSM = 1)
* IADREF    : The reference PAC is included in the output list

* ========
*. Output
* ========
* NEXC_PAC   : Number of excited PACs
* I_OC_EXPAC : Occupation of the excited configurations (IFLAG.NE.1)
*
*. Jeppe Olsen, Warwick, October 5, 2005
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER IOC_PAC(NOBAC), NOBPAT(NAT)
*. Output (if IFLAG.NE.1)
      INTEGER IOC_EXPAC(NOBAC,*)
*. Local scratch
      INTEGER ISCR(MXPORB), IATFORB(MXPORB)
*
      NTEST = 1000
*
      IF(NEXC.GT.2) THEN
        WRITE(6,*) ' EXC_PAC called for exc. level = ', NEXC
        WRITE(6,*) ' EXC_PAC is pt only programmed for NEXC upto',2
        STOP ' EXC_PAC called for too large exc. level '
      END IF
*. Total number of orbitals
      NOB_TOT = 0
      DO IAT = 1, NAT
       NOB_TOT = NOB_TOT + NOBPAT(IAT)
      END DO
      WRITE(6,*) ' EXC_PAC, NOB_TOT = ', NOB_TOT
*. Set up array orbital => atom
      DO IAT = 1, NAT
        IF(IAT.EQ.1) THEN
          ISTART = 1
        ELSE 
          ISTART = ISTART + NOBPAT(IAT-1)
        END IF
        ISTOP = ISTART + NOBPAT(IAT) - 1
        DO IORB = ISTART, ISTOP
         IATFORB(IORB) = IAT
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Orbital => atom map '
        CALL IWRTMA(IATFORB,1,NOB_TOT,1,NOB_TOT)
      END IF
*. Include reference PAC if required
      IF(IADREF.EQ.0) THEN
        NCONF = 0
      ELSE 
        NCONF = 1
      END IF
      IF(IADREF.EQ.1.AND.IFLAG.NE.1) THEN
        CALL ICOPVE(IOC_PAC,IOC_EXPAC(1,1),NOBAC)
      END IF
*
      DO IEXC = 1, NEXC
*
      IF(IEXC.EQ.1) THEN
*
* Single excitations
*
        DO IANNI = 1,  NOBAC
          DO ICREA = 1, NOBAC
            IF((IANNI.NE.ICREA).AND.
     &         (INTRA.EQ.1.OR.IATFORB(IANNI).NE.IATFORB(ICREA))) THEN
              CALL ICOPVE(IOC_PAC(1),ISCR(1),NOBAC)
              ISCR(IANNI) = ISCR(IANNI) - 1
              ISCR(ICREA) = ISCR(ICREA) + 1
*. Ensure that changed occupations are within bounds
              IF(ISCR(IANNI).GE.0.AND.ISCR(ICREA).LE.2) THEN
*. Occupation okay, check symmetry if required
                ISYMOK = 1
                IF(ICHKSM.EQ.1) THEN
                 ISM = ISYMCN2(ISCR,NOBAC)
                 IF(ISM.NE.IREFSM) ISYMOK = 0
                END IF
                IF(ISYMOK.EQ.1) THEN
*. Yet another allowed configuration has been found 
                  NCONF = NCONF + 1
                  IF(IFLAG.NE.1) THEN
                    CALL ICOPVE(ISCR,IOC_EXPAC(1,NCONF),NOBAC)
                  END IF
                END IF
              END IF
*             ^ End if occupations were inside bounds 
            END IF
*           ^ End if excitation was included
          END DO
        END DO
*       ^ Ends of loop over IANNI, ICREA
      ELSE IF(IEXC.EQ.2) THEN
*
*. Double excitations
*
        DO IANNI = 1,  NOBAC
         DO JANNI = 1, IANNI
          DO ICREA = 1, NOBAC
           DO JCREA = 1, ICREA
*. Test whether an excitation contains an intraatomic part
            I_AM_INTRA = 0
            IF(IATFORB(IANNI).EQ.IATFORB(ICREA).OR.
     &         IATFORB(IANNI).EQ.IATFORB(JCREA).OR.
     &         IATFORB(JANNI).EQ.IATFORB(ICREA).OR.
     &         IATFORB(JANNI).EQ.IATFORB(JCREA)    ) I_AM_INTRA = 1
            IF((IANNI.NE.ICREA.AND.JANNI.NE.JCREA).AND.
     &         (INTRA.EQ.1.OR.I_AM_INTRA.EQ.1)) THEN 
              CALL ICOPVE(IOC_PAC(1),ISCR(1),NOBAC)
              ISCR(IANNI) = ISCR(IANNI) - 1
              ISCR(ICREA) = ISCR(ICREA) + 1
              ISCR(JANNI) = ISCR(JANNI) - 1
              ISCR(JCREA) = ISCR(JCREA) + 1
*. Ensure that changed occupations are within bounds
              IN_BOUNDS = 0
              IF(ISCR(IANNI).GE.0.AND.ISCR(ICREA).LE.2.AND.
     &           ISCR(JANNI).GE.0.AND.ISCR(JCREA).LE.2) THEN
                 IN_BOUNDS = 1
              ELSE 
                 IN_BOUNDS = 0
              END IF
              IF(IN_BOUNDS.EQ.1) THEN
*. Occupation okay, check symmetry if required
                ISYMOK = 1
                IF(ICHKSM.EQ.1) THEN
                 ISM = ISYMCN2(ISCR,NOBAC)
                 IF(ISM.NE.IREFSM) ISYMOK = 0
                END IF
                IF(ISYMOK.EQ.1) THEN
*. Yet another allowed configuration has been found 
                  NCONF = NCONF + 1
                  IF(IFLAG.NE.1) THEN
                    CALL ICOPVE(ISCR,IOC_EXPAC(1,NCONF),NOBAC)
                  END IF
                END IF
              END IF
*             ^ End if occupations were inside bounds 
            END IF
*           ^ End if excitation was included
           END DO
          END DO
         END DO
        END DO
*       ^ Ends of loops over IANNI, JANNI, ICREA, JCREA
      END IF
*     ^ End of switch between different excitations levels
      END DO
*     ^ End of loop over excitation levels
      NEXC_PAC = NCONF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Number of configurations generated ', NCONF
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' List of generated configurations : '
          WRITE(6,*) ' ==================================='
          DO ICONF = 1, NCONF
            WRITE(6,*) ' Configuration ', ICONF
            WRITE(6,'(20I3)') (IOC_EXPAC(I,ICONF),I=1,NOBAC)
          END DO
        END IF
      END IF
*
      RETURN
      END
      FUNCTION ISYMCN2(IOCC,NORB)
*
* The occupation of a configuration is given as 
* a set of integers : 0 => unoccupied 
*                    -1,2 => doubly occupied
*                     1 => single occupied
* Find symmetry of configuration
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER IOCC(*)
*
      ISYM = 1
      DO IORB = 1, NORB
        IF(IOCC(IORB).EQ.1) ISYM = MULTD2H(ISYM,ISMFTO(IORB))
      END DO
      ISYMCN2 = ISYM
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation of configuration '
        CALL IWRTMA(IOCC,1,NORB,1,NORB)
        WRITE(6,*) ' Symmetry = ', ISYM 
      END IF
*
      RETURN
      END
      SUBROUTINE Z_HIGH_SPIN_PACS(NPAC,IPACOCC,NACOB,NOBAT,NAT,IPATSPIN)
*
* NPAC PAC's are given in IPAC
* Construct corresponding high-spin PAT's by giving each PAC
* high-spin
*
*. Jeppe Olsen, Oct. 5. 2005 in Warwick
*
      INCLUDE 'implicit.inc'
*.Input
      INTEGER IPACOCC(NACOB,NPAC), NOBAT(NAT)
*. Output
      INTEGER IPATSPIN(NAT,NPAC)
*
      DO IPAC = 1, NPAC
        DO IAT = 1, NAT
          NOPEN = 0
          IF(IAT.EQ.1) THEN
            IOB_OFF = 1
          ELSE 
            IOB_OFF = IOB_OFF + NOBAT(IAT-1)
          END IF
          DO IOB = 1, NOBAT(IAT)
            IF(IPACOCC(IOB+IOB_OFF-1,IPAC).EQ.1) NOPEN = NOPEN + 1
          END DO
          IPATSPIN(IAT,IPAC) = NOPEN + 1
        END DO
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' High-spin PAT''s generated '
        IONLYOCC = 0
        CALL WRITE_PATLIST(IPACOCC,IPATSPIN,NPAC,NOBAT,NAT,NACOB,
     &                     IONLYOCC)
      END IF
*
      RETURN
      END
      SUBROUTINE WRITE_PATLIST(IPACOCC,IPATSPIN,NPAT,NOBAT,NAT,NACOB,
     &                         IONLYOCC)
*
* Write configurations and multiplicites (if IONLYOCC = 0)
* for a list of PATS
*
* Jeppe Olsen, Oct. 6, Warwick
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IPACOCC(NACOB,NPAT), IPATSPIN(NAT,NPAT)
      INTEGER NOBAT(NAT)
*
      IF(IONLYOCC.EQ.0) THEN
        WRITE(6,*) ' ( Multiplicity) and occupation for each atom : '
      ELSE 
        WRITE(6,*) ' Occupation for each atom : '
      END IF
      DO IPAT = 1, NPAT
        WRITE(6,*) '   Product of atomic terms ', IPAT
        WRITE(6,*) '   ================================'
        WRITE(6,*) 
        DO IAT = 1, NAT
          IF(IAT.EQ.1) THEN
            IOB_OFF = 1
          ELSE 
            IOB_OFF = IOB_OFF + NOBAT(IAT-1)
          END IF
          NOB = NOBAT(IAT)
          IF(IONLYOCC.EQ.0 ) THEN
            WRITE(6,'(5X,A,I3,A,30I2)') '(',IPATSPIN(IAT,IPAT),')',
     &      (IPACOCC(IOB_OFF+IOB-1,IPAT),IOB= 1, NOB)
          ELSE
            WRITE(6,'(5X,30I2)') 
     &      (IPACOCC(IOB_OFF+IOB-1,IPAT),IOB= 1, NOB)
          END IF
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GEN_PAT_FROM_PAC(IPAC,NAT,NOBAT,IMULTP,NPAT,IMULTPAT,
     &                            IFLAG)
*
* Generate spins  for the various atoms for a given product of 
* atomic configurations. The spins are returned as multiplicities.
* Type of generated spins are defined by IMULTTP
*
* IMULTP = 0 : Generate only highspin for each atom
*        = N(>0) : Allow multiplicities Highspin - N -- Highspin
*        = N(<0) : Allow all spin-patterns
*
*. Jeppe Olsen, October 5, 2005, Warwick
*
* =======
*. Input
* =======
*
* IPAC : Occupation of the PAC  ( product atomic configuration)
* NAT :  Number of atoms
* NOBAT : Number of orbitals per atom
* IMULTP  : Specifies allowed spins as given above
* IFLAG : = 1 =>calculate only NPAT
*         ne.1 => Determine also IMULTPAT
* 
* ========
*  Output
* ========
* NPAT : Number of obtained Product atomic terms 
* IMULTPAT : Multiplicities of the various atoms for each PAT
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER IPAC(*)
      INTEGER NOBAT(NAT)
*. Output (if IFLAG.ne.1)
      INTEGER IMULTPAT(NAT,*)
*. Local scratch
      INTEGER IOPEN(MXPNGAS)
*
      DO IAT = 1, NAT
*. Number of singly occupied orbitals for this atom
        IF(IAT.EQ.1) THEN
          IOB_START = 1
        ELSE 
          IOB_START = IOB_START + NOBAT(IAT-1)
        END IF
        IOB_END = IOB_START + NOBAT(IAT) - 1
        NOPEN = 0
        DO IOB = IOB_START, IOB_END
          IF(IPAC(IOB).EQ.1)  NOPEN = NOPEN + 1
        END DO
        IOPEN(IAT) = NOPEN
      END DO
*
* Number of PATs
*
      IF(IMULTP.EQ.0) THEN
         NPAT = 1
      ELSE 
         NPAT = 1
         DO IAT = 1, NAT
           MULT_MAX = 2*IOPEN(IAT) + 1
           IF(IMULTP.GT.0) THEN
             MULT_MIN = MAX(0,MULT_MAX-IMULTP)
           ELSE 
             MULT_MIN = MOD(IOPEN(IAT),2)
           END IF
           NPAT = NPAT*(MULT_MAX-MULT_MIN+1)
         END DO
       END IF
*
*. and the actual PATS
*
      IF(IFLAG.NE.1) THEN
        IF(IMULTP.EQ.0) THEN
          DO IAT = 1, NAT
            MULT_MAX = 2*IOPEN(IAT)+1
            IMULTPAT(IAT,1) = MULT_MAX
          END DO
        ELSE 
          WRITE(6,*) 
     &    ' GEN_PAT_FROM_PAC has not been programmed for IMULTP.NE.0'
          STOP       
     &    ' GEN_PAT_FROM_PAC has not been programmed for IMULTP.NE.0'
        END IF
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GEN_PAT_FROM_PAC '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Number of generated PAT''s ', NPAT
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' And the multiplicities of the various pats '
          DO IPAT = 1, NPAT
            WRITE(6,*) ' PAT  ', IPAT
            CALL IWRTMA(IMULTPAT(I,IPAT),1,NAT,1,NAT)
          END DO
        END IF
      END IF
*
      RETURN
      END
      FUNCTION NPAS_FOR_PAT(IMULT,NAT,MS2TOT)
*
* The multiplicities of each atom in a PAT is given by IMULT
* Find the number of PAS for this PAT with total ms2 value MS2TOT
*
*. Jeppe Olsen, Warwick, Oct 5, 2005
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. input
      INTEGER IMULT(NAT)
*. Local scratch 
      INTEGER MS2_MINMAX(MXPNGAS,2)
      PARAMETER(NMXCOMP = 1000)
      INTEGER NCOMP_IN(NMXCOMP), NCOMP_OUT(NMXCOMP)
*
*. Set up Min and max arrays for MS2
*
      CALL MINMAX_MS2_FOR_PAT(IMULT,NAT,MS2TOT,MS2_MINMAX)
*
* We now have min and max of accumulated MS2TOT. Find number of components
*
*. Largest and smallest possible MSvalues
      MAX_MS2 = 0
      DO IAT = 1, NAT
       MAX_MS2 = MAX_MS2 + IMULT(IAT) -1
      END DO
      MIN_MS2 = -MAX_MS2
*. A given MS2 will in the following be stored as MS2-MIN_MS2 + 1
      IF(2*MAX_MS2+1 .GT. NMXCOMP) THEN
        WRITE(6,*) ' NPAS_FOR_PAT, NMXCOMP too small'
        WRITE(6,*) ' Current and required value ', NMXCOMP,2*MAX_MS2+1
        STOP       ' NPAS_FOR_PAT, NMXCOMP too small'
      END IF
      NCOMP = 2*MAX_MS2+1
      IZERO = 0
*. Allowed MS2 values for atom 1
      CALL ISETVC(NCOMP_IN,IZERO,NCOMP)
      DO MS2 = MS2_MINMAX(1,1),MS2_MINMAX(1,2),2
        NCOMP_IN(MS2-MIN_MS2+1) = 1
      END DO
      DO IAT = 1, NAT-1
*. Atom IAT  => Atom IAT +1 
        CALL ISETVC(NCOMP_OUT,IZERO,NCOMP)
        IDELTA = IMULT(IAT+1)-1
        DO MS2_IN = MS2_MINMAX(IAT,1),MS2_MINMAX(IAT,2),2
          DO IIDELTA = - IDELTA, IDELTA,2
            MS2_OUT = MS2_IN + IIDELTA
C?          WRITE(6,*) ' IAT, MS2_IN, IIDELTA, MS2_OUT ',
C?   &                   IAT, MS2_IN, IIDELTA, MS2_OUT 
            IF(MS2_MINMAX(IAT+1,1).LE.MS2_OUT.AND.
     &         MS2_MINMAX(IAT+1,2).GE.MS2_OUT     ) THEN
               NCOMP_OUT(MS2_OUT-MIN_MS2+1) 
     &       = NCOMP_OUT(MS2_OUT-MIN_MS2+1) + NCOMP_IN(MS2_IN-MIN_MS2+1)
C?           WRITE(6,*) ' Allowed ! '
            END IF
          END DO
        END DO
C?      WRITE(6,*) ' NCOMP_IN, NCOMP_OUT for atom ', IAT
C?      CALL IWRTMA(NCOMP_IN,1,NCOMP,1,NCOMP)
C?      CALL IWRTMA(NCOMP_OUT,1,NCOMP,1,NCOMP)
        CALL ICOPVE(NCOMP_OUT,NCOMP_IN,NCOMP)
      END DO
      NPAS = NCOMP_OUT(MS2TOT-MIN_MS2+1)
*. And calculate number of components
CE    NPAS = 1
CE    DO IAT = 1, NAT
CE      NPAS = NPAS*(MS2_MINMAX(IAT,2)-MS2_MINMAX(IAT,1)+1)
CE    END DO
*
      NPAS_FOR_PAT = NPAS
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of PAS for given PAT ', NPAS
      END IF
*
      RETURN
      END
      SUBROUTINE MINMAX_MS2_FOR_PAT(IMULT,NAT,MS2TOT,MS2_MINMAX)
*
* A PAT is defined by the multiplicities IMULT
* Obtained max/min for accumulated MS2 with the constraint that
* the total MS2 of the PAT is MS2TOT
*
*. Jeppe Olsen, Oct. 5, 2005, Warwick
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER IMULT(NAT)
*. Output 
      INTEGER MS2_MINMAX(MXPNGAS,2)
*
* built up from atom 1, no considerations that the total MS2 has to be 
* MS2TOT
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MINMAX_MS2_FOR_PAT speaking : '
        WRITE(6,*) ' =============================='
        WRITE(6,*) ' Multiplicities for each atom '
        CALL IWRTMA(IMULT,1,NAT,1,NAT)
        WRITE(6,*) ' Required MS2 = ',MS2TOT
      END IF
*
      MS2ACC_MIN = 0
      MS2ACC_MAX = 0
      DO IAT = 1, NAT
        IDMS2 = IMULT(IAT) - 1
        IF(IAT.EQ.1) THEN
          MS2_MINMAX(IAT,1) = -IDMS2
          MS2_MINMAX(IAT,2) =  IDMS2
        ELSE 
          MS2_MINMAX(IAT,1) = MS2_MINMAX(IAT-1,1)-IDMS2
          MS2_MINMAX(IAT,2) = MS2_MINMAX(IAT-1,2)+IDMS2
        END IF
      END DO
*
*. In the above, we skipped all considerations that the total MS2 should
*. be MS2TOT. Impose this constraint and modify MS2_MINMAX
      MS2_MINMAX(NAT,1) = MS2TOT
      MS2_MINMAX(NAT,2) = MS2TOT
      DO IAT = NAT-1,1,-1
*. We know the limits for MS2 for atom IAT+1, what are then the limits for
* for atom I
        IDMS2 = IMULT(IAT+1)-1
        MS2_MIN = MS2_MINMAX(IAT+1,1)-IDMS2
        MS2_MAX = MS2_MINMAX(IAT+1,2)+IDMS2
        MS2_MINMAX(IAT,1) = MAX( MS2_MINMAX(IAT,1),MS2_MIN)
        MS2_MINMAX(IAT,2) = MIN( MS2_MINMAX(IAT,2),MS2_MAX)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Smallest allowed values of acc. MS2 '
        CALL IWRTMA(MS2_MINMAX(1,1),1,NAT,1,NAT)
        WRITE(6,*) ' largest allowed values of acc. MS2 '
        CALL IWRTMA(MS2_MINMAX(1,2),1,NAT,1,NAT)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_PAC_FROM_INI_AC(IPAC,NAT,NOBAT,IACFAT)
*
* obtain PAC as product of initial AC's
*
*. Jeppe Olsen, Oct. 5, 2005
*
* ======
*. Input
* ======
*.
* NAT : Number of atoms
* IACFAT(I) : Select Initial AT IACFAT(I) for atom I
*
* =======
*. Output
* =======
* IPAC : The PAC for all atoms collected
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'prdwvf.inc'
*. Input
      INTEGER IACFAT(NAT), NOBAT(NAT)
*. Output
      INTEGER IPAC(*)
*
      DO IAT = 1, NAT
        IF(IAT.EQ.1) THEN
         IOB_OFF = 1
        ELSE 
         IOB_OFF = IOB_OFF + NOBAT(IAT-1)
        END IF
        NOB = NOBAT(IAT)
        DO IOB = 1, NOB
          IPAC(IOB_OFF-1+IOB) =  ISUBSPCWF_OCC(IOB,1,IAT)
        END DO
      END DO
      NOB_TOT = IOB_OFF -1 + NOB
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Occupation of initial PAC '
       CALL IWRTMA(IPAC,1,NOB_TOT,1,NOB_TOT)
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_FOR_PASEXPANSION(IPATSPIN,NPAT,
     &           NAT,MS2TOT,NPAS,XPAS,MXLENPAT,NPAS_PER_PAT)
*. A PAT expansion is defined by IPACOCC, IPATSPIN. Find 
*. the number of PAS for a given MS2 ( MS2TOT)
*
* MXLENPAT is the largest dimension of a single PAT
*
*. Jeppe Olsen, Oct. 6, Warwick
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IPATSPIN(NAT,NPAT)
*. Output
      INTEGER NPAS_PER_PAT(NPAT)
*
      NPAS = 0
      XPAS = 0
      MXLENPAT = 0
      DO IPAT = 1, NPAT
         NPAS_TERM =  NPAS_FOR_PAT(IPATSPIN(1,IPAT),NAT,MS2TOT)
C                     NPAS_FOR_PAT(IMULT,NAT,MS2TOT)
         NPAS_PER_PAT(IPAT) = NPAS_TERM
         NPAS = NPAS + NPAS_TERM
         MXLENPAT = MAX(MXLENPAT,NPAS_TERM)
         XPAS = XPAS + FLOAT(NPAS_TERM)
      END DO
*
C?    WRITE(6,*) ' Total number of parameters for MS2TOT = ', MS2TOT,
C?   &           ' is ', NPAS ,'(', XPAS, ')'
*
      RETURN
      END
      FUNCTION NOPEN_FOR_PAC(IPACOCC,NACOB)
*
* The occupation  of a PAC is given in IPACOCC, find
* number of open orbitals
*
*. Jeppe Olsen, Oct. 6, 2005 in Warwick
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IPACOCC(NACOB)
*
      NOPEN = 0
      DO IOB = 1, NACOB
        IF(IPACOCC(IOB).EQ.1) NOPEN = NOPEN + 1
      END DO
*
      NOPEN_FOR_PAC = NOPEN
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
       WRITE(6,'(A,30I3)') ' Occupation : ', (IPACOCC(I),I=1,NACOB)
       WRITE(6,'(A,I7)')   ' Number of singly occupied orbitals ',NOPEN
      END IF
*
      RETURN
      END 
      SUBROUTINE SIGMA_PRDWF1(C, SIGMA)
*
* Direct CI for product wavefunction. 
* outer routine adding 
* from /KPRDWVF/ and /PRDWVF/ and  MS2 value is obtained from 
* /CSTATE/
* Version -2
*
* Incore version 
* 
*. Jeppe Olsen, Warwick, Oct. 6, 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'kprdwvf.inc'
      INCLUDE 'prdwvf.inc'
      INCLUDE 'cstate.inc'
*. And get a slave to do the job
      CALL SIGMAS_PRDWF1(C,SIGMA,WORK(KPACS),WORK(KPACSPIN),NEXC_PAS,
     &                   NGAS,MS2)
*
      RETURN
      END 
      SUBROUTINE SIGMAS_PRDWF1(C,SIGMA,IPACOCC,IPACMULT,NPAT,NAT,MS2)
* 
* Initial slave routine for the direct CI part of product wave function 
*
*. Jeppe Olsen, Warwick Oct. 6, 2005
*
*. General input
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER IPACOCC(NACOB,NPAT),IPACMULT(NAT,NPAT)
      DIMENSION C(*)
*. Output
      DIMENSION SIGMA(*)
*. Local scratch
      INTEGER IDIFFAT(4)
*
*· Length of expansion 

*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'ADDL  ',IDUM,'SIGPR1')
*
*. Length of expansion 
*
C     DIM_FOR_PASEXPANSION(IPATSPIN,NPAT,
C    &           NAT,MS2TOT,NPAS,XPAS,MXLENPAT)
      CALL DIM_FOR_PASEXPANSION(IPACMULT,NPAC,NAT,MS2,NPAS,XPAS,
     &                          MXLENPAT)
*. (I assume that NPAS gives the number of terms ( no integer overflow) )
      ZERO = 0.0D0
      CALL SETVEC(SIGMA,ZERO,NPAS)
*. Loop over PATs in sigma and C and check if they connect
      ISOFF = 1
      DO ISPAT = 1, NPAT
        ICOFF = 1
        DO ICPAT = 1, NCPAT
*. Are these PATS connected by Hamiltonian , and if yes : what 
*. type of Hamiltonian connects
           CALL CONNECT_PAT(IPACOCC(1,ISPAT),IPACOCC(1,ICPAT),
     &                      IPACMULT(1,ISPAT),IPACMULT(1,ICPAT),
     &                      ICONNECT,
     &                      NAT,NOBPT,I1CREA,I2CREA,I1ANNI,I2ANNI,
     &                      NDIFF,NDIFFAT,IDIFFAT)
*. We know how the two PATS differ, go to the appropriate 
*. subroutine 
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE CONNECT_PAT(ILPAT,IRPAT,ILMULT,IRMULT,ICONNECT,NAT,
     %           NOBPT,I1CREA, I2CREA, I1ANNI,I2ANNI,NDIFF,NDIFFAT,
     &           IDIFFAT)
*
* Two PATs, ILPAT and IRPAT are given. Check whether they 
* are connected by atmost two-electron operator, and 
* find the minimal operator connecting the two PATs
*
* At the moment, the Spins of the various atoms are not checked !!!
*
*. Jeppe Olsen, Oct. 6, 2005 in Warwick
*
      INCLUDE 'implicit.inc'
*. General input
      INTEGER NOBPT(NAT)
*. specific input
      INTEGER ILPAT(*),IRPAT(*)
      INTEGER ILMULT(NAT),IRMULT(NAT)
*. Output
       INTEGER IDIFFAT(4)
*
      NDIFF = 0
      LDIFF = 0
      I1CREA = 0
      I2CREA = 0
      I1ANNI = 0
      I2ANNI = 0
      DO IAT = 1, NAT
        IF(IAT.EQ.1) THEN
          IOFF = 1
        ELSE 
          IOFF = IOFF + NOBPT(IAT-1)
        END IF
        NOB = NOBPT(IAT)
        DO IOB = 1, NOB
         IF(ILPAT(IOB+IOFF-1).NE.IRPAT(IOB+IOFF-1)) THEN
           NDIFF = NDIFF + IABS(ILPAT(IOB+IOFF-1)-IRPAT(IOB+IOFF-1))
           IF(LDIFF.EQ.0.OR.(LDIFF.GE.1.AND.IAT.NE.IDIFFAT(LDIFF)))THEN
             LDIFF = LDIFF + 1
             IDIFFAT(LDIFF) = IAT
           END IF
           IF(ILPAT(IOB+IOFF-1).EQ.IRPAT(IOB+IOFF-1)+2) THEN
*. Double creation required
             I1CREA = IOB
             I2CREA = IOB
           ELSE IF (ILPAT(IOB+IOFF-1).EQ.IRPAT(IOB+IOFF-1)+1) THEN
*. Single creation required
             IF(I1CREA.NE.0) THEN 
                I2CREA = IOB
             ELSE
                I1CREA = IOB
             END IF
           ELSE IF (ILPAT(IOB+IOFF-1).EQ.IRPAT(IOB+IOFF-1)-1) THEN
*. Single annihilation required
             IF(I1ANNI.NE.0) THEN 
                I2ANNI = IOB
             ELSE
                I1ANNI = IOB
             END IF
           ELSE IF (ILPAT(IOB+IOFF-1).EQ.IRPAT(IOB+IOFF-1)-2) THEN
*. Double annihilation required
             I1ANNI = IOB
             I2ANNI = IOB
           END IF
*          ^ End of switch between the various forms of differences
         END IF
*        ^ End of the two PATs differed here
        END DO
*       ^ End of loop over orbitals at atom
      END DO
*     ^ End of loop over atoms
      NDIFFAT = LDIFF
*
      IF(NDIFF.LE.4) THEN
*. There is connection 
       ICONNECT = 1
      ELSE
*. Well, no connection 
       ICONNECT = 0
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CONNECT_PAT '
        IF(ICONNECT.EQ.1) THEN
          WRITE(6,*) 
     &    ' Connection, number of different occupations ', NDIFF
        IF(NDIFF.EQ.2) THEN
          WRITE(6,*) ' I1CREA, I1ANNI = ', I1CREA, I1ANNI
        ELSE IF (NDIFF.EQ.4) THEN
          WRITE(6,*) ' I1CREA, I1ANNI, I2CREA, I2ANNI = ', 
     &                 I1CREA, I1ANNI, I2CREA, I2ANNI
        END IF
        WRITE(6,*) ' Number of atoms with different occupations',
     &             NDIFFAT
        WRITE(6,*) ' Atoms with differing occupations ',
     &             (IDIFFAT(I),I=1,NDIFFAT)
        ELSE
         WRITE(6,*) ' No connection '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE GENPAS_FOR_PAT(IMULT,NAT,MS2TOT,NPAS,IPAS,ISCR)
*
* The multiplicities of each atom in a PAT is given by IMULT and 
* the number of PAS'ses for this PAT is given by NPAS. Generate 
* the PAS'ses in the form of MS2-values for the various atoms.
*
*. Jeppe Olsen, November 05
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. input
      INTEGER IMULT(NAT)
*. Local scratch 
      INTEGER MS2_MINMAX(MXPNGAS,2)
      PARAMETER(NMXCOMP = 1000)
      INTEGER NCOMP_IN(NMXCOMP), NCOMP_OUT(NMXCOMP)
      INTEGER IOFF_IN(NMXCOMP), IOFF_OUT(NMXCOMP)
*. Scratch of the same dimension as IPAS
      INTEGER ISCR(NAT,NPAS)
*. Output
      INTEGER IPAS(NAT,NPAS)
*
      NTEST = 100 
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' All PAS''ses will be generated for PAT '
        CALL IWRTMA(IMULT,1,NAT,1,NAT)
        WRITE(6,*) ' Required MS2TOT = ', MS2TOT
        WRITE(6,*) ' Number of PAS''ses are (input) ', NPAS
      END IF
*
*. Zero so we can print and copy without problems
      IZERO = 0.0D0
      CALL ISETVC(IPAS,IZERO,NAT*NPAS)
      CALL ISETVC(ISCR,IZERO,NAT*NPAS)
*
* The Passes are first constructed in terms of accumulated MS2, and are then 
*.changed to individual MS2's
*
*. Set up Min and max arrays for MS2
*
      CALL MINMAX_MS2_FOR_PAT(IMULT,NAT,MS2TOT,MS2_MINMAX)
*
* We now have min and max of accumulated MS2TOT. Find number of components
*
*. Largest and smallest possible MSvalues
      MAX_MS2 = 0
      DO IAT = 1, NAT
       MAX_MS2 = MAX_MS2 + IMULT(IAT) -1
      END DO
      MIN_MS2 = -MAX_MS2
*. A given MS2 will in the following be stored as MS2-MIN_MS2 + 1
      IF(2*MAX_MS2+1 .GT. NMXCOMP) THEN
        WRITE(6,*) ' NPAS_FOR_PAT, NMXCOMP too small'
        WRITE(6,*) ' Current and required value ', NMXCOMP,2*MAX_MS2+1
        STOP       ' NPAS_FOR_PAT, NMXCOMP too small'
      END IF
      NCOMP = 2*MAX_MS2+1
      IZERO = 0
*. Allowed MS2 values for atom 1
      CALL ISETVC(NCOMP_IN,IZERO,NCOMP)
      DO MS2 = MS2_MINMAX(1,1),MS2_MINMAX(1,2),2
        NCOMP_IN(MS2-MIN_MS2+1) = 1
        ISCR(1,MS2-MS2_MINMAX(1,1)+1) = MS2
        IOFF_IN(MS2-MIN_MS2+1) = MS2-MS2_MINMAX(1,1)+1
      END DO
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' Initial ISCR '
        NINI = (MS2_MINMAX(1,2)-MS2_MINMAX(1,1))/2 + 1
        DO ICOMP = 1, NINI
          WRITE(6,*)  ISCR(1,ICOMP) 
        END DO
      END IF
*. At hand : In ISCR we have strings of IAT atoms, The strings with total MS2
* are stored in ISCR(*,IOFF_IN(MS2-MIN_MS2+1)) to 
* ISCR(*,IOFF_IN(MS2-MIN_MS2+1)-1+NCOMP_IN(MS2-MIN_MS2+1)
*. We can therefore proceed and built the occupations for atom IAT + 1
      DO IAT = 1, NAT-1
*. Atom IAT  => Atom IAT +1 
        IDELTA = IMULT(IAT+1)-1
        LCOMP_AT = 0
        IOFF_OUT(MS2_MINMAX(IAT+1,1)-MIN_MS2+1) = 1
        DO MS2_OUT = MS2_MINMAX(IAT+1,1),MS2_MINMAX(IAT+1,2),2
          LCOMP_AT_MS2 = 0
          DO IIDELTA = - IDELTA, IDELTA,2
            MS2_IN = MS2_OUT - IIDELTA
C?          WRITE(6,*) ' IAT,MS2_OUT,MS2_IN,IIDELTA = ',
C?   &                   IAT,MS2_OUT,MS2_IN,IIDELTA
            IF(MS2_MINMAX(IAT,1).LE.MS2_IN.AND.
     &         MS2_MINMAX(IAT,2).GE.MS2_IN     ) THEN
*. We can combine MS2_OUT for atom IAT + 1 with MS_IN for the previous IAT atoms
               IOFF_IAT_MS2 = IOFF_IN(MS2_IN-MIN_MS2+1)
C?             WRITE(6,*) ' IOFF_IAT_MS2', IOFF_IAT_MS2
               L_IAT_MS2 = NCOMP_IN(MS2_IN-MIN_MS2+1)
               DO ICOMP = 1, NCOMP_IN(MS2_IN-MIN_MS2+1)
                 LCOMP_AT_MS2 =  LCOMP_AT_MS2 + 1
                 LCOMP_AT = LCOMP_AT + 1
                 CALL ICOPVE(ISCR(1,IOFF_IAT_MS2-1+ICOMP),
     &                       IPAS(1,LCOMP_AT),IAT)
                 IPAS(IAT+1,LCOMP_AT) = MS2_OUT
               END DO
            END IF
          END DO
*. We have now constructed all PAS with IAT+1 atoms with acc with given MS2_OUT
          NCOMP_OUT(MS2_OUT-MIN_MS2+1) = LCOMP_AT_MS2
*. And offset for these
          IF(MS2_OUT.LT.MS2_MINMAX(IAT+1,2)) 
     &    IOFF_OUT(MS2_OUT+2-MIN_MS2+1) = LCOMP_AT+1
        END DO
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IPAS after atom (IAT) ', IAT
          DO JPAS = 1, LCOMP_AT
           WRITE(6,'(I6,2x,15(1X,I3))') JPAS, (IPAS(JAT,JPAS),JAT=1,NAT)
          END DO
        END IF
*. We have now constructed the PAS'ses for atom IAT + 1, save
        CALL ICOPVE(NCOMP_OUT,NCOMP_IN,NCOMP)
        CALL ICOPVE(IOFF_OUT,IOFF_IN,NCOMP)
        CALL ICOPVE(IPAS,ISCR,NPAS*NAT)
*
*
      END DO
*. We have now constructed the passes as accumulated MS2's, convert to individual MS2's
      DO JPAS = 1, NPAS
        DO IAT = 2, NAT
          IPAS(IAT,JPAS) = ISCR(IAT,JPAS) - ISCR(IAT-1,JPAS)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' List of Pas''es '
        WRITE(6,*) ' ================'
        DO JPAS = 1, NPAS
          WRITE(6,'(I6,2x,15(1X,I3))') JPAS, (IPAS(IAT,JPAS),IAT=1,NAT)
C         WRITE(6,*) JPAS, (IPAS(JAT,JPAS),JAT=1,NAT)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE CPAS_CSF_SD(CPAS_CSF,CPAS_SD,IWAY,CPACOCC,
     &           CPACMULTS,CPACMS2,NAT,ITRN_AT)
*
* A set of expansion coefficients are given for a PAS. Transform 
* between CSF and SD format
*
* IWAY = 1 : CSF to SD transformation
* IWAY = 1 : SD to CSF transformation 
*
      RETURN
      END
C      SPNCOM_LUCIA(NOPEN,MS2,NDET,IABDET,
C    &                  IABUPP,IFLAG,PSSIGN,IPRCSF)
*. Reform the prototype SD's so they are in lexical order
*. We know have the prototype SD's and the coupling schemes. 
      SUBROUTINE FIND_PROTO_AS(NPAT,IPATSPIN,IPATOCC,NAT,IFLAG,
     &                       NOBAC,NOBPAT,NAS_PROTO,IAS_PROTO)
*
* A set of atomic terms is given in the form of IPAT and IPAC
* Find the various atomic prototype states
* 
* An atomic prototype state is defined by 
* 1 : Number of unpaired electrons 
* 2 : Spinmultiplicity
* 3 : Ms2
* For given number of unpaired electrons and spinmultiplicities, 
* all MS2's are generated.
*
*. IFLAG = 1 => Only number of prototypes
*. IFLAG = 2 => Also the actual prototype
*
*. Jeppe Olsen, November 2005
*
      INCLUDE 'implicit.inc'
*. For getting MXPORB
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER IPATSPIN(NAT,NPAT),IPATOCC(NOBAC,NPAT)
      INTEGER NOBPAT(NAT)
*. Output
       INTEGER IAS_PROTO(3,*)
*. Local scratch 
       DIMENSION IMULTS(MXPORB+1)
*
       NTEST = 100
*. Find min and max number of unpaired electrons for any atom
       DO IPAT = 1, NPAT
         IOB = 0
C        WRITE(6,*) ' IPATOCC(*,IPAT) ',(IPATOCC(I,IPAT),I=1,NOBAC)
         DO IAT = 1, NAT
*. Number of unpaired electrons 
           IOP = 0
           DO IIOB = 1, NOBPAT(IAT)
             IOB = IOB + 1
             IF(IPATOCC(IOB,IPAT).EQ.1) IOP = IOP + 1
           END DO
           IF(IPAT.EQ.1.AND.IAT.EQ.1) THEN
            MIN_OP = IOP
            MAX_OP = IOP
           ELSE 
            MIN_OP = MIN(MIN_OP,IOP)
            MAX_OP = MAX(MAX_OP,IOP)
           END IF
         END DO
       END DO
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' MIN_OP, MAX_OP = ', MIN_OP, MAX_OP
       END IF
* 
       NAS_PROTO = 0
       DO IOP = MIN_OP, MAX_OP
*. Indicate active multiplicities in IMULTS
        IZERO = 0
        CALL ISETVC(IMULTS,IZERO,IOP+1)
        DO IPAT = 1, NPAT
          IOB = 0
          DO IAT = 1, NAT
            IIOP = 0
            DO IIOB = 1, NOBPAT(IAT)
              IOB = IOB + 1
              IF(IPATOCC(IOB,IPAT).EQ.1) IIOP = IIOP + 1
            END DO
            IF(IOP.EQ.IIOP) THEN
              IIMULTS =  IPATSPIN(IAT,IPAT)
              IMULTS(IIMULTS) = IMULTS(IIMULTS) + 1
            END IF
          END DO
        END DO
*. We have now marked the active multiplicities for a given number of 
*. unpaired orbitals, register the results
        DO IIMULTS = 1,  IOP + 1
          IF(IMULTS(IIMULTS).NE.0) THEN
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Active atomic term, IOP, IIMULTS = ',
     &        IOP,IIMULTS
            END IF
*. Include all MS2 components
            IF(IFLAG.EQ.2) THEN
*. Lowest MS2 for this multiplicity 
              MS2 = -(IIMULTS-1)
              DO IADD = 1, IIMULTS
                IAS_PROTO(1,NAS_PROTO+IADD) = IOP
                IAS_PROTO(2,NAS_PROTO+IADD) = IIMULTS
                IAS_PROTO(3,NAS_PROTO+IADD) = MS2
                MS2 = MS2 + 2
              END DO
            END IF
            NAS_PROTO = NAS_PROTO + IIMULTS
C?          WRITE(6,*) ' NAS_PROTO, IIMULTS = ', NAS_PROTO, IIMULTS 
          END IF
        END DO
      END DO
*     ^ End of loop over number of open orbitals
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output from FIND_PROTO '
        WRITE(6,*) ' Number of Prototype atomic states ', NAS_PROTO
        IF(IFLAG.EQ.2) THEN
          WRITE(6,*) ' The prototype atomic states : '
          WRITE(6,*) ' =============================='
          WRITE(6,*)
          WRITE(6,*)  '    I   Nopen   Mult. Ms2 '
          WRITE(6,*)  ' ========================'
          DO IPRO = 1, NAS_PROTO
             WRITE(6,'(I6,2X,I4,2X,I4,2X,I4)') 
     &       IPRO, (IAS_PROTO(I,IPRO),I=1,3)
          END DO
        END IF
*       ^ End if IFLAG = 2
      END IF
*      ^ End if NTEST. GE. 100
      RETURN
      END
      SUBROUTINE GEN_SPINFO_FOR_ALL_PROTO_AS(NPROTO_AS,IPROTO_AS,
     &           NSD_PROTO,NCSF_PROTO,KSD_PROTO,KCSF_PROTO,
     &           KCSDCSF_PROTO)
*
* Set up proto types information for all prototype atomic states
*
*. Jeppe Olsen, November 2005
*
*. Note : Memory is allocated in this routine and returned
*.
*. Input
* =======
*  NPROTO_AS : Number of prototype atomic states
*  IPROTO_AS : The actual prototype atomic states
*
*. Output
* =======
*  NSD_PROTO : Number of SD's per prototype atomic state
*  NCSF_PROTO : Number of CSF's per prototype atomic state
*  KSD_PROTO : Pointers to prototype determinants in work
*  KCSF_PROTO : Pointers to prototype CSF's in work
*  KCSDCSF_PROTO : Pointers to prototype SD-CSD transformation in work
*
*. Note : the 
*
      INCLUDE 'wrkspc.inc'
*. Input
      INTEGER IPROTO_AS(3,NPROTO_AS)
*. Output
      INTEGER NSD_PROTO(NPROTO_AS),NCSF_PROTO(NPROTO_AS)
      INTEGER KSD_PROTO(NPROTO_AS),KCSF_PROTO(NPROTO_AS)
      INTEGER KCSDCSF_PROTO(NPROTO_AS)
      DO JPROTO_AS = 1, NPROTO_AS
        WRITE(6,*) ' JPROTO_AS =', JPROTO_AS
        CALL MEMCHK2('CHECK1')
*
* ==========================
*. Number of SD's and CSF's
* ==========================
*
        NOPEN = IPROTO_AS(1,JPROTO_AS)
        MULTS = IPROTO_AS(2,JPROTO_AS)
        MS2   = IPROTO_AS(3,JPROTO_AS)
*. Number of alpha and beta electrons
        NAEL = (NOPEN + MS2 ) / 2
        NBEL = (NOPEN - MS2 ) / 2
        IF(NAEL+NBEL .EQ. NOPEN .AND. NAEL-NBEL .EQ. MS2 .AND.
     &            NAEL .GE. 0 .AND. NBEL .GE. 0) THEN
*. Allowed case
          NSD = IBION(NOPEN,NAEL)
          IF(NOPEN .GE. MULTS-1) THEN
            NCSF = IWEYLF(NOPEN,MULTS)
          ELSE
            NCSF = 0
          END IF
        ELSE
*. Case without SD's and CSF's
          NSD = 0
          NCSF = 0
        END IF
        CALL MEMCHK2('CHECK2')
        NSD_PROTO(JPROTO_AS) = NSD
        CALL MEMCHK2('CHEC2a')
        NCSF_PROTO(JPROTO_AS) = NCSF
        CALL MEMCHK2('CHECK3')
* 
* ===========================================================
* Allocate space for SD's, CSF's and transformation matrices
* ===========================================================
*
        LEN_SD = NOPEN*NSD
        LEN_CSF = NOPEN*NCSF
        LEN_CTRA = NSD*NCSF
*
        CALL MEMMAN(KSD_PROTO(JPROTO_AS),LEN_SD,'ADDL  ',2,'SD_PRO')
        CALL MEMMAN(KCSF_PROTO(JPROTO_AS),LEN_CSF,'ADDL  ',2,'CS_PRO')
        CALL MEMMAN(KCSDCSF_PROTO(JPROTO_AS),LEN_CTRA,'ADDL  ',2,
     &              'CT_PRO')
*. And the actual transformation matrices and prototype determinants
        CALL MEMCHK2('BE_GEN')
        CALL GEN_SPINFO_FOR_PROTO(MS2,MULTS,NOPEN,
     &       WORK(KSD_PROTO(JPROTO_AS)),WORK(KCSF_PROTO(JPROTO_AS)),
     &       WORK(KCSDCSF_PROTO(JPROTO_AS)) )
        CALL MEMCHK2('AF_GEN')
      END DO
      RETURN
      END
      SUBROUTINE GEN_SPINFO_FOR_PROTO(MS2,MULTS,NOPEN,ISD,ICSF,CSDCSF)
*
* A prototype is defined by M2, MULTS, NOPEN.
*.Construct the prototype determinants, spin-couplings and 
*.csf-sd transformations.
*
*. Output is
*
* ISD : List of prototype determinants
* ICSD : List of prototype CSF couplings
* CSDCSF : Transformation matrix as C(NSD,NCSF) matrix
*
* Jeppe Olsen, November 2005
*
* Note : The prototype determinants are not generated in lexical order
*        and are not resorted. 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cprnt.inc'
*. Output
      INTEGER ISD(NOPEN,*)
      INTEGER ICSF(NOPEN,*)
      DIMENSION CSDCSF(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GEN_SP')
*
      NTEST = 100 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' CSF-SD information will be generated for '
        WRITE(6,*) '      Number of unpaired orbitals ', NOPEN
        WRITE(6,*) '      Spin-multiplicity           ', MULTS
        WRITE(6,*) '      MS2 * 2                     ', MS2
      END IF
*. Number of alpha and beta electrons 
      NAEL = (NOPEN + MS2 ) / 2
      NBEL = (NOPEN - MS2 ) / 2
      IF(NAEL+NBEL .EQ. NOPEN .AND. NAEL-NBEL .EQ. MS2 .AND.
     &            NAEL .GE. 0 .AND. NBEL .GE. 0) THEN
*. Allowed case 
        NSD = IBION(NOPEN,NAEL)
        IF(NOPEN .GE. MULTS-1) THEN
          NCSF = IWEYLF(NOPEN,MULTS)
        ELSE
          NCSF = 0
        END IF
      ELSE
*. Case without SD's and CSF's
        NSD = 0
        NCSF = 0
      END IF
*
      IF(NSD.NE.0.AND.NCSF.NE.0) THEN
        ZERO = 0.0D0
*. Prototype determinants
        IFLAG = 1
        CALL SPNCOM_LUCIA(NOPEN,MS2,NNDET,ISD,ICSF,IFLAG,ZERO,
     &                    IPRCSF)
*. And proto-type couplings
        IFLAG = 3
        CALL SPNCOM_LUCIA(NOPEN,MULTS-1,NNDET,ISD,ICSF,IFLAG,ZERO,
     &                    IPRCSF)
*. And the CSF-SD transformation matrices
        LSCR = (NSD+1)*NOPEN
        WRITE(6,*) ' LSCR, NSD, NOPEN = ', LSCR,NSD, NOPEN
        CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'SCR   ')
        CALL CSFDET_LUCIA(NOPEN,ISD,NSD,ICSF,NCSF,CSDCSF,
     &                    WORK(KLSCR),ZERO,IPRCSF)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GEN_SP')
      RETURN
      END
      SUBROUTINE SIGDEN_PRWF(ISD,SIGMA,C,MS2TOT,NPATL,NPATR,
     &           IPATSPINL,IPATSPINR,
     &           IPATOCCL,IPATOCCR,NPAS_PER_PATL,NPAS_PER_PATR,NAT)
*
* Calculate one- and two-body density matrices or 
* sigma-vector for a product wave-function using 
* partial expansion to Slater-determiant basis
* 
*. Jeppe Olsen, December 2005
*
* ISD = 1 => sigmavector, Sigma = H C
* ISD = 2 => densitymatrices, <Sigma!E!sigma>, <sigma!e!sigma>
*
*
*. I am here using NACOB as the number of orbitals in the various atoms, 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'prdwvf.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*· Input : coefficients of the PASses
      DIMENSION C(*)
*. CI : Output, Densi : Input
      DIMENSION SIGMA(*)
*. The various PATs in expansion 
      INTEGER IPATSPINL(NAT,NPATL),IPATSPINR(NAT,NPATR)
      INTEGER IPATOCCL(NACOB,NPATL),IPATOCCR(NACOB,NPATR)
      INTEGER NPAS_PER_PATL(NPATL),NPAS_PER_PATR(NPATR)
*. Local scratch
      INTEGER IDIFFAT_PAT(100), IDIFFAT_PAS(100)
      INTEGER ILMS2(MXPORB),IRMS2(MXPORB)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,MARK,'SG_PRW')
*. For PAS expansion of given PAT
      CALL MEMMAN(KLPASL,MXLENPAT,'ADDL  ',2,'PAS_L ')
      CALL MEMMAN(KLPASR,MXLENPAT,'ADDL  ',2,'PAS_R ')
      CALL MEMMAN(KLPASSCR,MXLENPAT,'ADDL  ',2,'PASSCR')
*
      NTEST = 1000
*
      NAT = NGAS
*. Orbital space of first atom is set to 1
      IBGAS_AT = 1
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',ADDL,' SIDEPR')
*. Total number of PATS in L and R expansions
      NPASL = IELSUM(NPAS_PER_PATL,NPATL)
      NPASR = IELSUM(NPAS_PER_PART,NPATR)
      IF(NTEST.GE.100) THEN
       WRITE(6,*)' Number of PAS''ses in N and R expansions',NPATL,NPATR
      END IF
*. Zero sigma-vector for sigmavectorgeneration
      IF(ISIGDE.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(SIGMA,ZERO,NPASL)
      END IF
*
      DO IPATL = 1, NPATL
*. Offset to PASses for given PAT
        IF(IPATL.EQ.1) THEN
          IOFFL = 1
        ELSE
          IOFFL = IOFFL + NPAS_PER_PATL(IPATL-1)
        END IF
        DO IPATR = 1, NPATR
          IF(IPATR.EQ.1) THEN
            IOFFR = 1
          ELSE
            IOFFR = IOFFR + NPAS_PER_PATR(IPATR-1)
          END IF
*. Check if the two PAT's interact with atmost two-body  operator
          CALL CONNECT_PAT(IPATOCCL(1,IPATL),IPATOCCR(1,IPATR),
     &         IPATSPINL(1,IPATL),IPATOCCR(1,IPATR),ICONNECT,NAT,
     /         NOBPT(IBGAS_AT),I1CREA,I2CREA,I1ANNI,I2ANNI,NDIFF,
     &         NDIFFAT,IDIFFAT_PAT)
C     SUBROUTINE CONNECT_PAT(ILPAT,IRPAT,ILMULT,IRMULT,ICONNECT,NAT,
C    %           NOBPT,I1CREA, I2CREA, I1ANNI,I2ANNI,NDIFF,NDIFFAT,
C    &           IDIFFAT)
          IF(NDIFF.LE.4) THEN
*. 
*. We have that the two PATs interact, generate and loop over  the PAS'es of 
*. this PAT
            NPASL = NPAS_PER_PATL(IPATL)
            NPASR = NPAS_PER_PATR(IPATR)
C     GENPAS_FOR_PAT(IMULT,NAT,MS2TOT,NPAS,IPAS,ISCR)
*. Generate Passes
            CALL GENPAS_FOR_PAT(IPATSPINL(1,IPATL),NAT,MS2TOT,NPASL,
     &           IPATOCCL(1,IPATL),WORK(KLPASL),WORK(KLPASSCR))
            CALL GENPAS_FOR_PAT(IPATSPINR(1,IPATR),NAT,MS2TOT,NPASR,
     &           IPATOCCR(1,IPATR),WORK(KLPASR),WORK(KLPASSCR))
         
            DO IPASL = 1, NPASL  
             CALL ICOPVE2(WORK(KLPASL),(IPASL-1)*NAT+1,NAT,ILMS2)
C                  ICOPVE2(IIN,IOFF,NDIM,IOUT)
             DO IPASR = 1, NPASR
              CALL ICOPVE2(WORK(KLPASR),(IPASR-1)*NAT+1,NAT,IRMS2)
*. Check number of atoms for which the PAS'ses differs
              CALL CONNECT_PAS(ILMS2,IRMS2,IPATSPINL(1,IPATL),
     &        IPATSPINR(1,IPATR),NAT,NDIFFAT_PAT,IDIFFAT_PAT,
     &        NDIFFAT_PAS,IDIFFAT_PAS)
C     CONNECT_PAS(ILMS2,IRMS2,ILMULT,IRMULT,NAT,
C    %           NDIFFAT_PAT,IDIFFAT_PAT,NDIFFAT_PAS,IDIFFAT_PAS)
*
* In the following, matrix elements will be evaluated of we have 
* differences in PATs for atmost four orbitals and differences 
* in PASs at atmost four atoms. This screening of vanishing
* terms could be improved 
              IF(NDIFFAT_PAS.LE.4) THEN
*. Transform the PAS'ses from CSF's to SD's

*
*. 
         IF(NDIFF.EQ.4) THEN
*. Differences in the occupation of four orbitals. Only the part 
* of the Hamiltonian that connects these four orbitals and the corresponding 
* atoms gives nonvanishing contribution. 
*
*. The PAS'ses for the various states
C             GENPAS_FOR_PAT(IMULT,NAT,MS2TOT,NPAS,IPAS,ISCR)
         
*
         END IF
         END IF
* Expand the expansions in the active atoms to the SD basis.
         END DO
         END DO
        END IF
        END DO
        END DO
*
       RETURN
       END
      SUBROUTINE CONNECT_PAS(ILMS2,IRMS2,ILMULT,IRMULT,NAT,
     &           NDIFFAT_PAT,IDIFFAT_PAT,NDIFFAT_PAS,IDIFFAT_PAS)
*
* The Passes are given with MS2 projections and totals spins ILMS2, IRMS2
* and ILMULT, IRMULT. The corresponding PAT's have already been examined 
* and it has been determined that the PAT's differ in NDIFFAT_PAT atoms 
* IDIFFAT_PAT.
* Compare the PAS'ses and find the number of atoms for which MS2 or S
* differs
* 
* Jeppe Olsen, December 2005
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ILMS2(NAT),IRMS2(NAT),ILMULT(NAT),IRMULT(NAT)
      INTEGER IDIFFAT_PAT(*)
*. Output
      INTEGER IDIFFAT_PAS(*)
*
      NDIFFAT_PAS = 0
      DO IAT = 1, NAT
*. Is atom AT among the atoms where PAT's differ ?
        IPAT_DIFFER = 0
        DO JAT = 1, NDIFFAT_PAT
          IF(IDIFFAT_PAT(JAT).EQ.JAT) IPAT_DIFFER = 1
        END DO
        IF(IPAT_DIFFER.EQ.1) THEN
*. PAS'es differ as mother PAT's differ
          NDIFFAT_PAS = NDIFFAT_PAS + 1
          IDIFFAT_PAT(NIDFFAT_PAS) = IAT
        ELSE
*. PATS identical, check S and MS2
          IF(ILMS2(IAT).NE.IRMS2(IAT).OR.
     &       ILMULT(IAT).NE.IRMULT(IAT)  ) THEN
               NDIFFAT_PAS = NDIFFAT_PAS + 1
               IDIFFAT_PAT(NIDFFAT_PAS) = IAT
          END IF
        END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Connections between two PAS '
        WRITE(6,*) ' ================================='
        WRITE(6,*)
        WRITE(6,*) ' Info on PAS (input)'
        WRITE(6,*) ' MS2 of L and R '
        CALL IWRTMA(ILMS2,1,NAT,1,NAT)
        CALL IWRTMA(IRMS2,1,NAT,1,NAT)
        WRITE(6,*) ' MULTS of L and R '
        CALL IWRTMA(ILMULTS,1,NAT,1,NAT)
        CALL IWRTMA(IRMULTS,1,NAT,1,NAT)
        WRITE(6,*) ' PATs differs at atoms '
        CALL IWRTMA(IDIFFAT_PAT,1,NDIFFAT_PAT,1,NDIFFAT_PAT)
        WRITE(6,*)
        WRITE(6,*) ' Number of atoms for which the two PAS differs ',
     &  NDIFFAT_PAS
        WRITE(6,*) ' atoms for which the PATS differs '
        CALL IWRTMA(IDIFFAT_PAS,1,NDIFFAT_PAS,1,NDIFFAT_PAS)
      END IF
*
      RETURN
      END
       SUBROUTINE EXPND_PAS_PRWF_CSFSD(
     &            CCSF,CSD,NAT,IMULT,IMS2,IOCC,NAT_EXPN,IAT_EXPND,
     &            NOBPAT, IWAY)
*
* A PAS is given by IOCC, IMULT,IMS2, Transform between CSF and SD 
* forms. Transformations are only performed for the the NAT_EXPND
* atoms given in IAT_EXPND. A sign, IAPSIGN, giving the signchange 
* required to put the wf for the expanded atoms in front is also 
* obtained. 
*
* IWAY = 1 : CSF to SD transformation
* IWAY = 1 : SD to CSF transformation 
*
*. Jeppe Olsen, December 2005
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'prdwvf.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Input/ output
      DIMENSION CSD(*)
*. Input : PAS
      INTEGER IOCC(NACOB),IMULT(NAT),IMS2(NAT)
*. Atoms to be expanded
      INTEGER IAT_EXPND
      RETURN
      END
c $Id$
