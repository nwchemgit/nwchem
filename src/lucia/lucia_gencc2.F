* GENCC production codes (whatever that means)
      SUBROUTINE CI_TO_CC_REFRM(LUCC,LUCI,ISPC,ISM)
*
* A CI vector is given as the only vector in LUCI
* Rewrite this vector to a set of Coupled CLuster amplitudes so
*
* Exp T |Ref> = CI
*
* Jeppe Olsen, April 14, early in the morning 
*
* Reference space (CI space 1 ) is assumed to be a single det
*
* Output CC coefficients are put on FILE LUCC in current form 
* of CC coefficients.
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cprnt.inc'

      CHARACTER*6 CCTYPE
*
      NTEST = 00
*
      IDUM = 0
      CCTYPE(1:6) = 'GEN_CC'
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CICCRF')
      CALL LUCIAQENTER('CI_CC ')
*
      WRITE(6,*) ' CI => CC transformation of coefficients '
*. Space for CI behind the curtain
      CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
*
* 1 : Find coefficient of reference det
*
      IREFSPC = 1
      LBLK = -1
C     EXPCIV(ISM,ISPCIN,LUIN,ISPCUT,LUUT,LBLK,
C    &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
      CALL EXPCIV(ISM,ISPC,LUCI,IREFSPC,LUSC1,LBLK,
     &            LUHC,1,0,IDC,NTEST)
      WRITE(6,*) ' LBLK after EXPCIV = ', LBLK
*. Number of records in reference space
      IF(IDC.EQ.1.OR.ISM.EQ.1) THEN
        NREC = NSMST
      ELSE
        NREC = NSMST/2
      END IF
*. Read reference coefficient
C          FRMDSCN(VEC,NREC,LBLK,LU)
      CALL REWINO(LUSC1)
      CALL FRMDSCN(CREF,NREC,LBLK,LUSC1)
      WRITE(6,*) ' CREF = ', CREF
      WRITE(6,*) ' LBLK after FRMDSCN = ', LBLK
*
* 2 : Normalize CI vector so reference coef is one 
*
      IF(CREF.EQ.0.0D0) THEN 
        WRITE(6,*) ' CI_TO_CC_RF: Problems, norm of ref coef = 0'
        STOP        'CI_TO_CC_RF: Problems, norm of ref coef = 0'
      ELSE
        FACTOR = 1.0D0/CREF
*. And the scaling, result on LUHC
        IREW = 1
C            SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
        WRITE(6,*) ' Before call to SCLVCD '
        WRITE(6,*) ' LUCI, LUHC = ', LUCI, LUHC
        WRITE(6,*) ' KVEC1, LBLK = ', KVEC1, LBLK
        CALL SCLVCD(LUCI,LUHC,FACTOR,WORK(KVEC1),IREW,LBLK)
        WRITE(6,*) ' After call to SCLVCD '
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Scaled CI vector '
          CALL WRTVCD(WORK(KVEC1),LUHC,1,LBLK)
        END IF
      END IF
*. So now : Complete vector with reference coef = 1 on LUHC
*
* 3. Information about CC space 
*
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL = NAEL + NBEL
*
      IREFSPC = 1
*. Find the type of reference state 
      CALL CC_AC_SPACES(IREFSPC,IREFTYP)
*. Number of active orbital spaces
      NACT_SPC = 0
      IACT_SPC = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) THEN
          NACT_SPC = NACT_SPC + 1
          IACT_SPC = IGAS
        END IF
      END DO
*. Info on active-active excitation types
      CALL ACAC_EXC_TYP(IAAEXC_TYP,MX_AAEXC,IPRCC)
*. Number of occupation classes for actual space
      CALL OCCLSE(1,NOCCLS,IOCCLS,NEL,ISPC,0,0,NOBPT)
*. And the occupation classes of actual space
      CALL MEMMAN(KLOCCLS,NOCCLS*NEL,'ADDL  ',1,'OCCLS ')
      CALL OCCLSE(2,NOCCLS,WORK(KLOCCLS),NEL,ISPC,0,0,NOBPT)
*. Number of occupation classes for reference space
      CALL OCCLSE(1,NOCCLS_REF,IOCCLS,NEL,IREFSPC,0,0,NOBPT)
      IF(NOCCLS_REF.GT.1) THEN
        WRITE(6,*) ' Problem in general CC '
        WRITE(6,*) 
     &  ' Reference space is not a single occupation space'
        STOP 
     &  ' Reference space is not a single occupation space'
      END IF
*. and the occupation classes of reference space
      CALL MEMMAN(KLOCCLS_REF,NGAS,'ADDL  ',1,'OCC_RF')
      CALL OCCLSE(2,NOCCLS_REF,WORK(KLOCCLS_REF),NEL,IREFSPC,0,0,NOBPT)
*. Excitation type => Original occupation class
*. 
*. Number of excitation types 
      IFLAG = 1
      IDUM = 1
      CALL TP_OBEX2(NOCCLS,NEL,NGAS,WORK(IDUM),
     &             WORK(IDUM),WORK(IDUM),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(IDUM),MX_AAEXC,IFLAG,
     &             I_OOCC,NOBEX_TP,NOAAEX,IPRCC)
C?    WRITE(6,*) ' NOBEX_TP,MX_EXC_LEVEL = ', NOBEX_TP,MX_EXC_LEVEL
*. And the actual orbital excitaions
      NOBEX_TPE = NOBEX_TP + 1
      CALL MEMMAN(KLCOBEX_TP,NOBEX_TPE,'ADDL  ',1,'LCOBEX')
      CALL MEMMAN(KLAOBEX_TP,NOBEX_TPE,'ADDL  ',1,'LAOBEX')
      CALL MEMMAN(KOBEX_TP ,NOBEX_TPE*2*NGAS,'ADDL  ',1,'IOBE_X')
*. Excitation type => Original occupation class
      CALL MEMMAN(KEX_TO_OC,NOBEX_TPE,'ADDL  ',1,'EX__OC')
      IFLAG = 0
      CALL TP_OBEX2(NOCCLS,NEL,NGAS,WORK(KOBEX_TP),
     &             WORK(KLCOBEX_TP),WORK(KLAOBEX_TP),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(KEX_TO_OC),MX_AAEXC,IFLAG,
     &             I_OOCC,NOBEX_TP,NOAAEX,IPRCC)
*. Spinorbital excitations
*. Spin combinations of CC excitations : Currently we assume that 
*. The T-operator is a singlet, can 'easily' be changed 
      IF(PSSIGN.EQ.0.0D0) THEN
        MSCOMB_CC = 0
      ELSE IF(PSSIGN.EQ.1.0D0) THEN 
        MSCOMB_CC = 1  
      END IF
      MSCOMB_CC = 0
*. Number of spin-orbital excitations
      CALL OBEX_TO_SPOBEX(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TP,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,IPRCC,IDUMMY,
     &     MXSPOX,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),
     &     NAEL,NBEL,IREFSPC)
*. And the actual spinorbital excitation operators 
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. First SOX of given OX ( including zero operator )
      CALL MEMMAN(KIBSOX_FOR_OX,NOBEX_TP+1,'ADDL  ',1,'IBSOXF')
*. Number of SOX's for given OX
      CALL MEMMAN(KNSOX_FOR_OX,NOBEX_TP+1,'ADDL  ',1,'IBSOXF')
*. SOX for given OX
      CALL MEMMAN(KISOX_FOR_OX,NSPOBEX_TP+1,'ADDL  ',1,'IBSOXF')

*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,
     &     IPRCC,WORK(KLSOX_TO_OX),MXSPOX,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPC)
*. Alpha- and beta-excitations constituting the spinorbital excitations
*. Number 
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      CALL MEMMAN(KLAOBEX,LENA,'ADDL  ',2,'IAOBEX')
      CALL MEMMAN(KLBOBEX,LENB,'ADDL  ',2,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     0,NAOBEX_TP,NBOBEX_TP,WORK(KLAOBEX),WORK(KLBOBEX))
*. Max dimensions of CCOP !KSTR> = !ISTR> maps
*. For alpha excitations
      IATP = 1
      IOCTPA = IBSPGPFTP(IATP)
      NOCTPA = NSPGPFTP(IATP)
      CALL LEN_GENOP_STR_MAP(
     &     NAOBEX_TP,WORK(KLAOBEX),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IBTP = 2
      IOCTPB = IBSPGPFTP(IBTP)
      NOCTPB = NSPGPFTP(IBTP)
      CALL LEN_GENOP_STR_MAP(
     &     NBOBEX_TP,WORK(KLBOBEX),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1


* Dimension of spinorbital excitation operators
      ITOTSM = 1
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TP,'ADDL  ',1,'SPOBAC')
*
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ITOTSM,
     &          MX_ST_TSOSO_MX,MX_ST_TSOSO_BLK_MX,MX_TBLK_MX,
     &          WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &          MSCOMB_CC,MX_SBSTR,
     &          WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &          NTCONF,IPRCC)
     
      N_CC_AMP = LEN_T_VEC
      WRITE(6,*) 'N_CC_AMP = ', N_CC_AMP
*. Allocate three CC vectors 
      CALL MEMMAN(KCCF,N_CC_AMP,'ADDL  ',2,'CCF   ')
      CALL MEMMAN(KCC1,N_CC_AMP,'ADDL  ',2,'CC1   ')
      CALL MEMMAN(KCC2,N_CC_AMP,'ADDL  ',2,'CC2   ')
*
      CALL MEMMAN(KLLCC,NSPOBEX_TP,'ADDL  ',1,'LCC   ')
      CALL MEMMAN(KLICC,NSPOBEX_TP,'ADDL  ',1,'ICC   ')
      CALL MEMMAN(KLJCC,NSPOBEX_TP,'ADDL  ',1,'JCC   ')
*
* Now the rest of the show goes as 
* Vector |LUHC> starts as complete CI with coef of refeence = 1
*
* Loop over excitation levels IEXC
*  Reform CI vector !LUHC> to CC form
*  Extract coefficients of excitation level IEXC,
*  These are the CC coeffcients for this level
*  Calculate Exp(-T(iexc))!LUHC> and store on LUHC 
* End of loop over excitation levels.
* 
      ZERO = 0.0D0
      CALL SETVEC(WORK(KCCF),ZERO,N_CC_AMP)
      DO IEXC = 1, MX_EXC_LEVEL 
       IF(NTEST.GE.100) WRITE(6,*) ' Excitation level = ', IEXC
*. Reform current CI coefficient to CC form, and store in WORK(KCC1)
       IREW = 1
       I_DO_CC_INFO = 0
       CALL CC_CI_REORD(WORK(KCC1),LUHC,2,ISPC,ISM,IREW,I_DO_CC_INFO)  
*. Copy coefficients of excitation level to vector containing 
*. final CC amplitudes
*. Spinorbital excitation types corresponding to this excitation level
C     GET_SPOBTP_FOR_EXC_LEVEL(ILEVEL,ILEVEL_FOR_EXTP,
C     &           NEXTP,NEXTP_AC,IEXTP_AC,ISOX_TO_OX)
      CALL GET_SPOBTP_FOR_EXC_LEVEL(IEXC,WORK(KLCOBEX_TP),NSPOBEX_TP,
     &     NEXTP_AC,WORK(KLJCC),WORK(KLSOX_TO_OX))
*^ The active spinorbital excitation types are stored in WORK(KLJCC)
*. first gathering from KCC1 to KCC2
C     SCAGAT_CCVEC(CC_CMP,CC_EXP,ISG,NEXTP_SG,IEXTP_SG,
C    &           IBEXTP,LEXTP,LEXTP_SG)
       CALL  SCAGAT_CCVEC(WORK(KCC2),WORK(KCC1),2,NEXTP_AC,
     &              WORK(KLJCC),WORK(KLIBSOBEX),WORK(KLLSOBEX),
     &              WORK(KLLCC) )
*. Then scatter from KCC1 to KCCF 
       CALL  SCAGAT_CCVEC(WORK(KCC2),WORK(KCCF),1,NEXTP_AC,
     &              WORK(KLJCC),WORK(KLIBSOBEX),WORK(KLLSOBEX),
     &              WORK(KLLCC) )
*
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Updated list of final CC coefficients : '  
         CALL WRTMAT(WORK(KCCF),1,N_CC_AMP,1,N_CC_AMP)
       END IF
*. calculate Exp(-T(iexc)|LUHC> 
*. Make only excitations with excitation level IEXC level
       IZERO = 0
       CALL ISETVC(WORK(KLSPOBEX_AC),IZERO,NSPOBEX_TP)
       IONE = 1
       CALL ISCASET(WORK(KLSPOBEX_AC),IONE,WORK(KLJCC),NEXTP_AC)
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' List of active spobex fresh from ISCASET '
         CALL IWRTMA(WORK(KLSPOBEX_AC),1,NSPOBEX_TP,1,NSPOBEX_TP)
       END IF
*.  Exp(-t) !LUHC on LUSC35
       MX_TERM = 100
       ICC_EXC = 1
       XCONV = 1.0D-20
       CALL COPVEC(WORK(KCCF),WORK(KCC1),N_CC_AMP)
       ONEM = -1.0D0
       CALL SCALVE(WORK(KCC1),ONEM,N_CC_AMP)
       CALL EXPT_REF(LUHC,LUSC35,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &               WORK(KVEC1),WORK(KVEC2),CCTYPE)
*. And transfer to LUHC
       CALL COPVCD(LUSC35,LUHC,WORK(KVEC1),1,LBLK)
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Updated CI vector on LUHC '
         CALL WRTVCD(WORK(KVEC1),LUHC,1,LBLK)
       END IF
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CC coefficents obtained from CI coefficients'
       CALL WRT_CC_VEC2(WORK(KCCF),IDUMMY,CCTYPE)
      END IF
*. Dump to LU_CI_TO_CC
C     CALL REWINO(LU_CC_FROM_CI)
C     WRITE(LU_CC_FROM_CI,'(I9)') N_CC_AMP
C     DO I = 1, N_CC_AMP
C       WRITE(LU_CC_FROM_CI,'(E25.15)') WORK(KCCF-1+I)       
C     END DO 
*. Dump to LUCCAMP
      CALL REWINO(LU_CC_FROM_CI)
      I_FORMATTED = 0
      IF(I_FORMATTED.EQ.1) THEN
        WRITE(LU_CC_FROM_CI,'(I9)') N_CC_AMP
        DO I = 1, N_CC_AMP
        WRITE(LU_CC_FROM_CI,'(E25.15)') WORK(KCCF-1+I)       
        END DO 
      ELSE
        WRITE(LU_CC_FROM_CI) N_CC_AMP
        WRITE(LU_CC_FROM_CI) (WORK(KCCF-1+I),I=1, N_CC_AMP)       
      END IF
      CALL REWINO(LU_CC_FROM_CI)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CICCRF')
      CALL LUCIAQEXIT('CI_CC ')
      RETURN
      END
      SUBROUTINE ISCASET(IARRAY,IVAL,ISCA,NSCA)
*
* IARRAY(ISCA(I)) = IVAL
*
* Jeppe Olsen, Aril 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER ISCA(NSCA)
*. Output
      INTEGER IARRAY(*)
*
      DO I = 1, NSCA
       IARRAY(ISCA(I)) = IVAL
      END DO
*
      RETURN
      END
      SUBROUTINE CC_CI_REORD(CCVEC,LUCI,IWAY,ISPC,ISM,IREW,
     &           I_DO_CC_INFO)
*
* Convert between CI and CC organizations of coupled cluster 
* coefficients. Note that in this routine, the coefficients 
* are only reordered, no exponentations are involved here.
*
* Reference state is assumed to be a single Slaterdeterminant
*
* Input
* =====
*
* CCVEC : Amplitudes organized as a CC vector
* LUCI : File containing initial/final CI coefficients
* IWAY  = 1 => CC to CI
*       = 2 => CI to CC
* ISPC : Space of expansions
* ISM  : Symmetry of expansions
*
* Note : In core version in line with current assumption 
* that all coefs can be stored in core
*
* CI coefficients are initially/finally on file LUCI
* but are in the routine store in a single array
*
* Jeppe Olsen, Magistratsvaegen 37, March 25  2000
*              - in the kitchen, smelling Dittes cake and 
*                listening to Stones, Sticky  Fingers
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
*
      DIMENSION CCVEC(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CC_CI_')
*
      NTEST = 00
      LBLK = -1
      IF(IREW.EQ.1) THEN
        CALL REWINO(LUCI)
      END IF
*
* Info on occupation classes in expansion 
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL = NAEL + NBEL
*
      ICSPC = ISPC
      ISSPC = ISPC
*
      IREFSPC = 1
*
*. Number of occupation classes in CI and CC  expansion
      CALL OCCLSE(1,NOCCLS,IOCCLS,NEL,ISPC,0,0,NOBPT)
COLD  CALL OCCLS (1,NOCCLS,IOCCLS,NEL,NGAS,
COLD &           IGSOCCX(1,1,ISPC),IGSOCCX(1,2,ISPC),
COLD &           0,0,NOBPT)
*. And the occupation classes
      CALL MEMMAN(KLOCCLS,NOCCLS*NEL,'ADDL  ',1,'OCCLS ')
      CALL OCCLSE(2,NOCCLS,WORK(KLOCCLS),NEL,ISPC,0,0,NOBPT)
COLD  CALL OCCLS (2,NOCCLS,WORK(KLOCCLS),NEL,NGAS,
COLD &           IGSOCCX(1,1,ISPC),IGSOCCX(1,2,ISPC),
COLD &           0,0,NOBPT)
*
      IF(I_DO_CC_INFO.EQ.1) THEN
*
*              ==========================
*              Information about CC space 
*              ==========================
*
*
* Find the type of reference state 
*
      CALL CC_AC_SPACES(IREFSPC,IREFTYP)
*. Number of active orbital spaces
      NACT_SPC = 0
      IACT_SPC = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) THEN
          NACT_SPC = NACT_SPC + 1
          IACT_SPC = IGAS
        END IF
      END DO
*. Info on active-active excitation types
      CALL ACAC_EXC_TYP(IAAEXC_TYP,MX_AAEXC,IPRCC)
*. Number of occupation classes for reference space
      IREFSPC = 1
      CALL OCCLSE(1,NOCCLS_REF,IOCCLS,NEL,IREFSPC,0,0,NOBPT)
      IF(NOCCLS_REF.GT.1) THEN
        WRITE(6,*) ' Problem in general CC '
        WRITE(6,*) 
     &  ' Reference space is not a single occupation space'
        STOP 
     &  ' Reference space is not a single occupation space'
      END IF
*. and the reference occupation space 
      CALL MEMMAN(KLOCCLS_REF,NGAS,'ADDL  ',1,'OCC_RF')
      CALL OCCLSE(2,NOCCLS_REF,WORK(KLOCCLS_REF),NEL,IREFSPC,0,0,NOBPT)
COLD  CALL OCCLS (2,NOCCLS_REF,WORK(KLOCCLS_REF),NEL,NGAS,
COLD &           IGSOCCX(1,1,IREFSPC),IGSOCCX(1,2,IREFSPC),
COLD &           0,0,NOBPT)
*. Number of excitation types 
      IFLAG = 1
      IDUM = 1
      CALL TP_OBEX2(NOCCLS,NEL,NGAS,WORK(IDUM),
     &             WORK(IDUM),WORK(IDUM),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(IDUM),MX_AAEXC,IFLAG,
     &             I_OOCC,NOBEX_TP,NOAAEX,IPRCC)
      WRITE(6,*) ' NOBEX_TP,MX_EXC_LEVEL = ', NOBEX_TP,MX_EXC_LEVEL
      CALL MEMMAN(KLCOBEX_TP,NOBEX_TP,'ADDL  ',1,'LCOBEX')
      CALL MEMMAN(KLAOBEX_TP,NOBEX_TP,'ADDL  ',1,'LAOBEX')
      CALL MEMMAN(KOBEX_TP ,NOBEX_TP*2*NGAS,'ADDL  ',1,'IOBE_X')
*. Excitation type => Original occupation class
      CALL MEMMAN(KEX_TO_OC,NOBEX_TP,'ADDL  ',1,'EX__OC')
      IFLAG = 0
      CALL TP_OBEX2(NOCCLS,NEL,NGAS,WORK(KOBEX_TP),
     &             WORK(KLCOBEX_TP),WORK(KLAOBEX_TP),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(KEX_TO_OC),MX_AAEXC,IFLAG,
     &             I_OOCC,NOBEX_TP,NOAAEX,IPRCC)
*. Spinorbital excitation types
*. Spin combinations of CC excitations : Currently we assume that 
*. The T-operator is a singlet, can 'easily' be changed 
      IF(PSSIGN.EQ.0.0D0) THEN
        MSCOMB_CC = 0
      ELSE IF(PSSIGN.EQ.1.0D0) THEN 
        MSCOMB_CC = 1  
      END IF
      MSCOMB_CC = 0
*. Number of spinorbital excitation types 
      CALL OBEX_TO_SPOBEX(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TP,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,IPRCC,IDUMMY,
     &     NAEL,NBEL)
*. And the actual spinorbital excitation types
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,
     &     IPRCC,WORK(KLSOX_TO_OX),NAEL,NBEL)
* Dimension of spinorbital excitation operators
      ITOTSM = 1
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TP,'ADDL  ',1,'SPOBAC')
*
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ITOTSM,
     &              MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &              WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &              MSCOMB_CC,MX_SBSTR,IPRCC)
      N_CC_AMP = LEN_T_VEC
      END IF
*
*                     ==========================
*                      Info for CI coefficients 
*                     ==========================
*
* 
*. Information about block structure- needed by new PICO2 routine.
*. Memory for partitioning of C vector
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      NTTS = MXNTTS
C?    WRITE(6,*) ' GASCI : NTTS = ', NTTS
      CALL MEMMAN(KLCLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLCLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLCI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLCIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
      CALL MEMMAN(KLC2B  ,  NTTS,'ADDL  ',1,'C2BT  ')
*. Additional info required to construct partitioning 
*. Additional info required to construct partitioning 
*
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*
      CALL IAIBCM(ISPC,WORK(KLCIOIO))
*. option KSVST not active so
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLCBLTP),WORK(KSVST))
*
*. Batches  of C vector
      ICOMP = 1
      ISIMSYM = 0
*. Length of batch does not matter as we specified complete CI vector
      LBLOCK = 1     
      CALL PART_CIV2(IDC,WORK(KLCBLTP),WORK(KNSTSO(IATP)),
     &              WORK(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLCIOIO),
     &              ISMOST(1,ISM),
     &              NBATCH,WORK(KLCLBT),WORK(KLCLEBT),
     &              WORK(KLCI1BT),WORK(KLCIBT),ICOMP,ISIMSYM)
*. Number of BLOCKS
        NBLOCK = IFRMR(WORK(KLCI1BT),1,NBATCH)
     &         + IFRMR(WORK(KLCLBT),1,NBATCH) - 1
C?      WRITE(6,*) ' Number of blocks ', NBLOCK
*. Length of each block
      CALL EXTRROW(WORK(KLCIBT),8,8,NBLOCK,WORK(KLCI1BT))
*. Length of CI expansion 
      LENGTH_CI = IELSUM(WORK(KLCI1BT),NBLOCK)
*. alphasupergroup, betasupergroup=> class
        CALL MEMMAN(KLSPSPCL,NOCTPA*NOCTPB,'ADDL  ',1,'SPSPCL')
        CALL SPSPCLS(WORK(KLSPSPCL),WORK(KLOCCLS),NOCCLS)
*. Class of each block
        CALL MEMMAN(KLBLKCLS,NBLOCK,'ADDL  ',1,'BLKCLS')
        CALL MEMMAN(KLCLSL,NOCCLS,'ADDL  ',1,'CLSL  ')
        CALL MEMMAN(KLCLSLR,NOCCLS,'ADDL  ',2,'CLSL_R  ')
        CALL BLKCLS(WORK(KLCIBT),NBLOCK,WORK(KLBLKCLS),WORK(KLSPSPCL),
     &              NOCCLS,WORK(KLCLSL),NOCTPA,NOCTPB,WORK(KLCLSLR))
*
* The connection between the CI and CC coefficients are 
* the mappings to the the Occupation classes 
*
* KLBLKCLS : Occupation class for each CI block
* KEXTP_TO_OCCLS : Occupation type for each excitation type 
*
* Scratch vector for storing CI vector 
      CALL MEMMAN(KLCIVEC,LENGTH_CI,'ADDL  ',2,'CIVEC ')
      IF(IWAY.EQ.2) THEN
*. Collect CI vector in WORK(KLCIVEC)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Vector from LUCI '
          CALL WRTVCD(WORK(KLCIVEC),LUCI,1,-1)
          WRITE(6,*) ' LUCI,NBLOCK = ', LUCI,NBLOCK
        END IF
        CALL REWINO(LUCI)
        LBLK = -1
C       FRMDSCN(VEC,NREC,LBLK,LU)
        CALL FRMDSCN(WORK(KLCIVEC),NBLOCK,LBLK,LUCI)
      END IF
*. Four blocks of string occupations for creation/annihilation strings
      WRITE(6,*) ' MX_ST_TSOSO_BLK_MX = ', MX_ST_TSOSO_BLK_MX
      CALL MEMMAN(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
*. Space for string generation : Z matrices and strings
*. Also used to hold an NORB*NORB matrix  
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',2,'KLZSCR')
      CALL MEMMAN(KLZ1,LZ,'ADDL  ',1,'KLZ1  ')
      CALL MEMMAN(KLZ2,LZ,'ADDL  ',1,'KLZ2  ')
*. Occupation af alpha- and betastrings
      CALL MEMMAN(KLOCSTR1,MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS1')         
      CALL MEMMAN(KLOCSTR2,MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS2')         
*. Reorder arrays
      CALL MEMMAN(KLREO1,MAX_STR_SPGP,'ADDL  ',1,'KLREO1')
      CALL MEMMAN(KLREO2,MAX_STR_SPGP,'ADDL  ',1,'KLREO2')
*. An alpha and betastring
      CALL MEMMAN(KLSTRAL,NAEL,'ADDL  ',2,'STR_AL')
      CALL MEMMAN(KLSTRBE,NBEL,'ADDL  ',2,'STR_BE')
*
      CALL CC_CI_REORD_S(CCVEC,WORK(KLCIVEC),IWAY,ISPC,ISM,
     &     WORK(KLCIBT),NBLOCK,WORK(KLBLKCLS),
     &     WORK(KLSOBEX),NSPOBEX_TP,WORK(KLSOX_TO_OX),
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),WORK(KEX_TO_OC),
     &     WORK(KLSTR1_OCC), WORK(KLSTR2_OCC), WORK(KLSTR3_OCC),
     &     WORK(KLSTR4_OCC),WORK(KLZ1),WORK(KLZ2),
     &     WORK(KLREO1),WORK(KLREO2),WORK(KLOCSTR1),WORK(KLOCSTR2),
     &     WORK(KLZSCR),WORK(KLSTRAL),WORK(KLSTRBE),N_CC_AMP)

      IF(IWAY.EQ.1) THEN
*. Write resulting CI vector to DISC
C       TODSCN(VEC,NREC,LREC,LBLK,LU)
        CALL TODSCN(WORK(KLCIVEC),NBLOCK,WORK(KLCI1BT),LBLK,LUCI)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reordering between CI and CC order '
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' CC to CI reordering '
        ELSE IF ( IWAY.EQ.2) THEN
          WRITE(6,*) ' CI to CC reordering '
        END IF
        WRITE(6,*) ' Vector of CC coefficients '
        CALL WRTMAT(CCVEC,1,N_CC_AMP,1,N_CC_AMP)
        WRITE(6,*) ' Vector of CI coefficients '
        CALL WRTMAT(WORK(KLCIVEC),1,LENGTH_CI,1,LENGTH_CI)
      END IF
*
COLD  STOP ' Enforced stop at end of CC_CI_REORD '
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CC_CI_')
      RETURN
      END
      SUBROUTINE SCAGAT_CCVEC(CC_CMP,CC_EXP,ISG,NEXTP_SG,IEXTP_SG,
     &           IBEXTP,LEXTP,LEXTP_SG)
*
*. Scatter or gather blocks of CCVECTOR
*  ISG = 1 : Scatter from ICC_CMP to ICC_EXP
*  ISG = 2 : Gather  from ICC_EXP to ICC_CMP
*
* Blocks to be scattered/gathered are the NEXP_TP blocks in IEXP_TP
*
* Jeppe Olsen, April 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IEXTP_SG(NEXTP_SG)
*. Input : Offset and length of spinorbital excitation blocks
      INTEGER IBEXTP(*), LEXTP(*)
*. Input/Output
      DIMENSION CC_CMP(*),CC_EXP(*)
*. Output (from gather)
      INTEGER LEXTP_SG(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'SCAGAT_CCVEC speaking '
      END IF  
    
      IOFF_CMP = 1
      DO ITP = 1, NEXTP_SG
        IOFF_EXP = IBEXTP(IEXTP_SG(ITP))
        LEN = LEXTP(IEXTP_SG(ITP))
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' ITP, IOFF_EXP, LEN = ',ITP,IOFF_EXP,LEN
        END IF
        IF(ISG.EQ.1) THEN
*. Scatter
          CALL COPVEC(CC_CMP(IOFF_CMP),CC_EXP(IOFF_EXP),LEN)
        ELSE
*. Gather
          CALL COPVEC(CC_EXP(IOFF_EXP),CC_CMP(IOFF_CMP),LEN)
          LEXTP_SG(ITP) = LEN
        END IF
        IOFF_CMP = IOFF_CMP + LEN
      END DO
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Gathered list, Vector and offsets '
         LEN_CMP = IOFF_CMP-1
         CALL WRTMAT(CC_CMP,1,LEN_CMP,1,LEN_CMP)
         CALL IWRTMA(LEXTP_SG,1,NEXTP_SG,1,NEXTP_SG)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_SPOBTP_FOR_EXC_LEVEL(ILEVEL,ILEVEL_FOR_EXTP,
     &           NEXTP,NEXTP_AC,IEXTP_AC,ISOX_TO_OX)
*
* Total number and blocknumbers of spinorbital excitations with 
* excitation level ILEVEL
*
* Jeppe Olsen, April 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ILEVEL_FOR_EXTP(NEXTP), ISOX_TO_OX(*)
*. Output
      INTEGER IEXTP_AC(*)
*
      NEXTP_AC = 0
      DO IEXTP = 1, NEXTP
*. Excitation level for this spinexcitation type
        JLEVEL = ILEVEL_FOR_EXTP(ISOX_TO_OX(IEXTP))
        IF(JLEVEL.EQ.ILEVEL) THEN
          NEXTP_AC = NEXTP_AC + 1
          IEXTP_AC(NEXTP_AC) = IEXTP
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Spinorbital excitation blocks with excitation level',ILEVEL
        WRITE(6,*) 
     &  ' Number of obtained spin-orbital excitation types', NEXTP_AC 
        WRITE(6,*) ' And the corresponding blocknumbers : '
        CALL IWRTMA(IEXTP_AC,1,NEXTP_AC,1,NEXTP_AC)
      END IF
*
      RETURN
      END
      SUBROUTINE CC_CI_REORD_S(CCVEC,CIVEC,IWAY,ISPC,ISM,
     &           ICIBLK,NBLOCK_CI,ICIBLK_OCCLS,
     &           ISOBEX,NSOBEX_TP,ISOX_TO_OX,LSOBEX,IBSOBEX,
     &           IEX_TO_OC,
     &           ICA_STR,ICB_STR,IAA_STR,IAB_STR,
     &           IZA,IZB,IREOA,IREOB,IOCSTA,IOCSTB,IZSCR,
     &           ISTRAL,ISTRBE,N_CC_AMP)  
*
* Inner routine ( sounds nicer than slave routine )
* for reordering between CI and CC orders.
* 
* Only reordering is performed, no scaling
*
* Jeppe Olsen, March 28 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc' 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'newccp.inc'
C     I_USE_NEWCCP
*. Input and output 
      DIMENSION CCVEC(*),CIVEC(*)
*. Input
      INTEGER ICIBLK(8,NBLOCK_CI), ICIBLK_OCCLS(NBLOCK_CI)
      INTEGER ISOBEX(4*NGAS,NSOBEX_TP),ISOX_TO_OX(NSOBEX_TP)
      INTEGER LSOBEX(NSOBEX_TP),IBSOBEX(NSOBEX_TP) 
      INTEGER IEX_TO_OC(*)
*. Space for creation and annihilation strings of given symmetry
      INTEGER ICA_STR(*),ICB_STR(*),IAA_STR(*),IAB_STR(*)
*. Space for strings, reorder arrays, and Z matrices, and scratch for 
*. constructing Z
      INTEGER IZA(*),IZB(*),IREOA(*),IREOB(*),IOCSTA(*),IOCSTB(*)
      INTEGER IZSCR(*)
*. Space for a single alpha and beta string
      INTEGER ISTRAL(*), ISTRBE(*)
*. Local scratch : Occupation in Reference space
      INTEGER IREF_OCC_AL(MXPNGAS),IREF_OCC_BE(MXPNGAS)
*. Actual reference strings 
      INTEGER IREF_STR_AL(MXPNGAS),IREF_STR_BE(MXPNGAS)
*. General occupation of a pair of alpha- and beta-strings
      INTEGER IOCC_AL(MXPNGAS),IOCC_BE(MXPNGAS)
*. And the corresponding groups 
      INTEGER IGRP_AL(MXPNGAS),IGRP_BE(MXPNGAS)
*. Offsets to CI blocks, with given TT as a function of symmetry of
*. alpha strings
      INTEGER ICIBLK_OFF(MXPOBS)
*. For group notation of annihilation/creation strings
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS)
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*. For local testing
CTEST INTEGER ITOUCH(1000)
*
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST WRITE(6,*) ' Jeppe : Remember local tests are active '
CTEST IZERO = 0
CTEST CALL ISETVC(ITOUCH,IZERO,N_CC_AMP+1)
*
      
      ITP_AL = 1
      ITP_BE = 2
      NEL_AL = NELEC(ITP_AL)
      NEL_BE = NELEC(ITP_BE)
      I_CC = 0
*
      NTEST = 1000
*. Check sums for CI and CC adressing
      ICC_CHECK = 0
      ICI_CHECK = 0
*
C?    WRITE(6,*) ' Included Spinorbital excitations'
C?    CALL WRT_SPOX_TP(ISOBEX,NSOBEX_TP)
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CI_CIS')
C     GET_REF_ALBE_OCC(IREFSPC,IREF_AL,IREF_BE)
*. Obtain alpha and beta occupation of reference space
      IREFSPC = 1
      CALL GET_REF_ALBE_OCC(IREFSPC,IREF_OCC_AL,IREF_OCC_BE)
*. Obtain atual alpha and beta strings for reference space 
      CALL GET_REF_ALBE_STR(IREFSPC,IREF_STR_AL,IREF_STR_BE)
*. Symmetry of reference strings
      ISM_REF_AL = ISYMST(IREF_STR_AL,NEL_AL)
      ISM_REF_BE = ISYMST(IREF_STR_BE,NEL_BE)
*. Loop over spinorbital excitation blocks
      DO JSOBEX = 1, NSOBEX_TP
C?     WRITE(6,*) ' Output for JSOBEX = ', JSOBEX
*. Resulting occupation of alpha and beta strings
C     EXOCC_STROCC(ICR_OCC,IAN_OCC,ISTR_IN_OCC,
C    &           ISTR_OUT_OCC,NGAS,IZERO_STR)
*. Occupation of alpha string
        CALL EXOCC_STROCC(ISOBEX(1+0*NGAS,JSOBEX),
     &       ISOBEX(1+2*NGAS,JSOBEX),IREF_OCC_AL,
     &       IOCC_AL,NGAS,IZERO_ALSTR)
        CALL OCC_TO_GRP(IOCC_AL,IGRP_AL,1)
*. Occupation of betastring 
        CALL EXOCC_STROCC(ISOBEX(1+1*NGAS,JSOBEX),
     &       ISOBEX(1+3*NGAS,JSOBEX),IREF_OCC_BE,
     &       IOCC_BE,NGAS,IZERO_BESTR)
        CALL OCC_TO_GRP(IOCC_BE,IGRP_BE,1)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Occupation of resulting strings for JSOBEX=',
     &    JSOBEX
          CALL IWRTMA(IOCC_AL,1,NGAS,1,NGAS)     
          CALL IWRTMA(IOCC_BE,1,NGAS,1,NGAS)     
          WRITE(6,*) ' And the corresponding groups '
          CALL IWRTMA(IGRP_AL,1,NGAS,1,NGAS)
          CALL IWRTMA(IGRP_BE,1,NGAS,1,NGAS)
        END IF
*. Supergroups corresponding to these occupation 
        CALL FIND_SPGRP_FROM_OCC(IOCC_AL,ISPGRP_AL,ITP_AL)
        CALL FIND_SPGRP_FROM_OCC(IOCC_BE,ISPGRP_BE,ITP_BE)
*. Relative number of these supergroups 
        ISPGRP_AL_REL = ISPGRP_AL - IBSPGPFTP(ITP_AL) + 1
        ISPGRP_BE_REL = ISPGRP_BE - IBSPGPFTP(ITP_BE) + 1
*. And then the TTS blocks with these occupation
C            CIBLOCKS_FOR_TT(ICIBLK,NCIBLK,IATP,IBTP,IFORM,ITTBLK)
        CALL CIBLOCKS_FOR_TT(ICIBLK,NBLOCK_CI,ISPGRP_AL_REL,
     &       ISPGRP_BE_REL,2,ICIBLK_OFF)
*. Transform creation/annihilations type from occupation to group notation 
        CALL OCC_TO_GRP(ISOBEX(1+0*NGAS,JSOBEX),IGRP_CA,1)
        CALL OCC_TO_GRP(ISOBEX(1+1*NGAS,JSOBEX),IGRP_CB,1)
        CALL OCC_TO_GRP(ISOBEX(1+2*NGAS,JSOBEX),IGRP_AA,1)
        CALL OCC_TO_GRP(ISOBEX(1+3*NGAS,JSOBEX),IGRP_AB,1)
*
        NEL_CA = IELSUM(ISOBEX(1+0*NGAS,JSOBEX),NGAS)
        NEL_CB = IELSUM(ISOBEX(1+1*NGAS,JSOBEX),NGAS)
        NEL_AA = IELSUM(ISOBEX(1+2*NGAS,JSOBEX),NGAS)
        NEL_AB = IELSUM(ISOBEX(1+3*NGAS,JSOBEX),NGAS)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' NEL_CA, NEL_CB, NEL_AA, NEL_AB',
     &                 NEL_CA, NEL_CB, NEL_AA, NEL_AB
        END IF

*. Loop over symmetries of creation/annihilation strings
*. Symmetry of excitations is assumed to be 1 (total sym)
        ISM = 1
        DO ISM_C = 1, NSMST
          ISM_A = MULTD2H(ISM,ISM_C) 
          DO ISM_CA = 1, NSMST
            ISM_CB = MULTD2H(ISM_C,ISM_CA)
            DO ISM_AA = 1, NSMST
             ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*. Obtain creation/annihilation strings
             CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,
     &            NSTR_CA,ICA_STR, NTOOB,0,IDUM,IDUM)
             CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,
     &            NSTR_CB,ICB_STR, NTOOB,0,IDUM,IDUM)
             CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,
     &            NSTR_AA,IAA_STR, NTOOB,0,IDUM,IDUM)
             CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,
     &            NSTR_AB,IAB_STR, NTOOB,0,IDUM,IDUM)
*
C?       WRITE(6,*) ' Beta annihilation strings : '
C?       CALL IWRTMA(IAB_STR,NEL_AB,NSTR_AB,NEL_AB,NSTR_AB)
*. Corresponding symmetries of alpha and beta strings
             ISM_OP_AL = MULTD2H(ISM_CA,ISM_AA)
             ISM_OP_BE = MULTD2H(ISM_CB,ISM_AB)
*. Symmetry of alpha and beta strings
             ISM_STR_AL = MULTD2H(ISM_REF_AL,ISM_OP_AL)
             ISM_STR_BE = MULTD2H(ISM_REF_BE,ISM_OP_BE)
*. Obtain all alpha and beta strings of correct supergroup and sym
*. The mapping from occupation to address will be used on the following
*.. Generate information about IA strings
C                 WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,ISCR,NTEST)        
             NTESTL = 0
             CALL WEIGHT_SPGP(IZA,NGAS,IOCC_AL,NOBPT,IZSCR,NTESTL)
C            GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
             CALL GETSTR2_TOTSM_SPGP(IGRP_AL,NGAS,ISM_STR_AL,NEL_AL,
     &            NSTR_AL,IOCSTA,NOCOB,1,IZA,IREOA)
C?           WRITE(6,*) ' Reorder array for alpha strings '
C?           CALL IWRTMA(IREOA,1,NSTR_AL,1,NSTR_AL)
*. And about beta string 
             CALL WEIGHT_SPGP(IZB,NGAS,IOCC_BE,NOBPT,IZSCR,NTESTL)
             CALL GETSTR2_TOTSM_SPGP(IGRP_BE,NGAS,ISM_STR_BE,NEL_BE,
     &            NSTR_BE,IOCSTB,NOCOB,1,IZB,IREOB)
*. Loop over T elements as  matric T(I_CA, I_CB, IAA, I_AB)
             DO I_AB = 1, NSTR_AB
              DO I_AA = 1, NSTR_AA
               DO I_CB = 1, NSTR_CB
                DO I_CA = 1, NSTR_CA
*. Alpha string obtained by alpha crea alpha anni alpha refstring
                 IOFF_CA = 1 + (I_CA-1)*NEL_CA
                 IOFF_AA = 1 + (I_AA-1)*NEL_AA
                 CALL CRAN_STR(ICA_STR(IOFF_CA),IAA_STR(IOFF_AA),
     &                NEL_CA,NEL_AA,IREF_STR_AL,NEL_AL,
     &                ISTRAL,ISIGN_AL,IZERO_AL) 
*. And number of this string 
C       ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
                 IANUM = ISTRNM(ISTRAL,NOCOB,NEL_AL,IZA,IREOA,1)
                 IF(NTEST.GE.1000) WRITE(6,*) ' IANUM = ', IANUM
*. And for beta string 
                 IOFF_CB = 1 + (I_CB-1)*NEL_CB
                 IOFF_AB = 1 + (I_AB-1)*NEL_AB
                 IF(NTEST.GE.1000) THEN
                   WRITE(6,*) ' I_AB, IOFF_AB = ',
     &                          I_AB, IOFF_AB
                 END IF
                 CALL CRAN_STR(ICB_STR(IOFF_CB),IAB_STR(IOFF_AB),
     &                NEL_CB,NEL_AB,IREF_STR_BE,NEL_BE,
     &                ISTRBE,ISIGN_BE,IZERO_BE) 
*. And number of this string 
                 IBNUM = ISTRNM(ISTRBE,NOCOB,NEL_BE,IZB,IREOB,1)
                 IF(NTEST.GE.1000) WRITE(6,*) ' IBNUM = ', IBNUM
                 
*. Number in CC order 
                 I_CC = I_CC + 1
                 ICC_CHECK = ICC_CHECK + I_CC
*. Number in CI order 
                 I_CI = ICIBLK_OFF(ISM_STR_AL)-1+(IBNUM-1)*NSTR_AL+
     &                  IANUM
CTEST            ITOUCH(I_CI) = ITOUCH(I_CI) + 1
                 ICI_CHECK = ICI_CHECK + I_CI
                 IF(NTEST.GE.1000) THEN
                   WRITE(6,*) 'ICIBLK_OFF, NSTR_AL = ',
     &                          ICIBLK_OFF(ISM_STR_AL),NSTR_AL
                   WRITE(6,'(A,4I4)') ' I_AB, I_AA, I_CB, I_CA', 
     &                                  I_AB, I_AA, I_CB, I_CA 
                   WRITE(6,*) 'I_CC, I_CI = ',I_CC,I_CI

                 END IF
*
                 IF(I_USE_NEWCCP.EQ.0) THEN
                   SIGN = DBLE(ISIGN_AL*ISIGN_BE) 
                 ELSE 
                   IF(MOD(NEL_CB*NEL_AA,2).EQ.0) THEN
                     SIGN_CBAA = 1
                   ELSE 
                     SIGN_CBAA = -1
                   END IF
                   SIGN = SIGN_CBAA* DBLE(ISIGN_AL*ISIGN_BE)
                 END IF
*
                 IF(IWAY.EQ.1) THEN
                   CIVEC(I_CI) = CCVEC(I_CC)*SIGN
                 ELSE
                   CCVEC(I_CC) = CIVEC(I_CI)*SIGN
                 END IF
*
                 IF(I_CI.LE.0.OR.I_CI.GT.N_CC_AMP+1) THEN 
                   WRITE(6,*) ' I_CI out of range = ', I_CI 
                 END IF
*
                END DO
               END DO
             END DO
*       C    ^ End of loop over elements of block
            END DO
           END DO
          END DO
        END DO
*       ^ End of loop over symmetries of creation/annihilation strings
      END DO
*     ^ End of loop over types of CC excitations
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CC and CI check sums = ', ICC_CHECK,ICI_CHECK
      END IF
*
*. Local tests : Print the numbers for the CI coefficients that 
*. were not touched exactly once
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST WRITE(6,*) ' Local tests : '
CTEST DO ICI = 1, N_CC_AMP + 1
CTEST   IF(ITOUCH(ICI).NE.1) THEN
CTEST     WRITE(6,*) ' Mess : ICI, ITOUCH(ICI) = ',
CTEST&                        ICI, ITOUCH(ICI)
CTEST   END IF
CTEST END DO
*
*. The part below does not work as we do not know the number of 
*. the CI det corresponding to the HF reference
C     IF(ICI_CHECK.NE.ICC_CHECK+1) THEN 
C       WRITE(6,*) ' Problem in reord, inconsistent checksums'
C       WRITE(6,*) ' ICC_CHECK, ICI_CHECK = ', ICC_CHECK,ICI_CHECK
C       STOP       ' Problem in reord, inconsistent checksums'
C     END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CI_CIS')
      RETURN
      END
      SUBROUTINE GET_REF_ALBE_OCC(IREFSPC,IREF_AL,IREF_BE)
*
* Obtain alpha and beta occupations for reference space 
*
* Reference space is assumed to be a single pair of occupations of 
* alpha and beta strings ( this includes closed shell HF,            
* Highspin open shell and CAS reference)
*
* Only a single valence orbital space is assumed
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc' 
      INCLUDE 'strinp.inc'
*. Output : Alpha and beta occupations for each GAS space 
      INTEGER IREF_AL(NGAS),IREF_BE(NGAS)
*
*. Total number of Hole orbitals
      NHOLE = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) THEN
          NHOLE = NHOLE + NGSOBT(IGAS)
        END IF
      END DO
*. Number of orbitals in valence space
      NELEC_AL = NELEC(1)
      NELEC_BE = NELEC(2)
      NVAL_AL = NELEC_AL - NHOLE
      NVAL_BE = NELEC_BE - NHOLE
*
      DO IGAS = 1, NGAS
        NORB = NGSOBT(IGAS)
        IF(IHPVGAS(IGAS).EQ.1) THEN
          IREF_AL(IGAS) = NORB
          IREF_BE(IGAS) = NORB
        ELSE IF( IHPVGAS(IGAS).EQ.2) THEN
          IREF_AL(IGAS) = 0
          IREF_BE(IGAS) = 0
        ELSE IF( IHPVGAS(IGAS).EQ.3) THEN
          IREF_AL(IGAS) = NVAL_AL
          IREF_BE(IGAS) = NVAL_BE
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Alpha and Beta Occupation for reference space '
        WRITE(6,*)
        CALL IWRTMA(IREF_AL,1,NGAS,1,NGAS)
        CALL IWRTMA(IREF_BE,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_REF_ALBE_STR(IREFSPC,IREF_AL,IREF_BE)
*
* Obtain alpha and beta strings for reference space 
*
* Reference space is assumed to be a single pair of STRINGS of 
* alpha and beta strings ( this includes closed shell HF,            
* Highspin open shell, but nor CAS !!)
*
* Only a single valence orbital space is assumed
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc' 
      INCLUDE 'strinp.inc'
*. Output : Alpha and beta occupations for each GAS space 
      INTEGER IREF_AL(NGAS),IREF_BE(NGAS)
*
*. Total number of Hole orbitals
      NHOLE = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) THEN
          NHOLE = NHOLE + NGSOBT(IGAS)
        END IF
      END DO
*. Number of orbitals in valence space
      NELEC_AL = NELEC(1)
      NELEC_BE = NELEC(2)
      NVAL_AL = NELEC_AL - NHOLE
      NVAL_BE = NELEC_BE - NHOLE
*
      IOFF = 1
      IOFF_AL = 1
      IOFF_BE = 1
*
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          IOFF = 1
        ELSE 
          IOFF = IOFF + NGSOBT(IGAS-1)
        END IF
        NORB = NGSOBT(IGAS)
        IF(IHPVGAS(IGAS).EQ.1) THEN
          DO IORB = 1, NORB
            IREF_AL(IOFF_AL-1+IORB) = IOFF-1+IORB
            IREF_BE(IOFF_BE-1+IORB) = IOFF-1+IORB
          END DO
          IOFF_AL = IOFF_AL + NORB
          IOFF_BE = IOFF_BE + NORB
        ELSE IF( IHPVGAS(IGAS).EQ.3) THEN
          DO IORB = 1, NVAL_AL 
            IREF_AL(IOFF_AL-1+IORB) = IOFF-1+IORB
          END DO
          IOFF_AL = IOFF_AL + NVAL_AL
          DO IORB = 1, NVAL_BE 
            IREF_BE(IOFF_BE-1+IORB) = IOFF-1+IORB
          END DO
          IOFF_BE = IOFF_BE + NVAL_BE
        END IF
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Alpha and Beta strings for reference space '
        WRITE(6,*)
        CALL IWRTMA(IREF_AL,1,NELEC_AL,1,NELEC_AL)
        CALL IWRTMA(IREF_BE,1,NELEC_BE,1,NELEC_BE)
      END IF
*
      RETURN
      END
      SUBROUTINE EXOCC_STROCC(ICR_OCC,IAN_OCC,ISTR_IN_OCC,
     &           ISTR_OUT_OCC,NGAS,IZERO_STR)
*
* Occupation of excitaion op,  occupation of string =>
* Occupation of excited string
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICR_OCC(*),IAN_OCC(*),ISTR_IN_OCC(*)
*. Output
      INTEGER ISTR_OUT_OCC(*)
*. Annihilation
      IZERO_STR = 0
      DO IGAS = 1, NGAS
        ISTR_OUT_OCC(IGAS) = ISTR_IN_OCC(IGAS) - IAN_OCC(IGAS)
        IF(ISTR_OUT_OCC(IGAS).LT.0) IZERO_STR = 1
      END DO
*. Creation
      DO IGAS = 1, NGAS
        ISTR_OUT_OCC(IGAS) = ISTR_OUT_OCC(IGAS) + ICR_OCC(IGAS)
      END DO 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from EXOCC_STROCC '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        WRITE(6,*) ' Occ of crea string : '
        CALL IWRTMA(ICR_OCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occ of anni string '
        CALL IWRTMA(IAN_OCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occ of input string '
        CALL IWRTMA(ISTR_IN_OCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occ of output string '
        CALL IWRTMA(ISTR_OUT_OCC,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE CRAN_STR(ICR,IAN,NCR,NAN,ISTR_IN,NEL_IN,
     &                    ISTR_OUT,ISIGN,IZERO_STR)
*
* ISTR_OUT = ISIGN* ICR IAN ISTR_IN
*
* Where ICR is a string of creation operators and IAN is a string 
* of annihilation operators.
*
* Input string is assumed to be given in ascending order, 
* and output string will be delivered with orbitals in 
* ascending order
*
*. Initial version, I hope it is not for mission critical routines 
* (could be speeded up)
*
* Jeppe Olsen, March 2000
*
* Change of phase of annihilations strings, Oct2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'newccp.inc'
*. Input
      INTEGER ICR(NCR),IAN(NAN)
      INTEGER ISTR_IN(NEL_IN)
*. Output
      INTEGER ISTR_OUT(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CRAN_STR speaking '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' Input string '
        CALL IWRTMA(ISTR_IN,1,NEL_IN,1,NEL_IN) 
        WRITE(6,*) ' Annihilation string '
        CALL IWRTMA(IAN,1,NAN,1,NAN)
        WRITE(6,*) ' Creation string '
        CALL IWRTMA(ICR,1,NCR,1,NCR)
      END IF
*. Make sure that annihilation strings are properly increasing
C?    DO JAN = 2, NAN
C?      IF(IAN(JAN).LE.IAN(JAN-1)) THEN
C?        WRITE(6,*) ' CRAN confused, strange annihilation string :' 
C?        CALL IWRTMA(IAN,1,NAN,1,NAN)
C?        STOP ' CRAN confused, strange annihilation string'
C?      END IF
C?    END DO
*
      NEL_OUT = NEL_IN - NAN + NCR
*
      IZERO_STR = 0
      ISIGN = 1.0D0
      CALL ICOPVE(ISTR_IN,ISTR_OUT,NEL_IN)
*. Annihilate  : IAN(1) IAN(2) .... !STR_IN>
      DO IANNI = 1, NAN
        IFOUND = 0
        DO IEL = 1, NEL_IN-IANNI+1
C?        WRITE(6,*) ' CRAN : IANNI IEL ISTR IAN ',
C?   &    IANNI,IEL,ISTR_OUT(IEL),IAN(NAN-IANNI+1)
          IF(ISTR_OUT(IEL).EQ.IAN(NAN-IANNI+1)) THEN 
            ISIGN = ISIGN*(-1)**(IEL-1)
            IFOUND = 1
            DO JEL = IEL, NEL_IN-IANNI
              ISTR_OUT(JEL) = ISTR_OUT(JEL+1)
             END DO
          END IF
        END DO
        IF(IFOUND.EQ.0) THEN
*. orbital to be annihilated not found, output string is zero
          IZERO_STR = 1
          GOTO 1001 
        END  IF
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Annihilated string '
        CALL IWRTMA(ISTR_OUT,NEL_IN-NAN,1,NEL_IN-NAN,1)
      END IF
*. Creation part 
      DO ICREA = 1, NCR
*. Place to insert orbital 
        ICR_ORB = ICR(NCR-ICREA+1)
        IPLACE = 1
        DO IEL = 1, NEL_IN-NAN + ICREA-1
          IF(ISTR_OUT(IEL).EQ.ICR_ORB) THEN
*. Electron is already around, zero
            IZERO_STR = 1
            GOTO 1001 
          END IF
*
          IF(IEL.LT.NEL_IN-NAN + ICREA-1) THEN
            IF(ISTR_OUT(IEL).LT.ICR_ORB.AND.
     &         ISTR_OUT(IEL+1).GT.ICR_ORB   ) THEN
               IPLACE = IEL+1
            END IF
          ELSE IF(IEL.EQ.NEL_IN-NAN+ICREA-1) THEN
            IF(ISTR_OUT(IEL).LT. ICR_ORB   ) THEN
              IPLACE = IEL + 1
            END IF
          END IF
        END DO
        ISIGN  = ISIGN*(-1)**(IPLACE-1)
        DO IEL = NEL_IN-NAN+ICREA,IPLACE + 1, -1
          ISTR_OUT(IEL) = ISTR_OUT(IEL-1)
        END DO
        ISTR_OUT(IPLACE) = ICR_ORB
      END DO
*
 1001 CONTINUE
*. A bit on the sign : In LUCIA the order of the annihilation is actually 
* descending !, find permutation sign
      IF(I_USE_NEWCCP.EQ.0.AND.NAN.GT.1) THEN
        NPERM = NAN*(NAN-1)/2
        ISIGN = ISIGN*(-1)**NPERM
      END IF
*
      IF(NTEST.GE.100) THEN 
        IF(IZERO_STR.EQ.0) THEN
          WRITE(6,*) ' Output string '
          CALL IWRTMA(ISTR_OUT,1,NEL_OUT,1,NEL_OUT)
          WRITE(6,*) ' ISIGN = ', ISIGN
        ELSE 
          WRITE(6,*) ' Vanishing string '
        END IF
      END IF
*
      RETURN
      END
C      
      SUBROUTINE FIND_SPGRP_FROM_OCC(IOCC,ISPGRP_NUM,ITP)
*
* Find the number(ISPGRP_NUM) corresponding to supergroup of type ITP 
* with occupation IOCC. If ITP = 0, all supergroup types are checked
* Returned supergroup number is absolute supergroup number 
*
* If no supergroup is identified a zero is returned
*
* Jeppe Olsen, April 2000
* ITP = 0 option added, March 2007 - not tested..
*
      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
*. Specific input 
      INTEGER IOCC(*)
*
      ISPGRP_NUM = 0
      IF(ITP.EQ.0) THEN
        ITP_MIN = 1
        ITP_MAX = NTSPGP
      ELSE
        ITP_MIN = ITP
        ITP_MAX = ITP
      END IF
*
      DO IITP = ITP_MIN, ITP_MAX
       DO  JSPGP = IBSPGPFTP(IITP), IBSPGPFTP(IITP)+NSPGPFTP(IITP)-1
         IDENTICAL = 1
         DO IGAS = 1, NGAS
           IF(NELFSPGP(IGAS,JSPGP).NE.IOCC(IGAS)) IDENTICAL = 0
         END DO
         IF(IDENTICAL.EQ.1) ISPGRP_NUM = JSPGP
       END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation of supergroup : '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        IF(ISPGRP_NUM.EQ.0) THEN  
          WRITE(6,*) ' Not identified '
        ELSE 
          WRITE(6,*) ' Number of supergroup : ', ISPGRP_NUM  
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE CIBLOCKS_FOR_TT(ICIBLK,NCIBLK,IATP,IBTP,IFORM,ITTBLK)
*
* A set of CI blocks is specified through ICIBLK
* Find block with TYPES IATP, IBTP
*
*
* Output :  ITTBLK : The blocks : IFORM = 1 => The number of the block
*                                       = 2 => Offset of block
*
*. Jeppe Olsen, April 20000
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'csm.inc'
      INTEGER ICIBLK(8,NCIBLK)
*. Output
      INTEGER ITTBLK(*)
*
      IZERO = 0
      CALL ISETVC(ITTBLK,IZERO,NSMST)
      DO JCIBLK = 1, NCIBLK
C?      WRITE(6,*) ' JCIBLK, IA and IB : ',
C?   &  JCIBLK,ICIBLK(1,JCIBLK),ICIBLK(2,JCIBLK)
        IF(ICIBLK(1,JCIBLK).EQ.IATP.AND.ICIBLK(2,JCIBLK).EQ.IBTP) THEN
          IASM = ICIBLK(3,JCIBLK)
          IOFF = ICIBLK(5,JCIBLK)
          IF(IFORM.EQ.1) THEN
            ITTBLK(IASM) = JCIBLK
          ELSE IF(IFORM.EQ.2) THEN
            ITTBLK(IASM) = IOFF
          END IF
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,2I6)') ' Blocks with IATP, IBTP ', IATP,IBTP
        IF(IFORM.EQ.1) THEN 
          WRITE(6,*) ' Block numbers '
        ELSE
          WRITE(6,*) ' Offsets '
        END IF
        CALL IWRTMA(ITTBLK,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END 
      
c $Id$
