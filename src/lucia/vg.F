      SUBROUTINE COUNT_OCCURENCE(IVEC,IOCC,NELMNT,MAXVAL)
* A string of integers IVEC is given with elements 
* from 1 to MAXVAL. Find the number of times each integer 
* occurs and store this info in IOCC
*
* Jeppe Olsen, Sept 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IVEC(NELMNT)
*. Output
      INTEGER IOCC(MAXVAL)
*
      CALL ISETVC(IOCC,0,MAXVAL)
      DO I = 1, NELMNT
        IOCC(IVEC(I)) = IOCC(IVEC(I)) + 1
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Integer vector : '
        CALL IWRTMA(IVEC,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' Number of occurence for each integer '
        CALL IWRTMA(IOCC,1,MAXVAL,1,MAXVAL)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_OBEX(IFLAG,NCREA,NANNI,NOBEX_TP,IOBEX_TP)
*
* Generate all possible excitation types containing 
* NCREA excitation operators and NANNI annihilation operators
*
* IFLAG = 1 => Just number of excitation types
*       = 0 => Number and the actual types
*
* Output
* ======
* NOBEX_TP : Number of orbital excitation types
* IOBEX_TP : The actual orbital excitation types (if IFLAG = 0)  
*
* Jeppe Olsen, Helsinki, Sept 2001
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Output
      DIMENSION IOBEX_TP(2*NGAS,*)
*. Local scratch
      DIMENSION IORD_C(MXPLCCOP), IORD_A(MXPLCCOP)
      DIMENSION IOCC_C(MXPNGAS),IOCC_A(MXPNGAS)
*. Loop over possible creation strings
C  NXTORD_NS(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GEN_OBEX speaking '
        WRITE(6,*) ' NCREA, NANNI = ', NCREA, NANNI 
      END IF
      NOBEX_TP = 0 
      IFIRST_C = 1
      NONEW_C = 0
 1002 CONTINUE
      IF(IFIRST_C.EQ.1) THEN
        CALL ISETVC(IORD_C,1,NCREA)
        NONEW_C = 0
      ELSE
        WRITE(6,*) ' NXTORD_S will be called for CREA '
        CALL NXTORD_NS(IORD_C,NCREA,1,NGAS,NONEW_C)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' New IORD_C '
        CALL IWRTMA(IORD_C,NCREA,1,NCREA,1)
        END IF
      IFIRST_C = 0
      IF(NONEW_C.EQ.0) THEN
*. Creation operator in form with number of crea per gasspace
C            COUNT_OCCURENCE(IVEC,IOCC,NELMNT,MAXVAL)
        CALL COUNT_OCCURENCE(IORD_C,IOCC_C,NCREA,NGAS)
*. Is number of operators in each space smaller than 2*number 
*  orbitals in this space
        IOKAY_C = 1
        DO IGAS = 1, NGAS
          IF(IOCC_C(IGAS).GT.2*NOBPT(IGAS))IOKAY_C = 0
          WRITE(6,*) ' IGAS, IORD, NOBPT ',
     &    IGAS,IOCC_C(IGAS),NOBPT(IGAS)
        END DO
        WRITE(6,*) ' IOKAY_C = ', IOKAY_C
        IF(IOKAY_C .EQ. 1) THEN
*. Generate annihilation strings 
          IFIRST_A = 1
 1001     CONTINUE
          IF(IFIRST_A.EQ.1) THEN
            CALL ISETVC(IORD_A,1,NANNI)
            NONEW_A = 0
          ELSE
            WRITE(6,*) ' NXTORD_S will be called for ANNI '
            CALL NXTORD_NS(IORD_A,NANNI,1,NGAS,NONEW_A)
          END IF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' New IORD_A '
            CALL IWRTMA(IORD_A,NANNI,1,NANNI,1)
          END IF
          IFIRST_A = 0
          IF(NONEW_A.EQ.0) THEN
*. Annihilation operator in form with number of anni per gasspace
C                COUNT_OCCURENCE(IVEC,IOCC,NELMNT,MAXVAL)
            CALL COUNT_OCCURENCE(IORD_A,IOCC_A,NANNI,NGAS)
*. Is number of operators in each space smaller than 2*number 
*  orbitals in this space
            IOKAY_A = 1
            DO IGAS = 1, NGAS
              IF(IOCC_A(IGAS).GT.2*NOBPT(IGAS))IOKAY_A = 0
            END DO
            IF(IOKAY_A .EQ. 1) THEN
*. A new operator has been born
              NOBEX_TP = NOBEX_TP + 1
              IF(IFLAG.EQ.0) THEN
                CALL ICOPVE(IOCC_C,IOBEX_TP(1,NOBEX_TP)     ,NGAS)
                CALL ICOPVE(IOCC_A,IOBEX_TP(1+NGAS,NOBEX_TP),NGAS)
              END IF
            END IF
*           ^ End if IOKAY_A = 1
            GOTO 1001 
          END IF
*         ^ End of NONEW = 0 for anni string
        END IF
*       ^ End if IOKAY_C = 1
      GOTO 1002
      END IF
*     ^ End if NONEW = 0 for crea string
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of orbital excitation types generated ', 
     &               NOBEX_TP
        IF(IFLAG.EQ.0) THEN
          WRITE(6,*) ' Creation part,  Annihilation  part '
          WRITE(6,*) ' ==================================='
          DO IOBEX = 1, NOBEX_TP
            WRITE(6,'(16I4,16I4)')
     &      (IOBEX_TP(I,IOBEX),I=1, NGAS),
     &      (IOBEX_TP(NGAS+I,IOBEX),I=1, NGAS) 
          END DO
        END IF
      END IF
*
      RETURN 
      END
      SUBROUTINE VGSIGDEN_M(ISIGDEN)
*
* Master Routine for general density/sigma routine
*
* Calculates densities <L!Op!R>
*
* Where L and R not neccessarily have the same number of 
* electrons
*
* Right hand side is defined by Space in /CANDS/ whereas
* left hand side is defined by GS state in /CRUN/
*
* Results is pt returned in work(klt), where klt is determined 
* by this subroutine !!! 
*
* Jeppe Olsen, Sept 2001, for the Quantum Dot project
*
* 
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'glstate.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'ctcc.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'VGSI_M')
*
      NTEST = 00
*
      WRITE(6,*)
      WRITE(6,*) ' ***********************************'
      WRITE(6,*) ' Control transferred to VG routines '
      WRITE(6,*) ' ***********************************'
      WRITE(6,*)
*. Number of alpha- and beta-electrons  in left hand side state
*
      NELEC_L = IGST_OCC(NGAS,2)
      NAELEC_L = (NELEC_L+IGST_MS2)/2
      NBELEC_L = (NELEC_L-IGST_MS2)/2
      IF( NAELEC_L+ NBELEC_L .NE.  NELEC_L) THEN
        WRITE(6,*) ' Inconsistent NELEC and MS2 for LHS '
        WRITE(6,*) ' NELEC, MS2 = ', NELEC_L, IGST_MS2 
        STOP ' Inconsistent NELEC and MS2 for LHS '
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of electrons in LHS state ', NELEC_L
        WRITE(6,*) ' Number of alpha-electrons in LHS state ',NAELEC_L
      END IF
*. Required change in MS2
      MS2_RL = IGST_MS2 - MS2 
      WRITE(6,*) ' Required change in MS2 = ', MS2_RL
*.Corresponding types : First types with this number of electrons
      IATP_L = 0
      IBTP_L = 0
      DO ITP = 1, NSTTP
       IF(NELFTP(ITP).EQ.NAELEC_L.AND.IATP_L.EQ.0) IATP_L = ITP
       IF(NELFTP(ITP).EQ.NBELEC_L.AND.IBTP_L.EQ.0) IBTP_L = ITP
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' alpha and beta types of LHS ',
     &                 IATP_L,IBTP_L
*. Occupation change compared to reference and right state
      IATP_R = 1
      IBTP_R = 2
      NAELEC_R = NELFTP(IATP_R)
      NBELEC_R = NELFTP(IBTP_R)
      NELEC_R = NAELEC_R+NBELEC_R
      IDEL_RL = NELEC_R - NELEC_L
      IADEL_RL = NAELEC_R - NAELEC_L
      IBDEL_RL = NBELEC_R - NBELEC_L
      IF(NTEST.GE.100) THEN
       WRITE(6,*) '  IDEL_RL,  IADEL_RL,  IBDEL_RL = ',
     &               IDEL_RL,  IADEL_RL,  IBDEL_RL
      END IF
*. In this initial version the number of electrons in R state 
*. should be greater than the number of electrons in L state
      IF(IDEL_RL.LT.0) THEN
        WRITE(6,*) ' Problem in VG routines'
        WRITE(6,*) ' More electrons in L state than in R state '
        STOP ' More electrons in L state than in R state '
      END IF
*. In addition to the operators dictated by the difference 
*. Between L and R, crea-anni pairs may be added 
      IAD_PAIRS = 0
*. Total number of creation and annihilation operators
      IF(IDEL_RL.LE.0) THEN
        NCREA_RL = IAD_PAIRS + ABS(IDEL_RL)
        NANNI_RL = IAD_PAIRS
      ELSE 
        NCREA_RL = IAD_PAIRS
        NANNI_RL = IAD_PAIRS + IDEL_RL
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCREA_RL, NANNI_RL = ', NCREA_RL, NANNI_RL
      END IF
*. Allocate space  
       CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
*. Pointers to KVEC1 and KVEC2, transferred through GLBBAS 
       KVEC1P = KVEC1
       KVEC2P = KVEC2
*
* ======================================
*. Set up list of types of excitations
* ======================================
*
C     GEN_OBEX(IFLAG,NCREA,NANNI,NOBEX_TP,IOBEX_TP)
*. Number of orbital excitations 
      CALL GEN_OBEX(1,NCREA_RL,NANNI_RL,NOBEX_TP,IDUMMY)
*. And the actual orbital excitations
      LEN = NOBEX_TP*2*NGAS
      CALL MEMMAN(KOBEX_TP,LEN,'ADDL  ',1,'OBEXTP')
      CALL MEMMAN(KLCOBEX_TP,NOBEX_TP,'ADDL  ',1,'COBEX ')
      CALL MEMMAN(KLAOBEX_TP,NOBEX_TP,'ADDL  ',1,'AOBEX ')
      CALL GEN_OBEX(0,NCREA_RL,NANNI_RL,NOBEX_TP,WORK(KOBEX_TP))
      CALL ISETVC(WORK(KLCOBEX_TP),NCREA_RL,NOBEX_TP)
      CALL ISETVC(WORK(KLAOBEX_TP),NANNI_RL,NOBEX_TP)
*. Number of types of spin-orbital excitations 
      IZERO = 0
      IAAEXC_TYP = 0
      IACT_SPC = 0
      MXSPOX = 0
      IPRCIX = 100
      CALL OBEX_TO_SPOBEX(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TP,NGAS,
     &     NOBPT,MS2_RL,IZERO    ,IAAEXC_TYP,IACT_SPC,IPRCIX,IDUMMY,
     &     MXSPOX,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAELEC_R,
     &     NBELEC_R,0)
*. And the actual spinorbital excitation operators 
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TP,'ADDL  ',1,'SPOBEX')
*. First SOX of given OX ( including zero operator )
      CALL MEMMAN(KIBSOX_FOR_OX,NOBEX_TP,'ADDL  ',1,'IBSOXF')
*. Number of SOX's for given OX
      CALL MEMMAN(KNSOX_FOR_OX,NOBEX_TP,'ADDL  ',1,'IBSOXF')
*. SOX for given OX
      CALL MEMMAN(KISOX_FOR_OX,NSPOBEX_TP,'ADDL  ',1,'IBSOXF')
      WRITE(6,*) ' NSPOBEX_TP = ', NSPOBEX_TP
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     NOBPT,MS2_RL,IZERO    ,IAAEXC_TYP,IACT_SPC,IPRCIX,
     &     WORK(KISOX_FOR_OX),
     &     MXSPOX,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAELEC_R,
     &     NBELEC_R,0)
       CALL MEMCHK
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
      IOCTPA = IBSPGPFTP(IATP_R)
      NOCTPA = NSPGPFTP(IATP_R)
      CALL LEN_GENOP_STR_MAP(
     &     NAOBEX_TP,WORK(KLAOBEX),NOCTPA,NELFSPGP(1,IOCTPA),
     &     NOBPT,NGAS,MAXLENA)
      IOCTPB = IBSPGPFTP(IBTP_R)
      NOCTPB = NSPGPFTP(IBTP_R)
      CALL LEN_GENOP_STR_MAP(
     &     NBOBEX_TP,WORK(KLBOBEX),NOCTPB,NELFSPGP(1,IOCTPB),
     &     NOBPT,NGAS,MAXLENB)
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(NTEST.GE.5) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
* Max Dimension of spinorbital excitation operators
*
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TP,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TP,'ADDL  ',1,'SPOBAC')
*
      MX_ST_TSOSO_MX = 0
      MX_ST_TSOSO_BLK_MX = 0
      MX_TBLK_MX = 0
      MX_TBLK_AS_MX = 0
      LEN_T_VEC_MX = 0
*
      IOPSM = MULTD2H(IGST_SM,ISSM)
*. Arrays not used
      KISOX_FOR_OCCLS = 1
      KIBSOX_FOR_OCCLS = 1
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,IOPSM,
     &       MX_ST_TSOSO_MX,MX_ST_TSOSO_BLK_MX,MX_TBLK_MX,
     &       WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &       0,MX_TBLK_AS_MX,
     &       WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &       NTCONF,IPRCIX)
      NDENSELMNT = LEN_T_VEC
      N_CC_AMP = LEN_T_VEC
      WRITE(6,*) ' Number of elements in general density=',
     &             NDENSELMNT
*. In the following we will construct density matrices 
*. using the lhs as specified by ICSM, ICSPC and IATP = 1, IBTP = 2
*. and the RHS specified by IATP_L, IBTP_L, in glstate.inc
* and IGST_OCC(MXPNGAS,2),IGST_SM,IGST_MS2 specified in crun
*. Jeppe, this could be made cleaner !!
*. The LHS state is assummed to be on LU17
*. Space for density matrices 
      CALL MEMMAN(KLT,NDENSELMNT,'ADDL  ',2,'TDENSI')
*. And transfer control and responsibility to the next level
      LULHS = 17
*. All blocks of T are active 
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TP)
*. Initialize densities to zero 
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLT),ZERO,NDENSELMNT)
*
      CALL VGSIGDEN_CC(WORK(KVEC1),WORK(KVEC2),LUC,LULHS,
     &                 WORK(KLT),ISIGDEN)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' The final density '
        WRITE(6,*) ' ================= '
        WRITE(6,*)
C            WRT_TVEC(ITSS_TP,LTSS_TP,NTSS_TP,T,ISM)
        CALL WRT_TVEC(WORK(KLSOBEX),WORK(KLLSOBEX),NSPOBEX_TP,
     &                WORK(KLT),IOPSM)
      END IF
*
      RETURN
      END
      SUBROUTINE VGSIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
*
* Outer routine for SIGMA and Density calculation for 
* general CC code
*
* Very General version : C and S space may have different 
* number of electrons ( and holes ...)
*
* For Sigma T is the input set of coefficients
* For Densi T is the output set of coefficients
*
* Jeppe Olsen, Helsinki in Sept. 2001 
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*
* =====
*.Input
* =====
*
*.Definition of c (RHS)
      INCLUDE 'cands.inc'
*. Definition of RHS 
      INCLUDE 'glstate.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cc_exc.inc'
*
      DIMENSION T(*)
      CALL LUCIAQENTER('SIGDE')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGDEN')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' SIGDEN_CC in action '
        WRITE(6,*)' Symmetry of left and right state = ',IGST_SM, ICSM
        WRITE(6,*) ' LUC and LUHC = ', LUC, LUHC
      END IF
*
      IATP_R = 1
      IBTP_R = 2
*
      NOCTPA_R = NOCTYP(IATP_R)
      NOCTPB_R = NOCTYP(IBTP_R)
*
      IOCTPA_R = IBSPGPFTP(IATP_R)
      IOCTPB_R = IBSPGPFTP(IBTP_R)
*
      NAEL_R = NELEC(IATP_R)
      NBEL_R = NELEC(IBTP_R)
*
      NOCTPA_L = NOCTYP(IATP_L)
      NOCTPB_L = NOCTYP(IBTP_L)
*
      IOCTPA_L = IBSPGPFTP(IATP_L)
      IOCTPB_L = IBSPGPFTP(IBTP_L)
*
      WRITE(6,*) ' IATP_L, IBTP_L = ', IATP_L, IBTP_L
*. Arrays giving allowed combinations of types in LHS ( called Sigma)
      CALL MEMMAN(KSIOIO,NOCTPA_L*NOCTPB_L,'ADDL  ',2,'SIOIO ')
C     IAIBCM_GAS(LCMBSPC,ICMBSPC,
C    &                      MNMXOC,NOCTPA,NOCTPB,IOCA,IOCB,NELFTP,
C    &                      MXPNGAS,NGAS,IOCOC,IPRNT,I_RE_MS2_SPACE,
C    &                      I_RE_MS2_VALUE)
      WRITE(6,*) ' NOCTPA_L, NOCTPB_L = ', NOCTPA_L, NOCTPB_L
      CALL IAIBCM_GAS(1,1,IGST_OCC,NOCTPA_L,NOCTPB_L,
     &            ISPGPFTP(1,IOCTPA_L),ISPGPFTP(1,IOCTPB_L),
     &            NELFGP,MXPNGAS,NGAS,WORK(KSIOIO),IPRCIX,
     &            I_RE_MS2_SPACE, I_RE_MS2_VALUE)
*. Arrays giving block type
      KSVST = 1
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL ZBLTP(ISMOST(1,IGST_SM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
*. Arrays for partitioning of S
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'SLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'SLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'SI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'SIBT  ')
*. Batches  of S vector
      LBLOCK = MAX(MXSOOB_AS,LCSBLK)
*. Well, SIGDEN uses full symmetry blocks of CI vectors
      ISIMSYM_LOC = 1
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP_L)),
     &     WORK(KNSTSO(IBTP_L)),NOCTPA_L,NOCTPB_L,NSMST,LBLOCK,
     &     WORK(KSIOIO),ISMOST(1,IGST_SM),
     &     NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),0,ISIMSYM_LOC)
*. Number of BLOCKS
      NBLOCK = IFRMR(WORK(KLSI1BT),1,NBATCH)
     &       + IFRMR(WORK(KLSLBT),1,NBATCH) - 1
C?    WRITE(6,*) ' Number of blocks ', NBLOCK
*. Start, length, ...  of  cc ampitudes 
      IOPSM = MULTD2H(IGST_SM,ICSM)
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,IOPSM,    
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MSCOMB_CC,MX_SBSTR,
     &     WORK(KISOX_FOR_OCCLS),NTOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &     NTCONF,IPRCC)
*. If combinations are in use, renormalize T-coefficients to 
*. spin-orbital normalization 
      IF(ISIGDEN.EQ.1.AND.MSCOMB_CC.EQ.1) 
     &CALL RENORM_T(WORK(KLSOBEX),WORK(KLLSOBEX),NSPOBEX_TP,
     &              T,IOPSM,2)
*. Well, SIGDEN uses full symmetry blocks 
      ISIMSYM_LOC = 1
      WRITE(6,*) ' NSPOBEX_TP before call to VGSIGDEN_CC2', NSPOBEX_TP
      CALL VGSIGDEN_CC2(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),LUC,LUHC,T,ISIGDEN,
     &     WORK(KVEC3))
*. Renormalize to combination form
      IF(MSCOMB_CC.EQ.1) 
     &CALL RENORM_T(WORK(KLSOBEX),WORK(KLLSOBEX),NSPOBEX_TP,
     &              T,IOPSM,1)
*. Eliminate local memory
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'SIGDEN')
*
      CALL LUCIAQEXIT('SIGDE')
*
      RETURN
      END
      SUBROUTINE VGSIGDEN_CC2(CB,SB,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &           LUC,LUHC,T,ISIGDEN,C2)
*
* <GSTATE! OP ! C >   = < Sigma ! Op ! C >
*
* First inner routine for general CC sigma/densi
*
* Jeppe Olsen   Modified SIGDEN_CC2, Sept. 2001 
*
* =====
* Input
* =====
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'glstate.inc'
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
*.Scratch
      DIMENSION SB(*),CB(*)
*
      CALL LUCIAQENTER('SIGD2')
      NTEST = 00
C     NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' SIGDEN2 speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*) ' LUC and LUHC ', LUC, LUHC
*
        IF(ISIGDEN.EQ.2) THEN
          WRITE(6,*) ' LHS on LUHC '
          CALL WRTVCD(CB,LUHC,1,-1)
        END IF
      END IF
*. We call the Left state Sigma, giving the S on the parameters
      IB_ATP_S = IBSPGPFTP(IATP_L)
      IB_BTP_S = IBSPGPFTP(IBTP_L)
      NOCTPA_S = NSPGPFTP(IATP_L)
      NOCTPB_S = NSPGPFTP(IBTP_L)
*
      IF(LUHC.GT.0) CALL REWINO(LUHC)
* Loop over batches over sigma blocks
      DO JBATS = 1, NBATS
*
        IB_SBAT = I1BATS(JBATS)
        L_SBAT  = LBATS(JBATS)
*
* Initialize sigma
*
* ISIGDEN = 1 : Set to zero
* ISIGDEN = 2 : Read in batch of blocks of lhs vector
*. Unique blocks are read in, expanded and normalized to determinant form
        DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
*
          ISATP = IBATS(1,ISBLK)
          ISBTP = IBATS(2,ISBLK)
          ISASM = IBATS(3,ISBLK)
          ISBSM = IBATS(4,ISBLK)
          ISATP_ABS = ISATP + IB_ATP_S - 1
          ISBTP_ABS = ISBTP + IB_BTP_S - 1
*
          ISOFF = IBATS(5,ISBLK)
          LEN2  = IBATS(7,ISBLK)
          IF(ISIGDEN.EQ.1) THEN
            ZERO = 0.0D0
            CALL SETVEC(SB(ISOFF),ZERO,LEN2)
          ELSE 
            XDUM = 0.0D0
            ISCALE = 1
*. Signs taken from CSTATE
            CALL GSTTBLD(SB(ISOFF),ISATP,ISASM,ISBTP,ISBSM,
     &            WORK(KNSTSO(IATP_L)),WORK(KNSTSO(IBTP_L)), 
     &            PSSIGN,IDC,PLSIGN,LUHC,C2,
     &            NSMST,ISCALE,XSCALE)
          END IF
        END DO
*       ^ End of loop over S-blocks in batch
*. Obtain sigma/density for batch of blocks
        CALL VGSIGDEN3(LBATS(JBATS),IBATS(1,I1BATS(JBATS)),1,
     &       CB,SB,LUC,T,ISIGDEN)
*
        IF(ISIGDEN.EQ.1) THEN
          IF(IDC.EQ.2) THEN
*. Determinant => combination form and scale
*. reform 
           CALL RFTTS(SB,C2,IBATS(1,IB_SBAT),L_SBAT,
     &                1,NSMST,NOCTPA_S,NOCTPB_S,
     &                WORK(KNSTSO(IATP_L)), WORK(KNSTSO(IBTP_L)),
     &                IDC,PSSIGN,1,NTEST)
*. scale
           CALL SCDTTS(SB,IBATS(1,IB_SBAT),L_SBAT,NSMST,NOCTPA,NOCTPB,
     &                WORK(KNSTSO(IATP_L)), WORK(KNSTSO(IBTP_L)),
     &                IDC,1,NTEST)
          END IF
*. Transfer packed S block to permanent storage
          DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
            IOFF = IBATS(6,ISBLK)
            LEN  = IBATS(8,ISBLK)
            CALL ITODS(LEN,1,-1,LUHC)
            CALL TODSC(SB(IOFF),LEN,-1,LUHC)
          END DO
        END IF
*
      END DO
*     ^ End of loop over S-batches
*
      IF(ISIGDEN.EQ.1) CALL ITODS(-1,1,-1,LUHC)
*
      CALL LUCIAQEXIT('SIGD2')
      RETURN
      END
      SUBROUTINE VGSIGDEN3(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,T,ISIGDEN)
*
* Generate a set of sigma blocks, 
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order 
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      INCLUDE 'cands.inc'
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'glstate.inc'
      COMMON/CMXCJ/MXCJ
*
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      INCLUDE 'cintfo.inc'
*
      DIMENSION CB(*), HCB(*)
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGDE3')
*
      NTEST = 00
      IF(NTEST.GE.5)
     &WRITE(6,*) ' SIGDE3 : ISSPC,ICSPC ', ISSPC,ICSPC
*. MXINKA_CC is pt used for length of batch of strings and Toccs
*
*.Local scratch arrays for blocks of C and sigma
      LSCR1 = MXSOOB_AS
      LSCR1 = MAX(LSCR1,LCSBLK)
*.T-elements                       
      INTSCR = MIN(MXINKA_CC**4, MX_TBLK_MX)
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*
*. Offsets for alpha and beta supergroups in C
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM(ICSPC,WORK(KCIOIO))
*
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + MXINKA_CC**4
*
*. Arrays for storing NEL consecutive annihilations/creations
      LSCR3 = MAXLEN_I1
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
*. Arrays for storing occupations of T-operators
      LENNY = MX_ST_TSOSO_BLK_MX
      CALL MEMMAN(KTOCC1,LENNY,'ADDL  ',2,'TOCC1 ')
      CALL MEMMAN(KTOCC2,LENNY,'ADDL  ',2,'TOCC2 ')
      CALL MEMMAN(KTOCC3,LENNY,'ADDL  ',2,'TOCC3 ')
      CALL MEMMAN(KTOCC4,LENNY,'ADDL  ',2,'TOCC4 ')
*. Arrays for storing NEL consecutive annihilations/creations
      LSCR4 = MXINKA_CC*MXINKA_CC*MXINKA_CC
      CALL MEMMAN(KI1G  ,LSCR4,'ADDL  ',1,'I1G   ')
      CALL MEMMAN(KXI1G,LSCR4,'ADDL  ',2,'XIG  ')
*
      CALL MEMMAN(KI2G  ,LSCR4,'ADDL  ',1,'I2G   ')
      CALL MEMMAN(KXI2G,LSCR4,'ADDL  ',2,'XI2G  ')
*
      CALL MEMMAN(KI3G  ,LSCR4,'ADDL  ',1,'I3G   ')
      CALL MEMMAN(KXI3G,LSCR4,'ADDL  ',2,'XI3G  ')
*
      CALL MEMMAN(KI4G  ,LSCR4,'ADDL  ',1,'I4G   ')
      CALL MEMMAN(KXI4G,LSCR4,'ADDL  ',2,'XI4G  ')
*
      CALL MEMMAN(KI1GE  ,LSCR4,'ADDL  ',1,'I1GE  ')
      CALL MEMMAN(KXI1GE,LSCR4,'ADDL  ',2,'XIGE ')
*
      CALL MEMMAN(KI2GE  ,LSCR4,'ADDL  ',1,'I2GE  ')
      CALL MEMMAN(KXI2GE,LSCR4,'ADDL  ',2,'XI2GE ')
*
      CALL MEMMAN(KI3GE  ,LSCR4,'ADDL  ',1,'I3GE  ')
      CALL MEMMAN(KXI3GE,LSCR4,'ADDL  ',2,'XI3GE ')
*
      CALL MEMMAN(KI4GE  ,LSCR4,'ADDL  ',1,'I4GE  ')
      CALL MEMMAN(KXI4GE,LSCR4,'ADDL  ',2,'XI4GE ')
*.TTS arrays  for partitioning of vector 
      NTTS = MXNTTS
      CALL MEMMAN(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      CALL MEMMAN(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      CALL MEMMAN(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      CALL MEMMAN(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      CALL MEMMAN(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')
*. Arrays giving block type
COLD  CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      KSVST = 1
*. Find batches of C - strings
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*. Use all symmetryblocks of given TT
      ISIMSYM_LOC = 1
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &              WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &              WORK(KCIOIO),ISMOST(1,ICSM),NCBATCH,WORK(KLLBT),
     &              WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &              0,ISIMSYM_LOC)
*. Space for four blocks of string occupations and arrays of 
*. reordering arrays
*. Also used to hold an NORB*NORB matrix  
*. MAx for left and right ?
*
      NAEL_L = NELFTP(IATP_L)
      NBEL_L = NELFTP(IBTP_L)
      NAEL_MAX = MAX(NAEL,NAEL_L)
      NBEL_MAX = MAX(NBEL,NBEL_L) 
      NEL_MAX = MAX(NAEL_MAX,NBEL_MAX)
*
      LZSCR = (NEL_MAX+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (NEL_MAX+2) * NOCOB
*. Set up to two blocks for orbital conserving operator 
      DO I1234 = 1, 2
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')         
      END DO
   
      DO I1234 = 1, 2
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',2,'KLZSCR')
*
      IOCTPA_L = IBSPGPFTP(IATP_L)
      IOCTPB_L = IBSPGPFTP(IBTP_L)
      CALL VGSIGDEN4(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KVEC3),
     &     WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &     WORK(KNSTSO(IATP_L)),WORK(KNSTSO(IBTP_L)),
     &     IATP,IBTP,IOCTPA,IOCTPB,
     &     IATP_L,IBTP_L,IOCTPA_L,IOCTPB_L,
     &     MXINKA_CC,WORK(KINSCR),
     &     WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &     WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),
     &     IPRDIA,LUC,WORK(KCJRES),WORK(KSIRES),
     &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &     ISSM,ICSM,ISIGDEN,NCBATCH,T,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),
     &     WORK(KLSPOBEX_AC),WORK(KLLSOBEX),
     &     WORK(KTOCC1),WORK(KTOCC2), WORK(KTOCC3),WORK(KTOCC4),
     &     WORK(KI1G),WORK(KXI1G), WORK(KI2G),WORK(KXI2G),
     &     WORK(KI3G),WORK(KXI3G), WORK(KI4G),WORK(KXI4G),
     &     WORK(KI1GE),WORK(KXI1GE), WORK(KI2GE),WORK(KXI2GE),
     &     WORK(KI3GE),WORK(KXI3GE), WORK(KI4GE),WORK(KXI4GE),
     &     WORK(KLSPOBEX_FRZ))
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SIGDE3')
      RETURN
      END
      SUBROUTINE VGSIGDEN4(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           NSSOA_C,NSSOB_C,NSSOA_S,NSSOB_S,
     &           IAGRP_C,IBGRP_C,
     &           IOCTPA_C,IOCTPB_C,
     &           IAGRP_S,IBGRP_S,
     &           IOCTPA_S,IOCTPB_S,
     &           MAXK,
     &           XINT,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           IPRNT,LUC,CJRES,SIRES,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,
     &           ISSM,ICSM,ISIGDEN,NCBATCH,
     &           T,NSPOBEX_TP,ITSPOBEX_TP,IBSPOBEX_TP,
     &           ISPOBEX_AC,LSPOBEX_TP,
     &           TOCC1,TOCC2,TOCC3,TOCC4,
     &           I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &           I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,XI4GE,
     &           ISPOBEX_FRZ)
*. IAGRP_C : Type of alpha strings in C
*. IOCTPA_C : First supergroup of alpha strings in C
*
*
* Jeppe Olsen, September 2001 from GENSIG4
*                                 
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included 
*
* NSSOA : Number of strings per type and symmetry for alpha strings
* NSSOB : Number of strings per type and symmetry for beta strings
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
c      IMPLICIT REAL*8(A-H,O-Z)
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strbas.inc'
*. Specific input
      INTEGER ISBLOCK(8,*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*.General input
      INTEGER NSSOA_C(NSMST ,*), NSSOB_C(NSMST ,*)
      INTEGER NSSOA_S(NSMST ,*), NSSOB_S(NSMST ,*)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*)
      DIMENSION  XI1S(*),XI2S(*),XI3S(*)
*
      DIMENSION CJRES(*),SIRES(*)
*. T-coefficients
      DIMENSION T(*), ITSPOBEX_TP(4*NGAS,*), IBSPOBEX_TP(*)
      INTEGER  ISPOBEX_AC(*),ISPOBEX_FRZ(*), LSPOBEX_TP(*)
*
*
      CALL LUCIAQENTER('SIGD4')
      IPRNT_ORIG = IPRNT
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' VGSIGDEN4 speaking :'
        WRITE(6,*) ' ===================='
        WRITE(6,*)  
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS and IOFF for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0) 
     &      WRITE(6,'(10X,I6,4I3,I8)') 
     &      IBLOCK,(ISBLOCK(II,IBLOCK),II=1,5)
          END DO
        WRITE(6,*) ' List af active T blocks : '
        CALL IWRTMA(ISPOBEX_AC,1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
*
      NSTT_BLK = NSBLOCK/NSMST
* Loop over batches over C blocks      
      REWIND LUC
*
      DO 20000 JCBATCH = 1, NCBATCH             
*
*. Read C blocks into core
*
        NJBLOCK = LCBLOCK(JCBATCH)
        NCTT_BLK = NJBLOCK/NSMST
        I1C = I1CBLOCK(JCBATCH)
*
*. Loop over TT blocks of sigma and C in batches and 
*  obtain  contributions 
        DO 9000 ICTT_BLK = 1, NCTT_BLK
          IREAD = 0
*. first block of next TT block of C    
          ICBLK = I1C + (ICTT_BLK-1)*NSMST
          JJATP = ICBLOCK(1,ICBLK)
          JJBTP = ICBLOCK(2,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Next block of C, ICBLK,JJATP,JJBTP,ICOFF ',
     &      ICBLK,JJATP,JJBTP, ICOFF
          END IF
*. Number of permutations of C-block
          NCPERM = 1
          IF(IDC.EQ.2.AND.JJATP.NE.JJBTP) NCPERM = 2
          IF(NTEST.GE.100) WRITE(6,*) ' NCPERM = ', NCPERM
*. Diagonal block
          ICDIAG = 0
          IF(IDC.EQ.2.AND.JJATP.EQ.JJBTP) ICDIAG = 1
          DO 8999 ICPERM = 1, NCPERM
          IF(NTEST.GE.100) WRITE(6,*) ' ICPERM = ', ICPERM
*
          IF(ICPERM.EQ.1) THEN
            JATP = JJATP
            JBTP = JJBTP
          ELSE
            JATP = JJBTP
            JBTP = JJATP
          END IF
          ICTRANSPOSED = 0
*
          IRATP = JATP + IOCTPA_C - 1
          IRBTP = JBTP + IOCTPB_C - 1
*
          DO 10000 ISTT_BLK = 1, NSTT_BLK
*. first block of next TT block of Sigma
            ISBLK = (ISTT_BLK-1)*NSMST + 1
            ISOFF = ISBLOCK(5,ISBLK)
*
            IIATP = ISBLOCK(1,ISBLK)
            IIBTP = ISBLOCK(2,ISBLK)
*
            NSPERM = 1
            IF(ISIGDEN.EQ.2.AND.IDC.EQ.2.AND.IIATP.NE.IIBTP) NSPERM = 2
            DO 9999 ISPERM = 1, NSPERM
*
            IF(ISPERM.EQ.1) THEN
              IATP = IIATP
              IBTP = IIBTP
            ELSE
              IATP = IIBTP
              IBTP = IIATP
            END IF
            IF(ISPERM.EQ.2) THEN
*. Transpose S- blocks
               CALL TRP_CITT_BLK(SB(ISOFF),C2,ISSM,
     &              NSSOA_S(1,IIATP),NSSOB_S(1,IIBTP),NSMST,1)
            END IF
*
            IF(IDC.EQ.2.AND.IATP.EQ.IBTP) THEN 
              ISRESTRICT = 1
            ELSE 
              ISRESTRICT = 0
            END IF
            ILATP = IATP + IOCTPA_S - 1
            ILBTP = IBTP + IOCTPB_S - 1
*. Connections ?
             DO 8000 ITTP = 1, NSPOBEX_TP
*
              ICA_OFF = 1
              ICB_OFF = 1 +  NGAS
              IAA_OFF = 1 +2*NGAS
              IAB_OFF = 1 +3*NGAS
*
              IF(MSCOMB_CC.EQ.1) THEN
*. Check if given T-block corresponds to one or two blocks
                CALL DIAG_EXC_CC(ITSPOBEX_TP(1+0*NGAS,ITTP),
     &               ITSPOBEX_TP(1+1*NGAS,ITTP),
     &               ITSPOBEX_TP(1+2*NGAS,ITTP),
     &               ITSPOBEX_TP(1+3*NGAS,ITTP),NGAS,IDIAG)
                IF(IDIAG.EQ.1.OR.ISIGDEN.EQ.2) THEN
                 NPBLK = 1
                ELSE
                 NPBLK = 2
                END IF
              ELSE
                NPBLK = 1
                IDIAG = 0
              END IF
              DO IPBLK = 1, NPBLK
*
               IF(IPBLK.EQ.1) THEN
                 ITABTRNSP = 0
               ELSE
                 ITABTRNSP = 1
               END IF
*
               WRITE(6,*) ' Spinorbitalexcitationtype=',ITTP
               CALL WRT_TP_GENOP(ITSPOBEX_TP(ICA_OFF,ITTP),
     &                           ITSPOBEX_TP(ICB_OFF,ITTP),
     &                           ITSPOBEX_TP(IAA_OFF,ITTP),
     &                           ITSPOBEX_TP(IAB_OFF,ITTP),NGAS)
               IF(IPBLK.EQ.1) THEN
*.   untransposed block
                 LCA_OFF  = ICA_OFF
                 LAA_OFF  = IAA_OFF
                 LCB_OFF  = ICB_OFF
                 LAB_OFF  = IAB_OFF
               ELSE
*.  Alpha-beta transposed block
                 LCB_OFF  = ICA_OFF
                 LAB_OFF  = IAA_OFF
                 LCA_OFF  = ICB_OFF
                 LAA_OFF  = IAB_OFF
               END IF
*. Connections ?
               CALL GXFSTR(NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),
     &         NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &         ITSPOBEX_TP(LCA_OFF,ITTP), ITSPOBEX_TP(LAA_OFF,ITTP),
     &         ITSPOBEX_TP(LCB_OFF,ITTP), ITSPOBEX_TP(LAB_OFF,ITTP),
     &                     NGAS,ICON)
               IF(ICON.EQ.1.AND.ISPOBEX_AC(ITTP).EQ.1.AND.
     &            .NOT.(ISPOBEX_FRZ(ITTP).EQ.1.AND.ISIGDEN.EQ.2)) THEN
*. Read block in   
                 IF(IREAD.EQ.0) THEN
*. Number of blocks
                   NCBLK = NSMST
                   IF(IDC.EQ.2.AND.ICSM.NE.1.AND.JATP.EQ.JBTP) 
     &             NCBLK = NSMST/2
*. And fetch the c-blocks and transform into unpacked blocks in
*. determinant normalization
                   JATP_ABS = JATP +  IOCTPA_C - 1
                   JBTP_ABS = JBTP +  IOCTPB_C - 1
                   DO JBLK = ICBLK, ICBLK -1 + NCBLK            
                     JASM = ICBLOCK(3,JBLK)
                     JBSM = ICBLOCK(4,JBLK)
                     JOFF = ICBLOCK(5,JBLK)
                     ISCALE = 1
                     CALL GSTTBLD(CB(JOFF),JJATP,JASM,JJBTP,JBSM,
     &                    WORK(KNSTSO(IAGRP_C)),WORK(KNSTSO(IBGRP_C)), 
     &                    PSSIGN,IDC,PLSIGN,
     &                    LUC,C2,NSMST,ISCALE,XSCALE)
*
                     IF(NTEST.GE.200) THEN
                       WRITE(6,*) ' C block read in, JBLK, JOFF : ',
     &                 JBLK,JOFF
                       LEN = ICBLOCK(7,JBLK)
                       CALL WRTMAT(CB(JOFF),1,LEN,1,LEN)
                     END IF
*
                    END DO
                   IREAD = 1
                 END IF
*
                 IF(ICPERM.EQ.2.AND.ICTRANSPOSED.EQ.0) THEN
*. Transpose 
                    CALL TRP_CITT_BLK(CB(ICOFF),C2,ICSM,
     &                   NSSOA_C(1,JJATP),NSSOB_C(1,JJBTP),NSMST,1)
                    ICTRANSPOSED = 1
                 END IF
                 IF(IDC.EQ.2.AND.JJATP.EQ.JJBTP) THEN 
                   ICRESTRICT = 1
                 ELSE 
                   ICRESTRICT = 0
                 END IF
*
* ICPERM = 2 : All permutations have been obtained by redefining IRATP ...
                 ITB = IBSPOBEX_TP(ITTP)
                 CALL VGNSIDE(ISIGDEN,
     &           ITSPOBEX_TP(LCA_OFF,ITTP),ITSPOBEX_TP(LCB_OFF,ITTP),
     &           ITSPOBEX_TP(LAA_OFF,ITTP),ITSPOBEX_TP(LAB_OFF,ITTP),
     &           NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),            
     &           NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &           NSSOA_S(1,IATP),NSSOB_S(1,IBTP),
     &           NSSOA_C(1,JATP),NSSOB_C(1,JBTP),
     &           T(ITB), 
     &           SB(ISOFF),CB(ICOFF),ISSM,ICSM,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           XINT,CJRES,SIRES,MAXK,
     &           TOCC1,TOCC2,TOCC3,TOCC4,
     &           I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &           I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,XI4GE,
     &           IDIAG,ITABTRNSP,ICRESTRICT,ISRESTRICT)
               END IF
*              ^ End if connection 
              END DO 
*             ^ End of loop over permutation of T-block
              IF(NTEST.GE.200) THEN
                WRITE(6,*) ' Updated density block '
                WRITE(6,*) ' ===================== '
                L = LSPOBEX_TP(ITTP)
                WRITE(6,*) ' Offset and length ', ITB, L
                CALL WRTMAT(T(ITB),1,L,1,L)
              END IF
 8000       CONTINUE
*.          ^  End of loop over TT sigma blocks 
            IF(ISPERM.EQ.2) THEN
*. Back transpose S- blocks
               CALL TRP_CITT_BLK(SB(ISOFF),C2,ISSM,
     &              NSSOB_S(1,IIBTP),NSSOA_S(1,IIATP),NSMST,1)
            END IF
 9999     CONTINUE
*         ^ End of loop over S-permutations
          IF(NTEST.GE.200) THEN
            IF(ISIGDEN.EQ.1) THEN
              WRITE(6,*) ' Updated sigma blocks ' 
              CALL WRTVH1(SB(ISOFF),ISSM,NSSOA_S(1,IIATP),
     &                    NSSOB_S(1,IIBTP),NSMST,0) 
            END IF
          END IF
10000     CONTINUE
*         ^ End of loop over blocks of sigma 
 8999   CONTINUE
*.      ^ End of loop over Permutations of C-block
        IF(IREAD.EQ.0) THEN
*. These C blocks were not needed, skip them
          DO JBLK = 1, NSMST
            CALL IFRMDS(LBL,-1,1,LUC)
            CALL SKPRCD2(LBL,-1,LUC)
          END DO
        END IF
*
 9000   CONTINUE
*.      ^ End of loop over TT C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
*
      CALL LUCIAQEXIT('SIGD4')
      RETURN
      END
      SUBROUTINE VGNSIDE(ISD,ICA,ICB,IAA,IAB,
     &           IAOC,IBOC,JAOC,JBOC,
     &           NIA,NIB,NJA,NJB,
     &           TB,SB,CB,ISSM,ICSM,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           TBSUB,CJRES,SIRES,MAXLB,
     &           ICA_STR,ICB_STR,IAA_STR,IAB_STR,
     &           KJA,XKJA,KJB,XKJB,KIA,XKIA,KIB,XKIB,
     &           KJAD,XKJAD,KJBD,XKJBD,KIAD,XKIAD,KIBD,XKIBD,
     &           ITDIAG,ITABTRNSP,IC_RESTRICT,IS_RESTRICT)
*
* Sigma routine and density routine for for general operator 
*
* sum(Ica,Icb,Iaa,Iab) T(Ica,Icb,Iaa,Iab) Ica+ Icb+ Iaa Iab
*
* contribution from C vector with occupation JAOC,IBOC
* to Sigma vector with occupation IAOC, IBOC
*
* All symmetryblocks of C and Sigma are treated 
*
* ITDIAG = 1 => Diagonal T-block 
*               Only symmetryblocks with (ICASM,IAASM).GE.(ICBSM,IABSM)
*               are included, and blocks with (ICASM,IAASM) = (ICBSM,IABSM)
*               are given in diagonal form 
*               (Iaa.ge.Iab, if Iaa=Iab, Ica .ge. Icb)
* ITABTRSNP = 1 => alpha-beta transposed block, input block is 
*               T(ICB,ICA,IAB,IAA)
* IX_RESTRICT = 1: BLocks in X are restricted to those with IASM.GE.IBSM
* IX_RESTRICT =-1: BLocks in X are restricted to those with IASM.LT.IBSM
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
*
* Jeppe Olsen, From GNSIDE September 2001
*

*
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc' 
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'newccp.inc'
*     
*.Input
      DIMENSION CB(*)
      INTEGER  NIA(*), NIB(*), NJA(*), NJB(*)
      INTEGER IAOC(*),IBOC(*),JAOC(*),JBOC(*)
      INTEGER ICA(*), ICB(*),IAA(*),IAB(*)
*.Input or Output
      DIMENSION SB(*), TB(*)
*.Local Scratch
      INTEGER KAOCX(MXPNGAS),KBOC(MXPNGAS)
      INTEGER NIKAINTM(MXPOBS,MXPLCCOP),NIKBINTM(MXPOBS,MXPLCCOP)
      INTEGER NJKAINTM(MXPOBS,MXPLCCOP),NJKBINTM(MXPOBS,MXPLCCOP)
      INTEGER NKA(MXPOBS), NKB(MXPOBS)
*. 
      INTEGER NICA(MXPOBS),NICB(MXPOBS),NIAA(MXPOBS),NIAB(MXPOBS)
      INTEGER ICA_EXP(100),ICB_EXP(100),IAA_EXP(100),IAB_EXP(100)
      INTEGER ICAGP(MXPNGAS),ICBGP(MXPNGAS)
      INTEGER IAAGP(MXPNGAS),IABGP(MXPNGAS)
      INTEGER IAC_AR(100)
      INTEGER IB_CA(8,8,MXPLCCOP), IB_CB(8,8,MXPLCCOP)
      INTEGER IB_AA(8,8,MXPLCCOP), IB_AB(8,8,MXPLCCOP)
      INTEGER IB_T(8,8,8)
      INTEGER IB_C(8),IB_S(8)
*. Scratch through input 
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION TBSUB(*)
      INTEGER KJA(MAXLB),KJB(MAXLB),KIA(MAXLB),KIB(MAXLB)
      DIMENSION XKJA(MAXLB),XKJB(MAXLB),XKIA(MAXLB),XKIB(MAXLB)
      INTEGER KJAD(MAXLB),KJBD(MAXLB),KIAD(MAXLB),KIBD(MAXLB)
      DIMENSION XKJAD(MAXLB),XKJBD(MAXLB),XKIAD(MAXLB),XKIBD(MAXLB)
*. ^Must hold excitations for all intermediate strings of given sym 
*. and all orbitals of given type
      INTEGER ICA_STR(*),ICB_STR(*),IAA_STR(*),IAB_STR(*)
*     ^ Must hold strings of given sym for T ops.
      DIMENSION CJRES(*),SIRES(*)
*
      COMMON/KKKDUMMY/LEN_C2,LEN_S2
      COMMON/CMXCJ/MXCJ
     
      CALL LUCIAQENTER('GENSIG')
*
      NTEST = 500 
C     NTEST = MAX(NTESTL,NTESTG)
*
      ISKIP = 0
      LEN_C = LEN_TT_BLOCK(ICSM,NJA,NJB,NSMST)
      XCNORM = INPROD(CB,CB,LEN_C)
      IF(XCNORM.EQ.0.0D0) ISKIP = 1
      IF(ISKIP.EQ.1) GOTO 9999
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' GNSIDE speaking '
        WRITE(6,*) ' ================ '
        WRITE(6,*)
        WRITE(6,*) ' IAOC and IBOC in GENSIG'
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' JAOC and JBOC in GENSIG'
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Type of operator '
        CALL WRT_TP_GENOP(ICA,ICB,IAA,IAB,NGAS)
        IF(ISD.EQ.1) THEN
          WRITE(6,*) ' Sigma generation '
        ELSE 
          WRITE(6,*) ' Density generation '
       END IF
*
       WRITE(6,*) ' Input C vector '
C              WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
       CALL WRTVH1(CB,ICSM,NJA,NJB,NSMST,0)
       WRITE(6,*) ' Input S vector '
       CALL WRTVH1(SB,ISSM,NIA,NIB,NSMST,0)
*
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Symmetry of T
      IOPSM = MULTD2H(ICSM,ISSM)
*. Offsets for symmetryblocks of C and S
*. Total length of C and sigma blocks
      LEN_S = LEN_TT_BLOCK(ISSM,NIA,NIB,NSMST)
*
      ISKIP = 0
      IF(ISD.EQ.2) THEN
       XSNORM = INPROD(SB,SB,LEN_S)
       IF(XSNORM.EQ.0.0D0) ISKIP = 1
      END IF
      IF(ISKIP.EQ.1) GOTO 9999
*. Offset of C and sigma routines
      CALL IOFF_SYMBLK_MAT(NSMST,NJA,NJB,ICSM,IB_C,IC_RESTRICT)
      CALL IOFF_SYMBLK_MAT(NSMST,NIA,NIB,ISSM,IB_S,IS_RESTRICT)
*. types of strings in T in groupnotation
      CALL OCC_TO_GRP(ICA,ICAGP,1)
      CALL OCC_TO_GRP(ICB,ICBGP,1)
      CALL OCC_TO_GRP(IAA,IAAGP,1)
      CALL OCC_TO_GRP(IAB,IABGP,1)
C          OCC_TO_GRP(IOCC,IGRP,IWAY)
*
C     XCNORM = INPROD(CB,CB,LEN_C)
C     IF(XCNORM.EQ.0.0D0) GOTO 9999
*. Operators in T in expanded form
C          REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
      CALL REF_OP(ICA,ICA_EXP,NCA_OP,NGAS,1)
      CALL REF_OP(ICB,ICB_EXP,NCB_OP,NGAS,1)
      CALL REF_OP(IAA,IAA_EXP,NAA_OP,NGAS,1)
      CALL REF_OP(IAB,IAB_EXP,NAB_OP,NGAS,1)
*. Sign for bringing operators and amplitudes into new order
      NPERM = NAA_OP*NCB_OP + NAA_OP*(NAA_OP-1)/2+NAB_OP*(NAB_OP-1)/2
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNXXX  =  -1.0D0
      ELSE
        SIGNXXX = 1.0D0
      END IF
      IF(I_USE_NEWCCP .EQ. 0) SIGNXXX = 1.0D0
*. Type of Ka and Kb
      CALL CCEX_OCC_OCC(JAOC,KAOCX,NGAS,1,IAA,IKA_ZERO)
      CALL CCEX_OCC_OCC(JBOC,KBOC,NGAS,1,IAB,IKB_ZERO)
      IF(IKA_ZERO.EQ.0.AND.IKB_ZERO.EQ.0) THEN
*  NST_SPGP(IOCC,NSTFSM)
       CALL NST_SPGP(KAOCX,NKA)
       CALL NST_SPGP(KBOC,NKB)
*. Number of strings in excitation ops
       CALL NST_SPGP(ICA,NICA)
       CALL NST_SPGP(ICB,NICB)
       CALL NST_SPGP(IAA,NIAA)
       CALL NST_SPGP(IAB,NIAB)
*. Offset for symmetryblocks of T
      IF(ITABTRNSP.EQ.0) THEN
        CALL Z_TCC_OFF(IB_T,NICA,NICB,NIAA,NIAB,IOPSM,NSMST,ITDIAG)
      ELSE
        CALL Z_TCC_OFF(IB_T,NICB,NICA,NIAB,NIAA,IOPSM,NSMST,ITDIAG)
      END IF
C          Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
*. Creation path from C/S strings to K strings
*. Mappings from !I> to !K> are annihilations mappings 
       NMXOP_T = MAX(NCA_OP,NCB_OP,NAA_OP,NAB_OP)
       DO JOP = 1, NMXOP_T
         IAC_AR(JOP) = 1
       END DO
       ONE = 1.0D0
* CA : KA => IA mappings in 1
C           MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SX,NK,IB,SCLFAC)
       CALL MAP_EXSTR(ICA_EXP,IAC_AR,NCA_OP,IAOC,I1,XI1S,
     &               NIKAINTM,IB_CA,1.0D0)
* CB : KB => IB     
       CALL MAP_EXSTR(ICB_EXP,IAC_AR,NCB_OP,IBOC,I2,XI2S,
     &               NIKBINTM,IB_CB,ONE)
* AA : KA => JA
       CALL MAP_EXSTR(IAA_EXP,IAC_AR,NAA_OP,JAOC,I3,XI3S,
     &            NJKAINTM,IB_AA,ONE)    
* AB : KB => JB
       CALL MAP_EXSTR(IAB_EXP,IAC_AR,NAB_OP,JBOC,I4,XI4S,
     &            NJKBINTM,IB_AB,ONE)    
*. Loop over symmetry-blocks of C and Sigma
       DO IASM = 1, NSMST 
        IBSM = MULTD2H(ISSM,IASM)
*
        IS_INCLUDED = 1
        IF(IS_RESTRICT.EQ.1.AND.IASM.LT.IBSM) IS_INCLUDED = 0
        IF(IS_RESTRICT.EQ.-1.AND.IASM.GE.IBSM) IS_INCLUDED = 0
*
        IF(IS_INCLUDED.EQ.1) THEN
        DO JASM = 1, NSMST
         JBSM = MULTD2H(ICSM,JASM)
*
         IC_INCLUDED = 1
         IF(IC_RESTRICT.EQ.1.AND.JASM.LT.JBSM) IC_INCLUDED = 0
         IF(IC_RESTRICT.EQ.-1.AND.JASM.GE.JBSM) IC_INCLUDED = 0
*
         IF(IC_INCLUDED.EQ.1) THEN
         IF(NTEST.GE.1000) WRITE(6,'(A,4I4)') 
     &   ' IASM, IBSM, JASM, JBSM ', IASM,IBSM,JASM,JBSM
         DO KASM = 1, NSMST
*
          IAA_SM = MULTD2H(KASM,JASM)
          ICA_SM = MULTD2H(IASM,KASM)
          LIAA_SM = NIAA(IAA_SM)
          LICA_SM = NICA(ICA_SM)
          LKA_SM = NKA(KASM)
*
          DO KBSM = 1, NSMST
*
           IAB_SM = MULTD2H(KBSM,JBSM)
           ICB_SM = MULTD2H(IBSM,KBSM)
*
           LIAB_SM = NIAB(IAB_SM)
           LICB_SM = NICB(ICB_SM)
*
           LKB_SM = NKB(KBSM)
*
           ICBOFF = IB_C(JASM)
           ISBOFF = IB_S(IASM)
*
*
           IF(LIAA_SM*LICA_SM*LIAB_SM*LICB_SM*LKA_SM*LKB_SM.NE.0)
     &     THEN 
*
             IF(NTEST.GE.100) THEN
               WRITE(6,*) ' ICA_SM, ICB_SM, IAA_SM,IAB_SM =',
     &                      ICA_SM, ICB_SM, IAA_SM,IAB_SM
             END IF
*
             IF(NTEST.GE.1000) 
     &       WRITE(6,'(A,2I4)') ' KASM,KBSM ', KASM, KBSM
             IF(NTEST.GE.1000) 
     &       WRITE(6,*) ' LKA_SM, LKB_SM', LKA_SM, LKB_SM
             IF(NTEST.GE.1000) 
     &       WRITE(6,'(A,4I3)')      
     &       'LIAA_SM, LICA_SM,LIAB_SM,LICB_SM',
     &        LIAA_SM, LICA_SM,LIAB_SM,LICB_SM
*
*. The operator strings of ICA, ICB, IAA, IAB 
            CALL GETSTR2_TOTSM_SPGP(ICAGP,NGAS,ICA_SM,NCA_OP,
     &           LICA_SM2,ICA_STR,NORBT,0,0,0)
            CALL GETSTR2_TOTSM_SPGP(ICBGP,NGAS,ICB_SM,NCB_OP,
     &           LICB_SM2,ICB_STR,NORBT,0,0,0)
            CALL GETSTR2_TOTSM_SPGP(IAAGP,NGAS,IAA_SM,NAA_OP,
     &           LIAA_SM2,IAA_STR,NORBT,0,0,0)
            CALL GETSTR2_TOTSM_SPGP(IABGP,NGAS,IAB_SM,NAB_OP,
     &           LIAB_SM2,IAB_STR,NORBT,0,0,0)
* Batching                   
           LEN_KA_BAT = MAXLB
           LEN_KB_BAT = MAXLB
*
           IF(LKA_SM.LT.MAXLB) THEN
             LEN_KA_BAT = LKA_SM
             LEN_KB_BAT = MIN(LKB_SM,MAXLB*MAXLB/LKA_SM)
           END IF
*
           IF(LKB_SM.LT.MAXLB) THEN
             LEN_KB_BAT = LKB_SM
             LEN_KA_BAT = MIN(LKA_SM,MAXLB*MAXLB/LKB_SM)
           END IF
*
           NKA_BT = LKA_SM/LEN_KA_BAT
           IF(NKA_BT*LEN_KA_BAT. LT. LKA_SM) NKA_BT = NKA_BT + 1
           NKB_BT = LKB_SM/LEN_KB_BAT
           IF(NKB_BT*LEN_KB_BAT. LT. LKB_SM) NKB_BT = NKB_BT + 1
           
*
           NCA_BT = LICA_SM/MAXLB
           IF(NCA_BT*MAXLB.LT.LICA_SM) NCA_BT = NCA_BT + 1
           NCB_BT = LICB_SM/MAXLB
           IF(NCB_BT*MAXLB.LT.LICB_SM) NCB_BT = NCB_BT + 1
           NAA_BT = LIAA_SM/MAXLB
           IF(NAA_BT*MAXLB.LT.LIAA_SM) NAA_BT = NAA_BT + 1
           NAB_BT = LIAB_SM/MAXLB
           IF(NAB_BT*MAXLB.LT.LIAB_SM) NAB_BT = NAB_BT + 1
*. Loop over batches
           DO IKA_BT = 1, NKA_BT
           DO IKB_BT = 1, NKB_BT
            IKA_B = (IKA_BT-1)*LEN_KA_BAT + 1
            IKA_E = MIN(LKA_SM,IKA_B + LEN_KA_BAT - 1)
            NKA_B = IKA_E - IKA_B + 1
            NKA_B0 = NKA_B
            IKB_B = (IKB_BT-1)*LEN_KB_BAT + 1
            IKB_E = MIN(LKB_SM,IKB_B +  LEN_KB_BAT - 1)
            NKB_B = IKB_E - IKB_B + 1
            NKB_B0 = NKB_B
*. Loop over annihilation batches
            DO IAA_BT = 1, NAA_BT
             IAA_B = (IAA_BT-1)*MAXLB + 1
             IAA_E = MIN(LIAA_SM,IAA_B + MAXLB -1)
             NAA_B = IAA_E-IAA_B+1
*. KA => JA mapping
             CALL K_TO_J_TOT(
     &            KJAD,XKJAD,KASM,IKA_B,IKA_E,IAA_STR,IAA_SM,NAA_B,
     &            IAA_B,IAA_E,I3,XI3S,IB_AA,NJKAINTM,NAA_OP,IZEROKAJA)
            IF(IZEROKAJA.EQ.0) THEN
            DO IAB_BT = 1, NAB_BT
*B(eginning) and E(nd) of each Batch
             IAB_B = (IAB_BT-1)*MAXLB + 1
             IAB_E = MIN(LIAB_SM,IAB_B + MAXLB -1)
             NAB_B = IAB_E-IAB_B+1
*. KB => JB mapping
             CALL K_TO_J_TOT(
     &            KJBD,XKJBD,KBSM,IKB_B,IKB_E,IAB_STR,IAB_SM,NAB_B,
     &            IAB_B,IAB_E,I4,XI4S,IB_AB,NJKBINTM,NAB_OP,IZEROKBJB)
            IF(IZEROKBJB.EQ.0) THEN
*. Loop over creation batches 
            DO ICA_BT = 1, NCA_BT
             ICA_B = (ICA_BT-1)*MAXLB + 1
             ICA_E = MIN(LICA_SM,ICA_B + MAXLB -1)
             NCA_B = ICA_E-ICA_B+1
*. KA => IA mapping
             CALL K_TO_J_TOT(
     &            KIAD,XKIAD,KASM,IKA_B,IKA_E,ICA_STR,ICA_SM,NCA_B,
     &            ICA_B,ICA_E,I1,XI1S,IB_CA,NIKAINTM,NCA_OP,IZEROKAIA)
            DO ICB_BT = 1, NCB_BT
             ICB_B = (ICB_BT-1)*MAXLB + 1
             ICB_E = MIN(LICB_SM,ICB_B + MAXLB -1)
             NCB_B = ICB_E-ICB_B+1
*. KB => IB mapping
             CALL K_TO_J_TOT(
     &            KIBD,XKIBD,KBSM,IKB_B,IKB_E,ICB_STR,ICB_SM,NCB_B,
     &            ICB_B,ICB_E,I2,XI2S,IB_CB,NIKBINTM,NCB_OP,IZEROKBIB)
             IF(IZEROKBJB.EQ.0.AND.IZEROKAJA.EQ.0.AND.
     &          IZEROKBIB.EQ.0.AND.IZEROKAIA.EQ.0     ) THEN
*. Compress Ka and Kb strings to active subset 
C       COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
*. Ka
             CALL COMPRS2LST_B(KIAD,XKIAD,NCA_B,KJAD,XKJAD,NAA_B,
     &                       NKA_B0,NKA_B,KIA,XKIA,KJA,XKJA)
*. Kb
             CALL COMPRS2LST_B(KIBD,XKIBD,NCB_B,KJBD,XKJBD,NAB_B,
     &                       NKB_B0,NKB_B,KIB,XKIB,KJB,XKJB)
             IM_ACTIVE = 1
* Off set and form of T-coefficients
               IF(ITABTRNSP.EQ.0) THEN
                 ITOFF = IB_T(ICA_SM,ICB_SM,IAA_SM)
                 IF(ITOFF.LT.0) THEN
                   ITOFF = -ITOFF
                   ITRNSP = 1
                 ELSE
                   ITRNSP = 0
                 END IF
               ELSE
                 ITOFF = IB_T(ICB_SM,ICA_SM,IAB_SM)
                 IF(ITOFF.LT.0) THEN
                   ITOFF = -ITOFF
CE                 ITRNSP = 1
                   ITRNSP = 0
                 ELSE
CE                 ITRNSP = 0
                   ITRNSP = 1
                 END IF
               END IF
               IF(ITDIAG.EQ.1.AND.ICA_SM.EQ.ICB_SM.AND.
     &            IAA_SM.EQ.IAB_SM) THEN
                  IDIAG = 1
               ELSE
                  IDIAG = 0
               END IF
*
             IF(ISD.EQ.1) THEN
               CALL TCC_SUBBLK(TB(ITOFF),TBSUB,1,
     &         LICA_SM,NCA_B,ICA_B,LICB_SM,NCB_B,ICB_B,
     &         LIAA_SM,NAA_B,IAA_B,LIAB_SM,NAB_B,IAB_B,
     &         ITRNSP,IDIAG)
             ELSE
*. Zero density subblock
              ZERO = 0.0D0
              CALL SETVEC(TBSUB,ZERO,NCA_B*NCB_B*NAA_B*NAB_B)
             END IF
*  
* ==========================
*. Obtain C(Ka,Kb,Iaa,Iab)  
* ==========================
*
*. C(Ka,Kb,Iaa,Iab) = sum(Ja,Jb) <Ka!Iaa!Ja><Kb!Iab!Jb> C(Ja,Jb)
             CALL C_TO_CKK(CB(ICBOFF),
     &       CJRES,1,NJA(JASM),NJB(JBSM),NKA_B,NKB_B,NAA_B,NAB_B,
     &       KJA,XKJA,KJB,XKJB)
*
             IF(ISD.EQ.2) THEN
               CALL C_TO_CKK(SB(ISBOFF),
     &         SIRES,1,NIA(IASM),NIB(IBSM),NKA_B,NKB_B,NCA_B,NCB_B,
     &         KIA,XKIA,KIB,XKIB)
             END IF
* =================================================================
* For sigma :
* SIRES(Ka,Kb,Ica,Icb) = CJRES(Ka,Kb,Iaa,Iab)*T(Ica,Icb,Iaa,Iab)
* For densi :
* T(Ica,Icb,Iaa,Iab)   = SIRES(Ka,Kb,Ica,Icb)*CJRES(Ka,Kb,Ica,Icb)
* =================================================================
*
             LKAB = NKA_B*NKB_B
             LCAB = NCA_B*NCB_B 
             LAAB = NAA_B*NAB_B
*
             FACTORC = 0.0D0
             FACTORAB = 1.0D0*SIGNXXX
*
             IF(ISD.EQ.1) THEN
               CALL MATML7(SIRES,CJRES,TBSUB,
     &         LKAB,LCAB,LKAB,LAAB,LCAB,LAAB,FACTORC,FACTORAB,2)
*
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated S(Ka,Kb,*)  from MATML7 '
                 CALL WRTMAT(SIRES,LKAB,LCAB,LKAB,LCAB)
               END IF
             ELSE
               FACTORC = 0.0D0
               CALL MATML7(TBSUB,SIRES,CJRES,
     &         LCAB,LAAB,LKAB,LCAB,LKAB,LAAB,FACTORC,FACTORAB,1) 
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated TBSUB '
                 CALL WRTMAT(TBSUB,LCAB,LAAB,LCAB,LAAB)
               END IF
             END IF
* =================================================================
* For sigma :
* S(Ia,IB) = SI(Ia,Ib) + <Ka!O+_ca!Ia><Kb!O+cb!Ib>SIRES(Ka,Kb,Ica,Icb)
* For densi :
* Just scatter out to complete density block
* =================================================================
*
             IF(ISD.EQ.1) THEN
               CALL C_TO_CKK(SB(ISBOFF),
     &         SIRES,2,NIA(IASM),NIB(IBSM),NKA_B,NKB_B,NCA_B,NCB_B,
     &         KIA,XKIA,KIB,XKIB)
             ELSE IF(ITRNSP.EQ.0) THEN
C?             WRITE(6,*) ' Densi-call to TCC, ITOFF  = ', ITOFF
               CALL TCC_SUBBLK(TB(ITOFF),TBSUB,2,
     &         LICA_SM,NCA_B,ICA_B,LICB_SM,NCB_B,ICB_B,
     &         LIAA_SM,NAA_B,IAA_B,LIAB_SM,NAB_B,IAB_B,
     &         ITRNSP,IDIAG)
             END IF
            END IF
*           ^ End if all maps are nonvanishings
            END DO
            END DO
*           ^ End of loop over batches of creation   operators
            END IF
*           ^ End if batch of beta annihilation was nontrivial
            END DO
            END IF
*           ^ End if batch of alpha annihilation was nontrivial
            END DO
*           ^ End of loop over batches of annihilation operators
           END DO
           END DO
*          ^ End of loop over batches of Ka,Kb
           END IF
*          ^ End if symmetry combinations have nonvanishing dimensions
         END DO
         END DO
*        ^ End of loop over KASM, KBSM
        END IF
*       ^ End if symmetryblock of C should be included
        END DO
*       ^ End of loop over JASM
       END IF
*      ^ End if symmetryblock of S should be included
       END DO
*      ^ End of loop over IASM
      END IF
*     ^ End if Ka and Kb are nontrivial strings
      IF(NTEST.GE.100) THEN 
        IF(ISD.EQ.1) THEN
          WRITE(6,*) ' Updated sigma block '
C              WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
          CALL WRTVH1(SB,ISSM,NIA,NIB,NSMST,0)
        ELSE IF( ISD.EQ.2) THEN
          WRITE(6,*) ' Updated density block '
          CALL WRT_TCC_BLK(TB,1,NICA,NICB,NIAA,NIAB,NSMST)
C              WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
        END IF
      END IF
*
 9999 CONTINUE
*
      CALL LUCIAQEXIT('GENSIG')
      RETURN
      END
c $Id$
