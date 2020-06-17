      SUBROUTINE OFFSET_FOR_FTERM(IFTP,NSPOBEX,ISPOBEX,LSPOBEX,NGAS,
     &           IB_F,LEN_F,IONLY_LEN)
*
* Find offsets for the F-vector associated with 
* operator T(IFTP) 
* Only terms in T F that are connected 
* with the excitation manifold by single and  double 
* excitation are included
* LEN_F is total length of this F-vector
* IF IONLY_LEN = 1, only LEN_F is obtained 
*
* Jeppe Olsen, Nov. 2001, 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER ISPOBEX(4*NGAS,NSPOBEX), LSPOBEX(NSPOBEX)
*. Output, address of spinorbital type ITTP in F 
      INTEGER IB_F(NSPOBEX)
*. Local scratch
      INTEGER ITF_OCC(4*MXPNGAS)
*
      NTEST = 000
      IB = 1
      DO ITTP = 1, NSPOBEX
*. Occupation of ITTP*IFTP
        CALL OP_T_OCC(ISPOBEX(1,ITTP),ISPOBEX(1,IFTP),
     &                ITF_OCC,IMZERO)
        IADD = 0
        IF(IMZERO.EQ.0) THEN
*. Address of ITF, or number of excitations required to 
*. bring it into space
          CALL INUM_FOR_OCC2(ITF_OCC,INUM,NDIFF) 
*. Space will be allocated only if TF is not in space, 
*. and can be reached by atmost 4 operators
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IFTP, ITTP, NDIFF = ', IFTP, ITTP, NDIFF
          IF(0.LT.NDIFF.AND.NDIFF.LE.4) IADD = 1
        END IF
        IF(IADD.EQ.1) THEN
          IF(IONLY_LEN.NE.1) IB_F(ITTP) = IB
          IB = IB + LSPOBEX(ITTP)
        ELSE
          IF(IONLY_LEN.NE.1) IB_F(ITTP) = -1
        END IF
*
      END DO
*
      LEN_F = IB - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Type of T-operator : ', IFTP
        WRITE(6,*) ' Total length of F-terms = ', LEN_F
        IF(IONLY_LEN.NE.1) THEN
          WRITE(6,*) ' Array of offsets for F-term '
          WRITE(6,*) ' ============================ '
          CALL IWRTMA(IB_F,NSPOBEX,1,NSPOBEX,1)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE LEN_FOR_TF(NSPOBEX,ISPOBEX,LSPOBEX,NGAS,LEN_TF)
*
* Find largest Length of a T block and the associated F-vector
*
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER ISPOBEX(4*NGAS,NSPOBEX), LSPOBEX(NSPOBEX)
*. Local scratch
      INTEGER ITF_OCC(4*MXPNGAS)
*
      NTEST = 000
      LEN_TF = 0
      DO IFTP = 1, NSPOBEX+1
        LEN_F = 0
        DO ITTP = 1, NSPOBEX
*. Occupation of ITTP*IFTP
          CALL OP_T_OCC(ISPOBEX(1,ITTP),ISPOBEX(1,IFTP),
     &                  ITF_OCC,IMZERO)
          IADD = 0
          IF(IMZERO.EQ.0) THEN
*. Address of ITF, or number of excitations required to 
*. bring it into space
            CALL INUM_FOR_OCC2(ITF_OCC,INUM,NDIFF) 
*. Space will be allocated only if TF is not in space, 
*. and can be reached by atmost 4 operators
            IF(0.LT.NDIFF.AND.NDIFF.LE.4) IADD = 1
          END IF
          IF(IADD.EQ.1) LEN_F = LEN_F + LSPOBEX(ITTP)
        END DO
        LEN_TF = MAX(LEN_TF,LEN_F+LSPOBEX(IFTP))
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Largest length of an T-block and F-vector',
     &               LEN_TF
      END IF
*
      RETURN
      END
      SUBROUTINE COMP_CAAB_WLIST(ICAAB,ICAAB_LIST,NGAS,LCAAB_LIST,
     &                           INDEX)
*
* Compare CAAB operator ICAAB with list of CAAB operators (ICAAB_LIST)
* and find first occurance (INDEX) of ICAAB in ICAAB_LIST.
* If ICAAB is not found in ICAAB_LIST, INDEX is returned as 0
*
* Jeppe Olsen, Aug. 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICAAB(NGAS,4),ICAAB_LIST(NGAS,4,LCAAB_LIST)
*
      INDEX = 0
      DO JCAAB = 1, LCAAB_LIST
C            CMP_CAAB(ICAAB1,ICAAB2,IDENTICAL)
        CALL CMP_CAAB(ICAAB,ICAAB_LIST(1,1,JCAAB),IDENTICAL)
        IF(IDENTICAL.EQ.1) THEN
         INDEX = JCAAB 
         GOTO 1001
        END IF
      END DO
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CAAB operator : '
        CALL WRT_SPOX_TP(ICAAB,1)
        WRITE(6,*) ' Address in list = ', INDEX
      END IF
*
      RETURN
      END
*
      SUBROUTINE SPINFLIP_CAAB(ICAAB_IN,ICAAB_OUT,NGAS)
*
* Obtain ICAAB_OUT by spinflipping ICAAB_IN 
*
* Jeppe Olsen, August 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICAAB_IN(NGAS,4)
*. Output
      INTEGER ICAAB_OUT(NGAS,4)
*
      DO ICA = 1, 2
      DO IAB = 1, 2
       IBA = 2/IAB
       ICAAB = (ICA-1)*2 + IAB
       ICABA = (ICA-1)*2 + IBA
       DO IGAS = 1, NGAS
         ICAAB_OUT(IGAS,ICABA) = ICAAB_IN(IGAS,ICAAB)
       END DO
      END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CAAB and spin-flipped CAAB '
       CALL WRT_SPOX_TP(ICAAB_IN,1)
       CALL WRT_SPOX_TP(ICAAB_OUT,1)
      END IF
*
      RETURN
      END
      SUBROUTINE CMP_CAAB(ICAAB1,ICAAB2,IDENTICAL)
*
* Check to see if two excitation types are identical 
*
* Jeppe Olsen, August 2001
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*
      INTEGER ICAAB1(NGAS,4),ICAAB2(NGAS,4)
*
      NDIFF = 0
      DO JCAAB = 1, 4
      DO JGAS = 1, NGAS
        NDIFF = NDIFF + ABS(ICAAB1(JGAS,JCAAB)-ICAAB2(JGAS,JCAAB))
      END DO
      END DO
*
      IF(NDIFF.EQ.0) THEN
       IDENTICAL = 1
      ELSE
       IDENTICAL = 0
      END IF
*
      RETURN
      END

      SUBROUTINE H_TYPE_SYM12(IHCAAB,NGAS,ISYM12)
*
* Check if given type of Hamiltonian has symmetry between 
* particle one and two
*
* Jeppe Olsen, July 18, 2001
*
      INCLUDE 'implicit.inc'
      INTEGER IHCAAB(NGAS,4)
*. Symmetry between particle one and two requires 
*. that both creation operators belongs to the same type and
*. spin , and that annihilation operators belongs to the same 
*. type and spin
*
      ISYM12 = 1
      DO ICAAB= 1, 4
       DO IGAS = 1, NGAS
        IF( IHCAAB(IGAS,ICAAB).EQ.1) ISYM12 = 0
       END DO
      END DO
*
      RETURN
      END 
      SUBROUTINE SPINFLIP_CC_BLOCKS(CC,NSPOBEX_TP,ISPOBEX_PAIRS,
     &           IBSPOBEX,ISPOBEX,ISM,NSMST,NGAS)
*
* A CC-vector CC is given where the active CC blocks are known.
* Find the full CC-vector, by obtaining the passive blocks throgh 
* spin-flip
*
* Jeppe Olsen, July12, 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISPOBEX_PAIRS(*),IBSPOBEX(*)
      INTEGER ISPOBEX(4*NGAS,*)
*. Input and output
      DIMENSION CC(*)
*. Local scratch
      DIMENSION LCAAB(8,4)
*
      NTEST = 00
*
      DO ITP = 1, NSPOBEX_TP
        IF(ISPOBEX_PAIRS(ITP).LT.0.OR.
     &     ISPOBEX_PAIRS(ITP).EQ.ITP) THEN
*. Block not explicitly calculated 
          ITP_IN = ABS(ISPOBEX_PAIRS(ITP))
          ITP_OUT = ITP
*. Number of strings per symmetry for input block 
          CALL NST_CAAB(ISPOBEX(1,ITP_IN),LCAAB)
*. Do the spinflip for this block
          IB_IN = IBSPOBEX(ITP_IN)
          IB_OUT = IBSPOBEX(ITP_OUT)
          IF(ISPOBEX_PAIRS(ITP).LT.0) THEN
            IDIAG = 0
          ELSE 
            IDIAG = 1
          END IF
          CALL SPINFLIP_CC_BLOCK(CC(IB_IN),CC(IB_OUT),
     &    LCAAB(1,1),LCAAB(1,2),LCAAB(1,3),LCAAB(1,4),NSMST,ISM,IDIAG)
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
C  CALL WRT_CC_VEC2(CC_VEC,6,CCTYPE)
        WRITE(6,*) ' Complete CC_vector with passive blocks'
        CALL WRT_CC_VEC2(CC,6,'GEN_CC')
      END IF
*
      RETURN
      END
      SUBROUTINE SPINFLIP_CC_BLOCK(CC_IN,CC_OUT,NCA,NCB,NAA,NAB,NSMST,
     &                             ISM,IDIAG)
*
* A CC block CC_IN is given.
* Obtain CC block CC_OUT by spinflipping
*
* If IDIAG = 1, the block is diagonal and is symmetrized 
*
* Jeppe Olsen, July 12 2001
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      DIMENSION CC_IN(*)
*. Dimensions for strings in CC_IN
      INTEGER NCA(NSMST),NCB(NSMST),NAA(NSMST),NAB(NSMST)
*. Output
      DIMENSION CC_OUT(*)
*. Local scratch
      INTEGER IB(8,8,8)
*. Set up Offset array for CC_OUT
C Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST,IDIAG)
      CALL Z_TCC_OFF(IB,NCB,NCA,NAB,NAA,ISM,NSMST,0)
*. Loop over symmetry blocks in INPUT block
      I_IN = 0
      DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        DO ISM_CA = 1, NSMST
          ISM_CB = MULTD2H(ISM_C,ISM_CA)
          DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
            IB_OUT = IB(ISM_CB,ISM_CA,ISM_AB)
*. Dimensions in input block
            LCA = NCA(ISM_CA)
            LCB = NCB(ISM_CB)
            LAA = NAA(ISM_AA)
            LAB = NAB(ISM_AB)
*. Loop over elements in input block 
            DO IAB = 1, LAB
            DO IAA = 1, LAA
            DO ICB = 1, LCB
            DO ICA = 1, LCA
              I_IN = I_IN + 1
              I_OUT = IB_OUT -1 + (IAA-1)*LAB*LCA*LCB
     &              + (IAB-1)*LCA*LCB + (ICA-1)*LCB + ICB
*
              IF(IDIAG.EQ.0) THEN
                CC_OUT(I_OUT) = CC_IN(I_IN) 
              ELSE IF(I_IN .GT. I_OUT) THEN
                CC_OUT(I_IN) = 0.5D0*( CC_IN(I_OUT) +  CC_IN(I_IN) )
                CC_OUT(I_OUT) = CC_OUT(I_IN)
              END IF
*
            END DO
            END DO
            END DO
            END DO
*           ^ End of loops over Elements in symmetry subblock
          END DO
        END DO
      END DO
*     ^ End of loops over symmetry subblocks 
      NELMNT = I_IN
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input CC block '
        CALL WRTMAT(CC_IN,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' Output CC block '
        CALL WRTMAT(CC_OUT,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
      SUBROUTINE TEST_D1
*
* Looking  for error in SMD1 NMD1
*
* Jeppe, June 24
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
*
      COMMON/CC_SCR/ KCCF,KCCE,KCCVEC1,KCCVEC2,
     &              KISTST1,KXSTST1,KISTST2,KXSTST2,
     &              KISTST3,KXSTST3,KISTST4,KXSTST4,
     &              KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &              KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTREO,
     &              KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &              KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &              KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &              KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &              KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &              KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &              KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &              KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &              KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
     &              KLOPSCR,
     &              KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &              KLSMD1,KLSMD2,KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &              KLNMD1,KLNMD2,KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &              KLOCK1, KLOCK2, KLOCL1, KLOCOT1T2, KL_IBF, 
     &              KLEXEORD
*
      LEN = 4*MXTSOB_P**2 * MXTSOB_H**2
      CALL TEST_D1S(WORK(KLSMD1),WORK(KLNMD1),LEN)
*
      RETURN
      END
      SUBROUTINE TEST_D1S(ISMD1,INMD1,LEN)
*
      INCLUDE 'implicit.inc'
      INTEGER ISMD1(LEN),INMD1(LEN)
*
      IOFF = 4*20+1
      NERR = 0
      DO I = IOFF, LEN
        IF(ISMD1(I).NE.-42) THEN
          WRITE(6,*) ' ISMD1 for I = ', I, ' is ', ISMD1(I)
          NERR = NERR + 1
        ELSE IF (INMD1(I).NE.-42) THEN
          WRITE(6,*) ' INMD1 for I = ', I, ' is ', INMD1(I) 
          NERR = NERR + 1
        END IF
      END DO
*
      IF(NERR.NE.0) THEN
        WRITE(6,*) ' Number of overwrites = ', NERR
        WRITE(6,*)' Passive array of ISMD1/INMD1 overwritten '
        STOP ' Passive array of ISMD1/INMD1 overwritten '
      END IF 
*
      RETURN
      END
      
            
      

      SUBROUTINE TI_TO_TOKBN(NSMST,
     &           IOP_SM_C1,IOP_SM_CA1,IOP_SM_AA1,IOP_SM,
     &           IOP_NM_AB1,IOP_NM_AA1,IOP_NM_CB1,IOP_NM_CA1,
     &           NOP_BAT,NOP_CAAB,
     &           KK_SM_C1,KK_SM_CA1,KK_SM_AA1,KK_SM,
     &           KK_NM_AB1,KK_NM_AA1,KK_NM_CB1,KK_NM_CA1,
     &           NK_BAT,NKK_CAAB,
     &           IKJ_CA_MAP,SKJ_CA_MAP,IBKJ_CA,
     &           IKJ_CB_MAP,SKJ_CB_MAP,IBKJ_CB,
     &           IKJ_AA_MAP,SKJ_AA_MAP,IBKJ_AA,
     &           IKJ_AB_MAP,SKJ_AB_MAP,IBKJ_AB,
     &           NISFSM_CA,NISFSM_CB,NISFSM_AA,NISFSM_AB,
     &           TOPK,TI,IB_TI,ISG,LDUM,SIGNI,IOPDAG)
*
* T(Idum, I) <=> T(Idum,OP,K) formats
*
* Input :
* =======
*
* NSMST : Number of symmetries of strings
* IOP_SM_C1,IOP_SM_CA1,IOP_SM_AA1 : start values of symmetries of OP
* for C, CA and AA part of operator, used for batching of operators. 
* IOP_SM : Overall symmetry of IOP 
* IOP_NM_AB1,IOP_NM_AA1,IOP_NM_CB1,IOP_NM_CA1, : Initial numbers of CAAB
* parts of IOP 
* NOP_BAT : Number of operators OP in current batch
* NOP_CAAB : Number of operators per sym for CAAB parts of OP
* KK_SM_C1,KK_SM_CA1,KK_SM_AA1 : start values of symmetries of KK
* for C, CA, AA part of operator, used for batching of K
* KK_SM : Symmetry of K
* KK_NM_AB1,KK_NM_AA1,KK_NM_CB1,KK_NM_CA1 : Initial numbers of CAAB 
* parts of K
* NK_BAT : Number of K operators in current batch
* NKK_CAAB : Number of operators per sym for CAAB parts of K
* IKJ_CA_MAP :  CA part of OP, K => I map 
* SKJ_CA_MAP :  Sign part of above map
* IBKJ_CA : off sets to IKJ_CA_MAP for K and OP with given sym
* IKJ_CB_MAP, IKJ_AA_MAP, IKJ_AB_MAP : CB, AA, AB part of map 
* SKJ_CB_MAP, SKJ_AA_MAP, SKJ_AB_MAP : sign part of above maps
* IBKJ_CB, IBKJ_AA, IBKJ_AB : off sets to mappings for K and Op with given sym
* NISFSM_CA : Number of CA- strings in I with given sym 
* NISFSM_CB : Number of CB- strings in I with given sym 
* NISFSM_AA : Number of AA- strings in I with given sym 
* NISFSM_AB : Number of AB- strings in I with given sym 
* TOPK : The T(Idum,Op,K) array
* TI   : The T(Idum,I) array
* IB_TI : Off sets to T with given sym of CA, CB, AA
* ISG : =1 => Scatter, =2 => Gather  
* LDUM : Dimension of dummy index 
* SIGNI : Initial sign 
* IOPDAG  = 1 :  IOP_NM, IOP_SM contains info for daggered operator.

* 
* A set of operators OP_CA, OP_CB, OP_AA, OP_AB is given 
* OP_NM_CAAB, OP_SM_CAAB contains a number of such operator quadruplets
*   
*
* ISG = 1 : 
*   T(IDUM,OP,K) :=  Sign(OP,K)* T(IDUM,I)
* ISG = 2 :
*   T(IDUM,I) := T(IDUM,OP,K) + Sign(OP,K)*T(IDUM,I)
*
* If IOPDAG = 1 IOP_NM, IOP_SM contains info for daggered operator 
* instead of operator itself
*
* Version with complete precomputed maps OP * K => I
*
* Jeppe Olsen, October 2000 (TOKB version)
*              June 2001    TOKBN (Turbo version)
*              
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
*. Operator
      INTEGER NOP_CAAB(8,4)
*. Number of K-strings per symmetry
      INTEGER NKK_CAAB(8,4)
*. Offset to strings with given symmetry
C     INTEGER IBOP_CA(8),IBOP_CB(8),IBOP_AA(8)
*. Offset to block with given sym in TI
      INTEGER IB_TI(8,8,8)
*. Number of strings in I per symmetry
      INTEGER NISFSM_CA(8),NISFSM_CB(8),NISFSM_AA(8),NISFSM_AB(8)
*. Kstrings
COLD  INTEGER K_SM_CAAB(4,NK_BAT),K_NM_CAAB(4,NK_BAT)
*. Info on OP * K to I mappings
      INTEGER IKJ_CA_MAP(*),IBKJ_CA(8,8)
      INTEGER IKJ_CB_MAP(*),IBKJ_CB(8,8)
      INTEGER IKJ_AA_MAP(*),IBKJ_AA(8,8)
      INTEGER IKJ_AB_MAP(*),IBKJ_AB(8,8)
*
      DIMENSION SKJ_CA_MAP(*),SKJ_CB_MAP(*),SKJ_AA_MAP(*),SKJ_AB_MAP(*)
*. Input or output
C     DIMENSION TOPK(LDUM,NOP_BAT,NK_BAT)
C     DIMENSION TI(LDUM,*)
      DIMENSION TOPK(*),TI(*)
      COMMON/ROLLO/NINNER,NINNER1,NINNER2,NINNER3,NINNER4,NINNER5,
     &             NINNER6,NINNER7,NINNER8,NINNER9,NINNER10
*

      NTEST = 0
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' TI_TO_TOKBN : SIGNI, LDUM, NOP_BAT, NK_BAT  = ',
     &                           SIGNI, LDUM, NOP_BAT, NK_BAT
      MX_TI_COL = 0
      I_CAAB_TI_MX = 0
*.
      IF(ISG.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(TOPK,ZERO,LDUM*NOP_BAT*NK_BAT)
      END IF
*. Loop over I-operators
      INI_LOOP_OP = 1
      IOP = 0
      DO IOP_SM_C = IOP_SM_C1, NSMST
        IOP_SM_A = MULTD2H(IOP_SM,IOP_SM_C) 
        IF(INI_LOOP_OP .EQ. 1 ) THEN
          IOP_SM_CA2 = IOP_SM_CA1
        ELSE
          IOP_SM_CA2 = 1
        END IF
        DO IOP_SM_CA = IOP_SM_CA2, NSMST
          IOP_SM_CB = MULTD2H(IOP_SM_C,IOP_SM_CA)
          IF(INI_LOOP_OP .EQ. 1 ) THEN
            IOP_SM_AA2 = IOP_SM_AA1
          ELSE
            IOP_SM_AA2 = 1
          END IF
          DO IOP_SM_AA = IOP_SM_AA2, NSMST
            IOP_SM_AB =  MULTD2H(IOP_SM_A,IOP_SM_AA)
*. Loop over operators as  matrix (I_CA,I_CB,I_AA,I_AB)
            NOP_STR_CA = NOP_CAAB(IOP_SM_CA,1)
            NOP_STR_CB = NOP_CAAB(IOP_SM_CB,2)
            NOP_STR_AA = NOP_CAAB(IOP_SM_AA,3)
            NOP_STR_AB = NOP_CAAB(IOP_SM_AB,4)
C?          WRITE(6,*) ' IOP_SM_CA, IOP_SM_CB,IOP_SM_AA, IOP_SM_AB ',
C?   &                   IOP_SM_CA, IOP_SM_CB,IOP_SM_AA, IOP_SM_AB
C?          WRITE(6,*) 
C?   &      ' NOP_STR_CA, NOP_STR_CB, NOP_STR_AA,NOP_STR_AB' ,
C?   &        NOP_STR_CA, NOP_STR_CB, NOP_STR_AA,NOP_STR_AB
            IF(NOP_STR_CA*NOP_STR_CB*NOP_STR_AA*NOP_STR_AB.NE.0)THEN
            NINNER10 = NINNER10 + 1
*
            IF(INI_LOOP_OP.EQ.1) THEN
             IOP_NM_AB2 = IOP_NM_AB1
            ELSE
             IOP_NM_AB2 = 1
            END IF
            DO IOP_NM_AB = IOP_NM_AB2, NOP_STR_AB
C?       WRITE(7,*) ' IOP_NM_AB = ', IOP_NM_AB
                NINNER9 = NINNER9 + 1
             IF(INI_LOOP_OP.EQ.1) THEN
              IOP_NM_AA2 = IOP_NM_AA1
             ELSE
              IOP_NM_AA2 = 1
             END IF
             DO IOP_NM_AA = IOP_NM_AA2, NOP_STR_AA
                NINNER8 = NINNER8 + 1
              IF(INI_LOOP_OP.EQ.1) THEN
                IOP_NM_CB2 = IOP_NM_CB1
              ELSE 
                IOP_NM_CB2 = 1
              END IF
C?       WRITE(7,*) ' IOP_NM_CB2 = ', IOP_NM_CB2
              DO IOP_NM_CB = IOP_NM_CB2, NOP_STR_CB
C?       WRITE(7,*) ' IOP_NM_CB = ', IOP_NM_CB
                NINNER7 = NINNER7 + 1
               IF(INI_LOOP_OP.EQ.1) THEN
                 IOP_NM_CA2 = IOP_NM_CA1
               ELSE 
                 IOP_NM_CA2 = 0
               END IF
               DO IOP_NM_CA = IOP_NM_CA2 + 1, NOP_STR_CA
                NINNER6 = NINNER6 + 1
                IOP = IOP + 1
*. A major part of the following could be moved outside
                IF(IOPDAG.EQ.0) THEN
                  IIOP_SM_CA = IOP_SM_CA
                  IIOP_NM_CA = IOP_NM_CA
*
                  IIOP_SM_CB = IOP_SM_CB
                  IIOP_NM_CB = IOP_NM_CB 
*
                  IIOP_SM_AA = IOP_SM_AA
                  IIOP_NM_AA = IOP_NM_AA 
*
                  IIOP_SM_AB = IOP_SM_AB
                  IIOP_NM_AB = IOP_NM_AB
                ELSE
                  IIOP_SM_CA = IOP_SM_AA
                  IIOP_NM_CA = IOP_NM_AA
*
                  IIOP_SM_CB = IOP_SM_AB
                  IIOP_NM_CB = IOP_NM_AB 
*
                  IIOP_SM_AA = IOP_SM_CA
                  IIOP_NM_AA = IOP_NM_CA 
*
                  IIOP_SM_AB = IOP_SM_CB
                  IIOP_NM_AB = IOP_NM_CB
                END IF
*
                IF(IOPDAG.EQ.0) THEN
                  LOP_CA = NOP_CAAB(IIOP_SM_CA,1)
                  LOP_CB = NOP_CAAB(IIOP_SM_CB,2)
                  LOP_AA = NOP_CAAB(IIOP_SM_AA,3)
                  LOP_AB = NOP_CAAB(IIOP_SM_AB,4)
                ELSE 
                  LOP_CA = NOP_CAAB(IIOP_SM_CA,3)
                  LOP_CB = NOP_CAAB(IIOP_SM_CB,4)
                  LOP_AA = NOP_CAAB(IIOP_SM_AA,1)
                  LOP_AB = NOP_CAAB(IIOP_SM_AB,2)
                END IF
*. Loop over K-strings 
      INI_LOOP_KK = 1
      KSTR = 0
      DO KK_SM_C = KK_SM_C1, NSMST
        KK_SM_A = MULTD2H(KK_SM,KK_SM_C) 
        IF(INI_LOOP_KK .EQ. 1 ) THEN
          KK_SM_CA2 = KK_SM_CA1
        ELSE
          KK_SM_CA2 = 1
        END IF
        DO KK_SM_CA = KK_SM_CA2, NSMST
          KK_SM_CB = MULTD2H(KK_SM_C,KK_SM_CA)
*
          NKK_STR_CA = NKK_CAAB(KK_SM_CA,1)
          NKK_STR_CB = NKK_CAAB(KK_SM_CB,2)
C?        WRITE(6,*) ' KK_SM_CA, KK_SM_CB, NKK_STR_CA, NKK_STR_CB ',
C?   &                 KK_SM_CA, KK_SM_CB, NKK_STR_CA, NKK_STR_CB
*
          IBOPK_CA = IBKJ_CA(IIOP_SM_CA, KK_SM_CA)
          I_CA_SM = MULTD2H(IIOP_SM_CA,KK_SM_CA)
*
          IBOPK_CB = IBKJ_CB(IIOP_SM_CB, KK_SM_CB)
          I_CB_SM = MULTD2H(IIOP_SM_CB,KK_SM_CB)
*
          IF(INI_LOOP_KK .EQ. 1 ) THEN
            KK_SM_AA2 = KK_SM_AA1
          ELSE
            KK_SM_AA2 = 1
          END IF
          IF( NKK_STR_CA* NKK_STR_CB.EQ.0) KK_SM_AA2 = NSMST + 1
          DO KK_SM_AA = KK_SM_AA2, NSMST
            KK_SM_AB =  MULTD2H(KK_SM_A,KK_SM_AA)
*
            IBOPK_AA = IBKJ_AA(IIOP_SM_AA, KK_SM_AA)
            I_AA_SM = MULTD2H(IIOP_SM_AA,KK_SM_AA)
*
            IBOPK_AB = IBKJ_AB(IIOP_SM_AB, KK_SM_AB)
C?          WRITE(6,*) ' KK_SM_AA, KK_SM_AB = ',
C?   &                   KK_SM_AA, KK_SM_AB
C?          WRITE(6,*) ' IIOP_SM_AB, KK_SM_AB, IBOPK_AB =',
C?   &                   IIOP_SM_AB, KK_SM_AB, IBOPK_AB

*. Offset to block in TI with this symmetry combination
            IBTI_SSSS = IB_TI(I_CA_SM,I_CB_SM,I_AA_SM)
C?          WRITE(6,*) ' I_CA_SM, I_CB_SM, I_AA_SM, IBTI_SSSS = ',
C?   &                   I_CA_SM, I_CB_SM, I_AA_SM, IBTI_SSSS

*. Loop over operators as  matrix (CA,CB,AA,AB)
            NKK_STR_AA = NKK_CAAB(KK_SM_AA,3)
            NKK_STR_AB = NKK_CAAB(KK_SM_AB,4)
C?          WRITE(6,*) ' NKK_STR_AA, NKK_STR_AB ',
C?   &                   NKK_STR_AA, NKK_STR_AB
C           IF(NKK_STR_CA*NKK_STR_CB*NKK_STR_AA*NKK_STR_AB.NE.0)THEN
           
*
            IF(INI_LOOP_KK.EQ.1) THEN
             KK_NM_AB2 = KK_NM_AB1
            ELSE
             KK_NM_AB2 = 1
            END IF
*. Tired of IF/ END IF today
            IF(NKK_STR_AA*NKK_STR_AB.EQ.0) THEN 
              KK_NM_AB2 = NKK_STR_AB + 1
            ELSE 
              NINNER5 = NINNER5 + 1
            END IF
            DO KK_NM_AB = KK_NM_AB2, NKK_STR_AB
          NINNER4 = NINNER4 + 1
C?        WRITE(7,*) ' IBOPK_AB, KK_NM_AB, LOP_AB, IIOP_NM_AB ',
C?   &                 IBOPK_AB, KK_NM_AB, LOP_AB, IIOP_NM_AB
          I_AB = IKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
          S_AB = SKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
             IF(INI_LOOP_KK.EQ.1) THEN
              KK_NM_AA2 = KK_NM_AA1
             ELSE
              KK_NM_AA2 = 1
             END IF
             DO KK_NM_AA = KK_NM_AA2, NKK_STR_AA
          NINNER3 = NINNER3 + 1
          I_AA = IKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
          S_AA = SKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
              IF(INI_LOOP_KK.EQ.1) THEN
                KK_NM_CB2 = KK_NM_CB1
              ELSE 
                KK_NM_CB2 = 1
              END IF
*
              I_CAAB_TI00 = IBTI_SSSS-1+ 
     &        (I_AB-1)*
     &        NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)*NISFSM_AA(I_AA_SM)
     &       +(I_AA-1)*NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)

*
C?            WRITE(6,*) 'IBTI_SSSS, I_CA_SM, I_CB_SM, I_AA_SM ',
C?   &                    IBTI_SSSS, I_CA_SM, I_CB_SM, I_AA_SM 
C?            WRITE(6,*) ' I_AA, I_AB = ', I_AA, I_AB 
C?            WRITE(6,*) ' I_CAAB_TI00  ', I_CAAB_TI00
*
              IKJ_CB_ADR = IBOPK_CB-1+(KK_NM_CB2-2)*LOP_CB+IIOP_NM_CB
              DO KK_NM_CB = KK_NM_CB2, NKK_STR_CB
C?             WRITE(6,*) ' KK_NM_CB = ', KK_NM_CB
               NINNER2 = NINNER2 + 1
               IKJ_CB_ADR =  IKJ_CB_ADR + LOP_CB
*
               I_CB = IKJ_CB_MAP(IKJ_CB_ADR)
               S_CB = SKJ_CB_MAP(IKJ_CB_ADR)
               SIGN123 = S_CB*S_AA*S_AB*SIGNI
*
               IF(INI_LOOP_KK.EQ.1) THEN
                 KK_NM_CA2 = KK_NM_CA1
               ELSE 
                 KK_NM_CA2 = 0
               END IF
*. 
               IF(I_AA*I_AB*I_CB.EQ.0) THEN
                 KSTR = KSTR + NKK_STR_CA - KK_NM_CA2
                 IF(KSTR.GE.NK_BAT) GOTO 2001
               ELSE
                 I_CAAB_TI0 =  I_CAAB_TI00  
     &          +(I_CB-1)*NISFSM_CA(I_CA_SM) - 1
*
C?               WRITE(6,*) 
C?   &           ' I_CAAB_TI0, I_CAAB_TI00, I_CB, NISFSM_CA(I_CA_SM)',
C?   &             I_CAAB_TI0, I_CAAB_TI00, I_CB, NISFSM_CA(I_CA_SM)
*
                 KJ_ADR = IBOPK_CA-1 + (KK_NM_CA2-1)*LOP_CA+IIOP_NM_CA
                 IF(LDUM.EQ.1) THEN
*. Separate loop for LDUM = 1
                   IADR_TOPK0 = (KSTR-1)*NOP_BAT + IOP-1
                   N_CA_EFF = NKK_STR_CA - KK_NM_CA2
                   IF(KSTR+N_CA_EFF.GT.NK_BAT) N_CA_EFF = NK_BAT - KSTR
                   KSTR = KSTR + N_CA_EFF
*
                   IF(ISG.EQ.1) THEN
                    DO KK_NM_CA = 1,  N_CA_EFF
                     IADR_TOPK0 = IADR_TOPK0 + NOP_BAT
                     KJ_ADR = KJ_ADR + LOP_CA
                     I_CA = IKJ_CA_MAP(KJ_ADR)
                     IF(I_CA .NE. 0) THEN
                       SIGN = SKJ_CA_MAP(KJ_ADR)*SIGN123
                       TOPK(1+IADR_TOPK0) = SIGN*TI(1+I_CAAB_TI0+I_CA)
                     END IF
*                    ^ End if Istrings was nonvanishing
                    END DO
*                   ^ End of loop over KK_NM_CA
                    IF(KSTR.EQ.NK_BAT) GOTO 2001
                   ELSE  
*                  ^ ISG switch 
                    DO KK_NM_CA = 1,  N_CA_EFF
C?                   WRITE(6,*) ' KK_NM_CA (b) = ', KK_NM_CA
                     IADR_TOPK0 = IADR_TOPK0 + NOP_BAT
                     KJ_ADR = KJ_ADR + LOP_CA
                     I_CA = IKJ_CA_MAP(KJ_ADR)
*
                     IF(I_CA .NE. 0) THEN
                       SIGN = SKJ_CA_MAP(KJ_ADR)*SIGN123
                       TI(1 + I_CAAB_TI0 + I_CA) = 
     &                 TI(1 + I_CAAB_TI0 + I_CA)
     &                +SIGN*TOPK(1+IADR_TOPK0)
                       I_CAAB_TI_MX = 
     &                 MAX(I_CAAB_TI_MX,1+I_CAAB_TI0+I_CA)
C?                     WRITE(6,*) ' I_CA, I_CAAB_TI0', 
C?   &                              I_CA, I_CAAB_TI0
                     END IF
*                    ^ End if Istrings was nonvanishing
                    END DO
*                   ^ End of loop over KK_NM_CA
                    IF(KSTR.EQ.NK_BAT) GOTO 2001
                  END IF 
*                 ^ ISG switch 
                 ELSE 
*                ^ If LDUM .eq .1
                   IADR_TOPK0 = (KSTR-1)*NOP_BAT*LDUM + (IOP-1)*LDUM
                   DO KK_NM_CA = KK_NM_CA2 + 1, NKK_STR_CA
C?                  WRITE(6,*) ' KK_NM_CA = ', KK_NM_CA
                    NINNER1 = NINNER1 + 1
                    KSTR = KSTR + 1
                    IADR_TOPK0 = IADR_TOPK0 + NOP_BAT*LDUM
                    KJ_ADR = KJ_ADR + LOP_CA
                    I_CA = IKJ_CA_MAP(KJ_ADR)
C?                  WRITE(6,*) ' I_CA = ', I_CA
*
                    IF(I_CA .NE. 0) THEN
                      NINNER = NINNER + LDUM
*. Adress of I_CB, I_CA, I_AA, I_AB
                      SIGN = SKJ_CA_MAP(KJ_ADR)*SIGN123
                      IF(ISG.EQ.1) THEN
                       DO I = 1, LDUM
                         IADR_TI0 = (I_CAAB_TI0 + I_CA)*LDUM 
                         TOPK(I+IADR_TOPK0) = SIGN*TI(I+IADR_TI0)
                       END DO
                      ELSE IF (ISG.EQ.2) THEN
                        DO I = 1, LDUM
                         I_CAAB_TI_MX 
     &                 = MAX(I_CAAB_TI_MX,I_CAAB_TI0 + I_CA+1)
                         IADR_TI0 = (I_CAAB_TI0 + I_CA)*LDUM 
                         TI(I+IADR_TI0) = 
     &                   TI(I+IADR_TI0)+SIGN*TOPK(I+IADR_TOPK0)
                        END DO
                      END IF
*                     ^ End of scatter/gather switch 
                    END IF
*                   ^ End if Istring was nonvanishing
                    IF(KSTR.EQ.NK_BAT) GOTO 2001
                   END DO
*                  ^ End of loop over KK_NM_CA
                 END IF
*                ^ End if LDUM = 1 switch

               END IF
*              ^ End if I_AA*I_AB*I_CB .NE. 0
               INI_LOOP_KK = 0
              END DO
*              ^ End of loop over KK_NM_CB
              INI_LOOP_KK = 0
             END DO
             INI_LOOP_KK = 0
            END DO
*           ^ End of loop over elements of block
C           END IF
*           ^ End if number of K-strings was nonvanishing
            INI_LOOP_KK = 0
          END DO
          INI_LOOP_KK = 0
*         ^ End of loop over ISM_AA
        END DO
        INI_LOOP_KK = 0
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
 2001 CONTINUE
* ^ End of loop over K-strings
*
                IF(IOP.EQ.NOP_BAT) GOTO 1001
               END DO
               INI_LOOP_OP = 0
              END DO
              INI_LOOP_OP = 0
             END DO
             INI_LOOP_OP = 0
            END DO
            INI_LOOP_OP = 0
*           ^ End of loop over elements of block
            END IF
*           ^ End if number of I strings was nonvanishing
          INI_LOOP_OP = 0
          END DO
*         ^ End of loop over ISM_AA
        END DO
        INI_LOOP_OP = 0
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
 1001 CONTINUE
*
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' TI_TO_TOKB speaking : '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' NOP_BAT, NK_BAT, LDUM = ', NOP_BAT, NK_BAT, LDUM
        WRITE(6,*) ' Largest column number used = ',  I_CAAB_TI_MX
        WRITE(6,*) ' ISG = ', ISG
        WRITE(6,*) ' TI as TI(Idum,I)  '
        CALL WRTMAT(TI,LDUM,I_CAAB_TI_MX,LDUM,I_CAAB_TI_MX)
        WRITE(6,*) ' TOPK as TOPK(IdumIop,Kstr) '
        CALL WRTMAT(TOPK,LDUM*NOP_BAT,NK_BAT,LDUM*NOP_BAT,NK_BAT)
      END IF

*
      RETURN
      END
      SUBROUTINE HEXP_T(T,E,F,NOBEX_TP,ISOBEX_TP,NSOBEX_TP,
     &           LSOBEX_TP,IBSOBEX_TP,ISOX_TO_OX,
     &           CCVEC1,CCVEC2,
     &           ISTST1,XSTST1,ISTST2,XSTST2,
     &           ISTST3,XSTST3,ISTST4,XSTST4,
     &           IZ,IZSCR,ISTOCC1,ISTOCC2,ISTREO,
     &           I_DO_ZERO,I_DO_FTERMS,ISPOBEX_AC,
     &           IEXEORD,HEXPT,ISPOBEX_FRZ,ISOBEX_PAIRS,IB_F,
     &           LUT)
*
* routine for calculation of the part of H exp(T) relevant for the 
* evaluation of the CC vector function
*
* Exp(T) is produced as
*
* E0 +  sum(I) E(I) + sum(I) T(I) F(I), 
* where sum is over spin orbital types , and E0 is the initial value 
* in E(N_CC_AMP+1). 
*
* Version where F-terms are contracted with H    
*
* T is on discfile LUT
*
*. About the algorithm
* =====================
*
* Contributions from each spin-orbital excitation operator T(I) are
* evaluated separately. 
*
*. After K terms we have the form 
*  E0 + sum_I E(K,I) + sum_I T(I) F(K,I)
* and assume that E(K,I) is in E
*. and we must now calculate in CCVEC3( setting T' = T(K+1))
* (1 + T' + 1/2 T'^2 + ... 1/l! T'^l) (1 + sum_I E(K,I) + sum_I T(I) F(K,I))
*
*  N = 0 term : Current E-term
*  (N = 1)  
*           Copy E to CCVEC1
*           T' + T' CCVEC1 in CCVEC2
*           add  CCVEC2 to E
*           copy CCVEC2 to CCVEC1
*  (N > 1)  T' CCVEC1 in CCVEC2
*           copy CCVEC2 to CCVEC1
*           Add  1/n! CCVEC2 to CCVEC3
*           Goto to next N if norm of CCVEC2 is gt. 0
*
*   The order in which T(N)T(N-1) ... T(2)T(1) is defined by 
*   IEXEORD
*           
*
* Jeppe Olsen, May 2000
*              Jan 2001 : IEXEORD added
*              June 2001 : Version where integrals and F terms are 
*                          contracted directly after construction 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
      REAL*8 INPROD
*
*. Input
*
      DIMENSION T(N_CC_AMP)
      INTEGER ISOBEX_TP(4*NGAS,NSOBEX_TP)
      INTEGER LSOBEX_TP(NSOBEX_TP),IBSOBEX_TP(NSOBEX_TP)
      INTEGER ISOX_TO_OX(NSOBEX_TP)
      INTEGER ISPOBEX_AC(NSOBEX_TP), ISPOBEX_FRZ(NSOBEX_TP)
      INTEGER IB_F(NSOBEX_TP)
      INTEGER IEXEORD(NSOBEX_TP)
      INTEGER ISOBEX_PAIRS(NSOBEX_TP)
*
*. Output 
*
      DIMENSION E(N_CC_AMP), F(*)
*
*. Scratch through input list
*
      DIMENSION CCVEC1(*),CCVEC2(*)
*     ^ Must hold CC vectors
      DIMENSION ISTST1(*),XSTST1(*),ISTST2(*),XSTST2(*)
      DIMENSION ISTST3(*),XSTST3(*),ISTST4(*),XSTST4(*)
*     ^ Must hold ST*ST maps for individual strings, given types,
*     all symmetries. 
      INTEGER IZ(*),IZSCR(*)
      INTEGER ISTOCC1(*),ISTOCC2(*)
*     ^ Hold occupations af all strings of given CAAB, SPGP and sym
      INTEGER ISTREO(*)
*     ^ Must hold reordering for all strings of given CAAB and SPGP
*
*. Local scratch
*
      INTEGER IJOCC(4*MXPNGAS)
*
      CALL LUCIAQENTER('HEXP_T')
      NTEST = 000
      IF (NTEST.GE. 50) THEN
*
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' Welcome to HExp_T '
        WRITE(6,*) ' ================ '
*
        WRITE(6,*) ' I_DO_ZERO, I_DO_FTERMS = ',
     &               I_DO_ZERO, I_DO_FTERMS
        WRITE(6,*) ' NSOBEX_TP = ', NSOBEX_TP
      END IF
*
      CALL REWINO(LUT)
*
      N_CC_AMPP1 = N_CC_AMP + 1
      IF(I_DO_ZERO.EQ.1) THEN
*. set initial operator to the unit operator
        ZERO = 0.0D0
        CALL SETVEC(E,ZERO,N_CC_AMPP1)
        E(N_CC_AMPP1) = 1.0D0
      END IF
*     ^ End if zeroing was required
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial E-vector '
        CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
      END IF
    
*  N = 0 term : Current E-term
*  (N = 1)  
*           Copy E to CCVEC1
*           T' + T' CCVEC1 in CCVEC2
*           add  CCVEC2 to E
*           copy CCVEC2 to CCVEC1
*  (N > 1)  T' CCVEC1 in CCVEC2
*           copy CCVEC2 to CCVEC1
*           Add  1/n! CCVEC2 to CCVEC3
*           Goto to next N if norm of CCVEC2 is gt. 0
      DO ITP_ORIG = 1, NSOBEX_TP    
      ITP = IEXEORD(ITP_ORIG)
C     ITP = ITP_ORIG
*.   Trying with reordered T again 
C?    WRITE(6,*) ' ITP, ISPOBEX_AC(ITP) = ', 
C?   &             ITP, ISPOBEX_AC(ITP)
      IF(ISPOBEX_AC(ITP).EQ.1) THEN
        ZERO = 0.0D0
        IF(I_DO_FTERMS.EQ.1) THEN
*. Offsets for current F-vector 
          CALL OFFSET_FOR_FTERM(ITP,NSOBEX_TP,ISOBEX_TP,LSOBEX_TP,NGAS,
     &         IB_F,LEN_F,0)
C             OFFSET_FOR_FTERM(IFTP,NSPOBEX,ISPOBEX,LSPOBEX,NGAS,
C    &           IB_F,LEN_F,IONLY_LEN)
          CALL SETVEC(F,ZERO,LEN_F)        
        ELSE
          LEN_F = 0
        END IF
*. Read in T-block at end of F-vector ( Memory has been defined 
*. so F-vector and T may be stored together as large T-block
*. usually is associated with small F-vector
        LEN_ITP = LSOBEX_TP(ITP)
        IB_T = LEN_F + 1
        CALL IFRMDS(LEN,1,-1,LUT)
C            FRMDSCO(ARRAY,NDIM,MBLOCK,IFILE,IMZERO)
        CALL FRMDSCO(F(IB_T),LEN,-1,LUT,IMZERO)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*)
          WRITE(6,*) ' ===================='
          WRITE(6,*) ' Exp_T for ITP = ', ITP
          WRITE(6,*) ' ===================='
          WRITE(6,*)
        END IF
*
C       IB_ITP = IBSOBEX_TP(ITP)
        IF(NTEST.GE.1000) WRITE(6,*) ' LEN_ITP =',     LEN_ITP
*. Loop over N from  1
        N = 0 
        FACN = 0
 1000   CONTINUE
         N = N + 1
         IF(NTEST.GE.1000) WRITE(6,*) ' order N = ' , N
         IF(N.EQ.1) THEN
          FACN = 1.0D0
         ELSE
          FACN = FACN*DBLE(N)
         END IF
         FACNI = 1.0D0/FACN
* 
         IF(N.EQ.1) THEN
          CALL COPVEC(E,CCVEC1,N_CC_AMPP1)
         END IF
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' CCVEC1 '
           CALL WRTMAT(CCVEC1,1,N_CC_AMPP1,1,N_CC_AMPP1)
         END IF
*. T CCVEC1 in CCVEC2
         ZERO = 0.0D0
         CALL SETVEC(CCVEC2,ZERO,N_CC_AMPP1)
*. Loop over Blocks in CCVEC1
         DO JTP = 1, NSOBEX_TP+1
C!       IF(ISPOBEX_AC(JTP).EQ.1) THEN
           IF(NTEST.GE.1000) WRITE(6,*) ' JTP = ', JTP
*. Spin-orbital type corresponding to ITP*JTP
*. Is ITP * JTP in active 
           CALL PROD_SPOB_EX_TP(ITP,JTP,IJTP,IEXCLEVEL,IJ_IS_ZERO)
           IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' ITP, JTP, IJTP, IJ_IS_ZERO = ', 
     &                  ITP, JTP, IJTP, IJ_IS_ZERO
           END IF
           IB_JTP = IBSOBEX_TP(JTP)
           IF(IJTP.NE.0) THEN
*. Occupation of the four strings in ITP*JTP
C                 T1T2_TO_T12_OCC(I1,I2,I12,NGAS)
             CALL T1T2_TO_T12_OCC(ISOBEX_TP(1,ITP),ISOBEX_TP(1,JTP),
     &                            IJOCC,NGAS)
*. T(ITP) CCVEC1(JTP) to be stored in CCVEC2(IJTP)
             IB_IJTP = IBSOBEX_TP(IJTP)
             IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' E-vector before T1T2... '
               CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
             END IF
             IF(JTP.EQ.NSOBEX_TP+1) THEN
*. JTP is unit operator, so just do a DAXPY : T12 = T12 + T(ITP)*Constant
               FACTOR = CCVEC1(IB_JTP)
               ONE = 1.0D0
               LEN_IJTP = LSOBEX_TP(IJTP)
               CALL VECSUM(CCVEC2(IB_IJTP),CCVEC2(IB_IJTP),F(IB_T),
     &                     ONE, FACTOR, LEN_IJTP)
             ELSE 
               CALL T1T2_TO_T12N(ISOBEX_TP(1,ITP),1,F(IB_T),
     &                          ISOBEX_TP(1,JTP),1,CCVEC1(IB_JTP),
     &                          IJOCC,IT12SM,CCVEC2(IB_IJTP),
     &                          ISTOCC1,ISTOCC2,ISTREO,
     &                          ISTST1,XSTST1,ISTST2,XSTST2,
     &                          ISTST3,XSTST3,ISTST4,XSTST4,
     &                          IZ,IZSCR)
             END IF
             IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' E-vector after  T1T2... '
               CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
             END IF
           ELSE IF(IEXCLEVEL.LE.2 .AND. I_DO_FTERMS.EQ.1
     &             .AND. IJ_IS_ZERO.EQ.0                ) THEN
*. Add to F terms
             ONE = 1.0D0
             LEN_JTP = LSOBEX_TP(JTP)
*
             IB_JTP_F = IB_F(JTP)
             IF(IB_JTP_F.GT.0) 
     &       CALL VECSUM(F(IB_JTP_F),F(IB_JTP_F),CCVEC1(IB_JTP),
     &                   ONE,FACNI,LEN_JTP)
           END IF
*          ^ End of tests if terms should be included
C!       END IF
C!*.       ^ End if Block JTP is active
         END DO
*        ^ End of loop over blocks JTP in e(k)
         XNORM = INPROD(CCVEC2,CCVEC2,N_CC_AMPP1)
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' CCVEC2 : '
           CALL WRTMAT(CCVEC2,1,N_CC_AMPP1,1,N_CC_AMPP1)
         END IF
         ONE = 1.0D0
         CALL VECSUM(E,E,CCVEC2,ONE,FACNI,N_CC_AMPP1)
         CALL COPVEC(CCVEC2,CCVEC1,N_CC_AMPP1)
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Updated E vector for N = ', N
           CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
         END IF
*. Temporary break
         MAXTRM = 20
         XTEST = 0.0D0    
        IF(N.EQ.MAXTRM) THEN
          DO K = 1, 1000
           WRITE(6,*) ' Termination due to N=MAXTRM test '
          END DO
        END IF
        IF(XNORM.GT.XTEST.AND.N.LE.MAXTRM) GOTO 1000
*     ^ End of loop over operators T(ITP)
*. An F-vector has been constructed, contract with Hamiltonian
C     H_TF_TERM(NSPOBEX_TP,IBSPOBEX_TP,LSPOBEX_TP,
C    &           ISPOBEX_TP,N_CC_AMP,T,F,HEXPT,ISPOBEX_AC,
C    &           N_TDL_MAX,ISPOBEX_FRZ,ITTP)          
      XFNORM = INPROD(F,F,LEN_F)
      I_DO_UNITOP = 0
      IUSE_TR2 = 0
      IF(XFNORM.GT.0.0D0)     
     &CALL H_TF_TERM(NSOBEX_TP,IBSOBEX_TP,LSOBEX_TP,
     &          ISOBEX_TP,N_CC_AMP,F(IB_T),F,HEXPT,ISPOBEX_AC,
     &          ISPOBEX_FRZ,ITP,I_DO_UNITOP,ISOBEX_PAIRS,IUSE_TR,
     &          IUSE_TR2,IB_F)
*
      ELSE
*. T-operator was not active, skip T block
        CALL IFRMDS(LEN,1,-1,LUT)
        CALL FRMDSCO(F(1),LEN,-1,LUT,IMZERO)
      END IF
*.    ^ End if operator is active
      END DO
*. Contract H 1 (1+E) with Hamiltonian
*. Add Unit operator to T
      ITP_UNI = NSOBEX_TP + 1
      IB_UNI = IBSOBEX_TP(ITP_UNI)
C     T(IB_UNI) = 1.0D0
C     F(IB_UNI) = 1.0D0
      ONE = 1.0D0
      I_DO_UNITOP = 1
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' E-vector : '
        CALL WRTBLKN(E,NSOBEX_TP+1,LSOBEX_TP)
      END IF
C?    WRITE(6,*) ' H times E '
      IF(IUSE_TR.EQ.1) IUSE_TR2 = 1
      CALL H_TF_TERM(NSOBEX_TP,IBSOBEX_TP,LSOBEX_TP,
     &          ISOBEX_TP,N_CC_AMP,ONE,E,HEXPT,ISPOBEX_AC,
     &          ISPOBEX_FRZ,ITP_UNI,I_DO_UNITOP,ISOBEX_PAIRS,IUSE_TR,
     &          IUSE_TR2,IBSOBEX_TP)
C     F(IB_UNI) = 0.0D0
*. Add core energy
      ONE = 1.0D0
      CALL VECSUM(HEXPT,HEXPT,E,ONE,ECORE,N_CC_AMP+1)
*. Alternative evaluation of energy = <E|H|CC > / <E|CC>
CMO   E_CC2 = INPROD(E,HEXPT,N_CC_AMP+1)/INPROD(E,E,N_CC_AMP+1)
CMO   WRITE(6,'(A,F22.15)') 
CMO  &' Alternative evaluation of E_CC = ', E_CC2
*
      IF(NTEST.GE. 50) THEN
        WRITE(6,*) ' H  EXP_T : '
        CALL WRTBLKN(HEXPT,NSOBEX_TP+1,LSOBEX_TP)
      END IF
*
C     STOP ' Stop at end of H Exp_T '
      CALL LUCIAQEXIT('HEXP_T')
      RETURN
      END
      SUBROUTINE H_TF_TERM(NSPOBEX_TP,IBSPOBEX_TP,LSPOBEX_TP,
     &           ISPOBEX_TP,N_CC_AMP,T,F,HEXPT,ISPOBEX_AC,
     &           ISPOBEX_FRZ,ITTP,I_DO_UNITOP,ISPOBEX_PAIRS,IUSE_TR,
     &           IUSE_TR2,IB_F)          
*
* Evaluate H T(ITTP) F and ADD to HEXPT
*
* T block is assumed to be delivered directly 
*
* (Remember to initialize HEXPT outside !)
*
* Jeppe Olsen, June 2001
*
* IUSE_TR2 : Extended use of time-reversal symmetry- 
*            assumes that F also enjoys time reversal symmetry
*
*
*. Note Extended spinorbital list is used where spinorbital excitation
*. NSPOBEX_TP + 1 is unit operator 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
C     INCLUDE 'ctcc.inc'
*
*. Input through argument list
*
      INTEGER IBSPOBEX_TP(*),LSPOBEX_TP(*),ISPOBEX_TP(4*NGAS,*)
      INTEGER ISPOBEX_AC(*)
      INTEGER ISPOBEX_FRZ(*)
      INTEGER ISPOBEX_PAIRS(*)
      DIMENSION F(*),T(*)
      INTEGER IB_F(*)
*
*. Local scratch
*
      INTEGER IOCC_H_AR(4*MXPNGAS)
      INTEGER IHINDX(4)
C     INTEGER IOCC_HT_AR(4*MXPNGAS)
      INTEGER IOCC_HTF_AR(4*MXPNGAS)
      INTEGER IOCC_T1_AR(4*MXPNGAS)
      INTEGER IOCC_T2_AR(4*MXPNGAS)
      INTEGER IOCC_T12_AR(4*MXPNGAS)
*. Remaining scratch is delivered through pointers in common block
*
*. Output
*
      DIMENSION HEXPT(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'H_TF_TE') 
*
*. Use time-reversal symmetry to simplify calculations
      IDUM = 0
      IZERO = 0
*
      NTEST = 000
      IF(NTEST.GE.50) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' Welcome to H_TF_TERM '
        WRITE(6,*) ' ====================='
        WRITE(6,*)
*. \/ Cannot be used with contracted F-vector
C       WRITE(6,*) ' Initial F-vector'
C       CALL WRTBLKN(F,NSPOBEX_TP+1,LSPOBEX_TP)
*. \/ cannot be used as complete T-vector is not delivered
C       WRITE(6,*) ' Initial T-vector'
C       CALL WRTBLKN(F,NSPOBEX_TP+1,LSPOBEX_TP)
        WRITE(6,*) ' Type of T = ', ITTP
      END IF
*
*. Number and types of operators in one- and two-electron operators
*
*. (Could and should be moved outside)
      CALL H_TYPES(0,N1TP,N2TP,IDUM,IDUM,IDUM,
     &             N1OBTP,N2OBTP,IDUM,IDUM,
     &             IDUM,IDUM,IDUM,IDUM,IDUM)
      N12TP = N1TP + N2TP
      LEN = 4*NGAS*(N1TP+N2TP)
      CALL MEMMAN(KLHTP,LEN,'ADDL  ',1,'H_TYP ')
      CALL MEMMAN(KLHINDX,LEN,'ADDL ',1,'HINDX ')
      CALL MEMMAN(KLHSIGN,LEN,'ADDL ',1,'HSIGN ')
      LEN = (N1OBTP + N2OBTP)*2*NGAS
      CALL MEMMAN(KLHOBTP,LEN,'ADDL  ',1,'HOBTP ')
      N12OBTP = N1OBTP + N2OBTP
      CALL MEMMAN(KLNSOX_FOR_OX_H,N12OBTP,'ADDL  ',1,'NSOX_H')
      CALL MEMMAN(KLISOX_TO_OX_H,N12TP,'ADDL  ',1,'ISOXTH') 
      CALL MEMMAN(KLISOX_FOR_OX_H,N12TP,'ADDL  ',1,'ISOXFH') 
      CALL MEMMAN(KLIBSOX_FOR_OX_H,N12OBTP,'ADDL  ',1,'IBSOXH') 
      CALL MEMMAN(KLSOX_SPFLIP,N12TP,'ADDL  ',1,'HSPFLP') 
      CALL MEMMAN(KLH_EXC2,N12TP,'ADDL  ',1,'H_EXC2')
      CALL H_TYPES(1,N1TP,N2TP,WORK(KLHTP),
     &             WORK(KLHINDX),WORK(KLHSIGN),
     &             N1OBTP,N2OBTP,WORK(KLHOBTP),
     &             WORK(KLNSOX_FOR_OX_H),WORK(KLISOX_TO_OX_H),
     &             WORK(KLISOX_FOR_OX_H),WORK(KLIBSOX_FOR_OX_H),
     &             WORK(KLSOX_SPFLIP),WORK(KLH_EXC2) )
*
      DO IHTP = 1, N1TP+N2TP
*. occupation of IHTP 
        CALL ICOPVE2(WORK(KLHTP),(IHTP-1)*4*NGAS+1,4*NGAS,IOCC_H_AR)
*. sign and indeces 
        CALL ICOPVE2(WORK(KLHINDX),(IHTP-1)*4+1,4,IHINDX)
        CALL ICOPVE2(WORK(KLHSIGN),IHTP,1,IHSIGN)
*. Occupation of T
        CALL ICOPVE(ISPOBEX_TP(1,ITTP),IOCC_T1_AR,4*NGAS)
*. And types in F
        IF(I_DO_UNITOP .EQ. 0 ) THEN
          NOP_TOT =  NSPOBEX_TP
        ELSE 
          NOP_TOT = NSPOBEX_TP + 1
        END IF
        DO IFTP = 1, NOP_TOT
C        IB_FTP = IBSPOBEX_TP(IFTP) 
         IB_FTP = IB_F(IFTP)
         IF(IB_FTP.GT.0) THEN
*. This F(IFTP) is active
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IHTP, ITTP, IFTP = ', IHTP, ITTP,IFTP 
*. Occupation of F operator 
          CALL ICOPVE(ISPOBEX_TP(1,IFTP),IOCC_T2_AR,4*NGAS)
*. Occupation of TF operator
          CALL OP_T_OCC(IOCC_T1_AR,IOCC_T2_AR,IOCC_T12_AR,IMZERO_TF)
*. Occupation of HTF operator
          CALL OP_T_OCC(IOCC_H_AR,IOCC_T12_AR,IOCC_HTF_AR,IMZERO_HTF)
*. And type of HTF
          CALL INUM_FOR_OCC(IOCC_HTF_AR,IHTFTP)
       
          IF(NTEST.GE.1000)  THEN
            WRITE(6,*) ' IHTFTP  = ', IHTFTP 
            IF(IHTFTP.GT.0) THEN
              WRITE(6,*) ' IOCC_HTF_AR ' 
              CALL WRT_SPOX_TP(IOCC_HTF_AR,1)
            END IF
          END IF
*
          IF(IHTFTP.GT.0.AND.
     &       IUSE_TR.EQ.1.AND.ISPOBEX_PAIRS(IHTFTP).LT.0) THEN
C            WRITE(6,*) ' Loop skipped for IHTFTP = ', IHTFTP
            IHTFTP = 0
          END IF 
*
          HFAC = 1.0D0
          IF(IHTFTP.GT.0.AND. IUSE_TR.EQ.1 .AND. 
     &       ISPOBEX_PAIRS(IHTFTP).EQ.IHTFTP) THEN  
*. The HTF block goes into itself by spinflip,      
*. Include if : 1) type of Hamiltonian is smaller than or equal to type of 
*                  spin flipped Hamiltonian 
*E              2) Type of Hamiltonian equals type of spinflipped 
*E                 Hamiltonian, and type of TF is greater than 
*E                 type of spinflipped TF
*. Type 2 does not work in general as we will not have that 
*.  F(I,J) = F(SPFLIP(I),SPFLIP(J))
*. However, for the term H E, the E-term has TR symmetry and 
*. this simplification may then be used
*. Extract type of spinflipped Hamiltonian   
            IH_SPINFLIP_TP = IFRMR(WORK(KLSOX_SPFLIP),1,IHTP)
            IF(IHTP.LT.IH_SPINFLIP_TP) THEN
*. In the symmetrizer for the blocks that goes into themselve 
*. by spin-flip, there is a factor 0.5. For terms where 
*. only half of the elements are calculated, one must multiply 
* with 2.0 to counteract this factor
               HFAC = 2.0D0
            ELSE IF (IUSE_TR2.EQ.1.AND.IHTP .EQ. IH_SPINFLIP_TP) THEN
*. At the moment this part is only entered for the H 1 E term, and for this 
*. term it does not eliminate any computations
*
              IT_FLIP_TP = ABS(ISPOBEX_PAIRS(ITTP))
              IF_FLIP_TP = ABS(ISPOBEX_PAIRS(IFTP))
*
              IF( ITTP.GT.IT_FLIP_TP.OR.
     &           (ITTP.EQ.IT_FLIP_TP.AND.IFTP.GT.IF_FLIP_TP))THEN
                  HFAC = 2.0D0
              ELSE IF (ITTP.LT.IT_FLIP_TP.OR.
     &           (ITTP.EQ.IT_FLIP_TP.AND.IFTP.LT.IF_FLIP_TP))THEN
                 IHTFTP = 0
              END IF
            ELSE IF (IHTP .GT. IH_SPINFLIP_TP) THEN
               IHTFTP = 0
            END IF
*.          ^ End of the various relations between IHTP and IH_SPINFLIP_TP
          END IF
*.              
          IF(IHTFTP.GT.0.AND.IMZERO_TF.EQ.0.AND.IMZERO_HTF.EQ.0
     &       .AND. ISPOBEX_FRZ(IHTFTP).EQ.0) THEN
*. Inside bounds, calculate H T F 
*. offsets in TCC arrays
C           IB_TTP = IBSPOBEX_TP(ITTP)
            IB_TTP = 1
            IB_HTFTP = IBSPOBEX_TP(IHTFTP)
*
            LENGTH_T = LSPOBEX_TP(ITTP)
            LENGTH_F = LSPOBEX_TP(IFTP)
            LENGTH_HT1T2 = LSPOBEX_TP(IHTFTP)
*      
            XFNORM = INPROD(F(IB_FTP),F(IB_FTP),LENGTH_F)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' IFTP,LENGTH_F,IB_FTP,IB_TTP, XFNORM =',
     &                     IFTP,LENGTH_F,IB_FTP,IB_TTP, XFNORM
            END IF
*
            IF(XFNORM.GT.0.0D0 ) THEN
*. Note IHSIGN is not used !
              CALL H_T1T2(IOCC_H_AR,IHINDX,HFAC,1,
     &             IOCC_T1_AR,1,T(IB_TTP),
     &             IOCC_T2_AR,1,F(IB_FTP),HEXPT(IB_HTFTP),
     &                 LENGTH_T,LENGTH_F,LENGTH_HT1T2)
            END IF
*           ^ End if F blocks was nonvanishing
          END IF
*         ^ End if HTF was inside bounds
         END IF
*        ^ End if F type was active
        END DO
*.      ^ End of loop over types in F
        IF(NTEST.GE.50) THEN
          WRITE(6,*)
          WRITE(6,*) ' Updated HEXPT after H op : ',  IHTP 
          WRITE(6,*) ' ================================== '
          WRITE(6,*)
          CALL WRTBLKN(HEXPT,NSPOBEX_TP+1,LSPOBEX_TP)
        END IF
      END DO
*     ^ End of loop over H types
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' Updated HEXPT  in string form '
        WRITE(6,*) ' ============================= '
        WRITE(6,*)
        CALL WRTBLKN(HEXPT,NSPOBEX_TP+1,LSPOBEX_TP)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'H_TF_TE') 
      RETURN
      END
      SUBROUTINE DIM_CNTR(ICONT,NCONT,LDIM,NCNTR)
*
* Find dimension, NCNTR,  
* of contraction operator ICONT having NCONT operators 
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER ICONT(LDIM,4)
*. Local scratch, assuming contraction is atmost 2-e oprator 
      INTEGER IUSED(4)
*
      IZERO = 0
      CALL ISETVC(IUSED,IZERO,NCONT)
      N = 1
      DO IOP = 1, NCONT
      IF(IUSED(IOP).EQ.0) THEN
*. other operators of the same GAS, SPIN AND CA ?
        NOP = 1
        DO JOP = IOP+1,NCONT
          IF(ICONT(IOP,1).EQ.ICONT(JOP,1).AND.
     &       ICONT(IOP,2).EQ.ICONT(JOP,2).AND.
     &       ICONT(IOP,3).EQ.ICONT(JOP,3)     )THEN
               NOP = NOP + 1
               IUSED(JOP) = 1
          END IF
        END DO
        LOB = NOBPT(ICONT(IOP,1))
        IDIV = 1
        DO I = 1, NOP
          N = N*(LOB+1-I)
          IDIV = IDIV * I
        END DO
        N = N/IDIV
      END IF
      END DO
*
      NCNTR = N
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Gasspace  Spin    Cr/An '
        WRITE(6,*) ' ========================'
        DO JCONT = 1, NCONT
          WRITE(6,'(I4,6X,I2,7X,I2)') 
     &    ICONT(JCONT,1),ICONT(JCONT,2),ICONT(JCONT,3)
        END DO
        WRITE(6,*) ' Number of contractions = ', NCNTR
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_SPOBEX_FOR_OBEX(ISPOBEX,IOBEX,NSPOBEX)
*
* Get all spinorbital excitations for orbital excitation IOBEX
*
* Jeppe Olsen, Winter of 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
*. Output 
      INTEGER ISPOBEX(*)
*
      CALL GET_SPOBEX_FOR_OBEXS(ISPOBEX,IOBEX,NSPOBEX,NSPOBEX_TP,
     &                         WORK(KLSOX_TO_OX))
*
      RETURN
      END 
      SUBROUTINE GET_SPOBEX_FOR_OBEXS(ISPOBEX,IOBEX,NSPOBEX,
     &           NSPOBEX_TP,ISOX_TO_OX)
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISOX_TO_OX(NSPOBEX_TP)
*. Output
      INTEGER ISPOBEX(*)
*
      NSPOBEX = 0
      DO JSPOBEX = 1, NSPOBEX_TP
        IF(ISOX_TO_OX(JSPOBEX).EQ.IOBEX) THEN
          NSPOBEX = NSPOBEX + 1
          ISPOBEX(NSPOBEX) = JSPOBEX 
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
       WRITE(6,*) ' Input orbital type = ', IOBEX
       WRITE(6,*) 
     & ' Corresponding number of spinorbital-excit =',NSPOBEX
       WRITE(6,*) '  Corresponding spinorbital-excitations : '
       CALL IWRTMA(ISPOBEX,1,NSPOBEX,1,NSPOBEX)
      END IF
*
      RETURN
      END
      SUBROUTINE INUM_FOR_OBEXC(IOBEXC,INUM)
*
* Number for orbital excitation given by IOBOCC
*
* Jeppe Olsen, Winter of 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc' 
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
*. Specific input 
      INTEGER IOBEXC(2*NGAS)
*
      CALL INUM_FOR_OBEXCS(IOBEXC,INUM,WORK(KOBEX_TP),NOBEX_TP,NGAS) 
*
      RETURN
      END
      SUBROUTINE INUM_FOR_OBEXCS(IOBEX,INUM,IGOBEX,NOBEX,NGAS)
*
* Number for orbital excitation given by IOBEXC, slave part 
*
* Jeppe Olsen, Winter of 2000
*. All orbital excitations
      INTEGER IGOBEX(2*NGAS,NOBEX)
*. Orbital excitation of interest
      INTEGER IOBEX(2*NGAS)
*
      INUM = 0
C?    WRITE(6,*) ' NGAS,NOBEX = ', NGAS,NOBEX
      DO JOBEX = 1, NOBEX
C?     WRITE(6,*) ' Orbital excitation ', JOBEX
C?     CALL IWRTMA(IGOBEX(1,JOBEX),1,2*NGAS,1,2*NGAS)
       IEQUAL = 1
       DO JGAS_CA = 1, 2*NGAS
         IF(IOBEX(JGAS_CA).NE.IGOBEX(JGAS_CA,JOBEX)) IEQUAL = 0
C?       WRITE(6,*) ' IOBEX(JGAS_CA), IGOBEX(JGAS_CA,JOBEX) ',
C?   &                IOBEX(JGAS_CA), IGOBEX(JGAS_CA,JOBEX)
C?       WRITE(6,*) ' JGAS_CA, IEQUAL = ', JGAS_CA, IEQUAL
       END DO
*
       IF(IEQUAL.EQ.1) INUM = JOBEX
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Orbital excitation in CA form '
        CALL IWRTMA(IOBEX,1,2*NGAS,1,2*NGAS)
        IF(INUM.EQ.0) THEN
          WRITE(6,*) ' Not included '
        ELSE
          WRITE(6,*) ' Number in list = ', INUM
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE ORD_CCEXTP(ICC_EXE_ORD)
*
* Find order of CC-excitation operators in which 

* Exp(T1) exp(T2) exp(T3) exp(T4) ....  should be executed.
*
* The right order is important for performing efficient MRCC
* calculations 
*
* The order is 
* Loop over # of creation operators in last particle space
* Loop over # of annihilation operators in first hole space
*
* Loop over  # of creation operators in last particle space - 1
* Loop over # of annihilation operators in second hole space
*
* ....
*
* Loop over # of creation operators in first particle space
* Loop over # of annihilation in last hole space 
*
*
* ICC_EXE_ORD(I) : Type of CC excitation that should be 
*                  used as I'th operator in the expansion 
*                  of Exp(T) 
* Jeppe Olsen, at the end of the GENCC rush, winter 2000
*
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
*. Output 
      INTEGER ICC_EXE_ORD(NSPOBEX_TP)
*. Local scratch
      INTEGER IEXC(MXPNGAS), IEXC1(MXPNGAS), MAXVAL_AR(MXPNGAS)
      INTEGER IEXC_CA(2*MXPNGAS)
*
*. The orbital excitations are first generated in the above 
*. order and these are subsequently expanded to spinorbital 
*. excitations
*
* The above order is generated by looping over integer with NGAS 
* integers arranged as 
* H space NHOLE, P space 1, H space NHOLE-1, P Space 2, .... 
*
      NTEST = 00
*
*, The current code does not work for open shell cases. 
*. As the current code furthermore is unused in general p.t.
*. we just return  the original order if there are any 
*. open shells in system
*
*. TEST : Enforce standard order 
      NVAL = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) NVAL = NVAL + 1
      END DO
      IF(NVAL.NE.0) THEN
        DO I = 1,  NSPOBEX_TP
           ICC_EXE_ORD(I)  = I 
        END DO
        I_AM_OKAY = 1
        GOTO 1001 
      END IF
*
      NPAR_SPC = NGAS - NHOL_SPC
      NCOM_SPC = MIN(NHOL_SPC,NPAR_SPC)
      CALL ISETVC(MAXVAL_AR,MX_EXC_LEVEL,NGAS)
      IFIRST = 1
      IOFF = 1
      NLOOP = 0
 1000 CONTINUE
        IF(IFIRST.EQ.1) THEN
          IZERO = 0
          CALL ISETVC(IEXC1,IZERO,NGAS)
          NONEW = 0
          IFIRST = 0
        ELSE 
C  NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
          CALL NXTNUM2(IEXC1,NGAS,0,MAXVAL_AR,NONEW)
        END IF
        NLOOP = NLOOP + 1
        IF(NONEW.EQ.0) THEN
*. Reorder to standard order of GAS spaces
*. Common spaces
         DO ICOM = 1, NCOM_SPC
           IEXC(NHOL_SPC+1-ICOM) = IEXC1((ICOM-1)*2+1)
           IEXC(NHOL_SPC+ICOM) = IEXC1((ICOM-1)*2+2)
         END DO
*. Remaining spaces, either particle or hole 
         IF(NHOL_SPC.GT.NPAR_SPC) THEN
           DO IHOL = 1,NGAS-2*NCOM_SPC
             IEXC(NHOL_SPC+1-NCOM_SPC-IHOL) = IEXC1(2*NCOM_SPC+IHOL)
           END DO
         ELSE IF(NPAR_SPC.GT.NHOL_SPC) THEN
           DO IPAR = 1,NGAS-2*NCOM_SPC
             IEXC(NHOL_SPC+NCOM_SPC+IPAR) = IEXC1(2*NCOM_SPC + IPAR)
           END DO
         END IF
*. Divide into creation and annihilation part
*
         IZERO = 0
         CALL ISETVC(IEXC_CA,IZERO,2*NGAS)
         DO JGAS = 1, NGAS
           IF(IHPVGAS(JGAS).EQ.1) THEN
             IEXC_CA(NGAS+JGAS) = IEXC(JGAS)
           ELSE IF(IHPVGAS(JGAS).EQ.2) THEN
             IEXC_CA(JGAS) = IEXC(JGAS)
           END IF
         END DO
*
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) 
     &     ' Next orbital excitation : original, ordered, CA form'
           CALL IWRTMA(IEXC1,1,NGAS,1,NGAS)
           CALL IWRTMA(IEXC,1,NGAS,1,NGAS)
           CALL IWRTMA(IEXC_CA,1,2*NGAS,1,2*NGAS)
         END IF
*. Is this orbital excitation included ?
         CALL INUM_FOR_OBEXC(IEXC_CA,INUM)
         IF(INUM.NE.0) THEN
*. Included, find corresponding spinorbital excitations
C               GET_SPOBEX_FOR_OBEX(ISPOBEX,IOBEX,NSPOBEX)
           CALL GET_SPOBEX_FOR_OBEX(ICC_EXE_ORD(IOFF),INUM,NUM)
           IOFF = IOFF + NUM
         END IF
C        IF(NLOOP.LT.100000) GOTO 1000
         GOTO 1000
        END IF
*. End of loop over numbers
*. Test that all spinorbital excitations have been included once
      I_AM_OKAY = 1
      DO ITP = 1, NSPOBEX_TP 
        IHIT = 0
        DO JTP = 1, NSPOBEX_TP
          IF(ICC_EXE_ORD(JTP).EQ.ITP) IHIT = IHIT + 1
        END DO
        IF(IHIT.NE.1) THEN  
          WRITE(6,*) ' Problems with Reorganization '
          WRITE(6,*) ' Type = ', ITP,' Obtained ', IHIT, ' times '
        END IF
      END DO
*
 1001 CONTINUE
*
      IF(NTEST.GE.100.OR.I_AM_OKAY.EQ.0) THEN
        WRITE(6,*) ' Execution order of spinorbital excitations'
        CALL IWRTMA(ICC_EXE_ORD,1,NSPOBEX_TP,1,NSPOBEX_TP) 
*
        WRITE(6,*) ' Spinorbital excitations in execution order '
        DO ITP = 1, NSPOBEX_TP
          IOFF = KLSOBEX + (ICC_EXE_ORD(ITP)-1)*4*NGAS/2 
C              WRT_SPOX_TP(IEX_TP,NEX_TP)
          CALL WRT_SPOX_TP(WORK(IOFF),1)
        END DO
      END IF
*
      IF(I_AM_OKAY.EQ.0) STOP ' Problem in ORD_CCEXC '
*
      RETURN
      END
      SUBROUTINE INUM_FOR_OCC2(IOCC,INUM,NDIFF)
*
* An operator is defined by OCC(NGAS,4). 
* Obtain the type number of this operator in list of CC operators
* A -1 indicates a nontrivial excitation not included in the list
* NDIFF is the number of smallest number of operators needed 
* to bring IOCC into the space 
*
* Difference to INUM_FOR_OCC2 : NDIFF added
*
* Jeppe Olsen, November 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'    
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
*
      CALL INUM_FOR_OCC21(IOCC,INUM,NGAS,WORK(KLSOBEX),
     &                   NSPOBEX_TP,NDIFF)
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Coupled cluster operator '
        CALL WRT_SPOX_TP(IOCC,1)
        WRITE(6,*) ' Corresponding type number ', INUM 
        WRITE(6,*) ' NDIFF = ', NDIFF
      END IF
*
      RETURN
      END
      SUBROUTINE INUM_FOR_OCC21(IOCC,INUM,NGAS,ISPOBEX_TP,
     &                         NSPOBEX_TP,NDIFF)
*
* Type number of CAAB operator
*. ISPOBEX_TP is assumed to be extended list also containing 
*. zero-particle unit operator
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP+1)
      INTEGER IOCC(4*NGAS)
*
      INUM = -1
      DO ITP = 1, NSPOBEX_TP+1
        IDIFF = 0
        DO I = 1, 4*NGAS
         IDIFF = IDIFF + ABS(IOCC(I)-ISPOBEX_TP(I,ITP))
        END DO
        IF(ITP.EQ.1) THEN
          NDIFF = IDIFF
        ELSE 
          NDIFF = MIN(NDIFF,IDIFF)
        END IF
        IF(IDIFF.EQ.0) INUM = ITP
      END DO
*
      RETURN
      END
      SUBROUTINE OFFSET_FOR_FTERMS(NSPOBEX,ISPOBEX,LSPOBEX,NGAS,IB_F,
     &           LEN_F,IONLY_LEN)
*
* Find offsets for F-terms of EXP T Ref
* The F terms are sum_I T_I F_I, where I runs over 
* all types of spinorbital excitations. 
* However, only the terms in T_I F_I that is connected 
* with the excitation manifold by atmost an double 
* excitation contributes.
* LEN_F is total length of needed F-terms
* IF IONLY_LEN = 1, only LEN_F is obtained 
*
* Jeppe Olsen, Shaping the new CC code up, November 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER ISPOBEX(4*NGAS,NSPOBEX), LSPOBEX(NSPOBEX)
*. Output, address of spinorbital type ITTP in F_IFTP is given in 
*  IB_F(ITTP,IFTP)
      INTEGER IB_F(NSPOBEX,NSPOBEX)
*. Local scratch
      INTEGER ITF_OCC(4*MXPNGAS)
*
      NTEST = 000
C?    WRITE(6,*) ' NSPOBEX, NGAS = ', NSPOBEX, NGAS 
      IB = 1
      DO IFTP = 1, NSPOBEX
        DO ITTP = 1, NSPOBEX
C?        WRITE(6,*) ' IFTP, ITTP = ', IFTP, ITTP
*. Occupation of ITTP*IFTP
C              OP_T_OCC(IOPOCC,ITOCC,IOPTOCC,IMZERO)
          CALL OP_T_OCC(ISPOBEX(1,ITTP),ISPOBEX(1,IFTP),
     &                  ITF_OCC,IMZERO)
          IADD = 0
          IF(IMZERO.EQ.0) THEN
*. Address of ITF, or number of excitations required to 
*. bring it into space
C                INUM_FOR_OCC2(IOCC,INUM,NDIFF)
            CALL INUM_FOR_OCC2(ITF_OCC,INUM,NDIFF) 
*. Space will be allocated only if TF is not in space, 
*. and can be reached by atmost 4 operators
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' IFTP, ITTP, NDIFF = ', IFTP, ITTP, NDIFF
            IF(0.LT.NDIFF.AND.NDIFF.LE.4) IADD = 1
          END IF
          IF(IADD.EQ.1) THEN
            IF(IONLY_LEN.NE.1) IB_F(ITTP,IFTP) = IB
            IB = IB + LSPOBEX(ITTP)
          ELSE
            IF(IONLY_LEN.NE.1) IB_F(ITTP,IFTP) = -1
          END IF
*
        END DO
      END DO
*
      LEN_F = IB - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Total length of F-terms = ', LEN_F
        IF(IONLY_LEN.NE.1) THEN
          WRITE(6,*) ' Array of offsets for F-terms '
          WRITE(6,*) ' ============================ '
          CALL IWRTMA(IB_F,NSPOBEX,NSPOBEX,NSPOBEX,NSPOBEX)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE LEN_TCC(NCA,NCB,NAA,NAB,NCAAB,NSM)
*
* A TCC block is given with specified number of 
* CA, CB, AA, AB strings of each sym. 
*
* Jeppe Olsen, November 2000
*
* New version, reducing scaling from NSM**4 to NSM**2, Juy 2002, HNIE
*
      INCLUDE 'implicit.inc'
*. Input 
      INTEGER NCA(8),NCB(8),NAA(8),NAB(8)
      INCLUDE 'multd2h.inc'
*. Output
      INTEGER NCAAB(8)
*. Local scratch 
      INTEGER NC(8), NA(8)
*
      IZERO = 0
      CALL ISETVC(NCAAB,IZERO,NSM)
      CALL ISETVC(NC,IZERO,NSM)
      CALL ISETVC(NA,IZERO,NSM)
*. Number of creation and annihilation strings with different symmetry
      DO IA_SM = 1, NSM
        DO IB_SM = 1, NSM
          IAB_SM = MULTD2H(IA_SM,IB_SM) 
          NC(IAB_SM) = NC(IAB_SM) + NCA(IA_SM)*NCB(IB_SM)
          NA(IAB_SM) = NA(IAB_SM) + NAA(IA_SM)*NAB(IB_SM)
        END DO
      END DO
*. And then combination of C and A
      DO IC_SM = 1, NSM
        DO IA_SM = 1, NSM
          ISM = MULTD2H(IC_SM,IA_SM)
          NCAAB(ISM) = NCAAB(ISM) +  NC(IC_SM)*NA(IA_SM)
        END DO
      END DO
*     
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of CAAB amplitudes per sym '    
        CALL IWRTMA(NCAAB,1,NSM,1,NSM)  
      END IF
*
      RETURN 
      END
      SUBROUTINE DIM_TCC_OP(ICCOP,NTCC)
*
*. A CC operator is defined by ICCOP
*. Find number of amplitudes per symmetry
*
*. Jeppe Olsen, November 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER ICCOP(NGAS,4)
*. Output
      INTEGER NTCC(*)
*. Local scratch  
      INTEGER NSTROP(8,4)
*. Number of CA,CB,AA,AB strings per symmetry 
C  NST_SPGP(IOCC,NSTFSM)
      CALL NST_SPGP(ICCOP(1,1),NSTROP(1,1))
      CALL NST_SPGP(ICCOP(1,2),NSTROP(1,2))
      CALL NST_SPGP(ICCOP(1,3),NSTROP(1,3))
      CALL NST_SPGP(ICCOP(1,4),NSTROP(1,4))
*
C     LEN_TCC(NCA,NCB,NAA,NAB,NCAAB,NSM)
      CALL LEN_TCC(NSTROP(1,1),NSTROP(1,2),NSTROP(1,3),NSTROP(1,4),
     &             NTCC,NSMST)
*
      RETURN
      END
      SUBROUTINE WRT_EXP_CCOP(IEXOP,LEXOP)
*
* Write CC operator given in expanded form
*
      INCLUDE 'implicit.inc'
      INTEGER  IEXOP(3,LEXOP)

*
      WRITE(6,*) '    LOP  GAS  ICA  IAB '
      WRITE(6,*) ' ======================'
      DO IOP = 1, LEXOP
        WRITE(6,'(4(3X,I3))') 
     &  IOP, IEXOP(1,IOP), IEXOP(2,IOP), IEXOP(3,IOP)
      END DO
*
      RETURN
      END
      SUBROUTINE EXP_CC_OP(IOP,IEXOP,LEXOP,NGAS)
*
* Expand CC excitation operator 
*
* Jeppe Olsen, Late november 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOP(NGAS,4)
*. Ouput 
      INTEGER IEXOP(3,*)
* 
      ICAAB = 0
      LOP = 0
      DO ICA = 1, 2
        DO IAB = 1, 2
          ICAAB = ICAAB + 1
          DO IGAS = 1, NGAS
            NOP = IOP(IGAS,ICAAB)
            DO JOP = 1, NOP
              LOP = LOP + 1
              IEXOP(1,LOP) = IGAS
              IEXOP(2,LOP) = ICA
              IEXOP(3,LOP) = IAB
            END DO
          END DO
        END DO
      END DO
      LEXOP = LOP
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Operator in standard and expanded form '
        WRITE(6,*) ' ======================================='
        CALL WRT_SPOX_TP(IOP,1)
C            WRT_EXP_CCOP(IEXOP,LEXOP)
        CALL WRT_EXP_CCOP(IEXOP,LEXOP)
      END IF
*
      RETURN
      END
      SUBROUTINE T_DL_DIM(ISPOBEX,NSPOBEX,MAXOP,N_TDL_MAX,ND_BAT)
*
* Obtain Max dimension of T(D,L) array where D is 
* an excitation operator with 1-MAXOP operators, and DL 
* belongs to the set of operators. D is batched in batches of 
* length ND_BAT
*
* Jeppe Olsen, November 2
*
*. General Input 
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'csm.inc'
*. Specific input 
      INTEGER ISPOBEX(4*NGAS,NSPOBEX)
*. Local scratch 
      INTEGER IEXPOP(3,2*MXPLCCOP)
      INTEGER INUM(MAXOP)
      INTEGER ILOP(4*MXPNGAS)
      INTEGER IDOP(4*MXPNGAS)
      INTEGER NDOP(8),NLOP(8)
*
      N_TDL_MAX = 0
      DO JSPOBEX = 1, NSPOBEX
*. Expand operator to array giving GAS, ICA, IAB for each operator
C       EXP_CC_OP(IOP,IEXOP,LEXOP,NGAS)
        CALL EXP_CC_OP(ISPOBEX(1,JSPOBEX),IEXPOP,LOP,NGAS)
*. Loop over number of operators that can be removed
        MAXOP_EFF = MIN(MAXOP,LOP)
        DO LDOP = 0, MAXOP_EFF
         NEW = 1
*. Loop over different ways of extracting LDOP operators from the list 
 9999    CONTINUE
         IF(NEW.EQ.1) THEN
*. First number 
           DO I = 1, LDOP
            INUM(I) = I
           END DO
           NEW = 0
           NONEW = 0
         ELSE
C               NXTORD(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
           CALL NXTORD(INUM,LDOP,1,LOP,NONEW)
         END IF
         IF(NONEW.EQ.0) THEN
*. A new D operator was obtained, find corresponding L and D operator 
*. in compact form 
           CALL ICOPVE(ISPOBEX(1,JSPOBEX),ILOP,4*NGAS)
           IZERO = 0
           CALL ISETVC(IDOP,IZERO,4*NGAS)
           DO JOP = 1, LDOP
             JGAS = IEXPOP(1,INUM(JOP))
             JCA  = IEXPOP(2,INUM(JOP))
             JAB  = IEXPOP(3,INUM(JOP))
             JCAAB = (JCA-1)*2 + JAB
             IDOP(JGAS+(JCAAB-1)*NGAS) = IDOP(JGAS+(JCAAB-1)*NGAS) + 1
             ILOP(JGAS+(JCAAB-1)*NGAS) = ILOP(JGAS+(JCAAB-1)*NGAS) - 1
           END DO
*. Dimensions of DOP and LOP
C     DIM_TCC_OP(ICCOP,NTCC)
           CALL DIM_TCC_OP(IDOP,NDOP)
           CALL DIM_TCC_OP(ILOP,NLOP)
           MAX_D = IMNMX(NDOP,NSMST,2)
           MAX_D = MIN(MAX_D,ND_BAT)
COLD *. D does not need to be connected with above op, one is actually 
COLD *. setting up T(D2,L1)
C          MAX_D = ND_BAT
           MAX_L = IMNMX(NLOP,NSMST,2)
C?         WRITE(6,*) ' MAX_D, MAX_L = ', MAX_D, MAX_L
           LEN_DL= MAX_D*MAX_L
           N_TDL_MAX = MAX(N_TDL_MAX,LEN_DL)
           IF(LDOP.NE.0) GOTO 9999
*          ^ End of loop over different ways of removing NDOP operators
         END IF
*        ^ End if new operators were found
        END DO
*       ^ End of loop over NDOP
      END DO
*     ^ End of loop over JSPOBEX
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Max dimension of T(D,L) = ', N_TDL_MAX
      END IF
*
      RETURN
      END
      SUBROUTINE K_RES_DIM(ISPOBEX,NSPOBEX,MAXOP,NELMNT_MAX)
*
* A set of spinorbital excitations, ISPOBEX_TP are given 
*
* Find largest dimension of operators that can be obtained 
* by removing 0-MAXOP operators from these 
*
* Jeppe Olsen, Nov. 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER ISPOBEX(NGAS,4,NSPOBEX)
*. Local scratch
      INTEGER KSPOBEX(NGAS,4)
*
      NELMNT_MAX = 0
      DO JSPOBEX = 1, NSPOBEX
*. Number of elements in this TCC block
        NELMNT = 1
        DO ICAAB = 1,4
          DO IGAS = 1, NGAS
            NELMNT = 
     &      NELMNT*IBION(NOBPT(IGAS),ISPOBEX(IGAS,ICAAB,JSPOBEX))
          END DO
        END DO
        NELMNT_MAX = MAX(NELMNT,NELMNT_MAX)
*. This was the number of element in the unmodified type. 
*. Larger dimensions can be obtained only by removing operators from 
*. spaces that are more than half filled. 
*. Use this to locate the spaces where the largest gains are obtained 
*. from removing electrons, done one electron at a time 
        KELMNT_MAX = NELMNT
        CALL ICOPVE(ISPOBEX(1,1,JSPOBEX),KSPOBEX(1,1),4*NGAS)
        DO KOP = 1, MAXOP
*. Place where most will be gained by removing one operator
          KELMNTP = KELMNT_MAX
          JGAS = 0
          JCAAB = 0
          DO ICAAB = 1, 4
          DO IGAS = 1, NGAS
            IF(2*KSPOBEX(IGAS,ICAAB).GT.NOBPT(IGAS)) THEN
*. more than halffiled with operatos
              KELMNT = (KELMNTP*KSPOBEX(IGAS,ICAAB))
     &                /(NOBPT(IGAS)-KSPOBEX(IGAS,ICAAB)+1)
C             KELMNT = (KELMNTP*(NOBPT(IGAS)-KSPOBEX(IGAS,ICAAB)+1))/
C    &                  KSPOBEX(IGAS,ICAAB)
              IF(KELMNT.GT.KELMNT_MAX) THEN
                JGAS = IGAS
                JCAAB = ICAAB
                KELMNT_MAX  = KELMNT
              END IF
            END IF
          END DO
          END DO
          IF(JGAS.GT.0) THEN
            KSPOBEX(JGAS,JCAAB) = KSPOBEX(JGAS,JCAAB)-1
          ELSE
*. No more spaces more than halffilled with operators so
            GOTO 1001
          END IF
        END DO
*       ^ End of loop over KOP
 1001   CONTINUE
        NELMNT_MAX = MAX(NELMNT_MAX,KELMNT_MAX)
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Max number of operators to be removed ', MAXOP
        WRITE(6,*) ' Largest T block = ', NELMNT_MAX
      END IF
*
      RETURN
      END
      SUBROUTINE GET_CCSTR_FROM_LIST(ISTRING,INUM_CAAB,ISYM_CAAB,
     &           NOP_CAAB,LIST_CA,LIST_CB,LIST_AA,LIST_AB,
     &           IB_CA,IB_CB,IB_AA,IB_AB)
* Obtain CC string defined by INUM_CAAB, ISYM_CAAB
*
* Jeppe Olsen, May 2000 is running out
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER INUM_CAAB(4),ISYM_CAAB(4)
      INTEGER NOP_CAAB(4)
      INTEGER LIST_CA(*),LIST_CB(*),LIST_AA(*),LIST_AB(*)
      INTEGER IB_CA(*),IB_CB(*),IB_AA(*),IB_AB(*)
*. Output
      INTEGER ISTRING(*)
*. CA part
      IOFF = 1 
      IF(NOP_CAAB(1).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(1),ISYM_CAAB(1),
     &                       NOP_CAAB(1),LIST_CA,IB_CA)
      IOFF = IOFF + NOP_CAAB(1)
      END IF
*. CB 
      IF(NOP_CAAB(2).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(2),ISYM_CAAB(2),
     &                       NOP_CAAB(2),LIST_CB,IB_CB)
      IOFF = IOFF + NOP_CAAB(2)
      END IF
*. AA
      IF(NOP_CAAB(3).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(3),ISYM_CAAB(3),
     &                       NOP_CAAB(3),LIST_AA,IB_AA)
      IOFF = IOFF + NOP_CAAB(3)
      END IF
*. AA
      IF(NOP_CAAB(4).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(4),ISYM_CAAB(4),
     &                       NOP_CAAB(4),LIST_AB,IB_AB)
      IOFF = IOFF + NOP_CAAB(4)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        NOP_TOT = IOFF - 1
        WRITE(6,*) ' CAAB string from GET_CCSTR ..'
        CALL IWRTMA(ISTRING,1,NOP_TOT,1,NOP_TOT)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_STR_FROM_LIST(ISTRING,INUM,ISYM,NOP,LIST,IB)
*
* Obtain string INUM of sym ISYM from complete list of strings
*
* Jeppe Olsen, Still May 2000
*
      INCLUDE 'implicit.inc'
*.Input
      INTEGER LIST(NOP,*),IB(*)
*. Output 
      INTEGER ISTRING(NOP)
*
      INUM_ABS = IB(ISYM)-1+INUM
      CALL ICOPVE(LIST(1,INUM_ABS),ISTRING,NOP)
*
C?    IF(NOP.GT.0) WRITE(6,*) ' GET_STR... LIST(1,1) = ', LIST(1,1)
C?    IF(NOP.GT.0) WRITE(6,*) ' INUM_ABS = ', INUM_ABS
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' String number and sym = ', INUM,ISYM
        WRITE(6,*) ' Corresponding string : '
        CALL IWRTMA(ISTRING,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_OPINT(OPSCR,
     &IOD2X,LD2,INM_CAAB_D2,ISM_CAAB_D2,ISTR_D2,IBO2DX,
     &IOD1X,LD1,INM_CAAB_D1,ISM_CAAB_D1,ISTR_D1,IBO1DX,
     &IOEX ,LEX,INM_CAAB_EX,ISM_CAAB_EX,ISTR_EX,IBOEX,
     &IEXD1D2_INDX,FACX)
*
* Fetch batch of operator integrals O(D2,EX,D1)
*
* Connected with standard integral input
*
* Jeppe Olsen, May of 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
*. Input 
*.. Occupation of gas spaces         
      INTEGER IOD2X(NGAS,4),IOD1X(NGAS,4),IOEX(NGAS,4)
*. Index in EXD1D2 to original index in H
      INTEGER IEXD1D2_INDX(4)
*.. The sym and number of the included strings 
      INTEGER ISM_CAAB_EX(4,*),INM_CAAB_EX(4,*)
      INTEGER ISM_CAAB_D1(4,*),INM_CAAB_D1(4,*)
      INTEGER ISM_CAAB_D2(4,*),INM_CAAB_D2(4,*)
*.. The actual occupation of the various strings
      INTEGER ISTR_D1(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D2(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_EX(MX_ST_TSOSO_BLK_MX*NSMST,4)
*.. Start of strings with given symmetry
      INTEGER IBOEX(8,4),IBO1DX(8,4),IBO2DX(8,4)
* 
*. Local scratch
*
C     INTEGER IOCC(4),ICREA(4),IANNI(4)
      INTEGER JD1STR(4),JD2STR(4),JEXSTR(4)
      INTEGER IJKL_EDD(4), IJKL_ORIG(4)
      INTEGER NOP_D1_CAAB(4),NOP_D2_CAAB(4),NOP_EX_CAAB(4)
*. Output
      DIMENSION OPSCR(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GET_OPINT speaking : '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' IOEX, IOD2X and  IOD1X  '
        CALL WRT_SPOX_TP(IOEX,1)
        CALL WRT_SPOX_TP(IOD2X,1)
        CALL WRT_SPOX_TP(IOD1X,1)
        WRITE(6,*)
      END IF
*
*. Number of C/A A/B indeces 
      NOP_D1_CA = IELSUM(IOD1X(1,1),NGAS)
      NOP_D1_CB = IELSUM(IOD1X(1,2),NGAS)
      NOP_D1_AA = IELSUM(IOD1X(1,3),NGAS)
      NOP_D1_AB = IELSUM(IOD1X(1,4),NGAS)
      NOP_D1 = NOP_D1_CA+NOP_D1_CB+NOP_D1_AA+NOP_D1_AB
C?    WRITE(6,*) ' NOP_D1 = ', NOP_D1
      NOP_D1_CAAB(1) = NOP_D1_CA
      NOP_D1_CAAB(2) = NOP_D1_CB
      NOP_D1_CAAB(3) = NOP_D1_AA
      NOP_D1_CAAB(4) = NOP_D1_AB
*
      NOP_D2_CA = IELSUM(IOD2X(1,1),NGAS)
      NOP_D2_CB = IELSUM(IOD2X(1,2),NGAS)
      NOP_D2_AA = IELSUM(IOD2X(1,3),NGAS)
      NOP_D2_AB = IELSUM(IOD2X(1,4),NGAS)
      NOP_D2 = NOP_D2_CA+NOP_D2_CB+NOP_D2_AA+NOP_D2_AB
C?    WRITE(6,*) ' NOP_D2 = ', NOP_D2
      NOP_D2_CAAB(1) = NOP_D2_CA
      NOP_D2_CAAB(2) = NOP_D2_CB
      NOP_D2_CAAB(3) = NOP_D2_AA
      NOP_D2_CAAB(4) = NOP_D2_AB
*
      NOP_EX_CA = IELSUM(IOEX(1,1),NGAS)
      NOP_EX_CB = IELSUM(IOEX(1,2),NGAS)
      NOP_EX_AA = IELSUM(IOEX(1,3),NGAS)
      NOP_EX_AB = IELSUM(IOEX(1,4),NGAS)
      NOP_EX = NOP_EX_CA+NOP_EX_CB+NOP_EX_AA+NOP_EX_AB
C?    WRITE(6,*) ' NOP_EX = ', NOP_EX
      NOP_EX_CAAB(1) = NOP_EX_CA
      NOP_EX_CAAB(2) = NOP_EX_CB
      NOP_EX_CAAB(3) = NOP_EX_AA
      NOP_EX_CAAB(4) = NOP_EX_AB
*
      NCREA_ALPHA  = NOP_D1_CA + NOP_D2_CA + NOP_EX_CA
      NCREA_BETA   = NOP_D1_CB + NOP_D2_CB + NOP_EX_CB
      NANNI_ALPHA  = NOP_D1_AA + NOP_D2_AA + NOP_EX_AA
      NANNI_BETA   = NOP_D1_AB + NOP_D2_AB + NOP_EX_AB
*
C?    WRITE(6,*) ' NANNI_ALPHA, NANNI_BETA = ', 
C?   &             NANNI_ALPHA, NANNI_BETA
*
      NCREA = NCREA_ALPHA + NCREA_BETA
      NANNI = NANNI_ALPHA + NANNI_BETA
*
      NALPHA = NANNI_ALPHA + NCREA_ALPHA
      NBETA  = NANNI_BETA  + NCREA_BETA
*
      NOP = NCREA + NANNI
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original index '
        CALL IWRTMA(IEXD1D2_INDX,1,NOP,1,NOP)
      END IF
*
      IF(.NOT.((NCREA.EQ.1.AND.NANNI.EQ.1 ) .OR.
     &         (NCREA.EQ.2.AND.NANNI.EQ.2 )     )) THEN
        WRITE(6,*) ' Unknown operator in GET_OPINT '
        WRITE(6,*) ' NCREA, NANNI = ', NCREA, NANNI
        STOP       ' Unknown operator in GET_OPINT '
      END IF
*.
      IF(.NOT.(NCREA_ALPHA.EQ.NANNI_ALPHA.AND.
     &         NCREA_BETA.EQ.NANNI_BETA       ) ) THEN
        WRITE(6,*) ' Unknown operator in GET_OPINT '
        WRITE(6,*) 
     &  ' NCREA_ALPHA, NANNI_ALPHA, NCREA_BETA, NANNI_BETA  ', 
     &    NCREA_ALPHA, NANNI_ALPHA, NCREA_BETA, NANNI_BETA 
        STOP       ' Unknown operator in GET_OPINT '
      END IF
*. Now we know that a usual MS = 0, one- or two-electron 
*. operator is in operation
*. More complicated operators can be concocted, once I have
*. these routines debugged .... 
*
* We will constrict the D2 EX D1 operators as 
* CA_D2 CB_D2 AA_D2,AB_D2, CA_EX, CB_EX, AA_EX,AB_EX ...
*
      INT = 0
      DO JD1 = 1, LD1
*. Occupation of J1 string in JD1STR
      IF(NOP_D1.NE.0) CALL GET_CCSTR_FROM_LIST 
     &     (JD1STR,INM_CAAB_D1(1,JD1),ISM_CAAB_D1(1,JD1),
     &      NOP_D1_CAAB,ISTR_D1(1,1),ISTR_D1(1,2),ISTR_D1(1,3),
     &      ISTR_D1(1,4),IBO1DX(1,1),IBO1DX(1,2),IBO1DX(1,3),
     &      IBO1DX(1,4) )
      DO JEX = 1, LEX
      IF(NOP_EX.NE.0) CALL GET_CCSTR_FROM_LIST 
     &     (JEXSTR,INM_CAAB_EX(1,JEX),ISM_CAAB_EX(1,JEX),
     &      NOP_EX_CAAB,ISTR_EX(1,1),ISTR_EX(1,2),ISTR_EX(1,3),
     &      ISTR_EX(1,4),IBOEX(1,1),IBOEX(1,2),IBOEX(1,3),
     &      IBOEX(1,4) )
      DO JD2 = 1, LD2
      IF(NOP_D2.NE.0) CALL GET_CCSTR_FROM_LIST 
     &     (JD2STR,INM_CAAB_D2(1,JD2),ISM_CAAB_D2(1,JD2),
     &      NOP_D2_CAAB,ISTR_D2(1,1),ISTR_D2(1,2),ISTR_D2(1,3),
     &      ISTR_D2(1,4),IBO2DX(1,1),IBO2DX(1,2),IBO2DX(1,3),
     &      IBO2DX(1,4) )
*. Indeces of operator EXD1D2 
      DO IOP_EX = 1, NOP_EX
        IJKL_EDD(IOP_EX) = JEXSTR(IOP_EX)
      END DO
      DO IOP_D1 = 1, NOP_D1
        IJKL_EDD(NOP_EX+IOP_D1) = JD1STR(IOP_D1)
      END DO
      DO IOP_D2 = 1, NOP_D2
        IJKL_EDD(NOP_EX+NOP_D1+IOP_D2) = JD2STR(IOP_D2)
      END DO
*. Original order
      DO IOP = 1, NOP   
        IJKL_ORIG(IEXD1D2_INDX(IOP)) = IJKL_EDD(IOP)
      END DO
*
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' IJKL_EDD and IJKL_ORIG '
        CALL IWRTMA(IJKL_EDD,1,NOP,1,NOP)
        CALL IWRTMA(IJKL_ORIG,1,NOP,1,NOP)    
      END IF
*
      INT = INT + 1
      IF(NOP.EQ.2) THEN
*. One-electron integral 
        IF(IREFTYP.NE.2) THEN
*. Similarity transformed integrals in orbital basis 
          OPSCR(INT) = 
     &    GETH1_B(IJKL_ORIG(1),IJKL_ORIG(2))*FACX                 
        ELSE
*. Similarity transformed integrals in spinorbital basis 
          IF(NALPHA.EQ.2) THEN
           OPSCR(INT) = 
     &     GETH1_B2(IJKL_ORIG(1),IJKL_ORIG(2),WORK(KFI_AL))*FACX
          ELSE IF (NBETA.EQ.2) THEN
           OPSCR(INT) = 
     &     GETH1_B2(IJKL_ORIG(1),IJKL_ORIG(2),WORK(KFI_BE))*FACX
          END IF
        END IF
      ELSE IF(NOP.EQ.4) THEN
*. Two-electron integral
        IF(IREFTYP.NE.2) THEN
*. Similarity transformed integral in orbital basis
          IF(NALPHA.EQ.2.AND.NBETA.EQ.2) THEN
*. coulomb integral
            OPSCR(INT) = GTIJKL(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                          IJKL_ORIG(3))*FACX
          ELSE IF (NALPHA.EQ.4.OR.NBETA.EQ.4) THEN
*. coulomb - exchange  
            OPSCR(INT) = (GTIJKL(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                          IJKL_ORIG(3)) 
     &                 - GTIJKL(IJKL_ORIG(1),IJKL_ORIG(3),IJKL_ORIG(2),
     &                          IJKL_ORIG(4)))*FACX
          END IF
        ELSE IF (IREFTYP.EQ.2) THEN
           IF(NALPHA.EQ.2.AND.NBETA.EQ.2) THEN
*. A single two-electron integral 
            OPSCR(INT) 
     &    = GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                   IJKL_ORIG(3),2,2)*FACX
           ELSE IF (NALPHA.EQ.4) THEN
*. coulomb - exchange  
            OPSCR(INT) 
     &    = GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                   IJKL_ORIG(3),4,0) 
     &    - GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(3),IJKL_ORIG(2),
     &                   IJKL_ORIG(4),4,0)
            OPSCR(INT) = FACX*OPSCR(INT)
           ELSE IF (NBETA.EQ.4) THEN
*. coulomb - exchange  
            OPSCR(INT) 
     &    = (GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(4),IJKL_ORIG(2),
     &                   IJKL_ORIG(3),0,4) 
     &    - GTIJKL_SM_AB(IJKL_ORIG(1),IJKL_ORIG(3),IJKL_ORIG(2),
     &                   IJKL_ORIG(4),0,4)) * FACX
           END IF
*          ^ End of switch for NALPHA, NBETA 
         END IF
*        ^ End of IREFTYP switch 
      END IF
*     ^ End of NOP = 2 / 4 switch
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' INT, OPSCR(INT) = ', INT, OPSCR(INT)
      END IF
*
      END DO
      END DO
      END DO
*
      IF(NTEST.GE.100) THEN 
         WRITE(6,*)
         WRITE(6,*) ' Output matrix from GET_OPINT as X(D2EX,D1)' 
         WRITE(6,*) ' ========================================= '
         WRITE(6,*)
         CALL WRTMAT(OPSCR,LD2*LEX,LD1,LD2*LEX,LD1)
      END IF
*
      RETURN
      END 
      SUBROUTINE INUM_FOR_OCC(IOCC,INUM)
*
* An operator is defined by OCC(NGAS,4). 
* Obtain the type number of this operator in list of CC operators
* A vanishing type number indicates that IOCC  contains no operators 
* at all 
* A -1 indicates a nontrivial excitation not included in the list
*
* Jeppe Olsen, May 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'   
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
*
      CALL INUM_FOR_OCC1(IOCC,INUM,NGAS,WORK(KLSOBEX),
     &                   NSPOBEX_TP)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Coupled cluster operator '
        CALL WRT_SPOX_TP(IOCC,1)
        WRITE(6,*) ' Corresponding type number ', INUM 
      END IF
*
      RETURN
      END
      SUBROUTINE INUM_FOR_OCC1(IOCC,INUM,NGAS,ISPOBEX_TP,
     &                         NSPOBEX_TP)
*
* Type number of CAAB operator
*. ISPOBEX_TP is assumed to be extended list also containing 
*. zero-particle unit operator
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP+1)
      INTEGER IOCC(4*NGAS)
*
      INUM = -1
      DO ITP = 1, NSPOBEX_TP+1
        IDIFF = 0
        DO I = 1, 4*NGAS
         IDIFF = IDIFF + ABS(IOCC(I)-ISPOBEX_TP(I,ITP))
        END DO
        IF(IDIFF.EQ.0) INUM = ITP
      END DO
* 
COLD  IF(INUM.EQ.-1) THEN
*. Check to see if excitation is zero-particle unit operator 
COLD   IUNIT = 1
COLD   DO I = 1, 4*NGAS 
COLD     IF(IOCC(I).NE.0) IUNIT = 0
COLD   END DO
COLD   IF(IUNIT.EQ.1) INUM = 0
COLD  END IF
*
      RETURN
      END
      SUBROUTINE OP_T_OCC(IOPOCC,ITOCC,IOPTOCC,IMZERO)
*
* The occupation of  general operator OP and an CC cluster operator 
* are given. Find occupation of resulting operator after all 
* contractions have been performed. All deexcitation operators
* (anni of particles and creation of holes) are contracted. If 
*  this is not possible, Imzero = 1 is returned
*
* Jeppe Olsen, still May 2000
*
*. The only allowed open shell reference allowed is 
*  the high spin state with all singly occupied electrons having alpha spin
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER IOPOCC(NGAS,4),ITOCC(NGAS,4)
*. Output
      INTEGER IOPTOCC(NGAS,4)
*
      CALL ICOPVE(ITOCC,IOPTOCC,4*NGAS)
*. Loop over parts of operator
      IMZERO = 0
      DO IGAS = 1, NGAS
       DO IAB = 1, 2
        IF(IHPVGAS_AB(IGAS,IAB).EQ.1)THEN 
*. Contraction of hole creations in operator with annihilations in T
          IOPTOCC(IGAS,2+IAB) = IOPTOCC(IGAS,2+IAB)-IOPOCC(IGAS,IAB)
          IF( IOPTOCC(IGAS,2+IAB) .LT. 0 ) IMZERO = 1
*. Addition of hole annihilation operators in Op and T
          IOPTOCC(IGAS,2+IAB) = IOPTOCC(IGAS,2+IAB)+IOPOCC(IGAS,2+IAB)
        ELSE IF (IHPVGAS_AB(IGAS,IAB).EQ.2) THEN
*. Contraction of particle annihilations in operator with creations in T
          IOPTOCC(IGAS,IAB) = IOPTOCC(IGAS,IAB) - IOPOCC(IGAS,2+IAB)
          IF( IOPTOCC(IGAS,IAB) .LT. 0 ) IMZERO = 1
*. Addition of particle creation operators in Op and T
          IOPTOCC(IGAS,IAB) = IOPTOCC(IGAS,IAB) + IOPOCC(IGAS,IAB)
        END IF
       END DO 
      END DO
*. Ensure that number of creations/annihilations do not exceed 
*. number of orbitals
      DO ICAAB = 1, 4
        DO IGAS = 1, NGAS
          IF(IOPTOCC(IGAS,ICAAB).GT.NOBPT(IGAS)) IMZERO = 1 
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '----------------'
        WRITE(6,*) ' OP_T_OCC says:'
        WRITE(6,*) '----------------'
        WRITE(6,*) ' Combination of general operator and CC operator'
        WRITE(6,*) ' Input general and CC operators '
        CALL WRT_SPOX_TP(IOPOCC,1)
        CALL WRT_SPOX_TP(ITOCC,1)
        WRITE(6,*) ' Output operator '
        CALL WRT_SPOX_TP(IOPTOCC,1)
        IF(IMZERO.EQ.1) THEN
          WRITE(6,*) ' Vanishing operator '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE H_TYPES(IFLAG,N1TP,N2TP,IHTP,IHINDX,IHSGN,
     &                   N1OBTP,N2OBTP,IHOBTP,NSOX_FOR_OX_H,
     &                   ISOX_TO_OX_H,ISOX_FOR_OX_H,IBSOX_FOR_OX_H,
     &                   IHTP_SPINFLIP,IH_EXC2)
*
* Obtain number of spin-orbital types of one- and two-electron 
* operator. 
* If IFLAG = 1, the type  are constructed and saved in IHTP,IHINDX,IHSIGN
*
* IHTP : Types of the four operators 
* IHINDX : Place of each index in original Hamiltonian
* IHSGN  : Sign required to bring operator into  ICA ICB IAA IAB order
*
* NSOX_FOR_OX_H : Number of SOX for a given OX
* ISOX_TO_OX_H : Spinorbitalexcitation => orbitalexcitations
* ISOX_FOR_OX_H: SOX for a given OX
* IBSOX_FOR_OX_H: Base for SOX for a given OX 
* IHTP_SPINFLIP : Type 
* IH_EXC2 : Excitation level (times two) for each operator 
* 
*
* Jeppe Olsen, May 2000
*              Feb. 2001 : Orbital types added
*              June 2002 : IH_EXC2 added 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*
      INTEGER IHTP(4*NGAS,*), IHINDX(4,*), IHSGN(*)
*. Orbital excitations are given as C A
      INTEGER IHOBTP(2*NGAS,*)
*
      INTEGER NSOX_FOR_OX_H(*), ISOX_TO_OX_H(*) 
      INTEGER ISOX_FOR_OX_H(*), IBSOX_FOR_OX_H(*)
*. Spinflipped operator type
      INTEGER IHTP_SPINFLIP(*)
*. Excitation level times two ( Number of exci op - number of deexci op)
*. Scratch for one spinflipped op
      INTEGER IH_FLIP(4*MXPNGAS)
*. Array giving excitation level times two of the various types
      INTEGER IH_EXC2(*)
*
      IZERO = 0
*
* ============================================
*. 1 : Numbers of 1- and 2-electron operators
* ============================================
*
*. One-electron operators : alpha alpha , beta-beta for all pairs of gasspaces
      N1TP = 2 * NGAS*NGAS
*. One-electron orbital operatos
      N1OBTP = NGAS * NGAS
*. Two-electron operators : a+i alpha a+j alpha ak alpha al alpha 
*                           (i.ge.j, k.ge.l)
*                           a+i beta  a+j beta  ak beta  al beta  
*                           (i.ge.j, k.ge.l)
*                           a+i alpha a+j beta ak beta   al alpha 
*                           ( all i,j,k,l)
      N2TP = 2* (NGAS*(NGAS+1)/2)**2 + NGAS**4
      N2OBTP = NGAS**4
*
*
* ============================================
*. 2 : And the operators
* ============================================
*
      IF(IFLAG.EQ.1) THEN
*
        IZERO = 0
        CALL ISETVC(NSOX_FOR_OX_H,IZERO,N1OBTP+N2OBTP)
        INEXT = 1
        INEXT_OB = 1
*. One-electron operators
        DO IGAS = 1, NGAS
          DO JGAS = 1, NGAS
*. Orbital excitation :
          CALL ISETVC(IHOBTP(1,INEXT_OB),IZERO,2*NGAS)
          IHOBTP(IGAS+NGAS,INEXT_OB)  = 1
          IHOBTP(JGAS     ,INEXT_OB)  = 1
          IELMNT = 2
          NSOX_FOR_OX_H(INEXT_OB) = IELMNT
          IBSOX_FOR_OX_H(INEXT_OB) = INEXT
          DO J = INEXT, INEXT+IELMNT-1
            ISOX_FOR_OX_H(J) = J
            ISOX_TO_OX_H(J) = INEXT_OB
          END DO
          INEXT_OB = INEXT_OB + 1
 
*. A+ ia A ja
            CALL ISETVC(IHTP(1,INEXT),IZERO,4*NGAS)
            IHTP(IGAS,INEXT) = 1
            IHTP(2*NGAS+JGAS,INEXT) = 1
            IHINDX(1,INEXT) = 1
            IHINDX(2,INEXT) = 2 
            IHSGN(INEXT) = 1
            INEXT = INEXT + 1
*. a+ ib a jb
            CALL ISETVC(IHTP(1,INEXT),IZERO,4*NGAS)
            IHTP(IGAS+NGAS,INEXT) = 1
            IHTP(3*NGAS+JGAS,INEXT) = 1
            IHINDX(1,INEXT) = 1
            IHINDX(2,INEXT) = 2 
            IHSGN(INEXT) = 1
            INEXT = INEXT + 1
          END DO
        END DO
*. Two-electron operators
        DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
          DO KGAS = 1, NGAS
           DO LGAS = 1, NGAS
*. Orbital excitation 
              CALL ISETVC(IHOBTP(1,INEXT_OB),IZERO,2*NGAS)
              IHOBTP(IGAS+NGAS,INEXT_OB) = 1
              IHOBTP(KGAS+NGAS,INEXT_OB) = 
     &        IHOBTP(KGAS+NGAS,INEXT_OB) + 1
              IHOBTP(JGAS,INEXT_OB) = 1
              IHOBTP(LGAS,INEXT_OB) = 
     &        IHOBTP(LGAS,INEXT_OB) + 1
*
              IF(IGAS.GE.JGAS.AND.KGAS.GE.LGAS) THEN
                IELMNT = 3
              ELSE
                IELMNT = 1
              END IF
              NSOX_FOR_OX_H(INEXT_OB) = IELMNT
              IBSOX_FOR_OX_H(INEXT_OB) = INEXT
              DO J = INEXT, INEXT+IELMNT-1
                ISOX_FOR_OX_H(J) = J
                ISOX_TO_OX_H(J) = INEXT_OB
              END DO
              INEXT_OB = INEXT_OB + 1
*
*. a+i alpha a+j alpha ak alpha al alpha
            IF(IGAS.GE.JGAS.AND.KGAS.GE.LGAS) THEN
             CALL ISETVC(IHTP(1,INEXT),IZERO,4*NGAS)
             IHTP(IGAS,INEXT) = 1
             IHTP(JGAS,INEXT) = IHTP(JGAS,INEXT) + 1
             IHTP(2*NGAS+KGAS,INEXT) = 1
             IHTP(2*NGAS+LGAS,INEXT) = IHTP(2*NGAS+LGAS,INEXT) + 1
             IHINDX(1,INEXT) = 1
             IHINDX(2,INEXT) = 2 
             IHINDX(3,INEXT) = 3 
             IHINDX(4,INEXT) = 4 
             IHSGN(INEXT) = 1
             INEXT = INEXT + 1
*. a+i beta  a+j beta  ak beta  al beta 
             CALL ISETVC(IHTP(1,INEXT),IZERO,4*NGAS)
             IHTP(NGAS+IGAS,INEXT) = 1
             IHTP(NGAS+JGAS,INEXT) = IHTP(NGAS+JGAS,INEXT) + 1
             IHTP(3*NGAS+KGAS,INEXT) = 1
             IHTP(3*NGAS+LGAS,INEXT) = IHTP(3*NGAS+LGAS,INEXT) + 1
             IHINDX(1,INEXT) = 1
             IHINDX(2,INEXT) = 2 
             IHINDX(3,INEXT) = 3 
             IHINDX(4,INEXT) = 4 
             IHSGN(INEXT) = 1
             INEXT = INEXT + 1
            END IF
* a+i alpha a+j beta a k beta a l beta
            CALL ISETVC(IHTP(1,INEXT),IZERO,4*NGAS)
            IHTP(IGAS,INEXT) = 1
            IHTP(JGAS+NGAS,INEXT) = 1
            IHTP(KGAS+3*NGAS,INEXT) = 1
            IHTP(LGAS+2*NGAS,INEXT) = 1
             IHINDX(1,INEXT) = 1
             IHINDX(2,INEXT) = 2 
             IHINDX(3,INEXT) = 4 
             IHINDX(4,INEXT) = 3 
             IHSGN(INEXT) = -1
            INEXT = INEXT + 1
           END DO
          END DO
         END DO
        END DO
*       ^ End of loops over IGAS,JGAS,KGAS,LGAS
      END IF
*. Find pairs of operators matched by spinflips 
      IF(IFLAG.EQ.1) THEN
        DO ITP = 1, N1TP + N2TP
*. Spinflip
          CALL SPINFLIP_CAAB(IHTP(1,ITP),IH_FLIP,NGAS)
          CALL COMP_CAAB_WLIST(IH_FLIP,IHTP,NGAS,N1TP+N2TP,
     &                         ITP_SPFLIP)
           IHTP_SPINFLIP(ITP) = ITP_SPFLIP 
        END DO
* Excitation level ( times two ) for the operators 
        DO ITP = 1, N1TP + N2TP
          IH_EXC2(ITP) = IEXC2_LEVEL_FOR_CAAB(IHTP(1,ITP))
        END DO
*
      END IF
      
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*)  ' Output from H_TYPES '
        WRITE(6,*)  ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' Total number of H-types = ', N1TP + N2TP
        WRITE(6,*) 
     &  ' Number of types of one-electron operators = ', N1TP
        WRITE(6,*) 
     &  ' Number of types of two-electron operators = ', N2TP
*
        WRITE(6,*)
     &  ' Number of types of one-electron ORBITAL excitations',
     &  N1OBTP
        WRITE(6,*)
     &  ' Number of types of two-electron ORBITAL excitations',
     &  N2OBTP
        IF(IFLAG.EQ.1) THEN
          WRITE(6,*) ' CAAB occupations of the operators '
          WRITE(6,*)
          CALL WRT_SPOX_TP(IHTP,N1TP+N2TP)
*
          WRITE(6,*) ' C and A of orbital excitations '
          DO IOBTP = 1, N1OBTP + N2OBTP
            WRITE(6,'(3X,I3,20I3)')
     &      IOBTP, (IHOBTP(I,IOBTP),I=1, 2*NGAS)  
          END DO
*
          WRITE(6,*) ' Orbital => spin-orbital translations '
          WRITE(6,*) ' ===================================== '
          DO IOBTP = 1, N1OBTP + N2OBTP
            WRITE(6,*) ' Number of SOXs for OX ', IOBTP ,'is ', 
     &      NSOX_FOR_OX_H(IOBTP)
            WRITE(6,*) ' and the SOXs : '
            IB = IBSOX_FOR_OX_H(IOBTP)
            L  = NSOX_FOR_OX_H(IOBTP)
            CALL IWRTMA(ISOX_FOR_OX_H(IB),1,L,1,L) 
          END DO
*
          WRITE(6,*) ' Operator => Spinflipped operator '
          CALL IWRTMA(IHTP_SPINFLIP,1,N1TP+N2TP,1,N1TP+N2TP)
*
          WRITE(6,*) ' Excitation level times two for '
          CALL IWRTMA(IH_EXC2,1,N1TP+N2TP,1,N1TP+N2TP)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_EXPTOP(NSPOBEX_TP,LSPOBEX_TP,N_CC_AMP,E,F,IB_F)
*
* Write the part of Exp(T) needed for the CC vector function
*
* Jeppe Olsen, May 2000
*              IB_F added November 2000
*
      INCLUDE 'implicit.inc'
      INTEGER IB_F(NSPOBEX_TP,NSPOBEX_TP),LSPOBEX_TP(NSPOBEX_TP)
*
      DIMENSION E(N_CC_AMP),F(*)
*
      WRITE(6,*)
      WRITE(6,*) ' Exp T is written as 1 + E + sum(I) T(I) F(I) '
      WRITE(6,*) ' ============================================='
      WRITE(6,*)
      WRITE(6,*) ' E-operator : '
      WRITE(6,*)
      CALL WRTBLKN(E,NSPOBEX_TP,LSPOBEX_TP)
      DO I = 1, NSPOBEX_TP
       WRITE(6,*)  
       WRITE(6,*) 'F-operator for spinorbital excitation type ', I
       DO J = 1, NSPOBEX_TP  
         IF(IB_F(J,I).GT.0) THEN 
           WRITE(6,*) ' Block = ', J
           IB = IB_F(J,I)
           LEN = LSPOBEX_TP(J)
           CALL WRTMAT(F(IB),1,LEN,1,LEN)
         END IF
       END DO
      END DO
*
      RETURN 
      END
      SUBROUTINE STRING_IN_CCEXCIT(ISTRING,ICAAB,INSPC,
     &                             NSPOBEX_TP,ISPOBEX_TP)
*
* Is String ISTRING ( given as occ in each gas space)
* include as CAAB string ICAAB in any of the NSPOBEX_TP TCC excitations 
* blocks given in ISPOBEX_TP 
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER ISTRING(NGAS)
*
      INSPC = 0
      IB = 1 + (ICAAB-1)*NGAS
      DO KTP = 1, NSPOBEX_TP
        IDEL = 0
        DO IGAS = 1, NGAS
          IDEL = IDEL + ABS(ISPOBEX_TP(IB-1+IGAS,KTP)-ISTRING(IGAS))
        END DO
        IF(IDEL.EQ.0) THEN
          INSPC = 1
          GOTO 1001 
        END IF
      END DO
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ICAAB of string = ', ICAAB
        WRITE(6,*) ' Occupation : '
        CALL IWRTMA(ISTRING,1,NGAS,1,NGAS)
        IF(INSPC.EQ.1) THEN
           WRITE(6,*) ' String is in space '
        ELSE
           WRITE(6,*) ' String is out of space '
        END IF
      END IF
*
      RETURN 
      END
      SUBROUTINE DIM_T1T2_TO_T12_MAP(LEN_T1T2_STRING,LEN_T1T2_TCC)
*
* Dimension of T1*T2 => T12 mappings of strings
*
* LEN_T1T2_STRING : Max dimension of MAT(T1,T2) where T1 and 
*                   T2 are spin strings, and T1*T2 is allowed 
*                   spin strings
*
* Jeppe Olsen, May 2000 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cgas.inc'
*
      CALL DIM_T1T2_TO_T12_S(LEN_T1T2_STRING,LEN_T1T2_TCC,
     &                        NSPOBEX_TP,WORK(KLSOBEX),
     &                        NGAS)
      RETURN
      END
      SUBROUTINE DIM_T1T2_TO_T12_S(LEN_T1T2_STRING,LEN_T1T2_TCC,
     &           NSPOBEX_TP,ISPOBEX_TP,NGAS)
*
* LEN_T1T2_STRING : Max dimension of MAT(T1,T2) where T1 and 
*                   T2 are spin strings, and T1*T2 is an allowed 
*                   spin strings. T1 and T2 are assumed to 
*                   run over all symmetries 
*                   
*. Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'csm.inc'
*. Input 
      INTEGER ISPOBEX_TP(4*NGAS,*)
*. Local scratch
      INTEGER I12OCC(MXPNGAS)
      INTEGER NI(8)
      INTEGER NJ(8)
* 
      LEN_T1T2_TCC = 0
*     ^ Not implemented yet
*
      LEN_T1T2_STRING = 0
      IONE = 1
      DO ICAAB = 1, 4
       DO ITP = 1, NSPOBEX_TP
        DO JTP = 1, NSPOBEX_TP
C?        WRITE(6,*) ' ITP, JTP, ICAAB = ', ITP, JTP, ICAAB
*. ITP*JTP, component ICAAB
          IB = 1 + (ICAAB-1)*NGAS
          CALL IVCSUM(I12OCC,ISPOBEX_TP(IB,ITP),ISPOBEX_TP(IB,JTP),
     &                IONE,IONE,NGAS)
*. Is I12OCC included in any of the strings
C              STRING_IN_CCEXCIT(ISTRING,ICAAB,INSPC,
C    &                           NSPOBEX_TP,ISPOBEX_TP)
*
          CALL STRING_IN_CCEXCIT(I12OCC,ICAAB,INSPC,NSPOBEX_TP,
     &                           ISPOBEX_TP) 
*
          IF(INSPC.EQ.1) THEN 
*. String is in space, find Dimensions of ITP and JTP per sym
            CALL NST_SPGP(ISPOBEX_TP(IB,ITP),NI)
            CALL NST_SPGP(ISPOBEX_TP(IB,JTP),NJ)
            NI_TOT = IELSUM(NI,NSMST)
            NJ_TOT = IELSUM(NJ,NSMST)
            LEN = NI_TOT*NJ_TOT
            LEN_T1T2_STRING = MAX(LEN_T1T2_STRING,LEN)
          END IF
        END DO
       END DO
*      ^ End of loop over ITP, JTP
      END DO
*     ^ End of loop over ICAAB
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' LEN_T1T2_STRING = ', LEN_T1T2_STRING
      END IF
*
      RETURN
      END
      SUBROUTINE PROD_SPOB_EX_TP(ITP,JTP,IJTP,IEXCLEVEL,I_AM_ZERO)
*
* Two spin-orbital excitation types ITP and JTP are 
* given. Find product excitation type IJTP
*
* IF IJTP = 0, then the product is not inside the 
* excitation space and IEXCLEVEL gives the 
* smallest rank of an operator required to 
* bring ITP and JTP into the space
*
* Jeppe Olsen, May 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
*
      CALL PROD_SPOB_EX_TP1(ITP,JTP,IJTP,IEXCLEVEL,WORK(KLSOBEX),
     &                      NSPOBEX_TP,I_AM_ZERO)
*
      RETURN
      END
      SUBROUTINE PROD_SPOB_EX_TP1(ITP,JTP,IJTP,IEXCLEVEL_MIN,
     &                            ISPOBEX_TP,NSPOBEX_TP,I_AM_ZERO)
*
* Product of spinorbital types ITP and JTP
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,*)
*. Local scratch 
      INTEGER IJSPOBEX_TP(4*MXPNGAS)
*. Tired of GNU warnings ...
      IADEL_MIN = -1
      ICDEL_MIN = -1    
      IEXCLEVEL_MIN =  0 
*
      IJTP = 0
*. Occupation of ITP JTP
CERRORIOFF = (ITP-1)*4*NGAS+1
CERRORJOFF = (JTP-1)*4*NGAS+1
      IONE = 1
      CALL IVCSUM(IJSPOBEX_TP,ISPOBEX_TP(1,ITP),ISPOBEX_TP(1,JTP),
     &            IONE,IONE,4*NGAS)
*. Does this match any allowed orbital type, and 
*. if not what is lowest rank of excitation needed for connection 
*. Assumning excitation rank less than 10000
      DO KTP = 1, NSPOBEX_TP
*. Creation 
       ICDEL = 0
       DO IGASAB = 1, 2*NGAS
         ICDEL = ICDEL + 
     &   ABS(ISPOBEX_TP(IGASAB,KTP)-IJSPOBEX_TP(IGASAB))
       END DO
       IF(KTP.EQ.1) THEN
         ICDEL_MIN = ICDEL
       ELSE 
         ICDEL_MIN = MIN(ICDEL_MIN,ICDEL)
       END IF
*. Annihilation 
       IADEL = 0
       DO IGASAB = 2*NGAS+1, 4*NGAS
         IADEL = IADEL + 
     &   ABS(ISPOBEX_TP(IGASAB,KTP)-IJSPOBEX_TP(IGASAB))
       END DO
       IF(KTP.EQ.1) THEN
         IADEL_MIN = IADEL
       ELSE 
         IADEL_MIN = MIN(IADEL_MIN,IADEL)
       END IF
*. A hit 
       IF(IADEL.EQ.0.AND.ICDEL.EQ.0) THEN
         IJTP = KTP
       ELSE
         IEXCLEVEL = MAX(ICDEL,IADEL)
         IF(KTP.EQ.1) THEN
           IEXCLEVEL_MIN = IEXCLEVEL
         ELSE
           IEXCLEVEL_MIN = MIN(IEXCLEVEL,IEXCLEVEL_MIN)
         END IF
       END IF
      END DO
*     ^ End of loop over KTP
*. Is the excitation vanishing : More operators than orbitals in any space ?
      I_AM_ZERO = 0
      DO ICAAB = 1, 4
        DO IGAS = 1, NGAS
          IB = (ICAAB-1)*NGAS+1
          IF(IJSPOBEX_TP(IB-1+IGAS).GT.NOBPT(IGAS)) I_AM_ZERO = 1
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)' Product of spin-orbital excitation types ',ITP,JTP
        WRITE(6,*)
        IF(IJTP.NE.0) THEN
          WRITE(6,*) ' In space, IJTP = ', IJTP
        ELSE 
          IF(I_AM_ZERO.EQ.0) THEN 
            WRITE(6,*) 
     &      ' Out of space, excitation level for connection =',
     &      IEXCLEVEL_MIN
          ELSE
            WRITE(6,*) ' Vanishing CC excitation operator '
          END IF
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE H_EXP_T_MEM(I_NEW_OR_OLD)
*
* set up scratch space for direct calculation of H EXP T
*
* Pointers to scratch are returned in CC_SCR
*
*. Jeppe Olsen, May 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
*
      COMMON/CC_SCR/ KCCF,KCCE,KCCVEC1,KCCVEC2,
     &              KISTST1,KXSTST1,KISTST2,KXSTST2,
     &              KISTST3,KXSTST3,KISTST4,KXSTST4,
     &              KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &              KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTREO,
     &              KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &              KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &              KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &              KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &              KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &              KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &              KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &              KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &              KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
     &              KLOPSCR,
     &              KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &              KLSMD1,KLSMD2,KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &              KLNMD1,KLNMD2,KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &              KLOCK1, KLOCK2, KLOCL1, KLOCOT1T2, KL_IBF, 
     &              KLEXEORD

*
*.1 : KCCF : Intermediate F-vectors for needed out of space 
*         excitations in the form sum_I T(I) F(I) 
*
*. Largest length of F-vector and T- block for given type of T-block
*. ( It is convenient to store these together )
C            LEN_FOR_TF(NSPOBEX,ISPOBEX,LSPOBEX,NGAS,LEN_TF)
        CALL LEN_FOR_TF(NSPOBEX_TP,WORK(KLSOBEX),WORK(KLLSOBEX),
     &                  NGAS,LEN_TF)
C       WRITE(6,*) ' ... MEM, LEN_TF (1) = ', LEN_TF
*. It is convenient to be able to use F as a normal CC also so
        N_CC_AMPP1 = N_CC_AMP + 1
        LEN_TF = MAX(LEN_TF,N_CC_AMPP1)
        CALL MEMMAN(KCCF,LEN_TF,'ADDL  ',2,'CCF   ')
*
* For a given orbital type I the start of the corresponding 
* intermediate array is given as KCCF    + (I-1)*N_CC_AMP
*
*.2 : KCCVEC1,KCCVEC2,KCCE    : Three  CC vectors, may allowed to hold 
*.    component corresponding to unit operator 
      CALL MEMMAN(KCCE   ,N_CC_AMPP1,'ADDL  ',2,'CCE   ')
C     CALL MEMMAN(KCCVEC1,N_CC_AMPP1,'ADDL  ',2,'CCVEC1')
      KCCVEC1 = 0
      KCCVEC2 = 0
*.3 : Maps for T1T2 => T12 for individual strings types 
      CALL DIM_T1T2_TO_T12_MAP(LEN_T1T2_STRING,LENT_T1T2_TCC)
      LEN = LEN_T1T2_STRING
*. Arrays for holding elementary operators times strings
      CALL MEMMAN(KISTST1,LEN,'ADDL  ',1,'ISTST1')
      CALL MEMMAN(KXSTST1,LEN,'ADDL  ',2,'XSTST1')
*
      CALL MEMMAN(KISTST2,LEN,'ADDL  ',1,'ISTST2')
      CALL MEMMAN(KXSTST2,LEN,'ADDL  ',2,'XSTST2')
*
      CALL MEMMAN(KISTST3,LEN,'ADDL  ',1,'ISTST3')
      CALL MEMMAN(KXSTST3,LEN,'ADDL  ',2,'XSTST3')
*
      CALL MEMMAN(KISTST4,LEN,'ADDL  ',1,'ISTST4')
      CALL MEMMAN(KXSTST4,LEN,'ADDL  ',2,'XSTST4')
*
      CALL MEMMAN(KIX1_CA,LEN,'ADDL  ',1,'IX1_CA')
      CALL MEMMAN(KSX1_CA,LEN,'ADDL  ',2,'SX1_CA')
*
      CALL MEMMAN(KIX1_CB,LEN,'ADDL  ',1,'IX1_CB')
      CALL MEMMAN(KSX1_CB,LEN,'ADDL  ',2,'SX1_CB')
*
      CALL MEMMAN(KIX1_AA,LEN,'ADDL  ',1,'IX1_AA')
      CALL MEMMAN(KSX1_AA,LEN,'ADDL  ',2,'SX1_AA')
*
      CALL MEMMAN(KIX1_AB,LEN,'ADDL  ',1,'IX1_AB')
      CALL MEMMAN(KSX1_AB,LEN,'ADDL  ',2,'SX1_AB')
*
      CALL MEMMAN(KIX2_CA,LEN,'ADDL  ',1,'IX2_CA')
      CALL MEMMAN(KSX2_CA,LEN,'ADDL  ',2,'SX2_CA')
*
      CALL MEMMAN(KIX2_CB,LEN,'ADDL  ',1,'IX2_CB')
      CALL MEMMAN(KSX2_CB,LEN,'ADDL  ',2,'SX2_CB')
*
      CALL MEMMAN(KIX2_AA,LEN,'ADDL  ',1,'IX2_AA')
      CALL MEMMAN(KSX2_AA,LEN,'ADDL  ',2,'SX2_AA')
*
      CALL MEMMAN(KIX2_AB,LEN,'ADDL  ',1,'IX2_AB')
      CALL MEMMAN(KSX2_AB,LEN,'ADDL  ',2,'SX2_AB')
*
      CALL MEMMAN(KIX3_CA,LEN,'ADDL  ',1,'IX3_CA')
      CALL MEMMAN(KSX3_CA,LEN,'ADDL  ',2,'SX3_CA')
*
      CALL MEMMAN(KIX3_CB,LEN,'ADDL  ',1,'IX3_CB')
      CALL MEMMAN(KSX3_CB,LEN,'ADDL  ',2,'SX3_CB')
*
      CALL MEMMAN(KIX3_AA,LEN,'ADDL  ',1,'IX3_AA')
      CALL MEMMAN(KSX3_AA,LEN,'ADDL  ',2,'SX3_AA')
*
      CALL MEMMAN(KIX3_AB,LEN,'ADDL  ',1,'IX3_AB')
      CALL MEMMAN(KSX3_AB,LEN,'ADDL  ',2,'SX3_AB')
*
      CALL MEMMAN(KIX4_CA,LEN,'ADDL  ',1,'IX4_CA')
      CALL MEMMAN(KSX4_CA,LEN,'ADDL  ',2,'SX4_CA')
*
      CALL MEMMAN(KIX4_CB,LEN,'ADDL  ',1,'IX4_CB')
      CALL MEMMAN(KSX4_CB,LEN,'ADDL  ',2,'SX4_CB')
*
      CALL MEMMAN(KIX4_AA,LEN,'ADDL  ',1,'IX4_AA')
      CALL MEMMAN(KSX4_AA,LEN,'ADDL  ',2,'SX4_AA')
*
      CALL MEMMAN(KIX4_AB,LEN,'ADDL  ',1,'IX4_AB')
      CALL MEMMAN(KSX4_AB,LEN,'ADDL  ',2,'SX4_AB')
*
*.4 : KLZ,KLZSCR : Memory for a Z matrix and scratch for 
*     constructing Z
*     
      IATP = 1
      IBTP = 2
      NAEL = NELFTP(IATP)
      NBEL = NELFTP(IBTP)
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      CALL MEMMAN(KLZ,LZ,'ADDL  ',1,'Z     ')
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'ZSCR  ')
*. String occupations for given CAAB, all symmetris
      LEN = MX_ST_TSOSO_BLK_MX*NSMST*4
      CALL MEMMAN(KLSTOCC1,LEN,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTOCC2,LEN,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTOCC3,LEN,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTOCC4,LEN,'ADDL  ',1,'STOCC4')
      CALL MEMMAN(KLSTOCC5,LEN,'ADDL  ',1,'STOCC5')
      CALL MEMMAN(KLSTOCC6,LEN,'ADDL  ',1,'STOCC6')

*. Reorder array for given CAAB, all symmetries 
      CALL MEMMAN(KLSTREO,MX_ST_TSOSO_MX*NSMST,'ADDL  ',1,'STREO ')
*. Intermediate blocks, (LD12B,LD12B)
*. For the moment 
      LD12B = LCCBD12
      LB = LCCB
      LEN = LD12B**2
      LEN3 = LD12B**3
      CALL MEMMAN(KLTSCR1,LEN,'ADDL  ',2,'TSCR1 ')
      CALL MEMMAN(KLTSCR2,LEN3,'ADDL  ',2,'TSCR2 ')
      CALL MEMMAN(KLTSCR3,LEN,'ADDL  ',2,'TSCR3 ')
*
C          K_RES_DIM(ISPOBEX,NSPOBEX,MAXOP,NELMNT_MAX)
      CALL K_RES_DIM(WORK(KLSOBEX),NSPOBEX_TP,4,MX_KBLK)    
COLD  CALL T_DL_DIM(WORK(KLSOBEX),NSPOBEX_TP,4,N_TDL_MAX,LD12B)
C     WRITE(6,*) ' MX_KBLK in H_EXP ... = ', MX_KBLK
*. A big one for storing T(D2,L)
      LEN = N_TDL_MAX
      CALL MEMMAN(KLTSCR4,LEN,'ADDL  ',2,'TSCR4 ')
*. For batches of indeces of Hamiltionian operator
      LEN = 4*LD12B
      CALL MEMMAN(KLIOD1_ST,LEN,'ADDL  ',2,'IOD1_S')
      CALL MEMMAN(KLIOD2_ST,LEN,'ADDL  ',2,'IOD2_S')
      CALL MEMMAN(KLIOEX_ST,LEN,'ADDL  ',2,'IOEX_S')
*and for a batch of Hamiltonian
      LOPSCR = LD12B*LD12B*LB
      CALL MEMMAN(KLOPSCR,LOPSCR,'ADDL  ',2,'OPSCR ')
*. Number and symmetries of each substring for 6 complete T-blocks 
*. Largest number of strings in intermediate arrays
*
*. D1, D2, EX holds part of Hamiltonian so
*. MXTSOB_P is largest number of particle orbitals of given type, ALL SYM!
C     LEN = 4*MXTSOB_P**2 * MXTSOB_H**2
      LEN = 4*LB
      CALL MEMMAN(KLSMD1,LEN,'ADDL  ',1,'SM_D1 ')
      CALL MEMMAN(KLNMD1,LEN,'ADDL  ',1,'NM_D1 ')
*
C?    CALL ISETVC(WORK(KLSMD1),-42,LEN)
C?    CALL ISETVC(WORK(KLNMD1),-42,LEN)
*
C     LEN = 4*MXTSOB_P**2 * MXTSOB_H**2
      CALL MEMMAN(KLSMD2,LEN,'ADDL  ',1,'SM_D2 ')
      CALL MEMMAN(KLNMD2,LEN,'ADDL  ',1,'NM_D2 ') 
*
      CALL MEMMAN(KLSMEX,LEN,'ADDL  ',1,'SM_EX ')
      CALL MEMMAN(KLNMEX,LEN,'ADDL  ',1,'NM_EX ') 
*. K1, K2, L1 holds general strings
C     LEN = 4*MX_KBLK
      LEN = 4*LB
      CALL MEMMAN(KLSMK1,LEN,'ADDL  ',1,'SM_K1 ')
      CALL MEMMAN(KLNMK1,LEN,'ADDL  ',1,'NM_K1 ') 
*
*
      CALL MEMMAN(KLSMK2,LEN,'ADDL  ',1,'SM_K2 ')
      CALL MEMMAN(KLNMK2,LEN,'ADDL  ',1,'NM_K2 ') 
*
      CALL MEMMAN(KLSML1,LEN,'ADDL  ',1,'SM_L1 ')
      CALL MEMMAN(KLNML1,LEN,'ADDL  ',1,'NM_L1 ') 
*     KLOCK1, KLOCK2, KLOCL1
      LEN = 4*NGAS
      CALL MEMMAN(KLOCK1,LEN,'ADDL  ',1,'IOC_K1')
      CALL MEMMAN(KLOCK2,LEN,'ADDL  ',1,'IOC_K2')
      CALL MEMMAN(KLOCL1,LEN,'ADDL  ',1,'IOC_L1')
      CALL MEMMAN(KLOCOT1T2,LEN,'ADDL  ',1,'IOC_OT')
*. Array for offsets in F
      CALL MEMMAN(KL_IBF,NSPOBEX_TP,'ADDL  ',1,'IB_F  ')
*. Array for execution order of spinorbital types in Exp T
      CALL MEMMAN(KLEXEORD,NSPOBEX_TP+1,'ADDL  ',1,'EXEORD')
*
      RETURN
      END
      SUBROUTINE EXP_MT_H_EXP_T(T,HEXPT,CCVEC2,E_CC2)
*
* Calculate the exp-T H exp T in CC subspace
* All deexcitation terms in Hamiltonian are contracted with T-operators
*
* Jeppe Olsen, Another weekend at Chemistry, Aarhus, May 2000
*              June 2001 : Low memory route added 
*              July 2001 : Some spin-flipping added
*
*
* Exp T contains three parts
*
*    1 : Part that is inside CC excitation space
*    2 : Part that is outside CC excitation space but can be 
*        brought into CC space by deexcitations in H
*    3 : Part that is outside CC excitation space and cannot 
*        be brought into CC space by deexcitations in H
*. Part 1 may be constructed simple and explicit, part 3 does not 
*. contribute so the problematic gay is part 2.
*. We will construct part 2 in the form 
*  sum_(I) T(I) OP(I)
* where T(I) is a TCC block and OP(I) is inside the CC space.
* To write the part 2 operator in 
* this form, we assume that single exitations are treated 
* as a Hamiltonian transformation. If single excitations where 
* present, there would be terms like sum(IJ) S(I) S(J) OP(IJ).
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*.Specific input
      DIMENSION T(*)
*. output
      DIMENSION HEXPT(*)
*. General input
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Scratch through argument list 
      DIMENSION CCVEC2(*)
*. Scr for CC calculations
      COMMON/CC_SCR/ KCCF,KCCE,KCCVEC1,KCCVEC2,
     &              KISTST1,KXSTST1,KISTST2,KXSTST2,
     &              KISTST3,KXSTST3,KISTST4,KXSTST4,
     &              KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &              KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTREO,
     &              KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &              KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &              KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &              KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &              KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &              KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &              KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &              KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &              KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
     &              KLOPSCR,
     &              KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &              KLSMD1,KLSMD2,KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &              KLNMD1,KLNMD2,KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &              KLOCK1, KLOCK2, KLOCL1, KLOCOT1T2, KL_IBF, 
     &              KLEXEORD
*
      COMMON/ROLLO/NINNER,NINNER1,NINNER2,NINNER3,NINNER4,NINNER5,
     &             NINNER6,NINNER7,NINNER8,NINNER9,NINNER10
      COMMON/CNTDL_MAX/NTDL_MAX_ACT,ID2_MAX_ACT(MXPNGAS*4),
     &      IL1_MAX_ACT(MXPNGAS*4),IEX_MAX_ACT(MXPNGAS*4),
     &      IT1_MAX_ACT(MXPNGAS*4),IT2_MAX_ACT(MXPNGAS*4)
*. A bit of Local scratch
*. Local scratch
      INTEGER ISX(2*MXPNGAS*MXPNGAS)

*
      IDUM = 0
      CALL LUCIAQENTER('HEXPT ')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HEXPT ')
      LCCBD12 = LCCB
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' EXP_MT_H_EXP speaking '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,*) ' Input T vector '
        CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
      END IF
C     WRITE(6,*) ' IUSE_TR = ', IUSE_TR
*. For stats for TOKB
      NINNER  = 0
      NINNER1 = 0
      NINNER2 = 0
      NINNER3 = 0
      NINNER4 = 0
      NINNER5 = 0
      NINNER6 = 0
      NINNER7 = 0
      NINNER8 = 0
      NINNER9 = 0
      NINNER10 = 0
      NTDL_MAX_ACT = 0
*. Memory
      CALL H_EXP_T_MEM(I_NEW_OR_OLD)
*
*. Calculate H exp T |HF> as H(1 + sum_(I) E(I) + sum_I T(I) F(I)) 
*. Order of execution for spinorbital types in Exp T
      CALL ORD_CCEXTP(WORK(KLEXEORD))
*. Identify Single Excitations - to be used for alternative E_CC calculation
      CALL GET_SPOBTP_FOR_EXC_LEVEL(1,WORK(KLCOBEX_TP),NSPOBEX_TP,
     &     NSXTP,ISX,WORK(KLSOX_TO_OX))
*
* Save T vector in blocks on file 63  - use reordered vector
      CALL MEMMAN(KLLREO,NSPOBEX_TP+1,'ADDL  ',2,'LREO  ')
      CALL REO_BLK_VEC(T,WORK(KLLSOBEX),NSPOBEX_TP,WORK(KCCE),
     &                 WORK(KLLREO),WORK(KLEXEORD),1)
*. and add the unitoperator
      WORK(KCCE-1+N_CC_AMP+1) = T(N_CC_AMP+1)
C          ICOPVE3(IIN,IOFFIN,IOUT,IOFFOUT,NDIM)
      CALL ICOPVE3(1,1,WORK(KLLREO),NSPOBEX_TP+1,1)
      CALL ICOPVE3(NSPOBEX_TP+1,1,WORK(KLEXEORD),NSPOBEX_TP+1,1)
* 
      LU63 = 63
      CALL REWINO(LU63)
      CALL TODSCN(WORK(KCCE),NSPOBEX_TP+1,WORK(KLLREO),-1,LU63)
      CALL ITODS(-1,1,-1,LU63)
C?    WRITE(6,*) ' Vector read in from LU63'
C?    CALL WRTVCD(WORK(KCCE),LU63,1,-1)

*
* New version where Each batch of F is contracted immediately with H
*
      ZERO = 0.0D0
      N_CC_AMPP1 = N_CC_AMP + 1
      CALL SETVEC(HEXPT,ZERO,N_CC_AMPP1)
* H Exp T
* T vector is also used as CCVEC1
      I_DO_FTERMS = 1
      I_DO_ZERO = 1
      CALL HEXP_T(T,WORK(KCCE),WORK(KCCF),
     &     NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,WORK(KLLSOBEX),
     &     WORK(KLIBSOBEX),WORK(KLSOX_TO_OX),
     &     T,CCVEC2,
     &     WORK(KISTST1),WORK(KXSTST1),WORK(KISTST2),WORK(KXSTST2),
     &     WORK(KISTST3),WORK(KXSTST3),WORK(KISTST4),WORK(KXSTST4),
     &     WORK(KLZ),WORK(KLZSCR),WORK(KLSTOCC1),WORK(KLSTOCC2),
     &     WORK(KLSTREO),I_DO_ZERO,I_DO_FTERMS,WORK(KLSPOBEX_AC),
     &     WORK(KLEXEORD),HEXPT,WORK(KLSPOBEX_FRZ),
     &     WORK(KLSPOBEX_AB),WORK(KL_IBF),LU63)
      XNORM_HEXPT = INPROD(HEXPT,HEXPT,N_CC_AMP+1)
      WRITE(6,'(A,E22.13)') ' Norm of H exp(t) | hf > ', XNORM_HEXPT
      XNORM_EXPT = INPROD(WORK(KCCE),WORK(KCCE),N_CC_AMP+1)
      WRITE(6,'(A,E22.13)') ' Norm of exp(t) | hf > ', XNORM_EXPT
      IF(IUSE_TR.EQ.1) THEN
        CALL SPINFLIP_CC_BLOCKS(HEXPT,NSPOBEX_TP+1,WORK(KLSPOBEX_AB),
     &       WORK(KLIBSOBEX),WORK(KLSOBEX),1,NSMST,NGAS) 
      END IF
*
*. Alternative evaluation of energy = <E|H|CC > / <E|CC>
*. There is a bit of problem with single excitations.
*. What we have in HEXPT is Exp(-T1) H Exp(T) !HF>
*. and what we have in E is Exp(-T1) Exp(T) !HF> ( E without single excitations)
*. Calculate Exp(T1) !E> and Exp(T1) HEXPT
*. Change so ISPOBEX_AC so only single excitations are active
      I_DO_ALT = 0
      IF(I_DO_ALT.EQ.1) THEN
        IZERO = 0 
        CALL ISETVC(WORK(KLSPOBEX_AC),IZERO,NSPOBEX_TP)
        IONE = 1
        CALL ISCASET(WORK(KLSPOBEX_AC),IONE,ISX,NSXTP)
*. Retrieve T
        CALL REWINO(LU63)
        CALL FRMDSCN(CCVEC2,NSPOBEX_TP+1,-1,LU63)
*. Obtain in original order 
        CALL REO_BLK_VEC(T,WORK(KLLSOBEX),NSPOBEX_TP+1,CCVEC2,
     &                   WORK(KLLREO),WORK(KLEXEORD),-1)
        I_DO_ZERO = 0
        I_DO_FTERMS = 0
        XDUM = 0.0D0
*. Exp(T1) E
        CALL EXP_T(T,WORK(KCCE),XDUM,
     &       NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,WORK(KLLSOBEX),
     &       WORK(KLIBSOBEX),WORK(KLSOX_TO_OX),
     &       WORK(KCCF),CCVEC2,
     &       WORK(KISTST1),WORK(KXSTST1),WORK(KISTST2),WORK(KXSTST2),
     &       WORK(KISTST3),WORK(KXSTST3),WORK(KISTST4),WORK(KXSTST4),
     &       WORK(KLZ),WORK(KLZSCR),WORK(KLSTOCC1),WORK(KLSTOCC2),
     &       WORK(KLSTREO),I_DO_ZERO,I_DO_FTERMS,WORK(KLSPOBEX_AC),
     &       WORK(KL_IBF),LEN_F,WORK(KLEXEORD))
*. Exp(T1) HEXPT
        CALL EXP_T(T,HEXPT,XDUM,
     &       NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,WORK(KLLSOBEX),
     &       WORK(KLIBSOBEX),WORK(KLSOX_TO_OX),
     &       WORK(KCCF),CCVEC2,
     &       WORK(KISTST1),WORK(KXSTST1),WORK(KISTST2),WORK(KXSTST2),
     &       WORK(KISTST3),WORK(KXSTST3),WORK(KISTST4),WORK(KXSTST4),
     &       WORK(KLZ),WORK(KLZSCR),WORK(KLSTOCC1),WORK(KLSTOCC2),
     &       WORK(KLSTREO),I_DO_ZERO,I_DO_FTERMS,WORK(KLSPOBEX_AC),
     &       WORK(KL_IBF),LEN_F,WORK(KLEXEORD))
        E_CC2 = INPROD(WORK(KCCE),HEXPT,N_CC_AMP+1)/
     &          INPROD(WORK(KCCE),WORK(KCCE),N_CC_AMP+1)
        WRITE(6,'(A,F22.15)') 
     &  ' Alternative evaluation (2) of E_CC = ', E_CC2
         END IF
         IF(NTEST.GE.500) THEN
           WRITE(6,*) ' HEXPT : '
           CALL WRTMAT(HEXPT,1,N_CC_AMP+1,1,N_CC_AMP+1)
         END IF
*. Exp(-T) H Exp(T)
*. As we in the alternative E_CC evaluation multiplied 
*. Exp(-T1) HexpT with Exp(T1) then all T operators should be included 
*. in the following so
        IF(I_DO_ALT.EQ.1) THEN
          IONE = 1 
          CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TP)
        END IF
      CALL REWINO(LU63)
      CALL FRMDSCN(CCVEC2,NSPOBEX_TP+1,-1,LU63)
*. Obtain in original order 
      CALL REO_BLK_VEC(T,WORK(KLLSOBEX),NSPOBEX_TP+1,CCVEC2,
     &                 WORK(KLLREO),WORK(KLEXEORD),-1)
        ONEM = -1.0D0
        CALL SCALVE(T,ONEM,N_CC_AMP)
*. No F terms and no zeroing of E
        I_DO_ZERO = 0
        I_DO_FTERMS = 0
        CALL EXP_T(T,HEXPT,WORK(KCCF),      
     &       NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,WORK(KLLSOBEX),
     &       WORK(KLIBSOBEX),WORK(KLSOX_TO_OX),
     &       WORK(KCCE),CCVEC2,
     &       WORK(KISTST1),WORK(KXSTST1),WORK(KISTST2),WORK(KXSTST2),
     &       WORK(KISTST3),WORK(KXSTST3),WORK(KISTST4),WORK(KXSTST4),
     &       WORK(KLZ),WORK(KLZSCR),WORK(KLSTOCC1),WORK(KLSTOCC2),
     &       WORK(KLSTREO),I_DO_ZERO,I_DO_FTERMS,WORK(KLSPOBEX_AC),
     &       WORK(KL_IBF),LEN_F,WORK(KLEXEORD))
        CALL SCALVE(T,ONEM,N_CC_AMP)
*. For stability, ensure final vector is spinsymmetrized
      IF(IUSE_TR.EQ.1) THEN
        CALL SPINFLIP_CC_BLOCKS(HEXPT,NSPOBEX_TP+1,WORK(KLSPOBEX_AB),
     &       WORK(KLIBSOBEX),WORK(KLSOBEX),1,NSMST,NGAS) 
      END IF
*and we are finished !!
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from EXP_MT_H_EXP_T : '     
        CALL WRTBLKN(HEXPT,NSPOBEX_TP+1,WORK(KLLSOBEX))
      END IF
*. Stats for TOK
CT    WRITE(6,*) ' NINNER NINNER1 NINNER2 =', NINNER,NINNER1,NINNER2
CT    WRITE(6,*) ' NINNER3 NINNER4 NINNER5 =',NINNER3,NINNER4,NINNER5
CT    WRITE(6,*) ' NINNER6 NINNER7 =',NINNER6,NINNER7
CT    WRITE(6,*) ' NINNER8 NINNER9 =',NINNER8,NINNER9
CT    WRITE(6,*) ' NINNER10 =',NINNER10
*
C     WRITE(6,*) ' NTDL_MAX_ACT = ', NTDL_MAX_ACT
CT    WRITE(6,*) ' Corresponding L1 D2 Ex : '
CT    CALL WRT_SPOX_TP(IL1_MAX_ACT,1)
CT    CALL WRT_SPOX_TP(ID2_MAX_ACT,1)
CT    CALL WRT_SPOX_TP(IEX_MAX_ACT,1)
CT    WRITE(6,*) ' Corresponding T1 T2  : '
CT    CALL WRT_SPOX_TP(IT1_MAX_ACT,1)
CT    CALL WRT_SPOX_TP(IT2_MAX_ACT,1)
*
COLD  CALL TEST_D1
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HEXPT ')
      CALL LUCIAQEXIT('HEXPT ')
C!    STOP ' Enforced stop at end of EXP_MT ....'
      RETURN
      END
      SUBROUTINE EXP_T(T,E,F,NOBEX_TP,ISOBEX_TP,NSOBEX_TP,
     &           LSOBEX_TP,IBSOBEX_TP,ISOX_TO_OX,
     &           CCVEC1,CCVEC2,
     &           ISTST1,XSTST1,ISTST2,XSTST2,
     &           ISTST3,XSTST3,ISTST4,XSTST4,
     &           IZ,IZSCR,ISTOCC1,ISTOCC2,ISTREO,
     &           I_DO_ZERO,I_DO_FTERMS,ISPOBEX_AC,IB_F,LEN_F,
     &           IEXEORD)
*
* routine for calculation of the part of exp(T) relevant for the 
* evaluation of the CC vector function
*
* Exp(T) is produced as
*
* E0 + sum(I) E(I) + sum(I) T(I) F(I), 
* where sum is over spin orbital types , and E0 is the initial value 
* in E(N_CC_AMP+1)
*
*. About the algorithm
* =====================
*
* Contributions from each spin-orbital excitation operator T(I) are
* evaluated separately. 
*
*. After K terms we have the form 
*  E0 + sum_I E(K,I) + sum_I T(I) F(K,I)
* and assume that E(K,I) is in E
*. and we must now calculate in CCVEC3( setting T' = T(K+1))
* (1 + T' + 1/2 T'^2 + ... 1/l! T'^l) (1 + sum_I E(K,I) + sum_I T(I) F(K,I))
*
*  N = 0 term : Current E-term
*  (N = 1)  
*           Copy E to CCVEC1
*           T' + T' CCVEC1 in CCVEC2
*           add  CCVEC2 to E
*           copy CCVEC2 to CCVEC1
*  (N > 1)  T' CCVEC1 in CCVEC2
*           copy CCVEC2 to CCVEC1
*           Add  1/n! CCVEC2 to CCVEC3
*           Goto to next N if norm of CCVEC2 is gt. 0
*
*   The order in which T(N)T(N-1) ... T(2)T(1) is defined by 
*   IEXEORD
*           
*
* Jeppe Olsen, May 2000
*              Jan 2001 : IEXEORD added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
      REAL*8 INPROD
*
*. Input
*
      DIMENSION T(N_CC_AMP)
      INTEGER ISOBEX_TP(4*NGAS,NSOBEX_TP)
      INTEGER LSOBEX_TP(NSOBEX_TP),IBSOBEX_TP(NSOBEX_TP)
      INTEGER ISOX_TO_OX(NSOBEX_TP)
      INTEGER ISPOBEX_AC(NSOBEX_TP)
      INTEGER IB_F(NSOBEX_TP,NSOBEX_TP)
      INTEGER IEXEORD(NSOBEX_TP)
*
*. Output 
*
      DIMENSION E(N_CC_AMP), F(*)
*
*. Scratch through input list
*
      DIMENSION CCVEC1(*),CCVEC2(*)
*     ^ Must hold CC vectors
      DIMENSION ISTST1(*),XSTST1(*),ISTST2(*),XSTST2(*)
      DIMENSION ISTST3(*),XSTST3(*),ISTST4(*),XSTST4(*)
*     ^ Must hold ST*ST maps for individual strings, given types,
*     all symmetries. 
      INTEGER IZ(*),IZSCR(*)
      INTEGER ISTOCC1(*),ISTOCC2(*)
*     ^ Hold occupations af all strings of given CAAB, SPGP and sym
      INTEGER ISTREO(*)
*     ^ Must hold reordering for all strings of given CAAB and SPGP
*
*. Local scratch
*
      INTEGER IJOCC(4*MXPNGAS)
*
      CALL LUCIAQENTER('EXP_T ')
      NTEST = 00
      IF (NTEST.GE. 50) THEN
*
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' Welcome to Exp_T '
        WRITE(6,*) ' ================ '
*
        WRITE(6,*) ' I_DO_ZERO, I_DO_FTERMS = ',
     &               I_DO_ZERO, I_DO_FTERMS
        WRITE(6,*) ' NSOBEX_TP = ', NSOBEX_TP
      END IF
*
      N_CC_AMPP1 = N_CC_AMP + 1
      IF(I_DO_ZERO.EQ.1) THEN
*. set initial operator to the unit operator
        ZERO = 0.0D0
        CALL SETVEC(E,ZERO,N_CC_AMPP1)
        E(N_CC_AMPP1) = 1.0D0
*. And zero the intermediate arrays
        IF(I_DO_FTERMS.EQ.1) THEN
          CALL SETVEC(F,ZERO,LEN_F)        
        END IF
      END IF
*     ^ End if zeroing was required
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial E-vector '
        CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
      END IF
    
*  N = 0 term : Current E-term
*  (N = 1)  
*           Copy E to CCVEC1
*           T' + T' CCVEC1 in CCVEC2
*           add  CCVEC2 to E
*           copy CCVEC2 to CCVEC1
*  (N > 1)  T' CCVEC1 in CCVEC2
*           copy CCVEC2 to CCVEC1
*           Add  1/n! CCVEC2 to CCVEC3
*           Goto to next N if norm of CCVEC2 is gt. 0
      DO ITP_ORIG = 1, NSOBEX_TP    
      ITP = IEXEORD(ITP_ORIG)
      IF(ISPOBEX_AC(ITP).EQ.1) THEN
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*)
          WRITE(6,*) ' ===================='
          WRITE(6,*) ' Exp_T for ITP = ', ITP
          WRITE(6,*) ' ===================='
          WRITE(6,*)
        END IF
*
        IB_ITP = IBSOBEX_TP(ITP)
        LEN_ITP = LSOBEX_TP(ITP)
        IF(NTEST.GE.1000) WRITE(6,*) ' IB_ITP, LEN_ITP =',
     &               IB_ITP, LEN_ITP
*. Loop over N from  1
        N = 0 
        FACN = 0
 1000   CONTINUE
         N = N + 1
         IF(NTEST.GE.1000) WRITE(6,*) ' order N = ' , N
         IF(N.EQ.1) THEN
          FACN = 1.0D0
         ELSE
          FACN = FACN*DBLE(N)
         END IF
         FACNI = 1.0D0/FACN
* 
         IF(N.EQ.1) THEN
          CALL COPVEC(E,CCVEC1,N_CC_AMPP1)
         END IF
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' CCVEC1 '
           CALL WRTMAT(CCVEC1,1,N_CC_AMPP1,1,N_CC_AMPP1)
         END IF
*. T CCVEC1 in CCVEC2
         ZERO = 0.0D0
         CALL SETVEC(CCVEC2,ZERO,N_CC_AMPP1)
C?       IF(N.EQ.1) THEN
C?         ONE = 1.0D0
C?         CALL VECSUM(CCVEC2(IB_ITP),CCVEC2(IB_ITP),T(IB_ITP),
C?   &                 ONE,ONE,LEN_ITP) 
C?       END IF
*. Loop over Blocks in CCVEC1
         DO JTP = 1, NSOBEX_TP+1
C!       IF(ISPOBEX_AC(JTP).EQ.1) THEN
           IF(NTEST.GE.1000) WRITE(6,*) ' JTP = ', JTP
*. Spin-orbital type corresponding to ITP*JTP
*. Is ITP * JTP in active 
           CALL PROD_SPOB_EX_TP(ITP,JTP,IJTP,IEXCLEVEL,IJ_IS_ZERO)
           IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' ITP, JTP, IJTP, IJ_IS_ZERO = ', 
     &                  ITP, JTP, IJTP, IJ_IS_ZERO
           END IF
           IB_JTP = IBSOBEX_TP(JTP)
           IF(IJTP.NE.0) THEN
*. Occupation of the four strings in ITP*JTP
C                 T1T2_TO_T12_OCC(I1,I2,I12,NGAS)
             CALL T1T2_TO_T12_OCC(ISOBEX_TP(1,ITP),ISOBEX_TP(1,JTP),
     &                            IJOCC,NGAS)
*. T(ITP) CCVEC1(JTP) to be stored in CCVEC2(IJTP)
             IB_IJTP = IBSOBEX_TP(IJTP)
             IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' E-vector before T1T2... '
               CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
             END IF
*
             IF(JTP.EQ.NSOBEX_TP+1) THEN
*. JTP is unit operator, so just do a DAXPY : T12 = T12 + T(ITP)*Constant
               FACTOR = CCVEC1(IB_JTP)
               ONE = 1.0D0
               LEN_IJTP = LSOBEX_TP(IJTP)
               CALL VECSUM(CCVEC2(IB_IJTP),CCVEC2(IB_IJTP),T(IB_ITP),
     &                     ONE, FACTOR, LEN_IJTP)
             ELSE 
*. Product of two nontrivial operators, do the real job
               CALL T1T2_TO_T12N(ISOBEX_TP(1,ITP),1,T(IB_ITP),
     &                          ISOBEX_TP(1,JTP),1,CCVEC1(IB_JTP),
     &                          IJOCC,IT12SM,CCVEC2(IB_IJTP),
     &                          ISTOCC1,ISTOCC2,ISTREO,
     &                          ISTST1,XSTST1,ISTST2,XSTST2,
     &                          ISTST3,XSTST3,ISTST4,XSTST4,
     &                          IZ,IZSCR)
             END IF
             IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' E-vector after  T1T2... '
               CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
             END IF
           ELSE IF(IEXCLEVEL.LE.2 .AND. I_DO_FTERMS.EQ.1
     &             .AND. IJ_IS_ZERO.EQ.0                ) THEN
*. Add to F terms
             ONE = 1.0D0
             LEN_JTP = LSOBEX_TP(JTP)
*
             IB_EFF = IB_F(JTP,ITP)
             CALL VECSUM(F(IB_EFF),F(IB_EFF),CCVEC1(IB_JTP),
     &                   ONE,FACNI,LEN_JTP)
           END IF
*          ^ End of tests if terms should be included
C!       END IF
*.       ^ End if Block JTP is active
         END DO
*        ^ End of loop over blocks JTP in e(k)
         XNORM = INPROD(CCVEC2,CCVEC2,N_CC_AMPP1)
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' CCVEC2 : '
           CALL WRTMAT(CCVEC2,1,N_CC_AMPP1,1,N_CC_AMPP1)
         END IF
         ONE = 1.0D0
         CALL VECSUM(E,E,CCVEC2,ONE,FACNI,N_CC_AMPP1)
         CALL COPVEC(CCVEC2,CCVEC1,N_CC_AMPP1)
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Updated E vector for N = ', N
           CALL WRTMAT(E,1,N_CC_AMPP1,1,N_CC_AMPP1)
         END IF
*. Temporary break
         MAXTRM = 20
         XTEST = 0.0D0   
        IF(N.EQ.MAXTRM) THEN
          DO K = 1, 1000
           WRITE(6,*) ' Termination due to N=MAXTRM test '
          END DO
        END IF
        IF(XNORM.GT.XTEST.AND.N.LE.MAXTRM) GOTO 1000
      END IF
*.    ^ End if operator is active
      END DO
*     ^ End of loop over operators T(ITP)
*
      IF(NTEST.GE. 50) THEN
        WRITE(6,*) ' Results from EXP_T : '
        CALL WRT_EXPTOP(NSOBEX_TP,LSOBEX_TP,N_CC_AMP,E,F,IB_F)
      END IF
*
C     STOP ' Stop at end of Exp_T '
      CALL LUCIAQEXIT('EXP_T ')
      RETURN
      END
      SUBROUTINE TI_TO_TOK(IOP_CA,LOP_CA,IBOP_CA,IOP_CB,LOP_CB,IBOP_CB,
     &                     IOP_AA,LOP_AA,IBOP_AA,IOP_AB,LOP_AB,IBOP_AB,
     &                     IOP_NM_CAAB,IOP_SM_CAAB,NOP_CAAB,
     &                     K_NM_CAAB,K_SM_CAAB,NK_BAT,
     &                     IKJ_CA_MAP,SKJ_CA_MAP,IBKJ_CA,NK_CA,
     &                     IKJ_CB_MAP,SKJ_CB_MAP,IBKJ_CB,NK_CB,
     &                     IKJ_AA_MAP,SKJ_AA_MAP,IBKJ_AA,NK_AA,
     &                     IKJ_AB_MAP,SKJ_AB_MAP,IBKJ_AB,NK_AB,
     &                     NISFSM_CA,NISFSM_CB,NISFSM_AA,NISFSM_AB,
     &                     TOPK,TI,IB_TI,ISG,LDUM,IOPOFF,IKOFF)
* A set of operators OP_CA, OP_CB, OP_AA, OP_AB is given 
*
* OP_NM_CAAB, OP_SM_CAAB contains a number of such operator quadruplets
*
* Find Mapping OP K => I and use this to gather or scatter
*
* T(IDUM,OP,K) <=> T(IDUM,I)
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
*. Operator
      INTEGER IOP_CA(LOP_CA,*), IOP_CB(LOP_CB,*)
      INTEGER IOP_AA(LOP_AA,*), IOP_AB(LOP_AB,*)
      INTEGER IOP_NM_CAAB(4,NOP_CAAB),IOP_SM_CAAB(4,NOP_CAAB)
*. Offset to strings with given symmetry
      INTEGER IBOP_CA(8),IBOP_CB(8),IBOP_AA(8)
C     INTEGER  IBOP(8)
*. Offset to block with given sym in TI
      INTEGER IB_TI(8,8,8)
*. Number of strings in I per symmetry
      INTEGER NISFSM_CA(8),NISFSM_CB(8),NISFSM_AA(8),NISFSM_AB(8)
*. Kstrings
      INTEGER K_SM_CAAB(4,NK_BAT),K_NM_CAAB(4,NK_BAT)
*. Info on K to J mappings
      INTEGER IKJ_CA_MAP(*),SKJ_CA_MAP(*),IBKJ_CA(8,8,LOP_CA),NK_CA(8)
      INTEGER IKJ_CB_MAP(*),SKJ_CB_MAP(*),IBKJ_CB(8,8,LOP_CA),NK_CB(8)
      INTEGER IKJ_AA_MAP(*),SKJ_AA_MAP(*),IBKJ_AA(8,8,LOP_CA),NK_AA(8)
      INTEGER IKJ_AB_MAP(*),SKJ_AB_MAP(*),IBKJ_AB(8,8,LOP_CA),NK_AB(8)
*. Input or output
      DIMENSION TOPK(LDUM,NOP_CAAB,NK_BAT)
      DIMENSION TI(LDUM,*)
*
      NTEST = 00
      I_CAAB_TI_MX = 0
      DO IOP =1, NOP_CAAB
        IIOP_SM_CA = IOP_SM_CAAB(1,IOP)
        IIOP_NM_CA = IOP_NM_CAAB(1,IOP)
*
        IIOP_SM_CB = IOP_SM_CAAB(2,IOP)
        IIOP_NM_CB = IOP_NM_CAAB(2,IOP)
*
        IIOP_SM_AA = IOP_SM_CAAB(3,IOP)
        IIOP_NM_AA = IOP_NM_CAAB(3,IOP)
*
        IIOP_SM_AB = IOP_SM_CAAB(4,IOP)
        IIOP_NM_AB = IOP_NM_CAAB(4,IOP)
        DO KSTR = 1, NK_BAT
          KK_SM_CA = K_SM_CAAB(1,IOP)
          KKOP_NM_CA = K_NM_CAAB(1,IOP)
*
          KK_SM_CB = K_SM_CAAB(2,IOP)
          KK_NM_CB = K_NM_CAAB(2,IOP)
*
          KK_SM_AA = K_SM_CAAB(3,IOP)
          KK_NM_AA = K_NM_CAAB(3,IOP)
*
          KK_SM_AB = K_SM_CAAB(4,IOP)
          KK_NM_AB = K_NM_CAAB(4,IOP)
*. corresponding J strings
C     K_TO_J_TOT_SINGLE(IKJ,XKJ,KSM,IT,KSTR,
C    &                            IM,XM,IBM,NK,LTOP)
*. I_CA
          IOP_CA_NUM = IBOP_CA(IIOP_SM_CA)-1+IIOP_NM_CA
          CALL K_TO_J_TOT_SINGLE(I_CA,S_CA,K_SM_CA,
     &         IOP_CA(1,IOP_CA_NUM),KK_NUM_CA,
     &         IKJ_MAP_CA,SKJ_MAP_CA,IBKJ_CA,NK_CA,LOP_CA)
*. I_CB
          IOP_CB_NUM = IBOP_CB(IIOP_SM_CB)-1+IIOP_NM_CB
          CALL K_TO_J_TOT_SINGLE(I_CB,S_CB,K_SM_CB,
     &         IOP_CB(1,IOP_CB_NUM),KK_NUM_CB,
     &         IKJ_MAP_CB,SKJ_MAP_CB,IBKJ_CB,NK_CB,LOP_CB)
*. I_AA
          IOP_AA_NUM = IBOP_AA(IIOP_SM_AA)-1+IIOP_NM_AA
          CALL K_TO_J_TOT_SINGLE(I_AA,S_AA,K_SM_AA,
     &         IOP_AA(1,IOP_AA_NUM),KK_NUM_AA,
     &         IKJ_MAP_AA,SKJ_MAP_AA,IBKJ_AA,NK_AA,LOP_AA)
*. I_AB
          IOP_AB_NUM = IBOP_AB(IIOP_SM_AB)-1+IIOP_NM_AB
          CALL K_TO_J_TOT_SINGLE(I_AB,S_AB,K_SM_AB,
     &         IOP_AB(1,IOP_AB_NUM),KK_NUM_AB,
     &         IKJ_MAP_AB,SKJ_MAP_AB,IBKJ_AB,NK_AB,LOP_AB)
          IF(I_CA*I_AA*I_CB*I_AB.NE.0) THEN
*. Adress of I_CB, I_CA, I_AA, I_AB
*. Symmetry of I_ strings
            I_CA_SM = MULTD2H(IIOP_SM_CA,KK_SM_CA)
            I_CB_SM = MULTD2H(IIOP_SM_CB,KK_SM_CB)
            I_AA_SM = MULTD2H(IIOP_SM_AA,KK_SM_AA)
*. Offset to block in TI with this symmetry combination
            IBTI_SSSS = IB_TI(I_CA_SM,I_CB_SM,I_AA_SM)
            I_CAAB_TI = IBTI_SSSS-1+ 
     &      (I_AB-1)*
     &      NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)*NISFSM_AA(I_AA_SM)
     &     +(I_AA-1)*NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)
     &     +(I_CB-1)*NISFSM_CA(I_CA_SM) + I_CA
            SIGN = S_CA*S_CB*S_AA*S_AB
            SIGN = 1.0D0 
            I_CAAB_TI_MX = MAX(I_CAAB_TI_MX, I_CAAB_TI)
            IF(ISG.EQ.1) THEN
             DO IDUM = 1, LDUM
               TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1) 
     &       = SIGN*TI(IDUM,I_CAAB_TI)
             END DO
            ELSE IF (ISG.EQ.2) THEN
             DO IDUM = 1, LDUM
               TI(IDUM,I_CAAB_TI) = 
     &         TI(IDUM,I_CAAB_TI) 
     &       + SIGN*TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1)
             END DO
            END IF
*           ^ End of scatter/gather switch 
          END IF
*         ^ End if Kstrings was nonvanishing
        END DO
*       ^ End of loop over K-strings 
      END DO
*     ^ End of loop over IOP operators
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' TI_TO_TOK speaking : '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) ' Largest column number used = ',  I_CAAB_TI_MX
        WRITE(6,*) ' TI as TI(Idum,I)  '
        CALL WRTMAT(TI,LDUM,I_CAAB_TI_MX,LDUM,I_CAAB_TI_MX)
        WRITE(6,*) ' TOPK as TOPK(IdumIop,Kstr) '
        CALL WRTMAT(TOPK,LDUM*NOP_CAAB,NK_BAT,LDUM*NOP_CAAB,NK_BAT)
      END IF
*
      RETURN
      END
      SUBROUTINE TI_TO_TOKBX(
     &           IOP_NM_CAAB,IOP_SM_CAAB,NOP_BAT, 
     &           NOP_CA,NOP_CB,NOP_AA,NOP_AB,
     &           K_NM_CAAB,K_SM_CAAB,NK_BAT,
     &           IKJ_CA_MAP,SKJ_CA_MAP,IBKJ_CA,
     &           IKJ_CB_MAP,SKJ_CB_MAP,IBKJ_CB,
     &           IKJ_AA_MAP,SKJ_AA_MAP,IBKJ_AA,
     &           IKJ_AB_MAP,SKJ_AB_MAP,IBKJ_AB,
     &           NISFSM_CA,NISFSM_CB,NISFSM_AA,NISFSM_AB,
     &           TOPK,TI,IB_TI,ISG,LDUM,SIGNI,IOPDAG,IOPOFF,IKOFF)
*
* Version with increased print level, may be discarded
* A set of operators OP_CA, OP_CB, OP_AA, OP_AB is given 
* OP_NM_CAAB, OP_SM_CAAB contains a number of such operator quadruplets
*   
*
* ISG = 1 : 
*   T(IDUM,OP,K) :=  Sign(OP,K)* T(IDUM,I)
* ISG = 2 :
*   T(IDUM,I) := T(IDUM,OP,K) + Sign(OP,K)*T(IDUM,I)
*
* If IOPDAG = 1 IOP_NM, IOP_SM contains info for daggered operator 
* instead of operator itself
*
* Version with complete precomputed maps OP * K => I
*
* Jeppe Olsen, October 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
*. Operator
      INTEGER IOP_NM_CAAB(4,NOP_BAT),IOP_SM_CAAB(4,NOP_BAT)
*. Number of operators per symmetry 
      INTEGER NOP_CA(*), NOP_CB(*), NOP_AA(*), NOP_AB(*)
*. Offset to strings with given symmetry
      INTEGER IBOP_CA(8),IBOP_CB(8),IBOP_AA(8)
*. Offset to block with given sym in TI
      INTEGER IB_TI(8,8,8)
*. Number of strings in I per symmetry
      INTEGER NISFSM_CA(8),NISFSM_CB(8),NISFSM_AA(8),NISFSM_AB(8)
*. Kstrings
      INTEGER K_SM_CAAB(4,NK_BAT),K_NM_CAAB(4,NK_BAT)
*. Info on OP * K to I mappings
      INTEGER IKJ_CA_MAP(*),IBKJ_CA(8,8)
      INTEGER IKJ_CB_MAP(*),IBKJ_CB(8,8)
      INTEGER IKJ_AA_MAP(*),IBKJ_AA(8,8)
      INTEGER IKJ_AB_MAP(*),IBKJ_AB(8,8)
*
      DIMENSION SKJ_CA_MAP(*),SKJ_CB_MAP(*),SKJ_AA_MAP(*),SKJ_AB_MAP(*)
*. Input or output
      DIMENSION TOPK(LDUM,NOP_BAT,NK_BAT)
      DIMENSION TI(LDUM,*)
*. Test
C?    IF(SIGNI.NE.1.0D0.AND.SIGNI.NE.-1.0D0) THEN
C?      WRITE(6,*) ' nonunit SIGNI = ', SIGNI
C?      STOP       ' nonunit SIGNI   '
C?    END IF
*

      NTEST = 00
      MX_TI_COL = 0
      I_CAAB_TI_MX = 0
C?    WRITE(6,*) ' SIGNI(1) = ', SIGNI
*
      DO IOP =1, NOP_BAT
        IF(IOPDAG.EQ.0) THEN
          IIOP_SM_CA = IOP_SM_CAAB(1,IOP)
          IIOP_NM_CA = IOP_NM_CAAB(1,IOP)
*
          IIOP_SM_CB = IOP_SM_CAAB(2,IOP)
          IIOP_NM_CB = IOP_NM_CAAB(2,IOP)
*
          IIOP_SM_AA = IOP_SM_CAAB(3,IOP)
          IIOP_NM_AA = IOP_NM_CAAB(3,IOP)
*
          IIOP_SM_AB = IOP_SM_CAAB(4,IOP)
          IIOP_NM_AB = IOP_NM_CAAB(4,IOP)
        ELSE
          IIOP_SM_CA = IOP_SM_CAAB(3,IOP)
          IIOP_NM_CA = IOP_NM_CAAB(3,IOP)
*
          IIOP_SM_CB = IOP_SM_CAAB(4,IOP)
          IIOP_NM_CB = IOP_NM_CAAB(4,IOP)
*
          IIOP_SM_AA = IOP_SM_CAAB(1,IOP)
          IIOP_NM_AA = IOP_NM_CAAB(1,IOP)
*
          IIOP_SM_AB = IOP_SM_CAAB(2,IOP)
          IIOP_NM_AB = IOP_NM_CAAB(2,IOP)
        END IF
*
        IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' IOP = ', IOP
        WRITE(6,*) ' IIOP_SM_CA, IIOP_SM_CB, IIOP_SM_AA, IIOP_SM_AB ',
     &               IIOP_SM_CA, IIOP_SM_CB, IIOP_SM_AA, IIOP_SM_AB
        WRITE(6,*) ' IIOP_NM_CA, IIOP_NM_CB, IIOP_NM_AA, IIOP_NM_AB ',
     &               IIOP_NM_CA, IIOP_NM_CB, IIOP_NM_AA, IIOP_NM_AB
        END IF
*
        DO KSTR = 1, NK_BAT
          KK_SM_CA = K_SM_CAAB(1,KSTR)
          KK_NM_CA = K_NM_CAAB(1,KSTR)
*
          KK_SM_CB = K_SM_CAAB(2,KSTR)
          KK_NM_CB = K_NM_CAAB(2,KSTR)
*
          KK_SM_AA = K_SM_CAAB(3,KSTR)
          KK_NM_AA = K_NM_CAAB(3,KSTR)
*
          KK_SM_AB = K_SM_CAAB(4,KSTR)
          KK_NM_AB = K_NM_CAAB(4,KSTR)
          IF(NTEST.GE.1000) THEN
C23456789012345678901234567890123456789012345678901234567890123456789012
             WRITE(6,*) ' KSTR = ', KSTR
             WRITE(6,*) ' KK_SM_CA, KK_SM_CB, KK_SM_AA, KK_SM_AB ',
     &                    KK_SM_CA, KK_SM_CB, KK_SM_AA, KK_SM_AB
             WRITE(6,*) ' KK_NM_CA, KK_NM_CB, KK_NM_AA, KK_NM_AB ',
     &                    KK_NM_CA, KK_NM_CB, KK_NM_AA, KK_NM_AB
          END IF
* IOP * K => I
          IBOPK_CA = IBKJ_CA(IIOP_SM_CA, KK_SM_CA)
          LOP_CA = NOP_CA(IIOP_SM_CA)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IBOPK_CA, KK_NM_CA, LOP_CA, IIOP_NM_CA ',
     &                   IBOPK_CA, KK_NM_CA, LOP_CA, IIOP_NM_CA 
          END IF
          I_CA = IKJ_CA_MAP(IBOPK_CA-1+(KK_NM_CA-1)*LOP_CA+IIOP_NM_CA)
          S_CA = SKJ_CA_MAP(IBOPK_CA-1+(KK_NM_CA-1)*LOP_CA+IIOP_NM_CA)
*
          IBOPK_CB = IBKJ_CB(IIOP_SM_CB, KK_SM_CB)
          LOP_CB = NOP_CB(IIOP_SM_CB)
          
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IBOPK_CB, KK_NM_CB, LOP_CB, IIOP_NM_CB ',
     &                   IBOPK_CB, KK_NM_CB, LOP_CB, IIOP_NM_CB 
          END IF
          I_CB = IKJ_CB_MAP(IBOPK_CB-1+(KK_NM_CB-1)*LOP_CB+IIOP_NM_CB)
          S_CB = SKJ_CB_MAP(IBOPK_CB-1+(KK_NM_CB-1)*LOP_CB+IIOP_NM_CB)
*
          IBOPK_AA = IBKJ_AA(IIOP_SM_AA, KK_SM_AA)
          LOP_AA = NOP_AA(IIOP_SM_AA)
          IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IBOPK_AA, KK_NM_AA, LOP_AA, IIOP_NM_AA = ',
     &                 IBOPK_AA, KK_NM_AA, LOP_AA, IIOP_NM_AA
          END IF
          I_AA = IKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
          S_AA = SKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
*
          IBOPK_AB = IBKJ_AB(IIOP_SM_AB, KK_SM_AB)
          LOP_AB = NOP_AB(IIOP_SM_AB)
          I_AB = IKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
          S_AB = SKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
          IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOP, KSTR, I_CA, I_AA, I_CB, I_AB ',
     &                 IOP, KSTR, I_CA, I_AA, I_CB, I_AB
          END IF
          IF(I_CA*I_AA*I_CB*I_AB.NE.0) THEN
*. Adress of I_CB, I_CA, I_AA, I_AB
*. Symmetry of I_ strings
            I_CA_SM = MULTD2H(IIOP_SM_CA,KK_SM_CA)
            I_CB_SM = MULTD2H(IIOP_SM_CB,KK_SM_CB)
            I_AA_SM = MULTD2H(IIOP_SM_AA,KK_SM_AA)
*. Offset to block in TI with this symmetry combination
            IBTI_SSSS = IB_TI(I_CA_SM,I_CB_SM,I_AA_SM)
            IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' IBTI_SSSS = ', IBTI_SSSS                       
             WRITE(6,*) ' I_CA, I_CB, I_AA, I_AB = ',
     &                    I_CA, I_CB, I_AA, I_AB
             WRITE(6,*) 
     &      'NISFSM_CA(I_CA_SM),NISFSM_CB(I_CB_SM),NISFSM_AA(I_AA_SM)',
     &       NISFSM_CA(I_CA_SM),NISFSM_CB(I_CB_SM),NISFSM_AA(I_AA_SM)  
            END IF
     
            I_CAAB_TI = IBTI_SSSS-1+ 
     &      (I_AB-1)*
     &      NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)*NISFSM_AA(I_AA_SM)
     &     +(I_AA-1)*NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)
     &     +(I_CB-1)*NISFSM_CA(I_CA_SM) + I_CA
            IF(NTEST.GE.1000) WRITE(6,*) ' I_CAAB_TI = ', I_CAAB_TI 
            I_CAAB_TI_MX = MAX(I_CAAB_TI_MX, I_CAAB_TI)
*
            SIGN = S_CA*S_CB*S_AA*S_AB*SIGNI
C?          WRITE(6,*) ' SIGNI(2) = ', SIGNI
C?          WRITE(6,*) ' S_CA, S_CB, S_AA, S_AB = ',
C?   &                   S_CA, S_CB, S_AA, S_AB
            IF(ISG.EQ.1) THEN
             DO IDUM = 1, LDUM
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' IDUM, IOP, KSTR, I_CAAB_TI = ',
     &                        IDUM, IOP, KSTR, I_CAAB_TI 
               END IF
*
               TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1) 
     &       = SIGN*TI(IDUM,I_CAAB_TI)
*
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' I_CAAB_TI, TI(IDUM,I_CAAB_TI) : ',
     &                        I_CAAB_TI, TI(IDUM,I_CAAB_TI)
                 WRITE(6,*) ' updated TopK : Idum,Iop,Kstr, value : ',
     &           IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1,
     &           TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1)
                 WRITE(6,*) ' Sign = ', sign
               END IF
             END DO
            ELSE IF (ISG.EQ.2) THEN
             DO IDUM = 1, LDUM
               TI(IDUM,I_CAAB_TI) = 
     &         TI(IDUM,I_CAAB_TI) 
     &       + SIGN*TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1)
C?             WRITE(6,*) ' updated TI : Idum, I_CAAB_TI, Value :',
C?   &         IDUM,I_CAAB_TI,TI(IDUM,I_CAAB_TI)
             END DO
            END IF
*           ^ End of scatter/gather switch 
          ELSE IF( ISG.EQ.1) THEN 
             DO IDUM = 1, LDUM
               TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1) = 0.0D0                     
             END DO
          END IF
*         ^ End if Istrings was nonvanishing
        END DO
*       ^ End of loop over K-strings 
      END DO
*     ^ End of loop over IOP operators
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' TI_TO_TOKB speaking : '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' NOP_BAT, NK_BAT, LDUM = ', NOP_BAT, NK_BAT, LDUM
        WRITE(6,*) ' Largest column number used = ',  I_CAAB_TI_MX
        WRITE(6,*) ' TI as TI(Idum,I)  '
        CALL WRTMAT(TI,LDUM,I_CAAB_TI_MX,LDUM,I_CAAB_TI_MX)
        WRITE(6,*) ' TOPK as TOPK(IdumIop,Kstr) '
        CALL WRTMAT(TOPK,LDUM*NOP_BAT,NK_BAT,LDUM*NOP_BAT,NK_BAT)
      END IF

*
      RETURN
      END
      SUBROUTINE TI_TO_TOKB(
     &           IOP_NM_CAAB,IOP_SM_CAAB,NOP_BAT, 
     &           NOP_CA,NOP_CB,NOP_AA,NOP_AB,
     &           K_NM_CAAB,K_SM_CAAB,NK_BAT,
     &           IKJ_CA_MAP,SKJ_CA_MAP,IBKJ_CA,
     &           IKJ_CB_MAP,SKJ_CB_MAP,IBKJ_CB,
     &           IKJ_AA_MAP,SKJ_AA_MAP,IBKJ_AA,
     &           IKJ_AB_MAP,SKJ_AB_MAP,IBKJ_AB,
     &           NISFSM_CA,NISFSM_CB,NISFSM_AA,NISFSM_AB,
     &           TOPK,TI,IB_TI,ISG,LDUM,SIGNI,IOPDAG,IOPOFF,IKOFF)
* A set of operators OP_CA, OP_CB, OP_AA, OP_AB is given 
* OP_NM_CAAB, OP_SM_CAAB contains a number of such operator quadruplets
*   
*
* ISG = 1 : 
*   T(IDUM,OP,K) :=  Sign(OP,K)* T(IDUM,I)
* ISG = 2 :
*   T(IDUM,I) := T(IDUM,OP,K) + Sign(OP,K)*T(IDUM,I)
*
* If IOPDAG = 1 IOP_NM, IOP_SM contains info for daggered operator 
* instead of operator itself
*
* Version with complete precomputed maps OP * K => I
*
* Jeppe Olsen, October 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
*. Operator
      INTEGER IOP_NM_CAAB(4,NOP_BAT),IOP_SM_CAAB(4,NOP_BAT)
*. Number of operators per symmetry 
      INTEGER NOP_CA(*), NOP_CB(*), NOP_AA(*), NOP_AB(*)
*. Offset to strings with given symmetry
      INTEGER IBOP_CA(8),IBOP_CB(8),IBOP_AA(8)
*. Offset to block with given sym in TI
      INTEGER IB_TI(8,8,8)
*. Number of strings in I per symmetry
      INTEGER NISFSM_CA(8),NISFSM_CB(8),NISFSM_AA(8),NISFSM_AB(8)
*. Kstrings
      INTEGER K_SM_CAAB(4,NK_BAT),K_NM_CAAB(4,NK_BAT)
*. Info on OP * K to I mappings
      INTEGER IKJ_CA_MAP(*),IBKJ_CA(8,8)
      INTEGER IKJ_CB_MAP(*),IBKJ_CB(8,8)
      INTEGER IKJ_AA_MAP(*),IBKJ_AA(8,8)
      INTEGER IKJ_AB_MAP(*),IBKJ_AB(8,8)
*
      DIMENSION SKJ_CA_MAP(*),SKJ_CB_MAP(*),SKJ_AA_MAP(*),SKJ_AB_MAP(*)
*. Input or output
      DIMENSION TOPK(LDUM,NOP_BAT,NK_BAT)
      DIMENSION TI(LDUM,*)
*. Test
C?    IF(SIGNI.NE.1.0D0.AND.SIGNI.NE.-1.0D0) THEN
C?      WRITE(6,*) ' nonunit SIGNI = ', SIGNI
C?      STOP       ' nonunit SIGNI   '
C?    END IF
*

      NTEST = 000
      MX_TI_COL = 0
      I_CAAB_TI_MX = 0
C?    WRITE(6,*) ' SIGNI(1) = ', SIGNI
*
      DO IOP =1, NOP_BAT
        IF(IOPDAG.EQ.0) THEN
          IIOP_SM_CA = IOP_SM_CAAB(1,IOP)
          IIOP_NM_CA = IOP_NM_CAAB(1,IOP)
*
          IIOP_SM_CB = IOP_SM_CAAB(2,IOP)
          IIOP_NM_CB = IOP_NM_CAAB(2,IOP)
*
          IIOP_SM_AA = IOP_SM_CAAB(3,IOP)
          IIOP_NM_AA = IOP_NM_CAAB(3,IOP)
*
          IIOP_SM_AB = IOP_SM_CAAB(4,IOP)
          IIOP_NM_AB = IOP_NM_CAAB(4,IOP)
        ELSE
          IIOP_SM_CA = IOP_SM_CAAB(3,IOP)
          IIOP_NM_CA = IOP_NM_CAAB(3,IOP)
*
          IIOP_SM_CB = IOP_SM_CAAB(4,IOP)
          IIOP_NM_CB = IOP_NM_CAAB(4,IOP)
*
          IIOP_SM_AA = IOP_SM_CAAB(1,IOP)
          IIOP_NM_AA = IOP_NM_CAAB(1,IOP)
*
          IIOP_SM_AB = IOP_SM_CAAB(2,IOP)
          IIOP_NM_AB = IOP_NM_CAAB(2,IOP)
        END IF
*
        IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' IOP = ', IOP
        WRITE(6,*) ' IIOP_SM_CA, IIOP_SM_CB, IIOP_SM_AA, IIOP_SM_AB ',
     &               IIOP_SM_CA, IIOP_SM_CB, IIOP_SM_AA, IIOP_SM_AB
        WRITE(6,*) ' IIOP_NM_CA, IIOP_NM_CB, IIOP_NM_AA, IIOP_NM_AB ',
     &               IIOP_NM_CA, IIOP_NM_CB, IIOP_NM_AA, IIOP_NM_AB
        END IF
*
        DO KSTR = 1, NK_BAT
          KK_SM_CA = K_SM_CAAB(1,KSTR)
          KK_NM_CA = K_NM_CAAB(1,KSTR)
*
          KK_SM_CB = K_SM_CAAB(2,KSTR)
          KK_NM_CB = K_NM_CAAB(2,KSTR)
*
          KK_SM_AA = K_SM_CAAB(3,KSTR)
          KK_NM_AA = K_NM_CAAB(3,KSTR)
*
          KK_SM_AB = K_SM_CAAB(4,KSTR)
          KK_NM_AB = K_NM_CAAB(4,KSTR)
          IF(NTEST.GE.1000) THEN
C23456789012345678901234567890123456789012345678901234567890123456789012
             WRITE(6,*) ' KSTR = ', KSTR
             WRITE(6,*) ' KK_SM_CA, KK_SM_CB, KK_SM_AA, KK_SM_AB ',
     &                    KK_SM_CA, KK_SM_CB, KK_SM_AA, KK_SM_AB
             WRITE(6,*) ' KK_NM_CA, KK_NM_CB, KK_NM_AA, KK_NM_AB ',
     &                    KK_NM_CA, KK_NM_CB, KK_NM_AA, KK_NM_AB
          END IF
* IOP * K => I
          IBOPK_CA = IBKJ_CA(IIOP_SM_CA, KK_SM_CA)
          LOP_CA = NOP_CA(IIOP_SM_CA)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IBOPK_CA, KK_NM_CA, LOP_CA, IIOP_NM_CA ',
     &                   IBOPK_CA, KK_NM_CA, LOP_CA, IIOP_NM_CA 
          END IF
          I_CA = IKJ_CA_MAP(IBOPK_CA-1+(KK_NM_CA-1)*LOP_CA+IIOP_NM_CA)
          S_CA = SKJ_CA_MAP(IBOPK_CA-1+(KK_NM_CA-1)*LOP_CA+IIOP_NM_CA)
*
          IBOPK_CB = IBKJ_CB(IIOP_SM_CB, KK_SM_CB)
          LOP_CB = NOP_CB(IIOP_SM_CB)
          
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IBOPK_CB, KK_NM_CB, LOP_CB, IIOP_NM_CB ',
     &                   IBOPK_CB, KK_NM_CB, LOP_CB, IIOP_NM_CB 
          END IF
          I_CB = IKJ_CB_MAP(IBOPK_CB-1+(KK_NM_CB-1)*LOP_CB+IIOP_NM_CB)
          S_CB = SKJ_CB_MAP(IBOPK_CB-1+(KK_NM_CB-1)*LOP_CB+IIOP_NM_CB)
*
          IBOPK_AA = IBKJ_AA(IIOP_SM_AA, KK_SM_AA)
          LOP_AA = NOP_AA(IIOP_SM_AA)
          IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IBOPK_AA, KK_NM_AA, LOP_AA, IIOP_NM_AA = ',
     &                 IBOPK_AA, KK_NM_AA, LOP_AA, IIOP_NM_AA
          END IF
          I_AA = IKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
          S_AA = SKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
*
          IBOPK_AB = IBKJ_AB(IIOP_SM_AB, KK_SM_AB)
          LOP_AB = NOP_AB(IIOP_SM_AB)
          I_AB = IKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
          S_AB = SKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
          IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOP, KSTR, I_CA, I_AA, I_CB, I_AB ',
     &                 IOP, KSTR, I_CA, I_AA, I_CB, I_AB
          END IF
          IF(I_CA*I_AA*I_CB*I_AB.NE.0) THEN
*. Adress of I_CB, I_CA, I_AA, I_AB
*. Symmetry of I_ strings
            I_CA_SM = MULTD2H(IIOP_SM_CA,KK_SM_CA)
            I_CB_SM = MULTD2H(IIOP_SM_CB,KK_SM_CB)
            I_AA_SM = MULTD2H(IIOP_SM_AA,KK_SM_AA)
*. Offset to block in TI with this symmetry combination
            IBTI_SSSS = IB_TI(I_CA_SM,I_CB_SM,I_AA_SM)
            IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' IBTI_SSSS = ', IBTI_SSSS                       
             WRITE(6,*) ' I_CA, I_CB, I_AA, I_AB = ',
     &                    I_CA, I_CB, I_AA, I_AB
             WRITE(6,*) 
     &      'NISFSM_CA(I_CA_SM),NISFSM_CB(I_CB_SM),NISFSM_AA(I_AA_SM)',
     &       NISFSM_CA(I_CA_SM),NISFSM_CB(I_CB_SM),NISFSM_AA(I_AA_SM)  
            END IF
     
            I_CAAB_TI = IBTI_SSSS-1+ 
     &      (I_AB-1)*
     &      NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)*NISFSM_AA(I_AA_SM)
     &     +(I_AA-1)*NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)
     &     +(I_CB-1)*NISFSM_CA(I_CA_SM) + I_CA
            IF(NTEST.GE.1000) WRITE(6,*) ' I_CAAB_TI = ', I_CAAB_TI 
            I_CAAB_TI_MX = MAX(I_CAAB_TI_MX, I_CAAB_TI)
*
            SIGN = S_CA*S_CB*S_AA*S_AB*SIGNI
C?          WRITE(6,*) ' SIGNI(2) = ', SIGNI
C?          WRITE(6,*) ' S_CA, S_CB, S_AA, S_AB = ',
C?   &                   S_CA, S_CB, S_AA, S_AB
            IF(ISG.EQ.1) THEN
             DO IDUM = 1, LDUM
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' IDUM, IOP, KSTR, I_CAAB_TI = ',
     &                        IDUM, IOP, KSTR, I_CAAB_TI 
               END IF
*
               TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1) 
     &       = SIGN*TI(IDUM,I_CAAB_TI)
*
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' I_CAAB_TI, TI(IDUM,I_CAAB_TI) : ',
     &                        I_CAAB_TI, TI(IDUM,I_CAAB_TI)
                 WRITE(6,*) ' updated TopK : Idum,Iop,Kstr, value : ',
     &           IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1,
     &           TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1)
                 WRITE(6,*) ' Sign = ', sign
               END IF
             END DO
            ELSE IF (ISG.EQ.2) THEN
             DO IDUM = 1, LDUM
               TI(IDUM,I_CAAB_TI) = 
     &         TI(IDUM,I_CAAB_TI) 
     &       + SIGN*TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1)
C?             WRITE(6,*) ' updated TI : Idum, I_CAAB_TI, Value :',
C?   &         IDUM,I_CAAB_TI,TI(IDUM,I_CAAB_TI)
             END DO
            END IF
*           ^ End of scatter/gather switch 
          ELSE IF( ISG.EQ.1) THEN 
             DO IDUM = 1, LDUM
               TOPK(IDUM,IOP-IOPOFF+1,KSTR-IKOFF+1) = 0.0D0                     
             END DO
          END IF
*         ^ End if Istrings was nonvanishing
        END DO
*       ^ End of loop over K-strings 
      END DO
*     ^ End of loop over IOP operators
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' TI_TO_TOKB speaking : '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' NOP_BAT, NK_BAT, LDUM = ', NOP_BAT, NK_BAT, LDUM
        WRITE(6,*) ' IOPOFF, IKOFF = ', IOPOFF, IKOFF
        WRITE(6,*) ' Largest column number used = ',  I_CAAB_TI_MX
        WRITE(6,*) ' TI as TI(Idum,I)  '
        CALL WRTMAT(TI,LDUM,I_CAAB_TI_MX,LDUM,I_CAAB_TI_MX)
        WRITE(6,*) ' TOPK as TOPK(IdumIop,Kstr) '
        CALL WRTMAT(TOPK,LDUM*NOP_BAT,NK_BAT,LDUM*NOP_BAT,NK_BAT)
      END IF

*
      RETURN
      END
      SUBROUTINE NST_CAAB(IOC_CAAB,LST_CAAB)
* Number of strings per sym of TCC block defined by IOC_CAAB
*
* Jeppe Olsen, May 2000
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Input
      INTEGER IOC_CAAB(NGAS,4)
*. Output
      INTEGER LST_CAAB(8,4)
*
      DO JCAAB = 1, 4
        CALL NST_SPGP(IOC_CAAB(1,JCAAB),LST_CAAB(1,JCAAB))
      END DO
*
      RETURN
      END 
      SUBROUTINE STR_CAAB(ICAAB,ISTR_CAAB)
*
* Obtain strings of CAAB types given by ICAAB  for all symmetries 
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'csm.inc'

*. Input
      INTEGER ICAAB(NGAS,4)
*. Output
      INTEGER ISTR_CAAB(MX_ST_TSOSO_BLK_MX*NSMST,4)
C     INTEGER IBSTR_CAAB(NSMST,4)
*. Local scratch
      INTEGER IGRP_AR(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' STR_CAAB in action : '
        WRITE(6,*) ' ====================='
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Strings will be generated for the four CAAB substrings of '
        CALL WRT_SPOX_TP(ICAAB,1)
      END IF
*
      IDUM = 0
      IBSTR_CA = 1
      IBSTR_CB = 1
      IBSTR_AA = 1
      IBSTR_AB = 1
*CA strings
      DO ISTSM = 1, NSMST
        CALL OCC_TO_GRP(ICAAB(1,1),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,1),NGAS) 
C?      WRITE(6,*) ' NEL(CA) = ', NEL
        CALL GETSTR2_TOTSM_SPGP(IGRP_AR,NGAS,ISTSM,NEL,NSTR_CA,
     &      ISTR_CAAB((IBSTR_CA-1)*NEL+1,1), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,1) = IBSTR_CA
        IBSTR_CA = IBSTR_CA + NSTR_CA
*. CB strings
        CALL OCC_TO_GRP(ICAAB(1,2),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,2),NGAS) 
C?      WRITE(6,*) ' NEL(CB) = ', NEL
        CALL GETSTR2_TOTSM_SPGP(IGRP_AR,NGAS,ISTSM,NEL,NSTR_CB,
     &       ISTR_CAAB((IBSTR_CB-1)*NEL+1,2), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,2) = IBSTR_CB
        IBSTR_CB = IBSTR_CB + NSTR_CB
*. AA strings
        CALL OCC_TO_GRP(ICAAB(1,3),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,3),NGAS) 
C?      WRITE(6,*) ' NEL(AA) = ', NEL
        CALL GETSTR2_TOTSM_SPGP(IGRP_AR,NGAS,ISTSM,NEL,NSTR_AA,
     &       ISTR_CAAB((IBSTR_AA-1)*NEL+1,3), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,3) = IBSTR_AA
        IBSTR_AA = IBSTR_AA + NSTR_AA
*. AB strings
        CALL OCC_TO_GRP(ICAAB(1,4),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,4),NGAS) 
C?      WRITE(6,*) ' NEL(AB) = ', NEL
        CALL GETSTR2_TOTSM_SPGP(IGRP_AR,NGAS,ISTSM,NEL,NSTR_AB,
     &       ISTR_CAAB((IBSTR_AB-1)*NEL+1,4), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,4) = IBSTR_AB
        IBSTR_AB = IBSTR_AB + NSTR_AB
*
      END DO
*
      RETURN
      END
      SUBROUTINE NEW_CAAB_OC(IOCC_L,IOCC_R,IOC_OP,ICE,ILR,NGAS)
*
* ICE = 1 contraction map 
*        IOCC_L = IOC_OP contracted with IOCC_R
* ICE = 2 Excitation map 
*        IOCC_L = IOC_OP multiplied with IOCC_R 
*
* ILR = 1 : Find IOCC_L
* ILR = 2 : Find IOCC_R
*
      INCLUDE 'implicit.inc'
*. Input or output
      INTEGER IOCC_L(NGAS,4),IOC_OP(NGAS,4)
      INTEGER IOCC_R(NGAS,4)
*
      IONE = 1
      MONE = -1
*
* Contraction map
*
      IF(ICE.EQ.1) THEN
*. CA_L = CA_R - AA_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,1),IOCC_R(1,1),IOC_OP(1,3),IONE,MONE,
     &                NGAS)
        ELSE IF (ILR.EQ.2) THEN
          CALL IVCSUM(IOCC_R(1,1),IOCC_L(1,1),IOC_OP(1,3),IONE,IONE,
     &               NGAS)
        END IF
*. CB_OUT = CB_IN -AB_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,2),IOCC_R(1,2),IOC_OP(1,4),IONE,MONE,
     &                NGAS)
        ELSE 
          CALL IVCSUM(IOCC_R(1,2),IOCC_L(1,2),IOC_OP(1,4),IONE,IONE,
     &                NGAS)
        END IF
*. AA_OUT = AA_IN - CA_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,3),IOCC_R(1,3),IOC_OP(1,1),IONE,MONE,
     &                NGAS)
         ELSE 
          CALL IVCSUM(IOCC_R(1,3),IOCC_L(1,3),IOC_OP(1,1),IONE,IONE,
     &                NGAS)
        END IF
*. AB_OUT = AB_IN - CB_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,4),IOCC_R(1,4),IOC_OP(1,2),IONE,MONE,
     &                NGAS)
        ELSE 
          CALL IVCSUM(IOCC_R(1,4),IOCC_L(1,4),IOC_OP(1,2),IONE,IONE,
     &                NGAS)
        END IF
      ELSE IF (ICE.EQ.2) THEN
*
* Excitation map
*
* IL_CA = IR _CA + IOP_CA
        
        IF(ILR.EQ.1) THEN
          DO ICAAB = 1, 4
           CALL IVCSUM(IOCC_L(1,ICAAB),IOCC_R(1,ICAAB),IOC_OP(1,ICAAB),
     &          IONE,IONE,NGAS)
          END DO
        ELSE IF( ILR.EQ.2) THEN
          DO ICAAB = 1, 4
           CALL IVCSUM(IOCC_R(1,ICAAB),IOCC_L(1,ICAAB),IOC_OP(1,ICAAB),
     &          IONE,MONE,NGAS)
          END DO
        END IF
      END IF
*     ^ End of ICE switch 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'NEW_CAAB_OC : ICE, ILR = ', ICE,ILR
        WRITE(6,*) ' IOCC_L,IOCC_OP,IOCC_R : '
        CALL IWRTMA(IOCC_L,NGAS,4,NGAS,4)
        CALL IWRTMA(IOC_OP,NGAS,4,NGAS,4)
        CALL IWRTMA(IOCC_R,NGAS,4,NGAS,4)
      END IF
*
      RETURN
      END
      SUBROUTINE CONJ_OP(IGAS_IN, ICA_IN,NOP,IGAS_OUT,ICA_OUT)
*
* An operator is defined by IGAS_IN and ICA_IN
*
* Obtain conjugated operator 
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IGAS_IN(NOP),ICA_IN(NOP)
*. Output
      INTEGER IGAS_OUT(NOP),ICA_OUT(NOP)
*
      DO IOP = 1, NOP     
        IGAS_OUT(IOP) = IGAS_IN(NOP+1-IOP)
        ICA_OUT(IOP) = ICA_IN(NOP+1-IOP)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Gas and CA of conjugated operator '
        CALL IWRTMA(IGAS_OUT,1,NOP,1,NOP) 
        CALL IWRTMA(ICA_OUT,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
      SUBROUTINE DE_ST_CAAB_INDX(
     &           IDE_SP_CA, IDE_SP_CB, IDE_SP_AA, IDE_SP_AB,
     &           IDE_CA_CA, IDE_CA_CB, IDE_CA_AA, IDE_CA_AB,
     &           LDE_CA, LDE_CB, LDE_AA, LDE_AB,
     &           IST_CA, IST_CB, IST_AA, IST_AB,
     &           IBK_CA, IBK_CB, IBK_AA, IBK_AB,
     &           NK_CA,  NK_CB,  NK_AA,  NK_AB,
     &           IMAP_CA,IMAP_CB,IMAP_AA,IMAP_AB,
     &           SMAP_CA,SMAP_CB,SMAP_AA,SMAP_AB)
*
* A deexcitation/excitation operator IDX is given
* for each of the four CAAB types of a string
*
* Find mappings to resulting strings,
* Mappings are given as individual index mappings
*
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'  
      INCLUDE 'cgas.inc'
*. input : Gasspaces and C/A for each elementary operators
*.   
      INTEGER IDE_SP_CA(*),IDE_SP_CB(*),IDE_SP_AA(*),IDE_SP_AB(*)
      INTEGER IDE_CA_CA(*),IDE_CA_CB(*),IDE_CA_AA(*),IDE_CA_AB(*)
*. Occupations of each GAS space of Input Strings 
      INTEGER IST_CA(NGAS),IST_CB(NGAS),IST_AA(NGAS),IST_AB(NGAS)
*. Output
      INTEGER IBK_CA(8,8,MXPLCCOP), IBK_CB(8,8,MXPLCCOP)
      INTEGER IBK_AA(8,8,MXPLCCOP), IBK_AB(8,8,MXPLCCOP)
      INTEGER NK_CA(MXPOBS,MXPLCCOP), NK_CB(MXPOBS,MXPLCCOP)
      INTEGER NK_AA(MXPOBS,MXPLCCOP), NK_AB(MXPOBS,MXPLCCOP)
      INTEGER IMAP_CA(*),IMAP_CB(*),IMAP_AA(*),IMAP_AB(*)
      INTEGER SMAP_CA(*),SMAP_CB(*),SMAP_AA(*),SMAP_AB(*)
*. Local scratch 
      DIMENSION IOP_AR(MXPLCCOP),ICA_AR(MXPLCCOP)
*
       NTEST = 00
       IF(NTEST.GE.100) THEN
          WRITE(6,*)
          WRITE(6,*) ' DE_ST_CAAB_INDX in action '
          WRITE(6,*) ' ========================= '
          WRITE(6,*) 
          WRITE(6,*) ' Operator to be contracted with CA of string'
C              WRT_CNTR2(ICONT_SP,ICONT_CA,NCONT)
          CALL WRT_CNTR2(IDE_SP_CA,IDE_CA_CA,LDE_CA)
          WRITE(6,*) ' Operator to be contracted with CB of string'
          CALL WRT_CNTR2(IDE_SP_CB,IDE_CA_CB,LDE_CA)
          WRITE(6,*) ' Operator to be contracted with AA of string'
          CALL WRT_CNTR2(IDE_SP_AA,IDE_CA_AA,LDE_AA)
          WRITE(6,*) ' Operator to be contracted with AB of string'
          CALL WRT_CNTR2(IDE_SP_AB,IDE_CA_AB,LDE_AB)
*
          WRITE(6,*) ' CA, CB, AA and AB of string : '
          CALL IWRTMA(IST_CA,1,NGAS,1,NGAS)
          CALL IWRTMA(IST_AA,1,NGAS,1,NGAS)
          CALL IWRTMA(IST_AA,1,NGAS,1,NGAS)
          CALL IWRTMA(IST_AB,1,NGAS,1,NGAS)
       END IF
*. CA mapping
      ONE = 1.0D0
      CALL MAP_EXSTR(IDE_SP_CA,IDE_CA_CA,LDE_CA,IST_CA,
     &               IMAP_CA,SMAP_CA,NK_CA,IB_CA,ONE)
*. CB mapping
      CALL MAP_EXSTR(IDE_SP_CB,IDE_CA_CB,LDE_CB,IST_CB,
     &               IMAP_CB,SMAP_CB,NK_CB,IB_CB,ONE)
*. AA mapping
*. ( annihilations strings are obtained as standard ordered 
*    creation strings conjugated. Conjugate DE operator
*. ( a sign is also needed I guess...)
C     CONJ_OP(IGAS_IN, ICA_IN,NOP,IGAS_OUT,ICA_OUT)
      CALL CONJ_OP(IDE_SP_AA,IDE_CA_AA,LDE_AA,IOP_AR,ICA_AR)
      CALL MAP_EXSTR(IOP_AR,ICA_AR,LDE_AA,IST_AA, 
     &               IMAP_AA,SMAP_AA,NK_AA,IB_AA,ONE)
*. AB mapping ( A sign should be added)
      CALL CONJ_OP(IDE_SP_AB,IDE_CA_AB,LDE_AB,IOP_AR,ICA_AR)
      CALL MAP_EXSTR(IOP_AR,ICA_AR,LDE_AB,IST_AB, 
     &               IMAP_AB,SMAP_AB,NK_AB,IB_AB,ONE)
*
      RETURN
      END
      SUBROUTINE ISMNM_FOR_TCC_ALL(NS_CAAB,ISM_CAAB,INM_CAAB,ISM,
     &                             NELMNT)     
*
* A Coupled cluster block is given with dimensions N_CAAB(ISM,I_CAAB) for 
* the four substrings.
* Obtain the symmetries and string numbers of the four substrings for  
* all strings in cc block                     
* 
*
* Jeppe Olsen, May of 2000  
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER NS_CAAB(8,4)             
*. Output
      INTEGER ISM_CAAB(4,*), INM_CAAB(4,*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) '  ISMNM_FOR_TCC_ALL speaking : '
        WRITE(6,*)
C       WRITE(6,*) ' NS_CAAB matrix : '
C       CALL IWRTMA(NS_CAAB,NSMST,4,8,4)
      END IF
*
      IT =  0
      DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        DO ISM_CA = 1, NSMST
          ISM_CB = MULTD2H(ISM_C,ISM_CA)
          DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*. Loop over operators as  matrix (I_CA,I_CB,I_AA,I_AB)
            NSTR_CA = NS_CAAB(ISM_CA,1)
            NSTR_CB = NS_CAAB(ISM_CB,2)
            NSTR_AA = NS_CAAB(ISM_AA,3)
            NSTR_AB = NS_CAAB(ISM_AB,4)

            DO I_AB = 1, NSTR_AB
             DO I_AA = 1, NSTR_AA
              DO I_CB = 1, NSTR_CB
               DO I_CA = 1, NSTR_CA
                IT = IT + 1
*
                  INM_CAAB(1,IT) = I_CA
                  INM_CAAB(2,IT) = I_CB
                  INM_CAAB(3,IT) = I_AA
                  INM_CAAB(4,IT) = I_AB
*
                  ISM_CAAB(1,IT) = ISM_CA
                  ISM_CAAB(2,IT) = ISM_CB
                  ISM_CAAB(3,IT) = ISM_AA
                  ISM_CAAB(4,IT) = ISM_AB
*
               END DO
              END DO
             END DO
            END DO
*           ^ End of loop over elements of block
          END DO
*         ^ End of loop over ISM_AA
        END DO
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
      NELMNT = IT
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' String numbers for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(INM_CAAB,4,IT,4,IT)
        WRITE(6,*)
        WRITE(6,*) ' Symmetry for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(ISM_CAAB,4,IT,4,IT)
        WRITE(6,*)
      END IF
*
      RETURN
      END
      SUBROUTINE ISMNM_FOR_TCC_BAT(NS_CAAB,ISM_CAAB,INM_CAAB,ISM,
     &           LBAT,ISM_INI,
     &           ISM_C1,ISM_CA1,ISM_AA1,
     &           INM_AB1,INM_AA1,INM_CA1,INM_CB1,
     &           ISM_CINI,ISM_CAINI,ISM_AAINI,
     &           INM_ABINI,INM_AAINI,INM_CAINI,INM_CBINI,IONLY_LIM)     
*
* A Coupled cluster block is given with dimensions NS_CAAB(ISM,I_CAAB) for 
* the four substrings.
* Obtain symmetries and string numbers of the four substrings for  
* a batch of strings with total sym ISM. 
* The last string in previous batch is defined 
* by  ISM_C1,ISM_CA1,ISM_AA1,
*     INM_AB1,INM_AA1,INM_CA1,INM_CB1
*
* IF IONLY_LIM = 1, only the limits are returned
*
* If ISM_INI = 1, this is first batch with this symmetry
* 
*
*  Input / Argument list 
* ======================
* NS_CAAB : Number of strings per sym for each CAAB part of string
* ISM_CAAB, INM_CAAB (OUTPUT) : List of symmetries and strings for batch 
* LBAT : Number of strings required
* ISM_INI : = 1 => Initial batch 
* ISM_C1, ISM_CA1, ISM_AA1 : Values of ISM_C, ISM_CA, ISM_AA for last string 
*                            in previous batch
* INM_AB1, INM_AA1, INM_CA1, INM_CB1 : Values of INM_AB, INM_AA, INM_AB, INM_CB
*                                      for last string in previous batch
*
*  ISM_CINI,ISM_CAINI,ISM_AAINI,INM_ABINI,INM_AAINI,INM_CAINI,INM_CBINI (OUTPUT)
*  Initial values of ISM_C, ISM_CA, ISM_AA, INM_AB,INM_AA,INM_CA,INM_CB   
* IONLY_LIM : = 1 => Obtain only the limits (*1, *INI ) for batch
*
*
* Jeppe Olsen, May of 2001  
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER NS_CAAB(8,4)             
*. Output
      INTEGER ISM_CAAB(4,*), INM_CAAB(4,*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) '  ISMNM_FOR_TCC_BAT speaking : '
        WRITE(6,*)
      END IF
*
      IF(ISM_INI.EQ.1) THEN
*. Initialize start indeces
       ISM_C1 = 1
       ISM_CA1 = 1
       ISM_AA1 = 1
*
       INM_AB1 = 1
       INM_AA1 = 1
       INM_CB1 = 1
       INM_CA1 = 0
      END IF
*. Save initial indeces
      ISM_CINI  =  ISM_C1
      ISM_CAINI =  ISM_CA1
      ISM_AAINI =  ISM_AA1 
*
      INM_ABINI =  INM_AB1
      INM_AAINI =  INM_AA1 
      INM_CBINI =  INM_CB1 
      INM_CAINI =  INM_CA1
*
      IT =  0
      INI_LOOP = 1
      DO ISM_C = ISM_C1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        IF(INI_LOOP .EQ. 1 ) THEN
          ISM_CA2 = ISM_CA1
        ELSE
          ISM_CA2 = 1
        END IF
        DO ISM_CA = ISM_CA2, NSMST
          ISM_CB = MULTD2H(ISM_C,ISM_CA)
          IF(INI_LOOP .EQ. 1 ) THEN
            ISM_AA2 = ISM_AA1
          ELSE
            ISM_AA2 = 1
          END IF
          DO ISM_AA = ISM_AA2, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*. Loop over operators as  matrix (I_CA,I_CB,I_AA,I_AB)
            NSTR_CA = NS_CAAB(ISM_CA,1)
            NSTR_CB = NS_CAAB(ISM_CB,2)
            NSTR_AA = NS_CAAB(ISM_AA,3)
            NSTR_AB = NS_CAAB(ISM_AB,4)
*
            IF(INI_LOOP.EQ.1) THEN
             INM_AB2 = INM_AB1
            ELSE
             INM_AB2 = 1
            END IF
            DO I_AB = INM_AB2, NSTR_AB
             IF(INI_LOOP.EQ.1) THEN
              INM_AA2 = INM_AA1
             ELSE
              INM_AA2 = 1
             END IF
             DO I_AA = INM_AA2, NSTR_AA
              IF(INI_LOOP.EQ.1) THEN
                INM_CB2 = INM_CB1
              ELSE 
                INM_CB2 = 1
              END IF
              DO I_CB = INM_CB2, NSTR_CB
               IF(INI_LOOP.EQ.1) THEN
                 INM_CA2 = INM_CA1
               ELSE 
                 INM_CA2 = 0
               END IF
               DO I_CA = INM_CA2 + 1, NSTR_CA
                IT = IT + 1
*
                IF(IONLY_LIM.EQ.0) THEN
                  INM_CAAB(1,IT) = I_CA
                  INM_CAAB(2,IT) = I_CB
                  INM_CAAB(3,IT) = I_AA
                  INM_CAAB(4,IT) = I_AB
*
                  ISM_CAAB(1,IT) = ISM_CA
                  ISM_CAAB(2,IT) = ISM_CB
                  ISM_CAAB(3,IT) = ISM_AA
                  ISM_CAAB(4,IT) = ISM_AB
                END IF
                  IF(IT.EQ.LBAT) THEN
*. Save last indeces for me
                    ISM_C1  = ISM_C
                    ISM_CA1 = ISM_CA
                    ISM_AA1 = ISM_AA 
*
                    INM_AB1 = I_AB 
                    INM_AA1 = I_AA
                    INM_CB1 = I_CB
                    INM_CA1 = I_CA
                    GOTO 1001
                   END IF
*                  ^ End if last element
               END DO
               INI_LOOP = 0
              END DO
              INI_LOOP = 0
             END DO
             INI_LOOP = 0
            END DO
            INI_LOOP = 0
*           ^ End of loop over elements of block
          END DO
          INI_LOOP = 0
*         ^ End of loop over ISM_AA
        END DO
        INI_LOOP = 0
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' String numbers for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(INM_CAAB,4,LBAT,4,LBAT)
        WRITE(6,*)
        WRITE(6,*) ' Symmetry for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(ISM_CAAB,4,LBAT,4,LBAT)
        WRITE(6,*)
      END IF
*
      RETURN
      END
      SUBROUTINE ISMNM_FOR_TCC(NS_CAAB,ISM_CAAB,INM_CAAB,ISM,
     &                        ISTART,ISTOP) 
*
* A Coupled cluster block is given with dimensions N_CAAB(ISM,I_CAAB) for 
* the four substrings.
* Obtain the symmetries and string numbers of the four substrings for  
* elements  ISTART to ISTOP of this block
* 
*
* Jeppe Olsen, May of 2000  
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER NS_CAAB(NGAS,4)             
*. Output
      INTEGER ISM_CAAB(*), INM_CAAB(*)
      IT = - ISTART + 1
      DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        DO ISM_CA = 1, NSMST
          ISM_CB = MULTD2H(ISM_C,ISM_CA)
          DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*. Loop over operators as  matrix (I_CA,I_CB,I_AA,I_AB)
            NSTR_CA = NS_CAAB(ISM_CA,1)
            NSTR_CB = NS_CAAB(ISM_CB,2)
            NSTR_AA = NS_CAAB(ISM_AA,3)
            NSTR_AB = NS_CAAB(ISM_AB,4)

            DO I_AB = 1, NSTR_AB
             DO I_AA = 1, NSTR_AA
              DO I_CB = 1, NSTR_CB
               DO I_CA = 1, NSTR_CA
                IT = IT + 1
                IF(IT.GE.1) THEN
*
                  INM_CAAB((IT-1)+1) = I_CA
                  INM_CAAB((IT-1)+2) = I_CB
                  INM_CAAB((IT-1)+3) = I_AA
                  INM_CAAB((IT-1)+4) = I_AB
*
                  ISM_CAAB((IT-1)+1) = ISM_CA
                  ISM_CAAB((IT-1)+2) = ISM_CB
                  ISM_CAAB((IT-1)+3) = ISM_AA
                  ISM_CAAB((IT-1)+4) = ISM_AB

                END IF
                IF(IT.EQ.ISTOP-ISTART+1) GOTO 1001
               END DO
              END DO
             END DO
            END DO
*           ^ End of loop over elements of block
          END DO
*         ^ End of loop over ISM_AA
        END DO
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
*
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' String numbers for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(INM_CAAB,IT,4,IT,4)
        WRITE(6,*)
        WRITE(6,*) ' Symmetry for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(ISM_CAAB,IT,4,IT,4)
        WRITE(6,*)
      END IF
*
      RETURN
      END
      FUNCTION LEN_TCCBLK(NCA,NCB,NAA,NAB,ITSYM,NSMST)
*
* Find length of Coupled cluster excitation block,
* given types, total symmetry is ITSYM
*
* Jeppe Olsen, May of 2000
* Changed from NSMST**3 to NSMST**2 procedure, July 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NCA(*),NCB(*),NAA(*),NAB(*)
*. Local scratch 
      INTEGER NC(8), NA(8)
*
*. Number of creation and annihilation strings per symmetry
      IZERO = 0
      CALL ISETVC(NC,IZERO,NSMST)
      CALL ISETVC(NA,IZERO,NSMST)
      DO IA_SM = 1, NSMST
        DO IB_SM = 1, NSMST
          IAB_SM = MULTD2H(IA_SM,IB_SM) 
          NC(IAB_SM) = NC(IAB_SM) + NCA(IA_SM)*NCB(IB_SM)
          NA(IAB_SM) = NA(IAB_SM) + NAA(IA_SM)*NAB(IB_SM)
        END DO
      END DO
*. And number of strings with sym ITSYM
      LEN = 0
      DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM_C,ITSYM) 
*. Number of C and A strings with this sym
        LEN = LEN + NC(ISM_C)*NA(ISM_A)
      END DO
*
      LEN_TCCBLK = LEN
*
      RETURN
      END
      SUBROUTINE OPCT1T2(IOEX,IO1DX,IO2DX,IT1,IT2,T1,T2,OT1T2,
     &           IT1SM,IT2SM,IOPSM,
     &           LD12B,LB,
     &           IOD2_ST,IOD1_ST,IOEX_ST,
     &           IX1_CA,SX1_CA,IX1_CB,SX1_CB,
     &           IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &           IX2_CA,SX2_CA,IX2_CB,SX2_CB,
     &           IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &           IX3_CA,SX3_CA,IX3_CB,SX3_CB,
     &           IX3_AA,SX3_AA,IX3_AB,SX3_AB,
     &           IX4_CA,SX4_CA,IX4_CB,SX4_CB,
     &           IX4_AA,SX4_AA,IX4_AB,SX4_AB,
     &           ISTR_D1,ISTR_D2,ISTR_EX, ISTR_K2,ISTR_L1,ISTR_K1,
     &           TSCR1,TSCR2,TSCR3,
     &           TSCR4,OPSCR,
     &           ISM_CAAB_D1,ISM_CAAB_D2,ISM_CAAB_K1,
     &           ISM_CAAB_K2,ISM_CAAB_EX,ISM_CAAB_L1,
     &           INM_CAAB_D1,INM_CAAB_D2,INM_CAAB_K1,
     &           INM_CAAB_K2,INM_CAAB_EX,INM_CAAB_L1,IEXD1D2_INDX,
     &           IOC_K1,IOC_K2,IOC_L1,IOC_OT1T2,IZ,IZSCR,ISTREO,
     &           ISIGNG,FACX,N_TDL_MAX)
*
* Contract indeces of Operator O with indeces of 
* excitation operator T1 and T2
*
* Operator O is defined by an operator part IOEX, and 
* two deexcitation parts IO1DX,IO2DX.
*
* LD12B : Batch size for D1 and D2
* LB2   : Batch size for for remaining expansions
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'ctcc.inc'
*. Input   
      INTEGER IOEX(NGAS,4),IO1DX(NGAS,4),IO2DX(NGAS,4)
      INTEGER IT1(NGAS,4),IT2(NGAS,4)
      DIMENSION T1(*),T2(*)
*. Index of EXD1D2 operator to original index of Hamiltonian
      INTEGER IEXD1D2_INDX(4)
*. Output
      DIMENSION OT1T2(*)
*
*. Local scratch
*
*. Occupation of conjugated operators
      INTEGER IO1DX_DAG(4*MXPNGAS),IO2DX_DAG(4*NGAS)
*. Occupation of OP T1 T2
      INTEGER IOC_OT1T2(NGAS,4)
*
*. Occupation of gasspaces for various strings 
COLD  INTEGER IOC_K1(MXPNGAS,4),IOC_K2(MXPNGAS,4)
COLD  INTEGER IOC_L1(MXPNGAS,4)
*. Number of strings per sym of the various string supergroups
      INTEGER NOEX(8,4),NO1DX(8,4),NO2DX(8,4)
      INTEGER NT1(8,4), NT2(8,4)
      INTEGER NK1(8,4),NK2(8,4),NL1(8,4)
      INTEGER NOT1T2(8,4)
*. Offsets of strings with given sym for the strings of various CC ops
      INTEGER IBOEX(8,4),IBO1DX(8,4),IBO2DX(8,4)
      INTEGER IBT1(8,4), IBT2(8,4)
      INTEGER IBK1(8,4),IBK2(8,4),IBL1(8,4)
      INTEGER IBOT1T2(8,4)
*. offsets to start of TCC strings with given symmetry
      INTEGER IBK1_T(8),IBK2_T(8),IBD1_T(8),IBD2_T(8) 
      INTEGER IBEX_T(8),IBL1_T(8)
*
      INTEGER IBK1_CA(8,8,MXPLCCOP), IBK1_CB(8,8,MXPLCCOP)
      INTEGER IBK1_AA(8,8,MXPLCCOP), IBK1_AB(8,8,MXPLCCOP)
      INTEGER NK1_CA(MXPOBS,MXPLCCOP), NK1_CB(MXPOBS,MXPLCCOP)
      INTEGER NK1_AA(MXPOBS,MXPLCCOP), NK1_AB(MXPOBS,MXPLCCOP)
*
      INTEGER IBK2_CA(8,8,MXPLCCOP), IBK2_CB(8,8,MXPLCCOP)
      INTEGER IBK2_AA(8,8,MXPLCCOP), IBK2_AB(8,8,MXPLCCOP)
      INTEGER NK2_CA(MXPOBS,MXPLCCOP), NK2_CB(MXPOBS,MXPLCCOP)
      INTEGER NK2_AA(MXPOBS,MXPLCCOP), NK2_AB(MXPOBS,MXPLCCOP)
*
      INTEGER IBK2B_CA(8,8,MXPLCCOP), IBK2B_CB(8,8,MXPLCCOP)
      INTEGER IBK2B_AA(8,8,MXPLCCOP), IBK2B_AB(8,8,MXPLCCOP)
      INTEGER NK2B_CA(MXPOBS,MXPLCCOP), NK2B_CB(MXPOBS,MXPLCCOP)
      INTEGER NK2B_AA(MXPOBS,MXPLCCOP), NK2B_AB(MXPOBS,MXPLCCOP)
*
      INTEGER IBL1_CA(8,8,MXPLCCOP), IBL1_CB(8,8,MXPLCCOP)
      INTEGER IBL1_AA(8,8,MXPLCCOP), IBL1_AB(8,8,MXPLCCOP)
      INTEGER NL1_CA(MXPOBS,MXPLCCOP), NL1_CB(MXPOBS,MXPLCCOP)
      INTEGER NL1_AA(MXPOBS,MXPLCCOP), NL1_AB(MXPOBS,MXPLCCOP)
*
      INTEGER I_CA_EXP(MXPLCCOP,2),I_CB_EXP(MXPLCCOP,2)
      INTEGER I_AA_EXP(MXPLCCOP,2),I_AB_EXP(MXPLCCOP,2)
*
C     INTEGER IAC_AR(MXPLCCOP),IOP_AR(MXPLCCOP)
*. Offset in operators to strings with given sym
      INTEGER IBT1_TCC(8,8,8), IBL1_TCC(8,8,8), IBT2_TCC(8,8,8)
      INTEGER IBOT1T2_TCC(8,8,8)
*
      INTEGER IB_D1K1(8,8,4), IB_D2K2(8,8,4)
      INTEGER IB_EXK1(8,8,4), IB_L1K2(8,8,4)
*. Number of CA, CB, AA. AB strings per sym
C     INTEGER NST_CAAB_D1(MXPNSMST*4),NST_CAAB_D2(MXPNSMST*4)
C     INTEGER NST_CAAB_K1(MXPNSMST*4),NST_CAAB_K2(MXPNSMST*4)
C     INTEGER NST_CAAB_EX(MXPNSMST*4),NST_CAAB_L1(MXPNSMST*4)
*
*. Scratch through parameter list. 
*
      INTEGER IOC_K1(NGAS,4),IOC_K2(NGAS,4), IOC_L1(NGAS,4)
      INTEGER IOD1_ST(4,*),IOD2_ST(4,*),IOEX_ST(4,*)
*     ^ Dimension : 4 * LD12B
*. IX1_* : Number of operators in operator * Largest C or A string
      DIMENSION IX1_CA(*),SX1_CA(*),IX1_CB(*),SX1_CB(*)
      DIMENSION IX1_AA(*),SX1_AA(*),IX1_AB(*),SX1_AB(*)
*. IX2, SX2
      DIMENSION IX2_CA(*),SX2_CA(*),IX2_CB(*),SX2_CB(*)
      DIMENSION IX2_AA(*),SX2_AA(*),IX2_AB(*),SX2_AB(*)
*. IX3, SX3
      DIMENSION IX3_CA(*),SX3_CA(*),IX3_CB(*),SX3_CB(*)
      DIMENSION IX3_AA(*),SX3_AA(*),IX3_AB(*),SX3_AB(*)
*. IX4, SX4
      DIMENSION IX4_CA(*),SX4_CA(*),IX4_CB(*),SX4_CB(*)
      DIMENSION IX4_AA(*),SX4_AA(*),IX4_AB(*),SX4_AB(*)
*. ISTR_D1, ISTR_K2, ISTR_L1 : for occupations of strings  
      INTEGER ISTR_D1(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D2(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_EX(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_K1(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_K2(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_L1(MX_ST_TSOSO_BLK_MX*NSMST,4)
*. For intermediates with both strings batched 
      DIMENSION TSCR1(LD12B*LD12B),TSCR2(LD12B*LD12B)
      DIMENSION TSCR3(LD12B*LD12B)
*. For an intermediate with only one string batched 
      DIMENSION TSCR4(N_TDL_MAX)
*.    ^ This dimension is probably too large !!
*. For a batch of coefficients for Operator
*. Operator accessed as OP(EX,D2,D1) so
      DIMENSION OPSCR(LD12B*LD12B*LB)
*. Number of CA, CB, AA. AB strings per sym
C     INTEGER NST_CAAB_D1(MXPNSMST*4),NST_CAAB_D2(MXPNSMST*4)
C     INTEGER NST_CAAB_K1(MXPNSMST*4),NST_CAAB_K2(MXPNSMST*4)
C     INTEGER NST_CAAB_EX(MXPNSMST*4),NST_CAAB_L1(MXPNSMST*4)
*. For part of Hamiltonian
      INTEGER ISM_CAAB_D1(4,*)
      INTEGER INM_CAAB_D1(4,*)
      INTEGER ISM_CAAB_D2(4,*) 
      INTEGER INM_CAAB_D2(4,*)
      INTEGER ISM_CAAB_EX(4,*) 
      INTEGER INM_CAAB_EX(4,*)
*. For CC operators
      INTEGER ISM_CAAB_K1(4,*) 
      INTEGER INM_CAAB_K1(4,*)
      INTEGER ISM_CAAB_K2(4,*) 
      INTEGER INM_CAAB_K2(4,*)
      INTEGER ISM_CAAB_L1(4,*) 
      INTEGER INM_CAAB_L1(4,*)
*
      INTEGER IZ(*), IZSCR(*), ISTREO(*)
*. For testing generation of strings in batch
      INTEGER IISM_CAAB_D1(4,100),IINM_CAAB_D1(4,100)
*. Collecting info on largest T(D,L) Block
      COMMON/CNTDL_MAX/NTDL_MAX_ACT,ID2_MAX_ACT(MXPNGAS*4),
     &      IL1_MAX_ACT(MXPNGAS*4),IEX_MAX_ACT(MXPNGAS*4),
     &      IT1_MAX_ACT(MXPNGAS*4),IT2_MAX_ACT(MXPNGAS*4)
*
* The story goes as
*
* Loop over sym of D2 
*  Loop over batches of D2
*   Loop over symmetry of EX => Symmetry of D1 => Sym of K1 
*    Loop over batches of Ex
*     Loop over batches of D1
*      Loop over batches of resolution strings K1 for T1
*       Obtain T1(d1,k1) = T1(d1,Kca1,Kcb1,Kaa1,Kaa2) 
*       Obtain O(Kex,d2,d1)
*       OT1(Kex,Kd2,K1) = 
*       Sum(Kd1) O(Kex,d2,Kd1)*T(d1,K1)
*       Expand OT1(Kex,d2,K1) to OT1(d2,L1)
*      End of loop over batches of K1
*     End of loop over batches of D1
*    End of loop over batches of Ex
*   End of loop over sym of Ex  
*   We have now OT1(D2,L1) for all strings L1
*   Loop over batch of K2
*    Loop over Batches of L1
*     Obtain T2(D2,K2)
*     OT1T2(L1,K2) = sum(d2) T2(D2,K2)*OT1(D2,L1)
*     Reform OT1T2(L1,K2) to OT1T2(I)    
*    End of loop over batches of L1
*   End of loop over batches of K2
*  End of loop over batches of D2
* End of loop over symmetries of D2
*
* - That's all she wrote
*
* The batches over indeces that are summation indeces in 
* matrix multiplications, that is D1 and D2 should be as 
* large as possible, preferable 100 - 1000.
* The size of the remaining indeces is less important   
* although they combined should be large enough to 
* make the matrix multiplication efficient.
*
      IDUM = 0
*. New or old TOKB routine
      I_NEW_OR_OLD_REFORM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'OPCTT ')
      NTEST = 00 
      IF(NTEST.GE.100 )THEN
        WRITE(6,*) ' OPCT1T2 entered '
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Excitation op, contraction op1, contraction op2 '
        WRITE(6,*)
        CALL WRT_SPOX_TP(IOEX,1)
        CALL WRT_SPOX_TP(IO1DX,1)
        CALL WRT_SPOX_TP(IO2DX,1)
        WRITE(6,*) 
        WRITE(6,*) ' T1 and T2 '
        CALL WRT_SPOX_TP(IT1,1)
        CALL WRT_SPOX_TP(IT2,1)
*
      END IF
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' N_TDL_MAX in OPCT1T2 ', N_TDL_MAX
      SIGNG = DBLE(ISIGNG)
*. Number of operators in IT2 - simplifications if this is unit operator
      NOP_T2 = IELSUM(IT2,4*NGAS)
*. Symmetry of OT1T2 
      IT1T2SM = MULTD2H(IT1SM,IT2SM)
      IOT1T2SM = MULTD2H(IOPSM,IT1T2SM)
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Symmetry of OT1T2 = ', IOT1T2SM
      END IF
*
* =================================================
* Occupation of the various gas spaces for strings 
* =================================================
*
*. K1 strings 
*. K1 = O1DX  I1 
C          NEW_CAAB_OC(IOCC_L,IOCC_R,IOCC_OP,ICE,ILR,NGAS)
      CALL NEW_CAAB_OC(IOC_K1,IT1,IO1DX,1,1,NGAS)
*. Occupation of L1 strings
*. L1 = OEX K1 
      CALL NEW_CAAB_OC(IOC_L1,IOC_K1,IOEX,2,1,NGAS)
*. Occupation of K2 strings
*. K2 = O2DX I2
      CALL NEW_CAAB_OC(IOC_K2,IT2,IO2DX,1,1,NGAS)
*. Occupation of Op T1 T2 = L1 K2
      CALL NEW_CAAB_OC(IOC_OT1T2,IOC_K2,IOC_L1,2,1,NGAS)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation of K1, L1 and K2 and OT1T2 ' 
        CALL WRT_SPOX_TP(IOC_K1,1)
        CALL WRT_SPOX_TP(IOC_L1,1)
        CALL WRT_SPOX_TP(IOC_K2,1)
        CALL WRT_SPOX_TP(IOC_OT1T2,1)
      END IF
*
* =============================================================
*. Obtain symmetry-dimensions and -offsets for various strings
* =============================================================
*
      DO I_CAAB = 1, 4
*. OEX 
       CALL NST_SPGP(IOEX(1,I_CAAB),NOEX(1,I_CAAB))
       CALL ZBASE(NOEX(1,I_CAAB),IBOEX(1,I_CAAB),NSMST)
C      ZBASE(NVEC,IVEC,NCLASS)
*. D1
       CALL NST_SPGP(IO1DX(1,I_CAAB),NO1DX(1,I_CAAB))
       CALL ZBASE(NO1DX(1,I_CAAB),IBO1DX(1,I_CAAB),NSMST)
*. D2
       CALL NST_SPGP(IO2DX(1,I_CAAB),NO2DX(1,I_CAAB))
       CALL ZBASE(NO2DX(1,I_CAAB),IBO2DX(1,I_CAAB),NSMST)
*. T1
       CALL NST_SPGP(IT1(1,I_CAAB),NT1(1,I_CAAB))
       CALL ZBASE(NT1(1,I_CAAB),IBT1(1,I_CAAB),NSMST)
*. T2
       CALL NST_SPGP(IT2(1,I_CAAB),NT2(1,I_CAAB))
       CALL ZBASE(NT2(1,I_CAAB),IBT2(1,I_CAAB),NSMST)
*. K1
       CALL NST_SPGP(IOC_K1(1,I_CAAB),NK1(1,I_CAAB))
C?     WRITE(6,*) ' NK1 for I_CAAB = ', I_CAAB 
C?     CALL IWRTMA(NK1(1,I_CAAB),1,NSMST,1,NSMST)
       CALL ZBASE(NK1(1,I_CAAB),IBK1(1,I_CAAB),NSMST)
*. K2
       CALL NST_SPGP(IOC_K2(1,I_CAAB),NK2(1,I_CAAB))
       CALL ZBASE(NK2(1,I_CAAB),IBK2(1,I_CAAB),NSMST)
*. L1
       CALL NST_SPGP(IOC_L1(1,I_CAAB),NL1(1,I_CAAB))
       CALL ZBASE(NL1(1,I_CAAB),IBL1(1,I_CAAB),NSMST)
*. OT1T2
       CALL NST_SPGP(IOC_OT1T2(1,I_CAAB),NOT1T2(1,I_CAAB))
       CALL ZBASE(NOT1T2(1,I_CAAB),IBOT1T2(1,I_CAAB),NSMST)
*
      END DO
C?    WRITE(6,*) ' Fresh NOEX array '
C?    CALL IWRTMA(NOEX,8,4,8,4)
*. We now have the various dimensions, so we can write T1 and T2 if 
*. Required
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NT1 :     '
        CALL IWRTMA(NT1,NSMST,4,8,4) 
        WRITE(6,*) ' IT1SM = ', IT1SM
        WRITE(6,*) ' Input T1 block '
C            WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
        CALL WRT_TCC_BLK(T1,IT1SM,NT1(1,1),NT1(1,2),NT1(1,3),NT1(1,4),
     &                   NSMST)
        WRITE(6,*) ' Input T2 block '
        CALL WRT_TCC_BLK(T2,IT2SM,NT2(1,1),NT2(1,2),NT2(1,3),NT2(1,4),
     &                   NSMST)
      END IF
*
*
* =================================================
*. Offsets to T1(ICA,ICB,IAA,IAB)( and T2 OT1T2)  with given sym
* ==============================================================
*
*. T1
      CALL Z_TCC_OFF2(IBT1_TCC,LEN_T1,NT1(1,1),NT1(1,2),NT1(1,3),
     &                NT1(1,4),IT1SM,NSMST)
*. T2
      CALL Z_TCC_OFF2(IBT2_TCC,LEN_T2,NT2(1,1),NT2(1,2),NT2(1,3),
     &                NT2(1,4),IT2SM,NSMST)
*. OT1T2
      CALL Z_TCC_OFF2(IBOT1T2_TCC,LEN_OT1T2,NOT1T2(1,1),
     &     NOT1T2(1,2),NOT1T2(1,3),NOT1T2(1,4),IOT1T2SM,NSMST)
*
*
* ========================
*.  D1{\dagger} K1  => T1 mapping
* =========================
*
*. Obtain D1{\dagger}
      CALL CONJ_CAAB(IO1DX,IO1DX_DAG,NGAS,SP_D1K1)
      CALL T1T2_TO_T12_MAP(IO1DX_DAG,IOC_K1,IT1,
     &     IX1_CA,SX1_CA,IX1_CB,SX1_CB,IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &     IB_D1K1,ISTR_D1,ISTR_K1,ISTREO,IZ,IZSCR,SIGN_D1K1)
*
* ======================
* Ex K1 => L1 mapping
* ======================
*
      CALL T1T2_TO_T12_MAP(IOEX,IOC_K1,IOC_L1,
     &     IX2_CA,SX2_CA,IX2_CB,SX2_CB,IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &     IB_EXK1,ISTR_K1,ISTR_L1,ISTREO,IZ,IZSCR,SIGN_EXK1)
*
* ==============
*  DX2 K2 => I2
* ==============
*
* ========================
*.  D2{\dagger} K2  => T2 mapping
* =========================
*
*. Obtain D2{\dagger}
      CALL CONJ_CAAB(IO2DX,IO2DX_DAG,NGAS,SP_D2K2)
      CALL T1T2_TO_T12_MAP(IO2DX_DAG,IOC_K2,IT2,
     &     IX3_CA,SX3_CA,IX3_CB,SX3_CB,IX3_AA,SX3_AA,IX3_AB,SX3_AB,
     &     IB_D2K2,ISTR_D2,ISTR_K2,ISTREO,IZ,IZSCR,SIGN_D2K2)
*
* ==========
* L1 K2 => I
* ==========
*
      CALL T1T2_TO_T12_MAP(IOC_L1,IOC_K2,IOC_OT1T2,
     &     IX4_CA,SX4_CA,IX4_CB,SX4_CB,IX4_AA,SX4_AA,IX4_AB,SX4_AB,
     &     IB_L1K2,ISTR_L1,ISTR_K2,ISTREO,IZ,IZSCR,SIGN_L1K2)
*. And the individual strings : all symmetries constructed
C     STR_CAAB(ICAAB,ISTR_CAAB)
*. D1 strings 
      CALL STR_CAAB(IO1DX,ISTR_D1)
*. D2 strings 
      CALL STR_CAAB(IO2DX,ISTR_D2)
*. Ex strings 
      CALL STR_CAAB(IOEX,ISTR_EX)
      DO ID2SM = 1, NSMST
*. Symmetries of K2 and L1 is now well-defined
       IOPT1SM = MULTD2H(IOPSM,IT1SM)
       L1SM = MULTD2H(ID2SM,IOPT1SM)
       IL1SM_NEW = 1
       K2SM = MULTD2H(IT2SM,ID2SM)
       IK2SM_NEW = 1
       IF(NTEST.GE.100) 
     & WRITE(6,*) ' ID2SM, L1SM, K2SM = ',ID2SM, L1SM, K2SM
*. Number of D2 strings with actual symmetry
       ND2_TOT = LEN_TCCBLK(NO2DX(1,1),NO2DX(1,2),NO2DX(1,3),
     &               NO2DX(1,4),ID2SM,NSMST)
*. Number of K2 strings with actual symmetry
       NK2_TOT = LEN_TCCBLK(NK2(1,1),NK2(1,2),NK2(1,3),
     &               NK2(1,4),K2SM,NSMST)
*. Number of L1 strings with actual symmetry
       NL1_TOT = LEN_TCCBLK(NL1(1,1),NL1(1,2),NL1(1,3),
     &               NL1(1,4),L1SM,NSMST)
       IF(NL1_TOT*NK2_TOT*ND2_TOT.NE.0) THEN
*. Offset to given symsubblocks of L1
      ITDIAG = 0
      CALL Z_TCC_OFF(IBL1_TCC,NL1(1,1),NL1(1,2),NL1(1,3),NL1(1,4),
     &               L1SM,NSMST,ITDIAG)
*. Number of D2 batches 
       ND2_BAT = ND2_TOT/LD12B
       IF(ND2_BAT*LD12B.LT.ND2_TOT) ND2_BAT=ND2_BAT+1
*. And loop over batches of D2
       ID2SM_NEW = 1
       DO ID2_BAT = 1, ND2_BAT
        IF(NTEST.GE.100) WRITE(6,*) ' ID2_BAT = ', ID2_BAT
        ID2_START = (ID2_BAT-1)*LD12B + 1
        ID2_STOP  = MIN(ND2_TOT,ID2_START+LD12B-1)
        ID2_BATLEN = ID2_STOP-ID2_START+1
*. Generate D2 strings for given sym and batch
       CALL ISMNM_FOR_TCC_BAT(NO2DX,ISM_CAAB_D2,INM_CAAB_D2,ID2SM,
     &           ID2_BATLEN,ID2SM_NEW,
     &           ISM_C1_D2,ISM_CA1_D2,ISM_AA1_D2,
     &           INM_AB1_D2,INM_AA1_D2,INM_CA1_D2,INM_CB1_D2,
     &           ISM_CINI_D2,ISM_CAINI_D2,ISM_AAINI_D2,      
     &           INM_ABINI_D2,INM_AAINI_D2,INM_CAINI_D2,INM_CBINI_D2,
     &           0)
            ID2SM_NEW = 0
            ID2_OFF = 1
*. We will in the following lines construct OT1(D2,L1) for all L1 and
*. D2 in batch, TSCR4 is used for OT1 so clear
        ZERO = 0.0D0
*
        IF(ID2_BATLEN*NL1_TOT.GT.NTDL_MAX_ACT) THEN
          NTDL_MAX_ACT = ID2_BATLEN*NL1_TOT
          CALL ICOPVE(IOC_L1,IL1_MAX_ACT,4*NGAS)
          CALL ICOPVE(IOEX,IEX_MAX_ACT,4*NGAS)
          CALL ICOPVE(IO2DX,ID2_MAX_ACT,4*NGAS)
          CALL ICOPVE(IT1,IT1_MAX_ACT,4*NGAS)
          CALL ICOPVE(IT2,IT2_MAX_ACT,4*NGAS)
        END IF
*
        IF(ID2_BATLEN*NL1_TOT .GT. N_TDL_MAX ) THEN
*. Problem, too small length for TSCR4 
          WRITE(6,*) ' Dimension of TSCR4 too small '
          WRITE(6,*) ' Required and allocated : ', 
     &    ID2_BATLEN*NL1_TOT, N_TDL_MAX 
          WRITE(6,*) '  ID2_BATLEN, NL1_TOT = ',
     &                  ID2_BATLEN, NL1_TOT
          WRITE(6,*) ' D2 and L1 : '
          CALL WRT_SPOX_TP(IO2DX,1)
          CALL WRT_SPOX_TP(IOC_L1,1)
*
          WRITE(6,*) ' Excitation op, contraction op1, contraction op2'
          CALL WRT_SPOX_TP(IOEX,1)
          CALL WRT_SPOX_TP(IO1DX,1)
          CALL WRT_SPOX_TP(IO2DX,1)
          WRITE(6,*) ' T1 and T2 '
          CALL WRT_SPOX_TP(IT1,1)
          CALL WRT_SPOX_TP(IT2,1)
          STOP ' Dimension of TSCR4 too small '
        END IF
        CALL SETVEC(TSCR4,ZERO,ID2_BATLEN*NL1_TOT)
*
        DO IEXSM = 1, NSMST
         IEXSM_NEW = 1
         IEXD2SM = MULTD2H(ID2SM,IEXSM)
         ID1SM   = MULTD2H(IOPSM,IEXD2SM)
         ID1SM_NEW = 1
*  I1 = Kd1 K1
         K1SM = MULTD2H(ID1SM,IT1SM)
         IK1SM_NEW = 1
         IF(NTEST.GE.100) 
     & WRITE(6,*) ' K1SM, IEXSM =', K1SM, IEXSM
*. Total number of K1 strings with this symmetry
         NK1_TOT = LEN_TCCBLK(NK1(1,1),NK1(1,2),NK1(1,3),
     &         NK1(1,4),K1SM,NSMST)
         IF(NTEST.GE.100) WRITE(6,*) ' NK1_TOT = ', NK1_TOT
*. Number of excitations strings in O with this symmetry
         NEX_TOT = LEN_TCCBLK(NOEX(1,1),NOEX(1,2),NOEX(1,3),NOEX(1,4),
     &                   IEXSM,NSMST)
*. Number of batches of excitation part of O
         NEX_BAT = NEX_TOT/LB
         IF(NEX_BAT*LB.LT.NEX_TOT) NEX_BAT = NEX_BAT + 1
*. And loop over these batches
         IEXSM_NEW = 1
         DO IEX_BAT = 1, NEX_BAT
          IF(NTEST.GE.100) WRITE(6,*) ' IEX_BAT = ', IEX_BAT
          IEX_START = (IEX_BAT-1)*LB + 1
          IEX_STOP  = MIN(NEX_TOT,IEX_START+LB-1)
          IEX_BATLEN = IEX_STOP - IEX_START + 1
*. Generate IEX strings for given sym and batch
       CALL ISMNM_FOR_TCC_BAT(NOEX,ISM_CAAB_EX,INM_CAAB_EX,IEXSM,
     &           IEX_BATLEN,IEXSM_NEW,
     &           ISM_C1_EX,ISM_CA1_EX,ISM_AA1_EX,
     &           INM_AB1_EX,INM_AA1_EX,INM_CA1_EX,INM_CB1_EX,
     &           ISM_CINI_EX,ISM_CAINI_EX,ISM_AAINI_EX,      
     &           INM_ABINI_EX,INM_AAINI_EX,INM_CAINI_EX,INM_CBINI_EX,
     &           0)
            IEXSM_NEW = 0
            IEX_OFF = 1
*. Number of D1 strings with given symmetry
          ND1_TOT = LEN_TCCBLK(NO1DX(1,1),NO1DX(1,2),NO1DX(1,3),
     &          NO1DX(1,4),ID1SM,NSMST)
*. And number of batches
          ND1_BAT = ND1_TOT/LD12B
          IF(ND1_BAT*LD12B.LT.ND1_TOT) ND1_BAT=ND1_BAT+1
C?        WRITE(6,*) ' ND1_BAT = ', ND1_BAT
*. Loop over D1 batches
          ID1SM_NEW = 1
          DO ID1_BAT = 1, ND1_BAT
           IF(NTEST.GE.100) WRITE(6,*) ' ID1_BAT = ', ID1_BAT
           ID1_START = (ID1_BAT-1)*LD12B + 1
           ID1_STOP  = MIN(ND1_TOT,ID1_START+LD12B-1)
           ID1_BATLEN = ID1_STOP-ID1_START+1
*. Generate ID1 strings for given sym and batch
       CALL ISMNM_FOR_TCC_BAT(NO1DX,ISM_CAAB_D1,INM_CAAB_D1,ID1SM,
     &           ID1_BATLEN,ID1SM_NEW,
     &           ISM_C1_D1,ISM_CA1_D1,ISM_AA1_D1,
     &           INM_AB1_D1,INM_AA1_D1,INM_CA1_D1,INM_CB1_D1,
     &           ISM_CINI_D1,ISM_CAINI_D1,ISM_AAINI_D1,      
     &           INM_ABINI_D1,INM_AAINI_D1,INM_CAINI_D1,INM_CBINI_D1,
     &           0)
            ID1SM_NEW = 0
            ID1_OFF = 1
*. Obtain integrals OP(ID2,IEX,ID1)
      IF(NK1_TOT.GT.0) 
     &  CALL GET_OPINT(OPSCR,IO2DX,ID2_BATLEN,
     &  INM_CAAB_D2(1,ID2_OFF),ISM_CAAB_D2(1,ID2_OFF),
     &  ISTR_D2,IBO2DX,IO1DX,ID1_BATLEN,
     &  INM_CAAB_D1(1,ID1_OFF),ISM_CAAB_D1(1,ID1_OFF),
     &  ISTR_D1,IBO1DX,IOEX,IEX_BATLEN,
     &  INM_CAAB_EX(1,IEX_OFF),ISM_CAAB_EX(1,IEX_OFF),
     &  ISTR_EX,IBOEX,IEXD1D2_INDX,FACX)
*  I1 = Kd1 K1
*. And number of batches
           NK1_BAT = NK1_TOT/LB
           IF(NK1_BAT*LB.LT.NK1_TOT) NK1_BAT=NK1_BAT+1
*. Loop over K1 batches
           IK1SM_NEW = 1
           DO IK1_BAT = 1, NK1_BAT
           IF(NTEST.GE.100) WRITE(6,*) ' IK1_BAT = ', IK1_BAT
            IK1_START = (IK1_BAT-1)*LB + 1
            IK1_STOP  = MIN(NK1_TOT,IK1_START+LB-1)
            IK1_BATLEN = IK1_STOP-IK1_START+1 
*. Generate K1 strings for given sym and batch
       CALL ISMNM_FOR_TCC_BAT(NK1,ISM_CAAB_K1,INM_CAAB_K1,K1SM,
     &           IK1_BATLEN,IK1SM_NEW,
     &           ISM_C1_K1,ISM_CA1_K1,ISM_AA1_K1,
     &           INM_AB1_K1,INM_AA1_K1,INM_CA1_K1,INM_CB1_K1,
     &           ISM_CINI_K1,ISM_CAINI_K1,ISM_AAINI_K1,      
     &           INM_ABINI_K1,INM_AAINI_K1,INM_CAINI_K1,INM_CBINI_K1,
     &           1)
            IK1SM_NEW = 0
*. Obtain T(ID1,K1) in TSCR1 for operators in batches
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' T1(I) => T1(ID1,K1) reordering '
            END IF
            SIGNX = SP_D1K1*SIGN_D1K1 
*. Offsets for operators : relative adresses employed for operators, 
* so offsets are 1
            IK1_OFF = 1
*
C?          WRITE(6,*) ' New reformer '
            CALL LUCIAQENTER('TOK_A')
            CALL TI_TO_TOKBN(NSMST,
     &           ISM_CINI_D1,ISM_CAINI_D1,ISM_AAINI_D1,ID1SM,
     &           INM_ABINI_D1,INM_AAINI_D1,INM_CBINI_D1,INM_CAINI_D1,
     &           ID1_BATLEN,NO1DX,
     &           ISM_CINI_K1,ISM_CAINI_K1,ISM_AAINI_K1,K1SM,
     &           INM_ABINI_K1,INM_AAINI_K1,INM_CBINI_K1,INM_CAINI_K1,
     &           IK1_BATLEN,NK1,
     &           IX1_CA,SX1_CA,IB_D1K1(1,1,1),
     &           IX1_CB,SX1_CB,IB_D1K1(1,1,2),
     &           IX1_AA,SX1_AA,IB_D1K1(1,1,3),
     &           IX1_AB,SX1_AB,IB_D1K1(1,1,4),
     &           NT1(1,1),NT1(1,2),NT1(1,3),NT1(1,4),
     &           TSCR1,T1,IBT1_TCC,1,1,SIGNX,1)                  
            CALL LUCIAQEXIT('TOK_A')
*. Matrix multiply OT1(D2,EX,K1) = Sum(D1) OP(D2,EX,D1)T1(D1,K1)
*. result in TSCR2
            FACTORC = 0.0D0
            FACTORAB = 1.0D0
            NR_OT1 = ID2_BATLEN*IEX_BATLEN
            NC_OT1 = IK1_BATLEN
*
C           WRITE(6,*) ' OP*T1 => OT1 '
            CALL MATML7(TSCR2,OPSCR,TSCR1,NR_OT1,NC_OT1,NR_OT1,
     &           ID1_BATLEN,ID1_BATLEN,NC_OT1,FACTORC,FACTORAB,0)
*. Expand OT1(D2,EX,K1) to OT1(D2,L1) ( in TSCR4)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) '  OT1(D2,EX,K1) =>  OT1(D2,L1) reordering'
            END IF
            CALL LUCIAQENTER('TOK_B')
            CALL TI_TO_TOKBN(NSMST,
     &           ISM_CINI_EX,ISM_CAINI_EX,ISM_AAINI_EX,IEXSM,
     &           INM_ABINI_EX,INM_AAINI_EX,INM_CBINI_EX,INM_CAINI_EX,
     &           IEX_BATLEN,NOEX,
     &           ISM_CINI_K1,ISM_CAINI_K1,ISM_AAINI_K1,K1SM,
     &           INM_ABINI_K1,INM_AAINI_K1,INM_CBINI_K1,INM_CAINI_K1,
     &           IK1_BATLEN,NK1,
     &           IX2_CA,SX2_CA,IB_EXK1(1,1,1),
     &           IX2_CB,SX2_CB,IB_EXK1(1,1,2),
     &           IX2_AA,SX2_AA,IB_EXK1(1,1,3),
     &           IX2_AB,SX2_AB,IB_EXK1(1,1,4),
     &           NL1(1,1),NL1(1,2),NL1(1,3),NL1(1,4),
     &           TSCR2,TSCR4,IBL1_TCC,2,ID2_BATLEN,SIGN_EXK1,0)
            CALL LUCIAQEXIT('TOK_B')
           END DO
*          ^ End of loops over K1 batches
          END DO
*         ^ End of loop over D1 batch
         END DO
*        ^ End of loop over batches of excitation operators OEX
        END DO
*       ^ End of loop over symmetries of OEX
* ^ We have now OT2(d2,L1) for D2 in batch and all L1 of given symmetry
        IF(NOP_T2.EQ.0) THEN
*. Second operator is the unitoperator and L1 = I, so just add
*. TSCR4 to OPT1T2
          ONE = 1.0D0
          CALL VECSUM(OT1T2,OT1T2,TSCR4,ONE,SIGNG,NL1_TOT)
        ELSE
*. T2 is a nontrivial operator so ..
*. 
*. Number of batches of K2
         NK2_BAT = NK2_TOT/LB
         IF(NK2_BAT*LB.LT.NK2_TOT) NK2_BAT = NK2_BAT + 1
*. Number of batches of L1
         NL1_BAT = NL1_TOT/LB
         IF(NL1_BAT*LB.LT.NL1_TOT) NL1_BAT = NL1_BAT + 1
*. And loop over batches of K2
         IK2SM_NEW = 1
         DO IK2_BAT = 1, NK2_BAT
          IK2_START = (IK2_BAT-1)*LB + 1
          IK2_STOP  = MIN(NK2_TOT,IK2_START+LB-1)
          IK2_BATLEN = IK2_STOP-IK2_START+1 
*. Generate K2 strings for given sym and batch
       CALL ISMNM_FOR_TCC_BAT(NK2,ISM_CAAB_K2,INM_CAAB_K2,K2SM,
     &           IK2_BATLEN,IK2SM_NEW,
     &           ISM_C1_K2,ISM_CA1_K2,ISM_AA1_K2,
     &           INM_AB1_K2,INM_AA1_K2,INM_CA1_K2,INM_CB1_K2,     
     &           ISM_CINI_K2,ISM_CAINI_K2,ISM_AAINI_K2,      
     &           INM_ABINI_K2,INM_AAINI_K2,INM_CAINI_K2,INM_CBINI_K2,
     &           1)
            IK2SM_NEW = 0
*. Obtain T2(D2,K2) in TSCR1
            IF(NTEST.GE.100) THEN
              WRITE(6,*) '  T2(I2) => T2(D2,K2) reordering'
            END IF
          SIGNX = SIGN_D2K2*SP_D2K2
          IK2_OFF = 1
            CALL LUCIAQENTER('TOK_C')
            CALL TI_TO_TOKBN(NSMST,
     &           ISM_CINI_D2,ISM_CAINI_D2,ISM_AAINI_D2,ID2SM,
     &           INM_ABINI_D2,INM_AAINI_D2,INM_CBINI_D2,INM_CAINI_D2,
     &           ID2_BATLEN,NO2DX,
     &           ISM_CINI_K2,ISM_CAINI_K2,ISM_AAINI_K2,K2SM,
     &           INM_ABINI_K2,INM_AAINI_K2,INM_CBINI_K2,INM_CAINI_K2,
     &           IK2_BATLEN,NK2,
     &           IX3_CA,SX3_CA,IB_D2K2(1,1,1),
     &           IX3_CB,SX3_CB,IB_D2K2(1,1,2),
     &           IX3_AA,SX3_AA,IB_D2K2(1,1,3),
     &           IX3_AB,SX3_AB,IB_D2K2(1,1,4),
     &           NT2(1,1),NT2(1,2),NT2(1,3),NT2(1,4),
     &           TSCR1,T2,IBT2_TCC,1,1,SIGNX,1)
            CALL LUCIAQEXIT('TOK_C')
*. And loop over batches of L1
          IL1SM_NEW = 1
          DO IL1_BAT = 1, NL1_BAT
           IL1_START = (IL1_BAT-1)*LB + 1
           IL1_STOP  = MIN(NL1_TOT,IL1_START+LB-1)
           IL1_BATLEN = IL1_STOP-IL1_START+1 
*. Generate L1 strings for given sym and batch
       CALL ISMNM_FOR_TCC_BAT(NL1,ISM_CAAB_L1,INM_CAAB_L1,L1SM,
     &           IL1_BATLEN,IL1SM_NEW,
     &           ISM_C1_L1,ISM_CA1_L1,ISM_AA1_L1,
     &           INM_AB1_L1,INM_AA1_L1,INM_CA1_L1,INM_CB1_L1,     
     &           ISM_CINI_L1,ISM_CAINI_L1,ISM_AAINI_L1,      
     &           INM_ABINI_L1,INM_AAINI_L1,INM_CAINI_L1,INM_CBINI_L1,
     &           1)
            IL1SM_NEW = 0
*. Obtain OT1(D2,L1) for L1 in batch in TSCR2 
            CALL COPVEC(TSCR4((IL1_START-1)*ID2_BATLEN+1),TSCR2,
     &                  ID2_BATLEN*IL1_BATLEN) 
*. Sum OT1T2(L1,K2) = sum(D2) OT1(D2,L1)*T2(D2,K2) in TSCR3
*. Signg is introduced here !!
C?         WRITE(6,*) ' OT1*T2 => OT1T2 '
           CALL MATML7(TSCR3,TSCR2,TSCR1,IL1_BATLEN,IK2_BATLEN,
     &                 ID2_BATLEN,IL1_BATLEN,ID2_BATLEN,IK2_BATLEN,
     &                 FACTORC, SIGNG,1)
*. Expand OT1T2(L1,K2) to OT1T2(I) 
            IF(NTEST.GE.100) THEN
              WRITE(6,*) '   OT1T2(L1,K2) =>  OT1T2(I) reordering '
            END IF
C         IL1_OFF = IBL1_T(L1SM)-1+IL1_START
          IL1_OFF = 1
            CALL LUCIAQENTER('TOK_D')
            CALL TI_TO_TOKBN(NSMST,
     &           ISM_CINI_L1,ISM_CAINI_L1,ISM_AAINI_L1,L1SM,
     &           INM_ABINI_L1,INM_AAINI_L1,INM_CBINI_L1,INM_CAINI_L1,
     &           IL1_BATLEN,NL1,
     &           ISM_CINI_K2,ISM_CAINI_K2,ISM_AAINI_K2,K2SM,
     &           INM_ABINI_K2,INM_AAINI_K2,INM_CBINI_K2,INM_CAINI_K2,
     &           IK2_BATLEN,NK2,
     &           IX4_CA,SX4_CA,IB_L1K2(1,1,1),
     &           IX4_CB,SX4_CB,IB_L1K2(1,1,2),
     &           IX4_AA,SX4_AA,IB_L1K2(1,1,3),
     &           IX4_AB,SX4_AB,IB_L1K2(1,1,4),
     &           NOT1T2(1,1),NOT1T2(1,2),NOT1T2(1,3),NOT1T2(1,4),
     &           TSCR3,OT1T2,IBOT1T2_TCC,2,1,SIGN_L1K2,0)
            CALL LUCIAQEXIT('TOK_D')
          END DO
*         ^ End of loop over batches of L1
         END DO
*.       ^ End of loop over batches of K2
       END IF
*.     ^ End if T2 is unit operator (zero crea/anni ops)
       END DO
*      ^ End of loop over D2 batches
      END IF
*     ^ End if nonvanishing number of strings
      END DO
*     ^ End of loop over symmetries of D2 batches
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Updated HT1T2 block '
        CALL WRT_TCC_BLK(OT1T2,IT2SM,NOT1T2(1,1),NOT1T2(1,2),
     &       NOT1T2(1,3),NOT1T2(1,4),NSMST)
      END IF
*
C?    STOP ' Enforced stop at end of OPCT1T2'
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'OPCTT ')
      RETURN
      END
      SUBROUTINE REFRM_OP(IOPC,IOPE,NOP,NGAS,IWAY,LDIMC)
*
* Reform between contracted and expanded forms of 
* an operator string
*
* Contracted form : IOPC(NOP,3) for each op GAS, ICA, IAB given
* Expanded   form : IOPE(NGAS,4) each CA CB AA AB string given
*
* IWAY = 1 : C => E
*      = 2 : E => C
*
* Jeppe Olsen, May 2000
*
      INCLUDE 'implicit.inc'
*. Input/Output
      INTEGER IOPC(LDIMC,3), IOPE(NGAS,4)
*
      IF(IWAY.EQ.1) THEN
*Compressed to Expanded form
       IZERO = 0
       CALL ISETVC(IOPE,IZERO,4*NGAS)
       DO IOP = 1, NOP
         IGAS = IOPC(IOP,1)
         IAB  = IOPC(IOP,2)
         ICA  = IOPC(IOP,3)
         I_ABCA = (ICA-1)*2 + IAB
         IOPE(IGAS,I_ABCA) = IOPE(IGAS,I_ABCA) + 1
       END DO
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Expanded to Compressed form 
       IOP = 0
       I_ABCA = 0
       DO ICA = 1, 2
        DO IAB = 1, 2
          I_ABCA = I_ABCA + 1
          DO IGAS = 1, NGAS
           LOP = IOPE(IGAS,I_ABCA)
           DO JOP = 1, LOP
            IOP = IOP + 1
            IOPC(JOP,1) = IGAS
            IOPC(JOP,2) = IAB
            IOPC(JOP,3) = ICA
           END DO
          END DO
        END DO
       END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reforming between compact and expanded form of op'
        IF(IWAY.EQ.1) THEN
          WRITE(6,*)  ' Compressed => Expanded '
        ELSE IF (IWAY.EQ.2) THEN
          WRITE(6,*)  ' Expanded  => Compressed '
        END IF
*
        WRITE(6,*) ' Compressed operator'
C           WRT_CNTR(ICONT,NCONT,LDIM)
        CALL WRT_CNTR(IOPC,NOP,LDIMC)
        WRITE(6,*) ' Expanded operator '
        CALL WRT_SPOX_TP(IOPCE,1)
      END IF
*
      RETURN
      END
      SUBROUTINE CONTR_STR(ICONT,NCONT,JSTRING,ISTRING,IMZERO)
*
* ISTRING is result of contracting JSTRING with contraction 
* operator ICONT.
*
* Find ISTRING and determine if the contraction is vanishing
*
* Jeppe Olsen, May 2000
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER ICONT(4,3),JSTRING(NGAS,4)
*. Output
      INTEGER ISTRING(NGAS,4)
*
      CALL ICOPVE(JSTRING,ISTRING,4*NGAS)
      IMZERO = 0
      DO IOP = 1, NCONT
        IGAS = ICONT(IOP,1)
        IAB  = ICONT(IOP,2)
        ICA  = ICONT(IOP,3)
*. ICA was for Hamiltonian so the operator in ISTRING  should 
*  have the opposite C/A index
C?      WRITE(6,*) ' IGAS, ICA, IAB = ', IGAS,ICA,IAB
        ICA = 2/ICA
        ICA_AB = (ICA-1)*2 + IAB
C?      WRITE(6,*) ' ICA_AB = ', ICA_AB
        ISTRING(IGAS,ICA_AB) = ISTRING(IGAS,ICA_AB)-1
        IF(ISTRING(IGAS,ICA_AB).LT.0) IMZERO = 1
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Contraction operator, input and output strings'
C     WRT_CNTR(ICONT,NCONT,LDIM)
        CALL WRT_CNTR(ICONT,NCONT,4)
        CALL WRT_SPOX_TP(JSTRING,1)
        CALL WRT_SPOX_TP(ISTRING,1)
      END IF
*
      RETURN
      END
      SUBROUTINE CONTR_ORD(ICONT1,NCONT1,ICONT2,NCONT2,I12FIRST,
     &                     IT1OCC,IT2OCC)
*
* Hamiltonian will be contracted with the operators T1 and T2.
*
* Decide which operator should be contracted first
*
* 
*. The right order is of very significant 
*. importance of the present procedure. 
*
*     As first contraction choose the contraction with 
*     most elements as first contraction 
*
*     If the number of elements in the two contractions are 
*     identical, use as first contraction, the one with 
*     fewest elements ( reduces memory requirements)
*
*. Input
      INTEGER ICONT1(4,3),ICONT2(4,3)
      INTEGER IT1OCC(*), IT2OCC(*)
*. Local scratch 
      INTEGER LEN(8)
*
      CALL DIM_CNTR(ICONT1,NCONT1,4,LEN_CONT1)
      CALL DIM_CNTR(ICONT2,NCONT2,4,LEN_CONT2)
      IF(LEN_CONT1.GT.LEN_CONT2) THEN
        I12FIRST = 1
      ELSE IF (LEN_CONT1.LT.LEN_CONT2) THEN
        I12FIRST = 2
      ELSE IF (LEN_CONT1.EQ.LEN_CONT2) THEN
        CALL  DIM_TCC_OP(IT1OCC,LEN)
        LEN1 = LEN(1)
        CALL  DIM_TCC_OP(IT2OCC,LEN)
        LEN2 = LEN(1)
        IF(LEN1.GE.LEN2) THEN
          I12FIRST = 2
        ELSE 
          I12FIRST = 1
        END IF
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' First operator to be contracted = ', I12FIRST
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_CNTR2(ICONT_SPC,ICONT_CA,NCONT)
*
* Write contraction operator ICONT consisting just of 
* of orbital space and C/A
*
      INCLUDE 'implicit.inc'
      INTEGER ICONT_SPC(NCONT), ICONT_CA(NCONT)
*
      WRITE(6,*) ' Gasspace  Cr/An  '
      WRITE(6,*) ' ================='
      DO JCONT = 1, NCONT
        WRITE(6,'(I4,6X,I2,7X)') 
     &  ICONT_SPC(JCONT),ICONT_CA(JCONT)
      END DO
*
      RETURN
      END 
      SUBROUTINE WRT_CNTR(ICONT,NCONT,LDIM)
*
* Write contraction operator ICONT
*
      INCLUDE 'implicit.inc'
      INTEGER ICONT(LDIM,4)
*
      WRITE(6,*)
     &' Information about operators to be contracted'
      WRITE(6,*) ' Gasspace  Spin    Cr/An    Index'
      WRITE(6,*) ' ================================'
      DO JCONT = 1, NCONT
        WRITE(6,'(I4,6X,I2,7X,I2,7X,I2)') 
     &  ICONT(JCONT,1),ICONT(JCONT,2),ICONT(JCONT,3),ICONT(JCONT,4)
      END DO
*
      RETURN
      END 
      SUBROUTINE T1T2_TO_T12_OCC(I1,I2,I12,NGAS)
*
*  Occupation of two excitation strings T1, T2 are given 
*  Find occupation of T1T2 operator 
*
*. Input
      INTEGER I1(NGAS,4),I2(NGAS,4)
*. Output 
      INTEGER I12(NGAS,4)
*
* Assuming no overlap between creation and annihilation spaces 
*
* Jeppe Olsen, May 2000
*
*
      IONE = 1
      CALL IVCSUM(I12(1,1),I1(1,1),I2(1,1),IONE,IONE,NGAS)
      CALL IVCSUM(I12(1,2),I1(1,2),I2(1,2),IONE,IONE,NGAS)
      CALL IVCSUM(I12(1,3),I1(1,3),I2(1,3),IONE,IONE,NGAS)
      CALL IVCSUM(I12(1,4),I1(1,4),I2(1,4),IONE,IONE,NGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T1T2_TO_T12_OCC '
        WRITE(6,*) ' Input OCC1 OCC2 '
        WRITE(6,*)
        CALL WRT_SPOX_TP(I1,1)
        CALL WRT_SPOX_TP(I2,1)
        WRITE(6,*) ' Output '
        CALL WRT_SPOX_TP(I12,1)
      END IF
*
      RETURN
      END 
      SUBROUTINE H_T1T2(IH,IHINDX,HFAC,IHSM,
     &                  IT1OCC,IT1SM,T1,IT2OCC,IT2SM,T2,HT1T2,
     &                  LT1,LT2,LHT1T2)
*
* 
* Contract Hamiltonian operator defined by IH with T1 and T2
*
* Disconnected terms included
*
* Jeppe Olsen, May 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. Input  : Occupations in each GAS space of T1,T2 and H
      INTEGER IT1OCC(NGAS,4),IT2OCC(NGAS,4),IH(NGAS,4)
      INTEGER IHINDX(4)
*. Original indeces of each operator in H
      INTEGER IEXD1D2_INDX(4)
*. And the numerical values of T1, T2
      DIMENSION T1(*),T2(*)
*. Input common blocks
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Offsets for scratch used for HT1T2 evaluation
      COMMON/CC_SCR/ KCCF,KCCE,KCCVEC1,KCCVEC2,
     &              KISTST1,KXSTST1,KISTST2,KXSTST2,
     &              KISTST3,KXSTST3,KISTST4,KXSTST4,
     &              KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &              KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTREO,
     &              KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &              KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &              KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &              KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &              KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &              KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &              KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &              KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &              KLTSCR1,KLTSCR2,KLTSCR3,KLTSCR4,
     &              KLOPSCR,
     &              KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &              KLSMD1,KLSMD2,KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &              KLNMD1,KLNMD2,KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &              KLOCK1, KLOCK2, KLOCL1, KLOCOT1T2, KL_IBF, 
     &              KLEXEORD
*. Output
      DIMENSION HT1T2(*)
*. Local scratch, assuming atmost 4 operators in Hamiltonian 
      INTEGER IHDEEX_CR(4,3), IHEX_CR(4,3)
      INTEGER IHDEEX_AN(4,3), IHEX_AN(4,3)
C     INTEGER IKOCC(MXPNGAS*4) 
      INTEGER IK1OCC(MXPNGAS*4),IK2OCC(MXPNGAS*4)
      INTEGER IHD1OCC(MXPNGAS*4),IHD2OCC(MXPNGAS*4),IHEXOCC(MXPNGAS*4)
C??   INTEGER JEPTEST(16)
*
      INTEGER ICONT1(4,4), ICONT2(4,4)
C?    CALL LUCIAQENTER('HT1T2')
*
      NTEST = 00
*. Total number of operators in H
      NHOP = IELSUM(IH,4*NGAS)
C?    WRITE(6,*) ' NHOP = ', NHOP
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' H_T1T2 in action '
        WRITE(6,*) ' ================ '
        WRITE(6,*)
        WRITE(6,*) ' Occupation of H, T1 and T2 operator '
        WRITE(6,*)
        CALL WRT_SPOX_TP(IH,1)
        CALL WRT_SPOX_TP(IT1OCC,1)
        CALL WRT_SPOX_TP(IT2OCC,1)
        WRITE(6,*) ' Indeces of IH '
        CALL IWRTMA(IHINDX,1,NHOP,1,NHOP)
*
        WRITE(6,*) ' Block of T1 ' 
        CALL WRTMAT(T1,1,LT1,1,LT1)
        WRITE(6,*) ' Block of T2 ' 
        CALL WRTMAT(T2,1,LT2,1,LT2)
        WRITE(6,*) ' Input block of HT1T2 '
        CALL WRTMAT(HT1T2,1,LHT1T2,1,LHT1T2)
      END IF
*. Check to see if T1 and T2 are identical , or more exact if T2 = 0.5*T1
C     CMP_CAAB(ICAAB1,ICAAB2,IDENTICAL)
      CALL CMP_CAAB(IT1OCC,IT2OCC,IT1EQT2)
      IF(IT1EQT2.EQ.1) THEN
C?    WRITE(6,*) ' Identical occupations : '
C?    CALL WRT_SPOX_TP(IT1OCC,1)
C?    CALL WRT_SPOX_TP(IT2OCC,1)
C?    WRITE(6,*) ' Indeces of IH '
C?    CALL IWRTMA(IHINDX,1,NHOP,1,NHOP)
*. Check to see if T1 and T2 have identical elements 
        XDIFF = 0.0D0
        XT1TOT = 0.0D0 
        HALF = 0.5D0
        DO I = 1, LT1
          XDIFF = XDIFF + ABS(HALF*T1(I)-T2(I))
          XT1TOT = XT1TOT + ABS(T1(I))
        END DO
        IF(XDIFF/XT1TOT.GT.1.0D-13) IT1EQT2 = 0
C?      WRITE(6,*) ' LT1, LT2, XDIFF, XT1TOT, DIFF/TOT = ',
C?   &               LT1, LT2, XDIFF, XT1TOT, XDIFF/XT1TOT
C??     CALL ISETVC(JEPTEST,0,16)
      END IF
*. TEST TEST For checking of problem : enforce difference 
      IT1EQT2 = 0
*. Check to see if operator has symmetry between particle one and two
      I_USE_SYM12 = 1
      IF(I_USE_SYM12.EQ.1) THEN
         CALL H_TYPE_SYM12(IH,NGAS,ISYM12)
      ELSE 
         ISYM12 = 0
      END IF
*
* Divide operators in H into excitation and deexcitation operators 
*
* Excitation : Annihilation of hole + creation of particle
* deexcitation : Annihilation of particle + creation of hole
      NHDEEX_CR = 0
      NHDEEX_AN = 0
      NHEX_CR = 0
      NHEX_AN = 0
*. Loop over operators in HCA HCB HAA HAB
      IOPT = 0
      ICAAB = 0
      DO ICA = 1, 2
       DO IAB = 1, 2
        ICAAB = ICAAB + 1
        DO JOBTP = 1, NGAS
          JDEEX = 0
*. creation of hole is deexcitation
          IF(IHPVGAS_AB(JOBTP,IAB).EQ.1.AND.ICA.EQ.1) JDEEX = 1
*. Annihilation of particle is deexcitation 
          IF(IHPVGAS_AB(JOBTP,IAB).EQ.2.AND.ICA.EQ.2) JDEEX = 1
*
          NOP = IH(JOBTP,ICAAB)
          DO JOB = 1, NOP
            IOPT = IOPT + 1
            IF(ICA.EQ.1.AND.JDEEX.EQ.1) THEN
             NHDEEX_CR = NHDEEX_CR + 1
             IHDEEX_CR(NHDEEX_CR,1) = JOBTP
             IHDEEX_CR(NHDEEX_CR,2) = IAB    
             IHDEEX_CR(NHDEEX_CR,3) = IHINDX(IOPT)
            ELSE IF (ICA.EQ.2.AND.JDEEX.EQ.1) THEN
             NHDEEX_AN = NHDEEX_AN + 1
             IHDEEX_AN(NHDEEX_AN,1) = JOBTP
             IHDEEX_AN(NHDEEX_AN,2) = IAB    
             IHDEEX_AN(NHDEEX_AN,3) = IHINDX(IOPT)
            ELSE IF (ICA.EQ.1.AND.JDEEX.EQ.0) THEN
             NHEX_CR = NHEX_CR + 1
             IHEX_CR(NHEX_CR,1) = JOBTP
             IHEX_CR(NHEX_CR,2) = IAB   
             IHEX_CR(NHEX_CR,3) = IHINDX(IOPT)
            ELSE IF (ICA.EQ.2.AND.JDEEX.EQ.0) THEN
             NHEX_AN = NHEX_AN + 1
             IHEX_AN(NHEX_AN,1) = JOBTP
             IHEX_AN(NHEX_AN,2) = IAB    
             IHEX_AN(NHEX_AN,3) = IHINDX(IOPT)
            END IF
          END DO
        END DO
*       ^ End of loop over JOBTP
      END DO
*     ^ End of loop over IAB
      END DO
*     ^ End of loop over ICA
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Number of deexcitation creation operators ', NHDEEX_CR
        WRITE(6,*) 
     &  ' Deexcitation creation operators, GAS, spin and index :'
        CALL IWRTMA(IHDEEX_CR(1,1),1,NHDEEX_CR,1,NHDEEX_CR)
        CALL IWRTMA(IHDEEX_CR(1,2),1,NHDEEX_CR,1,NHDEEX_CR)
        CALL IWRTMA(IHDEEX_CR(1,3),1,NHDEEX_CR,1,NHDEEX_CR)
        WRITE(6,*) 
     &  ' Number of deexcitation annihilation operators ', NHDEEX_AN
        WRITE(6,*) 
     &  ' Deexcitation annihilation operators, GAS, spin and index  :'
        CALL IWRTMA(IHDEEX_AN(1,1),1,NHDEEX_AN,1,NHDEEX_AN)
        CALL IWRTMA(IHDEEX_AN(1,2),1,NHDEEX_AN,1,NHDEEX_AN)
        CALL IWRTMA(IHDEEX_AN(1,3),1,NHDEEX_AN,1,NHDEEX_AN)
        WRITE(6,*) 
     &  ' Number of excitation creation operators ', NHEX_CR
        WRITE(6,*) 
     &  ' Excitation creation operators, GAS and spin and index :'
        CALL IWRTMA(IHEX_CR(1,1),1,NHEX_CR,1,NHEX_CR)
        CALL IWRTMA(IHEX_CR(1,2),1,NHEX_CR,1,NHEX_CR)
        CALL IWRTMA(IHEX_CR(1,3),1,NHEX_CR,1,NHEX_CR)
        WRITE(6,*) 
     &  ' Number of excitation annihilation operators ', NHEX_AN
        WRITE(6,*) 
     &  ' Excitation annihilation operators, GAS, spin and index '
        CALL IWRTMA(IHEX_AN(1,1),1,NHEX_AN,1,NHEX_AN)
        CALL IWRTMA(IHEX_AN(1,2),1,NHEX_AN,1,NHEX_AN)
        CALL IWRTMA(IHEX_AN(1,3),1,NHEX_AN,1,NHEX_AN)
*.      
      END IF
*. The deexcitation must all be contracted with excitations operators
*. Consider f.ex the possible nontrivial deexcitation creation operators :
* 1 creation operator  => 2 possibilities 
*                      1 : contract with T1
*                      2 : contract with T2
* 2 creation operators => 4 possibilities 
*                 Contract with T1     Contract with T2
*             ===========================================
*              1 : 1,2              
*              2 :                        1,2 
*              3 : 1                      2
*              4 : 2                      1
* 
      NDEEX_CR_TP = 0
      NDEEX_AN_TP = 0
      IF(NHDEEX_CR.EQ.0) THEN
        NDEEX_CR_TP = 1
      ELSE IF (NHDEEX_CR.EQ.1) THEN
        NDEEX_CR_TP = 2
      ELSE IF (NHDEEX_CR.EQ.2) THEN
        NDEEX_CR_TP = 4
      END IF
      IF(NHDEEX_AN.EQ.0) THEN
        NDEEX_AN_TP = 1
      ELSE IF (NHDEEX_AN.EQ.1) THEN
        NDEEX_AN_TP = 2
      ELSE IF (NHDEEX_AN.EQ.2) THEN
        NDEEX_AN_TP = 4
      END IF
*. And loop over the different contraction possibilities
      DO IDEEX_CR_TP = 1, NDEEX_CR_TP
      DO IDEEX_AN_TP = 1, NDEEX_AN_TP
        IF(NTEST.GE.100) WRITE(6,*) ' IDEEX_CR_TP, IDEEX_AN_TP = ',
     &               IDEEX_CR_TP, IDEEX_AN_TP
         FACX = HFAC
*
*  Define contraction in terms of ICONT1, ICONT2
*
* ICONTX(I,1) : Type of GAS space to be contracted
* ICONTX(I,2) : Spin of index     to be contracted
* ICONTX(I,3) : Does this index correspond to creation or annihilation in H
* ICONTX(I,4) : Index in original Hamiltonian 
*
* FACX :
* The terms with 4 alpha or 4 beta reads
* sum(i.gt.k,j.gt.l) (ij ! kl) a+i a+k al aj
* when indeces ik and jl are split in different operators, 
* the restrictions are not enforced when the orbitals ik (jl)
* belong to the same type. Therefore, there is included a 
* factor for these terms
*
*. Indeces to be contracted with T1           
       NCONT1 = 0
       NCONT2 = 0
*. Creation deexcitations
       ICA = 1
       IF(IDEEX_CR_TP.EQ.1) THEN
         ICA = 1
         IF(NHDEEX_CR.EQ.1) THEN
           NCONT1=NCONT1 + 1
           ICONT1(NCONT1,1) = IHDEEX_CR(1,1)
           ICONT1(NCONT1,2) = IHDEEX_CR(1,2)
           ICONT1(NCONT1,3) = ICA
           ICONT1(NCONT1,4) = IHDEEX_CR(1,3)
         ELSE IF(NHDEEX_CR.EQ.2) THEN
           NCONT1=NCONT1 + 1
           ICONT1(NCONT1,1) = IHDEEX_CR(1,1)
           ICONT1(NCONT1,2) = IHDEEX_CR(1,2)
           ICONT1(NCONT1,3) = ICA
           ICONT1(NCONT1,4) = IHDEEX_CR(1,3)
           NCONT1 = NCONT1 + 1
           ICONT1(NCONT1,1) = IHDEEX_CR(2,1)
           ICONT1(NCONT1,2) = IHDEEX_CR(2,2)
           ICONT1(NCONT1,3) = ICA
           ICONT1(NCONT1,4) = IHDEEX_CR(2,3)
         END IF
       ELSE IF (IDEEX_CR_TP.EQ.2) THEN 
         ICA = 1
         IF(NHDEEX_CR.EQ.1) THEN
           NCONT2=NCONT2 + 1
           ICONT2(NCONT2,1) = IHDEEX_CR(1,1)
           ICONT2(NCONT2,2) = IHDEEX_CR(1,2)
           ICONT2(NCONT2,3) = ICA
           ICONT2(NCONT2,4) = IHDEEX_CR(1,3)
         ELSE IF(NHDEEX_CR.EQ.2) THEN
           NCONT2=NCONT2 + 1
           ICONT2(NCONT2,1) = IHDEEX_CR(1,1)
           ICONT2(NCONT2,2) = IHDEEX_CR(1,2)
           ICONT2(NCONT2,3) = ICA
           ICONT2(NCONT2,4) = IHDEEX_CR(1,3)
           NCONT2 = NCONT2 + 1
           ICONT2(NCONT2,1) = IHDEEX_CR(2,1)
           ICONT2(NCONT2,2) = IHDEEX_CR(2,2)
           ICONT2(NCONT2,3) = ICA
           ICONT2(NCONT2,4) = IHDEEX_CR(2,3)
         END IF
       ELSE IF (IDEEX_CR_TP.EQ.3) THEN
         NCONT1 = NCONT1 + 1
         NCONT2 = NCONT2 + 1
         IF( IHDEEX_CR(1,1).EQ.IHDEEX_CR(2,1).AND.
     &       IHDEEX_CR(1,2).EQ.IHDEEX_CR(2,2)) FACX = 0.5D0*FACX
         ICONT1(NCONT1,1) = IHDEEX_CR(1,1)
         ICONT1(NCONT1,2) = IHDEEX_CR(1,2)
         ICONT1(NCONT1,3) = ICA
         ICONT1(NCONT1,4) = IHDEEX_CR(1,3)
*
         ICONT2(NCONT2,1) = IHDEEX_CR(2,1)
         ICONT2(NCONT2,2) = IHDEEX_CR(2,2)
         ICONT2(NCONT2,3) = ICA
         ICONT2(NCONT1,4) = IHDEEX_CR(2,3)
       ELSE IF (IDEEX_CR_TP.EQ.4) THEN
         NCONT1 = NCONT1 + 1
         NCONT2 = NCONT2 + 1
C        IF( IHDEEX_CR(1,1).EQ.IHDEEX_CR(2,1)) FACX = 0.5D0*FACX
         IF( IHDEEX_CR(1,1).EQ.IHDEEX_CR(2,1).AND.
     &       IHDEEX_CR(1,2).EQ.IHDEEX_CR(2,2)) FACX = 0.5D0*FACX
         ICONT1(NCONT1,1) = IHDEEX_CR(2,1)
         ICONT1(NCONT1,2) = IHDEEX_CR(2,2)
         ICONT1(NCONT1,3) = ICA
         ICONT1(NCONT1,4) = IHDEEX_CR(2,3)
         ICONT2(NCONT2,1) = IHDEEX_CR(1,1)
         ICONT2(NCONT2,2) = IHDEEX_CR(1,2)
         ICONT2(NCONT2,3) = ICA
         ICONT2(NCONT2,4) = IHDEEX_CR(1,3)
       END IF
*. Annihilation deexcitations
       ICA = 2
       IF(IDEEX_AN_TP.EQ.1) THEN
         ICA = 2
         IF(NHDEEX_AN.EQ.1) THEN
           NCONT1=NCONT1 + 1
           ICONT1(NCONT1,1) = IHDEEX_AN(1,1)
           ICONT1(NCONT1,2) = IHDEEX_AN(1,2)
           ICONT1(NCONT1,3) = ICA
           ICONT1(NCONT1,4) = IHDEEX_AN(1,3)
         ELSE IF(NHDEEX_AN.EQ.2) THEN
           NCONT1=NCONT1 + 1
           ICONT1(NCONT1,1) = IHDEEX_AN(1,1)
           ICONT1(NCONT1,2) = IHDEEX_AN(1,2)
           ICONT1(NCONT1,3) = ICA
           ICONT1(NCONT1,4) = IHDEEX_AN(1,3)
           NCONT1 = NCONT1 + 1
           ICONT1(NCONT1,1) = IHDEEX_AN(2,1)
           ICONT1(NCONT1,2) = IHDEEX_AN(2,2)
           ICONT1(NCONT1,3) = ICA
           ICONT1(NCONT1,4) = IHDEEX_AN(2,3)
         END IF
       ELSE IF (IDEEX_AN_TP.EQ.2) THEN 
         ICA = 2
         IF(NHDEEX_AN.EQ.1) THEN
           NCONT2=NCONT2 + 1
           ICONT2(NCONT2,1) = IHDEEX_AN(1,1)
           ICONT2(NCONT2,2) = IHDEEX_AN(1,2)
           ICONT2(NCONT2,3) = ICA
           ICONT2(NCONT2,4) = IHDEEX_AN(1,3)
         ELSE IF(NHDEEX_AN.EQ.2) THEN
           NCONT2=NCONT2 + 1
           ICONT2(NCONT2,1) = IHDEEX_AN(1,1)
           ICONT2(NCONT2,2) = IHDEEX_AN(1,2)
           ICONT2(NCONT2,3) = ICA
           ICONT2(NCONT2,4) = IHDEEX_AN(1,3)
           NCONT2 = NCONT2 + 1
           ICONT2(NCONT2,1) = IHDEEX_AN(2,1)
           ICONT2(NCONT2,2) = IHDEEX_AN(2,2)
           ICONT2(NCONT2,3) = ICA
           ICONT2(NCONT2,4) = IHDEEX_AN(2,3)
         END IF
       ELSE IF (IDEEX_AN_TP.EQ.3) THEN
         ICA = 2
         NCONT1 = NCONT1 + 1
         NCONT2 = NCONT2 + 1
         IF( IHDEEX_AN(1,1).EQ.IHDEEX_AN(2,1).AND.
     &       IHDEEX_AN(1,2).EQ.IHDEEX_AN(2,2)     ) FACX = 0.5D0*FACX
         ICONT1(NCONT1,1) = IHDEEX_AN(1,1)
         ICONT1(NCONT1,2) = IHDEEX_AN(1,2)
         ICONT1(NCONT1,3) = ICA
         ICONT1(NCONT1,4) = IHDEEX_AN(1,3)
         ICONT2(NCONT2,1) = IHDEEX_AN(2,1)
         ICONT2(NCONT2,2) = IHDEEX_AN(2,2)
         ICONT2(NCONT2,3) = ICA
         ICONT2(NCONT2,4) = IHDEEX_AN(2,3)
       ELSE IF (IDEEX_AN_TP.EQ.4) THEN
         ICA = 2
         NCONT1 = NCONT1 + 1
         NCONT2 = NCONT2 + 1
C        IF( IHDEEX_AN(1,1).EQ.IHDEEX_AN(2,1)) FACX = 0.5D0*FACX
         IF( IHDEEX_AN(1,1).EQ.IHDEEX_AN(2,1).AND.
     &       IHDEEX_AN(1,2).EQ.IHDEEX_AN(2,2)     ) FACX = 0.5D0*FACX
         ICONT1(NCONT1,1) = IHDEEX_AN(2,1)
         ICONT1(NCONT1,2) = IHDEEX_AN(2,2)
         ICONT1(NCONT1,3) = ICA
         ICONT1(NCONT1,4) = IHDEEX_AN(2,3)
         ICONT2(NCONT2,1) = IHDEEX_AN(1,1)
         ICONT2(NCONT2,2) = IHDEEX_AN(1,2)
         ICONT2(NCONT2,3) = ICA
         ICONT2(NCONT2,4) = IHDEEX_AN(1,3)
       END IF
*
       IF(NTEST.GE.100) THEN
         WRITE(6,*)
     &   ' Information about operator to be contracted with T1'
         CALL WRT_CNTR(ICONT1,NCONT1,4)
         WRITE(6,*)
     &   ' Information about operator to be contracted with T2'
         CALL WRT_CNTR(ICONT2,NCONT2,4)
       END IF
*. Order of contraction
       CALL CONTR_ORD(ICONT1,NCONT1,ICONT2,NCONT2,I12FIRST,
     &                IT1OCC,IT2OCC)
*. Find strings resulting from contraction 
       CALL CONTR_STR(ICONT1,NCONT1,IT1OCC,IK1OCC,IMZERO1)
       CALL CONTR_STR(ICONT2,NCONT2,IT2OCC,IK2OCC,IMZERO2)
       IF(IMZERO1.EQ.0.AND.IMZERO2.EQ.0) THEN
C?     WRITE(6,*) ' Contraction is active '
*. Number of operators in excitation part
        NHEX = NHOP - NCONT1 - NCONT2
*. Contraction operators as CAAB arrays
        IZERO = 0
        CALL ISETVC(IHD1OCC,IZERO,4*NGAS)
        DO JCONT = 1, NCONT1
          IGAS = ICONT1(JCONT,1)
          IAB  = ICONT1(JCONT,2)
          ICA  = ICONT1(JCONT,3)
          ICAAB = (ICA-1)*2+IAB 
          IHD1OCC((ICAAB-1)*NGAS+IGAS) = 
     &    IHD1OCC((ICAAB-1)*NGAS+IGAS) + 1
          INDEX = ICONT1(JCONT,4)
          IEXD1D2_INDX(NHEX+JCONT) = INDEX
         END DO
*
        CALL ISETVC(IHD2OCC,IZERO,4*NGAS)
        DO JCONT = 1, NCONT2
          IGAS = ICONT2(JCONT,1)
          IAB  = ICONT2(JCONT,2)
          ICA  = ICONT2(JCONT,3)
          ICAAB = (ICA-1)*2+IAB 
          IHD2OCC((ICAAB-1)*NGAS+IGAS) = 
     &    IHD2OCC((ICAAB-1)*NGAS+IGAS) + 1
          INDEX = ICONT2(JCONT,4)
          IEXD1D2_INDX(NHEX+NCONT1+JCONT) = INDEX
        END DO
*. Excitation part of Hamilton operator 
        CALL ICOPVE(IH,IHEXOCC,4*NGAS)
        MONE = -1
        IONE =  1
        CALL IVCSUM(IHEXOCC,IHEXOCC,IHD1OCC,IONE,MONE,4*NGAS)
        CALL IVCSUM(IHEXOCC,IHEXOCC,IHD2OCC,IONE,MONE,4*NGAS)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' D1 contraction part of H in CAAB form '
          CALL WRT_SPOX_TP(IHD1OCC,1)
          WRITE(6,*) ' D2 contraction part of H in CAAB form '
          CALL WRT_SPOX_TP(IHD2OCC,1)
          WRITE(6,*) ' Excitation part of H in CAAB form '
          CALL WRT_SPOX_TP(IHEXOCC,1)
        END IF
*
*. Operators for the excitation part : Those left after decontracting
        IIEX = 0
        DO IORIG = 1, NHOP
*. Was original index IORIG used in decontracting
          IUSED = 0
          DO ICINDEX = 1, NCONT1+NCONT2
C           IF(IEXD1D2_INDX(NHEX+ICINDEX).EQ.IORIG)IUSED = 1
            IF(IEXD1D2_INDX(NHEX+ICINDEX).EQ.IHINDX(IORIG))IUSED = 1
          END DO
          IF(IUSED.EQ.0) THEN 
*. Operator in H that corresponds to original index IORIG
            IIEX = IIEX + 1
            IEXD1D2_INDX(IIEX) = IHINDX(IORIG)
          END IF
        END DO
*. Sign for bringing H into (ex d2 d1)
        ISIGN1 = IPERM_PARITY(IEXD1D2_INDX,NHOP)
*. Move indeces of T1 and T2 around if planned order is T2 T1
        IF(I12FIRST.EQ.2) THEN 
*
          DO JCONT = 1, NCONT2
            INDEX = ICONT2(JCONT,4)
            IEXD1D2_INDX(NHEX+JCONT) = INDEX
          END DO
*
          DO JCONT = 1, NCONT1
            INDEX = ICONT1(JCONT,4)
            IEXD1D2_INDX(NHEX+NCONT2+JCONT) = INDEX
          END DO
        END IF
*. if particlesym is present, combinations 4 3 and 44 are skipped 
*. we are dealing with terms where all operators have the same 
*. spin symmetry so 34 can also be skipped 
        ISKIP = 0
        IREDUCED = 0
        IF(ISYM12.EQ.1) THEN
          IF((IDEEX_CR_TP.EQ.4.AND.IDEEX_AN_TP.EQ.4).OR.
     &       (IDEEX_CR_TP.EQ.4.AND.IDEEX_AN_TP.EQ.3).OR. 
     &       (IDEEX_CR_TP.EQ.3.AND.IDEEX_AN_TP.EQ.4)   ) THEN 
             ISKIP = 1
             IREDUCED = 1
          ELSE IF (IDEEX_CR_TP.EQ.3.AND.IDEEX_AN_TP.EQ.3) THEN 
             ISKIP = 0
             IREDUCED = 1
             FACX = 4.0D0*FACX
          END IF
        END IF
*. Number of terms may also be reduced if T1 = T2  
 
*. Test
C       IT1EQT2 = 0
        IF(IREDUCED.EQ.0.AND.IT1EQT2.EQ.1) THEN 
*. Assign a unique number to the current contraction 
          IINUM = ( IDEEX_CR_TP - 1 ) * 4 +  IDEEX_AN_TP 
          JDEEX_CR_TP = 0
          JDEEX_AN_TP = 0
*. Interchange T1 and T2
          IF(IDEEX_CR_TP.EQ.1) THEN
             JDEEX_CR_TP = 2
          ELSE IF (IDEEX_CR_TP.EQ.2) THEN
                   JDEEX_CR_TP = 1
          ELSE IF (IDEEX_CR_TP.EQ.3) THEN
                   JDEEX_CR_TP = 4
          ELSE IF (IDEEX_CR_TP.EQ.4) THEN
                   JDEEX_CR_TP = 3
          END IF
          IF(IDEEX_AN_TP.EQ.1) THEN
             JDEEX_AN_TP = 2
          ELSE IF (IDEEX_AN_TP.EQ.2) THEN
                   JDEEX_AN_TP = 1
          ELSE IF (IDEEX_AN_TP.EQ.3) THEN
                   JDEEX_AN_TP = 4
          ELSE IF (IDEEX_AN_TP.EQ.4) THEN
                   JDEEX_AN_TP = 3
          END IF
*. If there are no annihilations or creations keep original numbers 
          IF(NHDEEX_CR .EQ.  0)  JDEEX_CR_TP = 1
          IF(NHDEEX_AN .EQ.  0)  JDEEX_AN_TP = 1
*
          JJNUM = ( JDEEX_CR_TP - 1 ) * 4 +  JDEEX_AN_TP 
          IF(JJNUM.GT.IINUM) THEN
            ISKIP = 1
C??         JEPTEST(IINUM) = JEPTEST(IINUM)-1
          ELSE IF (JJNUM.LT.IINUM) THEN
            ISKIP = 0
            FACX = 2.0D0*FACX
C??         JEPTEST(IINUM) = JEPTEST(IINUM)+1
          END IF
        END IF
*       ^ End if T1 = T2, so futher simplifications where possible
   

        IF(ISKIP.EQ.0) THEN
*. Calculate ((H contracted T1) contracted T2)
        IF(I12FIRST.EQ.1) THEN 
          CALL OPCT1T2(IHEXOCC,IHD1OCC,IHD2OCC,IT1OCC,IT2OCC,
     &    T1,T2,HT1T2,IT1SM,IT2SM,IHSM,LCCBD12,LCCB,
     &    WORK(KLIOD2_ST),WORK(KLIOD1_ST),WORK(KLIOEX_ST),
     &    WORK(KIX1_CA),WORK(KSX1_CA),WORK(KIX1_CB),WORK(KSX1_CB),
     &    WORK(KIX1_AA),WORK(KSX1_AA),WORK(KIX1_AB),WORK(KSX1_AB),
     &    WORK(KIX2_CA),WORK(KSX2_CA),WORK(KIX2_CB),WORK(KSX2_CB),
     &    WORK(KIX2_AA),WORK(KSX2_AA),WORK(KIX2_AB),WORK(KSX2_AB),
     &    WORK(KIX3_CA),WORK(KSX3_CA),WORK(KIX3_CB),WORK(KSX3_CB),
     &    WORK(KIX3_AA),WORK(KSX3_AA),WORK(KIX3_AB),WORK(KSX3_AB),
     &    WORK(KIX4_CA),WORK(KSX4_CA),WORK(KIX4_CB),WORK(KSX4_CB),
     &    WORK(KIX4_AA),WORK(KSX4_AA),WORK(KIX4_AB),WORK(KSX4_AB),
     &    WORK(KLSTOCC1),WORK(KLSTOCC2),WORK(KLSTOCC3),
     &    WORK(KLSTOCC4),WORK(KLSTOCC5),WORK(KLSTOCC6),
     &    WORK(KLTSCR1),WORK(KLTSCR2),WORK(KLTSCR3),
     &    WORK(KLTSCR4),WORK(KLOPSCR),
     &    WORK(KLSMD1),WORK(KLSMD2),WORK(KLSMK1),WORK(KLSMK2),
     &    WORK(KLSMEX),WORK(KLSML1),
     &    WORK(KLNMD1),WORK(KLNMD2),WORK(KLNMK1),WORK(KLNMK2),
     &    WORK(KLNMEX),WORK(KLNML1),IEXD1D2_INDX,
     &    WORK(KLOCK1),WORK(KLOCK2),WORK(KLOCL1),WORK(KLOCOT1T2),
     &    WORK(KLZ),WORK(KLZSCR),WORK(KLSTREO),ISIGN1,FACX,
     &    N_TDL_MAX)
        ELSE
          ISIGN1 = ISIGN1*(-1)**(NCONT2*NCONT1)
*. ( it is assumed that T1 and T2 contains an even number of ops)
*. Calculate ((H contracted T2) contracted T1)
          CALL OPCT1T2(IHEXOCC,IHD2OCC,IHD1OCC,IT2OCC,IT1OCC,
     &    T2,T1,HT1T2,IT2SM,IT1SM,IHSM,LCCBD12,LCCB,
     &    WORK(KLIOD2_ST),WORK(KLIOD1_ST),WORK(KLIOEX_ST),
     &    WORK(KIX1_CA),WORK(KSX1_CA),WORK(KIX1_CB),WORK(KSX1_CB),
     &    WORK(KIX1_AA),WORK(KSX1_AA),WORK(KIX1_AB),WORK(KSX1_AB),
     &    WORK(KIX2_CA),WORK(KSX2_CA),WORK(KIX2_CB),WORK(KSX2_CB),
     &    WORK(KIX2_AA),WORK(KSX2_AA),WORK(KIX2_AB),WORK(KSX2_AB),
     &    WORK(KIX3_CA),WORK(KSX3_CA),WORK(KIX3_CB),WORK(KSX3_CB),
     &    WORK(KIX3_AA),WORK(KSX3_AA),WORK(KIX3_AB),WORK(KSX3_AB),
     &    WORK(KIX4_CA),WORK(KSX4_CA),WORK(KIX4_CB),WORK(KSX4_CB),
     &    WORK(KIX4_AA),WORK(KSX4_AA),WORK(KIX4_AB),WORK(KSX4_AB),
     &    WORK(KLSTOCC1),WORK(KLSTOCC2),WORK(KLSTOCC3),
     &    WORK(KLSTOCC4),WORK(KLSTOCC5),WORK(KLSTOCC6),
     &    WORK(KLTSCR1),WORK(KLTSCR2),WORK(KLTSCR3),
     &    WORK(KLTSCR4),WORK(KLOPSCR),
     &    WORK(KLSMD1),WORK(KLSMD2),WORK(KLSMK1),WORK(KLSMK2),
     &    WORK(KLSMEX),WORK(KLSML1),
     &    WORK(KLNMD1),WORK(KLNMD2),WORK(KLNMK1),WORK(KLNMK2),
     &    WORK(KLNMEX),WORK(KLNML1),IEXD1D2_INDX,
     &    WORK(KLOCK1),WORK(KLOCK2),WORK(KLOCL1),WORK(KLOCOT1T2),
     &    WORK(KLZ),WORK(KLZSCR),WORK(KLSTREO),ISIGN1,FACX,
     &    N_TDL_MAX)
        END IF
*       ^ End of I12FIRST switch
       END IF 
*      ^ End if ISKIP = 0
       END IF
*      ^ End if contraction is nonvanishing
      END DO
      END DO
*.    ^ End of loop over possible divisions of crea and anni contractions
*
C??   IF(IT1EQT2.EQ.1) THEN
C??     WRITE(6,*) ' JEPTEST array '
C??     CALL IWRTMA(JEPTEST,4,4,4,4)
C??   END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated block of HT1T2 '
        CALL WRTMAT(HT1T2,1,LHT1T2,1,LHT1T2)
      END IF
*
C?    CALL LUCIAQEXIT('HT1T2')
      RETURN
      END
      SUBROUTINE T1T2_TO_T12_MAP(I1SPOBEX,I2SPOBEX,I12SPOBEX,
     &                       ICA_MAP,XCA_MAP,ICB_MAP,XCB_MAP,
     &                       IAA_MAP,XAA_MAP,IAB_MAP,XAB_MAP,
     &                       IB,I1OCC,I2OCC,I1REO,IZ,IZSCR,SIGNP) 
*
* Obtain mappings for T1T2 => T12 as 4 mappings for each CAAB component
*
* All creation and annihilation operators are assumed to commute.
*
* Jeppe Olsen, Oct 2000     
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
*
*. Input
*
      INTEGER I1SPOBEX(NGAS,4)
      INTEGER I2SPOBEX(NGAS,4)
      INTEGER I12SPOBEX(NGAS,4)
*
*. Output
*
*. T1*T2 => T12 mappings
      INTEGER ICA_MAP(*),ICB_MAP(*),IAA_MAP(*),IAB_MAP(*)               
      DIMENSION XCA_MAP(*),XCB_MAP(*),XAA_MAP(*),XAB_MAP(*)               
      INTEGER IB(8,8,4)
*
*. Local scratch 
*
      INTEGER NI1CA(MXPNGAS),NI1CB(MXPNGAS),
     &        NI1AA(MXPNGAS),NI1AB(MXPNGAS)
      INTEGER NI2CA(MXPNGAS),NI2CB(MXPNGAS),
     &        NI2AA(MXPNGAS),NI2AB(MXPNGAS)
      INTEGER NI12CA(MXPNGAS),NI12CB(MXPNGAS),
     &        NI12AA(MXPNGAS),NI12AB(MXPNGAS)
      INTEGER I1CA_EXP(100),I1CB_EXP(100),I1AA_EXP(100),I1AB_EXP(100)
      INTEGER I2CA_EXP(100),I2CB_EXP(100),I2AA_EXP(100),I2AB_EXP(100)
      INTEGER I12CA_EXP(100),I12CB_EXP(100),
     &        I12AA_EXP(100),I12AB_EXP(100)
C     INTEGER ICA(100)
      INTEGER I1OFF(8,8,8), I2OFF(8,8,8), I12OFF(8,8,8)
C     INTEGER IB_CA_SS(8,8),IB_CB_SS(8,8),IB_AA_SS(8,8),IB_AB_SS(8,8)
*
*. Scratch through parameter list
*
*. For occ of two sets of excitation strings and an reorder array 
      INTEGER I1OCC(*),I2OCC(*),I1REO(*)
*. For Z-matrix and its construction
      INTEGER IZ(*),IZSCR(*)
*
      IDUM = 0
C?    CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TT_TM ')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' Welcome to T1T2_TO_T12_MAP '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Excitation operator 1 :'
        WRITE(6,*)
C   WRT_SPOX_TP(IEX_TP,NEX_TP)
        CALL WRT_SPOX_TP(I1SPOBEX,1)
        WRITE(6,*) ' Excitation operator 2 :'
        WRITE(6,*)
        CALL WRT_SPOX_TP(I2SPOBEX,1)
        WRITE(6,*) ' Compound operator '
        CALL WRT_SPOX_TP(I12SPOBEX,1)
      END IF
*. Number of operators in each string 
      N1CA_OP = IELSUM(I1SPOBEX(1,1),NGAS)
      N1CB_OP = IELSUM(I1SPOBEX(1,2),NGAS)
      N1AA_OP = IELSUM(I1SPOBEX(1,3),NGAS)
      N1AB_OP = IELSUM(I1SPOBEX(1,4),NGAS)
*
      N2CA_OP = IELSUM(I2SPOBEX(1,1),NGAS)
      N2CB_OP = IELSUM(I2SPOBEX(1,2),NGAS)
      N2AA_OP = IELSUM(I2SPOBEX(1,3),NGAS)
      N2AB_OP = IELSUM(I2SPOBEX(1,4),NGAS)
*
      N12CA_OP = IELSUM(I12SPOBEX(1,1),NGAS)
      N12CB_OP = IELSUM(I12SPOBEX(1,2),NGAS)
      N12AA_OP = IELSUM(I12SPOBEX(1,3),NGAS)
      N12AB_OP = IELSUM(I12SPOBEX(1,4),NGAS)
*. Sign to bring O(ca1)O(cb1)O(aa1)O(ab1)
*.               O(ca2)O(cb2)O(aa2)O(ab2)
*. into          O(ca1)O(ca2)O(cb1)O(cb2)O(aa1)O(aa2)O(ab1)O(ab2)
*
      NPERM =  N1AB_OP*( N2CA_OP+N2CB_OP+N2AA_OP)
     &      +  N1AA_OP*( N2CA_OP+N2CB_OP )
     &      +  N1CB_OP*  N2CA_OP
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNP = -1.0D0
      ELSE
        SIGNP = 1.0D0
      END IF

*
      MX_NOP = MAX(N12CA_OP,N12CB_OP,N12AA_OP,N12AB_OP)
*
C     CALL REF_OP(ICA,ICA_EXP,NCA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,1),I1CA_EXP,N1CA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,2),I1CB_EXP,N1CB_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,3),I1AA_EXP,N1AA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,4),I1AB_EXP,N1AB_OP,NGAS,1)
*
      CALL REF_OP(I2SPOBEX(1,1),I2CA_EXP,N2CA_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,2),I2CB_EXP,N2CB_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,3),I2AA_EXP,N2AA_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,4),I2AB_EXP,N2AB_OP,NGAS,1)
*
      CALL REF_OP(I12SPOBEX(1,1),I12CA_EXP,N12CA_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,2),I12CB_EXP,N12CB_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,3),I12AA_EXP,N12AA_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,4),I12AB_EXP,N12AB_OP,NGAS,1)
*. Number of strings per sym
      CALL NST_SPGP(I1SPOBEX(1,1),NI1CA)
      CALL NST_SPGP(I1SPOBEX(1,2),NI1CB)
      CALL NST_SPGP(I1SPOBEX(1,3),NI1AA)
      CALL NST_SPGP(I1SPOBEX(1,4),NI1AB)
*
      CALL NST_SPGP(I2SPOBEX(1,1),NI2CA)
      CALL NST_SPGP(I2SPOBEX(1,2),NI2CB)
      CALL NST_SPGP(I2SPOBEX(1,3),NI2AA)
      CALL NST_SPGP(I2SPOBEX(1,4),NI2AB)
*
      CALL NST_SPGP(I12SPOBEX(1,1),NI12CA)
      CALL NST_SPGP(I12SPOBEX(1,2),NI12CB)
      CALL NST_SPGP(I12SPOBEX(1,3),NI12AA)
      CALL NST_SPGP(I12SPOBEX(1,4),NI12AB)
*. And offsets
C          Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
C     CALL Z_TCC_OFF2(I1OFF,LT1,I1SPOBEX(1,1),I1SPOBEX(1,2),
C    &     I1SPOBEX(1,3),I1SPOBEX(1,4),IT1SM,NSMST)   
C     CALL Z_TCC_OFF2(I2OFF,LT2,I2SPOBEX(1,1),I2SPOBEX(1,2),
C    &     I2SPOBEX(1,3),I2SPOBEX(1,4),IT2SM,NSMST)   
C     CALL Z_TCC_OFF2(I12OFF,LT11,I12SPOBEX(1,1),I12SPOBEX(1,2),
C    &     I12SPOBEX(1,3),I12SPOBEX(1,4),IT12SM,NSMST)   
*
*  ================
*. T1 * T2 mappings 
*  ================
*
*. CA
C?    WRITE(6,*) ' Memory check before STST_TO_ST_MAP'
C?    CALL MEMCHK
C?    WRITE(6,*) ' Memcheck passed '
*
C          STST_TO_ST_MAP(IS1OC,IS2OC,IS12OC,
C    &           IBS1S2,IS1S2_TO_S12,XS1S2_TO_S12,
C    &           IZ,IZSCR,IS1_STR,IS2_STR,IS12_REO)
      CALL STST_TO_ST_MAP(I1SPOBEX(1,1),I2SPOBEX(1,1),
     &     I12SPOBEX(1,1), IB(1,1,1),ICA_MAP,XCA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. CB
      CALL STST_TO_ST_MAP(I1SPOBEX(1,2),I2SPOBEX(1,2),
     &     I12SPOBEX(1,2), IB(1,1,2),ICB_MAP,XCB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. AA ( Jeppe Phase factor due to def of Annistrings must be added)
      CALL STST_TO_ST_MAP(I1SPOBEX(1,3),I2SPOBEX(1,3),
     &     I12SPOBEX(1,3), IB(1,1,3),IAA_MAP,XAA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. AB ( Jeppe Phase factor due to def of Annistrings must be added)
      CALL STST_TO_ST_MAP(I1SPOBEX(1,4),I2SPOBEX(1,4),
     &     I12SPOBEX(1,4), IB(1,1,4),IAB_MAP,XAB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*
C?    CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TT_TM ')
      RETURN 
      END
      SUBROUTINE T1T2_TO_T12(I1SPOBEX,IT1SM,T1,
     &                       I2SPOBEX,IT2SM,T2,
     &                       I12SPOBEX,IT12SM,T12,
     &                       I1OCC,I2OCC,I1REO,
     &                       ICA_MAP,XCA_MAP,ICB_MAP,XCB_MAP,
     &                       IAA_MAP,XAA_MAP,IAB_MAP,XAB_MAP,
     &                       IZ,IZSCR) 
*
* Obtain excitation operator T12 as a product of two excitation operators  
* T1 and T2. 
*
* All creation and annihilation operators are assumed to commute.
*
* Note : T12 is not set to zero at start  
* Jeppe Olsen, May 1, 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
*
*. Input
*
      INTEGER I1SPOBEX(NGAS,4)
      INTEGER I2SPOBEX(NGAS,4)
      DIMENSION T1(*),T2(*)
*
*. Output
*
      INTEGER I12SPOBEX(NGAS,4)
      DIMENSION T12(*)
*
*. Local scratch 
*
      INTEGER NI1CA(MXPNGAS),NI1CB(MXPNGAS),
     &        NI1AA(MXPNGAS),NI1AB(MXPNGAS)
      INTEGER NI2CA(MXPNGAS),NI2CB(MXPNGAS),
     &        NI2AA(MXPNGAS),NI2AB(MXPNGAS)
      INTEGER NI12CA(MXPNGAS),NI12CB(MXPNGAS),
     &        NI12AA(MXPNGAS),NI12AB(MXPNGAS)
      INTEGER I1CA_EXP(100),I1CB_EXP(100),I1AA_EXP(100),I1AB_EXP(100)
      INTEGER I2CA_EXP(100),I2CB_EXP(100),I2AA_EXP(100),I2AB_EXP(100)
      INTEGER I12CA_EXP(100),I12CB_EXP(100),
     &        I12AA_EXP(100),I12AB_EXP(100)
C     INTEGER ICA(100)
      INTEGER I1OFF(8,8,8), I2OFF(8,8,8), I12OFF(8,8,8)
      INTEGER IB_CA_SS(8,8),IB_CB_SS(8,8),IB_AA_SS(8,8),IB_AB_SS(8,8)
*
*. Scratch through parameter list
*
*. Space for String*String => String maps      
      INTEGER ICA_MAP(*),ICB_MAP(*),IAA_MAP(*),IAB_MAP(*)               
      DIMENSION XCA_MAP(*),XCB_MAP(*),XAA_MAP(*),XAB_MAP(*)               
*. For occ of two sets of excitation strings and an reorder array 
      INTEGER I1OCC(*),I2OCC(*),I1REO(*)
*. For Z-matrix and its construction
      INTEGER IZ(*),IZSCR(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'T1T2_T')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' Welcome to T1T2_TO_T12 '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Excitation operator 1 :'
        WRITE(6,*)
C   WRT_SPOX_TP(IEX_TP,NEX_TP)
        CALL WRT_SPOX_TP(I1SPOBEX,1)
        WRITE(6,*) ' Excitation operator 2 :'
        WRITE(6,*)
        CALL WRT_SPOX_TP(I2SPOBEX,1)
        WRITE(6,*) ' Symmetry of T1 and T2 : ', IT1SM,IT2SM
      END IF
*. Type and symmetry of compound excitation operators
      IONE = 1
      DO JPART = 1, 4
        CALL IVCSUM(I12SPOBEX(1,JPART),I1SPOBEX(1,JPART),
     &              I2SPOBEX(1,JPART),IONE,IONE,NGAS) 
      END DO 
      IT12SM = MULTD2H(IT1SM,IT2SM)
      IF(NTEST.GE.1000) THEN 
        WRITE(6,*)
        WRITE(6,*) ' Compound operator '
        CALL WRT_SPOX_TP(I12SPOBEX,1)
        WRITE(6,*)
        WRITE(6,*) ' Symmetry of compound operator', IT12SM
      END IF
*. Number of operators in each string 
      N1CA_OP = IELSUM(I1SPOBEX(1,1),NGAS)
      N1CB_OP = IELSUM(I1SPOBEX(1,2),NGAS)
      N1AA_OP = IELSUM(I1SPOBEX(1,3),NGAS)
      N1AB_OP = IELSUM(I1SPOBEX(1,4),NGAS)
*
      N2CA_OP = IELSUM(I2SPOBEX(1,1),NGAS)
      N2CB_OP = IELSUM(I2SPOBEX(1,2),NGAS)
      N2AA_OP = IELSUM(I2SPOBEX(1,3),NGAS)
      N2AB_OP = IELSUM(I2SPOBEX(1,4),NGAS)
*
      N12CA_OP = IELSUM(I12SPOBEX(1,1),NGAS)
      N12CB_OP = IELSUM(I12SPOBEX(1,2),NGAS)
      N12AA_OP = IELSUM(I12SPOBEX(1,3),NGAS)
      N12AB_OP = IELSUM(I12SPOBEX(1,4),NGAS)
*. Sign to bring O(ca1)O(cb1)O(aa1)O(ab1)
*.               O(ca2)O(cb2)O(aa2)O(ab2)
*. into          O(ca1)O(ca2)O(cb1)O(cb2)O(aa1)O(aa2)O(ab1)O(ab2)
*
      NPERM =  N1AB_OP*( N2CA_OP+N2CB_OP+N2AA_OP)
     &      +  N1AA_OP*( N2CA_OP+N2CB_OP )
     &      +  N1CB_OP*  N2CA_OP
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNP = -1.0D0
      ELSE
        SIGNP = 1.0D0
      END IF
*
      MX_NOP = MAX(N12CA_OP,N12CB_OP,N12AA_OP,N12AB_OP)
*
C     CALL REF_OP(ICA,ICA_EXP,NCA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,1),I1CA_EXP,N1CA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,2),I1CB_EXP,N1CB_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,3),I1AA_EXP,N1AA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,4),I1AB_EXP,N1AB_OP,NGAS,1)
*
      CALL REF_OP(I2SPOBEX(1,1),I2CA_EXP,N2CA_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,2),I2CB_EXP,N2CB_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,3),I2AA_EXP,N2AA_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,4),I2AB_EXP,N2AB_OP,NGAS,1)
*
      CALL REF_OP(I12SPOBEX(1,1),I12CA_EXP,N12CA_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,2),I12CB_EXP,N12CB_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,3),I12AA_EXP,N12AA_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,4),I12AB_EXP,N12AB_OP,NGAS,1)
*. Number of strings per sym
      CALL NST_SPGP(I1SPOBEX(1,1),NI1CA)
      CALL NST_SPGP(I1SPOBEX(1,2),NI1CB)
      CALL NST_SPGP(I1SPOBEX(1,3),NI1AA)
      CALL NST_SPGP(I1SPOBEX(1,4),NI1AB)
*
      CALL NST_SPGP(I2SPOBEX(1,1),NI2CA)
      CALL NST_SPGP(I2SPOBEX(1,2),NI2CB)
      CALL NST_SPGP(I2SPOBEX(1,3),NI2AA)
      CALL NST_SPGP(I2SPOBEX(1,4),NI2AB)
*
      CALL NST_SPGP(I12SPOBEX(1,1),NI12CA)
      CALL NST_SPGP(I12SPOBEX(1,2),NI12CB)
      CALL NST_SPGP(I12SPOBEX(1,3),NI12AA)
      CALL NST_SPGP(I12SPOBEX(1,4),NI12AB)
*. And offsets
C          Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
      CALL Z_TCC_OFF2(I1OFF,LT1,NI1CA,NI1CB,NI1AA,NI1AB,IT1SM,NSMST)   
      CALL Z_TCC_OFF2(I2OFF,LT2,NI2CA,NI2CB,NI2AA,NI2AB,IT2SM,NSMST)   
      CALL Z_TCC_OFF2(I12OFF,LT12,NI12CA,NI12CB,NI12AA,NI12AB,
     &                IT12SM,NSMST)   
*
*  ================
*. T1 * T2 mappings 
*  ================
*
*. CA
*
C          STST_TO_ST_MAP(IS1OC,IS2OC,IS12OC,
C    &           IBS1S2,IS1S2_TO_S12,XS1S2_TO_S12,
C    &           IZ,IZSCR,IS1_STR,IS2_STR,IS12_REO)
      CALL STST_TO_ST_MAP(I1SPOBEX(1,1),I2SPOBEX(1,1),
     &     I12SPOBEX(1,1), IB_CA_SS,ICA_MAP,XCA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. CB
      CALL STST_TO_ST_MAP(I1SPOBEX(1,2),I2SPOBEX(1,2),
     &     I12SPOBEX(1,2), IB_CB_SS,ICB_MAP,XCB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. AA 
      CALL STST_TO_ST_MAP(I1SPOBEX(1,3),I2SPOBEX(1,3),
     &     I12SPOBEX(1,3), IB_AA_SS,IAA_MAP,XAA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. AB 
      CALL STST_TO_ST_MAP(I1SPOBEX(1,4),I2SPOBEX(1,4),
     &     I12SPOBEX(1,4), IB_AB_SS,IAB_MAP,XAB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. And then the mapping
      DO I12SM_C = 1, NSMST
       I12SM_A = MULTD2H(IT12SM,I12SM_C) 
       DO I12SM_CA = 1, NSMST
        I12SM_CB = MULTD2H(I12SM_C,I12SM_CA)
        DO I12SM_AA = 1, NSMST
         I12SM_AB =  MULTD2H(I12SM_A,I12SM_AA)
         DO I1SM_C = 1, NSMST
          I1SM_A = MULTD2H(IT1SM,I1SM_C) 
          DO I1SM_CA = 1, NSMST
           I1SM_CB = MULTD2H(I1SM_C,I1SM_CA)
           DO I1SM_AA = 1, NSMST
            I1SM_AB =  MULTD2H(I1SM_A,I1SM_AA)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' I12SM_CA, I12SM_CB, I12SM_AA, I12SM_AB ',
     &                     I12SM_CA, I12SM_CB, I12SM_AA, I12SM_AB
              WRITE(6,*) ' I1SM_CA, I1SM_CB, I1SM_AA, I1SM_AB ',
     &                     I1SM_CA, I1SM_CB, I1SM_AA, I1SM_AB
            END IF
            I2SM_CA = MULTD2H(I1SM_CA,I12SM_CA)
            I2SM_CB = MULTD2H(I1SM_CB,I12SM_CB)
            I2SM_AA = MULTD2H(I1SM_AA,I12SM_AA)
            I2SM_AB = MULTD2H(I1SM_AB,I12SM_AB)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' I2SM_CA, I2SM_CB, I2SM_AA, I2SM_AB ',
     &                     I2SM_CA, I2SM_CB, I2SM_AA, I2SM_AB
            END IF
            NI1CA_FSM = NI1CA(I1SM_CA)
            NI1CB_FSM = NI1CB(I1SM_CB)
            NI1AA_FSM = NI1AA(I1SM_AA)
            NI1AB_FSM = NI1AB(I1SM_AB)
*
            NI2CA_FSM = NI2CA(I2SM_CA)
            NI2CB_FSM = NI2CB(I2SM_CB)
            NI2AA_FSM = NI2AA(I2SM_AA)
            NI2AB_FSM = NI2AB(I2SM_AB)
*
            NI12CA_FSM = NI12CA(I12SM_CA)
            NI12CB_FSM = NI12CB(I12SM_CB)
            NI12AA_FSM = NI12AA(I12SM_AA)
            NI12AB_FSM = NI12AB(I12SM_AB)
* . Offsets in st*st => st mappings
            IB_CA = IB_CA_SS(I1SM_CA, I2SM_CA)
            IB_CB = IB_CB_SS(I1SM_CB, I2SM_CB)
            IB_AA = IB_AA_SS(I1SM_AA, I2SM_AA)
            IB_AB = IB_AB_SS(I1SM_AB, I2SM_AB)
*. Offsets in T1,T2, T12
            J1OFF = I1OFF(I1SM_CA,I1SM_CB,I1SM_AA)
            J2OFF = I2OFF(I2SM_CA,I2SM_CB,I2SM_AA)
            J12OFF = I12OFF(I12SM_CA,I12SM_CB,I12SM_AA)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' J1OFF, J2OFF, J12OFF = ', 
     &                     J1OFF, J2OFF, J12OFF
              WRITE(6,*) 
     &        ' NI1CA_FSM,  NI2CA_FSM,  NI1CB_FSM, NI2CB_FSM = ',
     &          NI1CA_FSM,  NI2CA_FSM,  NI1CB_FSM, NI2CB_FSM 
            END IF
*. Loop over the four pair of operators ( 8 more indeces !)
            DO I1CA_ST = 1, NI1CA_FSM
            DO I2CA_ST = 1, NI2CA_FSM
             DO I1CB_ST = 1, NI1CB_FSM
             DO I2CB_ST = 1, NI2CB_FSM
              DO I1AA_ST = 1, NI1AA_FSM
              DO I2AA_ST = 1, NI2AA_FSM
               DO I1AB_ST = 1, NI1AB_FSM
               DO I2AB_ST = 1, NI2AB_FSM

                I1_CAAB = J1OFF - 1
     &          +(I1AB_ST-1)*NI1AA_FSM*NI1CA_FSM*NI1CB_FSM
     &          +(I1AA_ST-1)*NI1CA_FSM*NI1CB_FSM
     &          +(I1CB_ST-1)*NI1CA_FSM + I1CA_ST

                I2_CAAB = J2OFF - 1
     &          +(I2AB_ST-1)*NI2AA_FSM*NI2CA_FSM*NI2CB_FSM
     &          +(I2AA_ST-1)*NI2CA_FSM*NI2CB_FSM
     &          +(I2CB_ST-1)*NI2CA_FSM + I2CA_ST
*
                K12CA = IB_CA - 1 + (I2CA_ST-1)*NI1CA_FSM + I1CA_ST
                K12CB = IB_CB - 1 + (I2CB_ST-1)*NI1CB_FSM + I1CB_ST
                K12AA = IB_AA - 1 + (I2AA_ST-1)*NI1AA_FSM + I1AA_ST
                K12AB = IB_AB - 1 + (I2AB_ST-1)*NI1AB_FSM + I1AB_ST
*
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' IB_CA, I2CA_ST, I1CA_ST, NI1CA_FSM',
     &                         IB_CA, I2CA_ST, I1CA_ST, NI1CA_FSM
                  WRITE(6,*) ' K12CA = ', K12CA
                END IF
*
                I12ST_CA = ICA_MAP(K12CA)
C?              WRITE(6,*) ' ICA_MAP(1) = ', ICA_MAP(1)
C?              WRITE(6,*) ' Memory check '
C?              CALL MEMCHK
                I12ST_CB = ICB_MAP(K12CB)
                I12ST_AA = IAA_MAP(K12AA)
                I12ST_AB = IAB_MAP(K12AB)
                IF(I12ST_CA.NE.0.AND.I12ST_CB.NE.0.AND.
     &             I12ST_AA.NE.0.AND.I12ST_AB.NE.0     ) THEN
*
                 X12_CA = XCA_MAP(K12CA)
                 X12_CB = XCB_MAP(K12CB)
                 X12_AA = XAA_MAP(K12AA)
                 X12_AB = XAB_MAP(K12AB)
                 SIGN = X12_CA*X12_CB*X12_AA*X12_AB*SIGNP
*           
                 I12_CAAB = J12OFF-1 
     &          +(I12ST_AB-1)*NI12AA_FSM*NI12CB_FSM*NI12CA_FSM
     &          +(I12ST_AA-1)*NI12CB_FSM*NI12CA_FSM
     &          +(I12ST_CB-1)*NI12CA_FSM + I12ST_CA
*
                 IF(NTEST.GE.1000) THEN
                  WRITE(6,*) 'I12ST_AB, I12ST_AA, I12ST_CB, I12ST_CA=',
     &                        I12ST_AB, I12ST_AA, I12ST_CB, I12ST_CA
  
                  WRITE(6,*) 'I1_CAAB, I2_CAAB, I12_CAAB, sign = ',
     &                        I1_CAAB, I2_CAAB, I12_CAAB, SIGN 
                 END IF
                 T12(I12_CAAB) = 
     &           T12(I12_CAAB) + SIGN*T1(I1_CAAB)*T2(I2_CAAB)
                 IF(NTEST.GE.1000) 
     &           WRITE(6,*) ' T12(I12_CAAB) = ',  T12(I12_CAAB)
                END IF
               END DO
               END DO
*              ^ End of loop over AB strings
              END DO
              END DO
*             ^ End of loop over AA strings
             END DO
             END DO
*            ^ End of loop over CB strings
            END DO
            END DO
*           ^ End of loop over CA strings
           END DO
*          ^ End of loop over I1SM_AA
          END DO
*         ^ End of loop over I1SM_CA
         END DO
*        ^ End of loop over I1SM_C
        END DO
*       ^ End of loop over I1SM_AA
       END DO
*       ^ End of loop over I12SM_CA
      END DO
*     ^ End of loop over I12SM_C
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Product of two excitation operators (T1T2_TO_T12)'
        WRITE(6,*) ' ================================================='
        WRITE(6,*)  
        WRITE(6,*) ' T1 TCC block '
C            WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
        CALL WRT_TCC_BLK(T1,IT1SM,NI1CA,NI1CB,NI1AA,NI1AB,NSMST)
        WRITE(6,*) ' T2 TCC block '
        CALL WRT_TCC_BLK(T2,IT2SM,NI2CA,NI2CB,NI2AA,NI2AB,NSMST)
        WRITE(6,*) ' T12 TCC block '
        CALL WRT_TCC_BLK(T12,IT12SM,NI12CA,NI12CB,NI12AA,NI12AB,NSMST)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'T1T2_T')
      RETURN 
      END
      SUBROUTINE T1T2_TO_T12N(I1SPOBEX,IT1SM,T1,
     &                       I2SPOBEX,IT2SM,T2,
     &                       I12SPOBEX,IT12SM,T12,
     &                       I1OCC,I2OCC,I1REO,
     &                       ICA_MAP,XCA_MAP,ICB_MAP,XCB_MAP,
     &                       IAA_MAP,XAA_MAP,IAB_MAP,XAB_MAP,
     &                       IZ,IZSCR) 
*
* Obtain excitation operator T12 as a product of two excitation operators  
* T1 and T2. 
*
* All creation and annihilation operators are assumed to commute.
*
* Note : T12 is not set to zero at start  
* Jeppe Olsen, May 1, 2000
*              Speeded up august 2001
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
*
*. Input
*
      INTEGER I1SPOBEX(NGAS,4)
      INTEGER I2SPOBEX(NGAS,4)
      DIMENSION T1(*),T2(*)
*
*. Output
*
      INTEGER I12SPOBEX(NGAS,4)
      DIMENSION T12(*)
*
*. Local scratch 
*
      INTEGER NI1CA(MXPNGAS),NI1CB(MXPNGAS),
     &        NI1AA(MXPNGAS),NI1AB(MXPNGAS)
      INTEGER NI2CA(MXPNGAS),NI2CB(MXPNGAS),
     &        NI2AA(MXPNGAS),NI2AB(MXPNGAS)
      INTEGER NI12CA(MXPNGAS),NI12CB(MXPNGAS),
     &        NI12AA(MXPNGAS),NI12AB(MXPNGAS)
      INTEGER I1CA_EXP(100),I1CB_EXP(100),I1AA_EXP(100),I1AB_EXP(100)
      INTEGER I2CA_EXP(100),I2CB_EXP(100),I2AA_EXP(100),I2AB_EXP(100)
      INTEGER I12CA_EXP(100),I12CB_EXP(100),
     &        I12AA_EXP(100),I12AB_EXP(100)
C     INTEGER ICA(100)
      INTEGER I1OFF(8,8,8), I2OFF(8,8,8), I12OFF(8,8,8)
      INTEGER IB_CA_SS(8,8),IB_CB_SS(8,8),IB_AA_SS(8,8),IB_AB_SS(8,8)
*
*. Scratch through parameter list
*
*. Space for String*String => String maps      
      INTEGER ICA_MAP(*),ICB_MAP(*),IAA_MAP(*),IAB_MAP(*)               
      DIMENSION XCA_MAP(*),XCB_MAP(*),XAA_MAP(*),XAB_MAP(*)               
*. For occ of two sets of excitation strings and an reorder array 
      INTEGER I1OCC(*),I2OCC(*),I1REO(*)
*. For Z-matrix and its construction
      INTEGER IZ(*),IZSCR(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'T1T2_T')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' Welcome to T1T2_TO_T12 '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Excitation operator 1 :'
        WRITE(6,*)
C   WRT_SPOX_TP(IEX_TP,NEX_TP)
        CALL WRT_SPOX_TP(I1SPOBEX,1)
        WRITE(6,*) ' Excitation operator 2 :'
        WRITE(6,*)
        CALL WRT_SPOX_TP(I2SPOBEX,1)
        WRITE(6,*) ' Symmetry of T1 and T2 : ', IT1SM,IT2SM
      END IF
*. Type and symmetry of compound excitation operators
      IONE = 1
      DO JPART = 1, 4
        CALL IVCSUM(I12SPOBEX(1,JPART),I1SPOBEX(1,JPART),
     &              I2SPOBEX(1,JPART),IONE,IONE,NGAS) 
      END DO 
      IT12SM = MULTD2H(IT1SM,IT2SM)
      IF(NTEST.GE.1000) THEN 
        WRITE(6,*)
        WRITE(6,*) ' Compound operator '
        CALL WRT_SPOX_TP(I12SPOBEX,1)
        WRITE(6,*)
        WRITE(6,*) ' Symmetry of compound operator', IT12SM
      END IF
*. Number of operators in each string 
      N1CA_OP = IELSUM(I1SPOBEX(1,1),NGAS)
      N1CB_OP = IELSUM(I1SPOBEX(1,2),NGAS)
      N1AA_OP = IELSUM(I1SPOBEX(1,3),NGAS)
      N1AB_OP = IELSUM(I1SPOBEX(1,4),NGAS)
*
      N2CA_OP = IELSUM(I2SPOBEX(1,1),NGAS)
      N2CB_OP = IELSUM(I2SPOBEX(1,2),NGAS)
      N2AA_OP = IELSUM(I2SPOBEX(1,3),NGAS)
      N2AB_OP = IELSUM(I2SPOBEX(1,4),NGAS)
*
      N12CA_OP = IELSUM(I12SPOBEX(1,1),NGAS)
      N12CB_OP = IELSUM(I12SPOBEX(1,2),NGAS)
      N12AA_OP = IELSUM(I12SPOBEX(1,3),NGAS)
      N12AB_OP = IELSUM(I12SPOBEX(1,4),NGAS)
*. Sign to bring O(ca1)O(cb1)O(aa1)O(ab1)
*.               O(ca2)O(cb2)O(aa2)O(ab2)
*. into          O(ca1)O(ca2)O(cb1)O(cb2)O(aa1)O(aa2)O(ab1)O(ab2)
*
      NPERM =  N1AB_OP*( N2CA_OP+N2CB_OP+N2AA_OP)
     &      +  N1AA_OP*( N2CA_OP+N2CB_OP )
     &      +  N1CB_OP*  N2CA_OP
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNP = -1.0D0
      ELSE
        SIGNP = 1.0D0
      END IF
*
      MX_NOP = MAX(N12CA_OP,N12CB_OP,N12AA_OP,N12AB_OP)
*
C     CALL REF_OP(ICA,ICA_EXP,NCA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,1),I1CA_EXP,N1CA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,2),I1CB_EXP,N1CB_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,3),I1AA_EXP,N1AA_OP,NGAS,1)
      CALL REF_OP(I1SPOBEX(1,4),I1AB_EXP,N1AB_OP,NGAS,1)
*
      CALL REF_OP(I2SPOBEX(1,1),I2CA_EXP,N2CA_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,2),I2CB_EXP,N2CB_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,3),I2AA_EXP,N2AA_OP,NGAS,1)
      CALL REF_OP(I2SPOBEX(1,4),I2AB_EXP,N2AB_OP,NGAS,1)
*
      CALL REF_OP(I12SPOBEX(1,1),I12CA_EXP,N12CA_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,2),I12CB_EXP,N12CB_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,3),I12AA_EXP,N12AA_OP,NGAS,1)
      CALL REF_OP(I12SPOBEX(1,4),I12AB_EXP,N12AB_OP,NGAS,1)
*. Number of strings per sym
      CALL NST_SPGP(I1SPOBEX(1,1),NI1CA)
      CALL NST_SPGP(I1SPOBEX(1,2),NI1CB)
      CALL NST_SPGP(I1SPOBEX(1,3),NI1AA)
      CALL NST_SPGP(I1SPOBEX(1,4),NI1AB)
*
      CALL NST_SPGP(I2SPOBEX(1,1),NI2CA)
      CALL NST_SPGP(I2SPOBEX(1,2),NI2CB)
      CALL NST_SPGP(I2SPOBEX(1,3),NI2AA)
      CALL NST_SPGP(I2SPOBEX(1,4),NI2AB)
*
      CALL NST_SPGP(I12SPOBEX(1,1),NI12CA)
      CALL NST_SPGP(I12SPOBEX(1,2),NI12CB)
      CALL NST_SPGP(I12SPOBEX(1,3),NI12AA)
      CALL NST_SPGP(I12SPOBEX(1,4),NI12AB)
*. And offsets
C          Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
      CALL Z_TCC_OFF2(I1OFF,LT1,NI1CA,NI1CB,NI1AA,NI1AB,IT1SM,NSMST)   
      CALL Z_TCC_OFF2(I2OFF,LT2,NI2CA,NI2CB,NI2AA,NI2AB,IT2SM,NSMST)   
      CALL Z_TCC_OFF2(I12OFF,LT12,NI12CA,NI12CB,NI12AA,NI12AB,
     &                IT12SM,NSMST)   
*
*  ================
*. T1 * T2 mappings 
*  ================
*
*. CA
*
C          STST_TO_ST_MAP(IS1OC,IS2OC,IS12OC,
C    &           IBS1S2,IS1S2_TO_S12,XS1S2_TO_S12,
C    &           IZ,IZSCR,IS1_STR,IS2_STR,IS12_REO)
      CALL STST_TO_ST_MAP(I1SPOBEX(1,1),I2SPOBEX(1,1),
     &     I12SPOBEX(1,1), IB_CA_SS,ICA_MAP,XCA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. CB
      CALL STST_TO_ST_MAP(I1SPOBEX(1,2),I2SPOBEX(1,2),
     &     I12SPOBEX(1,2), IB_CB_SS,ICB_MAP,XCB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. AA 
      CALL STST_TO_ST_MAP(I1SPOBEX(1,3),I2SPOBEX(1,3),
     &     I12SPOBEX(1,3), IB_AA_SS,IAA_MAP,XAA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. AB 
      CALL STST_TO_ST_MAP(I1SPOBEX(1,4),I2SPOBEX(1,4),
     &     I12SPOBEX(1,4), IB_AB_SS,IAB_MAP,XAB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO)
*. And then the mapping
      DO I12SM_C = 1, NSMST
       I12SM_A = MULTD2H(IT12SM,I12SM_C) 
       DO I12SM_CA = 1, NSMST
        I12SM_CB = MULTD2H(I12SM_C,I12SM_CA)
        NI12CA_FSM = NI12CA(I12SM_CA)
        NI12CB_FSM = NI12CB(I12SM_CB)
        IF( NI12CA_FSM*NI12CB_FSM .NE. 0 ) THEN
        DO I12SM_AA = 1, NSMST
         I12SM_AB =  MULTD2H(I12SM_A,I12SM_AA)
         NI12AA_FSM = NI12AA(I12SM_AA)
         NI12AB_FSM = NI12AB(I12SM_AB)
         IF( NI12AA_FSM*NI12AB_FSM.NE.0 ) THEN
         DO I1SM_C = 1, NSMST
          I1SM_A = MULTD2H(IT1SM,I1SM_C) 
          DO I1SM_CA = 1, NSMST
           I1SM_CB = MULTD2H(I1SM_C,I1SM_CA)
           DO I1SM_AA = 1, NSMST
            I1SM_AB =  MULTD2H(I1SM_A,I1SM_AA)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' I12SM_CA, I12SM_CB, I12SM_AA, I12SM_AB ',
     &                     I12SM_CA, I12SM_CB, I12SM_AA, I12SM_AB
              WRITE(6,*) ' I1SM_CA, I1SM_CB, I1SM_AA, I1SM_AB ',
     &                     I1SM_CA, I1SM_CB, I1SM_AA, I1SM_AB
            END IF
            I2SM_CA = MULTD2H(I1SM_CA,I12SM_CA)
            I2SM_CB = MULTD2H(I1SM_CB,I12SM_CB)
            I2SM_AA = MULTD2H(I1SM_AA,I12SM_AA)
            I2SM_AB = MULTD2H(I1SM_AB,I12SM_AB)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' I2SM_CA, I2SM_CB, I2SM_AA, I2SM_AB ',
     &                     I2SM_CA, I2SM_CB, I2SM_AA, I2SM_AB
            END IF
            NI1CA_FSM = NI1CA(I1SM_CA)
            NI1CB_FSM = NI1CB(I1SM_CB)
            NI1AA_FSM = NI1AA(I1SM_AA)
            NI1AB_FSM = NI1AB(I1SM_AB)
            L1 =  NI1CA_FSM*NI1CB_FSM*NI1AA_FSM*NI1AB_FSM
*
            NI2CA_FSM = NI2CA(I2SM_CA)
            NI2CB_FSM = NI2CB(I2SM_CB)
            NI2AA_FSM = NI2AA(I2SM_AA)
            NI2AB_FSM = NI2AB(I2SM_AB)
            L2 =  NI2CA_FSM*NI2CB_FSM*NI2AA_FSM*NI2AB_FSM
*
            NI12CA_FSM = NI12CA(I12SM_CA)
            NI12CB_FSM = NI12CB(I12SM_CB)
            NI12AA_FSM = NI12AA(I12SM_AA)
            NI12AB_FSM = NI12AB(I12SM_AB)
            L12 =  NI12CA_FSM*NI12CB_FSM*NI12AA_FSM*NI12AB_FSM
            IF(L1*L2*L12.NE.0) THEN
* . Offsets in st*st => st mappings
            IB_CA = IB_CA_SS(I1SM_CA, I2SM_CA)
            IB_CB = IB_CB_SS(I1SM_CB, I2SM_CB)
            IB_AA = IB_AA_SS(I1SM_AA, I2SM_AA)
            IB_AB = IB_AB_SS(I1SM_AB, I2SM_AB)
*. Offsets in T1,T2, T12
            J1OFF = I1OFF(I1SM_CA,I1SM_CB,I1SM_AA)
            J2OFF = I2OFF(I2SM_CA,I2SM_CB,I2SM_AA)
            J12OFF = I12OFF(I12SM_CA,I12SM_CB,I12SM_AA)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' J1OFF, J2OFF, J12OFF = ', 
     &                     J1OFF, J2OFF, J12OFF
              WRITE(6,*) 
     &        ' NI1CA_FSM,  NI2CA_FSM,  NI1CB_FSM, NI2CB_FSM = ',
     &          NI1CA_FSM,  NI2CA_FSM,  NI1CB_FSM, NI2CB_FSM 
            END IF
*. Loop over the four pair of operators ( 8 more indeces !)
            DO I1AB_ST = 1, NI1AB_FSM
            DO I2AB_ST = 1, NI2AB_FSM
             K12AB = IB_AB - 1 + (I2AB_ST-1)*NI1AB_FSM + I1AB_ST
             I12ST_AB = IAB_MAP(K12AB)
             X12_AB = XAB_MAP(K12AB)
             DO I1AA_ST = 1, NI1AA_FSM
             DO I2AA_ST = 1, NI2AA_FSM
               K12AA = IB_AA - 1 + (I2AA_ST-1)*NI1AA_FSM + I1AA_ST
               I12ST_AA = IAA_MAP(K12AA)
               X12_AA = XAA_MAP(K12AA)
               DO I1CB_ST = 1, NI1CB_FSM
               I1_CAAB0 = J1OFF - 1
     &         +(I1AB_ST-1)*NI1AA_FSM*NI1CA_FSM*NI1CB_FSM
     &         +(I1AA_ST-1)*NI1CA_FSM*NI1CB_FSM
     &         +(I1CB_ST-1)*NI1CA_FSM 
               DO I2CB_ST = 1, NI2CB_FSM
*
                K12CB = IB_CB - 1 + (I2CB_ST-1)*NI1CB_FSM + I1CB_ST
                I12ST_CB = ICB_MAP(K12CB)
                X12_CB = XCB_MAP(K12CB)
                IF(I12ST_AA.NE.0.AND.I12ST_AB.NE.0.AND.
     &             I12ST_CB.NE.0) THEN
*
                I2_CAAB0 = J2OFF - 1
     &          +(I2AB_ST-1)*NI2AA_FSM*NI2CA_FSM*NI2CB_FSM
     &          +(I2AA_ST-1)*NI2CA_FSM*NI2CB_FSM
     &          +(I2CB_ST-1)*NI2CA_FSM 
*
                I12_CAAB0 = J12OFF-1 
     &          +(I12ST_AB-1)*NI12AA_FSM*NI12CB_FSM*NI12CA_FSM
     &          +(I12ST_AA-1)*NI12CB_FSM*NI12CA_FSM
     &          +(I12ST_CB-1)*NI12CA_FSM 
        
                DO I1CA_ST = 1, NI1CA_FSM
                 I1_CAAB = I1_CAAB0 + I1CA_ST
                 K12CA = IB_CA - 1 - NI1CA_FSM + I1CA_ST
                 CONST = X12_CB*X12_AA*X12_AB*SIGNP*T1(I1_CAAB)
                 
                DO I2CA_ST = 1, NI2CA_FSM
C                K12CA = IB_CA - 1 + (I2CA_ST-1)*NI1CA_FSM + I1CA_ST
                 K12CA = K12CA + NI1CA_FSM
                 I12ST_CA = ICA_MAP(K12CA)
                 IF(I12ST_CA.NE.0) THEN
*
                  X12_CA = XCA_MAP(K12CA)
                  I2_CAAB = I2_CAAB0 + I2CA_ST
                  I12_CAAB = I12_CAAB0 + I12ST_CA
C                 SIGN = X12_CA*X12_CB*X12_AA*X12_AB*SIGNP
*
                  T12(I12_CAAB) = 
     &            T12(I12_CAAB) + CONST*X12_CA*T2(I2_CAAB)
*
                 END IF
*               ^ End of CA map was nonvanishing
                END DO
                END DO
*               ^ End of loop over CA strings
              END IF 
*             ^ End if AA,AB,CB map was nonvanishing
              END DO
              END DO
*             ^ End of loop over CB strings
             END DO
             END DO
*            ^ End of loop over AA strings
            END DO
            END DO
*           ^ End of loop over AB strings
           END IF
*          ^ End if dimensions were nonvanishing
           END DO
*          ^ End of loop over I1SM_AA
          END DO
*         ^ End of loop over I1SM_CA
         END DO
*        ^ End of loop over I1SM_C
         END IF
*        ^ End if where was nonvanishing I12_A strings
        END DO
*       ^ End of loop over I12SM_AA
        END IF
*       ^ End of there was nonvanishing I12_C strings with these sym
       END DO
*       ^ End of loop over I12SM_CA
      END DO
*     ^ End of loop over I12SM_C
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Product of two excitation operators (T1T2_TO_T12)'
        WRITE(6,*) ' ================================================='
        WRITE(6,*)  
        WRITE(6,*) ' T1 TCC block '
C            WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
        CALL WRT_TCC_BLK(T1,IT1SM,NI1CA,NI1CB,NI1AA,NI1AB,NSMST)
        WRITE(6,*) ' T2 TCC block '
        CALL WRT_TCC_BLK(T2,IT2SM,NI2CA,NI2CB,NI2AA,NI2AB,NSMST)
        WRITE(6,*) ' T12 TCC block '
        CALL WRT_TCC_BLK(T12,IT12SM,NI12CA,NI12CB,NI12AA,NI12AB,NSMST)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'T1T2_T')
      RETURN 
      END
      SUBROUTINE STST_TO_ST_MAP(IS1OC,IS2OC,IS12OC,
     &           IBS1S2,IS1S2_TO_S12,XS1S2_TO_S12,
     &           IZ,IZSCR,IS1_STR,IS2_STR,IS12_REO)
*
* Consider string multiplication S1*S2 => S12 
* Find the mapping between strings             
* Maps are organized as matrices (S1,S2) for each symmetry of S1 and S2
*
* Jeppe Olsen, May 1, 2000
*
* New version reducing number of calls to GET_STR2_TOTSM..
* a rainy day at HNIE, July 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input :  number of operators in each gasspace
C     INTEGER IS1OP(NGAS),IS2OP(NGAS),IS12OP(NGAS)
      INTEGER IS1OC(*),IS2OC(*),IS12OC(*)
*
*. Local scratch
*
      INTEGER IS1GRP(NGAS),IS2GRP(NGAS),IS12GRP(NGAS)
      INTEGER ISTR_OUT(100)
*
*. Scratch space through parameter list
*
      INTEGER IZ(*),IZSCR(*)
      INTEGER IS1_STR(*),IS2_STR(*),IS12_REO(*)
*     ^ Should hold largest list of strings with given CAAB
*. Local scratch 
      INTEGER LS1_STR(8) , IBS1_STR(8)
      INTEGER LS2_STR(8) , IBS2_STR(8)
C     INTEGER LS12_STR(8), IBS12_STR(8)
*
*. Output :
*
*. Offset to mappings for given sym of S1 and S2
      INTEGER IBS1S2(8,8)
*. And the mappings 
      INTEGER IS1S2_TO_S12(*)
      DIMENSION XS1S2_TO_S12(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' STST_TO_ST_MAP entered '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of S1, S2 and S12 '
        CALL IWRTMA(IS1OC,1,NGAS,1,NGAS)
        CALL IWRTMA(IS2OC,1,NGAS,1,NGAS)
        CALL IWRTMA(IS12OC,1,NGAS,1,NGAS)
      END IF
*
*. Occupation to group notation
      CALL OCC_TO_GRP(IS1OC,IS1GRP,1)
      CALL OCC_TO_GRP(IS2OC,IS2GRP,1)
      CALL OCC_TO_GRP(IS12OC,IS12GRP,1)
*. Number of operators 
      NS1OP = IELSUM(IS1OC,NGAS)
      NS2OP = IELSUM(IS2OC,NGAS)
      NS12OP = IELSUM(IS12OC,NGAS)
*. Z array for 12 strings 
      CALL WEIGHT_SPGP(IZ,NGAS,IS12OC,NOBPT,IZSCR,0)
*. Set up I12 Reorder array for all symmetries
      IOFF = 1
      DO IS12SM = 1, NSMST
        CALL GETSTR2_TOTSM_SPGP(IS12GRP,NGAS,IS12SM,NS12OP,NS12STR,
     &       IS1_STR,NOCOB,1,IZ,IS12_REO)
      END DO
*. Set up all I1 strings
      IOFF = 1
      DO IS1SM = 1, NSMST
        IBS1_STR(IS1SM) = IOFF
        CALL GETSTR2_TOTSM_SPGP(IS1GRP,NGAS,IS1SM,NS1OP,NS1STR,
     &       IS1_STR(IOFF),NOCOB,0,0,0)              
        LS1_STR(IS1SM) = NS1STR
        IOFF = IOFF + NS1STR*NS1OP
      END DO
*. And I2 strings
      IOFF = 1
      DO IS2SM = 1, NSMST
        IBS2_STR(IS2SM) = IOFF
        CALL GETSTR2_TOTSM_SPGP(IS2GRP,NGAS,IS2SM,NS2OP,NS2STR,
     &       IS2_STR(IOFF),NOCOB,0,0,0)              
        LS2_STR(IS2SM) = NS2STR
        IOFF = IOFF + NS2STR*NS2OP
      END DO
*
      IBOFF = 1
      DO IS12SM = 1, NSMST
        DO IS1SM = 1, NSMST
          IS2SM = MULTD2H(IS1SM,IS12SM)
          IBS1S2(IS1SM,IS2SM) = IBOFF
          NS1STR = LS1_STR(IS1SM)
          NS2STR = LS2_STR(IS2SM)
          IBS1STR = IBS1_STR(IS1SM)
          IBS2STR = IBS2_STR(IS2SM)
          IF(NTEST.GE.100) 
     &    WRITE(6,*) ' NS1STR, NS2STR = ', NS1STR, NS2STR
          DO IS2 = 1, NS2STR
          DO IS1 = 1, NS1STR
* S1 * S2 => S12 ( all strings are considered to be creation strings
            IS1_OFF =  IBS1STR + (IS1-1)*NS1OP
            IS2_OFF =  IBS2STR + (IS2-1)*NS2OP
            CALL CRAN_STR(IS1_STR(IS1_OFF),IDUM,NS1OP,0,
     &                    IS2_STR(IS2_OFF),NS2OP,ISTR_OUT,
     &                    ISIGN,IZERO_STR)
            IJ = IBOFF -1 +(IS2-1)*NS1STR + IS1
            IF(IZERO_STR.NE.1) THEN
             INUM = ISTRNM(ISTR_OUT,NOCOB,NS12OP,IZ,IS12_REO,1)
             IS1S2_TO_S12(IJ) = INUM
             XS1S2_TO_S12(IJ) = DBLE(ISIGN)
            ELSE
             IS1S2_TO_S12(IJ) = 0
             XS1S2_TO_S12(IJ) = 0.0              
            END IF
          END DO
          END DO
*         ^ End of loop over IS1, IS2
          IBOFF = IBOFF + NS1STR*NS2STR
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' S1S2 => S12 mapping, adress and signs' 
        WRITE(6,*)  ' ====================================='
        DO IS1SM = 1, NSMST
        DO IS2SM = 1, NSMST
          NS1 = LS1_STR(IS1SM)
          NS2 = LS2_STR(IS2SM)
          IB = IBS1S2(IS1SM,IS2SM)
          WRITE(6,*) ' Sym of S1 and S2 : ', IS1SM, IS2SM
          WRITE(6,*)
          CALL IWRTMA(IS1S2_TO_S12(IB),NS1,NS2,NS1,NS2)
          CALL WRTMAT(XS1S2_TO_S12(IB),NS1,NS2,NS1,NS2)
        END DO
        END DO
      END IF
*
      RETURN
      END     
      SUBROUTINE STST_TO_ST_MAPO(IS1OC,IS2OC,IS12OC,
     &           IBS1S2,IS1S2_TO_S12,XS1S2_TO_S12,
     &           IZ,IZSCR,IS1_STR,IS2_STR,IS12_REO)
*
* Consider string multiplication S1*S2 => S12 
* Find the mapping between strings             
* Maps are organized as matrices (S1,S2) for each symmetry of S1 and S2
*
* Jeppe Olsen, May 1, 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input :  number of operators in each gasspace
C     INTEGER IS1OP(NGAS),IS2OP(NGAS),IS12OP(NGAS)
      INTEGER IS1OC(*),IS2OC(*),IS12OC(*)
*
*. Local scratch
*
      INTEGER IS1GRP(NGAS),IS2GRP(NGAS),IS12GRP(NGAS)
      INTEGER ISTR_OUT(100)
*
*. Scratch space through parameter list
*
      INTEGER IZ(*),IZSCR(*)
      INTEGER IS1_STR(*),IS2_STR(*),IS12_REO(*)
*     ^ Should hold largest list of strings with given CAAB
      INTEGER LS1_STR(8), LS2_STR(8)
*
*. Output
*
*. Offset to mappings for given sym of S1 and S2
      INTEGER IBS1S2(8,8)
*. And the mappings 
      INTEGER IS1S2_TO_S12(*)
      DIMENSION XS1S2_TO_S12(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' STST_TO_ST_MAP entered '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of S1, S2 and S12 '
        CALL IWRTMA(IS1OC,1,NGAS,1,NGAS)
        CALL IWRTMA(IS2OC,1,NGAS,1,NGAS)
        CALL IWRTMA(IS12OC,1,NGAS,1,NGAS)
      END IF
C?    WRITE(6,*) ' Memchk at start of ...MAP'
C?    CALL MEMCHK
C?    WRITE(6,*) ' Memcheck passed '
*
*. Occupation to group notation
C     CALL OCC_TO_GRP(ICA,ICAGP,1)
      CALL OCC_TO_GRP(IS1OC,IS1GRP,1)
      CALL OCC_TO_GRP(IS2OC,IS2GRP,1)
      CALL OCC_TO_GRP(IS12OC,IS12GRP,1)
*. Number of operators 
      NS1OP = IELSUM(IS1OC,NGAS)
      NS2OP = IELSUM(IS2OC,NGAS)
      NS12OP = IELSUM(IS12OC,NGAS)
*. Z array for 12 strings 
C          WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,ISCR,NTEST)        
      CALL WEIGHT_SPGP(IZ,NGAS,IS12OC,NOBPT,IZSCR,0)
*. Loop over symmetry of compound strings
      IBOFF = 1
      DO IS12SM = 1, NSMST
C?     WRITE(6,*) ' IS12SM(a) = ', IS12SM
*. Generate reorder array for strings of this symmetry
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
C?      WRITE(6,*) ' Generation of T12 strings '
        CALL GETSTR2_TOTSM_SPGP(IS12GRP,NGAS,IS12SM,NS12OP,NS12STR,
     &       IS1_STR,NOCOB,1,IZ,IS12_REO)
C?      WRITE(6,*) ' IREOS12 : '
C?      CALL IWRTMA(IREOS12,1,NS12STR,1,NS12STR)
C?     WRITE(6,*) ' IS12SM(b) = ', IS12SM
*. Generate reorder array for strings of this symmetry
        DO IS1SM = 1, NSMST
          IS2SM = MULTD2H(IS1SM,IS12SM)
C?        WRITE(6,*) ' IS12SM, IS1SM, IS2SM = ',
C?   &                 IS12SM, IS1SM, IS2SM
          IBS1S2(IS1SM,IS2SM) = IBOFF
*. Generate S1 and S2 strings of given symmetry
C?        WRITE(6,*) ' Generation of T1 strings '
          CALL GETSTR2_TOTSM_SPGP(IS1GRP,NGAS,IS1SM,NS1OP,NS1STR,
     &         IS1_STR,NOCOB,0,0,0)              
          LS1_STR(IS1SM) = NS1STR
C?        WRITE(6,*) ' Generation of T2 strings '
          CALL GETSTR2_TOTSM_SPGP(IS2GRP,NGAS,IS2SM,NS2OP,NS2STR,
     &         IS2_STR,NOCOB,0,0,0)              
          LS2_STR(IS2SM) = NS2STR
          IF(NTEST.GE.100) 
     &    WRITE(6,*) ' NS1STR, NS2STR = ', NS1STR, NS2STR
          DO IS2 = 1, NS2STR
          DO IS1 = 1, NS1STR
* S1 * S2 => S12 ( all strings are considered to be creation strings
C           CRAN_STR(ICR,IAN,NCR,NAN,ISTR_IN,NEL_IN,
C    &               ISTR_OUT,ISIGN,IZERO_STR)
            CALL CRAN_STR(IS1_STR((IS1-1)*NS1OP+1),IDUM,NS1OP,0,
     &                    IS2_STR((IS2-1)*NS2OP+1),NS2OP,ISTR_OUT,
     &                    ISIGN,IZERO_STR)
C?          WRITE(6,*) ' Output string from CRAN_STR '
C?          CALL IWRTMA(ISTR_OUT,1,NS12OP,1,NS12OP)
            IJ = IBOFF -1 +(IS2-1)*NS1STR + IS1
C?          WRITE(6,*) ' IJ = ', IJ
            IF(IZERO_STR.NE.1) THEN
C                   ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
             INUM = ISTRNM(ISTR_OUT,NOCOB,NS12OP,IZ,IS12_REO,1)
             IS1S2_TO_S12(IJ) = INUM
             XS1S2_TO_S12(IJ) = DBLE(ISIGN)
            ELSE
             IS1S2_TO_S12(IJ) = 0
             XS1S2_TO_S12(IJ) = 0.0              
            END IF
          END DO
          END DO
*         ^ End of loop over IS1, IS2
          IBOFF = IBOFF + NS1STR*NS2STR
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' S1S2 => S12 mapping, adress and signs' 
        WRITE(6,*)  ' ====================================='
        DO IS1SM = 1, NSMST
        DO IS2SM = 1, NSMST
          NS1 = LS1_STR(IS1SM)
          NS2 = LS2_STR(IS2SM)
          IB = IBS1S2(IS1SM,IS2SM)
          WRITE(6,*) ' Sym of S1 and S2 : ', IS1SM, IS2SM
          WRITE(6,*)
          CALL IWRTMA(IS1S2_TO_S12(IB),NS1,NS2,NS1,NS2)
          CALL WRTMAT(XS1S2_TO_S12(IB),NS1,NS2,NS1,NS2)
        END DO
        END DO
      END IF
*
      RETURN
      END     
      SUBROUTINE K_TO_J_TOT_SINGLE(IKJ,XKJ,KSM,IT,KSTR,
     &                            IM,XM,IBM,NK,LTOP)
*
* Obtain map !J> = T-Oper !K>
*
* Special version for a single T operator and a single K operator 
*
* !K> is of sym KSM and is string number KSTR of this symmetry
*
* Maps for each elementary operator is provided by IM with offset mat IBM
*
* Jeppe Olsen, May of 2000 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Input
      INTEGER IT(LTOP) 
      INTEGER IM(*),IBM(8,8,LTOP),NK(8,LTOP)
      DIMENSION XM(*)
*. Output
*
* IKJ Obtained strings 
* XKJ Sign
*
      NTEST = 00
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Information from K_TO_J.... '
        WRITE(6,*) ' =========================== '
        WRITE(6,*) ' LTOP  = ', LTOP
      END IF
*
      IF(LTOP.GT.0) THEN
C       DO ITOP = IT_B, IT_E
          KNSM = KSM
          SIGN = 1.0D0
          KNSTR = KSTR
          DO IOP = 1, LTOP
            JOB_ABS = IT(IOP)
            JOB_SM = ISMFTO(JOB_ABS)
            JOB_TP = ITPFTO(JOB_ABS)
            JOB_REL = JOB_ABS - IOBPTS(JOB_TP,JOB_SM) + 1
            LK = NK(KNSM,IOP)
            IADR = IBM(JOB_SM,KNSM,IOP) -1 + (JOB_REL-1)*LK+KNSTR
            KNSTR = IM(IADR)
            SIGN = SIGN*XM(IADR) 
            IF(KNSTR.EQ.0) GOTO 1001
            KNSM = MULTD2H(KNSM,JOB_SM)
          END DO
*
 1001     CONTINUE
          IKJ = KNSTR
          XKJ = SIGN
C       END DO
      ELSE
*. No K-operators => Identity map
        IKJ = KSTR
        XKJ = 1.0D0
      END IF
*
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Output from K_TO_J_TOT_SINGLE ' 
        WRITE(6,*) ' T-Operator '
        CALL IWRTMA(IT,LTOP,1,LTOP,1)
        WRITE(6,*) 'IKJ and XKJ = ', IKJ, XKJ
      END IF
*
      RETURN
      END
      SUBROUTINE Z_TCC_OFF2(IBT,LENGTH,NCA,NCB,NAA,NAB,ITSYM,NSMST)
*
* Offsets for symmetryblocks of TCC elements, sym of CA,CB,AA used 
*
* Jeppe Olsen, Summer of 99
*
* Compared to Z_TCC_OFF : Total length added to argument list
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NCA(*),NCB(*),NAA(*),NAB(*)
*. Output
      INTEGER IBT(8,8,8)
*
      IOFF = 1
      DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM_C,ITSYM) 
        DO ISM_CA = 1, NSMST
          ISM_CB = MULTD2H(ISM_C,ISM_CA)
          DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
            IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
            IOFF = IOFF + 
     &      NCA(ISM_CA)*NCB(ISM_CB)*NAA(ISM_AA)*NAB(ISM_AB)
          END DO
        END DO
      END DO
      LENGTH = IOFF - 1 
*
      RETURN
      END
      SUBROUTINE CONJ_CAAB(ICAAB_IN,ICAAB_OUT,NGAS,SIGN)
*
* Conjugate ICAAB_IN to obtain ICAAB_OUT
* 
*. Jeppe Olsen, Oct 2000
      INCLUDE 'implicit.inc'
*.Input 
      INTEGER ICAAB_IN(NGAS,4)
*. Output
      INTEGER ICAAB_OUT(NGAS,4)
*. AA_out is obtain by conjugating CA_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,1), ICAAB_OUT(1,3),NGAS)
      CALL ICOPVE(ICAAB_IN(1,1), ICAAB_OUT(1,3),NGAS)
*. AB_out is obtain by conjugating CB_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,2), ICAAB_OUT(1,4),NGAS)
      CALL ICOPVE(ICAAB_IN(1,2), ICAAB_OUT(1,4),NGAS)
*. CA_out is obtain by conjugating AA_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,3), ICAAB_OUT(1,1),NGAS)
      CALL ICOPVE(ICAAB_IN(1,3), ICAAB_OUT(1,1),NGAS)
*. CB_out is obtain by conjugating AB_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,4), ICAAB_OUT(1,2),NGAS)
      CALL ICOPVE(ICAAB_IN(1,4), ICAAB_OUT(1,2),NGAS)
*. Directly the conjugated operator is 
*  O(ab){\dag}O(aa){\dag}O(cb){\dag}O(ca){\dag}
*. Sign required to change  
*  O(aa){\dag}O(ab){\dag}O(ca){\dag}O(cb){\dag}
      NCA = IELSUM(ICAAB_IN(1,1),NGAS)
      NCB = IELSUM(ICAAB_IN(1,2),NGAS)
      NAA = IELSUM(ICAAB_IN(1,3),NGAS)
      NAB = IELSUM(ICAAB_IN(1,4),NGAS)
*
      NPERM = NCA*NCB + NAA*NAB
*. Sign required to bring the individual strings into ascending order
      NPERM = NPERM + 
     &        NCA*(NCA-1)/2+NCB*(NCB-1)/2+NAA*(NAA-1)/2+NAB*(NAB-1)/2
*
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGN = -1.0D0
      ELSE
        SIGN = 1.0D0
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CAAB operator and conjugated CAAB operator : '
        CALL WRT_SPOX_TP(ICAAB_IN,1)
        CALL WRT_SPOX_TP(ICAAB_OUT,1)
        WRITE(6,*) ' Sign = ', sign
      END IF
*
      RETURN
      END
      SUBROUTINE REV_ORD_IARR(IIN,IOUT,NELMNT)
*
* Reverse order of elements of input array IIN  to obtain IOUT
*
* Jeppe Olsen, October 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IIN(NELMNT)
*. Output 
      INTEGER IOUT(NELMNT)
*
      DO IELMNT = 1, NELMNT
        IOUT(NELMNT+1-IELMNT) = IIN(IELMNT)
      END DO
*
      RETURN
      END  
      FUNCTION IPERM_PARITY(IPERM,NELMNT)
*
* Find sign required to bring permutation IPERM 
* of the first NELMNT integers into order.
*
* KISS version
*
* Jeppe Olsen, Oct. 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IPERM(NELMNT)
*. Local scratch 
      PARAMETER(MAXELMNT=100)
      INTEGER ISCR(MAXELMNT)
*
      IF(NELMNT.GT.MAXELMNT) THEN
        WRITE(6,*) ' IPERM_PARITY in trouble, NELMNT > MAXELMNT'
        WRITE(6,*) ' NELMNT, MAXELMNT = ',  NELMNT, MAXELMNT
        STOP       ' IPERM_PARITY in trouble, NELMNT > MAXELMNT'
      END IF
*. Ensure that all elements from 1 to NELMNT are included
      I_AM_OKAY = 1
      DO I = 1, NELMNT
        IFOUND = 0
        DO IELMNT = 1, NELMNT
          IF(IPERM(IELMNT).EQ.I) IFOUND = 1
        END DO
        IF(IFOUND.EQ.0) I_AM_OKAY = 0
      END DO
*
      IF(I_AM_OKAY.EQ.0) THEN
        WRITE(6,*) ' Illegal input to IPERM_PARITY'
        CALL IWRTMA(IPERM,1,NELMNT,1,NELMNT)
        STOP       ' Illegal input to IPERM_PARITY'
      END IF
*
      CALL ICOPVE(IPERM,ISCR,NELMNT)
      LPERM = 0
      DO IELMNT = 1, NELMNT
*. Find IELMNT in ISCR (Cannot occur before IELMNT)
        KELMNT = 0
        DO JELMNT = IELMNT, NELMNT
          IF(ISCR(JELMNT).EQ.IELMNT) KELMNT = JELMNT
        END DO
        DO JELMNT = KELMNT, IELMNT+1, -1
          ISCR(JELMNT) = ISCR(JELMNT-1)
        END DO
        ISCR(IELMNT) = IELMNT
        LPERM = LPERM + KELMNT - IELMNT
C?      WRITE(6,*) ' Updated list for IELMNT =', IELMNT
C?      CALL IWRTMA(ISCR,1,NELMNT,1,NELMNT)
      END DO
      IF(MOD(LPERM,2).EQ.0) THEN
        ISIGN = 1
      ELSE
        ISIGN = -1
      END IF
*
      IPERM_PARITY= ISIGN
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Permutation : '
       CALL IWRTMA(IPERM,1,NELMNT,1,NELMNT)
       WRITE(6,*) ' Parity of permutation = ', ISIGN
      END IF
*
      RETURN
      END
      FUNCTION  IEXC2_LEVEL_FOR_CAAB(ICAAB)
*
* Find excitation level times two for elementary excitations
*
* Excitation level times 2 : Number of creation in sec 
*                          + number of annihilations of inac
*                          - Number of anni of sec 
*                          - Number of crea of inac
*                              
* Jeppe Olsen, June 2002
*
      INCLUDE 'implicit.inc'
*. General input 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER ICAAB(NGAS,4)
*
      IEXC_RANK = 0
*. Creation 
      DO IAB = 1, 2
        DO IGAS = 1, NGAS
          IF(IHPVGAS(IGAS).EQ.1) THEN
            IEXC_RANK = IEXC_RANK - ICAAB(IGAS,IAB)
          ELSE IF ( IHPVGAS(IGAS).EQ.2) THEN
            IEXC_RANK = IEXC_RANK + ICAAB(IGAS,IAB)
          END IF
        END DO
      END DO
*. Annihilation 
      DO IAB = 1, 2
        DO IGAS = 1, NGAS
          IF(IHPVGAS(IGAS).EQ.1) THEN
            IEXC_RANK = IEXC_RANK + ICAAB(IGAS,2+IAB)
          ELSE IF ( IHPVGAS(IGAS).EQ.2) THEN
            IEXC_RANK = IEXC_RANK - ICAAB(IGAS,2+IAB)
          END IF
        END DO
      END DO
*
      IEXC2_LEVEL_FOR_CAAB = IEXC_RANK
*
      NTEST = 0
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' Excitation rank ( 2 x exc level)', IEXC_RANK
      END IF
*
      RETURN
      END
      SUBROUTINE REO_BLK_VEC(VECIN,LBLKIN,NBLK,VECOUT,LBLKOUT,
     &                       INEW_TO_OLD,IWAY)
*
* Reorder a blocked vector using INEW_TO_OLD 
*
* output is reordered blocked vector, and number elements
* per block in reordered vector.
*
* As the offset to a given block pt is calculated on the flight
* there is an NBLK**2 procedure included, so it will be 
* timeconsuming when there is more than a few thousand blocks...
*
* Jeppe Olsen, August 2004
*
* IWAY = 1 : OLD => NEW, LBLKOUT is constructed
* IWAY =-1 : NEW => OLD, NO block info is constructed
*
* IWAY = -1 is written with the purpose of reverting reorder made using IWAY = 1
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VECIN(*),LBLKIN(NBLK)
*. Required reordering of blocks
      DIMENSION INEW_TO_OLD(NBLK)
*. Output
      DIMENSION VECOUT(*),LBLKOUT(*)
*.Loop over blocks in output order
      IOFFO = 1
      DO IBLKO = 1, NBLK
        IBLKI = INEW_TO_OLD(IBLKO)
        IOFFI = 1 + IELSUM(LBLKIN,IBLKI-1)
        LBLK = LBLKIN(IBLKI)
        IF(IWAY.EQ.1) LBLKOUT(IBLKO) = LBLK
        IF(IWAY.EQ.1) THEN
          CALL COPVEC(VECIN(IOFFI),VECOUT(IOFFO),LBLK)
        ELSE 
          CALL COPVEC(VECOUT(IOFFO),VECIN(IOFFI),LBLK)
        END IF
        IOFFO = IOFFO + LBLK
      END DO
*
      RETURN
      END 
      SUBROUTINE ICOPVE3(IIN,IOFFIN,IOUT,IOFFOUT,NDIM)
*
* IOUT(IOFFOUT-1+I) = IIN(IOFFIN-1+I),I = 1, NDIM
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION IIN(*)
*. Output
      DIMENSION IOUT(*)
*
      DO I = 1, NDIM
        IOUT(IOFFOUT-1+I) = IIN(IOFFIN-1+I)
      END DO
*
      RETURN
      END
      SUBROUTINE CHECK_HTF_APR(NSPOBEX,ISPOBEX,IOKAY)
*
* For the CC expansion given by ISPOBEX, check that
* the relevant parts of HEXP T can be written as
*
* H EXP T |HF > = H( sum_I E_I + sum_{IJ}T_I F_{IJ} ) | HF >
*
* The approach requires that no product T(I)T(J)T(K), 
* I.ge.J.ge.K(execution order) exists for  which
* 1 : T(J)T(K) is outside operator space
* 2 : T(I)T(J)T(K) can be connected to operator space by
*     excitation in H
*
* Jeppe Olsen, August 2001
*
*. WARNING (NOV. 2004) This does not catch all exceptions, 
*. even for rather simple expansions ....


c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER ISPOBEX(4*NGAS,*)
*. Local scratch
      PARAMETER(MXNTOP = 3)
      INTEGER ITTTOP(MXNTOP), ITTTOP2(MXNTOP)
      INTEGER IJKOCC(4*MXPNGAS),IIJKOCC(4*MXPNGAS)
*
      NTEST = 10  
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CH_HTF')  
*. local scratch for execution order
      CALL MEMMAN(KL_EXEORD,NSPOBEX,'ADDL  ',1,'EXE_OR')
*. Execution order
      I_DO_REO = 1
      IF(I_DO_REO.EQ.0) THEN
        CALL ISTVC2(WORK(KL_EXEORD),0,1,NSPOBEX)
C  ISTVC2(IVEC,IBASE,IFACT,NDIM)
      ELSE 
        CALL ORD_CCEXTP(WORK(KL_EXEORD))
      END IF
*. Loop over MXNTOP operators in nonrestrict ascending order
      IFIRST = 1
      NLOOP = 0
      NEXCP = 0
 1000 CONTINUE
      IF( IFIRST. EQ. 1 ) THEN
*. Initiate
       CALL ISETVC(ITTTOP2,1,MXNTOP)
       IFIRST = 0
       NONEW = 0
      ELSE
       CALL NXTORD_NS(ITTTOP2,MXNTOP,1,NSPOBEX,NONEW)
      END IF
* We obtained from NXTORD_NS a row of numbers (in execution order)
* in nonrestrict ascending  order. 
* But we wanted the T operators as T(I)T(K)T(K), I.ge.j ge.k -
* and we want the original t-operators
* but according to the introductory comments 
      DO IOP = 1, MXNTOP
        JOP = ITTTOP2(IOP)
C  IFRMR(WORK,IROFF,IELMNT)
        JJOP = IFRMR(WORK,KL_EXEORD,JOP)
        ITTTOP(MXNTOP+1-IOP) = JJOP
      END DO
*. Check if T(J)T(K) does not correspond to a CC type and 
*.          T(I)T(J)T(K) can be connected to type by H
*. Occupation of T(J)T(K)
C  OP_T_OCC(IOPOCC,ITOCC,IOPTOCC,IMZERO)
      CALL OP_T_OCC(ISPOBEX(1,ITTTOP(2)),ISPOBEX(1,ITTTOP(3)),
     &              IJKOCC,IJKZERO)
*. Well, single excitations are handled differently so check of any of the 
*. operators are single excitations
      NOP1 = IELSUM(ISPOBEX(1,ITTTOP(1)),4*NGAS)
      NOP2 = IELSUM(ISPOBEX(1,ITTTOP(2)),4*NGAS)
      NOP3 = IELSUM(ISPOBEX(1,ITTTOP(3)),4*NGAS)
      IF(NOP1.EQ.2.OR.NOP2.EQ.2.OR.NOP3.EQ.2) IJKZERO = 1
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Triplet of operators ', (ITTTOP(I),I=1,3)
        WRITE(6,*) ' IJKZERO = ', IJKZERO
      END IF
      IF(IJKZERO.EQ.0) THEN
*. Is T(J)T(K) operator in current list
        CALL INUM_FOR_OCC2(IJKOCC,IJKNUM,NDIFF)
        IF(NTEST.GE.1000) WRITE(6,*) ' IJKNUM = ', IJKNUM
*. Continue check only if TJK is out of space
        IF(IJKNUM.EQ.-1) THEN                      
*. Find occupation of T(I)T(J)T(K)
         CALL OP_T_OCC(ISPOBEX(1,ITTTOP(1)),IJKOCC,
     &               IIJKOCC,IIJKZERO)
        IF(NTEST.GE.1000) WRITE(6,*) ' IIJKZERO = ', IIJKZERO
*. Number of operators required to connect I321OCC with reference
         IF(IIJKZERO.EQ.0) THEN
          CALL INUM_FOR_OCC2(IIJKOCC,IIJKNUM,NIJKDIFF)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IIJKNUM,NIJKDIFF = ', IIJKNUM,NIJKDIFF
          END IF
          IF(NIJKDIFF.LE.4) THEN
*. An exception has been found
           NEXCP = NEXCP + 1
           IF(NTEST.GE.10) THEN
             WRITE(6,*) ' T(I) T(J) T(K) not included in H_TF approach'
             CALL WRT_SPOX_TP(ISPOBEX(1,ITTTOP(1)),1)
             CALL WRT_SPOX_TP(ISPOBEX(1,ITTTOP(2)),1)
             CALL WRT_SPOX_TP(ISPOBEX(1,ITTTOP(3)),1)
             WRITE(6,*) 
     &       ' Number of ops required for connection = ', NIJKDIFF
           END IF
          END IF
         END IF
        END IF
*       ^ End if IJKNUM = 0  
      END IF
*     ^ End if IJKZERO = 0
      NLOOP = NLOOP + 1
      IF(NONEW.EQ.0) GOTO 1000
*
C?    IF(NLOOP.EQ.10000) THEN
C?      WRITE(6,*) ' Check stopped due to infinite loop (?) '
C?      STOP ' Check stopped due to infinite loop (?) '
C?    END IF
*
      IF(NEXCP.EQ.0) THEN
*. No exceptions occured, so expansion is OKAY
       IOKAY = 1
      ELSE
       WRITE(6,*) ' Number of exceptions in H_TF approach = ', NEXCP
       STOP ' H_TF approach is not valid '
      END IF
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CH_HTF')  
*
      RETURN
      END
      SUBROUTINE NXTORD_NS(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An ordered set of numbers INUM(I),I=1,NELMNT is
* given in nonstrict ascending order
* (adjecent elements may be identical)
* Values of INUM(*) is
* restricted to the interval MINVAL,MAXVAL .
*
* Find next higher number in nonstrict ascending order.
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
* If the array is zero-dimensional (NELMNT = 0),
* the routine sends NONEW = 1 back
*
* Jeppe Olsen May 1989
*
      DIMENSION INUM(*)
*
       NTEST = 00
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTORD '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
        NONEW = 1
        GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF( IPLACE .LT. NELMNT .AND.
     &      INUM(IPLACE)+1 .LE. INUM(IPLACE+1)
     &  .OR.IPLACE.EQ. NELMNT .AND.
     &      INUM(IPLACE)+1.LE.MAXVAL) THEN
              INUM(IPLACE) = INUM(IPLACE) + 1
              NONEW = 0
              GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
              IF(IPLACE .EQ. 1 ) THEN
                INUM(IPLACE) = MINVAL
              ELSE
                INUM(IPLACE) = INUM(IPLACE-1) 
              END IF
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
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
*


c $Id$
