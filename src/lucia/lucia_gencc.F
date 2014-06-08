      SUBROUTINE PRINT_CAAB_LIST(ISM)
*
*. Print complete  list of CAAB operators in standard
*. order- zero-partice operator included
*
* Jeppe Olsen, Aug. 2009, ( modified ANA_GENCC)
*
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*
* =====
*.Input
* =====
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cc_exc.inc'
*
      CALL LUCIAQENTER('PRCAA')
      CALL MEMMAN(KLOFF,DUMMY,'MARK  ',DUMMY,'PRCAAB')
*
*. Four blocks of string occupations
      CALL MEMMAN(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
      IUSLAB = 0
*
*. Occupation of strings of given sym and supergroup
      CALL PRINT_CAAB_LIST_INNER(WORK(KLSOBEX),WORK(KLOBEX),
     &            NSPOBEX_TPE,NOBEX_TP,ISM,
     &            WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &            WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &            IUSLAB,
     &            IDUMMY,
     &            NTOOB,IPRNCIV,WORK(KLSOX_TO_OX),MSCOMB_CC)
   
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'PRCAAB')
      CALL LUCIAQEXIT('PRCAA')
*
      RETURN
      END
      SUBROUTINE PRINT_CAAB_LIST_INNER(ISPOBEX_TP,IOBEX_TP,
     &                 NSPOBEX_TP,NOBEX_TP,ISM,
     &                 IOCC_CA, IOCC_CB, IOCC_AA, IOCC_AB,
     &                 IUSLAB,
     &                 IOBLAB,
     &                 NORB,IPRNCIV,ISOX_TO_OX,MSCOMB_CC)
*
* Print list of CAAB operators in standard order
*
*
*
* Jeppe Olsen , Aug. 2009
*                                                  

*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
      DIMENSION ISOX_TO_OX(NSPOBEX_TP)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      CHARACTER*6 IOBLAB(*)
*
*
      WRITE(6,*)
      WRITE(6,*) ' Operators are written as : '
      WRITE(6,*)
      WRITE(6,*)   ' Creation of alpha '
      WRITE(6,*)   ' Creation of beta '
      WRITE(6,*)   ' Annihilation of alpha '
      WRITE(6,*)   ' Annihilation of beta '
      WRITE(6,*)
      NTVAR = 0
*
      IT = 0
      IIT = 0
      DO ITSS = 1, NSPOBEX_TP
*. Transform from occupations to groups
       CALL OCC_TO_GRP(ISPOBEX_TP(1+0*NGAS,ITSS),IGRP_CA,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+1*NGAS,ITSS),IGRP_CB,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+2*NGAS,ITSS),IGRP_AA,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+3*NGAS,ITSS),IGRP_AB,1      )
*
       NEL_CA = IELSUM(ISPOBEX_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ISPOBEX_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ISPOBEX_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ISPOBEX_TP(1+3*NGAS,ITSS),NGAS)
*. Diagonal restricted type of spinorbital excitation ?
       IF(MSCOMB_CC.EQ.1) THEN
         CALL DIAG_EXC_CC(
     &        ISPOBEX_TP(1+0*NGAS,ITSS),ISPOBEX_TP(1+1*NGAS,ITSS),
     &        ISPOBEX_TP(1+2*NGAS,ITSS),ISPOBEX_TP(1+3*NGAS,ITSS),
     &        NGAS,IDIAG)
       ELSE
         IDIAG = 0
       END IF
*
       IF(IPRNCIV.EQ.1) THEN
          WRITE(6,*) ' Type of spin-orbital-excitation '
       END IF
*
       DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        DO ISM_CA = 1, NSMST
         ISM_CB = MULTD2H(ISM_C,ISM_CA)
         DO ISM_AA = 1, NSMST
          ISM_AB =  MULTD2H(ISM_A,ISM_AA)
          ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
          ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
*
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.LT.ISM_BETA) GOTO 777
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN 
           IRESTRICT_LOOP = 1
          ELSE
           IRESTRICT_LOOP = 0
          END IF
*. obtain strings
          CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,NSTR_CA,
     &         IOCC_CA, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,NSTR_CB,
     &         IOCC_CB, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,NSTR_AA,
     &         IOCC_AA, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,NSTR_AB,
     &         IOCC_AB, NORB,0,IDUM,IDUM)
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
          DO I_AB = 1, NSTR_AB
           IF(IRESTRICT_LOOP.EQ.1) THEN
             I_AA_MIN = I_AB
           ELSE
             I_AA_MIN = 1
           END IF
           DO I_AA = I_AA_MIN, NSTR_AA
            DO I_CB = 1, NSTR_CB
             IF(IRESTRICT_LOOP.EQ.1.AND.I_AA.EQ.I_AB) THEN
               I_CA_MIN = I_CB
             ELSE
               I_CA_MIN = 1
             END IF
             DO I_CA = I_CA_MIN, NSTR_CA
              IT = IT + 1
C?            WRITE(6,*) ' IT, T(IT) = ', IT,T(IT)
*
                WRITE(6,'(A)')
                WRITE(6,'(A)')
     &          '                 =================== '
                WRITE(6,*)

                WRITE(6,'(A,2I8,2X,E14.8)')
     &          '  Type, number, CAAB: ',  ITSS, IT
                IF(IUSLAB.EQ.0) THEN
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CA(IEL+(I_CA-1)*NEL_CA),IEL = 1, NEL_CA)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CB(IEL+(I_CB-1)*NEL_CB),IEL = 1, NEL_CB)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AA(IEL+(I_AA-1)*NEL_AA),IEL = 1, NEL_AA)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AB(IEL+(I_AB-1)*NEL_AB),IEL = 1, NEL_AB)
                END IF
             END DO
*            ^ End of loop over alpha creation strings
            END DO
*           ^ End of loop over beta creation strings
           END DO
*          ^ End of loop over alpha annihilation 
          END DO 
*         ^ End of loop over beta annihilation 
  777    CONTINUE
         END DO
        END DO
       END DO
*      ^ End of loop over symmetry blocks
      END DO
*     ^ End of loop over over types of excitations
      RETURN
      END
      SUBROUTINE SELECT_AB_TYPES(NSPOBEX_TP,ISPOBEX_TP,
     &           ISPOBEX_PAIRS,NGAS)
*
* Select which of the paired orbital excitation types that 
* should be active 
*
* Jeppe Olsen, July12, 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)  
*. Input and output
      INTEGER ISPOBEX_PAIRS(*)
*. On output an passive type is flagged by a negative number 
*. in ISPOBEX_PAIRS
*
      DO ITP = 1, NSPOBEX_TP
        JTP = ISPOBEX_PAIRS(ITP)
c quick fix
        IF (JTP.EQ.0) CYCLE
        IF(ITP.GT.JTP) THEN
*. Select the type with most alphaoperators as the active
          NCREA_I = IELSUM(ISPOBEX_TP(1,ITP),NGAS)
          NCREA_J = IELSUM(ISPOBEX_TP(1,JTP),NGAS)
          IF(NCREA_I.GE.NCREA_J) THEN
            ISPOBEX_PAIRS(JTP) = - ISPOBEX_PAIRS(JTP)
          ELSE 
            ISPOBEX_PAIRS(ITP) = - ISPOBEX_PAIRS(ITP)
          END IF
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Excitation types divided into active and passive '
        CALL IWRTMA(ISPOBEX_PAIRS,NSPOBEX_TP,1,NSPOBEX_TP,1)
      END IF
*
      RETURN
      END
         
      SUBROUTINE ABFLIP_SPOXTP(ICAAB_IN,ICAAB_OUT,NGAS)
*
* Obtain spin-orbital excitation type ICAAB_OUT by 
* spinflipping ICAAB_IN
*
* Jeppe Olsen, July 11, 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICAAB_IN(NGAS,4)
*. Output
      INTEGER ICAAB_OUT(NGAS,4)
*
      CALL ICOPVE(ICAAB_IN(1,1),ICAAB_OUT(1,2),NGAS)
      CALL ICOPVE(ICAAB_IN(1,2),ICAAB_OUT(1,1),NGAS)
      CALL ICOPVE(ICAAB_IN(1,3),ICAAB_OUT(1,4),NGAS)
      CALL ICOPVE(ICAAB_IN(1,4),ICAAB_OUT(1,3),NGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Spin-orbital excitation type '
        CALL WRT_SPOX_TP(ICAAB_IN,1)
        WRITE(6,*) ' Spinflipped orbital excitation type '
        CALL WRT_SPOX_TP(ICAAB_OUT,1)
      END IF
*
      RETURN
      END
      SUBROUTINE SPOBEXTP_PAIRS(
     &           NSPOBEX_TP,ISPOBEX_TP,NGAS,ISPOBEX_PAIRS)
*
* Find pairs of spin-orbital excitaion types related to each other 
* by spin-flip
*
* Jeppe Olsen, July11 2001
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      DIMENSION ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
*. Output
      INTEGER ISPOBEX_PAIRS(NSPOBEX_TP)
*. Local scratch 
      INTEGER IABFLIP(MXPNGAS*4)
*
      IZERO = 0
      CALL ISETVC(ISPOBEX_PAIRS,IZERO,NSPOBEX_TP)
*
      DO ITP = 1, NSPOBEX_TP
      IF(ISPOBEX_PAIRS(ITP).EQ.0) THEN
*. Perform spinflip on type ITP
        CALL ABFLIP_SPOXTP(ISPOBEX_TP(1,ITP),IABFLIP,NGAS)
*. Find address of spinflipped type in input list
        DO JTP = 1, NSPOBEX_TP
        IF(ISPOBEX_PAIRS(JTP).EQ.0) THEN
          IDIFF = 0
          DO IOP = 1, 4*NGAS
            IF(IABFLIP(IOP).NE.ISPOBEX_TP(IOP,JTP)) IDIFF = IDIFF + 1
          END DO
          IF(IDIFF.EQ.0) THEN
           ISPOBEX_PAIRS(ITP) = JTP
           ISPOBEX_PAIRS(JTP) = ITP
          END IF
        END IF
        END DO
*       ^ End of loop over JTP
      END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Table of paired spin-orbital excitations '
        CALL IWRTMA(ISPOBEX_PAIRS,1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*
      RETURN
      END
      SUBROUTINE ISTVC3(IVEC,IOFF,IVAL,NDIM)
*
* IVEC(IOFF-1+I) = IVAL, I = 1, NDIM
*
* Jeppe Olsen, Oct 2000
*
      INCLUDE 'implicit.inc'
*. Output
      INTEGER IVEC(*)
      DO IELMNT = IOFF, IOFF-1+NDIM
        IVEC(IELMNT) = IVAL
      END DO
*
      RETURN
      END
      SUBROUTINE FRZ_SPOBEX(ISPOBEX_FRZ,IEXC_FOR_OBEX,NSPOBEX_TP,
     &                      ISOX_TO_OX,IFRZ_CC_AR,NFRZ_CC)
*
* Obtain the spinorbital excitations that should be frozen 
*
* Jeppe Olsen, September 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IFRZ_CC_AR(NFRZ_CC), IEXC_FOR_OBEX(NSPOBEX_TP)
      INTEGER ISOX_TO_OX(NSPOBEX_TP)
*. Output 
      INTEGER ISPOBEX_FRZ(NSPOBEX_TP)
*
      DO ISPOBEX = 1, NSPOBEX_TP
        IEXC =  IEXC_FOR_OBEX(ISOX_TO_OX(ISPOBEX))
*. Is this level frozen 
        IFREEZE = 0
        DO IFRZ_CC = 1, NFRZ_CC
         IF( IFRZ_CC_AR(IFRZ_CC).EQ.IEXC) IFREEZE = 1
        END DO
C?      WRITE(6,*) ' ISPOBEX, IEXC, IFREEZE =',
C?   &               ISPOBEX, IEXC, IFREEZE
        ISPOBEX_FRZ(ISPOBEX) = IFREEZE 
      END DO
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Excitation levels to be frozen ' 
       CALL IWRTMA(ISPOBEX_FRZ,1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*
      RETURN
      END
      SUBROUTINE READ_CITT_BLK(C,IBLOCK,NSMST,ISM,IDC,IDIAG,LUC)
*
* Read in TT block of determinant/combination coefficients
* (all symmetries read in )
*
* Elements are stored in C according to IBLOCK
* Jeppe Olsen, July 2000, HNIE
*
      INCLUDE 'implicit.inc'
*. Input 
      INTEGER IBLOCK(8,*)
*. Output 
      DIMENSION C(*)
*. Number of blocks
      NBLOCK = NSMST
      IF(IDC.EQ.2.AND.IDIAG.EQ.1.AND.ISM.NE.1) NBLOCK = NSMST/2
*
      DO JBLOCK = 1, NBLOCK
        JOFF = IBLOCK(6,JBLOCK)
        LEN  = IBLOCK(8,JBLOCK)
        CALL IFRMDS(LEN2,1,-1,LUC)
        CALL FRMDSC(C(JOFF),LEN,-1,LUC,
     &       I_AM_ZERO,I_AM_PACKED)
      END DO
*
      RETURN
      END
      SUBROUTINE CITT_BLK_REFRM(CPCK,CEXP,IDC,ISM,NASTR,NBSTR,
     &                          IDIAG,NSMST,IWAY,ICOPY)
*
* Reform CI type type block containing all symmetry blocks
*
* Jeppe Olsen, July 2000, HNIE
*
*. IWAY = 1 : Det to comb
*       = 2   Comb to det
*  ICOPY = 1 : Copy result back to incoming matrix
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*.Input
      INTEGER NASTR(NSMST), NBSTR(NSMST)
*. Output and input
      DIMENSION CPCK(*),CEXP(*)
*
      IF(IDC.EQ.1.OR.IDIAG.EQ.0.AND.ICOPY.EQ.0) THEN
*. Just copy 
        LEN = 0
        DO IASM = 1, NSMST
          IBSM = MULTD2H(ISM,IASM)
          LEN = LEN + NASTR(IASM)*NBSTR(IBSM)
        END DO
        IF(IWAY.EQ.1) THEN
          CALL COPVEC(CEXP,CPCK,LEN)
        ELSE
          CALL COPVEC(CPCK,CEXP,LEN)
        END IF
      ELSE 
*. IDC = 2 and IDIAG = 1 : Blocks are packed
        IOFF_PCK = 1
        IOFF_EXP = 1
        IF(ISM.EQ.1) THEN
*. Each block is lower triangular, expand to complete block
          DO IASM = 1, NSMST
            NSTR = NASTR(IASM)
            ONE = 1.0D0
            CALL TRIPK3(CEXP,CPCK,IWAY,NSTR,NSTR,ONE) 
            IOFF_PCK = IOFF_PCK + NSTR*(NSTR+1)/2
            IOFF_EXP = IOFF_EXP + NSTR*NSTR
          END DO
        ELSE 
*. Blocks with IASM.GT.IBSM are given in complete form in CPCK
          IOFF_EXP = 1
          IOFF_PCK = 1
          DO IASM = 1, NSMST
            IBSM = MULTD2H(IASM,ISM)
            NA = NASTR(IASM)
            NB = NBSTR(IBSM)
            IF(IASM.GT.IBSM) THEN
*. C(IASM,IBSM) is directly given in CPCK, just copy 
              LEN = NA*NB
              IF(IWAY.EQ.1) THEN
                CALL COPVEC(CEXP(IOFF_EXP),CPCK(IOFF_PCK),LEN)
                ELSE
                CALL COPVEC(CPCK(IOFF_PCK),CEXP(IOFF_EXP),LEN)
              END IF
              IOFF_PCK = IOFF_PCK + NA*NB
              IOFF_EXP = IOFF_EXP + NA*NB
            ELSE IF(IBSM.GT.IASM.AND.IWAY.EQ.2) THEN
*. Obtain C(IA,IB) by transposing C(IB,IA)
*. Offset to C(IBSM,IASM) in packed block
              JOFF_PCK = 1
              DO JASM = 1, IBSM-1
                JBSM = MULTD2H(ISM,JASM)
                IF(JASM.GT.JBSM)   
     &          JOFF_PCK = JOFF_PCK + NASTR(JASM)*NBSTR(JBSM)
              END DO
              CALL TRPMT3(CPCK(JOFF_PCK),NB,NA,CEXP(IOFF_EXP))
              IOFF_EXP = IOFF_EXP + NA*NB
            END IF
          END DO
*         ^ End of loop over symmetries
        END IF
*       ^ End if ISM = 1
        IF(ICOPY.EQ.1) THEN
*. Copy back 
          IF(IWAY.EQ.1) THEN
            LEN = IOFF_PCK-1
            CALL COPVEC(CPCK,CEXP,LEN)
          ELSE 
            LEN = IOFF_EXP - 1
            CALL COPVEC(CEXP,CPCK,LEN)
          END IF
        END IF
*       ^ End if copying is required
      END IF
*     ^ End if packed
*
      RETURN
      END
      SUBROUTINE TRP_CITT_BLK(CIN,COUT,ISM,NASTR,NBSTR,NSMST,ICOPY)
*
* Transpose all symmetry blocks of a TT block of a CI expansion
* C(IA,IB) => C(IB,IA)
*
* All blocks are assumed to be present in complete form in CIN
*
*. Order of blocks is dictated by symmetry of first index
*
* Jeppe Olsen, HNIE, July 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'multd2h.inc'
*. Input
      DIMENSION CIN(*) 
      INTEGER NASTR(NSMST),NBSTR(NSMST)
*. Output
       DIMENSION COUT(*)
*. Local scratch
      INTEGER IOUT_OFF(MXPNSMST)
*. Set up offsets for output block
      IOFF = 1
      DO IBSM = 1, NSMST
        IASM = MULTD2H(ISM,IBSM)
        IOUT_OFF(IBSM) = IOFF
        IOFF = IOFF + NBSTR(IBSM)*NASTR(IASM)
      END DO
*
      IOFF_IN = 1
      DO IASM = 1, NSMST
        IBSM = MULTD2H(IASM,ISM)
        NA = NASTR(IASM)
        NB = NBSTR(IBSM)
        CALL TRPMT3(CIN(IOFF_IN),NA,NB,COUT(IOUT_OFF(IBSM)))
        IOFF_IN = IOFF_IN + NA*NB
      END DO
*
      IF(ICOPY.EQ.1) THEN
        LEN = IOFF_IN - 1 
        CALL COPVEC(COUT,CIN,LEN)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_TVEC(ITSS_TP,LTSS_TP,NTSS_TP,T,ISM)
*
* Write T-vector 
*
* Jeppe Olsen, July 2000, HNIE 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER ITSS_TP(4*NGAS,NTSS_TP), LTSS_TP(*)
*. Input and Output
      DIMENSION T(*)
*. Local scratch
      INTEGER NCA(8),NCB(8),NAA(8),NAB(8)
*
      IOFF = 1
      DO ITSS = 1, NTSS_TP
*. Number of strings per symmetry
       CALL NST_SPGP(ITSS_TP(1+0*NGAS,ITSS),NCA)
       CALL NST_SPGP(ITSS_TP(1+1*NGAS,ITSS),NCB)
       CALL NST_SPGP(ITSS_TP(1+2*NGAS,ITSS),NAA)
       CALL NST_SPGP(ITSS_TP(1+3*NGAS,ITSS),NAB)
*. Diagonal block ?
       IF(MSCOMB_CC.EQ.1) THEN
         CALL DIAG_EXC_CC(ITSS_TP(1+0*NGAS,ITSS),
     &                    ITSS_TP(1+1*NGAS,ITSS),
     &                    ITSS_TP(1+2*NGAS,ITSS),
     &                    ITSS_TP(1+3*NGAS,ITSS),NGAS,IDIAG)
       ELSE 
        IDIAG = 0
       END IF
C           WRT_TCC_BLK2(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST,IDIAG)
       CALL WRT_TCC_BLK2(T(IOFF),ISM,NCA,NCB,NAA,NAB,NSMST,IDIAG)
       IOFF = IOFF + LTSS_TP(ITSS)
      END DO
*
      RETURN
      END
      SUBROUTINE WRT_TCC_BLK2(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST,IDIAG)
*
*. Write TCC block containing all symmetries, total sym of T is ITCC_SM
*. If IDIAG = 1, only the lower half of the block is assumes stored
*
* Jeppe Olsen, summer of 99
*              July 2000 : Idiag added 
*
       INCLUDE 'implicit.inc'
       INCLUDE 'multd2h.inc'      
*. Input
       DIMENSION TCC(*)
       INTEGER NCA(*), NCB(*), NAA(*), NAB(*)
*
       WRITE(6,*) ' Block of Coupled cluster vector '
       WRITE(6,*) ' ================================'
       IOFF = 1
       DO ISM_C = 1, NSMST
         ISM_A = MULTD2H(ITCC_SM,ISM_C) 
         DO ISM_CA = 1, NSMST
           ISM_CB = MULTD2H(ISM_C,ISM_CA)
           DO ISM_AA = 1, NSMST
             ISM_AB =  MULTD2H(ISM_A,ISM_AA)
C            WRITE(6,*) ' ISM_AB = ', ISM_AB
*
             ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
             ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
             LENGTH = 0
             LCA = NCA(ISM_CA)
             LCB = NCB(ISM_CB)
             LAA = NAA(ISM_AA)
             LAB = NAB(ISM_AB)
*
             IF(IDIAG.EQ.0.OR.ISM_ALPHA.GT.ISM_BETA) THEN
               LENGTH = LCA*LCB*LAA*LAB
             ELSE IF(IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN
               LEN_ALPHA = LCA*LAA
               LENGTH = LEN_ALPHA*(LEN_ALPHA+1)/2
             END IF
*
             IF(LENGTH.NE.0) THEN
               WRITE(6,'(A,4I4)') ' Sym of CA, CB, AA, AB ',
     &         ISM_CA, ISM_CB, ISM_AA, ISM_AB  
               CALL WRTMAT(TCC(IOFF),1,LENGTH,1,LENGTH)
             END IF
             IOFF = IOFF + LENGTH
           END DO
         END DO
      END DO
*
      RETURN
      END
      SUBROUTINE RENORM_T(ITSS_TP,LTSS_TP,NTSS_TP,T,ISM,IWAY)
*
* Renormalize T-coefficients
*
* IWAY = 1 : Spinorbital => Combination (multiply by sqrt(2))
* IWAY = 2 : Combination => Spinorbital (divide by sqrt(2))
*
* Jeppe Olsen, July 2000, HNIE 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER ITSS_TP(4*NGAS,NTSS_TP), LTSS_TP(NTSS_TP)
*. Input for output
      DIMENSION T(*)
*. Local scratch
      INTEGER NCA(8),NCB(8),NAA(8),NAB(8)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' In RENORM_T '
        WRITE(6,*) ' IWAY = ', IWAY
      END IF
*
      IF(IWAY.EQ.1) THEN
        FACTOR = SQRT(2.0D0)
        FACTORI = 1.0D0/FACTOR
      ELSE
        FACTOR = 1.0D0/SQRT(2.0D0)
        FACTORI = 1.0D0/FACTOR
      END IF
      IOFF = 1
      DO ITSS = 1, NTSS_TP
*. Number of strings per symmetry
       CALL NST_SPGP(ITSS_TP(1+0*NGAS,ITSS),NCA)
       CALL NST_SPGP(ITSS_TP(1+1*NGAS,ITSS),NCB)
       CALL NST_SPGP(ITSS_TP(1+2*NGAS,ITSS),NAA)
       CALL NST_SPGP(ITSS_TP(1+3*NGAS,ITSS),NAB)
*. Diagonal block ?
       CALL DIAG_EXC_CC(ITSS_TP(1+0*NGAS,ITSS),
     &                  ITSS_TP(1+1*NGAS,ITSS),
     &                  ITSS_TP(1+2*NGAS,ITSS),
     &                  ITSS_TP(1+3*NGAS,ITSS),NGAS,IDIAG)
       IF(IDIAG.EQ.0) THEN
         IRESTRICT = 0
       ELSE 
         IRESTRICT = 1
       END IF
       DO ISM_C = 1, NSMST
         ISM_A = MULTD2H(ISM,ISM_C) 
         DO ISM_CA = 1, NSMST
           ISM_CB = MULTD2H(ISM_C,ISM_CA)
           DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*
            ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
            ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
            IF(IRESTRICT.EQ.1.AND.ISM_BETA.GT.ISM_ALPHA) GOTO 777
            IF(IRESTRICT.EQ.0.OR.ISM_ALPHA.GT.ISM_BETA) THEN
             IRESTRICT_LOOP = 0
            ELSE
             IRESTRICT_LOOP = 1
            END IF
*
            LCA = NCA(ISM_CA)
            LCB = NCB(ISM_CB)
            LAA = NAA(ISM_AA)
            LAB = NAB(ISM_AB)
*
            IF(IRESTRICT_LOOP.EQ.0) THEN
*. Nondiagonal block
              LEN = LCA*LCB*LAA*LAB
              CALL SCALVE(T(IOFF),FACTOR,LEN)
              IOFF = IOFF + LEN
            ELSE
*. Diagonal block
              LEN_ALPHA = LCA*LAA
              LEN = LEN_ALPHA*(LEN_ALPHA+1)/2
              CALL SCALVE(T(IOFF),FACTOR,LEN)
*. Diagonal elements were also multiplied with FACTOR, 
*. countermultiply
              DO ICA = 1, LCA
                DO IAA = 1, LAA
C                        ITDIANUM(ICA,ICB,IAA,IAB,NC,NA)
                  IADR = ITDIANUM(ICA,ICA,IAA,IAA,LCA,LAA)
                  T(IOFF-1+IADR) = FACTORI*T(IOFF-1+IADR)
                END DO
              END DO
              IOFF = IOFF + LEN
            END IF
  777       CONTINUE
           END DO
         END DO
       END DO
*      ^ End of loop over symmetries of given TT block 
      END DO
*     ^ End of loop over TT blocks
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Renormalized T-vector '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
C     WRT_TVEC(ITSS_TP,LTSS_TP,NTSS_TP,T,ISM)
        CALL WRT_TVEC(ITSS_TP,LTSS_TP,NTSS_TP,T,ISM)
      END IF
*
      RETURN
      END
      FUNCTION NST_FOR_OCC(NELEC_PER_GAS,NORB_PER_GAS,NGAS)
*
* A supergroup is defined by NELEC_PER_GAS. Find the 
* Total number of strings of this supergroup
*
* Jeppe Olsen, July 2000 (HNIE)
*
      INCLUDE 'implicit.inc'
*
      INTEGER NELEC_PER_GAS(*),NORB_PER_GAS(*)
*
      NSTR = 1
      DO IGAS = 1, NGAS
       IF(NELEC_PER_GAS(IGAS).NE.0) THEN
        NSTR = NSTR*IBION(NORB_PER_GAS(IGAS),NELEC_PER_GAS(IGAS))
       END IF
      END DO
*
      NST_FOR_OCC = NSTR
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of electrons and orbitals per gasspace '
        CALL IWRTMA(NELEC_PER_GAS,1,NGAS,1,NGAS)
        CALL IWRTMA(NORB_PER_GAS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of strings = ', NSTR
      END IF
*
      RETURN
      END
      SUBROUTINE LEN_GENOP_STR_MAP(
     &           NGENOP,IGENOP,NSPGRP,ISPGRP,NOBPT,NGAS,
     &           MAXLEN)
*
* A set of general operators, IGENOP, and a set of supergroups, ISPGRP, 
* are defined. Obtain max.length of mappings from resolution strings 
* to istrings. 
*
* The resolution strings KSTR are assumed to be inserted between the 
* creation and annihilation operators 
*
*   <ISPGRP!ICREA!KSTR> <KSTR!IANNI!ISPGRP'>
* So the mappings from KSTR to ISTR are always creation mappings 
*
* Jeppe Olsen, July 2000 ( At summerschool, HNIE)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Specific input
      INTEGER IGENOP(2*NGAS,*), ISPGRP(MXPNGAS,*), NOBPT(NGAS)
*. Local scratch
      INTEGER KSPGRP(MXPNGAS), IGENOP_EXP(MXPLCCOP)
*
      NTEST = 00
      MAXLEN = 0
      IONE = 1
      IMONE = -1
      DO ICA = 1, 2
      DO JSPGRP = 1, NSPGRP
*. ICREA/IANNI(dag)*!KSTR> = !ISTR>
       DO JGENOP = 1, NGENOP
*. Occupation of KSPGRP
         IF(ICA.EQ.1) THEN
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(1,JGENOP),
     &                 IONE,IMONE,NGAS)
         ELSE
           CALL IVCSUM(KSPGRP,ISPGRP(1,JSPGRP),IGENOP(NGAS+1,JGENOP),
     &                 IONE,IMONE,NGAS)
         END IF
*. Is KSPGRP a correct supergroup( all occ larger than zero) 
         I_AM_OKAY = 1
         DO IGAS = 1, NGAS
           IF(KSPGRP(IGAS).LT.0) I_AM_OKAY = 0
           IF(KSPGRP(IGAS).GT.NOBPT(IGAS)) I_AM_OKAY = 0
         END DO
         IF(I_AM_OKAY.EQ.1) THEN
*
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Initial occ. of K supergroup '
             CALL IWRTMA(KSPGRP,1,NGAS,1,NGAS)
             WRITE(6,*) ' Active part of IGENOP'
             IF(ICA.EQ.1) THEN 
               CALL IWRTMA(IGENOP(1,JGENOP),1,NGAS,1,NGAS)
             ELSE 
               CALL IWRTMA(IGENOP(NGAS+1,JGENOP),1,NGAS,1,NGAS)
             END IF
           END IF
*
C  REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
          LEN = 0
          IF(ICA.EQ.1) THEN
            CALL REF_OP(IGENOP(1,JGENOP),IGENOP_EXP,NOP,NGAS,1)
          ELSE
            CALL REF_OP(IGENOP(NGAS+1,JGENOP),IGENOP_EXP,NOP,NGAS,1)
          END IF
          DO IOP = 1, NOP
            IOPTP = IGENOP_EXP(IOP)
C?          WRITE(6,*) ' IOP, IOPTP = ', IOP, IOPTP
*. Number of Kstrings
C                  NST_FOR_OCC(NELEC_PER_GAS,NORB_PER_GAS,NGAS)
            NSTR = NST_FOR_OCC(KSPGRP,NOBPT,NGAS)
            NORB = NOBPT(IOPTP)
            LEN = LEN + NSTR*NORB 
*. Update KSPGRP 
            KSPGRP(IOPTP) = KSPGRP(IOPTP) + 1
          END DO
          MAXLEN = MAX(MAXLEN,LEN)
         END IF
       END DO
      END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Max length of KSTR => ISTR map = ', MAXLEN
      END IF
*
      RETURN
      END
      FUNCTION IOP_SPSTR_ACTIVE(IOP,ICA,KSPGRP)
*
* Is operator IOP times KSPGRP an active supergroup
* IOP and KSPSGRP are both defined as occupation in each gas spaces
*
* IF ICA = 1 then IOP is creation string
* IF ICA = 2 then IOP is annihilation string
*
* Jeppe Olsen, July 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
*. Specific input
      INTEGER IOP(NGAS),KSPGRP(NGAS)
*. Local scratch
      INTEGER ISPGRP(MXPNGAS)
*. Occupation of IOP*KSPGRP
      IONE = 1
      IF(ICA.EQ.1) THEN
        IDELTA = 1
      ELSE
        IDELTA = -1
      END IF
      CALL IVCSUM(ISPGRP,IOP,KSPGRP,IDELTA,IONE,NGAS)
*. Number of electrons in ISPGRP
      NIEL = IELSUM(ISPGRP,NGAS)
*
      IFOUND = 0
      DO JSPGRP = 1, NTSPGP
        I_DENTICAL = 1
        DO IGAS = 1, NGAS
          IF(ISPGRP(IGAS).EQ.NELFSPGP(IGAS,JSPGRP)) I_DENTICAL = 0
        END DO
        IF(I_DENTICAL.EQ.1) IFOUND = 0
      END DO
*
      IOP_SPSTR_ACTIVE = IFOUND 
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' creation operator and K spgrp '
        CALL IWRTMA(IOP,1,NGAS,1,NGAS)
        CALL IWRTMA(KSPGRP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Resulting supergroup '
        CALL IWRTMA(ISPGRP,1,NGAS,1,NGAS)
        IF(IFOUND.EQ.1) THEN 
          WRITE(6,*)' Resulting supergroup is active '
        ELSE
          WRITE(6,*)' Resulting supergroup is not active '
        END IF
      END IF
*
      RETURN
      END 

      SUBROUTINE SPOBEX_TO_ABOBEX(ISPOBEX_TP,NSPOBEX_TP,NGAS,
     &           IFLAG,NAOBEX_TP,NBOBEX_TP,IAOBEX_TP,IBOBEX_TP)
*
* Split spin-orbital excitations into alpha and beta-orbital excitations
*
* IFLAG = 1 : Find only number of alpha- and beta- orbital excitations
*
* Jeppe Olsen, July 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
*. Output
      INTEGER IAOBEX_TP(2*NGAS,*),IBOBEX_TP(2*NGAS,*)
*
C?    WRITE(6,*) ' SPOBEX_TO_ABOBEX : NSPOBEX_TP = ', NSPOBEX_TP
      DO IAB = 1, 2
        LEX_TP = 0
        DO JSPOBEX_TP = 1, NSPOBEX_TP
*. Has this a or b excitation been observed before
          I_AM_OLD_HAT = 0
          DO KSPOBEX_TP = 1, JSPOBEX_TP-1
            I_DENTICAL = 1
            DO ICA = 1, 2
              IEXP_OFF = (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
              DO IGAS = 1, NGAS
                IF(ISPOBEX_TP(IGAS+IEXP_OFF-1,JSPOBEX_TP).NE.
     &             ISPOBEX_TP(IGAS+IEXP_OFF-1,KSPOBEX_TP)    )
     &             I_DENTICAL = 0
              END DO
            END DO
            IF(I_DENTICAL.EQ.1) I_AM_OLD_HAT = 1
          END DO
          IF(I_AM_OLD_HAT.EQ.0) THEN
            LEX_TP = LEX_TP + 1
            IF(IFLAG.EQ.0) THEN
              DO ICA = 1, 2
               IEXP_OFF =  (IAB-1)*NGAS + (ICA-1)*2*NGAS + 1
               IAB_OFF  =  (ICA-1)*NGAS + 1
               IF(IAB.EQ.1) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IAOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               ELSE IF (IAB.EQ.2) THEN
                 CALL ICOPVE(ISPOBEX_TP(IEXP_OFF,JSPOBEX_TP),
     &                        IBOBEX_TP(IAB_OFF,LEX_TP),NGAS)
               END IF
              END DO
            END IF
          END IF
*         ^ End of I_AM_OLD_HAT = 0
        END DO
*       ^ End of loop over JSPOBEX_TP
      IF(IAB.EQ.1) THEN
        NAOBEX_TP = LEX_TP
      ELSE
        NBOBEX_TP = LEX_TP
      END IF
      END DO
*     ^ End of loop over IAB
*
      NTEST = 00
      IF(NTEST.GE.3) THEN
       WRITE(6,*) ' Number of alpha-excitation operators ', NAOBEX_TP
       WRITE(6,*) ' Number of beta-excitation operators ',  NBOBEX_TP
      END IF
      IF(NTEST.GE.5.AND.IFLAG.EQ.0) THEN
        WRITE(6,*) ' Alpha-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NAOBEX_TP  
         WRITE(6,*)
         WRITE(6,*) ' alphaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)') 
     &   ' Creation      :',  (IAOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)') 
     &   ' Annihilation  :',  (IAOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' beta-excitation operators : '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        DO JEX_TP = 1, NBOBEX_TP  
         WRITE(6,*)
         WRITE(6,*) ' betaorbitalexcitation ', JEX_TP
         WRITE(6,'(A,16I4)') 
     &   ' Creation      :',  (IBOBEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
         WRITE(6,'(A,16I4)') 
     &   ' Annihilation  :',  (IBOBEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        END DO
      END IF
*     ^ End of print is active
*
      RETURN
      END
      SUBROUTINE CC3_XXXX(CC3JAC,NDIM,EIGVEC)
*
* Analyze eigenvector EIGVEC of CC3 Jacobian
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CC3JAC(*),EIGVEC(*)
*. Local scratch
      DIMENSION IVEC(100), IDEGPAIR(2,100)
      DIMENSION SUBJAC(100*100),X(100*100),SCR(100,100)
      DIMENSION T(100*100)
*
      WRITE(6,*) ' Welcome to CC3_XXXX  '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) ' Input eigenvector : '
      CALL WRTMAT(EIGVEC,NDIM,1,NDIM,1)
*
*. Elements of eigenvector that are nonvanishing
*
      N_NONZERO = 0
      TEST = 1.0D-10
      DO I = 1, NDIM
        IF(ABS(EIGVEC(I)).GT.TEST) THEN
          N_NONZERO = N_NONZERO + 1
          IVEC(N_NONZERO) = I
         END IF
      END DO
*
      WRITE(6,*) ' Number of nonvanishing components ', N_NONZERO
      WRITE(6,*) ' Nonvanishing components : '
      CALL IWRTMA(IVEC,1,N_NONZERO,1,N_NONZERO) 
*. Jacobian over nonvanishing elements 
      DO I = 1, N_NONZERO
       DO J = 1, N_NONZERO
        IEXP = IVEC(I)
        JEXP = IVEC(J)
        SUBJAC((J-1)*N_NONZERO + I ) = 
     &  CC3JAC((JEXP-1)*NDIM+IEXP)
       END DO
      END DO
*
      WRITE(6,*) ' Jacobian in active subspace '
      CALL WRTMAT(SUBJAC,N_NONZERO,N_NONZERO,N_NONZERO,N_NONZERO)
*. The elements are assumed to occur in degenerate pairs, find these
      SIGN = -1.0D0
      NDEG = 0 
      DO I = 1, N_NONZERO
        IF(IVEC(I).GT.0) THEN
          ELMNT = EIGVEC(IVEC(I))
          DO J = I+1,N_NONZERO 
            IF(IVEC(J).GT.0.AND.
     &         ABS(ELMNT - SIGN*EIGVEC((IVEC(J)))).LE.TEST) THEN
              NDEG = NDEG + 1
              IDEGPAIR(1,NDEG) = I
              IDEGPAIR(2,NDEG) = J
              IVEC(I) = - IVEC(I)
              IVEC(J) = - IVEC(J)
            END IF
          END DO
        END IF
      END DO
      WRITE(6,*) ' Degenerate pairs '
      CALL IWRTMA(IDEGPAIR,2,NDEG,2,NDEG)
      IF(NDEG.NE.N_NONZERO/2) THEN
        WRITE(6,*) ' Problems : not simple degeneracies '
        RETURN
C       STOP       ' Problems : not simple degeneracies '
      END IF
*. Clean up
      DO I = 1, N_NONZERO
        IVEC(I) = - IVEC(I)
      END DO
*. Set up transformations matrix 
      ZERO = 0.0D0
      CALL SETVEC(X,ZERO,N_NONZERO**2)
      SQRT2I = 1.0D0/SQRT(2.0D0)
      SQRT2IM = (-1.0D0)/SQRT(2.0D0)
      DO IDEG = 1, NDEG
*
        I1 = IDEGPAIR(1,IDEG)
        I2 = IDEGPAIR(2,IDEG)
*
        X((IDEG-1)*N_NONZERO + I1) = SQRT2I
        X((IDEG-1)*N_NONZERO + I2) = SQRT2IM
      END DO
*
      WRITE(6,*) ' Transformation matrix '
      CALL WRTMAT(X,N_NONZERO,N_NONZERO,N_NONZERO,N_NONZERO)
*. Transform subspace Jacobian 
C TRAN_SYM_BLOC_MAT2(AIN,X,NBLOCK,LBLOCK,AOUT,SCR,ISYM)
      CALL TRAN_SYM_BLOC_MAT2(SUBJAC,X,1,N_NONZERO,T,SCR,0)
*
      WRITE(6,*) ' Transformed active Jacobian '
      CALL WRTMAT(T,N_NONZERO,N_NONZERO,N_NONZERO,N_NONZERO)
*
      RETURN
      END
      SUBROUTINE CC3_VECFNC(CCVEC,CCDIA,CCAMP,VEC1,VEC2)
*  
* Triples part of CC3 vector function 
* 
* = Diag*T3 coef + <mu!exp(-T1) [H, T2] Exp(T1) |HF>
*
* Jeppe Olsen, March 2000
*
* Output vector CCVEC is assumed to be complete vector - 
* although only triples part is updated
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cc_exc.inc'
C     COMMON/CTCC/KLSOBEX,NSPOBEX_TP,KLLSOBEX,KLIBSOBEX,LEN_T_VEC,
C    &             MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
C    &             KLOBEX,NOBEX_TP,KLSOX_TO_OX,KLSPOBEX_AC
      CHARACTER*8 CCTYPE
*. Input CC amplitudes and diagonal
      DIMENSION CCAMP(*), CCDIA(*)
*. Updated CC vector function
      DIMENSION CCVEC(*)
*. Local scratch for saving single and doubles part
      DIMENSION SCRLOC(1000)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN  
        WRITE(6,*) ' Entering CC3_VECFNC '
        WRITE(6,*) ' =================== '
      END IF
*
      CCTYPE(1:6) = 'GEN_CC'
*. Offset and number of triples
C     GET_IN_OBOEX(IOBEXC,NDIM,IOFF,ISOX_TO_OX,LEN_SOX,NSOX)
      CALL GET_IN_OBOEX(1,NTRIPLE,IBTRIPLE,WORK(KLSOX_TO_OX),
     &     WORK(KLLSOBEX),NSPOBEX_TP)
*. Save Singles and doubles part
      NSINDOU = IBTRIPLE-1
      CALL COPVEC(CCVEC,SCRLOC,NSINDOU)
* Exp T_1 |HF >
      IZERO = 0
      CALL ISETVC(WORK(KLSPOBEX_AC),IZERO,NSPOBEX_TP)
*.Activate single = excitation class 3
      CALL SET_AC_EXTP(3,WORK(KLSPOBEX_AC),WORK(KLSOX_TO_OX),
     &                 NSPOBEX_TP)
*. 
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-20
*.  Exp(t1) !ref> on LUHC
      CALL EXPT_REF2(LUC,LUHC,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &               CCAMP,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
C?    WRITE(6,*) '  Exp(t1) !ref> on LUHC '
C?    CALL WRTVCD(VEC1,LUHC,1,-1)
*. T2 times Exp(t1) !ref> on LUSC1 
*.Activate Doubles = excitation class 2
      CALL ISETVC(WORK(KLSPOBEX_AC),IZERO,NSPOBEX_TP)
      CALL SET_AC_EXTP(2,WORK(KLSPOBEX_AC),WORK(KLSOX_TO_OX),
     &                 NSPOBEX_TP)
C     SIG_GCC(C,HC,LUC,LUHC,T)
      CALL SIG_GCC(VEC1,VEC2,LUHC,LUSC1,CCAMP)
C?    WRITE(6,*) '  T2 Exp(t1) !ref> on LUSC1 '
C?    CALL WRTVCD(VEC1,LUSC1,1,-1)
*.  H T2  exp (T1) |hf> on LUSC3
      ICC_EXC = 0
      CALL MV7(VEC1,VEC2,LUSC1,LUSC3,0,0)
C?    WRITE(6,*) '  H T2 Exp(t1) !ref> on LUSC3 '
C?    CALL WRTVCD(VEC1,LUSC3,1,-1)
*
*.  H Exp(T1) !HF> on LUSC2
      CALL MV7(VEC1,VEC2,LUHC,LUSC2,0,0)
*. T2 H Exp(T1) |HF> on LUSC1
      CALL SIG_GCC(VEC1,VEC2,LUSC2,LUSC1,CCAMP)
* [H,T2] Exp( T1) | HF > on LUSC2
C     VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      FAC1 = 1.0D0
      FAC2 =-1.0D0
      CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUSC3,LUSC1,LUSC2,1,-1) 
*.     Exp(-T1)  [H,T2] Exp( T1) | HF > on LUHC
C?    WRITE(6,*) '  [H,T2] Exp( T1) | HF > : ' 
C?    CALL WRTVCD(VEC1,LUSC2,1,-1)
      CALL ISETVC(WORK(KLSPOBEX_AC),IZERO,NSPOBEX_TP)
      CALL SET_AC_EXTP(3,WORK(KLSPOBEX_AC),WORK(KLSOX_TO_OX),
     &                 NSPOBEX_TP)
*
      ONEM = -1.0D0
      CALL SCALVE(CCAMP,ONEM,N_CC_AMP)
      CALL EXPT_REF2(LUSC2,LUHC,LUSC35,LUSC1,LUSC3,XCONV,MX_TERM,
     &               CCAMP,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
      CALL SCALVE(CCAMP,ONEM,N_CC_AMP)
*
C?    WRITE(6,*) ' Exp(-T1)  [H,T2] Exp( T1) | HF > : ' 
C?    CALL WRTVCD(VEC1,LUHC,1,-1)
*
*. And obtain <mu 3! 
      CALL ISETVC(WORK(KLSPOBEX_AC),IZERO,NSPOBEX_TP)
      CALL SET_AC_EXTP(1,WORK(KLSPOBEX_AC),WORK(KLSOX_TO_OX),
     &                 NSPOBEX_TP)
C            DEN_GCC(C,HC,LUC,LUHC,T)
C       DEN_GCC(VEC1,VEC2,LUC,LUHC,JAC_VEC) 
      CALL DEN_GCC(VEC1,VEC2,LUC,LUHC,CCVEC)
*. And add Diag*CCAMP  
*
      DO ITRIP = IBTRIPLE,IBTRIPLE+NTRIPLE-1
        CCVEC(ITRIP) = CCVEC(ITRIP) + CCDIA(ITRIP)*CCAMP(ITRIP)  
      END DO
*. And copy in single-doubles part
      CALL COPVEC(SCRLOC,CCVEC,NSINDOU)
*. Reset active array 
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TP)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input (CCAMP) and output(CCVEC) from CC3_VEC..'
        CALL WRTMAT(CCAMP,1,N_CC_AMP,1,N_CC_AMP)
        CALL WRTMAT(CCVEC,1,N_CC_AMP,1,N_CC_AMP)
      END IF
*
      RETURN
      END
      SUBROUTINE SET_AC_EXTP(IACTP,ISPOBEX_AC,ISOX_TO_OX,NSPOBEX_TP)
*
* Activate spinorbital excitations belonging to type IACTP
*
* Notice : The other orbital excitations are not zeroed,
*          so several orbital excitation types may be 
*          activated by several calls to this routine.
*
* Jeppe Olsen, March 2000
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISOX_TO_OX(*)
*. Output
      INTEGER ISPOBEX_AC(*) 
*
      DO JSPOBEX_TP = 1, NSPOBEX_TP 
        IF(ISOX_TO_OX(JSPOBEX_TP).EQ.IACTP) THEN
           ISPOBEX_AC(JSPOBEX_TP) = 1
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN  
        WRITE(6,*) ' Orbital excitation type activated = ', IACTP
        WRITE(6,*) ' ISPOBEX_AC :  ' 
        CALL IWRTMA(ISPOBEX_AC,1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*
      RETURN
      END
      SUBROUTINE ACAC_EXC_TYP(IAAEXC,MX_AAEXC,IPRNT)
*
* Information about active-active excitations 
*
* IAAEXC_TYP = 0 => No active-active excitations (Closed shell HF)
*              1 => Active alpha => Active beta (High-spin, Ms = Max)
*              2 => Active beta  => Active alpha(High-spin, Ms = Min)
*              3 => All types of active-active   (General CAS)
*
* Identification of double occupied/valence orbitals is based on 
* info in IHPVGAS
*
* Used for Coupled Cluster Calculations
*
* Find type of referene state 
*
* Jeppe Olsen, March of 2000 
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc' 
      INCLUDE 'strinp.inc'
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
*. Number of hole and valence orbitals
      NHOLE = 0
      NVAL  = 0
      DO IGAS = 1, NGAS
*
       IF(IHPVGAS(IGAS).EQ.1) THEN
*. hole space
          NHOLE = NHOLE + NGSOBT(IGAS)
       ELSE IF(IHPVGAS(IGAS).EQ.3) THEN
*. Valence space
          NVAL = NVAL + NGSOBT(IGAS)
       END IF
      END DO
*
      NEL_AL = NELEC(1)
      NEL_BE = NELEC(2)
      IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE) THEN
*. Closed shell Hartree-Fock
        IAAEXC = 0
        MX_AAEXC = 0
      ELSE IF(NEL_BE.EQ.NHOLE.AND.NEL_AL.EQ.NHOLE+NVAL) THEN
*. High spin open shell case with Max MS
        IAAEXC = 1
        MX_AAEXC = NVAL
      ELSE IF(NEL_AL.EQ.NHOLE.AND.NEL_BE.EQ.NHOLE+NVAL ) THEN 
*. High spin open shell case with Min MS
        IAAEXC = 2
        MX_AAEXC = NVAL
      ELSE 
*. More general, not analyzed in detail p.t. , assumed CAS
        IAAEXC = 3
        MX_AAEXC = NVAL
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
*
        WRITE(6,*) ' Allowed types of active-active excitations '
        WRITE(6,*) ' ========================================== ' 
        WRITE(6,*)
        IF(IAAEXC.EQ.0) THEN
          WRITE(6,*) ' No excitations in active orbital space '
        ELSE IF(IAAEXC.EQ.1) THEN
          WRITE(6,*) ' Alpha => Beta active excitations '
        ELSE IF (IAAEXC.EQ.2) THEN
          WRITE(6,*) ' Beta => alpha active excitations '
        ELSE IF (IAAEXC.EQ.3) THEN 
          WRITE(6,*) ' All active-active excitations '
        END IF
* 
        WRITE(6,'(A,I5)') 
     &  ' Largest excitation level for active-active excitations',
     &    MX_AAEXC
      END IF
*
      RETURN
      END
      SUBROUTINE GET_IN_OBOEX(IOBEXC,NDIM,IOFF,ISOX_TO_OX,LEN_SOX,NSOX)
*
* Number (NDIM) and offset(IOFF) for excitation type IOBEXC
*
      INCLUDE 'implicit.inc'
*. Number of spinorbital excitation per type and spinorbital exc => orb exc
*. INPUT 
      INTEGER LEN_SOX(NSOX),ISOX_TO_OX(NSOX)
*
C?    WRITE(6,*) ' GET_IN ... LEN_SOX,ISOX_TO_OX'
C?    CALL IWRTMA(LEN_SOX,1,NSOX,1,NSOX)
C?    CALL IWRTMA(ISOX_TO_OX,1,NSOX,1,NSOX)
      NDIM = 0
      NTOT = 0
      IOFF = 0
      DO JSOX = 1, NSOX
        LEN = LEN_SOX(JSOX)
        IF(ISOX_TO_OX(JSOX).EQ.IOBEXC) THEN
          IF (IOFF.EQ.0) IOFF = NTOT + 1     
          NDIM = NDIM + LEN
        END IF
        NTOT = NTOT + LEN
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ORBEX_TP, NDIM, IOFF = ', IOBEXC,NDIM,IOFF 
      END IF
*
      RETURN
      END
      SUBROUTINE ZERO_BLK(AMAT,NROWT,NCOLT,IROW_ZERO,ICOL_ZERO,
     &                    LROW_ZERO,LCOL_ZERO)
*
* A matrix of DIM NROWT*NCOLT is given. Zero subblock
* containing LROW_ZERO rows, LCOL_ZERO columns 
* starting from row IROW_ZERO, columns ICOL_ZERO
*
* Jeppe Olsen, Feb. 2000
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION AMAT(NROWT,NCOLT)
*
      DO ICOL = ICOL_ZERO,ICOL_ZERO+LCOL_ZERO-1
        DO IROW = IROW_ZERO,IROW_ZERO+LROW_ZERO-1
          AMAT(IROW,ICOL) = 0.0D0
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Matrix with zero block '
        CALL WRTMAT(AMAT,NROWT,NCOLT,NROWT,NCOLT)
      END IF
*
      RETURN
      END
      SUBROUTINE CC3_JACO(CC_AMP,VEC1,VEC2,CCVEC)
*
* Construct and diagonalize CC3 Jacobiant
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
*
C     COMMON/CTCC/KLSOBEX,NSPOBEX_TP,KLLSOBEX,KLIBSOBEX,LEN_T_VEC,
C    &             MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
C    &             KLOBEX,NOBEX_TP,KLSOX_TO_OX


*
      PARAMETER(MAXDIM = 100)
      DIMENSION VECIN(MAXDIM),VECOUT(MAXDIM)
      REAL*8 JACOR(MAXDIM*MAXDIM), AMAT(MAXDIM*MAXDIM) 
      REAL*8 JACO(MAXDIM*MAXDIM)
      REAL*8 BMAT(MAXDIM*MAXDIM), CMAT(MAXDIM*MAXDIM) 
      REAL*8 VECA(MAXDIM),VECB(MAXDIM),VECC(MAXDIM),VECD(4*MAXDIM)
*. Converged CC amplitudes
      DIMENSION CC_AMP(*)
*. Scratch vectors for CI
      DIMENSION VEC1(*), VEC2(*)
*. Scratch vector for cc amplitudes 
      DIMENSION CCVEC(*)
*
*
      NTEST = 00
      IF(N_CC_AMP.GT.MAXDIM) THEN
        WRITE(6,*) ' N_CC_AMP > MAXDIM , N_CC_AMP, MAXDIM= ',
     &   N_CC_AMP, MAXDIM
        STOP ' N_CC_AMP > MAXDIM'
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Set of CC amplitudes '
        CALL WRTMAT_EP(CC_AMP,1,N_CC_AMP,1,N_CC_AMP)
      END IF
*
*. Construct jacobiant from a sequence of finite difference 
*. right transformations 
*
*. zero triples contributions to reference function 
*. Offset and number of triples
      CALL GET_IN_OBOEX(1,NTRIPLE,IBTRIPLE,WORK(KLSOX_TO_OX),
     &     WORK(KLLSOBEX),NSPOBEX_TP)
*. Off set and number for singles
      CALL GET_IN_OBOEX(3,NSINGLE,IBSINGLE,WORK(KLSOX_TO_OX),
     &     WORK(KLLSOBEX),NSPOBEX_TP)
*
C?    ZERO = 0.0D0
C?    CALL SETVEC(CC_AMP(IBTRIPLE),ZERO,NTRIPLE)
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
C?    WRITE(6,*) ' Triples in reference Zeroed ' 
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Set of CC amplitudes '
        CALL WRTMAT_EP(CC_AMP,1,N_CC_AMP,1,N_CC_AMP)
      END IF
*
      ZERO = 0.D0
      DO IAMP = 1, N_CC_AMP 
        CALL SETVEC(VECIN,ZERO,N_CC_AMP)
        VECIN(IAMP) = 1.0D0
C       CALL JAC_T_VEC(2,CC_AMP,VECOUT,VECIN,VEC1,VEC2,CCVEC)
C            JAC_T_VEC_FUSK(CC_AMP,JACVEC,TVEC,VEC1,VEC2)
        CALL JAC_T_VEC_FUSK(CC_AMP,VECOUT,VECIN,VEC1,VEC2)
C?      WRITE(6,*) ' Enforced stop after first call to JAC_FUSK'
C?      STOP       ' Enforced stop after first call to JAC_FUSK'
        CALL COPVEC(VECOUT,JACOR(1+(IAMP-1)*N_CC_AMP),N_CC_AMP)
      END DO

*
      WRITE(6,*)
      WRITE(6,*) ' ============='
      WRITE(6,*) ' CC3 Jacobian '
      WRITE(6,*) ' ============='
      WRITE(6,*)
      CALL WRTMAT_EP(JACOR,N_CC_AMP,N_CC_AMP,N_CC_AMP,N_CC_AMP)
*. zero all blocks connected with triples, except triples-triplesS
      I_ZERO_TRIP = 0
      IF(I_ZERO_TRIP.EQ.1) THEN
*. For zeroing all triple coupling blocks
C       DO I = 1, IBTRIPLE-1
C       DO J = IBTRIPLE, N_CC_AMP     
*. For zeroing triple-single coupling blocks
C       WRITE(6,*) ' NSINGLE = ', NSINGLE
C       DO I = 1, NSINGLE
C       DO J = IBTRIPLE, N_CC_AMP     
*. For zeroing triple-double coupling blocks
        DO I = NSINGLE + 1, IBTRIPLE-1
        DO J = IBTRIPLE, N_CC_AMP
          IJ = (J-1)*N_CC_AMP + I
          JI = (I-1)*N_CC_AMP + J
C!        JACOR(IJ) = 0.0D0
          JACOR(JI) = 0.0D0
        END DO
        END DO
*
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        WRITE(6,*) ' CC3 jacobiant with S-T couplings zeroed '
        CALL WRTMAT(JACOR,N_CC_AMP,N_CC_AMP,N_CC_AMP,N_CC_AMP)
      END IF
      CALL COPVEC(JACOR,JACO,N_CC_AMP**2)
C     EIGGMTN(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,
C    &                   Z,W,SCR)
      CALL EIGGMTN(JACOR,N_CC_AMP,VECA,VECB,AMAT,BMAT,
     &            CMAT,VECC,VECD)
C     CC3_XXXX(CC3JAC,NDIM,EIGVEC)
*. Find part of Jacobian relevant for given eigenvector
      ISOL = 1
      IF(ISOL.NE.0) THEN
        DO IEIG = 1, N_CC_AMP
          WRITE(6,*) 
          WRITE(6,*) ' Analysis of Subspace Jac for eigenvalue=',
     &    VECA(IEIG)
          WRITE(6,*) 
     &    ' =================================================='
          CALL CC3_XXXX(JACO,N_CC_AMP,AMAT((IEIG-1)*N_CC_AMP+1))
       END DO
      END IF
    
*. Sort real parts of eigenvalues
*
      RETURN
      END
      SUBROUTINE CONJ_CCAMP(CCIN,ISM,CCOUT)
*
* Conjugate a set of CC amplitudes. The input amplitudes 
* are assumed defined by ctcc. So CONJ_CCAMP should be called
* before these are conjugated
*
* Jeppe Olsen, Febr. 11, 3 a'clock listening to Sticky Fingers.
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
C     COMMON/CTCC/KLSOBEX,NSPOBEX_TP,KLLSOBEX,KLIBSOBEX,LEN_T_VEC,
C    &             MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK
*. Input 
      DIMENSION CCIN(*)
*. Output
      DIMENSION CCOUT(*)
*
      CALL CONJ_CCAMP_S(WORK(KLSOBEX),NSPOBEX_TP,CCIN,ISM,CCOUT)
*
      RETURN
      END 
      SUBROUTINE CONJ_CCAMP_S(ITSS_TP,NTSS_TP,CCIN,ISM,CCOUT)
*
* Slave routine for conjugating CC amplitudes
*
* Jeppe Olsen, Febr.11 2000 
*
* Adapted for MSCOMB_CC.NE.0 (after a hard fight)
*
* ak, 10-03-2004
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER ITSS_TP(4*NGAS,NTSS_TP)
      DIMENSION CCIN(*)
*. Output
      DIMENSION CCOUT(*)
*. Local scratch
      INTEGER NCA(MXPNGAS),NCB(MXPNGAS) 
      INTEGER NAA(MXPNGAS),NAB(MXPNGAS)
      INTEGER IBT(8,8,8)
*
      IT_OUT = 0
      DO ITSS = 1, NTSS_TP
*. Number of strings per sym , C/A refers to input matrix
       CALL NST_SPGP(ITSS_TP(1+0*NGAS,ITSS),NCA)
       CALL NST_SPGP(ITSS_TP(1+1*NGAS,ITSS),NCB)
       CALL NST_SPGP(ITSS_TP(1+2*NGAS,ITSS),NAA)
       CALL NST_SPGP(ITSS_TP(1+3*NGAS,ITSS),NAB)
       IF(MSCOMB_CC.NE.0) THEN
        CALL DIAG_EXC_CC(ITSS_TP(1+0*NGAS,ITSS),
     &                   ITSS_TP(1+1*NGAS,ITSS),
     &                   ITSS_TP(1+2*NGAS,ITSS),
     &                   ITSS_TP(1+3*NGAS,ITSS),NGAS,IDIAG)
       ELSE
        IDIAG = 0
       END IF
*. Offsets for symmetryblocks of unconjugated CC amplitudes 
*. relative to start of block
C           Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
       CALL Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ISM,NSMST,IDIAG)
cnew       CALL Z_TCC_OFF(IBT,NAA,NAB,NCA,NCB,ISM,NSMST,IDIAG)
*. Start of this excitation block relative to start of CC amplitudes
       IOFF_BL = IT_OUT+1
*. Loop over symmetryblocks in OUTPUT order
       DO ISM_A = 1, NSMST
         ISM_C = MULTD2H(ISM,ISM_A) 
         DO ISM_AA = 1, NSMST
           ISM_AB = MULTD2H(ISM_A,ISM_AA)
           DO ISM_CA = 1, NSMST
            ISM_CB =  MULTD2H(ISM_C,ISM_CA)
            ISM_AL_O = (ISM_CA-1)*NSMST + ISM_AA
            ISM_BE_O = (ISM_CB-1)*NSMST + ISM_AB
            IF (IDIAG.EQ.1.AND.ISM_AL_O.LT.ISM_BE_O) CYCLE
            IF (IDIAG.EQ.1.AND.ISM_AL_O.EQ.ISM_BE_O) THEN
              IRESTR = 1
            ELSE
              IRESTR = 0
            END IF

            ISM_AL_I = (ISM_AA-1)*NSMST + ISM_CA
            ISM_BE_I = (ISM_AB-1)*NSMST + ISM_CB
            IF (IDIAG.EQ.0.OR.ISM_AL_I.GE.ISM_BE_I) THEN
*. Offset for this block in INPUT matrix
              IBT_OFF = IBT(ISM_CA,ISM_CB,ISM_AA)
              IBT_OFF2 = IBT(ISM_CA,ISM_CB,ISM_AA)
              ITRNSP = 0
            ELSE
              IBT_OFF = IBT(ISM_CB,ISM_CA,ISM_AB)
              IBT_OFF2 = IBT(ISM_CA,ISM_CB,ISM_AA)
              ITRNSP = 1
            END IF

            ! to be skipped?
            IF (IBT_OFF.LT.0) THEN
              WRITE (6,*) 'unexpected event in CONJ_CCAMP_S'
              STOP 'unexpected event in CONJ_CCAMP_S'
            END IF
            IOFF_IN = IBT_OFF + IOFF_BL - 1

*. Numbers of strings for this sym, input notation
            LAA = NAA(ISM_AA)
            LAB = NAB(ISM_AB)
            LCA = NCA(ISM_CA)
            LCB = NCB(ISM_CB)

*. depending on irestr and itrnsp, we have three cases to process:
            IF ( IRESTR.EQ.0 .AND. ITRNSP.EQ.0 ) THEN
*. Loop over T elements in output block as matrix T(I_AA,I_AB,I_CA,I_CB)
              DO I_CB = 1, LCB         
                DO I_CA = 1, LCA           
                  DO I_AB = 1, LAB    
                    DO I_AA = 1, LAA        
*. Output address 
                 IT_OUT = IT_OUT + 1
*. Address of the corresponding input element T(I_CA,I_CB,I_AA,I_AB)
                 IT_IN = IOFF_IN -1 +
     &           (I_AB-1)*LAA*LCB*LCA + (I_AA-1)*LCB*LCA +
     &           (I_CB-1)*LCA + I_CA
                 CCOUT(IT_OUT) = CCIN(IT_IN)

                    END DO
                  END DO
                END DO
              END DO
*             ^ End of loop over elements of block
            ELSE IF ( IRESTR.EQ.0 .AND. ITRNSP.EQ.1 ) THEN
*. Loop over T elements in output block as matrix T(I_AA,I_AB,I_CA,I_CB)
              DO I_CB = 1, LCB         
                DO I_CA = 1, LCA           
                  DO I_AB = 1, LAB    
                    DO I_AA = 1, LAA        
*. Output address 
                 IT_OUT = IT_OUT + 1
*. Address of the corresponding input element T(I_CB,I_CA,I_AB,I_AA)
                 IT_IN = IOFF_IN -1 +
     &           (I_AA-1)*LAB*LCA*LCB + (I_AB-1)*LCA*LCB +
     &           (I_CA-1)*LCB + I_CB
                 CCOUT(IT_OUT) = CCIN(IT_IN)

                    END DO
                  END DO
                END DO
              END DO
*             ^ End of loop over elements of block
            ELSE ! IRESTR.EQ.1 (ITRNSP does not matter)
*. Loop over T elements in output block as matrix T(I_AA,I_AB,I_CA,I_CB)
              DO I_CB = 1, LCB         
                DO I_CA = 1, LCA
                  DO I_AB = 1, LAB    
                    DO I_AA = 1, LAA
        
                 IF(I_CA.LT.I_CB .OR. (I_CA.EQ.I_CB.AND.I_AA.LT.I_AB))
     &                 CYCLE
*. Output address 
                 IT_OUT = IT_OUT + 1
*. Address of the corresponding input element T(I_CA,I_CB,I_AA,I_AB)
                 IF (I_AA.GT.I_AB.OR.
     &                (I_AA.EQ.I_AB.AND.I_CA.GE.I_CB)) THEN
                   IT_IN = IOFF_IN - 1 +
     &                ITDIANUM(I_CA,I_CB,I_AA,I_AB,LCA,LAA)
                 ELSE 
                   IT_IN = IOFF_IN - 1 +
     &                ITDIANUM(I_CB,I_CA,I_AB,I_AA,LCA,LAA)
                 END IF
                 CCOUT(IT_OUT) = CCIN(IT_IN)

                    END DO
                  END DO
                END DO
              END DO
*             ^ End of loop over elements of block
            END IF ! IRESTR
           END DO
*          ^ End of loop over ISM_AA
         END DO
*        ^ End of loop over ISM_CA
       END DO
*      ^ End of loop over ISM_C
      END DO
*     ^ End of loop over ITSS
*
      RETURN
      END
* END OF CLONE 2
      SUBROUTINE JAC_T_VEC_FUSK(CC_AMP,JACVEC,TVEC,VEC1,VEC2)
*
*. Fusk calculation of Jacobiant times vector : Finite difference approach
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
*. Input
      DIMENSION CC_AMP(*),TVEC(*)
*. Output
      REAL*8    JACVEC(*)
*. Scratch
      DIMENSION VEC1(*), VEC2(*)
*. Local scratch
      PARAMETER(MAXDIM = 1000)
      DIMENSION XX(MAXDIM)
*
      CHARACTER*6 CCTYPE 
*
      WRITE(6,*) ' Not working anymore, change call to CC_VEC'
      STOP       ' Not working anymore, change call to CC_VEC'
          
      CCTYPE(1:6) = 'GEN_CC'
*. Calculate CC-vector function with amplitudes CC_AMP+alpha*TVEC
*. Alpha
      ALPHA = 1.0D-3
      ONE = 1.0D0
      CALL VECSUM(CC_AMP,CC_AMP,TVEC,ONE,ALPHA,N_CC_AMP)
C          CC_VEC_FNC(CC_AMP,CC_VEC,E_CC,E_CC_A,VEC1,VEC2,IBIO,CCTYPE)
      CALL CC_VEC_FNC(CC_AMP,XX,E_CC,E_CC2,VEC1,VEC2,1,CCTYPE)
*. 2*alpha 
      FACTOR = ALPHA
      CALL VECSUM(CC_AMP,CC_AMP,TVEC,ONE,FACTOR,N_CC_AMP)
      CALL CC_VEC_FNC(CC_AMP,XX,E_CC,E_CC2,VEC1,VEC2,1,CCTYPE)
*
      FACTORA =  8.0D0/(12.0D0*ALPHA)
      FACTOR2A = (-1.0D0)/(12.0D0*ALPHA)
      CALL VECSUM(JACVEC,JACVEC,XX,FACTORA,FACTOR2A,N_CC_AMP)
*. -alpha 
      FACTOR = (-3.0D0)*ALPHA
      CALL VECSUM(CC_AMP,CC_AMP,TVEC,ONE,FACTOR,N_CC_AMP)
      CALL CC_VEC_FNC(CC_AMP,XX,E_CC,E_CC2,VEC1,VEC2,1,CCTYPE)
      FACTORAM = - FACTORA
      CALL VECSUM(JACVEC,JACVEC,XX,ONE,FACTORAM,N_CC_AMP)
*. -2*alpha
      FACTOR = -ALPHA
      CALL VECSUM(CC_AMP,CC_AMP,TVEC,ONE,FACTOR,N_CC_AMP)
      CALL CC_VEC_FNC(CC_AMP,XX,E_CC,E_CC2,VEC1,VEC2,1,CCTYPE)
      FACTOR2AM = 1.0D0/(12.0D0*ALPHA)
      CALL VECSUM(JACVEC,JACVEC,XX,ONE,FACTOR2AM,N_CC_AMP)
*. Clean up
      FACTOR = 2.0D0*ALPHA
      CALL VECSUM(CC_AMP,CC_AMP,TVEC,ONE,FACTOR,N_CC_AMP)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Trial vector and Jacobiant times trial vector '
        WRITE(6,*) ' (supplied from FUSK code) '
        CALL WRTMAT(TVEC,1,N_CC_AMP,1,1)
        CALL WRTMAT(JACVEC,1,N_CC_AMP,1,1)
      END IF
*
      RETURN
      END
      SUBROUTINE JAC_T_VEC_NUM(CC_AMP,CC_VEC,VEC1,VEC2,
     &                         LUCCAMP,LU_RAMP,LU_SIG,LU_SIG_NUM)
*
*. Numerical calculation of Jacobiant times vector : Finite difference approach
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
*. Input
      DIMENSION CC_AMP(*),CC_VEC(*)
*. Scratch
      DIMENSION VEC1(*), VEC2(*)
*
      CHARACTER*6 CCTYPE 
      REAL(8) INPROD
*
      CCTYPE(1:6) = 'GEN_CC'
*. Calculate CC-vector function with amplitudes CC_AMP+alpha*TVEC
*. Alpha
      LUSCRAMP = IOPEN_UUS()
      LUSCRAMP2 = IOPEN_UUS()
      LUSCRAMP3 = IOPEN_UUS()
      LUSCRAMP4 = IOPEN_UUS()
      LUSCRAMP5 = IOPEN_UUS()

      ALPHA = 1.0D-4
      WRITE(6,*) 'NUMERICAL JACOBIAN TIMES VECTOR PRODUCT:'
      WRITE(6,*) ' 2-POINT CENTRAL DIFFERENCE FORMULA, INCR.=',ALPHA
      WRITE(6,*) ' SCRATCH UNITS ARE: ',LUSCRAMP,LUSCRAMP2,LUSCRAMP3,
     &     LUSCRAMP4,LUSCRAMP5


      WRITE(6,*) 'OUTPUT FOR + INCREMENT:'
      ONE = 1.0D0
      CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,-1,LUCCAMP)
      CALL VEC_FROM_DISC(CC_VEC,N_CC_AMP,1,-1,LU_RAMP)

      XNORM = SQRT(INPROD(CC_VEC,CC_VEC,N_CC_AMP,1,-1))
      WRITE(6,*) 'NORM OF TRIAL VECTOR: ',XNORM

      CALL VECSUM(CC_AMP,CC_AMP,CC_VEC,ONE,ALPHA,N_CC_AMP)
      CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,-1,LUSCRAMP)
C          CC_VEC_FNC(CC_AMP,CC_VEC,E_CC,E_CC_A,VEC1,VEC2,IBIO,CCTYPE)
      CALL CC_VEC_FNC2(CC_AMP,CC_VEC,E_CC,E_CC2,
     &                 XDUM1,XDUM2,XDUM3,
     &                 VEC1,VEC2,1,CCTYPE,
     &                 DUM,
     &                 LUSCRAMP,LUSCRAMP2,-202,!LUDUM,
     &                 LUSCRAMP4,LUSCRAMP5,-101)!LUDUM2)

      WRITE(6,*) 'OUTPUT FOR - INCREMENT:'
      CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,-1,LUCCAMP)
      CALL VEC_FROM_DISC(CC_VEC,N_CC_AMP,1,-1,LU_RAMP)
      CALL VECSUM(CC_AMP,CC_AMP,CC_VEC,ONE,-ALPHA,N_CC_AMP)
      CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,-1,LUSCRAMP)
C          CC_VEC_FNC(CC_AMP,CC_VEC,E_CC,E_CC_A,VEC1,VEC2,IBIO,CCTYPE)
      CALL CC_VEC_FNC2(CC_AMP,CC_VEC,E_CC,E_CC2,
     &                 XDUM1,XDUM2,XDUM3,
     &                 VEC1,VEC2,1,CCTYPE,
     &                 DUM,
     &                 LUSCRAMP,LUSCRAMP3,LUDUM,
     &                 LUSCRAMP4,LUSCRAMP5,LUDUM2)

      CALL VEC_FROM_DISC(CC_AMP,N_CC_AMP,1,-1,LUSCRAMP2)
      CALL VEC_FROM_DISC(CC_VEC,N_CC_AMP,1,-1,LUSCRAMP3)

      FAC = 1d0/(2d0*ALPHA)
      CALL VECSUM(CC_AMP,CC_AMP,CC_VEC,FAC,-FAC,N_CC_AMP)

      CALL VEC_FROM_DISC(CC_VEC,N_CC_AMP,1,-1,LU_SIG)

      WRITE(6,*) 'COMPARING WITH RESULT ON UNIT ',LU_SIG
      WRITE(6,*) ' col 1 -- numerical'
      WRITE(6,*) ' col 2 -- analytical'

      CALL CMP2VC(CC_AMP,CC_VEC,N_CC_AMP,100d0*ALPHA**2)

      WRITE(6,*) 'SAVING NUMERICAL RESULT ON UNIT ',LU_SIG_NUM

      CALL VEC_TO_DISC(CC_AMP,N_CC_AMP,1,-1,LU_SIG)      

      CALL RELUNIT(LUSCRAMP, 'delete')
      CALL RELUNIT(LUSCRAMP2,'delete')
      CALL RELUNIT(LUSCRAMP3,'delete')
      CALL RELUNIT(LUSCRAMP4,'delete')
      CALL RELUNIT(LUSCRAMP5,'delete')
*
      RETURN
      END
      SUBROUTINE JACO_TEST(CC_AMP,VEC1,VEC2,CCVEC)
*
* Test of general jacobiant times vector
*
      INCLUDE 'implicit.inc'
*
      PARAMETER(MAXDIM = 100)
      DIMENSION VECIN(MAXDIM),VECOUT(MAXDIM)
      REAL*8 JACOR(MAXDIM*MAXDIM), JACOL(MAXDIM*MAXDIM) 
      REAL*8 JACOD(MAXDIM*MAXDIM)
*
      REAL*8 AMAT(MAXDIM*MAXDIM) 
      REAL*8 BMAT(MAXDIM*MAXDIM), CMAT(MAXDIM*MAXDIM) 
      REAL*8 VECA(MAXDIM),VECB(MAXDIM),VECC(MAXDIM),VECD(4*MAXDIM)

*. Converged CC amplitudes
      DIMENSION CC_AMP(*)
*. Scratch vectors for CI
      DIMENSION VEC1(*), VEC2(*)
*. Scratch vector for cc amplitudes 
      DIMENSION CCVEC(*)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
*
      NTEST = 100
      IF(N_CC_AMP.GT.MAXDIM) THEN
        WRITE(6,*) ' N_CC_AMP > MAXDIM , N_CC_AMP, MAXDIM= ',
     &   N_CC_AMP, MAXDIM
        STOP ' N_CC_AMP > MAXDIM'
      END IF
*
C?    ZERO = 0.0D0
C?    CALL SETVEC(CC_AMP,ZERO,N_CC_AMP)
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
C?    WRITE(6,*) ' CC_amp set to zero '
      
*
*. Construct jacobiant from a sequence of right transformations 
*
      ZERO = 0.0D0
      DO IAMP = 1, N_CC_AMP 
        CALL SETVEC(VECIN,ZERO,N_CC_AMP)
        VECIN(IAMP) = 1.0D0
C            JAC_T_VEC(L_OR_R,CC_AMP,JAC_VEC,TVEC,VEC1,VEC2,CCVEC)
        CALL JAC_T_VEC(2,CC_AMP,VECOUT,VECIN,VEC1,VEC2,CCVEC)
        CALL COPVEC(VECOUT,JACOR(1+(IAMP-1)*N_CC_AMP),N_CC_AMP)
      END DO
*
      IF(NTEST .GE.100 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================================='
        WRITE(6,*) ' CC Jacobiant from JACO_TEST, rhs trans '
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        CALL WRTMAT(JACOR,N_CC_AMP,N_CC_AMP,N_CC_AMP,N_CC_AMP)
      END IF
*

*
*. Construct jacobiant from a sequence of left transformations 
*
      ZERO = 0.0D0
      DO IAMP = 1, N_CC_AMP 
        CALL SETVEC(VECIN,ZERO,N_CC_AMP)
        VECIN(IAMP) = 1.0D0
        CALL JAC_T_VEC(1,CC_AMP,VECOUT,VECIN,VEC1,VEC2,CCVEC)
        DO J = 1, N_CC_AMP
          JACOL(IAMP + (J-1)*N_CC_AMP) = VECOUT(J)
        END DO
      END DO
*
      IF(NTEST .GE.100 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================================='
        WRITE(6,*) ' CC Jacobiant from JACO_TEST, lhs trans '
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        CALL WRTMAT(JACOL,N_CC_AMP,N_CC_AMP,N_CC_AMP,N_CC_AMP)
      END IF
*
*. Construct Jacobian from a sequence of Finite difference calculations
*
C?    ZERO = 0.0D0
C?    DO IAMP = 1, N_CC_AMP 
C?      CALL SETVEC(VECIN,ZERO,N_CC_AMP)
C?      VECIN(IAMP) = 1.0D0
C?      CALL JAC_T_VEC_FUSK(CC_AMP,VECOUT,VECIN,VEC1,VEC2)
C            JAC_T_VEC_FUSK(CC_AMP,JAC_VEC,TVEC,VEC1,VEC2)
C?      CALL COPVEC(VECOUT,JACOD(1+(IAMP-1)*N_CC_AMP),N_CC_AMP)
C?    END DO
*
C?    IF(NTEST .GE.100 ) THEN
C?      WRITE(6,*)
C?      WRITE(6,*) ' ========================================'
C?      WRITE(6,*) ' Jacobiant obtained as finite difference '
C?      WRITE(6,*) ' ========================================'
C?      WRITE(6,*)
C?      CALL WRTMAT(JACOD,N_CC_AMP,N_CC_AMP,N_CC_AMP,N_CC_AMP)
C?    END IF
*. 
*
C?    WRITE(6,*) ' Comparison of RHS and FD Jacobian '
C?    XDIFF = 1.0D-9  
C?    CALL CMP2VC(JACOR,JACOD,N_CC_AMP**2,XDIFF)
*
      CALL EIGGMTN(JACOR,N_CC_AMP,VECA,VECB,AMAT,BMAT,
     &            CMAT,VECC,VECD)

      RETURN
      END
      SUBROUTINE JAC_T_VEC(L_OR_R,CC_AMP,JAC_VEC,TVEC,VEC1,VEC2,
     &                     CCVEC)
*
* Jeppe + Jesper, Dec. 1 1999
*. Symmetry added to RHS, Summer of 2000
*
* L_OR_R Flags whether left- or right- transformation is carried out
*
* Right Transformation
*     <hf!tau+_{mu} exp -T [H, sum(nu) tau_nu V_nu] exp T |hf>
* Left transformation
*
* Symmetry of tau-operator is given by ITEX_SM in /CC_EXC/   
*
*
*. Input
*
*  TVEC   : CC amplitudes defining Linear transformation
*  CC_AMP : CC amplitudes of current CC state
*
*. Output
*
* JACVEC : Jacobiant times vector 
* 
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc' 
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cstate.inc'
*. Input
      DIMENSION CC_AMP(*), TVEC(*)
*. Output
      REAL*8  JAC_VEC(*)
*. Scratch
      DIMENSION VEC1(*), VEC2(*), CCVEC(*)
*
      CHARACTER*6 CCTYPE
*
      NTEST = 1
C?    WRITE(6,*) ' In JAC ..... '
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'JAC_VE')
*. Works only for generalized cc
      CCTYPE(1:6) = 'GEN_CC'
*
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-20
*
      IF(L_OR_R.EQ.2) THEN
*
*                  ======================
*                   Right transformation
*                  ======================
*
*     <hf!tau+_{mu} exp -T [H, sum(nu) tau_nu V_nu] exp T |hf>
*.  Exp(t) !ref> on LUHC
        ICSM = IREFSM
        ISSM = IREFSM
*
        ICSPC = IETSPC
        ISSPC = IETSPC
*
        CALL EXPT_REF2(LUC,LUHC,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &               CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
*
*. Term 1 : H Sum(nu) tau_nu V_nu on LUSC34 
*
*. a : tau_nu V_nu exp (T) |hf> on LUSC35
C       SIG_GCC(C,HC,LUC,LUHC,T)
        ICSM = IREFSM
        ISSM = MULTD2H(ITEX_SM,ICSM)
        CALL SIG_GCC(VEC1,VEC2,LUHC,LUSC35,TVEC)
*. b : H *  tau_nu V_nu exp (T) |hf>
        ICC_EXC = 0
        ICSM = ISSM
*
        ICSPC = IETSPC
        ISSPC = ITSPC
        CALL MV7(VEC1,VEC2,LUSC35,LUSC34,0,0)
C       WRITE(6,*) ' H V_op exp(T) !HF > '
C       CALL WRTVCD(VEC1,LUSC34,1,-1)
*
*. Term 2 : Sum(nu) tau_nu V_nu H exp(T) |hf> on LUSC36 
*
        ICC_EXC = 0
        ICSM = IREFSM
        ISSM = IREFSM
        ICSPC = IETSPC
        ISSPC = ITSPC
        CALL MV7(VEC1,VEC2,LUHC,LUSC51,0,0)
        ICC_EXC = 1
        ISSM = MULTD2H(ITEX_SM,ICSM)
*

        ICSPC = ITSPC
        ISSPC = ITSPC
        CALL SIG_GCC(VEC1,VEC2,LUSC51,LUSC36,TVEC)
C       WRITE(6,*) '  V_op H exp(T) !HF > '
C       CALL WRTVCD(VEC1,LUSC36,1,-1)
*.     (H Sum(nu) tau_nu V_nu -  sum(nu) tau_nu V_nu exp (T) |hf> on LUSC35
        FAC1 = 1.0D0
        FAC2 =-1.0D0
C           VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUSC34,LUSC36,LUSC35,1,-1)
C?    CALL WRTVCD(WORK(KVEC1),LUSC1,1,-1)
C?      WRITE(6,*) ' Input vector to exp(-T) '
C?      CALL WRTVCD(VEC1,LUSC35,1,-1) 
* Exp(-T) [H,sum_nu tau_nu V_nu] Exp(T) |hf>  
        ONEM = -1.0D0
        CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
        ICSM = ISSM
        ICSPC = ITSPC
        ISSPC = ITSPC
        CALL EXPT_REF2(LUSC35,LUHC,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &              CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0) 
        CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
*<hf|tau+_mu  Exp(-T) [H,sum_nu tau_nu V_nu] Exp(T) |hf> =
*<LUC| tau+mu |LUHC> = <LUHC| tau_mu |LU> as density matrix
C            DEN_GCC(C,HC,LUC,LUHC,T)
        ICSM = IREFSM
        ICSPC = IETSPC
        ISSPC = ITSPC
        CALL DEN_GCC(VEC1,VEC2,LUC,LUHC,JAC_VEC) 
      ELSE
*
*                  =====================
*                   Left transformation
*                  =====================
*
* Term 1 : 
* sum(mu)<hf|(tau_mu+ V_mu)exp -T H  tau_nu exp T|hf>
*
*  Calculate as <1| tau_nu exp T |hf>
*  with |1> = H (exp -T)+ (sum(mu) V_mu tau_mu)|hf>
*
*. A: (sum(mu) V_mu tau_mu)|hf> on LUSC34
*
         WRITE(6,*) ' Left Transformation '
C             SIG_GCC(C,HC,LUC,LUHC,T)
         call SIG_GCC(VEC1,VEC2,LUC,LUSC34,TVEC)
        ICSPC = IETSPC
        ISSPC = IETSPC
*. B: exp(-T)+ A on LUSC35
         ONEM = -1.0D0
         CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
*. Conjugate amplitudes (reorder)  and operators
C             CONJ_CCAMP(CCIN,ISM,CCOUT)
* we may change the next lines a little as EXPT_REF now accepts any CCAMP:
         CALL CONJ_CCAMP(CC_AMP,1,CCVEC)
         CALL COPVEC(CCVEC,CC_AMP,N_CC_AMP)
         CALL CONJ_T
*
         CALL EXPT_REF2(LUSC34,LUSC35,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &              CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
         CALL CONJ_CCAMP(CC_AMP,1,CCVEC)
         CALL COPVEC(CCVEC,CC_AMP,N_CC_AMP)
         CALL CONJ_T
         CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
*. C: H B on LUSC36
         ICC_EXC = 0
         CALL MV7(VEC1,VEC2,LUSC35,LUSC36,0,0)
*. D: exp(T)|HF> on LUHC
         ICC_EXC = 1
         CALL EXPT_REF2(LUC,LUHC,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &              CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
*. E: Complex konjugate of D
CJO      CALL CONJ_T
*. F: sum(mu)<hf|(tau_mu+ V_mu)exp -T H exp T tau_nu|hf>
CJO*. <LUHC|tau_nu|LUSC36> on LUSC34
CJO*. <LUSC36|tau_nu|LUHC> in CCVEC  
CJO   <LUSC36|tau_nu|LUHC> is a density, not a CI-vector and is therefore
CJO   returned in the T-array of DEN_GCC
C             DEN_GCC(C,HC,LUC,LUHC,T)
         CALL DEN_GCC(VEC1,VEC2,LUHC,LUSC36,CCVEC)
*. G: H exp(T)|HF> on LUSC36
         CALL EXPT_REF2(LUC,LUSC51,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &               CC_AMP,DUM,VEC1,VEC2,N_CC_AMP,CCTYPE,0)
         ICC_EXC = 0
         CALL MV7(VEC1,VEC2,LUSC51,LUSC36,0,0)
*. H: sum_mu<HF|(tau_mu+V_mu)exp(-T) tau_nu H exp(T)|HF>
*.    <LUSC35|tau_nu|LUSC36> on LUSC51
C        CALL CONJ_T
         CALL DEN_GCC(VEC1,VEC2,LUSC36,LUSC35,JAC_VEC)
*. F-H
         ONE = 1.0D0
         ONEM = -1.0D0
         CALL VECSUM(JAC_VEC,JAC_VEC,CCVEC,ONEM,ONE,N_CC_AMP)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '********************'
        WRITE(6,*) ' JAC_T_VEC speaking '
        WRITE(6,*) '********************'
        WRITE(6,*) ' Output vector '
        LU = 0
        CALL WRT_CC_VEC2(JAC_VEC,LU,CCTYPE)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'JAC_VE')
*
      RETURN
      END 
      SUBROUTINE JAC_T_VEC2(L_OR_R,IADD_RHS,
     &                      IWRMOD,IREWRL,IREWSIG,
     &                      CC_AMP,CC_VEC,VEC1,VEC2,
     &                      NAMP,NAMP0,
     &                      ECC,XLAMPNRM,XRESNRM,
     &                      LUCCAMP,LU_OMG,LU_RLAMP,LU_SIG,LU_RHS,
     &                      LU_EXRF,LUHEXRF,LULAMST)
*
* Jeppe + Jesper, Dec. 1 1999
*. Symmetry added to RHS, Summer of 2000
*
* L_OR_R Flags whether left- or right- transformation is carried out
*
* Right Transformation
*     <hf!tau+_{mu} exp -T [H, sum(nu) tau_nu V_nu] exp T |hf>
* Left transformation
*
*   ECC       :: will be incremented with <L|Omg>
*   XLAMPNRM  :: norm L
*   XRESNRM   :: norm (L A - RHS)
*
*
* Symmetry of tau-operator is given by ITEX_SM in /CC_EXC/   
*
*
*. Input
*
*  LU_RLAMP   : CC amplitudes defining Linear transformation
*  LUCCAMP    : CC amplitudes of current CC state
*  
*  LU_EXRF    : e^T|ref>
*  LUHEXRF    : H e^T|ref>
*  (unit==0 implies recalculation)
*
*  IWRMOD : use vec_to_disc to load R/L and store sigma vectors (0)
*           or use frmdsc,todsc
*. Output
*
* JACVEC : Jacobiant times vector 
* 
*  LULAMST    : the lambda state (note: HF contrib is still missing!)

c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc' 
      INCLUDE 'clunit.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cstate.inc'
*. Scratch
      DIMENSION VEC1(*), VEC2(*), CC_AMP(*), CC_VEC(*)
*
      CHARACTER*6 CCTYPE
      REAL(8), EXTERNAL :: INPROD, INPRDD
*
      CALL ATIM(CPU0,WALL0)

      NTEST = 5

      IF (NTEST.GE.10) THEN
        WRITE(6,*) 'JAC_T_VEC (MOD. VERSION) reports:'
        WRITE(6,*) ' L_OR_R : ', L_OR_R
        WRITE(6,*) ' LUCCAMP, LU_OMG, LU_RLAMP, LU_SIG: ',
     &               LUCCAMP, LU_OMG, LU_RLAMP, LU_SIG
        CALL UNIT_INFO(LUCCAMP)
        CALL UNIT_INFO(LU_OMG)
        CALL UNIT_INFO(LU_RLAMP)
        CALL UNIT_INFO(LU_SIG)
        WRITE(6,*) ' LU_EXRF, LUHEXRF, LULAMST: ',
     &               LU_EXRF, LUHEXRF, LULAMST
        CALL UNIT_INFO(IABS(LU_EXRF))
        CALL UNIT_INFO(IABS(LUHEXRF))
        CALL UNIT_INFO(IABS(LULAMST))
      END IF

C?    WRITE(6,*) ' In JAC ..... '
*
*. Works only for generalized cc
      CCTYPE(1:6) = 'GEN_CC'
*
      LBLK = -1
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-20
*
      IF(L_OR_R.EQ.2) THEN
*=====================================================================*
*
*                  ======================
*                   Right transformation
*                  ======================
*
*=====================================================================*

        IF (CCFORM(1:3).NE.'TCC'.AND.CCFORM(1:3).NE.'VCC') THEN
          WRITE(6,*) ' CCFORM = >',CCFORM(1:6),'<'
          WRITE(6,*) ' jac_t_vec2 is not prepared for this ...'
          STOP 'jac_t_vec2'
        END IF
*
* load amplitudes
        WRITE(6,*) ' Right Transformation '

        IF (CCFORM(1:3).EQ.'TCC') THEN

        call vec_from_disc(cc_amp,namp0,1,lblk,luccamp)

*     <hf!tau+_{mu} exp -T [H, sum(nu) tau_nu V_nu] exp T |hf>
*.  Exp(t) !ref> on LUHC
        ICSM = IREFSM
        ISSM = IREFSM
*
        ICSPC = IETSPC
        ISSPC = IETSPC
*
        IF (LU_EXRF.GT.0) THEN
          LUEXPTREF = LU_EXRF
        ELSE
          LUEXPTREF = -LU_EXRF
c          CALL ATIM(CPUX,WALLX)
          CALL EXPT_REF2(LUC,LUEXPTREF,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &         CC_AMP,DUM,VEC1,VEC2,NAMP0,CCTYPE,0)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'RJT: time in EXPT_REF (1)',
c     &         cpu-cpux,wall-wallx)
        END IF

*
*. Term 1 : H Sum(nu) tau_nu V_nu on LUSC34 
*
*. a : tau_nu V_nu exp (T) |hf> on LUSC35
C       SIG_GCC(C,HC,LUC,LUEXPTREF,T)
c load R
        IF (IWRMOD.EQ.0) THEN
          CALL VEC_FROM_DISC(CC_AMP,NAMP,IREWRL,LBLK,LU_RLAMP)
        ELSE
          CALL FRMDSC(CC_AMP,NAMP,LBLK,LU_RLAMP,IMZERO,IAMPACK)
        END IF
        IF (NTEST.GE.5) THEN
          XNRM = INPROD(CC_AMP,CC_AMP,NAMP)
          WRITE(6,*) 'NORM OF INPUT VECTOR: ',XNRM
        END IF
        ICSM = IREFSM
        ISSM = MULTD2H(ITEX_SM,ICSM)
c        CALL ATIM(CPUX,WALLX)
        CALL SIG_GCC(VEC1,VEC2,LUEXPTREF,LUSC35,CC_AMP)
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'RJT: time in SIG_GCC (1)',
c     &        cpu-cpux,wall-wallx)
*. b : H *  tau_nu V_nu exp (T) |hf>
        ICC_EXC = 0
        ICSM = ISSM
*
        ICSPC = IETSPC
        ISSPC = ITSPC
c        CALL ATIM(CPUX,WALLX)
        CALL MV7(VEC1,VEC2,LUSC35,LUSC34,0,0)
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'RJT: time in MV7 (1)',
c     &        cpu-cpux,wall-wallx)
C       WRITE(6,*) ' H V_op exp(T) !HF > '
C       CALL WRTVCD(VEC1,LUSC34,1,-1)
*
*. Term 2 : Sum(nu) tau_nu V_nu H exp(T) |hf> on LUSC36 
*
        ICC_EXC = 0
        ICSM = IREFSM
        ISSM = IREFSM
        ICSPC = IETSPC
        ISSPC = ITSPC
        IF (LUHEXRF.GT.0) THEN
          LUHEXPTREF = LUHEXRF
        ELSE
          ! LUSC51 -> LUHEXPTREF
          LUHEXPTREF = -LUHEXRF
c          CALL ATIM(CPUX,WALLX)
          CALL MV7(VEC1,VEC2,LUEXPTREF,LUHEXPTREF,0,0)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'RJT: time in MV7 (2)',
c     &        cpu-cpux,wall-wallx)
        END IF
        ICC_EXC = 1
        ISSM = MULTD2H(ITEX_SM,ICSM)
*
c CC_AMP still contains the R vector
        ICSPC = ITSPC
        ISSPC = ITSPC
c        CALL ATIM(CPUX,WALLX)
        CALL SIG_GCC(VEC1,VEC2,LUHEXPTREF,LUSC36,CC_AMP)
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'RJT: time in SIG_GCC (2)',
c     &        cpu-cpux,wall-wallx)
C       WRITE(6,*) '  V_op H exp(T) !HF > '
C       CALL WRTVCD(VEC1,LUSC36,1,-1)
*.     (H Sum(nu) tau_nu V_nu -  sum(nu) tau_nu V_nu exp (T) |hf> on LUSC35
        FAC1 = 1.0D0
        FAC2 =-1.0D0
C           VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUSC34,LUSC36,LUSC35,1,-1)
C?    CALL WRTVCD(WORK(KVEC1),LUSC1,1,-1)
C?      WRITE(6,*) ' Input vector to exp(-T) '
C?      CALL WRTVCD(VEC1,LUSC35,1,-1) 
* Exp(-T) [H,sum_nu tau_nu V_nu] Exp(T) |hf>  
c reload the T vector
        call vec_from_disc(cc_amp,namp0,1,lblk,luccamp)
        ONEM = -1.0D0
        CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
        ICSM = ISSM
        ICSPC = ITSPC
        ISSPC = ITSPC
c        CALL ATIM(CPUX,WALLX)
        CALL EXPT_REF2(LUSC35,LUHC,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &              CC_AMP,DUM,VEC1,VEC2,NAMP0,CCTYPE,0) 
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'RJT: time in EXPT_REF (2)',
c     &        cpu-cpux,wall-wallx)
c        CALL SCALVE(CC_AMP,ONEM,N_CC_AMP)
*<hf|tau+_mu  Exp(-T) [H,sum_nu tau_nu V_nu] Exp(T) |hf> =
*<LUC| tau+mu |LUHC> = <LUHC| tau_mu |LU> as density matrix
C            DEN_GCC(C,HC,LUC,LUHC,T)
        ICSM = IREFSM
        ICSPC = IETSPC
        ISSPC = ITSPC
c        CALL ATIM(CPUX,WALLX)
        CALL DEN_GCC(VEC1,VEC2,LUC,LUHC,CC_AMP)
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'RJT: time in DEN_GCC',
c     &        cpu-cpux,wall-wallx)

*---------------------------------------------------------------------*
* right-trafo for VCC:
*---------------------------------------------------------------------*
        ELSE IF (CCFORM(1:3).EQ.'VCC') THEN

          ICSM = IREFSM
          ISSM = IREFSM
          ICSPC = IETSPC
          ISSPC = IETSPC
          
          IF (LU_EXRF.LT.0.OR.LUHEXRF.LT.0)
     &         CALL VEC_FROM_DISC(CC_AMP,NAMP0,1,LBLK,LUCCAMP)

c unless passed, make: exp(T)|Ref> --> LUEXPTREF
          IF (LU_EXRF.GT.0) THEN
            LUEXPTREF = LU_EXRF
          ELSE
            LUEXPTREF = -LU_EXRF
c            CALL ATIM(CPUX,WALLX)
            CALL EXPT_REF2(LUC,LUEXPTREF,LUSC1,LUSC2,LUSC3,
     &                                                 XCONV,MX_TERM,
     &           CC_AMP,DUM,VEC1,VEC2,NAMP0,CCTYPE,0)
c            CALL ATIM(CPU,WALL)
c            CALL PRTIM(6,'RJT: time in EXPT_REF (1)',
c     &           cpu-cpux,wall-wallx)
          END IF
          
c unless passed, make: H exp(T)|Ref> --> LUHEXPTREF
          IF (LUHEXRF.GT.0) THEN
            LUHEXPTREF = LUHEXRF
          ELSE
          ! LUSC51 -> LUHEXPTREF
            LUHEXPTREF = -LUHEXRF
c            CALL ATIM(CPUX,WALLX) 
            ICC_EXC=0
            CALL MV7(VEC1,VEC2,LUEXPTREF,LUHEXPTREF,0,0)
            ICC_EXC=1
c            CALL ATIM(CPU,WALL)
c            CALL PRTIM(6,'RJT: time in MV7 (1)',
c     &           cpu-cpux,wall-wallx)
          END IF

c get exptval <R|exp(T^+)Hexp(T)|Ref>
          EXPH = INPRDD(VEC1,VEC2,LUEXPTREF,LUHEXPTREF,1,LBLK)

c get overlap <R|exp(T^+)exp(T)|Ref>
          OVL = INPRDD(VEC1,VEC2,LUEXPTREF,LUEXPTREF,1,LBLK)
          ENERG = EXPH/OVL

c make:  (H - E) exp(T)|Ref> -> LUSC1
          CALL VECSMD(VEC1,VEC2,1d0,-ENERG,
     &         LUHEXPTREF,LUEXPTREF,LUSC1,1,LBLK)

c load R
          IF (IWRMOD.EQ.0) THEN
            CALL VEC_FROM_DISC(CC_AMP,NAMP,IREWRL,LBLK,LU_RLAMP)
          ELSE
            CALL FRMDSC(CC_AMP,NAMP,LBLK,LU_RLAMP,IMZERO,IAMPACK)
          END IF
          IF (NTEST.GE.5) THEN
            XNRM = INPROD(CC_AMP,CC_AMP,NAMP)
            WRITE(6,*) 'NORM OF INPUT VECTOR: ',XNRM
          END IF

          ICSM = IREFSM
          ISSM = MULTD2H(ITEX_SM,IREFSM)
c make:  R exp(T)|Ref> --> LUSC34
c          CALL ATIM(CPUX,WALLX)
          CALL SIG_GCC(VEC1,VEC2,LUEXPTREF,LUSC34,CC_AMP)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'RJT: time in SIG_GCC (1)',
c     &         cpu-cpux,wall-wallx)

c make:  H R exp(T)|Ref> --> LUSC35
          ICSM = MULTD2H(ITEX_SM,IREFSM)
          ISSM = ICSM
c          CALL ATIM(CPUX,WALLX)
          ICC_EXC=0
          CALL MV7(VEC1,VEC2,LUSC34,LUSC35,0,0)
          ICC_EXC=1
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'RJT: time in MV7 (1)',
c     &         cpu-cpux,wall-wallx)

c make:  (H - E) R exp(T)|Ref> --> LUSC2
          CALL VECSMD(VEC1,VEC2,1d0,-ENERG,LUSC35,LUSC34,LUSC2,1,LBLK)

c get density <R|exp(T^+)R^+ tau_mu^+ (H - E) exp(T)|Ref>
          ICSM=IREFSM
          ISSM=MULTD2H(ITEX_SM,IREFSM)
c          CALL ATIM(CPUX,WALLX)
c          CALL DEN_GCC(VEC1,VEC2,LUSC1,LUSC34,CC_AMP)
          CALL DEN_GCC(VEC1,VEC2,LUSC34,LUSC1,CC_AMP)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'RJT: time in DEN_GCC (1)',
c     &         cpu-cpux,wall-wallx)

c get density <R|exp(T^+) tau_mu^+ (H - E) R exp(T)|Ref>
          ICSM=MULTD2H(ITEX_SM,IREFSM)
          ISSM=IREFSM
c          CALL ATIM(CPUX,WALLX)
c          CALL DEN_GCC(VEC1,VEC2,LUSC2,LUSC35,CC_VEC)
          CALL DEN_GCC(VEC1,VEC2,LUSC35,LUSC2,CC_VEC)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'RJT: time in DEN_GCC (2)',
c     &         cpu-cpux,wall-wallx)

          CALL VECSUM(CC_AMP,CC_AMP,CC_VEC,1D0,-1D0,NAMP)

c gradient passed and non-zero?
          IF (LU_OMG.GT.0.AND.NAMP0.EQ.NAMP) THEN
            CALL VEC_FROM_DISC(CC_VEC,NAMP0,1,LBLK,LU_OMG)
            XNRM = SQRT(INPROD(CC_VEC,CC_VEC,NAMP0))

c   get overlap <R|exp(T^+) R exp(T)|Ref>
            OVLR = INPRDD(VEC1,VEC2,LUEXPTREF,LUSC34,1,LBLK)
c   subtract g_mu times overlap
            CALL VECSUM(CC_AMP,CC_AMP,CC_VEC,1d0,-OVLR,NAMP)
          END IF

c normalize gradient by overlap
          CALL SCALVE(CC_AMP,1d0/OVL,NAMP)

        END IF

c save the sigma vector
        IF (IWRMOD.EQ.0) THEN
          CALL VEC_TO_DISC(CC_AMP,NAMP,IREWSIG,LBLK,LU_SIG)
        ELSE
          CALL TODSC(CC_AMP,NAMP,LBLK,LU_SIG)
        END IF

      ELSE
*=====================================================================*
*
*                  =====================
*                   Left transformation
*                  =====================
*
*=====================================================================*
*
* Term 1 : 
* sum(mu)<hf|(tau_mu+ V_mu)exp -T H  tau_nu exp T|hf>
*
*  Calculate as <1| tau_nu exp T |hf>
*  with |1> = H (exp -T)+ (sum(mu) V_mu tau_mu)|hf>
*
*. A: (sum(mu) V_mu tau_mu)|hf> on LUSC34
*
        WRITE(6,*) ' Left Transformation '
        WRITE(6,*) '  CC-FORM: >',CCFORM(1:6),'<'

c load L vector
        IF (CCFORM(1:3).EQ.'TCC') THEN
          IF (IWRMOD.EQ.0) THEN
            CALL VEC_FROM_DISC(CC_AMP,NAMP,IREWRL,LBLK,LU_RLAMP)
          ELSE
            CALL FRMDSC(CC_AMP,NAMP,LBLK,LU_RLAMP,IMZERO,IAMPACK)
          END IF
          IF (NTEST.GE.5) THEN
            XNRM = INPROD(CC_AMP,CC_AMP,NAMP)
            WRITE(6,*) 'NORM OF INPUT VECTOR: ',XNRM
          END IF
          
          IF (LU_OMG.GT.0) THEN
            CALL VEC_FROM_DISC(CC_VEC,NAMP,1,LBLK,LU_OMG)
          
            ECC = ECC + INPROD(CC_AMP,CC_VEC,NAMP)
          END IF

          XLAMPNRM = SQRT(INPROD(CC_AMP,CC_AMP,NAMP))

          ICSM = IREFSM
          ISSM = MULTD2H(ITEX_SM,IREFSM) 
         
          ICSPC = IETSPC
          ISSPC = ITSPC
c          CALL ATIM(CPUX,WALLX)
          CALL SIG_GCC(VEC1,VEC2,LUC,LUSC34,CC_AMP)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'LJT: time in SIG_GCC',
c     &         cpu-cpux,wall-wallx)
        ELSE IF (CCFORM(1:3).EQ.'ECC') THEN
          CALL COPVCD(LULAMST,LUSC34,VEC1,1,LBLK)
        ELSE
          WRITE(6,*) 'ILLEGAL CCFORM!'
          STOP 'JAC_T_VEC2'
        END IF

*. B: exp(-T)+ A on LULAMST
c load T vector
        CALL VEC_FROM_DISC(CC_AMP,NAMP0,1,LBLK,LUCCAMP)
        ONEM = -1.0D0
        CALL SCALVE(CC_AMP,ONEM,NAMP0)
*. Conjugate amplitudes (reorder)  and operators
C             CONJ_CCAMP(CCIN,ISM,CCOUT)
        CALL CONJ_CCAMP(CC_AMP,1,CC_VEC)
c         CALL COPVEC(CC_VEC,CC_AMP,N_CC_AMP)
        CALL CONJ_T
*
c        CALL ATIM(CPUX,WALLX)
        ICSM = MULTD2H(ITEX_SM,IREFSM)
        ISSM = ICSM
        IF (CCFORM(1:3).EQ.'TCC') THEN
          ICSPC = ITSPC
          ISSPC = ITSPC
        ELSE
          ! until we know better ...
          ICSPC = IETSPC
          ISSPC = IETSPC
        END IF
        CALL EXPT_REF2(LUSC34,LULAMST,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &              CC_VEC,DUM,VEC1,VEC2,NAMP0,CCTYPE,0)
         ! restore T if we need it later in exp(T)|ref>
        IF(LU_EXRF.LE.0)
     &       CALL CONJ_CCAMP(CC_VEC,1,CC_AMP)
        CALL CONJ_T
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'LJT: time in EXPT_REF (1)',
c     &        cpu-cpux,wall-wallx)
*. C: H B on LUSC36
        ICC_EXC = 0
c        CALL ATIM(CPUX,WALLX)
        ICSM = MULTD2H(ITEX_SM,IREFSM)
        ISSM = ICSM
        IF (CCFORM(1:3).EQ.'TCC') THEN
          ICSPC=ITSPC
          ISSPC=IETSPC
        ELSE
          ICSPC=IETSPC
          ISSPC=IETSPC
        END IF
        CALL MV7(VEC1,VEC2,LULAMST,LUSC36,0,0)
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'LJT: time in MV7 (1)',
c     &        cpu-cpux,wall-wallx)
*. D: exp(T)|HF> on LUEXPTREF
        IF (LU_EXRF.GT.0) THEN
          LUEXPTREF = LU_EXRF
        ELSE
          LUEXPTREF = -LU_EXRF
          ICC_EXC = 1
          CALL SCALVE(CC_AMP,ONEM,NAMP0)
          ICSM = IREFSM
          ISSM = IREFSM
          ICSPC = IETSPC
          ISSPC = IETSPC
c          CALL ATIM(CPUX,WALLX)
          CALL EXPT_REF2(LUC,LUEXPTREF,LUSC1,LUSC2,LUSC3,XCONV,MX_TERM,
     &         CC_AMP,DUM,VEC1,VEC2,NAMP0,CCTYPE,0)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'LJT: time in EXPT_REF (2)',
c     &          cpu-cpux,wall-wallx)
        END IF
*. F: sum_mu<HF|(tau_mu+V_mu)exp(-T) H tau_nu exp(T)|HF>
C             DEN_GCC(C,HC,LUC,LUHC,T)

        ICSM = IREFSM
        ISSM = MULTD2H(ITEX_SM,IREFSM)
        ICSPC = IETSPC
        ISSPC = IETSPC 
c        CALL ATIM(CPUX,WALLX)
        CALL DEN_GCC(VEC1,VEC2,LUEXPTREF,LUSC36,CC_VEC)
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'LJT: time in DEN_GCC (1)',
c     &        cpu-cpux,wall-wallx)
*. G: H exp(T)|HF> on LUHEXPTREF
        ICC_EXC = 0
        IF (LUHEXRF.GT.0) THEN
          LUHEXPTREF = LUHEXRF
        ELSE
          LUHEXPTREF = -LUHEXRF
          
c          CALL ATIM(CPUX,WALLX)
          ICSM = IREFSM
          ISSM = IREFSM
          IF(CCFORM(1:3).EQ.'TCC') THEN
            ICSPC = IETSPC
            ISSPC = ITSPC
          ELSE
            ICSPC = IETSPC
            ISSPC = IETSPC
          END IF
          CALL MV7(VEC1,VEC2,LUEXPTREF,LUHEXPTREF,0,0)
c          CALL ATIM(CPU,WALL)
c          CALL PRTIM(6,'LJT: time in MV7 (2)',
c     &         cpu-cpux,wall-wallx)
          
        END IF
*. H: sum_mu<HF|(tau_mu+V_mu)exp(-T) tau_nu H exp(T)|HF>
*.    <LUSC35|tau_nu|LUSC36>
c        CALL ATIM(CPUX,WALLX)
        ICSM = IREFSM
        ISSM = MULTD2H(ITEX_SM,IREFSM)
        IF(CCFORM(1:3).EQ.'TCC') THEN
          ICSPC = ITSPC
          ISSPC = ITSPC
        ELSE
          ICSPC = IETSPC
          ISSPC = IETSPC
        END IF
        CALL DEN_GCC(VEC1,VEC2,LUHEXPTREF,LULAMST,CC_AMP)
c        CALL ATIM(CPU,WALL)
c        CALL PRTIM(6,'LJT: time in DEN_GCC (2)',
c     &       cpu-cpux,wall-wallx)

*. F-H
        ONE = 1.0D0
        ONEM = -1.0D0
        CALL VECSUM(CC_AMP,CC_AMP,CC_VEC,ONEM,ONE,NAMP)

c add rhs
        IF (IADD_RHS.EQ.1) THEN
          CALL VEC_FROM_DISC(CC_VEC,NAMP,1,LBLK,LU_RHS)
          CALL VECSUM(CC_AMP,CC_AMP,CC_VEC,1D0,-1D0,N_CC_AMP)

        END IF

c report the residual norm
        XRESNRM = SQRT(INPROD(CC_AMP,CC_AMP,NAMP))

c save the sigma vector
        IF (IWRMOD.EQ.0) THEN
          CALL VEC_TO_DISC(CC_AMP,NAMP,IREWSIG,LBLK,LU_SIG)
        ELSE
          CALL TODSC(CC_AMP,NAMP,LBLK,LU_SIG)
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '********************'
        WRITE(6,*) ' JAC_T_VEC speaking '
        WRITE(6,*) '********************'
        WRITE(6,*) ' Output vector '
        LU = 0
        CALL WRT_CC_VEC2(CC_AMP,LU,CCTYPE)
      END IF

      IF (NTEST.GE.5) THEN
        XNRM = INPROD(CC_AMP,CC_AMP,NAMP)
        WRITE(6,*) 'NORM OF OUTPUT VECTOR: ',XNRM
      END IF
      
      CALL ATIM(CPU,WALL)
      CALL PRTIM(6,'time in Jacobian-vector contraction',
     &     cpu-cpu0,wall-wall0)
      
      RETURN
      END 
      SUBROUTINE WRT_SPOX_TP_JEPPE(IEX_TP,NEX_TP)
*
* Print types of spin-orbital excitations
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*
      INTEGER IEX_TP(4*NGAS,NEX_TP)
*
      WRITE(6,*) 
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*) ' Information about spinorbital excitations '
      WRITE(6,*) ' ***************************************** '
      WRITE(6,*)
*
      DO JEX_TP = 1, NEX_TP  
        WRITE(6,*)
        WRITE(6,*) ' Included spinorbitalexcitation ', JEX_TP
        WRITE(6,'(A,16I4)') 
     &  ' Creation of alpha     :', 
     &  (IEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)') 
     &  ' Creation of beta      :', 
     &  (IEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)') 
     &  ' Annihilation of alpha :', 
     &  (IEX_TP(I+2*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A,16I4)') 
     &  ' Annihilation of beta  :', 
     &  (IEX_TP(I+3*NGAS,JEX_TP),I=1,NGAS)
      END DO
*
      RETURN
      END
      SUBROUTINE WRT_SPOX_TP(IEX_TP,NEX_TP)
*
* Print types of spin-orbital excitations
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*
      INTEGER IEX_TP(4*NGAS,NEX_TP)
      CHARACTER LINE*80, FMT*16
*
*
      LINE(1:26)=' +------------------------'
      LEN = 26
      DO I = 1, MIN(NGAS,12)
        LINE(LEN+1:LEN+4)='----'
        LEN = LEN+4
      END DO
      LINE(LEN+1:LEN+2) = '-+'
      LEN = LEN+2

      WRITE(FMT,'(A,I2,A)') '(A,',NGAS,'I4," |")'

      DO JEX_TP = 1, NEX_TP  
        IF (NEX_TP.GT.1)
     &   WRITE(6,'(/A,I4)') ' Included spinorbitalexcitation ', JEX_TP
        WRITE(6,'(A)') LINE(1:LEN)
        WRITE(6,FMT) 
     &  ' | Creation of alpha     :', 
     &  (IEX_TP(I+0*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,FMT) 
     &  ' | Creation of beta      :', 
     &  (IEX_TP(I+1*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,FMT) 
     &  ' | Annihilation of alpha :', 
     &  (IEX_TP(I+2*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,FMT) 
     &  ' | Annihilation of beta  :', 
     &  (IEX_TP(I+3*NGAS,JEX_TP),I=1,NGAS)
        WRITE(6,'(A)') LINE(1:LEN)
      END DO
*
      RETURN
      END
      SUBROUTINE CONJ_T_PAIRS(ICTP,IERR,IEX_TP,NEX_TP,NGAS)
*
* Set up array ICTP which tells which T operators are conjugate
*
* a little contribution from A.K., 2004
*
      INCLUDE 'implicit.inc'

      PARAMETER(NTEST = 00)

      DIMENSION ICTP(NEX_TP)
      DIMENSION IEX_TP(NGAS,4,NEX_TP)
* will not matter much, but maybe this small scratch helps a bit:
      DIMENSION ISCR_TP(NGAS,4)
      LOGICAL FOUND, IDENT

      IERR = 0
      DO  ITP = 1, NEX_TP
        ISCR_TP(1:NGAS,1:4) = IEX_TP(1:NGAS,1:4,ITP)
        FOUND = .FALSE.
        DO JTP = 1, NEX_TP
          IDENT = .TRUE.
          DO IGAS = 1, NGAS
            IF (ISCR_TP(IGAS,3).NE.IEX_TP(IGAS,1,JTP) .OR.
     &          ISCR_TP(IGAS,4).NE.IEX_TP(IGAS,2,JTP) .OR.
     &          ISCR_TP(IGAS,1).NE.IEX_TP(IGAS,3,JTP) .OR.
     &          ISCR_TP(IGAS,2).NE.IEX_TP(IGAS,4,JTP) ) THEN
              IDENT = .FALSE.
              EXIT
            END IF
          END DO
          FOUND = IDENT
          IF (FOUND) THEN
            ICTP(ITP) = JTP
            EXIT
          END IF
        END DO
        IF (.NOT.FOUND) THEN
          WRITE(6,*) 'WARNING: No conjugate Operator found for ',ITP
          IERR = IERR+1
          ICTP(ITP) = 0
        END IF

      END DO

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'CONJ_T_PAIRS: the ICTP array:'
        CALL IWRTMA(ICTP,NEX_TP,1,NEX_TP,1)
      END IF

      RETURN
      END
      SUBROUTINE SYMMET_T(ISIGN,ISM,
     &                    VEC,VECSCR,
     &                    ICTP,IEX_TP,NEX_TP,ILEN_TP,IOFF_TP,NGAS)
*
* Symmetrize T as  T = 0.5 (T +/- T^+)  the sign depending on ISIGN
*
      INCLUDE 'implicit.inc'

      PARAMETER(NTEST = 00)

      DIMENSION VEC(*), VECSCR(*)
      DIMENSION ICTP(NEX_TP), IEX_TP(NGAS,4,NEX_TP),
     &     ILEN_TP(NEX_TP), IOFF_TP(NEX_TP)
      
      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'input amplitudes:'
        CALL WRT_CC_VEC2(VEC,6,'GEN_CC')
      END IF

      IF (ISIGN.EQ.+1) THEN
        FAC=+0.5d0
      ELSE IF (ISIGN.EQ.-1) THEN
        FAC=-0.5d0
      ELSE
        WRITE(6,*) 'Illegal value for ISIGN: ',ISIGN
        STOP 'SYMMET_T'
      END IF

* get reordered T on vecscr
      CALL CONJ_CCAMP_S(IEX_TP,NEX_TP,VEC,ISM,VECSCR)
* loop over operator types and use ICTP information to get the partners
      DO ITP = 1, NEX_TP
        JTP = ICTP(ITP)
        IF (JTP.LE.0.OR.JTP.GT.NEX_TP) THEN
          WRITE(6,*) 'Erroneous ICTP array in SYMMET_T!'
          STOP 'SYMMET_T'
        END IF
        IOFFI = IOFF_TP(ITP)
        IOFFJ = IOFF_TP(JTP)
        ILENI = ILEN_TP(ITP)
        ILENJ = ILEN_TP(JTP)

        IF (ILENI.NE.ILENJ) THEN
          WRITE(6,*) 'Lengthes not matching in SYMMET_T: ',ILENI, ILENJ
          STOP 'SYMMET_T'
        END IF

        CALL VECSUM(VEC(IOFFI),VEC(IOFFI),VECSCR(IOFFJ),0.5D0,FAC,ILENI)

      END DO

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'output amplitudes:'
        CALL WRT_CC_VEC2(VEC,6,'GEN_CC')
      END IF

      RETURN
      END
      SUBROUTINE CHKSYM_T(ISIGN,ISM,
     &                    VEC,VECSCR,
     &                    ICTP,IEX_TP,NEX_TP,ILEN_TP,IOFF_TP,NGAS)
*
* Check if T has the desired symmetry  T = +/- T^+  the sign depending on ISIGN
*
      INCLUDE 'implicit.inc'

      PARAMETER(NTEST = 00)

      DIMENSION VEC(*), VECSCR(*)
      DIMENSION ICTP(NEX_TP), IEX_TP(NGAS,4,NEX_TP),
     &     ILEN_TP(NEX_TP), IOFF_TP(NEX_TP)
      
      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'input amplitudes:'
        CALL WRT_CC_VEC2(VEC,6,'GEN_CC')
      END IF

      THRSH = 100D0*EPSILON(1.0d0)
      IF (.NOT.(ISIGN.EQ.+1.OR.ISIGN.EQ.-1)) THEN
        WRITE(6,*) 'Illegal value for ISIGN: ',ISIGN
        STOP 'SYMMET_T'
      END IF

      FAC = DBLE(ISIGN)

* get reordered T on vecscr
      CALL CONJ_CCAMP_S(IEX_TP,NEX_TP,VEC,ISM,VECSCR)
* loop over operator types and use ICTP information to get the partners
      DO ITP = 1, NEX_TP
        JTP = ICTP(ITP)
        IF (JTP.LE.0.OR.JTP.GT.NEX_TP) THEN
          WRITE(6,*) 'Erroneous ICTP array in SYMMET_T!'
          STOP 'SYMMET_T'
        END IF
        IOFFI = IOFF_TP(ITP)
        IOFFJ = IOFF_TP(JTP)
        ILENI = ILEN_TP(ITP)
        ILENJ = ILEN_TP(JTP)

        IF (ILENI.NE.ILENJ) THEN
          WRITE(6,*) 'Lengthes not matching in SYMMET_T: ',ILENI, ILENJ
          STOP 'SYMMET_T'
        END IF

        DO II = 1, ILENI
          DIFF = VEC(IOFFI-1+II)-FAC*VECSCR(IOFFJ-1+II)
          IF (ABS(DIFF).GT.THRSH)
     &         WRITE(6,'(X,A,E12.4,A,2I4,A,I8)')
     &         'WARNING: violation of symmetry by ',DIFF,
     &         ' op.-pair',ITP,JTP,'ampl. ',II
        END DO

      END DO

      IF (NTEST.GE.100) THEN
        WRITE(6,*) 'output amplitudes:'
        CALL WRT_CC_VEC2(VEC,6,'GEN_CC')
      END IF

      RETURN
      END
      SUBROUTINE CONJ_T
*
* Conjugate the types of the T-operators
*
* Jeppe Olsen, early 2000
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'ctcc.inc'
C     COMMON/CTCC/KLSOBEX,NSPOBEX_TP,KLLSOBEX,KLIBSOBEX,LEN_T_VEC,
C    &             MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK
      INCLUDE 'cgas.inc'
*
      CALL CONJ_T_S(WORK(KLSOBEX),NSPOBEX_TP,NGAS)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Set of conjugated spinorbital excitations'
        CALL WRT_SPOX_TP(WORK(KLSOBEX),NSPOBEX_TP)
C            WRT_SPOX_TP(IEX_TP,NEX_TP)
      END IF
*
      RETURN
      END
      SUBROUTINE CONJ_T_S(IEX_TP,NEX_TP,NGAS)
*
* Conjugate the set of T-operators
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input and output
      INTEGER IEX_TP(4*NGAS,NEX_TP)
*. Local scratch 
      INTEGER ISCR(MXPNGAS)
*
      DO  JTP = 1, NEX_TP
*. Order of excit is crea-alpha, crea-beta, anni-alpha, anni-beta
       CALL ICOPVE(IEX_TP(1+0*NGAS,JTP),ISCR,NGAS)
       CALL ICOPVE(IEX_TP(1+2*NGAS,JTP),IEX_TP(1+0*NGAS,JTP),NGAS)
       CALL ICOPVE(ISCR,IEX_TP(1+2*NGAS,JTP),NGAS)
*
       CALL ICOPVE(IEX_TP(1+1*NGAS,JTP),ISCR,NGAS)
       CALL ICOPVE(IEX_TP(1+3*NGAS,JTP),IEX_TP(1+1*NGAS,JTP),NGAS)
       CALL ICOPVE(ISCR,IEX_TP(1+3*NGAS,JTP),NGAS)
      END DO
*
      RETURN
      END
C     SUBROUTINE REF_TT_BLK(TBLKI,TBLKI,ICA,IAA,IWAY)
*
* Reform T block of cc coefs between packed and unpacked form.
* It is assumed that block is diagonal. i.e. has identical alpha and
* beta form.
*
* IWAY = 1 : Packed to unpacked form
* IWAY = 2 : Unpacked to packed form
*
* Jeppe Olsen, Jan 2000
*
C     INCLUDE 'implicit.inc'
C     INCLUDE 'mxpdim.inc'
C     INCLUDE 'cgas.inc'
C     INCLUDE 'csm.inc'
C     INCLUDE 'multd2h.inc'
*. Specific input : Occupations of the creation and annihilation strings
C     INTEGER ICA(NGAS),IAA(NGAS)
*. Input block of CC coefs
C     DIMENSION TBLKI(*)
*. Output block of CC coefs
C     DIMENSION TBLKO(*)
 
*. Local scratch
C     INTEGER LCA(MXPNGAS),LAA(MXPNGAS)
C     INTEGER IB_T(8,8,8)
*. Length of creation and annihilation strings
C     CALL NST_SPGP(ICA,LCA)
C     CALL NST_SPGP(IAA,LAA)
* Offsets for symmetryblocks in expanded form
C     ISM = 1
C     CALL Z_TCC_OFF(IB_T,LCA,LCA,LAA,LAA,ISM,NSMST)
*. Loop over symmetry blocks in packed matrix
C     DO ISM_C = 1, NSMST
C       ISM_A = MULTD2H(ISM,ISM_C) 
C        DO ISM_CA = 1, NSMST
C          ISM_CB = MULTD2H(ISM_C,ISM_CA)
C          DO ISM_AA = 1, NSMST
C           ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*
C           ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
C           ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
C           IF(ISM_BETA.GT.ISM_ALPHA) GOTO 777
C           IF(ISM_ALPHA.GT.ISM_BETA) THEN
C            IDIAGSM = 0
C           ELSE
C            IDIAGSM = 1
C           END IF
*  
C           NSTR_AA = LAA(ISM_AA)
C           NSTR_AB = LAA(ISM_AB)
C           NSTR_CA = LCA(ISM_CA)
C           NSTR_CB = LCA(ISM_CB)
*
C           IOFF_EXP = IB_T(ISM_CA,ISM_CB,ISM_AA)
C           IOFF_EXP_TRP = IB_T(ISM_CB,ISM_CA,ISM_AB)
C?          WRITE(6,'(A,4I5)') ' ISM_CA, ISM_CB, ISM_AA, ISM_AB',
C?   &                           ISM_CA, ISM_CB, ISM_AA, ISM_AB
*. Loop over T elements as  matric T(I_CA, I_CB, IAA, I_AB)
C            DO I_AB = 1, NSTR_AB
C             IF(IDIAGSM.EQ.1) THEN
C               I_AA_MIN = I_AB
C             ELSE
C               I_AA_MIN = 1
C             END IF
C             DO I_AA = I_AA_MIN, NSTR_AA
C              DO I_CB = 1, NSTR_CB
C               IF(IDIAGSM.EQ.1.AND.I_AB.EQ.I_AA) THEN
C                ICA_MIN = I_CB
C               ELSE
C                ICA_MIN = 1
C               END IF
C               DO I_CA = ICA_MIN, NSTR_CA
*. Several form of expansions :
*  IDIAGSM = 1
*  Expand lower packed symmetry block to complete symmetry block
*  IDIAGSM = 0
*  Transpose symmetry block ISM_CA, ISM_CB, ISM_AA, ISM_AB to 
*                           ISM_CB, ISM_CA, ISM_AB, ISM_AA
C                IT = IT + 1 
C                IADR_EXP = (I_AB-1)*NSTR_AA*NSTR_CB*NSTR_CA 
C    &                    + (I_AA-1)*NSTR_CB*NSTR_CA
C    &                    + (I_CB-1)*NSTR_CA 
C    &                    +  I_CA + IOFF_EXP -1
C                IF(IDIAGSM.EQ.0) THEN
C                 IADR_EXP_TRP = (I_AA-1)*NSTR_AB*NSTR_CA*NSTR_CB
C    &                         + (I_AB-1)*NSTR_CA*NSTR_CB
C    &                         + (I_CA-1)*NSTR_CB +IOFF_EXP_TRP -1
C                 IF(IWAY.EQ.1) THEN
C                   TBLKO(IADR_EXP)     = TBLKI(IT)
C                   TBLKO(IADR_EXP_TRP) = TBLKI(IT)
C                 ELSE IF (IWAY.EQ.2) THEN 
C                   TBLKO(IT) = TBLKI(IADR_EXP)
C                 END IF
C                ELSE IF(IDIAGSM.EQ.1) THEN
C                 IF(IWAY.EQ.1) THEN
C                  TBLKO(IADR_EXP) = TBLKI(IT)
C                 ELSE IF(IWAY.EQ.2) THEN
C                  TBLKO(IT) = TBLKI(IADR_EXP)
C                 END IF
C                END IF
C                 
C                 
C     
*
      SUBROUTINE DIAG_EXC_CC(NCA,NCB,NAA,NAB,NGAS,IDIAG)
*
* Check of alpha and betaexcitation parts are identical
*
*. Output
* IDIAG = 0 : Not diagonal in alpha, beta
* IDIAG = 1 : Diag in alpha, beta
*
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NCA(NGAS),NCB(NGAS),NAA(NGAS),NAB(NGAS)
*
      IDIAG = 1 
      DO IGAS = 1, NGAS
        IF(NCA(IGAS).NE.NCB(IGAS)) IDIAG = 0
        IF(NAA(IGAS).NE.NAB(IGAS)) IDIAG = 0
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCA, NCB, NAA, NAB = '
        CALL IWRTMA(NCA,1,NGAS,1,NGAS)
        CALL IWRTMA(NCB,1,NGAS,1,NGAS)
        CALL IWRTMA(NAA,1,NGAS,1,NGAS)
        CALL IWRTMA(NAB,1,NGAS,1,NGAS)
        WRITE(6,*) ' IDIAG = ', IDIAG
      END IF
*
      RETURN 
      END
      SUBROUTINE WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
*
*. Write TCC block containing all symmetries, total sym of T is ITCC_SM
*
* Jeppe Olsen, summer of 99
*
       INCLUDE 'implicit.inc'
       INCLUDE 'multd2h.inc'      
*. Input
       DIMENSION TCC(*)
       INTEGER NCA(*), NCB(*), NAA(*), NAB(*)
*
       WRITE(6,*) ' Block of Coupled cluster vector '
       WRITE(6,*) ' ================================'
       IOFF = 1
       DO ISM_C = 1, NSMST
         ISM_A = MULTD2H(ITCC_SM,ISM_C) 
         DO ISM_CA = 1, NSMST
           ISM_CB = MULTD2H(ISM_C,ISM_CA)
           DO ISM_AA = 1, NSMST
             ISM_AB =  MULTD2H(ISM_A,ISM_AA)
C?           WRITE(6,*) ' ISM_AB = ', ISM_AB
             LCA = NCA(ISM_CA)
             LCB = NCB(ISM_CB)
             LAA = NAA(ISM_AA)
             LAB = NAB(ISM_AB)
             LENGTH = LCA*LCB*LAA*LAB
*
             IF(LENGTH.NE.0) THEN
               WRITE(6,'(A,4I4)') ' Sym of CA, CB, AA, AB ',
     &         ISM_CA, ISM_CB, ISM_AA, ISM_AB  
               CALL WRTMAT(TCC(IOFF),1,LENGTH,1,LENGTH)
             END IF
             IOFF = IOFF + LENGTH
           END DO
         END DO
       END DO
*
       RETURN
       END
       SUBROUTINE LUCIA_GENCC(ISM,ISPC,ICISPC,IPRNT,II_RESTRT_CC,
     &                        I_TRANS_WF,II_RES_EXC,II_DO_CMPCCI,
     &                        LGRAD,
     &                        E_FINAL,ERROR_NORM_FINAL,CONV_F)
*
* Master routine for general coupled cluster calc with LUCIA
*
* Jeppe Olsen, August 99
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPRDD
      LOGICAL CONV_F
      INCLUDE 'strinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'csmprd.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'multd2h.inc'
*
      COMMON/CMXCJ/MXCJ
      CHARACTER*6 CCTYPE
      CHARACTER FIR2DEN*8, CCLABEL*80
      DIMENSION IOCCUN(100)
      DIMENSION NOOEXCL(2)
      LOGICAL LGRAD   ! input flag: dump info for subsequent 
                      ! geom. gradient evaluation
      LOGICAL LEXIST,LRSPCNT,LSVCMO
* some small scratch arrays:
      INTEGER NFRZ(8), NACT(8), NDEL(8)
      REAL(8) INPROD
*
      WRITE(6,*)
      WRITE(6,*) ' ========================================'
      WRITE(6,*) ' General Coupled Cluster section entered '
      WRITE(6,*) ' ========================================'
      WRITE(6,*)
      WRITE(6,*) ' Space defining T-operators ', ISPC
      IF(II_RESTRT_CC.EQ.1) THEN
        WRITE(6,*) ' Restarted calculation '
      ELSE 
        WRITE(6,*) ' No restart  '
      END IF
      WRITE(6,*) ' I_TRANS_WF  = ',  I_TRANS_WF
      WRITE(6,*) ' II_DO_CMPCCI = ', II_DO_CMPCCI
*
      NTEST = 5 
      IPRNT = MAX(IPRNT,NTEST)
*
      I12 = 2
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GENCC ')
*
      ITSPC  = ISPC
      IETSPC = ICISPC
*. The reference space is the first space 
      IREFSPC = 1
*
* Find the type of reference state 
* ================================
*
*. Divide orbitals into HOLE ( Double occupied in reference space, 
*  PARTICLES( unoccupied in reference) and Valence orbitals
*  (not complete occupied or completely unoccupied)
*
*   IREFTP = 1 => CLOSED Shell HARTRE-FOCK, 
*   IREFTP = 2 => High spin open shell single det state
*   IREFTP = 3 => Cas state or more general multireference state
      CALL CC_AC_SPACES(IREFSPC,IREFTYP)
      WRITE(6,*) ' IREFTYP after call to CC_AC_SPACES ', IREFTYP

      nspin = 1
      if (ireftyp.eq.2) nspin = 2
*. Number of active orbital spaces
      NACT_SPC = 0
      IACT_SPC = 0
*. Number of hole spaces
      NHOL_SPC = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) THEN
          NACT_SPC = NACT_SPC + 1
          IACT_SPC = IGAS
        END IF
        IF(IHPVGAS(IGAS).EQ.1) THEN
          NHOL_SPC = NHOL_SPC + 1
        END IF
      END DO
      IF(IPRNT.GE.5) 
     &WRITE(6,*) ' Number of hole spaces = ', NHOL_SPC
*
      IF(NACT_SPC.GT.1) THEN
        WRITE(6,*) ' GENCC  in problems '
        WRITE(6,*) ' More than one active orbital spaces '
        WRITE(6,*) ' NACT_SPC = ',  NACT_SPC 
        STOP ' GENCC :  More than one active orbital spaces '
      END IF
      IF(IPRNT.GE.5) 
     &WRITE(6,*) ' GENCC : IACT_SPC,NACT_SPC',IACT_SPC,NACT_SPC
*. Info on active-active excitation types
      CALL ACAC_EXC_TYP(IAAEXC_TYP,MX_AAEXC,IPRCC)
*. pt active rotations are eliminated for CAS
c      IF(IAAEXC_TYP.EQ.3.AND.NOAAEX.EQ.1) THEN
c        MX_AAEXC = 0
c        DO I = 1, 20
c         WRITE(6,*) ' Active-active rotations inactivated '
c        END DO
c      END IF
      IF(IPRNT.GE.5) 
     &WRITE(6,*) ' GENCC : MX_AAEXC ', MX_AAEXC

      if ((i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1)
     &     .and.isimtrh.eq.1) then
        do ii = 1, 10
          write(6,*) ' WARNING: ISIMTRH deactivated'
        end do
        isimtrh = 0
      end if

* if we tried to save memory up to now:
*   allocate and read the two-electron integrals here
      if (isvmem.eq.1.and.incore.eq.1) then
        call memman(kint2,nint2,'ADDL  ',2,'INT2 L') 
        if (isimtrh.eq.1.and.ireftyp.ne.2) then
          len =  nint2_no_ccsym
          call memman(kint2_simtrh,len,      'ADDL  ',2,'SIMTR2')
          call memman(kpint2_simtrh,nsmob**3,'ADDL  ',2,'PSMTR2')
        else if (isimtrh.eq.1.and.ireftyp.eq.2) then
          len =  nint2_no_ccsym
          call memman(kint2_simtrh_aa,len,      'ADDL  ',2,'SMTHAA')
          call memman(kint2_simtrh_bb,len,      'ADDL  ',2,'SMTHBB')
          len =  nint2_no_ccsym_no12sym
          call memman(kint2_simtrh_ab,len,      'ADDL  ',2,'SMTHAB')
          call memman(kpint2_simtrh,nsmob**3,'ADDL  ',2,'PSMTXX')
          call memman(kpint2_simtrh_ab,nsmob**3,'ADDL  ',2,'PSMTAB')
        end if

* the simplest possible hack:
        isvmem=0
        call intim(iprorb)
        isvmem=1

      end if

* get some additional memory for integrals (well, well, well, ...)
      if (ireftyp.eq.2.and.(i_obcc.eq.1.or.i_bcc.eq.1.or.i_oocc.eq.1)
     &     .and.incore.eq.1) then
        call memman(kint1b,nint1,'ADDL  ',2,'IN1B L') 
        call memman(kint2bb,nint2,'ADDL  ',2,'I2BB L') 
        call memman(kint2ab,nint2_no12sym,'ADDL  ',2,'I2AB L')
        call memman(kpint2ab,nsmob**3,'ADDL  ',2,'PSMTR2')
        
        i12loc = 1
        i34loc = 1
        i1234loc = 0
        call pnt4dm(nsmob,nsmsx,mxpobs,ntoobs,ntoobs,ntoobs,ntoobs,
     &         itsdx,adsxa,sxdxsx,i12loc,i34loc,i1234loc,
     &         work(kpint2ab),work(klsm2),adasx)

      end if

      ! get a new rho1 array
      if (nspin.eq.2) then
        krho1old = krho1
        lrho1 = nspin*ntoob**2
        call memman(krho1,lrho1,'ADDS  ',2,'RHO1 L')
      end if

      idensi_merk = idensi
      if (((i_obcc.eq.1.or.i_bcc.eq.1.or.i_oocc.eq.1)
     &     .and.idensi.lt.2).or.
     &     (idensi.eq.2.and.isvmem.eq.1) .or. lgrad ) then
        idensi = 2
        lrho2 = nspin*ntoob**2*(ntoob**2+1)/2 + (nspin-1)*ntoob**4
        call memman(krho2,lrho2,'ADDS  ',2,'RHO2 L')
      end if


      if (ireftyp.eq.2.and.
     &    (isimtrh.eq.1.or.i_obcc.eq.1.or.i_bcc.eq.1.or.i_oocc.eq.1))
     &     then
* use unrestricted orbitals
        i_unrorb = 1
      end if
*
* new preliminary position for saving integrals:
        ! save orig. one-electron integrals
      if (i_bcc.eq.1.or.i_oocc.eq.1.or.i_obcc.eq.1) then
        lu1into  = iopen_nus('LUC_1INT')
        lu2into  = iopen_nus('LUC_2INT')
        lblk = -1
        call vec_to_disc(work(kint1o),nint1,1,lblk,lu1into)

        ! save orig. two-electron integrals
        call vec_to_disc(work(kint2),nint2,1,lblk,lu2into)
      end if

*
* 1 : Find the set of T amplitudes. 
*     These are defined as the set of excitation operators needed 
*     to obtain the CI space ISPC by exciting determinants in 
*     the reference space
*
*. Orbital excitation operators
*. Number of types of orbital excitation operators and length to store 
*  these types
*
*. Number of orbital excitations equals the number of occupation 
* classes in ISPC ( minus the reference space, which should be 
* included )
*. General CC so
      CCTYPE(1:6) = 'GEN_CC'

      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL = NAEL + NBEL
*
      ICSPC = ICISPC
      ISSPC = ICISPC
*. Routines use complete TT blocks so
C     ISIMSYM = 1
* 
* ========================
* info for reference space
* ========================
*
*. Make sure that there is just a single occupation space
      CALL OCCLSE(1,NOCCLS_REF,IOCCLS,NEL,IREFSPC,0,0,NOBPT)
      IF(NOCCLS_REF.NE.1) THEN
        WRITE(6,*) ' Problem in general CC '
        WRITE(6,*) 
     &  ' Reference space is not a single occupation space'
        STOP 
     &  ' Reference space is not a single occupation space'
      END IF
*. and the reference occupation space 
      CALL MEMMAN(KLOCCLS_REF,NGAS,'ADDL  ',1,'OCC_RF')
      CALL OCCLSE(2,NOCCLS_REF,WORK(KLOCCLS_REF),NEL,IREFSPC,0,0,NOBPT)
*
* ====================================
* Info for space defining excitations
* ====================================
*
*. Number 
      CALL OCCLSE(1,NOCCLS,IOCCLS,NEL,ISPC,0,0,NOBPT)
*. And the occupation classes
      CALL MEMMAN(KLOCCLS,NOCCLS*NGAS,'ADDL  ',1,'OCCLS ')
      CALL OCCLSE(2,NOCCLS,WORK(KLOCCLS),NEL,ISPC,0,0,NOBPT)
*. Number of occupation classes for T-operators 
      NTOCCLS = NOCCLS
*. It could be an idea to check that reference space is included 
*
* ========================
* Orbital excitation types 
* ========================
*
*. Number of excitation types 
* only for pure OCC we turn off singles; else, we need the T1-residual
      IFLAG = 1
      IDUM = 1
      CALL TP_OBEX2(NOCCLS,NEL,NGAS,WORK(IDUM),
     &             WORK(IDUM),WORK(IDUM),
     &             WORK(KLOCCLS),WORK(KLOCCLS_REF),MX_NCREA,MX_NANNI,
     &             MX_EXC_LEVEL,WORK(IDUM),MX_AAEXC,IFLAG,
     &             I_OOCC,NOBEX_TP,NOAAEX,IPRCC)
C     TP_OBEX2(NOCCLS,NEL,NGAS,IOBEX_TP,LCOBEX_TP,LAOBEX_TP,
C    &                   IOCCLS,IOCCLS_REF,MX_NCREA,MX_NANNI,
C    &                   MX_EXC_LEVEL,IEXTP_TO_OCCLS,MX_AAEXC,IFLAG,
C    &                   NOBEX_TP,IPRCC)
      IF(IPRNT.GE.5) 
     &WRITE(6,*) ' NOBEX_TP,MX_EXC_LEVEL = ', NOBEX_TP,MX_EXC_LEVEL
*. And the actual orbital excitations
*.  An orbital excition operator is defined by 
*   1 : Number of creation operators
*   2 : Number of annihilation operators 
*   3 : The actual creation and annihilation operators
*. The number of orbital excitations is increased by one to include 
*. excitations within the reference space
      NOBEX_TPE = NOBEX_TP+1
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
*
* =======================
* Spinorbital excitations
* =======================
*
*. Spin combinations of CC excitations : Currently we assume that 
*. The T-operator is a singlet, can 'easily' be changed 
C?    WRITE(6,*) ' MSCOMB_CC in LUCIA_GENCC ', MSCOMB_CC
*
*. Notice : The first time in OBEX_TO_SPOBEX we always use MSCOMB_CC = 0.
*. This may lead to the allocation of too much space for 
*. spinorbital excitations, but MSCOMB_CC = 1, requires access 
*. to WORK(KLSOBEX) which has not been defined 
*
*. Largest spin-orbital excitation level 
      IF(MXSPOX.NE.0) THEN
        MXSPOX_L = MXSPOX
      ELSE 
        MXSPOX_L = MX_EXC_LEVEL
      END IF
      WRITE(6,*) ' MXSPOX, MXSPOX_L, MX_EXC_LEVEL = ',
     &             MXSPOX, MXSPOX_L, MX_EXC_LEVEL
      IZERO = 0
      CALL OBEX_TO_SPOBEX(1,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,IDUMMY,NSPOBEX_TP,NGAS,
     &     NOBPT,0,IZERO ,IAAEXC_TYP,IACT_SPC,IPRCC,IDUMMY,
     &     MXSPOX_L,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPC)
*. Extended number of spin-orbital excitations : Include 
*. unit operator as last spinorbital excitation operator 
      NSPOBEX_TPE = NSPOBEX_TP + 1
      IF(IPRNT.GE.5) WRITE(6,*) ' NSPOBEX_TP = ', NSPOBEX_TP 
*. And the actual spinorbital excitation operators 
      CALL MEMMAN(KLSOBEX,4*NGAS*NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. Map spin-orbital exc type => orbital exc type
      CALL MEMMAN(KLSOX_TO_OX,NSPOBEX_TPE,'ADDL  ',1,'SPOBEX')
*. First SOX of given OX ( including zero operator )
      CALL MEMMAN(KIBSOX_FOR_OX,NOBEX_TP+1,'ADDL  ',1,'IBSOXF')
*. Number of SOX's for given OX
      CALL MEMMAN(KNSOX_FOR_OX,NOBEX_TP+1,'ADDL  ',1,'IBSOXF')
*. SOX for given OX
      CALL MEMMAN(KISOX_FOR_OX,NSPOBEX_TP+1,'ADDL  ',1,'IBSOXF')
      CALL OBEX_TO_SPOBEX(2,WORK(KOBEX_TP),WORK(KLCOBEX_TP),
     &     WORK(KLAOBEX_TP),NOBEX_TP,WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &     NOBPT,0,MSCOMB_CC,IAAEXC_TYP,IACT_SPC,
     &     IPRCC,WORK(KLSOX_TO_OX),MXSPOX_L,WORK(KNSOX_FOR_OX),
     &     WORK(KIBSOX_FOR_OX),WORK(KISOX_FOR_OX),NAEL,NBEL,IREFSPC)
      NSPOBEX_TPE = NSPOBEX_TP + 1
*. Add unit-operator as last spinorbital excitation
C     ISTVC3(IVEC,IOFF,IVAL,NDIM)
      IZERO = 0
      CALL ISTVC3(WORK(KLSOBEX),(NSPOBEX_TPE-1)*4*NGAS+1,IZERO,4*NGAS)
      IF(IPRNT.GE.5) THEN
        WRITE(6,*) ' Extended list of spin-orbital excitations : '
        CALL WRT_SPOX_TP(WORK(KLSOBEX),NSPOBEX_TPE) 
      END IF
      CALL ISTVC3(WORK(KLSOX_TO_OX),NSPOBEX_TPE,NOBEX_TP+1,1)
*. Mapping spinorbital excitations => occupation classes 
      CALL MEMMAN(KIBSOX_FOR_OCCLS,NOCCLS,'ADDL  ',1,'IBSXOC')
      CALL MEMMAN(KNSOX_FOR_OCCLS,NOCCLS,'ADDL  ',1,' NSXOC')
      CALL MEMMAN(KISOX_FOR_OCCLS,NSPOBEX_TPE,'ADDL  ',1,' ISXOC')
C       SPOBEX_FOR_OCCLS(
C    &           IEXTP_TO_OCCLS,NOCCLS,ISOX_TO_OX,NSOX,
C    &           NSOX_FOR_OCCLS,ISOX_FOR_OCCLS,IBSOX_FOR_OCCLS)
      CALL SPOBEX_FOR_OCCLS(WORK(KEX_TO_OC),NOCCLS,WORK(KLSOX_TO_OX),
     &     NSPOBEX_TPE,WORK(KNSOX_FOR_OCCLS),WORK(KISOX_FOR_OCCLS),
     &     WORK(KIBSOX_FOR_OCCLS))
*
*. Frozen spin-orbital excitation types
      CALL MEMMAN(KLSPOBEX_FRZ, NSPOBEX_TP+1,'ADDL  ',1,'SPOBFR')
      CALL FRZ_SPOBEX(WORK(KLSPOBEX_FRZ),WORK(KLCOBEX_TP),NSPOBEX_TP,
     &                WORK(KLSOX_TO_OX),IFRZ_CC_AR,NFRZ_CC)
      IZERO = 0
      CALL ISTVC3(WORK(KLSPOBEX_FRZ),NSPOBEX_TP+1,IZERO,1)
*. Spin-orbital excitation types related by spin-flip
      CALL MEMMAN(KLSPOBEX_AB,NSPOBEX_TP+1,'ADDL  ',1,'SPOBAB')
      CALL SPOBEXTP_PAIRS(NSPOBEX_TP+1,WORK(KLSOBEX),NGAS,
     &                    WORK(KLSPOBEX_AB))
      IF(NTEST.GE.5) WRITE(6,*) ' After SPOBEXTP_PAIRS'
C          SPOBEXTP_PAIRS(NSPOBEX_TP,ISPOBEX,NGAS,ISPOBEX_PAIRS)
C     SELECT_AB_TYPES(NSPOBEX_TP,ISPOBEX_TP,
C    &           ISPOBEX_PAIRS,NGAS)
      CALL SELECT_AB_TYPES(NSPOBEX_TP+1,WORK(KLSOBEX),
     &                     WORK(KLSPOBEX_AB),NGAS)
      IF(NTEST.GE.5) WRITE(6,*) ' After SELECT_AB_TYPES '
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
      IF(NTEST.GE.5) WRITE(6,*) ' After SPOBEX_TO_ABOBEX '
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
      IF(NTEST.GE.10) WRITE(6,*) ' After LEN_GENOP_STR_MAP '
      MAXLEN_I1 = MAX(MAXLENA,MAXLENB)
      IF(IPRNT.GE.5) WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1

*
* Max Dimension of spinorbital excitation operators
*
      CALL MEMMAN(KLLSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLIBSOBEX,NSPOBEX_TPE,'ADDL  ',1,'LSPOBX')
      CALL MEMMAN(KLSPOBEX_AC,NSPOBEX_TPE,'ADDL  ',1,'SPOBAC')
*. ALl spinorbital excitations are initially active 
      IONE = 1
      CALL ISETVC(WORK(KLSPOBEX_AC),IONE,NSPOBEX_TPE)
*
      MX_ST_TSOSO_MX = 0
      MX_ST_TSOSO_BLK_MX = 0
      MX_TBLK_MX = 0
      MX_TBLK_AS_MX = 0
      LEN_T_VEC_MX = 0
*
      DO ICCAMP_SM = 1, NSMST
*
        CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ICCAMP_SM,
     &       MX_ST_TSOSOL,MX_ST_TSOSO_BLKL,MX_TBLKL,
     &       WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VECL,
     &       MSCOMB_CC,MX_TBLK_AS,
     &       WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &       NTCONF,IPRCC)
*
        MX_ST_TSOSO_MX = MAX(MX_ST_TSOSO_MX,MX_ST_TSOSOL)
        MX_ST_TSOSO_BLK_MX = MAX(MX_ST_TSOSO_BLK_MX,MX_ST_TSOSO_BLKL)
        MX_TBLK_MX = MAX(MX_TBLK_MX,MX_TBLKL)
        MX_TBLK_AS_MX = MAX(MX_TBLK_AS_MX,MX_TBLK_AS)
        LEN_T_VEC_MX = MAX(LEN_T_VEC_MX, LEN_T_VECL)
*
      END DO
      IF(IPRNT.GE.5) WRITE(6,*) ' MX_TBLK_AS_MX = ', MX_TBLK_AS_MX
      IF(IPRNT.GE.5) WRITE(6,*) ' LEN_T_VEC_MX = ', LEN_T_VEC_MX
*
*. And dimensions for symmetry 1
      ITOP_SM = 1
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,ITOP_SM,
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MSCOMB_CC,MX_SBSTR,
     &     WORK(KISOX_FOR_OCCLS),NOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &     NTCONF,IPRCC)
      IF(NTEST.GE.100) WRITE(6,*) ' After IDIM_TCC section '
      N_CC_AMP = LEN_T_VEC
      IF(IPRNT.GE.5) 
     &WRITE(6,*) ' LUCIA_GENCC : N_CC_AMP = ', N_CC_AMP
      IF(IPRNT.GE.5) THEN
        WRITE(6,*) ' Number of amplitudes per operator type: '
        CALL IWRTMA(WORK(KLLSOBEX),NSPOBEX_TP,1,NSPOBEX_TP,1)
      END IF
*. Hard wire info for unit operator stored as last spinorbital excitation
C  ISTVC2(IVEC,IBASE,IFACT,NDIM)
      IONE = 1
      CALL ISTVC3(WORK(KLLSOBEX),NSPOBEX_TPE,IONE,1)
      N_CC_AMPP1 = N_CC_AMP + 1
      CALL ISTVC3(WORK(KLIBSOBEX),NSPOBEX_TPE,N_CC_AMPP1,1)
      IF(I_DO_NEWCCV.GE.1) THEN
*. Space for largest T(D,L) block
        CALL DIM_TDL2(NSPOBEX_TP+1,WORK(KLSOBEX),N_TDL_MAX,NSMST)
        WRITE(6,*) ' Memory for storing T(D,L) = ',N_TDL_MAX
COLD    CALL T_DL_DIM(WORK(KLSOBEX),NSPOBEX_TP,4,N_TDL_MAXOLD,LCCB)
COLD    WRITE(6,*) ' N_TDL_MAXOLD = ', N_TDL_MAXOLD
      END IF
*
* ==========================================
* Prototype information for Spinadaptation 
* ==========================================
CSPIN_ADAPT_CC_OP(
C    &           NSPOBEX_TP,ISPOBEX_TP,ISPOBEX_TO_OCCLS,
C    &           IBSPOBEX_TO_OCCLS,NSPOBEX_FOR_OCCLS,NOCCLS,S)
      IF(ISPIN_RESTRICTED.EQ.1) THEN
*. Set up S2 matrix for the various occupation classes (configurations)
*. and divide eigenvectors of S2 into P and Q space according to          
*. Szalay and  Gauss
        XMULTS = MULTS
        S = (XMULTS-1)/2
        WRITE(6,*) ' S = ', S
        CALL SPIN_ADAPT_CC_OP(NSPOBEX_TP,WORK(KLSOBEX),
     &       WORK(KISOX_FOR_OCCLS),WORK(KIBSOX_FOR_OCCLS),
     &       WORK(KNSOX_FOR_OCCLS),NOCCLS,S,WORK(KLOCCLS))
      END IF

      ! set up I_ADX array using the T-space
      CALL GASSPC2(I_IADX,IGSOCCX(1,1,ITSPC))

* ========================================
*  some specials for OCC and Brueckner CC
* ========================================

      if (i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1) then
        if (ntest.ge.5) then
          write(6,*) '============================================='
          write(6,*) ' setting up non-redundant orbital rotations:' 
          write(6,*) '============================================='
       end if

        if (i_unrorb.eq.0) then
          ! Nonredundant type-type excitations
          call memman(klttact,ngas*ngas,'ADDL  ',1,'TTACT ')
c          call nonred_tt_exc(work(klttact),1)
          call nonred_tt_exc_cc(work(klttact),itspc,1)
          ! nonredundant orbital excitations
          call memman(klooexc,ntoob*ntoob,'ADDL  ',1,'OOEXC ')
          call memman(klooexcc,2*ntoob*ntoob,'ADDL  ',1,'OOEXCC')
          call nonred_oo_exc(nooexcl(1),work(klooexc),work(klooexcc),
     &                       1,work(klttact),2)
          nooexcl(2) = 0
        else if (i_unrorb.eq.1) then
          ! Nonredundant type-type excitations
          call memman(klooexcc,2*2*ntoob*ntoob,'ADDL  ',1,'OOEXCC')
          idum = 0
          call memman(klttact,2*ngas*ngas,'ADDL  ',1,'TTACT ')
          call memman(klooexc,2*ntoob*ntoob,'ADDL  ',1,'OOEXC ')

          call nonred_os(work(klttact),work(klooexc),work(klooexcc),
     &                   nooexcl,itspc,ntoob,ngas)
        else
          write(6,*) 'not prepared for this: i_unrorb = ',i_unrorb
          stop 'lucia_gencc'
        end if
      else
        klooexcc = 1
        nooexcl(1:2) = 0
      end if

*
* =============
* Scratch space 
* =============
*
*. Scratch space for CI - behind the curtain 
       IF(I_DO_NEWCCV.EQ.0.OR.I_DO_CC_EXP.EQ.1.OR.I_TRANS_WF.EQ.1) THEN
         CALL GET_3BLKS_GCC(KVEC1,KVEC2,KVEC3,MXCJ)
C      ELSE
C        CALL GET_3BLKS(KVEC1,KVEC2,KVEC3) 
       END IF
*. Space is needed for CI in reference space 
C        IREFSPC = 1
C        ICSPC = IREFSPC
C        ISSPC = IREFSPC
C        CALL GET_3BLKS(KVEC1,KVEC2,KVEC3) 
C      END IF
*. Pointers to KVEC1 and KVEC2, transferred through GLBBAS 
       KVEC1P = KVEC1
       KVEC2P = KVEC2
       WRITE(6,*) ' KVEC3 after GET_3BLKS.. ', KVEC3
*. and two CC vectors , extra element for unexcited SD at end of vectors
       N_SD_INT = 1
       LENNY = LEN_T_VEC_MX + N_SD_INT
       if (i_obcc.eq.1.or.i_oocc.eq.1.or.i_bcc.eq.1) then
         lenny = max(lenny,nooexcl(1)+nooexcl(2))
       end if
       CALL MEMMAN(KCC1,LENNY,'ADDL  ',2,'CC1_VE')
       CALL MEMMAN(KCC2,LENNY,'ADDL  ',2,'CC2_VE')
*. And the CC diagonal 
       CALL MEMMAN(KDIA,LENNY,'ADDL  ',2,'CC_DIA')
*. Check that H_TF space is allowed
       IF(I_DO_NEWCCV.EQ.1) THEN
*. Check that no triple exctations T(I)T(J)T(K) are left out in this 
*. approach
        CALL CHECK_HTF_APR(NSPOBEX_TP,WORK(KLSOBEX),I_HTF_OKAY) 
       END IF
*
* ============================================================
* Initialize : Restart, set coefs to zero, start from CI coefs
* ============================================================
*
      IF(II_RESTRT_CC.EQ.1) THEN
*. Restart from previous CC coefs
        IFORM = 2
      ELSE   
*. No restart 
          IFORM = 1
      END IF
*
      LU_1DEN = IOPEN_NUS('CC1DEN')
      LU_2DEN = IOPEN_NUS('CC2DEN')
      IF (I_OBCC.EQ.1.OR.I_OOCC.EQ.1.OR.I_BCC.EQ.1)
     &     LUCMO = IOPEN_NUS('CCNEWCMO')
      LU_OMG  = IOPEN_NUS('CC_OMG')
      LUCCAMP = IOPEN_NUS('CCTVEC')
      LU_LAMP = IOPEN_NUS('CCLVEC')

*
* ==========================================
*. Transfer control to optimization routine
* ==========================================
*
*. New optimizer
c new amplitudes will be returned on LUCCAMP
      CALL OPTIM_CC_NEW(II_RESTRT_CC,ECC,LGRAD,
     &                 WORK(KCC1),WORK(KCC2),WORK(KDIA),
     &                 WORK(KVEC1),WORK(KVEC2),
     &                 WORK(KLTTACT),WORK(KLOOEXCC),WORK(KLOOEXC),
     &                 NOOEXCL,CCTYPE,
     &                 LUCCAMP,LU_OMG,LU_LAMP,LU_1DEN,
     &                 LU_2DEN,
     &                 LUCMO,LU1INTO,LU2INTO,
     &                 E_FINAL,ERROR_NORM_FINAL,CONV_F)
      WRITE(6,*) ' KVEC3 after OPTIM_CC2 ', KVEC3

*
* ===========================
* Calculate density matrices 
* ===========================
      IF(IDENSI.NE.0.AND.I_DO_NEWCCV.EQ.0) THEN

c        LENRHO1 = NTOOB**2
c        LENRHO2 = NTOOB**2*(NTOOB**2+1)/2
c        ICORRDEN = IRELAX
        IRDEN = 0
        IFGEN = 0
        ! for a gradient calculation we have to supply (on files):
        IF (LGRAD) IRDEN = 1  ! irreducible 2-particle density
        IF (LGRAD) IFGEN = 1  ! generalized Fock-matrix
        IREFDEN = IRELAX

c note: IRELAX and LGRAD; we currently leave it open, whether LUCIA shall
c       provide relaxed densities or the subsequent program solves the
c       corresponding equations.
        
        ntbsq = ntoob*ntoob
        lrho1 = nspin*ntbsq
        lrho2 = nspin*ntbsq*(ntbsq+1)/2 + (nspin-1)*ntbsq**2

        IF (IRDEN.EQ.1) THEN
          FIR2DEN = 'CCIR2DEN'  ! MAX. 8 CHARACTERS !!
          LU_IR2DEN = IOPEN_NUS(FIR2DEN)
        END IF
        IF (IFGEN.EQ.1) THEN
          LU_FGEN = IOPEN_NUS('CCFGEN')
        END IF

        IF (IREFDEN.EQ.1) THEN
          LU_1DENR = IOPEN_NUS('RF1DEN')
          LU_2DENR = IOPEN_NUS('RF2DEN')
        END IF

*. Unless this was orbital-optimized CC, we need to solve some more
*. equations
        IF (I_OBCC.EQ.0.AND.I_OOCC.EQ.0.AND.CCFORM(1:3).EQ.'TCC')
     &    CALL CC_DENSI(ECC,IREFDEN,
     &        WORK(KCC1),WORK(KCC2),
     &        WORK(KVEC1),WORK(KVEC2),
     &        LUCCAMP,LU_OMG,LU_LAMP,LUC,LUDIA,
     &        LU_1DEN,LU_2DEN,
     &        LU_1DENR,LU_2DENR)

* preliminary: set our own arrays for frozen, active, deleted
        CALL SET_FAD(NFRZ,NACT,NDEL)

        IF (IRELAX.NE.0.OR.LGRAD) THEN

          IF (NSPIN.NE.1) STOP 'ADAPT FOR NSPIN.NE.1'

          LBLK = -1
          CALL VEC_FROM_DISC(WORK(KRHO1),LRHO1,1,LBLK,LU_1DEN)
          CALL VEC_FROM_DISC(WORK(KRHO2),LRHO2,1,LBLK,LU_2DEN)

          IF (IRDEN.EQ.1) THEN

            ! make irreducible 2-particle density
            IMODE = 0

            CALL MKIR2DENS(IMODE,WORK(KRHO1),WORK(KRHO2),NTOOB,NSPIN)

            ! test irr. 2-p. density
            icc_exc_merk = icc_exc
            icc_exc = 0
            CALL EN_FROM_DENS(ENER,2,1)
            icc_exc = icc_exc_merk

            print *,'Energy from EN_FROM_DENS(*): ',ENER
            print *,'CC-Energy :                ',ECC
            print *,'Difference:                ',ENER-ECC
            

            CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MBUF2D')
            LEN = 0
            DO IJSM = 1,NSMOB
              LENIJ = 0
              DO JSM = 1,NSMOB
                ISM = MULTD2H(JSM,IJSM)
                LENIJ = LENIJ+NACT(ISM)*NACT(JSM)
              END DO
              LEN = MAX(LEN,LENIJ)
            END DO
            CALL MEMMAN(KLBUF2D,LEN,'ADDL  ',2,'BUFF2D')

            IMODE = 0
            CALL WR_IR2DEN_BLK(IMODE,WORK(KRHO2),WORK(KLBUF2D),
     &           LU_IR2DEN,1,
     &           NTOOB,NSMOB,NFRZ,NACT,NDEL,
     &           IBSO,IREOST)

            CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'MBUF2D')

c            CALL VEC_TO_DISC(WORK(KRHO2),LRHO2,1,LBLK,LU_IR2DEN)
            
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Irreducible part of two-body CC density: ' 
              CALL PRINT_2DENS(WORK(KRHO2),NTOOB,0,1d-1)
            END IF

            ! Well, well, we reload ....
            CALL VEC_FROM_DISC(WORK(KRHO2),LRHO2,1,LBLK,LU_2DEN)

          END IF ! IF IRDEN

          icc_exc_merk = icc_exc
          icc_exc = 0
          CALL EN_FROM_DENS(ENER,2,0)
          print *,'Energy from EN_FROM_DENS: ',ENER
          print *,'CC-Energy :               ',ECC
          print *,'Difference:               ',ENER-ECC
          NTOOB2_BLK = 0
          DO ISM = 1, NSMOB
            NTOOB2_BLK = NTOOB2_BLK + NTOOBS(ISM)*NTOOBS(ISM)
          END DO

          CALL MEMMAN(KLFOO,NTOOB2_BLK,'ADDL  ',2,'FOO   ')
cccc eventuell: fock_mat that uses f0 and lam2
          CALL FOCK_MAT(WORK(KLFOO),2)

          LRSPCNT = .FALSE.
          IF (IRELAX.EQ.1.AND.
     &        I_OOCC.EQ.0.AND.I_OBCC.EQ.0.AND.I_BCC.EQ.0) THEN
            CALL MEMMAN(KLRESP_DEN,NTOOB**2,'ADDL  ',2,'RESP_D')
            IRELAX_LOC = 1
            LRSPCNT = .TRUE.
            CALL RESPDEN_FROM_F(WORK(KLFOO),WORK(KLRESP_DEN))
            ! response-density is in blocked form
            LU_RSPCNT = IOPEN_NUS('CCRSPCNT')
            CALL VEC_TO_DISC(WORK(KLRESP_DEN),NTOOB2_BLK,
     &                                       1,LBLK,LU_RSPCNT)
            IF (LGRAD) THEN

              ! add relaxation contribution to FGEN

              CALL FGEN_RELAX(WORK(KLFOO),WORK(KFIZ),
     &             WORK(KLRESP_DEN),IHPVGAS_AB,MXPNGAS)

            END IF

          ELSE
            IRELAX_LOC = 0
            KLRESP_DEN = 1
          END IF
          icc_exc = icc_exc_merk

          IF (IFGEN.EQ.1)
     &         CALL VEC_TO_DISC(WORK(KLFOO),NTOOB2_BLK,1,LBLK,LU_FGEN)

          IF (LGRAD) THEN

            ! Reform RHO1 to symmetry-blocked form
            ! KLFOO used for blocked 1-density
            CALL REORHO1(WORK(KRHO1),WORK(KLFOO),1,1)
            
            LU_1DEN_R = IOPEN_NUS('CC1DEN_SCR')
            CALL VEC_TO_DISC(WORK(KLFOO),NTOOB2_BLK,1,LBLK,LU_1DEN_R)

            CALL CCWFLABEL(CCFORM,I_OOCC,I_BCC,I_OBCC,CCLABEL)
            LSVCMO = I_OOCC.EQ.1.OR.I_OBCC.EQ.1.OR.I_BCC.EQ.1
            ! set IRELAX flag: 
            !  1 -- relaxed densities provided (orbital relaxed methods)
            !  2 -- orbital response density provided
            IRELAX_L = 2*IRELAX
            IF (I_OOCC.EQ.1.OR.I_OBCC.EQ.1) IRELAX_L = 1
            CALL WRTLUCIFC(ECC,CCLABEL,IRELAX_L,
     &                     LUCMO,LSVCMO,
     &                     LU_FGEN,.TRUE.,
     &                     LU_1DEN_R,.TRUE.,
     &                     LU_RSPCNT,LRSPCNT,
     &                     LU_IR2DEN,FIR2DEN,.TRUE.,
     &                     NSMOB,NSPIN,NTOOBS,NAOS_ENV,
     &                     NFRZ,NACT,NDEL)

            CALL RELUNIT(LU_1DEN_R,'delete')
            IF (LRSPCNT)
     &           CALL RELUNIT(LU_RSPCNT,'delete')
c ???? here ????
            CALL RELUNIT(LU_IR2DEN,'keep')
            CALL RELUNIT(LU_FGEN,'keep')

          END IF

          ! restore unrelaxed density
          CALL VEC_FROM_DISC(WORK(KRHO1),LRHO1,1,LBLK,LU_1DEN)

        ELSE
          IRELAX_LOC = 0
          KLRESP_DEN = 1
        END IF
      END IF
*
* =========================================
*. Coupled cluster first order properties
* =========================================
*
      IF(NPROP.GT.0) THEN
        LBLK = -1
        CALL VEC_FROM_DISC(WORK(KRHO1),LRHO1,1,LBLK,LU_1DEN)
*. Calculate properties 
c        KLDUM = 1
c        CALL ONE_EL_PROP(1,0,WORK(KLDUM))
        CALL ONE_EL_PROP(1,IRELAX_LOC,WORK(KLRESP_DEN))
      END IF

      idensi = idensi_merk
*
* ======================
*. Analyze T- ampitudes
* ======================
*
      LBLK = -1
      CALL VEC_FROM_DISC(WORK(KCC1),N_CC_AMP,1,LBLK,LUCCAMP)
      WRITE(6,*) ' KVEC3 before ANA_GCC ', KVEC3
      CALL ANA_GENCC(WORK(KCC1),1)
      WRITE(6,*) ' KVEC3 after ANA_GCC ', KVEC3

*
*
* =========================
* Expectation value energy 
* =========================
*
*.  Exp(t) !ref> 
      IF(I_DO_CC_EXP.EQ.1) THEN
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-20
      I_USE_SIMTRH = 0
*
      ICSPC = IETSPC
      ISSPC = IETSPC
*
      CALL EXPT_REF(LUC,LUSC1,LUHC,LUSC2,LUSC3,XCONV,MX_TERM,
     &             WORK(KVEC1),WORK(KVEC2),CCTYPE)
*. For Frank J : Print coeffiecients of expanded  wave-function
      WRITE(6,*) ' Note : CI expansion of expanded wf printed '

*. H Exp T !ref>
      ICC_EXC = 0
C?    WRITE(6,*) ' Input to MV7 '
C?    CALL WRTVCD(WORK(KVEC1),LUSC1,1,-1)
*. Regenerate inactive Fock matrix - different definition of 
*. inactive Fock matrix may have been used in CC calculation
      CALL COPVEC(WORK(KINT1O),WORK(KINT1),NINT1)
      I_USE_SIMTRH = 0
      CALL FI(WORK(KINT1),ECCX,1)
      ECORE = ECORE_INI

      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC,0,0)
C?    WRITE(6,*) ' Output from MV7 '
C?    CALL WRTVCD(WORK(KVEC1),LUHC,1,-1)
*. E = <ref! exp(T)+ H exp(T)!ref>/<ref! exp(T)+exp(T)!ref>
      LBLK = -1
      CHC  = INPRDD(WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC  ,1,LBLK)
      CHHC = INPRDD(WORK(KVEC1),WORK(KVEC2),LUHC,LUHC  ,1,LBLK)
      CC   = INPRDD(WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC1 ,1,LBLK)
      WRITE(6,*)
      WRITE(6,*) ' Energy as expectation value : ' 
      WRITE(6,*) ' ============================================='
      WRITE(6,*)
      WRITE(6,'(5X,A,E25.12)')
     &'  <ref! exp(T)+ H exp(T)!ref> (- Ecore)= ', CHC
      WRITE(6,'(5X,A,E25.12)')
     &'  <ref! exp(T)+   exp(T)!ref> = ', CC
      WRITE(6,'(5X,A,E25.12)')
     &'                        ECORE = ', ECORE
      WRITE(6,'(5X,A,E25.12)') 
     &' Expectation value coupled cluster energy   = ', CHC/CC + ECORE
      WRITE(6,*)
      WRITE(6,*) ' (<CC|H^2|CC> - <CC|H|CC>^2/<CC!CC>)/<CC!CC>  = ',
     &              CHHC/CC - (CHC/CC)**2 
*
      WRITE(6,'(/,2(A,/),/)')
     & ' NOTE: Make sure that this result was obtained with the GEN_CC',
     & '   option (! not CC !) using the FCI space as underlying GAS!'
*
      WRITE(6,*) 'I_TRANS_WF (2) = ', I_TRANS_WF
      IF(I_TRANS_WF.EQ.1) THEN
        WRITE(6,*) ' cc wf is normalized and transferred to LUC '
        XNORM = SQRT(CC)
        FACTOR = 1.0D0/XNORM
        CALL SCLVCD(LUSC1,LUC,FACTOR,WORK(KVEC1),1,-1)
      END IF
*
      IF(IDENSI.NE.0) THEN
*. Normalize wavefunction 
        WRITE(6,*) ' cc wf is normalized and transferred to LUC '
        XNORM = SQRT(CC)
        FACTOR = 1.0D0/XNORM
        CALL SCLVCD(LUSC1,LUSC2,FACTOR,WORK(KVEC1),1,-1)
*. Calculate density matrices in standard expectation value formulation 
        LBLK = -1
        XDUM = 0.0D0
        CALL COPVCD(LUSC2,LUSC1,WORK(KVEC1),1,LBLK)
        CALL DENSI2(IDENSI,WORK(KRHO1),WORK(KRHO2),
     &       WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC2,EXPS2,
     &       0,XDUM,XDUM,XDUM,XDUM,1)
      END IF
*.    ^ End if density matrices should be calculated
      IF(NPROP.GT.0) THEN
*. Calculate properties 
        KLDUM = 1
        CALL ONE_EL_PROP(1,0,WORK(KLDUM))
      END IF
*.    ^ End if properties should be calculated
      END IF
*     ^ End of expectation values are calculated 
      IF(I_DO_CC3.EQ.1) THEN
        CALL CC3_JACO(WORK(KCC1),WORK(KVEC1),WORK(KVEC2),      
     &                WORK(KCC2))
      END IF
*
*. Test Jacobian by explicit calculation
*
      I_TEST_JACO = 0
      IF(I_TEST_JACO.EQ.1) THEN
        WRITE(6,*) ' Control will be transferred to TEST_JACO'
        CALL JACO_TEST(WORK(KCC1),WORK(KVEC1),WORK(KVEC2),
     &                 WORK(KCC2))
      END IF
*
*
* =======================
* CC excitation energies
* =======================
*
      IF(I_DO_CC_EXC_E.EQ.1) THEN
*. Three more vectors
c        CALL MEMMAN(KCC3,LEN_T_VEC_MX,'ADDL  ',2,'CC3_VE')
c        CALL MEMMAN(KCC4,LEN_T_VEC_MX,'ADDL  ',2,'CC4_VE')
c        CALL MEMMAN(KCC5,LEN_T_VEC_MX,'ADDL  ',2,'CC5_VE')
*
        CALL CC_EXC_E(WORK(KCC1),WORK(KCC2),WORK(KVEC1),WORK(KVEC2),
     &                LUCCAMP,II_RES_EXC)
      END IF
*
* =====================================
* Reform CC to CI form and store on LUC
* =====================================
*
      IF(I_TRANS_WF.EQ.1) THEN
*. CI will be used in the following, 
*. Regenerate inactive Fock matrix - different definition of 
*. inactive Fock matrix may have been used in CC calculation
        CALL COPVEC(WORK(KINT1O),WORK(KINT1),NINT1)
        ICC_EXC = 0
        I_USE_SIMTRH = 0
        CALL FI(WORK(KINT1),ECCX,1)
        ECORE = ECORE_INI

        WRITE(6,*) ' cc wf is normalized and transferred to LUC '
        MX_TERM = 100
        ICC_EXC = 1
        XCONV = 1.0D-20
        I_USE_SIMTRH = 0
*
        ICSPC = IETSPC
        ISSPC = IETSPC
*
        CALL EXPT_REF(LUC,LUSC1,LUHC,LUSC2,LUSC3,XCONV,MX_TERM,
     &               WORK(KVEC1),WORK(KVEC2),CCTYPE)
        CC = INPRDD(WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC1 ,1,LBLK)
        XNORM = SQRT(CC)
        FACTOR = 1.0D0/XNORM
        CALL SCLVCD(LUSC1,LUC,FACTOR,WORK(KVEC1),1,-1)
      END IF
*
      IF(II_DO_CMPCCI.EQ.1) THEN
*. Compare the expanded CI vector with the CI vector on LU17
C ANADIFF(C1,C2,INSPC,LUC1,LUC2,ICISTR)
        LU17 = 17
        CALL ANADIFF(WORK(KVEC1),WORK(KVEC2),ICISPC,ISM,LU17,LUC)
      END IF
*
      ! restore integrals
      if (isvmem.eq.0.and.(i_bcc.eq.1.or.i_oocc.eq.1.or.i_obcc.eq.1))
     &     then
        call vec_from_disc(work(kint1),nint1,1,lblk,lu1into)

        call vec_from_disc(work(kint2),nint2,1,lblk,lu2into)
        CALL RELUNIT(LU1INTO ,'delete')
        CALL RELUNIT(LU2INTO ,'delete')
      end if

      WRITE(6,*)  ' Returning from LUCIA_GENCC '
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GENCC ')
*  
      ! reset entry on common oper
      i_unrorb = 0
      ispcas = 0

      ! reset krho1
      if (nspin.eq.2) krho1 = krho1old

c      if (i_bcc.eq.1) then
c        call save_wo_t1(lu_ccamp,n_cc_amp,work(kvec1))
c        call save_wo_t1(lu_lamp,n_cc_amp,work(kvec1))
c      end if

*
      CALL RELUNIT(LU_1DEN,'delete')
      CALL RELUNIT(LU_2DEN,'delete')
      CALL RELUNIT(LU_LAMP,'keep')
      CALL RELUNIT(LUCCAMP,'keep')
      CALL RELUNIT(LU_OMG ,'delete')

      RETURN
      END
      subroutine nonred_os(ittact,iooexc,iooexcc,
     &                     nooexcl,itspc,ntoob,ngas)
*     little slave routine to correctly address the integer arrays
      implicit none

      integer, intent(inout) ::
     &     ittact(*), iooexc(*), iooexcc(*), nooexcl(*)
      integer, intent(in) ::
     &     ngas, ntoob, itspc

      integer ::
     &     ispin

      do ispin = 1, 2
        call nonred_tt_exc_cc(ittact((ispin-1)*ngas*ngas+1),itspc,ispin)
        ! manipulate here ttact array for open-shell!
        ! nonredundant orbital excitations
        call nonred_oo_exc(nooexcl(ispin),
     &                     iooexc((ispin-1)*ntoob*ntoob+1),
     &                     iooexcc((ispin-1)*2*nooexcl(1)+1),
     &                     1,ittact((ispin-1)*ngas*ngas+1),2)
      end do

      return
      end

      SUBROUTINE GENCC_F_DIAG_M(I_RES_OR_UNRES,
     &           ITSS_TP,NTSS_TP,CCDIA,ISM,
     &           XOHSS,NOOEXCL,IOOEXCC,NORBHSS,
     &           VEC1,VEC2,MX_ST_TSOSO_MX,MX_ST_TSOSO_BLK_MX)
*
* Obtain diagonal <0![F,T]!0>
*
* Jeppe Olsen, Summer of 99, 
*              Modified to allow for alpha- and beta- dependent F matrices,
*              Oct 2002 ( to see if this improves converegence for Open shell)
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'oper.inc'
      DIMENSION XOHSS(*),NOOEXCL(2),IOOEXCC(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GCC_FD')
      CALL MEMMAN(KFDIA_AL,NTOOB,'ADDL  ',2,'FDIA_A')
      CALL MEMMAN(KFDIA_BE,NTOOB,'ADDL  ',2,'FDIA_B')
*. Four blocks of string occupations
      CALL MEMMAN(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
*. Four arrays ( will be used to store sums of orbital energies)
      CALL MEMMAN(KLSTR1,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR1 ')
      CALL MEMMAN(KLSTR2,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR2 ')
      CALL MEMMAN(KLSTR3,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR3 ')
      CALL MEMMAN(KLSTR4,MX_ST_TSOSO_MX,'ADDL  ',2,'STAR4 ')
*
* 1 : Obtain Inactive Fock matrix 
*
*. One-body density matrix
      LBLK = -1
        IFUSK = 1
      IF(IFUSK.EQ.1) THEN
c        DO I = 1, 100
c          WRITE(6,*) ' One-particle density will not be calculated !!! '
c        END DO
c        WRITE(6,*) ' Jeppe : I hope you know what you are doing ! '
      ELSE
        CALL COPVCD(LUC,LUSC1,VEC1,1,LBLK)
        XDUM = 0.0D0
        CALL DENSI2(1,WORK(KRHO1),WORK(KRHO2),VEC1,VEC2,
     &               LUSC1,LUC,EXPS2,0,XDUM,XDUM,XDUM,XDUM,1)
      END IF

*. And Fock matrix
      IF(I_RES_OR_UNRES.EQ.1) THEN
*. Standard spinrestricted Fock matrices
        CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
        CALL FIFAM(WORK(KFI))
        CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
        ECORE_H = 0.0D0
*. Extract diagonal of fock-matrix
        CALL GT1DIS(WORK(KFDIA_AL),IREOTS,WORK(KPINT1),
     &            WORK(KFI),ISMFTO,IBSO,NACOB)
        CALL COPVEC(WORK(KFDIA_AL),WORK(KFDIA_BE),NACOB)
        
        IF (NORBHSS.GT.0) THEN
          WRITE(6,*) 'NOT IMPLEMENTED (BUT EASY-PEASY)'
          STOP 'GENCC_F_DIAG_M'
        END IF

      ELSE 
*. Spin unrestricted Fock matrices 
        I_UNRORB_MERK = I_UNRORB
        I_UNRORB = 0
        IF (I_RES_OR_UNRES.EQ.2) THEN
          CALL FI_HS(WORK(KINT1O),WORK(KFI_AL),WORK(KFI_BE),ECORE_X,1)
        ELSE
          ! version for unrestricted orbitals
          CALL FI_HS_AB(WORK(KINT1O),
     &                  WORK(KFI_AL),WORK(KFI_BE),ECORE_X,1)
        END IF
*. Extract diagonal of fock-matrix
        CALL GT1DIS(WORK(KFDIA_AL),IREOTS,WORK(KPINT1),
     &            WORK(KFI_AL),ISMFTO,IBSO,NACOB)
        CALL GT1DIS(WORK(KFDIA_BE),IREOTS,WORK(KPINT1),
     &            WORK(KFI_BE),ISMFTO,IBSO,NACOB)
        I_UNRORB = I_UNRORB_MERK

*. Build the approx. orbital hessian from diagonal fock-matrix
        IF (NORBHSS.EQ.1) FAC = 4d0
        IF (NORBHSS.EQ.2) FAC = 2d0
        IF (NORBHSS.GT.0) THEN
          IDX = 0
          DO ISPIN = 1, NORBHSS
            IF (ISPIN.EQ.1) KOFF = KFDIA_AL-1
            IF (ISPIN.EQ.2) KOFF = KFDIA_BE-1
            IOFF = (ISPIN-1)*NOOEXCL(1)*2
            DO IEXC = 1, NOOEXCL(ISPIN)
              II = IOOEXCC(IOFF+(IEXC-1)*2+1)
              JJ = IOOEXCC(IOFF+(IEXC-1)*2+2)
              IDX = IDX+1
              XOHSS(IDX) = FAC*(WORK(KOFF+II) - WORK(KOFF+JJ))
            END DO
          END DO
          
        END IF


        NTEST = 00
        IF (NTEST.GE.100) THEN
          WRITE(6,*) 'Approx. orbital Hessian:'
          IDX = 0
          DO ISPIN = 1, NORBHSS
            IF (NSPIN.EQ.2.AND.ISPIN.EQ.1) WRITE(6,*) 'Alpha:'
            IF (NSPIN.EQ.2.AND.ISPIN.EQ.2) WRITE(6,*) 'Beta:'
            IOFF = (ISPIN-1)*NOOEXCL(1)*2
            DO IEXC = 1, NOOEXCL(ISPIN)
              IDX = IDX+1
              WRITE(6,'(x,2i5,g20.10)')
     &                    IOOEXCC(IOFF+(IEXC-1)*2+1),
     &                    IOOEXCC(IOFF+(IEXC-1)*2+2),
     &                    XOHSS(IDX)
            END DO
          END DO
        END IF

      END IF
*. And then the diagonal 
      CALL GENCC_F_DIAG(ITSS_TP,NTSS_TP,CCDIA,ISM,
     &     WORK(KFDIA_AL),WORK(KFDIA_BE),
     &     WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &     WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &     WORK(KLSTR1),WORK(KLSTR2),WORK(KLSTR3),WORK(KLSTR4))

*
C?    WRITE(6,*) ' Enforced stop in GENCC_F_DIAG_M '
C?    STOP '  Enforced stop in GENCC_F_DIAG_M'
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GCC_FD')
*
      RETURN
      END 
      SUBROUTINE GENCC_F_DIAG(ITSS_TP,NTSS_TP,CCDIA,ISM,FDIA_AL,FDIA_BE,
     &           IOCC_CA,IOCC_CB,IOCC_AA,IOCC_AB,
     &           E_CA,E_CB,E_AA,E_AB) 
*
* <0![F,T]!0)
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cc_exc.inc'
*. Specific input
      INTEGER ITSS_TP(4*NGAS,NTSS_TP)
      DIMENSION FDIA_AL(*), FDIA_BE(*)
*. Output
      DIMENSION CCDIA(*)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
      DIMENSION E_CA(*), E_CB(*), E_AA(*), E_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*
      REAL(8) INPROD
*
      IT = 0
      DO ITSS = 1, NTSS_TP
*. Transform from occupations to groups
       CALL OCC_TO_GRP(ITSS_TP(1+0*NGAS,ITSS),IGRP_CA,1      )
       CALL OCC_TO_GRP(ITSS_TP(1+1*NGAS,ITSS),IGRP_CB,1      )
       CALL OCC_TO_GRP(ITSS_TP(1+2*NGAS,ITSS),IGRP_AA,1      )
       CALL OCC_TO_GRP(ITSS_TP(1+3*NGAS,ITSS),IGRP_AB,1      )
*
       NEL_CA = IELSUM(ITSS_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ITSS_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ITSS_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ITSS_TP(1+3*NGAS,ITSS),NGAS)
*. Diagonal block ?
       CALL DIAG_EXC_CC(ITSS_TP(1+0*NGAS,ITSS),
     &                  ITSS_TP(1+1*NGAS,ITSS),
     &                  ITSS_TP(1+2*NGAS,ITSS),
     &                  ITSS_TP(1+3*NGAS,ITSS),NGAS,IDIAG)
       IF(MSCOMB_CC.EQ.0.OR.IDIAG.EQ.0) THEN
         IRESTRICT = 0
       ELSE 
         IRESTRICT = 1
       END IF
       DO ISM_C = 1, NSMST
         ISM_A = MULTD2H(ISM,ISM_C) 
         DO ISM_CA = 1, NSMST
           ISM_CB = MULTD2H(ISM_C,ISM_CA)
           DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*
            ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
            ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
            IF(IRESTRICT.EQ.1.AND.ISM_BETA.GT.ISM_ALPHA) GOTO 777
            IF(IRESTRICT.EQ.0.OR.ISM_ALPHA.GT.ISM_BETA) THEN
             IRESTRICT_LOOP = 0
            ELSE
             IRESTRICT_LOOP = 1
            END IF
C?          WRITE(6,'(A,4I5)') ' ISM_CA, ISM_CB, ISM_AA, ISM_AB',
C?   &                           ISM_CA, ISM_CB, ISM_AA, ISM_AB
*. obtain strings
            CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,NSTR_CA,
     &           IOCC_CA, NORBT,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,NSTR_CB,
     &           IOCC_CB, NORBT,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,NSTR_AA,
     &           IOCC_AA, NORBT,0,IDUM,IDUM)
            CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,NSTR_AB,
     &           IOCC_AB, NORBT,0,IDUM,IDUM)
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
*. Set up sums of F elements
             CALL SUM_FDIA_FOR_STR(E_CA,FDIA_AL,NSTR_CA,NEL_CA,IOCC_CA)
             CALL SUM_FDIA_FOR_STR(E_CB,FDIA_BE,NSTR_CB,NEL_CB,IOCC_CB)
             CALL SUM_FDIA_FOR_STR(E_AA,FDIA_AL,NSTR_AA,NEL_AA,IOCC_AA)
             CALL SUM_FDIA_FOR_STR(E_AB,FDIA_BE,NSTR_AB,NEL_AB,IOCC_AB)
*. Loop over T elements as  matric T(I_CA, I_CB, IAA, I_AB)
             DO I_AB = 1, NSTR_AB
              IF(IRESTRICT_LOOP.EQ.1) THEN
                I_AA_MIN = I_AB
              ELSE
                I_AA_MIN = 1
              END IF
              DO I_AA = I_AA_MIN, NSTR_AA
               DO I_CB = 1, NSTR_CB
                IF(IRESTRICT_LOOP.EQ.1.AND.I_AB.EQ.I_AA) THEN
                 ICA_MIN = I_CB
                ELSE
                 ICA_MIN = 1
                END IF
                DO I_CA = ICA_MIN, NSTR_CA
                 IT = IT + 1
                 CCDIA(IT) = 
     &           E_CA(I_CA) + E_CB(I_CB) - E_AA(I_AA) - E_AB(I_AB)
                 YTEST = 1.D-10
*. A temp measure for eliminating some redundant operations 
                 IF(ABS(CCDIA(IT)).LE.YTEST) THEN
                   CCDIA(IT) = 1.0D+10
                 END IF
C?               WRITE(6,'(A,4I4)') ' I_AB, I_AA, I_CB, I_CA', 
C?   &                                I_AB, I_AA, I_CB, I_CA 
C?               WRITE(6,*) ' IT and Diag(IT)', IT,CCDIA(IT)
                END DO
               END DO
              END DO
             END DO
*            ^ End of loop over elements of block
  777       CONTINUE
            END DO
*           ^ End of loop over ISM_AA
         END DO
*        ^ End of loop over ISM_CA
       END DO
*      ^ End of loop over ISM_C
      END DO
*     ^ End of loop over ITSS
*
      NTEST =  10
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' <0![F,T]!0> constructed '
        WRITE(6,*) ' Number of elements ', IT
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' The diagonal   '
        WRITE(6,*) ' |diag| = ',SQRT(INPROD(CCDIA,CCDIA,IT))
c        CALL WRTMAT(CCDIA,1,IT,1,IT)
        CALL WRT_CC_VEC2(CCDIA,6,'GEN_CC')
      END IF
*
      RETURN
      END
      SUBROUTINE SUM_FDIA_FOR_STR(FSUM,FDIA,NSTR,NEL,IOCC)
*
* F(ISTR) = SUM(IEL) F(IOCC(IEL,ISTR)
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION FDIA(*)
      INTEGER IOCC(NEL,NSTR)
*. Output
      DIMENSION FSUM(NSTR)
*
      DO ISTR = 1, NSTR
        X = 0.0D0
        DO IEL = 1, NEL
          X = X + FDIA(IOCC(IEL,ISTR))
        END DO
        FSUM(ISTR) = X
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' First two elements of Fdia in Orb basis '
        CALL WRTMAT(FDIA,1,2,1,2)
        WRITE(6,*) ' Occ of Strings '
        CALL IWRTMA(IOCC,NEL,NSTR,NEL,NSTR)
        WRITE(6,*) ' FSUM in Str basis '
        CALL WRTMAT(FSUM,1,NSTR,1,NSTR)
      END IF
        
*
      RETURN 
      END
      SUBROUTINE IDIM_TCC(ITSOSO_TP,NTSOSO_TP,ISYM,
     &           MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &           LTSOSO_TP,IBTSOSO_TP,IDIM_T,MSCOMB_CC,
     &           MX_TBLK_AS,ISPOX_FOR_OCCLS,NOCCLS,IBSPOX_FOR_OCCLS,
     &           NTCONF,IPRCC)
*
* Dimension of T operators in spin-orbital basis
*
* Largest number of strings of given sym in T(ICA,ICB,IAA,IAB)
* i.e. largest block of ICA, ICB,IAA,IAB of given sym 
*
* Size of block required to hold above blocks
*
* Number of T-configurations is also calculated
*
* Jeppe Olsen, Summer of 99
*              T-configurations added summer of 01
*              (Does not work, neither in theory or practive 
*               and should be removed !) 
*
c      INCLUDE 'implicit.inc'
*. General input
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
*. Specific input
      INTEGER ITSOSO_TP(4*NGAS,NTSOSO_TP)
      INTEGER ISPOX_FOR_OCCLS(*),IBSPOX_FOR_OCCLS(NOCCLS)
*. Local scratch 
      INTEGER ICA_GRP(MXPNGAS), ICB_GRP(MXPNGAS)
      INTEGER IAA_GRP(MXPNGAS), IAB_GRP(MXPNGAS)
*.Output : Length of each type and offset for each type
      INTEGER LTSOSO_TP(NTSOSO_TP), IBTSOSO_TP(NTSOSO_TP)
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRCC)
C?    WRITE(6,*) ' IDIM_TCC... : IPRCC = ', IPRCC
*
C?    WRITE(6,*) ' IDIM_TCC : NOCCLS = ', NOCCLS
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' IDIM_TCC : MSCOMB_CC, ISYM = ', MSCOMB_CC,ISYM
      END IF
*
      LENGTH = 0
      MX_ST_TSOSO = 0
      MX_ST_TSOSO_BLK = 0
      MX_TBLK = 0
      MX_SBSTR = 0
      MX_TBLK_AS = 0
      NTCONF = 0
*
      IF(NTEST.GE.100) WRITE(6,*) 'NTSOSO_TP =', NTSOSO_TP
      DO ITSS_TP = 1, NTSOSO_TP
       IF(NTEST.GE.100) WRITE(6,*) ' ITSS_TP = ', ITSS_TP 
*. Occupation to group translation 
       CALL OCC_TO_GRP(ITSOSO_TP(1+0*NGAS,ITSS_TP),ICA_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+1*NGAS,ITSS_TP),ICB_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+2*NGAS,ITSS_TP),IAA_GRP,1)
       CALL OCC_TO_GRP(ITSOSO_TP(1+3*NGAS,ITSS_TP),IAB_GRP,1)
C           OCC_TO_GRP(IOCC,IGRP,LGRP)
*
       NEL_CA = IELSUM(ITSOSO_TP(1+0*NGAS,ITSS_TP),NGAS)
       NEL_CB = IELSUM(ITSOSO_TP(1+1*NGAS,ITSS_TP),NGAS)
       NEL_AA = IELSUM(ITSOSO_TP(1+2*NGAS,ITSS_TP),NGAS)
       NEL_AB = IELSUM(ITSOSO_TP(1+3*NGAS,ITSS_TP),NGAS)
C?     WRITE(6,*) ' NEL_CA, NEL_CB, NEL_AA, NEL_AB=',
C?   &              NEL_CA, NEL_CB, NEL_AA, NEL_AB
C     
*. Check if the alpha- and beta part are identical 
*
C     DIAG_EXC_CC(NCA,NCB,NAA,NAB,NGAS,IDIAG)
       IF(MSCOMB_CC.NE.0) THEN
        CALL DIAG_EXC_CC(ITSOSO_TP(1+0*NGAS,ITSS_TP),
     &                   ITSOSO_TP(1+1*NGAS,ITSS_TP),
     &                   ITSOSO_TP(1+2*NGAS,ITSS_TP),
     &                   ITSOSO_TP(1+3*NGAS,ITSS_TP),NGAS,IDIAG)
       ELSE
        IDIAG = 0
       END IF
       LENGTH2 = 0
       IBTSOSO_TP(ITSS_TP) = LENGTH + 1
       DO I_CR_SM = 1, NSMST
        I_AN_SM = MULTD2H(ISYM,I_CR_SM)
        DO I_CR_AL_SM = 1, NSMST
          DO I_AN_AL_SM = 1, NSMST
*
            I_CR_BE_SM = MULTD2H(I_CR_AL_SM,I_CR_SM)
            I_AN_BE_SM = MULTD2H(I_AN_AL_SM,I_AN_SM)
*
            I_AL_SYM = (I_AN_AL_SM-1)*NSMST + I_CR_AL_SM
            I_BE_SYM = (I_AN_BE_SM-1)*NSMST + I_CR_BE_SM
            IF(IDIAG.EQ.0.OR.I_AL_SYM.GE.I_BE_SYM) THEN
*
             CALL NST_SPGRP2(NGAS   ,ICA_GRP,I_CR_AL_SM,NSMST,
     &                       LEN_CA,NDIST_CA)
             CALL NST_SPGRP2(NGAS   ,ICB_GRP,I_CR_BE_SM,NSMST,
     &                       LEN_CB,NDIST_CB)
             CALL NST_SPGRP2(NGAS   ,IAA_GRP,I_AN_AL_SM,NSMST,
     &                      LEN_AA,NDIST_AA)
             CALL NST_SPGRP2(NGAS   ,IAB_GRP,I_AN_BE_SM,NSMST,
     &                      LEN_AB,NDIST_AB)
*
c             IF(NTEST.GE.1000) THEN
c               WRITE(6,'(A,4I5)') ' sym of CA CB AA AB',
c     &         I_CR_AL_SM, I_CR_BE_SM, I_AN_AL_SM, I_AN_BE_SM
c               WRITE(6,'(A,4I5)') ' LEN_CA, LEN_CB, LEN_AA, LEN_AB ',
c     &                      LEN_CA, LEN_CB, LEN_AA, LEN_AB
c             END IF
*
             IF( IDIAG.EQ.0.OR.I_AL_SYM.GT.I_BE_SYM) THEN
                   LENGTH = LENGTH + LEN_CA*LEN_CB*LEN_AA*LEN_AB
                   LENGTH2 = LENGTH2 + LEN_CA*LEN_CB*LEN_AA*LEN_AB
c test
                   LENGTH3 = LEN_CA*LEN_CB*LEN_AA*LEN_AB
             ELSE IF(IDIAG.EQ.1.AND.I_AL_SYM.EQ.I_BE_SYM) THEN
                   LL = LEN_CA*LEN_AA
                   LENGTH = LENGTH + LL*(LL+1)/2
                   LENGTH2 = LENGTH2 + LL*(LL+1)/2
c test
                   LENGTH3 = LL*(LL+1)/2
             END IF

             IF(NTEST.GE.1000.AND.LENGTH3.GT.0) THEN
               WRITE(6,'(/A,4I5)') ' sym of CA CB AA AB',
     &         I_CR_AL_SM, I_CR_BE_SM, I_AN_AL_SM, I_AN_BE_SM
               WRITE(6,'(A,4I5/)') ' LEN_CA, LEN_CB, LEN_AA, LEN_AB ',
     &                      LEN_CA, LEN_CB, LEN_AA, LEN_AB,'>>',length3
             END IF
*. Assumes that combination are not used for spin-restricted OSCC
CERR         IF(INI_T.EQ.1) 
CERR &       NTCONF = NTCONF + LEN_CA*LEN_CB*LEN_AA*LEN_AB
*
             MX_ST_TSOSO = 
     &       MAX(MX_ST_TSOSO,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
             MX_ST_TSOSO_BLK = 
     &       MAX(MX_ST_TSOSO_BLK,LEN_CA*NEL_CA,LEN_CB*NEL_CB,
     &                           LEN_AA*NEL_AA,LEN_AB*NEL_AB)
             MX_TBLK = MAX(MX_TBLK,LEN_CA*LEN_CB*LEN_AA*LEN_AB)
             MX_SBSTR = MAX(MX_SBSTR,LEN_CA,LEN_CB,LEN_AA,LEN_AB)
            END IF
*           ^ End if symmetry combination should be included 
          END DO
        END DO
       END DO
      LTSOSO_TP(ITSS_TP) = LENGTH2
      MX_TBLK_AS = MAX(MX_TBLK_AS,LENGTH2)
      END DO
*
      IDIM_T = LENGTH
*
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' Number of T-coefficients ', LENGTH
CERR    WRITE(6,*) ' Number of T-configurations', NTCONF
        WRITE(6,*) ' Largest symmetry block of T ', MX_TBLK
        WRITE(6,*) ' Largest block of T ', MX_TBLK_AS
        WRITE(6,*) ' Largest substring : ', MX_SBSTR
    
      END IF
*
      IF(NTEST.GE.20) THEN
        WRITE(6,*) 
     &  ' Largest number of strings of given sym and type ',
     &               MX_ST_TSOSO
        WRITE(6,*) 
     & ' Largest block of strings of given sym and type ',
     &               MX_ST_TSOSO_BLK
        WRITE(6,*) ' Number of elements per block '
        CALL IWRTMA(LTSOSO_TP,1,NTSOSO_TP,1,NTSOSO_TP)
        WRITE(6,*) ' Offset for each block '
        CALL IWRTMA(IBTSOSO_TP,1,NTSOSO_TP,1,NTSOSO_TP)
      END IF
*
      RETURN
      END 
      SUBROUTINE TP_OBEX2(NOCCLS,NEL,NGASX,IOBEX_TP,LCOBEX_TP,LAOBEX_TP,
     &                   IOCCLS,IOCCLS_REF,MX_NCREA,MX_NANNI,
     &                   MX_EXC_LEVEL,IEXTP_TO_OCCLS,MX_AAEXC,IFLAG,
     &                   I_OOCC,NOBEX_TP,NOAAEX,IPRCC)
*
* Obtain the orbital excitation types needed to generate occupation classes 
* in IOCCLS from IOCCLS_REF
* 
* Jeppe Olsen, Updated version of TP_OBEX, March 2000, Still on the train
*              April 2001         unit operator added as excitation 
*                                 NOBEX_TP + 1. The arrays 
*                                 IOBEX_TP, LCOBEX, LAOBEX, IEXTP_TO_OCCLS
*                                 should this be dimensioned with NOBEX_TP+1
*                                 
*
* Allows active-active excitations. Active orbitals are assumed to 
* be in a single orbital space, 
*
* If IFLAG = 1, only the number of orbital excitation types is generated
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IOCCLS(NGAS,*), IOCCLS_REF(NGAS)
*. Output
      INTEGER IOBEX_TP(2*NGAS,*),  LCOBEX_TP(*), LAOBEX_TP(*)
      INTEGER IEXTP_TO_OCCLS(*)
*. Local scratch
      DIMENSION ICREA(MXPNGAS),IANNI(MXPNGAS)
*. Number of active orbital spaces
      NACT_SPC = 0
      IACT_SPC = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) THEN
          NACT_SPC = NACT_SPC + 1
          IACT_SPC = IGAS
        END IF
      END DO
*
      IF(NACT_SPC.GT.1) THEN
        WRITE(6,*) ' TP_OBEX2 in problems '
        WRITE(6,*) ' More than one active orbital spaces '
        WRITE(6,*) ' NACT_SPC = ',  NACT_SPC 
        STOP ' TP_OBEX2 :  More than one active orbital spaces '
      END IF
C?    WRITE(6,*) ' TP_OBEX2 : IACT_SPC,NACT_SPC',IACT_SPC,NACT_SPC
C?    WRITE(6,*) ' TP_OBEX2 : MX_AAEXC ', MX_AAEXC
*
* The orbital excitation operator IEXTP is  organized as 
*
* LCOBEX(IEXTP) : Number of creation operators    
* LAOBEX(IEXTP) : Number of annihilation operators
* IOBEX_TP(1 - NGAS, IEXTP) : Number of creation operators per gassspace
* IOBEX(NGAS+1  -  2*NGAS, , IEXTP) : Number of annihilation operators 
*
* IEXTP_TO_OCCLS is map from orbital excitation type to occupation 
* class for CI coefficients
      NTEST = 0000
      NTEST = MAX(NTEST,IPRCC)
      MX_NCREA = 0
      MX_NANNI = 0
      NOBEX_TP = 0
      MX_EXC_LEVEL = 0
      JREFCLS = 0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'TP_OBEX speaking '
        WRITE(6,*) ' Reference occupation class ' 
        CALL IWRTMA(IOCCLS_REF,1,NGAS,1,NGAS)
        WRITE(6,*) '  IPRCC = ', IPRCC
      END IF
*
      IZERO = 0
      DO JOCCLS = 1, NOCCLS 
        NANNI = 0
        NCREA = 0
        CALL ISETVC(ICREA,IZERO,NGAS)
        CALL ISETVC(IANNI,IZERO,NGAS)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Excited occupation class '
          CALL IWRTMA(IOCCLS(1,JOCCLS),1,NGAS,1,NGAS)
        END IF
        DO IGAS = 1, NGAS
          IF(IOCCLS(IGAS,JOCCLS).GT.IOCCLS_REF(IGAS)) THEN
            ICREA(IGAS) = IOCCLS(IGAS,JOCCLS) - IOCCLS_REF(IGAS)
            NCREA = NCREA + ICREA(IGAS)
          ELSE IF (IOCCLS(IGAS,JOCCLS).LT.IOCCLS_REF(IGAS)) THEN
            IANNI(IGAS) = -( IOCCLS(IGAS,JOCCLS) - IOCCLS_REF(IGAS))
            NANNI = NANNI + IANNI(IGAS)
          END IF
        END DO
        IF(NANNI.EQ.0.AND.NCREA.EQ.0) JREFCLS = JOCCLS
*. Add active -active excitation 
        DO IAA_EXC = 0, MX_AAEXC
         ITSOKAY = 1
         IF(IAA_EXC.GT.0) THEN
*. Can another active-active excitation be added ?
          IF(IANNI(IACT_SPC)+1.LE.MX_AAEXC.AND.
     &       ICREA(IACT_SPC)+1.LE.MX_AAEXC) THEN
              ITSOKAY = 1
          ELSE 
              ITSOKAY = 0
          END IF
         END IF
         IF(ITSOKAY.EQ.1) THEN
*
          IF(IAA_EXC.GT.0) THEN
            ICREA(IACT_SPC)  = ICREA(IACT_SPC) + 1
            IANNI(IACT_SPC)  = IANNI(IACT_SPC) + 1
            NCREA = NCREA + 1
            NANNI = NANNI + 1
          END IF
*     If we do not include pure active-active rotations:
          IF (NOAAEX.EQ.1) THEN
*     Test whether this is one:
            IPAA = 1
            DO JGAS = 1, NGAS
              IF (JGAS.NE.IACT_SPC.AND.ICREA(JGAS).NE.0) IPAA=0
              IF (JGAS.NE.IACT_SPC.AND.IANNI(JGAS).NE.0) IPAA=0
            END DO
            IF (IPAA.EQ.1) CYCLE
          END IF
*
          MX_EXC_LEVEL = MAX(MX_EXC_LEVEL,NCREA)
          IF(NCREA+NANNI.NE.0.AND.
     &         (I_OOCC.EQ.0.OR.(NCREA.NE.1.AND.NANNI.NE.1))) THEN 
            NOBEX_TP = NOBEX_TP + 1
            IF(IFLAG.NE.1) THEN
              LCOBEX_TP(NOBEX_TP) = NCREA
              LAOBEX_TP(NOBEX_TP) = NANNI
              IEXTP_TO_OCCLS(NOBEX_TP) = JOCCLS
              CALL ICOPVE(ICREA,IOBEX_TP(1,NOBEX_TP),NGAS )
              CALL ICOPVE(IANNI,IOBEX_TP(NGAS+1,NOBEX_TP),NGAS )
            END IF
          END IF
*
         END IF
       END DO
*      ^ End of loop over active-active excitations
      END DO
*
      IF(IFLAG.NE.1) THEN
*. Add unit operator as excition NOBEX_TP + 1
        LCOBEX_TP(NOBEX_TP + 1) = 0
        LAOBEX_TP(NOBEX_TP + 1) = 0
        IEXTP_TO_OCCLS(NOBEX_TP+1) = JREFCLS
        IZERO = 0
        CALL ISETVC(IOBEX_TP(1,NOBEX_TP+1),IZERO,NGAS)
        CALL ISETVC(IOBEX_TP(NGAS+1,NOBEX_TP+1),IZERO,NGAS)
      END IF
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Largest excitation level : ', MX_EXC_LEVEL
        WRITE(6,*)
        WRITE(6,*) ' Number of types of orbital excitations ', NOBEX_TP 
        WRITE(6,*)
*
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' Creation part,  Annihilation  part '
          WRITE(6,*) ' ==================================='
          DO IOBEX = 1, NOBEX_TP+1
            WRITE(6,'(16I4,16I4)')
     &      (IOBEX_TP(I,IOBEX),I=1, NGAS),
     &      (IOBEX_TP(NGAS+I,IOBEX),I=1, NGAS) 
          END DO
*
          WRITE(6,*) ' Orbital excitation type to occupation class '
          CALL IWRTMA(IEXTP_TO_OCCLS,1,NOBEX_TP+1,1,NOBEX_TP+1)
        END IF
*
      END IF
*
      RETURN
      END
      FUNCTION  NAB_COMP_FOR_OBOP(NCREA,NANNI,MS2TOT)
*
* Number of possible spincomponente of orbital operator containing 
* Ncrea excitation operators, Nanni annihilation operators
* and total spin projection MSTOT
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
*
      NUM = 0
*
      DO NCREA_AL = 1, NCREA
        NCREA_BE = NCREA-NCREA_AL
        MS2_CREA = NCREA_AL - NCREA_BE 
        MS2_ANNI = MS2TOT - MS2_CREA
*
        IF(MS2_ANNI.GT.NANNI) THEN
*. For annihilations it holds that Nalpha + Nbeta = Nanni
*                                  Nbeta  - Nalpha= MS2_ANNI
*
          NANNI_BE = (MS2_ANNI + NANNI)/2
          NANNI_AL = NANNI - NANNI_BE
*. Number of operators is given by number of ways one can 
*. distribute the alpha electrons
          NUM = NUM + IBION(NCREA,NCREA_AL)*IBION(NANNI,NANNI_AL)
        END IF
      END DO
*
      NAB_COMP_FOR_OBOP = NUM                  
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from NAB_COMP_FOR_OBOP '
        WRITE(6,*) ' Input : NCREA, NANNI, MS2TOT : ',
     &                       NCREA, NANNI, MS2TOT 
        WRITE(6,*) 
     &  ' Output : Number of spincomponents : ', NAB_COMP_FOR_OBOP
      END IF
*
      RETURN
      END
      SUBROUTINE AB_COMP_FOR_OBOP(NCREA,NANNI,MS2TOT,NCOMP,IABCOMP,
     &                            NAEL,NBEL,IFLAG)
*
* Spinprojection components of orbital operator containing 
* NCREA creation operators and NANNI annihilation operators. 
*
* (Modified SPNCOM routine )
*
* Jeppe Olsen,  Summer of 99
*
      INCLUDE 'implicit.inc'
      INTEGER ADD
      INCLUDE 'mxpdim.inc'
*. Output
      DIMENSION IABCOMP(NCREA+NANNI,*)
*. Local scratch
      DIMENSION IWORK(2*MXPORB)
*
      NTEST = 00
C     WRITE(6,*) ' IFLAG, NAEL, NBEL at start of AB... ', 
C    &             IFLAG, NAEL, NBEL
      NCOMP=0
*
* combinations are considered as binary numbers,1=alpha,0=beta
*
      NOPEN = NCREA+NANNI
      MX=2 ** NOPEN
      CALL ISETVC(IWORK,0,NOPEN)
* Loop over all possible binary numbers
      DO 200 I=1,MX
C.. 1 : NEXT BINARY NUMBER
        ADD=1
        J=0
  190   CONTINUE
        J=J+1
        IF(IWORK(J).EQ.1) THEN
          IWORK(J)=0
        ELSE
          IWORK(J)=1
          ADD=0
        END IF
        IF( ADD .EQ. 1 ) GOTO 190
C.. 2 :  CORRECT SPIN PROJECTION ?
        NUP=0
        DO J=1,NCREA
          NUP=NUP+IWORK(J)
        END DO   
        DO J = NCREA+1, NCREA+NANNI 
          IF(IWORK(J).EQ.0) NUP = NUP + 1
        END DO
        NDOWN = NOPEN - NUP
        MS2 = NUP - NDOWN
*. Number of alpha and beta annihilations must not be greater than 
*       the number of electrons with this spinprojection
        MALPHA = 0
        MBETA = 0
        DO J = NCREA + 1, NCREA + NANNI
          IF(IWORK(J).EQ.1) THEN 
            MALPHA = MALPHA + 1
          ELSE 
            MBETA = MBETA + 1
          END IF
        END DO
C
        IF(MS2.EQ.MS2TOT.AND.MALPHA.LE.NAEL.AND.MBETA.LE.NBEL) THEN
          NCOMP = NCOMP + 1
C         WRITE(6,*) ' ICOMP, NALPHA, NBETA, MALPHA, MBETA',
C    &                 ICOMP, NALPHA, NBETA, MALPHA, MBETA
          CALL ICOPVE(IWORK,IABCOMP(1,NCOMP),NOPEN)
        END IF
C
  200 CONTINUE
C
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Generation of spincomponents of orbitaloperator'
         WRITE(6,*)
         WRITE(6,*) ' NCREA, NANNI, MS2TOT', NCREA, NANNI,MS2TOT 
         WRITE(6,*) ' Number of terms generated ', NCOMP 
         WRITE(6,*)
         IF(IFLAG.NE.1) THEN
           WRITE(6,*) ' The operators (alpha = 1, beta = 0)'
           WRITE(6,*) ' ==================================='
           WRITE(6,*)
           DO 20 J=1,NCOMP
             WRITE(6,1020) J,(IABCOMP(K,J),K=1,NOPEN)
  20       CONTINUE
         END IF
 1020    FORMAT(1H0,I5,2X,30I2,/,(1H ,7X,30I2))
      END IF
*
      RETURN
      END
* note the preliminary copy obex_to_spobex2
      SUBROUTINE OBEX_TO_SPOBEX(IFLAG,IOBEX_TP,LCOBEX_TP,LAOBEX_TP,
     &                          NOBEX_TP,ISPOBEX_TP,NSPOBEX_TP,NGAS,
     &                          NOBPT,MS2TOT,MSCOMB_CC,IAAEXC,IACT_SPC,
     &                          IPRSTR,ISOX_TO_OX,MXSPOX,NSOX_FOR_OX,
     &                          IBSOX_FOR_OX,ISOX_FOR_OX,NAEL,NBEL,
     &                          IREFSPC)
*
* Orbital excitation types => Spin-orbital excitations
*
* IFLAG = 1 : Just Number of spinorbital excitations
* IFLAG = 2 : Also the actual excitations
*
*. Jeppe Olsen, Summer of 99
*               Updated with Active-Active excitations, March 2000
*               SOX<=>OX mapping added, Spring 2001
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. Input
      INTEGER IOBEX_TP(2*NGAS, NOBEX_TP), NOBPT(*)
      INTEGER LCOBEX_TP(NOBEX_TP), LAOBEX_TP(NOBEX_TP)
*. Output : Creation in alpha, Creation in beta, 
*           Annihilation in alpha, annihilation in beta
      INTEGER ISPOBEX_TP(4*NGAS,*)
*. Spin orbital type => orbital type
      INTEGER ISOX_TO_OX(*)
*. Number of spinorbital excitations for given orbital excitation 
      INTEGER NSOX_FOR_OX(*)
*. And the actual spinorbital excitations for given orbital excitation 
      INTEGER ISOX_FOR_OX(*)
*. Offset to given orbital type in ISOX_FOR_OX 
      INTEGER IBSOX_FOR_OX(*)
*. Local scratch for occ of reference alpha and beta
      INTEGER IREF_AL(MXPNGAS),IREF_BE(MXPNGAS)
*
C?    WRITE(6,*) ' IPRSTR in OBEX_TO_SPOBEX ', IPRSTR
      NTEST = 1000
      NTEST = MAX(NTEST,IPRSTR)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from OBEX_TO_SPOBEX...'
        WRITE(6,*) ' -----------------------------'
        WRITE(6,*) ' LCOBEX_TP, LAOBEX_TP: '
        CALL IWRTMA(LCOBEX_TP,1,NOBEX_TP,1,NOBEX_TP)
        CALL IWRTMA(LAOBEX_TP,1,NOBEX_TP,1,NOBEX_TP)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'OE_SOE')
*. Largest number of creation and annihilation operators
      MX_CREA = IMNMX(LCOBEX_TP,NOBEX_TP,2)
      MX_ANNI = IMNMX(LAOBEX_TP,NOBEX_TP,2)
*. Largest possible block of spincombination
*. Nup + Ndown = MX_OPS, Nup - Ndown = Ms2_TOT
      MX_OPS = MX_CREA + MX_ANNI
      MX_UP = (MX_OPS+MS2TOT)/2
      MX_OPSD2 = MX_OPS/2
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NOBEX_TP = ', NOBEX_TP
        WRITE(6,*) ' MX_CREA, MX_ANNI = ', MX_CREA,MX_ANNI
        WRITE(6,*) ' MX_OPS, MX_UP = ',  MX_OPS, MX_UP
        WRITE(6,*) ' MS2TOT, NAEL, NBEL = ', MS2TOT, NAEL, NBEL
      END IF 
*
C     MX_SOEX_BLK = IBION(MX_OPS,MX_OPSD2)
      MX_SOEX_BLK = IBION(MX_OPS,MX_UP)
      CALL MEMMAN(KLSOEX ,MX_SOEX_BLK*MX_OPS,'ADDL  ',2,'SOEXBL')
*.Alpha and beta-occupations for reference space 
*. Notice : IREFSPC is not used, info in IHPVGAS is used ! 
      IF(IREFSPC.NE.0) THEN
        CALL GET_REF_ALBE_OCC(IREFSPC,IREF_AL,IREF_BE)
      END IF
*. 
      NSPOBEX_TP = 0
      DO ICREA = 0, MX_CREA
        DO IANNI = 0, MX_ANNI
*. Any operators with this number  of creation and annihilation operators ?
          NCOMP = 0
          DO JOBEX_TP = 1, NOBEX_TP
            IF(LCOBEX_TP(JOBEX_TP).EQ.ICREA .AND.
     &         LAOBEX_TP(JOBEX_TP).EQ.IANNI ) NCOMP = NCOMP + 1
          END DO
*
          IF(NCOMP.NE.0) THEN 
*. ICREA creation ops and IANNI annihilation strings with total 
*. spin projection MSTOT
            CALL AB_COMP_FOR_OBOP(ICREA,IANNI,MS2TOT,NABCOMP,
     &           WORK(KLSOEX),NAEL,NBEL,IFLAG)
C?          WRITE(6,*) ' NABCOMP = ', NABCOMP
            DO JOBEX_TP = 1, NOBEX_TP
              IF(LCOBEX_TP(JOBEX_TP).EQ.ICREA .AND.
     &           LAOBEX_TP(JOBEX_TP).EQ.IANNI .AND. 
     &           (MXSPOX.EQ.0.OR.ICREA.LE.MXSPOX)  ) THEN
*. Find the number of the above spincombinations that can be included 
*  for this orbital excitaion 
                 NOPEN = ICREA + IANNI
*
                 CALL ACT_SPOBEX_OBEX(IFLAG,IOBEX_TP(1,JOBEX_TP),
     &                NGAS,NOPEN,NABCOMP,WORK(KLSOEX),NSPCOMP_ACT, 
     &                ISPOBEX_TP(1,NSPOBEX_TP+1),NOBPT,MSCOMB_CC,
     &                IACT_SPC,IAAEXC,IREF_AL,IREF_BE,  
     &                IREFSPC,NAEL,NBEL)
*
                 IF(IFLAG.EQ.2) THEN
C?    WRITE(6,*) ' First SPOBEX_TP (a) '
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)      
                   IOFF = NSPOBEX_TP+1
C?                 WRITE(6,*) ' IOFF, NSPCOMP_ACT = ',
C?   &                          IOFF, NSPCOMP_ACT
                   CALL ISETVC(ISOX_TO_OX(IOFF),JOBEX_TP,NSPCOMP_ACT)
C?    WRITE(6,*) ' First SPOBEX_TP (b) '
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)      
                 END IF
*
                 NSPOBEX_TP =  NSPOBEX_TP + NSPCOMP_ACT
               END IF
            END DO
          END IF
*         ^ End if nonvanishing number of operators 
        END DO
      END DO
*     ^ End of loop over creation and annihilation operators 
C?    WRITE(6,*) ' Memchk after construction of spobex '
C?    CALL MEMCHK
C?    WRITE(6,*) ' Memcheck passed '
*
C?    WRITE(6,*) ' First SPOBEX_TP'
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)      
*
      IF(IFLAG.EQ.2) THEN
*. OX => SOX lists
*. Number of sox per ox
        IZERO = 0
        CALL ISETVC(NSOX_FOR_OX,IZERO,NOBEX_TP)
        DO ISOX = 1, NSPOBEX_TP
          IOX = ISOX_TO_OX(ISOX)
          NSOX_FOR_OX(IOX) = NSOX_FOR_OX(IOX)+1
        END DO
*. offsets for sox's with given ox
        DO IOX = 1, NOBEX_TP
         IF(IOX.EQ.1) THEN
           IBSOX_FOR_OX(IOX) = 1
         ELSE
           IBSOX_FOR_OX(IOX) = IBSOX_FOR_OX(IOX-1) + NSOX_FOR_OX(IOX-1)
         END IF
        END DO
*. and the sox for a given ox
        CALL ISETVC(NSOX_FOR_OX,IZERO,NOBEX_TP)
        DO ISOX = 1, NSPOBEX_TP  
          IOX = ISOX_TO_OX(ISOX)
          NSOX_FOR_OX(IOX) = NSOX_FOR_OX(IOX)+1
          IADR = IBSOX_FOR_OX(IOX) -1 + NSOX_FOR_OX(IOX)
          ISOX_FOR_OX(IADR) = ISOX 
        END DO
      END IF
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ***************************************** '
        WRITE(6,*) ' Information about spinorbital excitations '
        WRITE(6,*) ' ***************************************** '
        WRITE(6,*)
        WRITE(6,*) ' Total number ', NSPOBEX_TP 
        IF(IFLAG.EQ.2) THEN
          DO JSPCOMP_ACT = 1, NSPOBEX_TP  
            WRITE(6,*)
            WRITE(6,*) ' Included spinorbitalexcitation ', JSPCOMP_ACT
            WRITE(6,'(A,16I4)') 
     &      ' Creation of alpha     :', 
     &      (ISPOBEX_TP(I+0*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Creation of beta      :', 
     &      (ISPOBEX_TP(I+1*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Annihilation of alpha :', 
     &      (ISPOBEX_TP(I+2*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Annihilation of beta  :', 
     &      (ISPOBEX_TP(I+3*NGAS,JSPCOMP_ACT),I=1,NGAS)
          END DO
*
          IF(NTEST.GE.5) THEN
           WRITE(6,*) ' OX => SOX mapping '
           DO IOX = 1, NOBEX_TP
             N = NSOX_FOR_OX(IOX)
             IB = IBSOX_FOR_OX(IOX)
             WRITE(6,*)
             WRITE(6,*) ' ===================='
             WRITE(6,*) ' Info for OX : ', IOX
             WRITE(6,*) ' ===================='
             WRITE(6,*) ' Number of SOXS = ', N
             WRITE(6,*)
             WRITE(6,*) ' SOXS : '
             CALL IWRTMA(ISOX_FOR_OX(IB),1,N,1,N)
           END DO
          END IF
*
          WRITE(6,*) ' ISOX_TO_OX  array '
          CALL IWRTMA(ISOX_TO_OX,1,NSPOBEX_TP,1,NSPOBEX_TP)
*
        END IF
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'OE_SOE')
      RETURN
      END 
* the preliminary copy:
      SUBROUTINE OBEX_TO_SPOBEX2(IFLAG,IOBEX_TP,LCOBEX_TP,LAOBEX_TP,
     &                          NOBEX_TP,ISPOBEX_TP,NSPOBEX_TP,NGAS,
     &                          NOBPT,MS2TOT,MSCOMB_CC,IAAEXC,IACT_SPC,
     &                          IPRDIA,ISOX_TO_OX,MXSPOX,NSOX_FOR_OX,
     &                          IBSOX_FOR_OX,ISOX_FOR_OX,NAEL,NBEL,
     &                          IREFSPC,MN_CREA,MN_ANNI)
*
* Orbital excitation types => Spin-orbital excitations
*
* IFLAG = 1 : Just Number of spinorbital excitations
* IFLAG = 2 : Also the actual excitations
*
*. Jeppe Olsen, Summer of 99
*               Updated with Active-Active excitations, March 2000
*               SOX<=>OX mapping added, Spring 2001
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. Input
      INTEGER IOBEX_TP(2*NGAS, NOBEX_TP), NOBPT(*)
      INTEGER LCOBEX_TP(NOBEX_TP), LAOBEX_TP(NOBEX_TP)
*. Output : Creation in alpha, Creation in beta, 
*           Annihilation in alpha, annihilation in beta
      INTEGER ISPOBEX_TP(4*NGAS,*)
*. Spin orbital type => orbital type
      INTEGER ISOX_TO_OX(*)
*. Number of spinorbital excitations for given orbital excitation 
      INTEGER NSOX_FOR_OX(*)
*. And the actual spinorbital excitations for given orbital excitation 
      INTEGER ISOX_FOR_OX(*)
*. Offset to given orbital type in ISOX_FOR_OX 
      INTEGER IBSOX_FOR_OX(*)
*. Local scratch for occ of reference alpha and beta
      INTEGER IREF_AL(MXPNGAS),IREF_BE(MXPNGAS)
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRDIA)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'OE_SOE')
*. Largest number of creation and annihilation operators
      MX_CREA = IMNMX(LCOBEX_TP,NOBEX_TP,2)
      MX_ANNI = IMNMX(LAOBEX_TP,NOBEX_TP,2)
*. Largest possible block of spincombination
*. Nup + Ndown = MX_OPS, Nup - Ndown = Ms2_TOT
      MX_OPS = MX_CREA + MX_ANNI
      MX_UP = (MX_OPS+MS2TOT)/2
      MX_OPSD2 = MX_OPS/2
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NOBEX_TP = ', NOBEX_TP
        WRITE(6,*) ' MX_CREA, MX_ANNI = ', MX_CREA,MX_ANNI
        WRITE(6,*) ' MX_OPS, MX_UP = ',  MX_OPS, MX_UP
        WRITE(6,*) ' MS2TOT, NAEL, NBEL = ', MS2TOT, NAEL, NBEL
      END IF 
*
C     MX_SOEX_BLK = IBION(MX_OPS,MX_OPSD2)
      MX_SOEX_BLK = IBION(MX_OPS,MX_UP)
      CALL MEMMAN(KLSOEX ,MX_SOEX_BLK*MX_OPS,'ADDL  ',2,'SOEXBL')
*.Alpha and beta-occupations for reference space 
*. Notice : IREFSPC is not used, info in IHPVGAS is used ! 
      IF(IREFSPC.NE.0) THEN
        CALL GET_REF_ALBE_OCC(IREFSPC,IREF_AL,IREF_BE)
      END IF
*. 
      NSPOBEX_TP = 0
      DO ICREA = MN_CREA, MX_CREA
        DO IANNI = MN_ANNI, MX_ANNI
*. Any operators with this number  of creation and annihilation operators ?
          NCOMP = 0
          DO JOBEX_TP = 1, NOBEX_TP
            IF(LCOBEX_TP(JOBEX_TP).EQ.ICREA .AND.
     &         LAOBEX_TP(JOBEX_TP).EQ.IANNI ) NCOMP = NCOMP + 1
          END DO
*
          IF(NCOMP.NE.0) THEN 
*. ICREA creation ops and IANNI annihilation strings with total 
*. spin projection MSTOT
            CALL AB_COMP_FOR_OBOP(ICREA,IANNI,MS2TOT,NABCOMP,
     &           WORK(KLSOEX),NAEL,NBEL,IFLAG)
C?          WRITE(6,*) ' NABCOMP = ', NABCOMP
            DO JOBEX_TP = 1, NOBEX_TP
              IF(LCOBEX_TP(JOBEX_TP).EQ.ICREA .AND.
     &           LAOBEX_TP(JOBEX_TP).EQ.IANNI .AND. 
     &           (MXSPOX.EQ.0.OR.ICREA.LE.MXSPOX)  ) THEN
*. Find the number of the above spincombinations that can be included 
*  for this orbital excitaion 
                 NOPEN = ICREA + IANNI
*
                 CALL ACT_SPOBEX_OBEX(IFLAG,IOBEX_TP(1,JOBEX_TP),
     &                NGAS,NOPEN,NABCOMP,WORK(KLSOEX),NSPCOMP_ACT, 
     &                ISPOBEX_TP(1,NSPOBEX_TP+1),NOBPT,MSCOMB_CC,
     &                IACT_SPC,IAAEXC,IREF_AL,IREF_BE,  
     &                IREFSPC,NAEL,NBEL)
*
                 IF(IFLAG.EQ.2) THEN
C?    WRITE(6,*) ' First SPOBEX_TP (a) '
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)      
                   IOFF = NSPOBEX_TP+1
                   WRITE(6,*) ' IOFF, NSPCOMP_ACT = ',
     &                          IOFF, NSPCOMP_ACT
                   CALL ISETVC(ISOX_TO_OX(IOFF),JOBEX_TP,NSPCOMP_ACT)
C?    WRITE(6,*) ' First SPOBEX_TP (b) '
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)      
                 END IF
*
                 NSPOBEX_TP =  NSPOBEX_TP + NSPCOMP_ACT
               END IF
            END DO
          END IF
*         ^ End if nonvanishing number of operators 
        END DO
      END DO
*     ^ End of loop over creation and annihilation operators 
C?    WRITE(6,*) ' Memchk after construction of spobex '
C?    CALL MEMCHK
C?    WRITE(6,*) ' Memcheck passed '
*
C?    WRITE(6,*) ' First SPOBEX_TP'
C?    CALL IWRTMA(ISPOBEX_TP(1,1),4*NGAS,1,4*NGAS)      
*
      IF(IFLAG.EQ.2) THEN
*. OX => SOX lists
*. Number of sox per ox
        IZERO = 0
        CALL ISETVC(NSOX_FOR_OX,IZERO,NOBEX_TP)
        DO ISOX = 1, NSPOBEX_TP
          IOX = ISOX_TO_OX(ISOX)
          NSOX_FOR_OX(IOX) = NSOX_FOR_OX(IOX)+1
        END DO
*. offsets for sox's with given ox
        DO IOX = 1, NOBEX_TP
         IF(IOX.EQ.1) THEN
           IBSOX_FOR_OX(IOX) = 1
         ELSE
           IBSOX_FOR_OX(IOX) = IBSOX_FOR_OX(IOX-1) + NSOX_FOR_OX(IOX-1)
         END IF
        END DO
*. and the sox for a given ox
        CALL ISETVC(NSOX_FOR_OX,IZERO,NOBEX_TP)
        DO ISOX = 1, NSPOBEX_TP  
          IOX = ISOX_TO_OX(ISOX)
          NSOX_FOR_OX(IOX) = NSOX_FOR_OX(IOX)+1
          IADR = IBSOX_FOR_OX(IOX) -1 + NSOX_FOR_OX(IOX)
          ISOX_FOR_OX(IADR) = ISOX 
        END DO
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ***************************************** '
        WRITE(6,*) ' Information about spinorbital excitations '
        WRITE(6,*) ' ***************************************** '
        WRITE(6,*)
        WRITE(6,*) ' Total number ', NSPOBEX_TP 
        IF(IFLAG.EQ.2) THEN
          DO JSPCOMP_ACT = 1, NSPOBEX_TP  
            WRITE(6,*)
            WRITE(6,*) ' Included spinorbitalexcitation ', JSPCOMP_ACT
            WRITE(6,'(A,16I4)') 
     &      ' Creation of alpha     :', 
     &      (ISPOBEX_TP(I+0*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Creation of beta      :', 
     &      (ISPOBEX_TP(I+1*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Annihilation of alpha :', 
     &      (ISPOBEX_TP(I+2*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Annihilation of beta  :', 
     &      (ISPOBEX_TP(I+3*NGAS,JSPCOMP_ACT),I=1,NGAS)
          END DO
*
          WRITE(6,*) ' OX => SOX mapping '
          DO IOX = 1, NOBEX_TP
            N = NSOX_FOR_OX(IOX)
            IB = IBSOX_FOR_OX(IOX)
            WRITE(6,*)
            WRITE(6,*) ' ===================='
            WRITE(6,*) ' Info for OX : ', IOX
            WRITE(6,*) ' ===================='
            WRITE(6,*) ' Number of SOXS = ', N
            WRITE(6,*)
            WRITE(6,*) ' SOXS : '
            CALL IWRTMA(ISOX_FOR_OX(IB),1,N,1,N)
          END DO
*
        END IF
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'OE_SOE')
      RETURN
      END 
      SUBROUTINE ACT_SPOBEX_OBEX(IFLAG,IOBEX_TP,NGAS,NOPEN,
     &           NSPCOMP_F,ISPCOMP_F,NSPCOMP_ACT,ISPCOMP_ACT,
     &           NOBPT,IMSCOMB,IACT_SPC,IAAEXC,IREF_AL,IREF_BE,
     &           IREFSPC,NAEL,NBEL)
*
* A set of spincomponents, ISCOMP_F  of an orbital excitation 
* is given. Find the number of spinorbital excitations consistent
* with the given orbital excitation IOBEX_TP
*
* Jeppe Olsen, Summer of 99
*
* Spin combinations added Jan 2. 2000 (Sitting in the kitchen,     
*                                      Jette preparing dinner )
*. Active-active excitations added March 2000, In the train
*. Test for excitation * Reference is nonvanishing, April 2001
*
* A spincomponent can be excluded of the following reasons :
* 1 : the number of alpha or beta operators in a given group is
*     larger than the number of orbitals
* 2 : When there are several creation (annihilation) operators 
*     in a given orbital group, each combination of alpha and 
*     beta operators is included only once. 
*     For example if there are three operators belonging to a
*     given orbital space we include : a a a, a a b, a b b, b b b
*     (and not b a a, a b a, b a b, b b a )
* 3 : If IMSCOMB = 1, then only only one of the two combinations 
*     related by spin flip are included 
* 4 : IF IAAEXC = 1 or 2, only selected active crea/anni of   
*     active orbitals are allowed
*     
*
* If IFLAG = 1, then only the number of spinorbital excitation operators is
*               returned, not the actual operators.
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Local scratch
      INTEGER NCA(MXPNGAS),NCB(MXPNGAS),NAA(MXPNGAS),NAB(MXPNGAS)
*.Input
      INTEGER IOBEX_TP(2*NGAS)
      INTEGER ISPCOMP_F(NOPEN,NSPCOMP_F)
      INTEGER NOBPT(*)
*. Alpha and beta occupation of reference 
      INTEGER IREF_AL(*),IREF_BE(*)
*.Output
      INTEGER ISPCOMP_ACT(4*NGAS,*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ACT_SPOEX....: Orbital excitation '
        CALL IWRTMA(IOBEX_TP,1,2*NGAS,1,2*NGAS)
      END IF
C?    WRITE(6,*) ' ACT.., alpha and beta of ref .. '
C?    CALL IWRTMA(IREF_AL,1,NGAS,1,NGAS)
C?    CALL IWRTMA(IREF_BE,1,NGAS,1,NGAS)
C?    WRITE(6,*) '  ACT_SPOBEX_OBEX, IREFSPC = ', IREFSPC 
*
      NSPCOMP_ACT = 0
      IZERO = 0
*
*- Total number of operators in IOBEX_TP
      N_OB_OPER = IELSUM(IOBEX_TP,2*NGAS)
      IF(N_OB_OPER.NE.NOPEN) GOTO 1001
      DO JSPCOMP = 1, NSPCOMP_F
*
        CALL ISETVC(NCA,IZERO,NGAS)
        CALL ISETVC(NCB,IZERO,NGAS)
        CALL ISETVC(NAA,IZERO,NGAS)
        CALL ISETVC(NAB,IZERO,NGAS)
*
        I_AM_OKAY = 1
        IOP = 1  
        DO ICA = 1, 2
          IF(ICA.EQ.1) THEN
            IB = 1
          ELSE
            IB = NGAS + 1
          END IF
*. Check on creation operators
          DO IGAS = 1, NGAS
C?          WRITE(6,*) ' IGAS = ', IGAS
            LOP = NOBPT(IGAS)
C?          WRITE(6,*) ' LOP = ', LOP
            IF(IOBEX_TP(IB-1+IGAS).GT.0) THEN
              NOPFGAS = IOBEX_TP(IB-1+IGAS)
              IBA = 0
              LALPHA = 0
              LBETA = 0
              DO JOP = IOP, IOP + NOPFGAS - 1
C?              WRITE(6,*) ' JOP = ', JOP
                IF(ISPCOMP_F(JOP,JSPCOMP).EQ.1) THEN
                  LALPHA = LALPHA + 1
                  IF(LBETA.NE.0) IBA = 1
                ELSE
                  LBETA = LBETA + 1
                END IF
              END DO
              IF(ICA.EQ.1) THEN
                NCA(IGAS) = LALPHA
                NCB(IGAS) = LBETA
              ELSE
                NAA(IGAS) = LALPHA
                NAB(IGAS) = LBETA
              END IF
              IF(LALPHA.GT.LOP.OR.LBETA.GT.LOP.OR.IBA.EQ.1) THEN
                I_AM_OKAY = 0
              END IF
C             IOP = IOP + LOP
              IOP = IOP + NOPFGAS
            END IF
*           ^ End of there is operators in this space
          END DO
*         ^ End of loop over Gas spaces
        END DO
*       ^ End of loop over crea/anni
*
        I_EXCLUDE = 0
        IF(I_AM_OKAY.EQ.1.AND.IMSCOMB.NE.0) THEN
*. Test with previous included Spincombinations
          DO KSPCOMP = 1, NSPCOMP_ACT
           IDENTICAL = 1
           DO IGAS = 1, NGAS
C?          WRITE(6,*) ' IGAS, KSPCOMP =', IGAS,KSPCOMP
            IF(NCA(IGAS).NE.ISPCOMP_ACT(IGAS+1*NGAS,KSPCOMP))
     &      IDENTICAL = 0
            IF(NCB(IGAS).NE.ISPCOMP_ACT(IGAS+0*NGAS,KSPCOMP))
     &      IDENTICAL = 0
            IF(NAA(IGAS).NE.ISPCOMP_ACT(IGAS+3*NGAS,KSPCOMP))
     &      IDENTICAL = 0
            IF(NAB(IGAS).NE.ISPCOMP_ACT(IGAS+2*NGAS,KSPCOMP))
     &      IDENTICAL = 0
           END DO
           IF(IDENTICAL.EQ.1) I_EXCLUDE = 1
          END DO
        END IF
        IF(IAAEXC.EQ.1.AND.
     &     (NCA(IACT_SPC).NE.0.OR.NAB(IACT_SPC).NE.0))I_EXCLUDE = 1
        IF(IAAEXC.EQ.2.AND.
     &     (NCB(IACT_SPC).NE.0.OR.NAA(IACT_SPC).NE.0))I_EXCLUDE = 1
*. Test that annihilation on reference state is not vanishing
        IF(IREFSPC.NE.0) THEN
*. There is a well defined reference space from which the 
*. excitation are applied, check for vanishing excitations
          DO IGAS = 1, NGAS
            IF(NAA(IGAS).GT.IREF_AL(IGAS) ) I_EXCLUDE = 1
            IF(NAB(IGAS).GT.IREF_BE(IGAS) ) I_EXCLUDE = 1
          END DO
          DO IGAS = 1, NGAS
            IF(NCA(IGAS)-NAA(IGAS)+IREF_AL(IGAS).GT.NOBPT(IGAS))
     &      I_EXCLUDE = 1
            IF(NCB(IGAS)-NAB(IGAS)+IREF_BE(IGAS).GT.NOBPT(IGAS))
     &      I_EXCLUDE = 1
          END DO
        ELSE
*. No well defined reference space, just check that excitation in 
*. principle is feasible
          DO IGAS = 1, NGAS
            IF(NAA(IGAS).GT.MIN(NOBPT(IGAS),NAEL) ) I_EXCLUDE = 1
            IF(NAB(IGAS).GT.MIN(NOBPT(IGAS),NBEL) ) I_EXCLUDE = 1
          END DO
          DO IGAS = 1, NGAS
            IF(NCA(IGAS).GT.NOBPT(IGAS))I_EXCLUDE = 1
            IF(NCB(IGAS).GT.NOBPT(IGAS))I_EXCLUDE = 1
          END DO
        END IF
*
        IF(I_AM_OKAY.EQ.1.AND.I_EXCLUDE.EQ.0) THEN
*. All test passed, welcome to the included 
          NSPCOMP_ACT = NSPCOMP_ACT+1
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Spincombination = ', JSPCOMP
            CALL IWRTMA(ISPCOMP_F(1,JSPCOMP),1,NOPEN,1,NOPEN)
            WRITE(6,*) ' I_AM_OKAY = ', I_AM_OKAY
            WRITE(6,*) ' NAA and NAB '
            CALL IWRTMA(NAA,1,NGAS,1,NGAS)
            CALL IWRTMA(NAB,1,NGAS,1,NGAS)
          END IF
          IF(IFLAG.EQ.2) THEN
            CALL ICOPVE(NCA,ISPCOMP_ACT(1+0*NGAS,NSPCOMP_ACT),NGAS)
            CALL ICOPVE(NCB,ISPCOMP_ACT(1+1*NGAS,NSPCOMP_ACT),NGAS)
            CALL ICOPVE(NAA,ISPCOMP_ACT(1+2*NGAS,NSPCOMP_ACT),NGAS)
            CALL ICOPVE(NAB,ISPCOMP_ACT(1+3*NGAS,NSPCOMP_ACT),NGAS)
          END IF
        END IF
*       ^ End of I_AM_OKAY = 1
      END DO
*     ^ End of loop over spin combinations
*
 1001 CONTINUE
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Generating allowed spincombinations '
        WRITE(6,*) 
     &  ' Input per GASpace : creation and annihilation'
        CALL IWRTMA(IOBEX_TP(1)     ,1,NGAS,1,NGAS)
        CALL IWRTMA(IOBEX_TP(1+NGAS),1,NGAS,1,NGAS)
        WRITE(6,*) 
     &  ' Number of included spincombinations ', NSPCOMP_ACT
        IF(IFLAG.EQ.2) THEN
          DO JSPCOMP_ACT = 1, NSPCOMP_ACT 
            WRITE(6,*)
            WRITE(6,*) ' Included spinorbitalexcitation ', JSPCOMP_ACT
            WRITE(6,'(A,16I4)') 
     &      ' Creation of alpha     :', 
     &      (ISPCOMP_ACT(I+0*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Creation of beta      :', 
     &      (ISPCOMP_ACT(I+1*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Annihilation of alpha :', 
     &      (ISPCOMP_ACT(I+2*NGAS,JSPCOMP_ACT),I=1,NGAS)
            WRITE(6,'(A,16I4)') 
     &      ' Annihilation of beta  :', 
     &      (ISPCOMP_ACT(I+3*NGAS,JSPCOMP_ACT),I=1,NGAS)
          END DO
        END IF
      END IF
*     ^ End if NTEST .ge. 100
*
      RETURN
      END
      SUBROUTINE GET_3BLKS_GCC(KVEC1,KVEC2,KC2,MXCJ)
*
* Allocate the three blocks VEC1, VEC2, C2 used in sigma, densi etc
*
* Version for general coupled cluster
*
* Size of Intermediate arrays are simply obtained from allowed max
* on resolution strings
*
* Jeppe Olsen, Summer of 99
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
*
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
*
      IDUM = 1
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
* 1 : Memory required for CI
*
      CALL GET_L3BLKS(LVEC1,LVEC2,LC2)
      WRITE(6,*) ' Length of vectors for CI  : ', LVEC1, LVEC2, LC2
*
* 2 : Length of vectors used in CC
      MXCJ = MXINKA_CC ** 4
      LBLOCK = MXSOOB_AS
*
      LBLOCK = MAX(LBLOCK,LCSBLK)
      LSCR12 = MAX(LBLOCK,2*MXCJ)  
*
      WRITE(6,*) ' Length of vectors for CC : ', LBLOCK,LSCR12
*
* 3 : Max lengths 
*
      LBLOCK = MAX(LBLOCK,LVEC1)
c      LBLOCK = 2*MAX(LBLOCK,LVEC1)
      LSCR12 = MAX(LSCR12,LC2)
*
      WRITE(6,*) ' Length of allocated buffers : ', LBLOCK,LBLOCK,LSCR12
*
      CALL MEMMAN(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
      CALL MEMMAN(KVEC2,LBLOCK,'ADDS  ',2,'VEC2  ')
      CALL MEMMAN(KC2,LSCR12,'ADDS  ',2,'C2    ')
*
      RETURN
      END
      SUBROUTINE IOFF_SYMBLK_MAT(NSMST,NA,NB,ITOTSM,IOFF,IRESTRICT)
*
* Offset for symmetry blocked matrix. 
* Row-symmetry used as primary index
*
* Jeppe Olsen, Summer of 99
*               July 2000 : IRESTRICT added
*
* IRESTRICT = 0 : No restrictions of symmetry blocks
*           = 1 : Only blocks with ASM .GE. BSM are included
*           =-1 : Only blocks with ASM .LT. BSM are included
      INCLUDE 'implicit.inc'
*. General input  
      INCLUDE 'multd2h.inc'
*. Specific input 
      INTEGER NA(NSMST),NB(NSMST)
*. Output
      INTEGER IOFF(NSMST)
*
      IB = 1
      DO IASM = 1, NSMST
        IBSM = MULTD2H(IASM,ITOTSM)
        INCLUDED = 0
        IF(IRESTRICT.EQ.0) INCLUDED = 1
        IF(IRESTRICT.EQ.1.AND.IASM.GE.IBSM) INCLUDED = 1
        IF(IRESTRICT.EQ.-1.AND.IASM.LT.IBSM) INCLUDED = -1
        IF(INCLUDED.EQ.1) THEN
          IOFF(IASM) = IB
          IB = IB + NA(IASM)*NB(IBSM)
        ELSE
          IOFF(IASM) = 0
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST,IDIAG)
*
* Offsets for symmetryblocks of TCC elements, sym of CA,CB,AA used 
* If Idiag.eq.1 only blocks (Icasm,Iaasm).ge.(Icbsm,Iabsm) are included
*
*
* Jeppe Olsen, Summer of 99
*              July 2000, HNIE : IDIAG added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NCA(*),NCB(*),NAA(*),NAB(*)
*. Output
      INTEGER IBT(8,8,8)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Z_TCC_OFF speaking '
        WRITE(6,*) ' NSMST, IDIAG = ', NSMST, IDIAG
      END IF
*
      IOFF = 1
      DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM_C,ITSYM) 
        DO ISM_CA = 1, NSMST
          ISM_CB = MULTD2H(ISM_C,ISM_CA)
          DO ISM_AA = 1, NSMST
            ISM_AB =  MULTD2H(ISM_A,ISM_AA)
*
            ISM_ALPHA = (ISM_AA-1)*NSMST+ISM_CA
            ISM_BETA  = (ISM_AB-1)*NSMST+ISM_CB
*
            IF((IDIAG.EQ.0).OR.
     &         (IDIAG.EQ.1.AND.ISM_ALPHA.GT.ISM_BETA)) THEN
              IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
              IF(IDIAG.EQ.1.AND.ISM_ALPHA.GT.ISM_BETA) 
     &        IBT(ISM_CB,ISM_CA,ISM_AB) = -IOFF
              IOFF = IOFF + 
     &        NCA(ISM_CA)*NCB(ISM_CB)*NAA(ISM_AA)*NAB(ISM_AB)
            ELSE IF (IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN
              IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
              LEN = NCA(ISM_CA)*NAA(ISM_AA)
              IOFF = IOFF + LEN*(LEN+1)/2
            END IF
*
          END DO
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE C_TO_CKK_MATRIX(C,CKK,IWAY,NJA,NJB,NKA,NKB,NTA,NTB,
     &                    IKJA,XKJA,IKJB,XKJB)
*
* IWAY = 1 : C(Ja,Jb) => C(Ka,Kb,ITa,ITb)
* IWAY = 2 : C(Ka,KB,Ita,Itb) => C(Ja,Jb) + C(Ka,Kb,Ita,Itb)
*
* Jeppe Olsen, Summer of 99
*
* Version where matrix notation is kept to the max
*
* Notice : K and T is in general batched. This is not visible here
*
*. Input
      INCLUDE 'implicit.inc'
      INTEGER    IKJA(NTA,NKA), IKJB(NTB,NKB)             
      DIMENSION  XKJA(NTA,NKA), XKJB(NTB,NKB)             
*. Input/Output 
      DIMENSION C(NJA,NJB),CKK(NKA,NKB,*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from C_TO_CKK '
        WRITE(6,*) '  IWAY =  ', IWAY
      END IF
*
      DO ITB = 1, NTB
      DO ITA = 1, NTA
       ITAB = (ITB-1)*NTA+ITA
       DO KB = 1, NKB         
       DO KA = 1, NKA           
         JA = IKJA(ITA,KA)
         JB = IKJB(ITB,KB)
         SIGN = XKJA(ITA,KA)*XKJB(ITB,KB)
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Sign = ', SIGN
           WRITE(6,'(A,6I4)') ' ITB, ITA, KB, KA, JB, JA',
     &                        ITB, ITA, KB, KA, JB, JA
           WRITE(6,*) ' XKJA, XKJB =', XKJA(ITA,KA),XKJB(ITB,KB)
           WRITE(6,*) ' ITAB = ', ITAB
         END IF
         IF(JA.NE.0.AND.JB.NE.0) THEN
           IF(IWAY.EQ.1) THEN
             CKK(KA,KB,ITAB) = SIGN*C(JA,JB)
           ELSE
             C(JA,JB) = C(JA,JB) + SIGN*CKK(KA,KB,ITAB)
           END IF
         ELSE IF(IWAY.EQ.1) THEN
            CKK(KA,KB,ITAB) = 0.0D0
         END IF
       END DO
       END DO
      END DO
      END DO
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*)
        WRITE(6,*) ' C(Ja,Jb) '
        CALL WRTMAT(C,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' C(Ka,Kb,ITa, Itb) '
        ITAB = 0 
        DO ITB = 1, NTB
        DO ITA = 1, NTA
          ITAB = ITAB + 1
          WRITE(6,*) ' C(Ka,Kb, Ita, Itb) for Ita,Itb =', Ita, Itb
          CALL WRTMAT(CKK(1,1,ITAB),NKA,NKB,NKA,NKB)
        END DO
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE C_TO_CKK(C,CKK,IWAY,NJA,NJB,NKA,NKB,NTA,NTB,
     &                    IKJA,XKJA,IKJB,XKJB)
*
* IWAY = 1 : C(Ja,Jb) => C(Ka,Kb,ITa,ITb)
* IWAY = 2 : C(Ka,KB,Ita,Itb) => C(Ja,Jb) + C(Ka,Kb,Ita,Itb)
*
* Jeppe Olsen, Summer of 99
*
* Version where matrix notation is eliminated
*
* Notice : K and T is in general batched. 
*
*. Input
      INCLUDE 'implicit.inc'
      INTEGER    IKJA(NTA*NKA), IKJB(NTB,NKB)             
      DIMENSION  XKJA(NTA*NKA), XKJB(NTB,NKB)             
*. Input/Output 
C     DIMENSION C(NJA,NJB),CKK(NKA,NKB,*)
      DIMENSION C(NJA*NJB),CKK(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Output from C_TO_CKK, IWAY =  ', IWAY
        WRITE(6,*) ' ======================================='
        WRITE(6,*) ' NJA, NJB, NKA, NKB, NTA, NTB = ',
     &               NJA, NJB, NKA, NKB, NTA, NTB 
      END IF
*
      DO ITB = 1, NTB
      DO ITA = 1, NTA
       ITAB = (ITB-1)*NTA+ITA-1
*
       KAKBITAB = ITAB*NKA*NKB
       IF(IWAY.EQ.1) THEN
*. Iway = 1 part
         DO KB = 1, NKB         
          JB = IKJB(ITB,KB)
C?        WRITE(6,*) ' ITB, KB, JB =', ITB, KB, JB
          IF(JB.NE.0) THEN
           SIGNB = XKJB(ITB,KB)
           JAJB0 = (JB-1)*NJA 
           ITAKA = ITA-NTA
           DO KA = 1, NKA           
            ITAKA = ITAKA + NTA
            JA = IKJA(ITAKA)
C?          WRITE(6,*) ' ITA, KA, JA =', ITA, KA, JA
            KAKBITAB = KAKBITAB + 1
            IF(JA.NE.0) THEN
                JAJB = JAJB0 + JA
                SIGN = XKJA(ITAKA)*SIGNB
                CKK(KAKBITAB) = SIGN*C(JAJB)
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Sign and C(JAJB) ', SIGN, C(JAJB)
                  WRITE(6,*) ' Adress of updated element ', KAKBITAB
                END IF
            ELSE 
               CKK(KAKBITAB) = 0.0D0
            END IF
           END DO
          ELSE IF (JB.EQ.0) THEN
           DO KAKBITAB_EFF =  KAKBITAB+1,  KAKBITAB + NKA  
               CKK(KAKBITAB_EFF) = 0.0D0
           END DO
           KAKBITAB  =  KAKBITAB + NKA
          END IF
*         ^ End if JB .NE. 0
         END DO
       ELSE 
*. Iway = 2 part
         DO KB = 1, NKB         
          JB = IKJB(ITB,KB)
          IF(JB.NE.0) THEN
           JAJB0 = (JB-1)*NJA
           SIGNB = XKJB(ITB,KB)
           ITAKA = ITA-NTA
           DO KA = 1, NKA           
            ITAKA = ITAKA + NTA
            JA = IKJA(ITAKA)
            KAKBITAB = KAKBITAB + 1
            IF(JA.NE.0) THEN
             JAJB = JAJB0 + JA
             SIGN = XKJA(ITAKA)*SIGNB
             C(JAJB) = C(JAJB) + SIGN*CKK(KAKBITAB)
            END IF
           END DO
          ELSE
            KAKBITAB =  KAKBITAB + NKA
          END IF
*         ^ End if JB = 0
         END DO
       END IF
*      ^ End of IWAY switch
      END DO
      END DO
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' C(Ja,Jb) '
        CALL WRTMAT(C,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' C(Ka,Kb,ITa, Itb) '
        ITAB = 0 
        DO ITB = 1, NTB
        DO ITA = 1, NTA
          ITAB = ITAB + 1
          WRITE(6,*) ' C(Ka,Kb, Ita, Itb) for Ita,Itb =', Ita, Itb
          CALL WRTMAT(CKK(1+(ITAB-1)*NKA*NKB),NKA,NKB,NKA,NKB)
        END DO
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE TCC_SUBBLK(TCC,TCC_SUB,IWAY,
     &           NCA,NCA_SUB,ICA_B,NCB,NCB_SUB,ICB_B,
     &           NAA,NAA_SUB,IAA_B,NAB,NAB_SUB,IAB_B,IAB_TRNSP,IDIAG)
*
* Extract - or add- subbblock of an TCC block
*
* IWAY = 1 Obtain TCC_SUB from TCC
* IWAY = 2 Add TCC_SUB to TCC
*
* Jeppe Olsen, Summer of 99
*              Updated July 2000, IAB_TRNSP, IDIAG added
*
* IAB_TRNSP = 1 : Input block is TCC(ICB,ICA,IAB,IAA)
* IDIAG     = 1 : Input block is lower diagonal
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION TCC(NCA*NCB*NAA*NAB)
      DIMENSION TCC_SUB(NCA_SUB,NCB_SUB,NAA_SUB,NAB_SUB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' Input to TCC_SUBBLK '
        WRITE(6,*) ' NCA, NCB, NAA, NAB', NCA,NCB,NAA,NAB
        WRITE(6,*) ' NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB',
     &               NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB
        WRITE(6,*) ' IAB_TRNSP, IDIAG = ', IAB_TRNSP, IDIAG
        WRITE(6,*) ' Iway = ', IWAY
      END IF
*
      IF(IAB_TRNSP.EQ.0.AND.IDIAG.EQ.0) THEN

       DO IAB = 1, NAB_SUB
       DO IAA = 1, NAA_SUB
       DO ICB = 1, NCB_SUB
       DO ICA = 1, NCA_SUB
         IADR_FULL = (IAB+IAB_B-2)*(NAA*NCB*NCA)+(IAA+IAA_B-2)*(NCB*NCA)
     +             + (ICB+ICB_B-2)*NCA + ICA+ICA_B-1
C?       WRITE(6,*) ' ICA, ICB, IAA, IAB =', ICA,ICB,IAA,IAB
         IF(IWAY.EQ.1) THEN
           TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)                                        
         ELSE
           TCC(IADR_FULL)  =  TCC(IADR_FULL) + TCC_SUB(ICA,ICB,IAA,IAB)
         END IF
       END DO
       END DO
       END DO
       END DO
      ELSE IF (IAB_TRNSP.EQ.1.AND.IDIAG.EQ.0) THEN
       DO IAA = 1, NAA_SUB
       DO IAB = 1, NAB_SUB
       DO ICA = 1, NCA_SUB
       DO ICB = 1, NCB_SUB
*. Obtain T(ICA,ICB,IAA,IAB) as T(ICB,ICA,IAB,IAA)
         IADR_FULL = (IAA+IAA_B-2)*NAB*NCA*NCB + (IAB+IAB_B-2)*NCA*NCB
     +             + (ICA+ICA_B-2)*NCB + ICB+ICB_B-1
C?       WRITE(6,*) ' ICA, ICB, IAA, IAB =', ICA,ICB,IAA,IAB
         IF(IWAY.EQ.1) THEN
           TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)                                        
         ELSE
           TCC(IADR_FULL)  =  TCC(IADR_FULL) + TCC_SUB(ICA,ICB,IAA,IAB)
         END IF
       END DO
       END DO
       END DO
       END DO
      ELSE IF (IDIAG.EQ.1) THEN
       DO ICA = 1, NCA_SUB
       DO ICB = 1, NCB_SUB
       DO IAA = 1, NAA_SUB
       DO IAB = 1, NAB_SUB
         ICA_ABS = ICA + ICA_B -1
         ICB_ABS = ICB + ICB_B -1
         IAA_ABS = IAA + IAA_B -1
         IAB_ABS = IAB + IAB_B -1
*
         IF((IAA_ABS.GT.IAB_ABS).OR.
     &      (IAA_ABS.EQ.IAB_ABS.AND.ICA_ABS.GE.ICB_ABS)) THEN
            IADR_FULL = ITDIANUM(ICA_ABS,ICB_ABS,IAA_ABS,IAB_ABS,
     &                  NCA,NAA)
            ILOW = 1
         ELSE 
            IADR_FULL = ITDIANUM(ICB_ABS,ICA_ABS,IAB_ABS,IAA_ABS,
     &                  NCA,NAA)
            ILOW = 0
         END IF
*
         IF(IWAY.EQ.1) THEN
           TCC_SUB(ICA,ICB,IAA,IAB) =  TCC(IADR_FULL)                                        
         ELSE IF(ILOW.EQ.1) THEN
           TCC(IADR_FULL)  =  TCC(IADR_FULL) + TCC_SUB(ICA,ICB,IAA,IAB)
         END IF
       END DO
       END DO
       END DO
       END DO
      END IF
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' TCC_SUBBLK : TCC_SUB '
        NCAB_SUB = NCA_SUB*NCB_SUB
        NAAB_SUB = NAA_SUB*NAB_SUB
        CALL WRTMAT(TCC_SUB,NCAB_SUB,NAAB_SUB,NCAB_SUB,NAAB_SUB)
      END IF
*  
CM    WRITE(6,*) ' Memcheck at end of TCC.... '
CM    CALL MEMCHK
CM    WRITE(6,*) ' Memcheck passed '
*
      RETURN
      END
      SUBROUTINE TCC_SUBBLK_ORIG(TCC,TCC_SUB,IWAY,
     &           NCA,NCA_SUB,ICA_B,NCB,NCB_SUB,ICB_B,
     &           NAA,NAA_SUB,IAA_B,NAB,NAB_SUB,IAB_B)
*
* Extract - or add- subbblock of an TCC block
*
* IWAY = 1 Obtain TCC_SUB from TCC
* IWAY = 2 Add TCC_SUB to TCC
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
*. Input and output
      DIMENSION TCC(NCA,NCB,NAA,NAB)
      DIMENSION TCC_SUB(NCA_SUB,NCB_SUB,NAA_SUB,NAB_SUB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' Output from TCC_SUBBLK '
        WRITE(6,*) ' NCA, NCB, NAA, NAB', NCA,NCB,NAA,NAB
        WRITE(6,*) ' NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB',
     &               NCA_SUB, NCB_SUB, NAA_SUB, NAB_SUB
      END IF
      DO ICA = 1, NCA_SUB
      DO ICB = 1, NCB_SUB
      DO IAA = 1, NAA_SUB
      DO IAB = 1, NAB_SUB
C?      WRITE(6,*) ' ICA, ICB, IAA, IAB =', ICA,ICB,IAA,IAB
        IF(IWAY.EQ.1) THEN
          TCC_SUB(ICA,ICB,IAA,IAB) = 
     &    TCC(ICA+ICA_B-1,ICB+ICB_B-1,IAA+IAA_B-1,IAB+IAB_B-1)
        ELSE
          TCC(ICA+ICA_B-1,ICB+ICB_B-1,IAA+IAA_B-1,IAB+IAB_B-1) = 
     &    TCC(ICA+ICA_B-1,ICB+ICB_B-1,IAA+IAA_B-1,IAB+IAB_B-1) +
     &    TCC_SUB(ICA,ICB,IAA,IAB)
        END IF
      END DO
      END DO
      END DO
      END DO
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' TCC_SUBBLK : TCC_SUB '
        NCAB_SUB = NCA_SUB*NCB_SUB
        NAAB_SUB = NAA_SUB*NAB_SUB
        CALL WRTMAT(TCC_SUB,NCAB_SUB,NAAB_SUB,NCAB_SUB,NAAB_SUB)
*
        WRITE(6,*) ' IWAY = ', IWAY
C?      WRITE(6,*) ' TCC(1,1,1,1) ', TCC(1,1,1,1)
      END IF
*
      RETURN
      END
      SUBROUTINE K_TO_J_TOT(IKJ,XKJ,KSM,IK_B,IK_E,IT,ITSM,NTOP,
     &                      IT_B,IT_E,IM,XM,IBM,NK,LTOP,IZERO_MAP)
*
* Obtain Total map !J> = T-Oper !K>
* T-oper : string of elementary operators 
* !K> String of sym KSM
*
* !K> are of sym KSM and are restricted to IK_B to IK_E ( within sym)
* T_operators are of sym ITSM and are restricted to IT_B to IT_E
*
* Maps for each elementary operator is provided by IM with offset mat IBM
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Input
      INTEGER IT(LTOP,*),IM(*),IBM(8,8,LTOP),NK(8,LTOP)
      DIMENSION XM(*)
*. Output
      INTEGER IKJ(NTOP,(IK_E-IK_B+1))
      DIMENSION XKJ(NTOP,(IK_E-IK_B+1))
*
      NTEST = 00
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Information from K_TO_J.... '
        WRITE(6,*) ' =========================== '
        WRITE(6,*) ' LTOP, NTOP = ', LTOP,NTOP
      END IF
C?    WRITE(6,*) ' Dim of IKJ, XKJ ',NTOP,IK_E-IK_B+1
*
      IZERO_MAP = 1
      IF(LTOP.GT.0) THEN
        DO KSTR = IK_B, IK_E
C?        WRITE(6,*) ' KSTR = ', KSTR
          DO ITOP = IT_B, IT_E
            KNSM = KSM
            SIGN = 1.0D0
            KNSTR = KSTR
C?          WRITE(6,'(A,10I4)') ' T-op : ', (IT(JOP,ITOP),JOP=1,LTOP)
            DO IOP = 1, LTOP
C?            WRITE(6,*) ' Info for IOP = ', IOP
              JOB_ABS = IT(IOP,ITOP)
              JOB_SM = ISMFTO(JOB_ABS)
              JOB_TP = ITPFTO(JOB_ABS)
C?            WRITE(6,*) ' JOB_ABS, IBTS ',
C?   &        JOB_ABS, JOB_TP, JOB_SM, IOBPTS(JOB_TP,JOB_SM)
              JOB_REL = JOB_ABS - IOBPTS(JOB_TP,JOB_SM) + 1
              LK = NK(KNSM,IOP)
C?            WRITE(6,*) ' IBM,JOB_REL,LK,KNSTR',
C?   &        IBM(JOB_SM,KNSM,IOP),JOB_REL,LK,KNSTR
              IADR = IBM(JOB_SM,KNSM,IOP) -1 + (JOB_REL-1)*LK+KNSTR
C?            WRITE(6,*) ' IADR = ', IADR
              KNSTR = IM(IADR)
C?            WRITE(6,*) ' IADR, KNSTR =  ', IADR, KNSTR
              SIGN = SIGN*XM(IADR) 
              IF(KNSTR.EQ.0) GOTO 1001
              KNSM = MULTD2H(KNSM,JOB_SM)
            END DO
*
 1001       CONTINUE
            IF(KNSTR.NE.0)  IZERO_MAP = 0
            IKJ(ITOP-IT_B+1,KSTR-IK_B+1) = KNSTR
            XKJ(ITOP-IT_B+1,KSTR-IK_B+1) = SIGN
C?          WRITE(6,*) ' K_TO_J Row and col : ', 
C?   &      ITOP-IT_B+1,KSTR-IK_B+1
          END DO
        END DO
      ELSE
*. No K-operators => Identity map
        IZERO_MAP = 0
        DO KSTR = IK_B,  IK_E
          IKJ(1,KSTR-IK_B+1) = KSTR
          XKJ(1,KSTR-IK_B+1) = 1.0D0
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Operators '
        CALL IWRTMA(IT(1,IT_B),LTOP,IT_E-IT_B+1,LTOP,IT_E-IT_B+1)
        WRITE(6,*) ' Output from K_TO_J_TOT ' 
        WRITE(6,*) 'IKJ and XKJ arrays '
        LT = IT_E-IT_B+1
        LK = IK_E-IK_B+1
        CALL IWRTMA(IKJ,LT,LK,LT,LK)
        CALL WRTMAT(XKJ,LT,LK,LT,LK)
      END IF
*
      RETURN
      END
* 
      SUBROUTINE REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
* An operatorstring may be specifed as
*
* IOPGAS : Number of operators per GASspace
* IOP    : GASpace of each operator
*
* Transform between these two form
*
* Iway = 1 : IOPGAS => IOP
* Iway = 2 : IOP    => IOPGAS
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
*. Input/Output
      INTEGER IOPGAS(NGAS),IOP(NOP)
*
      IF(IWAY.EQ.1) THEN
        JOP = 0
        DO JGAS = 1, NGAS
          LJGAS = IOPGAS(JGAS)
          DO JJOP = 1, LJGAS
            JOP = JOP + 1
            IOP(JOP) = JGAS
          END DO
        END DO
        NOP = JOP
      ELSE
        DO JGAS = 0, NGAS
          JOP = 0
          DO JJOP = 1, NOP
            IF(IOP(JJOP).EQ.JGAS) JOP = JOP +1 
          END DO
          IOPGAS(JGAS) = JOP
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' IOPGAS => IOP '
        ELSE
         WRITE(6,*) ' IOP => IOPGAS '
        END IF
        WRITE(6,*) ' IOPGAS and IOP '
        CALL IWRTMA(IOPGAS,1,NGAS,1,NGAS)
        CALL IWRTMA(IOP,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
      SUBROUTINE MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SX,NK,IBX,SCLFAC)

* Information for mapping |Kstr> = IOP |Irefoc>
*
* Mapping for each creation/annihilation operator is given
*
*     Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc' 
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'ctcc.inc'
*. Input
      INTEGER IOP(NGAS),IAC(NGAS), IREFOC(NGAS)
*. Output
      INTEGER IX(*),IBX(8,8,NOP)          
      INTEGER NK(8,NOP)
      DIMENSION SX(*)
*. Local Scratch
      INTEGER KOC(MXPNGAS), KGRP(MXPNGAS)
*. IBX(ISM,JSM,JOP) : Start of map in JX for orbitals of symmetry ISM,
*.                    input strings of symmetry JSM for operator JOP
*. Notice at each level N the mapping is from !K(N)> to !K(N-1)> 
*  where !K(0)> is input type and !K(NOP)> is output type
*
C     MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SX,NK,IB,SCLFAC)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' MAP_EXSTR speaking' 
         WRITE(6,*) ' =================='
        WRITE(6,*) ' Reference occupation '
        CALL IWRTMA(IREFOC,1,NGAS,1,NGAS) 
        WRITE(6,*) ' NOP = ', NOP
        WRITE(6,*) ' Operator string, Orbspaces and C/A '
        CALL IWRTMA(IOP,1,NOP,1,NOP)
        CALL IWRTMA(IAC,1,NOP,1,NOP)
      END IF
*
      CALL ICOPVE(IREFOC,KOC,NGAS)
      CALL OCC_TO_GRP(KOC,KGRP,1)
      IOFF = 1
      IF(NTEST.GE.1000) WRITE(6,*) ' NOP = ', NOP
      DO JOP = NOP,1,-1
       IF(NTEST.GE.1000) WRITE(6,*) ' JOP = ', JOP
*
       JGAS = IOP(JOP)
       JAC  = IAC(JOP)
*
       IF(JAC.EQ.1) THEN 
         JAC_ADJ = 2 
       ELSE 
         JAC_ADJ = 1 
       END IF
*. 
       DO JSMOB = 1, NSMOB
        DO JSMSPGP = 1, NSMST
         KSMSPGP = MULTD2H(JSMSPGP,JSMOB)
*
         IBX(JSMOB,KSMSPGP,JOP) = IOFF
         CALL ADAST_GAS(JSMOB,JGAS,NGAS,KGRP,JSMSPGP,
     &        IX(IOFF),SX(IOFF),NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC,
     &        JAC_ADJ)
C     ADAST_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
C    &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC,
C    &                    IAC)
         NK(KSMSPGP,JOP) = NKSTR
         IOFF = IOFF + NOBPTS(JGAS,JSMOB)*NKSTR
C?       WRITE(6,*) ' IOFF, MAXLEN_I1 = ', IOFF, MAXLEN_I1
         IF(IOFF-1.GT. MAXLEN_I1) THEN
           WRITE(6,*) ' MAP... MAXLEN_I1 too small '
           WRITE(6,*) ' IOFF and MAXLEN_I1 = ', IOFF, MAXLEN_I1
           STOP       ' MAP... MAXLEN_I1 too small '
         END IF
*
        END DO
       END DO
*. Updated Kgroup  
*. Info on !K(JOP)> 
       IF(JAC.EQ.1) THEN 
         JDELTA = -1 
       ELSE 
         JDELTA = +1
       END IF
       KOC(JGAS) = KOC(JGAS) + JDELTA
       CALL OCC_TO_GRP(KOC,KGRP,1)
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' KOC and KGRP '
         CALL IWRTMA(KOC,1,NGAS,1,NGAS)
         CALL IWRTMA(KGRP,1,NGAS,1,NGAS)
       END IF
*
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from MAP_EXSTR '
        WRITE(6,*) ' Mappings from a+/a !K(J)> to !K(J-1)> '
        
        DO JOP = 1, NOP
         IF(IAC(JOP).EQ.1) THEN
           WRITE(6,*) 
     &     ' Operator ', JOP, ' is annihilation operator of type',
     &     IOP(NOP) 
         ELSE
           WRITE(6,*) 
     &     ' Operator ', JOP, ' is creation     operator of type',
     &     IOP(NOP) 
         END IF
         DO JSM = 1, NSMOB
           DO KNSM = 1, NSMST
             WRITE(6,*) ' Sym of operator and !K(N)> ',
     &       JSM,KNSM 
             LK = NK(KNSM,JOP)
             WRITE(6,*) ' Number of strings !K(N)> ',LK
             LJ = NOBPTS(IOP(JOP),JSM) 
             DO JJ = 1, LJ
              WRITE(6,*) ' Info for orbital ', JJ
              WRITE(6,*) ' Excited strings and sign '
              IOFF2 = IBX(JSM,KNSM,JOP) + (JJ-1)*LK
C?            WRITE(6,*) ' JOP JSM KNSM LK ',JOP,JSM,KNSM,LK
C?            WRITE(6,*) ' IOFF2 ', IOFF2
              CALL IWRTMA(IX(IOFF2),1,LK,1,LK)        
              CALL WRTMAT(SX(IOFF2),1,LK,1,LK)        
             END DO
           END DO
         END DO
        END DO
      END IF
*
      RETURN
      END 
*
      SUBROUTINE WRT_TP_GENOP(ICA,ICB,IAA,IAB,NGAS)
*
* Write occupation in each GASspace for general operator
*
* Jeppe Olsen, August 1999
*
      INCLUDE 'implicit.inc'
*
      INTEGER ICA(NGAS),ICB(NGAS),IAA(NGAS),IAB(NGAS)
*
      WRITE(6,*) ' Occupation of alpha-creation string '
      CALL IWRTMA(ICA,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of beta -creation string '
      CALL IWRTMA(ICB,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of alpha-annihilation string '
      CALL IWRTMA(IAA,1,NGAS,1,NGAS)
      WRITE(6,*) ' Occupation of beta -annihilation string '
      CALL IWRTMA(IAB,1,NGAS,1,NGAS)
*
      RETURN
      END
      SUBROUTINE GNSIDE(ISD,ICA,ICB,IAA,IAB,
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
* Jeppe Olsen, August 1999    
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
      INTEGER KXAOC(MXPNGAS),KBOC(MXPNGAS)
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
      NTEST = 000 
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
       WRITE(6,*) ' NIA, NIB, NJA, NJB: '
       CALL IWRTMA(NIA,1,NSMST,1,NSMST)
       CALL IWRTMA(NIB,1,NSMST,1,NSMST)
       CALL IWRTMA(NJA,1,NSMST,1,NSMST)
       CALL IWRTMA(NJB,1,NSMST,1,NSMST)
*
       WRITE(6,*) ' Input C vector '
C              WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
       CALL WRTVH1(CB,ICSM,NJA,NJB,NSMST,0)
       WRITE(6,*) ' Input S vector '
       CALL WRTVH1(SB,ISSM,NIA,NIB,NSMST,0)
*
      END IF
*. Number of operators in IA, IB, JA, JB
      NIA_OP = IELSUM(IAOC,NGAS)
      NIB_OP = IELSUM(IBOC,NGAS)
      NJA_OP = IELSUM(JAOC,NGAS)
      NJB_OP = IELSUM(JBOC,NGAS)
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
*. Well, mappings are actually made as annihilation, mappings, so
*. it is the other way around.- aug 09
      NPERM = NAA_OP*NCB_OP + NCA_OP*(NCA_OP-1)/2+NCB_OP*(NCB_OP-1)/2
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNXXX  =  -1.0D0
      ELSE
        SIGNXXX = 1.0D0
      END IF
*. There is another sign term, for general operators-
*. introduced aug 2009!
      NPERM2 = NIB_OP*(NIA_OP+NCA_OP+NAA_OP) 
     &       + NJA_OP*(NIB_OP+NCB_OP+NAB_OP)
      IF(MOD(NPERM2,2).EQ.0) THEN
        SIGNXXX = SIGNXXX
      ELSE
        SIGNXXX = -SIGNXXX
      END IF
*
      IF(I_USE_NEWCCP .EQ. 0) SIGNXXX = 1.0D0
      IF(NTEST.GE.1000) THEN
      WRITE(6,*) ' I_USE_NEWCCP,NPERM2, SIGNXXX = ', 
     &             I_USE_NEWCCP,NPERM2, SIGNXXX
      WRITE(6,*) ' NIA_OP, NIB_OP, NJA_OP, NJB_OP',
     &             NIA_OP, NIB_OP, NJA_OP, NJB_OP
      WRITE(6,*) ' NCA_OP,NCB_OP,NAA_OP,NAB_OP ',
     &             NCA_OP,NCB_OP,NAA_OP,NAB_OP
      WRITE(6,*) ' NIB_OP*(NIA_OP+NCA_OP+NAA_OP) = ',
     &             NIB_OP*(NIA_OP+NCA_OP+NAA_OP)
      WRITE(6,*) ' NJA_OP*(NIB_OP+NCB_OP+NAB_OP) =',
     &             NJA_OP*(NIB_OP+NCB_OP+NAB_OP)
      END IF
*. Type of Ka and Kb
      CALL CCEX_OCC_OCC(JAOC,KXAOC,NGAS,1,IAA,IKA_ZERO)
      CALL CCEX_OCC_OCC(JBOC,KBOC,NGAS,1,IAB,IKB_ZERO)
      IF(IKA_ZERO.EQ.0.AND.IKB_ZERO.EQ.0) THEN
*  NST_SPGP(IOCC,NSTFSM)
       CALL NST_SPGP(KXAOC,NKA)
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
                 CALL WRTMAT2(SIRES,LKAB,LCAB,LKAB,LCAB)
               END IF
             ELSE
               FACTORC = 0.0D0
               CALL MATML7(TBSUB,SIRES,CJRES,
     &         LCAB,LAAB,LKAB,LCAB,LKAB,LAAB,FACTORC,FACTORAB,1) 
               IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' Updated TBSUB '
                 CALL WRTMAT2(TBSUB,LCAB,LAAB,LCAB,LAAB)
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
      SUBROUTINE CCEX_OCC_OCC(INOCC,IOUTOCC,NGAS,IAC,ICCEXC,IZERO)
*
* A set of occupations INOCC is given
*
* Find occupation generated by applying coupled cluster operator ICCEXC
* in either creation or annihilation form
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER INOCC(NGAS)
      INTEGER ICCEXC(NGAS)
*. Output
      INTEGER IOUTOCC(NGAS)
*
      IZERO = 0
      DO JGAS = 1, NGAS
        IF(IAC.EQ.1) THEN
          JDELTA = - ICCEXC(JGAS)
        ELSE
          JDELTA = ICCEXC(JGAS)
        END IF
        IOUTOCC(JGAS) = INOCC(JGAS) + JDELTA
        IF(IOUTOCC(JGAS).LT.0) IZERO = 1
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input occupation '   
        CALL IWRTMA(INOCC,1,NGAS,1,NGAS)
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' String of annihilation operators'
        ELSE
          WRITE(6,*) ' String of creation operators'
        END IF
        CALL IWRTMA(ICCEXC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Output string '
        CALL IWRTMA(IOUTOCC,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE NEWOCC(INOCC,IOPTP,IOPAC,NOP,IOUTOCC)
*
* Occupation of outputstring from Operator * Occupation of input string
*
* Jeppe Olsen, March 1999 
*
      include 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
       INTEGER INOCC(NGAS),IOPAC(NOP),IOPTP(NOP)
*. Output
      INTEGER IOUTOCC(NGAS)
*
      CALL ICOPVE(INOCC,IOUTOCC,NGAS)
*
      DO IOP = 1, NOP
        IOBTP = IOPTP(IOP)
        IF(IOPAC(IOP).EQ.1) THEN
         IDELTA = -1
        ELSE 
         IDELTA = 1
        END IF
        IOUTOCC(IOBTP) = IOUTOCC(IOBTP) + IDELTA
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' NEWOCC : '               
        WRITE(6,*) ' Output and input strings '
        CALL IWRTMA(IOUTOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(INOCC  ,1,NGAS,1,NGAS)
        WRITE(6,*) ' operator : type and a/c' 
        CALL IWRTMA(IOPTP,1,NOP,1,NOP)
        CALL IWRTMA(IOPAC,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
      SUBROUTINE OCC_TO_GRP(IOCC,IGRP,IWAY)
*. Translate between occupation and group labels of supergroup
*
* IWAY = 1 => OCC to group
*      = 2 => Group to occ
*
* Jeppe Olsen, March 1999
*
      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc' 
      INCLUDE 'gasstr.inc'
*. Local scratch
C     INTEGER JGRP(MXPNGAS)
*. Specific input/output
      INTEGER IOCC(NGAS),IGRP(NGAS)
*
      NTEST = 00
*
      IF(IWAY.EQ.1) THEN 
*. Occupation => Group number 
        DO IOBTP = 1, NGAS
          JJGRP = 0
          DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
            IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
          END DO
          IGRP(IOBTP) = JJGRP
*
          IF(JJGRP.EQ.0) THEN
           WRITE(6,*) ' Group not included in list '
           WRITE(6,*) ' GAS space with problem : ', IOBTP
           WRITE(6,*) ' Input occupations : '
           CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
           STOP       ' Group not included in list '
          END IF
        END DO
*
      ELSE
*. Group => Occupation 
        DO IOBTP = 1, NGAS
          IOCC(IOBTP) = NELFGP(IGRP(IOBTP))
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Occupation and corresponding group array '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS) 
      END IF
*
      RETURN
      END 
      SUBROUTINE NST_SPGP(IOCC,NSTFSM)
* Number of strings for given supergroup. 
*.Input supergroup is defined by  occupation in each orb space
*
* Jeppe Olsen , March 99   
* 
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IOCC(NGAS)
*. Output
      DIMENSION NSTFSM(NSMST)        
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc' 
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
*. Scratch 
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER IGRP(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NST_SPGP is speaking '
        WRITE(6,*) ' ======================='
*
        WRITE(6,*) ' Occupation of supergroup '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
      END IF
*. Occupation => groups
C  OCC_TO_GRP_(NSTR,IOCC,IGRP,IWAY)
      CALL OCC_TO_GRP(IOCC,IGRP,1)
      IZERO = 0
      CALL ISETVC(NSTFSM,IZERO,NSMST)
      
*. Max and Min for allowed symmetries
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        DO ISYM = 1, NSMST
          IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
        END DO
        MNSM(IGAS) = NSMST
        DO ISYM = NSMST,1, -1
          IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
        END DO
      END DO
*. Last space with more than one symmetry
      NGASL = 1
      DO IGAS = 1, NGAS
        IF(MXSM(IGAS).NE.MNSM(IGAS)) NGASL = IGAS
      END DO
*. First symmetry combination
      DO IGAS = 1, NGAS
         ISM(IGAS) = MNSM(IGAS)
      END DO
*. Loop over symmetries in each gas space
      IFIRST = 1
 1000 CONTINUE
      IF(IFIRST.EQ.1) THEN
        CALL ISETVC(ISM,1,NGAS)
        IFIRST = 0
        NONEW = 0
      ELSE
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in this supergroup
        ISMSPGP = ISM(1)
        NST = NSTFSMGP(ISM(1),IGRP(1))
        DO JGRP = 2, NGAS
          CALL SYMCOM(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
          ISMSPGP = ISMSPGPO
          NST = NST * NSTFSMGP(ISM(JGRP),IGRP(JGRP))
        END DO
        NSTFSM(ISMSPGP) =   NSTFSM(ISMSPGP) + NST
        GOTO 1000
      END IF         
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
     &  ' Number of strings per symmetry for supergroup'
        CALL IWRTMA10(NSTFSM(1),1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
      SUBROUTINE OPSTR_REST(LOP,IOPTP,IOPAC,IRESTRCT)
*
* An operator string consists top LOP operators, 
* type defined by IOPTP, Anni/Crea by IOPAC
*
*. Find strings of similar type (same orbital space and 
*. anni/crea), and flag that these operators should be restricted.
*
* Jeppe Olsen, March 1999 for the General coupled cluster code.
*
*. Output
*  IRESTRCT : IRESTRCT(I) = 0 => Operator I should not be restricted  
*  IRESTRCT : IRESTRCT(I) = J => Operator I should be restricted  to be less
*                                than index J
*
* Written so the IRESTRCT(I) .GT. I ( or 0) , i.e. the last index can never be
* restricted 
*
* Jeppe Olsen for the General Coupled Cluster Program
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOPTP(LOP),IOPAC(LOP)
*. Output
      INTEGER IRESTRCT(LOP)
*. At the moment I expect that operators belonging to similar
*  types are placed next to each other
*
      DO IOP = 1, LOP-1
        IF(IOPTP(IOP).EQ.IOPTP(IOP+1).AND.
     &     IOPAC(IOP).EQ.IOPAC(IOP+1)     ) THEN
          IRESTRCT(IOP) = IOP+1
        ELSE
          IRESTRCT(IOP) = 0
        END  IF 
        IRESTRCT(LOP) = 0
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Restrictions of operators '
        WRITE(6,*)
        WRITE(6,*) ' TYPE and AC of operators '
        CALL IWRTMA(IOPTP,1,LOP,1,LOP)
        CALL IWRTMA(IOPAC,1,LOP,1,LOP)
        WRITE(6,*)
        WRITE(6,*) ' IRESTRCT array'
        CALL IWRTMA(IRESTRCT,1,LOP,1,LOP)
      END IF
*
      RETURN
      END 
      SUBROUTINE SIG_GCC(C,HC,LUC,LUHC,T)
*
* Outer routine for general sigma
*
* LUHC = \hat T LUC
*
      INCLUDE 'implicit.inc'
      ISIGDEN = 1
      CALL SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
*
      RETURN
      END
      SUBROUTINE SIG_GCC_U(C,HC,LUC,LUHC,LUSCR1,LUSCR2,T,TSCR)
*
* Outer routine for general sigma with antihermitian T, TSCR is scratch
*
* LUHC = \hat T LUC
*
      INCLUDE 'implicit.inc'
      ISIGDEN = 1
      CALL CONJ_CCAMP(T,1,TSCR)
      CALL SIGDEN_CC(C,HC,LUC,LUSCR1,T,ISIGDEN)
      CALL CONJ_T
      CALL SIGDEN_CC(C,HC,LUC,LUSCR2,TSCR,ISIGDEN)
      CALL CONJ_T
      CALL VECSMD(C,HC,1D0,-1D0,LUSCR1,LUSCR2,LUHC,1,-1)
*
      RETURN
      END
      SUBROUTINE DEN_GCC(C,HC,LUC,LUHC,T)
*
* Outer routine for general densi with general coupled cluster
*
* t_i = <LUHC!\hat t_i |LUC>
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
*
      DIMENSION C(*),HC(*),T(*)    
      ISIGDEN = 2
      ZERO = 0.0D0
      
C?    WRITE(6,*) ' N_CC_AMP = ', N_CC_AMP
      CALL SETVEC(T,ZERO,N_CC_AMP)
*
      CALL SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
*
      RETURN
      END
      SUBROUTINE DEN_GCC_S(C,HC,LUC,LUHC,T,TSCR,ISIGN)
*
* Outer routine for general densi with general coupled cluster
*
*     special routine for symmetrized densities
*
*     <LUHC| tau |LUC> +/- <LUC| tau |LUHC>
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
*
      DIMENSION C(*),HC(*),T(*)    
      ISIGDEN = 2
      ZERO = 0.0D0

      CALL SETVEC(TSCR,ZERO,N_CC_AMP)
      CALL SIGDEN_CC(C,HC,LUC,LUHC,TSCR,ISIGDEN)

      CALL SETVEC(T,ZERO,N_CC_AMP)
*
      CALL SIGDEN_CC(C,HC,LUHC,LUC,T,ISIGDEN)

      CALL VECSUM(T,T,TSCR,1d0,dble(ISIGN),N_CC_AMP)
*
      RETURN
      END
      SUBROUTINE SIGDEN_CC(C,HC,LUC,LUHC,T,ISIGDEN)
*
* Outer routine for SIGMA and Density calculation for 
* general CC code
*
* For Sigma T is the input set of coefficients
* For Densi T is the output set of coefficients
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'wrkspc.inc'
*
* =====
*.Input
* =====
*
*.Definition of c and sigma
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
      INCLUDE 'multd2h.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cands.inc'
*
      DIMENSION T(*)
      CALL LUCIAQENTER('SIGDE')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGDEN')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' SIGDEN_CC in action '
        WRITE(6,*) ' ICSM, ISSM   = ', ICSM, ISSM
        WRITE(6,*) ' ICSPC, ISSPC = ', ICSPC, ISSPC
      END IF
*. Info for S-space
      IF(ISSPC.LE.NCMBSPC) THEN
        IATP = 1
        IBTP = 2
      ELSE
        IATP = IALTP_FOR_GAS(ISSPC)
        IBTP = IBETP_FOR_GAS(ISSPC)
C?      WRITE(6,*) ' Mod refspace ', IATP, IBTP
      END IF
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Arrays giving allowed type combinations 
COLD  CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL IAIBCM(ISSPC,WORK(KSIOIO))
*. Arrays giving block type
      KSVST = 1
COLD  CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
*. Arrays for partitioning of sigma  
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
*. Batches  of S vector
      LBLOCK = MAX(MXSOOB_AS,LCSBLK)
C?    WRITE(6,*) ' Info on Batches of Sigma : '
*. Well, SIGDEN uses full symmetry blocks 
      ISIMSYM_LOC = 1
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),0,ISIMSYM_LOC)
*. Number of BLOCKS
      NBLOCK = IFRMR(WORK(KLSI1BT),1,NBATCH)
     &       + IFRMR(WORK(KLSLBT),1,NBATCH) - 1
C?    WRITE(6,*) ' Number of blocks ', NBLOCK
*. Start, length, ...  of  cc ampitudes 
      IOPSM = MULTD2H(ISSM,ICSM)
      CALL IDIM_TCC(WORK(KLSOBEX),NSPOBEX_TP,IOPSM,    
     &     MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,
     &     WORK(KLLSOBEX),WORK(KLIBSOBEX),LEN_T_VEC,
     &     MSCOMB_CC,MX_SBSTR,
     &     WORK(KISOX_FOR_OCCLS),NTOCCLS,WORK(KIBSOX_FOR_OCCLS),
     &     NTCONF,IPRCC)
C    &           MX_TBLK_AS,ISPOX_FOR_CLS,NOCCLS,IBSPOX_FOR_CLS,
C    &           NTCONF,IPRCC)
*. If combinations are in use, renormalize T-coefficients to 
*. spin-orbital normalization 
C     RENORM_T(ITSS_TP,LTSS_TP,NTSS_TP,T,ISM,IWAY)
      IF(ISIGDEN.EQ.1.AND.MSCOMB_CC.EQ.1) 
     &CALL RENORM_T(WORK(KLSOBEX),WORK(KLLSOBEX),NSPOBEX_TP,
     &              T,IOPSM,2)
*
*. Well, SIGDEN uses full symmetry blocks 
      ISIMSYM_LOC = 1
CM    WRITE(6,*) ' Memcheck before call to SIGDEN '
CM    WRITE(6,*) ' Memcheck passed '
      CALL SIGDEN_CC2(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),LUC,LUHC,T,ISIGDEN,
     &     WORK(KVEC3))  
CM    WRITE(6,*) ' Memcheck after call to SIGDEN '
CM    CALL MEMCHK
CM    WRITE(6,*) ' Memcheck passed '
*. Renormalize to combination form
      IF(MSCOMB_CC.EQ.1) 
     &CALL RENORM_T(WORK(KLSOBEX),WORK(KLLSOBEX),NSPOBEX_TP,
     &              T,IOPSM,1)
CM    WRITE(6,*) ' Memcheck after call to RENORM '
CM    CALL MEMCHK
CM    WRITE(6,*) ' Memcheck passed '
*. Eliminate local memory
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'SIGDEN')
*
      CALL LUCIAQEXIT('SIGDE')
*
      RETURN
      END
      SUBROUTINE SIGDEN_CC2(CB,SB,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &           LUC,LUHC,T,ISIGDEN,C2)
*
* First inner routine for general CC sigma/densi
*
* Jeppe Olsen   Summer of 1999
*
* =====
* Input
* =====
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cgas.inc'
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
*
        IF(ISIGDEN.EQ.2) THEN
          WRITE(6,*) ' LHS on LUHC '
          CALL WRTVCD(CB,LUHC,1,-1)
        END IF
      END IF
*
      IF(ISSPC.LE.NCMBSPC) THEN
        IA_OCTP = 1
        IB_OCTP = 2
      ELSE 
        IA_OCTP = IALTP_FOR_GAS(ISSPC)
        IB_OCTP = IBETP_FOR_GAS(ISSPC)
C?      WRITE(6,*) ' Mod refspace ', IA_OCTP, IB_OCTP
      END IF
*
      IB_AOCTP = IBSPGPFTP(IA_OCTP)
      IB_BOCTP = IBSPGPFTP(IB_OCTP)
      NOCTPA = NSPGPFTP(IA_OCTP)
      NOCTPB = NSPGPFTP(IB_OCTP)
*
C?    WRITE(6,*) ' NSSOA and NSSOB '
C?    CALL IWRTMA(WORK(KNSTSO(1)),NSMST,NOCTPA,NSMST,NOCTPA)
C?    CALL IWRTMA(WORK(KNSTSO(2)),NSMST,NOCTPB,NSMST,NOCTPB)
C?    WRITE(6,*) ' NSMST = ', NSMST
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
* ISIGDEN = 1 
* ISIGDEN = 2 : Read in batch of blocks of lhs vector
*. Unique blocks are read in, expanded and normalized to determinant form
        DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
*
          ISATP = IBATS(1,ISBLK)
          ISBTP = IBATS(2,ISBLK)
          ISASM = IBATS(3,ISBLK)
          ISBSM = IBATS(4,ISBLK)
          ISATP_ABS = ISATP + IB_AOCTP - 1
          ISBTP_ABS = ISBTP + IB_BOCTP - 1
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
     &            WORK(KNSTSO(IA_OCTP)),WORK(KNSTSO(IB_OCTP)), 
     &            PSSIGN,IDC,PLSIGN,LUHC,C2,
     &            NSMST,ISCALE,XSCALE)
C     GSTTBLD(CTT,IATP,IASM,IBTP,IBSM,
C    &                  NSASO,NSBSO,PSSIGN,IDC,
C    &                  PLSIGN,LUC,SCR,NSMST,ISCALE,SCLFAC)
          END IF
        END DO
*       ^ End of loop over S-blocks in batch
*. Obtain sigma/density for batch of blocks
        CALL SIGDEN3(LBATS(JBATS),IBATS(1,I1BATS(JBATS)),1,
     &       CB,SB,LUC,T,ISIGDEN,IB_AOCTP,IB_BOCTP,
     &       IA_OCTP,IB_OCTP)
*
C?      WRITE(6,*) ' Elements 2 to 5 after call to SIGDEN3 '
C?      CALL WRTMAT(SB(2),1,4,1,4)

        IF(ISIGDEN.EQ.1) THEN
          IF(IDC.EQ.2) THEN
*. Determinant => combination form and scale
*. reform 
           CALL RFTTS(SB,C2,IBATS(1,IB_SBAT),L_SBAT,
     &                1,NSMST,NOCTPA,NOCTPB,
     &                WORK(KNSTSO(IA_OCTP)), WORK(KNSTSO(IB_OCTP)),
     &                IDC,PSSIGN,1,NTEST)
*. scale
           CALL SCDTTS(SB,IBATS(1,IB_SBAT),L_SBAT,NSMST,NOCTPA,NOCTPB,
     &                WORK(KNSTSO(IA_OCTP)), WORK(KNSTSO(IB_OCTP)),
     &                IDC,1,NTEST)
          END IF
*. Transfer packed S block to permanent storage
          DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
            IOFF = IBATS(6,ISBLK)
            LEN  = IBATS(8,ISBLK)
            CALL ITODS(LEN,1,-1,LUHC)
            CALL TODSC(SB(IOFF),LEN,-1,LUHC)
C?          WRITE(6,*) ' sigma block to disc : ISBLK,IOFF,LEN =',
C?   &                  ISBLK,IOFF,LEN
C?          WRITE(6,*) ' Elements transferrred : '
C?          CALL WRTMAT(SB(IOFF),1,LEN,1,LEN)
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
      SUBROUTINE SIGDEN3(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,T,ISIGDEN,
     &                   IB_AOCTP_S,IB_BOCTP_S,IA_OCTP,IB_OCTP)
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
      INCLUDE 'cands.inc'
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
C     IF(IPRCIX.GE.3)   WRITE(6,*) ' LSCR1 ',LSCR1
*.T-elements                       
CER   INTSCR = MIN(MXINKA_CC**4, MX_TBLK_MX)
      INTSCR = MAX(MXINKA_CC**4, MX_TBLK_MX)
C?    WRITE(6,*) ' MXINKA_CC**4, MX_TBLK_MX = ',
C?   &             MXINKA_CC**4, MX_TBLK_MX
C?    WRITE(6,*) ' INTSCR = ', INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*. Offsets for alpha and beta supergroups in C
      IF(NTEST.GE.100) WRITE(6,*) ' SIGDE3: ICSPC =', ICSPC
      IF(ICSPC.LE.NCMBSPC) THEN
        IATP = 1
        IBTP = 2
      ELSE 
        IATP = IALTP_FOR_GAS(ICSPC)
        IBTP = IBETP_FOR_GAS(ICSPC)
        IF(NTEST.GE.100) WRITE(6,*) ' Modified IATP, IBTP =', IATP, IBTP
      END IF
* It is a bit sloppy to call the occupatio types for C for IATP,IBTP,
* but my time is not for cleaning up pt..
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Arrays giving allowed type combinations '
COLD  CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM(ICSPC,WORK(KCIOIO))
*
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + MXINKA_CC**4
*
*. Arrays for storing NEL consecutive annihilations/creations
C     LSCR3 = MXNSTR*MXTOB*MAX(NAEL,NBEL)
C     LSCR3 = MXNSTR*MXTOB*MX_EXC_LEVEL*NSMST   
C     LSCR3 = MAX_STR_SPGP*MXTOB*MX_EXC_LEVEL
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
C?    WRITE(6,*) ' LENNY = ', LENNY
      CALL MEMMAN(KTOCC1,LENNY,'ADDL  ',2,'TOCC1 ')
      CALL MEMMAN(KTOCC2,LENNY,'ADDL  ',2,'TOCC2 ')
      CALL MEMMAN(KTOCC3,LENNY,'ADDL  ',2,'TOCC3 ')
      CALL MEMMAN(KTOCC4,LENNY,'ADDL  ',2,'TOCC4 ')
*. Arrays for storing NEL consecutive annihilations/creations
      LSCR4 = MXINKA_CC*MXINKA_CC*MXINKA_CC
      CALL MEMMAN(KI1G  ,LSCR4,'ADDL  ',1,'I1G   ')
      CALL MEMMAN(KXI1G,LSCR4,'ADDL  ',2, 'XIG   ')
*
      CALL MEMMAN(KI2G  ,LSCR4,'ADDL  ',1,'I2G   ')
      CALL MEMMAN(KXI2G,LSCR4,'ADDL  ',2, 'XI2G  ')
*
      CALL MEMMAN(KI3G  ,LSCR4,'ADDL  ',1,'I3G   ')
      CALL MEMMAN(KXI3G,LSCR4,'ADDL  ',2, 'XI3G  ')
*
      CALL MEMMAN(KI4G  ,LSCR4,'ADDL  ',1,'I4G   ')
      CALL MEMMAN(KXI4G,LSCR4,'ADDL  ',2, 'XI4G  ')
*
      CALL MEMMAN(KI1GE ,LSCR4,'ADDL  ',1,'I1GE  ')
      CALL MEMMAN(KXI1GE,LSCR4,'ADDL  ',2,'XIGE  ')
*
      CALL MEMMAN(KI2GE  ,LSCR4,'ADDL  ',1,'I2GE  ')
      CALL MEMMAN(KXI2GE,LSCR4,'ADDL  ',2, 'XI2GE ')
*
      CALL MEMMAN(KI3GE  ,LSCR4,'ADDL  ',1,'I3GE  ')
      CALL MEMMAN(KXI3GE,LSCR4,'ADDL  ',2, 'XI3GE ')
*
      CALL MEMMAN(KI4GE  ,LSCR4,'ADDL  ',1,'I4GE  ')
      CALL MEMMAN(KXI4GE,LSCR4,'ADDL  ',2, 'XI4GE ')
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
      CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
C?    WRITE(6,*) ' Info on Batches of C     : '
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
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
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
*. IA_OCTP, IB_OCTP: Main types for sigma
*  IATP, IBTP      : Main types for C     (sic)
C?    WRITE(6,*) ' NSPOBEX_TP before call to SIGDEN4', NSPOBEX_TP 
      CALL SIGDEN4(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KVEC3),
     &     WORK(KNSTSO(IA_OCTP)),WORK(KNSTSO(IB_OCTP)),
     &     WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &     IATP,IBTP,IOCTPA,IOCTPB,
     &     MXINKA_CC,WORK(KINSCR),
     &     WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &     WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),
     &     IPRDIA,LUC,WORK(KCJRES),WORK(KSIRES),
     &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &     ISSM,ICSM,ISIGDEN,NCBATCH,T,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),
     &     WORK(KLSPOBEX_AC),
     &     WORK(KTOCC1),WORK(KTOCC2), WORK(KTOCC3),WORK(KTOCC4),
     &     WORK(KI1G),WORK(KXI1G), WORK(KI2G),WORK(KXI2G),
     &     WORK(KI3G),WORK(KXI3G), WORK(KI4G),WORK(KXI4G),
     &     WORK(KI1GE),WORK(KXI1GE), WORK(KI2GE),WORK(KXI2GE),
     &     WORK(KI3GE),WORK(KXI3GE), WORK(KI4GE),WORK(KXI4GE),
     &     WORK(KLSPOBEX_FRZ),IB_AOCTP_S,IB_BOCTP_S,ICSPC,ISSPC)
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SIGDE3')
      RETURN
      END
      SUBROUTINE SIGDEN4(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           NSSOA,NSSOB,NCSOA,NCSOB,IAGRP,IBGRP,IOCTPA,IOCTPB,
     &           MAXK,
     &           XINT,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           IPRNT,LUC,CJRES,SIRES,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,
     &           ISSM,ICSM,ISIGDEN,NCBATCH,
     &           T,NSPOBEX_TP,ITSPOBEX_TP,IBSPOBEX_TP,
     &           ISPOBEX_AC,
     &           TOCC1,TOCC2,TOCC3,TOCC4,
     &           I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &           I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,XI4GE,
     &           ISPOBEX_FRZ,IB_AOCTP_S,IB_BOCTP_S,ICSPC,ISSPC)
*
*
* Jeppe Olsen, Summer of 1999
*              July 2000, Combinations of T-amplitudes and 
*                         determinants introduced
*              September 2000, ISPOBEX_FRZ added
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
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER NCSOA(NSMST ,*), NCSOB(NSMST ,*)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*)
      DIMENSION  XI1S(*),XI2S(*),XI3S(*)
*
      DIMENSION CJRES(*),SIRES(*)
*. T-coefficients
      DIMENSION T(*), ITSPOBEX_TP(4*NGAS,*), IBSPOBEX_TP(*)
      INTEGER  ISPOBEX_AC(*),ISPOBEX_FRZ(*)
*
*
      CALL LUCIAQENTER('SIGD4')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGDE4')
      IPRNT_ORIG = IPRNT
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' SIGDEN4 speaking :'
        WRITE(6,*) ' ================'
        WRITE(6,*)  
        WRITE(6,*) ' ISIGDEN = ', ISIGDEN
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
        WRITE(6,*) ' LUC = ', LUC
      END IF
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
*
*. Offsets for alpha and beta supergroups in C
      IF(ICSPC.LE.NCMBSPC) THEN
        IA_OCTP_C = 1
        IB_OCTP_C = 2
      ELSE 
C?      WRITE(6,*) ' ICSPC = ', ICSPC
        IA_OCTP_C = IALTP_FOR_GAS(ICSPC)
        IB_OCTP_C = IBETP_FOR_GAS(ICSPC)
        WRITE(6,*) ' Mod space ', IA_OCTP_C ,IB_OCTP_C
      END IF
*
      IB_AOCTP_C = IBSPGPFTP(IA_OCTP_C)
      IB_BOCTP_C = IBSPGPFTP(IB_OCTP_C)
C?    WRITE(6,*) ' IA_OCTP_C, IB_AOCTP_C = ',
C?   &             IA_OCTP_C, IB_AOCTP_C
*
      NSTT_BLK = NSBLOCK/NSMST
* Loop over batches over C blocks      
      REWIND LUC
*
      DO 20000 JCBATCH = 1, NCBATCH             
        IF(NTEST.GE.10000) WRITE(6,*) 'JCBATCH = ', JCBATCH
*
*. Read C blocks into core
*
        NJBLOCK = LCBLOCK(JCBATCH)
        NCTT_BLK = NJBLOCK/NSMST
        IF(NTEST.GE.10000) 
     &  WRITE(6,*) ' NJBLOCK, NCTT_BLK = ', NJBLOCK, NCTT_BLK      
        I1C = I1CBLOCK(JCBATCH)
C       DO JBLOCK = I1C, I1C+NJBLOCK-1          
C         JOFF = ICBLOCK(6,JBLOCK)
C         LEN  = ICBLOCK(8,JBLOCK)
C         CALL IFRMDS(LEN2,1,-1,LUC)
C         CALL FRMDSC(CB(JOFF),LEN,-1,LUC,I_AM_ZERO,I_AM_PACKED)
C       END DO
*
*. Loop over TT blocks of sigma and C in batches and 
*  obtain  contributions 
        DO 9000 ICTT_BLK = 1, NCTT_BLK
          IF(NTEST.GE.10000)
     &    WRITE(6,*) ' 9000, ICTT_BLK = ', ICTT_BLK
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
          IRATP = JATP + IB_AOCTP_C  - 1
          IRBTP = JBTP + IB_BOCTP_C - 1
C?        WRITE(6,*) ' IRATP, JATP, IB_AOCTP_C = ',
C?   &                 IRATP, JATP, IB_AOCTP_C
*
          DO 10000 ISTT_BLK = 1, NSTT_BLK
            IF(NTEST.GE.10000)
     &    WRITE(6,*) ' 10000 ISTT_BLK = ', ISTT_BLK
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
            IF(NTEST.GE.10000)
     &    WRITE(6,*) ' ISPERM = ', ISPERM
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
     &              NSSOA(1,IIATP),NSSOB(1,IIBTP),NSMST,1)
            END IF
*
            IF(IDC.EQ.2.AND.IATP.EQ.IBTP) THEN 
              ISRESTRICT = 1
            ELSE 
              ISRESTRICT = 0
            END IF
            ILATP = IATP + IB_AOCTP_S - 1
            ILBTP = IBTP + IB_BOCTP_S - 1
*. Connections ?
             DO 8000 ITTP = 1, NSPOBEX_TP
              IF(NTEST.GE.10000)
     &    WRITE(6,*) ' 8000, ITTP = ', ITTP
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
               IF(NTEST.GE.10000)
     &    WRITE(6,*) ' IPBLK = ', IPBLK
*
               IF(IPBLK.EQ.1) THEN
                 ITABTRNSP = 0
               ELSE
                 ITABTRNSP = 1
               END IF
*
C?             WRITE(6,*) ' Spinorbitalexcitationtype=',ITTP
C?             CALL WRT_TP_GENOP(ITSPOBEX_TP(ICA_OFF,ITTP),
C?   &                           ITSPOBEX_TP(ICB_OFF,ITTP),
C?   &                           ITSPOBEX_TP(IAA_OFF,ITTP),
C?   &                           ITSPOBEX_TP(IAB_OFF,ITTP),NGAS)
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
               IF(NTEST.GE.10000)
     &    WRITE(6,*) ' ICON, AC, FRZ = ',
     &         ICON, ISPOBEX_AC(ITTP),  ISPOBEX_FRZ(ITTP)
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
                   JATP_ABS = JATP + IB_AOCTP_C - 1 
                   JBTP_ABS = JBTP + IB_BOCTP_C - 1
                   DO JBLK = ICBLK, ICBLK -1 + NCBLK            
                     JASM = ICBLOCK(3,JBLK)
                     JBSM = ICBLOCK(4,JBLK)
                     JOFF = ICBLOCK(5,JBLK)
                     ISCALE = 1
                     CALL GSTTBLD(CB(JOFF),JJATP,JASM,JJBTP,JBSM,
     &                    WORK(KNSTSO(IA_OCTP_C)),
     &                    WORK(KNSTSO(IB_OCTP_C)), 
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
*. Should the block be expanded to full block form
C                  IF(ICDIAG.EQ.1) THEN
C                         CITT_BLK_REFRM(CPCK,CEXP,IDC,ISM,NASTR,NBSTR,
C    &                          IDIAG,NSMST,IWAY,ICOPY)
C                    CALL CITT_BLK_REFRM(CB(ICOFF),C2,ICSM,
C    &                    NSSOA(1,JJATP),NSSOB(1,JJBTP),
C    &                    ICDIAG,NSMST,2,1)
C                  END IF
                 END IF
*
                 IF(ICPERM.EQ.2.AND.ICTRANSPOSED.EQ.0) THEN
*. Transpose 
                    CALL TRP_CITT_BLK(CB(ICOFF),C2,ICSM,
     &                   NSSOA(1,JJATP),NSSOB(1,JJBTP),NSMST,1)
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
                 CALL GNSIDE(ISIGDEN,
     &           ITSPOBEX_TP(LCA_OFF,ITTP),ITSPOBEX_TP(LCB_OFF,ITTP),
     &           ITSPOBEX_TP(LAA_OFF,ITTP),ITSPOBEX_TP(LAB_OFF,ITTP),
     &           NELFSPGP(1,ILATP),NELFSPGP(1,ILBTP),            
     &           NELFSPGP(1,IRATP),NELFSPGP(1,IRBTP),
     &           NSSOA(1,IATP),NSSOB(1,IBTP),
     &           NCSOA(1,JATP),NCSOB(1,JBTP),
     &           T(ITB), 
     &           SB(ISOFF),CB(ICOFF),ISSM,ICSM,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           XINT,CJRES,SIRES,MAXK,
     &           TOCC1,TOCC2,TOCC3,TOCC4,
     &           I1G,XI1G,I2G,XI2G,I3G,XI3G,I4G,XI4G,
     &           I1GE,XI1GE,I2GE,XI2GE,I3GE,XI3GE,I4GE,XI4GE,
     &           IDIAG,ITABTRNSP,ICRESTRICT,ISRESTRICT)
*
                 IF(NTEST.GE.200.AND.ISIGDEN.EQ.2) THEN
                  WRITE(6,*) ' First element of updated density block ',
     &                     T(ITB)
                  WRITE(6,*) ' First 3 elements of T '
                  CALL WRTMAT(T,1,3,1,3)
                 END IF
*
               END IF
*              ^ End if connection 
              END DO 
*             ^ End of loop over permutation of T-block
 8000       CONTINUE
*.          ^  End of loop over TT sigma blocks 
            IF(ISPERM.EQ.2) THEN
*. Back transpose S- blocks
               CALL TRP_CITT_BLK(SB(ISOFF),C2,ISSM,
     &              NSSOB(1,IIBTP),NSSOA(1,IIATP),NSMST,1)
            END IF
 9999     CONTINUE
*         ^ End of loop over S-permutations
          IF(NTEST.GE.200.AND.ISIGDEN.EQ.1) THEN
              WRITE(6,*) ' Updated sigma blocks ' 
              CALL WRTVH1(SB(ISOFF),ISSM,NSSOA(1,IIATP),
     &                    NSSOB(1,IIBTP),NSMST,0) 
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
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SIGDE4')
      CALL LUCIAQEXIT('SIGD4')
      IF(NTEST.GE.10) WRITE(6,*) ' Leaving SIGDEN4 '
      RETURN
      END
      SUBROUTINE GXFSTR(ILA,ILB,IRA,IRB,ICA,IAA,ICB,IAB,NGAS,ICON)
*
* Left type defined by ILA, ILB
* Right type defined by IRA, IRB
* 
* Excitation operator defined by ICA, IAA, ICB,IAB
*
* Does this operator connect the two types ?
*
* ICON = 1 => Connection
* ICON = 0 => No connection
*
*
* Jeppe Olsen, Summer of 99
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ILA(NGAS),ILB(NGAS),IRA(NGAS),IRB(NGAS)
      INTEGER ICA(NGAS),ICB(NGAS),IAA(NGAS),IAB(NGAS)
*
      ICON = 1
      DO IGAS = 1, NGAS
        IF(.NOT.( ICA(IGAS) - IAA(IGAS) + IRA(IGAS).EQ.ILA(IGAS)))
     &  ICON = 0
        IF(.NOT.( ICB(IGAS) - IAB(IGAS) + IRB(IGAS).EQ.ILB(IGAS)))
     &  ICON = 0
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'GXFSTR : IRA, IRB, ILA, ILB '
        CALL IWRTMA(IRA,1,NGAS,1,NGAS)
        CALL IWRTMA(IRB,1,NGAS,1,NGAS)
        CALL IWRTMA(ILA,1,NGAS,1,NGAS)
        CALL IWRTMA(ILB,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' ICA, IAA, ICB, IAB '
        CALL IWRTMA(ICA,1,NGAS,1,NGAS)
        CALL IWRTMA(IAA,1,NGAS,1,NGAS)
        CALL IWRTMA(ICB,1,NGAS,1,NGAS)
        CALL IWRTMA(IAB,1,NGAS,1,NGAS)
        IF(ICON.EQ.1) THEN
          WRITE(6,*) ' We have contact '
        ELSE
          WRITE(6,*) ' No contact '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE ANA_GENCC(T,ISM)
*
* Analyze T-coefficients of general coupled cluster wavefunction
*
* Jeppe Olsen, September 1999 ( modified GASANA)
*
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
*. specific input Coefficients and types
      DIMENSION T(*)
*
* =====
*.Input
* =====
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'ctcc.inc'
      INCLUDE 'cc_exc.inc'
*
      CALL LUCIAQENTER('ANAGC')
      CALL MEMMAN(KLOFF,DUMMY,'MARK  ',DUMMY,'ANAGCC')
*
*. Four blocks of string occupations
      CALL MEMMAN(KLSTR1_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC1')
      CALL MEMMAN(KLSTR2_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC2')
      CALL MEMMAN(KLSTR3_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC3')
      CALL MEMMAN(KLSTR4_OCC,MX_ST_TSOSO_BLK_MX,'ADDL  ',1,'STOCC4')
*. Number of terms to be printed 
      IF(IPRNCIV.EQ.0) THEN
        THRES = 0.1  
        MAXTRM = 200
      ELSE
        THRES = 0.0D0
*. Well atmost 100000 coefs - to save disk ..
        MAXTRM = 100000
      END IF
*
*Numbers and weight of spinorbital excitation type 
      LENGTH_SOX = NSPOBEX_TPE*10
      CALL MEMMAN(KNCPMT_SOX,LENGTH_SOX,'ADDL  ',1,'KNCPMT')
      CALL MEMMAN(KWCPMT_SOX,LENGTH_SOX,'ADDL  ',2,'KWCPMT')
*. Number and weight of orbital excitation type
      LENGTH_OX = (NOBEX_TP+1)*10
      CALL MEMMAN(KNCPMT_OX,LENGTH_OX,'ADDL  ',1,'KNCPMT')
      CALL MEMMAN(KWCPMT_OX,LENGTH_OX,'ADDL  ',2,'KWCPMT')
      IUSLAB = 0
*
*. Occupation of strings of given sym and supergroup
C     WRITE(6,*) ' KISOX_TO_OX before ANA... ', KISOX_TO_OX
      CALL ANA_GENCCS(T,WORK(KLSOBEX),WORK(KLOBEX),
     &            NSPOBEX_TPE,NOBEX_TP+1,ISM,
     &            THRES,MAXTRM,
     &            WORK(KLSTR1_OCC),WORK(KLSTR2_OCC),
     &            WORK(KLSTR3_OCC),WORK(KLSTR4_OCC),
     &            IUSLAB,
     &            IDUMMY,WORK(KNCPMT_SOX),WORK(KWCPMT_SOX),      
     &            WORK(KNCPMT_OX),WORK(KWCPMT_OX),      
     &            NTOOB,IPRNCIV,WORK(KLSOX_TO_OX),MSCOMB_CC)
   
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ANAGCC')
      CALL LUCIAQEXIT('ANAGC')
*
      RETURN
      END
      SUBROUTINE ANA_GENCCS(T,ISPOBEX_TP,IOBEX_TP,
     &                 NSPOBEX_TP,NOBEX_TP,ISM,
     &                 THRES,MAXTRM,
     &                 IOCC_CA, IOCC_CB, IOCC_AA, IOCC_AB,
     &                 IUSLAB,
     &                 IOBLAB,NCPMT_SOX,WCPMT_SOX,
     &                        NCPMT_OX, WCPMT_OX,
     &                 NORB,IPRNCIV,ISOX_TO_OX,MSCOMB_CC)
*
* Analyze T-CC vector :
*
*      1) Print atmost MAXTRM  operators with coefficients
*         larger than THRES
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given 
*         excitation type, spin orbital types and orbital types
*
* Jeppe Olsen , September 1999 
*               July 2000, MSCOMB_CC added
*                                                  

*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
      DIMENSION T(*)
      DIMENSION ISOX_TO_OX(NSPOBEX_TP)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      CHARACTER*6 IOBLAB(*)
*. Output
      DIMENSION NCPMT_SOX(10,NSPOBEX_TP)                           
      DIMENSION WCPMT_SOX(10,NSPOBEX_TP)                          
      DIMENSION NCPMT_OX(10,NOBEX_TP)                           
      DIMENSION WCPMT_OX(10,NOBEX_TP)                          
*
      NTEST = 000
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) '------------------------'
        WRITE(6,*) ' Output from ANA_GENCCS '
        WRITE(6,*) '------------------------'
        WRITE(6,*)
        WRITE(6,*) ' NOBEX_TP = ', NOBEX_TP
        WRITE(6,*) ' ISOX_TO_OX: '
        CALL IWRTMA(ISOX_TO_OX,1,NSPOBEX_TP,1,NSPOBEX_TP)
      END IF
*
      CALL ISETVC(NCPMT_SOX,0    ,10*NSPOBEX_TP)
      CALL SETVEC(WCPMT_SOX,0.0D0,10*NSPOBEX_TP)
      CALL ISETVC(NCPMT_OX,0    ,10*NOBEX_TP)
      CALL SETVEC(WCPMT_OX,0.0D0,10*NOBEX_TP)
      CALL MEMCHK2('AFTZER')
*
* ===========================================================
*.1 : Printout of coefficients and largest occupation vectors       
*.2 : Group coefficients by type and size
* ===========================================================
*
      WRITE(6,*)
      WRITE(6,*) ' Operators are written as : '
      WRITE(6,*)
      WRITE(6,*)   ' Creation of alpha '
      WRITE(6,*)   ' Creation of beta '
      WRITE(6,*)   ' Annihilation of alpha '
      WRITE(6,*)   ' Annihilation of beta '
      WRITE(6,*)
      MINPRT = 0
      ITRM = 0
      ILOOP = 0
      NTVAR = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      TNORM = 0.0D0
      TTNORM = 0.0D0
2001  CONTINUE
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = 1.0D0
        XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNCIV.EQ.1) THEN
*. Print in one shot
       XMAX = 3006.1956
       XMIN =-3006.1956
      END IF
*
      WRITE(6,*)
      WRITE(6,'(A,E10.4,A,E10.4)')
     &'  Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(6,'(A)')
     &'  ========================================================='
      WRITE(6,*)
*
      IT = 0
      IIT = 0
      DO ITSS = 1, NSPOBEX_TP
       IF(NTEST.GE.100) WRITE(6,*) ' Spinorbital  type = ', ITSS
*
       IF(IPRNCIV.EQ.1.OR.NTEST.GE.100) THEN
          WRITE(6,*) ' Type of spin-orbital-excitation '
          CALL WRT_SPOX_TP_JEPPE(ISPOBEX_TP(1,ITSS),1)
       END IF
*. Transform from occupations to groups
       CALL OCC_TO_GRP(ISPOBEX_TP(1+0*NGAS,ITSS),IGRP_CA,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+1*NGAS,ITSS),IGRP_CB,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+2*NGAS,ITSS),IGRP_AA,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+3*NGAS,ITSS),IGRP_AB,1      )
*
       NEL_CA = IELSUM(ISPOBEX_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ISPOBEX_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ISPOBEX_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ISPOBEX_TP(1+3*NGAS,ITSS),NGAS)
*. Diagonal restricted type of spinorbital excitation ?
       IF(MSCOMB_CC.EQ.1) THEN
         CALL DIAG_EXC_CC(
     &        ISPOBEX_TP(1+0*NGAS,ITSS),ISPOBEX_TP(1+1*NGAS,ITSS),
     &        ISPOBEX_TP(1+2*NGAS,ITSS),ISPOBEX_TP(1+3*NGAS,ITSS),
     &        NGAS,IDIAG)
       ELSE
         IDIAG = 0
       END IF
*
*
       DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        DO ISM_CA = 1, NSMST
         ISM_CB = MULTD2H(ISM_C,ISM_CA)
         DO ISM_AA = 1, NSMST
          ISM_AB =  MULTD2H(ISM_A,ISM_AA)
          ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
          ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
*
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.LT.ISM_BETA) GOTO 777
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN 
           IRESTRICT_LOOP = 1
          ELSE
           IRESTRICT_LOOP = 0
          END IF
*. obtain strings
          CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,NSTR_CA,
     &         IOCC_CA, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,NSTR_CB,
     &         IOCC_CB, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,NSTR_AA,
     &         IOCC_AA, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,NSTR_AB,
     &         IOCC_AB, NORB,0,IDUM,IDUM)
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
          IF(NTEST.GE.100)  WRITE(6,*) ' NSTR_AB = ', NSTR_AB
          IF(NTEST.GE.100)  WRITE(6,*) ' NSTR_AA = ', NSTR_AA
          IF(NTEST.GE.100)  WRITE(6,*) ' NSTR_CB = ', NSTR_CB
          IF(NTEST.GE.100)  WRITE(6,*) ' NSTR_CA = ', NSTR_CA
          DO I_AB = 1, NSTR_AB
           IF(IRESTRICT_LOOP.EQ.1) THEN
             I_AA_MIN = I_AB
           ELSE
             I_AA_MIN = 1
           END IF
           DO I_AA = I_AA_MIN, NSTR_AA
            DO I_CB = 1, NSTR_CB
             IF(IRESTRICT_LOOP.EQ.1.AND.I_AA.EQ.I_AB) THEN
               I_CA_MIN = I_CB
             ELSE
               I_CA_MIN = 1
             END IF
             DO I_CA = I_CA_MIN, NSTR_CA
              IT = IT + 1
C?            WRITE(6,*) ' IT, T(IT) = ', IT,T(IT)
*
              IF(ILOOP .EQ. 1 ) THEN 
                NTVAR = NTVAR + 1
                TNORM = TNORM + T(IT)**2
*. Classify element according to size 
                DO IPOT = 1, 10
                  IF(10.0D0 ** (-IPOT+1).GE.ABS(T(IT)).AND.
     &            ABS(T(IT)).GT. 10.0D0 ** ( - IPOT )      )THEN
                    IOEX = ISOX_TO_OX(ITSS)
C?                  WRITE(6,*) ' IOEX, ITSS = ', IOEX, ITSS
                    NCPMT_OX(IPOT,IOEX)= NCPMT_OX(IPOT,IOEX)+1  
                    WCPMT_OX(IPOT,IOEX)= WCPMT_OX(IPOT,IOEX)+T(IT)**2
                  END IF
                END DO
*             ^ End of loop over powers of ten
              END IF
*             ^ end if we are in loop 1
              IF( XMAX .GE. ABS(T(IT)) .AND.
     &        ABS(T(IT)).GT. XMIN .AND. ITRM.LE.MAXTRM) THEN
                ITRM = ITRM + 1
                IIT = IIT + 1
                TTNORM = TTNORM + T(IT) ** 2
                WRITE(6,'(A)')
                WRITE(6,'(A)')
     &          '                 =================== '
                WRITE(6,*)

                WRITE(6,'(A,2I8,2X,E14.8)')
     &          '  Type, number, size of amplitude : ',
     &          ITSS, IT, T(IT)
                IF(IUSLAB.EQ.0) THEN
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CA(IEL+(I_CA-1)*NEL_CA),IEL = 1, NEL_CA)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_CB(IEL+(I_CB-1)*NEL_CB),IEL = 1, NEL_CB)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AA(IEL+(I_AA-1)*NEL_AA),IEL = 1, NEL_AA)
                  WRITE(6,'(4X,10I4)')
     &            (IOCC_AB(IEL+(I_AB-1)*NEL_AB),IEL = 1, NEL_AB)
                END IF
              END IF
*             ^ End if could and should be printed
             END DO
*            ^ End of loop over alpha creation strings
            END DO
*           ^ End of loop over beta creation strings
           END DO
*          ^ End of loop over alpha annihilation 
          END DO 
*         ^ End of loop over beta annihilation 
  777    CONTINUE
         END DO
        END DO
       END DO
*      ^ End of loop over symmetry blocks
      END DO
*     ^ End of loop over over types of excitations
       IF(IIT .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
       IF( XMIN .GT. THRES .AND. ILOOP .LE. 20 ) GOTO 2001
*
       WRITE(6,'(A,E15.8)')
     & '  Norm of printed coefficients .. ', TTNORM
       WRITE(6,'(A,E15.8)')
     & '  Norm of all     coefficients .. ',  TNORM
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') '   Magnitude of T coefficients '
      WRITE(6,'(A)') '  =============================='
      WRITE(6,'(A)')
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, 10
        W = 0.0D0
        N = 0
        DO 290 ITSS = 1, NOBEX_TP
            N = N + NCPMT_OX(IPOT,ITSS)                    
            W = W + WCPMT_OX(IPOT,ITSS)                    
  290   CONTINUE
        WACC = WACC + W
        NACC = NACC + N
        WRITE(6,'(A,I2,A,I2,3X,I9,X,E15.8,3X,E15.8)')
     &  '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
  300 CONTINUE
*
      WRITE(6,*) ' Number of coefficients less than  10-11',
     &           ' IS  ',NTVAR - NACC
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Magnitude of CI coefficients for each type of excitation '
      WRITE(6,'(A)') 
     & '  ========================================================='
      WRITE(6,'(A)')
      DO 400 ITSS   = 1, NOBEX_TP 
          N = 0
          DO 380 IPOT = 1, 10
            N = N + NCPMT_OX(IPOT,ITSS)                     
  380     CONTINUE
          IF(N .NE. 0 ) THEN
            WRITE(6,*) ' Orbital excitation type = ', ITSS
            WRITE(6,'(A,I9)')  
     &      '         Number of coefficients larger than 10-11 ', N
            WRITE(6,*)
            WACC = 0.0D0
            DO IPOT = 1, 10
              N =  NCPMT_OX(IPOT,ITSS)                    
              W =  WCPMT_OX(IPOT,ITSS)                    
              WACC = WACC + W
              WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
            END DO
          END IF 
  400 CONTINUE
*
*. Total weight and number of dets per excitation level
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Total weight and number of SD''s (> 10 ** -11 )  : '          
      WRITE(6,'(A)') 
     & '  ================================================='
      WRITE(6,'(A)')
      WRITE(6,*) '        N      Weight      Acc. Weight   Exc. type '
      WRITE(6,*) ' ==================================================='
      WACC = 0.0D0
      DO ITSS = 1, NOBEX_TP  
          N = 0
          W = 0.0D0
          DO IPOT = 1, 10
            N = N + NCPMT_OX(IPOT,ITSS)                   
            W = W + WCPMT_OX(IPOT,ITSS)                   
          END DO
          WACC = WACC + W
          IF(N .NE. 0 ) THEN
            WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,I3)') 
     &      N,W,WACC,ITSS                                   
          END IF
      END DO
*
      RETURN
      END
      SUBROUTINE EIGGMTN(AMAT,NDIM,ARVAL,AIVAL,ARVEC,AIVEC,
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
      NTEST = 00
      IF(NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Output from EIGGMT '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Real and imaginary parts of eigenvalues '
        CALL WRTMAT_EP(ARVAL,1,NDIM,1,NDIM)
        CALL WRTMAT_EP(AIVAL,1,NDIM,1,NDIM)
      END IF
*
      IF(NTEST.GE.100) THEN
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
      SUBROUTINE COMPRS2LST_B(I1I,XI1I,N1,I2I,XI2I,N2,NKIN,NKOUT,
     &                        I1O,XI1O,I2O,XI2O)
*
* Two lists of excitations/annihilations/creations are given.
* COmpress to common nonvanishing entries
*
* Jeppe Olsen, July 2000 from COMPRS2LST 
*
* Compared to COMPRS2LST : Order of arrys in I1 and I2 interchanged
*                          Output lists differs from input lists
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input lists
      DIMENSION I1I(N1,NKIN),XI1I(N1,NKIN)
      DIMENSION I2I(N2,NKIN),XI2I(N2,NKIN)
*. Output lists
      DIMENSION I1O(N1,NKIN),XI1O(N1,NKIN)
      DIMENSION I2O(N2,NKIN),XI2O(N2,NKIN)
*
      NKOUT = 0
      DO K = 1, NKIN
        I1ACT  = 0
        DO I = 1, N1
          IF(I1I(I,K).NE.0) I1ACT = 1
        END DO
        I2ACT = 0
        DO I = 1, N2
          IF(I2I(I,K).NE.0) I2ACT = 1
        END DO
        IF(I1ACT.EQ.1.AND.I2ACT.EQ.1) THEN
          NKOUT = NKOUT + 1
            DO I = 1, N1
               I1O(I,NKOUT) = I1I(I,K)
              XI1O(I,NKOUT) =XI1I(I,K)
            END DO
            DO I = 1, N2
               I2O(I,NKOUT) = I2I(I,K)
              XI2O(I,NKOUT) =XI2I(I,K)
            END DO
        END IF
      END DO
*
      RETURN
      END
      FUNCTION ITDIANUM(ICA,ICB,IAA,IAB,NC,NA)
*
* Find adress of element T(ICA,ICB,IAA,IAB) in diagonal 
* block
* It is assumed that (ICA,IAA) .ge. (ICB,IAB) has been checked outside
* (i.e. IAA.gt.IAB .or. IAA.eq.IAB.and.ICA.ge.ICB)
*
*
* Jeppe Olsen, July 2000, HNIE
*
      INCLUDE 'implicit.inc'
*. T(1,1,Iab,Iab)
      I11IABIAB = ((IAB-1)*NA - IAB*(IAB-1)/2)*NC*NC+
     &             (IAB-1)*NC*(NC+1)/2
*. T(1,1,Iaa,Iab)
      IF(IAA.GT.IAB) THEN
       I11IAAIAB =  I11IABIAB + NC*(NC+1)/2 + (IAA-IAB-1)*NC*NC
      ELSE
       I11IAAIAB =  I11IABIAB 
      END IF
*. T(Ica,Icb,Iaa,Iab)
      IF(IAA.GT.IAB) THEN
       IADR = I11IAAIAB  + (ICB-1)*NC + ICA
      ELSE
       IADR = I11IAAIAB  + (ICB-1)*NC + ICA - ICB*(ICB-1)/2
      END IF
*
      ITDIANUM = IADR
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Address of diagonal T-element : '
        WRITE(6,*) ' ICA, ICB, IAA, IAB, NC,NA, IADR = ',
     &               ICA, ICB, IAA, IAB, NC,NA, IADR
      END IF
*
      RETURN
      END
C                    CALL GSTTBLD(CB(JOFF),JJATP,JASM,JJBTP,JBSM,
C    &                    WORK(KNSTSO(IA_OCTP)),WORK(KNSTSO(IB_OCTP)),
C    &                    PSSIGN,IDC,PLSIGN,
C    &                    LUC,C2,NSMST,ISCALE,XSCALE)

      SUBROUTINE GSTTBLD(CTT,IATP,IASM,IBTP,IBSM,
     &                  NSASO,NSBSO,PSSIGN,IDC,
     &                  PLSIGN,LUC,SCR,NSMST,ISCALE,SCLFAC)
*
* obtain  determinant block (iatp iasm, ibtp ibsm )
* from vector packed in combination format according to IDC
*
*. If ISCALE = 1, the routine scales and returns the block
*  in determinant normalization, and SCLFAC = 1.0D0
*
* If ISCALE = 0, the routine does not perform any overall
* scaling, and a scale factor is returned in SCLFAC
*
* IF ISCALE = 0, zero blocks are not set explicitly to zero,
* instead  zero is returned in SCLFAC
*
* ISCALE, SCLFAC added May 97
*
* Simplified version working only for vectors on disc
*
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION CTT(*),NSASO(NSMST, *),NSBSO(NSMST, *)
      DIMENSION SCR(*)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN 
        write(6,*) ' GSTTBLD  ,IATP,IASM,IBTP,IBSM,ISCALE'
        write(6,*)            IATP,IASM,IBTP,IBSM,ISCALE     
        WRITE(6,*) ' LUC = ', LUC
      END IF
* =================
* Read in from disc
* =================
      CALL IFRMDS(LBL,1,-1,LUC)
      IF(NTEST.GE.100) write(6,*) ' LBL = ', LBL
      IF(ISCALE.EQ.1) THEN
        CALL FRMDSC(SCR,LBL,-1,LUC,IMZERO,IAMPACK)
      ELSE
        NO_ZEROING = 1
        CALL FRMDSC2(SCR,LBL,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
      END IF
*
      IF(IMZERO.EQ.1.AND.ISCALE.EQ.0) THEN
        SCLFAC = 0.0D0
      ELSE
        NAST = NSASO(IASM,IATP)
        NBST = NSBSO(IBSM,IBTP)
        IF(LBL.NE.0) THEN
          PLSIGN = 1.0D0
          ISGVST = 1
          CALL SDCMRF(CTT,SCR,2,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &         IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
        ELSE
          SCLFAC = 0.0D0
        END IF
      END IF
*
C?    WRITE(6,*) ' ISCALE and SCLFAC on return in GSTTBL',
C?   &ISCALE,SCLFAC
*
      RETURN
      END
      SUBROUTINE OCCLSE(IWAY,NOCCLS,IOCCLS,NEL,ICISPC,
     &                  I_DO_BASSPC,IBASSPC,NOBPT)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in 
* CI space ICISPC
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and 
* IOCCLS = the allowed distributions of electrons
*
* Extended OCCLS : Allows the use of CI spaces obtained 
*                  as combinations of occupation spaces
*
*
*
* Jeppe Olsen, October 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      DIMENSION NOBPT(NGAS)
*. Output
      DIMENSION IOCCLS(NGAS,*)
      DIMENSION IBASSPC(*)
*. Local scratch 
      INTEGER IOCMIN(MXPNGAS),IOCMAX(MXPNGAS)
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' ICISPC = ', ICISPC  
      END IF
*. Largest and smallest accumulated occupations 
      DO IGAS = 1, NGAS
        KMIN = IGSOCCX(IGAS,1,ICMBSPC(1,ICISPC))
        KMAX = IGSOCCX(IGAS,2,ICMBSPC(1,ICISPC))
        DO JSPC = 2, LCMBSPC(ICISPC)
          LSPC = ICMBSPC(JSPC,ICISPC)
          KMIN = MIN(KMIN, IGSOCCX(IGAS,1,LSPC))
          KMAX = MAX(KMAX, IGSOCCX(IGAS,2,LSPC))
        END DO
        IOCMIN(IGAS) = KMIN
        IOCMAX(IGAS) = KMAX
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IOCMIN and IOCMAX from OCCLSE '
        CALL IWRTMA(IOCMIN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCMAX,1,NGAS,1,NGAS)
      END IF
*
      NOCCLS = 0
*. start with smallest allowed number 
      DO IGAS = 1, NGAS
*. Smallest allowed occ in this GASspace
        IOCA(IGAS) =  IOCMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation 
          CALL NXTNUM3(IOCA,NGAS,IOCMIN,IOCMAX,NONEW)
        END IF
        IFIRST = 0
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IOCMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
C?      WRITE(6,*) ' Another accumulated occupation: ' 
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space 
        NEGA=0
        IM_TO_STUFFED = 0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1 
            IF(IOC(IGAS).GT.2*NOBPT(IGAS)) IM_TO_STUFFED = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: ' 
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
*. Correct number of electrons ?
        IEL_TOT = IELSUM(IOC,NGAS)
*. Belongs to one of the accumulated spaces ?
        IN_SPC = 0
        DO JSPC = 1,  LCMBSPC(ICISPC)
          JCISPC = ICMBSPC(JSPC,ICISPC)
C?        WRITE(6,*) ' JSPC, JCISPC = ', JSPC, JCISPC
          INSPC_LOC = 1
          DO IGAS = 1, NGAS 
            IEL = IOCA(IGAS)
            IF(IEL.LT. IGSOCCX(IGAS,1,JCISPC).OR.
     &         IEL.GT. IGSOCCX(IGAS,2,JCISPC)    ) INSPC_LOC = 0
          END DO
C?        WRITE(6,*) ' JSPC, INSPC_LOC =', JSPC, INSPC_LOC 
          IF(INSPC_LOC.EQ.1) IN_SPC = 1
        END DO
*
C?      WRITE(6,*) ' IEL_TOT, NEGA, IN_SPC = ',
C?   &               IEL_TOT, NEGA, IN_SPC 
        IF(IEL_TOT.EQ.NEL.AND.NEGA.EQ.0.AND.IM_TO_STUFFED.EQ.0.AND.
     &     IN_SPC.EQ.1) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Another allowed class : ' 
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
            CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
*
            IF(I_DO_BASSPC.EQ.1) THEN
              IBASSPC(NOCCLS) = IBASSPC_FOR_CLS(IOC)
            END IF
*
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes : '
           WRITE(6,*) ' ===================='
           WRITE(6,*)
           WRITE(6,*) ' Class    Occupation in GASpaces '
           WRITE(6,*) ' ================================'
           DO I = 1, NOCCLS
             WRITE(6,'(1H ,I5,3X,16I3)')
     &       I, (IOCCLS(IGAS,I),IGAS=1, NGAS)
           END DO
C          CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
*
      IF(I_DO_BASSPC.EQ.1) THEN
C       WRITE(6,*) ' Base CI spaces for the classes '
C       CALL IWRTMA(IBASSPC,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END 
      SUBROUTINE DIM_TDL2(NSPOBEX_TP,ISPOBEX_TP,
     &           N_TDL_MAX2,NSMST)          
*
* Find largest dimension of T_DL array by looping over 
* pairs of T-operators.
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cprnt.inc'
*. NO_SX = 1 eliminates single T excitations 
*
*. Input through argument list
*
      INTEGER ISPOBEX_TP(4*NGAS,*)
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
*
      INTEGER ID1OCC_MX(4*MXPNGAS),ID2OCC_MX(4*MXPNGAS) 
      INTEGER IEXOCC_MX(4*MXPNGAS),IT1OCC_MX(4*MXPNGAS) 
      INTEGER IT2OCC_MX(4*MXPNGAS)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIM_TD')  
*
      IZERO = 0
      NTEST = 000
      NTEST = MAX(NTEST,IPRCC)
      N_TDL_MAX2 = 0
*. T(D,L) dim for active Ex operators
      N_TDL_MAXD = 0
*
      NO_SX = 1
*. If NO_SX .ne. 0, single excitations are excluded 
*
*. Number and types of operators in one- and two-electron operators
*
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
      DO IHTP = 1, N1TP+N2TP
*. occupation of IHTP 
        CALL ICOPVE2(WORK(KLHTP),(IHTP-1)*4*NGAS+1,4*NGAS,IOCC_H_AR)
*. indeces 
        CALL ICOPVE2(WORK(KLHINDX),(IHTP-1)*4+1,4,IHINDX)
*. Operators T(I)
        DO ITTP = 1, NSPOBEX_TP
C       IF(ISPOBEX_AC(ITTP).EQ.1) THEN
*. Occupation of T
          CALL ICOPVE(ISPOBEX_TP(1,ITTP),IOCC_T1_AR,4*NGAS)
          NCREA_T = IELSUM(IOCC_T1_AR,2*NGAS)
*. And types in F
          DO IFTP = 1, NSPOBEX_TP
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' IHTP, ITTP, IFTP = ', IHTP, ITTP,IFTP 
*. Occupation of F operator 
            CALL ICOPVE(ISPOBEX_TP(1,IFTP),IOCC_T2_AR,4*NGAS)
*. Number of creations in F
            NCREA_F = IELSUM(IOCC_T2_AR,2*NGAS)
*. Occupation of TF operator
            CALL OP_T_OCC(IOCC_T1_AR,IOCC_T2_AR,IOCC_T12_AR,IMZERO_TF)
*. Occupation of HTF operator
            CALL OP_T_OCC(IOCC_H_AR,IOCC_T12_AR,IOCC_HTF_AR,IMZERO_HTF)
*. And type of HTF
            CALL INUM_FOR_OCC(IOCC_HTF_AR,IHTFTP)
C?          WRITE(6,*) ' NUM for HTF = ', IHTFTP
            ISKIP = 0
            IF((NCREA_T.EQ.1.OR.NCREA_F.EQ.1).AND.NO_SX.EQ.1) ISKIP = 1
            IF((IHTFTP.GT.0.AND.IMZERO_TF.EQ.0.AND.IMZERO_HTF.EQ.0).AND.
     &         ISKIP.EQ.0) THEN
*. Inside bounds, calculate H T F 
*. offsets in TCC arrays
                CALL DIM_TDL_FOR_HTT(IOCC_H_AR,IHINDX,
     &               IOCC_T1_AR,IOCC_T2_AR,NSMST,N_TDL_MAX2,
     &           ID1OCC_MX,ID2OCC_MX,IEXOCC_MX,IT1OCC_MX,IT2OCC_MX,
     &           N_TDL_MAXD)
            END IF
*           ^ End if HTF was inside bounds
          END DO
*         ^ End of loop over F types
C       END IF
*.      ^ End if T-type is active
        END DO
*       ^ End of loop over T types
      END DO
*     ^ End of loop over H types
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Required scratch space for T(D,L) = ',
     &  N_TDL_MAX2
      END IF
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Required scratch space for T(D,L)(active D) = ',
     &  N_TDL_MAXD
        WRITE(6,*) ' Corresponding D1,D2,Ex,I1,I2 : '
        CALL WRT_SPOX_TP(ID1OCC_MX,1)
        CALL WRT_SPOX_TP(ID2OCC_MX,1)
        CALL WRT_SPOX_TP(IEXOCC_MX,1)
        CALL WRT_SPOX_TP(IT1OCC_MX,1)
        CALL WRT_SPOX_TP(IT2OCC_MX,1)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DIM_TD')  
      RETURN
      END
      SUBROUTINE DIM_TDL_FOR_HTT(IH,IHINDX,IT1OCC,IT2OCC,
     &           NSMST,N_TDL_MAXL,
     &           ID1OCC_MX,ID2OCC_MX,IEXOCC_MX,IT1OCC_MX,IT2OCC_MX,
     &           N_TDL_MAXD)
*
* Type of T1, T2 and H are given. Find largest dimension of T(D,L)
*
* 
* Disconnected terms included
*
* Jeppe Olsen, July 2001 from H_T1T2
*
c      INCLUDE 'implicit.inc'
*. Input common blocks
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Input  : Occupations in each GAS space of T1,T2 and H
      INTEGER IT1OCC(NGAS,4),IT2OCC(NGAS,4),IH(NGAS,4)
      INTEGER IHINDX(4)
*. Original indeces of each operator in H
C     INTEGER IEXD1D2_INDX(4)
*. Output : The operator strings connected with the Max resolution dim
      INTEGER ID1OCC_MX(NGAS,4),ID2OCC_MX(NGAS,4),IEXOCC_MX(NGAS,4)
      INTEGER IT1OCC_MX(NGAS,4),IT2OCC_MX(NGAS,4)
*. Local scratch, assuming atmost 4 operators in Hamiltonian 
      INTEGER IHDEEX_CR(4,3), IHEX_CR(4,3)
      INTEGER IHDEEX_AN(4,3), IHEX_AN(4,3)
      INTEGER IK1OCC(MXPNGAS*4),IK2OCC(MXPNGAS*4),IL1OCC(MXPNGAS*4)
      INTEGER IHD1OCC(MXPNGAS*4),IHD2OCC(MXPNGAS*4),IHEXOCC(MXPNGAS*4)
*
      INTEGER ICONT1(4,4), ICONT2(4,4)
      INTEGER NL1(8), ND(8) 
*
      NTEST = 00
*. Total number of operators in H
      NHOP = IELSUM(IH,4*NGAS)
*
      IF(NTEST.GE.50) THEN
        WRITE(6,*)
        WRITE(6,*) ' Occupation of H, T1 and T2 operator '
        WRITE(6,*)
        CALL WRT_SPOX_TP(IH,1)
        CALL WRT_SPOX_TP(IT1OCC,1)
        CALL WRT_SPOX_TP(IT2OCC,1)
        WRITE(6,*) ' Indeces of IH '
        CALL IWRTMA(IHINDX,1,NHOP,1,NHOP)
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
          IF(ICA.EQ.1.AND.IHPVGAS_AB(JOBTP,IAB).EQ.1) JDEEX = 1
*. Annihilation of particle is deexcitation 
          IF(ICA.EQ.2.AND.IHPVGAS_AB(JOBTP,IAB).EQ.2 ) JDEEX = 1
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
         FACX = 1.0D0
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
         IF( IHDEEX_CR(1,1).EQ.IHDEEX_CR(2,1)) FACX = 0.5D0*FACX
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
         IF( IHDEEX_CR(1,1).EQ.IHDEEX_CR(2,1)) FACX = 0.5D0*FACX
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
         IF( IHDEEX_AN(1,1).EQ.IHDEEX_AN(2,1)) FACX = 0.5D0*FACX
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
         IF( IHDEEX_AN(1,1).EQ.IHDEEX_AN(2,1)) FACX = 0.5D0*FACX
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
        END DO
*. Excitation part of Hamilton operator 
        CALL ICOPVE(IH,IHEXOCC,4*NGAS)
        MONE = -1
        IONE =  1
        CALL IVCSUM(IHEXOCC,IHEXOCC,IHD1OCC,IONE,MONE,4*NGAS)
        CALL IVCSUM(IHEXOCC,IHEXOCC,IHD2OCC,IONE,MONE,4*NGAS)
*
        IF(NTEST.GE.50) THEN
          WRITE(6,*) ' D1 contraction part of H in CAAB form '
          CALL WRT_SPOX_TP(IHD1OCC,1)
          WRITE(6,*) ' D2 contraction part of H in CAAB form '
          CALL WRT_SPOX_TP(IHD2OCC,1)
          WRITE(6,*) ' Excitation part of H in CAAB form '
          CALL WRT_SPOX_TP(IHEXOCC,1)
        END IF
* IF(I12FIRST.EQ.1) : 
*.      K1 = D1 T1
*.      L  = EX K1
*.      T(D,L) = T(D2,L)
*.ELse 
*       K1 = D2T2
*       L  = EX K1
*       T(D,L) = D(D1,L)
*
        IF(I12FIRST.EQ.1) THEN
*. Occupation of K1 was determined above and stored in IK1OCC 
*. Determine L 
C T1T2_TO_T12_OCC(I1,I2,I12,NGAS)
          CALL T1T2_TO_T12_OCC(IHEXOCC,IK1OCC,IL1OCC,NGAS)
        ELSE 
*. Occupation was stored above in IK2OCC so
          CALL T1T2_TO_T12_OCC(IHEXOCC,IK2OCC,IL1OCC,NGAS)
        END IF
        IF(NTEST.GE.50) THEN
        WRITE(6,*) ' K1, K2 and L '
        CALL WRT_SPOX_TP(IK1OCC,1)
        CALL WRT_SPOX_TP(IK2OCC,1)
        CALL WRT_SPOX_TP(IL1OCC,1)
        END IF
*. Number of strings per type for K1
C DIM_TCC_OP(ICCOP,NTCC)
        CALL DIM_TCC_OP(IL1OCC,NL1)
*. Number of strings for type for D2
         IF(I12FIRST.EQ.1) THEN
           CALL DIM_TCC_OP(IHD2OCC,ND)
         ELSE
           CALL DIM_TCC_OP(IHD1OCC,ND)
         END IF
*. Largest dimension of L 
         NL_MAX = IMNMX(NL1,NSMST,2)
*. Largest dimension of D2
         ND_MAX = IMNMX(ND,NSMST,2)
         ND_MAX = MIN(ND_MAX,LCCB)
         LEN_TDL = ND_MAX*NL_MAX
         IF(LEN_TDL.GT.N_TDL_MAXL) THEN
           N_TDL_MAXL = LEN_TDL
           IF(NHEX.GT.0) N_TDL_MAXD = LEN_TDL
*. and save D1, D2, Ex, I1, I2
           IF(I12FIRST.EQ.1) THEN
             CALL ICOPVE(IHD1OCC,ID1OCC_MX,4*NGAS)
             CALL ICOPVE(IHD2OCC,ID2OCC_MX,4*NGAS)
             CALL ICOPVE(IHEXOCC,IEXOCC_MX,4*NGAS)
             CALL ICOPVE(IT1OCC,IT1OCC_MX,4*NGAS)
             CALL ICOPVE(IT2OCC,IT2OCC_MX,4*NGAS)
           ELSE 
             CALL ICOPVE(IHD2OCC,ID1OCC_MX,4*NGAS)
             CALL ICOPVE(IHD1OCC,ID2OCC_MX,4*NGAS)
             CALL ICOPVE(IHEXOCC,IEXOCC_MX,4*NGAS)
             CALL ICOPVE(IT2OCC,IT1OCC_MX,4*NGAS)
             CALL ICOPVE(IT1OCC,IT2OCC_MX,4*NGAS)
           END IF
         END IF
*       ^ End of new T(D,L) was the largest obtained
       END IF
*      ^ End if contraction is nonvanishing
      END DO
      END DO
*.    ^ End of loop over possible divisions of crea and anni contractions
*
      RETURN
      END
      SUBROUTINE ANADIFF(C1,C2,INSPC,INSM,LUC1,LUC2)
*
* Analyze differences between two CI expansions residing on LUC1, LUC2.
* The CI vectors are assumed to span the same space 
*
* Jeppe Olsen, September 2002 from GASANA
*
*
c      INCLUDE 'implicit.inc'
* =====
*.Input
* =====
*
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'crun.inc'
*. Local scratch
      DIMENSION C1(*),C2(*)
*
      CALL LUCIAQENTER('ANADI')
      CALL MEMMAN(KLOFF,DUMMY,'MARK  ',DUMMY,'ANADIF')
*
** Specifications of internal space
*
      NTEST = 10
      NTEST = MAX(NTEST,IPRDIA)
* Type of alpha and beta strings
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
*. Arrays giving allowed type combinations 
COLD  CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL IAIBCM(INSPC,WORK(KSIOIO))
*. Arrays giving block type
      KSVST = 1
COLD  CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
*. Assume total symmetric op 
      CALL ZBLTP(ISMOST(1,INSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
*. Arrays for partitioning of sigma  
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
*. Batches  of S vector
      LBLOCK =  MXSOOB
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KSIOIO),ISMOST(1,INSM),
     &     NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),0,0)
*. Number of BLOCKS
      NBLOCK = IFRMR(WORK(KLSI1BT),1,NBATCH)
     &       + IFRMR(WORK(KLSLBT),1,NBATCH) - 1
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' ANADIFF speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL '
        WRITE(6,*)   IATP,IBTP,NAEL,NBEL
        WRITE(6,*) ' NOCTPA NOCTPB ', NOCTPA,NOCTPB 
        WRITE(6,*) ' INSPC,NBLOCK = ', INSPC, NBLOCK
      END IF
*
      IF( ICISTR .NE. 1 ) THEN 
         CALL REWINO(LUC1)
         CALL REWINO(LUC2)
      END IF
*. Number of occupation classes 
      IWAY = 1
      NEL = NAEL + NBEL
      CALL OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),
     &           0,0,NOBPT)
*. and then the occupation classes 
      CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL')
      IWAY = 2
      CALL OCCLS(IWAY,NOCCLS,WORK(KLOCCLS),NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),
     &           0,0,NOBPT)
*
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
*
      LENGTH = NOCCLS*10
      CALL MEMMAN(KNCPMT,LENGTH,'ADDL  ',1,'KNCPMT')
      CALL MEMMAN(KWCPMT,LENGTH,'ADDL  ',2,'KWCPMT')
*
      CALL MEMMAN(KLDIFF,NOCCLS,'ADDL  ',2,'KLDIFF')
      CALL MEMMAN(KLOVLP,NOCCLS,'ADDL  ',2,'KLOVLP')
      CALL MEMMAN(KLNORM,NOCCLS,'ADDL  ',2,'KLNORM')
*
*. Occupation of strings of given sym and supergroup
      CALL ANADIFFS(C1,C2,LUC1,LUC2,
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &            NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &            NBLOCK,WORK(KLSIBT),
     &            NAEL,NBEL,
     &            WORK(KLASTR),WORK(KLBSTR),
     &            WORK(KSBLTP),NSMST,
     &            WORK(KNCPMT),WORK(KWCPMT),NELFSPGP,      
     &            NOCCLS,NGAS,WORK(KLOCCLS),ICISTR,NTOOB,
     &            WORK(KLDIFF),WORK(KLOVLP),WORK(KLNORM))
   
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'ANADIF')
      CALL LUCIAQEXIT('ANADI')
*
      RETURN
      END
      SUBROUTINE ANADIFFS(C1,C2,LUC1,LUC2,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                 MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                 NAEL,NBEL,
     &                 IASTR,IBSTR,IBLTP,NSMST,
     &                 NCPMT,WCPMT,NELFSPGP,NOCCLS,NGAS,
     &                 IOCCLS,ICISTR,NORB,
     &                 DIFF,OVLAP,XNORM)
*
* Analyze Difference between two vectors spanning the same CI space :
*
*      1) Number of differences in a given range
*
*      2) Number of differences in given range for each
*         occupation class         
*      3) Norm of difference vector for each occ class
*         Overlap of vectors for each occ class
*
* Jeppe Olsen , Sept 2002 from GASANA   
*               For analyzing differences between CI and CC vectors
*                                                  
      INCLUDE 'implicit.inc'
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)  
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IBLTP(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
*. Output
      DIMENSION NCPMT(10,NOCCLS)                           
      DIMENSION WCPMT(10,NOCCLS)                          
      DIMENSION DIFF(NOCCLS),OVLAP(NOCCLS),XNORM(NOCCLS)
*. Scratch
      DIMENSION C1(*), C2(*)

*
      WRITE(6,*)
      WRITE(6,*) ' ===================================================='
      WRITE(6,*) ' The difference between two vectors will be analyzed '
      WRITE(6,*) ' ===================================================='
      WRITE(6,*)
*
*.Size of differences for the various occ classes 
*
      ZERO = 0.0D0
      IZERO = 0
      CALL ISETVC(NCPMT,IZERO,10*NOCCLS)
      CALL SETVEC(WCPMT,ZERO, 10*NOCCLS)
      CALL SETVEC(DIFF,ZERO,NOCCLS)
      CALL SETVEC(XNORM,ZERO,NOCCLS)
      CALL SETVEC(OVLAP,ZERO,NOCCLS)
*
      IDET = 0
      NCIVAR = 0
      IF(ICISTR .GE. 2 ) THEN
        CALL REWINO(LUC1)
        CALL REWINO(LUC2)
      END IF
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
*. Occupation class corresponding to given occupation
        JOCCLS = 0
        DO JJOCCLS = 1, NOCCLS
          IM_THE_ONE = 1
          DO IGAS = 1, NGAS
            IF(NELFSPGP(IGAS,IATP-1+IOCTPA)+
     &         NELFSPGP(IGAS,IBTP-1+IOCTPB).NE.IOCCLS(IGAS,JJOCCLS))
     &         IM_THE_ONE = 0
          END DO
          IF(IM_THE_ONE .EQ. 1 ) JOCCLS = JJOCCLS
        END DO
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        IBBAS = 1
        IABAS = 1
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN 
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(IDET,1,-1,LUC1)
          CALL FRMDSC(C1,IDET,-1,LUC1,IMZERO,IAMPACK)
          CALL IFRMDS(IDET,1,-1,LUC2)
          CALL FRMDSC(C2,IDET,-1,LUC2,IMZERO,IAMPACK)
          IDET = 0
        END IF
        DO IB = IBBAS,IBBAS+NIB-1
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB - IBBAS + IABAS
          ELSE
            MINIA = IABAS
          END IF
          DO IA = MINIA,IABAS+NIA-1
            IDET = IDET + 1
            NCIVAR = NCIVAR + 1
*
            OVLAP(JOCCLS) = OVLAP(JOCCLS) + C1(IDET)*C2(IDET)
            DIFF(JOCCLS)  = DIFF(JOCCLS)  + (C1(IDET)-C2(IDET))**2
            XNORM(JOCCLS) = XNORM(JOCCLS) + C1(IDET)*C1(IDET)
*
            DIFFL = ABS(C1(IDET)-C2(IDET))
            DO IPOT = 1, 10
              IF(10.0D0 ** (-IPOT+1).GE.DIFFL.AND.
     &           DIFFL.GT. 10.0D0 ** ( - IPOT )) THEN
                 NCPMT(IPOT,JOCCLS)= NCPMT(IPOT,JOCCLS)+ 1  
                 WCPMT(IPOT,JOCCLS)= WCPMT(IPOT,JOCCLS)+ DIFFL**2
              END IF
            END DO
*           ^ End of loop over powers of ten
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
      END DO
*     ^ End of loop over blocks
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') '   Magnitude of differences '
      WRITE(6,'(A)') '  ==========================='
      WRITE(6,'(A)')
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, 10
        W = 0.0D0
        N = 0
        DO 290 JOCCLS = 1, NOCCLS 
            N = N + NCPMT(IPOT,JOCCLS)                    
            W = W + WCPMT(IPOT,JOCCLS)                    
  290   CONTINUE
        WACC = WACC + W
        NACC = NACC + N
        WRITE(6,'(A,I2,A,I2,3X,I9,X,E15.8,3X,E15.8)')
     &  '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
  300 CONTINUE
*
      WRITE(6,*) ' Number of coefficients less than  10-11',
     &           ' IS  ',NCIVAR - NACC
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Differences of CI coefficients for each excitation level '
      WRITE(6,'(A)') 
     & '  ========================================================='
      WRITE(6,'(A)')
      DO 400 JOCCLS = 1, NOCCLS  
          N = 0
          DO 380 IPOT = 1, 10
            N = N + NCPMT(IPOT,JOCCLS)                     
  380     CONTINUE
          IF(N .NE. 0 ) THEN
            WRITE(6,*)
            WRITE(6,'(A,15I3)')'       Occupation of active sets :',
     &      (IOCCLS(IGAS,JOCCLS),IGAS=1, NGAS)
            WRITE(6,'(A,I9)')  
     &      '         Number of coefficients larger than 10-11 ', N
            WRITE(6,*)
            WACC = 0.0D0
            DO 370 IPOT = 1, 10
              N =  NCPMT(IPOT,JOCCLS)                    
              W =  WCPMT(IPOT,JOCCLS)                    
              WACC = WACC + W
              WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
  370       CONTINUE
          END IF 
  400 CONTINUE
*
*. Total weight and number of dets per excitation level
*
      WRITE(6,*) ' Test of OVLAP '
      CALL WRTMAT_EP(OVLAP,1,NOCCLS,1,NOCCLS)
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Norm of difference and overlap, absolute '
      WRITE(6,'(A)') 
     & '  ================================================='
      WRITE(6,'(A)')
      WRITE(6,*) '     Difference   Overlap    Sqnorm     Occupation '
      WRITE(6,*) '      ||1>-|2>|     <1|2>     <1|1>   '
      WRITE(6,*) ' ==================================================='
      DO JOCCLS = 1, NOCCLS
          WRITE(6,'(6X,E10.4,2X,E10.4,1X,E10.5,2X,16(1X,I2))') 
     &    SQRT(DIFF(JOCCLS)),OVLAP(JOCCLS),XNORM(JOCCLS),
     &    (IOCCLS(IGAS,JOCCLS),IGAS=1,NGAS)
      END DO
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Norm of difference and overlap, relative '
      WRITE(6,'(A)') 
     & '  ================================================='
      WRITE(6,'(A)')
      WRITE(6,*) '       Difference       Overlap     Occupation    '
      WRITE(6,*) '   ||1>-|2>|/||1>|    <1|2>/<1|1>      '
      WRITE(6,*) ' ==================================================='
      DO JOCCLS = 1, NOCCLS
          IF(XNORM(JOCCLS).EQ.0.0D0) THEN
            WRITE(6,*) ' Zero norm => no relative numbers '
          ELSE
            WRITE(6,'(6X,E10.4,6X,E10.4,2X,16(1X,I2))') 
     &      SQRT(DIFF(JOCCLS))/SQRT(XNORM(JOCCLS)),
     &      OVLAP(JOCCLS)/XNORM(JOCCLS),
     &      (IOCCLS(IGAS,JOCCLS),IGAS=1,NGAS)
          END IF
      END DO
*
      WRITE(6,*) ' Info for total wavefunction : '
      WRITE(6,*) ' ============================= '
      WRITE(6,*)
*. norm of |1> - |2>
      XDIFFT_SQ = ELSUM(DIFF,NOCCLS)
      XDIFFT  = SQRT(XDIFFT_SQ)
*. Overlap <1|2>
      OVLAPT = ELSUM(OVLAP,NOCCLS)
*
      WRITE(6,*) '     Norm of difference |1> - |2> = ' , XDIFFT 
      WRITE(6,*) '     Overlap <1|2>                = ', OVLAPT
*
      RETURN
      END
      FUNCTION ELSUM(VEC,NELMNT)
*
* Sum elements of a vector
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION VEC(*)
*
      X = 0.0D0
      DO I = 1, NELMNT
        X = X + VEC(I) 
      END DO
*
      ELSUM = X
*
      RETURN
      END 
*----------------------------------------------------------------------*
      SUBROUTINE SET_FAD(NFRZ,NACT,NDEL)
*----------------------------------------------------------------------*
*     
*     work-around: set up arrays with info on frozen, active and
*                  deleted orbitals (until ORBINF works properly)
*
*----------------------------------------------------------------------*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'

      CALL ISETVC(NFRZ,0,NSMOB)
      CALL ISETVC(NACT,0,NSMOB)
      CALL ISETVC(NDEL,0,NSMOB)

      DO IGAS = 1, NGAS
        IF (I_IADX(IGAS).EQ.1)
     &       CALL IVCSUM(NFRZ,NFRZ,NGSOB(1,IGAS),1,1,NSMOB)
        IF (I_IADX(IGAS).EQ.2)
     &       CALL IVCSUM(NACT,NACT,NGSOB(1,IGAS),1,1,NSMOB)
        IF (I_IADX(IGAS).EQ.3)
     &       CALL IVCSUM(NDEL,NDEL,NGSOB(1,IGAS),1,1,NSMOB)
      END DO

      RETURN
      END
*----------------------------------------------------------------------*
      SUBROUTINE CCWFLABEL(CCFORM,I_OOCC,I_BCC,I_OBCC,CCLABEL)
*----------------------------------------------------------------------*
*     based on the entry in CCFROM, return a descriptive label for
*     more common use
*----------------------------------------------------------------------*
      IMPLICIT NONE

      INTEGER, PARAMETER :: NTEST = 00

      CHARACTER, INTENT(IN) ::
     &     CCFORM*6
      INTEGER, INTENT(IN) ::
     &     I_OOCC, I_BCC, I_OBCC
      CHARACTER, INTENT(OUT) ::
     &     CCLABEL*80

      IF (CCFORM(1:5).EQ.'TCC  '.AND.I_OOCC.EQ.1) THEN
        CCLABEL = 'Orbital-Optimized Coupled-Cluster'
      ELSE IF (CCFORM(1:5).EQ.'TCC  '.AND.I_BCC.EQ.1) THEN
        CCLABEL = 'Brueckner Coupled-Cluster'
      ELSE IF (CCFORM(1:5).EQ.'TCC  '.AND.I_OBCC.EQ.1) THEN
        CCLABEL = 'Orbital-Optimized/Brueckner Coupled-Cluster'
      ELSE IF (CCFORM(1:5).EQ.'TCC  ') THEN
        CCLABEL = 'Coupled-Cluster'
      ELSE IF (CCFORM(1:5).EQ.'VCC  ') THEN
        CCLABEL = 'Variational Coupled-Cluster'
      ELSE IF (CCFORM(1:5).EQ.'UCC  ') THEN
        CCLABEL = 'Unitary Coupled-Cluster'
      ELSE IF (CCFORM(1:5).EQ.'ECC  ') THEN
        CCLABEL = 'Extended Coupled-Cluster'
      ELSE IF (CCFORM(1:5).EQ.'ECC2 ') THEN
        CCLABEL = 'Quadratic Extended Coupled-Cluster'
      ELSE
        CCLABEL = 'A very special Coupled-Cluster model '//
     &            ' (see LUCIA output)'
      END IF

      IF (NTEST.GE.100) THEN
        WRITE(6,*) '========================='
        WRITE(6,*) 'Debug info from CCWFLABEL'
        WRITE(6,*) '========================='
        WRITE(6,*) ' CCFORM:  "',CCFORM(1:5),'"'
        WRITE(6,*) ' CCLABEL: "',CCLABEL(1:LEN_TRIM(CCLABEL)),'"'
      END IF

      RETURN
      END
*----------------------------------------------------------------------*
      SUBROUTINE WR_IR2DEN_BLK(IMOD,DEN,BUFF,LUIR2DEN,IREW,
     &                         NTOOB,NSMOB,NFRZ,NACT,NDEL,
     &                         IBSO,IREOST)
*----------------------------------------------------------------------*
*     store the 2-particle density on disk in symmetry-blocked form (but
*     redundant wrt. particle exchange), active orbitals only (all other
*     entries should be zero, as we deal with the irreducible 2-particle
*     density. the file structure is:
*
*     for each r,s of lambda(p,q;r,s) store
*
*     |(number of elements,zero flag)(elements of p,q-block)|
*
*----------------------------------------------------------------------*
      IMPLICIT NONE
      INCLUDE 'multd2h.inc'

      INTEGER, PARAMETER ::
     &     NTEST = 000

      INTEGER, INTENT(IN) ::
     &     IMOD, LUIR2DEN, IREW, NSMOB, NTOOB,
     &     NFRZ(NSMOB), NACT(NSMOB), NDEL(NSMOB),
     &     IBSO(*), IREOST(*)
      REAL(8), INTENT(IN) ::
     &     DEN(*)
      REAL(8), INTENT(INOUT) ::
     &     BUFF(*)

      INTEGER ::
     &     IOFF, IDX, JDX, NP, NQ, IZERO, NLEN,
     &     IPSYM, IQSYM, IRSYM, ISSYM, IRSSYM,
     &     IPORB, IQORB, IRORB, ISORB,
     &     IPTOT, IQTOT, IRTOT, ISTOT, IPQTOT, IRSTOT
      REAL(8) ::
     &     XNRM, FAC, XNRM_BLK, XNRM_TOT
      REAL(8), EXTERNAL ::
     &     INPROD


      IF (NTEST.GT.0) THEN
        WRITE(6,*) '====================='
        WRITE(6,*) 'entered WR_IR2DEN_BLK'
        WRITE(6,*) '====================='
        WRITE(6,*) ' IMOD,LUIR2DEN,IREW: ',IMOD,LUIR2DEN,IREW
        WRITE(6,*) ' NTOOB,NSMOB:       ',NTOOB,NSMOB
        WRITE(6,*) ' NFRZ:  ',NFRZ(1:NSMOB)
        WRITE(6,*) ' NACT:  ',NACT(1:NSMOB)
        WRITE(6,*) ' NDEL:  ',NDEL(1:NSMOB)
      END IF

      IF (IREW.EQ.1) THEN
        CALL REWINO(LUIR2DEN)
      END IF

      XNRM_TOT = 0D0
      DO ISSYM = 1, NSMOB
        DO ISORB = 1, NACT(ISSYM)
          ISTOT = IREOST(IBSO(ISSYM)+NFRZ(ISSYM)+ISORB-1)

          IF (NTEST.GE.1000) THEN
            WRITE(6,*) 'IBSO(ISSYM): ',IBSO(ISSYM)
            WRITE(6,*) 'NFRZ(ISSYM): ',NFRZ(ISSYM)
            WRITE(6,*) 'ISORB:       ',ISORB
            WRITE(6,*) 'INDEX (SO) : ',
     &           IBSO(ISSYM)+NFRZ(ISSYM)+ISORB-1
            WRITE(6,*) 'INDEX (TO) : ',ISTOT
          END IF

          DO IRSYM = 1, NSMOB
            DO IRORB = 1, NACT(IRSYM)
              XNRM_BLK = 0D0
              IRTOT = IREOST(IBSO(IRSYM)+NFRZ(IRSYM)+IRORB-1)
              IRSTOT = (ISTOT-1)*NTOOB + IRTOT
              
              IRSSYM = MULTD2H(ISSYM,IRSYM)

              IF (NTEST.GE.100) THEN
                WRITE(6,*) 'IRSYM,ISORB,ISSYM,IRORB: ',
     &                      IRSYM,ISORB,ISSYM,IRORB
                WRITE(6,*) 'IRTOT,ISTOT,IRSTOT:',IRTOT,ISTOT,IRSTOT
              END IF
              IF (NTEST.GE.1000) THEN
                WRITE(6,*) 'IBSO(IRSYM): ',IBSO(IRSYM)
                WRITE(6,*) 'NFRZ(IRSYM): ',NFRZ(IRSYM)
                WRITE(6,*) 'IRORB:       ',IRORB
                WRITE(6,*) 'INDEX (SO) : ',
     &               IBSO(IRSYM)+NFRZ(IRSYM)+IRORB-1
                WRITE(6,*) 'INDEX (TO) : ',IRTOT
              END IF

              IOFF = 0

              ! NOTE that loop is over IPSYM such that typically 
              ! NP >= NQ if we later restrict the loop
              DO IPSYM = 1, NSMOB
                IQSYM = MULTD2H(IRSSYM,IPSYM)
                
                NQ = NACT(IQSYM)
                NP = NACT(IPSYM)
                IF (IMOD.EQ.0) THEN
                  DO IQORB = 1, NQ
                    IQTOT = IREOST(IBSO(IQSYM)+NFRZ(IQSYM)+IQORB-1)
                    DO IPORB = 1, NP
                      IPTOT = IREOST(IBSO(IPSYM)+NFRZ(IPSYM)+IPORB-1)
                      IPQTOT = (IQTOT-1)*NTOOB + IPTOT
                      
c factors needed ?????
                      FAC =1D0
c 1:
c                      FAC = 0.25D0
c                      IF (IPORB.EQ.IQORB.OR.IRORB.EQ.ISORB.OR.
c     &                    IPORB.EQ.ISORB.OR.IQORB.EQ.IRORB) FAC = 0.5D0
c                      IF ((IPORB.EQ.IQORB.AND.IRORB.EQ.ISORB).OR.
c     &                    (IPORB.EQ.ISORB.AND.IQORB.EQ.IRORB)) FAC = 1D0
c 2:
c                      FAC = 0.5D0
c                      IF (IPORB.EQ.IQORB.OR.IRORB.EQ.ISORB.OR.
c     &                    IPORB.EQ.ISORB.OR.IQORB.EQ.IRORB) FAC = 0.5D0
c                      IF ((IPORB.EQ.IQORB.AND.IRORB.EQ.ISORB).OR.
c     &                    (IPORB.EQ.ISORB.AND.IQORB.EQ.IRORB)) FAC = 1D0
c 3:
c                      FAC = 0.5D0
c                      IF (IPORB.EQ.IQORB.OR.IRORB.EQ.ISORB.OR.
c     &                    IPORB.EQ.ISORB.OR.IQORB.EQ.IRORB.OR.
cc 3a:
c     &                    IPORB.EQ.IRORB.OR.IQORB.EQ.ISORB)  FAC = 1D0
c                      IF ((IPORB.EQ.IQORB.AND.IRORB.EQ.ISORB).OR.
c     &                    (IPORB.EQ.ISORB.AND.IQORB.EQ.IRORB)) FAC = 1D0
c 4:
c                      FAC = 0.5D0
c                      IF ((IPORB.EQ.IRORB.AND.IQORB.EQ.ISORB)) FAC = 1D0
cc 4a:
c                      IF ((IPORB.EQ.IQORB.AND.IRORB.EQ.ISORB).OR.
c     &                    (IPORB.EQ.ISORB.AND.IQORB.EQ.IRORB)) FAC = 1D0
c 5:
c                      FAC = 0.25D0
c                      IF (IPORB.EQ.IQORB.OR.IRORB.EQ.ISORB.OR.
c     &                    IPORB.EQ.ISORB.OR.IQORB.EQ.IRORB.OR.
c     &                    IPORB.EQ.IRORB.OR.IQORB.EQ.ISORB)  FAC = 0.5D0
c                      IF ((IPORB.EQ.IQORB.AND.IRORB.EQ.ISORB).OR.
c    &                    (IPORB.EQ.ISORB.AND.IQORB.EQ.IRORB)) FAC = 1D0
                      IDX = MAX(IPQTOT,IRSTOT)
                      JDX = MIN(IPQTOT,IRSTOT)
                      BUFF(IOFF+(IQORB-1)*NP+IPORB) =
     &                     FAC*DEN(IDX*(IDX-1)/2 + JDX)
                      
                    END DO
                  END DO
                ELSE
                  STOP 'do not enter this route ....'
                  DO IQORB = 1, NQ
                    IQTOT = IBSO(IQSYM)+NFRZ(IQSYM)+IQORB
                    DO IPORB = 1, NP
                      IPTOT = IBSO(IPSYM)+NFRZ(IPSYM)+IPORB
                      IPQTOT = (IQTOT-1)*NTOOB + IPTOT
                      
                      IDX = MAX(IPQTOT,IRSTOT)
                      JDX = MIN(IPQTOT,IRSTOT)
                      BUFF(IOFF+(IQORB-1)*NP+IPORB) =
     &                     DEN(IRSTOT*NTOOB*NTOOB + IPQTOT)
                      
                    END DO
                  END DO

                END IF

                IF (NTEST.GE.1000) THEN
                  WRITE(6,*) 'CURRENT BLOCK:',IRSTOT,IPSYM,IQSYM
                  CALL WRTMAT2(BUFF(IOFF+1),NACT(IPSYM),NACT(IQSYM),
     &                         NACT(IPSYM),NACT(IQSYM))
                END IF

                IOFF = IOFF + NQ*NP
              END DO ! IPQSYM

              NLEN = IOFF
              IF (NLEN.EQ.0) THEN
                IZERO = 1
              ELSE
                XNRM = SQRT(INPROD(BUFF,BUFF,NLEN))
                XNRM_BLK = XNRM_BLK + XNRM
                IZERO = 0
                IF (XNRM.LT.100d0*EPSILON(1d0)) IZERO = 1
              END IF
              WRITE(LUIR2DEN) NLEN,IZERO
              IF (NLEN.GT.0.AND.IZERO.EQ.0)
     &             WRITE(LUIR2DEN) BUFF(1:NLEN)

              IF (NTEST.GT.10) WRITE(6,*) 'Length and Norm of block ',
     &             IRSTOT,NLEN,XNRM_BLK
              XNRM_TOT = XNRM_TOT + XNRM_BLK

            END DO
          END DO 
        END DO
      END DO

      IF (NTEST.GT.10) WRITE(6,*) 'Norm of density: ', XNRM_TOT

      RETURN
      END
      SUBROUTINE ANA_CUMULANTS(T,ISPOBEX_TP,IOBEX_TP,
     &                 NSPOBEX_TP,NOBEX_TP,ISM,
     &                 THRES,MAXTRM,
     &                 IOCC_CA, IOCC_CB, IOCC_AA, IOCC_AB,
     &                 IUSLAB,
     &                 IOBLAB,NCPMT_SOX,WCPMT_SOX,
     &                        NCPMT_OX, WCPMT_OX,
     &                 NORB,IPRNCIV,ISOX_TO_OX,MSCOMB_CC)
*
* Analyze CUMUMANT matrices 
*
*      1) Print coefficients 
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given 
*         excitation type, spin orbital types and orbital types
*
* Jeppe Olsen , September 1999 
*               April 2005, in Bochum
*                                                  

*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX_TP)
      INTEGER IOBEX_TP(2*NGAS,NOBEX_TP)
      DIMENSION T(*)
      DIMENSION ISOX_TO_OX(NSPOBEX_TP)
*. Scratch
      INTEGER IOCC_CA(*),IOCC_CB(*),IOCC_AA(*),IOCC_AB(*)
*. Local scratch
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS) 
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
      CHARACTER*6 IOBLAB(*)
*. Output
      DIMENSION NCPMT_SOX(10,NSPOBEX_TP)                           
      DIMENSION WCPMT_SOX(10,NSPOBEX_TP)                          
      DIMENSION NCPMT_OX(10,NOBEX_TP)                           
      DIMENSION WCPMT_OX(10,NOBEX_TP)                          
*
      CALL ISETVC(NCPMT_SOX,0    ,10*NSPOBEX_TP)
      CALL SETVEC(WCPMT_SOX,0.0D0,10*NSPOBEX_TP)
      CALL ISETVC(NCPMT_OX,0    ,10*NOBEX_TP)
      CALL SETVEC(WCPMT_OX,0.0D0,10*NOBEX_TP)
*
* ===========================================================
*.1 : Printout of coefficients and largest occupation vectors       
*.2 : Group coefficients by type and size
* ===========================================================
*
      WRITE(6,*)
      WRITE(6,*) ' Operators are written as : '
      WRITE(6,*)
      WRITE(6,*)   ' Upper index alpha (UA)'
      WRITE(6,*)   ' Upper index beta  (UB)'
      WRITE(6,*)   ' Lower index alpha (LA)'
      WRITE(6,*)   ' Lower index beta  (LB)'
      WRITE(6,*)
      MINPRT = 0
      ITRM = 0
      ILOOP = 0
      NTVAR = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      TNORM = 0.0D0
      TTNORM = 0.0D0
2001  CONTINUE
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = 1.0D0
        XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNCIV.EQ.1) THEN
*. Print in one shot
       XMAX = 3006.1956
       XMIN =-3006.1956
      END IF
*
      WRITE(6,*)
      WRITE(6,'(A,E10.4,A,E10.4)')
     &'  Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(6,'(A)')
     &'  ========================================================='
      WRITE(6,*)
*
      IT = 0
      IIT = 0
      DO ITSS = 1, NSPOBEX_TP
C?     WRITE(6,*) ' NSPOBEX_TP = ', NSPOBEX_TP
*. Transform from occupations to groups
       CALL OCC_TO_GRP(ISPOBEX_TP(1+0*NGAS,ITSS),IGRP_CA,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+1*NGAS,ITSS),IGRP_CB,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+2*NGAS,ITSS),IGRP_AA,1      )
       CALL OCC_TO_GRP(ISPOBEX_TP(1+3*NGAS,ITSS),IGRP_AB,1      )
*
       NEL_CA = IELSUM(ISPOBEX_TP(1+0*NGAS,ITSS),NGAS)
       NEL_CB = IELSUM(ISPOBEX_TP(1+1*NGAS,ITSS),NGAS)
       NEL_AA = IELSUM(ISPOBEX_TP(1+2*NGAS,ITSS),NGAS)
       NEL_AB = IELSUM(ISPOBEX_TP(1+3*NGAS,ITSS),NGAS)
*. Diagonal restricted type of spinorbital excitation ?
       IF(MSCOMB_CC.EQ.1) THEN
         CALL DIAG_EXC_CC(
     &        ISPOBEX_TP(1+0*NGAS,ITSS),ISPOBEX_TP(1+1*NGAS,ITSS),
     &        ISPOBEX_TP(1+2*NGAS,ITSS),ISPOBEX_TP(1+3*NGAS,ITSS),
     &        NGAS,IDIAG)
       ELSE
         IDIAG = 0
       END IF
*
       IF(IPRNCIV.EQ.1) THEN
         WRITE(6,*) 
     &   ' Number of indeces in the various spaces for next type'
         WRITE(6,'(A,5I5)') 
     &   ' Upper alpha :', (ISPOBEX_TP(I+0*NGAS,ITSS),I=1,NGAS)
         WRITE(6,'(A,5I5)') 
     &   ' Upper beta : ', (ISPOBEX_TP(I+1*NGAS,ITSS),I=1,NGAS)
         WRITE(6,'(A,5I5)') 
     &   ' Lower alpha :', (ISPOBEX_TP(I+2*NGAS,ITSS),I=1,NGAS)
         WRITE(6,'(A,5I5)') 
     &   ' Lower beta  :', (ISPOBEX_TP(I+3*NGAS,ITSS),I=1,NGAS)
       END IF
*
       DO ISM_C = 1, NSMST
        ISM_A = MULTD2H(ISM,ISM_C) 
        DO ISM_CA = 1, NSMST
         ISM_CB = MULTD2H(ISM_C,ISM_CA)
         DO ISM_AA = 1, NSMST
          ISM_AB =  MULTD2H(ISM_A,ISM_AA)
          ISM_ALPHA = (ISM_AA-1)*NSMST + ISM_CA
          ISM_BETA  = (ISM_AB-1)*NSMST + ISM_CB
*
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.LT.ISM_BETA) GOTO 777
          IF(IDIAG.EQ.1.AND.ISM_ALPHA.EQ.ISM_BETA) THEN 
           IRESTRICT_LOOP = 1
          ELSE
           IRESTRICT_LOOP = 0
          END IF
*. obtain strings
          CALL GETSTR2_TOTSM_SPGP(IGRP_CA,NGAS,ISM_CA,NEL_CA,NSTR_CA,
     &         IOCC_CA, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_CB,NGAS,ISM_CB,NEL_CB,NSTR_CB,
     &         IOCC_CB, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_AA,NGAS,ISM_AA,NEL_AA,NSTR_AA,
     &         IOCC_AA, NORB,0,IDUM,IDUM)
          CALL GETSTR2_TOTSM_SPGP(IGRP_AB,NGAS,ISM_AB,NEL_AB,NSTR_AB,
     &         IOCC_AB, NORB,0,IDUM,IDUM)
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
          DO I_AB = 1, NSTR_AB
           IF(IRESTRICT_LOOP.EQ.1) THEN
             I_AA_MIN = I_AB
           ELSE
             I_AA_MIN = 1
           END IF
           DO I_AA = I_AA_MIN, NSTR_AA
            DO I_CB = 1, NSTR_CB
             IF(IRESTRICT_LOOP.EQ.1.AND.I_AA.EQ.I_AB) THEN
               I_CA_MIN = I_CB
             ELSE
               I_CA_MIN = 1
             END IF
             DO I_CA = I_CA_MIN, NSTR_CA
              IT = IT + 1
C?            WRITE(6,*) ' IT, T(IT) = ', IT,T(IT)
*
              IF(ILOOP .EQ. 1 ) THEN 
                NTVAR = NTVAR + 1
                TNORM = TNORM + T(IT)**2
*. Classify element according to size 
                DO IPOT = 1, 10
                  IF(10.0D0 ** (-IPOT+1).GE.ABS(T(IT)).AND.
     &            ABS(T(IT)).GT. 10.0D0 ** ( - IPOT )      )THEN
                    IOEX = ISOX_TO_OX(ITSS)
                    NCPMT_OX(IPOT,IOEX)= NCPMT_OX(IPOT,IOEX)+1  
                    WCPMT_OX(IPOT,IOEX)= WCPMT_OX(IPOT,IOEX)+T(IT)**2
                  END IF
                END DO
*             ^ End of loop over powers of ten
              END IF
*             ^ end if we are in loop 1
              IF( XMAX .GE. ABS(T(IT)) .AND.
     &        ABS(T(IT)).GT. XMIN .AND. ITRM.LE.MAXTRM) THEN
                ITRM = ITRM + 1
                IIT = IIT + 1
                TTNORM = TTNORM + T(IT) ** 2
                WRITE(6,'(A)')
                WRITE(6,'(A)')
     &          '                 =================== '
                WRITE(6,*)

                WRITE(6,'(A,2I8,2X,E14.8)')
     &          '  size of amplitude : ',T(IT)
                IF(IUSLAB.EQ.0) THEN
                  WRITE(6,'(A,4X,10I4)')
     &            'UA', (IOCC_CA(IEL+(I_CA-1)*NEL_CA),IEL = 1, NEL_CA)
                  WRITE(6,'(A,4X,10I4)')
     &            'UB', (IOCC_CB(IEL+(I_CB-1)*NEL_CB),IEL = 1, NEL_CB)
                  WRITE(6,'(A,4X,10I4)')
     &            'LA', (IOCC_AA(IEL+(I_AA-1)*NEL_AA),IEL = 1, NEL_AA)
                  WRITE(6,'(A,4X,10I4)')
     &            'LB', (IOCC_AB(IEL+(I_AB-1)*NEL_AB),IEL = 1, NEL_AB)
                END IF
              END IF
*             ^ End if could and should be printed
             END DO
*            ^ End of loop over alpha creation strings
            END DO
*           ^ End of loop over beta creation strings
           END DO
*          ^ End of loop over alpha annihilation 
          END DO 
*         ^ End of loop over beta annihilation 
  777    CONTINUE
         END DO
        END DO
       END DO
*      ^ End of loop over symmetry blocks
      END DO
*     ^ End of loop over over types of excitations
       IF(IIT .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
       IF( XMIN .GT. THRES .AND. ILOOP .LE. 20 ) GOTO 2001
*
       WRITE(6,'(A,E15.8)')
     & '  Norm of printed coefficients .. ', TTNORM
       WRITE(6,'(A,E15.8)')
     & '  Norm of all     coefficients .. ',  TNORM
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') '   Magnitude of T coefficients '
      WRITE(6,'(A)') '  =============================='
      WRITE(6,'(A)')
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, 10
        W = 0.0D0
        N = 0
        DO 290 ITSS = 1, NOBEX_TP
            N = N + NCPMT_OX(IPOT,ITSS)                    
            W = W + WCPMT_OX(IPOT,ITSS)                    
  290   CONTINUE
        WACC = WACC + W
        NACC = NACC + N
        WRITE(6,'(A,I2,A,I2,3X,I9,X,E15.8,3X,E15.8)')
     &  '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
  300 CONTINUE
*
      WRITE(6,*) ' Number of coefficients less than  10-11',
     &           ' IS  ',NTVAR - NACC
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Magnitude of CI coefficients for each type of excitation '
      WRITE(6,'(A)') 
     & '  ========================================================='
      WRITE(6,'(A)')
      DO 400 ITSS   = 1, NOBEX_TP 
          N = 0
          DO 380 IPOT = 1, 10
            N = N + NCPMT_OX(IPOT,ITSS)                     
  380     CONTINUE
          IF(N .NE. 0 ) THEN
            WRITE(6,*) ' Orbital excitation type = ', ITSS
            WRITE(6,'(A,I9)')  
     &      '         Number of coefficients larger than 10-11 ', N
            WRITE(6,*)
            WACC = 0.0D0
            DO IPOT = 1, 10
              N =  NCPMT_OX(IPOT,ITSS)                    
              W =  WCPMT_OX(IPOT,ITSS)                    
              WACC = WACC + W
              WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC           
            END DO
          END IF 
  400 CONTINUE
*
*. Total weight and number of dets per excitation level
*
      WRITE(6,'(A)')
      WRITE(6,'(A)') 
     & '   Total weight and number of SD''s (> 10 ** -11 )  : '          
      WRITE(6,'(A)') 
     & '  ================================================='
      WRITE(6,'(A)')
      WRITE(6,*) '        N      Weight      Acc. Weight   Exc. type '
      WRITE(6,*) ' ==================================================='
      WACC = 0.0D0
      DO ITSS = 1, NOBEX_TP  
          N = 0
          W = 0.0D0
          DO IPOT = 1, 10
            N = N + NCPMT_OX(IPOT,ITSS)                   
            W = W + WCPMT_OX(IPOT,ITSS)                   
          END DO
          WACC = WACC + W
          IF(N .NE. 0 ) THEN
            WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,I3)') 
     &      N,W,WACC,ITSS                                   
          END IF
      END DO
*
      RETURN
      END
      SUBROUTINE TP_OBEX3(NOCCLS,NEL,NGASX,IOBEX_TP,LCOBEX_TP,LAOBEX_TP,
     &                   IOCCLS,IOCCLS_REF,MX_NCREA,MX_NANNI,
     &                   MX_EXC_LEVEL,IEXTP_TO_OCCLS,MX_AAEXC,IFLAG,
     &                   I_OOCC,NOBEX_TP,NOAAEX,MX_EXT_IND,IPRCC)
*
* Slightly modified version of TP_OBEX2 : Different way of counting 
* allowed number active indeces : MX_AAEXC is now the allowed number 
*                                 of active indeces in the operator
*
* Obtain the orbital excitation types needed to generate occupation classes 
* in IOCCLS from IOCCLS_REF
* MX_EXT_IND( Largest number of inactive or active indeces added March
* 2010
* 
* Jeppe Olsen, October 2006
*              March 2010: MX_EXT_IND added
*
* If IFLAG = 1, only the number of orbital excitation types is generated
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IOCCLS(NGAS,*), IOCCLS_REF(NGAS)
*. Output
      INTEGER IOBEX_TP(2*NGAS,*),  LCOBEX_TP(*), LAOBEX_TP(*)
      INTEGER IEXTP_TO_OCCLS(*)
*. Local scratch
      DIMENSION ICREA(MXPNGAS),IANNI(MXPNGAS)
*. Number of active orbital spaces
      NACT_SPC = 0
      IACT_SPC = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3) THEN
          NACT_SPC = NACT_SPC + 1
          IACT_SPC = IGAS
        END IF
      END DO
*
      IF(NACT_SPC.GT.1) THEN
        WRITE(6,*) ' TP_OBEX3 in problems '
        WRITE(6,*) ' More than one active orbital spaces '
        WRITE(6,*) ' NACT_SPC = ',  NACT_SPC 
        STOP ' TP_OBEX3 :  More than one active orbital spaces '
      END IF
C?    WRITE(6,*) ' TP_OBEX3, MX_NCREA, MX_NANNI = ', 
C?   &                       MX_NCREA, MX_NANNI
*
* The orbital excitation operator IEXTP is  organized as 
*
* LCOBEX(IEXTP) : Number of creation operators    
* LAOBEX(IEXTP) : Number of annihilation operators
* IOBEX_TP(1 - NGAS, IEXTP) : Number of creation operators per gassspace
* IOBEX(NGAS+1  -  2*NGAS, , IEXTP) : Number of annihilation operators 
*
* IEXTP_TO_OCCLS is map from orbital excitation type to occupation 
* class for CI coefficients
      NTEST = 03
      NTEST = MAX(NTEST,IPRCC)
      NOBEX_TP = 0
      MX_EXC_LEVEL = 0
      JREFCLS = 0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'TP_OBEX speaking '
        WRITE(6,*) ' Reference occupation class ' 
        CALL IWRTMA(IOCCLS_REF,1,NGAS,1,NGAS)
        WRITE(6,*) '  IPRCC = ', IPRCC
        WRITE(6,*) ' MX_AAEXC = ', MX_AAEXC
      END IF
*
      IZERO = 0
      DO JOCCLS = 1, NOCCLS 
        NANNI = 0
        NCREA = 0
        CALL ISETVC(ICREA,IZERO,NGAS)
        CALL ISETVC(IANNI,IZERO,NGAS)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Excited occupation class '
          CALL IWRTMA(IOCCLS(1,JOCCLS),1,NGAS,1,NGAS)
        END IF
        DO IGAS = 1, NGAS
          IF(IOCCLS(IGAS,JOCCLS).GT.IOCCLS_REF(IGAS)) THEN
            ICREA(IGAS) = IOCCLS(IGAS,JOCCLS) - IOCCLS_REF(IGAS)
            NCREA = NCREA + ICREA(IGAS)
          ELSE IF (IOCCLS(IGAS,JOCCLS).LT.IOCCLS_REF(IGAS)) THEN
            IANNI(IGAS) = -( IOCCLS(IGAS,JOCCLS) - IOCCLS_REF(IGAS))
            NANNI = NANNI + IANNI(IGAS)
          END IF
        END DO
        IF(NANNI.EQ.0.AND.NCREA.EQ.0) JREFCLS = JOCCLS
*. Add active -active excitation 
        DO IAA_EXC = 0, MX_AAEXC
         ITSOKAY = 1
         IF(IAA_EXC.GT.0) THEN
C?        WRITE(6,*) ' IANNI(IACT_SPC), ICREA(IACT_SPC) ',
C?   &                 IANNI(IACT_SPC), ICREA(IACT_SPC) 
C?        WRITE(6,*) '  NANNI, NCREA,  MX_NANNI, MX_NCREA ',
C?   &                  NANNI, NCREA,  MX_NANNI, MX_NCREA
*. Can another active-active excitation be added ?
          IF(IANNI(IACT_SPC)+ICREA(IACT_SPC)+2.LE.MX_AAEXC.AND.
     &       (NANNI .LT. MX_NANNI. AND. NCREA .LT. MX_NCREA)) THEN
              ITSOKAY = 1
          ELSE 
              ITSOKAY = 0
          END IF
         END IF
         IF(ITSOKAY.EQ.1) THEN
*
          IF(IAA_EXC.GT.0) THEN
            ICREA(IACT_SPC)  = ICREA(IACT_SPC) + 1
            IANNI(IACT_SPC)  = IANNI(IACT_SPC) + 1
            NCREA = NCREA + 1
            NANNI = NANNI + 1
          END IF
*     If we do not include pure active-active rotations:
          IF (NOAAEX.EQ.1) THEN
*     Test whether this is one:
            IPAA = 1
            DO JGAS = 1, NGAS
              IF (JGAS.NE.IACT_SPC.AND.ICREA(JGAS).NE.0) IPAA=0
              IF (JGAS.NE.IACT_SPC.AND.IANNI(JGAS).NE.0) IPAA=0
            END DO
            IF (IPAA.EQ.1) CYCLE
          END IF
*- Number of external indeces
          N_EXT_IND = 0
          DO JGAS = 1, NGAS
            IF (JGAS.NE.IACT_SPC) 
     &      N_EXT_IND = N_EXT_IND +ICREA(JGAS) + IANNI(JGAS)
          END DO
C?        WRITE(6,*) ' N_EXT_IND, MX_EXT_IND = ', 
C?   &                 N_EXT_IND, MX_EXT_IND
*
          MX_EXC_LEVEL = MAX(MX_EXC_LEVEL,NCREA)
          IF(NCREA+NANNI.NE.0.AND.N_EXT_IND.LE.MX_EXT_IND.AND.
     &         (I_OOCC.EQ.0.OR.(NCREA.NE.1.AND.NANNI.NE.1))) THEN 
            NOBEX_TP = NOBEX_TP + 1
            IF(IFLAG.NE.1) THEN
              LCOBEX_TP(NOBEX_TP) = NCREA
              LAOBEX_TP(NOBEX_TP) = NANNI
              IEXTP_TO_OCCLS(NOBEX_TP) = JOCCLS
              CALL ICOPVE(ICREA,IOBEX_TP(1,NOBEX_TP),NGAS )
              CALL ICOPVE(IANNI,IOBEX_TP(NGAS+1,NOBEX_TP),NGAS )
            END IF
          END IF
*
         END IF
        END DO
*       ^ End of loop over active-active excitations
      END DO
*.    ^ End of loop over  occupation classes
*
      IF(IFLAG.NE.1) THEN
*. Add unit operator as excition NOBEX_TP + 1
        LCOBEX_TP(NOBEX_TP + 1) = 0
        LAOBEX_TP(NOBEX_TP + 1) = 0
        IEXTP_TO_OCCLS(NOBEX_TP+1) = JREFCLS
        IZERO = 0
        CALL ISETVC(IOBEX_TP(1,NOBEX_TP+1),IZERO,NGAS)
        CALL ISETVC(IOBEX_TP(NGAS+1,NOBEX_TP+1),IZERO,NGAS)
      END IF
*
      IF(NTEST.GE.3) THEN
        WRITE(6,*) ' Largest excitation level : ', MX_EXC_LEVEL
        WRITE(6,*)
        WRITE(6,*) ' Number of types of orbital excitations (w. unit) ',
     &  NOBEX_TP+1
        WRITE(6,*)
*
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' Creation part,  Annihilation  part '
          WRITE(6,*) ' ==================================='
          DO IOBEX = 1, NOBEX_TP+1
            WRITE(6,'(16I4,16I4)')
     &      (IOBEX_TP(I,IOBEX),I=1, NGAS),
     &      (IOBEX_TP(NGAS+I,IOBEX),I=1, NGAS) 
          END DO
*
          WRITE(6,*) ' Orbital excitation type to occupation class '
          CALL IWRTMA(IEXTP_TO_OCCLS,1,NOBEX_TP+1,1,NOBEX_TP+1)
        END IF
*
      END IF
*
      RETURN
      END
c $Id$
