      SUBROUTINE WR_CPKJJB_MS
     &           (CPKJJB_MS,NP,NK,ICSM,KSM,NOB,NJB_AC,NSM,IB)
*
* Print multiple symmetry matrix C(P,Ka,j,Jb_ac)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*
      INTEGER NP(NSM), NOB(NSM),NJB_AC(NSM) 
      INTEGER IB(20,20)
      DIMENSION CPKJJB_MS(*)
*
      DO JSM = 1, NSM
       DO JB_AC_SM = 1, NSM
         J_JB_SM = MULTD2H(JSM,JB_AC_SM)
         K_J_JB_SM = MULTD2H(KSM,J_JB_SM)
         IPSM = MULTD2H(ICSM, K_J_JB_SM) 
         NPK = NP(IPSM)*NK
         IF(NPK*NOB(JSM)*NJB_AC(JB_AC_SM).GT.0) THEN
*
         WRITE(6,*) 
     &   ' ================================================'
         WRITE(6,*) 
     &   ' C(PKa,j,Jb_ac) with sym of j, Jb_ac =',JSM,Jb_ac_sm 
         WRITE(6,*) 
     &   ' ================================================'
         WRITE(6,*)
         IOFF = IB(JSM,JB_AC_SM) 
C?       WRITE(6,*) 
C?   &   ' IOFF NPK NJ, NJB', IOFF,NPK,NOB(JSM),NJB_AC(JB_AC_SM)
         CALL WR_CKJJB(CPKJJB_MS(IOFF),NPK,NOB(JSM),NJB_AC(JB_AC_SM))
*
         END IF
       END DO
      END DO
*
      RETURN
      END
      SUBROUTINE WR_CKJJB(CKJJB,NK,NJ,NJB)
*
* Print three-dimensional matrix C(Ka,j,Jb)
*
      INCLUDE 'implicit.inc'
      DIMENSION CKJJB(NK,NJ,NJB)
*
C     WRITE(6,*) ' Matrix C(Ka,J,Jb) ' 
C     WRITE(6,*) ' =================='
      WRITE(6,*)
      DO JB = 1, NJB
        WRITE(6,*) ' JB = ', JB
        WRITE(6,*)
        CALL WRTMAT(CKJJB(1,1,JB),NK,NJ,NK,NJ)
      END DO
*
      RETURN
      END
      FUNCTION LEN_TT_BLOCK(ISM,NIA,NIB,NSMST)
*
* Length of TT block with total sym ISM
*
* Jeppe Olsen, May 99 in Aarhus
*
      INCLUDE 'implicit.inc'
*. Specific input
      INTEGER NIA(*), NIB(*)
*. General input
      INCLUDE 'multd2h.inc'
*
      LEN = 0
      DO IASM = 1, NSMST
        IBSM = MULTD2H(IASM,ISM)
        LEN = LEN + NIA(IASM)*NIB(IBSM)
      END DO
*
      LEN_TT_BLOCK = LEN
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' LEN_TT_BLOCK = ', LEN_TT_BLOCK 
      END IF
*
      RETURN
      END
      SUBROUTINE IOFF_C_P_KA_J_JB(NP,JTYP,NOBPTS,MXPNGAS,LKA,KASM,
     &           ICSM,NSMST,NJB,IOFFA,LEN)
*
* Offset array for C(P,Ka,J,JB)
*
* Jeppe Olsen, Last night of May 1999 , at dept. of chem. Aarhus 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      DIMENSION NP(*), NOBPTS(MXPNGAS,*),NJB(*) 
*. Output
      DIMENSION IOFFA(20,20)
*
      IOFF = 1
      DO JSM = 1, NSMST 
        NJ = NOBPTS(JTYP,JSM)
        DO JB_AC_SM = 1, NSMST
          JA_AC_SM  = MULTD2H(KASM,JSM)
          JAB_AC_SM = MULTD2H(JA_AC_SM,JB_AC_SM)
          JAB_PA_SM = MULTD2H(ICSM,JAB_AC_SM)
          NNP =    NP(JAB_PA_SM)
          NNJB = NJB(JB_AC_SM)
          IOFFA(JSM,JB_AC_SM)= IOFF
* C(P,Jaa,Jba) => C(P,Ka,j,Jba)
          IOFF = IOFF + NNP*NJ*LKA*NNJB
C?        WRITE(6,*) ' JSM, JB_AC_SM, JAB_PA_SM', 
C?   &                 JSM, JB_AC_SM, JAB_PA_SM
C?        WRITE(6,*)  ' IOFF, NNP,NJ,LKA,NNJB',IOFF,  NNP,NJ,LKA,NNJB
        END DO
      END DO
      LEN = IOFF - 1
C     WRITE(6,*) ' LEN = ', LEN
      NTEST = 0
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' Offset array for C(P,Ka,j,Jba), Kasm= ', KASM 
      CALL IWRTMA(IOFFA,NSMST,NSMST,20,20)
      END IF
*
      RETURN
      END
      SUBROUTINE IOFF_CJPA_MS(NSMST,NAA,NAB,NPA,NPB,NPAB,NAAB,IBCPA,
     &                       IBCPA2,ICSM)
*
* Obtain offsets for Multiple symmetry C(Jpab,Jaa,Jab) = C(Jpa,Jpb,Jaa,Jab) 
* blocked matrix
*
* Order of matrix is 
*
* Loop over symmetry of Jaa*Jab => Symmetry of Jpa*Jpb
*  Loop over Symmetry of Jaa => Symmetry of Jab
*    Loop over all Jab            
*     Loop over all Jaa
*       Loop over Symmetry of Jpa => Symmetry of Jpb
*        Loop over Jpb
*         Loop over Jpa
*         End of loop over Jpa
*        End of loop over Jpb
*       End of loop over symmetry of Jpa
*     End of Loop over Jaa
*    End of Loop over Jab
*  End of loop over symmetries of Jaa
* End of Loop over symmetry of Jaa*Jab
*
* All products of passive strings with given total symmetry are
* therefore grouped together, maximizing the number of rows,
* for a given column
*
* Symmetry restricted to D2H at the moment
*
* Jeppe Olsen, Late April 99 In Aarhus
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* =======
*  Input
* =======
*. Number of strings per symmetry
C     INTEGER NA(NSMST),NB(NSMST)
      INTEGER NPA(NSMST),NPB(NSMST),NAA(NSMST),NAB(NSMST)
      INCLUDE 'multd2h.inc'
      
*
* ========
*. Output 
* ========
* 
*. Offset in CPA of dets with given symmetry of active product
      INTEGER IBCPA(NSMST)
*. Offset in CPA of dets with given sym af active alpha and active beta
      INTEGER IBCPA2(20,20)          
*. Number of product strings with given symmetry
      INTEGER NPAB(NSMST), NAAB(NSMST)
*. Number of active strings JAA*JAB with given sym
      DO JA_SM = 1, NSMST
        LA = 0
        DO JAA_SM = 1, NSMST
          JAB_SM = MULTD2H(JA_SM,JAA_SM)
          LA = LA + NAA(JAA_SM)*NAB(JAB_SM)
        END DO
        NAAB(JA_SM) = LA
      END DO
*. Number of passive strings JPA*JPB with given sym
      DO JP_SM = 1, NSMST
        LP = 0
        DO JPA_SM = 1, NSMST
          JPB_SM = MULTD2H(JP_SM,JPA_SM)
          LP = LP + NPA(JPA_SM)*NPB(JPB_SM)
        END DO
        NPAB(JP_SM) = LP
      END DO
*. Offset to block with given sym of active blocks in PA matrix
      IOFF = 1
      DO IASM = 1, NSMST
        IBCPA(IASM) = IOFF
        IPSM = MULTD2H(IASM,ICSM)
        LBLOCK = NAAB(IASM)*NPAB(IPSM)
        IOFF = IOFF + LBLOCK
      END DO
*. Offset to block with given sym of Jaa and Jba
      IOFF = 1
      DO JAB_AC_SM = 1, NSMST
      DO JA_AC_SM = 1, NSMST
       JB_AC_SM = MULTD2H(JAB_AC_SM,JA_AC_SM)
       IBCPA2(JA_AC_SM,JB_AC_SM) = IOFF
*
       JAB_PA_SM = MULTD2H(ICSM,JAB_AC_SM)
       NAB_PA =    NPAB(JAB_PA_SM)
       IOFF = IOFF + NAB_PA*NAA(JA_AC_SM)*NAB(JB_AC_SM)
      END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Offsets for C(Pab,Jaa,Jab) array '
        WRITE(6,*) 
        WRITE(6,*) ' Offsets for active strings with given sym'
        CALL IWRTMA(IBCPA,1,NSMST,1,NSMST)
        WRITE(6,*)
        WRITE(6,*) ' Offsets for Jaa, Jab with given sym '
        CALL IWRTMA(IBCPA2,NSMST,NSMST,20,20)       
      END IF
*
      RETURN
      END
      SUBROUTINE CJPA_MS(C,CPA,IWAY,NSMST,NA,NB,NAA,NAB,NPA,NPB,
     &                 IAREO,IBREO,IBAREO,IBBREO,ICSM,SIGN,
     &                 NPAB,NAAB,IBCPA,IBCPA2,LEN,I12ORD,NOREO)
*
* Standard <=> Passive/active division of Blocks of coefficients.
*
* IWAY = 1 : Standard => Active/Passive form
* IWAY = 2 : Active/Passive => Standard form
*
* IF NOREO = 1, no reordering is done, the addressing arrays are
*               however set up
*
* Version were all symmetryblocks belonging to a given occupation
* are included.
*
* PA form is C(Jpa,Jpb,Jaa,Jab) with  ordering blocks as
*
* Loop over symmetry of Jaa*Jab => Symmetry of Jpa*Jpb
*  Loop over Symmetry of Jaa => Symmetry of Jab
*    Loop over all Jab            
*     Loop over all Jaa
*       Loop over Symmetry of Jpa => Symmetry of Jpb
*        Loop over Jpb
*         Loop over Jpa
*         End of loop over Jpa
*        End of loop over Jpb
*       End of loop over symmetry of Jpa
*     End of Loop over Jaa
*    End of Loop over Jab
*  End of loop over symmetries of Jaa
* End of Loop over symmetry of Jaa*Jab
*
* All products of passive strings with given total symmetry are
* therefore grouped together, maximizing the number of rows,
* for a given column
*
* Symmetry restricted to D2H at the moment
*
* Jeppe Olsen, Late April 98 In Aarhus
*              Finished at Hampton Inn, Atlanta, Georgia, March 99
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* =======
*  Input
* =======
*. Symmetry of CBLOCK
      INTEGER ICSM
      
*. Number of strings per symmetry
      INTEGER NA(NSMST),NB(NSMST)
      INTEGER NPA(NSMST),NPB(NSMST),NAA(NSMST),NAB(NSMST)
*. Reorder arrays, Passive/Active => Standard form
      INTEGER IAREO(*),IBREO(*)
*. Offset for given symmetry block in standard format
      INTEGER IBAREO(8,8)        ,IBBREO(8,8)        
*. Note : The reorder arrays contain info for all symmetries.
      INCLUDE 'multd2h.inc'
      
*
* ========
*. Output 
* ========
* 
*. Offset in CPA of dets with given symmetry of active product
      INTEGER IBCPA(NSMST)
*. Offset in CPA of dets with given sym af active alpha and active beta
*
      INTEGER IBCPA2(20,20)        
*. Number of product strings with given symmetry
      INTEGER NPAB(NSMST), NAAB(NSMST)
*  
* ==================
*. Input and Output
* ==================
      DIMENSION C(*),CPA(*)
*
* =========
*. Scratch 
* =========
*
*. offset to strings of given sym
      INTEGER IAOFF(20),IBOFF(20)
*. Offset to block in standard order
      INTEGER IBCST(20)
*. Offset to passive(active) ab strings with given totsym and given  
*  sym of alpha. Offset is defined with respect to start of 
*  given totsym
      INTEGER IBPA2(20,20)
      INTEGER IBAC2(20,20)

*
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CJPA : Initial C matrix '
        CALL WRTVH1(C,ICSM,NA,NB,NSMST,0)
        WRITE(6,*) ' Sign = ', SIGN
      END IF
*
*. Arrays giving symmetry info on PA ordered strings and matrices
*. Number of active strings JAA*JAB with given sym
      DO JA_SM = 1, NSMST
        LA = 0
        DO JAA_SM = 1, NSMST
          JAB_SM = MULTD2H(JA_SM,JAA_SM)
          LA = LA + NAA(JAA_SM)*NAB(JAB_SM)
        END DO
        NAAB(JA_SM) = LA
      END DO
*. Number and offset of passive strings JPA*JPB with given sym
      DO JP_SM = 1, NSMST
        LP = 0
        IBPA2(1,JP_SM) = 1
        DO JPA_SM = 1, NSMST
          JPB_SM = MULTD2H(JP_SM,JPA_SM)
          IF(JPA_SM.NE.1) 
     &    IBPA2(JPA_SM,JP_SM) = IBPA2(1,JP_SM)+LP
          LP = LP + NPA(JPA_SM)*NPB(JPB_SM)
        END DO
        NPAB(JP_SM) = LP
      END DO
*. Offset for active ab strings with given a and totsym
      DO J_SM = 1, NSMST
        L = 0
        IBAC2(1,J_SM) = 1
        DO JA_SM = 1, NSMST
          JB_SM = MULTD2H(J_SM,JA_SM)
          IF(JA_SM.NE.1) 
     &    IBAC2(JA_SM,J_SM) = IBAC2(1,J_SM)+L
          L = L + NAA(JA_SM)*NPB(JB_SM)
        END DO
      END DO
*. Offset to block with given sym of active blocks in PA matrix
      IOFF = 1
      DO IASM = 1, NSMST
        IBCPA(IASM) = IOFF
        IPSM = MULTD2H(IASM,ICSM)
        LBLOCK = NAAB(IASM)*NPAB(IPSM)
        IOFF = IOFF + LBLOCK
      END DO
*. Offsets to strings of given sym
      IBA = 1
      IBB = 1
      DO ISM = 1, NSMST
       IAOFF(ISM) = IBA
       IBOFF(ISM) = IBB
       IBA = IBA + NA(ISM)
       IBB = IBB + NB(ISM)
      END DO
*. Offset to sym block, standard order
      IBASE = 1
C     I12ORD = 1
      DO IASM = 1, NSMST
        IBSM = MULTD2H(IASM,ICSM)
        IF(I12ORD.EQ.1) THEN
          IBCST(IBSM) = IBASE
        ELSE
          IBCST(IASM) = IBASE
        END IF
C       IBCST(IASM) = IBASE
        IF(I12ORD.EQ.1) THEN
          IBASE = IBASE + NB(IASM)*NA(IBSM)
        ELSE
          IBASE = IBASE + NA(IASM)*NB(IBSM)
        END IF
C       IBASE = IBASE + NA(IASM)*NB(IBSM)
      END DO
*. Offset to ci block with given sym of active alpha and active beta
      IBASE = 1
      DO J_AC_AB_SM = 1, NSMST
       DO J_AC_A_SM = 1, NSMST
         J_AC_B_SM = MULTD2H(J_AC_A_SM,J_AC_AB_SM)
         IBCPA2(J_AC_A_SM,J_AC_B_SM)=IBASE
         J_PA_AB_SM = MULTD2H(ICSM,J_AC_AB_SM)
         IBASE = 
     &   IBASE + NPAB(J_PA_AB_SM)*NAA(J_AC_A_SM)*NAB(J_AC_B_SM)
       END DO
      END DO
      IF(NOREO.EQ.0) THEN
*
*. Loop over blocks in passive/active order
C     J_AB_PA = 0
      DO J_AC_AB_SM = 1, NSMST
       J_PA_AB_SM = MULTD2H(ICSM,J_AC_AB_SM)
       DO J_AC_A_SM = 1, NSMST
         J_AC_B_SM = MULTD2H(J_AC_A_SM,J_AC_AB_SM)
         IF(NTEST.GE.1000) 
     &   WRITE(6,*) ' J_AC_A_SM, J_AC_B_SM', J_AC_A_SM, J_AC_B_SM
C        IBCPA2(J_AC_A_SM,J_AC_B_SM)=J_AB_PA+1
         J_AB_PA_OF = IBCPA2(J_AC_A_SM,J_AC_B_SM)
*. And loop over the actual active strings
         L_AC_A_SM = NAA(J_AC_A_SM)
         L_AC_B_SM = NAB(J_AC_B_SM)
*
           DO J_PA_A_SM = 1, NSMST
         J_ACAB_PA = 0
C        DO J_AC_A = 1, L_AC_A_SM
C        DO J_AC_B = 1, L_AC_B_SM
*. And a new Jaa Jpb String has been born 
            J_PAAB_PA_OF = IBPA2(J_PA_A_SM,J_PA_AB_SM)
            J_PA_B_SM = MULTD2H(J_PA_A_SM,J_PA_AB_SM)
            L_PA_A_SM = NPA(J_PA_A_SM)
            L_PA_B_SM = NPB(J_PA_B_SM)
            IF(L_PA_A_SM *  L_PA_B_SM .NE. 0 ) THEN
*. Symmetry of alpha and betastrings
            J_A_SM = MULTD2H(J_AC_A_SM,J_PA_A_SM)
            J_B_SM = MULTD2H(J_AC_B_SM,J_PA_B_SM)
*. Length of strings
            L_A_SM = NA(J_A_SM)
            L_B_SM = NB(J_B_SM)
*. Offset of these strings from start of stringtype, reoordered
            JA_OF = IBAREO(J_AC_A_SM,J_A_SM)+IAOFF(J_A_SM) - 1
            JB_OF = IBBREO(J_AC_B_SM,J_B_SM)+ IBOFF(J_B_SM) - 1
*. Start of this block of CI coef, standard order
            J_AB_ST_OF = IBCST(J_A_SM)
*
         DO J_AC_A = 1, L_AC_A_SM
         DO J_AC_B = 1, L_AC_B_SM
           J_ACAB_PA = J_ACAB_PA + 1
            J_A_ST_PA_0   = JA_OF + (J_AC_A-1)*L_PA_A_SM - 1 
            J_B_ST_PA     = JB_OF + (J_AC_B-1)*L_PA_B_SM - 1
            J_AB_PA2_0 = IBCPA2(J_AC_A_SM,J_AC_B_SM) - 1 + 
     &                (J_ACAB_PA-1)*NPAB(J_PA_AB_SM)
            DO J_PA_B = 1, L_PA_B_SM
              J_B_ST_PA  = J_B_ST_PA + 1
              J_B_ST = IBREO(J_B_ST_PA  )
              J_A_ST_PA = J_A_ST_PA_0
              J_AB_ST_0 =  J_AB_ST_OF - 1 -L_B_SM + J_B_ST
              J_PAAB_PA_OF2 = J_PAAB_PA_OF -1 + (J_PA_B-1)*L_PA_A_SM
*
            IF(IWAY.EQ.1) THEN
             J_AB_PA2 =  J_AB_PA2_0 +  J_PAAB_PA_OF2 
C            DO J_PA_A = 1, L_PA_A_SM
             DO JJ_A_ST_PA = J_A_ST_PA+1,J_A_ST_PA+ L_PA_A_SM
*. Strings in standard format
C             J_A_ST_PA =  J_A_ST_PA + 1
              J_A_ST = IAREO(JJ_A_ST_PA  )
*. Address in standard order ( arrays transposed )
              J_AB_ST =  J_AB_ST_0 + J_A_ST*L_B_SM
C             J_AB_PA2 =  J_AB_PA2_0 +  J_PAAB_PA_OF2 + J_PA_A
              J_AB_PA2 =  J_AB_PA2 + 1
*. Standard => passive/active form
               CPA(J_AB_PA2) = SIGN*C(J_AB_ST)
             END DO
*            ^ End of loop over passive alpha strings
            ELSE 
*. Passive/Active => Standard form
             J_AB_PA2 =  J_AB_PA2_0 +  J_PAAB_PA_OF2 
C            DO J_PA_A = 1, L_PA_A_SM
             DO JJ_A_ST_PA = J_A_ST_PA+1,J_A_ST_PA+ L_PA_A_SM
*. Strings in standard format
C             J_A_ST_PA =  J_A_ST_PA + 1
              J_A_ST = IAREO(JJ_A_ST_PA  )
*. Address in standard order ( arrays transposed )
              J_AB_ST =  J_AB_ST_0 + J_A_ST*L_B_SM
C             J_AB_PA2 =  J_AB_PA2_0 +  J_PAAB_PA_OF2 + J_PA_A
              J_AB_PA2 = J_AB_PA2  + 1
*. Passive/Active => Standard form
C             C(J_AB_ST) = C(J_AB_ST) +  SIGN*CPA(J_AB_PA2)
              C(J_AB_ST) =  SIGN*CPA(J_AB_PA2)
             END DO
*            ^ End of loop over passive alpha strings
            END IF 
*             ^ End of switch passive/active <=  => Standard
            END DO
*           ^ End of loop over passive beta strings 
         END DO
         END DO
*        ^ End of loop over active strings of given sym
           END IF
*          ^ End if passive strings have nonvanishing dimension
           END DO
*          ^ End of loop over symmetry of passive alpha strings
       END DO
*      ^ End of loop over symmetry of Active alpha strings
      END DO
*     ^ End of loop over symmetry of product of alpha strings
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NSMST,ICSM = ', NSMST, ICSM
        WRITE(6,*) ' Number of passive strings per sym '
        CALL IWRTMA(NPAB,1,NSMST,1,NSMST) 
        WRITE(6,*) ' Number of active strings per sym '
        CALL IWRTMA(NAAB,1,NSMST,1,NSMST)
        IF( IWAY.EQ.1) THEN 
           WRITE(6,*) ' Standard => Active/Passive form '
        ELSE 
           WRITE(6,*) ' Active/Passive => Standard  form '
        END IF
        WRITE(6,*) ' Block in standard form '
        CALL WRTVH1(C,ICSM,NA,NB,NSMST,0)
        WRITE(6,*) ' Block in Passive/Active order '
        CALL WRTVH1(CPA,ICSM,NPAB,NAAB,NSMST,0)
C            WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
        WRITE(6,*) ' IBCPA2 array '
        CALL IWRTMA(IBCPA2,NSMST,NSMST,20,20)
      END IF
*
      END IF
*     ^ End of noreo = 0
*
      RETURN
      END
      SUBROUTINE RSBB2BEN(
     &           IATP,IBTP,NIA,NIB,JATP,JBTP,NJA,NJB,
     &           IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &           SB,CB,NOBPTS,MAXKOLD,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           XINT,NSMOB,NSMST,IUSEAB,CJRES,SIRES,SCLFAC,NTESTG,
     &           IUSE_PH,IPHGAS,ICSM,ISSM,C2,
     &           IREO_IA,IREO_IB,IREO_JA,IREO_JB,IUSE_HW,I12ORD)
*
* Version after Very New : Extremely New
* all symmetryblocks belonging to given TTS block treated
* simultaneously
*
*
* alpha-beta double excitation
* contribution from C blocks to S blocks
*
* Several C and S blocks with the same occupations, but
* different symmetries can be included
*
*
*. If IUSAB only half the terms are constructed
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
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
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT  : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + CJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks 
* October 96 : New routines for accessing annihilation information
*             Cleaned and shaved, only IROUTE = 3 option active
* October   97 : allowing for N-1/N+1 switch
* March 98  : Allows for splitting of strings into active and passive groups
* April 98 : Reodering loops for active/passive, allowing active and passive
*-April 99   division of alpha and beta strings, and simultaneous 
*            treatment of several symmetry blocks with the same type

*
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
*
      INTEGER NOBPTS(MXPNGAS,*)
*     
      INTEGER  NIA(*), NIB(*) 
      INTEGER  NJA(*), NJB(*)
*
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
*. Must hold excitations for all intermediate strings of given sym 
*. and all orbitals of given type
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION H(MXPTSOB*MXPTSOB)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*. Arrays for reorganization 
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
*. Arrays for active/passive division
      INTEGER IACIA(20),IPAIA(20),IACIB(20),IPAIB(20)
      INTEGER IACJA(20),IPAJA(20),IACJB(20),IPAJB(20)
*. Total number of passive/active strings per sym (ab combined)
      INTEGER NJPAAB(20),NJACAB(20) 
      INTEGER NIPAAB(20),NIACAB(20) 
*. Offset for given sym of active strings in PA block
      INTEGER IBCPA(20),IBCPA2(20,20),IBSPA(20),IBSPA2(20,20)
      INTEGER NOBPTS2(20)
*
      INTEGER NIBAC_S(8),NIBPA_S(8),IBREO_IB(8,8)
      INTEGER NJBAC_S(8),NJBPA_S(8),IBREO_JB(8,8)
      INTEGER NIAAC_S(8),NIAPA_S(8),IBREO_IA(8,8)
      INTEGER NJAAC_S(8),NJAPA_S(8),IBREO_JA(8,8)
C     INTEGER ICOFF3(8,8), ISOFF3(8,8)
*
      INTEGER KACGRP(20),IB_C_P_KA_J_JBA(20,20),IB_S_P_KA_I_IBA(20,20)
*. Dimension ^ : Number of string symmetries
      INTEGER IREO_IB(*),IREO_JB(*),IREO_IA(*),IREO_JA(*)
*.    ^ Dimension  : Largest number of strings of given type, all sym
      INCLUDE 'multd2h.inc'
      COMMON/KKKDUMMY/LEN_C2,LEN_S2
      COMMON/CMXCJ/MXCJ
     
      CALL QENTER('RS2B ')
      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' RSBB2BEN speaking '
        WRITE(6,*) ' ================= '
        WRITE(6,*)
        WRITE(6,*) ' IAOC and IBOC in RSBB2BEN'
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' JAOC and JBOC in RSBB2BEN'
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' TT-block of C   '
        CALL WRTVH1(CB,ICSM,NJB,NJA,NSMST,0)
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*
      IOPSM = MULTD2H(ICSM,ISSM)
C?    WRITE(6,*) ' IOPSM = ', IOPSM
*. Total length of C and sigma blocks
C     LEN_TT_BLOCK(ISM,NIA,NIB,NSMST)
      LEN_C = LEN_TT_BLOCK(ICSM,NJA,NJB,NSMST)
      LEN_S = LEN_TT_BLOCK(ISSM,NIA,NIB,NSMST)
*
      XCNORM = INPROD(CB,CB,LEN_C)
      IF(XCNORM.EQ.0.0D0) GOTO 9999
      IF(NTEST.GE.1000) WRITE(6,*) ' 1 : LEN_S = ', LEN_S
*     ^ I have been having problems with the value of LEN_S when  
*      using optimized code without this statement !!!!
*. Use passive/active splitting ?
      IUSE_PA = 1
*. i.e works only when passive/active seperation is on
*. Groups defining each supergroup
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
*.Types of SX that connects the two strings
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)           
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        IF(NTEST.GE.1000)
     &  WRITE(6,*) ' IJTYP, ITYP, JTYP =', IJTYP,ITYP,JTYP
*. Should N-1 or N+1 projection be used for alpha strings
        IJ_TYP(1) = ITYP
        IJ_TYP(2) = JTYP
        IJ_AC(1)  = 2
        IJ_AC(2) =  1
        NOP = 2
        CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,SIGNIJ)
        SIGNIJ2 = SIGNIJ*SCLFAC
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' SIGNIJ2 = ', SIGNIJ2
        IJAC = IJ_REO(2)
C?      WRITE(6,*) ' IJAC =', IJAC
*
        IF(IJ_REO(1).EQ.1) THEN
*. No reo
         IAC = 2
         JAC = 1
        ELSE
*. reordering
         IAC = 1
         JAC = 2
         IITYP = ITYP
         ITYP = JTYP
         JTYP = IITYP
*        ^ ITYP, JTYP reordered so ITYP is type of first index, 
*                                  JTYP is type of second index
        END IF
*
        DO 2000 KLTYP = 1, NKLTYP
           KTYP = KTP(KLTYP)
           LTYP = LTP(KLTYP)
*. Timings for different types of integrals, space 2 is assumed hole space
           IHSPC = 2
           NHSPC = 0
           IF(ITYP.EQ.IHSPC) NHSPC = 1
           IF(JTYP.EQ.IHSPC) NHSPC = NHSPC + 1
           IF(KTYP.EQ.IHSPC) NHSPC = NHSPC + 1
           IF(LTYP.EQ.IHSPC) NHSPC = NHSPC + 1
           IF(NHSPC.EQ.0) THEN 
             CALL QENTER('PPPP ')
           ELSE IF(NHSPC.EQ.1) THEN
             CALL QENTER('PPPH ')
           ELSE IF(NHSPC.EQ.2) THEN
             IF(ITYP.EQ.JTYP) CALL QENTER('PPHH ')
             IF(ITYP.NE.JTYP) THEN 
               IF(ITYP.EQ.KTYP) THEN 
                 CALL QENTER('PHPH ')
               ELSE
                 CALL QENTER('PHHP ')
               END IF
             END IF
           ELSE IF(NHSPC.EQ.3) THEN
             CALL QENTER('PHHH ')
           ELSE IF(NHSPC.EQ.4) THEN
             CALL QENTER('HHHH ')
           END IF
*

           IF(NTEST.GE.1000)
     &     WRITE(6,*) ' KLTYP, KTYP, LTYP =', KLTYP,KTYP,LTYP
*
           KL_TYP(1) = KTYP
           KL_TYP(2) = LTYP
           KL_AC(1)  = 2
           KL_AC(2) =  1
           NOP = 2
           CALL ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,SIGNKL)
           IF(KL_REO(1).EQ.1) THEN
*. No reo
            KAC = 2
            LAC = 1
           ELSE
*. reordering
            KAC = 1
            LAC = 2
            KKTYP = KTYP
            KTYP = LTYP
            LTYP = KKTYP
*          ^ KTYP, LTYP reordered so KTYP is type of first index, 
*                                    LTYP is type of second index
           END IF
           KLAC = KL_REO(2)
           IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' SIGNKL = ', SIGNKL
*
*. Division of alpha- and beta strings of C into active and passive parts 
*
           NOREOJA = 1
           NOREOJB = 1
           NOREOIA = 1
           NOREOIB = 1
*
           DO ISTSM = 1, NSMST
             IF(ISTSM.EQ.1) THEN
              IJAOFF = 1
              IJBOFF = 1
              IIAOFF = 1
              IIBOFF = 1
             ELSE
              IJAOFF = IJAOFF + NJA(ISTSM-1)
              IJBOFF = IJBOFF + NJB(ISTSM-1)
              IIAOFF = IIAOFF + NIA(ISTSM-1)
              IIBOFF = IIBOFF + NIB(ISTSM-1)
             END IF
           IF(JAC.NE.1.AND.JAC.NE.2) THEN
             WRITE(6,*) ' JAC output of range, 1 = ', JAC
             STOP 
           END IF
*. Split JA strings into active/passive part ( done for all symmetries )
*. JA strings of  sym ISTSM required
             KBSM = MULTD2H(ICSM,ISTSM)
             IF(NJB(KBSM).GT.0) THEN
               CALL REO_STR_SPGRP3(JASPGP,NGAS,ISTSM,NJA(ISTSM),2,
     &              IJ_TYP,NACJA,IACJA,NJAAC_S,NJAPA_S,
     &              IBREO_JA(1,ISTSM),IREO_JA(IJAOFF),SIGNJA,NACACJA,
     &              NOREOJA_L)
               IF(NOREOJA_L.EQ.0) NOREOJA = 0
             END IF

*. Split JB strings into active/passive part
             KASM = MULTD2H(ICSM,ISTSM)
             IF(NJA(KASM).GT.0) THEN
               CALL REO_STR_SPGRP3(JBSPGP,NGAS,ISTSM,NJB(ISTSM),2,
     &              KL_TYP,NACJB,IACJB,NJBAC_S,NJBPA_S,
     &              IBREO_JB(1,ISTSM),IREO_JB(IJBOFF),SIGNJB,NACACJB,
     &              NOREOJB_L)
               IF(NOREOJB_L.EQ.0) NOREOJB = 0
             END IF
*. Split Ia strings into active/passive part
             KBSM = MULTD2H(ISSM,ISTSM)
             IF(NIB(KBSM).GT.0) THEN
               CALL REO_STR_SPGRP3(IASPGP,NGAS,ISTSM,NIA(ISTSM),2,
     &              IJ_TYP,NACIA,IACIA,NIAAC_S,NIAPA_S,
     &              IBREO_IA(1,ISTSM),IREO_IA(IIAOFF),SIGNIA,NACACIA,
     &              NOREOIA_L)
               IF(NOREOIA_L.EQ.0) NOREOIA = 0
             END IF
*. Split Ib strings into active/passive part
             KASM = MULTD2H(ISSM,ISTSM)
             IF(NIA(KASM).GT.0) THEN
             CALL REO_STR_SPGRP3(IBSPGP,NGAS,ISTSM,NIB(ISTSM),2,
     &            KL_TYP,NACIB,IACIB,NIBAC_S,NIBPA_S,
     &            IBREO_IB(1,ISTSM),IREO_IB(IIBOFF),SIGNIB,NACACIB,
     &            NOREOIB_L)
               IF(NOREOIB_L.EQ.0) NOREOIB = 0
             END IF
           END DO
           NOREOC = NOREOJA*NOREOJB
           NOREOS = NOREOIA*NOREOIB
*
           CALL ACOP_SPGRP(NACJA,IACJA,JAC,JTYP,KACGRP,IJCODE)
C          IF(IJCODE.EQ.0) GOTO 2001
           IF(IJCODE.EQ.0) GOTO 1999
*          ^ Temp. changed for qenter/qexit timings
       IF(NTEST.GE.1000) WRITE(6,*) ' 2 : LEN_S = ', LEN_S
           SIGNJAB = SIGNJA*SIGNJB
           SIGNIAB = SIGNIA*SIGNIB
*. Reorganize C(Ja,Jb) as C(Ja_pa,Jb_pa,Ja_ac,Jb_ac)
           XDUM = 0.0D0
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' C(A,B) => C(P,A) '
           END IF
           IF(JAC.NE.1.AND.JAC.NE.2) THEN
             WRITE(6,*) ' JAC output of range, 2 = ', JAC
             STOP 
           END IF
           CALL CJPA_MS(CB,C2,1,NSMST,NJA,NJB,NJAAC_S,
     &                 NJBAC_S,NJAPA_S,NJBPA_S,IREO_JA,
     &                 IREO_JB,IBREO_JA,IBREO_JB,ICSM,SIGNJAB,
     &                 NJPAAB,NJACAB,IBCPA,IBCPA2,LEN_C,I12ORD,
     &                 NOREOC)
             IF(NOREOC.EQ.0) CALL COPVEC(C2,CB,LEN_C)
*. Reorganize S(Ia,IB) as C(Ia_pa,Ib_pa,Ia_ac,Ia_pa)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' S(A,B) => S(P,A) '
           END IF
             CALL CJPA_MS(SB,C2,1,NSMST,NIA,NIB,NIAAC_S,
     &                 NIBAC_S,NIAPA_S,NIBPA_S,IREO_IA,
     &                 IREO_IB,IBREO_IA,IBREO_IB,ISSM,SIGNIAB,
     &                 NIPAAB,NIACAB,IBSPA,IBSPA2,LEN_S,I12ORD,
     &                 NOREOS)
             IF(NOREOS.EQ.0) CALL COPVEC(C2,SB,LEN_S)
*
* Check if hardwired routines can be used for this sigma segment 
*
C?        WRITE(6,*) ' IUSE_HW = '  , IUSE_HW
C?        WRITE(6,*) ' NACACIA, NACACIB ', NACACIA, NACACIB
C?        WRITE(6,*) ' NACACJA, NACACJB ', NACACJA, NACACJB
          IF(IUSE_HW.EQ.1) THEN 
*
*. one electron in each active string
*
            IF(NACACIA.EQ.1.AND.NACACIB.EQ.1.AND.NACACJA.EQ.1
     &         .AND.NACACJB.EQ.1.AND.IJAC.EQ.2.AND.KLAC.EQ.2) THEN
              CALL SIGMA_AB_1111(SB,CB,ITYP,JTYP,KTYP,LTYP,
     &             NIPAAB,ICSM,ISSM,IBSPA,IBCPA,IBSPA2,IBCPA2,XINT) 
C     SIGMA_AB_1111(SB,CB,ITYP,JTYP,KTYP,LTYP,
C    &           NPA,ICSM,ISSM,IBS,IBC,IBS2,IBC2,XINT)
              GOTO 1950 
            END IF 
*
*. Atmost two electrons in each active strings
*
            IF(NACACIA.LE.2.AND.NACACIB.LE.2.AND.NACACJA.LE.2
     &         .AND.NACACJB.LE.2.AND.IJAC.EQ.2.AND.KLAC.EQ.2) THEN
               CALL SIGMA_AB_2222(SB,CB,
     &         NACIA,IACIA,NACIB,IACIB,NACJA,IACJA,NACJB,IACJB,
     &         ITYP,JTYP,KTYP,LTYP,NIAAC_S,NIBAC_S,NJAAC_S,NJBAC_S,
     &         NIPAAB,ICSM,ISSM,IBS,IBC,IBSPA2,IBCPA2,
     &         XINT,SIRES,CJRES,MXCJ)
               GOTO 1950
C      SIGMA_AB_2222(SB,CB,
C    &           NACIA,IACIA,NACIB,IACIB,NACJA,IACJA,NACJB,IACJB,
C    &           ITYP,JTYP,KTYP,LTYP,LACIA,LACIB,LACJA,LACJB,
C    &           NPA,ICSM,ISSM,IBS,IBC,IBS2,IBC2,XINT,SSCR,CSCR,LSCR)
            END IF
*
          END IF
*
*. Ka strings will be obtained as <Ia! a+/a ia |Ka><Ka!a/a+ ja!Ja>
*  Obtain Groups in Ka strings
           IF(JAC.NE.1.AND.JAC.NE.2) THEN
             WRITE(6,*) ' JAC output of range = ,3', JAC
             STOP 
           END IF
*. Notice Loop 1940 has been changed to loop over symmetry of Kastring,
*  several symmetries of ISM, JSM will be used ( as the Cblocks 
*  contain several symmetries
           DO 1940 KASM = 1, NSMST 
C                 NST_SPGRP2(NIGRP,IGRP,ISM_TOT,NSMST,NSTRIN,NDIST)
             CALL NST_SPGRP2(NACJA,KACGRP,KASM,NSMST,NKASTR,NKADIST)
             NKAEFF = NKASTR
             IF(NKAEFF.EQ.0) GOTO 1940
             IF(NTEST.GE.1000) 
     &       WRITE(6,*) ' KASM, NKASTR = ', KASM, NKASTR
*
*. Generate creation- or annihilation- mappings for all Ka strings
*
             IBLIST = 1
             JBLIST = 1
             DO IJSM = 1, NSMOB
              KAFRST = 1
              NI = NOBPTS(ITYP,IJSM)
              NJ = NOBPTS(JTYP,IJSM)
*
              JASM = MULTD2H(KASM,IJSM)
              IASM = MULTD2H(KASM,IJSM)
*. For operator connecting to |Ka_ac> and |Ja_ac> i.e. operator 2
C?            WRITE(6,*) ' IJSM and start of ADinfo' ,IJSM,JBLIST 
              CALL ADAST_GAS(IJSM,JTYP,NACJA,IACJA,JASM,
     &             I1(JBLIST),XI1S(JBLIST),NKASTR,IEND,IFRST,KFRST,
     &             KACT,SIGNIJ2,IJAC)
*. For operator connecting |Ka_ac> and |Ia_ac>, i.e. operator 1
              CALL ADAST_GAS(IJSM,ITYP,NACIA,IACIA,IASM,
     &             I3(IBLIST),XI3S(IBLIST),NKASTR,IEND,IFRST,KFRST,
     &             KACT,ONE,IJAC)
*
               IBLIST = IBLIST + NI*NKAEFF
               JBLIST = JBLIST + NJ*NKAEFF
             END DO
           IF(JAC.NE.1.AND.JAC.NE.2) THEN
             WRITE(6,*) ' JAC output of range = ,4', JAC
             STOP 
           END IF
*            ^ End of loop over IJSM for generating connections witj Ka
*
* Determine Length of Ka strings 
*
               L_C_P_J_JBA = 0
               DO JSM = 1, NSMOB
                NJ = NOBPTS(JTYP,JSM)
                DO JB_AC_SM = 1, NSMST
                 JA_AC_SM = MULTD2H(KASM,JSM)
                 JAB_AC_SM = MULTD2H(JA_AC_SM,JB_AC_SM)
                 JAB_PA_SM = MULTD2H(ICSM,JAB_AC_SM)
                 NNP =    NJPAAB(JAB_PA_SM)
                 NNJB = NJBAC_S(JB_AC_SM)
                 L_C_P_J_JBA = L_C_P_J_JBA + NJ*NNP*NNJB
                END DO
               END DO
               L_S_P_I_IBA = 0
               DO ISM = 1, NSMOB
                NI = NOBPTS(ITYP,ISM)
                DO IB_AC_SM = 1, NSMST
                 IA_AC_SM = MULTD2H(KASM,ISM)
                 IAB_AC_SM = MULTD2H(IA_AC_SM,IB_AC_SM)
                 IAB_PA_SM = MULTD2H(ISSM,IAB_AC_SM)
                 NNP =    NJPAAB(IAB_PA_SM)
                 NNIB = NIBAC_S(IB_AC_SM)
                 L_S_P_I_IBA = L_S_P_I_IBA + NI*NNP*NNIB
                END DO
               END DO
               LSCMAX = MAX(L_C_P_J_JBA,  L_S_P_I_IBA)
               IF(LSCMAX.EQ.0) GOTO 1940
               MAXK1 = MXCJ/LSCMAX 
               IF(MAXK1.EQ.0) THEN
                 WRITE(6,*) ' RSBB2BEN : MAXK1 = 0 '
                 STOP       ' RSBB2BEN : MAXK1 = 0 '
               END IF
               MAXK = MIN(NKAEFF,MAXK1)
*
*
*. Loop over batches of KA strings
*
             NKABTC = NKAEFF/MAXK   
             IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
             DO 1801 IKABTC = 1, NKABTC
               KABOT = (IKABTC-1)*MAXK + 1
               KATOP = MIN(KABOT+MAXK-1,NKAEFF)
               LKABTC = KATOP-KABOT+1
*. The game to play is 
* Sigma(P,Ka,i,Iab) = sum(jkl,Jab)<Iab!a+kb alb!Jab> (ij!kl)
*                                 C(P,Ka,j,Jab)
*. Obtain C(P,Ka,j,Jb) for Ka in batch, all symmetries of j, all sym of Jb
*. C(P,Ka,j,Jba) will be organized as a set of blocks, each with 
*  well-defined sym of j,Jb
               IIB_C_P_KA_J_JBA = 1
               DO JSM = 1, NSMOB
*               \/ Offset to start of excitations with sym JSM in I1
                IF(JSM.EQ.1) THEN
                  IBJORB = 1
                ELSE
                  IBJORB = IBJORB + NKAEFF*NOBPTS(JTYP,JSM-1)
                END IF
                NJ = NOBPTS(JTYP,JSM)
                DO JB_AC_SM = 1, NSMST
*. Info on input C(P,Jaa,Jab) block
                 JA_AC_SM = MULTD2H(KASM,JSM)
                 JAB_AC_SM = MULTD2H(JA_AC_SM,JB_AC_SM)
                 JAB_PA_SM = MULTD2H(ICSM,JAB_AC_SM)
*. Offset for C(P,Jaa,Jab) for Jaa, Jab with given sym
                 IIBCPA = IBCPA2(JA_AC_SM,JB_AC_SM)
                 NNP =    NJPAAB(JAB_PA_SM)
                 NNJB = NJBAC_S(JB_AC_SM)
*. Offset to start of blocks C(P,Ka,J,Jba)       
                 IB_C_P_KA_J_JBA(JSM,JB_AC_SM)= IIB_C_P_KA_J_JBA
* set up C(P,Ka,j,Jab) = sum(Jaa) <Ka!a ja!Jaa> C(P,Jab,Jaa)
*     for all symmetries of j and Ka in batch 
* C(P,Jaa,Jba) => C(P,Ka,j,Jba)
                 IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' JSM, JA_AC_SM, JB_AC_SM', 
     &                        JSM, JA_AC_SM, JB_AC_SM
                 END IF
C?               WRITE(6,*) ' IIBCPA, IIB_C_P_KA_J_JBA',
C?   &                        IIBCPA, IIB_C_P_KA_J_JBA
C?               WRITE(6,*) ' JSM and start of ADinfo',JSM,IBJORB
                 DO JJ = 1, NJ             
                   CALL GET_CPKAJJB(CB(IIBCPA),NJ,NJA,
     &                  CJRES(IIB_C_P_KA_J_JBA),LKABTC,NNJB,
     &                  JJ,I1(IBJORB-1 + (JJ-1)*NKASTR + KABOT),
     &                  XI1S (IBJORB-1 + (JJ-1)*NKASTR + KABOT),NNP)
                 END DO
                 IIB_C_P_KA_J_JBA = 
     &           IIB_C_P_KA_J_JBA + NNP*NJ*LKABTC*NNJB
                END DO
               END DO
               LEN_CJRES = IIB_C_P_KA_J_JBA - 1
*
               IF(LEN_CJRES.GT.MXCJ) THEN
                 WRITE(6,*) ' RSBB2BEN : LEN_CJRES > MXCJ '
                 WRITE(6,*) ' LEN_CJRES, MXCJ ', LEN_CJRES,MXCJ
                 WRITE(6,*) ' JAOC, JBOC '
                 CALL IWRTMA(JAOC,1,NGAS,1,NGAS) 
                 CALL IWRTMA(JBOC,1,NGAS,1,NGAS) 
                 WRITE(6,*) ' JTYP,IJAC ', JTYP , IJAC 
                 WRITE(6,*) ' LKABTC = ', LKABTC
C                IB_C_P_KA_J_JBA(JSM,JB_AC_SM)= IIB_C_P_KA_J_JBA
                 WRITE(6,*) '  IB_C_P_KA_J_JBA '
                 CALL IWRTMA(IB_C_P_KA_J_JBA,NSMST,NSMST,20,20)
                 STOP ' RSBB2BEN : LEN_CJRES > MXCJ '
               END IF
*
               IF(NTEST.GE.1000) THEN
                DO JJJSM = 1, NSMST
                  NOBPTS2(JJJSM) = NOBPTS(JTYP,JJJSM)
                END DO
                WRITE(6,*) ' The C(P,Ka,j,Jb_ac) array, KASM= ', KASM
                CALL WR_CPKJJB_MS(CJRES,NJPAAB,LKABTC,ICSM,KASM,
     &                           NOBPTS2,NJBAC_S,NSMST,IB_C_P_KA_J_JBA)
               END IF
* For sigma we will be using arrays S(P,Ka,I,Ib_ac) and S(PA,AC)
*Obtain offsets for these array 
*
               CALL IOFF_C_P_KA_J_JB(NIPAAB,ITYP,NOBPTS,MXPNGAS,
     &              LKABTC,KASM,ISSM,NSMOB,NIBAC_S,IB_S_P_KA_I_IBA,
     &              LEN_SIRES)
*
               IF(LEN_SIRES.GT.MXCJ) THEN
                 WRITE(6,*) ' RSBB2BEN : LEN_SIRES > MXCJ '
                 WRITE(6,*) ' LEN_SIRES, MXCJ ', LEN_SIRES,MXCJ
                 WRITE(6,*) ' IAOC, IBOC '
                 CALL IWRTMA(IAOC,1,NGAS,1,NGAS) 
                 CALL IWRTMA(IBOC,1,NGAS,1,NGAS) 
                 WRITE(6,*) ' ITYP,IJAC ', ITYP , IJAC 
                 STOP ' RSBB2BEN : LEN_SIRES > MXCJ '
               END IF
*
C?             WRITE(6,*) ' LEN_SIRES = ', LEN_SIRES 
               CALL SETVEC(SIRES,ZERO,LEN_SIRES)              
               FACS = 1.0D0
*
               DO 1931 KLSM = 1, NSMOB
               DO 1930 KSM = 1, NSMOB
*
                 IFIRST = 1
C                LSM = ADSXA(KSM,KLSM)
                 LSM = MULTD2H(KSM,KLSM)
                 IF(NTEST.GE.100) THEN
                 WRITE(6,*) ' KSM, LSM', KSM, LSM 
                 END IF
                 IF(LSM.EQ.0) GOTO 1930
                 NK = NOBPTS(KTYP,KSM)
                 NL = NOBPTS(LTYP,LSM)
*
                 KLAC = KL_REO(2)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                 IKORD = 0  
C                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
C                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
C    &           GOTO 1930
C                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
C    &           IKORD = 1
*
                 IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Loop over symmetries of active strings
                 JBSMMN = 1
                 JBSMMX = NSMST
                 DO JB_AC_SM = JBSMMN,JBSMMX
*. Is this symmetry active
                   II_ACTIVE = 0
                   DO JSM = 1, NSMST
                     JKLSM = MULTD2H(JSM,KLSM)
                     ISM = MULTD2H(JKLSM,IOPSM)
                     NI = NOBPTS(ITYP,ISM)
                     NJ = NOBPTS(JTYP,JSM)
                     J_JBACSM = MULTD2H(JSM,JB_AC_SM) 
                     KA_J_JBACSM = MULTD2H(KASM,J_JBACSM) 
                     IPSM = MULTD2H(ICSM,KA_J_JBACSM)
                     IF(NOBPTS(ITYP,ISM)*NOBPTS(JTYP,JSM).NE.0.AND.   
     &                  NJPAAB(IPSM).NE.0) II_ACTIVE = 1
                   END DO
*
                   IF(NTEST.GE.1000) WRITE(6,*) 'JB_AC_SM=',JB_AC_SM
*. Symmetry of active part of I_b string
                   IB_AC_SM = MULTD2H(KLSM,JB_AC_SM)
                   IF(NIBAC_S(IB_AC_SM)*NJBAC_S(JB_AC_SM).NE.0.AND.
     &                II_ACTIVE.EQ.1) THEN
*
                   CALL ADAST_GAS(LSM,LTYP,NACJB,IACJB,
     &                  JB_AC_SM,I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,
     &                  SIGNKL,KLAC)
                   IF(NKBSTR.NE.0) THEN
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                   CALL ADAST_GAS(KSM,KTYP,NACIB,IACIB,
     &                  IB_AC_SM,I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,
     &                  ONE,KLAC)
                   DO JSM = 1, NSMST
                     JKLSM = MULTD2H(JSM,KLSM)
                     ISM = MULTD2H(JKLSM,IOPSM)
                     NI = NOBPTS(ITYP,ISM)
                     NJ = NOBPTS(JTYP,JSM)
C?                   WRITE(6,*) ' ISM, JSM = ', ISM, JSM
                     J_JBACSM = MULTD2H(JSM,JB_AC_SM) 
                     KA_J_JBACSM = MULTD2H(KASM,J_JBACSM) 
                     IPSM = MULTD2H(ICSM,KA_J_JBACSM)
                     IF(NOBPTS(ITYP,ISM)*NOBPTS(JTYP,JSM).GT.0) THEN
                     IF(NJPAAB(IPSM).NE.0) THEN
* Fetch Integrals as (iop2 iop1 |  k l )
                     IXCHNG = 0
                     ICOUL = 1
                     ONE = 1.0D0
                     CALL GETINT(XINT,JTYP,JSM,ITYP,ISM,
     &                    KTYP,KSM,LTYP,LSM,IXCHNG,0,0,ICOUL,
     &                    1.0D0,1.0D0)
* S(P,Ka,j,Ib) = sum(k,l,Jb)<Ib!a+kba lb!Jb>C(P,Ka,j,Jb)*(ji!kl)
                     ISK_OFF = IB_S_P_KA_I_IBA(ISM,IB_AC_SM)
                     ICK_OFF = IB_C_P_KA_J_JBA(JSM,JB_AC_SM)
                     IF(NTEST.GE.1000) THEN
                      WRITE(6,*) ' ISK_OFF, ICK_OFF', ISK_OFF,ICK_OFF
                      WRITE(6,*) ' ISM,IB_AC_SM', ISM,IB_AC_SM 
                     END IF
C?                     WRITE(6,*) ' IPSM = ', IPSM
*
                       IROUTE = 3
                       CALL SKICKJ(SIRES(ISK_OFF),CJRES(ICK_OFF),
     &                      LKABTC*NJPAAB(IPSM),
     &                      NIBAC_S(IB_AC_SM),NJBAC_S(JB_AC_SM),
     &                      NKBSTR,XINT,NI,NJ,NK,NL,          
     &                      NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &                      FACS,IROUTE )
                     END IF
                   END IF
                   END DO
*                  ^ End of loop over JSM
                   END IF
*                  ^ End is NKBSTR.NE.0 
                   END IF
*                  ^ End of both active symmetries have nonvanishing dims.
C                  END DO
*                  ^ End of loop over IJSM
                 END DO
* .              ^ End of loop over JB_AC_SM
 1930          CONTINUE
*              ^ End of loop over KSM
 1931          CONTINUE 
*              ^ End of loop over KLSM
* Add contributions from S(P,Ka,i,Iba) to S(Ipa,Ipb)
               DO ISM = 1, NSMOB
*               \/ Offset to start of excitations with sym ISM in I2
                IF(ISM.EQ.1) THEN
                  IBIORB = 1
                ELSE
                  IBIORB = IBIORB + NKASTR*NOBPTS(ITYP,ISM-1)
                END IF
                NI = NOBPTS(ITYP,ISM)
                DO IB_AC_SM = 1, NSMST
                 IA_AC_SM = MULTD2H(KASM,ISM)
                 IAB_AC_SM = MULTD2H(IA_AC_SM,IB_AC_SM)
                 IAB_PA_SM = MULTD2H(ISSM,IAB_AC_SM)
*. Offset for S(P,Iaa,Iab) for Jaa, Jab with given sym
                 IIBSPA = IBSPA2(IA_AC_SM,IB_AC_SM)
                 NNP =    NIPAAB(IAB_PA_SM)
                 NNIB = NIBAC_S(IB_AC_SM)
C?               WRITE(6,*) ' NNP, NNIB ', NNP,NNIB
*. Scatter out from s(P,Ka,i,Ib)
* S(P,Jaa,Jba) <= S(P,Ka,j,Jba)
                 IF(NTEST.GE.1000) THEN
                  DO IIISM = 1, NSMST
                    NOBPTS2(IIISM) = NOBPTS(ITYP,IIISM)
                  END DO
                  IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' The S(P,Ka,i,Ib_ac) array '
                  CALL WR_CPKJJB_MS(SIRES,NIPAAB,LKABTC,ISSM,KASM,
     &                           NOBPTS2,NIBAC_S,NSMST,IB_S_P_KA_I_IBA)
                  END IF
                 END IF
                 DO II = 1, NI            
                 IIB_S_P_KA_I_IBA = IB_S_P_KA_I_IBA(ISM,IB_AC_SM)
C?                 WRITE(6,*) ' IIB_S_P_KA_I_IBA, IIBSPA',
C?   &                          IIB_S_P_KA_I_IBA, IIBSPA
                   CALL ADD_SPKAIIB(NNP,SB(IIBSPA),NI,NIA,
     &                  SIRES(IIB_S_P_KA_I_IBA),LKABTC,NNIB,
     &                  II,I3(IBIORB-1 + (II-1)*NKASTR + KABOT),
     &                  XI3S (IBIORB-1 + (II-1)*NKASTR + KABOT))
C                    ADD_SPKAIIB(NP,SB,NI,NIA,SKAIIB,NKA,NIB,I,ISCA,SSCA)
                 END DO
                END DO
               END DO
*              ^ End of loops over ISM, IB_AC_SM
 1801        CONTINUE 
*            ^ End of loop over Batches of Ka strings
 1940      CONTINUE
*          ^ End of loop over KASM
* S(Pa,Ac) => S(Ia,Ib)
 1950       CONTINUE
            SIGNPA = SIGNJA*SIGNJB*SIGNIA*SIGNIB
COLD        ZERO = 0.0D0
COLD        CALL SETVEC(C2,ZERO,LEN_S)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' S(P,A) => S(A,B) '
           END IF
            CALL CJPA_MS(C2,SB,2,NSMST,NIA,NIB,NIAAC_S,NIBAC_S,
     &           NIAPA_S,NIBPA_S,IREO_IA,IREO_IB,IBREO_IA,IBREO_IB,
     &           ISSM,SIGNIAB,NIPAAB,NIACAB,IBSPA,IBSPA2,LEN_S,I12ORD) 
            CALL COPVEC(C2,SB,LEN_S)
*. Reform also C from passive/active form back to Ja,Jb
COLD        ZERO = 0.0D0
COLD        CALL SETVEC(C2,ZERO,LEN_C)
            ONE = 1.0D0
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' C(P,A) => C(A,B) '
           END IF
            CALL CJPA_MS(C2,CB,2,NSMST,NJA,NJB,NJAAC_S,NJBAC_S,
     &           NJAPA_S,NJBPA_S,IREO_JA,IREO_JB,IBREO_JA,IBREO_JB,
     &           ICSM,SIGNJAB,NJPAAB,NJACAB,IBCPA,IBCPA2,LEN_C,I12ORD) 
            CALL COPVEC(C2,CB,LEN_C)
*
C1940      CONTINUE
*          ^ End of loop over KASM
 1999      CONTINUE
           IF(NHSPC.EQ.0) THEN 
             CALL QEXIT('PPPP ')
           ELSE IF(NHSPC.EQ.1) THEN
             CALL QEXIT('PPPH ')
           ELSE IF(NHSPC.EQ.2) THEN
             IF(ITYP.EQ.JTYP) CALL QEXIT('PPHH ')
             IF(ITYP.NE.JTYP) THEN 
                IF(ITYP.EQ.KTYP) THEN 
                  CALL QEXIT('PHPH ')
                ELSE
                  CALL QEXIT('PHHP ')
                END IF
             END IF
           ELSE IF(NHSPC.EQ.3) THEN
             CALL QEXIT('PHHH ')
           ELSE IF(NHSPC.EQ.4) THEN
             CALL QEXIT('HHHH ')
           END IF
 2000   CONTINUE
*           ^ End of loop over KLTYP
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Updated TT block of sigma '
        CALL WRTVH1(SB,ISSM,NIB,NIA,NSMST,0)
      END IF
      IF(NTEST.GE.200) WRITE(6,*) ' Return from RSBB2BEN '
CM    WRITE(6,*) ' Memcheck at end of RSBB2BEN '
CM    CALL MEMCHK
      CALL QEXIT('RS2B ')
      RETURN
      END
C     SUBROUTINE WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
*
* Write one-electron integrals with symmetry IVSM
* ISYM = 1 => Only lower triangular matrix included
*
* Jeppe Olsen, Jan. 1999
*
C     IMPLICIT REAL*8(A-H,O-Z)
*. General input
C     INTEGER NRPSM(NSMOB),NCPSM(NSMOB)
C     INCLUDE 'multd2h.inc'
*. Specific input
C     DIMENSION H(*)
*
C     IOFF = 1
C     DO ISM = 1, NSMOB
C       JSM = MULTD2H(ISM,IHSM)
C       NI = NRPSM(ISM)
C       NJ = NCPSM(JSM)
C       IF(ISYM.EQ.0.OR.ISM.GT.JSM) THEN
*. Complete block
C         WRITE(6,*) ' Block with symmetry ISM, JSM '
C         CALL WRTMAT(H(IOFF),NI,NJ,NI,NJ)
C         IOFF = IOFF + NI*NJ
C       ELSE IF (ISYM.EQ.1.AND.ISM.EQ.JSM) THEN
C         CALL PRSYM(H(IOFF),NI)
C         IOFF = IOFF * NI*(NI+1)/2
C       END IF
C     END DO
*
C     RETURN
C     END
      SUBROUTINE GET_CPKAJJB(CB,NJ,NJA,CPKAJJB,NKA,NJB,J,ISCA,SSCA,NP)
*
* Obtain for given orbital index j the gathered matrix
*
* C(P,Ka,j,Jb) = SSCA(Ka)C(P,Jb,ISCA(Ka))
*
* Note the occurance of a passive index
* Atlanta, Georgia, March24 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
       DIMENSION CB(NP,NJB,NJA), SSCA(*),ISCA(*)
*. Output
       DIMENSION CPKAJJB(NP,*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
         WRITE(6,*) ' From GET_CPKAJJB'
         WRITE(6,*)  ' ISCA AND SSCA for J=', J
         CALL WRTMAT(SSCA,1,NKA,1,NKA)
         CALL IWRTMA(ISCA,1,NKA,1,NKA)
      END IF
C     LBLK = 100
      LBLK = 40
      NBLK = NJB/LBLK
      IF(LBLK*NBLK.LT.NJB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NJB)
        ICONST = NKA*NJ 
        IADR0 =  (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
CT-commented out : Commented out when inserting passive/active
CT-                may be reinstated
CT      IF(ICEND.GT.ICOFF) THEN
*. Inner loop over JB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
*. Adress of C(Ka,j,Jb)
                IADR = IADR + ICONST
                DO IP = 1, NP
                  CPKAJJB(IP, IADR) = S*CB(IP, JB,IROW)
                END DO
              END DO
            ELSE  
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
                IADR = IADR + ICONST
                DO IP = 1, NP
                  CPKAJJB(IP,IADR) = 0.0D0          
                END DO
              END DO
            END IF
          END DO
CT      ELSE
*. No inner loop over JB
CT        DO KA  = 1, NKA
CT          IF(ISCA(KA).NE.0) THEN
CT            S = SSCA(KA)
CT            IROW = ISCA(KA)
CT            IADR = IADR0 + KA
*. Adress of C(Ka,j,Jb)
CT              IADR = IADR + ICONST
CT              CKAJJB(IADR) = S*CB(ICOFF,IROW)
CT          ELSE  
CT            IADR = IADR0 + KA
CT              IADR = IADR + ICONST
CT              CKAJJB(IADR) = 0.0D0          
CT          END IF
CT        END DO
CT      END IF
*       ^ End of test ICEND,ICOFF
      END DO
*
      RETURN
      END
      SUBROUTINE ADD_SPKAIIB(NP,SB,NI,NIA,SKAIIB,NKA,NIB,I,ISCA,SSCA)
*
* Update Transposed sigma block with contributions for given orbital index j 
* from the matrix S(Ka,i,Ib)
*
* S(P,Ib,Isca(Ka)) =  S(P,Ib,Isca(Ka)) + Ssca(Ka)*S(P,Ka,I,Ib)

* For efficient processing of alpha-beta loop
* Version containing passive index P
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
       DIMENSION SKAIIB(NP,*),SSCA(*),ISCA(*)
*. Input and Output
       DIMENSION SB(NP,NIB,NIA)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' From ADD_SPKAIIB '
         WRITE(6,*)  ' ISCA AND SSCA for I=', I
         CALL WRTMAT(SSCA,1,NKA,1,NKA)
         CALL IWRTMA(ISCA,1,NKA,1,NKA)
      END IF
*
C     LBLK = 100
      LBLK = 40
      NBLK = NIB/LBLK
      IF(LBLK*NBLK.LT.NIB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NIB)
        ICONST = NKA*NI
        IADR0 =  (I-1)*NKA+(ICOFF-1-1)*NKA*NI
*       \/ Switch has been commented out, second part not programmed for
*       passive strings
C       IF(ICEND.GT.ICOFF) THEN
*. Use form with Inner loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (I-1)*NKA+(ICOFF-1-1)*NKA*NI
              IADR = IADR0 + KA
              DO IB = ICOFF,ICEND
*. Adress of S(Ka,i,Ib)
                IADR = IADR + ICONST
                DO IP = 1, NP
                  SB(IP,Ib,IROW) = SB(IP,Ib,IROW)+S*SKAIIB(IP,IADR)
                END DO
              END DO
            END IF
          END DO
C       ELSE
*. Form with no loop over IB
C         DO KA  = 1, NKA
C           IF(ISCA(KA).NE.0) THEN
C             S = SSCA(KA)
C             IROW = ISCA(KA)
C             IADR = IADR0 + KA + ICONST
C             DO IB = ICOFF,ICEND
*. Adress of S(Ka,i,Ib)
C               IADR = IADR + ICONST
C               SB(ICOFF,IROW) = SB(ICOFF,IROW)+S*SKAIIB(IADR)
C             END DO
C           END IF
C         END DO
C       END IF
*       ^ End of test of ICOFF=ICEND
      END DO
*
      RETURN
      END
      SUBROUTINE ACOP_SPGRP(NNGRP,IGRP_IN,IAC,IGAS,IGRP_OUT,ICODE)
*
* A supergroup of NNGRP groups are given, IGRP
*
* Find supergroup obtained by annihilating or creating an      
* electron in GAS space IGAS
*
* Jeppe Olsen, Easter at Mag. vaegen 37D, April 99
*
* If the outgroup contains zero strings ICODE = 0 is returned
* If the outgroup is nonvanishing but cannot be generated 
* ICODE = -1  is returned
* A nonvanishing output group that could be generated is 
* indicated by a return code ICODE = 1
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      DIMENSION IGRP_IN(NNGRP)
*. Output
      DIMENSION IGRP_OUT(NNGRP)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  'ACOP_SPGRP speaking '
        WRITE(6,*)  '===================='
        WRITE(6,*)
        WRITE(6,*) ' IGAS, IAC, NNGRP ', IGAS,IAC, NNGRP
        WRITE(6,*) ' Input supergroup '
        CALL IWRTMA(IGRP_IN,1,NNGRP,1,NNGRP) 
      END IF
*. Group corresponding to GAS space IGAS: Must be present
      ICODE = 1
      IACGRP = 0
      DO JGRP = 1, NNGRP
       IF(IGSFGP(IGRP_IN(JGRP)).EQ.IGAS) IACGRP = JGRP
      END DO
*
      IF(IACGRP.EQ.0) THEN
        WRITE(6,*) ' ACOP_SPGRP in problems' 
        WRITE(6,*) ' Active gasspace not included '  
        WRITE(6,*) ' Active GASSPACE = ', IGAS
        STOP       ' ACOP_SPGRP in problems' 
      END IF
*. Number of electrons in active group, in and out
      IEL_IN = NELFGP(IGRP_IN(IACGRP))
      IF(IAC.EQ.1) THEN
        IEL_OUT = IEL_IN - 1
      ELSE IF(IAC.EQ.2)  THEN
        IEL_OUT = IEL_IN + 1
      ELSE
        WRITE(6,*) ' Unknown value of IAC = ', IAC
        STOP 'ACOP : Unknown value of IAC ' 
      END IF
C?    WRITE(6,*) ' ACOP : IACGRP, IAC,IEL_IN, IEL_OUT',
C?   &                    IACGRP, IAC,IEL_IN, IEL_OUT
*. Trivial vanishing group
      IF(IEL_OUT.LT.0.OR.IEL_OUT.GT.NGSOBT(IGAS)) THEN
        ICODE = 0
        GOTO 1000
      END IF
*. Find output group included 
      JGRP_OUT = 0
      DO JGRP = IBGPSTR(IGAS),IBGPSTR(IGAS)+NGPSTR(IGAS)-1
        IF(NELFGP(JGRP).EQ.IEL_OUT) JGRP_OUT = JGRP
      END DO
*
      IF(JGRP_OUT.EQ.0) THEN
        WRITE(6,*) ' ACOP_SPGP in problems '
        WRITE(6,*) ' Required output group not included '
        WRITE(6,*) ' IGAS and IEL_OUT ', IGAS,IEL_OUT
        WRITE(6,*) '  IAC, NNGRP ', IAC, NNGRP
        WRITE(6,*) ' Input supergroup '
        CALL IWRTMA(IGRP_IN,1,NNGRP,1,NNGRP) 
        STOP       ' ACOP_SPGP in problems '
      END IF
*
 1000 CONTINUE
*. Output supergroup
      CALL ICOPVE(IGRP_IN,IGRP_OUT,NNGRP)
      IGRP_OUT(IACGRP) = JGRP_OUT
*
      IF(NTEST.GE.100) THEN
        IF(ICODE.EQ.0) THEN
          WRITE(6,*) ' Trivial vanishAing output group'
        ELSE
          WRITE(6,*) ' output supergroup'
          CALL IWRTMA(IGRP_OUT,1,NNGRP,1,NNGRP)
        END IF
      END IF
*
      RETURN
      END
C     SIGMA_AB_1111(SB,CB,ITYP,JTYP,KTYP,LTYP,NIPAAB,XINT) 
      SUBROUTINE SIGMA_AB_1111(SB,CB,ITYP,JTYP,KTYP,LTYP,
     &           NPA,ICSM,ISSM,IBS,IBC,IBS2,IBC2,XINT)
*
* alpha- beta interaction between 4 singly occupied orbital spaces
*
* Jeppe Olsen, Magistratsvaegen 37D, July 6 - Nicklas watching Asterix 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'csm.inc'
*
* Sigma(P,ik) = Sum(jl)C(P,jl)[jl,ik]
*
* 
* ======
* Input
* ======
*
* C block in active/passive form
       DIMENSION CB(*)
*. Number of passive strings per sym
       INTEGER NPA(*)
*. Offsets to active alpha and betastrings with given sym
       INTEGER IBS(NSMST),IBC(NSMST)
       INTEGER IBS2(20,20)      ,IBC2(20,20)      
*
* ======
* Output
* ======
*. Updated sigma block in active/passive form
       DIMENSION SB(*)
       CALL QENTER('S1111 ')
*
* Sigma(P,ik) = Sum(jl)C(P,jl)[jl,ik]
      DO JLSM = 1, NSMST
        IPASM = MULTD2H(JLSM,ICSM)
        IKSM = MULTD2H(IPASM,ISSM)
        DO JSM = 1, NSMST
          LSM = MULTD2H(JLSM,JSM)
          DO ISM = 1, NSMST
            KSM = MULTD2H(IKSM,ISM)
*. Fetch integrals <jsm,lsm!ism,ksm> A
            IXCHNG = 0
            ICOUL = 0
            ONE = 1.0D0
            CALL GETINT(XINT,LTYP,LSM,KTYP,KSM,JTYP,JSM,ITYP,ISM,
     &                  IXCHNG,0,0,ICOUL,ONE,ONE) 
            ICOFF = IBC2(JSM,LSM)
            ISOFF = IBS2(ISM,KSM)
            NIK = NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
            NJL = NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM) 
            NP = NPA(IPASM)
            ONE = 1.0D0
            IF(NP*NIK*NJL.NE.0) THEN
              CALL MATML7(SB(ISOFF),CB(ICOFF),XINT,NP,NIK,
     &                    NP,NJL,NJL,NIK,ONE,ONE,0)
            END IF
          END DO
        END DO
      END DO
*
       CALL QEXIT('S1111 ')
      RETURN
      END
      SUBROUTINE SIGMA_AB_2222(SB,CB,
     &           NACIA,IACIA,NACIB,IACIB,NACJA,IACJA,NACJB,IACJB,
     &           ITYP,JTYP,KTYP,LTYP,LACIA,LACIB,LACJA,LACJB,
     &           NPA,ICSM,ISSM,IBS,IBC,IBS2,IBC2,XINT,SSCR,CSCR,LSCR)
*
* Alpha - beta contribution to sigma loop.
* Hardwired codes for case where there is atmost 2 electrons in 
* each of the active strings
*
* Jeppe Olsen, Magistratsvaegen 37 D, July 8-15 1999  
*              - Nicklas playing Age of Empires
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'stinf.inc'
*. Specific input 
      INTEGER NPA(*), IBS(*),IBC(*)
      INTEGER IBS2(20,20),IBC2(20,20)
      INTEGER IACIA(NACIA),IACIB(NACIB),IACJA(NACJA),IACJB(NACJB)
      INTEGER LACIA(NSMST),LACIB(NSMST),LACJA(NSMST),LACJB(NSMST)
      DIMENSION CB(*)
*
      INTEGER IACKA(2),IACKB(2)
*. Input and output  
      DIMENSION SB(*)
*. Scratch via argument list
      DIMENSION CSCR(LSCR),SSCR(LSCR),XINT(*)
*. Local scratch
C     INTEGER IB_IA(MXPNSMST,MXPNSMST), IB_JA(MXPNSMST,MXPNSMST)
C     INTEGER IB_IB(MXPNSMST,MXPNSMST) ,IB_JB(MXPNSMST,MXPNSMST)
*
      COMMON/SOMESCR/SCR(MXPTSOB*MXPTSOB*MXPTSOB*MXPTSOB) 
*
*. The show goes on as
*     C(P,Ja,Jb) => C(P,Ka,Kb,jl) 
*     S(P,Ka,Kb,ik) = C(P,Ka,Kb,jl)<jl!ik>
*     S(P,Ka,Kb,ik> => S(P,Ia,Ib)
*
* As we restrict ourselve to active strings with atmost one electron, 
* there is only one electron in the spectator strings and the 
*
* There are three relations between strings in C, S and K strings
*
* 1 : 1 elec in IJ string => K is vacuum string
* 2 : 2 elecs in different orbspaces => K contains on elec in one space
* 3 : 2 elecs in same orbspace      =>  K contains on elec in one space
*
*.  Offsets to I/J A/B strings with given sym in each type
*. Type of alpha and beta spectator strings
*
      CALL QENTER('S2222')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN  
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' SIGMA_AB_2222 entered '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,*) ' LACIA, LACIB : '
        CALL IWRTMA(LACIA,1,NSMST,1,NSMST)  
        CALL IWRTMA(LACIB,1,NSMST,1,NSMST)  
        WRITE(6,*) ' IBS2 '
        CALL IWRTMA(IBS2,NSMST,NSMST,20,20)
        WRITE(6,*) ' LACJA, LACJB : '
        CALL IWRTMA(LACJA,1,NSMST,1,NSMST)  
        CALL IWRTMA(LACJB,1,NSMST,1,NSMST)  
        WRITE(6,*) ' IBC2 '
        CALL IWRTMA(IBC2,NSMST,NSMST,20,20)
        WRITE(6,*) ' NPA  '
        CALL IWRTMA(NPA,1,NSMST,1,NSMST)
      END IF
*
      IOPSM = MULTD2H(ICSM,ISSM)
      CALL ACOP_SPGRP(NACJA,IACJA,1,JTYP,IACKA,IJCODE)
      CALL ACOP_SPGRP(NACJB,IACJB,1,LTYP,IACKB,KLCODE)
      IF(IJCODE.EQ.0.OR.KLCODE.EQ.0) GOTO 9999
      DO JLSM = 1, NSMST
       IF(NTEST.GE.100) WRITE(6,*) ' JLSM = ', JLSM
       IKSM = MULTD2H(IOPSM,JLSM)
       NJL = NOBPAIR(JLSM,JTYP,LTYP,0)
       NIK = NOBPAIR(IKSM,ITYP,KTYP,0)
       LSCR = MXPTSOB ** 4
*
       IF(NJL*NIK.GT. LSCR) THEN
         WRITE(6,*) 
     &   ' SIGMA_AB_2222 in problems : NJL*NIK > MXPTSOB ** 4'
         WRITE(6,*) ' Increase MXPTSOB '
         STOP
     &   ' SIGMA_AB_2222 in problems : NJL*NIK > MXPTSOB ** 4'
       END IF
*
*. Fetch all integrals (JL!IK) 
*
       IB_JL = 1
       DO JSM = 1, NSMST
        LSM = MULTD2H(JLSM,JSM)
        NJLS = NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM)
        IB_IK = 1
        DO ISM = 1, NSMST
         KSM = MULTD2H(ISM,IKSM)
         NIKS = NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
*. Fetch (I J K L) as <J L ! I K >
         IXCHNG = 0
         ICOUL = 0
         ONE = 1.0D0
         CALL GETINT(XINT,JTYP,JSM,ITYP,ISM,LTYP,LSM,KTYP,KSM,
     &                  IXCHNG,0,0,ICOUL,ONE,ONE) 
*
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Integral list as delivered '
           CALL WRTMAT(XINT,NIKS,NJLS,NIKS,NJLS) 
         END IF
*
         DO IK = 1, NIKS
          DO JL = 1, NJLS
            SCR((IK+IB_IK-1-1)*NJL + JL + IB_JL-1) = 
     &      XINT((IK-1)*NJLS + JL)
          END DO
         END DO
*        ^ End of loops over JL, IK
         IB_IK = IB_IK + NIKS
        END DO
*       ^ End of loop over ISM
        IB_JL = IB_JL + NJLS
       END DO
*      ^ End of loop over JSM
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Expanded Integral list '
         CALL WRTMAT(SCR,NIK,NJL,NIK,NJL) 
       END IF
*
       DO KASM = 1, NSMST
        DO KBSM = 1, NSMST
          KSM = MULTD2H(KASM,KBSM)
          JSM = MULTD2H(JLSM,KSM)
          IPSM = MULTD2H(ICSM,JSM)
          NP = NPA(IPSM)
*
          CALL NST_SPGRP2(NACJA,IACKA,KASM,NSMST,NKASTR,NKADIST)
          CALL NST_SPGRP2(NACJB,IACKB,KBSM,NSMST,NKBSTR,NKADIST)
C              NST_SPGRP2(NIGRP,IGRP,ISM_TOT,NSMST,NSTRIN,NDIST)
*
          LSCR_USE = MAX(NJL,NIK)*NKASTR*NKBSTR*NP
          IF(LSCR_USE.GT.LSCR) THEN
           WRITE(6,*) ' SIGMA_AB_2222 in trouble, LSCR_USE > LSCR '
           WRITE(6,*) ' LSCR_USE, LSCR = ', LSCR_USE,LSCR
           STOP       ' SIGMA_AB_2222 in trouble, LSCR_USE > LSCR ' 
          END IF
*
*. C(P,Kb,Ka,JL)
*
          IWAY = 1
          CALL CPA_CPKKJL_SP(CB,CSCR,IWAY,JLSM,NJL,JTYP,LTYP,
     &         KASM,KBSM,NKASTR,NKBSTR,NP,IBC2,
     &         NACJA,IACJA,NACJB,IACJB,LACJA,LACJB)
C         CPA_CPKKJL_SP(CPA,CPKKJL,IWAY,JLSM,NJL,JTYP,LTYP,
C    &        KASM,KBSM,NKASTR,NKBSTR,NP,IBC2,
C    &          NACJA,IACJA,NACJB,IACJB,LACJA,LACJB)
*
* S(P,Kb, Ka,IK) = C(P,Kb,Ka,JL)(JL!IK) 
*
          LPKAKB = NP*NKASTR*NKBSTR
          FACTORAB = 1.0D0
          FACTORC = 0.0D0
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Input C(PKaKb,jl) to matmult '
            CALL WRTMAT(CSCR,LPKAKB,NJL,LPKAKB,NJL)
          END IF
          CALL MATML7(SSCR,CSCR,SCR,LPKAKB,NIK,LPKAKB,NJL,
     &                NJL,NIK,FACTORC,FACTORAB,0) 
          IF(NTEST.GE.100) THEN 
            WRITE(6,*) ' Result of matrix multiply as S(PKaKb,ik)'
            CALL WRTMAT(SSCR,LPKAKB,NIK,LPKAKB,NIK)
          END IF
*
* S(P,Kb,Ka,IK) => S(P,Jb,Ja)
*
          IWAY = 2
C?        WRITE(6,*) ' SB(1) before call ', SB(1)
          CALL CPA_CPKKJL_SP(SB,SSCR,IWAY,IKSM,NIK,ITYP,KTYP,
     &         KASM,KBSM,NKASTR,NKBSTR,NP,IBS2,
     &         NACIA,IACIA,NACIB,IACIB,LACIA,LACIB)
C?        WRITE(6,*) ' SB(1) after call ', SB(1)
        END DO
       END DO
*      ^ End of Loop over KASM,KBSM
      END DO
*     ^ End of loop over JLSM
 9999 CONTINUE
*
      CALL QEXIT('S2222')
*
      RETURN
      END
      FUNCTION NOBPAIR(ITOTSM,ITYP,JTYP,IRES)
*
* Number of orbital pairs of total sym ITOTSM
* Type of orbitals is ITYP, JTYP
*
* If IRES = 1, then the orbitals are restricted so IORB.GE.JORB
*
* Jeppe Olsen, July 1999
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'multd2h.inc'
*
      NN = 0
      DO ISYM = 1, NSMOB
        JSYM = MULTD2H(ITOTSM,ISYM)
        NI = NOBPTS(ITYP,ISYM)
        NJ = NOBPTS(JTYP,JSYM)
        IF(IRES.EQ.0.OR.
     &     (IRES.EQ.1.AND.
     &      (ITYP.GT.JTYP.OR.(ITYP.EQ.JTYP.AND.ISYM.GT.JSYM)))) THEN
         NN = NN + NI*NJ
        ELSE IF(IRES.EQ.1.AND.ITYP.EQ.JTYP.AND.ISYM.EQ.JSYM) THEN
         NN = NN + NI*NJ
        END IF
      END DO
*
      NOBPAIR = NN
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Symmetry of orbital pairs  ', ITOTSM
       WRITE(6,*) ' Types of orbitals ', ITYP,JTYP
       WRITE(6,*) ' Dimension = ', NOBPAIR
      END IF
*
      RETURN 
      END
      SUBROUTINE IB_SPGRP_TT(IB_TT,L1,L2,MXPNSMST,NSMST)
*
* Offsets to types with given sym for supergroup 
* containing two types
*
* Jeppe Olsen, July 99
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER L1(NSMST),L2(NSMST)
      INCLUDE 'multd2h.inc'
*. Output
      INTEGER IB_TT(MXPNSMST,MXPNSMST) 
*. IB_TT(I1SM,I2SM) will give offset to block with type 1 strings of sym1
*                   and type 2 strings of type 2 with respect to start of 
*                   of this symmetry
*
      DO ITOTSM = 1, NSMST
        IB = 1
        DO I1SM = 1, NSMST
          I2SM = MULTD2H(I1SM,ITOTSM)
          IB_TT(I1SM,I2SM) = IB
          IB = IB + L1(I1SM)*L2(I2SM)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' IB_SPGRP_TT speaking, offset array for two type supergroup'
        WRITE(6,*)
        CALL IWRTMA(IB_TT,NSMST,NSMST,MXPNSMST,MXPNSMST)
      END IF
*
      RETURN
      END
      SUBROUTINE CPA_CPKKJL_SP(CPA,CPKKJL,IWAY,JLSM,NJL,JTYP,LTYP,
     &           KASM,KBSM,NKASTR,NKBSTR,NP,IBC2,
     &           NACJA,IACJA,NACJB,IACJB,LACJA,LACJB)
*
* C(P,Ja,Jb) <=> C(P,Ka,Kb,JL) for given sym of P, Ka, Kb, j*l
*
* Special case routine for J-strings containing atmost two electrons
* (may be used/extended to k-strings defined by just two types)
*
*
* Iway = 1
* C(P,Kb,Ka,jl) = sum(Ja,Jb) <Ja!a+ja!Ka> <Jb!a+lb!Kb> C(P,Jb,Ja)
*
* IWay = 2
* C(P,Jb,Ja) = C(P,Jb,Ja) + sum(Ka,Kb) <Ja!a+ja!Ka><Jb!a+lb!Kb> C(P,Ka,Kb,jl)
*
* Jeppe Olsen, July 1999
*
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
*. Specific input 
      INTEGER  IBC2(20,20)
      INTEGER IACJA(NACJA),IACJB(NACJB)
      INTEGER LACJA(NSMST),LACJB(NSMST)
*. Local scratch
      INTEGER IACKA(2),IACKB(2)
*. input and output  
      DIMENSION CPKKJL(*)
      DIMENSION CPA(*)
*. Local scratch
      INTEGER IB_JA(MXPNSMST,MXPNSMST), IB_JB(MXPNSMST,MXPNSMST)
      INTEGER IB_KA(MXPNSMST,MXPNSMST), IB_KB(MXPNSMST,MXPNSMST)
*
*
* As we restrict ourselves to active strings with atmost two electrons, 
* there is only one electron in the spectator strings  
*
* There are three relations between strings in C, S and K strings
*
* 1 : 1 elec in IJ string => K is vacuum string
* 2 : 2 elecs in different orbspaces => K contains on elec in one space
* 3 : 2 elecs in same orbspace      =>  K contains on elec in one space
*
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) '  CPA_CPKKJL_SP entered '
       WRITE(6,*) ' ======================='
       WRITE(6,*)
       WRITE(6,*) ' IACJA, IACJB :'
       CALL IWRTMA(IACJA,1,NACJA,1,NACJA)
       CALL IWRTMA(IACJB,1,NACJB,1,NACJB)
       WRITE(6,*) ' JTYP,LTYP =', JTYP,LTYP
       WRITE(6,*) ' KASM, KBSM' , KASM, KBSM
      END IF
*
      CALL ACOP_SPGRP(NACJA,IACJA,1,JTYP,IACKA,JCODE)
      CALL ACOP_SPGRP(NACJB,IACJB,1,LTYP,IACKB,LCODE)
*. Sign to bring operators in place
      IF(NACJA.EQ.1) THEN
        SIGNJL = 1.0D0
      ELSE
        IF(JTYP.EQ.IGSFGP(IACKA(1))) THEN
         SIGNJL = 1.0D0
        ELSE
         SIGNJL = (-1)**(NELFGP(IACKA(1)))
        END IF
      END IF
      IF(NACJB.EQ.2) THEN
        IF(LTYP.EQ.IGSFGP(IACKB(2))) THEN
         SIGNJL = SIGNJL*(-1)**(NELFGP(IACKB(1)))
        END IF
      END IF
         
*. Mapping to start wrt to start of sym for given symmetrypairs
      IF(NACJA.EQ.2) THEN
        KA1TP = IACKA(1)
        KA2TP = IACKA(2)
        CALL IB_SPGRP_TT(IB_KA,NSTFSMGP(1,KA1TP),NSTFSMGP(1,KA2TP),
     &                   MXPNSMST,NSMST)
      END IF
*
      IF(NACJB.EQ.2) THEN
        KB1TP = IACKB(1)
        KB2TP = IACKB(2)
        CALL IB_SPGRP_TT(IB_KB,NSTFSMGP(1,KB1TP),NSTFSMGP(1,KB2TP),
     &                   MXPNSMST,NSMST)
      END IF
*
*. Creation mappings assumed
      JAC = 2
      LAC = 2
*. space in active strings corresponding to J,L
      IF(IGSFGP(IACJA(1)).EQ.JTYP) THEN
       JACT_12 = 1
      ELSE
       JACT_12 = 2
      END IF
      JPAS_12 = 0
      IF(NACJA.EQ.2) THEN
        IF(JACT_12.EQ.1) JPAS_12 = 2
        IF(JACT_12.EQ.2) JPAS_12 = 1
      END IF
      IF(IGSFGP(IACJB(1)).EQ.LTYP) THEN
       LACT_12 = 1
      ELSE
       LACT_12 = 2
      END IF
      LPAS_12 = 0
      IF(NACJB.EQ.2) THEN
        IF(LACT_12.EQ.1) LPAS_12 = 2
        IF(LACT_12.EQ.2) LPAS_12 = 1
      END IF
* The string types modified under process
      KAGRP_ACT = IACKA(JACT_12)
      JAGRP_ACT = IACJA(JACT_12)
      KBGRP_ACT = IACKB(LACT_12)
      JBGRP_ACT = IACJB(LACT_12)
*. The passive groups
      KAGRP_PAS = 0
      IF(JPAS_12.NE.0) KAGRP_PAS = IACJA(JPAS_12)
      JAGRP_PAS = KAGRP_PAS
      KBGRP_PAS = 0
      IF(LPAS_12.NE.0) KBGRP_PAS = IACJB(LPAS_12)
C     WRITE(6,*) ' LPAS_12, KBGRP_PAS', LPAS_12,KBGRP_PAS
      JBGRP_PAS = KBGRP_PAS
*
      NEL_KAGRP_ACT = NELFGP(KAGRP_ACT)
      NEL_KBGRP_ACT = NELFGP(KBGRP_ACT)
*. Offsets to types of given sym in J A/B
      JA1TP = IACJA(1)
      IF(NACJA.EQ.1) THEN
        JA2TP = 0 
      ELSE
        JA2TP = IACJA(2)
        CALL IB_SPGRP_TT(IB_JA,NSTFSMGP(1,JA1TP),NSTFSMGP(1,JA2TP),
     &                   MXPNSMST,NSMST)
C            IB_SPGRP_TT(IB_TT,L1,L2,MXPNSMST,NSMST)
      END IF
*
      JB1TP = IACJB(1)
      IF(NACJB.EQ.1) THEN
        JB2TP = 0 
      ELSE
        JB2TP = IACJB(2)
        CALL IB_SPGRP_TT(IB_JB,NSTFSMGP(1,JB1TP),NSTFSMGP(1,JB2TP),
     &                   MXPNSMST,NSMST)
      END IF
*. Are group mappings in expanded or compact form 
      IF(JAC.EQ.1.AND.ISTAC(KAGRP_ACT,2).EQ.0) THEN
        JEC = 2
        LR_J = NEL_KAGRP_ACT
      ELSE 
        JEC = 1
        LR_J = NOBPT(JTYP)
      END IF
*
      IF(LAC.EQ.1.AND.ISTAC(KBGRP_ACT,2).EQ.0) THEN
        LEC = 2
        LR_L = NEL_KBGRP_ACT
      ELSE 
        LEC = 1
        LR_L = NOBPT(LTYP)
      END IF
*. Number and symmetry of active and passive strings 
      IF(NACJA.EQ.1) THEN
        KA_ACT_SYM = KASM
        KA_PAS_SYM = 0
      ELSE
        IF(NEL_KAGRP_ACT.NE.0) THEN
          KA_ACT_SYM = KASM
          KA_PAS_SYM = 1
        ELSE 
          KA_ACT_SYM = 1
          KA_PAS_SYM = KASM
        END IF
      END IF
*
      IF(NACJB.EQ.1) THEN      
        KB_ACT_SYM = KBSM
        KB_PAS_SYM = 0
      ELSE
        IF(NEL_KBGRP_ACT.NE.0) THEN
          KB_ACT_SYM = KBSM
          KB_PAS_SYM = 1
        ELSE 
          KB_ACT_SYM = 1
          KB_PAS_SYM = KBSM
        END IF
      END IF
*
      N_KA_ACT = NSTFSMGP(KA_ACT_SYM,KAGRP_ACT)
      IF(KAGRP_PAS.EQ.0) THEN
        N_KA_PAS = 1
      ELSE
        N_KA_PAS = NSTFSMGP(KA_PAS_SYM,KAGRP_PAS)
      END IF
*
      N_KB_ACT = NSTFSMGP(KB_ACT_SYM,KBGRP_ACT)
      IF(KBGRP_PAS.EQ.0) THEN
        N_KB_PAS = 1
      ELSE
        N_KB_PAS = NSTFSMGP(KB_PAS_SYM,KBGRP_PAS)
      END IF
*. Offsets of K strings with respect to start of type
      IB_KAGRP_ACT = ISTFSMGP(KA_ACT_SYM,KAGRP_ACT)
      IB_KBGRP_ACT = ISTFSMGP(KB_ACT_SYM,KBGRP_ACT)
*
* Addresses of K strings with respect to start of sym from parts 
*
* When going to compound index, remember that the rightmost string  
* changes faster. 
* The index of a string (IST1, IST2) is therefore (IST1-1)*NST2 + IST2
* the strin
*  KA = IKA_ACT*KA_ACT + IKA_PAS*KA_PAS + IKA0 -1 
      IF(NACJA.EQ.1) THEN
*. No passive types
        IKA_PAS = 0
        IKA_ACT = 1
        IKA0    = 1
C     ELSE IF(JACT_12.EQ.2) THEN
      ELSE IF(JACT_12.EQ.1) THEN
* ( strings run with rightmost index as the inner index )
*. (KACT,KPAS)
        IKA_PAS = 1
        IKA_ACT = N_KA_PAS
C       IKA0 = IB_KA(KA_PAS_SYM,KA_ACT_SYM)-N_KA_PAS
        IKA0 = IB_KA(KA_ACT_SYM,KA_PAS_SYM)-N_KA_PAS
      ELSE
*. (KPAS, KACT)
        IKA_PAS = N_KA_ACT
        IKA_ACT = 1
        IKA0 = IB_KA(KA_PAS_SYM,KA_ACT_SYM)-N_KA_ACT
      END IF
*  KB = IKB_ACT*KB_ACT + IKB_PAS*KB_PAS + IKB0 -1 
      IF(NACJB.EQ.1) THEN
*. No passive types
        IKB_PAS = 0
        IKB_ACT = 1
        IKB0    = 1
      ELSE IF(LACT_12.EQ.1) THEN
*. (KACT, KPAS) runs as normal matrix (KPAS,KACT)
        IKB_PAS = 1
        IKB_ACT = N_KB_PAS
C       IKB0 = IB_KB(KB_PAS_SYM,KB_ACT_SYM)-N_KB_PAS
        IKB0 = IB_KB(KB_ACT_SYM,KB_PAS_SYM)-N_KB_PAS
      ELSE
*. (KPAS, KACT)
        IKB_PAS = N_KB_ACT
        IKB_ACT = 1
C       IKB0 = IB_KB(KB_ACT_SYM,KB_PAS_SYM)-N_KB_ACT
        IKB0 = IB_KB(KB_PAS_SYM,KB_ACT_SYM)-N_KB_ACT
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IKB_PAS, IKB_ACT, IKB0', IKB_PAS, IKB_ACT, IKB0
      END IF
*
      IBJL = 1
      DO JSM = 1, NSMST
       LSM = MULTD2H(JLSM,JSM) 
* 
       NJORB = NOBPTS(JTYP,JSM)
       NLORB = NOBPTS(LTYP,LSM)
       IBJORB = IOBPTS(JTYP,JSM)-IOBPTS(JTYP,1)+1
       IBLORB = IOBPTS(LTYP,LSM)-IOBPTS(LTYP,1)+1
C?     WRITE(6,*) ' IBJORB, IBLORB = ', IBJORB,IBLORB
*
       JASM = MULTD2H(KASM,JSM)
       JBSM = MULTD2H(KBSM,LSM)
       N_JA_TOT = LACJA(JASM)
       N_JB_TOT = LACJB(JBSM)
*. Symmetry of active and passive part of JA, JB 
       JA_ACT_SYM =  MULTD2H(JSM,KA_ACT_SYM)
       IF(NACJA.EQ.2) JA_PAS_SYM =  KA_PAS_SYM
       JB_ACT_SYM = MULTD2H(LSM, KB_ACT_SYM)
       IF(NACJB.EQ.2) JB_PAS_SYM  = KB_PAS_SYM
*. Dimensions of these J blocks 
       N_JA_ACT  = NSTFSMGP(JA_ACT_SYM,JAGRP_ACT)
       N_JA_PAS  = N_KA_PAS
       N_JB_ACT  = NSTFSMGP(JB_ACT_SYM,JBGRP_ACT)
       N_JB_PAS = N_KB_PAS
*
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' KAGRP_ACT, N_KA_ACT = ',  KAGRP_ACT, N_KA_ACT
         WRITE(6,*) ' KBGRP_ACT, N_KB_ACT = ',  KBGRP_ACT, N_KB_ACT
         WRITE(6,*) ' KAGRP_PAS, N_KA_PAS = ',  KAGRP_PAS, N_KA_PAS
         WRITE(6,*) ' KBGRP_PAS, N_KB_PAS = ',  KBGRP_PAS, N_KB_PAS
      END IF
*  JA = IJA_ACT*JA_ACT + IJA_PAS*JA_PAS + IJA0 - 1 
       IF(NACJA.EQ.1) THEN
*. No passive types
         IJA_PAS = 0
         IJA_ACT = 1
         IJA0    = 1
C      ELSE IF (JACT_12.EQ.2) THEN
       ELSE IF (JACT_12.EQ.1) THEN
*. (JACT, JPAS)
         IJA_PAS = 1
         IJA_ACT = N_JA_PAS
C        IJA0    = IB_JA(JA_PAS_SYM,JA_ACT_SYM)-N_JA_PAS
         IJA0    = IB_JA(JA_ACT_SYM,JA_PAS_SYM)-N_JA_PAS
       ELSE
*. (JPAS, JACT)
         IJA_PAS = N_JA_ACT
         IJA_ACT = 1
C        IJA0    = IB_JA(JA_ACT_SYM,JA_PAS_SYM)-N_JA_ACT
         IJA0    = IB_JA(JA_PAS_SYM,JA_ACT_SYM)-N_JA_ACT
       END IF
*  JB = IJB_ACT*JB_ACT + IJB_PAS*JB_PAS + IJB0 - 1 
       IF(NACJB.EQ.1) THEN
*. No passive types
         IJB_PAS = 0
         IJB_ACT = 1
         IJB0    = 1
       ELSE IF(LACT_12.EQ.1) THEN
*. (JACT, JPAS)
         IJB_PAS = 1
         IJB_ACT = N_JB_PAS
C        IJB0    = IB_JB(JB_PAS_SYM,JB_ACT_SYM)-N_JB_PAS
         IJB0    = IB_JB(JB_ACT_SYM,JB_PAS_SYM)-N_JB_PAS
       ELSE
*. (JPAS, JACT)
         IJB_PAS = N_JB_ACT
         IJB_ACT = 1
C        IJB0    = IB_JB(JB_ACT_SYM,JB_PAS_SYM)-N_JB_ACT
         IJB0    = IB_JB(JB_PAS_SYM,JB_ACT_SYM)-N_JB_ACT
       END IF
*. Offset for active J strings w.r.t start of supergroup
       IBB_JA = ISTFSMGP(JA_ACT_SYM,JAGRP_ACT)
       IBB_JB = ISTFSMGP(JB_ACT_SYM,JBGRP_ACT)
C      IBB_JA = ISTFSMGP(JASM,JAGRP_ACT)
C      IBB_JB = ISTFSMGP(JBSM,JBGRP_ACT)
*. CPA : Offset to block with given sym of P, Ja, Jb
       IB_CPA = IBC2(JASM,JBSM)
*. CPKAKBJL : Start of given sym of J and L
C      IB_CPKKJL = (IBJL-1)*NKASTR*NKBSTR+1
       IB_CPKKJL = (IBJL-1)*NP*NKASTR*NKBSTR+1
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' IB_CPA, IB_CPKKJL', IB_CPA, IB_CPKKJL
         WRITE(6,*) ' JASM, JBSM, N_JA_TOT, N_JB_TOT ', 
     &                JASM, JBSM, N_JA_TOT, N_JB_TOT 
       END IF
       SIGN0 = SIGNJL
*. And call routine to do the expansion 
       CALL 
     & CPA_CPKKJL_SP2(CPA(IB_CPA),CPKKJL(IB_CPKKJL),SIGN0,IWAY,
     & WORK(KSTSTM(KAGRP_ACT,1)),WORK(KSTSTM(KAGRP_ACT,2)),LR_J,
     & 2,N_JA_TOT,
     & WORK(KSTSTM(KBGRP_ACT,1)),WORK(KSTSTM(KBGRP_ACT,2)),LR_L,
     & 2,N_JB_TOT,NJORB,NLORB,IBJORB,IBLORB,NP,
     & IB_KAGRP_ACT,N_KA_ACT,IB_KBGRP_ACT,N_KB_ACT,
     & N_KA_PAS,N_KB_PAS,
     & IBB_JA,IBB_JB,
     & IKA_ACT,IKA_PAS,IKA0,IKB_ACT,IKB_PAS,IKB0,
     & IJA_ACT,IJA_PAS,IJA0,IJB_ACT,IJB_PAS,IJB0)
       IBJL = IBJL + NJORB*NLORB 
      END DO
*     ^ End of loop over JSM
      RETURN
      END
*
      SUBROUTINE CPA_CPKKJL_SP2(CPA,CPKKJL,SIGN0,IWAY,
     &           IKJAAO,IKJAAS,LKAA,IKJAAC,NJA,
     &           IKJBAO,IKJBAS,LKBA,IKJBAC,NJB,
     &           NJ,NL,IBJ,IBL,NP,
     &           IBKAA,NKAA,IBKBA,NKBA,NKAP,NKBP,
     &           IBJAA,IBJBA,
     &           IKA_ACT,IKA_PAS,IKA0,IKB_ACT,IKB_PAS,IKB0,  
     &           IJA_ACT,IJA_PAS,IJA0,IJB_ACT,IJB_PAS,IJB0 )
*
* C(P,Jb,Ja) <=> C(P,Kb,Ka,jl) Given sym of Ja, Jb, Ka,Kb,j,l
*
* Hardwired case for K strings only containing a single active type
*
* NOTE : 
*. CPA should be called with start of given symmetryblock as first element
*. CPKKJL should be called with start of given sym of JL as first element.
*
* Iway = 1
* C(P,Kb,Ka,jl) = sum(Ja,Jb) <Ja!a+ja!Ka> <Jb!a+lb!Kb> C(P,Jb,Ja)
*
* IWay = 2
* C(P,Jb,Ja) = C(P,Jb,Ja) + sum(Ka,Kb) <Ja!a+ja!Ka><Jb!a+lb!Kb> C(P,Ka,Kb,jl)
*
* The active Ja and Jb strings can consist of upto two stringtypes 
* an active and an passive part.
*
* Only creation mappings implemented 
*
* Jeppe Olsen, July 1999
*
      INCLUDE 'implicit.inc' 
*.Kaa => Jaa mappings
      INTEGER IKJAAO(LKAA,*), IKJAAS(LKAA,*)
*.Kba => Jba mappings
      INTEGER IKJBAO(LKBA,*), IKJBAS(LKBA,*)
*. Input and output    
      DIMENSION CPA(*),CPKKJL(*)
*
C     WRITE(6,*) ' Memtest at start of SP2 '
C     CALL MEMCHK
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CPA_CPKKJL_SP2 entered '
        WRITE(6,*) ' NKAA, NKBA : ', NKAA, NKBA
        WRITE(6,*) ' NKAP, NKBP : ', NKAP, NKBP
        WRITE(6,*) ' IKJAAC = ', IKJAAC
        WRITE(6,*) ' CPA(1) = ', CPA(1)
        WRITE(6,*) ' 1 : IKB_ACT, IKB_PAS, IKB0',
     &                   IKB_ACT, IKB_PAS, IKB0
        WRITE(6,*) ' NJA, NJB ', NJA, NJB
*
      END IF
*
      NKA = NKAA * NKAP
      NKB = NKBA * NKBP
*
      DO KAASTR = 1,NKAA
       KAASTR_ABS = KAASTR + IBKAA - 1
       DO JORB = 1, NJ 
        JORB_ABS = JORB + IBJ - 1
        KA_ACTIVE = 0
C?      WRITE(6,*) ' JORB, JORB_ABS = ', JORB, JORB_ABS
        IF(IKJAAC.EQ.2.AND.IKJAAO(JORB_ABS,KAASTR_ABS).GT.0) THEN
*. Creation is nonvanishing
          KA_ACTIVE = 1
          IF(IKJAAS(JORB_ABS,KAASTR_ABS) .GT. 0 ) THEN
            SIGNJ = 1.0D0
            JAASTR_ABS = IKJAAS(JORB_ABS,KAASTR_ABS)
          ELSE
            SIGNJ = -1.0D0
            JAASTR_ABS = -IKJAAS(JORB_ABS,KAASTR_ABS)
          END IF
        END IF
*. Jaa string relative to start of symmetry
C?      WRITE(6,*) ' KAASTR, JORB, KA_ACTIVE ', KAASTR,JORB,KA_ACTIVE
        JAASTR = JAASTR_ABS - IBJAA + 1
C?      IF(KA_ACTIVE.EQ.1) WRITE(6,*) ' JAASTR_ABS, JAASTR',
C?   &  JAASTR_ABS,JAASTR
        IF(KA_ACTIVE.EQ.1.AND.NTEST.GE.100) 
     &  WRITE(6,*) ' JORB,KAASTR, JAASTR', JORB,KAASTR,JAASTR
        DO KBASTR = 1,NKBA
         KBASTR_ABS = KBASTR + IBKBA - 1
         DO LORB = 1, NL 
          LORB_ABS = LORB + IBL - 1
          KB_ACTIVE = 0
          IF(IKJBAC.EQ.2.AND.IKJBAO(LORB_ABS,KBASTR_ABS).GT.0) THEN
*. Creation is nonvanishing
            KB_ACTIVE = 1
            IF(IKJBAS(LORB_ABS,KBASTR_ABS) .GT. 0 ) THEN
              SIGNL = 1.0D0
              JBASTR_ABS = IKJBAS(LORB_ABS,KBASTR_ABS)
            ELSE
              SIGNL = -1.0D0
              JBASTR_ABS =-IKJBAS(LORB_ABS,KBASTR_ABS)
            END IF
          END IF
*. Jba string relative to start of symmetry
          IF(KB_ACTIVE.EQ.1) 
     &    JBASTR = JBASTR_ABS - IBJBA + 1
          IF(KB_ACTIVE.EQ.1.AND.NTEST.GE.100) 
     &    WRITE(6,*) ' LORB,KBASTR,JBASTR_ABS,JBASTR,IBJBA',
     &                 LORB,KBASTR,JBASTR_ABS,JBASTR,IBJBA
C         WRITE(6,*) ' KBASTR, KBASTR_ABS, LORB', 
C    &                 KBASTR, KBASTR_ABS, LORB
C         IF(KB_ACTIVE.EQ.1) WRITE(6,*) ' JBASTR, JBASTR_ABS',
C    &    JBASTR, JBASTR_ABS
*. Loop over passive parts of Ka and Kb
          K_ACTIVE = KA_ACTIVE*KB_ACTIVE
          IF(K_ACTIVE.EQ.1) SIGN = SIGNJ*SIGNL*SIGN0
          DO KAP = 1, NKAP
          DO KBP = 1, NKBP
*  JA = IJA_ACT*JA_ACT + IJA_PAS*JA_PAS + IJA0 - 1 
           KASTR = IKA_ACT*KAASTR + IKA_PAS*KAP + IKA0 - 1
           KBSTR = IKB_ACT*KBASTR + IKB_PAS*KBP + IKB0 - 1
           IF(NTEST.GE.100) WRITE(6,*) ' KASTR, KBSTR', KASTR, KBSTR
           IF(NTEST.GE.100) WRITE(6,*) ' IKB_ACT, IKB_PAS, IKB0',
     &                                   IKB_ACT, IKB_PAS, IKB0
           IF(NTEST.GE.100) WRITE(6,*) ' KBASTR, KBP =', KBASTR,KBP
           IF(K_ACTIVE.EQ.1) 
     &     JASTR = IJA_ACT*JAASTR + IJA_PAS*KAP + IJA0 - 1
           IF(K_ACTIVE.EQ.1) 
     &     JBSTR = IJB_ACT*JBASTR + IJB_PAS*KBP + IJB0 - 1
           IF(NTEST.GE.100.AND.K_ACTIVE.EQ.1) 
     &     WRITE(6,*) ' JASTR, JBSTR', JASTR,JBSTR
           IF(K_ACTIVE.EQ.1.AND.(JASTR.LE.0.OR.JBSTR.LE.0)) THEN
             WRITE(6,*) ' Problem : JASTR, JBSTR =', JASTR,JBSTR 
             WRITE(6,*)
             WRITE(6,*) ' IJA_ACT, IJA_PAS, IJA0 = ',
     &                    IJA_ACT, IJA_PAS, IJA0
             WRITE(6,*) ' JAASTR, KAP = ', JAASTR, KAP
             WRITE(6,*) ' IJB_ACT, IJB_PAS, IJB0 = ',
     &                    IJB_ACT, IJB_PAS, IJB0
             WRITE(6,*) ' JBASTR, KBP = ', JBASTR, KBP
             STOP '..SP2 : Negative JASTR or JBSTR '
           END IF
*. Address C(0,Jb,Ja)
           IA_PJBJA = (JASTR-1)*NJB*NP + (JBSTR-1)*NP
           IF(K_ACTIVE.EQ.1.AND.IA_PJBJA.LT.0) THEN
             WRITE(6,*) ' Problem IA_PJBJA < 0 ', IA_PJBJA
             STOP ' ..SP2  Problem IA_PJBJA < 0 '
           END IF
*. Address C(0,Kb,Ka,jl)
           IA_PKBKAJL = ((LORB-1)*NJ+JORB-1)*NKA*NKB*NP+
     &                   (KASTR-1)*NKB*NP+(KBSTR-1)*NP
C          WRITE(6,*) ' SP2 : LORB,JORB,NJ,NKA,NKB,NP',
C    &                        LORB,JORB,NJ,NKA,NKB,NP
           IF(K_ACTIVE.EQ.0) THEN
            IF(IWAY.EQ.1) THEN
             DO IP = 1,NP
              IA_PKBKAJL = IA_PKBKAJL + 1
              IF(NTEST.GE.100) 
     &        WRITE(6,*)' IA_PKBKAJL set to zero ',  IA_PKBKAJL
              CPKKJL( IA_PKBKAJL ) = 0.0D0
             END DO
            END IF
           ELSE
            IF(IWAY.EQ.1) THEN
             DO IP = 1, NP
              IA_PKBKAJL = IA_PKBKAJL + 1
              IA_PJBJA   = IA_PJBJA + 1
              IF(NTEST.GE.100) 
     &        WRITE(6,*)' IA_PKBKAJL, IA_PJBJA',  IA_PKBKAJL, IA_PJBJA
              CPKKJL( IA_PKBKAJL ) = SIGN*CPA(IA_PJBJA)
             END DO
            ELSE IF(IWAY.EQ.2) THEN
             DO IP = 1, NP
              IA_PKBKAJL = IA_PKBKAJL + 1
              IA_PJBJA   = IA_PJBJA + 1
              CPA(IA_PJBJA) =  CPA(IA_PJBJA) + SIGN*CPKKJL( IA_PKBKAJL ) 
              IF(NTEST.GE.100) 
     &        WRITE(6,*)' IA_PKBKAJL, IA_PJBJA',  IA_PKBKAJL, IA_PJBJA
             END DO
            END IF
*           End IWAY switch
           END IF
*          ^ End of K_ACTIVE switch 
          END DO
          END DO
*         ^ End of loop over passive strings
         END DO
*        ^ End of loop over L
        END DO
*       ^ End of loop over KBASTR
       END DO
*      ^ end of loop over J
      END DO
*     ^ end of loop over KAASTR
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from  CPA_TO_CPKKJL_SP2 '
        JL = 0
        DO L = 1, NL
         DO J = 1, NJ
          JL = JL + 1
          IA_PKBKAJL = (JL-1)*NKA*NKB*NP+1
          WRITE(6,*) ' C(P,Ka,Kb,jl) as C(P,KaKb) for j,l =',J,L
          CALL WRTMAT(CPKKJL(IA_PKBKAJL),NP,NKA*NKB,NP,NKA*NKB) 
*
         END DO
        END DO
        WRITE(6,*) ' C(P,Jb,Ja) as C(PJb,Ja) '
        CALL WRTMAT(CPA,NP*NJB,NJA,NP*NJB,NJA)
      END IF
*
C     WRITE(6,*) ' Memtest at End of SP2 '
C     CALL MEMCHK
C     WRITE(6,*) ' Memchk passed '
      RETURN
      END
