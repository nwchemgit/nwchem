      SUBROUTINE TRP_TT_BLK(C,NA,NB,ISM,NSMST,IFB,SCR)
*
* Transpose or back transpose TT block 
*
* Jeppe Olsen, August 99
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NA(*),NB(*)
*. Input and output
      DIMENSION C(*)
*. Scratch  should hold largest TTS block)
      DIMENSION SCR(*)
*. Symmetry info
      INCLUDE 'multd2h.inc'
*
      IOFF = 1
      DO IASM = 1, NSMST
        IBSM = MULTD2H(IASM,ISM)
        NIA = NA(IASM)
        NIB = NB(IBSM)
        IF(IFB.EQ.1) THEN
          CALL TRPMT3(C(IOFF),NIA,NIB,SCR)
        ELSE
          CALL TRPMT3(C(IOFF),NIB,NIA,SCR)
        END IF
        CALL COPVEC(SCR,C(IOFF),NIA*NIB)
        IOFF = IOFF + NIA*NIB
      END DO
*
      RETURN
      END
      SUBROUTINE ADVICE_SIGMA2(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,IIADVICE,
     &           IJKLTP,NTERM)
*
* Advice Sigma routine about best route to take
*
* ITERM : Term  to be studied :  
*         =1 alpha-beta term 
*         ....... ( to be continued )
*
* LADVICE : ADVICE given ( short, an integer !!)
*
* For ITERM = 1 : 
*           LADVICE = 1 : Business as usual, no transpose of matrix
*                         (resolution on alpha strings, direct exc on beta)
*           LADVICE = 2 = Transpose matrices
*                         (resolution on beta strings, direct exc on alpha)
*
* Version returning all types of excitations, as well as advice for each
* excitation
*
*      Jeppe Olsen, August 99
*              
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Specific input
      INTEGER IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
*. Output
      INTEGER  IJKLTP(4,256),IIADVICE(256)
*. Local scratch
      INTEGER ILTP(16),JLTP(16),KLTP(16),LLTP(16)
*
      NTEST = 00
      NTERM = 0
      IF(ITERM.EQ.1) THEN
*.
*. sigma(i,Ka,Ib) = sum(i,kl)<Ib!Eb_kl!Jb>(ij!kl)C(j,Ka,Jb)
*
* Number of ops : Number of sx(kl) N_i*N_j_dimension of C(j,Ka,Jb)
*.No absolute calc of flops is made, only a relative measure
*
* Single excitations connecting the two types
*
C            SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ILTP,JLTP,NGAS,IAOCC,JAOCC,IPHGAS)
      CALL SXTYP2_GAS(NKLTYP,KLTP,LLTP,NGAS,IBOCC,JBOCC,IPHGAS)
*
      DO LIJTP = 1, NIJTYP
      DO LKLTP = 1, NKLTYP
*
        ITP = ILTP(LIJTP)
        JTP = JLTP(LIJTP)
        KTP = KLTP(LKLTP)
        LTP = LLTP(LKLTP)
*
        NTERM = NTERM + 1
        IJKLTP(1,NTERM) = ITP
        IJKLTP(2,NTERM) = JTP
        IJKLTP(3,NTERM) = KTP
        IJKLTP(4,NTERM) = LTP
*. Resolution : Particle or hole resolution
        IF(IPHGAS(ITP).EQ.2.AND.IPHGAS(JTP).EQ.2) THEN
          IPHIJ = 2
        ELSE
          IPHIJ = 1
        END IF
*
        IF(IPHGAS(KTP).EQ.2.AND.IPHGAS(LTP).EQ.2) THEN
          IPHKL = 2
        ELSE
          IPHKL = 1
        END IF
*
        IF(IADVICE.EQ.0) THEN
*. ph modifications or no advice asked
          LADVICE = 1
        ELSE
* =========================================
*.. Index for flops along C(j,Ka,Jb) route
* =========================================
*.Dim of C(j,Ka,Jb) relative to C(Ja,Jb)
          XNJOB = FLOAT(NOBPT(JTP))
          XNJEL = FLOAT(JAOCC(JTP))
          IF(IPHIJ.EQ.1) THEN
*. going from Ja to Ka reduces occ by one elec, changes dim by n/(N-n+1)
            XCJKAJB = XNJOB*XNJEL/(XNJOB-XNJEL+1)
          ELSE
*. going from Ja to Ka increases occ by one elec, changes dim by (N-n)/(n+1)
            XCJKAJB = XNJOB*(XNJOB-XNJEL)/(XNJEL+1)
          END IF
*. Number of kl excitations per beta string : 
          XNKLSX = FLOAT((NOBPT(KTP)-JBOCC(KTP))*JBOCC(LTP))
*. Number of ops (relative to dim of C)
          XNIOB = FLOAT(NOBPT(ITP))
          XFLOPA = XCJKAJB*XNKLSX*XNIOB
* =========================================
*.. Index for flops along C(l,Ja,Kb) route
* =========================================
*.Dim of C(l,Ja,Kb) relative to C(Ja,Jb)
          XNLOB = FLOAT(NOBPT(LTP))
          XNLEL = FLOAT(JBOCC(LTP))
          IF(IPHKL.EQ.1) THEN
            XCLJAKB = XNLOB*XNLEL/(XNLOB-XNLEL+1)
          ELSE
            XCLJAKB = XNLOB*(XNLOB-XNLEL)/(XNLEL+1)
          END IF
*. Number of ij excitations per alpha string : 
          XNIJSX = FLOAT((NOBPT(ITP)-JAOCC(ITP))*JAOCC(JTP))
*. Number of ops (relative to dim of C)
          XNKOB = FLOAT(NOBPT(KTP))
          XFLOPB = XCLJAKB*XNIJSX*XNKOB
*. Switch to second route if atleast 20 percent less work
          IF(XFLOPB.LE.0.8*XFLOPA) THEN
            LADVICE = 2
          ELSE
            LADVICE = 1
          END IF
*. Well, an additional consideration :
* If the C block involes the smallest allowed number of elecs in hole space,
* and the annihilation is in hole space
* then we do the annihilation in the space with the smallest number of 
* hole electrons.
          LHOLEA =0
          LHOLEB =0
          DO IGAS = 1, NGAS
            IF(IPHGAS(IGAS).EQ.2) THEN
              LHOLEA = LHOLEA + JAOCC(IGAS)
              LHOLEB = LHOLEB + JBOCC(IGAS)
            END IF
          END DO
*
          IF(IPHIJ.EQ.1.AND.IPHKL.EQ.1.AND.LHOLEA+LHOLEB.EQ.MNHL.AND.
     &      (IPHGAS(JTP).EQ.2.OR.IPHGAS(LTP).EQ.2))  THEN
*
            IF(IPHGAS(JTP).EQ.2) THEN
             KHOLEA = LHOLEA-1
             KHOLEB = LHOLEB 
            ELSE 
             KHOLEA = LHOLEA
             KHOLEB = LHOLEB - 1
            END IF
*
            IF(KHOLEA.EQ.KHOLEB) THEN
              LLADVICE = LADVICE
            ELSE IF(KHOLEA.LT.KHOLEB) THEN
              LLADVICE= 1
            ELSE
              LLADVICE = 2
            END IF
            IF(NTEST.GE.100.AND.LADVICE.NE.LLADVICE) THEN
              WRITE(6,*) ' Advice changed by hole considetions'
              WRITE(6,*) ' LADVICE, LLADVICE', LADVICE,LLADVICE
            END IF
            LADVICE = LLADVICE  
          END IF
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' ADVICE active '
            WRITE(6,*) ' IAOCC JAOCC IBOCC JBOCC'
            CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
            WRITE(6,*) ' ITP JTP KTP LTP ',ITP,JTP,KTP,LTP
            WRITE(6,*) ' XFLOPA,XFLOPB', XFLOPA,XFLOPB
            WRITE(6,*) ' ADVICE given : ', LADVICE
          END IF
         END IF 
*        ^ End if advice is sought
         IIADVICE(NTERM) = LADVICE
        END DO
        END DO
*       ^ End of loops over ij,kl excitations
      END IF
*     ^ End if ITERM test ( type of excitation)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ADVICE routine '
        WRITE(6,*)
*
        WRITE(6,*) ' IAOCC JAOCC IBOCC JBOCC'
        CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
        
        WRITE(6,*) ' Number of connecing excitations : ', NTERM
        WRITE(6,*)
        WRITE(6,*) ' ITP  JTP  KTP  LTP  ADVICE '
        WRITE(6,*) ' ==========================='
        DO KTERM = 1, NTERM
          WRITE(6,'(5I5)') (IJKLTP(II,KTERM),II=1,4),IIADVICE(KTERM)
        END DO
      END IF
*

      RETURN
      END
c $Id$
