      SUBROUTINE RSBB2B_M(IATP,IBTP,NIA,NIB,JATP,JBTP,NJA,NJB,
     &               IAGRP,IBGRP,ISBLOCK,ICBLOCK,
     &               NGAS,IAOCC,IBOCC,JAOCC,JBOCC,
     &               SB,CB,NOBPTS,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &               XINT,NSMOB,NSMST,IUSEAB,
     &               CJRES,SIRES,SCLFACX,NTEST,IUSE_PH,IPHGAS,
     &               ICSM,ISSM,C2,IREO_IA,IREO_IB,IREO_JA,IREO_JB,
     &               IUSE_HW)
*
* Master routine for alpha-beta loop for path using complete TT block
*
* ISBLOCK, ICBLOCK, SB,CB shold be positioned at start of TT block
*                                 
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*.input
      INTEGER NIA(*), NIB(*), NJA(*), NJB(*)
      INTEGER NOBPTS(MXPNGAS,*)
      INTEGER ISBLOCK(8,*),ICBLOCK(8,*)
*.Input   
      DIMENSION CB(*)
*.Output 
      DIMENSION SB(*)
*. Scratch 
      DIMENSION C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),XI1S(*),XI2S(*),XI3S(*)
      DIMENSION CJRES(*),SIRES(*)
*. Local scratch 
      DIMENSION IADVICE(256),IJKLTP(4,256)
*
C     DIMENSION SCLFAC(*)
*
*. Should be reinstated 
C     SCLFACX = SCLFAC(ICBLK)
      SCLFACX = 1.0D0
*. Call advice routine
      CALL ADVICE_SIGMA2(IAOCC,IBOCC,JAOCC,JBOCC,1,
     &                   IADVICE,IJKLTP,NTERM)
C     ADVICE_SIGMA2(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,IADVICE,
C    &           IJKLTP,NTERM)
*
*. First all terms without transposing
*
      DO ITERM = 1, NTERM
        IF(IADVICE(ITERM).EQ.1) THEN
*. Symmetryblockcks in TT block are ordered according to alpha sym, ie. 
*. second index 
          I12ORD = 2
          CALL RSBB2BEN2(IATP,IBTP,NIA,NIB,JATP,JBTP,NJA,NJB,                           
     &         IAGRP,IBGRP,NGAS,IAOCC,IBOCC,JAOCC,JBOCC,
     &         IJKLTP(1,ITERM),IJKLTP(2,ITERM),IJKLTP(3,ITERM),
     &         IJKLTP(4,ITERM),
     &         SB,CB,NOBPTS,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &         XINT,NSMOB,NSMST,IUSEAB,
     &         CJRES,SIRES,SCLFACX,NTEST,IUSE_PH,IPHGAS,
     &         ICSM,ISSM,C2,IREO_IA,IREO_IB,IREO_JA,IREO_JB,
     &         IUSE_HW,I12ORD)
        END IF
      END DO
*. Are there terms requiring additional transposing
      NTRNSP = 0
      DO ITERM = 1, NTERM
        IF(IADVICE(ITERM).EQ.2) THEN
          NTRNSP = NTRNSP + 1
        END IF
      END DO
*
      IF(NTRNSP.GT.0) THEN
*. Use transposed matrices, i.e transpose the transposed matrices to 
*. get untransposed matrices 
*
C            TRP_TT_BLK(C,NA,NB,ISM,NSMST,IFB,SCR)
        CALL TRP_TT_BLK(SB,NIA,NIB,ISSM,NSMST,2,C2)
        CALL TRP_TT_BLK(CB,NJA,NJB,ICSM,NSMST,2,C2)
*. Symmetryblockcs in TT block are ordered according to alpha sym, ie. 
*. first index 
        I12ORD = 1
        DO ITERM = 1, NTERM
          IF(IADVICE(ITERM).EQ.2) THEN
            CALL RSBB2BEN2(IBTP,IATP,NIB,NIA,JBTP,JATP,NJB,NJA,   
     &           IBGRP,IAGRP,NGAS,IBOCC,IAOCC,JBOCC,JAOCC,  
     &           IJKLTP(3,ITERM),IJKLTP(4,ITERM),IJKLTP(1,ITERM),
     &           IJKLTP(2,ITERM),
     &           SB,CB,NOBPTS,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           XINT,NSMOB,NSMST,IUSEAB,
     &           CJRES,SIRES,SCLFACX,NTEST,IUSE_PH,IPHGAS,
     &           ICSM,ISSM,C2,IREO_IB,IREO_IA,IREO_JB,IREO_JA,
     &           IUSE_HW,I12ORD)
          END IF
        END DO
*
        CALL TRP_TT_BLK(SB,NIA,NIB,ISSM,NSMST,1,C2)
        CALL TRP_TT_BLK(CB,NJA,NJB,ICSM,NSMST,1,C2)
      END IF
*     ^ End if there are terms requiring additional transposing 
      RETURN
      END
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
      SUBROUTINE RSBB2BEN2(
     &           IATP,IBTP,NIA,NIB,JATP,JBTP,NJA,NJB,
     &           IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &           ITYP,JTYP,KTYP,LTYP,
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
* Version where loop over active ITP,JTP,KTP,LTP has been moved outside
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
* April 99   division of alpha and beta strings, and simultaneous 
*            treatment of several symmetry blocks with the same type
* August 99 : ITYP,JTYP,KTYP,LTYP moved outside to improve efficiency 
*             (allows transpose for selected excitations )

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
      COMMON/CMXCJ/MXCJ,MAXK1_MX,LSCMAX_MX
     
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
        WRITE(6,'(A,4I4)') 'ITYP,JTYP,KTYP,LTYP', ITYP,JTYP,KTYP,LTYP
        WRITE(6,*) ' TT-block of C   '
        CALL WRTVH1(CB,ICSM,NJB,NJA,NSMST,0)
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*
      IOPSM = MULTD2H(ICSM,ISSM)
*. Total length of C and sigma blocks
      LEN_C = LEN_TT_BLOCK(ICSM,NJA,NJB,NSMST)
      LEN_S = LEN_TT_BLOCK(ISSM,NIA,NIB,NSMST)
*
C     XCNORM = INPROD(CB,CB,LEN_C)
C     IF(XCNORM.EQ.0.0D0) GOTO 9999
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
C     CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
C     CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)           
C     IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
C     DO 2001 IJTYP = 1, NIJTYP
C       ITYP = ITP(IJTYP)
C       JTYP = JTP(IJTYP)
        IF(NTEST.GE.1000)
     &  WRITE(6,*) ' ITYP, JTYP =', ITYP,JTYP
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
C       DO 2000 KLTYP = 1, NKLTYP
C          KTYP = KTP(KLTYP)
C          LTYP = LTP(KLTYP)
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
             IF(ITYP.NE.JTYP) CALL QENTER('PHPH ')
           ELSE IF(NHSPC.EQ.3) THEN
             CALL QENTER('PHHH ')
           ELSE IF(NHSPC.EQ.4) THEN
             CALL QENTER('HHHH ')
           END IF
*

           IF(NTEST.GE.1000)
     &     WRITE(6,*) ' KTYP, LTYP =', KTYP,LTYP
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
C?         WRITE(6,*) ' NSMST = ', NSMST
           NOREO_JA = 1
           NOREO_JB = 1
           NOREO_IA = 1
           NOREO_IB = 1
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
     &              NOREO_JAL)
               IF(NOREO_JAL.EQ.0) NOREO_JA = 0
             END IF

*. Split JB strings into active/passive part
             KASM = MULTD2H(ICSM,ISTSM)
             IF(NJA(KASM).GT.0) THEN
               CALL REO_STR_SPGRP3(JBSPGP,NGAS,ISTSM,NJB(ISTSM),2,
     &              KL_TYP,NACJB,IACJB,NJBAC_S,NJBPA_S,
     &              IBREO_JB(1,ISTSM),IREO_JB(IJBOFF),SIGNJB,NACACJB,
     &              NOREO_JBL)
               IF(NOREO_JBL.EQ.0) NOREO_JB = 0
             END IF
*. Split Ia strings into active/passive part
             KBSM = MULTD2H(ISSM,ISTSM)
             IF(NIB(KBSM).GT.0) THEN
               CALL REO_STR_SPGRP3(IASPGP,NGAS,ISTSM,NIA(ISTSM),2,
     &              IJ_TYP,NACIA,IACIA,NIAAC_S,NIAPA_S,
     &              IBREO_IA(1,ISTSM),IREO_IA(IIAOFF),SIGNIA,NACACIA,
     &              NOREO_IAL)
               IF(NOREO_IAL.EQ.0) NOREO_IA = 0
             END IF
*. Split Ib strings into active/passive part
             KASM = MULTD2H(ISSM,ISTSM)
             IF(NIA(KASM).GT.0) THEN
             CALL REO_STR_SPGRP3(IBSPGP,NGAS,ISTSM,NIB(ISTSM),2,
     &            KL_TYP,NACIB,IACIB,NIBAC_S,NIBPA_S,
     &            IBREO_IB(1,ISTSM),IREO_IB(IIBOFF),SIGNIB,NACACIB,
     &              NOREO_IBL)
               IF(NOREO_IBL.EQ.0) NOREO_IB = 0
             END IF
           END DO
           NOREO_C = NOREO_JA*NOREO_JB
           NOREO_S = NOREO_IA*NOREO_IB
*. NOREO_C business is not working so 
           NOREO_C = 0
           NOREO_S = 0
*
           CALL ACOP_SPGRP(NACJA,IACJA,JAC,JTYP,KACGRP,IJCODE)
C          IF(IJCODE.EQ.0) GOTO 2001
           IF(IJCODE.EQ.0) GOTO 1999
*          ^ Temp. changed for qenter/qexit timings
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
     &                 NOREO_C)
           IF(NOREO_C.EQ.0) CALL COPVEC(C2,CB,LEN_C)
*. Reorganize S(Ia,IB) as C(Ia_pa,Ib_pa,Ia_ac,Ia_pa)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' S(A,B) => S(P,A) '
           END IF
           CALL CJPA_MS(SB,C2,1,NSMST,NIA,NIB,NIAAC_S,
     &                 NIBAC_S,NIAPA_S,NIBPA_S,IREO_IA,
     &                 IREO_IB,IBREO_IA,IBREO_IB,ISSM,SIGNIAB,
     &                 NIPAAB,NIACAB,IBSPA,IBSPA2,LEN_S,I12ORD,NOREO_S)
           IF(NOREO_S.EQ.0) CALL COPVEC(C2,SB,LEN_S)
C?           WRITE(6,*) ' NACACIB 3 ', NACACIB
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
               LSCMAX_MX = MAX(LSCMAX,LSCMAX_MX)
               MAXK1_MX = MAX(MAXK1,MAXK1_MX)
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
               CALL SETVEC(SIRES,ZERO,LEN_SIRES)              
               FACS = 1.0D0
*
               DO 1931 KLSM = 1, NSMOB
               DO 1930 KSM = 1, NSMOB
*
                 IFIRST = 1
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
     &                    ONE,ONE)
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
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' S(P,A) => S(A,B) '
           END IF
            CALL CJPA_MS(C2,SB,2,NSMST,NIA,NIB,NIAAC_S,NIBAC_S,
     &           NIAPA_S,NIBPA_S,IREO_IA,IREO_IB,IBREO_IA,IBREO_IB,
     &           ISSM,SIGNIAB,NIPAAB,NIACAB,IBSPA,IBSPA2,LEN_S,I12ORD,
     &           NOREO_S) 
            IF(NOREO_S.EQ.0) CALL COPVEC(C2,SB,LEN_S)
*. Reform also C from passive/active form back to Ja,Jb
            ONE = 1.0D0
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' C(P,A) => C(A,B) '
           END IF
            CALL CJPA_MS(C2,CB,2,NSMST,NJA,NJB,NJAAC_S,NJBAC_S,
     &           NJAPA_S,NJBPA_S,IREO_JA,IREO_JB,IBREO_JA,IBREO_JB,
     &           ICSM,SIGNJAB,NJPAAB,NJACAB,IBCPA,IBCPA2,LEN_C,I12ORD,
     &           NOREO_C) 
            IF(NOREO_C.EQ.0) CALL COPVEC(C2,CB,LEN_C)
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
             IF(ITYP.NE.JTYP) CALL QEXIT('PHPH ')
           ELSE IF(NHSPC.EQ.3) THEN
             CALL QEXIT('PHHH ')
           ELSE IF(NHSPC.EQ.4) THEN
             CALL QEXIT('HHHH ')
           END IF
C2000   CONTINUE
*           ^ End of loop over KLTYP
C2001 CONTINUE
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
         
