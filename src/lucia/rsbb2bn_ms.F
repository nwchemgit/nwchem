      SUBROUTINE RSBB2BN_MS(ICSM,ISSM,IATP,IBTP,NIA_MS,NIB_MS,
     &                   JATP,JBTP,NJA_MS,NJB_MS,
     &                   IAGRP,IBGRP,IB_IAGRP,IB_IBGRP,NGAS,IAOC,IBOC,
     &                   JAOC,JBOC,
     &                   SB,CB,ADSXA,STSTSX,MXPNGASX,
     &                   NOBPTS,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                   XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBSX,IUSEAB,
     &                   CJRES,SIRES,SCLFAC,NTESTG,
     &                   NSEL2E,ISEL2E,IUSE_PH,IPHGAS,XINT2,
     &                   LEN_MAP)
*.
*. On input  the symmetry-blocks of CB are assumed transposed
*. On output the symmetry-blocks of SB are assumed transposed
*. LEN_MAP is length of map arrays I1, XI1S ... for given symmetry
*
* Combined alpha-beta double excitation, ALL SYMMETRIES
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
* =====
* Input
* =====
*
* ICSM,ISSM : Symmetry of C and S, respectively,
* IASM,IATP : Symmetry and type(rel) of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type(rel) of beta   strings in sigma
* JASM,JATP : Symmetry and type(rel) of alpha  strings in C
* JBSM,JBTP : Symmetry and type(rel) of beta   strings in C
* NIA_MS,NIB_MS : Number of alpha-(beta-) strings in sigma, all symmetries
* NJA_MS,NJB_MS : Number of alpha-(beta-) strings in C, all symmetries
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
*
* Last change : Aug 2000
* June 2003 : Adding possibility for treating all symmetries at the 
* same time 
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'multd2h.inc'
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
      REAL*8 INPROD
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
      DIMENSION NIA_MS(NSMST),NIB_MS(NSMST),NJA_MS(NSMST),NJB_MS(NSMST)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(LEN_MAP,NSMST),XI1S(LEN_MAP,NSMST)
      DIMENSION I2(LEN_MAP,NSMST),XI2S(LEN_MAP,NSMST)
      DIMENSION I3(LEN_MAP,NSMST),XI3S(LEN_MAP,NSMST)
      DIMENSION I4(LEN_MAP,NSMST),XI4S(LEN_MAP,NSMST)
      DIMENSION XINT(*), XINT2(*)
      DIMENSION CJRES(*),SIRES(*)
*
C     DIMENSION H(MXPTSOB*MXPTSOB)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
C     DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_REO(2),IJ_AC(2)
      DIMENSION KL_TYP(2),KL_REO(2),KL_AC(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
      INTEGER KASPGP(20), NKA_MS(8)
      INTEGER IB_C(8), IB_S(8)
*. Arrays for reorganization 
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C    &          SIGNREO,NADOP,NADDEL,IADDEL,ADSIGN)
* 
C     INTEGER JASM_MS(8),JBSM_MS(8) 
      INTEGER NOP1_MS(8),NOP2_MS(8)
      INTEGER NOP3_MS(8),NOP4_MS(8)
      INTEGER IB_CKAJJB(8), IB_SKAIIB(8)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'oper.inc'
      CALL LUCIAQENTER('RS2Bs')
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.200) THEN
*
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' RSBB2BN_MS speaking '
        WRITE(6,*) ' =============== '
*
        WRITE(6,*) ' Occupation of IA '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of IB '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JA '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JB '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        
*
        WRITE(6,*) ' Initial TT block of C (transposed)'
C            WRT_TT_BLK(C,NA,NB,ITRP,ISM,NSMST)
        CALL WRT_TT_BLK(CB,NJA_MS,NJB_MS,1,ICSM,NSMST)

        WRITE(6,*) ' Memcheck at start of RSBB2BN '
        CALL MEMCHK
        WRITE(6,*) ' Memory check passed '
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Offsets for symmetryblocks in C and S - remember that the blocks 
* of C and S are transposed 
* have been transposed 
C     GET_OFFSET_FOR_SYMBLK(NROW,NCOL,IRC,ISYM,NSMST,IB)
      CALL GET_OFFSET_FOR_SYMBLK(NJB_MS,NJA_MS,2,ICSM,NSMST,IB_C)
      CALL GET_OFFSET_FOR_SYMBLK(NIB_MS,NIA_MS,2,ISSM,NSMST,IB_S)
*. occupation of GAS spaces for each supergroup ( stored in **SPGP)
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
*. Symmetry of operator connecting C and S
      IHSM = MULTD2H(ICSM,ISSM)
*.Types of SX that connects the two strings
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)           
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
*
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
*. Should N-1 or N+1 projection be used for alpha strings
        IJ_TYP(1) = ITYP
        IJ_TYP(2) = JTYP
        IJ_AC(1)  = 2
        IJ_AC(2) =  1
        NOP = 2
        IF(IUSE_PH.EQ.1) THEN
          CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,
     &         SIGNIJ)
        ELSE
*. Enforced a+ a
          IJ_REO(1) = 1
          IJ_REO(2) = 2
          SIGNIJ = 1.0D0
        END IF
*. Two choices here :
*  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> ( good old creation mapping)
*  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka>  + delta(i,j)                   
        IF(IJ_REO(1).EQ.1.AND.IJ_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
          IJAC = 2
          SIGNIJ2 = SCLFAC
*
          IJ_TYP(1) = ITYP
          IJ_TYP(2) = JTYP
*
          IOP1TP = ITYP
          IOP2TP = JTYP
        ELSE
*. Terra Nova, annihilation map 
          IJAC = 1
          SIGNIJ2 = -SCLFAC
*
          IJ_TYP(1) = JTYP
          IJ_TYP(2) = ITYP
*
          IOP1TP = JTYP
          IOP2TP = ITYP
        END IF
        IF(IJAC.EQ.1) THEN
           IOP1AC = 1
           IOP2AC = 2 
        ELSE
           IOP1AC = 2
           IOP2AC = 1
        END IF
*
        DO IOPSM = 1, NSMST
          NOP1_MS(IOPSM) = NOBPTS(IOP1TP,IOPSM)
          NOP2_MS(IOPSM) = NOBPTS(IOP2TP,IOPSM)
        END DO
*. Ka strings are resolution strings <Ia! Iop1 !Ka><Ka! Iop2 !Ja> 
*. Obtain type of Ka strings as occupation in each GAS space 
C  NEWTYP(IATPABS,IAORC,IOBTP,1,KATP)
*. 
C       JASPGP_ABS = JASPGP + IB_IAGRP - 1
C       CALL NEWTYP(JASPGP_ABS,IOP2AC,IOP2TP,1,KASPGP)
C       WRITE(6,*) ' NSMST = ', NSMST
*. Number of Ka strings per symmetry
C       CALL NST_SPGP(KASPGP,NKA_MS)
*. Loop over symmetries of Ka-strings 
        DO 1940 KASM = 1, NSMST
C         NKA = NKA_MS(KASM)
C         NKAEFF = NKA
*. Set up all mappings from Ka-strings of this symmetry to Ja and Ia. 
* I1, XI1S : <Ja! Op2 + ! Ka>, stored for all symmetries of Op2
* I3, XI3S : <Ia! Op1   ! Ka>, stored for all symmetries of Op1
          DO IOPSM = 1, NSMST
            IJASM = MULTD2H(KASM,IOPSM)
*. The folowing are not active ..
            IFRST = 1
            KFRST = 1
            CALL ADAST_GAS(IOPSM,IJ_TYP(2),NGAS,JASPGP,IJASM,
     &           I1(1,IOPSM),XI1S(1,IOPSM),NKASTR,IEND,IFRST,KFRST,KACT,
     &           SIGNIJ2,IJAC)
*
            CALL ADAST_GAS(IOPSM,IJ_TYP(1),NGAS,IASPGP,IJASM,
     &           I3(1,IOPSM),XI3S(1,IOPSM),NKASTR,IEND,IFRST,KFRST,KACT,
     &           ONE,IJAC)
          END DO
          NKA = NKASTR
          NKAEFF = NKASTR
C?        WRITE(6,*) ' KASM, NKA = ', KASM, NKA
*. Compress list to common nonvanishing elements, not active (?)
C         IDOCOMP = 0
C         IF(IDOCOMP.EQ.1) THEN
C             CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),
C    &                        NKASTR,NKAEFF)
C         ELSE 
C             NKAEFF = NKASTR
C         END IF
*. Loop over batches of KA strings
          NKABTC = NKA/MAXK   
          IF(NKABTC*MAXK.LT.NKA) NKABTC = NKABTC + 1
*
          DO 1801 IKABTC = 1, NKABTC
C?          WRITE(6,*) ' 1801, IKABTC = ', IKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKA)
            LKABTC = KATOP-KABOT+1
*
*. Obtain C(ka,J,JB) for Ka in batch
*
C     GET_CKAJJB_MS(CB,ICSM,NJA,NJB,IB_C,KABOT,LKABTC,KASM,
C    &                         NJORB,NSMST,IKJMAP,XIKJMAP,LKJMAP,
C    &                         CJRES,IB_CJRES)
*. Determine offsets for symmetry blocks of C(Ka,j,Jb)
            CALL IOFF_CKAJJB_MS(ICSM,NJA_MS,NJB_MS,LKABTC,KASM,
     &           NOP2_MS,NSMST,IB_CKAJJB,LEN_CKAJJB)
*
            CALL GET_CKAJJB_MS(
     &      CB,ICSM,NJA_MS,NJB_MS,IB_C,KABOT,LKABTC,KASM,NKASTR,
     &      NOP2_MS,NSMST,I1,XI1S,LEN_MAP,CJRES,IB_CKAJJB) 
*
*. Offsets and length of matrix S(Ka,I,IB)
C     IOFF_CKAJJB_MS(ICSM,NJA_MS,NJB_MS,LKABTC,KASM,
C    &                         NJORB,NSMST,IB_CKAJJB,LEN_CKAJJB)
       
            CALL IOFF_CKAJJB_MS(ISSM,NIA_MS,NIB_MS,LKABTC,KASM,
     &                          NOP1_MS,NSMST,IB_SKAIIB,LEN_SKAIIB)
            CALL SETVEC(SIRES,ZERO,LEN_SKAIIB)
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*. Allowed double excitation ?
              IJKL_ACT = I_DX_ACT(ITYP,KTYP,LTYP,JTYP)
              IF(IJKL_ACT.EQ.0) GOTO 2000
              IF(NTEST.GE.100) THEN
                WRITE(6,*) ' KTYP, LTYP', KTYP, LTYP 
              END IF
*. Should this group of excitations be included 
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              KL_TYP(1) = KTYP
              KL_TYP(2) = LTYP
              KL_AC(1)  = 2
              KL_AC(2) =  1
              NOP = 2
              IF(IUSE_PH.EQ.1) THEN
                CALL ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,
     &               SIGNKL)
              ELSE
*. Enforced a+ a
                KL_REO(1) = 1
                KL_REO(2) = 2
                SIGNKL = 1.0D0
              END IF
*
              IFIRST = 1
*
              IF(KL_REO(1).EQ.1.AND.KL_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
                KLAC = 2
                KL_TYP(1) = KTYP
                KL_TYP(2) = LTYP
                IOP3TP = KTYP
                IOP4TP = LTYP
                ELSE
*. Terra Nova, annihilation map 
                KLAC = 1
                KL_TYP(1) = LTYP
                KL_TYP(2) = KTYP
                IOP3TP = LTYP
                IOP4TP = KTYP
              END IF
              DO IOPSM = 1, NSMST
                NOP3_MS(IOPSM) = NOBPTS(IOPSM,IOP3TP)
                NOP4_MS(IOPSM) = NOBPTS(IOPSM,IOP4TP)
              END DO
*. Loop over symmetry of formal resolution string Kb. These
*. resolution strings are not used for actual resolution. but 
*. only for obtaining adressing of single excitations.
              DO KBSM = 1, NSMST
* set up connections between !Ia> and !Ja> as 
* <Ia!Iop2!Ka><Ka!Iop1!Ka> for all symmetries of IOP2 and IOP1
                DO IOPSM = 1, NSMST
                 IJBSM = MULTD2H(IOPSM,KBSM)
*
                 CALL ADAST_GAS(IOPSM,KL_TYP(2),NGAS,JBSPGP,IJBSM,
     &           I2(1,IOPSM),XI2S(1,IOPSM),NKBSTR,IEND,IFRST,KFRST,KACT,
     &           SIGNKL,KLAC)
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                 CALL ADAST_GAS(IOPSM,KL_TYP(1),NGAS,IBSPGP,IJBSM,
     &           I4(1,IOPSM),XI4S(1,IOPSM),NKBSTR,IEND,IFRST,KFRST,KACT,
     &           ONE,KLAC)
                END DO
C               WRITE(6,*) ' KBSM, NKBSTR = ', KBSM, NKBSTR
*. It is now time to loop over SYM of I,J,K,L, and get some work 
*. done ...
                DO ISM = 1, NSMST 
                DO JSM = 1, NSMST
                DO KSM = 1, NSMST
                  IJSM = MULTD2H(ISM,JSM)
                  IJKSM = MULTD2H(IJSM,KSM)
                  LSM = MULTD2H(IJKSM,IHSM)
C?                WRITE(6,*) ' ISM,JSM,KSM,LSM = ',
C?   &                         ISM,JSM,KSM,LSM
*
                  NI = NOBPTS(IOP1TP,ISM)
                  NJ = NOBPTS(IOP2TP,JSM)
                  NK = NOBPTS(IOP3TP,KSM)
                  NL = NOBPTS(IOP4TP,LSM)
*
                  JASM = MULTD2H(KASM,JSM)
                  IASM = MULTD2H(KASM,ISM)
                  JBSM = MULTD2H(KBSM,LSM)
                  IBSM = MULTD2H(KBSM,KSM)
                  IF(MULTD2H(IASM,IBSM).NE.ISSM.OR.
     &               MULTD2H(JASM,JBSM).NE.ICSM    ) GOTO 1920
C?                WRITE(6,*) ' IASM, IBSM, JASM, JBSM = ',
C?   &                         IASM, IBSM, JASM, JBSM
*
                  NJA = NJA_MS(JASM)
                  NJB = NJB_MS(JBSM)
                  NIA = NIA_MS(IASM)
                  NIB = NIB_MS(IBSM)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                  IKORD = 0  
                  IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                  IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &            GOTO 1930
                  IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &            IKORD = 1
*
                  IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*
* Fetch Integrals as (iop2 iop1 |  k l )
*
                  IXCHNG = 0
                  ICOUL = 1
                  ONE = 1.0D0
                  IF(I_USE_SIMTRH .EQ.0 ) THEN
*. Normal integrals with conjugation symmetry
                    CALL LUCIAGETINT(XINT,IOP2TP,JSM,IOP1TP,ISM,
     &                   IOP3TP,KSM,IOP4TP,LSM,IXCHNG,0,0,ICOUL,
     &                   ONE,ONE)
                    IF(NTEST.GE.500) THEN
                      WRITE(6,*) ' Integrals as (ji,kl) '
                      CALL WRTMAT(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
                    END IF
                  ELSE IF (I_USE_SIMTRH.EQ.1) THEN
*. Integrals does not have conjugation symmetry so be careful...
*. The following is not enough is particle hole symmetry is encountered
*. Obtain ( i j ! k l )
                    CALL LUCIAGETINT(XINT,IOP1TP,ISM,IOP2TP,JSM,
     &                               IOP3TP,KSM,IOP4TP,LSM,
     &                          IXCHNG,0,0,ICOUL,
     &                   ONE,ONE)
                    IF(KLAC.EQ.2.AND.IJAC.EQ.2) THEN
*. Transpose to obtain ( j i ! k l )
                      CALL TRP_H2_BLK(XINT,12,NI,NJ,NK,NL,XINT2)
                    ELSE IF(KLAC.EQ.1.AND.IJAC.EQ.2) THEN  
*. Transpose to obtain (j i | l k)
                      CALL TRP_H2_BLK(XINT,46,NI,NJ,NK,NL,XINT2)
                    ELSE IF (KLAC.EQ.1.AND. IJAC .EQ. 1 ) THEN
*. Transpose to obtai (i j | l k)
                      CALL TRP_H2_BLK(XINT,34,NI,NJ,NK,NL,XINT2)
                    END IF
                  END IF
*
* S(Ka,i,Ib) = sum(j,k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                  IJKL_DIM = NI*NJ*NK*NL
                  IF(INPROD(XINT,XINT,IJKL_DIM).NE.0.0D0) THEN
                    IROUTE = 3
                    CALL SKICKJ(SIRES(IB_SKAIIB(ISM)),
     &                   CJRES(IB_CKAJJB(JSM)),LKABTC,NIB,NJB,
     &                   NKBSTR,XINT,NI,NJ,NK,NL,
     &                   NKBSTR,I4(1,KSM),XI4S(1,KSM),
     &                   I2(1,LSM),XI2S(1,LSM),IKORD,
     &                     FACS,IROUTE )
                  END IF
*
                  IF(NTEST.GE.500) THEN
                    WRITE(6,*) ' Updated Sires as S(Kai,Ib) for ISM=',
     &              ISM
                    CALL WRTMAT(
     &              SIRES(IB_SKAIIB(ISM)),LKABTC*NI,NIB,LKABTC*NI,NIB)
                  END IF
*
 1930           CONTINUE
 1920           CONTINUE
                END DO
                END DO
                END DO
*               ^ End of loop over ISM, JSM, KSM
              END DO
*             ^ End of loop over KBSM
 2000       CONTINUE
*           ^ End of loop over KLTYP
*
*. Scatter out from s(Ka,Ib,i)
*
*
C     ADD_SKAIIB_MS(SB,ISSM,NIA_MS,NIB_MS,IB_S,KABOT,
C    &                         LKABTC,KASM,
C    &                         NIORB,NSMST,IKIMAP,XIKIMAP,LKIMAP,
C    &                         SIRES,IB_SKAIIB)
            CALL ADD_SKAIIB_MS(SB,ISSM,NIA_MS,NIB_MS,IB_S,KABOT,
     &           LKABTC,KASM,NKASTR,NOP1_MS,NSMST,I3,XI3S,LEN_MAP,
     &           SIRES,IB_SKAIIB)
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 1940   CONTINUE
*       ^ End of loop over KASM
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL LUCIAQEXIT('RS2Bs')
      RETURN
      END
      SUBROUTINE GET_CKAJJB_MS(CB,ICSM,NJA_MS,NJB_MS,IB_C,KABOT,
     &                         LKABTC,KASM,NKASTR,
     &                         NJORB,NSMST,IKJMAP,XIKJMAP,LKJMAP,
     &                         CJRES,IB_CJRES)
*
* Obtain C(Ka,j,JB) for all orbitals j of given type. The 
* matrix is ordered using symmetry of orbital j
*
* Jeppe Olsen, June 2003 in Aarhus 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NJA_MS(NSMST),NJB_MS(NSMST)
      INTEGER NJORB(NSMST)
      DIMENSION IKJMAP(LKJMAP,NSMST),XIKJMAP(LKJMAP,NSMST)
*. Offset to block in C with given sym of JA
      INTEGER IB_C(NSMST)
*. as usual C is transposed on input
      DIMENSION CB(*) 
*. Output
      DIMENSION CJRES(*), IB_CJRES(NSMST)
*
*. Set up IB_CJRES : IB_CJRES(ISM) gives first element of C(Ka,j,Jb)
*. where j has given sym 
*
      IOFF = 1
      DO ISM_J = 1, NSMST
*
       IB_CJRES(ISM_J) = IOFF
*
       ISM_JA = MULTD2H(ISM_J,KASM)
       ISM_JB = MULTD2H(ICSM,ISM_JA)
       IOFF = IOFF + NJORB(ISM_J)*NJB_MS(ISM_JB)*LKABTC
      END DO
*
*. Obtain C(ka,J,JB) for Ka in batch
*
      DO ISM_J = 1, NSMST
        ISM_JA = MULTD2H(ISM_J,KASM)
        ISM_JB = MULTD2H(ICSM,ISM_JA)
        IB = IB_C(ISM_JA)
        DO J = 1, NJORB(ISM_J)
          CALL GET_CKAJJB(CB(IB),NJORB(ISM_J),NJA_MS(ISM_JA),
     &    CJRES(IB_CJRES(ISM_J)),LKABTC,NJB_MS(ISM_JB),
     &    J,IKJMAP(KABOT+(J-1)*NKASTR,ISM_J),
     &    XIKJMAP(KABOT+(J-1)*NKASTR,ISM_J))
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' C(Ka,j,Jb) matrix as C(Kaj,Jb)'
        DO ISM_J = 1, NSMST
          WRITE(6,*) ' For symmetry of orbital = ', ISM_J
          ISM_JA = MULTD2H(ISM_J,KASM)
          ISM_JB = MULTD2H(ICSM,ISM_JA)
          NJ = NJORB(ISM_J)
          NJB = NJB_MS(ISM_JB)
          LEN_ROW = NJ*LKABTC
          CALL WRTMAT(CJRES(IB_CJRES(ISM_J)),
     &    LEN_ROW,NJB,LEN_ROW,NJB)
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE GET_OFFSET_FOR_SYMBLK(NROW,NCOL,IRC,ISYM,NSMST,IB)
*
* Obtain offsets for symmetryblocked matrix. 
*.If IRC = 1, matrix is ordered according to symmetry of row 
*.If IRC = 2, matrix is ordered according to symmetry of column
*
* Jeppe Olsen, June 2003
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NROW(NSMST),NCOL(NSMST)
*. Output
      INTEGER IB(NSMST)
*
      IOFF = 1
      IF(IRC.EQ.1) THEN
        DO ISM_ROW = 1, NSMST
          IB(ISM_ROW) = IOFF
          ISM_COL = MULTD2H(ISYM,ISM_ROW)
          IOFF = IOFF + NROW(ISM_ROW)*NCOL(ISM_COL)
        END DO
      ELSE
        DO ISM_COL = 1, NSMST
          IB(ISM_COL) = IOFF
          ISM_ROW = MULTD2H(ISYM,ISM_COL)
          IOFF = IOFF + NROW(ISM_ROW)*NCOL(ISM_COL)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE IOFF_CKAJJB_MS(ICSM,NJA_MS,NJB_MS,LKABTC,KASM,
     &                         NJORB,NSMST,IB_CKAJJB,LEN_CKAJJB)
*
* Obtain offsets for C(Ka,j,Jb) matrix containing all symmetries
*
* Jeppe Olsen, June 2003 on the train to Odense an early morning ...
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NJA_MS(NSMST),NJB_MS(NSMST)
      INTEGER NJORB(NSMST)
*. Output
      DIMENSION IB_CKAJJB(NSMST)
*
*. IB_CKAJJB(ISM) gives first element of C(Ka,j,Jb) where j has sym ISM
*
      IOFF = 1
      DO ISM_J = 1, NSMST
*
       IB_CKAJJB(ISM_J) = IOFF
*
       ISM_JA = MULTD2H(ISM_J,KASM)
       ISM_JB = MULTD2H(ICSM,ISM_JA)
       IOFF = IOFF + NJORB(ISM_J)*NJB_MS(ISM_JB)*LKABTC
      END DO
      LEN_CKAJJB = IOFF -1 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Offset array for C(Ka,j,Jb) matrix for given sym of j'
        CALL IWRTMA(IB_CKAJJB,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
      SUBROUTINE ADD_SKAIIB_MS(SB,ISSM,NIA_MS,NIB_MS,IB_S,KABOT,
     &                         LKABTC,KASM,NKASTR,
     &                         NIORB,NSMST,IKIMAP,XIKIMAP,LKIMAP,
     &                         SIRES,IB_SKAIIB)
*
*
* Add sigma contributions in S(Ka,i,IB) to sigma vector,
* for form where SIRES is calculated for all symmetries
*
*
* Jeppe Olsen, June 2003 on the way to Odense
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER NIA_MS(NSMST),NIB_MS(NSMST)
      INTEGER NIORB(NSMST)
      DIMENSION IKIMAP(LKIMAP,NSMST),XIKIMAP(LKIMAP,NSMST)
*. Offset to block in S with given sym of IA
      INTEGER IB_S(NSMST)
*. as usual S is transposed on input
      DIMENSION SB(*) 
*. Output
      DIMENSION SIRES(*), IB_SKAIIB(NSMST)
*
      DO ISM_I = 1, NSMST
        ISM_IA = MULTD2H(ISM_I,KASM)
        ISM_IB = MULTD2H(ISSM,ISM_IA)
        IB = IB_S(ISM_IA)
        DO I = 1, NIORB(ISM_I)
          CALL ADD_SKAIIB(SB(IB),NIORB(ISM_I),NIA_MS(ISM_IA),
     &    SIRES(IB_SKAIIB(ISM_I)),LKABTC,NIB_MS(ISM_IB),
     &    I,IKIMAP(KABOT+(I-1)*NKASTR,ISM_I),
     &    XIKIMAP(KABOT+(I-1)*NKASTR,ISM_I))
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' S(Ka,i,Ib) matrix as S(Kai,Ib)'
        DO ISM_I = 1, NSMST
          WRITE(6,*) ' For symmetry of orbital = ', ISM_I
          ISM_IA = MULTD2H(ISM_I,KASM)
          ISM_IB = MULTD2H(ISSM,ISM_IA)
          NI = NIORB(ISM_I)
          NIB = NIB_MS(ISM_IB)
          WRITE(6,*) ' NI, NIB = ', NI, NIB
          LEN_ROW = NI*LKABTC
          CALL WRTMAT(SIRES(IB_SKAIIB(ISM_I)),
     &    LEN_ROW,NIB,LEN_ROW,NIB)
        END DO
*
        WRITE(6,*) ' Updated block of sigma : '
        WRITE(6,*) ' ========================'
        DO IASM = 1, NSMST
          IBSM = MULTD2H(IASM,ISSM)
          NIA = NIA_MS(IASM)
          NIB = NIB_MS(IBSM)
          IB = IB_S(IASM)
          CALL WRTMAT(SB(IB),NIB,NIA,NIB,NIA)
        END DO
*
      END IF
*
      RETURN
      END 
      SUBROUTINE ENFORCE_SIMSYM_ACT_BLK(NBLK,IBLK,XSCA)
*
* NBLK blocks are defined by IBLK. The zero blocks 
* are flagged by a zero in XSCA. Modify XSCA, so all
* blocks with the same TT have the same activity
*
* Jeppe Olsen, July 1, 2003 -trying to get on coding mode again
*
      INCLUDE 'implicit.inc'
*
      INTEGER IBLK(8,NBLK)
      DIMENSION XSCA(NBLK)
*
      IB_TT = 1
      DO JBLK = 1, NBLK
C?     WRITE(6,*) ' IBLK(1,JBLK), IBLK(2,JBLK) = ',
C?   &              IBLK(1,JBLK), IBLK(2,JBLK)
*. Does this block have the same TT as the next 
       IEND_TT = 0
       IF(JBLK.EQ.NBLK) THEN
         IEND_TT = 1
       ELSE 
         IF(IBLK(1,JBLK).NE.IBLK(1,JBLK+1).OR.
     &      IBLK(2,JBLK).NE.IBLK(2,JBLK+1)    ) IEND_TT = 1
       END IF
*
       IF(IEND_TT.EQ.1) THEN
C?       WRITE(6,*) ' IEND_TT = 1 for JBLK = ', JBLK
C?       WRITE(6,*) ' IB_TT = ', IB_TT
*. Block JBLK is the end of a TT block, ensure that 
*. all blocks have the same activity index 
         I_AM_ACTIVE = 0
         DO KBLK = IB_TT,JBLK
           IF(XSCA(KBLK).EQ.1.0D0) I_AM_ACTIVE = 1
         END DO
C?       WRITE(6,*) ' I_AM_ACTIVE = ', I_AM_ACTIVE
         IF(I_AM_ACTIVE.EQ.1) THEN
           DO KBLK = IB_TT,JBLK
             XSCA(KBLK) = 1.0D0
           END DO
         END IF
         IB_TT = JBLK + 1
       END IF
*      ^ End of block was active
      END DO
*     ^ End of loop over blocks
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated activity index XSCA '
        CALL WRTMAT(XSCA,1,NBLK,1,NBLK)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_TT_BLK(C,NA,NB,ITRP,ISM,NSMST)
*
* Write all S blocks of a given TT block
*
*. Jeppe Olsen, July 2003
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INTEGER NA(NSMST),NB(NSMST)
      DIMENSION C(*)
      INCLUDE 'multd2h.inc'
*
      IOFF = 1
      DO ISM_A = 1, NSMST
        ISM_B = MULTD2H(ISM,ISM_A)
*
        LA = NA(ISM_A)
        LB = NB(ISM_B)
        WRITE(6,*) ' Symmetry of alpha strings = ', ISM_A
        IF(ITRP.NE.1) THEN
           CALL WRTMAT(C(IOFF),LA,LB,LA,LB)
        ELSE 
           CALL WRTMAT(C(IOFF),LB,LA,LB,LA)
        END IF
        IOFF = IOFF + LA*LB 
      END DO
*
      RETURN
      END
      
c $Id$
