C Some routines added for configurations, spin and CSF's 
* Storing of configurations and CSF's for CI
* The code is generating the configurations seperately for each 
* occupation class. I have, however, decided for the this time around,
* that the configurations are stored according to number of open 
* orbitals without regard for the occupation class.
*. The reordered configurations are thus ordered according to number
*. of open orbitals, 
* The lexical ordering goes as follows
* 1) Lexical addressing arrays are obtained for each occupation class
* 2) A base array, IB_OCCLS is obtained giving start of LEXICAL addresses
*    for given occupation class
* A configuration ICONF of class ICLS has thus the lexical adress
* IADR = IB_OCCLS(ICLS) -1 + IZNUM(ICONF,IOCLS) 
* The reorder array at adress IADT then givees the actual address
* (Dec. 2011)
*. Notice that the orbitals in the configurations are numbered without adding
* inactive orbitals- orbital 1 is first active orbital
      SUBROUTINE INFO_CONF_LIST(ISYM,LENGTH_LIST,NCONF_TOT)
*
* The number of configurations as a function of number of open orbitals
* is given in NCONF_PER_OPEN(*,ISYM). Determine offsets to various arrays
* and store in arrays in spinfo
*
*
* Jeppe Olsen, November 2001
*              Dec. 2011, some modifications
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
*. Output:  
*
*    INTEGER IB_CN_OPEN(MAXOP+1)
*    Offset in list of configurations for configurations 
*    with a given number of open orbitals -in configuration order
*
*    INTEGER IB_CNOCC_OPEN(MAXOP+1)
*    Offset in list of configuration occupations for configurations 
*    with a given number of  open orbitals - in configuration order
*
*    INTEGER IB_SD_OPEN(MAXOP+1)
*    Offset in list of SD's  for SD
*    with a given number of  open orbitals - in configuration order
*
*    INTEGER IB_CS_OPEN(MAXOP+1)
*    Offset in list of CSF's  for CSF's
*    with a given number of  open orbitals - in configuration order
*
*    INTEGER IB_CM_OPEN(MAXOP+1)
*    Offset in list of CM's for CM's (combinations..)
*    with a given number of  open orbitals - in configuration order

* 
*    Offset for determinants with a given number of 
*    open 
*. 
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
       WRITE(6,*)
       WRITE(6,*) ' =========================='
       WRITE(6,*) ' Output from INFO_CONF_LIST'
       WRITE(6,*) ' =========================='
       WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NACTEL: MAXOP, NACTEL = ', MAXOP, NACTEL
      END IF
*
      LENGTH = 0
      JB_CN = 1
      JB_OCC = 1 
      JB_SD = 1
      JB_CS = 1
      JB_CM = 1
      DO IOP = 0, MAXOP
        IB_CN_OPEN(IOP+1) = JB_CN
        IB_SD_OPEN(IOP+1) = JB_SD
        IB_CS_OPEN(IOP+1) = JB_CS  
        IB_CM_OPEN(IOP+1) = JB_CM  
        IB_CNOCC_OPEN(IOP+1) =  JB_OCC
*
        NCNFOP = NCONF_PER_OPEN(IOP+1,ISYM)
        IF(MOD(NACTEL-IOP,2).EQ.0) THEN
          NOCOB = IOP + (NACTEL - IOP)/2
          JB_OCC = JB_OCC + NOCOB*NCNFOP
          JB_CN = JB_CN + NCNFOP
          JB_SD = JB_SD + NPDTCNF(IOP+1)*NCNFOP
          JB_CS = JB_CS + NPCSCNF(IOP+1)*NCNFOP
          JB_CM = JB_CM + NPCMCNF(IOP+1)*NCNFOP
        END IF
      END DO
*
      LENGTH_LIST = JB_OCC-1
      NCONF_TOT = JB_CN - 1
* 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCONF_PER_OPEN list (input)'
        CALL IWRTMA(NCONF_PER_OPEN(1,ISYM),1,MAXOP+1,1,MAXOP+1) 
        WRITE(6,*) ' Length of configuration list :', LENGTH_LIST 
        WRITE(6,*) ' Total number of configurations : ', NCONF_TOT
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  'IB_CN_OPEN:'
        CALL IWRTMA(IB_CN_OPEN,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*)  'IB_SD_OPEN:'
        CALL IWRTMA(IB_SD_OPEN,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*)  'IB_CS_OPEN:'
        CALL IWRTMA(IB_CS_OPEN,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*)  'IB_CM_OPEN:'
        CALL IWRTMA(IB_CM_OPEN,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*)  'IB_CNOCC_OPEN:'
        CALL IWRTMA(IB_CNOCC_OPEN,1,MAXOP+1,1,MAXOP+1)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_CONF_FOR_CISPC(IOCCLS,NOCCLS,ISYM,IOCCLS_LIST)
*
* Generate configuration information for CI space with symmetry ISYM 
* defined by  the NOCCLS class occupation spaces IOCCLS. 
* The numbering of the occuation classes refers to the occupation classes
* in IOCCLS_LIST
*
* The configuration are ordered according to the number of 
* open orbitals.
*
* It is assumed that GET_DIM_MINMAX_SPACE has been called to set
* the NCONF_PER_OPEN(1,ISYM) arrays - and others
*
* Jeppe Olsen, November 2001
*              Cleaned and shaved, Dec. 2011
* Last modification; Sub configuration approach added; Jeppe Olsen; Apr. 23. 2013
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'crun.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*. Input
      INTEGER IOCCLS(NOCCLS), IOCCLS_LIST(NGAS,*)
*. Local scratch 
      INTEGER NCONF_OPP(MXPOPORB)
*. number of electrons
*. Total number of electrons _ I hope all occ classes have the same 
      NELEC = IELSUM(IOCCLS_LIST(1,IOCCLS(1)),NGAS)
*
      CALL LUCIAQENTER('GENCN')
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================='
        WRITE(6,*)  'GEN_CONF_FOR_CISPC in action '
        WRITE(6,*) ' ============================='
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' NOCCLS, ISYM = ', NOCCLS, ISYM
        WRITE(6,*) ' List of occupation classes '
        CALL IWRTMA(IOCCLS,1,NOCCLS,1,NOCCLS)
      END IF
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' IB_CN_OPEN'
        CALL IWRTMA(IB_CN_OPEN,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*) ' IB_CNOCC_OPEN'
        CALL IWRTMA(IB_CNOCC_OPEN,1,MAXOP+1,1,MAXOP+1)
      END IF
      CALL MEMCHK2('GNCN01')
*
*. memory for storing configuration info 
*
*. At the moment we are 
*      1)  Storing configurations belonging to the 
*          various occupation classes together
*      2)  We are not storing info for different CI spaces together
* We are therefore just using
      ISPC = 1
*. 3 : Array giving start of each occupation class 
*. Scratch memory for setting up configurations 
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GEN_CO')
C?    WRITE(6,*) ' Test: NELEC, NOCOB = ', NELEC, NOCOB
      CALL MEMMAN(KLZSCR,(NOCOB+1)*(NELEC+1),'ADDL  ',1,'ZSCR  ')  !done
C     CALL MEMMAN(KLZ,NOCOB*NELEC*2,'ADDL  ',1,'Z     ')
      CALL MEMMAN(KLOCMIN,NOCOB,'ADDL  ',1,'OCMIN ')  !done
      CALL MEMMAN(KLOCMAX,NOCOB,'ADDL  ',1,'OCMAX ')  !done
*. The KICONF_REO is from lexical to actual order and is therefore
*. common for all symmetries, start by setting to -1
      NCONF_AS = NCONF_ALL_SYM_GN(ISPC)
      IF(NTEST.GE.1000) 
     & WRITE(6,*) ' NCONF_AS in GEN_CONF... =', NCONF_AS
      M1 = -1
C?    WRITE(6,*) ' KICONF_REO(1) = ', KICONF_REO(1)
      CALL ISETVC(int_mb(KICONF_REO(ISPC)),M1,NCONF_AS)
*
* Set up the configurations 
*
      IB_OCCLS = 1
      NCONF_CISPC = 0
      DO JJOCCLS = 1, NOCCLS
        JOCCLS = IOCCLS(JJOCCLS)
*. Save offset to current occupation class 
C             ITOR(WORK,IROFF,IVAL,IELMNT)
c.. dongxia: now we use int_mb, does not switch pointers.
c       CALL ITOR(WORK(KIB_OCCLS(ISYM)),1,IB_OCCLS,JOCCLS)
        int_mb(kib_occls(isym)+joccls-1)=ib_occls
*.Max and min arrays for strings
C            MXMNOC_OCCLS(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,MINOP,NTESTG)
        CALL MXMNOC_OCCLS(int_mb(KLOCMIN),int_mb(KLOCMAX),NGAS,
     &      NOBPT,IOCCLS_LIST(1,JOCCLS),MINOP,NTEST)
C            CONF_GRAPH(IOCC_MIN,IOCC_MAX,NORB,NEL,IARCW,NCONF,ISCR)
*. the arcweights
        IF(NTEST.GE.10000) THEN
          WRITE(6,*) ' KLOCMIN, KLOCMAX, KZCONF, KLZSCR',
     &                 KLOCMIN, KLOCMAX, KZCONF, KLZSCR
        END IF
        CALL CONF_GRAPH(int_mb(KLOCMIN),int_mb(KLOCMAX),NACOB,
     &        NELEC,int_mb(KZCONF),NCONF_P,int_mb(KLZSCR))
*
         IF(JJOCCLS.EQ.1) THEN
            INITIALIZE_CONF_COUNTERS = 1
         ELSE 
            INITIALIZE_CONF_COUNTERS = 0
         END IF
         IDOREO = 1
*. Offset in global orbital list of orbitals in configurations
         IB_ORB = NINOB + 1
*. Lexical addressing for configurations of this type 
         IF(I_DO_SBCNF.EQ.0) THEN
*
*. Set up the good old way
*
C          GEN_CONF_FOR_OCCLS(
C    &         IOCCLS,IB_OCCLS,INITIALIZE_CONF_COUNTERS,
C    &         NGAS,ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,NTORB,NOBPT,
C    &         NCONF_OP,NCONF,IBCONF_REO,IBCONF_OCC,ICONF,
C    &         IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM)
C?         WRITE(6,*) ' IB_OCCLS before GEN_CONF', IB_OCCLS
           CALL GEN_CONF_FOR_OCCLS(IOCCLS_LIST(1,JOCCLS),
     &        IB_OCCLS,INITIALIZE_CONF_COUNTERS,
     &        NGAS,ISYM,MINOP,MAXOP,NSMST,0,NACOB,
     &        NOBPT,NCONF_PER_OPEN(1,ISYM),NCONF_OCCLSL,
     &        IB_CN_OPEN,IB_CNOCC_OPEN,
     &        int_mb(KICONF_OCC(ISYM)),IDOREO,int_mb(KZCONF),
     &        int_mb(KICONF_REO(1)),NCONF_ALL_SYM,IB_ORB,
     &        NCONF_OCCLS_ALLSYM)
          IB_OCCLS = IB_OCCLS + NCONF_OCCLS_ALLSYM
          NCONF_CISPC = NCONF_CISPC + NCONF_OCCLSL
        ELSE
*
* Subconfiguration approach
*
*. Initialize pointers for CI-space
          IF(INITIALIZE_CONF_COUNTERS.EQ.1) THEN
           IZERO = 0
           CALL ISETVC(NCONF_PER_OPEN(1,ISYM),IZERO,MAXOP+1)
           NCONF_ALL_SYM = 0
          END IF 
*. Number of configurations generated before current occupation class in NCONF_OPP
          IF(JOCCLS.EQ.1) THEN
            IZERO = 0
            CALL ISETVC(NCONF_OPP,IZERO,MAXOP+1)
          ELSE
            CALL ICOPVE(NCONF_PER_OPEN(1,ISYM),NCONF_OPP,MAXOP+1)
          END IF
*
          CALL GEN_OCCONF_FOR_OCCLS_FROM_OCSBCLS(
     &         IOCCLS_LIST(1,JOCCLS),ISYM,int_mb(KICONF_OCC(ISYM)),
     &         IB_CNOCC_OPEN,IB_CN_OPEN,NCONF_PER_OPEN(1,ISYM),
     &         MINOP, MAXOP,int_mb(KKOCSBCNF),
     &         int_mb(KNSBCNF),int_mb(KIBSBCNF))
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' NCONF_PER_OPEN after GEN_OCCONF'
            CALL IWRTMA(NCONF_PER_OPEN(1,ISYM),1,MAXOP+1,1,MAXOP+1)
          END IF
*. And the mapping from lexical number (+offset) to actual number
          CALL REO_FOR_CONFS(int_mb(KICONF_OCC(ISYM)),
     &         MAXOP,NSMOB,NCONF_OPP,NCONF_PER_OPEN(1,ISYM),
     &         IB_CN_OPEN,IB_CNOCC_OPEN,IB_OCCLS,int_mb(KZCONF),
     &         NACOB,NELEC,
     &         int_mb(KICONF_REO(1)))
C              REO_FOR_CONFS(IOCC,IBCONF,MAXOP,NSMOB,
C    &             NCONF_OP1,NCONF_OP2,IBCONF_OP,IBOCCCONF_OP,IB_OCCLS,
C    &             IZCONF,NORBT,NELEC,IREO)
*. Update offset to start of configurations of given occupation class
*. Total number of configurations of current occupation class
           CALL NCONF_OCCLS(NOBPT,IOCCLS_LIST(1,JOCCLS),NGAS,0,
     &                      NCONF_OCCLS_ALLSYM,LOCCL)
C               NCONF_OCCLS(NOBPSP,NELPSP,NSPC,MINOP,NCONF,LOCC)
           IB_OCCLS = IB_OCCLS + NCONF_OCCLS_ALLSYM
           NCONF_ALL_SYM = NCONF_ALL_SYM  + NCONF_OCCLS_ALLSYM
        END IF !Switch to subconfiguration approach
*
      END DO
*
      NCONF_CISPC = IELSUM(NCONF_PER_OPEN(1,ISYM),MAXOP+1)
C?    WRITE(6,*) ' WORK(KICONF_REO(1)) as integer '
C?    CALL IWRTMA(WORK(KICONF_REO(1)),1,1,1,1)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*)  ' ============================================ '
        WRITE(6,*)  ' Final results from configuration generator : '
        WRITE(6,*)  ' ============================================ '
        WRITE(6,*)
        WRITE(6,*) ' Number of configurations of correct symmetry ',
     &       NCONF_CISPC
        WRITE(6,*) ' Number of configurations of all symmetries   ',
     &       NCONF_ALL_SYM
        WRITE(6,*) 
     &  ' Number of configurations for various number of open orbs'
        CALL IWRTMA(NCONF_PER_OPEN(1,ISYM),1,MAXOP+1,1,MAXOP+1)
      END IF
      IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Final list of configurations: '
          CALL WRT_CONF_LIST2(
     &    int_mb(KICONF_OCC(ISYM)),IB_CNOCC_OPEN,
     &    NCONF_PER_OPEN(1,ISYM),MAXOP,NCONF_CISPC,NELEC)
C         WRT_CONF_LIST2
C    &    (ICONF,IB_CONF_OCC,NCONF_FOR_OPEN,MAXOP,NCONF,NELEC)
          WRITE(6,*) ' Final reordering of confs, Lex=>Act:'
          CALL IWRTMA(int_mb(KICONF_REO(1)),
     &         1,NCONF_ALL_SYM,1,NCONF_ALL_SYM)
      END IF

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GEN_CO')
      CALL LUCIAQEXIT('GENCN')
      RETURN
      END
      SUBROUTINE MAX_OPEN_ORB(MAXOP,IOCLS,NGAS,NOCLS,NOBPT)
*
* Max number of open orbitals in occupation classes 
*
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
*. Input 
      INTEGER IOCLS(NGAS,NOCLS)
      INTEGER NOBPT(NGAS)
*
      MAXOP = 0
C?    WRITE(6,*) ' NOCLS, NGAS = ', NOCLS, NGAS
      DO JOCLS = 1, NOCLS
        MAXOP_J = 0
        DO IGAS = 1, NGAS
          NEL = IOCLS(IGAS,JOCLS)
          NORB = NOBPT(IGAS)
C?        WRITE(6,*) ' IGAS, NEL, NORB = ', IGAS, NEL, NORB
          MAXOP_IGAS = MIN(NEL,2*NORB-NEL)
          MAXOP_J = MAXOP_J + MAXOP_IGAS
        END DO
        MAXOP = MAX(MAXOP,MAXOP_J)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Max number of unpaired orbitals = ', MAXOP
      END IF
*
      RETURN
      END
      FUNCTION ILEX_FOR_CONF3(ICONF,IEL1,NOCC_ORB,NORB,NEL,IARCW)
*
* Calculate contribution to lexcical address for part of configuration
* The configuration is given by ICONF and starts at electron IEL1
*
* The configuration ICONF of NOCC_ORB orbitals is given as
* ICONF(I) = IORB implies  IORB is singly occupied
* ICONF(I) = -IORB  implies that IORB is doubly occupied 
* 
*
*
      INCLUDE 'implicit.inc'
*. Arcweights for single and doubly occupied arcs
      INTEGER IARCW(NORB,NEL,2)
*. Configuration 
      INTEGER ICONF(NOCC_ORB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ILEX_FOR_CONF3'
        WRITE(6,*) ' Number of occupied orbitals' , NOCC_ORB
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
        WRITE(6,*) ' NORB, NEL, IEL1 = ', NORB, NEL, IEL1
      END IF
*
      IEL = IEL1-1
      ILEX = 1
   
      DO IOCC = 1, NOCC_ORB
       IF(ICONF(IOCC).GT.0) THEN
         IEL = IEL + 1
         ILEX = ILEX + IARCW(ICONF(IOCC),IEL,1)
       ELSE IF(ICONF(IOCC).LT.0) THEN
         IEL = IEL + 2
         ILEX = ILEX + IARCW(-ICONF(IOCC),IEL,2)
       END IF
      END DO
*
       ILEX_FOR_CONF3 = ILEX
*
      IF(NTEST.GE.100) THEN
C?      WRITE(6,*) ' Configuration '
C?      CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
        WRITE(6,*) ' Lexical number = ', ILEX
      END IF
*
      RETURN
      END    
      FUNCTION ILEX_FOR_CONF2(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO,
     &                        IB_LEXADD)
*
* A configuration ICONF of NOCC_ORB orbitals is given
* ICONF(I) = IORB implies  IORB is singly occupied
* ICONF(I) = -IORB  implies that IORB is doubly occupied 
* 
* Find lexical address by adding IB_LEXADD-1 to address obtained from 
* arcweights, and find reordere address if requestes
*
* Differs from CONF (by the inclusion of IB_LEXADD)
*
* IF IDOREO .ne. 0, IREO is used to reorder lexical number 
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
*. Arcweights for single and doubly occupied arcs
      INTEGER IARCW(NORB,NEL,2)
*. Reorder array
      INTEGER IREO(*)
*. Configuration 
      INTEGER ICONF(NOCC_ORB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
C?      WRITE(6,*) ' Number of occupied orbitals' , NOCC_ORB
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
        WRITE(6,*) ' NORB, NEL = ', NORB, NEL
        WRITE(6,*) ' IB_LEXADD = ' , IB_LEXADD 
      END IF
*
      IEL = 0
      ILEX = 1
   
      DO IOCC = 1, NOCC_ORB
       IF(ICONF(IOCC).GT.0) THEN
         IEL = IEL + 1
         ILEX = ILEX + IARCW(ICONF(IOCC),IEL,1)
       ELSE IF(ICONF(IOCC).LT.0) THEN
         IEL = IEL + 2
         ILEX = ILEX + IARCW(-ICONF(IOCC),IEL,2)
       END IF
      END DO
      ILEX = ILEX + IB_LEXADD - 1
*
      IF(IDOREO.NE.0) THEN
       ILEX_FOR_CONF2 = IREO(ILEX)
      ELSE 
       ILEX_FOR_CONF2 = ILEX
      END IF
*
      IF(NTEST.GE.100) THEN
C?      WRITE(6,*) ' Configuration '
C?      CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
        WRITE(6,*) ' Lexical number = ', ILEX
        IF(IDOREO.NE.0)   
     &  WRITE(6,*) ' Reordered number = ', ILEX_FOR_CONF2
      END IF
*
      RETURN
      END    
      FUNCTION ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
*
* A configuration ICONF of NOCC_ORB orbitals are given
* ICONF(I) = IORB implies  IORB is singly occupied
* ICONF(I) = -IORB  implies that IORB is doubly occupied 
* 
* Find lexical address
*
* IF IDOREO .ne. 0, IREO is used to reorder lexical number 
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
*. Arcweights for single and doubly occupied arcs
      INTEGER IARCW(NORB,NEL,2)
*. Reorder array
      INTEGER IREO(*)
*. Configuration 
      INTEGER ICONF(NOCC_ORB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ILEX_FOR_CONF: '
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NORB, NEL = ', NORB, NEL
        WRITE(6,*) ' Number of occupied orbitals' , NOCC_ORB
      END IF
*
      IEL = 0
      ILEX = 1
   
      DO IOCC = 1, NOCC_ORB
       IF(ICONF(IOCC).GT.0) THEN
         IEL = IEL + 1
         ILEX = ILEX + IARCW(ICONF(IOCC),IEL,1)
       ELSE IF(ICONF(IOCC).LT.0) THEN
         IEL = IEL + 2
         ILEX = ILEX + IARCW(-ICONF(IOCC),IEL,2)
       END IF
      END DO
*
      IF(IDOREO.NE.0) THEN
       ILEX_FOR_CONF = IREO(ILEX)
      ELSE 
       ILEX_FOR_CONF = ILEX
      END IF
*
      IF(NTEST.GE.100) THEN
C?      WRITE(6,*) ' Configuration '
C?      CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
        WRITE(6,'(A,I9)') ' Lexical address= ', ILEX
        IF(IDOREO.NE.0)   
     &  WRITE(6,'(A,I9)') ' Reordered address = ', ILEX_FOR_CONF
      END IF
*
      RETURN
      END    
      SUBROUTINE GEN_CONF_FOR_OCCLS(
     &    IOCCLS,IB_OCCLS,INITIALIZE_CONF_COUNTERS,
     &    NGAS,ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,NTORB,NOBPT,
     &    NCONF_OP,NCONF,IBCONF_REO,IBCONF_OCC,ICONF,
     &    IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM,IB_ORB,
     &    NCONF_OCCLS_ALLSYM)
*
* IONLY_NCONF = 1 :
*
* Generate number of configurations of occclass IOCCLS and sym ISYM
*
* IONLY_NCONF = 0 :
*
* Generate number and actual configurations of occclass IOCCLS 
* and sym ISYM
*
*
* Jeppe Olsen, Nov. 2001
*       
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
*
*.. Input
*
*. Number of electrons per gas space 
      INTEGER IOCCLS(NGAS)  
*. Number of orbitals per gasspace 
      INTEGER NOBPT(NGAS)
*. Arc weights for configurations
      INTEGER IZ_CONF(*)
*. Offset for reordering array and occupation array
      INTEGER IBCONF_REO(*), IBCONF_OCC(*)
*
*.. Output
*
*. Number of configurations per number of open shells with given symmetry
      INTEGER NCONF_OP(MAXOP+1)
*. And the actual configurations
      INTEGER ICONF(*)
*. Reorder array : Lex number => Actual number 
      INTEGER IREO(*)
*. Local scratch
      INTEGER JCONF(2*MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GEN_CONF_FOR_OCCLS '
        WRITE(6,*) ' ============================ '
        WRITE(6,*) 
        WRITE(6,*) ' Occupation class in action  '
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Symmetry ', ISYM
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NGAS, MAXOP = ', NGAS, MAXOP
        WRITE(6,*) ' IB_ORB = ', IB_ORB
        WRITE(6,*) ' NOBPT: '
        CALL IWRTMA(NOBPT,1,NGAS,1,NGAS)
        WRITE(6,*) ' IONLY_NCONF = ', IONLY_NCONF
        WRITE(6,*) ' NCONF_ALL_SYM = ', NCONF_ALL_SYM
C?      WRITE(6,*) ' IBCONF_REO(1) = ', IBCONF_REO(1)
C?      WRITE(6,*) ' IBCONF_OCC(1) = ', IBCONF_OCC(1)
        WRITE(6,*) ' IB_OCCLS = ', IB_OCCLS
        WRITE(6,*) ' INITIALIZE_CONF_COUNTERS = ',
     &               INITIALIZE_CONF_COUNTERS
      END IF
*. Total number of electrons 
       NEL = IELSUM(IOCCLS,NGAS)
       IF(NTEST.GE.10000) WRITE(6,*) ' NEL = ', NEL
       IF(INITIALIZE_CONF_COUNTERS.EQ.1) THEN
         IZERO = 0
         CALL ISETVC(NCONF_OP,IZERO,MAXOP+1)
         NCONF_ALL_SYM = 0
       END IF
*. Loop over configurations 
       INI = 1
       NCONF = 0
       ISUM = 0
       NCONF_OCCLS_ALLSYM = 0
 1000  CONTINUE
         IF(NTEST.GE.10000)
     &   WRITE(6,*) ' NEXT_CONF_FOR_OCCLS will be called '
         CALL NEXT_CONF_FOR_OCCLS
     &         (JCONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
         ISUM = ISUM + 1
*
         INI = 0
         IF(NONEW.EQ.0) THEN
*. Check symmetry and number of open orbitals for this space

           IADD = IB_ORB - 1
C               IADCONST(IVEC,IADD, NDIM)
           CALL IADCONST(JCONF,IADD,NEL)
           ISYM_CONF = ISYMST(JCONF,NEL)
           IADD = - IADD
           CALL IADCONST(JCONF,IADD,NEL)
           NOPEN     = NOP_FOR_CONF(JCONF,NEL)
C?         WRITE(6,*) ' Number of open shells ', NOPEN
           NOCOB =  NOPEN + (NEL-NOPEN)/2
           NCONF_ALL_SYM = NCONF_ALL_SYM + 1 
           NCONF_OCCLS_ALLSYM = NCONF_OCCLS_ALLSYM + 1 
           IF(ISYM_CONF.EQ.ISYM.AND.NOPEN.GE.MINOP) THEN
*. A new configuration to be included, reform and save in packed form
             NCONF = NCONF + 1
             NCONF_OP(NOPEN+1) = NCONF_OP(NOPEN+1) + 1
             IF(IONLY_NCONF .EQ. 0 ) THEN
*. Actual and lexical address of this configuration 
               IB_OCC = IBCONF_OCC(NOPEN+1) 
     &                + (NCONF_OP(NOPEN+1)-1)*NOCOB
               IF(NTEST.GE.10000)
     &         WRITE(6,*) ' Offfset for storing conf =', IB_OCC
               CALL REFORM_CONF_OCC(JCONF,ICONF(IB_OCC),NEL,NOCOB,1)
               IF(IDOREO.NE.0) THEN 
C                           ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,
C                                         IARCW,IDOREO,IREO)
                 ILEXNUM =  ILEX_FOR_CONF(ICONF(IB_OCC),NOCOB,NTORB,
     &                      NEL,IZ_CONF,0,IDUM)
                 IF(NTEST.GE.10000) THEN
                   WRITE(6,*) ' Next configuration : '
                   CALL IWRTMA(ICONF(IB_OCC),1,NOCOB,1,NOCOB)
                 END IF
*. Actual address of this configuration 
                 JREO = IBCONF_REO(NOPEN+1) -1 + NCONF_OP(NOPEN+1)
                 IREO(IB_OCCLS-1+ILEXNUM) = JREO
                 IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' LEXCONF: JREO, ILEXNUM = ', JREO, ILEXNUM
                  WRITE(6,*) ' IB_OCCLS, IB_OCCLS-1+ILEXNUM, JREO= ',
     &                         IB_OCCLS, IB_OCCLS-1+ILEXNUM, JREO 
                 END IF
               END IF
             END IF
           END IF
C      IF(ISUM.LE.10) GOTO 1000
           GOTO 1000
         END IF
*        ^ End if nonew = 0
* 
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*)  ' ====================================== '
        WRITE(6,*)  ' Results from configuration generator : '
        WRITE(6,*)  ' ====================================== '
        WRITE(6,*)
        WRITE(6,*) ' Occupation class in action : '
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of configurations of correct symmetry ',
     &       NCONF
        WRITE(6,*) ' Number of configurations of all symmetries   ',
     &       NCONF_ALL_SYM
        WRITE(6,*) 
     &  ' Number of configurations for various number of open orbs'
        CALL IWRTMA(NCONF_OP,1,MAXOP+1,1,MAXOP+1)
        IF(IONLY_NCONF.EQ.0.AND.NTEST.GE.1000) THEN
          WRITE(6,*) 
     &   ' Updated list of configurations (may not be the final...)'
          CALL WRT_CONF_LIST2
     &    (ICONF,IBCONF_OCC,NCONF_OP,MAXOP,NCONF,NEL)
C         WRT_CONF_LIST2
C    &           (ICONF,IB_CONF_OCC,NCONF_FOR_OPEN,MAXOP,NCONF,NELEC)
          WRITE(6,*) 
     &   ' Updated reordering of conf, Lex=>Act (may not be the final'
          CALL IWRTMA(IREO,1,NCONF_ALL_SYM,1,NCONF_ALL_SYM)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE NEXT_CONF_FOR_OCCLS
     &           (ICONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
*
* Obtain next configuration for occupation class
*
* Jeppe Olsen, Nov. 2001
*       
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
*
*. Input
*
*. Number of electrons per gas space 
      INTEGER IOCCLS(NGAS)  
*. Number of orbitals per gasspace 
      INTEGER NOBPT(NGAS)
*. Input and output
      INTEGER ICONF(*)
*. Local scratch
      INTEGER IBORB(MXPNGAS), ICONF_GAS(MXPORB)
      INTEGER IBEL(MXPNGAS)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' NEXT_CONF_FOR_OCCLS entered '
       WRITE(6,*) ' Occupation class in action '
       CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
      END IF
*. Total number of electrons 
      NEL = IELSUM(IOCCLS,NGAS)
      IF(NTEST.GE.1000) 
     &WRITE(6,*) ' NEXT_CONF ... NEL, NGAS = ', NEL, NGAS
*. Offset for orbitals and electrons
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          IBORB(IGAS) = 1
          IBEL(IGAS)  = 1
        ELSE
          IBORB(IGAS) = IBORB(IGAS-1)+NOBPT(IGAS-1)
          IBEL(IGAS)  = IBEL(IGAS-1) + IOCCLS(IGAS-1)
        END IF
      END DO
*
      NONEW = 1
      
      IF(INI.EQ.1) THEN 
*
*. Initial configuration 
*
        NONEW = 0
        INI_L = 1
        NONEW_L = 0
        DO IGAS = 1, NGAS
*. Initial configuration for this GASSPACE 
          NEL_GAS = IOCCLS(IGAS)
          NORB_GAS = NOBPT(IGAS)
C?        WRITE(6,*) ' IGAS, NEL_GAS, NORB_GAS = ',
C?   &                 IGAS, NEL_GAS, NORB_GAS
          CALL NXT_CONF(ICONF_GAS,NEL_GAS,NORB_GAS,INI_L,NONEW_L)
          IF(NONEW_L.EQ.1) THEN
             NONEW = 1
             GOTO 1001
          ELSE 
             JBEL   = IBEL(IGAS)
             JBORB  = IBORB(IGAS)
             JEL = IOCCLS(IGAS)
             JORB  = NOBPT(IGAS)
             CALL REFORM_CONF_FOR_GAS
     &            (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,2)
C                 (ICONF_GAS,ICONF,IBORB,IBEL,NORB,NEL,IWAY)
          END IF
        END DO
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Initial configuration '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        END IF
*
      ELSE   
*
*. Next configuration 
*
*. Loop over GAS spaces and find first GASspace where a new configuration 
*. could be obtained
        DO IGAS = 1, NGAS
C?        WRITE(6,*) ' IGAS = ', IGAS
*. Remove the offsets for this space 
          JBEL   = IBEL(IGAS)
          JBORB  = IBORB(IGAS) 
          JEL = IOCCLS(IGAS)
          JORB = NOBPT(IGAS)
C?        WRITE(6,*) ' JBEL, JBORB, JEL, JORB = ',
C?   &                 JBEL, JBORB, JEL, JORB
          CALL REFORM_CONF_FOR_GAS
     &         (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,1)
*. Generate next configuration for this space 
          INI_L = 0
          CALL NXT_CONF(ICONF_GAS,JEL,JORB,INI_L,NONEW_L)
          IF(NONEW_L.EQ.0) THEN
            NONEW = 0
*. Configuration in space IGAS, was increased. Copy this and reset configurations
*. in previous gasspaces to initial form 
            CALL REFORM_CONF_FOR_GAS
     &           (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,2)
*
            DO JGAS = 1, IGAS-1
              JBEL   = IBEL(JGAS)
              JBORB  = IBORB(JGAS) 
              JEL = IOCCLS(JGAS)
              JORB = NOBPT(JGAS)
              CALL REFORM_CONF_FOR_GAS
     &             (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,1)
              INI_L = 1
              CALL NXT_CONF(ICONF_GAS,JEL,JORB,INI_L,NONEW_L)
C                                     NEL_GAS,NORB_GAS,INI_L,NONEW_L)
              CALL REFORM_CONF_FOR_GAS
     &             (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,2)
            END DO
*. Get out of the loop 
            GOTO 1001
          END IF
        END DO
*       ^ End of loop over gasspaces
      END IF
*     ^ End if swith between initialization/next
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
          WRITE(6,*) ' No new configuration '
        ELSE 
          WRITE(6,*) ' New configuration '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE NXT_CONF(ICONF,NEL,NORB,INI,NONEW)
*
* Next configuration of NEL electrons distributed in NORB orbitals
*
* A configuration is stored as the occupied orbitals 
* in nonstrict ascending order - two consecutivw orbitals are allowed 
* to be identical
* allowing two 
*
* IF INI = 1 : Generate initial configuration
*    NONOEW = 1 : No new configuration could be generated
*
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
      INTEGER ICONF(NEL)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input configuration to NXT_CONF '
        CALL IWRTMA(ICONF,1,NEL,1,NEL)
        WRITE(6,*) ' NEL, NORB = ',  NEL, NORB
      END IF
      IF(INI.EQ.1) THEN
*. Check that NEL electrons can be distributed in NORB orbitals 
        IF(NEL.LE.2*NORB) THEN
          NONEW = 0
          N_DOUBLE = NEL/2
          DO I = 1, N_DOUBLE
            ICONF(2*I-1) = I
            ICONF(2*I)   = I
          END DO
          IF(2*N_DOUBLE.NE.NEL) ICONF(2*N_DOUBLE+1) = N_DOUBLE + 1
        ELSE 
          NONEW = 1
        END IF
*
      ELSE IF(INI.EQ.0) THEN
*
        IADD  = 1
        IEL = 0
*. Increase orbital number of next electron
 1000   CONTINUE
          IEL = IEL + 1
*. Can orbital number be increased for electron IEL ?
          INCREASE = 0
          IF(IEL.LT.NEL) THEN
            IF(ICONF(IEL).LT.ICONF(IEL+1)-1)  INCREASE = 1
            IF(ICONF(IEL).EQ.ICONF(IEL+1)-1) THEN
*. If ICONF(IEL) is increased, ICONF(IEL) = ICONF(IEL+1), check if this is ok
              IF(IEL.EQ.NEL-1) THEN 
                 INCREASE = 1
              ELSE IF (ICONF(IEL+1).NE.ICONF(IEL+2)) THEN
                 INCREASE = 1
              END IF
            END IF
          ELSE 
C-jwk new            IF(ICONF(IEL).LT.NORB) THEN 
            IF(IEL .EQ. NEL .AND. ICONF(IEL) .LT. NORB) THEN
              INCREASE = 1
            ELSE
*. Nothing more to do
              NONEW = 1
              GOTO 1001
            END IF
          END IF
*
          IF(INCREASE.EQ.1) THEN
*. Increase orbital for elec IEL
            NONEW = 0
            ICONF(IEL) = ICONF(IEL)+1
*. Minimize orbital occupations
            NDOUBLE = (IEL-1)/2
            DO JORB = 1, NDOUBLE
              ICONF(2*JORB-1) = JORB
              ICONF(2*JORB  ) = JORB
            END DO
            IF(2*NDOUBLE.LT.IEL-1) ICONF(IEL-1) = NDOUBLE+1
            IADD = 0
          END IF
        IF(IADD.EQ.1)  GOTO 1000
      END IF
*     ^ End if INI = 0
*
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
          WRITE(6,*) ' No new configurations '
          WRITE(6,*) ' Input configuration '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        ELSE 
          WRITE(6,*) ' Next configurations '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_CONF_FOR_GAS(ICONF_GAS,ICONF,IBORB,IBEL,
     &                                NORB,NEL,IWAY)
*
* Reform between local and global numbering of 
* configuration for given GAS space
*
* IWAY = 1 : Global => Local
* IWAY = 2 : Local => GLobal
*
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
*
      INTEGER ICONF_GAS(NORB)
      INTEGER ICONF(*)
*
      IF(IWAY.EQ.1) THEN
        DO IEL = 1, NEL
          ICONF_GAS(IEL) = ICONF(IBEL-1+IEL)  - IBORB + 1
        END DO
      ELSE IF (IWAY.EQ.2) THEN
        DO IEL = 1, NEL
          ICONF(IBEL-1+IEL) = ICONF_GAS(IEL) + IBORB - 1
        END DO
      ELSE 
        WRITE(6,*) ' Problem in REFORM_CONF ... , IWAY = ', IWAY
        STOP       ' Problem in REFORM_CONF ... , IWAY = '
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Global => Local reform of conf '
        ELSE
          WRITE(6,*) ' Local => Global reform of conf '
        END IF
        WRITE(6,*) ' ICONF_GAS : '
        CALL IWRTMA(ICONF_GAS,1,NEL,1,NEL) 
        WRITE(6,*) ' Accessed part of ICONF ' 
        CALL IWRTMA(ICONF,1,IBEL-1+NEL,1,IBEL-1+NEL)
      END IF
*
      RETURN
      END
      FUNCTION NCL_FOR_CONF(ICONF,NEL)
*
* A configuration is given as a nonstrict ascending sequence of occupied 
* occupied orbitals. Find number of Doubly occupied orbitals
*
* Jeppe Olsen, Aug. 2004, NOP_FOR_CONF with a twist
*
      INCLUDE 'implicit.inc'
      INTEGER ICONF(NEL)
*. Loop over electrons 
      NOPEN = 0
      NCL = 0
      IEL = 1
 1000 CONTINUE
        IF(IEL.LT.NEL) THEN
          IF(ICONF(IEL).NE.ICONF(IEL+1)) THEN
           NOPEN = NOPEN + 1
           IEL = IEL + 1
          ELSE IF (ICONF(IEL).EQ.ICONF(IEL+1) ) THEN
           NCL = NCL + 1
           IEL = IEL + 2
          END IF
        END IF
*
        IF(IEL.EQ.NEL) THEN
*. The last orbital is not identical to any later orbitals so
         NOPEN = NOPEN+1
         IEL = IEL + 1
        END IF
      IF(IEL.LT.NEL) GOTO 1000
*
      NCL_FOR_CONF = NCL
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NEL,1,NEL)
        WRITE(6,*) ' Number of closed orbitals = ', NCL_FOR_CONF
      END IF  
*
      RETURN
      END
      FUNCTION NOP_FOR_CONF(ICONF,NEL)
*
* A configuration is given as a nonstrict ascending sequence of occupied 
* occupied orbitals. Find number of single occupied orbitals
*
* Jeppe Olsen, Nov. 2001
*
      INCLUDE 'implicit.inc'
      INTEGER ICONF(NEL)
*. Loop over electrons 
      NOPEN = 0
      NCL = 0
      IEL = 1
 1000 CONTINUE
        IF(IEL.LT.NEL) THEN
          IF(ICONF(IEL).NE.ICONF(IEL+1)) THEN
           NOPEN = NOPEN + 1
           IEL = IEL + 1
          ELSE IF (ICONF(IEL).EQ.ICONF(IEL+1) ) THEN
           NCL = NCL + 2
           IEL = IEL + 2
          END IF
        END IF
*
        IF(IEL.EQ.NEL) THEN
*. The last orbital is not identical to any later orbitals so
         NOPEN = NOPEN+1
         IEL = IEL + 1
        END IF
      IF(IEL.LT.NEL) GOTO 1000
*
      NOP_FOR_CONF = NOPEN
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NEL,1,NEL)
        WRITE(6,*) ' Number of open orbitals = ', NOP_FOR_CONF
      END IF  
*
      RETURN
      END
      SUBROUTINE REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY)
*
* Reform between two ways of writing occupations
*
* IOCC_EXP : Occupation in expanded form, i.e. the orbital for each 
*            electron is given
*
* IOCC_PCK  : Occupation is given in packed form, i.e. each occupied 
*             orbitals is given once, and a negative index indicates
*             a double occupattion 
*
* IWAY = 1 Expanded to Packed form
* IWAY = 2 Packed to expanded form
*
* Jeppe Olsen, Nov. 2001
*
      INCLUDE 'implicit.inc'
*. Input/Output
      INTEGER IOCC_EXP(NEL),IOCC_PCK(NOCOB)
*
      IF(IWAY.EQ.1) THEN
*
*. Expanded => Packed form 
*
*. Loop over electrons
        IEL = 1
        IOCC = 0
 1000   CONTINUE
C         IEL = IEL + 1
          IF(IEL.LT.NEL) THEN
            IF(IOCC_EXP(IEL).EQ.IOCC_EXP(IEL+1)) THEN
              IOCC = IOCC + 1
              IOCC_PCK(IOCC) = -IOCC_EXP(IEL)
              IEL = IEL + 2
            ELSE 
              IOCC = IOCC + 1
              IOCC_PCK(IOCC) =  IOCC_EXP(IEL)
              IEL = IEL + 1
            END IF
          ELSE 
*. Last occupation was not identical to previous, so single occupied
            IOCC = IOCC + 1
            IOCC_PCK(IOCC) =  IOCC_EXP(IEL)
            IEL = IEL + 1
          END IF
        IF(IEL.LE.NEL) GOTO 1000
        NOCOB = IOCC
*. Zero electrons gives problems in above loop, so
        IF(NEL.EQ.0) NOCOB = 0

*
      ELSE IF( IWAY.EQ.2) THEN
*
* Packed to expanded form 
*
        IEL = 0
        DO IORB = 1, NOCOB
          IF(IOCC_PCK(IORB).LT.0) THEN
            JORB = - IOCC_PCK(IORB)
            IEL = IEL +1 
            IOCC_EXP(IEL) = JORB
            IEL = IEL + 1
            IOCC_EXP(IEL) = JORB
          END IF
        END DO
      ELSE 
        WRITE(6,*) ' REFORM_CONF... in error, IWAY = ', IWAY
        STOP       ' REFORM_CONF... in error, IWAY '
      END IF
*     ^ End of IWAY switch 
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reforming form of configuration '
        IF(IWAY.EQ.1) THEN
           WRITE(6,*) ' Expanded to packed form '
        ELSE 
           WRITE(6,*) ' Packed to expanded form '
        END IF
        WRITE(6,*) ' IOCC_EXP : '
        CALL IWRTMA(IOCC_EXP,1,NEL,1,NEL)
        WRITE(6,*) ' IOCC_PCK : '
        CALL IWRTMA(IOCC_PCK,1,NOCOB,1,NOCOB)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_CONF_LIST(ICONF,NCONF_FOR_OPEN,MAXOP,NCONF,NELEC)
*
* Write list of configurations, given in packed form
*
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
*   
      INTEGER ICONF(*), NCONF_FOR_OPEN(MAXOP+1)
*
      IB = 1
      DO IOPEN = 0, MAXOP
        NCONF_OP = NCONF_FOR_OPEN(IOPEN+1)
        IF(NCONF_OP.NE.0) THEN
          WRITE(6,'(A,I3,A,I6)') 
     &    ' Number of configurations with ', IOPEN, 
     &               ' open orbitals is ', NCONF_OP
*
          NOCC_ORB = IOPEN + (NELEC-IOPEN)/2
          DO JCONF = 1, NCONF_OP
            CALL IWRTMA(ICONF(IB),1,NOCC_ORB,1,NOCC_ORB)
            IB = IB + NOCC_ORB
          END DO
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE WRT_CONF(ICONF,L)
*
* Write configuration ICONF given in packed form 
*
* Jeppe Olsen, Nov. 2001
*
      INCLUDE 'implicit.inc'
      INTEGER ICONF(L)
      WRITE(6,*)
      WRITE(6,'(A,20(1X,I4),/(8X,20(1X,I4)))') 
     &' Orbitals : ',(ICONF(I), I=1,L)
*. Add somwthing about occupation 
C     WRITE(6,'(A,20I3,/(8X,20I3))') 
C     ' Occ      : ',(ICONF(I), I=1,L)
      RETURN
      END
      SUBROUTINE ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,
     &           IDET_OC,IDET_SP,ISIGN)
*
* An alpha string (IA) and a betastring (IB) is given. 
* Combine these two strings to give an determinant with 
* orbitals in ascending order. For doubly occupied orbitals
* the alphaorbital is given first. 
* The output is given as IDET_OC : Orbital occupation (configuration ) 
*                        IDET_SP : Spin projections

* The phase required to change IA IB into IDET is computes as ISIGN
*
* Jeppe Olsen, November 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IA_OC(NAEL),IB_OC(NBEL)
*. Output
      INTEGER IDET_OC(NAEL+NBEL)
      INTEGER IDET_SP(NAEL+NBEL)
*  
      NEXT_AL = 1
      NEXT_BE = 1
      NEXT_EL = 0
      ISIGN = 1
*. Loop over next electron in outputstring
      DO NEXT_EL = 1, NAEL+NBEL
       IF(NEXT_AL.LE.NAEL.AND.NEXT_BE.LE.NBEL) THEN
*
         IF(IA_OC(NEXT_AL).LE.IB_OC(NEXT_BE)) THEN
*. Next electron is alpha electron
           IDET_OC(NEXT_EL) = IA_OC(NEXT_AL)
           IDET_SP(NEXT_EL) = +1
           NEXT_AL = NEXT_AL + 1
         ELSE
*. Next electron is beta electron
           IDET_OC(NEXT_EL) = IB_OC(NEXT_BE)
           IDET_SP(NEXT_EL) = -1
           NEXT_BE = NEXT_BE + 1
           ISIGN = ISIGN*(-1)**(NAEL-NEXT_AL+1) 
         END IF
       ELSE IF(NEXT_BE.GT.NBEL) THEN
*. Next electron is alpha electron
           IDET_OC(NEXT_EL) = IA_OC(NEXT_AL)
           IDET_SP(NEXT_EL) = +1
           NEXT_AL = NEXT_AL + 1
       ELSE IF(NEXT_AL.GT.NAEL) THEN
*. Next electron is beta electron
           IDET_OC(NEXT_EL) = IB_OC(NEXT_BE)
           IDET_SP(NEXT_EL) = -1
           NEXT_BE = NEXT_BE + 1
           ISIGN = ISIGN*(-1)**(NAEL-NEXT_AL+1) 
       END IF
      END DO
*    ^ End of loop over electrons in outputlist
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ABSTR to ORDSTR : '
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' Input alpha and beta strings '
        CALL IWRTMA(IA_OC,1,NAEL,1,NAEL)
        CALL IWRTMA(IB_OC,1,NBEL,1,NBEL)
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(IDET_OC,1,NAEL+NBEL,1,NAEL+NBEL)
        WRITE(6,*) ' Spin projections '
        CALL IWRTMA(IDET_SP,1,NAEL+NBEL,1,NAEL+NBEL)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_W_FOR_PTYPES(MINOP,MAXOP,NEL,MS2,KLZ_FOR_TYPES)
*
* Generate arcweights for addressing the various determinants of a 
* given configuration
*
* A determinant for a given configuraration is defined by the 
* open alpha-electrons, and a standard reverse lexical ordering 
* index is set up for addressing these. 
*
* Separate indeces are set up for the various number of open shells
* Info is not created for numbers of open electrons  that are 
* inconsistent with the total number of electrons
*
* The arrays are delivered in KLZ_FOR_TYPES, such that 
* the determinants for a configuration with IOPEN shells 
* is delivered in WORK(KLZ_FOR_TYPES(IOPEN))
* The pointers as well as the addresses are determined in the routine
*
* Jeppe Olsen, November 2001
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INTEGER KLZ_FOR_TYPES(MAXOP)
*
      NTEST = 00
C?    WRITE(6,*) ' MINOP, MAXOP entering GEN_W ...', MINOP, MAXOP
*. Set up pointers KLZ_FOR_OPEN
      IONE = 1
      CALL ISETVC(KLZ_FOR_TYPES,IONE,MINOP-1)
      DO IOPEN = MINOP, MAXOP
        IF(MOD(IOPEN-MS2,2).EQ.0) THEN
          NAEL = (IOPEN+MS2)/2
          NBEL = (IOPEN-MS2)/2
          L = NAEL*IOPEN
          CALL MEMMAN(KLZ_FOR_TYPES(IOPEN+1),L,'ADDL  ',1,'Z_F_TY')
        ELSE
          KLZ_FOR_TYPES(IOPEN) = 1
        END IF
      END DO
*. Local scratch for setting up arrays
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SCR_GE')
      CALL MEMMAN(KLMAX,MAXOP,'ADDL  ',1,'L_MAXO')
      CALL MEMMAN(KLMIN,MAXOP,'ADDL  ',1,'L_MINO') 
      L = (MAXOP+1)*(MAXOP+1)
      CALL MEMMAN(KLY,L,'ADDL  ',1,'L_Y   ')
*. And then the arrays
      DO IOPEN = MINOP, MAXOP
        IF(MOD(IOPEN-MS2,2).EQ.0) THEN
          NAEL = (IOPEN+MS2)/2
          NBEL = (IOPEN-MS2)/2
          CALL MXMNOC_SPGP(int_mb(KLMIN),int_mb(KLMAX),1,IOPEN,NAEL,
     &                     NAEL,NINOB,NTEST)
          CALL GRAPW(int_mb(KLY),int_mb(KLZ_FOR_TYPES(IOPEN+1)),
     &               int_mb(KLMIN),
     &               int_mb(KLMAX),IOPEN,NAEL,NINOB,NTEST)
C     MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),NORBTP,NORBFTP,NELFTP,
C    &                 NEL,NINOB,NTEST)
C              CALL GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL,NTEST)
        END IF
      END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SCR_GE')
*
      RETURN
      END
      SUBROUTINE ITOR(WORK,IROFF,IVAL,IELMNT)
*
* An integer array is stored in real array WORK,
* starting from WORK(IROFF). Give IELEMNT  value IVAL
* IELMNT of this array
*
      INTEGER WORK(*)
*
      INCLUDE 'irat.inc'
*. offset when work is integer array
      IIOFF = 1 + IRAT * (IROFF-1)
      WORK(IIOFF-1+IELMNT) = IVAL
C     IFRMR = WORK(IIOFF-1+IELMNT)
*
      RETURN
      END
      SUBROUTINE MXMNOC_OCCLS(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,
     &                  MINOP,NTESTG)
*
* Construct accumulated MAX and MIN arrays for an occupation class
*
* MINOP ( Smallest allowed number of open orbitals) added 
* April2, 2003, JO (modified by JWK, April - June 2003)
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Output
      DIMENSION  MINEL(*),MAXEL(*)
*. Input
      INTEGER NORBFTP(*),NELFTP(*)
*. Local scratch added April 2, 2003
*
      INCLUDE 'mxpdim.inc'
C     INTEGER MINOP_ORB(MXPORB), MINOP_GAS(MXPNGAS), MAXOP_GAS(MXPNGAS)
      INTEGER MINOP_GAS(MXPNGAS), MAXOP_GAS(MXPNGAS)
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.10000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============'
        WRITE(6,*) ' MXMNOC_OCCLS'
        WRITE(6,*) ' ============'
        WRITE(6,*)
C?      WRITE(6,*) ' MINOP  = ', MINOP
C?      WRITE(6,*) ' NORBTP = ', NORBTP
C?      WRITE(6,*) ' NORBFTP : '
C?      CALL IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
      END IF
*. Well
      NGAS = NORBTP
*
*. Largest number of unpaired electrons in each gas space 
*
      DO IGAS = 1, NGAS
        MAXOP_GAS(IGAS) = MIN(NELFTP(IGAS),2*NORBFTP(IGAS)-NELFTP(IGAS))
      END DO
*
*. Smallest number of electrons in each GAS space 
*
*. 1 : Just based on number of electrons in each space 
      DO IGAS = 1, NGAS
        IF(MOD(NELFTP(IGAS),2).EQ.1) THEN
          MINOP_GAS(IGAS) = 1
        ELSE
          MINOP_GAS(IGAS) = 0
        END IF
      END DO
*. 2 : the total number of open orbitals should be MINOP, this puts 
*. also a constraint on the number of open orbitals
*
*. The largest number of open orbitals, all spaces
      MAXOP_T = IELSUM(MAXOP_GAS,NGAS)
      DO IGAS = 1, NGAS
*. Max number of open orbitals in all spaces except IGAS
       MAXOP_EXL = MAXOP_T - MAXOP_GAS(IGAS)
       MINOP_GAS(IGAS) = MAX(MINOP_GAS(IGAS),MINOP-MAXOP_EXL)
       IF (MOD(NELFTP(IGAS)-MINOP_GAS(IGAS),2) .EQ. 1) THEN
          MINOP_GAS(IGAS) = MINOP_GAS(IGAS) + 1
       ENDIF
      END DO
*. We now have the min and max number of open shells per occls,
*. Find the corresponding min and max number accumulated electrons, 
*
* The Max occupation is obtained by occupying in max in the 
* first orbitals 
* The Min occupation is obtained by occopying max in the 
* last orbitals.
*
      NEL_INI = 0
      IBORB   = 1
      DO IGAS = 1, NGAS
        NELEC = NELFTP(IGAS)
        MAX_DOUBLE = (NELEC-MINOP_GAS(IGAS))/2
*
* If you are in a situation with no electrons to spare
*
        IF (NELEC .EQ. 0) THEN
           DO IORB = 1,NORBFTP(IGAS)
              IF (IORB+IBORB-1 .EQ. 1) THEN
                 MINEL(IORB+IBORB-1) = 0
                 MAXEL(IORB+IBORB-1) = 0
              ELSE
                 MINEL(IORB+IBORB-1) = MINEL(IORB+IBORB-2)
                 MAXEL(IORB+IBORB-1) = MAXEL(IORB+IBORB-2)
              END IF
           END DO
           GOTO 10
        END IF
*
* The min number of electrons
*
*. Doubly occupy the last MAX_DOUBLE orbitals
C Start Jesper !!!
C       IF (NORBFTP(IGAS)-MAX_DOUBLE .LE. 0
C    &        .AND. MINOP_GAS(IGAS) .GT. 0) CALL QUIT(9999)
C End Jesper !!!
        IORB_START = MAX(1,NORBFTP(IGAS)-MAX_DOUBLE)
        DO IORB = IORB_START,NORBFTP(IGAS)
           MINEL(IORB+IBORB-1) = 
     &           NEL_INI + NELEC - 2*(NORBFTP(IGAS)-IORB)
C?        write(6,*) ' 1 IORB+IBORB-1, MINEL() ',
C?   &    IORB+IBORB-1,  MINEL(IORB+IBORB-1) 
        END DO
*. Singly occupy 
        DO IORB = NORBFTP(IGAS)-MAX_DOUBLE-1,1,-1
           MINEL(IORB+IBORB-1) = MAX(NEL_INI,MINEL(IORB+IBORB-1+1)-1)
C?        write(6,*) ' 2 IORB+IBORB-1, MINEL() ',
C?   &    IORB+IBORB-1,  MINEL(IORB+IBORB-1)
        END DO
*
*. The max number of electrons 
*
       DO IORB = 1, MAX_DOUBLE
         MAXEL(IORB+IBORB-1) = NEL_INI + 2*IORB 
       END DO 
       DO IORB = MAX_DOUBLE+1, NORBFTP(IGAS)  
         IF (IORB+IBORB-1 .EQ. 1) THEN
           MAXEL(IORB+IBORB-1) = 1
         ELSE 
           MAXEL(IORB+IBORB-1)=MIN(NEL_INI+NELEC,MAXEL(IORB+IBORB-2)+1)
         ENDIF
       END DO
  10   CONTINUE
       NEL_INI = NEL_INI + NELFTP(IGAS)
       IBORB = IBORB + NORBFTP(IGAS)
      END DO
*
      IF( NTEST .GE. 10000 ) THEN
        NORB = IELSUM(NORBFTP,NORBTP)
        WRITE(6,*) ' MINEL and MAXEL: '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
* routines for spinadapting CC amplitudes 
*
      SUBROUTINE S2_MAT_FOR_DETS(IDET,NDET,NOP,MS2,IZ,IREO,S2MAT,
     &                           IWORK,XWORK)
*
* Set up S^2 matrix in basis defined by determinants IDET
*
* Jeppe Olsen, Jan. 2002
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IDET(NOP,NDET), IREO(*), IZ(*)  
*. Output 
      DIMENSION S2MAT(NDET,NDET)
*. Local scratch  : Should hold info for all dets in configuration.
      INTEGER IWORK(NOP,*)  
      DIMENSION XWORK(*)

*
      NAL = (NOP + MS2)/2
      ZERO = 0.0D0
      CALL SETVEC(S2MAT,ZERO,NDET**2)
      DO JDET = 1, NDET
*. S2 times determinant
C     S2_DET(IDET_IN,NOPEN,IDET,NDET,CONST,IFLAG)
        CALL S2_DET(IDET(1,JDET),NOP,IWORK,NDET_OUT,XWORK,0)
        DO KDET = 1, NDET_OUT
C                IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
          KNUM = IZNUM_PTDT(IWORK(1,KDET),NOP,NAL,IZ,IREO,1)
          IF(KNUM.GT.0) S2MAT(KNUM,JDET) = XWORK(KDET)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The S2 matrix '
        CALL WRTMAT(S2MAT,NDET,NDET,NDET,NDET)
      END IF
*
      RETURN
      END
C       REO_SD(WORK(KLSD),NOP,NOP_AEL,NDET,NDET_CNF,WORK(KLZ),WORK(KLREO))
      SUBROUTINE REO_SD(IDET,NOP,NOP_AL,NDET,NDET_CNF,IZ,IREO)
*
* A set of open orbitals with various spinprojections is given
* Obtain reorder array Lexical to actual
*. A nonpresent SD is flagged by a zero
*
*. Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IDET(NOP,NDET), IZ(*)
*. Output
      INTEGER IREO(NDET_CNF)
*
      NTEST = 00
*
      IDUM = 0
      IZERO = 0
      CALL ISETVC(IREO,IZERO,NDET_CNF)
      DO JDET = 1, NDET
C     IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
*. Lexical number of det without reordering
        JZ = IZNUM_PTDT(IDET(1,JDET),NOP,NOP_AL,IZ,IDUM,0)
        IREO(JZ) = JDET
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Determinant : '
          CALL IWRTMA(IDET(1,JDET),1,NOP,1,NOP)
          WRITE(6,*) ' Lex and actual address ', JZ, JDET
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) '  Determinants : '
        CALL IWRTMA(IDET,NOP,NDET,NOP,NDET)
        WRITE(6,*) ' Reorder array : Lex order => actual order '
        CALL IWRTMA(IREO,NDET_CNF,1,NDET_CNF,1)
      END IF
*
      RETURN 
      END 
      SUBROUTINE IALBE_TO_MAXOP(IAL,IBE,IMAXOP,IMS2)
*
* An spinorbital excitation is defined as an alpha and an beta string
*
* Obtain Max number of open orbitals per gas and MS2 in each gas
*
* Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IAL(NGAS), IBE(NGAS)
*. Output
      INTEGER IMAXOP(NGAS), IMS2(NGAS)
*
      DO IGAS = 1, NGAS
        IMS2(IGAS) = IAL(IGAS)-IBE(IGAS)
        NABEL = IAL(IGAS) + IBE(IGAS)
        IMAXOP(IGAS) = 
     &  MIN(NABEL,(2*NOBPT(IGAS)-NABEL)) 
      END DO
*
      NTEST =  00
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' IALBE_TO_MAX .... '
        WRITE(6,*) ' Input IAL and IBE '
        CALL IWRTMA(IAL,1,NGAS,1,NGAS)
        CALL IWRTMA(IBE,1,NGAS,1,NGAS)
        WRITE(6,*) ' Output IMAXOP, IMS2 : '
        CALL IWRTMA(IMAXOP,1,NGAS,1,NGAS)
        CALL IWRTMA(IMS2,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END 
      SUBROUTINE S2_DET(IDET_IN,NOPEN,IDET,NDET,CONST,IFLAG)
*
* A determinant is given as NOPEN open shells in IDET_IN
*
* Apply S2 on this determinant
* S2 = S+S- + Sz(Sz-1) 
*
* IFLAG = 1 => Only number of determinants
*       = 0 => the actual determinants
*
* Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input 
      INTEGER IDET_IN(NOPEN)
*. Output
      INTEGER IDET(NOPEN,*) 
      DIMENSION CONST(*)
*. Local scratch : 
      INTEGER IORB_AL(MXPORB),IORB_BE(MXPORB)
      
*. Number of alpha electrons 
      NALPHA = 0
      DO IOP = 1, NOPEN
        IF(IDET_IN(IOP).EQ.1) NALPHA = NALPHA + 1
      END DO
      NBETA = NOPEN - NALPHA 
*. Number of terms : Spinflips between alpha and beta and a term times original
      NDET = 1 + NALPHA*NBETA
*
      IF(IFLAG.EQ.0) THEN
*. Set up the actual combinations 
*1 : The input det, constant is MS*(MS-1) + NALPHA
       NDET = 1
       CALL ICOPVE(IDET_IN,IDET(1,1),NOPEN)
       XMS = 0.5D0*(NALPHA-NBETA)
C?     WRITE(6,*) ' Nalpha Nbeta XMS ', NALPHA, NBETA, XMS
       CONST(1)  = XMS*(XMS-1.0D0) + DBLE(NALPHA)
*. Find the alpha and beta orbitals 
       IALPHA = 0
       IBETA = 0
       DO JORB = 1, NOPEN
         IF(IDET_IN(JORB).EQ.1) THEN
           IALPHA = IALPHA + 1
           IORB_AL(IALPHA) = JORB
         ELSE
           IBETA = IBETA + 1
           IORB_BE(IBETA) = JORB
         END IF
       END DO
*
       DO IALPHA = 1, NALPHA
         DO IBETA = 1, NBETA
           NDET = NDET + 1
           CALL ICOPVE(IDET_IN,IDET(1,NDET),NOPEN)
           IDET(IORB_AL(IALPHA),NDET) = 0
           IDET(IORB_BE(IBETA),NDET)  = 1
           CONST(NDET) = 1.0D0
         END DO
       END DO
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Output from S2_DET '
         WRITE(6,*) ' Number of determinants ', NDET
         IF(IFLAG.EQ.0) THEN
           WRITE(6,*) ' And the actual determinants '
           DO JDET = 1, NDET
             WRITE(6,*) CONST(JDET) 
             CALL IWRTMA(IDET(1,JDET),1,NOPEN,1,NOPEN)
           END DO
         END IF
      END IF
*
      RETURN
      END 
      SUBROUTINE DET_FOR_SPOC(NOP_GAS,IMS2_GAS,NGAS,NDET,IDET,
     &                        IBDET,IFLAG)
*
* A given spinorbital occupation is given in the form of NOP_GAS(IGAS)
* open orbitals in gasspace IGAS. The spinprojection in GASpace IGAS
* must be IMS_GAS. Obtain determinants with these spinprojections.
*
* An determinant is an array of 1 and 0,  1 => alpha 0=>  beta
*
* IFLAG = 1 => Only number of determinants
* IFLAG = 0 => Also the actual determinants 
*
* Jeppe Olsen, December 2001
*
* 
      INCLUDE 'implicit.inc'
      INTEGER ADD
      INCLUDE 'mxpdim.inc'
*. input
      INTEGER NOP_GAS(NGAS), IMS2_GAS(NGAS)
*. Output
      INTEGER IDET(*)
*. local scratch
      INTEGER IWORK(MXPORB)
*. Total number of open orbitals
      NTEST = 00
      NOPEN = IELSUM(NOP_GAS,NGAS)
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' Number of open orbitals = ', NOPEN
      MX=2 ** NOPEN
      NDET = 0
      IFIRST = 1
* Loop over all possible binary numbers
      DO 200 I=1, MX
        IF(IFIRST.EQ.1) THEN
*. Initial number 
          IZERO = 0
          CALL ISETVC(IWORK,IZERO,NOPEN)
          IFIRST = 0
        ELSE
*. Next number
          ADD=1
          J=0
  190     CONTINUE
          J=J+1
          IF(IWORK(J).EQ.1) THEN
            IWORK(J)=0
          ELSE
            IWORK(J)=1
            ADD=0
          END IF
          IF( ADD .EQ. 1 ) GOTO 190
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Next possible det '
          CALL IWRTMA(IWORK,1,NOPEN,1,NOPEN)
        END IF
*. Correct spinprojections in the various gasspaces
        IBORB = 1
        IMOKAY = 1
        DO IGAS = 1, NGAS
          NALPHA = 0
          DO IORB = 1, NOP_GAS(IGAS)
            IF(IWORK(IBORB-1+IORB).EQ.1) NALPHA = NALPHA + 1
          END DO
          NBETA = NOP_GAS(IGAS) - NALPHA
          JMS2 = NALPHA - NBETA
          IF(JMS2.NE.IMS2_GAS(IGAS))IMOKAY = 0
C?        WRITE(6,*) ' NALPHA, NBETA, JMS2, IMOKAY = ',
C?   &                 NALPHA, NBETA, JMS2, IMOKAY
          IBORB = IBORB + NOP_GAS(IGAS)
        END DO
        IF(IMOKAY.EQ.1) THEN
          NDET = NDET + 1
          IF(IFLAG.NE.1) THEN
            CALL ICOPVE(IWORK,IDET((IBDET-1+NDET-1)*NOPEN+1),NOPEN)
          END IF
        END IF
  200 CONTINUE
*
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Number of open shells per GAS : '
        CALL IWRTMA(NOP_GAS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Spin-projection per GAS '
        CALL IWRTMA(IMS2_GAS,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' Number of generated determinants ', NDET
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' Generated dets, 1 => alpha, 0 => beta '
          WRITE(6,*) ' ===================================== '
          CALL IWRTMA(IDET((IBDET-1)*NOPEN+1),NOPEN,NDET,NOPEN,NDET)
        END IF
      END IF
*
      RETURN 
      END
      SUBROUTINE SPIN_ADAPT_CC_OP(
     &           NSPOBEX_TP,ISPOBEX,ISPOBEX_TO_OCCLS,
     &           IBSPOBEX_TO_OCCLS,NSPOBEX_FOR_OCCLS,NOCCLS,S)
*
* Master routine for spinadapting CC amplitudes. Required spin 
* is S
*
* Jeppe Olsen, December 2001
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

*
      INTEGER ISPOBEX(4*NGAS,*)
      INTEGER ISPOBEX_TO_OCCLS(*), IBSPOBEX_TO_OCCLS(*)
      INTEGER NSPOBEX_FOR_OCCLS(*)
      INTEGER KLZ_TYPES(NACTEL)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================== '
        WRITE(6,*) ' SPIN_ADAPT_CC_OP entered '
        WRITE(6,*) ' ======================== '
        WRITE(6,*)
      END IF
*
* Local memory take 1
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'ADDL  ',IDUM,'SPIN_A')
*. Alpha and beta occupation of reference space
      CALL MEMMAN(KLALREF_OC,NGAS,'ADDL  ',1,'ALRFOC')
      CALL MEMMAN(KLBEREF_OC,NGAS,'ADDL  ',1,'BERFOC')
*. And of general sd
      CALL MEMMAN(KLAL_OC,NGAS,'ADDL  ',1,'ALOC  ')
      CALL MEMMAN(KLBE_OC,NGAS,'ADDL  ',1,'BEOC  ')
*. MS2 and MAXOP per space
      CALL MEMMAN(KLMS2,NGAS,'ADDL  ',1,'MS2   ')
      CALL MEMMAN(KLMAXOP,NGAS,'ADDL  ',1,'MAXOP ')
*          GET_REF_ALBE_OCC(IREFSPC,IREF_AL,IREF_BE)
      CALL GET_REF_ALBE_OCC(IREFSPC,int_mb(KLALREF_OC),
     &                      int_mb(KLBEREF_OC))
*
*. Largest number of determinants belonging to a given occls
*
*. NMX_DET is largest number of determinants occuring for a given 
*. configuration occuring in actual CC expansion
      NMX_DET = 0 
*. NMX_DET_DNF is largest number of dets for a given configuration.
*. May be larger than NMX_DET as not all dets are neccessarily included.
      NMX_DET_CNF = 0 
      NMX_OP = 0
*. Loop over occupation classes 
      DO JOCCLS = 1, NOCCLS
C?      WRITE(6,*) ' JOCCLS = ', JOCCLS
        NDET = 0
        JB = IBSPOBEX_TO_OCCLS(JOCCLS)
*. Loop over spinorbitalexcitations belonging to this class
        DO JSPOX = 1,  NSPOBEX_FOR_OCCLS(JOCCLS)
         ISPOX = ISPOBEX_TO_OCCLS(JB-1+JSPOX)
C?       WRITE(6,*) ' Spin-orbital excitation : ', ISPOX 
*. Determinant obtained by applying this operator to reference
*. ( Jeppe : Give me a sign !)
C     EXOCC_STROCC(ICR_OCC,IAN_OCC,ISTR_IN_OCC,
C    &           ISTR_OUT_OCC,NGAS,IZERO_STR)
*. Excited alpha string
         CALL EXOCC_STROCC(
     &        ISPOBEX((1-1)*NGAS+1,ISPOX),
     &        ISPOBEX((3-1)*NGAS+1,ISPOX),
     &        int_mb(KLALREF_OC),int_mb(KLAL_OC),NGAS,IZERO)
*. Excited beta string
         CALL EXOCC_STROCC(
     &        ISPOBEX((2-1)*NGAS+1,ISPOX),
     &        ISPOBEX((4-1)*NGAS+1,ISPOX),
     &        int_mb(KLBEREF_OC),int_mb(KLBE_OC),NGAS,IZERO)
*. MS2 and number of unpaired electrons per GAS space 
*     IALBE_TO_MAXOP(IAL,IBE,IMAXOP,IMS2)
         CALL IALBE_TO_MAXOP(int_mb(KLAL_OC),int_mb(KLBE_OC),
     &        int_mb(KLMAXOP),int_mb(KLMS2) )
*. Number of open orbitals 
         NOP = IELSUM(int_mb(KLMAXOP),NGAS)
         NMX_OP = MAX(NOP,NMX_OP)
*. Number of SD's with this MAXOP and MS2
         IFLAG = 1
C             DET_FOR_SPOC(NOP_GAS,IMS_GAS,NGAS,NDET,IDET,IFLAG)
         CALL DET_FOR_SPOC(int_mb(KLMAXOP),int_mb(KLMS2),NGAS,NDETL,
     &                     IDUM,IDUM,IFLAG)
         NDET = NDET + NDETL
        END DO
*.      ^ End of loop over spinorbital excitations belonging to 
*.        a given class
         NMX_DET = MAX(NMX_DET,NDET)  
      END DO
*
      NOP_AEL = (NMX_OP+MS2)/2
      NMX_DET_CNF = IBION(NMX_OP,NOP_AEL)
      WRITE(6,*) ' Largest number of open orbitals ', NMX_OP
      WRITE(6,*) 
     &' Largest number of OCCURING dets belonging to a given class', 
     &  NMX_DET
      WRITE(6,*) 
     &' Largest number of POSSIBLE dets belonging to a given class', 
     &  NMX_DET_CNF
*
* Memory allocation II
*
*. Space for holding all occuring SD's belonging to given occ class
      L = NMX_DET*NMX_OP
      CALL MEMMAN(KLSD,L,'ADDL  ',1,'SDs   ')
*. Space for holding S2 matrix 
      L = NMX_DET ** 2
      CALL MEMMAN(KLS2,L,'ADDL  ',2,'S2MAT ')
*. and in packed form 
      L = NMX_DET*(NMX_DET+1)/2
      CALL MEMMAN(KLS2P,L,'ADDL  ',2,'S2P   ')
*. Pointers for matrices of Lexical addressing
C     CALL MEMMAN(KLZ,NMX_OP,'ADDL  ',1,'Z_SD  ')
*. Not appropriate as we then would need to access arrays 
*. as work(work(klz))!
*. Lexical address to actual adress 
       L = NMX_DET_CNF
       CALL MEMMAN(KLREO,L,'ADDL   ',1,'REO   ')
*. Matrix holding occupation all dets belonging to a given conf 
       CALL MEMMAN(KLSD_OCC_F,NMX_DET_CNF*NMX_OP,'ADDL  ',1,'SD_OCF')
*. Vector for all dets belonging to given conf 
       CALL MEMMAN(KLSD_F,NMX_DET_CNF,'ADDL  ',2,'SD_F  ')
*.

*
* Loop over occupation classes, generate SD's,  obtain and 
* diagonalize S2 matrix 
*
*. Loop over occupation classes 
      DO JOCCLS = 1, NOCCLS
        IBDET = 1
        WRITE(6,*) ' JOCCLS = ', JOCCLS
        NDET = 0
        JB = IBSPOBEX_TO_OCCLS(JOCCLS)
*. Loop over spinorbitalexcitations belonging to this class
        DO JSPOX = 1,  NSPOBEX_FOR_OCCLS(JOCCLS)
         ISPOX = ISPOBEX_TO_OCCLS(JB-1+JSPOX)
*. Determinant obtained by applying this operator to reference
*. Excited alpha string
         CALL EXOCC_STROCC(
     &        ISPOBEX((1-1)*NGAS+1,ISPOX),
     &        ISPOBEX((3-1)*NGAS+1,ISPOX),
     &        int_mb(KLALREF_OC),int_mb(KLAL_OC),NGAS,IZERO)
*. Excited beta string
         CALL EXOCC_STROCC(
     &        ISPOBEX((2-1)*NGAS+1,ISPOX),
     &        ISPOBEX((4-1)*NGAS+1,ISPOX),
     &        int_mb(KLBEREF_OC),int_mb(KLBE_OC),NGAS,IZERO)
*. MS2 and number of unpaired electrons per GAS space 
*     IALBE_TO_MAXOP(IAL,IBE,IMAXOP,IMS2)
         CALL IALBE_TO_MAXOP(int_mb(KLAL_OC),int_mb(KLBE_OC),
     &        int_mb(KLMAXOP),dbl_mb(KLMS2) )
         NOP = IELSUM(int_mb(KLMAXOP),NGAS)
*
*. SD's of this type
*
         IFLAG = 0
         CALL DET_FOR_SPOC(int_mb(KLMAXOP),dbl_mb(KLMS2),NGAS,NDETL,
     &                     int_mb(KLSD),IBDET,IFLAG)
C?       CALL MEMCHK2('DET_F1')
         IBDET = IBDET + NDETL
         NDET = NDET + NDETL
        END DO
*.      ^ End of loop over spinorbital excitations 
*. Number of possible Dets belonging to this type 
         NOP_AL = (NOP+MS2)/2
         NDET_CNF = IBION(NOP,NOP_AL)
*
*. We have now all the SD's belonging to this occlass.
*. Generate reverse lexical order array and mapping 
*. array Lexical order => actual order
*
*  GEN_W_FOR_PTYPES(MINOP,MAXOP,NEL,MS2,KLZ_FOR_PTYPES)
*. KLZ is obtained in GEN_W_FOR_TYPES
*. Notice KLZ is here obtained as a single number. When we will 
*. study several number of open shells for a given occ class 
*. some changes must be made here !
         CALL GEN_W_FOR_PTYPES(NOP,NOP,NOP,MS2,KLZ_TYPES)
*
*. Collect the number of spinorbitalexcitations included for this excitation
*. type.
*
* Obtain reorder array Lexical to actual number 
*
         CALL REO_SD(int_mb(KLSD),NOP,NOP_AL,NDET,NDET_CNF,
     &              int_mb(KLZ_TYPES(NOP+1)),int_mb(KLREO))
*. Set up S2 matrix for this set of determinants
C     S2_MAT_FOR_DETS(IDET,NDET,NOP,MS2,IZ,IREO,S2MAT,
C    &                           IWORK,XWORK)
         CALL S2_MAT_FOR_DETS(int_mb(KLSD),NDET,NOP,MS2,
     &      int_mb(KLZ_TYPES(NOP+1)),
     &      int_mb(KLREO),dbl_mb(KLS2),int_mb(KLSD_OCC_F),
     &      dbl_mb(KLSD_F))
*. Diagonalize S2 matrix 
C DIAG_SYM_MAT(A,X,SCR,NDIM,ISYM)
         CALL DIAG_SYM_MAT(dbl_mb(KLS2),dbl_mb(KLS2),dbl_mb(KLS2P),
     &                     NDET,0)
         WRITE(6,*) ' Eigenvalues of S2 '
         CALL WRTMAT(dbl_mb(KLS2P),NDET,1,NDET,1)
*. . Find the eigensolutions with spin S
         THRES = 1.0D-10
         DO IEIG = 1, NDET
           IF(ABS(dbl_mb(KLS2P-1+IEIG)-S*(S+1.0D0)).LE.THRES) THEN
             dbl_mb(KLSD_F-1+IEIG) = 1.0D0
           ELSE
             dbl_mb(KLSD_F-1+IEIG) = -1.0D0
           END IF
         END DO
*
         IF(NTEST.GE.100) THEN
           WRITE(6,*) ' Vector defining P/Q partitioning '
           CALL WRTMAT(dbl_mb(KLSD_F),1,NDET,1,NDET)
         END IF
*
      END DO
*     ^ End of loop over occupation classes 
*
      RETURN
      END
      FUNCTION IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
*
* Adress of prototype determinant IAB
* alpha occupation is used to define lex address
*
* Jeppe Olsen, Dec. 2001  
*
      INCLUDE 'implicit.inc'
      INTEGER Z(NOPEN,NALPHA)
      DIMENSION IAB(*),NEWORD(*)
*
      NTEST = 000 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IZNUM_PTDT, NOPEN, NALPHA = ', NOPEN,NALPHA
        WRITE(6,*) ' Input Z- matrix '
        CALL IWRTMA(Z,NOPEN,NALPHA,NOPEN,NALPHA)
C?      IF(IREORD.NE.0) THEN
C?        WRITE(6,*) ' First 20 elements of reorder '
C?        CALL IWRTMA3(NEWORD,1,20,1,20)
C?      END IF
      END IF
*
      IZ = 1
      IALPHA = 0
      DO I = 1,NOPEN
        IF(IAB(I).GT.0) THEN
          IALPHA = IALPHA + 1
          IZ = IZ + Z(I,IALPHA)
        END IF
      END DO
*
C?    WRITE(6,*) ' IZ = ', IZ
      IF(IREORD.EQ.0) THEN
        IZNUM_PTDT = IZ
      ELSE
        IZNUM_PTDT = NEWORD(IZ)
      END IF
*
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from IZNUM_PTDT '
        WRITE(6,*) ' Prototype determinant '
        CALL IWRTMA(IAB,1,NOPEN,1,NOPEN)
        WRITE(6,*) ' IZ and Address = ', IZ, IZNUM_PTDT
      END IF
*
      RETURN
      END
      FUNCTION IZNUM_PTDT2(IAB,NOPEN_DIM,NOPEN,IOPEN1,IALPHA1,IZ1,
     &         NALPHA,Z,NEWORD,IREORD)
*
* Adress of prototype determinant IAB
* alpha occupation is used to define lex address
*
* Jeppe Olsen, May 2013; Modified form to allow calculation of partial address
*
      INCLUDE 'implicit.inc'
      INTEGER Z(NOPEN_DIM,NALPHA)
      DIMENSION IAB(*),NEWORD(*)
*
      NTEST = 000 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IZNUM_PTDT2: NOPEN, NALPHA = ', NOPEN,NALPHA
        WRITE(6,*) ' IOPEN1, IALPHA1, IZ1 = ', IOPEN1, IALPHA1, IZ1
 
        WRITE(6,*) ' Input Z- matrix '
        CALL IWRTMA(Z,NOPEN_DIM,NALPHA,NOPEN_DIM,NALPHA)
C?      IF(IREORD.NE.0) THEN
C?        WRITE(6,*) ' First 20 elements of reorder '
C?        CALL IWRTMA3(NEWORD,1,20,1,20)
C?      END IF
      END IF
*
C     IZ = 1 
      IZ = IZ1
C?    write(6,*) ' IZ, 0 = ', IZ
      IALPHA = IALPHA1
      DO I = IOPEN1,NOPEN
        IF(IAB(I-IOPEN1+1).GT.0) THEN
          IALPHA = IALPHA + 1
          IZ = IZ + Z(I,IALPHA)
C?        WRITE(6,*) ' Updated IZ = ', IZ
        END IF
      END DO
*
      IF(NTEST.GE.100) WRITE(6,*) ' IZ = ', IZ
      IF(IREORD.EQ.0) THEN
        IZNUM_PTDT2 = IZ
      ELSE
        IZNUM_PTDT2 = NEWORD(IZ)
      END IF
*
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from IZNUM_PTDT2 '
        WRITE(6,*) ' Prototype determinant '
        CALL IWRTMA(IAB,1,NOPEN-IOPEN1+1,1,NOPEN-IOPEN1+1)
        WRITE(6,*) ' IZ and Address = ', IZ, IZNUM_PTDT2
      END IF
*
      RETURN
      END
      SUBROUTINE REO_PTDET(NOPEN,NALPHA,IZ_PTDET,IREO_PTDET,ILIST_PTDET,
     &                     NLIST_PTDET,ISCR)
*
* A list(ILIST_PTDET) of prototype determinants with NOPEN unpaired electrons and 
* NALPHA alpha electrons is given. 
*
* Obtain 1 : Z matrix(IZ_PTDET) for this set of prototype dets
*        2 : Reorder array(IREO_PTDET) going from lexical order to 
*            the order specified by the prototype dets 
*            given in ILIST_PTDET. The reordering goes 
*            from lexical order to actual order. 
*            Prototype determinants not included in 
*            ILIST_PTDET are given zero address
*          
*. Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ILIST_PTDET(NOPEN,NLIST_PTDET)
*. Output
      INTEGER IZ_PTDET(NALPHA,NOPEN), IREO_PTDET(*)
*. Local scratch : Min length : 2*NOPEN + (NALPHA+1)*(NOPEN+1)
      INTEGER ISCR(*)
      NTEST = 000
*
* 1 : Set up lexical order array for prototype determinants
*     (alpha considered as occupied electron)
      KLMIN = 1
      KLMAX = KLMIN + NOPEN
      KLW   = KLMAX + NOPEN
*. No consideration of inactive orbitals
      NINOBL = 0
      CALL MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),1,NOPEN,NALPHA,
     &                 NINOBL,NTEST)
*
*. Arc weights
*
      CALL GRAPW(ISCR(KLW),IZ_PTDET,ISCR(KLMIN),ISCR(KLMAX),NOPEN,
     &           NALPHA,NTEST)
*
*. Reorder array
*
*. Total number of prototype determinants
      NTOT_PTDET = 0
      IF (NALPHA .GE. 0 .AND. NALPHA .LE. NOPEN) THEN
         NTOT_PTDET = IBION(NOPEN,NALPHA)
      ELSE
         NTOT_PTDET = 0
      ENDIF
      IZERO = 0
      CALL ISETVC(IREO_PTDET,IZERO,NTOT_PTDET)
*
      DO JPTDT = 1, NLIST_PTDET
C?     WRITE(6,*) ' JPTDT = ', JPTDT
*. Lexical address of prototype determiant JPTDT
       ILEX = IZNUM_PTDT(ILIST_PTDET(1,JPTDT),NOPEN,NALPHA,IZ_PTDET,
     &                   IDUM,0)
C             IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
       IREO_PTDET(ILEX) = JPTDT 
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' REO_PTDET: List of prototype dets '
        CALL IWRTMA3(ILIST_PTDET,NOPEN,NLIST_PTDET,
     &                           NOPEN,NLIST_PTDET)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reorder array for prototype determinants '
        CALL IWRTMA(IREO_PTDET,1,NTOT_PTDET,1,NTOT_PTDET)
      END IF
*
      RETURN 
      END 
      SUBROUTINE NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN,NCOMP_PER_OPEN,
     &                         NCOMP)
*
* Number of configurations per number of open orbitals is given
* Find total number of some components, defined by number
* of components per open
*
* In practice : components are SD's, CSF's or CMB's
*
* Jeppe Olsen, Dec. 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NCONF_PER_OPEN(*), NCOMP_PER_OPEN(*)
*
      NCOMP = 0
      DO IOPEN = 0, MAXOP  
        NCOMP = 
     &  NCOMP + NCONF_PER_OPEN(IOPEN+1)*NCOMP_PER_OPEN(IOPEN+1)
      END DO
*
      RETURN
      END
      SUBROUTINE CSDTMT_GAS(IPDTCNF,IPCSCNF,DTOC,IPRCSF)
*
* Construct in IPDTCNF list of proto type combinations in IDFTP
* Construct in IPCSCNF list of proto type CSF's in ICFTP
* Construct in DTOC matrix expanding proto type CSF's in terms of
* prototype combinations in DTOC
*
* Also:
* KZ_PTDET(IOPEN+1): lexical array for prototype dets
* KREO_PTDET(IOPEN+1): Reorder array from lexical to actual address for protodets
*
* where MXPDBL is size of largest prototype combination block .
*
* Jeppe Olsen
*
* Changed to Combination form, June 92
* Adapted to LUCIA December 2001
* Cleaned, corrected and shaved a bit in June 2011 
* Some further modifications, Dec. 2011
*.
*. Last modification; Oct. 13, 2012; Jeppe Olsen, correcting for lowspin
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      INTEGER IPDTCNF(*),IPCSCNF(*)
      DIMENSION DTOC(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CSDTMT')
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRCSF)
*. 
*. Size of largest scratch block for CSF-SD routine: (NDET+1)*NOPEN
*
      MAX_SCR = 0
      DO IOPEN = 0, MAXOP
        L = (NPCMCNF(IOPEN+1)+1)*IOPEN
        MAX_SCR = MAX(MAX_SCR,L)
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' Size of largest scratch block ',
     &MAX_SCR
      LSCR = MAX_SCR
      CALL MEMMAN(KLSCR1,LSCR,'ADDL  ',2,'SCR_CS')  !done
*
*
* .. Set up combinations and upper determinants
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' **************************************'
        WRITE(6,*) ' Generation of proto type determinants '
        WRITE(6,*) ' **************************************'
        WRITE(6,*)
      END IF
*. Still tired of stupid compiler warnings
      IDTBS = 0
      ICSBS = 0
C     DO IOPEN = 0, MAXOP
      DO IOPEN = MINOP, MAXOP
        ITP = IOPEN + 1 
        IF( NTEST .GE. 5 ) THEN
          WRITE(6,*)
          WRITE(6,'(A,I3,A)')
     &    '       Type with ',IOPEN,' open orbitals '
          WRITE(6,'(A)')
     &    '       **********************************'
          WRITE(6,*)
        END IF
        IF( ITP .EQ. MINOP+1) THEN
          IDTBS = 1
          ICSBS = 1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NPCSCNF(ITP-1)
        END IF
COLD    WRITE(6,*) ' a: IOPEN, IDTBS, ICSBS = ', IOPEN,IDTBS,ICSBS
C
        IF( IOPEN .NE. 0 ) THEN
*. Proto type combinations and branching diagram for
*  proto type combinations
          IF( MS2+1 .EQ. MULTS ) THEN
            IFLAG = 2
            CALL SPNCOM_LUCIA(IOPEN,MS2,NNDET,IPDTCNF(IDTBS),
     &                  IPCSCNF(ICSBS),IFLAG,PSSIGN,IPRCSF)
C                SPNCOM(NOPEN,MS2,NDET,IABDET,
C    &                  IABUPP,IFLAG,PSSIGN,IPRCSF)
          ELSE
            IFLAG = 1
            CALL SPNCOM_LUCIA(IOPEN,MS2,NNDET,IPDTCNF(IDTBS),
     &                  IPCSCNF(ICSBS),IFLAG,PSSIGN,IPRCSF)
            IFLAG = 3
            CALL SPNCOM_LUCIA(IOPEN,MULTS-1,NNDET,IPDTCNF(IDTBS),
     &                  IPCSCNF(ICSBS),IFLAG,PSSIGN,IPRCSF)
           END IF
         END IF
      END DO
*     ^ End of loop over number of open orbitals
*
* Set up z-matrices for addressing prototype determinants with 
* a given number of open orbitals, and for readdressing to 
* the order given in IDCNF
*. Scr : largest block of 2*NOPEN + (NALPHA+1)*(NOPEN+1)
      LSCR = 0
      DO IOPEN = MINOP, MAXOP
        IF(MOD(IOPEN-MS2,2).EQ.0) THEN
          IALPHA = (IOPEN+MS2)/2
          L = 2*IOPEN + (IALPHA+1)*(IOPEN+1)
          LSCR = MAX(L,LSCR)
        END IF
      END DO
      CALL MEMMAN(KLSCR2,LSCR,'ADDL  ',1,'KLSCR2')  !done
*
      IDTBS = 1
      DO IOPEN = MINOP, MAXOP
        ITP = IOPEN + 1
*. Correction, oct. 13, 2012, JO
*. Was
CERR    IF( ITP .EQ. 1) THEN
*.should be
        IF( ITP .EQ. MINOP+1 ) THEN
          IDTBS = 1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
        END IF
COLD    WRITE(6,*) ' b: IOPEN, IDTBS = ', IOPEN,IDTBS
COLD    IB_SD_OPEN(IOPEN+1) = IDTBS

        IALPHA = (IOPEN+MS2)/2
        IDET = NPCMCNF(IOPEN+1)
C?      WRITE(6,*) ' IOPEN, IDET = ', IOPEN, IDET
        CALL REO_PTDET(IOPEN,IALPHA,int_mb(KZ_PTDT(ITP)),
     &                 int_mb(KREO_PTDT(ITP)),IPDTCNF(IDTBS),
     &                 IDET,int_mb(KLSCR2) )
C?     WRITE(6,*) ' Z array as delivered by REO_PTDET'
C?     CALL IWRTMA(WORK(KZ_PTDT(ITP)),IOPEN,IALPHA,IOPEN,IALPHA)
C     REO_PTDET(NOPEN,NALPHA,IZ_PTDET,IREO_PTDET,ILIST_PTDET,
C    &                     NLIST_PTDET,ISCR)
      END DO

*
*. matrix expressing csf's in terms of determinants
*
*
*. Tired of compiler warnings
      IDTBS = 0
      ICSBS = 0
      ICDCBS = 0
C     DO IOPEN = 0, MAXOP
      DO IOPEN = MINOP, MAXOP
        ITP = IOPEN + 1
        IF( ITP .EQ. MINOP+1 ) THEN
          IDTBS = 1
          ICSBS = 1
          ICDCBS =1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NPCSCNF(ITP-1)
          ICDCBS = ICDCBS + NPCMCNF(ITP-1)*NPCSCNF(ITP-1)
        END IF
C       IF(NPCMCNF(ITP)*NPCSCNF(ITP).EQ.0) GOTO 30
        IF( NTEST .GE. 5 ) THEN
          WRITE(6,*)
          WRITE(6,*) ' ************************************'
          WRITE(6,*) ' CSF - SD/COMB transformation matrix '
          WRITE(6,*) ' ************************************'
          WRITE(6,'(A)')
          WRITE(6,'(A,I3,A)')
     &    '  Type with ',IOPEN,' open orbitals '
          WRITE(6,'(A,I10)')
     &    '  Dimension ',ICDCBS
          WRITE(6,'(A)')
     &    '  ************************************'
          WRITE(6,*)
        END IF
        IF(IOPEN .EQ. 0 ) THEN
           DTOC(ICDCBS) = 1.0D0
        ELSE
          CALL CSFDET_LUCIA(IOPEN,IPDTCNF(IDTBS),NPCMCNF(ITP),
     &               IPCSCNF(ICSBS),NPCSCNF(ITP),DTOC(ICDCBS),
     &               dbl_mb(KLSCR1),PSSIGN,IPRCSF)
C              CSFDET(NOPEN,IDET,NDET,ICSF,NCSF,CDC,WORK,PSSIGN,
C    &                IPRCSF)
        END IF
      END DO
*     ^ End of loop over number of open shells
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CSDTMT')
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)  ' List of CSF-SD transformation matrices '
        WRITE(6,*)  ' ======================================='
        WRITE(6,*)
        IDTBS = 1
        ICSBS = 1
        ICDCBS = 1
C       DO IOPEN = 0, MAXOP
        DO IOPEN = MINOP, MAXOP
          ITP = IOPEN + 1
          IF( ITP .EQ. MINOP + 1 ) THEN
            IDTBS = 1
            ICSBS = 1
            ICDCBS =1
          ELSE
            IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
            ICSBS = ICSBS + (IOPEN-1)*NPCSCNF(ITP-1)
            ICDCBS = ICDCBS + NPCMCNF(ITP-1)*NPCSCNF(ITP-1)
          END IF
          NNCS = NPCSCNF(ITP)
          NNCM = NPCMCNF(ITP)
          IF(NNCS.GT.0.AND.NNCM.GT.0) THEN
            WRITE(6,*) ' Number of open shells : ', IOPEN
            WRITE(6,*) ' Number of combinations per conf ', NNCM
            WRITE(6,*) ' Number of CSFs per conf         ', NNCS
            CALL WRTMAT(DTOC(ICDCBS),NNCM,NNCS,NNCM,NNCS)
          END IF
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE CNFORD_GAS(IOCCLS,NOCCLS,IOCCLS_LIST,ISYM,PSSIGN,
     &           IPRCSF,ICONF_OCC,ICONF_REO,ICTSDT,SGNCTS,
     &           IBLTP,IBLOCK,NBLOCK) 
*
*
* Generate determinants in configuration order and obtain
* sign array for switching between the two formats.
*
*
* It is assumed that CSFDIM has been called 
*
* Jeppe Olsen Dec. 2001 from CNFORD
*             Dec. 2011, cleaned and shaved
*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
       INTEGER IOCCLS(NOCCLS)
       INTEGER IOCCLS_LIST(NGAS,*)
       INTEGER IBLOCK(8,NBLOCK), IBLTP(*)
*. Output 
      INTEGER ICONF_OCC(*),ICONF_REO(*)
      DIMENSION ICTSDT(*),SGNCTS(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' Output from CNFORD_GAS '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
      END IF
*
      NELEC = IELSUM(IOCCLS_LIST(1,IOCCLS(1)),NGAS)
C 
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CNFORD')
c..dongxia: the following arrays are not used. comment off
c     CALL MEMMAN(KLZSCR,(NOCOB+1)*(NELEC+1),'ADDL  ',1,'ZSCR  ')
c     CALL MEMMAN(KLZ,NOCOB*NELEC*2,'ADDL  ',1,'Z     ')
c     CALL MEMMAN(KLOCMIN,NOCOB,'ADDL  ',1,'OCMIN ')
c     CALL MEMMAN(KLOCMAX,NOCOB,'ADDL  ',1,'OCMAX ')
c..
*
* Reorder Determinants from configuration order to ab-order
*
C          REO_GASDET(INBLOCK,NBLOCK,ISYM,IREO,SREO )
      CALL REO_GASDET(IBLTP,IBLOCK,NBLOCK,ISYM,ICTSDT,SGNCTS)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CNFORD')
      RETURN
      END
      SUBROUTINE SPNCOM_LUCIA(NOPEN,MS2,NDET,IABDET,
     &                  IABUPP,IFLAG,PSSIGN,IPRCSF)
*
* Combinations of nopen unpaired electrons.Required
* spin projection MS2/2.
*
      INCLUDE 'implicit.inc'
*. MXPDIM is included to have access to MXPORB 
      INCLUDE 'mxpdim.inc'
      INTEGER ADD
      DIMENSION IABDET(NOPEN,*),IABUPP(NOPEN,*)
*. Should have length of max number of open orbitals
      DIMENSION IWORK(MXPORB)
*
* LENGTH OF IWORK MUST BE AT LEAST NOPEN
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRCSF)
      NDET=0
      NUPPER = 0
*
* Determinants are considered as binary numbers,1=alpha,0=beta
*
      MX=2 ** NOPEN
      CALL ISETVC(IWORK,0,NOPEN)
      IFIRST = 1
* Loop over all possible binary numbers
      DO 200 I=1, MX
        IF(IFIRST.EQ.1) THEN
*. Initial number 
          IZERO = 0
          CALL ISETVC(IWORK,IZERO,NOPEN)
          IFIRST = 0
        ELSE
*. Next number 
          ADD=1
          J=0
  190     CONTINUE
          J=J+1
          IF(IWORK(J).EQ.1) THEN
            IWORK(J)=0
          ELSE
            IWORK(J)=1
            ADD=0
          END IF
          IF( ADD .EQ. 1 ) GOTO 190
        END IF
C
C.. 2 :  CORRECT SPIN PROJECTION ?
        NALPHA=0
        DO J=1,NOPEN
          NALPHA=NALPHA+IWORK(J)
        END DO
C
        IF(2*NALPHA-NOPEN.EQ.MS2.AND. 
     &    .NOT.(PSSIGN.NE.0.AND. IWORK(1).EQ.0)) THEN
          IF (IFLAG .LT. 3 ) THEN
            NDET=NDET+1
            CALL ICOPVE(IWORK,IABDET(1,NDET),NOPEN)
          END IF
*
          IF (IFLAG .GT. 1 ) THEN
C UPPER DET ?
            MS2L = 0
            LUPPER = 1
C
            DO IEL = 1,NOPEN
              IF (IWORK(IEL).EQ.1) THEN
                 MS2L = MS2L + 1
              ELSE
                 MS2L = MS2L - 1
              END IF
              IF( MS2L .LT. 0 ) LUPPER = 0
            END DO
*
            IF( LUPPER .EQ. 1 ) THEN
              NUPPER = NUPPER + 1
              CALL ICOPVE(IWORK,IABUPP(1,NUPPER),NOPEN)
            END IF
          END IF
        END  IF
C
  200 CONTINUE
C
      XMSD2=DBLE(MS2)/2
C
      IF(NTEST.GE.5.AND.IFLAG .NE.3) THEN
         WRITE(6,1010) NOPEN,NDET,XMSD2
 1010    FORMAT(1H0,2X,I3,' Unpaired electrons give ',I5,/,
     +'           combinations with spin projection ',F12.7)
         WRITE(6,*)
         WRITE(6,'(A)') '  Combinations : '
         WRITE(6,'(A)') '  ============== '
         WRITE(6,*)
         DO 20 J=1,NDET
           WRITE(6,1020) J,(IABDET(K,J),K=1,NOPEN)
  20     CONTINUE
 1020    FORMAT(1H0,I5,2X,30I2,/,(1H ,7X,30I2))
      END IF
C
      IF( IFLAG.GT.1.AND.NTEST.GE.5) THEN
         WRITE(6,*)
         WRITE(6,'(A)') ' Upper determinants '
         WRITE(6,'(A)') ' ================== '
         WRITE(6,*)
         DO 22 J=1,NUPPER
           WRITE(6,1020) J,(IABUPP(K,J),K=1,NOPEN)
  22     CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE CSFDET_LUCIA(NOPEN,IDET,NDET,ICSF,NCSF,CDC,WORK,PSSIGN,
     &                  IPRCSF)
*
* Expand csf's in terms of combinations with
* the use of the Graebenstetter method ( I.J.Q.C.10,P142(1976) )
*
* Input :
*         NOPEN : NUMBER OF OPEN ORBITALS
*         IDET  : OCCUPATION OF combinations
*         NDET  : NUMBER OF combinations
*         ICSF  : INTERMEDIATE SPIN COUPLINGS OF
*                 CSF'S IN BRANCHING DIAGRAM
* Output :
*         CDC :  NDET X NCSF MATRIX
*                GIVING EXPANSION FROM COMB'S TO CSF,S
*                CSF BASIS = Comb basis *CDC
* Scratch :
*          WORK SHOULD AT LEAST BE (NDET+1)*NOPEN
*
* If combinations are use ( signaled by PSSIGN .ne. 0 )
* the factors are multiplies with sqrt(2), corresponding to 
* a combination being 1/sqrt(2) times the sum or difference of two
* determinants
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IDET(NOPEN,NDET),ICSF(NOPEN,NCSF)
      DIMENSION CDC(NDET,NCSF)
      DIMENSION WORK(*)

      NTEST = 0
      NTEST = MAX(IPRCSF,NTEST)
      IF(PSSIGN.EQ.0.0D0) THEN
       CMBFAC = 1.0D0
      ELSE 
       CMBFAC = SQRT(2.0D0)
      END IF
C
      KLFREE = 1
      KLMDET = KLFREE
      KLFREE = KLMDET + NDET * NOPEN
      KLSCSF =  KLFREE
      KLFREE = KLSCSF + NOPEN
C
C.. OBTAIN INTERMEDIATE VALUES OF MS FOR ALL DETERMINANTS
      DO 10 JDET = 1, NDET
        CALL MSSTRN_LUCIA(IDET(1,JDET),WORK(KLMDET+(JDET-1)*NOPEN),
     &        NOPEN,IPRCSF)
   10 CONTINUE
C
      DO 1000 JCSF = 1, NCSF
       IF( NTEST .GE. 105 ) WRITE(6,*) ' ....Output for CSF ',JCSF
C
C OBTAIN INTERMEDIATE COUPLINGS FOR CSF
      CALL MSSTRN_LUCIA(ICSF(1,JCSF),WORK(KLSCSF),NOPEN,IPRCSF)
C
      DO 900 JDET = 1, NDET
C EXPANSION COEFFICIENT OF DETERMINANT JDET FOR CSF JCSF
      COEF = 1.0D0
      SIGN = 1.0D0
      JDADD = (JDET-1)*NOPEN
      DO 700 IOPEN = 1, NOPEN
C
C + + CASE
        IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.1) THEN
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.1.AND.IDET(IOPEN,JDET).EQ.0) THEN
C + - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) )
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN) )
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.1) THEN
C - + CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)-WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
          SIGN  = - SIGN
        ELSE IF(ICSF(IOPEN,JCSF).EQ.0.AND.IDET(IOPEN,JDET).EQ.0) THEN
C - - CASE
          COEF = COEF *
     &    (WORK(KLSCSF-1+IOPEN)+WORK(KLMDET-1+JDADD+IOPEN) +1.0D0)
     &    /(2.0D0*WORK(KLSCSF-1+IOPEN)+2.0D0 )
        END IF
  700 CONTINUE
       CDC(JDET,JCSF) = SIGN * CMBFAC * SQRT(COEF)
  900 CONTINUE
 1000 CONTINUE
C
      IF( NTEST .GE. 5) THEN
        WRITE(6,*)
        WRITE(6,'(A,2I2)')
     &  '  The CDC array for  NOPEN ',NOPEN
        WRITE(6,*) ' NDET, NCSF = ', NDET,NCSF
        WRITE(6,*)
        CALL WRTMAT(CDC,NDET,NCSF,NDET,NCSF)
       END IF
C
      RETURN
      END
      SUBROUTINE MSSTRN_LUCIA(INSTRN,UTSTRN,NOPEN,IPRCSF)
C
C A STRING IS GIVEN IN FORM A SEQUENCE OF ZEROES
C AND ONE ' S
C
C REINTERPRET THIS AS :
C
C 1 : THE INPUT STRING IS A DETERMINANT AND THE
C     1'S INDICATE ALPHA ELECTRONS AND THE
C     0'S INDICATE BETA ELECTRONS .
C     UTSTRN IS THE MS-VALUES ATE EACH VERTEX
C
C 2 : THE INPUT STRING IS A CSF GIVEN IN A
C     BRANCHING DIAGRAM, WHERE
C     1'S INDICATE UPWARDS SPIN COUPLEING
C     WHILE THE 0'S INDICATES DOWNWARDS SPIN COUPLING ,
C     REEXPRESS THIS AS S VALUES OF ALL COUPLINGS
C
C THE TWO PROCEDURES ARE IDENTICAL .
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION INSTRN(NOPEN),UTSTRN(NOPEN)
C
      UTSTRN(1) = DBLE(INSTRN(1)) - 0.5D0
      DO 10 IOPEN = 2, NOPEN
        UTSTRN(IOPEN) = UTSTRN(IOPEN-1) +DBLE(INSTRN(IOPEN))-0.5D0
   10 CONTINUE
C
      NTEST = 0
      NTEST = MAX(NTEST,IPRCSF)
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ... Output from MSSTRN '
         WRITE(6,*) ' INSTRN AND UTSTRN'
         CALL IWRTMA(INSTRN,1,NOPEN,1,NOPEN)
         CALL WRTMAT(UTSTRN,1,NOPEN,1,NOPEN)
       END IF
C
      RETURN
      END
      SUBROUTINE CONF_GRAPH(IOCC_MIN,IOCC_MAX,NORB,NEL,IARCW,NCONF,
     &                      ISCR)
*
* A group of configurations is described by the 
* accumulated min and max, IOCC_MIN and IOCC_MAX.
*
* Find arcweights of corresponding graph and total number 
* of configurations ( all symmetries) 
*
* Jeppe Olsen, Oct. 2001
*
*. NORB should be the number of active orbitals
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC_MIN(NORB), IOCC_MAX(NORB)
*. Output
      INTEGER IARCW(NORB,NEL,2)
* IARCW(I,J,K) gives weight of arc ending at vertex (I,J) 
* with occupation K (=1,2)
*. Local scratch : Length should be (NORB+1)*(NEL+1)
      INTEGER ISCR(NORB+1,NEL+1)
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from CONF_GRAPH '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' IOCMIN and IOCMAX for orbitals '
        CALL IWRTMA(IOCC_MIN,1,NORB,1,NORB)
        CALL IWRTMA(IOCC_MAX,1,NORB,1,NORB)
      END IF
*. Set up vertex weights
      CALL CONF_VERTEX_W(IOCC_MIN,IOCC_MAX,NORB,NEL,ISCR)
      NCONF = ISCR(NORB+1,NEL+1)
*. Obtain arcweights from vertex weights
C?    WRITE(6,*) ' CONF_GRAPH, NORB, NEL = ', NORB, NEL
      CALL CONF_ARC_W(IOCC_MIN,IOCC_MAX,NORB,NEL,ISCR,IARCW)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Arcweights for single occupied arcs '
        CALL IWRTMA(IARCW(1,1,1),NORB,NEL,NORB,NEL)
        WRITE(6,*) ' Arcweights for double occupied arcs '
        CALL IWRTMA(IARCW(1,1,2),NORB,NEL,NORB,NEL)
        WRITE(6,*) ' Total number of configurations ', NCONF
      END IF
*
      RETURN
      END
*
      SUBROUTINE CONF_VERTEX_W(IOCC_MIN,IOCC_MAX,NORB,NEL,IVERTEXW)
*
* Obtain vertex weights for configuration graph
*
* Jeppe Olsen, October 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC_MIN(NORB),IOCC_MAX(NORB)
*. Output
      INTEGER IVERTEXW(NORB+1,NEL+1)
*
      IZERO = 0
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CONF_VERTEX speaking '
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' NORB, NEL = ', NORB, NEL
        WRITE(6,*) ' IOCC_MIN, IOCC_MAX: '
        CALL IWRTMA(IOCC_MIN,1,NORB,1,NORB)
        CALL IWRTMA(IOCC_MAX,1,NORB,1,NORB)
      END IF
     
      CALL ISETVC(IVERTEXW,IZERO,(NORB+1)*(NEL+1)) 
*
      IVERTEXW(0+1,0+1) = 1
      DO IORB = 1, NORB
        DO IEL  = IOCC_MIN(IORB), IOCC_MAX(IORB)
*
          IF(IEL.EQ. 0 ) 
     &    IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)
*
          IF(IEL.EQ. 1 ) 
     &    IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)
     &                           + IVERTEXW(IORB-1+1,IEL+1-1) 
*
          IF(IEL.GE. 2 ) 
     &    IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)
     &                           + IVERTEXW(IORB-1+1,IEL+1-1) 
     &                           + IVERTEXW(IORB-1+1,IEL+1-2) 
*
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' Vertex weights as an (NORB+1)*(NEL+1) matrix '
        CALL IWRTMA(IVERTEXW,NORB+1,NEL+1,NORB+1,NEL+1)
      END IF
*
      RETURN
      END
      SUBROUTINE CONF_ARC_W(IOCC_MIN,IOCC_MAX,NORB,NEL,IVERTEXW,IARCW)
*
* Obtain arcweights for single and double occupied arcs 
* from vertex weights
*
* Jeppe Olsen, October 2001
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IVERTEXW(NORB+1,NEL+1)
      INTEGER IOCC_MIN(NORB),IOCC_MAX(NORB)
*. Output
      INTEGER IARCW(NORB,NEL,2)
      IZERO = 0
      CALL ISETVC(IARCW,IZERO,2*NORB*NEL)
* IARCW(I,J,K) is weight of arc with occupation K ending at (I,J) 
* IARCW(I,J,K) = Sum(J-K < L <= J)   IVERTEXW(I-1,L)
      DO I = 1, NORB
       DO J = 1, NEL
        IF(IOCC_MIN(I).LE.J .AND. J.LE.IOCC_MAX(I)) THEN
          DO K = 1, NEL
            IF(K.EQ.1) IARCW(I,J,K) = IVERTEXW(I-1+1,J+1)
            IF(K.EQ.2.AND.J.GE.2) IARCW(I,J,K) = IVERTEXW(I-1+1,J+1)
     &                                         + IVERTEXW(I-1+1,J-1+1)
          END DO
        END IF
       END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Arc weights for single occupied arcs '
        CALL IWRTMA(IARCW(1,1,1),NORB,NEL,NORB,NEL)
        WRITE(6,*) ' Arc weights for double occupied arcs '
        CALL IWRTMA(IARCW(1,1,2),NORB,NEL,NORB,NEL)
      END IF
*
      RETURN
      END
      SUBROUTINE EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_OPEN_MS,NEL)
*
* A determinant IDET_OC, IDET_MS is given. Extract spinprojections 
* for open shells 
*
* Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
*.input
      INTEGER IDET_OC(NEL),IDET_MS(NEL)
*. Output
      INTEGER IDET_OPEN_MS(*)
*
      IEL = 1
      IOPEN = 0
*. Loop over electrons
 1000 CONTINUE
       IF(IEL.LT.NEL) THEN
         IF(IDET_OC(IEL).NE.IDET_OC(IEL+1)) THEN
*. Single occupied orbital so
            IOPEN = IOPEN + 1
            IDET_OPEN_MS(IOPEN) = IDET_MS(IEL) 
            IEL = IEL + 1
          ELSE 
            IEL = IEL + 2
          END IF
       ELSE 
*. Last electron was not identical to previous, so 
*. neccessarily single occupied.
          IOPEN = IOPEN + 1
          IDET_OPEN_MS(IOPEN) = IDET_MS(IEL)
          IEL = IEL + 1
       END IF
      IF(IEL.LE.NEL) GOTO 1000
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' EXTRT_MS_OP... Input det, occ and ms '
        CALL IWRTMA(IDET_OC,1,NEL,1,NEL)
        CALL IWRTMA(IDET_MS,1,NEL,1,NEL)
        WRITE(6,*) ' Number of open orbitals = ', IOPEN
        WRITE(6,*) ' Output det : ms of open orbitals '
        CALL IWRTMA(IDET_OPEN_MS,1,IOPEN,1,IOPEN)
      END IF
*
      RETURN
      END
      SUBROUTINE CMP_IVEC_ILIST(IVEC,ILIST,LLIST,NLIST,INUM)
* An integer IVEC of LLIST entries are given.
* compare with list of vectors in ILIST and find first
* vector in LLIST that is identical to IVEC.
*
* If INUM = 0, the list was not found
*
*  Jeppe Olsen, December 2001
*
      INCLUDE 'implicit.inc'
*. General input
      INTEGER ILIST(LLIST,NLIST)
*. Specific input
      INTEGER IVEC(LLIST)
*
      INUM = 0
      DO JLIST = 1, NLIST
        IFOUND = 1
        DO IELMNT = 1, LLIST
          IF(IVEC(IELMNT).NE.ILIST(IELMNT,JLIST))IFOUND = 0
        END DO
        IF(IFOUND.EQ.1) INUM = JLIST
        IF(INUM.NE.0) GOTO 1001
      END DO
*
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input list : '
        CALL IWRTMA(IVEC,1,LLIST,1,LLIST)
        WRITE(6,*) ' Address of list : ', INUM
      END IF
*
      RETURN
      END
      SUBROUTINE IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
*
* Find address of the occupation class corresponding to given types 
* of alpha and beta strings 
*
* Jeppe Olsen, December 2001
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Local scratch
      INTEGER IABOCC(MXPNGAS)
*
      IATP_ABS = IATP + IBSPGPFTP(IAGRP) - 1
      IBTP_ABS = IBTP + IBSPGPFTP(IBGRP) - 1
C?    WRITE(6,*) ' IATP, IBTP, IAGRP, IBGRP = ',
C?   &             IATP, IBTP, IAGRP, IBGRP 
C?    WRITE(6,*) ' IATP_ABS, IBTP_ABS ', IATP_ABS, IBTP_ABS
*
C  IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
      IONE = 1
      CALL IVCSUM(IABOCC,NELFSPGP(1,IATP_ABS),NELFSPGP(1,IBTP_ABS),
     &            IONE,IONE,NGAS)
*. And the address of this occupation class 
      CALL CMP_IVEC_ILIST(IABOCC,int_mb(KIOCCLS),NGAS,NMXOCCLS,INUM)
*
      IOC = INUM
*
      IF(INUM.EQ.0) THEN
        WRITE(6,*) 
     &  ' Combination of alpha and beta string not found as occ-class'
        WRITE(6,*) ' Occ of alpha, Occ of beta, Occ of alpha+beta '
        CALL IWRTMA(NELFSPGP(1,IATP_ABS),1,NGAS,1,NGAS)
        CALL IWRTMA(NELFSPGP(1,IBTP_ABS),1,NGAS,1,NGAS)
        CALL IWRTMA(IABOCC,1,NGAS,1,NGAS)
        STOP 
     &  ' Combination of alpha and beta string not found as occ-class'
      END IF
*
      RETURN 
      END
      SUBROUTINE REO_GASDET(IBLTP,IBLOCK,NBLOCK,ISYM,IREO,SREO)
*
* Create reorder array for determinants : configuration order => Ab order                            order
*
*
* Jeppe Olsen, November 2001, from GASANA
*
*
c      INCLUDE 'implicit.inc'
*
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
      INCLUDE 'spinfo.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      DIMENSION IBLOCK(8,NBLOCK),IBLTP(*)
*
*. Output 
      INTEGER IREO(*), SREO(*)
      CALL LUCIAQENTER('REOGA')
      CALL MEMMAN(KLOFF,DUMMY,'MARK  ',DUMMY,'REO_GA')
*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
* Type of alpha and beta strings
      IATP = 1             
      IBTP = 2              
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL= NAEL+NBEL
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)

*
*
**. Info on block structure of space
*
*
*Space for alpha and beta strings
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')  !done
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')  !done
*. Space for constructing arc weights for configurations 
      CALL MEMMAN(KLZSCR,(NOCOB+1)*(NEL+1),'ADDL  ',1,'ZSCR  ') !done
      CALL MEMMAN(KLZ,NOCOB*NEL*2,'ADDL  ',1,'Z     ')  !done
      CALL MEMMAN(KLOCMIN,NOCOB,'ADDL  ',1,'OCMIN ')  !done
      CALL MEMMAN(KLOCMAX,NOCOB,'ADDL  ',1,'OCMAX ')  !done
*. Occupation and projections of a given determinant
      CALL MEMMAN(KLDET_OC,NAEL+NBEL,'ADDL  ',1,'CONF_O') !done
      CALL MEMMAN(KLDET_MS,NAEL+NBEL,'ADDL  ',1,'CONF_M') !done
      CALL MEMMAN(KLDET_VC,NOCOB,'ADDL  ',1,'CONF_M')  !done

*. Configurations are constructed and stored with orbitals
*. labels starting with 1, whereas strings are stored
*. with orbitals starting with NINOB+1. Provide offset  
      IB_ORB = NINOB+1
      CALL REO_GASDET_S(IREO,int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &            NOCTPA,NOCTPB,IOCTPA,IOCTPB,
     &            NBLOCK,IBLOCK,
     &            NAEL,NBEL,int_mb(KLASTR),int_mb(KLBSTR),IBLTP,NSMST,
     &            NELFSPGP,NMXOCCLS,NGAS,      
     &            int_mb(KIOCCLS),NTOOB,NOBPT,
     &            int_mb(KZCONF),int_mb(KDFTP),
     &            int_mb(KIB_OCCLS(ISYM)),IB_CN_OPEN, 
     &            int_mb(KICONF_REO(1)),IB_CM_OPEN,
     &            int_mb(KLZSCR),int_mb(KLZ),int_mb(KLOCMIN),
     &            int_mb(KLOCMAX),int_mb(KLDET_OC),int_mb(KLDET_MS),
     &            int_mb(KLDET_VC),
     &            KZ_PTDT,KREO_PTDT,MINOP,IBCONF_ALL_SYM_FOR_OCCLS,
     &            IB_ORB,NACOB,NOCOB,PSSIGN,NPCMCNF)
C IB_CN_OPEN
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'REO_GA')
      CALL LUCIAQEXIT('REOGA')
*
      RETURN
      END
      SUBROUTINE REO_GASDET_S(IREO,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                 IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                 NAEL,NBEL,
     &                 IASTR,IBSTR,IBLTP,NSMST,
     &                 NELFSPGP,NOCCLS,NGAS,
     &                 IOCCLS,NORB,NOBPT,
     &                 IZCONF,DFTP,
     &                 IB_CONF_OCCLS,IB_CONF_OPEN, ICONF_REO,
     &                 IB_CM_OPEN,
     &                 IZSCR,IZ,IOCMIN,IOCMAX,IDET_OC,IDET_MS,
     &                 IDET_VC,KZ_PTDT,KREO_PTDT,MINOP,
     &                 IBCONF_ALL_SYM_FOR_OCCLS,IB_ORB,
     &                 NACOB,NOCOB,PSSIGN,NPCMCNF)
*
* Reorder determinants in GAS space from det to configuration order
* Reorder array created is Conf-order => AB-order 
*
      include 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)  
      DIMENSION IBLTP(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
      INTEGER NOBPT(*)
      INTEGER DFTP(*) 
      INTEGER NPCMCNF(*)
*. IB_CONF_OPEN, IB_CM_OPEN(IOPEN+1) Gives start of confs/CM's
*. with given symmetry and number of open orbitals
      INTEGER IB_CONF_OPEN(*), IB_CONF_OCCLS(*), ICONF_REO(*)
      INTEGER IB_CM_OPEN(*)
*. Offset to start of configurations of given occls in list containing all symmetries
      INTEGER IBCONF_ALL_SYM_FOR_OCCLS(NOCCLS)
*. WORK(KZ_PTDT(IOPEN+1) gives Z  array for prototype dets with IOPEN 
*. WORK(KREO_PTDT(IOPEN+1) gives the corresponding reorder array
*. open orbitals
      INTEGER KZ_PTDT(*), KREO_PTDT(*)
*. The work array containing used for WORK(KZ_PTDET()),WORK(KREO_PTDT())
c     DIMENSION WORK(*)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
*. Scratch space 
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      INTEGER IZSCR(*),IZ(*),IOCMIN(*),IOCMAX(*)
      INTEGER IDET_OC(*), IDET_MS(*) , IDET_VC(*)
*. Output
      INTEGER IREO(*)
C     DIMENSION SREO(*)
*
      NTEST = 000
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' output from REO_GASDET_S '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' PSSIGN = ', PSSIGN
      END IF
      IF(NTEST.GE.2000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' IB_CONF_OCCLS(1) = ', IB_CONF_OCCLS(1)
        WRITE(6,*) ' IB_CONF_OPEN(1) = ', IB_CONF_OPEN(1)
        WRITE(6,*) ' IB_CM_OPEN(1) = ', IB_CM_OPEN(1)
        WRITE(6,*) ' IB_CONF_OPEN(3) = ', IB_CONF_OPEN(3)
        WRITE(6,*) ' IB_CM_OPEN(3) = ', IB_CM_OPEN(3)
        WRITE(6,*) ' ICONF_REO(1) = ', ICONF_REO(1)
        WRITE(6,*) ' MINOP = ', MINOP
        WRITE(6,*) ' IB_ORB = ', IB_ORB
      END IF
*
      IAGRP = 1
      IBGRP = 2
*
      NEL = NAEL + NBEL
*
      IDET = 0
      DO JBLOCK = 1, NBLOCK
        CALL MEMCHK2('CHECKA')
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
        IF(NTEST.GE.10000.and.ga_nodeid().eq.0) THEN
          WRITE(6,'(A,4(2X,I6))') ' IATP, IBTP, IASM, IBSM = ',
     &                              IATP, IBTP, IASM, IBSM
        END IF
*. Occupation class of this combination of string 
C            IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
        CALL IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
        IF(NTEST.GE.10000.and.ga_nodeid().eq.0) WRITE(6,'(A,I6)')
     &  ' Corresponding occupation class number = ', IOC
*. Arcweights for this occupation class
        CALL MXMNOC_OCCLS(IOCMIN,IOCMAX,NGAS,NOBPT,IOCCLS(1,IOC),
     &       MINOP,NTEST)
C     MXMNOC_OCCLS(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,NTESTG)
*. the arcweights
         CALL CONF_GRAPH(IOCMIN,IOCMAX,NACOB,NEL,IZ,NCONF_P,IZSCR)
C        CONF_GRAPH(IOCC_MIN,IOCC_MAX,NORB,NEL,IARCW,NCONF,ISCR)
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                           NORB,0,IDUM,IDUM)
*. Offset to this occupation class in occupation class ordered cnf list
        IB_OCCLS = IB_CONF_OCCLS(IOC)
C       IB_OCCLS = IBCONF_ALL_SYM_FOR_OCCLS(IOC)
C?      WRITE(6,*) ' IOC, IB_OCCLS = ', IOC, IB_OCCLS
*. Info for this occupation class :
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        DO  IB = 1,NIB
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB 
            IRESTR2 = 1
          ELSE
            MINIA = 1
            IRESTR2 = 0
          END IF
          DO  IA = MINIA,NIA
            IDET = IDET + 1
C                ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,IDET_OC,IDET_SP,ISIGN)
            CALL ABSTR_TO_ORDSTR(IASTR(1,IA),IBSTR(1,IB),NAEL,NBEL,
     &                           IDET_OC,IDET_MS,ISIGN)
*. Orbital numbers in strings are absolute, but relative in confs.
*. subtract offset(number of inactive orbitals)
            ISUB = 1-IB_ORB
            CALL IADCONST(IDET_OC,ISUB,NAEL+NBEL)
C     IADCONST(IVEC,IADD, NDIM)
*. Number of open orbitals in this configuration 
            NOPEN = NOP_FOR_CONF(IDET_OC,NEL)
            IF(NOPEN.GE.MINOP) THEN
C                    NOP_FOR_CONF(ICONF,NEL)
             NDOUBLE = (NEL-NOPEN)/2
             NOCOB_L = NOPEN + NDOUBLE
             NOPEN_AL = NAEL - NDOUBLE
C?           WRITE(6,*) ' NOPEN, NOPEN_AL = ', NOPEN,NOPEN_AL
             NPTDT = NPCMCNF(NOPEN+1)
*. Packed form of this configuration 
C                 REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY) 
             CALL REFORM_CONF_OCC(IDET_OC,IDET_VC,NEL,NOCOB_L,1)
*. Address of this configuration 
*. Offset to configurations with this number of open orbitals in 
*. reordered cnf list
             IF(NTEST.GE.10000.and.ga_nodeid().eq.0) THEN
               WRITE(6,*) ' IASTR, IBSTR, ICONF:'
               CALL IWRTMA(IASTR(1,IA),1,NAEL,1,NAEL)
               CALL IWRTMA(IBSTR(1,IB),1,NBEL,1,NBEL)
               CALL IWRTMA(IDET_VC,1,NOCOB_L,1,NOCOB_L)
             END IF
*
             ICNF_OUT = ILEX_FOR_CONF2(IDET_VC,NOCOB_L,NACOB,NEL,IZ,1,
     &                  ICONF_REO(1),IB_OCCLS)
C                      ILEX_FOR_CONF2(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,
C                      IREO,IB_OCCLS)
C?           WRITE(6,*) ' number of configuration in output list',
C?   &       ICNF_OUT
*. Spinprojections of open orbitals
             CALL EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_VC,NEL)
C                 EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_OPEN_MS,NEL)
*. Address of this spinprojection pattern   
C  IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
             IPTDT = IZNUM_PTDT(IDET_VC,NOPEN,NOPEN_AL,
     &               int_mb(KZ_PTDT(NOPEN+1)),
     &               int_mb(KREO_PTDT(NOPEN+1)),1)
             ISIGNP = 1
             IF(IPTDT.EQ.0) THEN
C?            IF(IRESTR2 .EQ. 1) THEN
              IF(PSSIGN.NE.0) THEN
*. The determinant was not found among the list of prototype dets. For combinations
*. this should be due to the prototype determinant is the MS- switched determinant, so find
*. address of this and remember sign
                M1 = -1
                CALL ABSTR_TO_ORDSTR(IBSTR(1,IB),IASTR(1,IA),NBEL,NAEL,
     &                            IDET_OC,IDET_MS,ISIGN)
*. Spinprojections of open orbitals
                CALL EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_VC,NEL)
                IPTDT = IZNUM_PTDT(IDET_VC,NOPEN,NOPEN_AL,
     &                  int_mb(KZ_PTDT(NOPEN+1)),
     &                  int_mb(KREO_PTDT(NOPEN+1)),1)
                IF(PSSIGN.EQ.-1.0D0) ISIGNP = -1
              ELSE 
*. Prototype determinant was not found in list
               WRITE(6,*) 
     &         ' Error: Determinant not found in list of protodets'
               WRITE(6,*) ' Detected in REO_GASDET_S '
               WRITE(6,*) ' AB pattern of det '
               CALL IWRTMA3(IDET_VC,1,NOPEN,1,NOPEN)
              END IF
             END IF
*
             IBCNF_OUT = IB_CONF_OPEN(NOPEN+1)
             IF(NTEST.GE.10000.and.ga_nodeid().eq.0) THEN
              WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
              WRITE(6,*) ' IB_CM_OPEN(NOPEN+1) = ',
     &                     IB_CM_OPEN(NOPEN+1)
              WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
              WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
             END IF
             IADR_SD_CONF_ORDER = IB_CM_OPEN(NOPEN+1) - 1
     &                          + (ICNF_OUT-IBCNF_OUT)*NPTDT + IPTDT
             IF(IADR_SD_CONF_ORDER.LE.0.and.ga_nodeid().eq.0) THEN
               WRITE(6,*) ' Problemo, IADR_SD_CONF_ORDER < 0 '
               WRITE(6,*) ' IADR_SD_CONF_ORDER = ', IADR_SD_CONF_ORDER
               WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
               WRITE(6,*) ' IB_CM_OPEN(NOPEN+1) = ',
     &                     IB_CM_OPEN(NOPEN+1)
               WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
               WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
C?             CALL XFLUSH(6)
             END IF
             IF(NTEST.GE.10000.and.ga_nodeid().eq.0) THEN
               WRITE(6,*) ' IADR_SD_CONF_ORDER, ISIGN, IDET = ',
     &                      IADR_SD_CONF_ORDER, ISIGN, IDET
             END IF
             IREO(IADR_SD_CONF_ORDER) = ISIGN*IDET*ISIGNP
             IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
               WRITE(6,*) ' IDET, IADR_SD_CONF_ORDER ',
     &                      IDET, IADR_SD_CONF_ORDER
             END IF
            END IF! Nopen .ge. MINOP
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
      END DO
*     ^ End of loop over blocks
*
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Reorder array for SDs, CONF order => string order'
        WRITE(6,*) ' ================================================='
        CALL IWRTMA(IREO,1,IDET,1,IDET)
      END IF
*
*. Check sum of reordering array
      I_DO_CHECKSUM = 0
      IF(I_DO_CHECKSUM.EQ.1) THEN
        ISUM = 0
        DO JDET = 1, IDET
          ISUM = ISUM + ABS(IREO(JDET))
        END DO
        IF(ISUM.NE.IDET*(IDET+1)/2) THEN
          if (ga_nodeid().eq.0) then
          WRITE(6,*) ' Problem with sumcheck in REO_GASDET'
          WRITE(6,'(A,2I9)') 
     &    'Expected and actual value ', ISUM, IDET*(IDET+1)/2
          endif
          STOP       ' Problem with sumcheck in REO_GASDET'
        ELSE
          if (ga_nodeid().eq.0) 
     &       WRITE(6,*) ' Sumcheck in REO_GASDET passed '
        END IF
      END IF !checksum is invoked
*
      RETURN
      END
      SUBROUTINE WRITE_TO_LEA_FORMAT(NOP,XPORQ,COEFS,IOPEN,NOPEN,
     &           IREF_AL,IREF_BE,NOPEN_PER_GAS,IOC_PER_GAS,LU_LEA)
*
* Write results of Spin analysis of spinorbital excitations to 
* disc in form readable by Lea's program
*
* Jeppe + Lea, Feb 19, 2002
*
* NOP : Number of operators/determinants
* XPORQ : IPORQ(IOP) = 1 => Operator is P-operator
*         IPORQ(IOP) =-1 => Operator is Q-operator
* COEFS : Coefficients for the various operators
* IOPEN  :Spinprojections of the open orbitals 
* NOPEN : Number of open(single occupied) orbitals in this 
*         configuration.
* NOPEN_PER_GAS : Number of open orbitals per GAS space
* IOC_PER_GAS : Number of electrons per GAS in actual configurations
* ISPOBEX : The spin-orbital excitations 
* ISPOBEX_FOR_IOPEN : The spinorbital excitations for the 
*                     various IOPEN patterns
*
*
*Loop over Occupationsklasser
* Skriv - med havelaager foran lidt info om occ class
* Loop over P eller Q
*   Loop over P operatorer for hver okkupationsklasse 
*     Skriv : # P/Q ( for at indikere at naeste operator er en P/Q op)
*     Skriv %L "name"   (Som A1 I2 ( "name" laeses ikke ind))
*     Skriv name, f.eks 3a+(v 1x)a+(v 2x) .....  
*     ( foerste tal er nummer for denne op som er I2)
*     Loop over Spinorbital excitationer for denne P operator
*       Skriv : * 2  tstring ( hvor 2 (I2) er antallet af crea/anni-ops  i denne 
*                            og string er for vores skyld)
*        Skriv coefficient
*        a+(v2a)a (o1a)            (tal er pt I2) (hver op a3a1I1A1a1) ingen mell
*        Skriv * 2 tstring
*        Skriv koefficient
*        a+(v1b)a (o1b) 
*     End of loop over spinorbital excitations
*    End of loop over operatorer
*   End of loop over P eller Q
* End of loop over Occupationsclasser
*       
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Input 
      DIMENSION XPORQ(NOP)
      DIMENSION COEF(NOP,NOP)
      INTEGER IOPEN(NOPEN,NOP)
      INTEGER IREF_AL(NGAS), IREF_BE(NGAS)
      INTEGER NOPEN_PER_GAS(NGAS)
      INTEGER IOC_PER_GAS(NGAS)
C     INTEGER ISPOBEX(4*NGAS,*)
C     INTEGER ISPOBEX_FOR_IOPEN(NOP)
*. Local scratch 
      CHARACTER*1 CDENT
*. (assumes atmost 100-fold excitations)
      CHARACTER*8 CDEL_OC(100)
      CHARACTER*8 CSPOBEX(100)
*
      INTEGER IREF_OC(MXPNGAS), IDEL_OC(MXPNGAS)
      INTEGER ICA(MXPORB),ICB(MXPORB),IAA(MXPORB),IAB(MXPORB)
*. Occupation of reference state 
      IONE = 1
      CALL IVCSUM(IREF_OC,IREF_AL, IREF_BE, IONE, IONE, NGAS)
C          IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
*
      WRITE(LU_LEA,'(1A)') '#'
      WRITE(LU_LEA,'(1A, 15I3)') '#', (IOC_PER_GAS(IGAS),IGAS=1,NGAS)
      WRITE(LU_LEA,'(1A)') '#'
*. Changes of orbital occupations, reference => Excited conf
      IONEM = -1
      CALL IVCSUM(IDEL_OC, IOC_PER_GAS, IREF_OC, IONE, IONEM,NGAS)
      NEXOP = 0
      DO IGAS = 1, NGAS
       NEXOP = NEXOP + ABS(IDEL_OC(IGAS))
      END DO
      WRITE(6,*) ' Number of changes of orbital occupations ', NEXOP
*. Character form of orbital excitation 
      CALL CHAR_FORM_OF_ORBOP(IDEL_OC,CDEL_OC,ILEN_ORBOP)
*
      
      INUM = 1
      DO IPQ = 1, 2
        IF(IPQ.EQ.1) THEN
          XMARK = 1.0D0
          CDENT = 'P'
        ELSE 
          XMARK = -1.0D0
          CDENT = 'Q'
        END IF
    
        WRITE(LU_LEA,'(A,A)') '# ', CDENT
*. Loop over P/Q operators
        LOP = 0
        DO KOP = 1, NOP
         IF(XPORQ(KOP).EQ.XMARK) THEN
           WRITE(LU_LEA,'(A)') '% "name"'
*. Number of indeces in name 
           WRITE(LU_LEA,*) ILEN_ORBOP
*. Number of this operator 
           LOP = LOP + 1
           WRITE(LU_LEA,*) LOP
*. And the orbital operator
           WRITE(LU_LEA,'(A)') (CDEL_OC(IOP),IOP=1,ILEN_ORBOP)
*. Loop over spin-orbital excitations 
           DO LSPOBEX = 1, NOP
*. The corresponding spinorbital excitation
             CALL SPOBEXTP_ORDER(IOPEN(1,LSPOBEX),NOPEN,
     &            NOPEN_PER_GAS,NOBPT,
     &            ICA,NCA,ICB,NCB,IAA,NAA,IAB,NAB)
C     SPOBEXTP_ORDER(IOPEN,NOPEN,NOPEN_PER_GAS,
C    &           NOB_PER_GAS,
C    &           ICA,NCA,ICB,NCB,IAA,NAA,IAB,NAB)
*. Change into character form 
C     CHAR_FOR_ELOP(ICA,IAB,IHPV,INUM,CHAROP)
           END DO
         END IF
        END DO
      END DO
*
      RETURN
      END 
      SUBROUTINE CHAR_FORM_OF_ORBOP(IDEL_OC,CDEL_OC,NELM_OP)
*
* Character form of orbital excitation operator
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input 
      INTEGER IDEL_OC(NGAS)
*. Output 
      CHARACTER*8 CDEL_OC(*)
*. Local scratch
      CHARACTER*1 CHPVGAS(3)
      DATA CHPVGAS/'I','V','T'/
*. Creation, virtual 
      NV_CREA = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.2) NV_CREA = NV_CREA + IDEL_OC(NGAS)
      END DO
*. Creation, Valence 
      NT_CREA = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3.AND.IDEL_OC(IGAS).GT.0) 
     &  NT_CREA = NT_CREA + IDEL_OC(NGAS)
      END DO
*. Annihilation, valence 
      NT_ANNI = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.3.AND.IDEL_OC(IGAS).LT.0) 
     &  NT_ANNI = NT_ANNI + ABS(IDEL_OC(NGAS))
      END DO
*. Annihilation, inactive 
      NI_ANNI = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) NI_ANNI = NI_ANNI + ABS(IDEL_OC(NGAS))
      END DO
      NELM_OP = NV_CREA+NT_CREA+NI_ANNI+NT_ANNI
*
      WRITE(6,*) ' NV_CREA, NT_CREA, NI_ANNI, NT_ANNI = ',
     &             NV_CREA, NT_CREA, NI_ANNI, NT_ANNI  
*
      IOP = 1
*. Creation of Valence 
      DO I = 1, NT_CREA
        CALL CHAR_FOR_ELOP(1,0,3,I,CDEL_OC(IOP))
        IOP = IOP + 1
      END DO
*. Creation of Virtual 
      DO I = 1, NV_CREA
        CALL CHAR_FOR_ELOP(1,0,2,I,CDEL_OC(IOP))
        IOP = IOP + 1
      END DO
*. Annihilation of valence 
      DO I = 1, NT_ANNI
        CALL CHAR_FOR_ELOP(2,0,3,I,CDEL_OC(IOP))
        IOP = IOP + 1
      END DO
*. Annihilation of inactive 
      DO I = 1, NI_ANNI
        CALL CHAR_FOR_ELOP(2,0,1,I,CDEL_OC(IOP))
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Form of orbital excitation operator (IDEL_OC) '
        CALL IWRTMA(IDEL_OC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of elementary operators : ', NELM_OP
        WRITE(6,*) ' The operator in character form '
        DO I = 1, NELM_OP
          WRITE(6,'(A)') CDEL_OC(I) 
        END DO
      END IF
*
      RETURN
      END 
      SUBROUTINE CHAR_FOR_ELOP(ICA,IAB,IHPV,INUM,CHAROP)
*
* Generate character representation of elementary operator
*
* Jeppe Olsen, Feb. 2002
*
      INCLUDE 'implicit.inc'
*. Character for various operators  in format for Leas program
* Inactive : I
* Secondary : V
* Valence  : T
*. Output 
      CHARACTER*8 CHAROP
      CHARACTER*1 CORBTYP(3)
      DATA  CORBTYP/'O','V','T'/
      CHARACTER*1 CHARNUM(10)
      DATA CHARNUM/'0','1','2','3','4','5','6','7','8','9'/
*
      CHAROP(1:1) = 'a'
      CHAROP(3:3) = '('
      CHAROP(8:8) = ')'
      IF(ICA.EQ.1) THEN
        CHAROP(2:2) = '+'
      ELSE 
        CHAROP(2:2) = ' '
      END IF
      CHAROP(4:4) = CORBTYP(IHPV)
*. Character rep of number (less or equal to 
      CALL INTEGER_TO_CHAR2(INUM,CHAROP(5:5))
      IF(IAB.EQ.1) THEN
        CHAROP(7:7) = 'A'
      ELSE IF(IAB.EQ.2.OR.IAB.EQ.-1) THEN
        CHAROP(7:7) = 'B'
      ELSE IF (IAB.EQ.0) THEN
         CHAROP(7:7) ='X'
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ICA, IAB,  INUM for operator ',
     &               ICA, IAB,  INUM
        WRITE(6,'(A,A)') ' Character for operator : ', CHAROP
      END IF
*
      RETURN
      END 
      SUBROUTINE INTEGER_TO_CHAR2(I,C)
*
*. Convert integer less than 100 to character*2
*
      INCLUDE 'implicit.inc'
      INTEGER I
      CHARACTER*2 C
      CHARACTER*1 CHARNUM(10)
      DATA CHARNUM/'0','1','2','3','4','5','6','7','8','9'/
*
      IONES = MOD(I,10)
      C(2:2) = CHARNUM(IONES)
      ITEENS = (I-IONES)/10
      IF(ITEENS.EQ.0) THEN
        C(1:1) = ' '
      ELSE
        C(1:1) = CHARNUM(ITEENS+1)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,I10,A)') ' Character and integer', I, C
      END IF
*
      RETURN
      END
      SUBROUTINE SPOBEXTP_ORDER(IOPEN,NOPEN,NOPEN_PER_GAS,
     &           NOB_PER_GAS,
     &           ICA,NCA,ICB,NCB,IAA,NAA,IAB,NAB)
*
* A pattern of open  orbitals (IOPEN) is given
*.Find the indeces of the spinorbitals in the spinorbitalexcitations
*
* IT IS PT ASSUMED THAT REFERENCE IS HIGH SPIN OPEN REFERENCE STATE,
* IE OPEN ORBITALS ARE ALPHA ORBITALS
*
* Jeppe Olsen, Feb. 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
*. Input
      INTEGER IOPEN(NOPEN), NOPEN_PER_GAS(NGAS)
      INTEGER NOB_PER_GAS(NGAS)
*. Local scratch
      INTEGER IA(MXPORB), IB(MXPORB)
*. Output
      INTEGER ICA(*), ICB(*), IAA(*), IAB(*)
*
CE    IB = 1
      NCA_TOT = 0
      NCB_TOT = 0
      NAA_TOT = 0
      NAB_TOT = 0
      DO IGAS = 1, NGAS 
        IF(IGAS.EQ.1) THEN
          IIB = 1
        ELSE 
          IIB = IIB + NOPEN_PER_GAS(IGAS-1) 
        END IF
*. Obtain open open orbitals for this gasspace 
        NOP = NOPEN_PER_GAS(IGAS)
C              AB_FOR_OPENST(IA,NA,IB,IB,IAB,NAB,ISIGN)
        CALL AB_FOR_OPENST(IA,NA,IB,NB,IOPEN(IIB),NOP,ISIGN)
*   
        IF(IHPVGAS(IGAS).EQ.2) THEN
*. Particle space, Copy created orbitals to CA and CB
          DO JCA = 1, NA
            ICA(NCA_TOT-1+JCA) = IA(JCA)
          END DO
          NCA_TOT = NCA_TOT + NA
          DO JCB = 1, NB
            ICB(NCB_TOT-1+JCB) = IB(JCB)
          END DO
          NCB_TOT = NCB_TOT + NB
        ELSE IF (IHPVGAS(IGAS).EQ.1) THEN
*. Hole space, annihilated alpha corresponds to open beta and vice versa
          DO JAB = 1, NB
            IAA(NAA_TOT-1+JAB) = IB(JAB)
          END DO
          NAA_TOT = NAA_TOT + NB
          DO JAA = 1, NA
            IAB(NAA_TOT-1+JAA) = IA(JAA)
          END DO
          NAB_TOT = NAB_TOT + NA
        ELSE IF (IHPVGAS(IGAS).EQ.3) THEN
*. Valence space. IT IS HERE THAT THE ASSUMPTION OF HIGH SPIN OS IS USED
*. Open beta corresponds to created beta 
          DO JCB = 1, NB
            ICB(NCB-1+JCB) = IB(JCB)
          END DO
          NCB_TOT = NCB_TOT + NB
*.  Lacking open alpha corresponds to annihilated alpha
          II = 1
          DO JORB = 1, NOB_PER_GAS(IGAS)
            IF(IA(II).NE.JORB) THEN
              IAA(NAA_TOT)=JORB
              II = II + 1
              NAA_TOT = NAA_TOT
              IF(II.GT.NA) GOTO 101
            END IF
          END DO
  101     CONTINUE
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from SPOBEXTP_ORDER'
        WRITE(6,*) ' MS projection of open orbitals '
        CALL IWRTMA(IOPEN,1,NOPEN,1,NOPEN)
        WRITE(6,*) ' Number of open orbital per GAS space '
        CALL IWRTMA(NOPEN_PER_GAS,1,NGAS,1,NGAS)
        WRITE(6,*) ' orbital indeces of ICA '
        CALL IWRTMA(ICA,1,NCA_TOT,1,NCA_TOT)
        WRITE(6,*) ' orbital indeces of ICB '
        CALL IWRTMA(ICB,1,NCB_TOT,1,NCB_TOT)
        WRITE(6,*) ' orbital indeces of IAA '
        CALL IWRTMA(IAA,1,NAA_TOT,1,NAA_TOT)
        WRITE(6,*) ' orbital indeces of IAB '
        CALL IWRTMA(IAB,1,NAB_TOT,1,NAB_TOT)
      END IF
*
      RETURN
      END
      SUBROUTINE AB_FOR_OPENST(IA,NA,IB,NB,IAB,NAB,ISIGN)
*
* A string of open orbitals IAB is given. Extract      
* the alpha and beta orbitals
*
* Jeppe Olsen, March 2002
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IAB(NAB)
*. Output
      INTEGER IA(*), IB(*)
*
      NA = 0
      NB = 0
      DO I = 1, NAB
        IF(IAB(I).EQ.1) THEN
          NA = NA + 1
          IA(NA) = I
        ELSE
          NB = NB + 1
          IB(NB) = I
        END IF
      END DO
*
      NTEST = 00 
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' AB string : '
        CALL IWRTMA(IAB,1,NAB,1,NAB)
        WRITE(6,*) ' Index of alphaspinorbitals'
        CALL IWRTMA(IA,1,NA,1,NA)
        WRITE(6,*) ' Index of betaspinorbitals'
        CALL IWRTMA(IB,1,NB,1,NB)
      END IF
*
      RETURN
      END
C     SUBROUTINE CHAR_FORM_OF_SPOBEX(ICA,NCA,ICB,ICB,IAA,NAA,IAB,NAB,
C                CHAR_SPOBEX,IHPV)
*
* A spinorbital excitation is given in the form of 
* orbital numbers for the various parts. 
*
* Obtain character form of spinorbital excitation 
*
* Jeppe Olsen, April 1, 2002
*
C     INCLUDE 'implicit.inc'
C     INTEGER IHPV(*)
*. Input
C     INTEGER ICA(NCA),ICB(NCB),IAA(NAA),IAB(NAB)
*. Output
C     CHARACTER*8 CHAR_SPOBEX(100)
*. Creation of alpha
C     CHAR_FOR_ELOP(ICA,IAB,IHPV,INUM,CHAROP)
C     IBOP = 0
C     DO JCA = 1, NCA
C       IBOP = IBOP + 1
C       CALL CHAR_FOR_ELOP(1,1,IHPVGAS(
C      
C
      SUBROUTINE GET_NORB_ITV(NI,NV,NA)
*. Obtain number of Inactive, Valence and Secondary orbitals
*. for case with only one inactive space (besides frozen)
*. and one secondary space (besides deleted)
*
* Inactive spaces are assumed to be the first hole spaces
* deleted spaces are assumed to be the last particle spaces
*
* Jeppe Olsen, April 2002
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*
      IP = 0
      IH = 0
      IV = 0
      DO IGAS = 1, NGAS
        IF(IHPVGAS(IGAS).EQ.1) IH = IGAS  
        IF(IHPVGAS(IGAS).EQ.2.AND.IP.EQ.0) IP = IGAS
        IF(IHPVGAS(IGAS).EQ.3) IV = IGAS
      END DO
*
      NI = NOBPT(IH)
      NA = NOBPT(IP)
      NV = NOBPT(IV)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Number of inactive orbitals ( beside frozen) ', NI
        WRITE(6,*) 
     &  ' Number of valence orbitals                  ', NV
        WRITE(6,*) 
     &  ' Number of secondary orbitals ( beside deleted) ', NA
      END IF
*
      RETURN
      END 
      SUBROUTINE FI_HS(H,FI_AL,FI_BE,ECC,IDOH2)
*
* Obtain core energy and alpha and beta part of inactive Fock-matrix 
* / effective 1- electron operator for high spin openshell 
* case
*
* Jeppe Olsen, July 3, 2002
*
* The formulae goes 
*
* FI_ALPHA(I,J) = H(I,J) 
*               + SUM(K) (RHO_ALPHA(K,K) + RHO_BETA(K,K)) G(IJKK)
*               - SUM(K)  RHO_ALPHA(K,K) * G(IKKJ)
*
* ECORE = sum(i) ((RHO_ALPHA(I,I) + RHO_BETA (I,I)) H(II)
*       + 1/2 sum(i,j) (RHO_ALPHA(I,I)*RHO_ALPHA(J,J)) ( G(IIJJ) - G(IJJI) )
*       + 1/2 sum(i,j) (RHO_BETA (I,I)*RHO_BETA (J,J)) ( G(IIJJ) - G(IJJI) )
*       +     sum(i,j) RHO_ALPHA(I,I)*RHO_BETA(J,J)*G(IIJJ)
*
* The alpha and beta parts of the density matrices are not accessed, 
* instead is IHPVGAS_AB used to obtain the occupied spaces
*
      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
*. Standard 1 electron integrals 
cnw   DIMENSION H(*)
      integer H
*. Output 
cnw   DIMENSION FI_AL(*), FI_BE(*)
      integer FI_AL, FI_BE
*
* Is matrix symmetrically packed or ??? 
      ISYMPACK = 1
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' FI_HS : Initial 1 e integral list '
cnw     CALL APRBLM2(H,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
        call ga_print(h)

        WRITE(6,*) ' FI_HS: IHPVGAS_AB and ITPFSO '
        CALL IWRTMA(IHPVGAS_AB,NGAS,2,MXPNGAS,2)
        CALL IWRTMA(ITPFSO,1,NTOOB,1,NTOOB)
      END IF
*
* Core energy 
*
      ECC = 0.0D0
      IJSM = 1
      IIOFF = 0

*. One-electron part 
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IF(ISYMPACK.EQ.0) THEN
            IIOFF = IIOFF + NTOOBS(ISM-1)** 2 
          ELSE
            IIOFF = IIOFF + NTOOBS(ISM-1)*(NTOOBS(ISM-1)+1)/2
          END IF
        END IF
        LEN = NTOOBS(ISM)
        DO I = 1, LEN                           
          IF(ISYMPACK.EQ.0) THEN
            II = IIOFF -1 + (I-1)*LEN + I 
          ELSE 
            II = IIOFF -1 + I*(I+1)/2
          END IF
          FACTOR = 0.0D0
          IGAS = ITPFSO(I+IBSO(ISM)-1)
          IF(IHPVGAS_AB(IGAS,1).EQ.1) FACTOR = 1.0D0
          IF(IHPVGAS_AB(IGAS,2).EQ.1) FACTOR = FACTOR + 1.0D0
          call ga_get(h,ii,ii,1,1,hval,1)
          ECC = ECC + FACTOR*hval
cnw       ECC = ECC + FACTOR*H(II)
        END DO
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' one-electron part to ECC ', ECC
*. Two-electron part
      IF(IDOH2.NE.0) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IBSO(ISM), IBSO(ISM) + NTOOBS(ISM)-1
          DO J = IBSO(JSM), IBSO(JSM) + NTOOBS(JSM)-1
              IP = IREOST(I)
              JP = IREOST(J)
              IGAS = ITPFSO(I)
              JGAS = ITPFSO(J)
*. Factor for Coulomb and exchange
              FACTOR_C = 0.0D0
              FACTOR_E = 0.0D0
              IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,1).EQ.1) 
     &        THEN
               FACTOR_C = 0.5D0
               FACTOR_E = -0.5D0
              END IF
              IF(IHPVGAS_AB(IGAS,2).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1)
     &        THEN
               FACTOR_C = FACTOR_C + 0.5D0
               FACTOR_E = FACTOR_E -0.5D0
              END IF
              IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1) 
     &        THEN
               FACTOR_C = FACTOR_C + 1.0D0
              END IF
              IF(FACTOR_C.NE.0.0D0) 
     &        ECC = ECC + FACTOR_C*GTIJKL(IP,IP,JP,JP)
              IF(FACTOR_E.NE.0.0D0)
     &        ECC = ECC + FACTOR_E*GTIJKL(IP,JP,JP,IP)
          END DO
          END DO
        END DO
        END DO
      END IF
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
      END IF
*
*.  Inactive Fock matrix
*
* Original one-electron matrix
      IF(IDOH2.NE.0) THEN 
*. Number of 1-e integrals
C  NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
        NINT = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
cnw     CALL COPVEC(H,FI_AL,NINT)
cnw     CALL COPVEC(H,FI_BE,NINT)
        call ga_copy(h,fi_al)
        call ga_copy(h,fi_be)
        IJSM = 1
*. ( works for lower half packed only if operator is total sym) 
        IJ = 0
        DO ISM = 1, NSMOB
          CALL SYMCOM(2,6,ISM,JSM,IJSM)
          IJOFF = IJ + 1
          IF(JSM.NE.0) THEN
            DO J = IBSO(JSM),IBSO(JSM) + NTOOBS(JSM) - 1
              IF(ISYMPACK.EQ.0) THEN
                IMIN = IBSO(ISM)
              ELSE
                IMIN = J
              END IF
              DO I = IMIN,IBSO(ISM) + NTOOBS(ISM) - 1
                IP = IREOST(I)
                JP = IREOST(J)
                IF(ISYMPACK.EQ.0) THEN
                  IJ= IJ + 1
                ELSE 
                  IREL = I - IBSO(ISM) + 1
                  JREL = J - IBSO(JSM) + 1
                  IJ = IJOFF - 1 + IREL*(IREL-1)/2 + JREL
                END IF
                al_val = 0.0d0
                be_val = 0.0d0
                DO KSYM = 1, NSMOB
                  DO K = IBSO(KSYM),IBSO(KSYM)-1+NTOOBS(KSYM)
                    KP = IREOST(K)
                    KGAS = ITPFSO(K)
                    FACTOR_C    = 0.0D0
                    FACTOR_E_AL = 0.0D0
                    FACTOR_E_BE = 0.0D0
                    IF(IHPVGAS_AB(KGAS,1).EQ.1) THEN
                      FACTOR_C = 1.0D0
                      FACTOR_E_AL = -1.0D0
                    END IF
                    IF(IHPVGAS_AB(KGAS,2).EQ.1) THEN
                      FACTOR_C = FACTOR_C + 1.0D0
                      FACTOR_E_BE = -1.0D0
                    END IF
                    IF(FACTOR_C.NE.0.0D0) THEN 
                     al_val = al_val + FACTOR_C*GTIJKL(IP,JP,KP,KP)
                     be_val = be_val + FACTOR_C*GTIJKL(IP,JP,KP,KP)
cnw                  FI_AL(IJ) = 
cnw  &               FI_AL(IJ) + FACTOR_C*GTIJKL(IP,JP,KP,KP)
cnw                  FI_BE(IJ) = 
cnw  &               FI_BE(IJ) + FACTOR_C*GTIJKL(IP,JP,KP,KP)
                    END IF
                    IF(FACTOR_E_AL.NE.0.0D0) THEN 
                     al_val = al_val + FACTOR_E_AL*GTIJKL(IP,KP,KP,JP)
cnw                  FI_AL(IJ) = 
cnw  &               FI_AL(IJ) + FACTOR_E_AL*GTIJKL(IP,KP,KP,JP)
                    END IF
                    IF(FACTOR_E_BE.NE.0.0D0) THEN
                     be_val = be_val + FACTOR_E_BE*GTIJKL(IP,KP,KP,JP)
cnw                  FI_BE(IJ) = 
cnw  &               FI_BE(IJ) + FACTOR_E_BE*GTIJKL(IP,KP,KP,JP)
                    END IF
                  END DO
                END DO
                if (ga_nodeid().eq.0) then
                call ga_acc(fi_al,ij,ij,1,1,al_val,1,1.0d0)
                call ga_acc(fi_be,ij,ij,1,1,be_val,1,1.0d0)
                endif
              END DO
            END DO
          END IF
        END DO
      END IF
      call ga_sync()
*
      IF(NTEST.NE.0.AND.IDOH2.NE.0) THEN
       WRITE(6,*) ' FI_alpha in Symmetry blocked form '
       WRITE(6,*) ' =================================='
       WRITE(6,*) 
       ISYM = 0
cnw    CALL APRBLM2(FI_AL,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
       call ga_print(fi_al)
       WRITE(6,*) ' FI_beta in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       ISYM = 0
cnw    CALL APRBLM2(FI_BE,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
       call ga_print(fi_be)
      END IF
*
      RETURN
      END
      SUBROUTINE FI_HS_AB(HA,HB,FI_AL,FI_BE,ECC,IDOH2)
*
* Obtain core energy and alpha and beta part of inactive Fock-matrix 
* / effective 1- electron operator for high spin openshell 
* for case where integrals are similarity transformed integrals
* and the alpha and beta spin orbitals differs
*
* Jeppe Olsen, July 5, 2002
*
* The formulae goes 
*
* FI_ALPHA(I,J) = HA(I,J) 
*               + SUM(K) RHO_ALPHA(K,K)*(G_AA(IJKK) - G_AA(IJKK))
*               + SUM(K) RHO_BETA(K,K))*GAB(IJKK)
*
* FI_BETA(I,J) = HB(I,J) 
*               + SUM(K) RHO_BETA(K,K)*(G_BB(IJKK) - G_BB(IJKK))
*               + SUM(K) RHO_ALPHA(K,K))*GAB(KKIJ)
*
* ECORE = sum(i) ((RHO_ALPHA(I,I) H_ALPHA(II) + RHO_BETA (I,I)) H(II)
*       + 1/2 sum(i,j) (RHO_ALPHA(I,I)*RHO_ALPHA(J,J)) ( GAA(IIJJ) - G(IJJI) )
*       + 1/2 sum(i,j) (RHO_BETA (I,I)*RHO_BETA (J,J)) ( GBB(IIJJ) - G(IJJI) )
*       +     sum(i,j) RHO_ALPHA(I,I)*RHO_BETA(J,J)*GAB(IIJJ)
*
* The alpha and beta parts of the density matrices are not accessed, 
* instead is IHPVGAS_AB used to obtain the occupied spaces
*
      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'oper.inc'
*. Unrestricted 1e integrals
      DIMENSION HA(*),HB(*)
*. Output 
      DIMENSION FI_AL(*), FI_BE(*)
*
* Matrices are lower triangular matrices so 
      ISYMPACK = 1
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' FI_HS_AB : Initial alpha 1 e integral list '
        CALL APRBLM2(HA,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
        WRITE(6,*) ' FI_HS_AB : Initial beta  1 e integral list '
        CALL APRBLM2(HB,NTOOBS,NTOOBS,NSMOB,ISYMPACK)

        WRITE(6,*) ' FI_HS_AB: IHPVGAS_AB and ITPFSO '
        CALL IWRTMA(IHPVGAS_AB,NGAS,2,MXPNGAS,2)
        CALL IWRTMA(ITPFSO,1,NTOOB,1,NTOOB)
*
      END IF
*
* Core energy 
*
      ECC = 0.0D0
      IJSM = 1
      IIOFF = 0

*. One-electron part 
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IF(ISYMPACK.EQ.0) THEN
            IIOFF = IIOFF + NTOOBS(ISM-1)** 2 
          ELSE
            IIOFF = IIOFF + NTOOBS(ISM-1)*(NTOOBS(ISM-1)+1)/2
          END IF
        END IF
        LEN = NTOOBS(ISM)
        DO I = 1, LEN                           
          IF(ISYMPACK.EQ.0) THEN
            II = IIOFF -1 + (I-1)*LEN + I 
          ELSE 
            II = IIOFF -1 + I*(I+1)/2
          END IF
          IGAS = ITPFSO(I+IBSO(ISM)-1)
          IF(IHPVGAS_AB(IGAS,1).EQ.1) ECC = ECC + HA(II)
          IF(IHPVGAS_AB(IGAS,2).EQ.1) ECC = ECC + HB(II)
        END DO
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' one-electron part to ECC ', ECC
*. Two-electron part
      IF(IDOH2.NE.0) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IBSO(ISM), IBSO(ISM) + NTOOBS(ISM)-1
          DO J = IBSO(JSM), IBSO(JSM) + NTOOBS(JSM)-1
              IP = IREOST(I)
              JP = IREOST(J)
              IGAS = ITPFSO(I)
              JGAS = ITPFSO(J)
              IF(IP.NE.JP) THEN
                IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,1).EQ.1) 
     &          THEN
                  ISPCAS = 1
                  ECC = ECC + 
     &               0.5D0*(  GTIJKL(IP,IP,JP,JP) 
     &                      - GTIJKL(IP,JP,JP,IP) )
                END IF
                IF(IHPVGAS_AB(IGAS,2).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1)
     &          THEN
                  ISPCAS = 2
                  ECC = ECC + 
     &               0.5D0*(  GTIJKL(IP,IP,JP,JP) 
     &                      - GTIJKL(IP,JP,JP,IP) )
                END IF
              END IF
              IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1) 
     &        THEN
                ISPCAS = 3
                ECC = ECC + GTIJKL(IP,IP,JP,JP)
              END IF
          END DO
          END DO
        END DO
        END DO
      END IF
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
      END IF
*
*.  Inactive Fock matrix
*
* Original one-electron matrix
      IF(IDOH2.NE.0) THEN 
*. Number of 1-e integrals
C  NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
        NINT = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
        CALL COPVEC(HA,FI_AL,NINT)
        CALL COPVEC(HB,FI_BE,NINT)
        IJSM = 1
        IJ = 0
        DO ISM = 1, NSMOB
          CALL SYMCOM(2,6,ISM,JSM,IJSM)
          IJOFF = IJ + 1
          IF(JSM.NE.0) THEN
            DO J = IBSO(JSM),IBSO(JSM) + NTOOBS(JSM) - 1
              IF(ISYMPACK.EQ.0) THEN
                IMIN = IBSO(ISM)
              ELSE
                IMIN = J
              END IF
              DO I = IMIN,IBSO(ISM) + NTOOBS(ISM) - 1
                IP = IREOST(I)
                JP = IREOST(J)
                IF(ISYMPACK.EQ.0) THEN
                  IJ= IJ + 1
                ELSE 
                  IREL = I - IBSO(ISM) + 1
                  JREL = J - IBSO(JSM) + 1
                  IJ = IJOFF - 1 + IREL*(IREL-1)/2 + JREL
                END IF
                DO KSYM = 1, NSMOB
                  DO K = IBSO(KSYM),IBSO(KSYM)-1+NTOOBS(KSYM)
                    KP = IREOST(K)
                    KGAS = ITPFSO(K)
                    IF(IHPVGAS_AB(KGAS,1).EQ.1) THEN
                      ISPCAS = 1
                      FI_AL(IJ) = FI_AL(IJ)  
     &                          + GTIJKL(IP,JP,KP,KP)
     &                          - GTIJKL(IP,KP,KP,JP)
                      ISPCAS = 3
                      FI_BE(IJ) = FI_BE(IJ) 
     &                          + GTIJKL(KP,KP,IP,JP)
                    END IF
                    IF(IHPVGAS_AB(KGAS,2).EQ.1) THEN
                      ISPCAS = 2
                      FI_BE(IJ) = FI_BE(IJ)  
     &                          + GTIJKL(IP,JP,KP,KP)
     &                          - GTIJKL(IP,KP,KP,JP)
                      ISPCAS = 3
                      FI_AL(IJ) = FI_AL(IJ) 
     &                          + GTIJKL(IP,JP,KP,KP)
                    END IF
                  END DO
*                 ^ End of loop over K
                END DO
*               ^ End of Loop over KSYM
              END DO
            END DO
          END IF
        END DO
      END IF
*
      IF(NTEST.NE.0.AND.IDOH2.NE.0) THEN
       WRITE(6,*) ' FI_alpha in Symmetry blocked form '
       WRITE(6,*) ' =================================='
       WRITE(6,*) 
       ISYM = 0
       CALL APRBLM2(FI_AL,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
       WRITE(6,*) ' FI_beta in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       ISYM = 0
       CALL APRBLM2(FI_BE,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
      END IF
*
      RETURN
      END
      SUBROUTINE FI_HS_SM_AB(HA,HB,FI_AL,FI_BE,ECC,IDOH2)
*
* Obtain core energy and alpha and beta part of inactive Fock-matrix 
* / effective 1- electron operator for high spin openshell 
* for case where integrals are similarity transformed integrals
* and the alpha and beta spin orbitals differs
*
* Jeppe Olsen, July 5, 2002
*
* The formulae goes 
*
* FI_ALPHA(I,J) = HA(I,J) 
*               + SUM(K) RHO_ALPHA(K,K)*(G_AA(IJKK) - G_AA(IJKK))
*               + SUM(K) RHO_BETA(K,K))*GAB(IJKK)
*
* FI_BETA(I,J) = HB(I,J) 
*               + SUM(K) RHO_BETA(K,K)*(G_BB(IJKK) - G_BB(IJKK))
*               + SUM(K) RHO_ALPHA(K,K))*GAB(KKIJ)
*
* ECORE = sum(i) ((RHO_ALPHA(I,I) H_ALPHA(II) + RHO_BETA (I,I)) H(II)
*       + 1/2 sum(i,j) (RHO_ALPHA(I,I)*RHO_ALPHA(J,J)) ( GAA(IIJJ) - G(IJJI) )
*       + 1/2 sum(i,j) (RHO_BETA (I,I)*RHO_BETA (J,J)) ( GBB(IIJJ) - G(IJJI) )
*       +     sum(i,j) RHO_ALPHA(I,I)*RHO_BETA(J,J)*GAB(IIJJ)
*
* The alpha and beta parts of the density matrices are not accessed, 
* instead is IHPVGAS_AB used to obtain the occupied spaces
*
      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
*. Similarity transformed 1e integrals
      DIMENSION HA(*),HB(*)
*. Output 
      DIMENSION FI_AL(*), FI_BE(*)
*
* Matrices are complete matrices so 
      ISYMPACK = 0
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' FI_HS_AB : Initial alpha 1 e integral list '
        CALL APRBLM2(HA,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
        WRITE(6,*) ' FI_HS_AB : Initial beta  1 e integral list '
        CALL APRBLM2(HB,NTOOBS,NTOOBS,NSMOB,ISYMPACK)

        WRITE(6,*) ' FI_HS_AB: IHPVGAS_AB and ITPFSO '
        CALL IWRTMA(IHPVGAS_AB,NGAS,2,MXPNGAS,2)
        CALL IWRTMA(ITPFSO,1,NTOOB,1,NTOOB)
*
      END IF
*
* Core energy 
*
      ECC = 0.0D0
      IJSM = 1
      IIOFF = 0

*. One-electron part 
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IIOFF = 1
        ELSE 
          IF(ISYMPACK.EQ.0) THEN
            IIOFF = IIOFF + NTOOBS(ISM-1)** 2 
          ELSE
            IIOFF = IIOFF + NTOOBS(ISM-1)*(NTOOBS(ISM-1)+1)/2
          END IF
        END IF
        LEN = NTOOBS(ISM)
        DO I = 1, LEN                           
          IF(ISYMPACK.EQ.0) THEN
            II = IIOFF -1 + (I-1)*LEN + I 
          ELSE 
            II = IIOFF -1 + I*(I+1)/2
          END IF
          IGAS = ITPFSO(I+IBSO(ISM)-1)
          IF(IHPVGAS_AB(IGAS,1).EQ.1) ECC = ECC + HA(II)
          IF(IHPVGAS_AB(IGAS,2).EQ.1) ECC = ECC + HB(II)
        END DO
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' one-electron part to ECC ', ECC
*. Two-electron part
      IF(IDOH2.NE.0) THEN
        DO ISM = 1, NSMOB
        DO JSM = 1, NSMOB
          DO I = IBSO(ISM), IBSO(ISM) + NTOOBS(ISM)-1
          DO J = IBSO(JSM), IBSO(JSM) + NTOOBS(JSM)-1
              IP = IREOST(I)
              JP = IREOST(J)
              IGAS = ITPFSO(I)
              JGAS = ITPFSO(J)
              IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,1).EQ.1) 
     &        THEN
               ECC = ECC + 
     &         0.5D0*(  GTIJKL_SM_AB(IP,IP,JP,JP,4,0) 
     &                - GTIJKL_SM_AB(IP,JP,JP,IP,4,0) )
C                       GTIJKL_SM_AB(I,J,K,L,NA,NB)
              END IF
              IF(IHPVGAS_AB(IGAS,2).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1)
     &        THEN
               ECC = ECC + 
     &         0.5D0*(  GTIJKL_SM_AB(IP,IP,JP,JP,0,4) 
     &                - GTIJKL_SM_AB(IP,JP,JP,IP,0,4) )
              END IF
              IF(IHPVGAS_AB(IGAS,1).EQ.1.AND.IHPVGAS_AB(JGAS,2).EQ.1) 
     &        THEN
               ECC = ECC + GTIJKL_SM_AB(IP,IP,JP,JP,2,2)
              END IF
          END DO
          END DO
        END DO
        END DO
      END IF
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Core-Core interaction energy ', ECC
      END IF
*
*.  Inactive Fock matrix
*
* Original one-electron matrix
      IF(IDOH2.NE.0) THEN 
*. Number of 1-e integrals
C  NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
        NINT = NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
        CALL COPVEC(HA,FI_AL,NINT)
        CALL COPVEC(HB,FI_BE,NINT)
        IJSM = 1
        IJ = 0
        DO ISM = 1, NSMOB
          CALL SYMCOM(2,6,ISM,JSM,IJSM)
          IJOFF = IJ + 1
          IF(JSM.NE.0) THEN
            DO J = IBSO(JSM),IBSO(JSM) + NTOOBS(JSM) - 1
              IF(ISYMPACK.EQ.0) THEN
                IMIN = IBSO(ISM)
              ELSE
                IMIN = J
              END IF
              DO I = IMIN,IBSO(ISM) + NTOOBS(ISM) - 1
                IP = IREOST(I)
                JP = IREOST(J)
                IF(ISYMPACK.EQ.0) THEN
                  IJ= IJ + 1
                ELSE 
                  IREL = I - IBSO(ISM) + 1
                  JREL = J - IBSO(JSM) + 1
                  IJ = IJOFF - 1 + IREL*(IREL-1)/2 + JREL
                END IF
                DO KSYM = 1, NSMOB
                  DO K = IBSO(KSYM),IBSO(KSYM)-1+NTOOBS(KSYM)
                    KP = IREOST(K)
                    KGAS = ITPFSO(K)
                    IF(IHPVGAS_AB(KGAS,1).EQ.1) THEN
                      FI_AL(IJ) = FI_AL(IJ)  
     &                          + GTIJKL_SM_AB(IP,JP,KP,KP,4,0)
     &                          - GTIJKL_SM_AB(IP,KP,KP,JP,4,0)
                      FI_BE(IJ) = FI_BE(IJ) 
     &                          + GTIJKL_SM_AB(KP,KP,IP,JP,2,2)
                    END IF
                    IF(IHPVGAS_AB(KGAS,2).EQ.1) THEN
                      FI_BE(IJ) = FI_BE(IJ)  
     &                          + GTIJKL_SM_AB(IP,JP,KP,KP,0,4)
     &                          - GTIJKL_SM_AB(IP,KP,KP,JP,0,4)
                      FI_AL(IJ) = FI_AL(IJ) 
     &                          + GTIJKL_SM_AB(IP,JP,KP,KP,2,2)
                    END IF
                  END DO
*                 ^ End of loop over K
                END DO
*               ^ End of Loop over KSYM
              END DO
            END DO
          END IF
        END DO
      END IF
*
      IF(NTEST.NE.0.AND.IDOH2.NE.0) THEN
       WRITE(6,*) ' FI_alpha in Symmetry blocked form '
       WRITE(6,*) ' =================================='
       WRITE(6,*) 
       ISYM = 0
       CALL APRBLM2(FI_AL,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
       WRITE(6,*) ' FI_beta in Symmetry blocked form '
       WRITE(6,*) ' ================================='
       WRITE(6,*) 
       ISYM = 0
       CALL APRBLM2(FI_BE,NTOOBS,NTOOBS,NSMOB,ISYMPACK)
      END IF
*
      RETURN
      END
      SUBROUTINE MXMNOC_OCCLS2(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,
     &                  MINOP,NTESTG)
*
* Construct accumulated MAX and MIN arrays for an occupation class
*
* MINOP ( Smallest allowed number of open orbitals) added 
* April2, 2003, JO
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Output
      DIMENSION  MINEL(*),MAXEL(*)
*. Input
      INTEGER NORBFTP(*),NELFTP(*)
*. Local scratch added April 2, 2003
*
      INCLUDE 'mxpdim.inc'
C     INTEGER MINOP_ORB(MXPORB), MINOP_GAS(MXPNGAS), MAXOP_GAS(MXPNGAS)
      INTEGER MINOP_GAS(MXPNGAS), MAXOP_GAS(MXPNGAS)
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==========='
        WRITE(6,*) ' MXMNOC_OCCLS'
        WRITE(6,*) ' ==========='
        WRITE(6,*)
        WRITE(6,*) ' NORBTP = ', NORBTP
        WRITE(6,*) ' NORBFTP : '
        CALL IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
      END IF
*
*. Largest number of unpaired electrons in each gas space 
*
      DO IGAS = 1, NORBTP
        MAXOP_GAS(IGAS) = MIN(NELFTP(IGAS),2*NORBFTP(IGAS)-NELFTP(IGAS))
      END DO
*
*. Smallest number of electrons in each GAS space 
*
*. 1 : Just based on number of electrons in each space 
      DO IGAS = 1, NORBTP
        IF(MOD(NELFTP(IGAS),2).EQ.1) THEN
          MINOP_GAS(IGAS) = 1
        ELSE
          MINOP_GAS(IGAS) = 0
        END IF
      END DO
*. 2 : the total number of open orbitals should be MINOP, this puts 
*. also a constraint on the number of open orbitals
*
*. The largest number of open orbitals, all spaces
      MAXOP_T = IELSUM(MAXOP_GAS,NORBTP)
      DO IGAS = 1, NORBTP
*. Max number of open orbitals in all spaces except IGAS
       MAXOP_EXL = MAXOP_T - MAXOP_GAS(IGAS)
       MINOP_GAS(IGAS) = MAX(MINOP_GAS(IGAS),MINOP-MAXOP_EXL)
      END DO
*. We now have the min and max number of open shells per occls,
*. Find the corresponding min and max number accumulated electrons, 
*
* The Max occupation is obtained by occupying in max in the 
* first orbitals 
* The Min occupation is obtained by occopying max in the 
* last orbitals.
*
      DO IGAS = 1, NORBTP
        MAX_DOUBLE = (NELFTP(IGAS)-MINOP_GAS(IGAS))/2
        IF(IGAS.EQ.1) THEN
          NEL_INI = 0
          IBORB = 1
        ELSE 
          NEL_INI = NEL_INI + NELFTP(IGAS)
          IBORB = IBORB + NORBFTP(IGAS)
        END IF
        NELEC = NELFTP(IGAS)
*
* The min number of electrons
*
*. Doubly occupy the last MAX_DOUBLE orbitals
        DO IORB = NORBFTP(IGAS)-MAX_DOUBLE,NORBFTP(IGAS)
           MINEL(IORB+IBORB-1) = 
     &     NEL_INI + NELEC - 2*(NORBFTP(IGAS)-IORB)
        END DO
*. Singly occupy 
        DO IORB = NORBFTP(IGAS)-MAX_DOUBLE-1,1,-1
           MINEL(IORB+IBORB-1) = MAX(NEL_INI,MINEL(IORB+IBORB-1+1)-1)
        END DO
*
*. The max number of electrons 
*
       DO IORB = 1, MAX_DOUBLE
         MAXEL(IORB+IBORB-1) = NEL_INI + 2*IORB 
       END DO 
       DO IORB = MAX_DOUBLE+1, NORBFTP(IGAS)  
         MAXEL(IORB+IBORB-1) = MIN(NEL_INI+NELEC,MAXEL(IORB+IBORB-2)+1)
       END DO
      END DO
*
C     IORB_START = 0
C     DO IORBTP = 1, NORBTP
*. Max and min at start of this type and at end of this type
C       IF(IORBTP.EQ.1) THEN
C         IORB_START = 1
C         IORB_END = NORBFTP(1)
C         NEL_START = 0
C         NEL_END   = NELFTP(1)
C       ELSE
C         IORB_START =  IORB_START + NORBFTP(IORBTP-1)
C         IORB_END   =  IORB_START + NORBFTP(IORBTP)-1
C         NEL_START = NEL_END
C         NEL_END   = NEL_START + NELFTP(IORBTP)
C       END IF
C       IF(NTEST.GE.1000) THEN
C         WRITE(6,*) ' IORBTP,IORB_START-IORB_END,NEL_START,NEL_END '
C         WRITE(6,*)   IORBTP,IORB_START-IORB_END,NEL_START,NEL_END
C       END IF
*
C       DO IORB = IORB_START, IORB_END
C         MAXEL(IORB) = MIN(2*IORB,NEL_END,MAXEL(IORB))
C         MINEL(IORB) = NEL_START
C         IF(NEL_END-MINEL(IORB).GT. 2*(IORB_END-IORB))
C    &    MINEL(IORB) = NEL_END - 2*( IORB_END - IORB )

C       END DO
C     END DO
*
      IF( NTEST .GE. 1000 ) THEN
        NORB = IELSUM(NORBFTP,NORBTP)
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_CONF_FOR_MAXMIN_OCC(IOCC_MIN,IOCC_MAX,NORBL,
     &    IORB_OFF,
     &    INITIALIZE_CONF_COUNTERS,
     &    ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,
     &    NCONF_OP,NCONF,IBCONF_REO,IBCONF_OCC,ICONF,
     &    IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM,IREO_MNMX_OB_NO)
*
* IONLY_NCONF = 1 :
*
* Generate number of configurations with occupation between IOCC_MIN, IOCC_MAX
* and sym ISYM
*
* IONLY_NCONF = 0 :
*
* Generate number and actual configurations with occupation 
* between IOCC_MIN, IOCC_MAX and sym ISYM
*
*
* Jeppe Olsen, June 2011
* Last modification; Jeppe Olsen; May 30 2013; IREO_MNMX_OB_NO added
*       
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
*
*.. Input
*
*. Min and max number of accumulated electrons
      INTEGER IOCC_MIN(NORBL),IOCC_MAX(NORBL)
*. Arc weights for configurations
      INTEGER IZ_CONF(*)
*. Offset for reordering array and occupation array
      INTEGER IBCONF_REO(*), IBCONF_OCC(*)
*. Reorder array of orbitals
      INTEGER IREO_MNMX_OB_NO(NORBL)
*
*.. Output
*
*. Number of configurations per number of open shells, all symmetries
      INTEGER NCONF_OP(MAXOP+1)
*. And the actual configurations
      INTEGER ICONF(*)
*. Reorder array : Lex number => Actual number 
      INTEGER IREO(*)
*. Local scratch
      INTEGER JOCC(2*MXPORB), JACOCC(MXPORB), JACOCC2(MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================================='
        WRITE(6,*) ' Entering GEN_CONF_FOR_MAXMIN_OCC '
        WRITE(6,*) ' ================================='
        WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NORBL, IDOREO = ', NORBL, IDOREO
        WRITE(6,*) ' INITIALIZE_CONF_COUNTERS = ',
     &               INITIALIZE_CONF_COUNTERS
        WRITE(6,*) ' Min and Max accumulated occupations '
        WRITE(6,*)
        WRITE(6,*) ' Orbital Min. occ Max. occ '
        WRITE(6,*) ' =========================='
        DO IORB = 1, NORBL
          WRITE(6,'(3X,I4,2I3)') 
     &    IORB, IOCC_MIN(IORB), IOCC_MAX(IORB)
        END DO
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,2I3)') 
     &  ' Smallest and largest number of singly occupied orbitals ',
     &  MINOP, MAXOP
        WRITE(6,*) ' IONLY_NCONF, IORB_OFF = ', IONLY_NCONF, IORB_OFF
      END IF
*
*. Total number of electrons 
      NEL = IOCC_MAX(NORBL)
      IF(NTEST.GE.1000) WRITE(6,*) ' NEL = ', NEL
      IF(INITIALIZE_CONF_COUNTERS.EQ.1) THEN
        IZERO = 0
        CALL ISETVC(NCONF_OP,IZERO,MAXOP+1)
        NCONF_ALL_SYM = 0
      END IF
*. Loop over configurations in the form of accumulated occupations
      INI = 1
      NCONF = 0
 1000 CONTINUE
        CALL NEXT_CONF_FROM_MINMAX_OCC(JACOCC,
     &             IOCC_MIN,IOCC_MAX,INI,NONEW,NORBL)
        INI = 0
        IF(NONEW.EQ.0) THEN
*. Reform from accumulated to occupation form
          CALL REFORM_CONF_ACCOCC(JACOCC,JOCC,1,NORBL)
C              REFORM_CONF_ACCOCC(IACCOCC,IOCC,IWAY,NORB)
*. Reform to actual ordering of orbitals in JACOCC2
C              REO_OB_CONFE(ICONFP_IN, ICONFP_UT,IREO_NO,NOB)
          CALL REO_OB_CONFE(JOCC,JACOCC2,IREO_MNMX_OB_NO,NORBL)
          CALL ICOPVE(JACOCC2,JOCC,NORBL)
C?        WRITE(6,*) ' Next conf in JOCC '
C?        CALL IWRTMA(JOCC,1,NORBL,1,NORBL)
*
*. Check symmetry and number of open orbitals for this space
          JSYM = ISYM_CONF(JOCC,NORBL,IORB_OFF)
C                ISYM_CONF(IOCC,NORBL,IORB_OFF)
          NOPEN = NOP_FOR_CONF_OCC(JOCC,NORBL)
          NOCOBL = NOPEN + (NEL-NOPEN)/2
          IF(NTEST.GE.10000)
     &    WRITE(6,*) ' Number of open and occupied orbitals ', 
     &    NOPEN, NOCOBL
*
          
          NCONF_ALL_SYM = NCONF_ALL_SYM + 1 
          IF(JSYM.EQ.ISYM.AND.NOPEN.GE.MINOP) THEN
*. A new configuration to be included, reform and save in packed form
            NCONF = NCONF + 1
            NCONF_OP(NOPEN+1) = NCONF_OP(NOPEN+1) + 1
            IF(IONLY_NCONF .EQ. 0 ) THEN
*. Lexical number of this configuration 
              IB_OCC = IBCONF_OCC(NOPEN+1) 
     &               + (NCONF_OP(NOPEN+1)-1)*NOCOBL
              IF(NTEST.GE.10000) 
     &                     WRITE(6,*) ' IBCONF_OCC,NCONF_OP,IB_OCC =',
     &                     IBCONF_OCC(NOPEN+1),NCONF_OP(NOPEN+1),
     &                     IB_OCC
C                  REFORM_CONF_OCC2(ICONF_EXP,ICONF_PACK,NORBL,NOCOBL,IWAY)
              CALL REFORM_CONF_OCC2(JOCC,ICONF(IB_OCC),NORBL,NOCOBL,1)
              IF(IDOREO.NE.0) THEN 
C                          ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,
C                                        IARCW,IDOREO,IREO)
*. Okay, for graphs we are using the original order, so we return to this
              CALL REFORM_CONF_ACCOCC(JACOCC,JOCC,1,NORBL)
              CALL REFORM_CONF_OCC2(JOCC,JACOCC2,NORBL,NOCOBL,1)
*
                ILEXNUM =  ILEX_FOR_CONF(JACOCC2,NOCOBL,NORBL,
     &                     NEL,IZ_CONF,0,IDUM)
                JREO = IBCONF_REO(NOPEN+1) -1 + NCONF_OP(NOPEN+1)
                IREO(ILEXNUM) = JREO
                IF(NTEST.GE.10000) THEN
                 WRITE(6,*) ' Next configuration: '
                 CALL IWRTMA(ICONF(IB_OCC),1,NOCOBL,1,NOCOBL)
                 WRITE(6,*) ' LEXCONF : JREO, ILEXNUM = ', JREO,ILEXNUM
                 WRITE(6,*) ' ILEXNUM, JREO= ', ILEXNUM, JREO 
                END IF
              END IF ! End if IDOREO
            END IF ! End if correct sym and number of open orbitals
          END IF ! End of IONLY_NCONF = 0
      GOTO 1000
        END IF !End if nonew = 0
* 
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*)  ' ====================================== '
        WRITE(6,*)  ' Results from configuration generator : '
        WRITE(6,*)  ' ====================================== '
        WRITE(6,*)
        WRITE(6,*) ' Number of configurations of correct symmetry ',
     &       NCONF
        WRITE(6,*) ' Number of configurations of all symmetries   ',
     &       NCONF_ALL_SYM
        WRITE(6,*) 
     &  ' Number of configurations for various number of open orbs'
        CALL IWRTMA(NCONF_OP,1,MAXOP+1,1,MAXOP+1)
      END IF
*
      IF(NTEST.GE.1000) THEN
        IF(IONLY_NCONF.EQ.0) THEN
          WRITE(6,*) 
     &   ' Updated list of configurations (may not be the final...)'
          CALL WRT_CONF_LIST(ICONF,NCONF_OP,MAXOP,NCONF,NEL)
          WRITE(6,*) 
     &   ' Updated reordering of conf, Lex=>Act (may not be the final'
          CALL IWRTMA(IREO,1,NCONF_ALL_SYM,1,NCONF_ALL_SYM)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE NEXT_CONF_FROM_MINMAX_OCC(IACOCC,
     &           IACOCC_MIN,IACOCC_MAX,INI,NO_NEW,NORB)
*
*. Next accumulated occupation in CI space with accumulated
*   occupations between IOCC_MIN,IOCC_MAX
*  IF INI = 1, generate first configuration
*  NO_NEW = 1 indicates no new configurations
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IACOCC_MIN(NORB),IACOCC_MAX(NORB)
*. Output
      INTEGER IACOCC(NORB)
*
      NTEST = 000
C?    WRITE(6,*) ' TEST: NORB = ', NORB
      IF(NTEST.GE.1000.AND.INI.EQ.0) THEN
        WRITE(6,*) ' Input configuration to NEXT_CONF..'
        CALL IWRTMA(IACOCC,1,NORB,1,NORB)
      END IF
*
*. A configuration is considered as the integer with digit I being
*  IACCOCC(I). The initial configuration corresponds to the lowest
*  possible number and the following configurations are obtained as 
*  increasing integers
      NO_NEW = 1
      IF(INI.EQ.1) THEN
*
* ==============================
* Generate initial configuration
* ==============================
*
*. The lowest possible integer has the lowest possible digits in all positions,
*  i.e. the min occupation
        DO IORB = 1, NORB
          IACOCC(IORB) = IACOCC_MIN(IORB)
        END DO
        NO_NEW = 0
      ELSE
*
* ==============================
* Generate NEXT configuration
* ==============================
*
*. Find the first orbital where the accumulated number of electrons 
*. may be increased
       IMOD = 0
       DO IORB = 1, NORB-1
         IF(IACOCC(IORB).LT.IACOCC(IORB+1) .AND.
     &      IACOCC(IORB).LT.IACOCC_MAX(IORB)) THEN
            IMOD = IORB
            GOTO 101
         END IF
       END DO
 101   CONTINUE
       IF(IMOD.NE.0) THEN
         NO_NEW = 0
*. Increase occupation in orbital IMOD
         IACOCC(IMOD) = IACOCC(IMOD)+1
*. Minimize occupation in all previous orbitals
         DO IORB = IMOD-1,1,-1
           IACOCC(IORB) = MAX(IACOCC(IORB+1)-2,IACOCC_MIN(IORB))
         END DO
       END IF
      END IF ! end of INI switch 
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output from NEXT_CONF_FROM_MINMAX_OCC: '
        IF(NO_NEW.EQ.1) THEN
          WRITE(6,*) ' No new configurations '
        ELSE
          WRITE(6,*) ' Accumulated occupation of next configuration '
          CALL IWRTMA(IACOCC,1,NORB,1,NORB)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_CONF_ACCOCC(IACOCC,IOCC,IWAY,NORB)
*
* Reform between accumulated and actual occupation of configuration
* IWAY = 1: Accumulated => Actual in expanded form (occnumber) 
* IWAY = 2: Actual => Accumulated
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER IACOCC(NORB),IOCC(NORB)
*
      NTEST = 000
*
      IF(IWAY.EQ.1) THEN 
*. Accumulated => Actual
        IOCC(1) = IACOCC(1)
        DO IORB = 2, NORB
          IOCC(IORB) = IACOCC(IORB)-IACOCC(IORB-1)
        END DO
      ELSE
*. Actual => Accumulated
        IACOCC(1) = IOCC(1)
        DO IORB = 2, NORB
          IACOCC(IORB) = IACOCC(IORB-1)+IOCC(IORB)
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
       IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' Accumulated to occupation form of config'
       ELSE
         WRITE(6,*) ' Occupation to accumulated form of config'
       END IF
*
       WRITE(6,*) ' Configuration in accumulated form '
       CALL IWRTMA(IACOCC,1,NORB,1,NORB)
       WRITE(6,*) ' Configuration in occupation form '
       CALL IWRTMA(IOCC,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
      FUNCTION ISYM_CONF(IOCC,NORBL,IORB_OFF)
*
* Obtain symmetry of configuration IOCC, specified as 
* occupations of each orbital - i.e as occupation number vector
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. General input
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
* Specific input
      DIMENSION IOCC(NORBL)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ISYM_CONF:'
        WRITE(6,*) ' NORBL, IORB_OFF =', NORBL,IORB_OFF
        WRITE(6,*) '  IOCC: '
        CALL IWRTMA3(IOCC,1,NORBL,1,NORBL)
      END IF
*
      ISYM = 1
      DO IORB = 1, NORBL
        IF(IOCC(IORB).EQ.1)
     &  ISYM = MULTD2H(ISYM,ISMFTO(IORB-1+IORB_OFF))
      END DO
*
       ISYM_CONF = ISYM
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from ISYM_CONF '
       WRITE(6,*) ' ======================'
       WRITE(6,*) ' Configuration: '
       CALL IWRTMA(IOCC,1,NORBL,1,NORBL)
       WRITE(6,*) ' Symmetry = ', ISYM
      END IF
*
      RETURN
      END
   
      FUNCTION NOP_FOR_CONF_OCC(IOCC,NORBL)
*
* Number of singly occupied orbitals in configuration IOCC
* given as orbital occuopation numbers
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INTEGER IOCC(NORBL)
*
      NOP = 0
      DO IORB = 1, NORBL
        IF(IOCC(IORB).EQ.1) NOP = NOP + 1
      END DO
*
      NOP_FOR_CONF_OCC = NOP
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output from NOP_FOR_CONF_OCC '
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(IOCC,1,NORBL,1,NORBL)
        WRITE(6,*) ' Number of singly occupied orbitals ', NOP
      END IF
*
      RETURN
      END
      SUBROUTINE REFORM_CONF_OCC2(ICONF_EXP,ICONF_PACK,NORBL,NOCOBL,
     &                           IWAY)
*
* Another routine for packing and unpacking configuration
*
* The two forms:
* ICONF_EXP: Occupation numbers, 0,1,2 for each orbital
* ICONF_PACK: The occupied orbitals, negative value for doubly occ
*
* IWAY = 1: Expanded => Packed 
* IWAY = 2: Packed => Expanded
*
* Jeppe Olsen, June 2011  
*
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER ICONF_EXP(NORBL), ICONF_PACK(NOCOBL)
*
      IF(IWAY.EQ.1) THEN
*. Expanded to packed
       IOCOB = 0
       DO IORB = 1, NORBL
         IF(ICONF_EXP(IORB).EQ.1) THEN
           IOCOB = IOCOB + 1
           ICONF_PACK(IOCOB) = IORB
         ELSE IF (ICONF_EXP(IORB).EQ.2) THEN
           IOCOB = IOCOB + 1
           ICONF_PACK(IOCOB) = -IORB
         END IF
       END DO
      ELSE
*. Packed to expanded form
      IZERO = 0
      CALL ISETVC(ICONF_EXP,IZERO,NORBL)
      DO IOCOB = 1, NOCOBL
        IF(ICONF_PACK(IOCOB).GT.0) THEN
          IORB = ICONF_PACK(IOCOB)
          ICONF_EXP(IORB) = 1
        ELSE
          IORB = -ICONF_PACK(IOCOB)
          ICONF_EXP(IORB) = 2
        END IF
       END DO
      END IF ! End of IWAY switch
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from REFORM_CONF_OCC2 '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Expanded => Packed '
        ELSE
          WRITE(6,*) ' Packed => Expanded '
        END IF
*
        WRITE(6,*) ' Expanded form of configuration: '
        CALL IWRTMA(ICONF_EXP,1,NORBL,1,NORBL)
        WRITE(6,*) ' Packed form of configuration: '
        CALL IWRTMA(ICONF_PACK,1,NOCOBL,1,NOCOBL)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_CONF_FOR_MINMAX_SPC(IOCC_MIN,IOCC_MAX,NORBL,
     &           ISYM,IORB_IB,ISPC)
*
* Generate configuration information for CI space with symmetry ISYM 
* defined by  min and max accumulated occupations
* The orbitals are occupied in the order specified by IREO_MNMX_OB_NO
*
* Information is stored in arrays adressed by pointers ISPC
*
* Jeppe Olsen, June 2011
* Last modification; May 30, 2013; Jeppe Olsen; IREO_MNMX_OB_NO added
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*. Input
      INTEGER IOCC_MIN(NORBL),IOCC_MAX(NORBL)
*. Local scratch 
C     INTEGER ICONF(2*MXPORB)
*. Local scratch for testing 
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 
        WRITE(6,*) ' ===================================='
        WRITE(6,*) ' Output from GEN_CONF_FOR_MINMAX_SPC '
        WRITE(6,*) ' ===================================='
        WRITE(6,*) 
      END IF
*
*. Total number of electrons 
      NELEC = IOCC_MIN(NORBL)
*. Min and max number of open orbitals - read from SPINFO
      MINOP_GN(ISPC) = MINOP
*. Max number of open orbitals - just upper limit, not memory critical
      MAXOP_GN(ISPC) = MAXOP
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Max number of open orbitals ', MAXOP
        WRITE(6,*) ' Min number of open orbitals ', MINOP
      END IF
*
*. Number of configurations and length of configuration list
*. for all occupation classes
*
      IDUM = 0
      INITIALIZE_CONF_COUNTERS = 1
*
      IDOREO = 0
      IONLY_NCONF = 1
C     GEN_CONF_FOR_MAXMIN_OCC(IOCC_MIN,IOCC_MAX,NORBL,
C    &    INITIALIZE_CONF_COUNTERS,
C    &    ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,
C    &    NCONF_OP,NCONF,IBCONF_REO,IBCONF_OCC,ICONF,
C    &    IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM,IREO_MNMX_OB_NO)
      CALL GEN_CONF_FOR_MAXMIN_OCC(IOCC_MIN,IOCC_MAX,NORBL,IORB_IB,
     &     INITIALIZE_CONF_COUNTERS,
     &     ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,
     &     NCONF_PER_OPEN_GN(1,ISYM,ISPC),NCONF,
     &     IB_CONF_REO_GN(1,ISYM,ISPC),IB_CONF_OCC_GN(1,ISYM,ISPC),
     &     IDUM,IDOREO,KZ_GN(ISPC),IDUM,NCONF_ALL_SYM_GN(ISPC),
     &     IREO_MNMX_OB_NO)
*. Total number of configurations and length of configuration list
      CALL ICOPVE(NCONF_PER_OPEN_GN(1,ISYM,ISPC),
     &            NCONF_PER_OPEN(1,ISYM),MAXOP+1)
*. In INFO_CONF_LIST, the number of active electrons is used. 
*. Remove tempoary those connected with CORE spaces
      NACTEL = NACTEL - 2*(IORB_IB-NINOB-1)
      CALL INFO_CONF_LIST(ISYM,LENGTH_LIST,NCONF_TOT)
      NACTEL = NACTEL + 2*(IORB_IB-NINOB-1)
C    &     IB_CONF_REO_GN(1,ISYM,ISPC),IB_CONF_OCC_GN(1,ISYM,ISPC))
      CALL ICOPVE(IB_CNOCC_OPEN,IB_CONF_OCC_GN(1,ISYM,ISPC),MAXOP+1)
C?    WRITE(6,*) ' IB_CONF_OCC_GN(1,ISYM,ISPC) after ICOPVE '
C?    CALL IWRTMA(IB_CONF_OCC_GN(1,ISYM,ISPC),1,MAXOP+1,1,MAXOP+1)
      CALL ICOPVE(IB_CN_OPEN,IB_CONF_REO_GN(1,ISYM,ISPC),MAXOP+1)
*. Add offsets to start in CSF CI vector of CSF's  with a given
* number of unpaired electrons
*. ISYM is pt implied
      CALL OFFSETS_CSF_FOR_NOPEN(IB_OPEN_CSF(1,ISYM,ISPC),
     &                           NCONF_PER_OPEN_GN(1,ISYM,ISPC),
     &                           MAXOP, NPCSCNF)
C     OFFSETS_CSF_FOR_NOPEN(IB_OPEN_CSF,NCONF_OPEN,MAXOP,
C    &                                 NPCSCNF)
*
*. memory for storing configuration info 
*
*. 1 : Occupation of configurations 
      CALL MEMMAN
     &(KICONF_OCC_GN(ISYM,ISPC),LENGTH_LIST,'ADDL  ',1,'I_CONF')
*. 2 : Reordering of configurations of correct sym : lex. => act.
      CALL MEMMAN
     &(KICONF_REO_GN(ISYM,ISPC),NCONF_ALL_SYM_GN(ISPC),
     & 'ADDL  ',1,'REOCON')
*. 3 : Array giving start of each occupation class 
*. 4: Arc- weights for setting lexical addressing of space
      CALL MEMMAN(KZ_GN(ISPC),NORBL*NELEC*2,'ADDL  ',1,'Z     ')
*. Scratch memory for setting up configurations 
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GEN_CO')
      CALL MEMMAN(KLZSCR,(NORBL+1)*(NELEC+1),'ADDL  ',1,'ZSCR  ')  !done
*
* Set up the configurations 
*
*. Arcweights
      CALL CONF_GRAPH(IOCC_MIN,IOCC_MAX,NORBL,
     &     NELEC,int_mb(KZ_GN(ISPC)),NCONF_P,int_mb(KLZSCR))
*
      INITIALIZE_CONF_COUNTERS = 1
      IDOREO = 1
      IONLY_NCONF = 0
*. Generate configurations and lexical addressing arrays 
*. for configurations of this type 
      IZERO = 0
      CALL ISETVC(int_mb(KICONF_REO_GN(ISYM,ISPC)),IZERO,
     &            NCONF_ALL_SYM_GN(ISPC))
      CALL GEN_CONF_FOR_MAXMIN_OCC(IOCC_MIN,IOCC_MAX,NORBL,IORB_IB,
     &     INITIALIZE_CONF_COUNTERS,
     &     ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,
     &     NCONF_PER_OPEN_GN(1,ISYM,ISPC),NCONF,
     &     IB_CONF_REO_GN(1,ISYM,ISPC),IB_CONF_OCC_GN(1,ISYM,ISPC),
     &     int_mb(KICONF_OCC_GN(ISYM,ISPC)),IDOREO,int_mb(KZ_GN(ISPC)),
     &     int_mb(KICONF_REO_GN(ISYM,ISPC)),
     &     NCONF_ALL_SYM_GN(ISPC),IREO_MNMX_OB_NO)
*
C?         WRITE(6,*) ' KICONF_REO_GN(ISYM,ISPC) = ', 
C?   &     KICONF_REO_GN(ISYM,ISPC)
C?         CALL IWRTMA(WORK(KICONF_REO_GN(ISYM,ISPC)),1,1,1,1)
         
C     GEN_CONF_FOR_MAXMIN_OCC(IOCC_MIN,IOCC_MAX,NORBL,
C    &    IORB_OFF,
C    &    INITIALIZE_CONF_COUNTERS,
C    &    ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,
C    &    NCONF_OP,NCONF,IBCONF_REO,IBCONF_OCC,ICONF,
C    &    IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM,IREO_MNMX_OB_NO)

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GEN_CO')
      RETURN
      END
      SUBROUTINE PROTO_CONF_DIM(NELEC)
*
* Find dimensions for prototype configurations for system with 
* given MULTS, MS2, MINOP, MAXOP (read from common files)
*
* I.e. number of determinants, combinations and CSF's per prototype 
*
*. Output is in SPINFO: NPDTCNF, NPCMCNF, NPDSCNF, IB_PDTCNF, IB_PCSCNF,IB_PCMCNF
*
*. Jeppe Olsen, June 2011, Modifications, Dec. 2011
*
*. Last modification;  Oct. 2012; Jeppe Olsen, Dimensions for IOPEN < MINOP set to zero
 
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cprnt.inc'
      NTEST = 00
      NTEST = MAX(IPRCSF,NTEST)
      IF(NTEST.GE.100) WRITE(6,*) ' MULTS, MS2 = ', MULTS,MS2
C
C.. Number of prototype sd's and csf's per configuration prototype
C
*
      IZERO = 0
      CALL ISETVC(NPDTCNF,IZERO,MAXOP)
      CALL ISETVC(NPCMCNF,IZERO,MAXOP)
      CALL ISETVC(NPCSCNF,IZERO,MAXOP)
*
      ITP = 0
      DO IOPEN = MINOP, MAXOP
        ITP = IOPEN + 1
*. Unpaired electrons :
        IAEL = (IOPEN + MS2 ) / 2
        IBEL = (IOPEN - MS2 ) / 2
        IF(IAEL+IBEL .EQ. IOPEN .AND. IAEL-IBEL .EQ. MS2 .AND.
     &            IAEL .GE. 0 .AND. IBEL .GE. 0) THEN
          NPDTCNF(ITP) = IBION(IOPEN,IAEL)
          IF(PSSIGN.EQ. 0.0D0 .OR. IOPEN .EQ. 0 ) THEN
            NPCMCNF(ITP) = NPDTCNF(ITP)
          ELSE   
            NPCMCNF(ITP) = NPDTCNF(ITP)/2
          END IF
          IF(IOPEN .GE. MULTS-1) THEN
            NPCSCNF(ITP) = IWEYLF(IOPEN,MULTS)
          ELSE
            NPCSCNF(ITP) = 0
          END IF
        ELSE
          NPDTCNF(ITP) = 0 
          NPCMCNF(ITP) = 0 
          NPCSCNF(ITP) = 0 
        END IF
      END DO
*. Construct the corresponding offset arrays
      CALL ZBASE(NPDTCNF,IB_PDTCNF,MAXOP+1)
      CALL ZBASE(NPCSCNF,IB_PCSCNF,MAXOP+1)
      CALL ZBASE(NPCMCNF,IB_PCMCNF,MAXOP+1)
*
      IF(NTEST.GE.10) THEN
      IF(PSSIGN .EQ. 0 ) THEN
        WRITE(6,*) '  (Combinations = Determinants ) '
      ELSE
        WRITE(6,*) '  (Spin combinations in use ) '
      END IF
      WRITE(6,'(/A)') ' Information about prototype configurations '
      WRITE(6,'( A)') ' ========================================== '
      WRITE(6,'(/A)')
     &'  Open orbitals   Combinations    CSFs '
      DO IOPEN = MINOP,MAXOP,2
        WRITE(6,'(5X,I3,10X,I6,7X,I6)')
     &  IOPEN,NPCMCNF(IOPEN+1),NPCSCNF(IOPEN+1)
      END DO
*
      END IF
*
      RETURN
      END
      SUBROUTINE NPARA_FOR_MINMAX_SPC(NCONF_OP,NCSF,NSD,NCMB,NCNF)
*
* Number of CSF's, SD's, combinations for system with given NCONF_OP
* Using prototype info in SPINFO 
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
*
      NTEST = 000
*. Number of CSF's in expansion 
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_OP,NPCSCNF,NCSF)
*. Number of SD's in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_OP,NPDTCNF,NSD)
*. Number of combinations in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_OP,NPCMCNF,NCMB) 
*. Number of configurations
      NCNF = IELSUM(NCONF_OP,MAXOP+1)
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Number of SDs   ', NSD
        WRITE(6,*) ' Number of CMBs  ', NCMB
        WRITE(6,*) ' Number of CSFs  ', NCSF
        WRITE(6,*) ' Number of Confs ', NCNF
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_PROTO_ARRAYS
*
* Length of arrays for defining occupation of prototype determinants
* and CSF's and transformation between these
*
* Using prototype info in SPINFO 
*
* Output data are store in csfbas
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'csfbas.inc'
*
      LPDT_OCC = 0
      LPCS_OCC = 0
      LPDTOC = 0
      MXPTDT = 0
      MXPOCCBL = 0
      DO IOPEN = 0, MAXOP
        ITP = IOPEN + 1
        LPDT_OCC = LPDT_OCC + NPCMCNF(ITP) * IOPEN
        LPCS_OCC = LPCS_OCC + NPCSCNF(ITP) * IOPEN
        LPDTOC= LPDTOC + NPCSCNF(ITP)*NPCMCNF(ITP)
        MXPTDT =   MAX(MXPTDT,NPCMCNF(ITP) )
        MXPOCCBL = MAX(NPCMCNF(ITP)*IOPEN,MXPOCCBL)
      END DO
*
      NTEST = 0000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 
        WRITE(6,*) ' Dimension for prototype arrays:'
        WRITE(6,*) ' Combinations = ', LPDT_OCC
        WRITE(6,*) ' CSF"s        = ', LPCS_OCC
        WRITE(6,*) ' CSF <-> Comb = ', LPDTOC
      END IF
*
      RETURN
      END
      FUNCTION LCONF_OCC(NCONF_PER_OP,ISYM,NELEC)
* 
* Length of array storing occupation of configurations in packed form
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
*. Specific input
      INTEGER NCONF_PER_OP(*)
*
      NTEST = 000
*
      LOCC = 0
      DO IOPEN = 0, MAXOP
        ITYP = IOPEN + 1
        ICL = ( NELEC-IOPEN)/2
        LOCC = LOCC + NCONF_PER_OP(ITYP)*(IOPEN+ICL)
      END DO
*
      LCONF_OCC = LOCC
*
      IF(NTEST.GE.100) THEN
      WRITE(6,'(/A,I8)')
     &'  Memory for holding list of configurations ',LOCC
      END IF
*
      RETURN
      END
      SUBROUTINE CNFORD2(ISM,ICTSDT,ICONF_OCC,NCONF_PER_OP,
     &           IDFTP,ICONF_ORBSPC,ICISPC)
C
C
C Generate determinants in configuration order and obtain
C sign array for switching between the two formats and store info in ICTSDT
C
C Jeppe Olsen June 2011, slight modification of CNFORD for nonorthogonal CI
C
*. NORBL is the number of orbitals in space defining configurations
*  NEL_CONF is the number of electrons in configurations
*
* The configurations are assumed to be defined in a single GASpace, IACTSPC
* and the occupations in ICONF_OCC refers to these orbitals, with first
* orbital having index 1.
* The orbitals in the preceeding GASpaces are assumed to be doubly occupied
* and this is preceeded by IB_ORB-1 orbitals, all 
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cstate.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER ICONF_OCC(*), NCONF_PER_OP(MAXOP+1)
      INTEGER IDFTP(*)
*. Output
      INTEGER ICTSDT(*)
*. Local scratch
      INTEGER IOC_AL(MXPNGAS),IOC_BE(MXPNGAS), IB_FOR_ALSM(8)
      INTEGER IOC_STRING(MXPORB)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CNFORD')
*
      NTEST = 0000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from CNFORD2 '
        WRITE(6,*) ' =================== '
        IF(NTEST.GE.10000) 
     &  WRITE(6,*) ' ICONF_ORBSPC = ', ICONF_ORBSPC
      END IF

*
*. Number of occupation classes
      IATP = 1
      IBTP = 2
*
      NAEL = NELFTP(IATP)
      NBEL = NELFTP(IBTP)
*
* ===================================
* Info on configuration orbital space
* ===================================
*
*. First orbital 
      IORB_IB= NINOB + 1
      DO IOBSPC = 1, ICONF_ORBSPC-1
        IORB_IB = IORB_IB + NOBPT(IOBSPC)
      END DO
      NORB_CORE = IORB_IB - 1
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' IORB_IB, NORB_CORE = ', IORB_IB, NORB_CORE
      END IF
*. Number of electrons in Configuration orbital space
      NAEL_CONF = NAEL - (NORB_CORE-NINOB)
      NBEL_CONF = NBEL - (NORB_CORE-NINOB)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NAEL, NBEL = ', NAEL, NBEL
        WRITE(6,*) ' NAEL_CONF, NBEL_CONF = ',
     &               NAEL_CONF, NBEL_CONF
      END IF
*. The corresponding string-groups
      IGRP_AL = ISTR_GRP(ICONF_ORBSPC,NAEL_CONF)
      IGRP_BE = ISTR_GRP(ICONF_ORBSPC,NBEL_CONF)
C               ISTR_GRP(IGAS,IEL)
*
      NEL_CONF = NAEL_CONF + NBEL_CONF
*. Occupation of alpha and beta-strings
      IZERO = 0
      CALL ISETVC(IOC_AL,IZERO,NGAS)
      CALL ISETVC(IOC_BE,IZERO,NGAS)
*. Orbital spaces before ICONF_ORBSPC is completely filled
      DO IGAS = 1, ICONF_ORBSPC
        IOC_AL(IGAS) = NOBPT(IGAS)
        IOC_BE(IGAS) = NOBPT(IGAS)
      END DO
      IOC_AL(ICONF_ORBSPC) = NAEL_CONF
      IOC_BE(ICONF_ORBSPC) = NBEL_CONF
*. supergroup types corresponding to these occupation
C     GET_SPGP_FROM_OCC(ISPGP,IGASOCC)
      CALL GET_SPGP_FROM_OCC(ISPGP_AL,IOC_AL,IATP)
      CALL GET_SPGP_FROM_OCC(ISPGP_BE,IOC_BE,IBTP)
      WRITE(6,*) ' Supergroup numbers ', 
     &ISPGP_AL, ISPGP_BE
*. Symmetry of string with constant occupation preceeding Conf. orb. space
C     SET_IARR_LIN_FUNC(IARR,IA,IB,NDIM)
      CALL SET_IARR_LIN_FUNC(IOC_STRING,NINOB,1,NORB_CORE-NINOB)
      ISYM_CORE_STRING = ISYMST(IOC_STRING,NORB_CORE-NINOB)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Symmetry of core string = ', ISYM_CORE_STRING
      END IF
*
*. ================================================
*. Generate information about string-ordered space
*. ================================================
*
*
      NEL = NAEL + NBEL
      IWAY = 1
      CALL OCCLS(1,NOCCLS,IOCCLS,NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),0,0,NOBPT)
*. and the occupation classes
      CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL')  !done
      CALL MEMMAN(KLBASSPC,NOCCLS,'ADDL  ',1,'BASSPC')  !done
      IWAY = 2
      CALL OCCLS(2,NOCCLS,int_mb(KLOCCLS),NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),1,int_mb(KLBASSPC),NOBPT)
*. String-block structure of the CI expansion
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Allocate space for largest encountered number of TTSS blocks
      NTTS = MXNTTS
C?    WRITE(6,*) ' GASCI : NTTS = ', NTTS
*. 
      CALL MEMMAN(KLCLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')  !done
      CALL MEMMAN(KLCLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')  !done
      CALL MEMMAN(KLCI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')  !done
      CALL MEMMAN(KLCIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')  !done
c.. dongxia: klc2b is not used. comment off
c     CALL MEMMAN(KLC2B  ,  NTTS,'ADDL  ',1,'C2BT  ')
c..
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')  !done
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')  !done
*. Matrix giving allowed combination of alpha- and beta-strings
      CALL IAIBCM(ICISPC,dbl_mb(KLCIOIO))
*. option KSVST not active so
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,dbl_mb(KLCBLTP),int_mb(KSVST))
*. Blocks of  CI vector, using a single batch for complete  expansion
      ICOMP = 1
      ISIMSYM = 1
      CALL PART_CIV2(IDC,dbl_mb(KLCBLTP),int_mb(KNSTSO(IATP)),
     &              int_mb(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,dbl_mb(KLCIOIO),
     &              ISMOST(1,ISM),
     &              NBATCH,int_mb(KLCLBT),int_mb(KLCLEBT),
     &              int_mb(KLCI1BT),int_mb(KLCIBT),ICOMP,ISIMSYM)
*. Number of BLOCKS
        NBLOCK = int_mb(KLCI1BT+nbatch-1)+int_mb(KLCLBT+nbatch-1)-1
        IF(NTEST.GE.1000) WRITE(6,*) ' Number of blocks ', NBLOCK
* and back to the configurations: 
*.the blocks in the CI vector with the given supergroups
C     OFFSETS_FOR_ALBESPGP(IBLOCK,NBLOCK,IALSPGP,IBESPGP,
C    &                                  IBLOCK_FOR_SPGP)
      CALL OFFSETS_FOR_ALBESPGP(int_mb(KLCIBT),NBLOCK,
     &     ISPGP_AL,ISPGP_BE,IB_FOR_ALSM)
*
* Loop over configurations and generate determinants, split these to
* strings and determine their symmetry and type
*
*. Space for storing determinants of a configuration 
      CALL MEMMAN(KLDET1,NEL_CONF,'ADDL  ',1,'CNFDET')  !done
      CALL MEMMAN(KLAL,NEL_CONF,'ADDL  ',1,'CNFAL ')  !done
      CALL MEMMAN(KLBE,NEL_CONF,'ADDL  ',1,'CNFBE ')  !done
      CALL MEMMAN(KLSCR,NEL_CONF,'ADDL  ',1,'CNFSCR')  !done
*
      IB_OCC = 1
      IDET_CONF = 0
      IB_PDT = 1
      DO IOPEN = 0, MAXOP
        ITP = IOPEN + 1
        NOCOBL = IOPEN + (NEL_CONF-IOPEN)/2
        NPDT = NPDTCNF(ITP)
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' IOPEN, NPDT, NCONF.. = ',
     &               IOPEN, NPDT, NCONF_PER_OP(ITP)
        DO ICONF = 1, NCONF_PER_OP(ITP)
          
*. The configuration is given in compact form with double occupied 
*. orbitals flagged by negative number. 
          DO JPDT = 1, NPDT
*. Obtain determinant in WORK(KLDET1)
C               GETDET_FOR_CONF(IOCC,NOCOB,NOPEN,NELEC,JPDET,IPDET,IDET,
C               IDOREO,IREO)
           CALL GETDET_FOR_CONF(ICONF_OCC(IB_OCC),NOCOBL,
     &          IOPEN,NEL_CONF,JPDT,IDFTP(IB_PDT),int_mb(KLDET1),1,
     &          IREO_SPCP_OB_ON)
*. Divide string into alpha- and beta-parts
C              DETSTR(IDET,IASTR,IBSTR,NAEL,NBEL,ISIGN,IWORK)
           CALL DETSTR(int_mb(KLDET1),int_mb(KLAL),int_mb(KLBE),
     &                NAEL_CONF,NBEL_CONF,ISIGN,int_mb(KLSCR),
     &                NTEST)
*. Symmetry of strings in configuration
*. Modify orbital indeces so they correspond to absolute orbital numbers
           IADD = IORB_IB - 1
C               ADD_C_TO_INTVEC(INTVEC,IC,NDIM)
           CALL ADD_C_TO_INTVEC(int_mb(KLAL),IADD,NAEL_CONF)
           CALL ADD_C_TO_INTVEC(int_mb(KLBE),IADD,NBEL_CONF)
           IAL_SYM_CONF = ISYMST(int_mb(KLAL),NAEL_CONF)
           IBE_SYM_CONF = ISYMST(int_mb(KLBE),NBEL_CONF)
*. Number of alpha-strings with this symmetry in group
           NAL_STR_CONF = NSTFSMGP(IAL_SYM_CONF,IGRP_AL)
*. Address of alpha- and beta- configuration strings
C                  ISTRNM2(IOCC,NORB,NEL,Z,NEWORD,IOFFSETS,IREORD,IRELNUM)
           IALNM = ISTRNM2(int_mb(KLAL),NOCOB,NAEL_CONF,
     &             int_mb(KZ(IGRP_AL)),int_mb(KSTREO(IGRP_AL)),
     &             ISTFSMGP(1,IGRP_AL),1,1)
           IBENM = ISTRNM2(int_mb(KLBE),NOCOB,NBEL_CONF,
     &             int_mb(KZ(IGRP_BE)),int_mb(KSTREO(IGRP_BE)),
     &             ISTFSMGP(1,IGRP_BE),1,1)
*. Symmetry of complete strings including core
           IAL_SYM = MULTD2H(IAL_SYM_CONF,ISYM_CORE_STRING)
*. Offset to symmetryblock
           IF(IB_FOR_ALSM(IAL_SYM).EQ.0) THEN
             WRITE(6,*) 
     &       ' Problem unknown symmetry-block in CI-expansion',
     &       IAL_SYM
             STOP
     &       ' Problem unknown symmetry-block in CI-expansion' 
           ELSE
             IB = IB_FOR_ALSM(IAL_SYM)
           END IF
*
           IDET_CONF = IDET_CONF + 1
           IDET_STR =  IB-1 + (IBENM-1)*NAL_STR_CONF+ IALNM
           IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' IB, NAL_STR_CONF = ',
     &                    IB, NAL_STR_CONF
             WRITE(6,*) '  IALNM, IBENM, IDET_STR = ',
     &                     IALNM, IBENM, IDET_STR
           END IF
*
           ICTSDT(IDET_CONF) = ISIGN*IDET_STR
*
          END DO! End of loop over prototype dets
          IB_OCC = IB_OCC + NOCOBL
        END DO! End of loop over Configs
        IB_PDT = IB_PDT + IOPEN*NPDT
      END DO ! End of loop over number of open orbitals
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Reorder array for determinants: conf => string order'
        CALL IWRTMA(ICTSDT,1,IDET_CONF,1,IDET_CONF)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CNFORD')
      RETURN
      END
      SUBROUTINE GETDET_FOR_CONF(IOCC,NOCOB,NOPEN,NELEC,JPDET,IPDET,
     &           IDET,I_DO_REO,IREO)
*
* The occupation IOCC of a configuration is given together with 
* NPDET prototype-determinants IPDET. Obtain determinant  corresponding 
* to prototype-determinant JPDET
* with positive values flagging alpha spin and negative values flagging 
* betaspin
*
* If I_DO_REO. ne. 0, then the proto-type det is supposed to be 
* the one with the orbitals specified by IREO, and these are 
* reordered accordingly
*
*. Jeppe Olsen, June 2011
*. Last modification; Jeppe Olsen; June 2, 2013; Lugano, reordering added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. input
      INTEGER IOCC(NOCOB),IPDET(NOPEN,*)
      INTEGER IREO(*)
*. Output
      INTEGER IDET(NELEC)
*. Scratch
      INTEGER IOP(MXPOPORB),IREOP(MXPOPORB),IPLACE(MXPOPORB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GETDET_FOR_CONF '
        WRITE(6,*) ' =========================='
        WRITE(6,*)
        WRITE(6,*) ' Input configuration '
        CALL IWRTMA(IOCC,1,NOCOB,1,NOCOB)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NOPEN, JPDET = ', NOPEN, JPDET
        WRITE(6,*) ' I_DO_REO = ', I_DO_REO
      END IF
*
      IF(I_DO_REO.EQ.1) THEN
*. The coupling of the unpaired orbitals is specified by the
*. IREO_SPCP_OB* arrays. Obtain first the open orbitals and 
*. then the order in which they are coupled
*
        NOP = 0
        DO IOB = 1, NOCOB
          IF(IOCC(IOB).GT.0) THEN
           NOP = NOP + 1
           IOP(NOP) = IOCC(IOB)
          END IF
        END DO
*
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Unpaired electrons '
         CALL IWRTMA(IOP,1,NOPEN,1,NOPEN)
        END IF
*. Obtain the order in which these orbitals should be coupled
C       ORDSTR_GEN(IINST,NELMNT,IORD,ISIGN,IOUTST,
C    &             IPLACE,IPRNT)
        CALL ORDSTR_GEN(IOP,NOPEN,IREO,ISIGN,IREOP,IPLACE,NTEST)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Order in which the open orbitals are coupled'
          CALL IWRTMA(IPLACE,1,NOPEN,1,NOPEN)
        END IF
      END IF ! Reordering is in action
*
      IEL = 0
      IOPEN = 0
      DO JJORB = 1, NOCOB
        JORB = IOCC(JJORB)
        IF(JORB.LT.0) THEN
*. Orbital JORB is double occupied in determinant
*. Alpha
          IEL = IEL + 1
          IDET(IEL) = -JORB
*. Beta
          IEL = IEL + 1
          IDET(IEL) =  JORB
        ELSE
*. Orbital JORB is single occupied, read projection from prototype det
          IOPEN = IOPEN + 1
          IIOPEN = IOPEN
          IF(I_DO_REO.EQ.1) THEN
*. Determine what orbital corresponds to IOPEN
            DO IIOP = 1, NOPEN
              IF(IPLACE(IIOP).EQ.IOPEN) IIOPEN = IIOP
            END DO
          END IF
          IEL = IEL + 1
*. Alpha spin
          IF(IPDET(IIOPEN,JPDET).EQ.1) IDET(IEL) = JORB
*. Beta  spin
          IF(IPDET(IIOPEN,JPDET).EQ.0) IDET(IEL) =-JORB
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GETDET_FOR_CONF '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Configuration in compact form: '
        CALL IWRTMA(IOCC,1,NOCOB,1,NOCOB)
        WRITE(6,*) ' Number of prototype determinant ', JPDET
        WRITE(6,*) ' Proto type determinant ', JPDET
        CALL IWRTMA(IPDET(1,JPDET),1,NOPEN,1,NOPEN)
        WRITE(6,*) ' Output determinant: '
        NEL = IEL
        CALL IWRTMA(IDET,1,NEL,1,NEL)
      END IF
*
      RETURN
      END
      SUBROUTINE ADD_C_TO_INTVEC(INTVEC,IC,NDIM)
*
* Add constant IC to elements of integer vector INTVEC
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER INTVEC(NDIM)
*
      DO I = 1, NDIM
        INTVEC(I) = INTVEC(I) + IC
      END DO
*
      RETURN
      END
      SUBROUTINE GET_SPGP_FROM_OCC(ISPGP,IGASOCC,ITP)
*
* Obtain supergroup from occupation for given type of supergroup
* Number is relative to start of supergroup
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. General input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER IGASOCC(NGAS)
*
      ISPGP = 0
      DO JSPGP = 1 , NSPGPFTP(ITP)
        JJSPGP = JSPGP + IBSPGPFTP(ITP)-1
        IFOUND = 1
        DO IGAS = 1, NGAS
          IF(NELFSPGP(IGAS,JJSPGP).NE.IGASOCC(IGAS)) IFOUND = 0
        END DO
        IF(IFOUND.EQ.1) ISPGP = JSPGP
      END DO
*
      IF(ISPGP.EQ.0) THEN
        WRITE(6,*) ' Supergroup not found '
        WRITE(6,*) ' Input type: ', ITP
        WRITE(6,*) ' Input occupation: '
        CALL IWRTMA(IGASOCC,1,NGAS,1,NGAS)
        STOP ' Supergroup not found'
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input type =', ITP
        WRITE(6,*) ' Input occupation of supergroup '
        CALL IWRTMA(IGASOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Corresponding supergroup number ', ISPGP
      END IF
*
      RETURN
      END
      SUBROUTINE OFFSETS_FOR_ALBESPGP(IBLOCK,NBLOCK,IALSPGP,IBESPGP,
     &                                  IBLOCK_FOR_SPGP)
*
* The blocks of a CI expansion is defined in IBLOCK
* Find offsets for the blocks with given alpha and beta-supergroups
* and various alpha-symmetries
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'csm.inc'
*. Input
      INTEGER IBLOCK(8,*)
*. Output: Block number for various alpha-symmetries
      INTEGER IBLOCK_FOR_SPGP(NSMST)
*
      IZERO = 0
      CALL ISETVC(IBLOCK_FOR_SPGP,IZERO,NSMST)
      NBLOCK_FOR_SPGP = 0
      DO JBLOCK = 1, NBLOCK
        IF(IBLOCK(1,JBLOCK).EQ.IALSPGP.AND.
     &     IBLOCK(2,JBLOCK).EQ.IBESPGP) THEN
*. Correct supergroups, this block belongs to the chosen ones, 
*. record offset wrt start of batch - in unpacked form
           IALSM = IBLOCK(3,JBLOCK)
           IBLOCK_FOR_SPGP(IALSM) = IBLOCK(5,JBLOCK)
        END IF
      END DO
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output from OFFSETS_FOR_ALBESPGP '
        WRITE(6,*) ' Requested alpha- and beta-supergroups =',
     &               IALSPGP, IBESPGP
        WRITE(6,*) ' Offsets for blocks with various alpha-symmetries'
        CALL IWRTMA(IBLOCK_FOR_SPGP,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
      SUBROUTINE SET_IARR_LIN_FUNC(IARR,IA,IB,NDIM)
*
* IARR(I) = IA + IB*I
*
      INCLUDE 'implicit.inc'
*. Output
      INTEGER IARR(NDIM)
*
      DO I = 1, NDIM
       IARR(I) = IA + IB*I
      END DO
*
      RETURN
      END  
      FUNCTION ISTR_GRP(IGAS,IEL)
*
* Obtain address/number of string group with IEL electron in GASpace IGAS
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
*
      IGRP = 0
      DO JGRP = IBGPSTR(IGAS), IBGPSTR(IGAS) + NGPSTR(IGAS) -1
        IF(NELFGP(JGRP).EQ.IEL) IGRP = JGRP
      END DO
*
      ISTR_GRP = IGRP
*
      IF(IGRP.EQ.0) THEN
        WRITE(6,*) ' ISTR_GRP: Required IGAS, IEL = ', IGAS, IEL
        WRITE(6,*) ' Requested group not found '
        STOP ' ISTR_GRM: Requested group not found '
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ISTR_GRP: Required IGAS, IEL = ', IGAS, IEL
        WRITE(6,*) ' Corresponding group = ', IGRP
      END IF
*
      RETURN
      END
      SUBROUTINE CSDTVCM(CSFVEC,DETVEC,SCR,IWAY,ICOPY,ISYM,ISPC,
     &           IMAXMIN_OR_GAS)
*
* Outer routine for transformation between CSF and CM forms of 
* vectors. The CM's may either by spin-combinations or determinants
*
* IWAY = 1 => CSF to CM
* IWAY = 2 => CM to CSF
*
* For GAS-expansions (where reorganization of SD's are required),
* the scratch vector SCR of length NCM_STRING is used
*
* Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cprnt.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      integer CSFVEC
      DIMENSION SCR(*)
      CALL LUCIAQENTER('CSDTV')
      NTESTL = 000
      NTEST  = MAX(NTESTL,IPRCSF) 
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' CSDTVCM speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*)
        WRITE(6,*) ' ICOPY, ISYM, ISPC,IWAY = ', 
     &               ICOPY, ISYM, ISPC, IWAY
        WRITE(6,*) ' IMAXMIN_OR_GAS = ', IMAXMIN_OR_GAS
      END IF
*
      NCMB = NCM_PER_SYM_GN(ISYM,ISPC)
      NCSF = NCSF_PER_SYM_GN(ISYM,ISPC) 
*
      IF(IMAXMIN_OR_GAS.EQ.1) THEN
*
*, The expansion is a MAXMIN space 
*
       IF(NTEST.GE.10) WRITE(6,*) ' NCMB, NCSF = ', NCMB, NCSF
*. No reordering
       IDO_REO = 0
       CALL CSDTVCB(CSFVEC,DETVEC,IWAY,dbl_mb(KDTOC),
     &      int_mb(KICONF_REO_GN(ISYM,ISPC)),NCMB,NCMB,NCSF,
     &      NCONF_PER_OPEN_GN(1,ISYM,ISPC),
     &      ICOPY,IDO_REO,NTEST)
C     CSDTVCB(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,
C    &                  NDET,NCSF,NCNFTP,
C    &                  ICOPY,IDO_REO,NTEST)
      ELSE 
*
*. The CI expansion is a standard GAS expansion, requiring reordering
*. Note that the reorder array is assumed in KSDREO_I(ISYM) irrespectively
* of space
       IDO_REO = 1
       IF(IWAY.EQ.1) THEN
*. CSF => SD transformation
        CALL COPVEC(CSFVEC,SCR,NCSF)
        CALL CSDTVCB(SCR,DETVEC,IWAY,dbl_mb(KDTOC),
     &       int_mb(KSDREO_I(ISYM)),NCCM_CONF,NCCM_STRING,NCSF,
     &       NCONF_PER_OPEN(1,ISYM),
     &       ICOPY,IDO_REO,NTEST)
       ELSE
*. SD => CSF transformation
        CALL CSDTVCB(SCR,DETVEC,IWAY,dbl_mb(KDTOC),
     &       int_mb(KSDREO_I(ISYM)),NCCM_CONF,NCCM_STRING,NCSF,
     &       NCONF_PER_OPEN(1,ISYM),
     &       ICOPY,IDO_REO,NTEST)
        CALL COPVEC(SCR,CSFVEC,NCSF)
       END IF
      END IF
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Output CSF and CM vectors from CSDTVCM' 
       WRITE(6,*) ' ======================================'
       WRITE(6,*) 
       CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF)
       CALL WRTMAT(DETVEC,1,NCCM_STRING,1,NCCM_STRING)
      END IF
*
      CALL LUCIAQEXIT('CSDTV')
      RETURN
      END 
      SUBROUTINE CSDTVCB(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,
     &                  NCM_CONF,NCM_STRING,NCSF,NCNFTP,
     &                  ICOPY,IDO_REO,NTEST)
C
C IWAY = 1 : CSF to DETERMINANT TRANSFORMATION
C IWAY = 2 : DETERMINANT TO CSF TRANSFORMATION
*
*. It is actually transformation between CSF's and combinations...
C
C ICOPY .NE. 0 : COPY OUTPUT INTO INPUT
C                SO INPUT BECOMES OUTPUT WHILE
C                OUTPUT REMAINS OUTPUT ( FOR THE MOMENT )
* IF IDO_REO = 0, no reordering of determinants are performed
* IF IDO_REO = 1, the determinants are reordered, and input vector
*                 is overwritten and must have the dimension of 
*                 the SD's
C
C
* Slight modification of CSDTVC routine, Jeppe Olsen, June 2011
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      DIMENSION DTOCMT(*),ICTSDT(*)
      DIMENSION NCNFTP(*)
*. Specific Input and output
      DIMENSION CSFVEC(*),DETVEC(*)
*. Note: Both CSFVEC and DETVEC should be able to hold DET expansion
*. and input vector is in  general destroyed
C
      LUPRI = 6
      IF( NTEST .GE. 20 ) WRITE(LUPRI,*) ' >>> CSDTVC SPEAKING <<<'
C
      IF(NTEST.GE.20) WRITE(LUPRI,*) ' NCM_STRING,NCM_CONF, NCSF ',
     &                                 NCM_STRING,NCM_CONF, NCSF
      IF(NTEST.GE.20) WRITE(LUPRI,*) ' ICOPY = ', ICOPY
*
      NTYP = MAXOP + 1
*. Minimize compiler warnings
      IOFFCS = 0
      IOFFDT = 0
      IOFFCD = 0

      IF(IWAY .EQ. 1 ) THEN
C
         IF(NTEST .GE. 40) THEN
            WRITE(LUPRI,*) ' INPUT VECTOR IN CSF BASIS '
            CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF)
         END IF
C.. CSF to DET transformation
C
C. Multiply with  expansion matrix
        DO 100 ITYP = MINOP+1,NTYP
          IDET = NPCMCNF(ITYP)
          ICSF = NPCSCNF(ITYP)
          ICNF = NCNFTP(ITYP)
          IF(NTEST.GE.100) THEN
            WRITE(6,'(a,4i8)') 
     &      ' ITYP, IDET, ICSF, ICNF = ', ITYP,IDET,ICSF,ICNF
          END IF
*
          IF(ITYP .EQ. MINOP+1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNFTP(ITYP-1)*NPCSCNF(ITYP-1)
            IOFFDT = IOFFDT+NCNFTP(ITYP-1)*NPCMCNF(ITYP-1)
            IOFFCD = IOFFCD + NPCMCNF(ITYP-1)*NPCSCNF(ITYP-1)
          END IF
          CALL MATML4(DETVEC(IOFFDT),DTOCMT(IOFFCD),CSFVEC(IOFFCS),
     &                IDET,ICNF,IDET,ICSF,ICSF,ICNF,0)
  100   CONTINUE
C. Change from csf to string ordering with sign changes
        IF(IDO_REO.EQ.1) THEN
          ZERO = 0.0D0
          CALL SETVEC(CSFVEC,ZERO,NCM_STRING)
          CALL SCACSF(NCM_CONF,CSFVEC,DETVEC,ICTSDT)
          IF(NTEST.GE.1000) THEN
             WRITE(6,*) 'ICTSDT array '
             CALL IWRTMA(ICTSDT,1,NCM_CONF,1,NCM_CONF)
          END IF
          CALL COPVEC(CSFVEC,DETVEC,NCM_STRING)
        END IF
        IF(NTEST .GE. 40) THEN
           WRITE(LUPRI,*) ' OUTPUT VECTOR IN DET BASIS '
           CALL WRTMAT(DETVEC,1,NCM_STRING,1,NCM_STRING,0)
        END IF
      ELSE
         IF( NTEST .GE. 40 ) THEN
            WRITE(LUPRI,*) ' INPUT VECTOR IN DETERMINANT BASIS '
            CALL WRTMAT(DETVEC,1,NCM_STRING,1,NCM_STRING,0)
         END IF
C
C.. Determinant to csf transformation
C
C. To CSF ordering
        IF(IDO_REO.EQ.1) THEN
          CALL GATCSF(NCM_CONF,DETVEC,CSFVEC,ICTSDT)
          CALL COPVEC(CSFVEC,DETVEC,NCM_CONF)
        END IF
C. Multiply with CIND expansion matrix
        DO 200 ITYP = 1,NTYP
          IDET = NPCMCNF(ITYP)
          ICSF = NPCSCNF(ITYP)
          ICNF = NCNFTP(ITYP)
          IF(ITYP .EQ. 1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNFTP(ITYP-1)*NPCSCNF(ITYP-1)
            IOFFDT = IOFFDT+NCNFTP(ITYP-1)*NPCMCNF(ITYP-1)
            IOFFCD = IOFFCD + NPCMCNF(ITYP-1)*NPCSCNF(ITYP-1)
          END IF
          CALL MATML4(CSFVEC(IOFFCS),DTOCMT(IOFFCD),DETVEC(IOFFDT),
     &                ICSF,ICNF,IDET,ICSF,IDET,ICNF,1)
  200   CONTINUE
        IF( ICOPY .NE. 0 ) CALL COPVEC(CSFVEC,DETVEC,NCSF)
        IF( NTEST .GE. 15 ) THEN
          WRITE(LUPRI,*) ' OUTPUT VECTOR IN CSF BASIS '
          CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF)
        END IF
       END IF
C
      RETURN
      END
      SUBROUTINE GATCSF(NDET,ADET,BCSF,IORD)
C
C Gather determinant vector in CSF order
C from determinant vector in string order,
C with sign changes caused by switch from
C string order to configuration order
C
      INCLUDE 'implicit.inc'
*
      DIMENSION ADET(NDET), BCSF(NDET)
      INTEGER   IORD(NDET)
      DO I = 1,NDET
        IF(IORD(I).GT.0) THEN
          BCSF(I) = ADET(IORD(I))
        ELSE IF (IORD(I).LT.0) THEN
          BCSF(I) = -ADET(-IORD(I))
        END IF
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from GATCSF: '
        WRITE(6,*) ' Scatter array '
        CALL IWRTMA(IORD,1,NDET,1,NDET)
      END IF
*
      RETURN
      END
      SUBROUTINE SCACSF(NDET,ADET,BCSF,IORD)
C
C Scatter determinant vector in CSF order
C to determinant vector in string order,
C with sign changes caused by switch from
C string order to configuration order
C
      INCLUDE 'implicit.inc'
*
      DIMENSION ADET(NDET), BCSF(NDET)
      INTEGER   IORD(NDET)
      DO I = 1,NDET
        IF (IORD(I) .LT. 0) THEN 
          ADET(-IORD(I)) = -BCSF(I)
        ELSE IF (IORD(I).GT.0) THEN
          ADET(IORD(I)) = BCSF(I)
        END IF
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SCACSF: '
        WRITE(6,*) ' Scatter array '
        CALL IWRTMA(IORD,1,NDET,1,NDET)
      END IF
*
      RETURN
      END
      FUNCTION ILEX_FOR_CONF_G(ICONF,NOCC_ORB,ICONF_SPC,IDOREO)
*
* Obtain Lexical order in CONFSPACE, perhaps reordered, 
* for configuration ICONF
*
* The configuration ICONF is given in packed form
*
*. Jeppe Olsen, July 2011, generalization to multiple spaces
*  
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'spinfo.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific Input
      INTEGER ICONF(NOCC_ORB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Address is requested for configuration '
       CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
       WRITE(6,*) ' In configuration space ', ICONF_SPC
      END IF
*. Symmetry
C     ISYM_CONF(IOCC,NORBL,IORB_OFF)
      IDUM = 0
C     ISYM_CONF_G(IOCC,NOCCL,NORBL,IORB_OFF,ICONF_FORM)
      ISYM = ISYM_CONF_G(ICONF,NOCC_ORB,IDUM,IB_ORB_CONF,2)      
*
      IF(NTEST.GE.1000) WRITE(6,*) ' IB_ORB_CONF, ISYM = ', 
     &                  IB_ORB_CONF, ISYM
*
      ILEX = ILEX_FOR_CONF(ICONF,NOCC_ORB,N_ORB_CONF,N_EL_CONF,
     &       int_mb(KZ_GN(ICONF_SPC)),
     &       IDOREO,int_mb(KICONF_REO_GN(ISYM,ICONF_SPC)))
*
C?    WRITE(6,*) ' Test: first three elements of KICONF_REO..'
C?    CALL IWRTMA(WORK(KICONF_REO_GN(ISYM,ICONF_SPC)),1,3,1,3)
*
C     ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
*
      ILEX_FOR_CONF_G = ILEX
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Configuration address = ', ILEX
      END IF
*
      RETURN
      END
      SUBROUTINE OFFSETS_CSF_FOR_NOPEN(IB_OPEN_CSF,NCONF_OPEN,MAXOP,
     &                                 NPCSCNF)
*
* Generate offsets to start of CSF's with a given number of 
* open orbitals
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*
      INTEGER NPCSCNF(MAXOP+1)
*.Specific input: Number of configurations per number of open orbitals
      INTEGER NCONF_OPEN(MAXOP+1)
*. Output
      INTEGER IB_OPEN_CSF(MAXOP+1)
*
      NTEST = 00
*
      IB_OPEN_CSF(1) = 1
      DO IOPEN = 0, MAXOP-1
       ITP = IOPEN + 1
       NCSF_PER_CONF = NPCSCNF(ITP)
       NCSF_PER_OPEN = NCSF_PER_CONF*NCONF_OPEN(ITP)
       IB_OPEN_CSF(ITP+1) = IB_OPEN_CSF(ITP) + NCSF_PER_OPEN
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) 
     & ' Offsets to CSFs with given number of open orbitals'
       CALL IWRTMA(IB_OPEN_CSF,1,MAXOP+1,1,MAXOP+1)
      END IF
*
      RETURN
      END
      SUBROUTINE GETDETS_FOR_CONF(IOCC,NOCOB,NOPEN,NELEC,NPDET,IPDET,
     &           IDET)
*
* The occupation IOCC_CONF of a configuration is given together with 
* NPDET prototype-determinants IPDET. Obtain the determinants corresponding 
* to the prototype-determinants
* with positive values flagging alpha spin and negative values flagging 
* betaspin
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC_CONF(NOCOB),IPDET(NOPEN,NPDET)
*. Output
      INTEGER IOCC_DET(NELEC,NPDET)
*
      NTEST = 000
*
      DO IDET = 1, NPDET
        IEL = 0
        IOPEN = 0
        DO JJORB = 1, NOCOB
         JORB = IOCC_CONF(JJORB)
         IF(JORB.LT.0) THEN
*. Orbital JORB is double occupied in determinant
*. Alpha
           IEL = IEL + 1
           IOCC_DET(IEL,IDET) = -JORB
*. Beta
           IEL = IEL + 1
           IOCC_DET(IEL,IDET) =  JORB
         ELSE
*. Orbital JORB is single occupied, read projection from prototype det
           IOPEN = IOPEN + 1
           IEL = IEL + 1
*. Alpha spin
           IF(IPDET(IOPEN,IDET).EQ.1) IOCC_DET(IEL,IDET) = JORB
*. Beta  spin
           IF(IPDET(IOPEN,IDET).EQ.0) IOCC_DET(IEL,IDET) =-JORB
         END IF
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output from GETDETS_FOR_CONF '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Configuration in compact form: '
        CALL IWRTMA(IOCC_CONF,1,NOCOB,1,NOCOB)
        WRITE(6,*) ' Determinants: det(column), electron(row)'
        CALL IWRTMA(IOCC_DET,NELEC,NPDET,NELEC,NPDET)
      END IF
*
      RETURN
      END
      SUBROUTINE MXMNOC_GAS(MINEL_ORB,MAXEL_ORB,NGAS,NOBPT,
     &                  MINEL_GAS,MAXEL_GAS,
     &                  NTESTG)
*
* Construct accumulated MAX and MIN arrays for orbitals 
* for given accumulated  occupations of electrons for GAS-paces
*
* 
* Dec 2011 - Jeppe Olsen - programmers deja vu: I have written this routine 
*                          before
*       
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Output
      DIMENSION  MINEL_ORB(*),MAXEL_ORB(*)
*. Input
      INTEGER NOBPT(NGAS),MINEL_GAS(NGAS),MAXEL_GAS(NGAS)
*
      INCLUDE 'mxpdim.inc'
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==========='
        WRITE(6,*) ' MXMNOC_GAS '
        WRITE(6,*) ' ==========='
        WRITE(6,*)
        WRITE(6,*) ' NOBPT : '
        CALL IWRTMA(NOBPT,1,NGAS,1,NGAS)
      END IF
*. Minimize compiler warnings
      MINEL = 0
      MAXEL = 0
      IB_ORB = 0
      LORB = 0
*
      NELEC = MAXEL_GAS(NGAS)
      DO IGAS = 1, NGAS
*. Max and min number of electrons in this space
       IF(IGAS.EQ.1) THEN
         MINEL = MINEL_GAS(1)
         MAXEL = MAXEL_GAS(1)
         IB_ORB = 1
         MIN_START = 0
         MAX_START = 0
         MIN_END = MINEL_GAS(1)
         MAX_END = MAXEL_GAS(1)
       ELSE
         MINEL = MAX(0,MINEL_GAS(IGAS)-MAXEL_GAS(IGAS-1))
         MAXEL = MIN(2*NOBPT(IGAS),MAXEL_GAS(IGAS)-MINEL_GAS(IGAS-1))
         IB_ORB = IB_ORB + NOBPT(IGAS-1)
         MIN_START = MINEL_GAS(IGAS-1)
         MAX_START = MAXEL_GAS(IGAS-1)
         MIN_END = MINEL_GAS(IGAS)
         MAX_END = MAXEL_GAS(IGAS)
       END IF
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' IGAS, IB_ORB, MINEL, MAXEL = ',
     &                IGAS, IB_ORB, MINEL, MAXEL
       END IF
       LORB = NOBPT(IGAS)
       DO IORB = IB_ORB, IB_ORB-1+LORB
         MINEL_ORB(IORB) = 
     &   MAX(MIN_START,MIN_END-2*(LORB-(IORB-IB_ORB+1)))
         MAXEL_ORB(IORB) = 
     &   MIN(MAX_START+MAXEL, MAX_START + 2*(IORB-IB_ORB+1),NELEC,
     &        MAX_END)
       END DO !orbitals
      END DO! gaspaces
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*)  'Output from MXMNOC_GAS '
        WRITE(6,*) ' ======================='
        WRITE(6,*) 
        WRITE(6,*) 
        WRITE(6,*) 
        NORB = IB_ORB + LORB - 1
        CALL WRT_MINMAX_OCC(MINEL_ORB,MAXEL_ORB,NORB)
      END IF
*
      RETURN
      END
      SUBROUTINE CSDIAG(CSFDIA,DETDIA,NCNFTP,MAXOP,ISM,
     &           ICTSDT,NPCMCNF,NPCSCNF,IPRCSF,
     &           ICNFBAT,NOCCLS_ACT,IOCCLS_ACT,
     &           LBLOCK,LUDIA_DET,LUDIA_CSF)
*
* Obtain averaged CI diagonal in CSF basis from 
* CI diagonal in SD basis
*
* ICNFBAT = 1 => ICISTR = 1 storage mode, i.e. diagonal is in core
* ICNFBAT = 2 => ICISTR = 2 storage mode, i.e. diagonal written to disc
* are in core. It is assumed that each occupation block defines a 
* batch.
*
*. Jeppe Olsen, Jan 2012, form CSDIAG in lucas
*               Febr. 2012: Updated with ICNFBAT option
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      integer DETDIA
CNW   DIMENSION DETDIA(*)
      DIMENSION NCNFTP(MAXOP+1),NPCMCNF(MAXOP+1),NPCSCNF(MAXOP+1)
      INTEGER ICTSDT(*)
*. (If ICNFBAT = 2)
      INTEGER IOCCLS_ACT(NOCCLS_ACT)
      INTEGER LBLOCK(NOCCLS_ACT)
*. Output
      integer CSFDIA
CNW   DIMENSION CSFDIA(*)
  
      NTEST = 000
      NTEST = MAX(NTEST,IPRCSF)
      IF(NTEST.GE.10) THEN
       WRITE(6,*)
       WRITE(6,*) ' ================== '
       WRITE(6,*) ' Wellcome to CSDIAG '
       WRITE(6,*) ' ================== '
       WRITE(6,*)
       WRITE(6,*) ' ISM, ICNFBAT = ', ISM, ICNFBAT 
       WRITE(6,*) ' LUDIA_DET, LUDIA_CSF = ', 
     &              LUDIA_DET, LUDIA_CSF
      END IF
*
      ICSOFF = 1
      IDTOFF = 1
      JCNABS = 0
*
      IF(ICNFBAT.EQ.1) THEN
       DO ITYP = 1, MAXOP + 1
         IDET = NPCMCNF(ITYP)
         ICSF = NPCSCNF(ITYP)
         ICNF = NCNFTP(ITYP)
         IF(NTEST.GE.1000) THEN
           WRITE(6,'(A,4I7)') ' ITYP, IDET, ICSF,ICNF =', 
     &                          ITYP, IDET, ICSF, ICNF
         END IF
         DO JCNF = 1, ICNF
           JCNABS = JCNABS + 1
           EAVER = 0.0D0
           DO JDET = 1, IDET
             IADD = ABS(ICTSDT(IDTOFF-1+JDET))
             call ga_get(DETDIA,IADD,IADD,1,1,EVAL,1)
             EAVER = EAVER + EVAL
CNW          EAVER = EAVER +DETDIA(ABS(ICTSDT(IDTOFF-1+JDET)))
           END DO
           IF( IDET .NE. 0 )EAVER = EAVER/IDET
           call ga_fill_patch(CSFDIA,ICSOFF,ICSOFF+ICSF-1,1,1,EAVER)
CNW        CALL SETVEC(CSFDIA(ICSOFF),EAVER,ICSF)
           ICSOFF = ICSOFF + ICSF
           IDTOFF = IDTOFF + IDET
         END DO !loop over confs
       END DO !loop over ITYP
*
      ELSE
*. Process each occupation classes individually
       LBLK = -1
       CALL REWINO(LUDIA_DET)
       CALL REWINO(LUDIA_CSF)
       DO IIOCLS = 1, NOCCLS_ACT
*. Read in the diagonal in det basis for this batch
        NBLK = LBLOCK(IIOCLS)
CBERT   CALL FRMDSCN(DETDIA,NBLK,LBLK,LUDIA_DET)
C       FRMDSCN(VEC,NREC,LBLK,LU)
        IOCCLS = IOCCLS_ACT(IIOCLS)
*. Obtain configurations for this occupation class 
C           GEN_CNF_INFO_FOR_OCCLS(IOCCLS_NUM,IDOSDREO,ISYM)
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCCLS,1,ISM)
        ICSOFF = 1
        IDTOFF = 1
        JCNABS = 0
        DO ITYP = 1, MAXOP + 1
         IDET = NPCMCNF(ITYP)
         ICSF = NPCSCNF(ITYP)
         ICNF = NCNFTP(ITYP)
         IF(NTEST.GE.1000) THEN
           WRITE(6,'(A,4I7)') ' ITYP, IDET, ICSF,ICNF =', 
     &                          ITYP, IDET, ICSF, ICNF
         END IF
         DO JCNF = 1, ICNF
           JCNABS = JCNABS + 1
           EAVER = 0.0D0
           DO JDET = 1, IDET
CBERT        EAVER = EAVER +DETDIA(ABS(ICTSDT(IDTOFF-1+JDET)))
           END DO
           IF( IDET .NE. 0 )EAVER = EAVER/IDET
CBERT      CALL SETVEC(CSFDIA(ICSOFF),EAVER,ICSF)
           ICSOFF = ICSOFF + ICSF
           IDTOFF = IDTOFF + IDET
         END DO !loop over confs
        END DO !loop over ITYP
*. Write to Disc
        NCSF_OCC = ICSOFF - 1
        CALL ITODS(NCSF_OCC,1,-1,LUDIA_CSF)
CBERT   CALL TODSC(CSFDIA,NCSF_OCC,-1,LUDIA_CSF)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Averaged CSF  diagonal for occls '
CBERT     CALL WRTMAT(CSFDIA,1,NCSF_OCC,1,NCSF_OCC)
        END IF
       END DO! over active occupation classes
*This is the end...
       CALL ITODS(-1,1,-1,LUDIA_CSF)
       IF(NTEST.GE.1000) WRITE(6,*) ' EOF written '
      END IF ! ICNFBAT switch
*
      IF( NTEST .GE. 1000 ) THEN
        WRITE(6,*) 'CI diagonal in DET and CSF basis '
        IF(ICNFBAT.EQ.1) THEN
          NCSTOT = ICSOFF-1
          NDTTOT = IDTOFF-1
CNW       CALL WRTMAT(DETDIA,1,NDTTOT,1,NDTTOT)
          call ga_print_patch(DETDIA,1,NCSTOT,1,NCSTOT,1)
          WRITE(6,*)
CNW       CALL WRTMAT(CSFDIA,1,NCSTOT,1,NCSTOT) 
          call ga_print_patch(CSFDIA,1,NCSTOT,1,NCSTOT,1)
        ELSE
CNW       CALL WRTVCD(DETDIA,LUDIA_DET,1,LBLK)
          WRITE(6,*)
CNW       CALL WRTVCD(DETDIA,LUDIA_CSF,1,LBLK)
C               WRTVCD(SEGMNT,LU,IREW,LBLK)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE NXT_OCCLS_IN_MINMAX_SPC(IOCC,IGSOCC_MIN,IGSOCC_MAX,
     &           NGAS,NOBPT,INI,NONEW)
*. Obtain next occupation class in MINMAX space
*
*. Jeppe Olsen, Jan. 2011
*
      INCLUDE 'implicit.inc' 
* Input
      INTEGER IGSOCC_MIN(NGAS),IGSOCC_MAX(NGAS),NOBPT(NGAS)
*. Output
      INTEGER IOCC(NGAS)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from NXT_OCCLS_IN_MINMAX_SPC'
        WRITE(6,*) ' ==================================='
C?      WRITE(6,*) ' INI, NGAS = ', INI, NGAS
      END IF
*
      IF(INI.EQ.1) THEN
*. Initialize
        DO  IGAS = 1, NGAS
          IOCC(IGAS) = IGSOCC_MIN(IGAS)
        END DO
        NONEW = 0
      ELSE
*. Next occupation class.
*. Find first GASpace, where the number of electrons may be increased  
        KGAS = 0
        DO IGAS = 1, NGAS-1
          IF((IOCC(IGAS).LT.IGSOCC_MAX(IGAS)). AND.
     &       (IOCC(IGAS).LT.IOCC(IGAS+1)) ) THEN
            KGAS = IGAS
            GOTO 101
          END IF
        END DO
  101   CONTINUE
        IF(KGAS.EQ.0) THEN
          NONEW = 1
        ELSE 
          IOCC(KGAS) = IOCC(KGAS) + 1
          NONEW = 0
*. Obtain smallest possible occupations in 1 - (KGAS-1) i.e. occupy the 
*. highest Gaspaces first...
*. Min number of electrons in 1 - (KGAS-1)
          DO IGAS = KGAS,2,-1
*. Max number of electrons in IGAS => occupation of IGAS-1
            NELEC_I_MAX = 
     &      MIN(2*NOBPT(IGAS),IOCC(IGAS)-IGSOCC_MIN(IGAS-1))
            IOCC(IGAS-1) = IOCC(IGAS)-NELEC_I_MAX
          END DO
*
CE        DO IGAS = 1, KGAS - 1
CE          IF(IGAS.EQ.1) THEN
CE           NELEC_I_MAX = IGSOCC_MAX(1)
CE          ELSE
CE           NELEC_I_MAX = 
CE   &       MIN(2*NOBPT(IGAS),IGSOCC_MAX(IGAS)-IOCC(IGAS-1))
CE          END IF
CE          IF(IGAS.EQ.1) THEN
CE            IOCC(1) = MIN(IGSOCC_MIN(1)+NADD,NELEC_I_MAX)
CE          ELSE
CE            IOCC(IGAS) = MIN(IGSOCC_MIN(IGAS)+NADD,IGSOCC_MAX(IGAS),
CE   &                         IOCC(IGAS-1) + NELEC_I_MAX)   
CE            IOCC(IGAS) = MAX(IOCC(IGAS),IOCC(IGAS-1))
CE          END IF
CE          NADD = NADD0 -  (IOCC(IGAS)-IGSOCC_MIN(IGAS))
CE        END DO
        END IF !NONEW
      END IF !Initialize
*
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
          WRITE(6,*) ' No new occupation classes'
        ELSE 
          WRITE(6,*) ' Accumulated occupation for next occ. class '
          CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE MAX_NOPEN_NOCCLS_FOR_CISPAC
     &(IGSOCC_MIN, IGSOCC_MAX,NOPEN_MAX,NOCCLS)
*
* Find max number of open orbitals and number of occupation
* classes for MINMAX space defined by IGSOCC
*
*. Jeppe Olsen, Jan. 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER IGSOCC_MIN(MXPNGAS), IGSOCC_MAX(MXPNGAS)
*. Local scratch
      DIMENSION IOCC_A(MXPNGAS),IOCC(MXPNGAS)

*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) '  MAX_NOPEN_NOCCLS... in action'
         WRITE(6,*) ' ======================'
         WRITE(6,*) ' NGAS ', NGAS
         WRITE(6,*) ' IGSOCC_MIN, IGSOCC_MAX: '
         CALL WRT_MINMAX_GASOCC(IGSOCC_MIN,IGSOCC_MAX,NGAS)
      END IF
*   
      MAXLOOP = 10000000
cVOG  MAXLOOP = 10000
      NLOOP = 0
      NOCCLS = 0
      NOPEN_MAX = 0
*. Loop over occupation classes
      INI  = 1
 1000 CONTINUE
C       NXT_OCCLS_IN_MINMAX_SPC(IOCC,IGSOCC_MIN,IGSOCC_MAX,
C    &             NGAS,NOBPT,INI,NONEW)
        CALL NXT_OCCLS_IN_MINMAX_SPC(IOCC_A,IGSOCC_MIN,IGSOCC_MAX,NGAS,
     &       NOBPT,INI,NONEW)
        INI = 0
        IF(NONEW.EQ.0) THEN
         NLOOP =  NLOOP + 1
*. Change occupation class from accumulated to actual occupation
         CALL REF_OCCLS(IOCC_A,IOCC,1,NGAS)
*. If requested, check for occupations in the compound GAS space
         IM_IN = 1
         IF(I_CHECK_ENSGS.EQ.1)
     &   CALL  CHECK_IS_OCC_IN_ENGSOCC(IOCC,-1,IM_IN)
C              CHECK_IS_OCC_IN_ENGSOCC(IGSOCCL,ISPC,IM_IN)
         IF(IM_IN.EQ.1) THEN
          NOCCLS = NOCCLS + 1
*. Largest number of open orbitals
          MAXOP_L = 0
          DO IGAS = 1, NGAS
            NEL = IOCC(IGAS)
            NORB = NOBPT(IGAS)
C?          WRITE(6,*) ' IGAS, NEL, NORB = ', IGAS, NEL, NORB
            MAXOP_IGAS = MIN(NEL,2*NORB-NEL)
            MAXOP_L = MAXOP_L + MAXOP_IGAS
          END DO
          NOPEN_MAX = MAX(NOPEN_MAX,MAXOP_L)
         END IF! IM_IN .eq. 1 
        END IF! NO_NEW = 0
      IF(NONEW.EQ.0.AND.NLOOP.LT.MAXLOOP) GOTO 1000
      IF(NLOOP.EQ.MAXLOOP) THEN
       if (ga_nodeid().eq.0) WRITE(6,*) ' Forced MAXLOOP exit in MAXOP.'
       if (ga_nodeid().eq.0)
     &    WRITE(6,*) ' Bug or test to be removed, Jeppe!!! '
cVOG    STOP ' Forced MAXLOOP exit in MAXOP... '
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Number of occupation classes ', NOCCLS
        WRITE(6,*) ' Largest number of unpaired electrons ', NOPEN_MAX
      END IF
*
      RETURN
      END
      SUBROUTINE REF_OCCLS(IOCC_A,IOCC,IWAY,NGAS)
*. Reform between accumulated and actual occupation of occupation class
*
*. IWAY = 1: Accumulated to actual
*. IWAY = 2: Actual to accumulated
*. Jeppe Olsen, Jan. 2011
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER IOCC_A(NGAS), IOCC(NGAS)
*. 
      IF(IWAY.EQ.1) THEN
        IOCC(1) = IOCC_A(1)
        DO IGAS = 2, NGAS
          IOCC(IGAS) = IOCC_A(IGAS)-IOCC_A(IGAS-1)
        END DO
      ELSE IF (IWAY.EQ.2) THEN
        IOCC_A(1) = IOCC(1)
        DO IGAS = 2, NGAS
          IOCC_A(IGAS) = IOCC(IGAS) + IOCC_A(IGAS-1)
        END DO
      ELSE 
        WRITE(6,*) ' REF_OCCLS, illegal IWAY = ', IWAY
        STOP       ' REF_OCCLS, illegal IWAY = '
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from REF_OCCLS '
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Accumulated to actual '
        ELSE
          WRITE(6,*) ' Actual to accumulated '
        END IF
        WRITE(6,*) ' Accumulated occupation '
        CALL IWRTMA(IOCC_A,1,NGAS,1,NGAS)
        WRITE(6,*) ' Actual occupation '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_MINMAX_GASOCC(IOCC_MIN,IOCC_MAX,NGAS)
*
* Write min and max accumulated occupation arrays for GASpaces
*
*. Jeppe Olsen, June 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC_MIN(NGAS),IOCC_MAX(NGAS)
*
      WRITE(6,*) ' Min and Max accumulated occupations: '
      WRITE(6,*)
      WRITE(6,*) ' GASpace Min. occ Max. occ '
      WRITE(6,*) ' =========================='
      DO IGAS = 1, NGAS
        WRITE(6,'(3X,I4,2(4X,I4))')
     &  IGAS, IOCC_MIN(IGAS), IOCC_MAX(IGAS)
      END DO
*
      RETURN
      END
      SUBROUTINE IADCONST(IVEC,IADD, NDIM)
*
* Add IADD to integer array IVEC
*
      INCLUDE 'implicit.inc'
      INTEGER IVEC(NDIM)
*
      DO I = 1, NDIM
        IVEC(I) = IVEC(I) + IADD
      END DO
*
      RETURN
      END
      SUBROUTINE GET_CSF_H_PRECOND(NCONF_FOR_OPEN,ICONF_OCC,H0,
     &           ECORE,LUDIA,NOCCLS_SPC,IOCCLS_SPC,ISYM)
*
* Obtain preconditioner H0 for the Hamiltonian matrix in the CSF basis
*
* IH0_CSF = 2: Diagonal of CSF's
*         = 3: Block diagonal with block consisting of CSF's belonging to given
*            configurations
* (IH0_CSF = 1 is reserved for something else)
*
* Storage mode for configurations are determined by ICNFBAT:
*
* ICNFBAT = 1 => All info is in core and diagonal is returned in H0 and 
*                not written to disc
*         = 2 => Info is constructed for each occ class and 
*                diagonal is written to disc with each occlass as a record
*
*. Jeppe Olsen, Jan. 2012, some speed up in June 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'gasstr.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input (or scratch)
      INTEGER ICONF_OCC(*), NCONF_FOR_OPEN(*)
*  Input for ICNFBAT = 2
      INTEGER IOCCLS_SPC(NOCCLS_SPC)
*. Output
      DIMENSION H0(*)
* 
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CSFH0 ')
      CALL LUCIAQENTER('CSFH0 ')
*
*. Scratch space for calculations
*
*. Length: INTEGER: (NDET_C + NDET_S)*N_EL_CONF + NDET_C + 6*NORB *. (NELEC .le. 2*NORB)
*          REAL   : 2 NDETL x NDETR matrix (could be reduced)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from GET_CSF_H_PRECOND'
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' IH0_CSF, ICNFBAT, ISYM = ', IH0_CSF, ICNFBAT, ISYM
        WRITE(6,*) ' MAXOP = ', MAXOP
      END IF
        
*. largest number of protype dets: that of maxop
      NPDT_MAX = NPCMCNF(MAXOP+1)
      NPCS_MAX = NPCSCNF(MAXOP+1)
*
      IATP = 1
      IBTP = 2
      NEL = NELFTP(IATP) + NELFTP(2)
*
      IF(NTEST.GE.10) 
     &WRITE(6,*) ' NPDT_NAX, NPCS_MAX = ', NPDT_MAX, NPCS_MAX
      LISCR = 2*NPDT_MAX*NEL + NPDT_MAX + 6*NACOB
      LRSCR = 2*NPDT_MAX**2
      CALL MEMMAN(KLISCR,LISCR,'ADDL  ',1,'IS_CHC')  !done
      CALL MEMMAN(KLRSCR,LRSCR,'ADDL  ',2,'RS_CHC')  !done
      LCSFHCSF = NPCS_MAX**2
      CALL MEMMAN(KLCSFHCSF,LCSFHCSF,'ADDL  ',2,'CSHCS ')  !done
*
      CALL MEMMAN(KLJ,NTOOB**2,'ADDL  ',2,'IIJJ  ')  !done
      CALL MEMMAN(KLK,NTOOB**2,'ADDL  ',2,'IJJI  ')  !done
*
*. Obtain J(I,J) = (II!JJ), K(I,J) = (IJ!JI)
*
      CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,XDUM,IREOTS)
*
      IF(ICNFBAT.EQ.1) THEN 
*
* All Conformation is in core, calculate diagonal and save in core
* 
        IB_OCC = 1
        IB_H0 = 1
        NCSF_TOT = 0
        IF(IH0_CSF.EQ.2) THEN
         IONLY_DIAG = 1
         ISYMG = 1
        ELSE
         IONLY_DIAG = 0
         ISYMG = 0
        END IF
*
C?      IONLY_DIAG = 0
C?      WRITE(6,*) ' IONLY_DIAG set to 0 '
        DO IOPEN = 0, MAXOP
          ITYP =  IOPEN + 1
          IORB = (NEL+IOPEN)/2
          NCONF = NCONF_FOR_OPEN(IOPEN+1)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' IOPEN, NCONF = ', IOPEN, NCONF
          NCSF = NPCSCNF(IOPEN+1)
          DO  ICONF = 1, NCONF
            IF(NTEST.GE.1000) 
     &      WRITE(6,*) ' Infor for CONF = ', ICONF
* Construct Hamiltonian matrix over CSF's in conf
C                CNHCN_CSF_BLK(ICNL,IOPL,ICNR,IOPR,CNHCNM,IADOB,
C    &                         IPRODT,DTOC,I12OP,ISCR,SCR,ECORE)
            CALL CNHCN_CSF_BLK(ICONF_OCC(IB_OCC),IOPEN,
     &      ICONF_OCC(IB_OCC),IOPEN,dbl_mb(KLCSFHCSF),NINOB,
     &      int_mb(KDFTP),dbl_mb(KDTOC),2,int_mb(KLISCR),
     &      dbl_mb(KLRSCR),ECORE,IONLY_DIAG,
     &      ISYMG,dbl_mb(KLJ),dbl_mb(KLK))
*. Transfer block
            IF(IH0_CSF.EQ.2) THEN
*.Extract diagonal
              CALL COPDIA(dbl_mb(KLCSFHCSF),H0(IB_H0),NCSF,0)
C                  COPDIA(A,VEC,NDIM,IPACK)
              IB_H0 = IB_H0 + NCSF
            ELSE
*. Keep complete block
              CALL COPVEC(dbl_mb(KLCSFHCSF),H0(IB_H0),NCSF**2)
              IB_H0 = IB_H0 + NCSF**2
            END IF
            IB_OCC = IB_OCC + IORB
            NCSF_TOT = NCSF_TOT + NCSF
          END DO! loop over confs for given IOPEN
        END DO! IOPEN
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 
          WRITE(6,*) ' Preconditioner '
          IF(IH0_CSF.EQ.2) THEN
            CALL WRTMAT(H0,1,NCSF_TOT,1,NCSF_TOT)
          ELSE
            IB_H0 = 1
            DO IOPEN = 0, MAXOP
             ITYP =  IOPEN + 1
             NCONF = NCONF_FOR_OPEN(IOPEN+1)
             NCSF = NPCSCNF(IOPEN+1)
             DO  ICONF = 1, NCONF
              CALL WRTMAT(H0(IB_H0),NCSF,NCSF,NCSF,NCSF)
              IB_H0 = IB_H0 + NCSF**2
             END DO! loop over confs for given IOPEN
            END DO! IOPEN
          END IF ! switch IH0_CSF
        END IF ! NTEST
*
      ELSE IF(ICNFBAT.EQ.2) THEN 
*
* Conformation is constructed separately for each occupation class. 
* Preconditioner is stored in records over occupationclasses
* 
       IF(IH0_CSF.EQ.2) THEN
        IONLY_DIAG = 1
         ISYMG = 1
       ELSE
        IONLY_DIAG = 0
         ISYMG = 0
       END IF
       CALL REWINO(LUDIA)
       IF(NTEST.GE.10) WRITE(6,*) ' NOCCLS_SPC = ', NOCCLS_SPC
       DO IIOCLS = 1, NOCCLS_SPC
        IB_OCC = 1
        IB_H0 = 1
        IDO_REO = 1
        IB_BLK = 1
        NCSF_TOT = 0
        IOCCLS = IOCCLS_SPC(IIOCLS)
        IF(NTEST.GE.10) WRITE(6,*) ' Output from IOCCLS = ', IOCCLS
*. Generate Conformation (only configurations are needed)  
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCCLS,0,ISYM)
        NCSF_OCCLS = IELSUM(NCS_FOR_OC_OP_ACT,MAXOP+1)
        NCM_OCCLS = IELSUM(NCM_FOR_OC_OP_ACT,MAXOP+1)
*
        IF(NTEST.GE.10) THEN
           WRITE(6,*) ' NCSF_OCCLS, NCM_OCCLS = ',
     &                  NCSF_OCCLS, NCM_OCCLS
        END IF
*
        DO IOPEN = 0, MAXOP
          ITYP =  IOPEN + 1
          IORB = (NEL+IOPEN)/2
          NCONF = NCONF_FOR_OPEN(IOPEN+1)
          NCSF = NPCSCNF(IOPEN+1)
          DO  ICONF = 1, NCONF
* Construct Hamiltonian matrix over CSF's in conf
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Next conf: '
              CALL IWRTMA(ICONF_OCC(IB_OCC),1,IORB,1,IORB)
            END IF
            CALL CNHCN_CSF_BLK(ICONF_OCC(IB_OCC),IOPEN,
     &      ICONF_OCC(IB_OCC),IOPEN,dbl_mb(KLCSFHCSF),NINOB,
     &      int_mb(KDFTP),
     &      dbl_mb(KDTOC),2,int_mb(KLISCR),dbl_mb(KLRSCR),ECORE,
     &      IONLY_DIAG,ISYMG,dbl_mb(KLJ),dbl_mb(KLK))
C              CNHCN_CSF_BLK(ICNL,IOPL,ICNR,IOPR,CNHCNM,IADOB,
C    &           IPRODT,DTOC,I12OP,ISCR,SCR,ECORE,IONLY_DIAG,ISYMG,
C    &           RJ, RK)
*
*
*. Transfer block
            IF(IH0_CSF.EQ.2) THEN
*.Extract diagonal
              CALL COPDIA(dbl_mb(KLCSFHCSF),H0(IB_H0),NCSF,0)
C                  COPDIA(A,VEC,NDIM,IPACK)
C?            WRITE(6,*) ' Added elements, IB_H0, NCSF = ', IB_H0, NCSF
C?            CALL WRTMAT(H0(IB_H0),1,NSCF,1,NCSF)
              IB_H0 = IB_H0 + NCSF
            ELSE
*. Keep complete block
              CALL COPVEC(dbl_mb(KLCSFHCSF),H0(IB_H0),NCSF**2)
              IB_H0 = IB_H0 + NCSF**2
            END IF
            CALL MEMCHK2('AFCOPX')
            IB_OCC = IB_OCC + IORB
            NCSF_TOT = NCSF_TOT + NCSF
          END DO! loop over confs for given IOPEN
        END DO! IOPEN
*
*. Write preconditioner block to disc
*
        LENGTH_H0 = IB_H0 - 1
        CALL TODSCN(H0,1,LENGTH_H0,-1,LUDIA)
        
        IF(NTEST.GE.100) THEN
          WRITE(6,*) 
          WRITE(6,*) ' Preconditioner for occupation class ', IOCCLS
          IF(IH0_CSF.EQ.2) THEN
            CALL WRTMAT(H0,1,NCSF_TOT,1,NCSF_TOT)
          ELSE
            IB_H0 = 1
            DO IOPEN = 0, MAXOP
             ITYP =  IOPEN + 1
             NCONF = NCONF_FOR_OPEN(IOPEN+1)
             NCSF = NPCSCNF(IOPEN+1)
             DO  ICONF = 1, NCONF
              CALL WRTMAT(H0(IB_H0),NCSF,NCSF,NCSF,NCSF)
              IB_H0 = IB_H0 + NCSF**2
             END DO! loop over confs for given IOPEN
            END DO! IOPEN
          END IF ! switch IH0_CSF
        END IF ! NTEST
*
       END DO! IIOCCLS
* And EOV 
       CALL ITODS(-1,1,-1,LUDIA)
*
       IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Complete preconditioner from Disc: '
        CALL WRTVCD(H0,LUDIA,1,-1)
       END IF
*
      END IF !ICNFBAT switch
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CSFH0 ')
      CALL LUCIAQEXIT('CSFH0 ')
      RETURN
      END
      SUBROUTINE ANACSF(CIVEC,ICONF_OCC,NCONF_FOR_OPEN,IPROCS,THRES,
     &           MAXTRM,IOUT)
*
* Analyze CI vector in CSF basis
*. Jeppe Olsen, Jan 2012 (after 25 years rest of code, not Jeppe)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      REAL*8 INPROD
*. Input
      DIMENSION ICONF_OCC(*),NCONF_FOR_OPEN(*)
      INTEGER IPROCS(1)
      DIMENSION CIVEC(*)
*. Local scratch
      CHARACTER*3 PAT(MXPORB)
      INTEGER IOCCL(MXPORB)
*
* The number of CSF's
*
      NCSF = 0
      DO IOPEN = 0, MAXOP
        NNCNF = NCONF_FOR_OPEN(IOPEN+1)
        NNCSF = NPCSCNF(IOPEN+1)
        NCSF = NCSF + NNCSF*NNCNF
      END DO
      WRITE(6,*) ' Number of CSF''s in expansion: ', NCSF
* Norm of totalwave function
      X2 = INPROD(CIVEC,CIVEC,NCSF)
      XNORM = SQRT(X2)
      WRITE(6,*) ' Norm of input vector = ', XNORM
       

* =================================
* Printout of configuration weights
* =================================
*
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(1H ,A)') ' ====================== '
      WRITE(IOUT,'(1H ,A)') ' Info on configurations '
      WRITE(IOUT,'(1H ,A)') ' ====================== '
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(1H ,A)') 
     &' (Negative orbital index implies doubly occupied orbital)'
*
      ITRM = 0
      ILOOP = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
*
3001  CONTINUE
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = XNORM
        XMIN = XMAX/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
C
      WRITE(IOUT,'(//A,E10.4,A,E10.4)')
     &'  Printout of contributions in interval  ',XMIN,' to ',XMAX
      WRITE(IOUT,'(A)')
     &'  =============================================================='
*. Loop over configurations and CSF's for given configuration
      ICNF = 0
      ICSF = 0
      INRANG = 0
      NOPRT = 0
      DO IOPEN = 0, MAXOP
        ITYP = IOPEN + 1
        ICL = (NACTEL - IOPEN) / 2
        IOCC = IOPEN + ICL
        IF( ITYP .EQ. 1 ) THEN
          ICNBS0 = 1
          IPBAS = 1
        ELSE
          ICNBS0 = ICNBS0 + NCONF_FOR_OPEN(IOPEN+1-1)*(NACTEL+IOPEN-1)/2
          IPBAS = IPBAS + NPCSCNF(IOPEN+1-1)*(IOPEN-1)
        END IF
*. Configurations of this type
        NNCNF = NCONF_FOR_OPEN(IOPEN+1)
        NNCSF = NPCSCNF(IOPEN+1)
        DO IC = 1, NNCNF
          ICNF = ICNF + 1
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
*. Weight of CSF's in this configuration
          W = 0.0D0
          DO IICSF = 1, NNCSF
            ICSF = ICSF+1
            W = W +  CIVEC(ICSF)*CIVEC(ICSF)
          END DO
          SQ_W = SQRT(W)
          IF(XMAX.GE.SQ_W .AND. SQ_W.GT.XMIN ) THEN
            ITRM  = ITRM + 1
            INRANG = INRANG + 1
            IF( ITRM .LE. MAXTRM ) THEN
              WRITE(IOUT,'(A,E10.5)') 
     &        '  Square root of weight for conf ',SQ_W
              CALL ICOPVE(ICONF_OCC(ICNBS),IOCCL,IOCC)
              WRITE(IOUT,*) ' Occupation of configuration: '
              IF(NINOB.NE.0) THEN
                CALL IADD_SIGN_CONST(IOCCL,NINOB,IOCC)
C                    IADD_SIGN_CONST(IVEC,IADD,NELMNT)
                WRITE(IOUT,'(4X,A,10(2X,I3),
     &                       (/,16X,10(2X,I3)))')
     &          ' (Inactive) ', (IOCCL(II),II = 1,IOCC)
              ELSE
                WRITE(IOUT,'(4X,10(2X,I3))')
     &          (IOCCL(II),II = 1,IOCC)
              END IF !NINOB switch
            ELSE
              NOPRT = NOPRT + 1
            END IF! MAXTRM check
          END IF! XMAX/XMIN check
        END DO ! Loop over configs
      END DO ! Loop over IOPEN
      IF(XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 3001
*
* ==========================
* Printout of CSF-weights
* ==========================
*
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(1H ,A)') ' ============ '
      WRITE(IOUT,'(1H ,A)') ' Info on CSFs '
      WRITE(IOUT,'(1H ,A)') ' ============ '
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      ITRM = 0
      ILOOP = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
2001  CONTINUE
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = XNORM
        XMIN = XMAX/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
C
      WRITE(IOUT,'(//A,E10.4,A,E10.4)')
     &'  Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(IOUT,'(A)')
     &'  =============================================================='
*. Loop over configurations and CSF's for given configuration
      ICNF = 0
      ICSF = 0
      INRANG = 0
      NOPRT = 0
      DO 1000 IOPEN = 0, MAXOP
        ITYP = IOPEN + 1
        ICL = (NACTEL - IOPEN) / 2
        IOCC = IOPEN + ICL
        IF( ITYP .EQ. 1 ) THEN
          ICNBS0 = 1
          IPBAS = 1
        ELSE
          ICNBS0 = ICNBS0 + NCONF_FOR_OPEN(IOPEN+1-1)*(NACTEL+IOPEN-1)/2
          IPBAS = IPBAS + NPCSCNF(IOPEN+1-1)*(IOPEN-1)
        END IF
*. Configurations of this type
        NNCNF = NCONF_FOR_OPEN(IOPEN+1)
        NNCSF = NPCSCNF(IOPEN+1)
        DO 900  IC = 1, NNCNF
          ICNF = ICNF + 1
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
*. CSF's in this configuration
          DO 800 IICSF = 1, NNCSF
            ICSF = ICSF+1
            IF( XMAX .GE. ABS(CIVEC(ICSF)) .AND.
     &         ABS(CIVEC(ICSF)).GT. XMIN ) THEN
              ITRM  = ITRM + 1
              INRANG = INRANG + 1
              IF( ITRM .LE. MAXTRM ) THEN
C                    SPIN_COUPLING_PATTERN(IOCC,ISPIN,IPAT,NCL,NOP)
                CALL SPIN_COUPLING_PATTERN(ICONF_OCC(ICNBS),
     &               IPROCS(IPBAS+(IICSF-1)*IOPEN),PAT,ICL,IOPEN)
*
                CNORM = CNORM + CIVEC(ICSF) ** 2
                WRITE(IOUT,*) ' Coefficient of CSF ',ICSF,CIVEC(ICSF)
                WRITE(IOUT,*) ' Occupation and spin coupling '
                IF(NINOB.NE.0) THEN
                  WRITE(IOUT,'(4X,A,10(2X,I3),/,16X,10(2X,I3))') 
     &            ' (Inactive) ', 
     &            (ABS(ICONF_OCC(ICNBS-1+II))+NINOB,II = 1,IOCC)
                  WRITE(IOUT,'(5X,A,10(2X,A3),/,17X,10(2X,A3))') 
     &            '            ', (PAT(II),II=1, IOCC) 
                ELSE
                  WRITE(IOUT,'(4X,A,10(2X,I3))') ' ',
     &            (ABS(ICONF_OCC(ICNBS-1+II)),II = 1,IOCC)
                  WRITE(IOUT,'(5X,A,10(2X,A3))') ' ',
     &            (PAT(II),II=1, IOCC) 
                END IF! NINOB switch
                NOPRT = NOPRT + 1
              END IF
            END IF
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
      NCIVAR = ICSF
      IF(INRANG .EQ. 0 ) WRITE(IOUT,*) '   ( no coefficients )'
      IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
      IF(NOPRT.NE.0) WRITE(IOUT,*)
     &' Number of coefficients not printed ', NOPRT
       WRITE(IOUT,'(//A,E15.8)')
     & '  Norm of printed CI vector .. ', CNORM
*
       WRITE(IOUT,'(/A)') '   Magnitude of CI coefficients '
       WRITE(IOUT,'(A/)') '  =============================='
*
       CNORM = 0.0D0
       ISUM = 0
       XMIN = 1.0D0
       DO 200 IPOT = 0, 10
         CLNORM = 0.0D0
         INRANG = 0
         XMAX = XMIN
         XMIN = XMIN * 0.1D0
C
         DO 180 IDET = 1, ICSF
           IF( ABS(CIVEC(IDET)) .LE. XMAX  .AND.
     &         ABS(CIVEC(IDET)) .GT. XMIN ) THEN
                 INRANG = INRANG + 1
                 CLNORM = CLNORM + CIVEC(IDET) ** 2
           END IF
  180    CONTINUE
         CNORM = CNORM + CLNORM
C
         IF (INRANG .GT. 0)
     &     WRITE(IOUT,'(A,I2,A,I2,3X,I7,3X,E15.8,3X,E15.8)')
     &     '  10-',IPOT+1,' to 10-',IPOT,INRANG,CLNORM,CNORM
C
         ISUM = ISUM + INRANG
  200 CONTINUE
C
      WRITE(IOUT,*) ' Number of coefficients less than  10-11',
     &           ' is  ',NCIVAR - ISUM
C
      RETURN
      END
      SUBROUTINE SPIN_COUPLING_PATTERN(IOCC,ISPIN,IPAT,NCL,NOP)
*
* A configuration IOCC and a coupling of the open orbitals ISPIN is given.
* Find spin-coupling pattern of this and save in IPAT
*
*. Jeppe Olsen, Jan 2012
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC(NCL+NOP),ISPIN(NOP)
*. Output
      CHARACTER*3 IPAT(*)
*
      NTEST = 00
*
      IOP = 0
      DO IORB = 1, NCL + NOP
        IF(IOCC(IORB).LT.0) THEN
*. Double occupied
          IPAT(IORB) = ' d '
        ELSE 
*. Single occupied
         IOP = IOP + 1
         IF(ISPIN(IOP).EQ.1) THEN
           IPAT(IORB) = ' + '
         ELSE
           IPAT(IORB) = ' - '
         END IF
        END IF 
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SPIN_COUPLING_PATTERN '
        WRITE(6,*)
        WRITE(6,*) ' Input configuration: '
        CALL IWRTMA(IOCC,1,NCL+NOP,1,NCL+NOP)
        WRITE(6,*) ' Input spincoupling pattern of open orbitals '
        CALL IWRTMA(ISPIN,1,NOP,1,NOP)
        WRITE(6,*) ' Output spin coupling pattern '
        WRITE(6,'(1H , 20A3)') (IPAT(IORB),IORB = 1, NOP)
      END IF
*
      RETURN
      END
      SUBROUTINE IADD_SIGN_CONST(IVEC,IADD,NELMNT)
*
* IVEC(I) = IVEC(I) + SIGN(IVEC(I))*IADD
*
*PS: sign of zero is defined positive...
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER IVEC(NELMNT)
*
      DO I = 1, NELMNT
        IF(IVEC(I).GE.0) THEN
          IVEC(I) = IVEC(I) + IADD
        ELSE
          IVEC(I) = IVEC(I) - IADD
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE HCONF_DIA(HCSF_DIA,HCNF_AVE_DIA,HCNF_MIN_DIA,
     &           NCONF_FOR_OPEN,NPCSCNF,MAXOP)
*
* Obtain diagonals of Hamilton over configurations:
*   HCNF_MIN_DIA(ICONF): Lowest element in ICONF
*   HCNF_AVE_DIA(ICONF): Average of CSF's in ICONF
*
*. Jeppe Olsen, Jan. 2012
*
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION HCSF_DIA(*)
      INTEGER NCONF_FOR_OPEN(*), NPCSCNF(*)
*. Output
      DIMENSION HCNF_AVE_DIA(*), HCNF_MIN_DIA(*)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' HCONF_DIA speaking '
        WRITE(6,*) ' ================== '
        WRITE(6,*)
      END IF
*
      IB_CNF = 0
      IB_CSF = 0
      DO IOPEN = 0, MAXOP
        NNCSF = NPCSCNF(IOPEN+1)
        DO IICNF = 1, NCONF_FOR_OPEN(IOPEN+1)
*. Average energy
         ESUM = ELSUM(HCSF_DIA(IB_CSF),NNCSF)
         IF(NNCSF.NE.0) THEN
          EAVE = ESUM/FLOAT(NNCSF)
         ELSE
          EAVE = 123456789.0D0
         END IF
*. Lowest energy
         EMIN = FNDMNX(HCSF_DIA(IB_CSF),NNCSF,1)
*. Save info
         HCNF_MIN_DIA(IB_CNF) = EMIN
         HCNF_AVE_DIA(IB_CNF) = EAVE
*. Update pointers
         IB_CNF = IB_CNF + 1
         IB_CSF = IB_CSF + NNCSF
        END DO! over CSF's
      END DO !over IOPEN
      NCONF_TOT = IB_CNF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) 
     & 'Average and minimal energy diagonal over configurations'
       WRITE(6,*)
       CALL WRTMAT(HCNF_AVE_DIA,1,NCONF_TOT,1,NCONF_TOT)
       CALL WRTMAT(HCNF_MIN_DIA,1,NCONF_TOT,1,NCONF_TOT)
      END IF
*
      RETURN
      END
      SUBROUTINE ISCLVEC(IVEC,ISCAL,NDIM)
*
* IVEC(I) = ISCAL*IVEC(I)
*
*
      INCLUDE 'implicit.inc'
      INTEGER IVEC(*)
*
      DO I = 1, NDIM
        IVEC(I) = ISCAL*IVEC(I)
      END DO
*
      RETURN
      END
      SUBROUTINE WRT_CONF_LIST2
     &           (ICONF,IB_CONF_OCC,NCONF_FOR_OPEN,MAXOP,NCONF,NELEC)
*
* Write list of configurations, given in packed form
*
*
* Jeppe Olsen, January 2011 
*              IB_CONF_OCC added (compared to WRT_CONF_LIST)
*
      INCLUDE 'implicit.inc'
*   
      INTEGER ICONF(*), NCONF_FOR_OPEN(MAXOP+1)
      INTEGER IB_CONF_OCC(MAXOP+1)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' WRT_CONF_LIST2 in action '
        WRITE(6,*) ' ========================='
        WRITE(6,*)
        WRITE(6,*) ' NCONF_FOR_OPEN and IB_CONF_OCC:'
        CALL IWRTMA(NCONF_FOR_OPEN,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*)
        CALL IWRTMA(IB_CONF_OCC,1,MAXOP+1,1,MAXOP+1)
      END IF
      DO IOPEN = 0, MAXOP
        NCONF_OP = NCONF_FOR_OPEN(IOPEN+1)
        IF(NTEST.GE.100)
     &  WRITE(6,*) ' Test: IOPEN, NCONF_OP = ', IOPEN, NCONF_OP  
        IF(NCONF_OP.NE.0) THEN
          WRITE(6,'(A,I3,A,I6)') 
     &    ' Number of configurations with ', IOPEN, 
     &               ' open orbitals is ', NCONF_OP
*
          NOCC_ORB = IOPEN + (NELEC-IOPEN)/2
          IB = IB_CONF_OCC(IOPEN+1)
          IF(NTEST.GE.100) WRITE(6,*) '  IB = ', IB
          DO JCONF = 1, NCONF_OP
            CALL IWRTMA(ICONF(IB),1,NOCC_ORB,1,NOCC_ORB)
            IB = IB + NOCC_ORB
          END DO
        END IF
      END DO
*
      RETURN
      END
      SUBROUTINE GEN_OCCLS_FOR_CISPAC
     &(IGSOCC_MIN, IGSOCC_MAX,ISPC,IOCCLS,NOCCLS)
*
* Obtain and store in IOCCLS the NOCCLS occupation classes for the MINMAX space defined by 
* IGSOCC_MIN, IGS_OCC_MAX.
*
*. Jeppe Olsen, Jan. 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
*. Input
      INTEGER IGSOCC_MIN(MXPNGAS), IGSOCC_MAX(MXPNGAS)
*. Local scratch
      DIMENSION IOCC_A(MXPNGAS), IOCC(MXPNGAS)
*. Output
      INTEGER IOCCLS(NGAS,NOCCLS)
*
      CALL LUCIAQENTER('GNOCCL')
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) '  GEN_OCCLS_FOR_CISPAC in action '
         WRITE(6,*) ' ================================'
C?       WRITE(6,*) ' NGAS ', NGAS
         WRITE(6,*) ' IGSOCC_MIN, IGSOCC_MAX: '
         CALL WRT_MINMAX_GASOCC(IGSOCC_MIN,IGSOCC_MAX,NGAS)
      END IF
*   
      MAXLOOP = 1000000
cVOG  MAXLOOP = 10000
      NLOOP = 0
      NOCCLSL = 0
*. Loop over occupation classes
      INI  = 1
 1000 CONTINUE
        CALL NXT_OCCLS_IN_MINMAX_SPC(IOCC_A,IGSOCC_MIN,IGSOCC_MAX,NGAS,
     &       NOBPT,INI,NONEW)
        INI = 0
        IF(NONEW.EQ.0) THEN
          NLOOP =  NLOOP + 1
*. Change occupation class from accumulated to actual occupation and save
          CALL REF_OCCLS(IOCC_A,IOCC,1,NGAS)
*. If requested, check for occupations in the compound GAS space
          IM_IN = 1
          IF(I_CHECK_ENSGS.EQ.1)
     &    CALL  CHECK_IS_OCC_IN_ENGSOCC(IOCC,-1,IM_IN)
C               CHECK_IS_OCC_IN_ENGSOCC(IGSOCCL,ISPC,IM_IN)
          IF(IM_IN.EQ.1) THEN
           NOCCLSL = NOCCLSL + 1
           CALL ICOPVE(IOCC,IOCCLS(1,NOCCLSL),NGAS)
          END IF
        END IF! NONEW .eq. 0
      IF(NONEW.EQ.0.AND.NLOOP.LT.MAXLOOP) GOTO 1000
      IF(NLOOP.EQ.MAXLOOP) THEN
        WRITE(6,*) ' Forced MAXLOOP exit in MAXOP...'
        WRITE(6,*) ' Bug or test to be removed, Jeppe!!! '
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of occupation generated', NOCCLSL
      END IF
      IF(NTEST.GE.1000) THEN
        CALL WRITE_OCCLS(IOCCLS,NOCCLS,NGAS)
      END IF
*
      CALL LUCIAQEXIT('GNOCCL')
      RETURN
      END
      SUBROUTINE WRITE_OCCLS(IOCCLS,NOCCLS,NGAS)
*
* Print the NOCCLS classes IOCCLS
*
*. Jeppe Olsen, Jan. 2012
*
      INCLUDE 'implicit.inc'
      INTEGER IOCCLS(NGAS,NOCCLS)
*
      WRITE(6,*)
      WRITE(6,*) ' ======================='
      WRITE(6,*) ' The occupation classes '  
      WRITE(6,*) ' ======================='
      WRITE(6,*)
*
      WRITE(6,'(A, 16(2X,I4))')   '  Gas: ', (I, I = 1, NGAS)
      WRITE(6,*)                  ' -----'
      DO KOCCLS = 1, NOCCLS
        WRITE(6,'(A, 16(2X,I4))') '       ', (IOCCLS(I,KOCCLS),I=1,NGAS)
      END DO
*
      RETURN
      END
      SUBROUTINE GEN_INFO_FOR_ALL_OCCLS(I_DO_SBCNF)
*
* Set up the information for the NOCCLS occupation classes defined in 
* WORK(KIOCCLS). Only information  about the number of configurations
* and not the actual configurations are generated
*
* Information about dimensions of CI-spaces are also determined
*
*. Jeppe Olsen, Jan. 2012
*  April 2, 2013; Jeppe Olsen; Subconf approach added
*  Latest modification; May 17 2013; Jeppe Olsen; Info on AB SDs and AB CMs added
*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      CALL LUCIAQENTER('GNINCL')
*
      NTEST = 0
C?    WRITE(6,*) ' GEN_INFO_FOR_ALL.., I_DO_SBCNF = ', I_DO_SBCNF
*
*. Call the slaves to do the work (avoid real/integer problems)
*
      CALL GEN_INFO_FOR_ALL_OCCLS_S(int_mb(KIOCCLS),NOCCLS_MAX,
     &     MINOP,MAXOP,NIRREP,NINOB,NGAS,NOBPT,
     &     int_mb(KNCN_PER_OP_SM),int_mb(KNCN_ALLSYM_FOR_OCCLS),
     &     int_mb(KNCN_FOR_OCCLS),
     &     NCN_ALLSYM_FOR_OCCLS_MAX, NCN_FOR_OCCLS_MAX,
     &     LEN_OCC_FOR_OCCLS_MAX,I_DO_SBCNF)
C          GEN_INFO_FOR_ALL_OCCLS_S(IOCCLS,NOCCLS_MAX,
C    &     MINOP,MAXTOP,NIRREP,NINOB,NGAS,NOBPT,
C    &     NCN_PER_OP_SM,NCN_ALLSYM_FOR_OCCLS,
C    &     NCN_FOR_OCCLS,
C    &     NCN_ALLSYM_FOR_OCCLS_MAX, NCN_FOR_OCCLS_MAX,
C    &     LEN_OCC_FOR_OCCLS_MAX,I_DO_SBCNF)
     
*
*. Number of CSFs, CMs, SDs per occupation class
*
      CALL GEN_NCSF_ETC_FOR_ALL_OCCLS(int_mb(KNCN_PER_OP_SM),
     &     MAXOP,NIRREP,NOCCLS_MAX,NPCSCNF,NPCMCNF,NPDTCNF,
     &     int_mb(KNCS_FOR_OCCLS),
     &     int_mb(KNCM_FOR_OCCLS), int_mb(KNSD_FOR_OCCLS),
     &     NCS_FOR_OCCLS_MAX,NCM_FOR_OCCLS_MAX,
     &     NSD_FOR_OCCLS_MAX                     )
C          GEN_NCSF_ETC_FOR_ALL_OCCLS(NCNF_PER_OP_SM,
C    &     MAXOP,NIRREP,NOCCLS_MAX,
C    &     NPCSCNF,NPCMCNF,NPDTCNF,
C    &     NCN_FOR_OCCLS, NCS_FOR_OCCLS,
C    &     NCM_FOR_OCCLS, NSD_FOR_OCCLS)
*
*. Number of AB SDs and CMs per occupation class
*
C     GEN_NSDAB_FOR_ALL_OCCLS(
C    &           NABSPGP_PER_OCCLS,IABSPGP_PER_OCCLS,
C    &           N_SDAB_PER_OCCLS,N_CMAB_PER_OCCLS,
C                N_SDAB_PER_OCCLS_MAX,N_CMAB_PER_OCCLS_MAX)
      CALL GEN_NSDAB_FOR_ALL_OCCLS(
     &     int_mb(KNABSPGP_FOR_OCCLS),int_mb(KIABSPGP_FOR_OCCLS),
     &     int_mb(KNSDAB_FOR_OCCLS),int_mb(KNCMAB_FOR_OCCLS),
     &     N_SDAB_PER_OCCLS_MAX, N_CMAB_PER_OCCLS_MAX) 
*
*. Number of CSF's, CM's, SD's and Confs' per CISPACE
*
C          OCCLS_TO_SPACE_DIM(NCN_FOR_OCCLS,NCS_FOR_OCCLS,
C    &           NCM_FOR_OCCLS,NSD_FOR_OCCLS,NCN_ALLSYM_FOR_OCCLS,
C    &           IOCCLS_ACT)
      CALL OCCLS_TO_SPACE_DIM(int_mb(KNCN_FOR_OCCLS),
     &     int_mb(KNCS_FOR_OCCLS),int_mb(KNCM_FOR_OCCLS),
     &     int_mb(KNSD_FOR_OCCLS),
     &     int_mb(KNCN_ALLSYM_FOR_OCCLS),int_mb(KSIOCCLS_ACT))
*
      WRITE(6,*) ' First 4 elements of WORK(KNCMAB_FOR_OCCLS) '
      CALL IWRTMA(int_mb(KNCMAB_FOR_OCCLS),1,3,1,3)
*
*. Length of configuration expansions
*
      NELL = IELSUM(int_mb(KIOCCLS),NGAS)
C     FUNCTION LEN_OCCLIST(NCONF_PER_OPEN, MAXOP, NELEC)
      DO ISPC = 1, NCISPC
       DO ISM = 1, NSMOB
        LCONFOCC_PER_SYM_GN(ISM,ISPC) = 
     &  LEN_OCCLIST(NCONF_PER_OPEN_GN(1,ISM,ISPC),MAXOP,NELL)
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Length of confocc for each sym and space '
        CALL IWRTMA(LCONFOCC_PER_SYM_GN,NSMOB,NCISPC,MXPCSM,
     &               MXPICI)
      END IF
*
*. Obtain various max dims
*
      CALL GET_MAX_CONF_DIMS
*
      CALL LUCIAQEXIT('GNINCL')
      RETURN
      END
      SUBROUTINE GEN_INFO_FOR_ALL_OCCLS_S(IOCCLS,NOCCLS_MAX,
     &     MINOP,MAXOP,NIRREP,NINOB,NGAS,NOBPT,
     &     NCN_PER_OP_SM,NCN_ALLSYM_FOR_OCCLS,
     &     NCN_FOR_OCCLS,
     &     NCN_ALLSYM_FOR_OCCLS_MAX, NCN_FOR_OCCLS_MAX,
     &     LEN_OCC_FOR_OCCLS_MAX,I_DO_SBCNF)
*
* Generate information on number of configurations: NCN_FOR_OP_SP
* various occupation classes and LEN_OCC_OR_OCCLS_MAX,
* NCN_ALLSYM_FOR_OCCLS_MAX, NCN_FOR_OCCLS_MAX
*
*. Jeppe Olsen, Jan. 2012
*  Latest modification; April 2, 2013; Jeppe Olsen; subconf approach added
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCCLS(NGAS,NOCCLS_MAX)
      INTEGER NOBPT(NGAS)
*. Output
      INTEGER NCN_PER_OP_SM(MAXOP+1,NIRREP,NOCCLS_MAX)
      INTEGER NCN_ALLSYM_FOR_OCCLS(NOCCLS_MAX)
      INTEGER NCN_FOR_OCCLS(NIRREP,NOCCLS_MAX) 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GEN_INFO_FOR_ALL_OCCLS_S '
        WRITE(6,*) ' ================================== '
        WRITE(6,*) ' NOCCLS_MAX = ', NOCCLS_MAX
        WRITE(6,*) ' I_DO_SBCNF = ', I_DO_SBCNF
      END IF
*
      LEN_OCC_FOR_OCCLS_MAX = 0
      NELEC = IELSUM(IOCCLS(1,1),NGAS)
*
      DO JOCCLS = 1, NOCCLS_MAX
        IF(I_DO_SBCNF.EQ.0) THEN
          CALL GEN_NCONF_FOR_OCCLSN(IOCCLS(1,JOCCLS),NGAS,
     &         MINOP,MAXOP,NIRREP,NINOB+1,NOBPT,
     &         NCN_ALLSYM_FOR_OCCLS(JOCCLS),
     &         NCN_PER_OP_SM(1,1,JOCCLS),NCN_FOR_OCCLS(1,JOCCLS))
        ELSE
C              OCCLSDIM_FROM_OCSBCLSDIM(IOCCLS,IOCCLSDIM)
          CALL OCCLSDIM_FROM_OCSBCLSDIM(IOCCLS(1,JOCCLS),
     &         NCN_PER_OP_SM(1,1,JOCCLS),MAXOP)
*
*. Number of confs per sym
*
          DO ISM = 1, NIRREP
            NCN_FOR_OCCLS(ISM,JOCCLS) = 
     &      IELSUM(NCN_PER_OP_SM(1,ISM,JOCCLS),MAXOP+1)
          END DO
*
*. Number of confs per occupation class without restriction on number of open orbitals
*
C              NCONF_OCCLS(NOBPSP,NELPSP,NSPC,MINOP,NNCONF,LLOCC)
          CALL NCONF_OCCLS(NOBPT,IOCCLS(1,JOCCLS),NGAS,0,NNCONF,LLOCC)
          NCN_ALLSYM_FOR_OCCLS(JOCCLS) = NNCONF
*
*. Largest block of occupation for any symmetry
*
          DO ISM = 1, NIRREP
           LENNY = LEN_OCCLIST(NCN_PER_OP_SM(1,ISM,JOCCLS),
     &             MAXOP, NELEC)
C                  LEN_OCCLIST(NCONF_PER_OPEN, MAXOP, NELEC)
           LEN_OCC_FOR_OCCLS_MAX =
     &     MAX(LEN_OCC_FOR_OCCLS_MAX, LENNY)
          END DO
       END IF ! Subconfigurations should be generated
      END DO! loop over JOCCLS
*
      NCN_ALLSYM_FOR_OCCLS_MAX =IMNMX(NCN_ALLSYM_FOR_OCCLS,NOCCLS_MAX,2)
      NCN_FOR_OCCLS_MAX = IMNMX(NCN_FOR_OCCLS,NIRREP*NOCCLS_MAX,2)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Total number of configs per occupation class '
        CALL IWRTMA10(NCN_ALLSYM_FOR_OCCLS,1,NOCCLS_MAX,1,NOCCLS_MAX)
*
        WRITE(6,*) ' LEN_OCC_FOR_OCCLS_MAX = ', LEN_OCC_FOR_OCCLS_MAX
        WRITE(6,*) ' NCN_ALLSYM_FOR_OCCLS_MAX, NCN_FOR_OCCLS_MAX  = ',
     &               NCN_ALLSYM_FOR_OCCLS_MAX, NCN_FOR_OCCLS_MAX
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_NCONF_FOR_OCCLS(IOCCLS,NGAS,MINOP,MAXOP,NIRREP,
     &           IB_ORB,NOBPT,NCONF_OCCLS_ALLSYM,NCONF_PER_OP_SM,
     &           NCONF_OCCLS)
*
*
* Generate 
* NCONF_PER_OP_SM: Number of configurations of occclass IOCCLS 
*                  for all symmetries and atleast MINOP open orbitals 
* NCONF_OCCLS_ALLSYM: Number of configurations of occclass IOCCLS 
*                     and any symmetry and any number of open orbitals
* NCONF_OCCLS: Number of configurations of occlass IOCCLS af given symmetry
*              and atleast MINOP open orbitals
*
* Jeppe Olsen, Jan. 2012
*       
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
*
*.. Input
*
*. Number of electrons per gas space 
      INTEGER IOCCLS(NGAS)  
*. Number of orbitals per gasspace 
      INTEGER NOBPT(NGAS)
*
*.. Output
*
*. Number of configurations per number of open shells, all symmetries
      INTEGER NCONF_PER_OP_SM(MAXOP+1,NIRREP)
      INTEGER NCONF_OCCLS(NIRREP)
*. Local scratch
      INTEGER JCONF(2*MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GEN_NCONF_FOR_OCCLS'
        WRITE(6,*) ' ============================ '
        WRITE(6,*) 
        WRITE(6,*) ' Occupation class in action '    
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NGAS, MAXOP = ', NGAS, MAXOP
        WRITE(6,*) ' IB_ORB = ', IB_ORB
        WRITE(6,*) ' NOBPT: '
        CALL IWRTMA(NOBPT,1,NGAS,1,NGAS)
      END IF
*. Total number of electrons 
      NEL = IELSUM(IOCCLS,NGAS)
      IF(NTEST.GE.1000) WRITE(6,*) ' NEL = ', NEL
*
      IZERO = 0
      CALL ISETVC(NCONF_PER_OP_SM,IZERO,(MAXOP+1)*NIRREP)
      CALL ISETVC(NCONF_OCCLS,IZERO,NIRREP)
*. Loop over configurations 
      INI = 1
      NCONF = 0
      ISUM = 0
      NCONF_OCCLS_ALLSYM = 0
*. Loop over configurations
 1000 CONTINUE
        IF(NTEST.GE.10000)
     &  WRITE(6,*) ' NEXT_CONF_FOR_OCCLS will be called '
        CALL NEXT_CONF_FOR_OCCLS
     &        (JCONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
        ISUM = ISUM + 1
        INI = 0
*
        IF(NONEW.EQ.0) THEN
*. Check symmetry and number of open orbitals for this space

          IADD = IB_ORB - 1
C              IADCONST(IVEC,IADD, NDIM)
          CALL IADCONST(JCONF,IADD,NEL)
          ISYM_CONF = ISYMST(JCONF,NEL)
          IADD = - IADD
          CALL IADCONST(JCONF,IADD,NEL)
          NOPEN     = NOP_FOR_CONF(JCONF,NEL) 
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Number of open shells and SYM', 
     &      NOPEN, ISYM_CONF
          END IF
          NOCOB =  NOPEN + (NEL-NOPEN)/2
          NCONF_OCCLS_ALLSYM = NCONF_OCCLS_ALLSYM + 1 
          IF(NOPEN.GE.MINOP) THEN
*. A new configuration to be included, reform and save in packed form
            NCONF_PER_OP_SM(NOPEN+1,ISYM_CONF) =
     &      NCONF_PER_OP_SM(NOPEN+1,ISYM_CONF) + 1
            NCONF_OCCLS(ISYM_CONF) = NCONF_OCCLS(ISYM_CONF) + 1
          END IF
      GOTO 1000
        END IF !End if nonew = 0
* 
C     WRITE(6,*) ' TEST,  NCONF_OCCLS_ALLSYM = ',  NCONF_OCCLS_ALLSYM
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*)  ' =============================================== '
        WRITE(6,*)  ' Information on number of configuration for occls'
        WRITE(6,*)  ' =============================================== '
        WRITE(6,*)
        WRITE(6,*) ' Occupation class in action: '
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of included configurations per symmetry'
        CALL IWRTMA(NCONF_OCCLS,1,NIRREP,1,NIRREP)
        WRITE(6,*) 
     &  ' Numbers:  Open orbitals +1 (ROW) and  sym (COLUMN)  '
          CALL IWRTMA(NCONF_PER_OP_SM,MAXOP+1,NIRREP,MAXOP+1,NIRREP)
      END IF

*
      RETURN
      END
      SUBROUTINE GEN_NCSF_ETC_FOR_ALL_OCCLS(NCNF_PER_OP_SM,
     &           MAXOP,NIRREP,NOCCLS_MAX,
     &           NPCSCNF,NPCMCNF,NPDTCNF,
     &           NCS_FOR_OCCLS,
     &           NCM_FOR_OCCLS, NSD_FOR_OCCLS,
     &           NCS_FOR_OCCLS_MAX,NCM_FOR_OCCLS_MAX,
     &           NSD_FOR_OCCLS_MAX                     )
*
* Generate number of CSFs, SDs CMs for the occupation classes
*
*. Jeppe Olsen, Jan. 2012
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NCNF_PER_OP_SM(MAXOP+1,NIRREP,NOCCLS_MAX)
      INTEGER NPCSCNF(MAXOP+1), NPCMCNF(MAXOP+1),
     &        NPDTCNF(MAXOP+1)
*. Output
      INTEGER NCS_FOR_OCCLS(NIRREP,NOCCLS_MAX),
     &        NCM_FOR_OCCLS(NIRREP,NOCCLS_MAX),
     &        NSD_FOR_OCCLS(NIRREP,NOCCLS_MAX) 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  Output from GEN_NCSF_ETC_FOR_ALL_OCCLS '
        WRITE(6,*) '  ======================================='
        WRITE(6,*) 
        WRITE(6,*) ' NOCCLS_MAX, NIRREP = ', NOCCLS_MAX, NIRREP
      END IF
*
      DO JOCCLS = 1, NOCCLS_MAX
       DO ISM = 1, NIRREP
*
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' JOCCLS, ISM = ', JOCCLS, ISM
           WRITE(6,*) ' NCNF_PER_OP_SM: '
           CALL IWRTMA(NCNF_PER_OP_SM(1,ISM,JOCCLS),
     &                 1,MAXOP+1,1,MAXOP+1)
         END IF
*
        CALL NCNF_TO_NCOMP(
     &       MAXOP,NCNF_PER_OP_SM(1,ISM,JOCCLS),NPCSCNF,
     &       NCS_FOR_OCCLS(ISM,JOCCLS))
        CALL NCNF_TO_NCOMP(
     &       MAXOP,NCNF_PER_OP_SM(1,ISM,JOCCLS),NPCMCNF,
     &       NCM_FOR_OCCLS(ISM,JOCCLS))
        CALL NCNF_TO_NCOMP(
     &       MAXOP,NCNF_PER_OP_SM(1,ISM,JOCCLS),NPDTCNF,
     &       NSD_FOR_OCCLS(ISM,JOCCLS))
       END DO
      END DO
*
*. Max dimensions
*
      NCS_FOR_OCCLS_MAX = IMNMX(NCS_FOR_OCCLS,NIRREP*NOCCLS_MAX,2)
      NCM_FOR_OCCLS_MAX = IMNMX(NCM_FOR_OCCLS,NIRREP*NOCCLS_MAX,2)
      NSD_FOR_OCCLS_MAX = IMNMX(NSD_FOR_OCCLS,NIRREP*NOCCLS_MAX,2)
*
      IF(NTEST.GE.100) THEN 
       DO ISM = 1, NIRREP
        WRITE(6,'(A,I4)') 
     &  ' Number of CSFs, CMs, SDs per occupation class for sym', ISM
        WRITE(6,*) 
     &  ' ==========================================================='
        WRITE(6,*)
     &    '   Occ class      CSFs       CMs        SDs '
        WRITE(6,*) 
     &    '   =========================================='
        DO JOCCLS = 1, NOCCLS_MAX
          WRITE(6,'(1H ,I6, 5X, 3(I9,2X))')
     &    JOCCLS, NCS_FOR_OCCLS(ISM,JOCCLS),NCM_FOR_OCCLS(ISM,JOCCLS),
     &    NSD_FOR_OCCLS(ISM,JOCCLS)
        END DO
       END DO
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,'(A,I9)') 
     &  ' Largest number of CSFs in a occupation class ',
     &  NCS_FOR_OCCLS_MAX
        WRITE(6,'(A,I9)') 
     &  ' Largest number of CMs in a occupation class ',
     &  NCM_FOR_OCCLS_MAX
        WRITE(6,'(A,I9)') 
     &  ' Largest number of SDs in a occupation class ',
     &  NSD_FOR_OCCLS_MAX
      END IF
*
      RETURN
      END
      FUNCTION LEN_OCCLIST(NCONF_PER_OPEN, MAXOP, NELEC)
*
* Obtain length of list of occupations for configurationlist 
* specified by NCONF_PER_OPEN, MAXOP, NELEC
*
*. Jeppe Olsen, February 2012
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NCONF_PER_OPEN(MAXOP+1)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Entering LEN_OCCLIST '
      END IF
*
      LEN_CONFOCC = 0
      DO NOP = 0, MAXOP
        NCL = (NELEC - NOP)/2
        NOB =  NOP + NCL
        LEN_CONFOCC = LEN_CONFOCC + NOB*NCONF_PER_OPEN(NOP+1)
      END DO
*
      LEN_OCCLIST = LEN_CONFOCC
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Length of list of occupations ', LEN_CONFOCC
      END IF
*
      RETURN
      END
      SUBROUTINE ABSPGP_TO_OCCLS(
     &           IASPGP,NASPGP,IBSPGP,NBSPGP,
     &           NGAS,IGSOCC_MNMX,NELFGP,
     &           NOCCLS,IOCCLS,
     &           N_ABSPGP_TOT,I_ABSPGP_FOR_OCCLS,N_ABSPGP_FOR_OCCLS,
     &           IB_ABSPGP_FOR_OCCLS,IONLY_NABSPGP,NNODES,N_ABSPGP_MAX,
     &           I_ABSPGP_A,I_ABSPGP_B,N_ABSPGP)
*
* Obtain the combination of A- and B- supergroups that
* belong to a given occupation class.
*
* IF IONLY_NABPSGP = 1, Then only the arrays NABSPGP_FOR_OCCLS,
* IB_ABSPGP_FOR_OCCLS and NABSP_TOT are set up
*
*
* IGSOCC_MNMX is compound space and is used for pre-screening
* 
*
*. Jeppe Olsen, Feb. 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IASPGP(MXPNGAS,NASPGP),IBSPGP(MXPNGAS,NBSPGP)
      INTEGER IGSOCC_MNMX(MXPNGAS,2)
      INTEGER NELFGP(*)
      INTEGER IOCCLS(NGAS,*)
      INTEGER I_ABSPGP_A(N_ABSPGP_MAX), I_ABSPGP_B(N_ABSPGP_MAX),
     &        N_ABSPGP(NNODES)
*. Output
      INTEGER N_ABSPGP_FOR_OCCLS(NOCCLS)
      INTEGER IB_ABSPGP_FOR_OCCLS(NOCCLS)
      INTEGER I_ABSPGP_FOR_OCCLS(2,*)
*. Local scratch
      INTEGER JAJBOC(MXPNGAS), N_ABSPGP_VOG(NOCCLS) !VOG
*. 
      CALL LUCIAQENTER('ABSPG')
      NTEST = 0
      IF(ga_nodeid().eq.0.AND.NTEST.GE.10) THEN
        WRITE(6,*) ' Information from ABSPGP_TO_OCCLS '
        WRITE(6,*) '==================================='
        WRITE(6,*) 
        IF(IONLY_NABSPGP.NE.0) THEN
          WRITE(6,*) ' Only construction of dimension arrays'
        ELSE
          WRITE(6,*) ' AB supergroups for each occlass constructed'
        END IF
        write(*,*)
        write(*,*) 'Max number of AB supergroups:', N_ABSPGP_MAX 
      END IF
      IF(ga_nodeid().eq.0.AND.NTEST.GE.1000) THEN
        WRITE(6,*) ' Input occupation classes '
        CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
      END IF
*
      N_ABSPGP_TOT = 0
      IZERO = 0
      JJOCCLS = 0
      IF(ga_nodeid().eq.0.AND.NTEST.GE.10)
     &       WRITE(6,*) ' NOCCLS in AB.... ', NOCCLS
      CALL ISETVC(N_ABSPGP_FOR_OCCLS,IZERO,NOCCLS)
*
*. VOG Parallelize outer loop (JASPGP)
*
      !
      ! Initialize
      !
      call dfill(NNODES, 0.0d0, N_ABSPGP, 1)
      !
      ! Organize distribution of alpha supergroups
      ! OCCLMAP is a GA that holds the occupation map
      !
      call ga_sync()
      NUMBLOCK = INT(NASPGP / ga_nnodes())
      NREM     = MOD(NASPGP,ga_nnodes())
      ihi = 0
      DO I = 1, ga_nodeid() + 1
       IF (I.LE.NREM) THEN
        NSIZE = NUMBLOCK + 1
       ELSE
        NSIZE = NUMBLOCK
       END IF
       ihi = ihi + NSIZE
       ilo = ihi - NSIZE + 1
      END DO
      if ( .not. ga_create(mt_int, NOCCLS, NNODES, 'OCCLMAP', NOCCLS,
     &     1, IOCCLMAP))
     &     call errquit('OCCLMAP ga_create fail', NOCCLS*(NNODES+1),
     &                  GA_ERR)
      call ga_zero(IOCCLMAP)
      call ga_sync()
      !
      ! Print some output (for debugging)
      !
      IF(NTEST.GE.10.and.ga_nodeid().eq.0)
     &                     write(*,*) 'NASPGP,NUMBLOCK,NREM',
     &                                 NASPGP,NUMBLOCK,NREM
      IF(NTEST.GE.10)
     &  write(*,*) 'Distribution', ga_nodeid(), ilo, ihi
*
*. STEP 1
*. Construct occupation map
*
      DO JASPGP = ilo,ihi
cVOG  DO JASPGP = 1, NASPGP
       DO JBSPGP = 1, NBSPGP
*. Is this combination included in combination
        I_AM_IN = IS_ABSPGP_IN_CI_SPACE(
     &            IASPGP(1,JASPGP),IBSPGP(1,JBSPGP),NGAS,
     &            1,1,IGSOCC_MNMX,NELFGP)
        IF(I_AM_IN.EQ.1) THEN
*. Find the occupation class housing the combination. 
*. First version is based on searching all occupation classes...
          DO IGAS = 1, NGAS
            JAJBOC(IGAS) = 
     &      NELFGP(IASPGP(IGAS,JASPGP)) + NELFGP(IBSPGP(IGAS,JBSPGP))
          END DO
          IM_IT = 3006
          DO JOCCLS = 1, NOCCLS
           IM_IT = 1
           DO IGAS = 1, NGAS
             IF(JAJBOC(IGAS).NE.IOCCLS(IGAS,JOCCLS))IM_IT = 0
           END DO
           IF(IM_IT.EQ.1) THEN
            JJOCCLS = JOCCLS
            GOTO 101
           END IF
          END DO! JOCCLS
  101     CONTINUE
*
          IF(IM_IT.EQ.1) THEN
           N_ABSPGP_FOR_OCCLS(JJOCCLS) = N_ABSPGP_FOR_OCCLS(JJOCCLS) + 1
           N_ABSPGP(ga_nodeid()+1) = N_ABSPGP(ga_nodeid()+1) + 1
cVOG       N_ABSPGP_TOT = N_ABSPGP_TOT + 1
cVOG       IF(IONLY_NABSPGP.EQ.0) THEN
cVOG         IAB_COMB = IB_ABSPGP_FOR_OCCLS(JJOCCLS)-1+
cVOG &                   N_ABSPGP_FOR_OCCLS(JJOCCLS)
cVOG         I_ABSPGP_FOR_OCCLS(1,IAB_COMB) = JASPGP
cVOG         I_ABSPGP_FOR_OCCLS(2,IAB_COMB) = JBSPGP
cVOG       END IF
          ELSE
*. Supergroup not belonging to an occupation class- may happen
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' STEP 1 IN ABSPGP_TO_OCCLS '
             WRITE(6,*) ' Occupation class not found for supergroup '
             WRITE(6,*) ' JASPGP, JBSPGP = ', JASPGP, JBSPGP
             WRITE(6,*) ' occupation: '
             CALL IWRTMA(JAJBOC,1,NGAS,1,NGAS)
           END IF
          END IF! IM_IT
        END IF ! I_AM_IN
       END DO
      END DO
*
*. Gather IOCCLMAP (occupation map), N_ABSPGP, N_ABSPGP_FOR_OCCLS
*. from all processors
*
      call ga_sync()
      IF ((ga_nodeid()+1).ne.NNODES) THEN
!       DO IOCCL = 1, NOCCLS
!        call ga_put(IOCCLMAP,IOCCL,IOCCL,ga_nodeid()+2,ga_nodeid()+2,
!     &              N_ABSPGP_FOR_OCCLS(IOCCL),1)
!       END DO
       call ga_put(IOCCLMAP,1,NOCCLS,ga_nodeid()+2,ga_nodeid()+2,
     &             N_ABSPGP_FOR_OCCLS(1),NOCCLS)
      END IF
      call ga_igop(1,N_ABSPGP,NNODES,'+')
      call ga_igop(1,N_ABSPGP_FOR_OCCLS,NOCCLS,'+')
      if (ga_nodeid().eq.0.and.NTEST.GE.100)
     &    write(*,*) 'N_ABSPGP_FOR_OCCLS:', N_ABSPGP_FOR_OCCLS
*
*. Update the columns with the elements of the previous columns
*
      call ga_sync()
! BERT THIS IS A GATHER THAT CAN BE DONE MORE EFFICIENTLY!!!!
      DO ICOL = 1, NNODES-1
       call ga_add_patch(1, IOCCLMAP, 1, NOCCLS, ICOL,   ICOL,
     &                   1, IOCCLMAP, 1, NOCCLS, ICOL+1, ICOL+1,
     &                      IOCCLMAP, 1, NOCCLS, ICOL+1, ICOL+1)
      END DO
*
      call ga_sync()
      if (NTEST.GE.100) then
       call ga_print(IOCCLMAP)
       call ga_print_distribution(IOCCLMAP)
      end if
*
*. Calculate total number of AB supergroups
*
      DO I = 1, NNODES
       N_ABSPGP_TOT = N_ABSPGP_TOT + N_ABSPGP(I)
      END DO
*
*. Distribute data before constructing the AB combinations
*
      call dfill(NOCCLS, 0.0d0, N_ABSPGP_VOG, 1)
      DO I = 1, N_ABSPGP_TOT
       I_ABSPGP_A(I) = 0
       I_ABSPGP_B(I) = 0
      END DO
      call ga_get(IOCCLMAP,1,NOCCLS,ga_nodeid()+1,ga_nodeid()+1,
     &            N_ABSPGP_VOG,1)
      call ga_sync()
*
*. STEP 2
*. Loop for constructing I_ABSPGP_A/I_ABSPGP_B --> I_ABSPGP_FOR_OCCLS
*
      DO JASPGP = ilo,ihi
       DO JBSPGP = 1, NBSPGP
        I_AM_IN = IS_ABSPGP_IN_CI_SPACE(
     &            IASPGP(1,JASPGP),IBSPGP(1,JBSPGP),NGAS,
     &            1,1,IGSOCC_MNMX,NELFGP)
        IF(I_AM_IN.EQ.1) THEN
          DO IGAS = 1, NGAS
            JAJBOC(IGAS) =
     &      NELFGP(IASPGP(IGAS,JASPGP)) + NELFGP(IBSPGP(IGAS,JBSPGP))
          END DO
          IM_IT = 3006
          DO JOCCLS = 1, NOCCLS
           IM_IT = 1
           DO IGAS = 1, NGAS
             IF(JAJBOC(IGAS).NE.IOCCLS(IGAS,JOCCLS))IM_IT = 0
           END DO
           IF(IM_IT.EQ.1) THEN
            JJOCCLS = JOCCLS
            GOTO 102
           END IF
          END DO! JOCCLS
  102     CONTINUE
*
          IF(IM_IT.EQ.1) THEN
           N_ABSPGP_VOG(JJOCCLS) = N_ABSPGP_VOG(JJOCCLS) + 1
            
           IF(IONLY_NABSPGP.EQ.0) THEN
             IAB_COMB = IB_ABSPGP_FOR_OCCLS(JJOCCLS)-1+
     &                   N_ABSPGP_VOG(JJOCCLS)
             I_ABSPGP_A(IAB_COMB) = JASPGP
             I_ABSPGP_B(IAB_COMB) = JBSPGP
cVOG         I_ABSPGP_FOR_OCCLS(1,IAB_COMB) = JASPGP
cVOG         I_ABSPGP_FOR_OCCLS(2,IAB_COMB) = JBSPGP
           END IF
          ELSE
*. Supergroup not belonging to an occupation class- may happen
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' STEP 2 IN ABSPGP_TO_OCCLS '
             WRITE(6,*) ' Occupation class not found for supergroup '
             WRITE(6,*) ' JASPGP, JBSPGP = ', JASPGP, JBSPGP
             WRITE(6,*) ' occupation: '
             CALL IWRTMA(JAJBOC,1,NGAS,1,NGAS)
           END IF
          END IF! IM_IT
        END IF ! I_AM_IN
       END DO
      END DO
*
*. Collect I_ABSPGP_A/I_ABSPGP_B data from all processors
*. and store them in I_ABSPGP_FOR_OCCLS
*
      call ga_sync()
      call ga_igop(1,I_ABSPGP_A,N_ABSPGP_TOT,'+')
      call ga_igop(1,I_ABSPGP_B,N_ABSPGP_TOT,'+')
      DO I = 1, N_ABSPGP_TOT
       I_ABSPGP_FOR_OCCLS(1,I) = I_ABSPGP_A(I)
       I_ABSPGP_FOR_OCCLS(2,I) = I_ABSPGP_B(I)
      END DO
*
      if (.not. ga_destroy(IOCCLMAP)) call
     &     errquit('ABSPG: ga_destroy IOCCLMAP failed', 0 , GA_ERR)
*
      IF(ga_nodeid().eq.0.and.NTEST.GE.1000) THEN
       write(*,*) 'NOCCLS', NOCCLS
       write(*,*) 'I_ABSPGP_FOR_OCCLS(1,*)'       
       write(*,*)  (I_ABSPGP_FOR_OCCLS(1,I), I = 1, N_ABSPGP_TOT)
       write(*,*) 'I_ABSPGP_FOR_OCCLS(2,*)'
       write(*,*)  (I_ABSPGP_FOR_OCCLS(2,I), I = 1, N_ABSPGP_TOT)
       write(*,*) 'N_ABSPGP_FOR_OCCLS'
       write(*,*)  N_ABSPGP_FOR_OCCLS 
       write(*,*) 'N_ABSPGP_TOT'
       write(*,*)  N_ABSPGP_TOT
       IF(IONLY_NABSPGP.EQ.0) THEN
        write(*,*) 'AB COMBINATIONS:'
        do i=1, N_ABSPGP_TOT
         write(*,*) I, I_ABSPGP_FOR_OCCLS(1,I),I_ABSPGP_FOR_OCCLS(2,I)
        end do
       END IF
      END IF
*
      IF(IONLY_NABSPGP.EQ.1) THEN
* Construct pointer array
        CALL ZBASE(N_ABSPGP_FOR_OCCLS,IB_ABSPGP_FOR_OCCLS,NOCCLS)
      END IF
*
      IF(ga_nodeid().eq.0.AND.NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' AB supergroups for the occupation classes '
        WRITE(6,*) ' =========================================='
        WRITE(6,*)
        WRITE(6,*) ' Number of AB supergroups in expansion ', 
     &  N_ABSPGP_TOT
        DO JOCCLS = 1, NOCCLS
          WRITE(6,*) ' Info on occupation class: ', JOCCLS
          NAB = N_ABSPGP_FOR_OCCLS(JOCCLS)
          WRITE(6,*) ' Number of AB supergroups: ', NAB
          IF(NTEST.GE.100.AND.IONLY_NABSPGP.EQ.0) THEN
           IB_AB = IB_ABSPGP_FOR_OCCLS(JOCCLS)
           DO IAB = IB_AB, IB_AB - 1 + NAB
            JA = I_ABSPGP_FOR_OCCLS(1,IAB)
            JB = I_ABSPGP_FOR_OCCLS(2,IAB)
            WRITE(6,'(A,2I9)')  '   ', JA,  JB
           END DO
          END IF
        END DO
      END IF !NTEST
      call ga_sync()
*
      CALL LUCIAQEXIT('ABSPG')
      RETURN
      END
      FUNCTION IS_ABSPGP_IN_CI_SPACE(IASPGP,IBSPGP,NGAS,
     &         NMNMX_SPC,IMNMX_SPC,MNMX_OCC,NELFGP)
*
* Is combination of IASPGP and IBSPGP in CI space ?
*
*. Jeppe Olsen, Feb. 2012
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
#include "global.fh"
*. Input
      INTEGER IASPGP(NGAS),IBSPGP(NGAS)
      INTEGER NELFGP(*)
*. And the CI space
      INTEGER MNMX_OCC(MXPNGAS,2,*), IMNMX_SPC(NMNMX_SPC)
*
      NTEST = 000
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Output from IS_ABSPGP_IN_CI_SPACE '
        WRITE(6,*) ' =================================='
      END IF
*
      INCLUDE = 0
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) WRITE(6,*) ' MNMX_SPC '
      DO ISPC = 1, NMNMX_SPC
        IISPC = IMNMX_SPC(ISPC)
        IF(NTEST.GE.1000.and.ga_nodeid().eq.0)
     &   WRITE(6,*) ' IISPC = ', IISPC
        IEL = 0
        IAMOKAY = 1
        DO IGAS = 1, NGAS
          IEL = IEL + NELFGP(IASPGP(IGAS))+NELFGP(IBSPGP(IGAS))
          IF(NTEST.GE.1000.and.ga_nodeid().eq.0) 
     &    WRITE(6,*) ' IGAS, IEL = ', IGAS,IEL
          IF(IEL.LT.MNMX_OCC(IGAS,1,IISPC).OR.
     &       IEL.GT.MNMX_OCC(IGAS,2,IISPC)    )  IAMOKAY = 0
        END DO! Loop over IGAS
        IF(IAMOKAY.EQ.1) INCLUDE = 1
      END DO! Loop over ISPC
*
      IS_ABSPGP_IN_CI_SPACE = INCLUDE
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' IASPGP and IBSPGP '
        CALL IWRTMA(IASPGP,1,NGAS,1,NGAS)
        CALL IWRTMA(IBSPGP,1,NGAS,1,NGAS)
        IF(INCLUDE.EQ.1) THEN
          WRITE(6,*) ' Combination is included '
        ELSE       
          WRITE(6,*) ' Combination is not included '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE ABSPGP_TO_CISPACE(
     &           IASPGP,NASPGP,IBSPGP,NBSPGP,
     &           NGAS,IGSOCC_MNMX,NELFGP,
     &           N_ABSPGP,I_ABSPGP,
     &           IONLY_NABSPGP)
*
* Obtain the combination of A- and B- supergroups that
* belong to a CI space. 
*
* IF IONLY_NABPSGP = 1, Then only the number N_ABSPGP
* (the number of TT combinations in the CI space) is calculated.
*
*
* IGSOCC_MNMX is compound space and is used for pre-screening
* 
*
*. Jeppe Olsen, Feb. 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IASPGP(MXPNGAS,NASPGP),IBSPGP(MXPNGAS,NBSPGP)
      INTEGER IGSOCC_MNMX(MXPNGAS,2)
      INTEGER NELFGP(*)
*. Output
      INTEGER I_ABSPGP(2,*)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN 
        WRITE(6,*) ' ABSPGP_TO_CISPACE reporting '
        WRITE(6,*) ' ============================='
        WRITE(6,*) 
        IF(IONLY_NABSPGP.EQ.1) THEN
          WRITE(6,*) ' Only dimension determined '
        ELSE
          WRITE(6,*) ' Actual mappings determined '
        END IF
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input alpha supergrops '
        DO ISPGP = 1, NASPGP   
          WRITE(6,'(A,I5)') ' Supergroup ', ISPGP
          WRITE(6,'(11X,15(2X,I4))') (IASPGP(IGAS,ISPGP),IGAS = 1, NGAS)
        END DO
        WRITE(6,*) ' Input beta supergrops '
        DO ISPGP = 1, NBSPGP   
          WRITE(6,'(A,I5)') ' Supergroup ', ISPGP
          WRITE(6,'(11X,15(2X,I4))') (IBSPGP(IGAS,ISPGP),IGAS = 1, NGAS)
        END DO
      END IF
*. 
      N_ABSPGP = 0
      DO JASPGP = 1, NASPGP
       DO JBSPGP = 1, NBSPGP
*. Is this combination included in combination
        I_AM_IN = IS_ABSPGP_IN_CI_SPACE(
     &            IASPGP(1,JASPGP),IBSPGP(1,JBSPGP),NGAS,
     &            1,1,IGSOCC_MNMX,NELFGP)
        IF(I_AM_IN.EQ.1) THEN
          N_ABSPGP = N_ABSPGP + 1
          IF(IONLY_NABSPGP.EQ.0) THEN
            I_ABSPGP(1,N_ABSPGP) = JASPGP
            I_ABSPGP(2,N_ABSPGP) = JBSPGP
          END IF
        END IF !I_AM_IN
       END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from ABSPGP_TO_CISPACE '
        WRITE(6,*) ' =============================='
        WRITE(6,*)
        WRITE(6,*) ' Number of AB supergroups in CI-space', 
     &  N_ABSPGP 
      END IF
*
      IF(NTEST.GE.100.AND.IONLY_NABSPGP.EQ.0) THEN
        WRITE(6,*) ' The AB supergroups in the expansion'
        WRITE(6,*) ' ==================================='
        WRITE(6,*)
        WRITE(6,*) 
     &  ' AB number  A-supergroup    B-supergroup '
        WRITE(6,*) 
     &  '*****************************************'
        DO IAB = 1, N_ABSPGP
          WRITE(6,'(1H , I8, 2(5X,I8))')
     &    IAB, I_ABSPGP(1,IAB),I_ABSPGP(2,IAB)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE OCCLS_TO_SPACE_DIM(NCN_FOR_OCCLS,NCS_FOR_OCCLS,
     &           NCM_FOR_OCCLS,NSD_FOR_OCCLS,NCN_ALLSYM_FOR_OCCLS,
     &           IOCCLS_ACT)
*
* Obtain the dimension of the CI spaces from the
* dimension of the occupation classes.
*
*. Jeppe Olsen, Feb. 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER NCN_FOR_OCCLS(NIRREP,*),NCS_FOR_OCCLS(NIRREP,*),
     &        NCM_FOR_OCCLS(NIRREP,*),NSD_FOR_OCCLS(NIRREP,*),
     &        NCN_ALLSYM_FOR_OCCLS(*)
*. Scratch
      INTEGER IOCCLS_ACT(*)
*
      NTEST = 00
*. Output is the N*_PER_SYM_GN arrays in spinfo
*
      DO ICISPC = 1, NCMBSPC
*. Set up in IOCCLS_ACT an array giving allowed occ classes for this CI space
          CALL OCCLS_IN_CISPACE(NOCCLS_ACT,IOCCLS_ACT,
     &         NOCCLS_MAX,int_mb(KIOCCLS),NGAS,
     &         LCMBSPC(ICISPC),ICMBSPC(1,ICISPC),IGSOCCX,ICISPC)
       DO ISM = 1, NIRREP
        NCM_PER_SYM_GN(ISM,ICISPC) = 0
        NSD_PER_SYM_GN(ISM,ICISPC) = 0
*. This is line 250.000 in LUCIA!!(Feb. 7, Jyvaeskylae, couch in Roberts office)
        NCSF_PER_SYM_GN(ISM,ICISPC) = 0
        NCONF_PER_SYM_GN(ISM,ICISPC) = 0 
        IF(ISM.EQ.1) NCONF_ALL_SYM_GN(ICISPC) = 0
        DO IIOCCLS = 1, NOCCLS_ACT
*. Construct in KSIOCCLS_ACT the active classes for this CI space
          IOCCLS = IOCCLS_ACT(IIOCCLS)
*
          NCM_PER_SYM_GN(ISM,ICISPC) = NCM_PER_SYM_GN(ISM,ICISPC)  +
     &    NCM_FOR_OCCLS(ISM,IOCCLS)
          NSD_PER_SYM_GN(ISM,ICISPC) =   NSD_PER_SYM_GN(ISM,ICISPC) +
     &    NSD_FOR_OCCLS(ISM,IOCCLS)
          NCSF_PER_SYM_GN(ISM,ICISPC) = NCSF_PER_SYM_GN(ISM,ICISPC) +
     &    NCS_FOR_OCCLS(ISM,IOCCLS)  
          NCONF_PER_SYM_GN(ISM,ICISPC) = NCONF_PER_SYM_GN(ISM,ICISPC) +
     &    NCN_FOR_OCCLS(ISM,IOCCLS) 
          IF(ISM.EQ.1) NCONF_ALL_SYM_GN(ICISPC) = 
     &    NCONF_ALL_SYM_GN(ICISPC) + NCN_ALLSYM_FOR_OCCLS(IOCCLS)
        END DO !over occlases
       END DO! over symmetries
      END DO! over CI spaces
*
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Information on dimensions per CI space '
       DO ICISPC = 1, NCMBSPC
         WRITE(6,*) 
         WRITE(6,'(A,I3)') ' Info for CIspace ', ICISPC
         WRITE(6,*) ' ============================ '
         WRITE(6,'(A,8(2X,I9))') 
     &   ' Confs: ', (NCONF_PER_SYM_GN(ISM,ICISPC),ISM = 1, NIRREP) 
         WRITE(6,'(A,8(1X,I10))') 
     &   ' CSFs   ', (NCSF_PER_SYM_GN(ISM,ICISPC),ISM = 1, NIRREP) 
         WRITE(6,'(A,8(1X,I10))') 
     &   ' CMs    ', (NCM_PER_SYM_GN(ISM,ICISPC),ISM = 1, NIRREP) 
         WRITE(6,'(A,8(1X,I10))') 
     &   ' SDs    ', (NSD_PER_SYM_GN(ISM,ICISPC),ISM = 1, NIRREP) 
         WRITE(6,'(A,1X,I10)')
     &   ' Configurations for all sym ', NCONF_ALL_SYM_GN(ICISPC)
       END DO
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_CNF_INFO_FOR_OCCLS(IOCCLS_NUM,
     &           IDOSDREO,ISYM)
* 
*. Generate information about configurations and SD reorder 
*. array for occupation class number IOCCLS_NUM
*
*. Jeppe Olsen, Feb. 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
c     INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'spinfo.inc'  
      INCLUDE 'lucinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Scratch
      DIMENSION  IOCCLS_OCC(MXPNGAS)
*
      CALL LUCIAQENTER('GTCNF')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' Information about CONformation for occlass number ',
     & IOCCLS_NUM
       WRITE(6,*) ' ================================================= '
       WRITE(6,*)
       WRITE(6,*) ' Symmetry = ', ISYM
      END IF
*
*. Obtain occupation of occupation class
* ======================================
      CALL MEMCHK2('BEFGOO')
      CALL GET_OCCLS_OCC_FOR_NUMB(IOCCLS_OCC,IOCCLS_NUM)
C          GET_OCCLS_OCC_FOR_NUMB(IOCCLS_OCC,IOCCLS_NUM)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Occupations of occupation class '
        CALL IWRTMA(IOCCLS_OCC,1,NGAS,1,NGAS)
      END IF
      CALL MEMCHK2('AFTGOO')
*
*. Copy Dimension and offsets arrays to active arrays
* =====================================================
      IOCCLS_ADD = (IOCCLS_NUM-1)*NIRREP*(MAXOP+1)
     &            +(ISYM-1)*(MAXOP+1) + 1
C?    WRITE(6,*) ' IOCCLS_ADD = ', IOCCLS_ADD
      CALL ICOPVE2(int_mb(KNCN_PER_OP_SM),IOCCLS_ADD,MAXOP+1,
     &             NCN_FOR_OC_OP_ACT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' NCN_FOR_OC_OP_ACT '
        CALL IWRTMA(NCN_FOR_OC_OP_ACT,1,MAXOP+1,1,MAXOP+1)
      END IF
C?    WRITE(6,*) ' First 3 elements of WORK(KNCN_PER_OP_SM)'
C?    CALL IWRTMA(WORK(KNCN_PER_OP_SM),1,3,1,3)
C     NCNF_TO_NCOMP_PER_OP(MAXOP,NCONF_PER_OPEN,NCOMP_PER_OPEN,NCOMP)
      CALL NCNF_TO_NCOMP_PER_OP(MAXOP,NCN_FOR_OC_OP_ACT,NPCSCNF,
     &                   NCS_FOR_OC_OP_ACT)
      CALL NCNF_TO_NCOMP_PER_OP(MAXOP,NCN_FOR_OC_OP_ACT,NPCMCNF,
     &                   NCM_FOR_OC_OP_ACT)
      CALL NCNF_TO_NCOMP_PER_OP(MAXOP,NCN_FOR_OC_OP_ACT,NPDTCNF,
     &                   NSD_FOR_OC_OP_ACT)
      CALL ZBASE(NCN_FOR_OC_OP_ACT,IBCN_FOR_OC_OP_ACT,MAXOP+1)
      CALL ZBASE(NCM_FOR_OC_OP_ACT,IBCM_FOR_OC_OP_ACT,MAXOP+1)
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' NCM_FOR_OC_OP, IBCM_FOR_OC... '
        CALL IWRTMA(NCM_FOR_OC_OP_ACT,1,MAXOP+1,1,MAXOP+1)
        CALL IWRTMA(IBCM_FOR_OC_OP_ACT,1,MAXOP+1,1,MAXOP+1)
      END IF
*.The following contains some overwriting that I am not fond of...
      CALL ICOPVE(NCN_FOR_OC_OP_ACT,NCONF_PER_OPEN(1,ISYM),
     &            MAXOP+1)
      CALL ICOPVE(IBCN_FOR_OC_OP_ACT,IB_CN_OPEN,MAXOP+1)
      CALL ICOPVE(IBCM_FOR_OC_OP_ACT,IB_CM_OPEN,MAXOP+1)
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' IB_CM_OPEN '
        CALL IWRTMA(IB_CM_OPEN,1,MAXOP+1,1,MAXOP+1)
      END IF
      CALL ICOPVE2(int_mb(KNCN_ALLSYM_FOR_OCCLS),IOCCLS_NUM,1,
     &     NCONF_ALL_SYM_GN(1))
*
      IATP = 1
      IBTP = 2
      NEL = NELFTP(IATP) + NELFTP(2)
C?    WRITE(6,*) ' NEL = ', NEL
C     ZIB_CONFOCC(N_CONF_FOR_OP,NEL,IB_CONF_FOR_OP,MAXOP)
      CALL ZIB_CONFOCC(NCN_FOR_OC_OP_ACT,NEL,IB_CNOCC_OPEN,
     &                 MAXOP)

*  
*  Configuration information
* ===========================
*. Storage of information:
* Arcweight matrix Z for configurations: KZCONF
* Occupation of Configuration: KICONF_OCC(ISYM)
* Reorder array of configurations: KICONF_REO
C          GEN_CONF_FOR_CISPC(IOCCLS,NOCCLS,ISYM)
      CALL GEN_CONF_FOR_CISPC(IOCCLS_NUM,1,ISYM,int_mb(KIOCCLS))
*  
*. Reorder array of determinants betweeen conf and string order
* =============================================================
      IF(IDOSDREO.EQ.1) THEN
        CALL REO_SD_FOR_OCCLS(IOCCLS_NUM,ISYM,int_mb(KSDREO_I(ISYM)))   
      END IF
*
      CALL LUCIAQEXIT('GTCNF')
*
      RETURN
      END
      SUBROUTINE REO_SD_FOR_OCCLS(IOCCLS_NUM,ISYM,IREO)
*
*
* For the determinants of occupation class defined by occupation class number 
* IOCCLS, determine the reorder array of going between the CONF and 
* AB ordering.
*
* It is assumed that all information about prototype confs
* and the actual occupation class has been obtained
*
* Jeppe Olsen Feb. 2012, from CNFORD_GAS
*             Completed in Geneva, Feb. 2012
*
* Last Modification: Jeppe Olsen; April 28, 2013, Minneapolis
*                    Improved reorder inserted
*                    Jeppe Olsen; April 29, 2013; Minneapolis
*                    Even more improved version added
*
*
* =====
*.Input
* =====
*
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
      INCLUDE 'spinfo.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOCCLS(NGAS)
*. Output 
      DIMENSION IREO(*)
*
      CALL LUCIAQENTER('REOOC')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'REO_OC')
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' Output from REO_SD_FOR_OCCLS'
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Number of occupation class in action ', 
     &             IOCCLS_NUM
      END IF
*
* Initial preparations
*
* Type of alpha and beta strings
      IATP = 1             
      IBTP = 2              
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL= NAEL+NBEL
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
*. Configurations are constructed and stored with orbitals
*. labels starting with 1, whereas strings are stored
*. with orbitals starting with NINOB+1. Provide offset  
      IB_ORB = NINOB+1
*
*. Approach
*
      INEW_OR_OLD = 1
      IF(NGAS.EQ.1) THEN
        INEW_OR_SUPERNEW = 1
      ELSE
        INEW_OR_SUPERNEW = 2
      END IF
*
* Scratch space 
*
*. For alpha and beta strings
      IF(INEW_OR_SUPERNEW.NE.2) THEN
*. Allocate strings of given sym and type
        CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')  !done
        CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')  !done
      ELSE
*
*. Allocate arrays for occupations of symmetries for P and L strings,
*. Could be reduced a lot- as it is not complete strings
*. we need to store
       CALL MEMMAN(KLASTRP,MXNSTRP_AS*NAEL,'ADDL  ',1,'ASTP  ')  !done
       CALL MEMMAN(KLBSTRP,MXNSTRP_AS*NBEL,'ADDL  ',1,'BSTP  ')  !done
       CALL MEMMAN(KLASTRL,MXNSTRL_AS*NAEL,'ADDL  ',1,'ASTL  ')  !nu
       CALL MEMMAN(KLBSTRL,MXNSTRL_AS*NBEL,'ADDL  ',1,'BSTL  ')  !nu
*. Contribution to adress of P-configurations
       CALL MEMMAN(KLCNSTSTP,MXNSTRSTRP_AS,'ADDL  ',1,'SSAP  ')  !done
*. Occupation of L-strings
       CALL MEMMAN(KLCNOCSTSTL,MXNSTRSTRLOC_AS,'ADDL  ',1,'LCNOC ') !d.
*. Contribution to adress of L-configurations
       CALL MEMMAN(KLCNSTSTL,MXNSTRSTRL_AS,'ADDL  ',1,'SSAL  ') !done
*. Number of open orbitals in P and L-strings
       CALL MEMMAN(KLNOPSTSTP,MXNSTRSTRP_AS,'ADDL  ',1,'SSNOPP') !done
       CALL MEMMAN(KLNOPSTSTL,MXNSTRSTRL_AS,'ADDL  ',1,'SSNOPL') !done
*. Number of open alpha-electrons in P-and L-strings
       CALL MEMMAN(KLNALSTSTP,MXNSTRSTRP_AS,'ADDL  ',1,'SSNALP') !done
       CALL MEMMAN(KLNALSTSTL,MXNSTRSTRL_AS,'ADDL  ',1,'SSNALP') !done
*. P-contribution to address of spin-combination
       CALL MEMMAN(KLIADOPSTSTP,MXNSTRSTRP_AS,'ADDL  ',1,'SSIOPP')
*. Are P-combinations AB switched
       IF(PSSIGN.NE.0.0D0) THEN
         CALL MEMMAN(KLABSWITCHP,MXNSTRSTRP_AS,'ADDL  ',1,'SSAPSW') !d
         CALL MEMMAN(KLABSWITCHL,MXNSTRSTRL_AS,'ADDL  ',1,'SSALSW') !d
         CALL MEMMAN(KLNOPSTSTL_AB,MXNSTRSTRL_AS,'ADDL  ',1,'SSNOPL')!nu
         CALL MEMMAN(KLCNOCSTSTL_AB,MXNSTRSTRLOC_AS,'ADDL  ',1,'LCNOC ')
         CALL MEMMAN(KLCNSTSTL_AB,MXNSTRSTRL_AS,'ADDL  ',1,'SSAL  ') !d
         CALL MEMMAN(KLNALSTSTL_AB,MXNSTRSTRL_AS,'ADDL  ',1,'SSNALP') !d
       ELSE
         KLABSWITCHP = 1
         KLABSWITCHL = 1
         KLNOPSTSTL_AB = 1
         KLCNOCSTSTL_AB = 1
         KLCNSTSTL_AB = 1
         KLNALSTSTL_AB = 1
       END IF
*. Offsets in STRSTR arrays to combinations with given symmetries
       CALL MEMMAN(KLSMSMP,NSMST**2,'ADDL  ',1,'SMSMP ')  !done
       CALL MEMMAN(KLSMSML,NSMST**2,'ADDL  ',1,'SMSML ')  !done
      END IF
*. For Occupation and projections of a given determinant
      CALL MEMMAN(KLDET_OC,NAEL+NBEL,'ADDL  ',1,'CONF_O')  !done
      CALL MEMMAN(KLDET_MS,NAEL+NBEL,'ADDL  ',1,'CONF_M')  !done
      CALL MEMMAN(KLDET_MS_AB,NAEL+NBEL,'ADDL  ',1,'CONF_M')  !done
      CALL MEMMAN(KLDET_VC,NOCOB,'ADDL  ',1,'CONF_M')  !done
      CALL MEMMAN(KLBLTP,NSMST,'ADDL   ',1,'BLTP  ')  !done
*. 
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,int_mb(KLBLTP),int_mb(KSVST))
      NCMB_CN = int_mb(KNCM_FOR_OCCLS+(IOCCLS_NUM-1)*NSMST+ISYM-1)
C?    WRITE(6,*) ' TESTY, NCMB_CN = ', NCMB_CN
      ILCHK = -2303
      IF(INEW_OR_OLD.EQ.2) THEN
        CALL REO_SD_FOR_OCCLS_S(IOCCLS_NUM,ISYM,IREO,
     &       NOCTPA,NOCTPB,IOCTPA,IOCTPB,NAEL,NBEL,
     &       NSMST,NGAS,IB_ORB,NACOB,NOCOB,PSSIGN,MINOP,
     &       NTOOB,NOBPT,int_mb(KLBLTP),ISMOST(1,ISYM),
     &       int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &       int_mb(KICONF_REO(1)),int_mb(KZCONF),
     &       NPCMCNF,int_mb(KDFTP),
     &       KZ_PTDT,KREO_PTDT,
     &       IB_CN_OPEN, IB_CM_OPEN,
     &       int_mb(KNABSPGP_FOR_OCCLS),int_mb(KIBABSPGP_FOR_OCCLS),
     &       int_mb(KIABSPGP_FOR_OCCLS),
     &       int_mb(KLASTR),int_mb(KLBSTR),
     &       int_mb(KLDET_OC),int_mb(KLDET_MS),int_mb(KLDET_VC),ILCHK)
      ELSE
       IF(INEW_OR_SUPERNEW.EQ.1) THEN
       CALL REO_SD_FOR_OCCLS_SN(IOCCLS_NUM,ISYM,IREO,
     &        NOCTPA,NOCTPB,IOCTPA,IOCTPB,NAEL,NBEL,
     &        NSMST,NGAS,IB_ORB,NACOB,NOCOB,PSSIGN,MINOP,
     &        NTOOB,NOBPT,int_mb(KLBLTP),ISMOST(1,ISYM),
     &        int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &        int_mb(KICONF_REO(1)),int_mb(KZCONF),
     &        NPCMCNF,int_mb(KDFTP),
     &        KZ_PTDT,KREO_PTDT,
     &        IB_CN_OPEN, IB_CM_OPEN,
     &        int_mb(KNABSPGP_FOR_OCCLS),int_mb(KIBABSPGP_FOR_OCCLS),
     &        int_mb(KIABSPGP_FOR_OCCLS),
     &        int_mb(KLASTR),int_mb(KLBSTR),
     &        int_mb(KLDET_OC),int_mb(KLDET_MS),int_mb(KLDET_VC),ILCHK)
       ELSE 
*. The supernew approach
       CALL REO_SD_FOR_OCCLS_SSN(IOCCLS_NUM,ISYM,IREO,
     &        NOCTPA,NOCTPB,IOCTPA,IOCTPB,NAEL,NBEL,
     &        NSMST,NGAS,IB_ORB,NACOB,NOCOB,PSSIGN,MINOP,
     &        NTOOB,NOBPT,int_mb(KLBLTP),ISMOST(1,ISYM),
     &        int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &        int_mb(KISTSO(IATP)),int_mb(KISTSO(IBTP)),
     &        int_mb(KICONF_REO(1)),int_mb(KZCONF),
     &        NPCMCNF,int_mb(KDFTP),
     &        KZ_PTDT,KREO_PTDT,
     &        IB_CN_OPEN, IB_CM_OPEN,
     &        int_mb(KNABSPGP_FOR_OCCLS),int_mb(KIBABSPGP_FOR_OCCLS),
     &        int_mb(KIABSPGP_FOR_OCCLS),
     &        int_mb(KLASTRP),int_mb(KLBSTRP),
     &        int_mb(KLDET_OC),int_mb(KLDET_MS),int_mb(KLDET_VC),
     &        int_mb(KLDET_MS_AB),
     &        ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &        int_mb(KNSTSGP(1)),int_mb(KISTSGP(1)),
     &        int_mb(KLCNSTSTP), int_mb(KLCNSTSTL), int_mb(KLSMSMP), 
     &        int_mb(KLSMSML),int_mb(KLCNOCSTSTL),
     &        int_mb(KLCNOCSTSTL_AB),
     &        int_mb(KLCNSTSTL_AB), int_mb(KLNALSTSTL_AB),
     &        int_mb(KLNOPSTSTP),int_mb(KLNOPSTSTL),
     &        int_mb(KLNALSTSTP),int_mb(KLNALSTSTL),
     &        int_mb(KLIADOPSTSTP),int_mb(KLABSWITCHP),
     &        int_mb(KLABSWITCHL),
     &        MAXOP,NCMB_CN,IDC,ILCHK)
       END IF ! Switch between new and supernew
      END IF !Switch between new or old
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'REO_OC')
      CALL LUCIAQEXIT('REOOC')
*
      RETURN
      END
      SUBROUTINE REO_SD_FOR_OCCLS_S(IOCCLS_NUM,ISYM,IREO,
     &           NOCTPA,NOCTPB,IOCTPA,IOCTPB,
     &           NAEL,NBEL,
     &           NSMST,NGAS,IB_ORB,NACOB,NOCOB,PSSIGN,MINOP,
     &           NTOOB,NOBPT,IBLTP,ISMOST,
     &           NSSOA,NSSOB,
     &           ICONF_REO,IZCONF,
     &           NPCMCNF,DFTP,
     &           KZ_PTDT,KREO_PTDT,
     &           IB_CN_OPEN, IB_CM_OPEN,
     &           NABSPGP_FOR_OCCLS,IBABSPGP_FOR_OCCLS,
     &           IABSPGP_FOR_OCCLS,
     &           IASTR,IBSTR,
     &           IDET_OC,IDET_MS,IDET_VC,
     &           IASPGPOC,IBSPGPOC,ILCHK)
*
* Reorder determinants in GAS space from det to configuration order
* for determinants in OCClass IOCCLS
* Reorder array created is Conf-order => AB-order 
*
c     IMPLICIT REAL*8(A-H,O-Z)
      include 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. General input
*.---------------
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)  
      INTEGER IBLTP(*), ISMOST(*)
      INTEGER NOBPT(*)
      INTEGER DFTP(*) 
      INTEGER NPCMCNF(*)
*. Info on the AB supergroups
      INTEGER NABSPGP_FOR_OCCLS(*), IABSPGP_FOR_OCCLS(2,*),
     &        IBABSPGP_FOR_OCCLS(*)
*. IB_CN_OPEN, IB_CM_OPEN(IOPEN+1) Gives start of confs/CM's
*. with given symmetry and number of open orbitals
      INTEGER IB_CN_OPEN(*), IB_CM_OPEN(*)
*. Info for the lexical adressing of configurations 
      INTEGER ICONF_REO(*), IZCONF(*)
*. WORK(KZ_PTDT(IOPEN+1) gives Z  array for prototype dets with IOPEN 
*. WORK(KREO_PTDT(IOPEN+1) gives the corresponding reorder array
*. open orbitals
      INTEGER KZ_PTDT(*), KREO_PTDT(*)
*. The occupations of the a- and b-super groups
      INTEGER IASPGPOC(MXPNGAS,*), IBSPGPOC(MXPNGAS,*)
*
*. The work array used for WORK(KZ_PTDET()),WORK(KREO_PTDT())
*. Scratch space 
*. --------------
      INTEGER IASTR(NAEL,*),IBSTR(NBEL,*)
      INTEGER IDET_OC(*), IDET_MS(*) , IDET_VC(*)
*
*. Output
*. ------
      INTEGER IREO(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output from REO_SD_OCCLS_S '
        WRITE(6,*) ' ========================='
        WRITE(6,*) 
        WRITE(6,*) ' Number of occupation class in action: ', IOCCLS_NUM
      END IF
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' ILCHK = ', ILCHK
        WRITE(6,*) ' PSSIGN = ', PSSIGN
        WRITE(6,*) ' IB_ORB = ', IB_ORB
      END IF
*
      IAGRP = 1
      IBGRP = 2
*
      NEL = NAEL + NBEL
*
      IDET = 0
      N_AB = NABSPGP_FOR_OCCLS(IOCCLS_NUM)
      IB_AB = IBABSPGP_FOR_OCCLS(IOCCLS_NUM)
      DO I_AB = IB_AB, IB_AB + N_AB - 1
       IATP = IABSPGP_FOR_OCCLS(1,I_AB)
       IBTP = IABSPGP_FOR_OCCLS(2,I_AB)
       DO IASM = 1, NSMST
        IBSM = ISMOST(IASM)
        IF(IBLTP(IASM).EQ.1.OR.(IBLTP(IASM).EQ.2.AND.IATP.GE.IBTP)) THEN
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         IF(NTEST.GE.10000) THEN
           WRITE(6,'(A,5(2X,I6))')
     &     ' I_AB, IATP, IBTP, IASM, IBSM = ',
     &       I_AB, IATP, IBTP, IASM, IBSM 
         END IF
*
*. Obtain alpha strings of sym IASM and type IATP
*
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                            NORB,0,IDUM,IDUM)
*
*. Obtain Beta  strings of sym IBSM and type IBTP
*
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                            NORB,0,IDUM,IDUM)
*. Offset to this occupation class in occupation class ordered cnf list
         IB_OCCLS = 1
         IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
         ELSE
          IRESTR = 0
         END IF
*
         DO  IB = 1,NIB
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB 
            IRESTR2 = 1
          ELSE
            MINIA = 1
            IRESTR2 = 0
          END IF
          DO  IA = MINIA,NIA
            IDET = IDET + 1
C                ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,IDET_OC,IDET_SP,ISIGN)
            CALL ABSTR_TO_ORDSTR(IASTR(1,IA),IBSTR(1,IB),NAEL,NBEL,
     &                           IDET_OC,IDET_MS,ISIGN)
*. Orbital numbers in strings are absolute, but relative in confs.
*. subtract offset(number of inactive orbitals)
            ISUB = 1-IB_ORB
            CALL IADCONST(IDET_OC,ISUB,NAEL+NBEL)
C     IADCONST(IVEC,IADD, NDIM)
*. Number of open orbitals in this configuration 
            NOPEN = NOP_FOR_CONF(IDET_OC,NEL)
            IF(NOPEN.GE.MINOP) THEN
C                    NOP_FOR_CONF(ICONF,NEL)
             NDOUBLE = (NEL-NOPEN)/2
             NOCOB_L = NOPEN + NDOUBLE
             NOPEN_AL = NAEL - NDOUBLE
C?           WRITE(6,*) ' NOPEN, NOPEN_AL = ', NOPEN,NOPEN_AL
             NPTDT = NPCMCNF(NOPEN+1)
*. Packed form of this configuration 
C                REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY) 
             CALL REFORM_CONF_OCC(IDET_OC,IDET_VC,NEL,NOCOB_L,1)
*. Address of this configuration 
*. Offset to configurations with this number of open orbitals in 
*. reordered cnf list
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' IASTR, IBSTR, ICONF:'
               CALL IWRTMA(IASTR(1,IA),1,NAEL,1,NAEL)
               CALL IWRTMA(IBSTR(1,IB),1,NBEL,1,NBEL)
               CALL IWRTMA(IDET_VC,1,NOCOB_L,1,NOCOB_L)
             END IF
*
             ICNF_OUT = ILEX_FOR_CONF2(IDET_VC,NOCOB_L,NACOB,NEL,IZCONF,
     &                  1,ICONF_REO,1)
C                      ILEX_FOR_CONF2(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,
C                      IREO,IB_OCCLS)
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' Configuration: '
               CALL IWRTMA(IDET_VC,1,NOCOB_L,1,NOCOB_L)
               WRITE(6,*) ' number of configuration in output list',
     &         ICNF_OUT
             END IF
*. Spinprojections of open orbitals
             CALL EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_VC,NEL)
C                 EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_OPEN_MS,NEL)
*. Address of this spinprojection pattern   
C  IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
             IPTDT = IZNUM_PTDT(IDET_VC,NOPEN,NOPEN_AL,
     &              int_mb(KZ_PTDT(NOPEN+1)),int_mb(KREO_PTDT(NOPEN+1)),
     &              1)
             ISIGNP = 1
             IF(IPTDT.EQ.0) THEN
              IF(PSSIGN.NE.0) THEN
*. The determinant was not found among the list of prototype dets. For combinations
*. this should be due to the prototype determinant is the MS- switched determinant, so find
*. address of this and remember sign
                M1 = -1
                CALL ABSTR_TO_ORDSTR(IBSTR(1,IB),IASTR(1,IA),NBEL,NAEL,
     &                            IDET_OC,IDET_MS,ISIGN)
*. Spinprojections of open orbitals
                CALL EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_VC,NEL)
                IPTDT = IZNUM_PTDT(IDET_VC,NOPEN,NOPEN_AL,
     &                  int_mb(KZ_PTDT(NOPEN+1)),
     &                  int_mb(KREO_PTDT(NOPEN+1)),1)
                IF(PSSIGN.EQ.-1.0D0) ISIGNP = -1
              ELSE 
*. Prototype determinant was not found in list
               WRITE(6,*) 
     &         ' Error: Determinant not found in list of protodets'
               WRITE(6,*) 
     &         ' Detected in REO_SD_FOR_OCCLS_S'
              END IF
             END IF
*
             IBCNF_OUT = IB_CN_OPEN(NOPEN+1)
             IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
              WRITE(6,*) ' IB_CM_OPEN(NOPEN+1) = ',
     &                     IB_CM_OPEN(NOPEN+1)
              WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
              WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
             END IF
             IADR_SD_CONF_ORDER = IB_CM_OPEN(NOPEN+1) - 1
     &                          + (ICNF_OUT-IBCNF_OUT)*NPTDT + IPTDT
             IF(IADR_SD_CONF_ORDER.LE.0) THEN
               WRITE(6,*) ' Problemo, IADR_SD_CONF_ORDER < 0 '
               WRITE(6,*) ' IADR_SD_CONF_ORDER = ', IADR_SD_CONF_ORDER
               WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
               WRITE(6,*) ' IB_CM_OPEN(NOPEN+1) = ',
     &                     IB_CM_OPEN(NOPEN+1)
               WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
               WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
C?             CALL XFLUSH(6)
             END IF
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' IADR_SD_CONF_ORDER, ISIGN, IDET = ',
     &                      IADR_SD_CONF_ORDER, ISIGN, IDET
             END IF
             IREO(IADR_SD_CONF_ORDER) = ISIGN*IDET*ISIGNP
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' IDET, IADR_SD_CONF_ORDER ',
     &                      IDET, IADR_SD_CONF_ORDER
             END IF
            END IF! Nopen .ge. MINOP
          END DO
*         ^ End of loop over alpha strings
         END DO
*        ^ End of loop over beta strings
        END IF! Block should be included
       END DO ! Loop over IASM
      END DO! Loop over AB blocks in occ. class
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Reorder array for SDs, CONF order => string order'
        WRITE(6,*) ' ================================================='
        CALL IWRTMA(IREO,1,IDET,1,IDET)
      END IF
*
      CALL MEMCHK2('AFTSSN')
     
*
*. Check sum of reordering array
      I_DO_CHECKSUM = 0
      IF(I_DO_CHECKSUM.EQ.1) THEN
        ISUM = 0
        DO JDET = 1, IDET
          ISUM = ISUM + ABS(IREO(JDET))
        END DO
        IF(ISUM.NE.IDET*(IDET+1)/2) THEN
          WRITE(6,*) ' Problem with sumcheck in REO_GASDET'
          WRITE(6,'(A,2I9)') 
     &    'Expected and actual value ', ISUM, IDET*(IDET+1)/2
          STOP       ' Problem with sumcheck in REO_GASDET'
        ELSE
          WRITE(6,*) ' Sumcheck in REO_GASDET passed '
        END IF
      END IF !checksum is invoked
*
      RETURN
      END
      SUBROUTINE GET_OCCLS_OCC_FOR_NUMB(IOCCLS_OCC,IOCCLS_NUM)
*
* Obtain the occupation of occupation class with (global) number
* IOCCLS_NUM
*
*. Jeppe Olsen, Febr. 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
c     INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      NTEST = 00
*
      CALL ICOPVE2(int_mb(KIOCCLS),(IOCCLS_NUM-1)*NGAS+1,
     $             NGAS,IOCCLS_OCC)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Occupation of occupation class ', IOCCLS_NUM  
       CALL IWRTMA(IOCCLS_OCC,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE TEST_CNF_INFO_FOR_OCCLS
*
* Loop over all occupation classes and generate all conf info
*
*. Jeppe Olsen, testing large scale CSF routines, Geneva, Febr. 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cstate.inc'
*
      DO IOCCLS = 1, NOCCLS_MAX
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCCLS,1,IREFSM)
      END DO
*
      RETURN
      END
      SUBROUTINE ZIB_CONFOCC(N_CONF_FOR_OP,NEL,IB_CONFOCC_FOR_OP,MAXOP)
*
* An expansion defined by N_CONF_FOR_OP is given
* Obtain offsets to lists giving the occupations 
*
*. Jeppe Olsen, Feb. 2012
*
      INCLUDE 'implicit.inc' 
*. Input
      INTEGER N_CONF_FOR_OP(MAXOP+1)
*. Output
      INTEGER IB_CONFOCC_FOR_OP(MAXOP+1)
*
      NTEST = 00
*
      IB_CONFOCC_FOR_OP(1) = 1
      DO NOP = 0, MAXOP 
        NOC = NOP + (NEL-NOP)/2
        IF(NOP.NE.MAXOP) THEN
          IB_CONFOCC_FOR_OP(NOP+2) = 
     &    IB_CONFOCC_FOR_OP(NOP+1)+NOC*N_CONF_FOR_OP(NOP+1)
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' N_CONF_FOR_OP (input) '
        CALL IWRTMA(N_CONF_FOR_OP,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*) ' IB_CONF_FOR_OP (output) '
        CALL IWRTMA(IB_CONFOCC_FOR_OP,1,MAXOP+1,1,MAXOP+1)
      END IF
*
      RETURN
      END
      SUBROUTINE NCNF_TO_NCOMP_PER_OP
     &           (MAXOP,NCONF_PER_OPEN,NCOMP_PER_OPEN,NCOMP_CONF_OPEN)
*
* Number of configurations per number of open orbitals is given
* Find number of some components, defined by NCOMP_PER_OPEN
* for the configuration list per open, 
*
* In practice : components are SD's, CSF's or CMB's
*
* Jeppe Olsen, Feb. 2012  
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NCONF_PER_OPEN(*), NCOMP_PER_OPEN(*)
*. Output
      INTEGER NCOMP_CONF_OPEN(*)
*
      DO IOPEN = 0, MAXOP  
        NCOMP_CONF_OPEN(IOPEN+1) = 
     &  NCONF_PER_OPEN(IOPEN+1)*NCOMP_PER_OPEN(IOPEN+1)
      END DO
*
      RETURN
      END
      SUBROUTINE Z_BASSPC_FOR_ALL_OCCLS(IOCCLS,NOCCLS,IBASSPC)
*
* A set of occupation classes IOCCLS is given. Obtain base space
* for each occupation class
*
*. Jeppe Olsen, Feb. 2012, Geneva
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Specific input
      INTEGER IOCCLS(NGAS,NOCCLS)
*. Output
      INTEGER IBASSPC(NOCCLS)
*
      DO JOCCLS = 1, NOCCLS
        IBASSPC(JOCCLS) = IBASSPC_FOR_CLS(IOCCLS(1,JOCCLS))
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The basespace of the various occupation classes'
        CALL IWRTMA(IBASSPC,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END
      SUBROUTINE CSDTVCMN(CSFVEC,DETVEC,SCR,IWAY,ICOPY,ISYM,ISPC,
     &           IMAXMIN_OR_GAS,ICNFBAT,LU_DET,LU_CSF,NOCCLS_ACT,
     &           IOCCLS_ACT,IBLOCK,NBLK_PER_BATCH)
*
* Outer routine for transformation between CSF and CM forms of 
* vectors. The CM's may either by spin-combinations or determinants
*
* IWAY = 1 => CSF to CM
* IWAY = 2 => CM to CSF
*
* For GAS-expansions (where reorganization of SD's are required),
* the scratch vector SCR of length NCM_STRING is used
*
* Jeppe Olsen, June 2011
*              ICNFBAT option added, Feb. 2012, Geneva
*              1 batch per occupation class assumed
*. Last modification; May 17, 2013; Jeppe Olsen;  Distinction between NCM_CN, NCM_ST
*                     for INFBAT = 2
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'csm.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. For ICNFBAT = 2:
      INTEGER IOCCLS_ACT(NOCCLS_ACT), IBLOCK(8,*), NBLK_PER_BATCH(*)
*. Input and output
      DIMENSION CSFVEC(*), DETVEC(*)
*. Scratch
      DIMENSION SCR(*)
      CALL LUCIAQENTER('CSDTV')
      NTESTL = 000
      NTEST  = MAX(NTESTL,IPRCSF) 
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' CSDTVCMN speaking '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' ICOPY, ISYM, ISPC,IWAY = ', 
     &               ICOPY, ISYM, ISPC, IWAY
        WRITE(6,*) ' IMAXMIN_OR_GAS = ', IMAXMIN_OR_GAS
        WRITE(6,*) ' ICNFBAT = ', ICNFBAT
        WRITE(6,*) ' NOCCLS_ACT = ', NOCCLS_ACT
        WRITE(6,*) ' NBLK_PER_BATCH '
        CALL IWRTMA(NBLK_PER_BATCH,1,NOCCLS_ACT,1,NOCCLS_ACT)
      END IF
*
      IF(ICNFBAT.EQ.1) THEN
*. Proceed as usual- all info and data in core
       NCMB_CM = NCM_PER_SYM_GN(ISYM,ISPC)
       NCSF = NCSF_PER_SYM_GN(ISYM,ISPC) 
*
       IF(IMAXMIN_OR_GAS.EQ.1) THEN
*
*, The expansion is a MAXMIN space 
*
        IF(NTEST.GE.100) WRITE(6,*) ' NCMB, NCSF = ', NCMB, NCSF
*. No reordering
        IDO_REO = 0
        CALL CSDTVCB(CSFVEC,DETVEC,IWAY,dbl_mb(KDTOC),
     &       int_mb(KICONF_REO_GN(ISYM,ISPC)),NCMB,NCMB,NCSF,
     &       NCONF_PER_OPEN_GN(1,ISYM,ISPC),
     &       ICOPY,IDO_REO,NTEST)
C     CSDTVCB(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,
C    &                  NDET,NCSF,NCNFTP,
C    &                  ICOPY,IDO_REO,NTEST)
       ELSE 
*
*. The CI expansion is a standard GAS expansion, requiring reordering
*. Note that the reorder array is assumed in KSDREO_I(ISYM) irrespectively
* of space
        IDO_REO = 1
        IF(IWAY.EQ.1) THEN
*. CSF => SD transformation
         CALL COPVEC(CSFVEC,SCR,NCSF)
         CALL CSDTVCB(SCR,DETVEC,IWAY,dbl_mb(KDTOC),
     &        int_mb(KSDREO_I(ISYM)),NCCM_CONF,NCCM_STRING,NCSF,
     &        NCONF_PER_OPEN(1,ISYM),
     &        ICOPY,IDO_REO,NTEST)
        ELSE
*. SD => CSF transformation
         CALL CSDTVCB(SCR,DETVEC,IWAY,dbl_mb(KDTOC),
     &        int_mb(KSDREO_I(ISYM)),NCCM_CONF,NCCM_STRING,NCSF,
     &        NCONF_PER_OPEN(1,ISYM),
     &        ICOPY,IDO_REO,NTEST)
         CALL COPVEC(SCR,CSFVEC,NCSF)
        END IF! IWAY switch
       END IF! Maxmin switch
      ELSE
*
*. Batching of all info and data
* ==============================
*
       IDO_REO = 1
       CALL REWINO(LU_DET)
       CALL REWINO(LU_CSF)
       IB_BLK = 1
       DO IIOCLS = 1, NOCCLS_ACT
        IOCCLS = IOCCLS_ACT(IIOCLS)
        IF(NTEST.GE.1000) WRITE(6,*) ' Output from IOCCLS = ', IOCCLS
        NBLK = NBLK_PER_BATCH(IIOCLS)
        IF(NTEST.GE.1000) WRITE(6,*) ' NBLK = ', NBLK
*. Generate Conformation
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCCLS,1,ISYM)
        NCSF_OCCLS = IELSUM(NCS_FOR_OC_OP_ACT,MAXOP+1)
        NCM_OCCLS = IELSUM(NCM_FOR_OC_OP_ACT,MAXOP+1)
        ncm_occls_ab=int_mb(kncmab_for_occls+(ioccls-1)*nsmst+isym-1)
c       NCM_OCCLS_AB = 
c    &  IFRMR(WORK,KNCMAB_FOR_OCCLS,(IOCCLS-1)*NSMST+ISYM)
C?      WRITE(6,*) ' TESTY, NCM_OCCLS_AB = ', NCM_OCCLS_AB
        IF(IWAY.EQ.1) THEN
*. CSF => SD transformation
*. Obtain CSF- coefficients
         CALL FRMDSCN(SCR,1,-1,LU_CSF)
         CALL CSDTVCB(SCR,DETVEC,IWAY,dbl_mb(KDTOC),
     &        int_mb(KSDREO_I(ISYM)),NCM_OCCLS,NCM_OCCLS_AB,
     &        NCSF_OCCLS,NCN_FOR_OC_OP_ACT,
     &        ICOPY,IDO_REO,NTEST)
         IF(NTEST.GE.1000) WRITE(6,*) ' Home from CSDTVCB '
*. Write determinants to disc blockwise
         IOFF = 1
         DO IBLK = IB_BLK, IB_BLK+NBLK -1
           LENP = IBLOCK(8,IBLK)
           IF(NTEST.GE.1000) 
     &     WRITE(6,*) ' IBLK, LENP = ', IBLK, LENP
           CALL TODSCN(DETVEC(IOFF),1,LENP,-1,LU_DET)
           IOFF = IOFF + LENP
         END DO
         IF(NTEST.GE.1000) WRITE(6,*) ' After TODSCN '
        ELSE
*. SD => CSF transformation
*. Read determinant blocks in
         CALL FRMDSCN(DETVEC,NBLK,-1,LU_DET)
         CALL CSDTVCB(SCR,DETVEC,IWAY,dbl_mb(KDTOC),
     &        int_mb(KSDREO_I(ISYM)),NCM_OCCLS,NCM_OCCLS_AB,
     &        NCSF_OCCLS,NCN_FOR_OC_OP_ACT,
     &        ICOPY,IDO_REO,NTEST)
*. Write CSF coefs to disc
         CALL TODSCN(SCR,1,NCSF_OCCLS,-1,LU_CSF)
        END IF! IWAY switch
        IB_BLK = IB_BLK + NBLK
       END DO! loop over occupation classes
* EOF marks
       IF(IWAY.EQ.1) THEN
         CALL ITODS(-1,1,-1,LU_DET)
       ELSE
         CALL ITODS(-1,1,-1,LU_CSF)
       END IF
      END IF! ICNFBAT switch
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Output CSF and CM vectors from CSDTVCMN' 
       WRITE(6,*) ' ======================================='
       WRITE(6,*) 
       IF(ICNFBAT.EQ.1) THEN
         CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF)
         WRITE(6,*)
         CALL WRTMAT(DETVEC,1,NCCM_STRING,1,NCCM_STRING)
       ELSE
         LBLK = -1
         CALL WRTVCD(DETVEC,LU_CSF,1,LBLK)
         WRITE(6,*)
         CALL WRTVCD(DETVEC,LU_DET,1,LBLK)
       END IF! ICNFBAT switch
      END IF! NTEST switch
*
      CALL LUCIAQEXIT('CSDTV')
      RETURN
      END
      SUBROUTINE GEN_NCONF_FOR_OCCLSN(IOCCLS,NGAS,MINOP,MAXOP,NIRREP,
     &           IB_ORB,NOBPT,NCONF_OCCLS_ALLSYM,NCONF_PER_OP_SM,
     &           NCONF_OCCLS)
*
*
* Generate 
* NCONF_PER_OP_SM: Number of configurations of occclass IOCCLS 
*                  for all symmetries and atleast MINOP open orbitals 
* NCONF_OCCLS_ALLSYM: Number of configurations of occclass IOCCLS 
*                     and any symmetry and any number of open orbitals
* NCONF_OCCLS: Number of configurations of occlass IOCCLS af given symmetry
*              and atleast MINOP open orbitals
*
* Jeppe Olsen, Febr. 2012, Geneva, speeded up version 
*       
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
*
*.. Input
*
*. Number of electrons per gas space 
      INTEGER IOCCLS(NGAS)  
*. Number of orbitals per gasspace 
      INTEGER NOBPT(NGAS)
*
*.. Output
*
*. Number of configurations per number of open shells, all symmetries
      INTEGER NCONF_PER_OP_SM(MAXOP+1,NIRREP)
      INTEGER NCONF_OCCLS(NIRREP)
*. Local scratch
      INTEGER JCONF(2*MXPORB)
      INTEGER MIN_OCC(MXPORB),MAX_OCC(MXPORB)
      INTEGER IOCC(2*MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GEN_NCONF_FOR_OCCLS'
        WRITE(6,*) ' ============================ '
        WRITE(6,*) 
        WRITE(6,*) ' Occupation class in action '    
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NGAS, MAXOP = ', NGAS, MAXOP
        WRITE(6,*) ' IB_ORB = ', IB_ORB
        WRITE(6,*) ' NOBPT: '
        CALL IWRTMA(NOBPT,1,NGAS,1,NGAS)
      END IF
*
*. Set up min and max arrays of occupations over orbitals
*
      JCONF(1) = IOCCLS(1)
      DO IGAS = 2, NGAS
        JCONF(IGAS) = JCONF(IGAS-1)+IOCCLS(IGAS)
      END DO
      CALL MXMNOC_GAS(MIN_OCC,MAX_OCC,NGAS,NOBPT,
     &     JCONF,JCONF,NTEST)
C     MXMNOC_GAS(MINEL_ORB,MAXEL_ORB,NGAS,NOBPT,
C    &                  MINEL_GAS,MAXEL_GAS,
C    &                  NTESTG)
*. Total number of electrons 
      NEL = IELSUM(IOCCLS,NGAS)
      IF(NTEST.GE.1000) WRITE(6,*) ' NEL = ', NEL
*. And total number of orbitals in GAS
      NORB = IELSUM(NOBPT,NGAS)
*
      IZERO = 0
      CALL ISETVC(NCONF_PER_OP_SM,IZERO,(MAXOP+1)*NIRREP)
      CALL ISETVC(NCONF_OCCLS,IZERO,NIRREP)
*. Loop over configurations 
      INI = 1
      NCONF = 0
      ISUM = 0
      NCONF_OCCLS_ALLSYM = 0
*. Loop over configurations
 1000 CONTINUE
        IF(NTEST.GE.10000)
     &  WRITE(6,*) ' NEXT_CONF_FOR_OCCLS will be called '
COLD    CALL NEXT_CONF_FOR_OCCLS
COLD &        (JCONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
C            NEXT_CONF_FROM_MINMAX_OCC(IACOCC,
C    &       IACOCC_MIN,IACOCC_MAX,INI,NO_NEW,NORB)
        CALL NEXT_CONF_FROM_MINMAX_OCC(IOCC,MIN_OCC,MAX_OCC,
     &       INI,NONEW,NORB)
        ISUM = ISUM + 1
        INI = 0
*. The configuration is returned in JCONF as an accumulated occupation
* over all orbitals. Change into list giving the occupied orbitals
        IEL = 0
        DO IORB = 1, NORB
          IF(IORB.EQ.1) THEN
            IIOCC = IOCC(1)
          ELSE
            IIOCC = IOCC(IORB)-IOCC(IORB-1)
          END IF
          IF(IIOCC.EQ.1) THEN
           IEL = IEL + 1
           JCONF(IEL) = IORB
          ELSE IF(IIOCC.EQ.2) THEN
           IEL = IEL + 1
           JCONF(IEL) = IORB
           IEL = IEL + 1
           JCONF(IEL) = IORB
          END IF
        END DO
*
        IF(NONEW.EQ.0) THEN
*. Check symmetry and number of open orbitals for this space
          IADD = IB_ORB - 1
C              IADCONST(IVEC,IADD, NDIM)
          CALL IADCONST(JCONF,IADD,NEL)
          ISYM_CONF = ISYMST(JCONF,NEL)
          IADD = - IADD
          CALL IADCONST(JCONF,IADD,NEL)
          NOPEN     = NOP_FOR_CONF(JCONF,NEL) 
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Number of open shells and SYM', 
     &      NOPEN, ISYM_CONF
          END IF
          NOCOB =  NOPEN + (NEL-NOPEN)/2
          NCONF_OCCLS_ALLSYM = NCONF_OCCLS_ALLSYM + 1 
          IF(NOPEN.GE.MINOP) THEN
*. A new configuration to be included, reform and save in packed form
            NCONF_PER_OP_SM(NOPEN+1,ISYM_CONF) =
     &      NCONF_PER_OP_SM(NOPEN+1,ISYM_CONF) + 1
            NCONF_OCCLS(ISYM_CONF) = NCONF_OCCLS(ISYM_CONF) + 1
          END IF
      GOTO 1000
        END IF !End if nonew = 0
* 
C     WRITE(6,*) ' TEST,  NCONF_OCCLS_ALLSYM = ',  NCONF_OCCLS_ALLSYM
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*)  ' =============================================== '
        WRITE(6,*)  ' Information on number of configuration for occls'
        WRITE(6,*)  ' =============================================== '
        WRITE(6,*)
        WRITE(6,*) ' Occupation class in action: '
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
        WRITE(6,*) ' Number of included configurations per symmetry'
        CALL IWRTMA(NCONF_OCCLS,1,NIRREP,1,NIRREP)
        WRITE(6,*) 
     &  ' Numbers:  Open orbitals +1 (ROW) and  sym (COLUMN)  '
          CALL IWRTMA(NCONF_PER_OP_SM,MAXOP+1,NIRREP,MAXOP+1,NIRREP)
      END IF

*
      RETURN
      END
      SUBROUTINE EXP_BLKVEC(LU_IN,NBLK_IN, IBLK_IN,
     &                      LU_OUT,NBLK_OUT,IBLK_OUT,
     &                      LBLK,ITASK,VEC,VEC_OUT,IREW,ICISTR,
     &                      INCORE)
*
* A blocking of some space with length LBLK(IBLK) for blocks IBLK 
* is given.
* 
* ITASK = 1:
* =========
* Input vector LU_IN contains the NBLK_IN blocks in IBLK_IN.
* A partition IBLK_OUT is also given. Both IBLK_IN and IBLK_OUT
* is assumed to be in increasing order
* Copy the blocks IBLK_OUT from LU_IN to LU_OUT.
* If a block in IBLK_OUT is not in IBLK_IN then 
* a zero block is written
*
* ITASK = 2:
* =========
* The file LU_IN contains the records given by IBLK_OUT
* copy to LU_OUT those record on LU_IN that are in IBLK_OUT..
*             but not in LU_IN (yes, is needed, sometimes...)
*
* If INCORE = 1, then vectors are in core at entry and are 
* not saved at exit
*
* Jeppe Olsen, Febr. 2012, Geneva
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IBLK_IN(NBLK_IN),IBLK_OUT(NBLK_OUT)
      INTEGER LBLK(*)
*. Scratch: Each holding one block(ICISTR>1) or one vector(ICISTR=1)
      REAL*8 VEC(*),VEC_OUT(*)
*
      NTEST = 1000
      IF(NTEST.GE.10) THEN
        WRITE(6,*)  ' Output from EXP_BLKVEC '
        WRITE(6,*)  ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' ITASK, ICISTR = ', ITASK, ICISTR
        WRITE(6,*) ' LU_IN, LU_OUT = ', LU_IN, LU_OUT
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' IBLK_IN, IBLK_OUT = '
        CALL IWRTMA(IBLK_IN,1,NBLK_IN,1,NBLK_IN)
        WRITE(6,*)
        CALL IWRTMA(IBLK_OUT,1,NBLK_OUT,1,NBLK_OUT)
      END IF
*
      IF(IREW.EQ.1.AND.INCORE.EQ.0) THEN
        CALL REWINO(LU_IN)
        CALL REWINO(LU_OUT)
      END IF
*
      IF(INCORE.EQ.1.AND.ICISTR.NE.1) THEN
        WRITE(6,*) ' Conflicting parameters: INCORE, ICISTR = ',
     &                                       INCORE, ICISTR
        STOP       ' Conflicting parameters: INCORE, ICISTR  '
      END IF
*
      ZERO = 0.0D0
      IIIBLK_IN = 1
*
      IF(ICISTR.EQ.1) THEN
*
* One vector in core mode
*
        LENGTH_IN = 0
        DO IBLK = 1, NBLK_IN
          LENGTH_IN = LENGTH_IN + LBLK(IBLK_IN(IBLK))
        END DO
*
        LENGTH_OUT = 0
        DO IBLK = 1, NBLK_OUT
          LENGTH_OUT = LENGTH_OUT + LBLK(IBLK_OUT(IBLK))
        END DO
*
        IF(INCORE.EQ.0) 
     &  CALL FRMDSC(VEC,LENGTH_IN,-1,LU_IN,IMZERO,IAMPACK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Input vector '
          CALL WRTMAT(VEC,1,LENGTH_IN,1,LENGTH_IN)
        END IF
*
        IF(ITASK.EQ.1) THEN
          IOFF_IN = 1
          IOFF_OUT = 1
*. Loop over occ classes
          IIBLK_IN = 1
          DO IIBLK_OUT = 1, NBLK_OUT
*
  999       CONTINUE
            IF(IBLK_IN(IIBLK_IN).LT.IBLK_OUT(IIBLK_OUT).AND.
     &         IIBLK_IN.LT.NBLK_IN) THEN
*. Input block is lower than outputblock, output block may arise later
              IOFF_IN = IOFF_IN +  LBLK(IBLK_IN(IIBLK_IN))
              IIBLK_IN = IIBLK_IN + 1
              GOTO 999
            ELSE IF(IBLK_IN(IIBLK_IN).EQ.IBLK_OUT(IIBLK_OUT)) THEN
*. Match between input and output blocks
              LENGTH = LBLK(IBLK_OUT(IIBLK_OUT))
              CALL COPVEC(VEC(IOFF_IN),VEC_OUT(IOFF_OUT),LENGTH)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' Match: IIBLK_IN, IIBLK_OUT = ',
     &                              IIBLK_IN, IIBLK_OUT
                WRITE(6,*) ' IOFF_IN, IOFF_OUT, LENGTH = ',
     &                       IOFF_IN, IOFF_OUT, LENGTH
                WRITE(6,*) ' VEC(IOFF_IN),VEC_OUT(IOFF_OUT) = ',
     &                       VEC(IOFF_IN),VEC_OUT(IOFF_OUT)
                WRITE(6,*) ' Updated block '
                CALL WRTMAT(VEC_OUT(IOFF_OUT),1,LENGTH,1,LENGTH)
              END IF
              IOFF_IN = IOFF_IN + LENGTH
              IOFF_OUT = IOFF_OUT + LENGTH
              IIBLK_IN = IIBLK_IN + 1
            ELSE
*. Outputblock is not in input, so set to zero
              LENGTH = LBLK(IBLK_OUT(IIBLK_OUT))
              CALL SETVEC(VEC_OUT(IOFF_OUT),ZERO,LENGTH)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' Block zeroed '
                WRITE(6,*) ' IIBLK_OUT, IOFF_OUT, LENGTH = ',
     &                       IIBLK_OUT, IOFF_OUT, LENGTH
              END IF
              IOFF_OUT = IOFF_OUT + LENGTH
            END IF
          END DO! loop over output blocks
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Input vector(again) '
            CALL WRTMAT(VEC,1,LENGTH_IN,1,LENGTH_IN)
            WRITE(6,*) ' Output vector '
            CALL WRTMAT(VEC_OUT,1,LENGTH_OUT,1,LENGTH_OUT)
          END IF
*. And save on Disc
          IF(INCORE.EQ.0) 
     &    CALL TODSCN(VEC_OUT,1,LENGTH_OUT,-1,LU_OUT)
*
        ELSE IF (ITASK.EQ.2) THEN
          WRITE(6,*) ' ITASK = 2, ICISTR = 1 not programmed '
          STOP       ' ITASK = 2, ICISTR = 1 not programmed '
        END IF ! ITASK switch
*
      ELSE IF (ICISTR.GE.2) THEN
*
* One block in core mode
*
        IF(ITASK.EQ.1) THEN
*
*
          DO IIIBLK_OUT = 1, NBLK_OUT
            IIBLK_OUT = IBLK_OUT(IIIBLK_OUT)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' IIIBLK_OUT, IIBLK_OUT = ',
     &                     IIIBLK_OUT, IIBLK_OUT
            END IF
 1001       CONTINUE
            IIBLK_IN = IBLK_IN(IIIBLK_IN)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' IIBLK_IN = ', IIBLK_IN
            END IF
            LL = LBLK(IIBLK_OUT)
*. Compare with next input block
            IF(IIBLK_IN.GT.IIBLK_OUT) THEN
*. The next input block is higher than the current output, so 
*. current output is vanishing
              IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' BLock is not in output '
              END IF
              CALL SETVEC(VEC,ZERO,LL)
              CALL TODSCN(VEC,1,LL,-1,LU_OUT)
            ELSE IF (IIBLK_IN.EQ.IIBLK_OUT) THEN
*. We gave a match
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' We have a match '
              END IF
              CALL FRMDSCN(VEC,1,-1,LU_IN)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' Block read in: '
                CALL WRTMAT(VEC,1,LL,1,LL)
              END IF
              CALL TODSCN(VEC,1,LL,-1,LU_OUT)   
              IIIBLK_IN = IIIBLK_IN + 1
            ELSE
*. Current output block is higher than current input block. step
*. input block and take another look
              READ(LU_IN)
              READ(LU_IN)
              IIIBLK_IN = IIIBLK_IN + 1
              GOTO 1001
            END IF! switch IIBLK_IN .GT. IIBLK_OUT
          END DO! loop over IIIBLK_OUT
*
*
        ELSE IF (ITASK.EQ.2) THEN
*
*
         IIIBLK_IN = 1
         DO IIIBLK_OUT = 1, NBLK_OUT
           IIBLK_OUT = IBLK_OUT(IIIBLK_OUT)
           LL = LBLK(IIBLK_OUT)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' IIIBLK_OUT, IIBLK_OUT, LL = ', 
     &                    IIIBLK_OUT, IIBLK_OUT, LL
           END IF
 2001      CONTINUE
           IIBLK_IN = IBLK_IN(IIIBLK_IN)
           IF(NTEST.GE.100) 
     &     WRITE(6,*) ' IIIBLK_IN, IIBLK_IN = ', IIIBLK_IN, IIBLK_IN
           IF(IIBLK_IN.GT.IIBLK_OUT) THEN
*. Not excluded
             IF(NTEST.GE.100)
     &       WRITE(6,*) ' Copied from in to out '
             CALL FRMDSCN(VEC,1,-1,LU_IN)
C                      (VEC,NREC,LBLK,LU)
C (VEC,NREC,LBLK,LU)
             CALL TODSCN(VEC,1,LL,-1,LU_OUT)
C                 TODSCN(VEC,NREC,LREC,LBLK,LU)
           ELSE IF(IIBLK_IN.EQ.IIBLK_OUT) THEN
*. Excluded
             IF(NTEST.GE.100)
     &       WRITE(6,*) ' Excluded '
             CALL SETVEC(VEC,ZERO,LL)
             CALL TODSCN(VEC,1,LL,-1,LU_OUT)
             READ(LU_IN)
             READ(LU_IN)
           ELSE
*. Input block was lower than output block, so no verdict yet
             IIIBLK_IN = IIIBLK_IN + 1
             GOTO 2001
           END IF! switch IIBLK_IN .GT. IIBLK_OUT
         END DO
        END IF! ITASK switch
*. Remember to put EOV on LU_OUT
        CALL ITODS(-1,1,-1,LU_OUT)
      END IF! ICISTR switch
*
      RETURN
      END
      SUBROUTINE GET_NCSF_PER_OCCLS_FOR_CISPACE(ISYM,IOCCLS_ACT,
     &           NOCCLS_ACT,NCS_FOR_OCCLS,NCS_FOR_OCCLS_ACT)
*
* A CI space is defined by the NOCCLS_ACT occ classes in IOCCLS_ACT
* Collect number of CSF's per occlass from the general NCS_FOR_OCCLS
* array
*
*. Jeppe Olsen, Febr. 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOCCLS_ACT(NOCCLS_ACT)
      INTEGER NCS_FOR_OCCLS(NIRREP,*)
*. Output
      INTEGER NCS_FOR_OCCLS_ACT(NOCCLS_ACT)
*
      NTEST = 0
*
      DO IIOCCLS = 1, NOCCLS_ACT
        IOCCLS = IOCCLS_ACT(IIOCCLS)
        NCS_FOR_OCCLS_ACT(IIOCCLS) = NCS_FOR_OCCLS(ISYM,IOCCLS)
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Symmetry: ', ISYM
       WRITE(6,*) ' Number of CSFs per occlass in CI expansion '
       CALL IWRTMA(NCS_FOR_OCCLS_ACT,1,NOCCLS_ACT,1,NOCCLS_ACT)
      END IF
*
      RETURN
      END
      FUNCTION ISYM_CONF_G(IOCC,NOCCL,NORBL,IORB_OFF,ICONF_FORM)
*
* Obtain symmetry of configuration IOCC 
*
* If ICONF_FORM = 1, the occupation is specified as 
* occupations of each orbital - i.e as occupation number vector
* (and NOCCL may be a dummy integer)
*
* If ICONF_FORM = 2, the occupation is specified in
* packed form - as NOCCL occupied orbitals with minus indicating
* a doubly occupied orbital
*
* Jeppe Olsen, June 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. General input
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
* Specific input
      DIMENSION IOCC(*)
*
      IF(ICONF_FORM.EQ.1) THEN
        NDIM = NORBL
      ELSE
        NDIM = NOCCL
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from ISYM_CONF:'
        WRITE(6,*) ' ==================='
        IF(ICONF_FORM.EQ.1) THEN
          WRITE(6,*) ' Configuration is given in ONV form '
        ELSE IF(ICONF_FORM.EQ.2) THEN
          WRITE(6,*) ' Configuration is given in PACKED form '
        END IF
*
        WRITE(6,*) ' NORBL, IORB_OFF =', NORBL,IORB_OFF
        WRITE(6,*) '  IOCC: '
        CALL IWRTMA3(IOCC,1,NDIM,1,NDIM)
      END IF
*
      IF(1.GT.ICONF_FORM.OR.ICONF_FORM.GT.2) THEN
        WRITE(6,*) ' Illegal value of ICONF_FORM in ISYM_CONF_G, ',
     &               ICONF_FORM
        STOP       ' Illegal value of ICONF_FORM in ISYM_CONF_G, '
      END IF
*
      ISYM = 1
      IF(ICONF_FORM.EQ.1) THEN
*. Occupation number representation
        DO IORB = 1, NORBL
          IF(IOCC(IORB).EQ.1)
     &    ISYM = MULTD2H(ISYM,ISMFTO(IORB-1+IORB_OFF))
        END DO
      ELSE IF (ICONF_FORM.EQ.2) THEN
*. Packed representation
       DO IORB = 1, NOCCL
         IF(IOCC(IORB).GT.0) THEN
           IIORB = IOCC(IORB)
C?         WRITE(6,*) ' IORB, IIORB = ', IORB, IIORB
           ISYM =  MULTD2H(ISYM,ISMFTO(IIORB-1+IORB_OFF))
         END IF
       END DO
      END IF! switch between the various representations of occupations
*
      ISYM_CONF_G = ISYM
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from ISYM_CONF_G '
       WRITE(6,*) ' ======================'
       WRITE(6,*) ' Configuration: '
       CALL IWRTMA(IOCC,1,NDIM,1,NDIM)
       WRITE(6,*) ' Symmetry = ', ISYM
      END IF
*
      RETURN
      END
      SUBROUTINE DIHDJ2_LUCIA_CONF
     &(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,IADOB,NORB,
     & IHORS,HAMIL,C,SIGMA,IWORK,ISYM,ECORE,ICOMBI,PSIGN,
     & NTERMS,NDIF0,NDIF1,NDIF2,I12OP,I_DO_ORBTRA,IORBTRA,
     & NTOOB,RJ,RK)
C
C A SET OF DETERMINANTS IA DEFINED BY ALPHA AND BETASTRINGS
C IASTR,IBSTR AND ANOTHER SET OF DETERMINATS DEFINED BY STRINGS
C JASTR AND JBSTR ARE GIVEN. 
*
* IHORS = 1: OBTAIN CORRESPONDING HAMILTONIAN MATRIX: <IDET!H!JDET>
* IHORS = 2:  Multiply elements on C and store in sigma
*            Sigma(I) = Sum_J <IDET!H!JDET>C_JDET
C
C IF ICOMBI .NE. 0 COMBINATIONS ARE USED FOR ALPHA AND BETA STRING
C THAT DIFFERS :
C   1/SQRT(2) * ( !I1A I2B! + PSIGN * !I2A I1B! )
C
* For IHORS = 1
C IF ISYM .EQ. 0 FULL HAMILTONIAN IS CONSTRUCTED
C IF ISYM .NE. 0 LOWER HALF OF HAMILTONIAN IS CONSTRUCTED
C
C JEPPE OLSEN JANUARY 1989, Slight changes July 2011
*. Last Modification; May 19, 2013; Jeppe Olsen, corrected an 
*                     error for combinations
*
* In the current version, it is assumed and required that the 
* left and right vectors belong to given configurations.
* This may be used to simplify the integral handling
*
* The orbitals in the strings may have a different offset than 
* the global addressing of the orbitals. A constant IADOB may
* therefore be added.
*
* IF I12OP = 1, then only one-electron operator is included
* IF I_DO_ORBTRA = 1, then orbital transformation is in progress, 
* with  orbital IORBTRA being transformed
*
*. Note: use of RJ, RK commented temporary out, June 4, 2013
C
      INCLUDE 'implicit.inc'
*. Input
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION JASTR(NAEL,*),JBSTR(NBEL,*)
      DIMENSION C(*)
*. Exchange and coulomb integrals
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*. Scratch: Length: 4*NORB  ( + NIDET if ICOMBI ne 0)
      INTEGER IWORK(*)
*. Output
      DIMENSION HAMIL(*), SIGMA(*)
*
      NTEST = 000
      IPRT = NTEST
*. To eliminate compiler warnings
      KLIAB = 0
      IAEQIB = 0
      JAEQJB = 0
      IEL1 = 0
      JEL1 = 0
      SIGNA = 0.0D0
      SIGNB = 0.0D0
      IPERM = 0   
      JPERM = 0
      SIGN = 0.0D0
      IUSE_JK = 1
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Welcome to DIHDJ2_LUCIA '
       WRITE(6,*) ' ======================== '
       WRITE(6,*) ' Number of Determinants in L and R ',
     &              NIDET, NJDET
       WRITE(6,*) ' ECORE, ISYM = ', ECORE, ISYM
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input IA and IB strings'
        CALL IWRTMA(IASTR,NAEL,NIDET,NAEL,NIDET)
        CALL IWRTMA(IBSTR,NBEL,NIDET,NBEL,NIDET)
        WRITE(6,*) ' Input JA and JB strings'
        CALL IWRTMA(JASTR,NAEL,NJDET,NAEL,NJDET)
        CALL IWRTMA(JBSTR,NBEL,NJDET,NBEL,NJDET)
        WRITE(6,*) ' IADOB = ', IADOB
      END IF
*
      IF(IHORS.EQ.2.AND.NTEST.GE.1000) THEN
        WRITE(6,*) ' C- vector '
        CALL WRTMAT(C,1,NJDET,1,NJDET)
        WRITE(6,*) ' Initial Sigma vector '
        CALL WRTMAT(SIGMA,1,NIDET,1,NIDET)
      END IF
*
C?    IF(NTEST.GE.10000) THEN
C?      WRITE(6,*) ' RJ and RK matrices '
C?      CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
C?      WRITE(6,*)
C?      CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
C?    END IF
*
*. Scratch space 
*
* .. 1 : EXPANSION OF ALPHA AND BETA STRINGS OF TYPE I
      KLIAE = 1
      KLIBE = KLIAE + NORB
      KLJAE = KLIBE + NORB
      KLJBE = KLJAE + NORB
      KLFREE = KLJBE + NORB
      IF( ICOMBI .NE. 0 ) THEN
        KLIAB  = KLFREE
        KLFREE = KLFREE + NIDET
      END IF
*
      IF( ICOMBI .NE. 0 ) THEN
* SET UP ARRAY COMBARING ALPHA AND BETA STRINGS IN IDET LIST
        DO IDET = 1, NIDET
          IAEQIB = 1
          DO IEL = 1, NAEL
            IF(IASTR(IEL,IDET) .NE. IBSTR(IEL,IDET))IAEQIB = 0
          END DO  
          IWORK(KLIAB-1+IDET) = IAEQIB
        END DO  
      END IF ! ICOMBI.ne.0
*
      IF(IHORS.EQ.1) THEN
        IF( ISYM .EQ. 0 ) THEN
          LHAMIL = NIDET*NJDET
        ELSE
          LHAMIL = NIDET*(NIDET+1) / 2
        END IF
        CALL SETVEC(HAMIL,0.0D0,LHAMIL)
      END IF
*
      NTERMS= 0
      NDIF0 = 0
      NDIF1 = 0
      NDIF2 = 0
C.. LOOP OVER J DETERMINANTS
C
      DO JDET = 1,NJDET
        IF( NTEST .GE. 2000 ) WRITE(6,*) '  ****** JDET ', JDET
*
* EXPAND JDET
*
        IF( ICOMBI .NE. 0 ) THEN
          JAEQJB = 1
          DO IEL = 1, NAEL
            IF(JASTR(IEL,JDET) .NE. JBSTR(IEL,JDET))JAEQJB = 0
          END DO   
C?        WRITE(6,*) ' JAEQJB ', JAEQJB
        END IF
*
        CALL ISETVC(IWORK(KLJAE),0,NORB)
        CALL ISETVC(IWORK(KLJBE),0,NORB)
*
        DO IAEL = 1, NAEL
          IWORK(KLJAE-1+JASTR(IAEL,JDET) ) = 1
        END DO   
*
        DO IBEL = 1, NBEL
          IWORK(KLJBE-1+JBSTR(IBEL,JDET) ) = 1
        END DO   
*
        IF( NTEST .GE. 2000 ) THEN
          WRITE(6,*) ' JDET =  ',JDET
          WRITE(6,*) ' JASTR AND JBSTR '
          CALL IWRTMA(JASTR(1,JDET),1,NAEL,1,NAEL)
          CALL IWRTMA(JBSTR(1,JDET),1,NBEL,1,NBEL)
          WRITE(6,*) ' EXPANDED ALPHA AND BETA STRING '
          CALL IWRTMA(IWORK(KLJAE),1,NORB,1,NORB)
          CALL IWRTMA(IWORK(KLJBE),1,NORB,1,NORB)
        END IF
C
        IF( ISYM .EQ. 0 ) THEN
          MINI = 1
        ELSE
          MINI = JDET
        END IF
*
        DO IDET = MINI, NIDET
C?      WRITE(6,*) '   IDET .... ',IDET
*
        IF( ICOMBI .EQ. 0 ) THEN
            NLOOP = 1
        ELSE
          IAEQIB = IWORK(KLIAB-1+IDET)
          IF(IAEQIB+JAEQJB .EQ. 0 ) THEN
            NLOOP = 2
          ELSE
            NLOOP = 1
          END IF
        END IF
*
        DO 899 ILOOP = 1, NLOOP
         NTERMS = NTERMS + 1
C?       WRITE(6,*) '   899 : ILOOP ' , ILOOP
*
*.. COMPARE DETERMINANTS
*
* SWAP IA AND IB FOR SECOND PART OF COMBINATIONS
        IF( ILOOP .EQ. 2 )
     &  CALL ISWPVE(IASTR(1,IDET),IBSTR(1,IDET),NAEL)
*. Number of common alpha occupations
        NACM = 0
        DO  IAEL = 1, NAEL
          NACM = NACM + IWORK(KLJAE-1+IASTR(IAEL,IDET))
        END DO   
*. Number of common beta occupations
        NBCM = 0
        DO IBEL = 1, NBEL
          NBCM = NBCM + IWORK(KLJBE-1+IBSTR(IBEL,IDET))
        END DO   
*. And number of differences..
        NADIF = NAEL-NACM
        NBDIF = NBEL-NBCM
*
        INTERACT = 1
        IF(NADIF+NBDIF.GT.I12OP) INTERACT = 0
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IDET, JDET, ILOOP, INTERACT = ',
     &                 IDET, JDET, ILOOP, INTERACT
          WRITE(6,*) ' COMPARISON , NADIF , NBDIF ', NADIF,NBDIF
        END IF
        IF(INTERACT.EQ.0) GOTO 898
*  FACTOR FOR COMBINATIONS
          CONST = 1.0D0
          IF( ICOMBI .EQ. 0 ) THEN
            CONST = 1.0D0
          ELSE
           IF((JAEQJB +IAEQIB) .EQ.2 ) THEN
             CONST = 1.0D0
           ELSE IF( (JAEQJB+IAEQIB) .EQ. 1 ) THEN
             CONST = 1.0D0/SQRT(2.0D0)*(1.0D0+PSIGN)
            ELSE IF( (JAEQJB+IAEQIB) .EQ. 0 ) THEN
             IF( ILOOP .EQ. 1)  THEN
               CONST = 1.0D0
             ELSE
               CONST = PSIGN
             END IF
            END IF
          END IF
*
* ======================================================================
*.. Obtain the orbitals differring in the two dets and sign to obtain 
*   max juxtaposition
* ======================================================================
*
*
* EXPAND IDET
          CALL ISETVC(IWORK(KLIAE),0,NORB)
          CALL ISETVC(IWORK(KLIBE),0,NORB)
*
          DO  IAEL = 1, NAEL
            IWORK(KLIAE-1+IASTR(IAEL,IDET) ) = 1
          END DO   
*
          DO IBEL = 1, NBEL
            IWORK(KLIBE-1+IBSTR(IBEL,IDET) ) = 1
          END DO   
*
          IF(NADIF .EQ. 1 ) THEN
            DO IAEL = 1,NAEL
              IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
                IA = IASTR(IAEL,IDET)
                IEL1 = IAEL
                GOTO 121
               END IF
            END DO
  121       CONTINUE
C
            DO JAEL = 1,NAEL
              IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
                JA = JASTR(JAEL,JDET)
                JEL1 = JAEL
                GOTO 131
              END IF
            END DO    
  131       CONTINUE
            SIGNA = (-1)**(JEL1+IEL1)
C?          WRITE(6,*) ' IA JA SIGNA... ',IA,JA,SIGNA
          END IF ! NADIF = 1
 
          IF(NBDIF .EQ. 1 ) THEN
            DO IBEL = 1,NBEL
              IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
                IB = IBSTR(IBEL,IDET)
                IEL1 = IBEL
                GOTO 221
               END IF
            END DO
  221       CONTINUE
*
            DO JBEL = 1,NBEL
              IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
                JB = JBSTR(JBEL,JDET)
                JEL1 = JBEL
                GOTO 231
               END IF
            END DO   
  231       CONTINUE
            SIGNB = (-1)**(JEL1+IEL1)
C?          WRITE(6,*) ' 1B: IB JB SIGNB... ',IB,JB,SIGNB
          END IF ! NBDIF = 1
*
          IF(NADIF .EQ. 2 ) THEN
            IDIFF = 0
            DO IAEL = 1,NAEL
              IF(IWORK(KLJAE-1+IASTR(IAEL,IDET)).EQ.0) THEN
                IF( IDIFF .EQ. 0 ) THEN
                  IDIFF = 1
                  I1 = IASTR(IAEL,IDET)
                  IPERM = IAEL
                ELSE
                I2 = IASTR(IAEL,IDET)
                  IPERM = IAEL + IPERM
                  GOTO 321
                END IF
              END IF
            END DO   
  321       CONTINUE
*
            JDIFF = 0
            DO JAEL = 1,NAEL
              IF(IWORK(KLIAE-1+JASTR(JAEL,JDET)).EQ.0) THEN
                IF( JDIFF .EQ. 0 ) THEN
                  JDIFF = 1
                  J1 = JASTR(JAEL,JDET)
                  JPERM = JAEL
                ELSE
                  J2 = JASTR(JAEL,JDET)
                  JPERM = JAEL + JPERM
                  GOTO 331
                END IF
              END IF
            END DO  
  331       CONTINUE
            SIGN = (-1)**(IPERM+JPERM)
          END IF ! NADIF = 2
C
          IF(NBDIF .EQ. 2 ) THEN
            IDIFF = 0
            DO IBEL = 1,NBEL
              IF(IWORK(KLJBE-1+IBSTR(IBEL,IDET)).EQ.0) THEN
                IF( IDIFF .EQ. 0 ) THEN
                  IDIFF = 1
                  I1 = IBSTR(IBEL,IDET)
                  IPERM = IBEL
                ELSE
                  I2 = IBSTR(IBEL,IDET)
                  IPERM = IBEL + IPERM
                  GOTO 421
                 END IF
              END IF
            END DO   
  421       CONTINUE
*
            JDIFF = 0
            DO JBEL = 1,NBEL
              IF(IWORK(KLIBE-1+JBSTR(JBEL,JDET)).EQ.0) THEN
                IF( JDIFF .EQ. 0 ) THEN
                  JDIFF = 1
                  J1 = JBSTR(JBEL,JDET)
                  JPERM = JBEL
                ELSE
                  J2 = JBSTR(JBEL,JDET)
                  JPERM = JBEL + JPERM
                  GOTO 431
                END IF
              END IF
            END DO
  431       CONTINUE
            SIGN = (-1)**(IPERM+JPERM)
          END IF ! NBDIF = 2
*
* =====================================
*  Obtain Value of Hamiltonian element
* =====================================
*
          XVAL = 0.0D0
          IF( NADIF .EQ. 2 .OR. NBDIF .EQ. 2 ) THEN
*
* Difference in occupation of two alpha or of two beta orbitals
*
* Value:    SIGN * (I1 J1 ! I2 J2 ) - ( I1 J2 ! I2 J1 )
            NDIF2 = NDIF2 + 1
            I1 = I1 + IADOB
            I2 = I2 + IADOB
            J1 = J1 + IADOB
            J2 = J2 + IADOB
            XVAL = 
     &      SIGN*( GTIJKL_GN(I1,J1,I2,J2)-GTIJKL_GN(I1,J2,I2,J1) )
          ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 1 ) THEN
*
* Difference in occupation of one alpha and one beta orbitals
*
* Value:    SIGN * (IA JA ! IB JB )
            NDIF2 = NDIF2 + 1
            IA = IA + IADOB
            IB = IB + IADOB
            JA = JA + IADOB
            JB = JB + IADOB
            IF(NTEST.GE.10000) WRITE(6,*) ' IA IB JA JB ', IA,IB,JA,JB
            IF(IUSE_JK.EQ.1.AND.IA.EQ.JB.AND.IB.EQ.JA) THEN
              XVAL = SIGNA*SIGNB*RK(IA,IB)
            ELSE
              XVAL = SIGNA*SIGNB* GTIJKL_GN(IA,JA,IB,JB)
            END IF
            IF(NTEST.GE.10000) 
     &      WRITE(6,*) ' SIGNA SIGNB XVAL ',SIGNA,SIGNB,XVAL
          ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 0 .OR.
     &             NADIF .EQ. 0 .AND. NBDIF .EQ. 1 )THEN
*
* Difference in occupation of one alpha or beta orbital
*
            NDIF1 = NDIF1 + 1
C SIGN *
C(  H(I1 J1 ) +
C  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 ! JORB JORB )
C -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB ! JORB J1 ) )
            IF( NADIF .EQ. 1 ) THEN
              I1 = IA + IADOB
              J1 = JA + IADOB
              SIGN = SIGNA
            ELSE
              I1 = IB + IADOB
              J1 = JB + IADOB
              SIGN = SIGNB
            END IF
C?          WRITE(6,*) ' ONE DIFF I1 J1 SIGN : ',I1,J1,SIGN
C
            XVAL = GETH1I(I1,J1)
            IF(I12OP.EQ.2) THEN
              DO JAEL = 1, NAEL
                JORB = JASTR(JAEL,JDET)+IADOB
C?              WRITE(6,*) ' A I1 = ', I1
                XVAL = XVAL + GTIJKL_GN(I1,J1,JORB,JORB)
              END DO
*
              DO JBEL = 1, NBEL
                JORB = JBSTR(JBEL,JDET)+IADOB
C?              WRITE(6,*) ' B I1 = ', I1
                XVAL = XVAL + GTIJKL_GN(I1,J1,JORB,JORB)
              END DO
*
              IF( NADIF .EQ. 1 ) THEN
                DO JAEL = 1, NAEL
                  JORB = JASTR(JAEL,JDET)+IADOB
                  XVAL = XVAL - GTIJKL_GN(I1,JORB,JORB,J1)
                END DO
              ELSE
                DO JBEL = 1, NBEL
                  JORB = JBSTR(JBEL,JDET)+IADOB
C?                WRITE(6,*) ' C I1 = ', I1
                  XVAL = XVAL - GTIJKL_GN(I1,JORB,JORB,J1)
                END DO
              END IF ! NADIF = 1
            END IF ! I12OP = 2
            XVAL = XVAL * SIGN
          ELSE IF( NADIF .EQ. 0 .AND. NBDIF .EQ. 0 ) THEN
*
* Diagonal elements
*
C SUM(I,J OF JDET) H(I,J) + (I I ! J J ) - (I J ! J I )
            NDIF0 = NDIF0 + 1
            XVAL = ECORE
            DO IAB = 1,2
              IF(IAB .EQ. 1 ) THEN
                NIABEL = NAEL
              ELSE
                NIABEL = NBEL
              END IF
              DO JAB = 1, 2
                IF(JAB .EQ. 1 ) THEN
                  NJABEL = NAEL
                ELSE
                  NJABEL = NBEL
                END IF
                DO IEL = 1, NIABEL
                  IF( IAB .EQ. 1 ) THEN
                    IORB = IASTR(IEL,IDET) + IADOB
                  ELSE
                    IORB = IBSTR(IEL,IDET) + IADOB
                  END IF
*
                  IF(IAB .EQ. JAB ) XVAL = XVAL + GETH1I(IORB,IORB)
                  IF(I12OP.EQ.2) THEN
C?                  write(6,*)  ' XVAL after one body term ', XVAL
                    DO JEL = 1, NJABEL
                      IF( JAB .EQ. 1 ) THEN
                        JORB = IASTR(JEL,IDET)+IADOB
                      ELSE
                        JORB = IBSTR(JEL,IDET)+IADOB
                      END IF
*
                      IF(NTEST.GE.1000) THEN
                        WRITE(6,*) ' IAB, JAB, IORB, JORB = ',
     &                               IAB, JAB, IORB, JORB 
                      END IF
*
                      IF(IUSE_JK.EQ.0) THEN
                        XVAL = XVAL+0.5D0*GTIJKL_GN(IORB,IORB,JORB,JORB)
                        IF( IAB . EQ. JAB )
     &                  XVAL = XVAL-0.5D0*GTIJKL_GN(IORB,JORB,JORB,IORB)
                      ELSE
                        XVAL = XVAL + 0.5D0*RJ(IORB,JORB)
                        IF( IAB . EQ. JAB )
     &                  XVAL = XVAL - 0.5D0*RK(IORB,JORB)
                      END IF
C?                    write(6,*) ' XVAL with coulomb', XVAL
COLD                  IF( IAB . EQ. JAB )
COLD &                XVAL = XVAL - 0.5D0*GTIJKL_GN(IORB,JORB,JORB,IORB)
COLD &                XVAL = XVAL - 0.5D0*RK(IORB,JORB)
                      IF(NTEST.GE.1000) 
     &                WRITE(6,*) ' Updated XVAL = ', XVAL
                    END DO   
                  END IF ! I12OP = 2
                END DO   
              END DO   
            END DO   
          END IF ! End of switch between different differences
*. And save
          IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' IDET, JDET, ILOOP = ', 
     &                    IDET, JDET, ILOOP
             WRITE(6,*) ' CONST XVAL  ', CONST ,XVAL
          END IF
*
          IF(IHORS.EQ.1) THEN
            IF( ISYM .EQ. 0 ) THEN
              HAMIL((JDET-1)*NIDET+IDET) =
     &        HAMIL((JDET-1)*NIDET+IDET) + CONST * XVAL
            ELSE
              HAMIL((IDET-1)*IDET/2 + JDET ) =
     &        HAMIL((IDET-1)*IDET/2 + JDET ) + CONST * XVAL
            END IF
          ELSE
            SIGMA(IDET) = SIGMA(IDET) + CONST*XVAL*C(JDET)
            IF(ISYM.EQ.1) 
     &      SIGMA(JDET) = SIGMA(JDET) + CONST*XVAL*C(IDET)
            IF(NTEST.GE.100) THEN
              WRITE(6,'(A,2I4,2E13.6)') 
     &      ' IDET, JDET, C(JDET), SIGMA(IDET) = ',
     &        IDET, JDET, C(JDET), SIGMA(IDET)
              WRITE(6,'(A,2E13.6)') ' CONST, XVAL = ',
     &        CONST,XVAL
            END IF
          END IF
  898     CONTINUE
C RESTORE ORDER !!!
          IF( ILOOP .EQ. 2 )
     &    CALL ISWPVE(IASTR(1,IDET),IBSTR(1,IDET),NAEL)
  899 CONTINUE ! Anchor for nonvanishing interaction
        END DO ! Loop over IDET
      END DO ! loop over JDET
*
 
      IF( IPRT .GT. 100 ) THEN
        WRITE(6,*)
     &' Number of elements differing by 0 excitation.. ',NDIF0
        WRITE(6,*)
     &' Number of elements differing by 1 excitation.. ',NDIF1
        WRITE(6,*)
     &' Number of elements differing by 2 excitation.. ',NDIF2
        WRITE(6,*)
     &' Number of vanishing elments                    ',
     &  NTERMS - NDIF0 - NDIF1 - NDIF2
      END IF
*
      IF( IPRT .GE. 100 ) THEN
        IF(IHORS.EQ.1) THEN
          WRITE(6,*) '  HAMILTONIAN MATRIX '
          IF( ISYM .EQ. 0 ) THEN
            CALL WRTMAT(HAMIL,NIDET,NJDET,NIDET,NJDET)
          ELSE
            CALL PRSYM(HAMIL,NIDET)
          END IF
        ELSE
          WRITE(6,*) ' Updated Sigma vector '
          CALL WRTMAT(SIGMA,1,NIDET,1,NIDET)
        END IF! IHORS = 1,2
      END IF
*
      RETURN
      END
      SUBROUTINE CONF_EXP_LEN_LIST(ILEN,NCONF_PER_OPEN,
     &           NELMNT_PER_OPEN,MAXOP)
*
* Obtain ILEN: Length of expansion of each configuration. 
*              Number of elements (CSF's or SD'S) per conftype is given in
*              NELMNT_PER_OPEN
* 
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NCONF_PER_OPEN(MAXOP+1), NELMNT_PER_OPEN(MAXOP+1)
*. Output
      INTEGER ILEN(*)
*
      IB_CONF = 1
      DO IOPEN = 0, MAXOP
        N = NCONF_PER_OPEN(IOPEN+1)
        L = NELMNT_PER_OPEN(IOPEN+1)
        CALL ISETVC(ILEN(IB_CONF),L,N)
        IB_CONF = IB_CONF + N
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        NCONF = IB_CONF -1
        WRITE(6,*) 'Length of expansion of each configuration '
        WRITE(6,*) ' =========================================='
        CALL IWRTMA(ILEN,1,NCONF,1,NCONF)
      END IF
*
      RETURN
      END
      SUBROUTINE NEXT_CONF_IN_CONFSPC(IOCC,IOPEN,INUM_OP,INI,
     &           ISYM,ISPC,NEW)
*
* Obtain occupation of next configuration of sym ISPC in configuration list ISPC
* (with current being defined by IOPEN,INUM_OP)
* INI => Initialize
*
* IOPEN is number of open orbitals in output configuration
* INUMOP is address of configuration wrt start of confs with IOPEN open orbs.
*
* Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      INTEGER IOCC(*)
*
      NTEST = 000
*
      IOPEN_INI = IOPEN
      INUM_OP_INI = INUM_OP
*
      IF(INI.EQ.1) THEN
        IOPEN = -1
        INUM_OP = 1
      END IF
*
      NEXTOP = 1 
      NEW = 0
      IF(INI.EQ.0) THEN
*. Are there more configurations with the current number of open orbitals?
       IF(INUM_OP+1.LE.NCONF_PER_OPEN_GN(IOPEN+1,ISYM,ISPC)) THEN
         INUM_OP = INUM_OP+1
         NEXTOP = 0
         NEW = 1
       END IF
      END IF
*
      IF(NEXTOP.EQ.1) THEN
*. Find next number of open orbitals with a nonvanishing number of configurations
       JOPEN_MIN = IOPEN+1
       DO JOPEN = JOPEN_MIN,MAXOP
         IF(NCONF_PER_OPEN_GN(JOPEN+1,ISYM,ISPC).NE.0) THEN
           IOPEN = JOPEN
           INUM_OP = 1
           NEW = 1
           GOTO 1001
         END IF
       END DO
 1001  CONTINUE
      END IF ! look for next open 
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) '  Output from NEXT_CONF_IN_CONFSPC'
        WRITE(6,*) '  Output conf: IOPEN, INUM_OP ', IOPEN, INUM_OP
      END IF
*. Stop(!) is no new configurations were generated 
      IF(NEW.EQ.0) THEN
       WRITE(6,*) ' No new configurations in NEXT_CONF_IN_CONFSPC'
       WRITE(6,'(A,2I8)') 
     & ' Input IOPEN, INUM_OP = ', IOPEN_INI,INUM_OP_INI 
       STOP  ' No new configurations in NEXT_CONF_IN_CONFSPC'
      END IF
*. Find occupation of the new configuration
      IF(NEW.EQ.1) THEN
        NOCOB_L = (N_EL_CONF+IOPEN)/2
        IB_CONFLIST = KICONF_OCC_GN(ISYM,ISPC)
        IB_INLIST = IB_CONF_OCC_GN(IOPEN+1,ISYM,ISPC)
        IADR_INLIST = IB_INLIST + (INUM_OP-1)*NOCOB_L 
C            ICOPVE2(IIN,IOFF,NDIM,IOUT)
        CALL ICOPVE2(int_mb(IB_CONFLIST),IADR_INLIST,NOCOB_L,IOCC)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Next configuration '
          CALL IWRTMA(IOCC,1,NOCOB_L,1,NOCOB_L)
        END IF
      END IF ! NEW
*
      RETURN
      END
      SUBROUTINE GET_LBLK_CONF_BATCH(ICNF_INI,NCNF,LBLK_BAT,ISYM,ISPC,
     &           NSD_BAT_TOT,NCSF_BAT_TOT) 
*
* Obtain length of csf-expansion for configurations ICNF_INI - ICNF_INI+NCNF-1 
* in COnf space ISPC and SYM ISYM
*
* The total number of CSF's and SD's are also obtained
*
* Note: Offset (ICNF_INI is not changed in routine)
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'spinfo.inc'
*. Output
      INTEGER LBLK_BAT(*)
*. Start, absolute relative to beginning of batch
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from GET_LBLK_CONF'
      END IF
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' ICNF_INI = ',ICNF_INI
       WRITE(6,*) ' NCNF = ', NCNF
      END IF
*
      ICNF_ABS = ICNF_INI
      NCNF_LFT = NCNF
*
      NSD_BAT_TOT = 0
      NCSF_BAT_TOT = 0
*
      IF(NTEST.GE.1000) THEN
      WRITE(6,*) ' TEST: ISYM, ISPC = ',
     &                   ISYM, ISPC
      END IF
      DO IOPEN = 0, MAXOP
*.Current offset relative to start of block
        IB_REL = ICNF_ABS - IB_CONF_REO_GN(IOPEN+1,ISYM,ISPC) + 1
C?      WRITE(6,*) ' IOPEN, IB_REL = ', IOPEN, IB_REL
        IF(IB_REL.GE.1) THEN
         N_IN_IOP = 
     &   MIN(NCONF_PER_OPEN_GN(IOPEN+1,ISYM,ISPC)-IB_REL+1,NCNF_LFT)
         IF(N_IN_IOP.GE.1) THEN
           NCSF = NPCSCNF(IOPEN+1)
C?         WRITE(6,*) ' NCONF_PER..., IB_REL = ',
C?   &     NCONF_PER_OPEN_GN(IOPEN+1,ISYM,ISPC), IB_REL
*
           IB_BATCH = ICNF_ABS - ICNF_INI + 1
           IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' ICNF_ABS, ICNF_INI, IB_BATCH = ', 
     &                    ICNF_ABS, ICNF_INI, IB_BATCH 
           END IF
*
           CALL ISETVC(LBLK_BAT(IB_BATCH),NCSF,N_IN_IOP)
           ICNF_ABS = ICNF_ABS + N_IN_IOP
           NCNF_LFT = NCNF_LFT - N_IN_IOP
*
           NSD_BAT_TOT  = NSD_BAT_TOT  + N_IN_IOP*NPDTCNF(IOPEN+1)
           NCSF_BAT_TOT = NCSF_BAT_TOT + N_IN_IOP*NPCSCNF(IOPEN+1)

           IF(NCNF_LFT.EQ.0) GOTO 1001
         END IF ! There was something left for this IOPEN
        END IF! We are at or below the relevant IOPEN
      END DO ! loop over IOPEN
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Number of CSFs per conf for batch'
       CALL IWRTMA(LBLK_BAT,1,NCNF,1,NCNF)
*
       WRITE(6,'(A,2I9)') ' Total number of CSFs and SDs in batch ',
     & NCSF_BAT_TOT, NSD_BAT_TOT
      END IF
*
      RETURN
      END
      SUBROUTINE CSDTVC_CONF(C_SD,C_CSF,NOPEN,ISIGN,IAC,IWAY)
*
* Transform between CSF and SD form of the coefficients for a given 
* configuration with IOPEN open orbitals
*
* IWAY = 1 => CSF to SD
* IWAY = 2 => SD to CSF
* 
* IAC = 1 => Result is added to outputvector
*     = 2 => Result is copied to outputvector
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input or output
      DIMENSION C_SD(*),C_CSF(*)
*. Input signs for going to alpha-beta ordering
      INTEGER ISIGN(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' CSDTVC_CONF in action '
      END IF
*. Offset to given transformation matrix
      IOFFCD = 1
      DO IOPEN = 0, NOPEN-1
       IOFFCD = IOFFCD + NPCSCNF(IOPEN+1)*NPDTCNF(IOPEN+1)
      END DO
      NPCSF = NPCSCNF(NOPEN+1)
      NPSD  = NPDTCNF(NOPEN+1)
*
      IF(IAC.EQ.1) THEN
       FACTORC = 1.0D0
      ELSE
       FACTORC = 0.0D0
      END IF
*
      IF(IWAY.EQ.1) THEN
C            MATVCC2(A,VIN,VOUT,NROW,NCOL,ITRNS,FACIN)
        CALL MATVCC2(dbl_mb(KDTOC + IOFFCD -1),C_CSF,C_SD,NPSD,NPCSF,
     &                0,FACTORC)
*
        IF(NTEST.GE.10000) THEN
          WRITE(6,*) ' Result before MATVCC2: '
          CALL WRTMAT(C_SD,1,NPSD,1,NPSD)
        END IF
*
        CALL ISIGN_TIMES_REAL(ISIGN,C_SD,NPSD)
C            ISIGN_TIMES_REAL(ISIGN,VEC,NDIM)
      ELSE
        CALL ISIGN_TIMES_REAL(ISIGN,C_SD,NPSD)
C?      WRITE(6,*) ' ISIGN: '
C?      CALL IWRTMA(ISIGN,1,NPSD,1,NPSD)
C?      WRITE(6,*) ' Result after ISIGN '
C?      CALL WRTMAT(C_SD,1,NPSD,1,NPSD)
        CALL MATVCC2(dbl_mb(KDTOC + IOFFCD -1),C_SD,C_CSF,NPSD,NPCSF,
     &               1,FACTORC)
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from CSDTVC_CONF '
       WRITE(6,*) ' ======================= '
       WRITE(6,*) ' NOPEN, NPCSF, NPSD = ', NOPEN, NPCSF, NPSD
       IF(IWAY.EQ.1) THEN
        WRITE(6,*) ' CSF => SD transformation '
       ELSE 
        WRITE(6,*) ' SD => CSF transformation '
       END IF
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' CSF-vector '
        CALL WRTMAT(C_CSF,NPCSF,1,NPCSF,1)
        WRITE(6,*) ' SD-vector '
        CALL WRTMAT(C_SD,NPSD,1,NPSD)
        WRITE(6,*) ' The ISIGN array '
        CALL IWRTMA(ISIGN, 1, NPSD,1, NPSD)
        WRITE(6,*) 
        WRITE(6,*) ' IOPEN, IOFFCD = ',IOPEN, IOFFCD
        WRITE(6,*) ' The DTOC array '
        CALL WRTMAT(dbl_mb(KDTOC + IOFFCD -1),NPSD,NPCSF,
     &  NPSD, NPCSF)
      END IF
*
      RETURN
      END
      SUBROUTINE CNHCN_CSF_BLK(ICNL,IOPL,ICNR,IOPR,CNHCNM,IADOB,
     &           IPRODT,DTOC,I12OP,ISCR,SCR,ECORE,IONLY_DIAG,ISYMG,
     &           RJ, RK)
*
*. Obtain Hamiltonian matrix between configurations ICNL and ICNR
*  If IONLY_DIAG = 1, then only the diagonal CSF elements are 
*  calculated.
*
*  If ISYMG = 1, then only the lower half of the matrix is calculated
*  (ISYMG = 1 is pt only active in connection with the IONLY_DIAG option)
*
*. Jeppe Olsen, Jan 2012, updated June 2012
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'priunit.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
*. Coulomb and exchange integrals
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*
* IADOB is start of orbitals in configurations - correcting
* for that inactive orbitals may preceed these
*
* I12OP = 1  => Only one-electron operator
* I12OP = 2  => One- + two-electron operator
*. Input
*. Complete prototype det info and CSF-SD transformation
      DIMENSION IPRODT(*), DTOC(*)
*. Occupation of configurations
      INTEGER ICNL(*),ICNR(*)
*. Output: 
      DIMENSION CNHCNM(*)
*. Scratch: Length: INTEGER: (NDET_C + NDET_S)*N_EL_CONF + NDET_C + 6*NORB *. (NELEC .le. 2*NORB)
*                   REAL   : 2 NDETL x NDETR matrix (could be reduced)
*                    
      DIMENSION ISCR(*)
      DIMENSION SCR(*)
*
      REAL*8 INPROD
      ITIME_DETAILS = 1
*
      IF(ITIME_DETAILS.EQ.1) CALL LUCIAQENTER('CNHCNM')
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from  CNHCN_CSF_BLK'
        WRITE(6,*) '=========================='
        WRITE(6,*)
        WRITE(6,*) ' IOPL, IOPR = ', IOPL, IOPR
        WRITE(6,*) ' ECORE = ', ECORE
        WRITE(6,*) ' ISYMG, IONLY_DIAG = ', ISYMG, IONLY_DIAG
        WRITE(6,*) ' N_EL_CONF = ', N_EL_CONF
*
        ICLL = (N_EL_CONF-IOPL)/2
        ICLR = (N_EL_CONF-IOPR)/2
*
        NOBL = IOPL + ICLL
        NOBR = IOPR + ICLR
*
        WRITE(6,*) ' Left and right configurations '
        CALL IWRTMA(ICNL,1,NOBL,1,NOBL)
        CALL IWRTMA(ICNR,1,NOBR,1,NOBR)
      END IF
*
      ISYM = ISYMG
      IF(ISYMG.EQ.1.AND.IONLY_DIAG.EQ.0) THEN
        ISYM = 0
        WRITE(6,*) ' ISYM set to 0 as IONLY_DIAG = 0 '
      END IF
*
*. A bit of info on the configurations
      ITPL = IOPL + 1
      ITPR = IOPR + 1
*
      ICLL = (N_EL_CONF-IOPL)/2
      ICLR = (N_EL_CONF-IOPR)/2
*
      NOBL = IOPL + ICLL
      NOBR = IOPR + ICLR
*. Is really combinations, so
      NDETL = NPCMCNF(ITPL)
      NDETR = NPCMCNF(ITPR)
*
      NCSFL = NPCSCNF(ITPL)
      NCSFR = NPCSCNF(ITPR)
*
      NAEL = (N_EL_CONF + MS2)/2
      NBEL = (N_EL_CONF - MS2)/2
*. Offsets in CSF-SD transformation matrices
      IB_CDL = 1
      DO IOP = 0, IOPL-1
       IB_CDL = IB_CDL + NPCSCNF(IOP+1)*NPCMCNF(IOP+1)
      END DO
      IB_CDR = 1
      DO IOP = 0, IOPR-1
       IB_CDR = IB_CDR + NPCSCNF(IOP+1)*NPCMCNF(IOP+1)
      END DO
*
*. ======================================================
*. Construct in SCR the Hamiltonian block in the SD basis
*. ======================================================
*
      KLSDHSD = 1
      KLCSFHSD = KLSDHSD + NDETL*NDETR
      XDUM = 0.0D0
      IF(ITIME_DETAILS.EQ.1) CALL LUCIAQENTER('CNHCNL')
      CALL CNHCN_LUCIA(ICNL,IOPL,ICNR,IOPR,XDUM,SCR(KLSDHSD),XDUM,
     &     IADOB,IPRODT,I12OP,0,IDUM,ECORE,1,ISYM,RJ,RK,ISCR)
C     CNHCN_LUCIA(ICNL,IOPL,ICNR,IOPR,C,CNHCNM,SIGMA,
C    &           IADOB,IPRODT,I12OP,I_DO_ORBTRA,IORBTRA,
C    &           ECORE,IHORS,ISYM,RJ,RK,ISCR)
      IF(ITIME_DETAILS.EQ.1) CALL LUCIAQEXIT('CNHCNL')
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Hamiltonian block in SD basis before sign changes '
        IF(ISYM.EQ.0) THEN
          CALL WRTMAT(SCR(KLSDHSD),NDETL,NDETR,NDETL,NDETR)
        ELSE
          CALL PRSYM(SCR(KLSDHSD),NDETL)
        END IF
      END IF
*
* ========================================================
* Sign changes for going from String to Conf order of SDs
* ========================================================
* 
*. Obtain sign-change for going between SD and CSF-order
*
      KLIL = 1
      KLIR = KLIL + NDETL
      KLIFREE = KLIR + NDETR
           
C     SIGN_CONF_SD(ICONF,NOB_CONF,IOP,ISGN,IPDET_LIST,ISCR)
      CALL SIGN_CONF_SD(ICNL,NOBL,IOPL,ISCR(KLIL),IPRODT,ISCR(KLIFREE))
      CALL SIGN_CONF_SD(ICNR,NOBR,IOPR,ISCR(KLIR),IPRODT,ISCR(KLIFREE))
      IF(NTEST.GE.1000) THEN
      WRITE(6,*) ' IL and IR sign arrays '
        CALL IWRTMA(ISCR(KLIL),1,NDETL,1,NDETL)
        CALL IWRTMA(ISCR(KLIR),1,NDETR,1,NDETR)
      END IF
      IF(ISYM.EQ.0) THEN
*. full matrix form
       IB_HSD = 1
       DO J = 1, NDETR
         IF(J.GT.1) IB_HSD = IB_HSD + NDETL
         CALL ISIGN_TIMES_REAL(ISCR(KLIL),SCR(KLSDHSD-1+IB_HSD),NDETL)
         IF(ISCR(KLIR-1+J).EQ.1) THEN
           SIGN = 1.0D0
         ELSE
           SIGN = -1.0D0
         END IF
         CALL SCALVE(SCR(KLSDHSD-1+IB_HSD),SIGN,NDETL)
       END DO
      ELSE
*. matrix packed rowwise
       DO I = 1, NDETL
         IB_HSD = I*(I-1)/2 + 1
         CALL ISIGN_TIMES_REAL(ISCR(KLIR),SCR(KLSDHSD-1+IB_HSD),I)
         IF(ISCR(KLIL-1+I).EQ.1) THEN
           SIGN = 1.0D0
         ELSE
           SIGN = -1.0D0
         END IF
         CALL SCALVE(SCR(KLSDHSD-1+IB_HSD),SIGN,I)
       END DO
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Hamiltonian block in SD basis with sign changes '
        IF(ISYM.EQ.0) THEN
         CALL WRTMAT(SCR(KLSDHSD),NDETL,NDETR,NDETL,NDETR)
        ELSE
         CALL PRSYM(SCR(KLSDHSD),NDETL)
        END IF
      END IF
*
* ===============================
* Transform from CSF to SD basis
* ===============================
* 
      IF(IONLY_DIAG.EQ.0) THEN
*
*. Complete matrix in CSF basis
*
*. Transform first index to CSF basis
*. Offset to transformation matrix for IOPL/IOR open orbitals
        FACTORC = 0.0D0
        FACTORAB = 1.0D0
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' CSF-SD transformation '
          CALL WRTMAT(DTOC(IB_CDL),NDETL,NCSFL,NDETL,NCSFL)
        END IF
        CALL MATML7(SCR(KLCSFHSD),DTOC(IB_CDL),SCR(KLSDHSD),
     &       NCSFL,NDETR,NDETL,NCSFL,NDETL,NDETR,
     &       FACTORC,FACTORAB,1)
*. Transform second index to CSF basis
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' CSF-SD transformation '
          CALL WRTMAT(DTOC(IB_CDR),NDETR,NCSFR,NDETR,NCSFR)
        END IF
        CALL MATML7(CNHCNM,SCR(KLCSFHSD),DTOC(IB_CDR),
     &       NCSFL,NCSFR,NCSFL,NDETR,NDETR,NCSFR,
     &       FACTORC,FACTORAB,0)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) 
     &    ' Hamiltonian matrix betweeen two confs in CSF basis'
          CALL WRTMAT(CNHCNM,NCSFL,NCSFR,NCSFL,NCSFR)
        END IF
      ELSE
*. Obtain just the diagonal values - and store as diagonal in CNHCNM
*. DTOC(IB_CDL) SDHSD as sparse matrix multiply
*  For the diagonal to be meaningful it is assumed that the 
*  left and right occupation are identical 
*. Transpose the DTOC matrix for improved matrix multiply performance
C           TRPMT3(XIN,NROW,NCOL,XOUT)
       CALL TRPMT3(DTOC(IB_CDL),NDETL,NCSFL,SCR(KLCSFHSD))
       CALL COPVEC(SCR(KLCSFHSD),DTOC(IB_CDL),NDETL*NCSFL)
C           MULT_MAT_SPMAT_MAT(AOUT,AIN,X,NAOUT_R,NAOUT_C,NX_C,IAINPAK )
       IF(ITIME_DETAILS.EQ.1) CALL LUCIAQENTER('SPAMMU')
       CALL MULT_MAT_SPMAT_MAT(SCR(KLCSFHSD),SCR(KLSDHSD),DTOC(IB_CDL),
     &      NCSFL,NDETR,NDETL,1)
       IF(ITIME_DETAILS.EQ.1) CALL LUCIAQEXIT('SPAMMU')
       CALL TRPMT3(SCR(KLCSFHSD),NCSFL,NDETL,SCR(KLSDHSD))
*. And clean up
       CALL TRPMT3(DTOC(IB_CDL),NCSFL,NDETL,SCR(KLCSFHSD))
       CALL COPVEC(SCR(KLCSFHSD),DTOC(IB_CDL),NDETL*NCSFL)
       DO I = 1, NCSFL
         CNHCNM((I-1)*NCSFL+I) = INPROD(DTOC(IB_CDL+(I-1)*NDETL),
     &   SCR(KLSDHSD+(I-1)*NDETL),NDETL)
       END DO
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Diagonal elements created '
         CALL WRTDIA(CNHCNM,NCSFL,1)
       END IF
      END IF! only diagonal should be calculated
      IF(ITIME_DETAILS.EQ.1) CALL LUCIAQEXIT('CNHCNM')
*
      RETURN
      END
      SUBROUTINE CNHCN_LUCIA(ICNL,IOPL,ICNR,IOPR,C,CNHCNM,SIGMA,
     &           IADOB,IPRODT,I12OP,I_DO_ORBTRA,IORBTRA,
     &           ECORE,IHORS,ISYM,RJ,RK,ISCR)
*
* Obtain Hamiltonian matrix elements between L and R configurations 
* and multiply  elements of C to update SIGMA (if IHORS = 2)
* - evrything in SD basis
*
* The Phase of the input C is supposed to be in String def, 
* and Sigma is returned with the same phase convention
* Jeppe Olsen, Summer of '89, Revival in Summer of '11
*
* IF ISYM = 1, only the lower half of the matrix is explicitly 
* calculated
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'priunit.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
* IADOB is start of orbitals in configurations - correcting
* for that inactive orbitals may preceed these
*
* IHORS = 1: Construct Block of Hamilton matrix
* IHORS = 2: Construct block of Hamilton matrix times vector
*
* I12OP = 1  => Only one-electron operator
* I12OP = 2  => One- + two-electron operator
* I_DO_ORBTRA = 1: Orbital transformation with orbital IORBTRA
*. Input
*. Prototype info
      DIMENSION IPRODT(*)
*. Occupation of configurations
      INTEGER ICNL(*),ICNR(*)
      DIMENSION C(*)
*. Coulomb and exchange integrals (used for diagonal elements)
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*. Output: one of the following depending on IHORS
      DIMENSION SIGMA(*)
      DIMENSION CNHCNM(*)
*. Scratch: Length: INTEGER: (NDET_C + NDET_S)*N_EL_CONF + NDET_C + 6*NORB
*. (Remembering that NELEC .le. 2*NORB)
*                    
      DIMENSION ISCR(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CNHCN_LUCIA '
        WRITE(6,*) '=========================='
        WRITE(6,*)
        WRITE(6,*) ' IOPL, IOPR = ', IOPL, IOPR
        WRITE(6,*) ' ECORE, IHORS = ', ECORE, IHORS
      END IF
*
* 1 : Obtain determinants of configurations
*
      ITPL = IOPL + 1
      ITPR = IOPR + 1
*
      ICLL = (N_EL_CONF-IOPL)/2
      ICLR = (N_EL_CONF-IOPR)/2
*
      NOBL = IOPL + ICLL
      NOBR = IOPR + ICLR
*
      NDETL = NPCMCNF(ITPL)
      NDETR = NPCMCNF(ITPR)
*
      NCSFL = NPCSCNF(ITPL)
      NCSFR = NPCSCNF(ITPR)
*
      NAEL = (N_EL_CONF + MS2)/2
      NBEL = (N_EL_CONF - MS2)/2
*
      KLFREE = 1
*
      KLDTLA = KLFREE
      KLFREE = KLFREE + NDETL*NAEL
*
      KLDTLB = KLFREE
      KLFREE = KLFREE + NDETL*NBEL
*
      KLDTRA = KLFREE
      KLFREE = KLFREE + NDETR*NAEL
*
      KLDTRB = KLFREE
      KLFREE = KLFREE + NDETR*NBEL
*
*
*. Obtain determinants of Left and Right configurations
*
C     CNFSTRN(ICONF,NOB_CONF,IOP,IASTR,IBSTR,IPDET_LIST,ISCR)
CTIME CALL LUCIAQENTER('CNFSTR')
* 3* NELEC is used for scratch in CNFSTRN
      CALL CNFSTRN(ICNL,NOBL,IOPL,ISCR(KLDTLA),ISCR(KLDTLB),
     &     int_mb(KDFTP),ISCR(KLFREE))
*
      CALL CNFSTRN(ICNR,NOBR,IOPR,ISCR(KLDTRA),ISCR(KLDTRB),
     &     int_mb(KDFTP),ISCR(KLFREE))
CTIME CALL LUCIAQEXIT('CNFSTR')
*
*
* Hamiltonian matrix over determinants times C
*
COLD  IHORS = 2
      XDUM = 0.0D0
C     DIHDJ2_LUCIA_CONF
C    &(IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,IADOB,
C    & IHORS,HAMIL,C,SIGMA,IWORK,ISYM,ECORE,ICOMBI,PSIGN,
C    & NTERMS,NDIF0,NDIF1,NDIF2,IORBTRA,IORB)
      IF(PSSIGN.EQ.0) THEN
        ICOMBI_L = 0
      ELSE
        ICOMBI_L = 1
      END IF
      CALL DIHDJ2_LUCIA_CONF(ISCR(KLDTLA),ISCR(KLDTLB),NDETL,
     &                 ISCR(KLDTRA),ISCR(KLDTRB),NDETR,
     &                 NAEL,NBEL,IADOB,N_ORB_CONF,
     &                 IHORS,CNHCNM,C,SIGMA,
     &                 ISCR(KLFREE),ISYM,ECORE,ICOMBI_L,   
     &                 PSSIGN,NTERMS,NDIF0,NDIF1,NDIF2,
     &                 I12OP, I_DO_ORBTRA,IORBTRA,NTOOB,RJ,RK)
*
      IF( NTEST .GE. 1000 ) THEN
        IF(IHORS.EQ.1) THEN
          LUPRI = 6
          WRITE(LUPRI,*)
     &    ' SD-Hamiltonian matrix between two configurations'
          IF(ISYM.EQ.0) THEN
            CALL WRTMAT(CNHCNM,NDETL,NDETR,NDETL,NDETR)
          ELSE
            CALL PRSYM(CNHCNM,NDETL)
          END IF
        ELSE
          WRITE(6,*) ' Input C vector for conf'
          CALL WRTMAT(C,1,NDETR,1,NDETR)
          WRITE(6,*) ' Updated sigma-vector for conf'
          CALL WRTMAT(SIGMA,1,NDETL,1,NDETL)
        END IF
      END IF
c
      RETURN
      END
      SUBROUTINE CNFSTRN(ICONF,NOB_CONF,IOP,IASTR,IBSTR,IPDET_LIST,
     &ISCR)
*
* A configuration ICONF with NOB_CONF occupied orbitals is given
* (in compact form). 
* Obtain the alpha- and beta-strings of determinants of the 
* configuration.
*
*. Note: Sign changes for switching from configuration to alpha-beta
*        order are not stored.
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER ICONF(NOB_CONF)
*. The complete list of prototype determinants
      INTEGER IPDET_LIST(*)
*. Output
      INTEGER IASTR(*), IBSTR(*)
*. Scratch: 3*NELEC 
      INTEGER ISCR(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from CNFSTRN '
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' NOB_CONF, IOP = ', NOB_CONF, IOP 
      END IF
*
      NAEL = (N_EL_CONF + MS2)/2
      NBEL = N_EL_CONF - NAEL
      IF(NTEST.GE.1000) WRITE(6,*) ' NAEL, NBEL = ',NAEL, NBEL
*
*. Division of scratch
      KLSCR1 = 1
      KLSCR2 = KLSCR1 + (NAEL+NBEL)
      KLFREE = KLSCR2 + (NAEL+NBEL)
*
C     NPDET = NPDTCNF(IOP + 1)
      NPDET = NPCMCNF(IOP + 1)
      IBPDET = 1
      DO IOPEN = 0, IOP - 1
C      IBPDET = IBPDET + NPDTCNF(IOPEN+1)*IOPEN
       IBPDET = IBPDET + NPCMCNF(IOPEN+1)*IOPEN
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NPDET, IBPDET = ',
     &                               NPDET, IBPDET
*
      DO IPDET = 1, NPDET
*. Obtain occupation corresponding to prototype determinant IPDET
C     GETDET_FOR_CONF(IOCC,NOCOB,NOPEN,NELEC,JPDET,IPDET,
C    &           IDET,I_DO_REO,IREO)
*
        CALL GETDET_FOR_CONF(ICONF,NOB_CONF,IOP,N_EL_CONF,IPDET,
     &       IPDET_LIST(IBPDET),ISCR(KLSCR1),1,IREO_SPCP_OB_ON)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*)  ' Determinant from GETDET_FOR_CONF'
          CALL IWRTMA(ISCR(KLSCR1),1,N_EL_CONF,1,N_EL_CONF)
        END IF
* Obtain the alpha- and beta-strings of the configuration
C            DETSTR(IDET,IASTR,IBSTR,NAEL,NBEL,ISIGN,IWORK,IPRNT)
        CALL DETSTR(ISCR(KLSCR1),
     &      IASTR((IPDET-1)*NAEL+1),IBSTR((IPDET-1)*NBEL+1),
     &      NAEL,NBEL,IDUM,ISCR(KLSCR2),NTEST)
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Determinants in AB form for configuration '
       DO IPDET = 1, NPDET
         WRITE(6,*)
         WRITE(6,'(A,I7)') ' Determinant number ', IPDET
         CALL WRT_ABSTR(IASTR((IPDET-1)*NAEL+1),IBSTR((IPDET-1)*NBEL+1),
     &        NAEL,NBEL)
       END DO
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_ABSTR(IASTR,IBSTR,NAEL,NBEL)
*
* Print alpha- and beta-strings
*
*. Jeppe Olsen, July 2011 (sic)
*
      INCLUDE 'implicit.inc'
      INTEGER IASTR(NAEL), IBSTR(NBEL)
*
      WRITE(6,'(A, 20(1X,I3),/,(15X, 20(1X,I3)))')
     &' Alpha-string: ', (IASTR(IEL),IEL=1, NAEL)
      WRITE(6,'(A, 20(1X,I3),/,(15X, 20(1X,I3)))')
     &' Beta-string:  ', (IBSTR(IEL),IEL=1, NBEL)
*
      RETURN
      END
      SUBROUTINE SIGN_CONF_SD
     &(ICONF,NOB_CONF,IOP,ISGN,IPDET_LIST,ISCR)
*
* A configuration ICONF with NOB_CONF occupied orbitals is given
* (in compact form). 
*
*. Obtain the phase required  for changing between configuration 
*. and alpha-beta order of strings
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'orbinp.inc'
*. Specific input
      INTEGER ICONF(NOB_CONF)
*. The complete list of prototype determinants
      INTEGER IPDET_LIST(*)
*. Output
      DIMENSION ISGN(*)
*. Scratch: NELEC + NAEL+NBEL
      INTEGER ISCR(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from SIGN_CONF_SD '
        WRITE(6,*) ' ======================== '
        WRITE(6,*) ' NOB_CONF, IOP = ', NOB_CONF, IOP 
        WRITE(6,*) ' Input configuration: '
        CALL IWRTMA(ICONF,1,NOB_CONF,1,NOB_CONF)
        WRITE(6,*) ' First 4 elements of IPDET_LIST '
        CALL IWRTMA(IPDET_LIST,1,4,1,4)
      END IF
*
      NAEL = (N_EL_CONF + MS2)/2
      NBEL = N_EL_CONF - NAEL
      IF(NTEST.GE.1000) WRITE(6,*) ' NAEL, NBEL = ',NAEL, NBEL
*
*. Division of scratch
      KLSCR1 = 1
      KLSCR2  = KLSCR1 + (NAEL+NBEL)
      KLSCR3A = KLSCR2 + (NAEL+NBEL)
      KLSCR3B = KLSCR3A + NAEL
*
C     NPDET = NPDTCNF(IOP + 1)
      NPDET = NPCMCNF(IOP + 1)
      IBPDET = 1
      DO IOPEN = 0, IOP - 1
C      IBPDET = IBPDET + NPDTCNF(IOPEN+1)*IOPEN
       IBPDET = IBPDET + NPCMCNF(IOPEN+1)*IOPEN
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NPDET, IBPDET = ',
     &                               NPDET, IBPDET
*
      DO IPDET = 1, NPDET
         IF(NTEST.GE.1000) WRITE(6,*) ' Info for IPDET = ', IPDET
*. Obtain spin-projections corresponding to prototype determinant IPDET
        CALL GETDET_FOR_CONF(ICONF,NOB_CONF,IOP,N_EL_CONF,IPDET,
     &       IPDET_LIST(IBPDET),ISCR(KLSCR1),1,IREO_SPCP_OB_ON)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*)  ' Determinant from GETDET_FOR_CONF'
          CALL IWRTMA(ISCR(KLSCR1),1,N_EL_CONF,1,N_EL_CONF)
        END IF
* Obtain the sign change for going to a b order
C            DETSTR(IDET,IASTR,IBSTR,NAEL,NBEL,ISIGN,IWORK,IPRNT)
        CALL DETSTR(ISCR(KLSCR1),
     &      ISCR(KLSCR3A),ISCR(KLSCR3B),
     &      NAEL,NBEL,ISGN(IPDET),ISCR(KLSCR2),NTEST)
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Sign changes for switching between SD and conf'
       CALL IWRTMA(ISGN,1,NPDET,1,NPDET)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_DIM_MINMAX_SPACE(MIN_OCC,MAX_OCC,MINMAX_ORB,NORB,
     &           ISYM,NCONF,NCSF,NSD,NCM,LCONFOCC,NCONF_AS)
* 
* Find Number of configurations, CSF's, SDs, CMs for a MINMAX expansion
* defined by MIN_OCC, MAX_OCC
*
* It is assumed that the Prototype info has been set up
*
* The number of configurations for each number of open orbitals is 
* returned in NCONF_PER_OPEN(*,ISYM)
*
* Jeppe Olsen, July 16, 2011
*
* Last modification; July 2013; Jeppe Olsen; MINMAX_ORB added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
*
      INTEGER MIN_OCC(NORB),MAX_OCC(NORB),MINMAX_ORB(NORB)
*
      IDUM = 0  
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTDMMM')
*
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' GET_DIM_MINMAX in action '
       WRITE(6,*) ' ======================= '
      END IF
*
C     GEN_CONF_FOR_MAXMIN_OCC(IOCC_MIN,IOCC_MAX,NORBL,
C    &    IORB_OFF,
C    &    INITIALIZE_CONF_COUNTERS,
C    &    ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,
C    &    NCONF_OP,NCONF,IBCONF_REO,IBCONF_OCC,ICONF,
C    &    IDOREO,IZ_CONF,IREO,NCONF_ALL_SYM,IREO_MNMX_OB_NO)

*
      CALL GEN_CONF_FOR_MAXMIN_OCC(MIN_OCC,MAX_OCC,NORB,
     &     IB_ORB_CONF,1,ISYM,MINOP,MAXOP,NSMST,1,
     &     NCONF_PER_OPEN(1,ISYM),NCONF,IDUM,IDUM,IDUM,0,IDUM,IDUM,
     &     NCONF_AS,MINMAX_ORB)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) 
     & ' Number of configurations for various number of open orbitals'
       CALL IWRTMA(NCONF_PER_OPEN(1,ISYM),MAXOP+1,1,MXPORB+1,MXPCSM)
      END IF
* We now have in NCONF_PER_OPEN the number of configurations for the 
* various number of open orbitals. Determine the number of 
* CSF's and SD's
*
      NCSF = IINPROD(NCONF_PER_OPEN(1,ISYM),NPCSCNF,MAXOP+1)
      NSD  = IINPROD(NCONF_PER_OPEN(1,ISYM),NPDTCNF,MAXOP+1)
      NCM  = IINPROD(NCONF_PER_OPEN(1,ISYM),NPCMCNF,MAXOP+1)
*. Various offsets
      CALL INFO_CONF_LIST(ISYM,LCONFOCC,NCONF_TOT)

      IF(NTEST.GE.10) THEN
       WRITE(6,'(A,4(2X,I9))') 
     & ' Number of CONFs, CSFs, CMs and SDs in MINMAX space ',
     &   NCONF,NCSF, NCM,NSD
       WRITE(6,'(A,I9)') 
     & ' Number of CONFs, all symmetries: ', NCONF_AS
       WRITE(6,'(A,I9)') ' Length of list of conf occs ', LCONFOCC
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTDMMM')
      RETURN
      END
      SUBROUTINE MINMAX_PER_SYM(MIN_OCC,MAX_OCC,MIN_PER_SYM,MAX_PER_SYM)
*
* Min and Max for each symmetry for a MINMAX space 
*
*. It is assumed that the configuration space is defined by a single
*. orbital space, so orbitals with the same symmetry are grouped
* together.
*
*. Jeppe Olsen, July 2011
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'  
      INCLUDE 'orbinp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
*.Input
      INTEGER MIN_OCC(*), MAX_OCC(*)
*.Output
      INTEGER MIN_PER_SYM(NSMOB), MAX_PER_SYM(NSMOB)
*
      DO ISYM = 1, NSMOB
*. First and number of orbital of this symmetry
        NFOUND = 0
        IFIRST = 0
        DO IORB = 1, N_ORB_CONF 
          IF(ISMFTO(IB_ORB_CONF-1+IORB).EQ.ISYM) THEN
            IF(NFOUND.EQ.0) IFIRST = IORB
            NFOUND = NFOUND + 1
          END IF
        END DO
*
        IF(NFOUND.EQ.0) THEN
          MIN_PER_SYM(ISYM) = 0
          MAX_PER_SYM(ISYM) = 0
        ELSE
          IF(IFIRST.EQ.1) THEN
*. First symmetry occuring 
            MIN_PER_SYM(ISYM) = MIN_OCC(NFOUND)
            MAX_PER_SYM(ISYM) = MAX_OCC(NFOUND)
          ELSE
            MIN_PER_SYM(ISYM) = MAX(0,MIN_OCC(NFOUND)-MAX_OCC(IFIRST-1))
            MAX_PER_SYM(ISYM) = 
     &      MIN(2*NFOUND,MAX_OCC(IFIRST+NFOUND-1)-MAX_OCC(IFIRST-1))
          END IF ! IFIRST = 1
        ENDIF !NFOUND = 0
      END DO ! loop over ISYM
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Min and max per symmetry '
        WRITE(6,*) ' ======================== '
        WRITE(6,*)
        CALL IWRTMA(MIN_PER_SYM,1,NSMOB,1,NSMOB)
        CALL IWRTMA(MAX_PER_SYM,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
      SUBROUTINE PROTO_CSF_DIM
* 
* Obtain information about the memory requirements for the prototype
* CSF info
*
*. Jeppe Olsen, Dec. 2011 (on train from Aarhus to Odense)
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'orbinp.inc'   
      INCLUDE 'csfbas.inc'
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from PROTO_CSF_DIM '
        WRITE(6,*) ' =========================='
        WRITE(6,*)
      END IF
*
      NELEC = NACTEL
*
*. Max. and min. number of open orbitals -as we want to treat 
*. several spins, detailed selection based on these have been removed.
* and electrons
C?    IF(MOD(MULTS-1,2).EQ.1) THEN
C?     MINOP = 1
C?    ELSE
C?     MINOP = 0
C?    END IF
*. Well, use MS2
      MINOP = ABS(MS2)
*. If only a single spin is used, we may use this
      MINOP = MULTS - 1
      IF(NTEST.GE.10)
     &WRITE(6,*) ' Min number of open orbitals ', MINOP
*. Max number of open orbitals - Depends on the actual form 
* of the CI spaces etc, is therefore read in from SPINFO from value
* determine by explicit check of occ spaces
      MAXOP = NOPEN_MAX
      IF(NTEST.GE.10)
     &WRITE(6,*) ' Max number of open orbitals ', MAXOP
*
* Number of Prototype SD's, CSF's ....
*
      CALL PROTO_CONF_DIM(NELEC)
*. Length of arrays for prototype info: Output is stored in CSFBAS
      CALL DIM_PROTO_ARRAYS
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from PROTO_CSF_DIM'
        WRITE(6,*) ' LPDT_OCC,LPCS_OCC,LPDTOC,MXPTDT,MXPOCCBL ',
     &  LPDT_OCC,LPCS_OCC,LPDTOC,MXPTDT,MXPOCCBL
      END IF
*
      RETURN
      END 
      SUBROUTINE PART_CIV_OCC(IDC,IBLTP,ISM,
     &           NSMST,MXLNG,ISMOST,NSSOA,NSSOB,
     &           NOCCLS_ACT,IOCCLS_ACT,
     &           N_ABSPGP_FOR_OCCLS,IB_ABSPGP_FOR_OCCLS,
     &           I_ABSPGP_FOR_OCCLS,
     &           NSD_PER_OCCLS,
     &           NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &           NOCCLS_BAT,IBOCCLS_BAT)
*
* Partition a CI vector into batches of blocks. 
* The length of a batch must be atmost MXLNG 
*
* Occupation class driven version
*
* IF ICOMP. eq. 1 the complete CI vector is constructed in one batch
*
*. Output 
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start 
*          of batch
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta sym
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta 
*   IBATCH(5,*) : Offset of block with respect to start of batch in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of batch in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form                   
*   IBATCH(8,*) : Length of block, packed form 
*
*  NOCCLS_BAT(IBAT): Number of occupation classes in batch IBAT
*  IBOCCLS_BAT(IBAT): First occupation block (IN IOCCLS_ACT) of batch
*
*
* Jeppe Olsen, Feb. 2012
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*. Info on the occupation classes
      INTEGER NSD_PER_OCCLS(NSMST,*)
*. Occupation => supergroup mappings
      INTEGER N_ABSPGP_FOR_OCCLS(*), IB_ABSPGP_FOR_OCCLS(*)
      INTEGER I_ABSPGP_FOR_OCCLS(2,*)
* The Active supergroups
      INTEGER IOCCLS_ACT(NOCCLS_ACT)
    
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      INTEGER NOCCLS_BAT(*)
      INTEGER IBOCCLS_BAT(*)
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================'
        WRITE(6,*) '   PART_CIV_OCC  '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*) ' MXLNG = ', MXLNG
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ISMOST array '
        CALL IWRTMA(ISMOST,1,NSMST,1,NSMST) 
        WRITE(6,*) ' IBLTP array '
        CALL IWRTMA(IBLTP,1,NSMST,1,NSMST) 
        WRITE(6,*) ' ICOMP = ', ICOMP
        WRITE(6,*) ' IOCCLS_ACT(NOCCLS_ACT): '
        CALL IWRTMA(IOCCLS_ACT,1,NOCCLS_ACT,1,NOCCLS_ACT)
      END IF
*
*. Loop over occ classes, generate  TTSS blocks and batch them
*
      LBLOCKP = 0
      LENGTHP_BAT = 0
      LENGTH_BAT = 0
      NNOCCLS_BAT = 0
      NBATCH = 0
      NBLK = 0
      NBLK_BAT = 0   
      I1BATCH(1) = 1
      IBOCCLS_BAT(1) = 1
*. Force one occupation class per batch?
      I_ONE_OCCLASS_PER_BATCH = 1
      DO IIOCCLS = 1, NOCCLS_ACT
       IOCCLS = IOCCLS_ACT(IIOCCLS)
       N_AB  =  N_ABSPGP_FOR_OCCLS(IOCCLS)
       IB_AB = IB_ABSPGP_FOR_OCCLS(IOCCLS)
       LENGTH_OCC = NSD_PER_OCCLS(ISM,IOCCLS)
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,'(A,5(2X,I5))') 
     &              ' IIOCCLS, IOCCLS, N_AB, IB_AB, LENGTH_OCC',
     &                IIOCCLS, IOCCLS, N_AB, IB_AB, LENGTH_OCC
         WRITE(6,'(A,3(2X,I6))') ' LENGTH_OCC, LENGTHP_BAT, MXLNG = ',
     &                             LENGTH_OCC, LENGTHP_BAT, MXLNG
       END IF
       IF(LENGTH_OCC + LENGTHP_BAT.GT.MXLNG.AND.ICOMP.NE.1.OR.
     &    (I_ONE_OCCLASS_PER_BATCH.EQ.1.AND.IIOCCLS.NE.1)) THEN
*. Not enough memory for adding occ block
        IF(NNOCCLS_BAT.EQ.0) THEN
         WRITE(6,*) 
     &   ' Error: A single occ class cannot be stored in one batch'
         WRITE(6,*) 
     &   ' Number and length of occ class: ', IOCCLS,LENGTH_OCC
         WRITE(6,*) ' Insufficient space detected in PART_CIV'
         CALL MEMCHK
         STOP
     &   ' Error: A single occ class cannot be stored in one batch'
        ELSE
*. Finalize previous batch and reset batch pointers
         NBATCH = NBATCH + 1
         LBATCH(NBATCH) = NBLK_BAT
         LEBATCH(NBATCH) = LENGTHP_BAT       
         I1BATCH(NBATCH+1) = NBLK + 1
         NOCCLS_BAT(NBATCH) = NNOCCLS_BAT
         IBOCCLS_BAT(NBATCH+1) = 
     &   IBOCCLS_BAT(NBATCH) + NNOCCLS_BAT
C?       WRITE(6,*) 
C?   &   ' NBATCH,IBOCCLS_BAT(NBATCH+1),IBOCCLS_BAT(NBATCH),NNOCLS_BAT'
C?       WRITE(6,*) 
C?   &     NBATCH, IBOCCLS_BAT(NBATCH+1),IBOCCLS_BAT(NBATCH),NNOCLS_BAT
     
*
         LENGTH_BAT = 0
         LENGTHP_BAT= 0
         NBLOCK = 0
         NNOCCLS_BAT = 0
         NBLK_BAT = 0
        END IF !NNOCLS_BAT = 0
       END IF !New batch required
*. Add occupation class to batch
       NNOCCLS_BAT = NNOCCLS_BAT + 1
       DO IAB = 1, N_AB
        IA = I_ABSPGP_FOR_OCCLS(1,IAB+IB_AB-1)
        IB = I_ABSPGP_FOR_OCCLS(2,IAB+IB_AB-1)
        IF(NTEST.GE.10000)
     &  WRITE(6,*) '  IAB, IA, IB = ',  IAB, IA, IB
        DO IASM = 1, NSMST
         IBSM = ISMOST(IASM)
CERR MAY13 IF(IBLTP(IASM).EQ.1.OR.(IBLTP(IASM).EQ.2.AND.IA.GE.IB)) THEN
         IF(IDC.EQ.2) THEN
            IF(IA.LT.IB) GOTO 1000
            IF(IA.EQ.IB.AND.IASM.LT.IBSM) GOTO 1000
         END IF
*. Block should be enrolled
           NBLK = NBLK + 1
           NBLK_BAT = NBLK_BAT + 1
*. Length
           NSTA = NSSOA(IASM,IA)
           NSTB = NSSOB(IBSM,IB)
           LBLOCK= NSTA*NSTB
           IF(IDC.EQ.1.OR.IA.GT.IB.OR.(IA.EQ.IB.AND.IASM.GT.IBSM)) THEN
            LBLOCKP = NSTA*NSTB
           ELSE IF(IDC.EQ.2.AND.IA.EQ.IB.AND.IASM.EQ.IBSM) THEN
            LBLOCKP = NSTA*(NSTA+1)/2
           END IF
*
           IF(NTEST.GE.10000) THEN
             WRITE(6,*) ' NBLK, IA, IB, IASM = ',
     &                    NBLK, IA, IB, IASM
             WRITE(6,*) ' NSTA, NSTB, LBLOCKP = ',
     &                    NSTA, NSTB, LBLOCKP
           END IF
           IBATCH(1,NBLK) = IA
           IBATCH(2,NBLK) = IB
           IBATCH(3,NBLK) = IASM
           IBATCH(4,NBLK) = IBSM
           IBATCH(5,NBLK) = LENGTH_BAT+1
           IBATCH(6,NBLK) = LENGTHP_BAT+1
           IBATCH(7,NBLK) = LBLOCK     
           IBATCH(8,NBLK) = LBLOCKP     
*
           LENGTHP_BAT = LENGTHP_BAT + LBLOCKP
           LENGTH_BAT = LENGTH_BAT + LBLOCK
C        END IF !Block should be included
 1000    CONTINUE
        END DO !Loop over IASM
       END DO! Loop over ABspgp in OCCLS
      END DO ! Loop over occls
*
*. Save  last batch 
*
      IF(NBLK_BAT.NE.0) THEN
         NBATCH = NBATCH + 1
         LBATCH(NBATCH) = NBLK_BAT
         LEBATCH(NBATCH) = LENGTHP_BAT       
C?       I1BATCH(NBATCH+1) = NBLK + 1
         NOCCLS_BAT(NBATCH) = NNOCCLS_BAT
      END IF
*
      IF(NTEST.GE.100) THEN
C?      WRITE(6,*) 'Output from PART_CIV'
C?      WRITE(6,*) '====================='
        WRITE(6,*)
        WRITE(6,*) ' Number of batches ', NBATCH    
        IBLOCKT = 0
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lengths of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            IBLOCKT = IBLOCKT + 1
            WRITE(6,'(10X,I5,2X,4I3,4I8)') 
     &      IBLOCKT,(IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
*
        WRITE(6,*)
        WRITE(6,*) ' Info on the batching of occ classes '
        WRITE(6,*) ' ==================================== '
        WRITE(6,*)
        WRITE(6,*) ' Batch   Number  Offset   '
        WRITE(6,*) ' ========================='
        DO JBATCH = 1, NBATCH
         WRITE(6,'(1H , I5,2X,I5,3X,I5)')
     &   JBATCH, NOCCLS_BAT(JBATCH), IBOCCLS_BAT(JBATCH)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE OCCLS_IN_CISPACE(NOCCLS_ACT,IOCCLS_ACT,
     &           NOCCLS,IOCCLS_OCC,NGAS,
     &           NMNMX_SPC,IMNMX_SPC,MNMX_OCC,ISPC)
*
* Obtain the occupation classes that is in CI space
* defined by NMNMX, IMNMX, MNMX_OCC, ISPC
*
*. Jeppe Olsen, Feb. 2012
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
#include "global.fh"
*. Input
      INTEGER IMNMX_SPC(NMNMX_SPC),MNMX_OCC(MXPNGAS,*)
      INTEGER IOCCLS_OCC(NGAS,NOCCLS)
      
*. Output
      INTEGER IOCCLS_ACT(*)
*
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' OCCLS_IN_CISPACE reporting'
      END IF
*
      NOCCLS_ACT = 0
      DO IOCCLS = 1, NOCCLS
        IM_IN = IS_OCCLS_IN_CISPACE(IOCCLS_OCC(1,IOCCLS),
     &          NMNMX_SPC,IMNMX_SPC,MNMX_OCC,ISPC)
        IF(IM_IN.EQ.1) THEN
          NOCCLS_ACT = NOCCLS_ACT + 1
          IOCCLS_ACT(NOCCLS_ACT) = IOCCLS
        END IF
      END DO
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Output from OCCLS_IN_CISPACE:'
        WRITE(6,*) 
        WRITE(6,*) ' Number of occlasses in CI-space ', NOCCLS_ACT
        WRITE(6,*) ' The active occupation classes'
        CALL IWRTMA(IOCCLS_ACT,1,NOCCLS_ACT,1,NOCCLS_ACT)
      END IF
*
      RETURN
      END
      FUNCTION IS_OCCLS_IN_CISPACE(IOCCLS,
     &         NMNMX_SPC,IMNMX_SPC,MNMX_OCC,ISPCSPC)
*
* Is occupation class IOCCLS in CI space defined by NMNMX,IMNMX,MNMX_OCC,
* ISPC
* 
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER IMNMX_SPC(NMNMX_SPC),MNMX_OCC(MXPNGAS,2,*)
      INTEGER IOCCLS(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IS_OCCLS_IN_CISPACE reporting '
        WRITE(6,*) ' The MINMAX spaces: '
        DO ISPC = 1, NMNMX_SPC
          WRITE(6,*) ' MINMAX space ', ISPC
          CALL IWRTMA(MNMX_OCC(1,1,ISPC),NGAS,2,MXPNGAS,2)
        END DO
      END IF
*
      INCLUDE = 0
      DO ISPC = 1, NMNMX_SPC
        IISPC = IMNMX_SPC(ISPC)
C?      WRITE(6,*) ' IISPC, ISPC = ', IISPC, ISPC
        IAMOKAY = 1
        IEL = 0
        DO IGAS = 1, NGAS
          IEL = IEL + IOCCLS(IGAS)
          IF(NTEST.GE.1000) WRITE(6,*) ' IGAS, IEL = ', IGAS,IEL,
C?        WRITE(6,*) ' MIN, MAX = ',
     &    MNMX_OCC(IGAS,1,IISPC), MNMX_OCC(IGAS,2,IISPC)
          IF(IEL.LT.MNMX_OCC(IGAS,1,IISPC).OR.
     &       IEL.GT.MNMX_OCC(IGAS,2,IISPC)    )  IAMOKAY = 0
C?        WRITE(6,*) ' IAMOKAY = ', IAMOKAY
        END DO! Loop over IGAS
        IF(IAMOKAY.EQ.1) INCLUDE = 1
      END DO! Loop over ISPC
*
      IF(INCLUDE.EQ.1.AND.I_CHECK_ENSGS.EQ.1) THEN
        CALL CHECK_IS_OCC_IN_ENGSOCC(IOCCLS,ISPCSPC,IM_IN)
        IF(IM_IN.EQ.1) INCLUDE = 1
      END IF
*
      IS_OCCLS_IN_CISPACE = INCLUDE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from IS_OCCLS_IN_CI_SPACE '
        WRITE(6,*) ' Occ lass in action: '
        CALL IWRTMA(IOCCLS,1,NGAS,1,NGAS)
        IF(INCLUDE.EQ.1) THEN
          WRITE(6,*) ' Is in space!'
        ELSE
          WRITE(6,*) ' Is not in space!'
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE Z_BLKFO_FOR_CISPACE(ISPC,ISM,LBLOCK,ICOMP,
     &           IPRNT,NBLOCK,NBATCH,
     &           IOIO,IBLTP,NOCCLS_ACT,IOCCLS_ACT,
     &           NBLK_PER_BATCH,NELMNT_PER_BATCH,LEN_PER_BLK,
     &           IB_FOR_BATCH, IBLOCKFO,NOCCLS_BAT,
     &           IBOCCLS_BAT,ILTEST)


*
* Construct information on allowed combination of blocks 
* and the batches for CI-space ISPC and symmetry ISM
*
*
* Jeppe Olsen, Cleaning a bit up in Feb. 2012, Geneva
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'csm.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
*. Output
*
*.Allowed combination of alpha- and beta-supergroups
      INTEGER IOIO(*)
*. Types of the various symmetry blocks
      INTEGER IBLTP(NSMST)
*. Info on the active occupation classes
      INTEGER IOCCLS_ACT(*)
*. Number of blocks per batch (CLBT)
      INTEGER NBLK_PER_BATCH(*)
*. Number of elements per batch (CLEBT)
      INTEGER NELMNT_PER_BATCH(*)
*. First block in a batch
      INTEGER IB_FOR_BATCH(*)
*. And the collected information about the blocks 
      INTEGER IBLOCKFO(8,*)
*. Length of each block (packed)
      INTEGER LEN_PER_BLK(*)
*. Number of occupation classes per batch and offset
      INTEGER NOCCLS_BAT(*), IBOCCLS_BAT(*)
*
      NTEST = 000
!      IPRNT = 0
      NTEST = MAX(NTEST,IPRNT)
*
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Output from Z_BLKFO_FOR_CISPACE '
        WRITE(6,*) ' ================================'
        WRITE(6,*)
        WRITE(6,*) ' Generation of blocks and batches for:'
        WRITE(6,'(A,6X,I2)') ' Space of state:', ISPC
        WRITE(6,'(A,6X,I2)') ' Symmetry of state:', ISM
        WRITE(6,*)
C       WRITE(6,*) ' ILTEST, ISIMSYM =  ', ILTEST, ISIMSYM
        WRITE(6,*) ' NOCSF = ', NOCSF
      END IF
*
      IF(ISPC.LE.NCMBSPC) THEN
        IATP = 1
        IBTP = 2
      ELSE
        IATP = IALTP_FOR_GAS(ISPC)
        IBTP = IBETP_FOR_GAS(ISPC)
      END IF
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0)
     &    WRITE(6,*) 'ZTEST: IATP, IBTP = ', IATP, IBTP

      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      IB_A = IBSPGPFTP(IATP)
      IB_B = IBSPGPFTP(IBTP)
* 
*. Information about block structure- needed by new PICO2 routine.
*. Memory for partitioning of C vector

      CALL IAIBCM(ISPC,IOIO)
*. option KSVST not active so
      IDUM = 1
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,IBLTP,IDUM)
*. The active occupation classes
      IF(ISPC.LE.NCMBSPC) THEN
        CALL OCCLS_IN_CISPACE(NOCCLS_ACT,IOCCLS_ACT,
     &       NOCCLS_MAX,int_mb(KIOCCLS),NGAS,
     &       LCMBSPC(ISPC),ICMBSPC(1,ISPC),IGSOCCX,ISPC)
*. Number of configurations per number of open orbitals
        IF(NOCSF.EQ.0) THEN
          CALL GET_NCONF_PER_OPEN_FOR_SUM_OCCLS(NCONF_PER_OPEN(1,ISM),
     &         MAXOP,NOCCLS_ACT,IOCCLS_ACT,ISM,int_mb(KNCN_PER_OP_SM),
     &         NIRREP)
*. And offsets
          CALL INFO_CONF_LIST(ISM,LENGTH_LIST_LOC, NCONF_TOT_LOC)
        END IF
      END IF
*
*. Batches of blocks
*. =================
*
      I_DO_BLOCKS = 1
      IF(NOCSF.EQ.0.AND.ICNFBAT.EQ.2.AND.ICOMP.EQ.1) I_DO_BLOCKS = 0
      CALL MEMCHK2('Z_BEPA')
      IF(NOCSF.EQ.0.AND.I_DO_BLOCKS.EQ.1) THEN
        CALL PART_CIV_OCC(IDC, IBLTP,ISM,NSMST,
     &       LBLOCK,ISMOST(1,ISM),int_mb(KNSTSO(IATP)),
     &       int_mb(KNSTSO(IBTP)),
     &       NOCCLS_ACT,IOCCLS_ACT,
     &       int_mb(KNABSPGP_FOR_OCCLS),int_mb(KIBABSPGP_FOR_OCCLS),
     &       int_mb(KIABSPGP_FOR_OCCLS),
     %       int_mb(KNSD_FOR_OCCLS),
     &       NBATCH,NBLK_PER_BATCH,NELMNT_PER_BATCH,
     &       IB_FOR_BATCH, IBLOCKFO,ICOMP,
     &       NOCCLS_BAT,IBOCCLS_BAT)
C    PART_CIV_OCC(IDC,IBLTP,ISM,
C    &           NSMST,MXLNG,ISMOST,NSSOA,NSSOB,
C    &           NOCCLS_ACT,IOCCLS_ACT,
C    &           N_ABSPGP_FOR_OCCLS,IB_ABSPGP_FOR_OCCLS,
C    &           I_ABSPGP_FOR_OCCLS,
C    &           NCM_PER_OCCLS,
C    &           NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
C    &           NOCCLS_BAT,IB_OCCLS_BAT)
      ELSE
        CALL PART_CIV2(IDC,IBLTP,int_mb(KNSTSO(IATP)),
     &       int_mb(KNSTSO(IBTP)),
     &       NOCTPA,NOCTPB,NSMST,LBLOCK,IOIO,
     &       ISMOST(1,ISM),
     &       NBATCH,NBLK_PER_BATCH,NELMNT_PER_BATCH,
     &       IB_FOR_BATCH, IBLOCKFO,ICOMP,ISIMSYM)
C            PART_CIV2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
C    &                  NSMST,MXLNG,IOCOC,ISMOST,
C    &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
C    &                  ISIMSYM)

 
      END IF
      CALL MEMCHK2('Z_AFPA')
*. Number of BLOCKS
      IF(I_DO_BLOCKS.EQ.1) THEN
        NBLOCK = IB_FOR_BATCH(NBATCH)+NBLK_PER_BATCH(NBATCH)-1
*. Extract length of each block
        CALL EXTRROW(IBLOCKFO,8,8,NBLOCK,LEN_PER_BLK)
*
        IF(IPRNT.GT.10.and.ga_nodeid().eq.0) THEN 
           WRITE(6,'(A,I9)') ' Number of blocks ', NBLOCK
           WRITE(6,'(A,I9)') ' Number of batches ', NBATCH
        END IF
        IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
          WRITE(6,*) ' Number of blocks per batch: '
          CALL IWRTMA(NBLK_PER_BATCH,1,NBATCH,1,NBATCH)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE NCONF_OCCLS(NOBPSP,NELPSP,NSPC,MINOP,
     &           NCONF,LOCC)
*
* An occupation class is given. 
* Determine number of configurations for this class with atleast MINOP open orbitals
* (Realizing that Weyl did this in a smarter way, but that is no surprise...)
*
*. Jeppe Olsen, March 2013
*  Last revision, April 2013; Jeppe Olsen; LOCC added 
*
* The routine goes as 
* 1) calculate the number of configurations for each orbital space as function of 
*    number of unpaired electrons
* 2) Combine the results for each orbital space to obtain the number of confs as 
*    a function of number of open electrons
* 3) sum these
*. Step 1) and 2) could be merged.
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER NELPSP(NSPC),NOBPSP(NSPC)
*. local scratch
      INTEGER NCNPOPSP(MXPNGAS*(MXPNEL+1))
      INTEGER NCNPOP(MXPNEL+1), NCNPOP2(MXPNEL+1)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from NCONF_OCCLS '
        WRITE(6,*) ' ========================'
        WRITE(6,*) 
        WRITE(6,*) ' Number of orbitals per subspace: '
        CALL IWRTMA3(NOBPSP,1,NSPC,1,NSPC)
        WRITE(6,*) ' Number of electrons per subspace '
        CALL IWRTMA3(NELPSP,1,NSPC,1,NSPC)
      END IF
*
      NELTOT = IELSUM(NELPSP,NSPC)
      NOBTOT = IELSUM(NOBPSP,NSPC)
*. Max for number of unpaired electrons
      MXOPGL = MIN(NELTOT,2*NOBTOT-NELTOT)
*
      IZERO = 0
      CALL ISETVC(NCNPOPSP,IZERO,(MXOPGL+1)*NSPC)
      
*. Generate for each orbitals-space, the number of configurations with given number of open electrons
      DO ISPC = 1, NSPC
        NEL = NELPSP(ISPC)
        NOB = NOBPSP(ISPC)
*. Largest number of singly occupied orbitals
        MXOP = MIN(NEL,2*NOB-NEL)
        DO IOP = 0, MXOP
          IF(MOD(NEL-IOP,2).EQ.0) THEN
            ICL = (NEL - IOP)/2
*. Number of ways of obtaining IOP/ICL open/occupied orbitals
            NCONF = IBION(NOB,IOP)*IBION(NOB-IOP,ICL)
            NCNPOPSP((ISPC-1)*(MXOPGL+1)+IOP+1) = NCONF
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) 
     &  ' Number of conf for given number of open orbitals and orbspace'
        CALL IWRTMA(NCNPOPSP,MXOPGL+1,NSPC,MXOPGL+1,NSPC)
      END IF
*
*. Combine the results for the various orbital spaces
*
      IF(NSPC.EQ.1) THEN
        CALL ICOPVE(NCNPOPSP,NCNPOP2,MXOPGL+1)
      END IF
        
      DO ISPC = 2, NSPC
        IF(ISPC.EQ.2) THEN
          CALL ICOPVE(NCNPOPSP,NCNPOP,MXOPGL+1)
        ELSE
          CALL ICOPVE(NCNPOP2,NCNPOP,MXOPGL+1)
        END IF
        CALL ISETVC(NCNPOP2,IZERO,MXOPGL+1)
        DO IOPCM = 0, MXOPGL
         DO IOPTO = IOPCM, MXOPGL
           IOPSP = IOPTO - IOPCM
           NCNPOP2(IOPTO+1) = NCNPOP2(IOPTO+1) +
     &     NCNPOP(IOPCM+1)*NCNPOPSP((ISPC-1)*(MXOPGL+1)+IOPSP+1)
         END DO
        END DO
      END DO
*
*. Length of occupations
*
      LOCC = 0
      DO IOPCM = 0, MXOPGL
        IF(MOD(NELTOT-IOPCM,2).EQ.0) THEN
          IOCCM = IOPCM + (NELTOT-IOPCM)/2
          LOCC = LOCC + IOCCM*NCNPOP2(IOPCM+1)
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) 
     & ' Number of configurations for various numbers of unpaired elec.'
       WRITE(6,*) 
     & ' =============================================================='
       CALL IWRTMA(NCNPOP2,MXOPGL+1,1,MXOPGL+1)
      END IF
*
*. And sum the configuration count
*
      NCONF = 0
      DO IOP = MINOP, MXOPGL
       NCONF = NCONF + NCNPOP2(IOP+1)
      END DO
*
      IF(NTEST.GE.10) THEN
       WRITE(6,*) ' Number of configurations ', NCONF
       WRITE(6,*) ' Length of occ list ',       LOCC
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' Leaving NCONF_OCCLS '
      RETURN
      END 
      SUBROUTINE OCCLS_TO_OCSBCLS(IOCCLS,NOCCLS,NOBSPC,
     &           NOCSBCLS,MINOCFSPC,MAXOCSPC,NOCFSPC,
     &           IFLAG)
*
* A set of occupation classes IOCCLS are given. Obtain
* the corresponding occupations of the occupation subclasses
*
* A occupation sub class is a given occupation in a given 
* orbital space
*
* if IFLAG = 1, then only the number of occupation sub classes 
* are determined
*
* 
*
*. Jeppe Olsen, March 20, 2013
*
      INCLUDE 'implicit.inc'
*. input
      INTEGER IOCCLS(NOBSPC, NOCCLS)
*. Output
      INTEGER MINOCFSPC(NOBSPC),MAXOCFSPC(NOBSPC)
      INTEGER NOCFSPC(NOBSPC)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================'
        WRITE(6,*) ' OCCLS_TO_OCSBCLS '
        WRITE(6,*) ' ================'
        WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The occupation classes (as columns) '
        CALL IWRTMA3(IOCCLS,NOBSPC,NOCCLS,NOBSPC,NOCCLS)
      END IF
*
      NOCSBCLS = 0
      DO IOBSPC = 1, NOBSPC
        MAXO = IOCCLS(IOBSPC,1)
        MINO = IOCCLS(IOBSPC,1)
        DO ICLS = 2, NOCCLS
         IF(IOCCLS(IOBSPC,ICLS).GT.MAXO) 
     &   MAXO = IOCCLS(IOBSPC,ICLS)
         IF(IOCCLS(IOBSPC,ICLS).LT.MINO) 
     &   MINO = IOCCLS(IOBSPC,ICLS)
        END DO
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IOBSPC, MINO, MAXO = ',
     &                 IOBSPC, MINO, MAXO
        END IF
        IF(IFLAG.NE.1) THEN
          MINOCFSPC(IOBSPC) = MINO
          MAXOCFSPC(IOBSPC) = MAXO
          NOCFSPC(IOBSPC) = MAXO- MINO
        END IF
*
        NOCSBCLS = NOCSBCLS+(MAXO-MINO+1)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info for occpation of individual orbital spaces '
        WRITE(6,*) ' ================================================'
        WRITE(6,*)
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' Orb. space  Min Occ   Max occ '
          WRITE(6,*) ' =============================='
          DO IOBSPC = 1, NOBSPC
            WRITE(6,'(3(3X,I5))')
     &      IOBSPC, MINOCFSPC(IOBSPC),MAXOCFSPC(IOBSPC)
          END DO
        END IF
        WRITE(6,*)' Total number of occupation sub classes ',
     &  NOCSBCLS
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_SUBCNF(NEL,NOB,MINOP,NSBCNF,LOCSBCNF)
*
* A class of subconfigurations is defined by NEL electrons and NOB 
* orbitals with atleast MINOP open orbitals. Find number of 
* subconfigurations and length of occupation list
*
*. Jeppe Olsen, March 20, 2013
*
      INCLUDE 'implicit.inc'
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from DIM_SUBCNF = '
        WRITE(6,*) ' NEL, NOB, MINOP = ', NEL, NOB,MINOP
      END IF
*
      MAXOP = MIN(NEL,NOB-NEL)
      NSBCNF = 0
      LOCSBCNF = 0
      DO IOP = MINOP, MAXOP
        IF(MOD(NEL-IOP,2).EQ.0) THEN
          ICL = (NEL-IOP)/2
          NSBCNF = NSBCNF + IBION(NOB,IOP)*IBION(NOB-IOP,ICL)
          LOCSBCNF = LOCSBCNF 
     &  + IBION(NOB,IOP)*IBION(NOB-IOP,ICL)*(ICL+IOP)
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of subconfigurations ', NSBCNF
        WRITE(6,*) ' Length of occupation list ', LOCSBCNF 
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_SUBCONF(NEL,NOB,ISMFOB,MINOP,NSBCNF_FOR_OP_SM,
     &           NSMOB,
     &           IBSBCNF_FOR_OP_SM,LSBCNF,IOCSBCNF,IFLAG,MAXOP_TOT,
     &           IBORB)
*
* Generate the subconfigurations defined by NEL electron in 
* NOB orbitals for each symmetry and number of occupied orbitals
*
* IFLAG = 0: Construct NSBCNF_FOR_OP_SM, IBSBCNF_FOR_OP_SM,LSBCNF  and not the actual occupations
* IFLAG = 1: Construct also IOCSBCNF
*
*. Jeppe Olsen, March 21, 2013
*  Last revision; April 26 2013; Jeppe Olsen; IBORB added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*
*. Input (IFLAG = 1), output (IFLAG = 0)
      INTEGER NSBCNF_FOR_OP_SM(MAXOP_TOT+1,NSMOB)
*. IBSBCNF_FOR_OP will be offset for the OCCUPATION with given
*  nopen and sym.
      INTEGER IBSBCNF_FOR_OP_SM(MAXOP_TOT+1,NSMOB)
*. Symmetry of orbitals
      INTEGER ISMFOB(NOB)
*. Output
      INTEGER IOCSBCNF(*)
*. Local scratch
      INTEGER JOP(MXPNEL),JJCL(MXPNEL),JCL(MXPNEL),JOB(MXPNEL)
*
      IZERO = 0
      CALL ISETVC(NSBCNF_FOR_OP_SM,IZERO,(MAXOP_TOT+1)*NSMOB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GEN_SUBCONF '
        WRITE(6,*) ' ====================== '
        WRITE(6,*) 
        WRITE(6,*) ' NEL, NOB = ', NEL, NOB
        WRITE(6,*) ' IFLAG = ', IFLAG
        WRITE(6,*) ' MINOP, MAXOP_TOT = ', MINOP, MAXOP_TOT
      END IF
*
      MAXOP = MIN(NEL,2*NOB-NEL)
 
      DO IOP = MINOP, MAXOP
         IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' Info for IOP = ', IOP
         END IF
         IF(MOD(NEL-IOP,2).EQ.0) THEN
           ICL = (NEL-IOP)/2
*. Loop over the various ways of distributing IOP electrons in NOB orbitals
           INI_OP = 1
 1001      CONTINUE
           IF(INI_OP.EQ.1) THEN
C                 ISTVC2(IVEC,IBASE,IFACT,NDIM)
             CALL ISTVC2(JOP,0,1,IOP)
             NONEW_OP = 0
             INI_OP = 0
           ELSE
C            NXTORD(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
             CALL NXTORD(JOP,IOP,1,NOB,NONEW_OP)
           END IF
           IF(NTEST.GE.10000) THEN
             WRITE(6,*) ' Next open orb config '
             CALL IWRTMA3(JOP,1,IOP,1,IOP)
           END IF
*
           IF(NONEW_OP.EQ.0)  THEN
*. A new string of unpaired electrons have been obtained, generate various closed shells
*. Symmetry of open string
C                         ISYMSTB(ISTR,ISYM,NEL)
             IOPSTR_SYM = ISYMSTB(JOP,ISMFOB,IOP)
             IF(NTEST.GE.10000) WRITE(6,*) ' IOPSTR_SYM = ', IOPSTR_SYM
             IF(IFLAG.EQ.0) THEN
*. Number of closed orbital contributions
               NCL_TERM = IBION(NOB-IOP,ICL)
               IF(NTEST.GE.10000) WRITE(6,*) ' NCL_TERM = ', NCL_TERM
               MCNF = NCL_TERM
               NSBCNF_FOR_OP_SM(IOP+1,IOPSTR_SYM) =
     &         NSBCNF_FOR_OP_SM(IOP+1,IOPSTR_SYM) + MCNF
             ELSE
*. Determine the various closed orbitals terms
*. The orbitals not single occupied
C                   INT_NOT_IN_STRING(ISTRING,NEL,NOB,I0STRING)
               CALL INT_NOT_IN_STRING(JOP,IOP,NOB,JOB)
               INI_CL = 1
  901          CONTINUE
                 IF(INI_CL.EQ.1) THEN
                   CALL ISTVC2(JJCL,0,1,ICL)
                   INI_CL = 0
                   NONEW_CL = 0
                 ELSE
                   CALL NXTORD(JJCL,ICL,1,NOB-IOP,NONEW_CL)
                 END IF
                 IF(NTEST.GE.10000) THEN
                   WRITE(6,*) ' Next closed config as JJCL '
                   CALL IWRTMA3(JJCL,1,ICL,1,ICL)
                 END IF
                 IF(NONEW_CL.EQ.0) THEN
*. We now have a new subconfiguration
                   DO KCL = 1, ICL
                     JCL(KCL) = JOB(JJCL(KCL))
                   END DO
                 IF(NTEST.GE.10000) THEN
                   WRITE(6,*) ' Next closed config as JCL '
                   CALL IWRTMA3(JCL,1,ICL,1,ICL)
                 END IF
*. Merge the closed and open orbital parts of config and save 
                   NSBCNF_FOR_OP_SM(IOP+1,IOPSTR_SYM) =
     &             NSBCNF_FOR_OP_SM(IOP+1,IOPSTR_SYM) + 1
                   NN = NSBCNF_FOR_OP_SM(IOP+1,IOPSTR_SYM)
                   IOC = IOP + ICL
                   IB = IBSBCNF_FOR_OP_SM(IOP+1,IOPSTR_SYM) + (NN-1)*IOC
                   IF(NTEST.GE.10000) 
     &             WRITE(6,*) ' Conf will be copied to IB = ', IB
C                  MERGE_CLOP_CONF(IOP,NOP,ICL,NCL,ICONF)
                   CALL MERGE_CLOP_CONF(JOP,IOP,JCL,ICL,
     &                  IOCSBCNF(IB),IBORB)
               GOTO 901 ! Take another turn
                 END IF ! NONEW_CL = 0
             END IF ! IFLAG switched
          GOTO 1001 ! Next open part
           END IF ! NONEW_OP = 0
        END IF ! MOD test passed
      END DO! Loop over IOP
*
      IF(IFLAG.EQ.0) THEN
*. Offset to the occupations of the subconfigurations
        LOC = 0
        IBSBCNF_FOR_OP_SM(1,1) = 1
        DO ISM = 1, NSMOB
          DO IOP = 0, MAXOP_TOT
            IOC = IOP + (NEL-IOP)/2
            NSBCNF = NSBCNF_FOR_OP_SM(IOP+1,ISM)
            LOC = IOC*NSBCNF
            IF(IOP.LT.MAXOP_TOT) THEN
             IBSBCNF_FOR_OP_SM(IOP+2,ISM) = 
     &       IBSBCNF_FOR_OP_SM(IOP+1,ISM) + LOC
            ELSE IF (ISM.NE.NSMOB) THEN
             IBSBCNF_FOR_OP_SM(1,ISM+1) = 
     &       IBSBCNF_FOR_OP_SM(IOP+1,ISM) + LOC
            ELSE
*. IOP = MXOP_TOT, ISM = NSMOB, we are at the end so
             LSBCNF = IBSBCNF_FOR_OP_SM(IOP+1,ISM) + LOC - 1
            END IF
          END DO
        END DO
      END IF ! IFLAG = 0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Number of subconfs for each number of open orbs(r) and sym(c)'
         CALL IWRTMA(NSBCNF_FOR_OP_SM,
     &        MAXOP_TOT+1,NSMOB,MAXOP_TOT+1,NSMOB)
        WRITE(6,*) 
        IF(IFLAG.EQ.0) 
     &  WRITE(6,*) ' Length of occupation list for subconfigurations ',
     &  LSBCNF
       END IF
*
       IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The corresponding offset array '
         CALL IWRTMA(IBSBCNF_FOR_OP_SM,
     &        MAXOP_TOT+1,NSMOB,MAXOP_TOT+1,NSMOB)
       END IF
*
       IF(IFLAG.EQ.1.AND.NTEST.GE.1000) THEN
         WRITE(6,*) ' Occupations of subconfigurations '
         DO ISM = 1, NSMOB
          DO IOP = MINOP, MAXOP
           IB = IBSBCNF_FOR_OP_SM(IOP+1,ISM)
           N  = NSBCNF_FOR_OP_SM(IOP+1,ISM)
           IOC = IOP + (NEL-IOP)/2
           IF(N.NE.0) 
     &     WRITE(6,'(A,2I4)') ' Subconfs with nopen and sym = ', IOP,ISM
           DO ISBCNF = 1, N
             CALL WRT_CONF(IOCSBCNF(IB+(ISBCNF-1)*IOC),IOC)
           END DO
          END DO
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE INT_NOT_IN_STRING(ISTRING,NEL,NOB,I0STRING)
*
* A string of NEL integers is given in ISTRING.
* Determine the integers not in ISTRING
*
*. Jeppe Olsen, March 20, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISTRING(NEL)
*. Output
      INTEGER I0STRING(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from INT_NOT_IN_STRING '
        WRITE(6,*) ' NEL, NOB = ', NEL, NOB
        WRITE(6,*) ' Input string: '
        CALL IWRTMA3(ISTRING,1,NEL,1,NEL)
      END IF
*
      IEL = 1
      I0EL = 1
      DO IOB = 1, NOB
        IF(IEL.GT.NEL) THEN
          I0STRING(I0EL) = IOB
          I0EL = I0EL + 1
        ELSE 
          IF(IOB.LT.ISTRING(IEL)) THEN
            I0STRING(I0EL) = IOB
            I0EL = I0EL + 1
          ELSE IF (IOB.EQ.ISTRING(IEL)) THEN
            IEL = IEL + 1
          END IF
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from INT_NOT_IN_STRING '
        WRITE(6,*) ' Input string: '
        CALL IWRTMA3(ISTRING,1,NEL,1,NEL)
        WRITE(6,*) ' Integers not in string '
        CALL IWRTMA3(I0STRING,NOB-NEL,1,NOB-NEL,1)
      END IF
*
      RETURN
      END
      FUNCTION ISYMSTB(ISTR,ISYM,NEL)
*
* Determined symmetry of string ISTR when the symmetry of each orbital is given by ISYM
*
*. Jeppe Olsen, March 20, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*. Specific input
      INTEGER ISTR(NEL)
*. general input
      INTEGER ISYM(*)
*
      IISYM = 1
      DO JEL = 1, NEL
        IISYM = MULTD2H(IISYM,ISYM(ISTR(JEL)))
      END DO
*
      ISYMSTB = IISYM
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ISYMSTB: symmetry of string ', IISYM
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Corresponding string '
        CALL IWRTMA3(ISTR,1,NEL,1,NEL)
      END IF
*
      RETURN
      END
      SUBROUTINE MERGE_CLOP_CONF(IOP,NOP,ICL,NCL,ICONF,IB)
*
* Open and closed orbital parts of a configuration is given in IOP, ICL
* Merge to obtain a configuration in compact form with doubly 
* occupied orbitals given by a -
* Offset to orbitals is IB, so IB-1 is added to orbital numbers
*
*. Jeppe Olsen, March 20, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOP(NOP), ICL(NCL)
*. Output
      INTEGER ICONF(NOP+NCL)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MERGE_CLOP_CONF in action '
        WRITE(6,*) ' =========================='
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input IOP and ICL: '
        CALL IWRTMA3(IOP,1,NOP,1,NOP)
        WRITE(6,*)
        CALL IWRTMA3(ICL,1,NCL,1,NCL)
      END IF
*
      IIOP = 1
      IICL = 1
      DO IEL = 1, NOP + NCL
        IF(IIOP.LE.NOP.AND.IICL.LE.NCL) THEN
          IF(IOP(IIOP).LT.ICL(IICL)) THEN
            ICONF(IEL) = IOP(IIOP) + IB -1
            IIOP = IIOP + 1
           ELSE 
            ICONF(IEL) = -(ICL(IICL) + IB - 1)
            IICL = IICL + 1
           END IF
        ELSE IF (IIOP.LE.NOP) THEN
           ICONF(IEL) = IOP(IIOP) + IB - 1
           IIOP = IIOP + 1
        ELSE IF (IICL.LE.NCL) THEN
            ICONF(IEL) = -(ICL(IICL) + IB - 1)
            IICL = IICL + 1
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
       NOC = NCL + NOP
       WRITE(6,*) ' Resulting configuration '
       CALL IWRTMA3(ICONF,1,NOC,1,NOC)
      END IF
*
      RETURN
      END
      SUBROUTINE INFO_OCSBCLS
*
* Generate the occupation subsclasses 
*
*. Jeppe Olsen, March 22, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' Info from INFO_OCSBCLS '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
      END IF
*
* Number and offset for occupation subclasses for given GASpace
*
      NSBCLS = 0
      DO IGAS = 1, NGAS
        IBOCSBCLS(IGAS) = NSBCLS + 1
        NOCSBSLS(IGAS) = MXGSOC(IGAS) - MNGSOC(IGAS) + 1
        NSBCLS = NSBCLS + NOCSBSLS(IGAS)
      END DO
      if (ga_nodeid().eq.0.and.NTEST.GE.100) WRITE(6,*) ' NSBCLS == ', NSBCLS
*
*. The actual occupation subclasses
*
      CALL GEN_OCSBCLS(MNGSOC, MXGSOC, NGAS, int_mb(KOGOCSBCLS))
C     GEN_OCSBCLS(MNGSOC, MXGSOC, NGAS, IOGSBCLS)
*
*. The occupation sub classes of each occupation class
      CALL OCSBCLS_OF_OCCLS(int_mb(KIOCCLS),NOCCLS_MAX, NGAS,
     &     IBOCSBCLS,MNGSOC,int_mb(KOCSBCLS_OF_OCCLS),
     &     int_mb(KMINOPGAS_FOR_OCCLS),MINOP,NOBPT)
C          OCSBCLS_OF_OCCLS(IOCCLS,NOCCLS,NGAS,
C    &           IBOCSBCLS,MNGSOC,IOCSBCLS_OF_OCCLS,
C    &           MINOPGAS_FOR_OCCLS,MINOP,NOBPT)
*
* The minimum number of open electrons for each occupation 
* sub classes
*
      CALL MINOP_FOR_OCSBCLS(int_mb(KMNOPOCSBCL),NGAS,NOCCLS_MAX,
     &     NOCSBCLST,int_mb(KMINOPGAS_FOR_OCCLS),
     &     int_mb(KOCSBCLS_OF_OCCLS)) 
C     MINOP_FOR_OCSBCLS(MNOPOCSBCLS,NGAS,NOCCLS,
C    &           NOCSBCLST,MINOPGAS_FOR_OCCLS,IOCSBCLS_OF_OCCLS)
*
      RETURN
      END
      SUBROUTINE OCSBCLS_OF_OCCLS(IOCCLS,NOCCLS,NGAS,
     &           IBOCSBCLS,MNGSOC,IOCSBCLS_OF_OCCLS,
     &           MINOPGAS_FOR_OCCLS,MINOP,NOBPT)
*
* A set of occupation classes are given. Find the occupation sub classes
* of each occupation class and the minumum number of open shells in each
* occupation sub class for a given occupation class
*
*. Jeppe Olsen, March 22, 2012
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOCCLS(NGAS,NOCCLS)
      INTEGER IBOCSBCLS(NGAS),MNGSOC(NGAS)
      INTEGER NOBPT(NGAS)
*. Output
      INTEGER IOCSBCLS_OF_OCCLS(NGAS,NOCCLS)
      INTEGER MINOPGAS_FOR_OCCLS(NGAS,NOCCLS)
*
      NTEST = 000
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Info from OCSBCLS_OF_OCCLS '
        WRITE(6,*) ' ========================== '
      END IF
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Input occupation classes '
        DO ICLS = 1, NOCCLS
          WRITE(6,*) ICLS, '     ', (IOCCLS(IGAS,ICLS),IGAS = 1, NGAS)
        END DO
        WRITE(6,*) ' MINOP = ', MINOP
      END IF
*
      DO ICLS = 1, NOCCLS
       DO IGAS = 1, NGAS
         ISUB = IBOCSBCLS(IGAS)-1+IOCCLS(IGAS,ICLS)-MNGSOC(IGAS)+1
C?       WRITE(6,*) ' IBOCSBCLS(IGAS), MNGSOC(IGAS), = ',
C?   &                IBOCSBCLS(IGAS), MNGSOC(IGAS)
         IOCSBCLS_OF_OCCLS(IGAS,ICLS) = ISUB
       END DO
      END DO
*
      DO ICLS = 1, NOCCLS
*. Max number of open shells in all GASpaces
         IF(NTEST.GE.1000.and.ga_nodeid().eq.0) 
     &      WRITE(6,*) ' ICLS = ', ICLS
         MAXOP_L= 0
         DO IGAS = 1, NGAS
           NELG = IOCCLS(IGAS,ICLS)
           NOBG = NOBPT(IGAS)
           MAXOP_L = MAXOP_L + MIN(NELG,2*NOBG-NELG)
         END DO
         IF(NTEST.GE.1000.and.ga_nodeid().eq.0)
     &   WRITE(6,*) ' ICLS, MAXOP_L = ', ICLS, MAXOP_L
         DO IGAS = 1, NGAS
*. Max. number of open orbitals in spaces different from IGAS
           NELL = IOCCLS(IGAS,ICLS)
           NOBL = NOBPT(IGAS)
           MAXOP_LO = MAXOP_L - MIN(NELL,2*NOBL-NELL)
           IF(NTEST.GE.1000.and.ga_nodeid().eq.0) 
     &     WRITE(6,*) ' IGAS, MAXOP_LO = ', IGAS,MAXOP_LO
           MINOPGAS_FOR_OCCLS(IGAS,ICLS) = MAX(0,MINOP-MAXOP_LO)
         END DO
      END DO
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) 
     &  ' Occupation sub classes for the various occupation classes '
        WRITE(6,*) 
     &  ' ========================================================= '
        DO ICLS = 1, NOCCLS
          WRITE(6,'(A, I4,A, 30(1X,I3))')
     &    ' Occupation class ', ICLS, ' has subclasses ', 
     &    (IOCSBCLS_OF_OCCLS(IGAS,ICLS),IGAS=1, NGAS)
        END DO
*
        WRITE(6,*) ' Min number of open orbitals for each GAS space'
        WRITE(6,*) ' =============================================='
        WRITE(6,*)
        WRITE(6,*) ' Occ. class   Min open per gas space '
        WRITE(6,*) ' ===================================='
        DO ICLS = 1, NOCCLS
         WRITE(6,'(3X,I5,6X,20(1X,I3))') 
     &   ICLS, (MINOPGAS_FOR_OCCLS(IGAS,ICLS),IGAS = 1, NGAS)
        END DO
*
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_OCSBCLS(MNGSOC, MXGSOC, NGAS, IOGSBCLS)
*
* Occupation and Gaspace of the various occupation sub classes
*
*. Jeppe Olsen, March 22, 2013
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER MNGSOC(NGAS), MXGSOC(NGAS)
*. Output
      INTEGER IOGSBCLS(2,*)
*
      NOCSBCLS = 0 
      DO IGAS = 1, NGAS
        DO IOCC = MNGSOC(IGAS), MXGSOC(IGAS)
          NOCSBCLS = NOCSBCLS +1
          IOGSBCLS(1,NOCSBCLS) = IOCC
          IOGSBCLS(2,NOCSBCLS) = IGAS
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' The generated occupation sub classes '
        WRITE(6,*) ' ====================================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation sub class    Occupation   Gaspace '
        WRITE(6,*) ' ============================================='
        DO MOCSBCLS = 1, NOCSBCLS
          WRITE(6,'(6X, I5, 16X, I3, 11X, I2)')
     &    MOCSBCLS, IOGSBCLS(1,MOCSBCLS), IOGSBCLS(2,MOCSBCLS)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE MINOP_FOR_OCSBCLS(MNOPOCSBCLS,NGAS,NOCCLS,
     &           NOCSBCLST,MINOPGAS_FOR_OCCLS,IOCSBCLS_OF_OCCLS)
*
* Obtain the minimum number of open orbitals per occupation sub class
*
*. Jeppe Olsen, March 2013
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER MINOPGAS_FOR_OCCLS(NGAS,NOCCLS)
      INTEGER IOCSBCLS_OF_OCCLS(NGAS,NOCCLS)
*. Output
      INTEGER MNOPOCSBCLS(NOCSBCLST)
*
      NTEST = 000
*
      MONE = -1
      CALL ISETVC(MNOPOCSBCLS,MONE,NOCSBCLST)
*
      DO IOCCLS = 1, NOCCLS
       DO IGAS = 1, NGAS
        IOCSBCLS = IOCSBCLS_OF_OCCLS(IGAS,IOCCLS)
        IF(MNOPOCSBCLS(IOCSBCLS).EQ.MONE) THEN
          MNOPOCSBCLS(IOCSBCLS) = 
     &    MINOPGAS_FOR_OCCLS(IGAS,IOCCLS)
        ELSE
          MNOPOCSBCLS(IOCSBCLS) = 
     &    MIN(MINOPGAS_FOR_OCCLS(IGAS,IOCCLS),MNOPOCSBCLS(IOCSBCLS))
        END IF
       END DO
      END DO
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Minimum number of open orbitals per oc sub class '
        WRITE(6,*) ' ================================================= '
        CALL IWRTMA(MNOPOCSBCLS,1,NOCSBCLST,1,NOCSBCLST)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_DIM_SBCNF(NSBCNF_FOR_OP_SM,IBSBCNF_FOR_OP_SM,
     &           LSBCNF,IOGOCSBCLS,MINOPFSBCLS)
*
* Generate info on the number of subconfigurations for each occupation sub class
*
* NSBCNF_FOR_OP_SM: Number of subconfs with given iopen, sym and ocsbcls
* IBSBCNF_FOR_OP_SM: Offset to occupations of confs with given iopen, sym and type
* LSBCNF: Length of occupations arrays for subconfigurations
*
*. Jeppe Olsen, March 25, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'orbinp.inc'
c      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOGOCSBCLS(2,NOCSBCLST)
      INTEGER MINOPFSBCLS(NOCSBCLST)
*. Output
      INTEGER NSBCNF_FOR_OP_SM(NOPEN_MAX+1,NSMOB,NOCSBCLST)
      INTEGER IBSBCNF_FOR_OP_SM(NOPEN_MAX+1,NSMOB,NOCSBCLST)
      INTEGER LSBCNF(NOCSBCLST)
*
      NTEST = 000
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Entering GEN_DIM_SBCNF '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
      END IF
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' MINOPFSBCLS array: '
        CALL IWRTMA3(MINOPFSBCLS,1,NOCSBCLST,1,NOCSBCLST)
      END IF
*
      DO IOCSBCLS = 1, NOCSBCLST
         NEL = IOGOCSBCLS(1, IOCSBCLS)
         IGAS = IOGOCSBCLS(2, IOCSBCLS)
         NOB = NOBPT(IGAS)
         MINOPL = MINOPFSBCLS(IOCSBCLS)
         IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
         WRITE(6,'(A,5I4)') ' IOCSBCLS, NEL, IGAS, NOB, MINOPL = ',
     &                IOCSBCLS, NEL, IGAS, NOB, MINOPL
         END IF
C?       WRITE(6,*) ' MINOPFSBCLS(1) = ', MINOPFSBCLS(1)
         IBOB = 1
         DO JGAS = 0, IGAS -1 
          DO JSM = 1, NSMOB
           IBOB = IBOB + NOBPTS_GN(JGAS,JSM)
          END DO
         END DO
C?       WRITE(6,*) ' IBOB = ' , IBOB
         IDUM = 0
         CALL GEN_SUBCONF(NEL,NOB,ISMFTO(IBOB),MINOPL,
     &        NSBCNF_FOR_OP_SM(1,1,IOCSBCLS),NSMOB,
     &        IBSBCNF_FOR_OP_SM(1,1,IOCSBCLS),LSBCNF(IOCSBCLS),
     &        IDUM,0,NOPEN_MAX,IBOB)
      END DO 
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' The LSBCNF array: length of occ per ocsblcs '
        CALL IWRTMA(LSBCNF,1,NOCSBCLST,1,NOCSBCLST)
*
        DO IOCSBCLS = 1, NOCSBCLST
          WRITE(6,*)
          WRITE(6,'(A,I4)') ' Info for occupation sub class', IOCSBCLS
          WRITE(6,'(A)')    ' =================================='
          WRITE(6,*)
          WRITE(6,*) ' Number of subconfs per open and sym:'
          CALL IWRTMA(NSBCNF_FOR_OP_SM(1,1,IOCSBCLS),
     &         NOPEN_MAX+1,NSMOB,NOPEN_MAX+1,NSMOB)
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE OCCLSDIM_FROM_OCSBCLSDIM(IOCCLS,IOCCLSDIM,
     &           IOPEN_OUT_DIM)
*
* Obtain the dimension of occupation class IOCCLS
*
*. Jeppe Olsen, Mar. 26, 2013
*
*. General input
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input: Occupation of the occupation class
      INTEGER IOCCLS(NGAS)
*. Output
      INTEGER IOCCLSDIM(IOPEN_OUT_DIM,NSMOB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Entering OCCLSDIM_FROM_OCSBCLSDIM '
C?      WRITE(6,*) ' IOPEN_OUT_DIM = ', IOPEN_OUT_DIM
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM, 'OCSBDM')
*. Two local arrays holding symmetry and number of open orbitals
      LDIM = (IOPEN_OUT_DIM+1)*NSMOB
      CALL MEMMAN(KLOS1,LDIM,'ADDL  ',1,'LO_S1 ')  !done
      CALL MEMMAN(KLOCSBCLS,NGAS,'ADDL  ',1,'OCSBCL')  !done
*. The occupation sub classes corresponding to the occupation class
      CALL OCC_TO_SBCLS_FOR_OCCLS(IOCCLS,int_mb(KLOCSBCLS))
*. 
      CALL OCCLSDIM_FROM_OCSBCLSDIM_IN(int_mb(KLOCSBCLS),
     &     int_mb(KNSBCNF),
     &     int_mb(KLOS1),IOCCLSDIM,IOPEN_OUT_DIM)

      CALL MEMMAN(IDUM,IDUM,'FLUSM ', IDUM, 'OCSBDM')
      RETURN
      END
      SUBROUTINE OCCLSDIM_FROM_OCSBCLSDIM_IN(IOCSBCLS,NSBCNF, ILOS1,
     &     IOCCLSDIM,IOPEN_OUT_DIM)
*
*. Dimension of occupation class from dimensions of occupation sub classes'
*
*. Jeppe Olsen, March 26, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER IOCSBCLS(NGAS)
      INTEGER NSBCNF(NOPEN_MAX+1,NSMOB,*)
*. Output
      INTEGER IOCCLSDIM(IOPEN_OUT_DIM+1,NSMOB)
*. Scratch
      INTEGER ILOS1(IOPEN_OUT_DIM+1,NSMOB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from OCCLSDIM_FROM_OCSBCLSDIM_IN'
        WRITE(6,*) ' ===================================='
        WRITE(6,*)
        WRITE(6,*) ' Input occupation sub classes: '
        CALL IWRTMA(IOCSBCLS,1,NGAS,1,NGAS)
      END IF
*
      IZERO = 0
      CALL ISETVC(ILOS1,IZERO, (IOPEN_OUT_DIM+1)*NSMOB)
*
      IF(NGAS.EQ.1) THEN
*. Just copy 
        J1OCSBCLS = IOCSBCLS(1)
        DO ISM = 1, NSMOB
          CALL ICOPVE(NSBCNF(1,ISM,J1OCSBCLS),IOCCLSDIM(1,ISM),
     &           IOPEN_OUT_DIM+1)
          END DO
      END IF
*
      DO IGAS = 2, NGAS
        JOCSBCLS = IOCSBCLS(IGAS)
*. Dimensions for GAS = 1, 2... IGAS in ILOS1
        IF(IGAS.EQ.2) THEN
          J1OCSBCLS = IOCSBCLS(1)
          DO ISM = 1, NSMOB
            CALL ICOPVE(NSBCNF(1,ISM,J1OCSBCLS),ILOS1(1,ISM),
     &           NOPEN_MAX+1)
          END DO
        ELSE 
          CALL ICOPVE(IOCCLSDIM, ILOS1, (IOPEN_OUT_DIM+1)*NSMOB)
        END IF
*. Multiply with subconfigurations in GAS = IGAS and save in IOCCLSDIM
C            DIM_PROD_OCSBCNF(I12DIM,I1DIM,I2DIM)
        CALL DIM_PROD_OCSBCNF(IOCCLSDIM,ILOS1, NSBCNF(1,1,JOCSBCLS),
     &       NOPEN_MAX+1,IOPEN_OUT_DIM)
      END DO
*. Zero terms with less than MINOP open orbitals
      IF(MINOP.GT.0) THEN
        DO ISM = 1, NSMOB
         DO IOP = 0, MINOP-1
          IOCCLSDIM(IOP+1,ISM) = 0
C         IOCCLSDIM(IOPEN_OUT_DIM+1,NSMOB)
         END DO
        END DO
      END IF
         
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Number of confs for occupation class: row: nopen, col: sym'
        CALL IWRTMA(IOCCLSDIM,NOPEN_MAX+1,NSMOB,IOPEN_OUT_DIM+1,NSMOB)
      END IF
*
      RETURN
      END
      SUBROUTINE OCC_TO_SBCLS_FOR_OCCLS(IOCCLS, IOCSBCLS)
*
* Obtain occupation subclass for occupation class defined by IOCCLS
*
*. Jeppe Olsen, March 26, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
*. Input
      INTEGER IOCCLS(NGAS)
*. Output
      INTEGER IOCSBCLS(NGAS)
*
      NTEST = 00
*
      DO IGAS = 1, NGAS
        ISUB = IBOCSBCLS(IGAS)-1+IOCCLS(IGAS)-MNGSOC(IGAS)+1
        IOCSBCLS(IGAS) = ISUB
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from OCC_TO_SBCLS_FOR_OCCLS '
        WRITE(6,*) ' Input occupation class and output occsubclasses '
        WRITE(6,*)
        CALL IWRTMA3(IOCCLS,1,NGAS,1,NGAS)
        CALL IWRTMA3(IOCSBCLS,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_PROD_OCSBCNF(I12DIM,I1DIM,I2DIM,NROW2,
     &           IOPEN_OUT_DIM)
*
* Two sets of occupation sub configurations are given with dimension I1DIM, I2DIM, respectively
* Obtain dimension of product 
*
*. Jeppe Olsen, March 26, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'multd2h.inc'
*. Input
      INTEGER I1DIM(IOPEN_OUT_DIM+1,NSMOB),I2DIM(NROW2,NSMOB)
*. Output
      INTEGER I12DIM(IOPEN_OUT_DIM+1,NSMOB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from DIM_PROD_OCSBCNF '
        WRITE(6,*) ' ========================== '
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input I1DIM and I2DIM '
        CALL IWRTMA(I1DIM,NOPEN_MAX+1, NSMOB, IOPEN_OUT_DIM+1,NSMOB)
        WRITE(6,*)
        CALL IWRTMA(I2DIM,NOPEN_MAX+1, NSMOB, NROW2,NSMOB)
      END IF
*
      DO I12OP = 0, NOPEN_MAX
       DO I12SM = 1, NSMOB
         ISUM = 0
         DO I1OP = 0, I12OP
          DO I1SM = 1, NSMOB
            I2OP = I12OP - I1OP
            I2SM = MULTD2H(I12SM,I1SM)
            ISUM = ISUM + I1DIM(I1OP+1,I1SM)*I2DIM(I2OP+1,I2SM)
          END DO
         END DO
         I12DIM(I12OP+1,I12SM) = ISUM
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Dimensions for product: nopen(row), sym(col) '
        CALL IWRTMA(I12DIM,NOPEN_MAX+1,NSMOB,IOPEN_OUT_DIM+1,NSMOB)
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_CISPACE_FROM_SBCNF(ICISPC,NCNFOPSM)
*
* Obtain dimensions of CI space ICISPC from subconfiguration info
*
*  Info is stores in NCNFOPSM
*
*. Jeppe Olsen, March 26, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      INTEGER NCNFOPSM(MXPOPORB+1,MXPCSM)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from DIM_CISPACE_FROM_SBCNF '
        WRITE(6,*) ' ================================ '
        WRITE(6,*) 
        WRITE(6,*) ' Info for CI-space ', ICISPC
        WRITE(6,*) ' ============================ '
        WRITE(6,*)
      END IF 
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DMCISB')
*
      WRITE(6,*) ' NOCCLS_MAX = ', NOCCLS_MAX
      CALL MEMMAN(KLOCCLS_ACT,NOCCLS_MAX,'ADDL  ',1,'OCCLAC')  !done
      CALL MEMMAN(KLNCNFOPSM2,(MXPOPORB+1)*MXPCSM,'ADDL  ',1,'CNOPSM')!d
*
*. Find the occupation classes that are active for this CI space
C          OCCLS_IN_CI(NOCCLS,IOCCLS,ICISPC,NINCCLS,INCCLS)
      CALL OCCLS_IN_CI(NOCCLS_MAX,int_mb(KIOCCLS),ICISPC,NACTCLS,
     &     int_mb(KLOCCLS_ACT))
      CALL DIM_CISPACE_FROM_SBCNF_IN(int_mb(KLOCCLS_ACT),
     &     NOCCLS_MAX,
     &     int_mb(KIOCCLS),NCNFOPSM, int_mb(KLNCNFOPSM2),NGAS,
     &     int_mb(KNSBCNF),
     &     MINOP,NOPEN_MAX,MXPOPORB,NSMOB)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Number of confs in CI-space per iopen(row) and sym(col): '
        CALL IWRTMA10(NCNFOPSM,NOPEN_MAX+1,NSMOB,MXPOPORB+1,MXPCSM)
      END IF
*  
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DMCISB')
*
      RETURN
      END
      SUBROUTINE DIM_CISPACE_FROM_SBCNF_IN(IOCCLS_ACT,NOCCLS_MAX,IOCCLS,
     &           NCNFOPSM, NCNFOPSMS,NGAS,NSBCNF,MINOP,NOPEN_MAX,
     &           MXPOPORB, NSMOB)
*
* Obtain dimension for the CI space defined by IOCCLS_ACT, IOOCLS
* Inner routine
*
*. Jeppe Olsen, March 26, 2013
*
      INCLUDE 'implicit.inc'
*. input
      INTEGER IOCCLS(NGAS,NOCCLS_MAX), IOCCLS_ACT(NOCCLS_MAX)
      INTEGER NSBCNF(NOPEN_MAX+1,NSMOB,*)
*. Output
      INTEGER NCNFOPSM(MXPOPORB+1,NSMOB)
*. Scratch
      INTEGER NCNFOPSMS(MXPOPORB+1,NSMOB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from DIM_CISPACE_FROM_SBCNF_IN '
        WRITE(6,*) ' =================================== '
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Array giving active occupation classes '
        CALL IWRTMA3(IOCCLS_ACT,1,NOCCLS_MAX,1,NOCCLS_MAX)
      END IF
*
      IZERO = 0
      CALL ISETVC(NCNFOPSM,IZERO,(MXPOPORB+1)*NSMOB)
*
      DO JOCCLS = 1, NOCCLS_MAX
        IF(IOCCLS_ACT(JOCCLS).EQ.1) THEN
*. Occupation class is active in this CI-space
C              OCCLSDIM_FROM_OCSBCLSDIM(IOCCLS, IOCCLSDIM)
          CALL OCCLSDIM_FROM_OCSBCLSDIM(IOCCLS(1,JOCCLS),NCNFOPSMS,
     &         MXPOPORB)
          IDIM = (MXPOPORB+1)*NSMOB
          IONE = 1
          CALL IVCSUM(NCNFOPSM,NCNFOPSM,NCNFOPSMS,IONE,IONE,IDIM)
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
     &  ' Number of confs IN CI SPACE per iopen(row), sym(col) '
        WRITE(6,*)
     &  ' ==================================================== '
        CALL IWRTMA(NCNFOPSM,NOPEN_MAX+1,NSMOB,MXPOPORB+1,NSMOB)
      END IF
*
       RETURN
       END
       SUBROUTINE GEN_OCCONF_FOR_OCCLS_FROM_OCSBCLS(IOCCLS,ISM,IOCC,
     &            IBOCC_OP,IB_OP,NCONF_OP,MINOP,MAXOP,IBOCSBCNF,
     &            NSBCNF, IBSBCNF)
*
*. Generate occupations of configurations of occupation class IOCCLS
*. and add these to IOCC. The pointers are 
*. not reset, so info for several occupation classed may be stored 
*. consecutively- giving a CI space
*
* Input:
*   IOCCLS: Occupation class as occupations for each gaspace
*   ISM   : Symmetry of configurations
*   IBOCC_OP: start of occupations with given number of open orbitals: 
*             all occupation classes of CI space
*   IB_OP: is start of confs with given number of open orbitals: 
*         all occupation classes of CI space
*   NCONF_OP: Number of confs with given number of open orbitals  
*             for this occupation class
*   MAXOP: Largest number of open orbitals
*   IBOCSBCNF: Offset for occupation of the various subconfiguration lists
*
*
*. Jeppe Olsen, April3, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOCCLS(NGAS)
      INTEGER IBOCC_OP(*), IB_OP(*)
      INTEGER IBOCSBCNF(*)
      INTEGER NSBCNF(MAXOP+1,NSMOB,*)
      INTEGER IBSBCNF(MAXOP+1,NSMOB,*)
*. Input /Output
      INTEGER IOCC(*)
      INTEGER NCONF_OP(*)
*. Local scratch
      INTEGER ISBCLS(MXPNGAS), MINOPGASL(MXPNGAS)
      INTEGER IBOPSML((MXPOPORB+1)*NSMOB),NOPSML((MXPOPORB+1)*NSMOB)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CNFRSB')
      IZERO = 0
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from GEN_OCCONF_FOR_OCCLS_FROM_OCSBCLS'
        WRITE(6,*) ' ==========================================='
        WRITE(6,*)
        WRITE(6,*) ' The occupation class (as occupations )'
        CALL IWRTMA3(IOCCLS,1,NGAS,1,NGAS)
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' MINOP, MAXOP = ', MINOP, MAXOP
        WRITE(6,*) ' Initial IBOCC_OP,IB_OP '
        CALL IWRTMA(IBOCC_OP,MAXOP+1,1,MAXOP+1,1)
        CALL IWRTMA(IB_OP   ,MAXOP+1,1,MAXOP+1,1)
*
C?      WRITE(6,*) ' NSBCNF for the subconfclasses'
C?      DO JSBCLS = 1, NOCSBCLST
C?        WRITE(6,*) ' Subconfiguration class ', JSBCLS
C?        CALL IWRTMA(NSBCNF(1,1,JSBCLS),MAXOP+1,NSMOB,MAXOP+1,NSMOB)
C?      END DO
C?      WRITE(6,*) ' WORK(KLSBCNF) '
C?      CALL IWRTMA(WORK(KLSBCNF),1,NOCSBCLST,1,NOCSBCLST)
      END IF
*
*. The input array IBOCC_OP, NCONF_OP refers to specific symmetry,
*. copy to generel symmetry arrays
*
      CALL ISETVC(IBOPSML,IZERO,(MAXOP+1)*NSMOB)
      CALL ICOPVE(IBOCC_OP,IBOPSML(1+(MAXOP+1)*(ISM-1)),MAXOP+1)
      CALL ISETVC(NOPSML,IZERO,(MAXOP+1)*NSMOB)
      CALL ICOPVE(NCONF_OP,NOPSML(1+(MAXOP+1)*(ISM-1)),MAXOP+1)
*
*. Obtain the occupation subclasses for the occupation class
*
C     CALL OCSBCLS_OF_OCCLS(WORK(KIOCCLS),NOCCLS_MAX, NGAS,
C    &     IBOCSBCLS,MNGSOC,WORK(KOCSBCLS_OF_OCCLS),
C    &     WORK(KMINOPGAS_FOR_OCCLS),MINOP,NOBPT)
C          OCSBCLS_OF_OCCLS(IOCCLS,NOCCLS,NGAS,
C    &             IBOCSBCLS,MNGSOC,IOCSBCLS_OF_OCCLS,
C    &             MINOPGAS_FOR_OCCLS,MINOP,NOBPT)
      CALL OCSBCLS_OF_OCCLS(IOCCLS,1, NGAS,
     &     IBOCSBCLS,MNGSOC,ISBCLS,
     &     MINOPGASL,MINOP,NOBPT)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Occupation sub classes of occupation class: '
        CALL IWRTMA3(ISBCLS,1,NGAS,1,NGAS)
        WRITE(6,*) ' MIN op per gas '
        CALL IWRTMA3(MINOPGASL,1,NGAS,1,NGAS)
      END IF
      NEL = IELSUM(IOCCLS,NGAS)
*
      IF(NGAS.EQ.1) THEN
*. A single GAS space, just add configurations to output list. 
*. The output list may contain several occupation classes, so it is 
* not just copy
*. The routines in charge use arrays for all symmetries, so adopt
        ISBCLS1 = ISBCLS(1)
*
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Initial IBOPSML '
          CALL IWRTMA(IBOPSML,MAXOP+1,1,MAXOP+1,1)
        END IF
*
        CALL ADD_CNFS_TO_LIST(IOCC,NOPSML,IBOPSML,
     &       int_mb(IBOCSBCNF(ISBCLS1)),
     &       NSBCNF(1,1,ISBCLS1),IBSBCNF(1,1,ISBCLS1),
     &       MAXOP,ISM,NSMOB,NEL)
C     ADD_CNFS_TO_LIST(IOCCOUT,NCONFOUT,IBCONFOUT, 
C    &           IOCCADD, NCONFADD,IBCONFADD,MAXOP,ISM,NSMOB,NEL)
*. And copy arrays back to single symmetry arrays
        CALL ICOPVE(IBOPSML(1+(MAXOP+1)*(ISM-1)),IBOCC_OP,MAXOP+1)
        CALL ICOPVE(NOPSML(1+(MAXOP+1)*(ISM-1)),NCONF_OP,MAXOP+1)
      ELSE
*
*. We have several GAS spaces which configurations should be merged
*
*. The configurations are generated as product of subconfigurations
*. Determine size of largest expansion that is not complete, this is the NGAS -1 configs
C       NCONF_OCCLS(NOBPSP,NELPSP,NSPC,MINOP,NCONF,LOCC)
        CALL NCONF_OCCLS(NOBPT,IOCCLS,NGAS-1,MINOPGASL,NCONFM1,LOCCM1)
*. Space for two such lists of occupations 
        CALL MEMMAN(KLOCCSB1,LOCCM1,'ADDL  ',1,'OCCM11')  !done
        CALL MEMMAN(KLOCCSB2,LOCCM1,'ADDL  ',1,'OCCM12')  !done
*. And for two sets of dimension and offsets
        CALL MEMMAN(KLIBSB1,(MAXOP+1)*NSMOB,'ADDL  ',1,'IBSB1 ') !d
        CALL MEMMAN(KLIBSB2,(MAXOP+1)*NSMOB,'ADDL  ',1,'IBSB2 ') !d
        CALL MEMMAN(KLNSB1 ,(MAXOP+1)*NSMOB,'ADDL  ',1,'NSB1  ')  !d
        CALL MEMMAN(KLNSB2 ,(MAXOP+1)*NSMOB,'ADDL  ',1,'NSB2  ')  !d
*
        NEL_IGM1 = IOCCLS(1)
        DO IGAS = 2, NGAS
          IF(NTEST.GE.1000) WRITE(6,*) ' IGAS = ', IGAS
*. Find product of configurations in GAS 1.. IGAS-1, and those of IGAS
          IF(IGAS.EQ.2) THEN
*. The configurations of GAS1... IGAS-1 are those of GAS1, 
            ISBCLS1 = ISBCLS(1)
            CALL ICOPVE(NSBCNF(1,1,ISBCLS1),int_mb(KLNSB1),
     &           (MAXOP+1)*NSMOB)
            CALL ICOPVE(IBSBCNF(1,1,ISBCLS1),int_mb(KLIBSB1),
     &           (MAXOP+1)*NSMOB)
C IFRMR(WORK,IROFF,IELMNT)
c.. dongxia replaced ifrmr
c           LLOCC = IFRMR(WORK,KLSBCNF,ISBCLS1)
            llocc = int_mb(klsbcnf + isbcls1 -1)
            CALL ICOPVE(int_mb(IBOCSBCNF(ISBCLS1)),
     &                  int_mb(KLOCCSB1),LLOCC)
          END IF
*. Info for GAS IGAS
          NEL_IG = IOCCLS(IGAS)
          ISBCLSI = ISBCLS(IGAS)
          IF(IGAS.NE.NGAS) THEN
*. Prepare dimension array for configs for space 1 - IGAS, in KLNSB2
            CALL DIM_PROD_OCSBCNF(int_mb(KLNSB2),int_mb(KLNSB1),
     &           NSBCNF(1,1,ISBCLSI),MAXOP+1,MAXOP)
C     DIM_PROD_OCSBCNF(I12DIM,I1DIM,I2DIM,NROW2,IOPEN_OUT_DIM)
*. Offset array for configurations for space 1 - IGAS, in KLIBSB2
            CALL GEN_IBOCC_SBCNF(int_mb(KLNSB2),int_mb(KLIBSB2),MAXOP,
     &         NSMOB,NEL_IGM1+NEL_IG)
*. And obtain the configurations in space 1 - IGAS in KLOCCSB2
*. Zero dimension for products - reconstructed in PROD_SBCNF
            IZERO = 0
            CALL ISETVC(int_mb(KLNSB2),IZERO,(MAXOP+1)*NSMOB)
            IF(NTEST.GE.1000) WRITE(6,*) ' NEL_IGM1, NEL_IG = ',
     &                   NEL_IGM1, NEL_IG
            CALL PROD_SBCNF(int_mb(KLOCCSB2),int_mb(KLNSB2),
     &                      int_mb(KLIBSB2),int_mb(KLNSB1),
     &                      NSBCNF(1,1,ISBCLSI),
     &                      int_mb(KLIBSB1), IBSBCNF(1,1,ISBCLSI),
     &                      int_mb(KLOCCSB1),
     &                      int_mb(IBOCSBCNF(ISBCLSI)),0,MAXOP,0,
     &                      NSMOB,NEL_IGM1,NEL_IG)
C                PROD_SBCNF(IOCC12,NCONF12,IBCONF12, NCONF1, NCONF2,
C    &           IBCONF1, IBCONF2, IOCC1, IOCC2, MINOP, MAXOP,ISM12,
C    &           NSMOB,NEL1,NEL2)
*. Copy the new configurations to *1 arrays
            LOCC2 = 
     &      LOCC_SBCNF(int_mb(KLNSB2),MAXOP,NSMOB,NEL_IGM1 + NEL_IG)
C           LOCC_SBCNF(NOPSM,NOPEN_MAX,NSMOB,NELEC)
*. Could be done simpler by just switching pointers, but for initial debugging
            CALL ICOPVE(int_mb(KLOCCSB2),int_mb(KLOCCSB1),LOCC2)
            CALL ICOPVE(int_mb(KLNSB2),int_mb(KLNSB1),(MAXOP+1)*NSMOB)
            CALL ICOPVE(int_mb(KLIBSB2),int_mb(KLIBSB1),(MAXOP+1)*NSMOB)
          ELSE
* IGAS  = NGAS, the final time around, obtain only configurations of the specified symmetry
            CALL PROD_SBCNF(IOCC,NOPSML,IBOPSML,
     &                      int_mb(KLNSB1),NSBCNF(1,1,ISBCLSI),
     &                      int_mb(KLIBSB1),IBSBCNF(1,1,ISBCLSI),
     &                      int_mb(KLOCCSB1),
     &                      int_mb(IBOCSBCNF(ISBCLSI)),MINOP,MAXOP,ISM,
     &                      NSMOB,NEL_IGM1,NEL_IG)
          END IF
          NEL_IGM1 = NEL_IGM1 + IOCCLS(IGAS)
        END DO ! loop over IGAS 
      END IF ! several gas spaces
*. And copy total number of confs back
      CALL ICOPVE(NOPSML(1+(MAXOP+1)*(ISM-1)),NCONF_OP,MAXOP+1)
C?    WRITE(6,*) ' Return from GEN_OCCONF_FOR.. '

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CNFRSB')
      RETURN
      END
      SUBROUTINE PROD_SBCNF(IOCC12,NCONF12,IBCONF12, NCONF1, NCONF2,
     &           IBCONF1, IBCONF2, IOCC1, IOCC2, MINOP, MAXOP,ISM12,
     &           NSMOB,NEL1,NEL2)
*
* Two list of subconfiguration occupation are given in IOCC1, IOCC2
* Obtain the occupations of the product configurations.
* The configurations are sorted as OCC(ICONF,IOP,ISM)
* 
* If ISM12 = 0, then all symmetries are generated, else only ISM12
*
*. Jeppe Olsen, April 3, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC1(*), IOCC2(*)
      INTEGER NCONF1(MAXOP+1,NSMOB), NCONF2(MAXOP+1,NSMOB)
      INTEGER IBCONF1(MAXOP+1,NSMOB), IBCONF2(MAXOP+1,NSMOB)
      INTEGER IBCONF12(MAXOP+1,NSMOB)
*. Input and output
      INTEGER IOCC12(*), NCONF12(MAXOP+1,NSMOB)
*. Symmetry
      INCLUDE 'multd2h.inc'
*. 
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*)  ' Info from PROD_SBCNF' 
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' NEL1, NEL2 = ', NEL1, NEL2
        WRITE(6,*) ' ISM12 =  ', ISM12
        WRITE(6,*) ' MINOP, MAXOP = ', MINOP, MAXOP
      END IF
*
       IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The NCONF1, NCONF2 arrays '
        CALL IWRTMA(NCONF1,MAXOP+1,NSMOB,MAXOP+1,NSMOB)
        WRITE(6,*)
        CALL IWRTMA(NCONF2,MAXOP+1,NSMOB,MAXOP+1,NSMOB)
        WRITE(6,*) ' The initial NCONF12 array '
        CALL IWRTMA(NCONF12,MAXOP+1,NSMOB,MAXOP+1,NSMOB)
        WRITE(6,*) ' The IBCONF12 array '
        CALL IWRTMA(IBCONF12,MAXOP+1,NSMOB,MAXOP+1,NSMOB)
      END IF
*
      IF(ISM12.EQ.0) THEN
        MINSM12 = 1
        MAXSM12 = NSMOB
      ELSE
        MINSM12 = ISM12
        MAXSM12 = ISM12
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' MINSM12, MAXSM12 = ', MINSM12, MAXSM12
      END IF
*
      NEL12 = NEL1 + NEL2
      DO JSM12 = MINSM12, MAXSM12
        DO JOP12 = MINOP, MAXOP
C?        WRITE(6,*) ' JSM12, JOP12 = ', JSM12, JOP12
          IF(MOD(NEL12-JOP12,2).EQ.0) THEN
            JCL12 = (NEL12-JOP12)/2 
            DO JOP1 = 0, JOP12
              JOP2 = JOP12-JOP1
C?            WRITE(6,*) ' JOP12, JOP1, JOP2 = ',
C?   &                     JOP12, JOP1, JOP2 
              JCL1 = (NEL1-JOP1)/2
              JCL2 = (NEL2-JOP2)/2
              JOC1 = JOP1 + JCL1 
              JOC2 = JOP2 + JCL2 
              JOC12 =  JOP12 + JCL12
              DO JSM1 =1, NSMOB
                JSM2 = MULTD2H(JSM12,JSM1)
C?              WRITE(6,*) ' JSM1, JSM2, JSM12 = ', 
C?   &          JSM1, JSM2, JSM12
C?              WRITE(6,*) ' NCONF1(JOP1+1, JSM1) =',
C?   &          NCONF1(JOP1+1, JSM1)
C?              WRITE(6,*) ' NCONF2(JOP2+1,JSM2) =',
C?   &          NCONF2(JOP2+1,JSM2)
    
                DO JCNF1 = 1, NCONF1(JOP1+1, JSM1)
                  DO JCNF2 = 1, NCONF2(JOP2+1,JSM2)
                    IB1 = IBCONF1(JOP1+1,JSM1)+(JCNF1-1)*JOC1
                    IB2 = IBCONF2(JOP2+1,JSM2)+(JCNF2-1)*JOC2
                    IF(NTEST.GE.1000) THEN
                      WRITE(6,*) ' Input confs 1, 2 ' 
                      CALL IWRTMA(IOCC1(IB1),1,JOC1,1,JOC1)
                      CALL IWRTMA(IOCC2(IB2),1,JOC2,1,JOC2)
                    END IF
                    NCONF12(JOP12+1,JSM12) = NCONF12(JOP12+1,JSM12) + 1
                    IB12 = IBCONF12(JOP12+1,JSM12) 
     &                   + (NCONF12(JOP12+1,JSM12)-1)*JOC12
                    DO JOB1 = 1, JOC1
                      IOCC12(IB12-1+JOB1) = IOCC1(IB1-1+JOB1)
                    END DO
                    DO JOB2 = 1, JOC2
                      IOCC12(IB12-1+JOC1+JOB2) = IOCC2(IB2-1+JOB2)
                    END DO
                    IF(NTEST.GE.1000) THEN
                      WRITE(6,*) ' Merged configuration '
                      CALL IWRTMA(IOCC12(IB12),1,JOC12,1,JOC12)
                      WRITE(6,'(A,4I2)') ' JOP1, JOP2, JSM1, JSM2 =',
     &                JOP1, JOP1, JSM1,  JSM2
                      WRITE(6,*) ' JCNF1, JCNF2 = ', JCNF1, JCNF2

                    END IF
*. Check number of electrons
                    CALL CHECK_NEL_IN_CONF(IOCC12(IB12),JOC12,NEL12)
C                   CHECK_NEL_IN_CONF(ICONF,NOC,NEL)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input subconf list 1: '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        CALL WRT_SBCNF_LIST(IOCC1, NCONF1, IBCONF1, MAXOP, NSMOB,
     &                      1, NSMOB,NEL1)
        WRITE(6,*)
        WRITE(6,*) ' Input subconf list 2: '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        CALL WRT_SBCNF_LIST(IOCC2, NCONF2, IBCONF2, MAXOP, NSMOB,
     &                      1, NSMOB,NEL2)
        WRITE(6,*)
        WRITE(6,*) ' Output subconf list: '
        WRITE(6,*) ' ====================='
        WRITE(6,*)
        CALL WRT_SBCNF_LIST(IOCC12, NCONF12, IBCONF12, MAXOP, NSMOB,
     &                      MINSM12, MAXSM12,NEL12)
      END IF
*
      RETURN
      END
      SUBROUTINE WRT_SBCNF_LIST(ICONF,NCONF_OP_SM,IBCONF_OP_SM,
     &           MAXOP,NSMOB,MINSM,MAXSM,NELEC)
*
* Write list of subconfigurations, given in packed form 
*
* The observant reader may wonder: what is the difference between 
* configuration and subconfiguation lists. Well, a subconfiguration
* list contains in general several symmetries. I hope that this should 
* answer the question of the curious 
* (or more likely, forgetsome programmer...)
*
* Jeppe Olsen, April 4, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICONF(*), NCONF_OP_SM(MAXOP+1,NSMOB)
      INTEGER IBCONF_OP_SM(MAXOP+1,NSMOB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from WRT_SBCNF_LIST '
        WRITE(6,*) ' MINSM, MAXSM = ', MINSM, MAXSM
        WRITE(6,*) ' MAXOP = ', MAXOP
      END IF
*
      DO ISM = MINSM, MAXSM
       DO IOPEN = 0, MAXOP
        NCONF_OPSM = NCONF_OP_SM(IOPEN+1,ISM)
        IF(NCONF_OPSM.NE.0) THEN
          WRITE(6,'(A,I3,A,I3,A,I6)') 
     &    ' Number of configurations with symmetry ', ISM, ' and ',
     &     IOPEN, ' open orbitals is ', NCONF_OPSM
          NOCC_ORB = IOPEN + (NELEC-IOPEN)/2
          DO JCONF = 1, NCONF_OPSM
            IB = IBCONF_OP_SM(IOPEN+1,ISM) + (JCONF-1)*NOCC_ORB
            CALL IWRTMA(ICONF(IB),1,NOCC_ORB,1,NOCC_ORB)
          END DO
        END IF
       END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GEN_OCC_SBCNF(NSBCNF_FOR_OP_SM,IBSBCNF_FOR_OP_SM,
     &                         IOGOCSBCLS,MINOPFSBCLS,KOCSBCNF)
*
*. Generate the occupations of the subconfigurations
*
*. Jeppe Olsen, April 4, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cgas.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IOGOCSBCLS(2,NOCSBCLST)
      INTEGER MINOPFSBCLS(NOCSBCLST)
*
      INTEGER NSBCNF_FOR_OP_SM(NOPEN_MAX+1,NSMOB,NOCSBCLST)
      INTEGER IBSBCNF_FOR_OP_SM(NOPEN_MAX+1,NSMOB,NOCSBCLST)
*. Pointer to occupation of subclass configs
      INTEGER KOCSBCNF(NOCSBCLST)
*
      CALL LUCIAQENTER('GNOCSB')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Entering GEN_OCC_SBCNF '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' The MINOPFSBCLS array '
        CALL IWRTMA(MINOPFSBCLS,1,NOCSBCLST,1,NOCSBCLST)
        WRITE(6,*) ' The IOGOCSBCLS array '
        CALL IWRTMA(IOGOCSBCLS,2,NOCSBCLST,2,NOCSBCLST)
      END IF
*
      DO IOCSBCLS = 1, NOCSBCLST
         NEL = IOGOCSBCLS(1, IOCSBCLS)
         IGAS = IOGOCSBCLS(2, IOCSBCLS)
         NOB = NOBPT(IGAS)
         MINOPL = MINOPFSBCLS(IOCSBCLS)
         IF(NTEST.GE.1000) THEN
          WRITE(6,'(A,5I4)') ' IOCSBCLS, NEL, IGAS, NOB, MINOPL = ',
     &    IOCSBCLS, NEL, IGAS, NOB, MINOPL
         END IF
         IBOB = 1
         DO JGAS = 0, IGAS -1 
          DO JSM = 1, NSMOB
           IBOB = IBOB + NOBPTS_GN(JGAS,JSM)
          END DO
         END DO
         IDUM = 0
         IBOB_MINA = IBOB-NINOB
         CALL GEN_SUBCONF(NEL,NOB,ISMFTO(IBOB),MINOPL,
     &        NSBCNF_FOR_OP_SM(1,1,IOCSBCLS),NSMOB,
     &        IBSBCNF_FOR_OP_SM(1,1,IOCSBCLS),IDUM,
     &        int_mb(KOCSBCNF(IOCSBCLS)),1,NOPEN_MAX,IBOB_MINA)
      END DO 
*
      CALL LUCIAQEXIT('GNOCSB')
      RETURN
      END
      SUBROUTINE ADD_CNFS_TO_LIST(IOCCOUT,NCONFOUT,IBCONFOUT, 
     &           IOCCADD, NCONFADD,IBCONFADD,MAXOP,ISM,NSMOB,NEL)
*
* 
*
* A list of configuration occupations are given BY IOCCOUT, NCONFOUT, IBCONFOUT.
* Add configurations defined by IOCCADD, NCONFADD, IBCONFADD
*
* 
* If ISM12 = 0, then all symmetries are generated, else only ISM12
*
*. Jeppe Olsen, April 10, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCCADD(*)  
      INTEGER NCONFADD(MAXOP+1,NSMOB),IBCONFADD(MAXOP+1,NSMOB)
      INTEGER IBCONFOUT(MAXOP+1,NSMOB), IBCONF2(MAXOP+1,NSMOB)
*. Input and output
      INTEGER IOCCOUT(*), NCONFOUT(MAXOP+1,NSMOB)
*. Symmetry
      INCLUDE 'multd2h.inc'
*. 
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)  ' Info from ADD_CNFS_TO_LIST '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' ISM,  NEL = ', ISM, NEL
      END IF
*
      IF(ISM.EQ.0) THEN
        MINSM = 1
        MAXSM = NSMOB
      ELSE
        MINSM = ISM
        MAXSM = ISM
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' MINSM, MAXSM = ', MINSM, MAXSM
      DO JSM = MINSM, MAXSM
        DO JOP = 0, MAXOP
          IF(NTEST.GE.1000) WRITE(6,*) ' JOP = ', JOP
          IF(MOD(NEL-JOP,2).EQ.0) THEN
            JCL = (NEL-JOP)/2 
            JOC = JOP + JCL
            DO JCNFADD = 1, NCONFADD(JOP+1, JSM)
              IBADD = IBCONFADD(JOP+1,JSM)+(JCNFADD-1)*JOC
              NCONFOUT(JOP+1,JSM) = NCONFOUT(JOP+1,JSM) + 1
              IBOUT = IBCONFOUT(JOP+1,JSM) 
     &              + (NCONFOUT(JOP+1,JSM)-1)*JOC
*
              IF(NTEST.GE.10000) THEN
                WRITE(6,*) ' IBCONFADD(JOP+1,JSM) =  ',
     &          IBCONFADD(JOP+1,JSM)
                WRITE(6,*) ' NCONFOUT(JOP+1,JSM) = ',
     &          NCONFOUT(JOP+1,JSM)
                WRITE(6,*) ' IBCONFOUT(JOP+1,JSM) = ',
     &          IBCONFOUT(JOP+1,JSM)
                WRITE(6,*) ' IBADD, IBOUT  = ', IBADD, IBOUT
                WRITE(6,*) ' JOC = ', JOC
              END IF
*
              DO JOB = 1, JOC
                IOCCOUT(IBOUT-1+JOB) = IOCCADD(IBADD-1+JOB)
              END DO
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Updated conf list from ADD_CNFS_TO_LIST'
        CALL WRT_SBCNF_LIST(IOCCOUT, NCONFOUT, IBCONFOUT, MAXOP, NSMOB,
     &                      MINSM, MAXSM,NEL)
C            WRT_SBCNF_LIST(ICONF,NCONF_OP_SM,IBCONF_OP_SM,
C    &           MAXOP,NSMOB,MINSM,MAXSM,NELEC)
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_IB_SBCNF(NOPSM,IBOPSM,NOPEN_MAX,NSMOB)
*
* Offset array for subconfs, all syms
*
*. Jeppe Olsen, April 11, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NOPSM(NOPEN_MAX+1,NSMOB)
*. Output
      INTEGER IBOPSM(NOPEN_MAX+1,NSMOB)
*
* Configurations are ordered as 
*. Loop over symmetries
*.  Loop over number of open orbtals
*
      NTEST = 00
*
      IB = 1
      DO ISM = 1, NSMOB
        DO IOP = 0, NOPEN_MAX
          IBOPSM(IOP+1,ISM) = IB
          IB = IB + NOPSM(IOP+1,ISM)
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GEN_IB_SBCNF reporting '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' NOPSM array (input) '
        WRITE(6,*)
        CALL IWRTMA(NOPSM,NOPEN_MAX+1,NSMOB,NOPEN_MAX+1,NSMOB)
        WRITE(6,*)
        WRITE(6,*) ' IBOPSM array (output) '
        WRITE(6,*)
        CALL IWRTMA(IBOPSM,NOPEN_MAX+1,NSMOB,NOPEN_MAX+1,NSMOB)
      END IF
*
      RETURN
      END
      FUNCTION LOCC_SBCNF(NOPSM,NOPEN_MAX,NSMOB,NELEC)
*
* A subconfiguration list (contains all symmetries) is defined by dimension array NOPSM
* Obtain length of the corresponding occupation list
*
*. Jeppe Olsen, April 23 (Hanging over the Atlantic on the way to Minneapolis)
*
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NOPSM(NOPEN_MAX+1,NSMOB)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from LOCC_SBCNF: Input NOPSM: '
        CALL IWRTMA(NOPSM,NOPEN_MAX+1,NSMOB,NOPEN_MAX+1,NSMOB)
      END IF
*
      LOCC = 0
      DO IOP = 0, NOPEN_MAX
       DO ISM = 1, NSMOB
         IOC = IOP + (NELEC-IOP)/2
         LOCC = LOCC + IOC*NOPSM(IOP+1,ISM)
       END DO
      END DO
*
      LOCC_SBCNF = LOCC
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Size of subconfiguration occupations = ', LOCC
      END IF
*
      RETURN
      END
      SUBROUTINE REO_FOR_CONFS(IOCC,MAXOP,NSMOB,
     &           NCONF_OP1,NCONF_OP2,IBCONF_OP,IBOCCCONF_OP,IB_OCCLS,
     &           IZCONF,NORBT,NELEC,IREO)
*
* Obtain reordering array for a set of configurations 
*
* Input
* IOCC: occupation of configurations
* NCONF_OP1: Start of configurations to be adressed for a given
*            number of open electrons
* NCONF_OP2: End of configurations to be addressed for a given number of 
*            open orbitals
* IBCONF_OP: Offset to configurations with a given number of open orbitals for CI_SPACE
* IBOCCCONF_OP: Offset to occupation of configurations with a given number of 
*               open orbitals for CI-SPACE
* NTORB: Total number of (active) orbitals
* NELEC: Number of active electrons
* IB_OCCLS: Overall offset to configurations of this occupation class
* IZCONF: Lexical adressing scheme for configurations
*
*. Output
*  IREO: Updated reorder array
*
*. Jeppe Olsen, April 2013, Commons Hotel, Minneapolis
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IBCONF(MAXOP+1),NCONF_OP1(MAXOP+1),NCONF_OP2(MAXOP+1)
      INTEGER IBCONF_OP(MAXOP+1),IBOCCCONF_OP(MAXOP+1)
      INTEGER IOCC(*)
      INTEGER IZCONF(*)
*. Output (updated) 
      INTEGER IREO(*)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Entering REO_FOR_CONFS '
        WRITE(6,*) ' ======================='
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NCONF_OP1(Prev. number confs per open)'
        CALL IWRTMA(NCONF_OP1,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*) ' NCONF_OP2(Current number of confs per open)'
        CALL IWRTMA(NCONF_OP2,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*) ' The IBOCCCONF_OP array '
        CALL IWRTMA(IBOCCCONF_OP,1,MAXOP+1,1,MAXOP+1)
      END IF
*
      ILEXNUM_MAX = 0
      DO IOP = 0, MAXOP
        N1OP = NCONF_OP1(IOP+1)
        N2OP = NCONF_OP2(IOP+1)
        NOC = IOP + (NELEC-IOP)/2
        IF(NTEST.GE.1000) 
     &  WRITE(6,*) ' IOP, N1OP, N2OP = ', IOP, N1OP, N2OP
        DO ICONF = N1OP+1, N2OP
*. Actual address
          IACT = IBCONF_OP(IOP+1) - 1 + ICONF
          IBOCC = IBOCCCONF_OP(IOP+1) + (ICONF-1)*NOC
*. Check: Number of electrons
          NELACT = NELEC_IN_CONF(IOCC(IBOCC),NOC)
C         NELEC_IN_CONF(ICONF,NOC)
*
          IF(NELACT.NE.NELEC) THEN
           WRITE(6,*) ' Configuration with wrong number of electrons'
           WRITE(6,*) ' Configuration: ',IBOCC,NELACT,NELEC
           CALL IWRTMA(IOCC(IBOCC),1,NOC,1,NOC)
           WRITE(6,*) ' IOP, NOC = ', IOP, NOC
           STOP ' Configuration with wrong number of electrons'
          END IF
*. Lexical address
          ILEXNUM = ILEX_FOR_CONF(IOCC(IBOCC),NOC,NORBT,NELEC,
     &    IZCONF,0,IDUM)
C         ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
          ILEXNUM_MAX = MAX(ILEXNUM,ILEXNUM_MAX)
          IREO(IB_OCCLS-1+ILEXNUM) = IACT
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IACT, ILEXNUM = ', IACT, ILEXNUM
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Offset and largest lexical number = ', 
     &               IB_OCCLS, ILEXNUM_MAX
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Updated reorder list actual from lexical to actual address'
        WRITE(6,*)
     &  ' =========================================================='
        WRITE(6,*)
        IDIM = IB_OCCLS - 1 + ILEXNUM_MAX
        CALL IWRTMA(IREO,1,IDIM,1,IDIM)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_MAX_CONF_DIMS
*
* Obtain various max dims for configurations
*
*. Jeppe Olsen, April 2013, Minneapolis, for the subconf codes
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
*
      NTEST = 000
*
      NCONF_MAX = IMNMX_IMAT(NCONF_PER_SYM_GN,
     &            NSMOB,NCISPC,MXPCSM,MXPICI,2)
      NCSF_MAX = IMNMX_IMAT(NCSF_PER_SYM_GN,
     &           NSMOB,NCISPC,MXPCSM,MXPICI,2)
      NSD_MAX = IMNMX_IMAT(NSD_PER_SYM_GN,
     &          NSMOB,NCISPC,MXPCSM,MXPICI,2)
      NCM_MAX = IMNMX_IMAT(NSD_PER_SYM_GN,
     &          NSMOB,NCISPC,MXPCSM,MXPICI,2)
      LCONFOCC_MAX = IMNMX_IMAT(LCONFOCC_PER_SYM_GN,
     &          NSMOB,NCISPC,MXPCSM,MXPICI,2)
      NCONF_AS_MAX = IMNMX(NCONF_ALL_SYM_GN,NCISPC,2)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCONF_MAX = ', NCONF_MAX
        WRITE(6,*) ' NCSF_MAX = ', NCSF_MAX
        WRITE(6,*) ' NSD_MAX = ', NSD_MAX
        WRITE(6,*) ' NCM_MAX = ', NCM_MAX
        WRITE(6,*) ' LCONFOCC_MAX = ', LCONFOCC_MAX
        WRITE(6,*) ' NCONF_AS_MAX = ', NCONF_AS_MAX
      END IF
*
      RETURN
      END
      FUNCTION IMNMX_IMAT(IMAT,NRA,NCA,NR,NC,IMNMX)
*
* Find Min or Max elements of a integer matrix
*
*. Jeppe Olsen, April 2013, Minneapolis 
*
      INCLUDE 'implicit.inc'
      INTEGER IMAT(NR,NC)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Info from MNMX_IMAT: '
      END IF
*
      IEXTR = IMAT(1,1)
      DO IC = 1, NCA
       DO IR = 1, NRA
         IF(IMNMX.EQ.1) THEN
           IEXTR = MIN(IEXTR,IMAT(IR,IC))
         ELSE
           IEXTR = MAX(IEXTR,IMAT(IR,IC))
         END IF
       END DO
      END DO
*
      IMNMX_IMAT = IEXTR
*
      IF(NTEST.GE.10) THEN
       IF(IMNMX.EQ.1) THEN
         WRITE(6,*) ' Smallest value ', IEXTR
       ELSE 
         WRITE(6,*) ' Largest value ', IEXTR
       END IF
      END IF
*
      RETURN
      END
      SUBROUTINE GEN_IBOCC_SBCNF(NOPSM,IBOCCOPSM,NOPEN_MAX,NSMOB,NEL)
*
* Offset array for OCCUPATIONS of subconfs, all syms
*
*. Jeppe Olsen, April 27, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NOPSM(NOPEN_MAX+1,NSMOB)
*. Output
      INTEGER IBOCCOPSM(NOPEN_MAX+1,NSMOB)
*
* Configurations are ordered as 
*. Loop over symmetries
*.  Loop over number of open orbtals
*
      NTEST = 00
*
      IB = 1
      DO ISM = 1, NSMOB
        DO IOP = 0, NOPEN_MAX
          IBOCCOPSM(IOP+1,ISM) = IB
          IF(MOD(NEL-IOP,2).EQ.0) THEN
            IOC = IOP + (NEL-IOP)/2
            IB = IB + NOPSM(IOP+1,ISM)*IOC
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GEN_IBOCC_SBCNF reporting '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' NOPSM array (input) '
        WRITE(6,*)
        CALL IWRTMA(NOPSM,NOPEN_MAX+1,NSMOB,NOPEN_MAX+1,NSMOB)
        WRITE(6,*)
        WRITE(6,*) ' IBOCCOPSM array (output) '
        WRITE(6,*)
        CALL IWRTMA(IBOCCOPSM,NOPEN_MAX+1,NSMOB,NOPEN_MAX+1,NSMOB)
      END IF
*
      RETURN
      END
      FUNCTION NELEC_IN_CONF(ICONF,NOC)
*
* A configuration is given in compact form
* Determine number of electrons
* 
*. Jeppe Olsen, April 27, 2013, Still debugging after all these years
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICONF(NOC)
*
      NEL = NOC
      DO IOC = 1, NOC
       IF(ICONF(IOC).LT.0) NEL = NEL + 1
      END DO
*
      NELEC_IN_CONF = NEL
*
      RETURN
      END
      SUBROUTINE CHECK_NEL_IN_CONF(ICONF,NOC,NEL)
*
* Check to see whether configuration ICONF has NEL electrons
*
*. Jeppe Olsen, April 27, Minneapolis
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICONF(NOC)
*
      NELA = NELEC_IN_CONF(ICONF,NOC)
*
      IF(NELA.NE.NEL) THEN
        WRITE(6,*) ' Configuration with wrong number of electrons'
        CALL IWRTMA(ICONF,1,NOC,1,NOC)
        WRITE(6,*) ' Expected and actual number of electrons: ',
     &  NEL, NELA
        STOP ' Configuration with wrong number of electrons'
      END IF
*
      RETURN
      END
      SUBROUTINE REO_SD_FOR_OCCLS_SN(IOCCLS_NUM,ISYM,IREO,
     &           NOCTPA,NOCTPB,IOCTPA,IOCTPB,
     &           NAEL,NBEL,
     &           NSMST,NGAS,IB_ORB,NACOB,NOCOB,PSSIGN,MINOP,
     &           NTOOB,NOBPT,IBLTP,ISMOST,
     &           NSSOA,NSSOB,
     &           ICONF_REO,IZCONF,
     &           NPCMCNF,DFTP,
     &           KZ_PTDT,KREO_PTDT,
     &           IB_CN_OPEN, IB_CM_OPEN,
     &           NABSPGP_FOR_OCCLS,IBABSPGP_FOR_OCCLS,
     &           IABSPGP_FOR_OCCLS,
     &           IASTR,IBSTR,
     &           IDET_OC,IDET_MS,IDET_VC,ILCHK)
*
* Jeppe Olsen, Minneapolis, Apr. 28, 2013: New version with 
*                           hopefully improved performance

*
* Reorder determinants in GAS space from det to configuration order
* for determinants in OCClass IOCCLS
* Reorder array created is Conf-order => AB-order 
*
c     IMPLICIT REAL*8(A-H,O-Z)
      include 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. General input
*.---------------
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)  
      INTEGER IBLTP(*), ISMOST(*)
      INTEGER NOBPT(*)
      INTEGER DFTP(*) 
      INTEGER NPCMCNF(*)
*. Info on the AB supergroups
      INTEGER NABSPGP_FOR_OCCLS(*), IABSPGP_FOR_OCCLS(2,*),
     &        IBABSPGP_FOR_OCCLS(*)
*. IB_CN_OPEN, IB_CM_OPEN(IOPEN+1) Gives start of confs/CM's
*. with given symmetry and number of open orbitals
      INTEGER IB_CN_OPEN(*), IB_CM_OPEN(*)
*. Info for the lexical adressing of configurations 
      INTEGER ICONF_REO(*), IZCONF(*)
*. WORK(KZ_PTDT(IOPEN+1) gives Z  array for prototype dets with IOPEN 
*. WORK(KREO_PTDT(IOPEN+1) gives the corresponding reorder array
*. open orbitals
      INTEGER KZ_PTDT(*), KREO_PTDT(*)
*. The work array used for WORK(KZ_PTDET()),WORK(KREO_PTDT())
c     DIMENSION WORK(*)
*. Scratch space 
*. --------------
      INTEGER IASTR(NAEL,*),IBSTR(NBEL,*)
      INTEGER IDET_OC(*), IDET_MS(*) , IDET_VC(*)
*. Output
*. ------
      INTEGER IREO(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output from REO_SD_OCCLS_SN '
        WRITE(6,*) ' ============================'
        WRITE(6,*) 
        WRITE(6,*) ' Number of occupation class in action: ', IOCCLS_NUM
      END IF
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' ILCHK = ', ILCHK
        WRITE(6,*) ' PSSIGN = ', PSSIGN
        WRITE(6,*) ' IB_ORB = ', IB_ORB
      END IF
*
      IAGRP = 1
      IBGRP = 2
*
      NEL = NAEL + NBEL
*
      IDET = 0
      N_AB = NABSPGP_FOR_OCCLS(IOCCLS_NUM)
      IB_AB = IBABSPGP_FOR_OCCLS(IOCCLS_NUM)
      DO I_AB = IB_AB, IB_AB + N_AB - 1
       IATP = IABSPGP_FOR_OCCLS(1,I_AB)
       IBTP = IABSPGP_FOR_OCCLS(2,I_AB)
       DO IASM = 1, NSMST
        IBSM = ISMOST(IASM)
        IF(IBLTP(IASM).EQ.1.OR.(IBLTP(IASM).EQ.2.AND.IATP.GE.IBTP)) THEN
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         IF(NTEST.GE.10000) THEN
           WRITE(6,'(A,5(2X,I6))')
     &     ' I_AB, IATP, IBTP, IASM, IBSM = ',
     &       I_AB, IATP, IBTP, IASM, IBSM 
         END IF
*
*. Obtain alpha strings of sym IASM and type IATP
*
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                            NORB,0,IDUM,IDUM)
*
*. Obtain Beta  strings of sym IBSM and type IBTP
*
         IDUM = 0
         CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                            NORB,0,IDUM,IDUM)
*. Offset to this occupation class in occupation class ordered cnf list
         IB_OCCLS = 1
         IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
         ELSE
          IRESTR = 0
         END IF
*
         DO  IB = 1,NIB
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB 
            IRESTR2 = 1
          ELSE
            MINIA = 1
            IRESTR2 = 0
          END IF
          DO  IA = MINIA,NIA
            IDET = IDET + 1
*
*. Reform from AB to CONF form
*
*. subtract offset(number of inactive orbitals)
            ISUB = 1-IB_ORB
C                 ABSTR_TO_CNFSTR(IA_OC,IB_OC,NAEL,NBEL,
C    &           ICONF,IOP_SP,ISIGN,IADD,NOP,NOC)
            CALL ABSTR_TO_CNFSTR(IASTR(1,IA),IBSTR(1,IB),NAEL,NBEL,
     &                           IDET_OC,IDET_MS,ISIGN,
     &                           ISUB,NOP,NOP_AL,NOC,IOP1)
*
* On output IDET_OC is conf in compact form, IDET_MS is projection of 
* open orbitals
*
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' IASTR, IBSTR, ICONF:'
              CALL IWRTMA(IASTR(1,IA),1,NAEL,1,NAEL)
              CALL IWRTMA(IBSTR(1,IB),1,NBEL,1,NBEL)
              CALL IWRTMA(IDET_OC,1,NOC,1,NOC)
            END IF
* 
            IF(NOP.GE.MINOP) THEN
             NPTDT = NPCMCNF(NOP+1)
*. Address of this configuration 
             ICNF_OUT = ILEX_FOR_CONF2(IDET_OC,NOC,NACOB,NEL,IZCONF,
     &                  1,ICONF_REO,1)
C                      ILEX_FOR_CONF2(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,
C                      IREO,IB_OCCLS)
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' Configuration: '
               CALL IWRTMA(IDET_OC,1,NOC,1,NOC)
               WRITE(6,*) ' Address of configuration in output list',
     &         ICNF_OUT
             END IF
*. Address of spinprojection pattern   
C                    IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
             IPTDT = IZNUM_PTDT(IDET_MS,NOP,NOP_AL,
     &                int_mb(KZ_PTDT(NOP+1)),int_mb(KREO_PTDT(NOP+1)),
     &              1)
             ISIGNP = 1
             IF(IPTDT.EQ.0) THEN
              IF(PSSIGN.NE.0) THEN
*. The determinant was not found among the list of prototype dets. For combinations
*. this should be due to the prototype determinant is the MS- switched determinant, so find
*. address of this and remember sign
                M1 = -1
                CALL ABSTR_TO_CNFSTR(IBSTR(1,IB),IASTR(1,IA),NBEL,NAEL,
     &                           IDET_OC,IDET_MS,ISIGN,
     &                           ISUB,NOP,NOP_AL,NOC,IOP1)
                IPTDT = IZNUM_PTDT(IDET_MS,NOP,NOP_AL,
     &                  int_mb(KZ_PTDT(NOP+1)),int_mb(KREO_PTDT(NOP+1)),
     &               1)
                IF(PSSIGN.EQ.-1.0D0) ISIGNP = -1
              ELSE 
*. Prototype determinant was not found in list
               WRITE(6,*) 
     &         ' Error: Determinant not found in list of protodets'
               WRITE(6,*) 
     &         ' Detected in REO_SD_FOR_OCCLS_S'
              END IF
             END IF
*
             IBCNF_OUT = IB_CN_OPEN(NOP+1)
*
             IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
              WRITE(6,*) ' IB_CM_OPEN(NOP+1) = ',
     &                     IB_CM_OPEN(NOP+1)
              WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
              WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
             END IF
*
             IADR_SD_CONF_ORDER = IB_CM_OPEN(NOP+1) - 1
*
     &                          + (ICNF_OUT-IBCNF_OUT)*NPTDT + IPTDT
             IF(IADR_SD_CONF_ORDER.LE.0) THEN
               WRITE(6,*) ' Problemo, IADR_SD_CONF_ORDER < 0 '
               WRITE(6,*) ' IADR_SD_CONF_ORDER = ', IADR_SD_CONF_ORDER
               WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
               WRITE(6,*) ' IB_CM_OPEN(NOP+1) = ',
     &                     IB_CM_OPEN(NOP+1)
               WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
               WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
C?             CALL XFLUSH(6)
             END IF
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' IADR_SD_CONF_ORDER, ISIGN, IDET = ',
     &                      IADR_SD_CONF_ORDER, ISIGN, IDET
             END IF
             IREO(IADR_SD_CONF_ORDER) = ISIGN*IDET*ISIGNP
             IF(NTEST.GE.10000) THEN
               WRITE(6,*) ' IDET, IADR_SD_CONF_ORDER ',
     &                      IDET, IADR_SD_CONF_ORDER
             END IF
            END IF! Nop .ge. MINOP
          END DO
*         ^ End of loop over alpha strings
         END DO
*        ^ End of loop over beta strings
        END IF! Block should be included
       END DO ! Loop over IASM
      END DO! Loop over AB blocks in occ. class
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Reorder array for SDs, CONF order => string order'
        WRITE(6,*) ' ================================================='
        CALL IWRTMA(IREO,1,IDET,1,IDET)
      END IF
*
*. Check sum of reordering array
      I_DO_CHECKSUM = 0
      IF(I_DO_CHECKSUM.EQ.1) THEN
        ISUM = 0
        DO JDET = 1, IDET
          ISUM = ISUM + ABS(IREO(JDET))
        END DO
        IF(ISUM.NE.IDET*(IDET+1)/2) THEN
          WRITE(6,*) ' Problem with sumcheck in REO_GASDET'
          WRITE(6,'(A,2I9)') 
     &    'Expected and actual value ', ISUM, IDET*(IDET+1)/2
          STOP       ' Problem with sumcheck in REO_GASDET'
        ELSE
          WRITE(6,*) ' Sumcheck in REO_GASDET passed '
        END IF
      END IF !checksum is invoked
*
      RETURN
      END
      SUBROUTINE ABSTR_TO_CNFSTR(IA_OC,IB_OC,NAEL,NBEL,
     &           ICONF,IOP_SP,ISIGN,IADD,NOP,NOP_AL,NOC,
     &           IOP1)
*
* An alpha string (IA) and a betastring (IB) is given. 
* Combine these two strings to 
* 1) The configuration corresponding to this string in ICNF, compact form
* 2) The spin-projection pattern of the open orbitals in IOP_SP
* 3) The number of open orbitals
* On output
* IOP1 = 1 => first open is alpha
* IOP2 = -1 => first open is beta
* IOP1 = 0 => No first open (closed)
*
* A constant IADD is added to the configuration to allow for
* different offsets
* orbitals in ascending order. 
* For doubly occupied orbitals
* the alphaorbital is given first. 

* The phase required to change IA IB into configuration order
* is computes as ISIGN
*
* Jeppe Olsen, April 28, 2013 in Minneapolis
*              Derived from ABSTR_TO_ORDSTR to obtain
*              improved reordering of orbitals
* Last Modification; Jeppe Olsen; May 15 2013;  IOP1 added
* 
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IA_OC(NAEL),IB_OC(NBEL)
*. Output
      INTEGER ICONF(NAEL+NBEL)
      INTEGER IOP_SP(NAEL+NBEL)
*  
      NEXT_AL = 1
      NEXT_BE = 1
      NOP = 0
      NOP_AL = 0
      NOC = 0
      ISIGN = 1
      NPERM = 0
      IOP1 = 0
*. Loop over next electron in outputstring
      IF(NAEL+NBEL.EQ.0) GOTO 1001
 1000 CONTINUE
       NOC = NOC + 1
       IF(NEXT_AL.LE.NAEL.AND.NEXT_BE.LE.NBEL) THEN
*
         IF(IA_OC(NEXT_AL).EQ.IB_OC(NEXT_BE)) THEN
*. Next orbital is doubly occupied
           ICONF(NOC) = -IA_OC(NEXT_AL)-IADD
           NEXT_AL = NEXT_AL + 1
           NEXT_BE = NEXT_BE + 1
C          ISIGN = ISIGN*(-1)**(NAEL-NEXT_AL+1) 
           NPERM = NPERM + NAEL-NEXT_AL+1
         ELSE IF(IA_OC(NEXT_AL).LE.IB_OC(NEXT_BE)) THEN
*. Next electron is alpha electron
           ICONF(NOC) = IA_OC(NEXT_AL) + IADD
           NEXT_AL = NEXT_AL + 1
           NOP = NOP + 1
           NOP_AL = NOP_AL + 1
           IOP_SP(NOP) = +1
           IF(IOP1.EQ.0) IOP1 = 1
         ELSE
*. Next electron is beta electron
           ICONF(NOC) = IB_OC(NEXT_BE) + IADD
           NEXT_BE = NEXT_BE + 1
           NOP = NOP + 1
           IOP_SP(NOP) = -1
C          ISIGN = ISIGN*(-1)**(NAEL-NEXT_AL+1) 
           NPERM = NPERM + NAEL-NEXT_AL+1
           IF(IOP1.EQ.0) IOP1 = -1
         END IF
       ELSE IF(NEXT_BE.GT.NBEL) THEN
*. Next electron is alpha electron
           ICONF(NOC) = IA_OC(NEXT_AL) + IADD
           NEXT_AL = NEXT_AL + 1
           NOP = NOP + 1
           NOP_AL = NOP_AL + 1
           IOP_SP(NOP) = +1
           IF(IOP1.EQ.0) IOP1 = 1
       ELSE IF(NEXT_AL.GT.NAEL) THEN
*. Next electron is beta electron
           ICONF(NOC) = IB_OC(NEXT_BE) + IADD
           NEXT_BE = NEXT_BE + 1
           NOP = NOP + 1
           IOP_SP(NOP) = -1
C          ISIGN = ISIGN*(-1)**(NAEL-NEXT_AL+1) 
           NPERM = NPERM + NAEL-NEXT_AL+1
           IF(IOP1.EQ.0) IOP1 = -1
       END IF
      IF(NEXT_AL.LE.NAEL.OR.NEXT_BE.LE.NBEL) GOTO 1000
*     ^ End of loop over orbital in outputlist
*
      IF(MOD(NPERM,2).EQ.0) THEN
       ISIGN = 1
      ELSE
       ISIGN = -1
      END IF
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' ABSTR to CNFSTR : '
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' Input alpha and beta strings '
        CALL IWRTMA(IA_OC,1,NAEL,1,NAEL)
        CALL IWRTMA(IB_OC,1,NBEL,1,NBEL)
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NOC,1,NOC)
        WRITE(6,*) ' Spin projections of open orbitals'
        CALL IWRTMA(IOP_SP,1,NOP,1,NOP)
        WRITE(6,*) ' IOP1 = ', IOP1
      END IF
*
      RETURN
      END
      SUBROUTINE REO_SD_FOR_OCCLS_SSN(IOCCLS_NUM,ISYM,IREO,
     &           NOCTPA,NOCTPB,IOCTPA,IOCTPB,
     &           NAEL,NBEL,
     &           NSMST,NGAS,IB_ORB,NACOB,NOCOB,PSSIGN,MINOP,
     &           NTOOB,NOBPT,IIBLTP,ISMOST,
     &           NSSOA,NSSOB,IBSSOA,IBSSOB,
     &           ICONF_REO,IZCONF,
     &           NPCMCNF,DFTP,
     &           KZ_PTDT,KREO_PTDT,
     &           IB_CN_OPEN, IB_CM_OPEN,
     &           NABSPGP_FOR_OCCLS,IBABSPGP_FOR_OCCLS,
     &           IABSPGP_FOR_OCCLS,
     &           IAPSTR,IBPSTR,
     &           IDET_OC,IDET_MS,IDET_VC,IDET_MS_AB,
     &           IASPGPGP,IBSPGPGP,
     &           NSTSGP, IBSTSGP,
     &           IADSTSTP,IADSTSTL,IB_SMSMP,IB_SMSML,
     &           ICNOCSTSTL,ICNOCSTSTL_AB,
     &           IADSTSTL_AB,NALSTSTL_AB,
     &           NOPSTSTP,NOPSTSTL,NALSTSTP,NALSTSTL,
     &           IADOPSTSTP,IABSWITCHP,IABSWITCHL,MAXOP,
     &           NCMB_CN,IDC,ILCHK)
*
* Jeppe Olsen, Minneapolis, Apr. 30, 2013: Super new version with 
*                           absolutely improved performance
*                           - the new version of two days ago
*                           disappointed
*
* Purpose:
* =========
*
* Reorder determinants in GAS space from det to configuration order
* for determinants in OCClass IOCCLS
* Reorder array created is Conf-order => AB-order 
*
* Approach: 
* =========
* The alpha- and beta-strings are each divided into 
* two groups, and as much calculation as possible 
* is performed for the separate groups
*
* So an alpha string is alpha = alphap * alphal
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      include 'wrkspc-static.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'multd2h.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
C   NGPSTR(IGAS) IBGPSTR(IGAS) 
*
*. General input
*.---------------
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)  
      INTEGER IIBLTP(*), ISMOST(*)
      INTEGER NOBPT(*)
      INTEGER DFTP(*) 
      INTEGER NPCMCNF(*)
*. Info on the AB supergroups
      INTEGER NABSPGP_FOR_OCCLS(*), IABSPGP_FOR_OCCLS(2,*),
     &        IBABSPGP_FOR_OCCLS(*)
*. IB_CN_OPEN, IB_CM_OPEN(IOPEN+1) Gives start of confs/CM's
*. with given symmetry and number of open orbitals
      INTEGER IB_CN_OPEN(*), IB_CM_OPEN(*)
*. Info for the lexical adressing of configurations 
      INTEGER ICONF_REO(*), IZCONF(*)
*. WORK(KZ_PTDT(IOPEN+1) gives Z  array for prototype dets with IOPEN 
*. WORK(KREO_PTDT(IOPEN+1) gives the corresponding reorder array
*. open orbitals
      INTEGER KZ_PTDT(*), KREO_PTDT(*)
*. The work array used for WORK(KZ_PTDET()),WORK(KREO_PTDT())
c     DIMENSION WORK(*)
*. Groups of the supergroups
      INTEGER IASPGPGP(MXPNGAS,*),IBSPGPGP(MXPNGAS,*)
*. Number of strings per symmetry and group and their offsets
      INTEGER NSTSGP(NSMST,*), IBSTSGP(NSMST,*)

*. Scratch through input
*. ----------------------
*. Number and offsets to p-string 
*. The occupation of p-strings
      INTEGER IAPSTR(*), IBPSTR(*)
      INTEGER IDET_OC(*), IDET_MS(*) , IDET_VC(*), IDET_MS_AB(*)
*. Contribution to configuration lexical address from pairs of P-strings and L-strings
      INTEGER IADSTSTP(*), IADSTSTL(*)
*. And offsets to the above 
      INTEGER IB_SMSMP(NSMST,NSMST), IB_SMSML(NSMST,NSMST)
*. Occupations of L-configuration
      INTEGER ICNOCSTSTL(*), ICNOCSTSTL_AB(*)
*. Number of open orbitals per stst for P and L
      INTEGER NOPSTSTP(*),  NOPSTSTL(*)
*. Number of alpha electrons per stst for p and L
      INTEGER NALSTSTP(*), NALSTSTL(*)
      INTEGER IADSTSTL_AB(*),NALSTSTL_AB(*)
*. P-contribution to address of open proto-type
      INTEGER IADOPSTSTP(*)
*. Are P-or L-part of strings AB switched?
       INTEGER IABSWITCHP(*), IABSWITCHL(*)
*. Output
*. ------
      INTEGER IREO(*)
*
*. Local (?) scratch space
*. ============================
      INTEGER IASTRF(MXPORB), IBSTRF(MXPORB)
*
      INTEGER NAPSTR(MXPNSMST),NBPSTR(MXPNSMST)
      INTEGER IBAPSTR(MXPNSMST),IBBPSTR(MXPNSMST)
*. Run also old route for check
      I_DO_ALSO_OLD = 0
*. If PSSIGN, only prototype dets with first open orbital being alpha
*. are stored. This may require AB switch to match det 
      I_CHECK_OP1 = 0
      IF(PSSIGN.NE.0.0D0) THEN
        I_CHECK_OP1 = 1
      END IF
      IPSSIGN = PSSIGN
*
      NTEST = 0000
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' output from REO_SD_OCCLS_SSN '
        WRITE(6,*) ' ============================'
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' Number of occupation class in action: ', IOCCLS_NUM
        WRITE(6,*) ' I_CHECK_OP1 =', I_CHECK_OP1  
      END IF
      IF(NTEST.GE.10000) THEN
        WRITE(6,*) ' PSSIGN = ', PSSIGN
        WRITE(6,*) ' IB_ORB = ', IB_ORB
        WRITE(6,*) ' NCMB_CN = ', NCMB_CN
        WRITE(6,*) ' ILCHK = ', ILCHK
      END IF
*
      IF(NGAS.EQ.1) THEN
       WRITE(6,*) ' Error: SSN reo should be called only if NGAS>1'
       STOP       ' Error: SSN reo should be called only if NGAS>1'
      END IF
*
*
      IAGRP = 1
      IBGRP = 2
*
      NEL = NAEL + NBEL
*. Initial:
      MINOPP = 0
      MINOPL = 0
      ISUB = 1-IB_ORB
*
      IDET = 0
      N_AB = NABSPGP_FOR_OCCLS(IOCCLS_NUM)
      IB_AB = IBABSPGP_FOR_OCCLS(IOCCLS_NUM)
      IF(NTEST.GE.100)
     &WRITE(6,*) ' N_AB, IB_AB = ', N_AB, IB_AB
      DO I_AB = IB_AB, IB_AB + N_AB - 1
       IATP = IABSPGP_FOR_OCCLS(1,I_AB)
       IBTP = IABSPGP_FOR_OCCLS(2,I_AB)
       IF(IDC.EQ.2.AND.IATP.LT.IBTP) GOTO 9999
       IF(NTEST.GE.100) WRITE(6,*) 'I_AB, IATP, IBTP = ',
     &             I_AB, IATP, IBTP
*. Determine the last space with nonvanishing number of 
*. electrons for the IATP and IBTP
C      LAST_OCCGAS_FOR_SUBGP(IOCSPGP,NGAS,LASTOCC,NLASTEL)
C?     WRITE(6,*) ' IASPGPGP(1,IATP): '
C?     CALL IWRTMA3(IASPGPGP(1,IATP),1,NGAS,1,NGAS)
       CALL LAST_OCCGAS_FOR_SUPGP(IASPGPGP(1,IATP),NGAS,
     &      IALGS,NALEL)
       CALL LAST_OCCGAS_FOR_SUPGP(IBSPGPGP(1,IBTP),NGAS,
     &      IBLGS,NBLEL)
*. At the moment, the last GASpaces should be identical- Else
*. the number of open electrons will not be the one in P + L
*. space
*. We must have one space in P pt
       IF(IALGS.EQ.1) THEN
         NALEL = NALEL - NELFGP(IASPGPGP(IALGS,IATP))
         IALGS = 2
       END IF
       IF(IBLGS.EQ.1) THEN
         NBLEL = NBLEL - NELFGP(IBSPGPGP(IBLGS,IBTP))
         IBLGS = 2
       END IF
*
       IF(IALGS.NE.IBLGS) THEN
        IALGS = MAX(IALGS,IBLGS)
        IBLGS = IALGS
        NALEL = 0
        NBLEL = 0
        DO IGAS = IALGS, NGAS
         NALEL = NALEL + NELFGP(IASPGPGP(IGAS,IATP))
         NBLEL = NBLEL + NELFGP(IBSPGPGP(IGAS,IBTP))
        END DO
       END IF
*
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' IALGS, IBLGS, NALEL, NBLEL = ',
     &                IALGS, IBLGS, NALEL, NBLEL
       END IF
       NLEL = NALEL + NBLEL
    
*. The types of the last GAS spaces
*. Obtain alpha and beta
       IALTP =  IASPGPGP(IALGS,IATP)
       IBLTP =  IBSPGPGP(IBLGS,IBTP)
       DO ITP = IBGPSTR(IBLGS),IBGPSTR(IBLGS)+NGPSTR(IBLGS)-1
        IF(NELFGP(ITP).EQ.NBLEL) IBLTP = ITP
       END DO
*
       IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Types of the last a,b-spaces =',
     &               IALTP, IBLTP 
       END IF
*. Number of electrons in active spaces
       NAPEL = NAEL - NALEL
       NBPEL = NBEL - NBLEL
       NPEL = NAPEL + NBPEL
       IF(NTEST.GE.100) THEN
         WRITE(6,*) ' NAPEL, NAEL, NALEL = ',
     &                NAPEL, NAEL, NALEL
         WRITE(6,*) ' NBPEL, NBEL, NBLEL = ',
     &                NBPEL, NBEL, NBLEL
       END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' =============================='
         WRITE(6,*) ' Assembling of PL information:'
         WRITE(6,*) ' =============================='
         WRITE(6,*)
      END IF
*
*. And then set up the strings in the P spaces, which exclude
*. the last nontrivial space
C      GETSTR_ALLSM_GNSPGP(ISTRTP,NGRPA,IGRPA,NSTRPSM,IBSTRPSM,NEL,ISTR)
       IF(NTEST.GE.1000) WRITE(6,*) ' Info for alpha p strings: '
       CALL GETSTR_ALLSM_GNSPGP(1,IALGS-1,IASPGPGP(1,IATP),
     &      NAPSTR,IBAPSTR,NAPEL,IAPSTR)
       IF(NTEST.GE.10000) WRITE(6,*) ' Info for beta p strings: '
       CALL GETSTR_ALLSM_GNSPGP(1,IBLGS-1,IBSPGPGP(1,IBTP),
     &      NBPSTR,IBBPSTR,NBPEL,IBPSTR)
*
*. Find the contributions to the configuration lexical address of P-strings
*
*. NEW START
       IPDET = 0
       DO IASM = 1, NSMST
        DO IBSM = 1, NSMST
          IB_SMSMP(IASM,IBSM) = IPDET + 1
          DO IBP = 1, NBPSTR(IBSM)
            DO IAP = 1, NAPSTR(IASM)
              IPDET = IPDET + 1
              IF(NTEST.GE.10000) WRITE(6,*) ' Info for IPDET: ',IPDET
              IIPA = (IBAPSTR(IASM)-1 + IAP-1)*NAPEL + 1
              IIPB = (IBBPSTR(IBSM)-1 + IBP-1)*NBPEL + 1
              IF(NTEST.GE.10000) THEN
                WRITE(6,*) ' IAP, IBP, IIPA, IIPB = ',
     &                       IAP, IBP, IIPA, IIPB
              END IF
*
*. Reform P strings from AB to CONF form
*
*. subtract offset(number of inactive orbitals)
              CALL ABSTR_TO_CNFSTR(IAPSTR(IIPA),IBPSTR(IIPB),
     &             NAPEL,NBPEL,IDET_OC,IDET_MS,ISIGNP,
     &             ISUB,NOP,NOP_AL,NOC,IOPP)
              SIGNAB = 1.0D0
              IF(NTEST.GE.10000) WRITE(6,*) ' IOPP = ', IOPP
              IF(I_CHECK_OP1.EQ.1) IABSWITCHP(IPDET) = IOPP
              IF(I_CHECK_OP1.EQ.1.AND.IOPP.EQ.-1)THEN
*. First open orbital is beta, interchange alpha and beta 
                CALL ABSTR_TO_CNFSTR(IBPSTR(IIPB),IAPSTR(IIPA),
     &               NBPEL,NAPEL,IDET_OC,IDET_MS,ISIGNP,
     &               ISUB,NOP,NOP_AL,NOC,IOPP2)
                SIGNAB = PSSIGN
              END IF
*

              IF(NOP.GE.MINOPP) THEN
*. Contribution to Address from P-configuration 
                ICNPLEX = ILEX_FOR_CONF3(IDET_OC,1,NOC,NACOB,NEL,
     &                    IZCONF)
                IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' Strings: '
                  CALL IWRTMA(IAPSTR(IIPA),1,NAPEL,1,NAPEL,1)
                  CALL IWRTMA(IBPSTR(IIPB),1,NBPEL,1,NBPEL,1)
                  WRITE(6,*) ' Configuration: '
                  CALL IWRTMA(IDET_OC,1,NOC,1,NOC)
                  WRITE(6,*) 
     &            ' P-term to lexical address of configuration',
     &            ICNPLEX
*. Add on
C?                WRITE(6,*) ' IPDET, NOP, NOP_AL ',
C?   &                         IPDET, NOP, NOP_AL
C?                WRITE(6,*) ' NAPSTR(IASM), NBPSTR(IBSM) = ',
C?   &                         NAPSTR(IASM), NBPSTR(IBSM)
C?                WRITE(6,*) ' IB_SMSMP(IASM,IBSM) = ',
C?   &                         IB_SMSMP(IASM,IBSM)
                END IF
                IADSTSTP(IPDET) = ISIGNP*ICNPLEX
*. Info on open orbitals
                NOPSTSTP(IPDET) = NOP
                NALSTSTP(IPDET) = NOP_AL
*. contribution to lexical address of proto-type det
                 IADOPSTSTP(IPDET) = IZNUM_PTDT2(
     &           IDET_MS,MAXOP,NOP,1,0,1,NOP_AL,
     &           int_mb(KZ_PTDT(MAXOP+1)),IDUM,0)
C              IZNUM_PTDT2(IAB,NOPEN_DIM,NOPEN,IOPEN1,IALPHA1,IZ1,
C    &         NALPHA,Z,NEWORD,IREORD)
C               IPTDT = IZNUM_PTDT(IDET_MS,NOP,NOP_AL,
C    &                  WORK(KZ_PTDT(NOP+1)),WORK(KREO_PTDT(NOP+1)),1)
              END IF ! number of open in P in range
            END DO !IPA
          END DO !IPB
        END DO !IBSM
       END DO ! IASM
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Contributions to conf-addresses from P '
        WRITE(6,*) ' ====================================== '
        WRITE(6,*)
        CALL IWRTMATMAT(IADSTSTP,NSMST,NSMST,
     &       NAPSTR,NBPSTR,IB_SMSMP)
C            IWRTMATMAT(IAA,NBLR,NBLC,LBLR,LBLC,IB)
      END IF
      CALL MEMCHK2('AFTPIN')
*
*. Find the contributions to the configuration of L-strings
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' >>> Info for L-strings: '
         WRITE(6,*)
       END IF
*
       ILDET = 0
       DO IASM = 1, NSMST
        DO IBSM = 1, NSMST
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IASM, IBSM = ', IASM, IBSM
          END IF
          IB_SMSML(IASM,IBSM) = ILDET + 1
          DO IBL = 1, NSTSGP(IBSM,IBLTP)
            DO IAL = 1, NSTSGP(IASM,IALTP)
              ILDET = ILDET + 1
              IF(NTEST.GE.10000) WRITE(6,*) ' ILDET = ', ILDET
              IBLA = (IBSTSGP(IASM,IALTP)-1+IAL-1)*NALEL+1
              IBLB = (IBSTSGP(IBSM,IBLTP)-1+IBL-1)*NBLEL+1
*
*. L-strings in integer arrays
*
              CALL ICOPVE2(int_mb(KOCSTR(IALTP)),IBLA,
     &             NALEL,IASTRF) 
              CALL ICOPVE2(int_mb(KOCSTR(IBLTP)),IBLB,
     &             NBLEL,IBSTRF) 
*
*. Reform L strings from AB to CONF form
*
              CALL ABSTR_TO_CNFSTR(IASTRF,IBSTRF,
     &             NALEL,NBLEL,IDET_OC,IDET_MS,ISIGNL,
     &             ISUB,NOP,NOP_AL,NOC,IOPL)
              CALL ICOPVE
     &        (IDET_MS,ICNOCSTSTL((ILDET-1)*NLEL+1),NOP)
              IF(I_CHECK_OP1.EQ.1) IABSWITCHL(ILDET) = IOPL
*
              IF(I_CHECK_OP1.NE.0) THEN
*. Info for the AB-switched string
                CALL ABSTR_TO_CNFSTR(IBSTRF,IASTRF,
     &               NBLEL,NALEL,IDET_OC,IDET_MS_AB,ISIGNL_AB,
     &               ISUB,NOP,NOP_AL_AB,NOC,IOPL2)
                CALL ICOPVE
     &          (IDET_MS_AB,ICNOCSTSTL_AB((ILDET-1)*NLEL+1),NOP)
              END IF
*
              IF(NOP.GE.MINOPL) THEN
*. Save the L-occupation
*. And obtain the lexical address 

C     ILEX_FOR_CONF3(ICONF,IEL1,NOCC_ORB,NORB,NEL,IARCW)
                ICNLLEX = ILEX_FOR_CONF3(IDET_OC,NPEL+1,
     &                    NOC,NACOB,NEL,IZCONF) 
                IF(I_CHECK_OP1.NE.0) THEN
*. Info for the AB-switched string
                  IADSTSTL_AB(ILDET) = ISIGNL_AB*ICNLLEX
                  NALSTSTL_AB(ILDET) = NOP_AL_AB
                END IF
*.
                IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' L- strings and configuration:'
                  CALL IWRTMA(IDET_OC,1,NOC,1,NOC)
                  WRITE(6,*) 
                  CALL IWRTMA(IASTRF, 1, NALEL,1, NALEL)
                  CALL IWRTMA(IBSTRF, 1, NBLEL,1, NBLEL)
                  WRITE(6,*) 
     &            ' L-term to address of configuration ',ICNLLEX
                  WRITE(6,*) ' ILDET, ISIGNL, ICNLLEX = ', 
     &                         ILDET, ISIGNL, ICNLLEX
                END IF
                IADSTSTL(ILDET) = ISIGNL*ICNLLEX
                NOPSTSTL(ILDET) = NOP
                NALSTSTL(ILDET) = NOP_AL
*
                IF(I_CHECK_OP1.EQ.1) THEN
                  IADSTSTL_AB(ILDET) = ISIGNL_AB*ICNLLEX
                  NALSTSTL_AB(ILDET) = NOP_AL_AB
                END IF ! I_CHECK_OP1 = 1
              END IF ! number of open in L in range
            END DO !ILA
          END DO !ILB
        END DO !IBSM
       END DO !IASM
*
       IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Contributions to conf-addresses from L '
        WRITE(6,*) ' ====================================== '
        WRITE(6,*)
        CALL IWRTMATMAT(IADSTSTL,NSMST,NSMST,
     &       NSTSGP(1,IALTP),NSTSGP(1,IALTP),IB_SMSML)
C            IWRTMATMAT(IAA,NBLR,NBLC,LBLR,LBLC,IB)
       END IF
       
*. END NEW STUFF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' =============================='
         WRITE(6,*) ' Assembling of AB information:'
         WRITE(6,*) ' =============================='
         WRITE(6,*)
      END IF

* The strings of a given group is stored at WORK(KOCSTR(IGRP))
*. Now: loop over the ab strings so dets come in standard order
       IB = 0
       DO IASM = 1, NSMST
         IBSM = ISMOST(IASM)
         IACT = 1
         IF(IDC.EQ.2.AND.IATP.EQ.IBTP.AND.IASM.LT.IBSM) IACT = 0
         IF(IACT.EQ.1) THEN
          IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
           IRESTR = 1
          ELSE
           IRESTR = 0
         END IF
*. Loop over betastrings !
         IB = 0

         DO IBLSM = 1, NSMST
           IBPSM = MULTD2H(IBLSM,IBSM)
           NBP = NBPSTR(IBPSM)
           NBL = NSTSGP(IBLSM,IBLTP)
           IF(NTEST.GE.1000) THEN
             WRITE(6,*) ' IBSM IBPSM IBLSM = ',
     &       IBSM,IBPSM,IBLSM
           END IF
           DO IBP = 1, NBP
            DO IBL = 1, NBL
             IB = IB + 1
             IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
               MINIA = IB
               IRESTR2 = 1
             ELSE
               MINIA = 1
               IRESTR2 = 0
             END IF
             IA = 0
             DO IALSM = 1, NSMST
              IAPSM = MULTD2H(IALSM,IASM)
              NAP = NAPSTR(IAPSM)
              NAL = NSTSGP(IALSM,IALTP)
              IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' IASM IAPSM IALSM = ',
     &                      IASM,IAPSM,IALSM
              END IF
              DO IAP = 1, NAP
               DO IAL = 1, NAL
                IA = IA + 1
                IF(IA.GE.MINIA) THEN
                 IDET = IDET + 1
*
*. New
*
*. Address of configuration
                 IF(NTEST.GE.1000) WRITE(6,*) ' >>> NEW '
                 IABL = IB_SMSML(IALSM,IBLSM)-1 + (IBL-1)*NAL+IAL
                 IABP = IB_SMSMP(IAPSM,IBPSM)-1 + (IBP-1)*NAP+IAP
                 IF(NTEST.GE.10000) THEN
                  WRITE(6,'(A,3I4)') 'IA, IB, IDET = ', IA, IB, IDET
                  WRITE(6,'(A,3I4)') ' IAL, IBL, IABL = ', IAL,IBL,IABL
                  WRITE(6,'(A,2I4)') 'IABL, IABP = ', IABL, IABP
                  WRITE(6,'(A,2I4)') 'NAP, NAL = ', NAP, NAL
                  WRITE(6,*) ' IADSTSTP(IABP), IADSTSTL(IABL) = ',
     &                         IADSTSTP(IABP), IADSTSTL(IABL) 
                 END IF
*. Should normal or AB transformed DET be used
                 I_DO_AB_SWITCH = 0
                 IF(PSSIGN.NE.0.0D0) THEN
                  IF(IABSWITCHP(IABP).EQ.-1.OR.
     &            (IABSWITCHP(IABP).EQ.0.AND.IABSWITCHL(IABL).EQ.-1))
     &            THEN
                     I_DO_AB_SWITCH = 1
                  END IF
                 END IF
                 IF(NTEST.GE.10000) WRITE(6,*)
     &           'I_DO_AB_SWITCH = ', I_DO_AB_SWITCH
*. Standard lexical address - independent of possible AB switch
                 ICONF_LEXN = 
     &           ABS(IADSTSTP(IABP)) + ABS(IADSTSTL(IABL)) - 1
                 IF(NTEST.GE.10000) WRITE(6,*) 
     &           ' ICONF_LEXN = ', ICONF_LEXN
                 ICNF_OUTN = ICONF_REO(ICONF_LEXN)
                 IF(NTEST.GE.10000) 
     &           WRITE(6,*) ' ICONF_LEXN,ICNF_OUTN = ', 
     &           ICONF_LEXN,ICNF_OUTN
                 IF(I_DO_AB_SWITCH.EQ.0) THEN
                  ISIGNN = SIGN(1,IADSTSTP(IABP))*SIGN(1,IADSTSTL(IABL))
                 ELSE
                  ISIGNN = 
     &            SIGN(1,IADSTSTP(IABP))*SIGN(1,IADSTSTL_AB(IABL))*
     &            IPSSIGN
                 END IF
                 IF(NTEST.GE.10000) THEN
                  WRITE(6,*) ' IADSTSTP(IABP),IADSTSTL(IABL),ISIGNN = ',
     &                         IADSTSTP(IABP),IADSTSTL(IABL),ISIGNN
                 END IF
*. Info on open orbitals
                 NOPL = NOPSTSTL(IABL)
                 NOPP = NOPSTSTP(IABP) 
                 NOPT = NOPP + NOPL
                 NPTDTN = NPCMCNF(NOPT+1)
                 IF(NTEST.GE.1000)
     &           WRITE(6,'(A,4I4)') ' NOPL, NOPP, NOPT, NPTDTN = ',
     &                                NOPL, NOPP, NOPT, NPTDTN
                 IBCNF_OUTN = IB_CN_OPEN(NOPT+1)
*
                 IF(NOPT.GE.MINOP) THEN
                  IF(I_DO_AB_SWITCH.EQ.0) THEN
                   NALPHAP = NALSTSTP(IABP) 
                   IF(NTEST.GE.10000) THEN
                     WRITE(6,*) ' NOPP, NALPHAP = ',
     &                            NOPP, NALPHAP
                   END IF
*. contribution to lexical address of proto-type det
                   IZP = IADOPSTSTP(IABP) 
                   NALPHAL = NALSTSTL(IABL)
*. There is a contribution to the sign from going from pa pb la lb TO pa la pb lb
                   IF(MOD(NALEL*NBPEL,2).NE.0) 
     &             ISIGNN = -ISIGNN
                   IPTDTN = IZNUM_PTDT2(ICNOCSTSTL((IABL-1)*NLEL+1),
     &             MAXOP,NOPT,NOPP+1,NALPHAP,IZP,NALPHAL,
     &             int_mb(KZ_PTDT(MAXOP+1)),int_mb(KREO_PTDT(NOPT+1)),1)
C                  IZNUM_PTDT2(IAB,NOPEN_DIM,NOPEN,IOPEN1,IALPHA1,IZ1,
C    &             NALPHA,Z,NEWORD,IREORD)
                  ELSE
*. WE have AB switch
                   NALPHAP = NALSTSTP(IABP) 
                   IF(NTEST.GE.10000) THEN
                     WRITE(6,*) ' NOPP, NALPHAP = ',
     &                            NOPP, NALPHAP
                   END IF
*. contribution to lexical address of proto-type det
                   IZP = IADOPSTSTP(IABP) 
                   NALPHAL = NALSTSTL_AB(IABL)
*. There is a contribution to the sign from going from pa pb la lb TO pa la pb lb
                   IF(MOD(NBLEL*NAPEL,2).NE.0) 
     &             ISIGNN = -ISIGNN
                   IPTDTN = IZNUM_PTDT2(ICNOCSTSTL_AB((IABL-1)*NLEL+1),
     &             MAXOP,NOPT,NOPP+1,NALPHAP,IZP,NALPHAL,
     &             int_mb(KZ_PTDT(MAXOP+1)),int_mb(KREO_PTDT(NOPT+1)),1)
                  END IF ! I_DO_AB_SWITCH
                  IF(NTEST.GE.10000) THEN
                    WRITE(6,*) ' New: '
                    WRITE(6,*) ' Number of det in list of PTDT ', IPTDTN
                    WRITE(6,*) ' IB_CM_OPEN(NOPT+1) = ',
     &                           IB_CM_OPEN(NOPT+1)
                    WRITE(6,*) ' ICNF_OUTN, NPTDTN ', ICNF_OUTN, NPTDTN
                    WRITE(6,*) ' IBCNF_OUTN = ', IBCNF_OUTN
                  END IF
*
                  IADR_SD_CONF_ORDERN = IB_CM_OPEN(NOPT+1) - 1
     &          + (ICNF_OUTN-IBCNF_OUTN)*NPTDTN+IPTDTN
*
                  IF(NTEST.GE.10000) THEN
                    WRITE(6,*) ' IDET, IADR_SD_CONF_ORDERN ',
     &                           IDET, IADR_SD_CONF_ORDERN
                  END IF
*. 
                  IREO(IADR_SD_CONF_ORDERN) = ISIGNN*IDET
                 ENDIF !number of open in range
*. End new
*
*
*. Old:
*
                IF(I_DO_ALSO_OLD .EQ. 1) THEN
* We now have an alpha and a beta-string, merge parts of alpha and beta
*. We now have a full alpha string, merge 
                 IF(NTEST.GE.1000) WRITE(6,*) ' >>> OLD '
                 IADR = 1+(IBAPSTR(IAPSM)+IAP-1-1)*NAPEL
                 CALL ICOPVE(IAPSTR(IADR),IASTRF,NAPEL)
                 IF(NTEST.GE.10000) WRITE(6,*) ' IADR, IAPSTR(1) = ',
     &           IADR, IAPSTR(1)
                 CALL ICOPVE2(int_mb(KOCSTR(IALTP)),
     &           (IBSTSGP(IALSM,IALTP)+IAL-1-1)*NALEL+1,
     &           NALEL,IASTRF(1+NAPEL))
                 IF(NTEST.GE.1000) THEN
                   WRITE(6,*) ' IA, IAP, IAL, IASTRF =', IA,IAP, IAL
                   CALL IWRTMA3(IASTRF,1,NAEL,1,NAEL)
                 END IF
                 IF(NTEST.GE.1000) THEN
                   WRITE(6,*) ' Off and IBPSTR(1) = ',
     &             1+(IBBPSTR(IBPSM)+IBP-1-1)*NBPEL, IBPSTR(1)
                 END IF
                 CALL ICOPVE(IBPSTR(1+(IBBPSTR(IBPSM)+IBP-1-1)*NBPEL),
     &                       IBSTRF,NBPEL)
C ICOPVE2(IIN,IOFF,NDIM,IOUT)
                 CALL ICOPVE2(int_mb(KOCSTR(IBLTP)),
     &                (IBSTSGP(IBLSM,IBLTP)+IBL-1-1)*NBLEL+1,
     &                NBLEL,IBSTRF(1+NBPEL))
                 IF(NTEST.GE.1000) THEN
                   WRITE(6,*) 
     &             ' IB, IBP, IBL, IBSTRF =', IB,IBP, IBL
                   CALL IWRTMA3(IBSTRF,1,NBEL,1,NBEL)
                 END IF
*
*. Reform from AB to CONF form
*
*. subtract offset(number of inactive orbitals)
                 ISUB = 1-IB_ORB
C                ABSTR_TO_CNFSTR(IA_OC,IB_OC,NAEL,NBEL,
C    &           ICONF,IOP_SP,ISIGN,IADD,NOP,NOC)
                 CALL ABSTR_TO_CNFSTR(IASTRF,IBSTRF,NAEL,NBEL,
     &                           IDET_OC,IDET_MS,ISIGN,
     &                           ISUB,NOP,NOP_AL,NOC,IOP1)
*
* On output IDET_OC is conf in compact form, IDET_MS is projection of 
* open orbitals
*
                 IF(NTEST.GE.200) THEN
                   WRITE(6,'(A,I10)') ' IDET =', IDET
                   WRITE(6,*) ' IASTR, IBSTR, ICONF:'
                   CALL IWRTMA(IASTRF(1),1,NAEL,1,NAEL)
                   CALL IWRTMA(IBSTRF(1),1,NBEL,1,NBEL)
                   CALL IWRTMA(IDET_OC,1,NOC,1,NOC)
                 END IF
* 
                 IF(NOP.GE.MINOP) THEN
                  NPTDT = NPCMCNF(NOP+1)
*. Address of this configuration 
                  ICNF_OUT = ILEX_FOR_CONF2(IDET_OC,NOC,NACOB,NEL,
     &                       IZCONF,1,ICONF_REO,1)
C                            ILEX_FOR_CONF2(ICONF,NOCC_ORB,NORB,NEL,IARCW,
C                            IDOREO,IREO,IB_OCCLS)
                  IF(NTEST.GE.10000) THEN
                    WRITE(6,*) ' Configuration: '
                    CALL IWRTMA(IDET_OC,1,NOC,1,NOC)
                    WRITE(6,*) 
     &              ' Address of configuration in output list',ICNF_OUT
                  END IF
*. Address of spinprojection pattern   
C                         IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
                  IPTDT = IZNUM_PTDT(IDET_MS,NOP,NOP_AL,
     &                    int_mb(KZ_PTDT(NOP+1)),
     &                    int_mb(KREO_PTDT(NOP+1)),1)
                  ISIGNC = 1
                  IF(IPTDT.EQ.0) THEN
                   IF(PSSIGN.NE.0) THEN
*. The determinant was not found among the list of prototype dets. 
*. For combinations this should be due to the prototype determinant 
*  is the MS- switched determinant, so find
*. address of this and remember sign
                     M1 = -1
                     CALL ABSTR_TO_CNFSTR(
     &                    IBSTRF,IASTRF,NBEL,NAEL,
     &                    IDET_OC,IDET_MS,ISIGN,
     &                    ISUB,NOP,NOP_AL,NOC,IOP1)
                     IPTDT = IZNUM_PTDT(IDET_MS,NOP,NOP_AL,
     &                       int_mb(KZ_PTDT(NOP+1)),
     &                       int_mb(KREO_PTDT(NOP+1)),1)
                     IF(PSSIGN.EQ.-1.0D0) ISIGNC = -1
                   ELSE 
*. Prototype determinant was not found in list
                    WRITE(6,*) 
     &              ' Error: Determinant not found in list of protodets'
                    WRITE(6,*) 
     &              ' Detected in REO_SD_FOR_OCCLS_S'
                   END IF
                  END IF
*
                  IBCNF_OUT = IB_CN_OPEN(NOP+1)
*
                  IF(NTEST.GE.10000) THEN
                   WRITE(6,*) ' Number of det in list of PTDT ', IPTDT
                   WRITE(6,*) ' IB_CM_OPEN(NOP+1) = ',
     &                          IB_CM_OPEN(NOP+1)
                   WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
                   WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
                  END IF
*
                  IADR_SD_CONF_ORDER = IB_CM_OPEN(NOP+1) - 1
     &                               + (ICNF_OUT-IBCNF_OUT)*NPTDT+IPTDT
*. Consistency, old and new
                  IF(IADR_SD_CONF_ORDER.NE.IADR_SD_CONF_ORDERN.OR.
     &              ISIGN.NE.ISIGNN) THEN
                    WRITE(6,'(A,3I6)')
     &              ' CMP:IDET,IADR_SD_CONF_ORDER,IADR_SD_CONF_ORDERN:',
     &                    IDET,IADR_SD_CONF_ORDER,IADR_SD_CONF_ORDERN 
                    WRITE(6,'(A,2I2)')
     &              ' CMP:ISIGN, ISIGNN = ', ISIGN, ISIGNN
*                 
                    WRITE(6,'(A,2I6)')
     &              ' CMP:ICNF_OUT, ICNF_OUTN =', ICNF_OUT, ICNF_OUTN
                    WRITE(6,'(A,2I6)')
     &              ' CMP:IPTDT, IPTDTN =', IPTDT, IPTDTN
                    WRITE(6,'(A,2I6)')
     &              ' CMP:IBCNF_OUT, IBCNF_OUTN = ',
     &                     IBCNF_OUT, IBCNF_OUTN 
                    WRITE(6,*) ' NOPP, NOPL, NOP = ',
     &                           NOPP, NOPL, NOP
                   WRITE(6,'(A,4I4)') ' IAL, IBL, IAP, IBP = ',
     &                                  IAL, IBL, IAP, IBP
                   WRITE(6,'(A,4I4)') ' IALSM, IBLSM, IAPSM, IBPSM = ',
     &                                  IALSM, IBLSM, IAPSM, IBPSM
                   WRITE(6,*) ' (IBP-1)*NAP+IAP = ', (IBP-1)*NAP+IAP
                   WRITE(6,*) 'IB_SMSMP(IAPSM,IBPSM) = ', 
     %                         IB_SMSMP(IAPSM,IBPSM)
                   WRITE(6,*) ' NAP, NBP = ', NAP, NBP
                   WRITE(6,*) ' IABP = ', IABP
                   WRITE(6,*) ' IASTR, IBSTR (from OLD):'
                   CALL IWRTMA(IASTRF(1),1,NAEL,1,NAEL)
                   CALL IWRTMA(IBSTRF(1),1,NBEL,1,NBEL)
*
                    WRITE(6,*) ' Configuration (from OLD): '
                    CALL IWRTMA(IDET_OC,1,NOC,1,NOC)
                    WRITE(6,*) 'PROBLEM IN SSN reorder '
                    CALL MEMCHK2('SSNPRO')
                    STOP 'PROBLEM IN SSN reorder'
                  END IF
*
                  IF(IADR_SD_CONF_ORDER.LE.0) THEN
                    WRITE(6,*) 
     &              ' Problemo, IADR_SD_CONF_ORDER < 0 '
                    WRITE(6,*) 
     &              ' IADR_SD_CONF_ORDER = ', IADR_SD_CONF_ORDER
                    WRITE(6,*) 
     &              ' Number of det in list of PTDT ', IPTDT
                    WRITE(6,*) 
     &              ' IB_CM_OPEN(NOP+1) = ', IB_CM_OPEN(NOP+1)
                    WRITE(6,*) ' ICNF_OUT, NPTDT ', ICNF_OUT, NPTDT
                    WRITE(6,*) ' IBCNF_OUT = ', IBCNF_OUT
C?                  CALL XFLUSH(6)
                  END IF
                  IF(NTEST.GE.10000) THEN
                    WRITE(6,*) ' IADR_SD_CONF_ORDER, ISIGN, IDET = ',
     &                           IADR_SD_CONF_ORDER, ISIGN, IDET
                  END IF
C                 IREO(IADR_SD_CONF_ORDER) = ISIGN*IDET*ISIGNC
                  IF(NTEST.GE.10000) THEN
                    WRITE(6,*) ' IDET, IADR_SD_CONF_ORDER ',
     &                           IDET, IADR_SD_CONF_ORDER
                  END IF
                 END IF! Nop .ge. MINOP
                 END IF! I_DO_ALSO_OLD
                END IF ! IA geq MINIA
               END DO ! IAL
              END DO ! IAP
             END DO ! IALSM
            END DO ! IBL
           END DO ! IBP
         END DO! IBLSM
       END IF! AB block is active
  999  CONTINUE
       END DO ! IASM
 9999  CONTINUE
      END DO! Loop over AB blocks in occ. class
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' Reorder array for SDs, CONF order => string order'
        WRITE(6,*) ' ================================================='
        CALL IWRTMA(IREO,1,NCMB_CN,1,NCMB_CN)
      END IF
*
*
      I_DO_CHECKSUM = 0
      IF(I_DO_CHECKSUM.EQ.1) THEN
        ISUM = 0
        DO JDET = 1, IDET
          ISUM = ISUM + ABS(IREO(JDET))
        END DO
        IF(ISUM.NE.IDET*(IDET+1)/2) THEN
          WRITE(6,*) ' Problem with sumcheck in REO_GASDET'
          WRITE(6,'(A,2I9)') 
     &    'Expected and actual value ', ISUM, IDET*(IDET+1)/2
          STOP       ' Problem with sumcheck in REO_GASDET'
        ELSE
          WRITE(6,*) ' Sumcheck in REO_GASDET passed '
        END IF
      END IF !checksum is invoked
*
      RETURN
      END
      SUBROUTINE LAST_OCCGAS_FOR_SUPGP(IGPSPGP,NGAS,LASTOCC,NLASTEL)
*
* A supergroup is given in the form of groups in IGPSPGP
* Find the last orbital space with a nonvanishing number
* of electrons and the number of electrons in this space
*
*. Jeppe Olsen, Apr. 30, 2013; Minnepolis
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
*. Input
      INTEGER IGPSPGP(NGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Info from LAST_OCCGAS_FOR_SUPGP '
       WRITE(6,*) ' ================================'
       WRITE(6,*)
      END IF
*
      LASTOCC = -1
      DO IGAS = 1, NGAS
        IF(NELFGP(IGPSPGP(IGAS)).GT.0) LASTOCC = IGAS
      END DO
      NLASTEL = NELFGP(IGPSPGP(LASTOCC))
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Supergroup, as groups: '
       CALL IWRTMA3(IGPSPGP,1,NGAS,1,NGAS)
       WRITE(6,*) ' Large nontrivial active space: ', LASTOCC
       WRITE(6,*) ' Number of electrons in this space: ', NLASTEL
      END IF
*
      RETURN
      END
      SUBROUTINE GET_DIM_PLSTRINGS
*
* Obtain some information about dimensions for PL approach
* (PL : Strings are splitted in two parts, so last GAspace with
* nonvanishing number of electrons is separated out
*
* The allowed combinations of strings are defined from 
* K*ABSPGP_FOR_OCCLS
*
*. Output is  MXNSTRP,MXNSTRP_AS,MXNSTRL,
*    &        MXNSTRL_AS,MXNSTRSTRP_AS,MXNSTRSTRL_AS 
*
* Information is returned in GASSTR
*
*
*. Jeppe Olsen, May 7, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*
      IATP = 1
      IBTP = 2
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      IB_A = IBSPGPFTP(IATP)
      IB_B = IBSPGPFTP(IBTP)
*
      CALL GET_DIM_PLSTRING_S(
     &     NOCCLS_MAX,int_mb(KNABSPGP_FOR_OCCLS),
     &     int_mb(KIABSPGP_FOR_OCCLS),
     &     NELFSPGP(1,IB_A),NELFSPGP(1,IB_B))
*
      RETURN
      END
      SUBROUTINE GET_DIM_PLSTRING_S(
     &     NOCCLS_MAX,NABSPGP_FOR_OCCLS,
     &     IABSPGP_FOR_OCCLS,
     &     NELFSPGP_A,NELFSPGP_B)
*
* Inner routine for determing various dimensions for PL approach
*
*. Jeppe Olsen, May 7, 2013
*
*
* Output is returned in GASSTR
*. Output:
*. ======
*
* MXNSTRP: Max number of Pstrings for given type and single sym
* MXNSTRP_AS: Max number of Pstrings for given type and all sym
* MXNSTRL: Max number of Lstrings for given type and single sym
* MXNSTRSTRP_AS: Max number of combinations of P alpha and P beta
*                strings for all sym
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
*
*. Input
*
      INTEGER NABSPGP_FOR_OCCLS(NOCCLS_MAX)
      INTEGER IABSPGP_FOR_OCCLS(2,NOCCLS_MAX)
      INTEGER NELFSPGP_A(MXPNGAS,*), NELFSPGP_B(MXPNGAS,*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_DIM_PLSTRING_S' 
        WRITE(6,*) ' ============================'
      END IF
*
      MXNSTRP = 0
      MXNSTRP_AS = 0
      MXNSTRSTRP_AS = 0
      MXNSTRL = 0
      MXNSTRL_AS = 0
      MXNSTRSTRL_AS = 0
      MXNSTRSTRLOC_AS = 0
*
      IATP = 1
      IBTP = 2
*
      IB_A = IBSPGPFTP(IATP)
      IB_B = IBSPGPFTP(IBTP)
*
      IABTP = 0
      DO IOCCLS = 1, NOCCLS_MAX
       NABTP = NABSPGP_FOR_OCCLS(IOCCLS)
       IF(NTEST.GE.1000) WRITE(6,*) ' IOCCLS, NABTP = ', IOCCLS, NABTP
       DO IIABTP = 1, NABTP
        IABTP = IABTP + 1
        IATP = IABSPGP_FOR_OCCLS(1,IABTP) + IB_A - 1
        IBTP = IABSPGP_FOR_OCCLS(2,IABTP) + IB_B - 1
        CALL DIM_ABPL(ISPGPFTP(1,IATP),ISPGPFTP(1,IBTP),
     &           LMXNSTP, LMXNSTP_AS,LMXNSTL,LMXNSTL_AS,
     &           LNSTSTP_AS, LNSTSTL_AS,LNSTSTLOC_AS)
*
        MXNSTRP = MAX(MXNSTRP,LMXNSTP)
        MXNSTRP_AS  = MAX(MXNSTRP_AS,LMXNSTP_AS)
        MXNSTRL = MAX(MXNSTRL,LMXNSTL)
        MXNSTRL_AS  = MAX(MXNSTRL_AS,LMXNSTL_AS)
        MXNSTRSTRP_AS = MAX(MXNSTRSTRP_AS, LNSTSTP_AS)
        MXNSTRSTRL_AS = MAX(MXNSTRSTRL_AS, LNSTSTL_AS)
        MXNSTRSTRLOC_AS  = MAX(MXNSTRSTRLOC_AS, LNSTSTLOC_AS)
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,'(A,2I6)')' MXNSTRP, MXNSTRP_AS  = ', MXNSTRP, MXNSTRP_AS
       WRITE(6,'(A,2I6)')' MXNSTRL, MXNSTRL_AS  = ', MXNSTRL, MXNSTRL_AS
       WRITE(6,'(A,2I6)')' MXNSTRSTRP_AS, MXNSTRSTRL_AS  = ',
     &                     MXNSTRSTRP_AS, MXNSTRSTRL_AS
       WRITE(6,'(A,I7)')' MXNSTRSTRLOC_AS  = ',  MXNSTRSTRLOC_AS
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_ABPL(ISPGPA,ISPGPB,
     &           MXNSTP, MXNSTP_AS,MXNSTL,MXNSTL_AS,
     &           NSTSTP_AS, NSTSTL_AS,NSTSTLOC_AS)

*
*. Alpha- and beta-supergroups are given in ISPGPA, ISPGPB
*. Determine dimension of P-strings and products of these
*
*. a P-string is a string where the last occupied space have been
*. removed
*
*. Jeppe Olsen, May 7, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Input
      INTEGER ISPGPA(NGAS),ISPGPB(NGAS)
*. Local scratch
      INTEGER NST_AP(MXPNSMST),NST_AL(MXPNSMST)
      INTEGER NST_BP(MXPNSMST),NST_BL(MXPNSMST)
      INTEGER IDUM(MXPNSMST)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from DIM_ABPL '
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input a- and b-groups '
        CALL IWRTMA(ISPGPA,NGAS,1,NGAS,1)
        CALL IWRTMA(ISPGPB,NGAS,1,NGAS,1)
      END IF
*
      IZERO = 0
* 
*. Find the last occupied space of the alpha- and beta-strings
C    LAST_OCCGAS_FOR_SUPGP(IGPSPGP,NGAS,LASTOCC,NLASTEL)
      CALL LAST_OCCGAS_FOR_SUPGP(ISPGPA,NGAS,LASTOCCA,NLASTELA)
      CALL LAST_OCCGAS_FOR_SUPGP(ISPGPB,NGAS,LASTOCCB,NLASTELB)
*. The split in gaspaces should be identical for a and b
      LASTOCC = MAX(LASTOCCA, LASTOCCB)
      LASTOCCA = LASTOCC
      LASTOCCB = LASTOCC
*
      NGAS_AP = LASTOCCA -1
      NGAS_BP = LASTOCCB -1
      
      IF(NGAS_AP.LE.0) THEN
        NGAS_AP = 1
      END IF
      IF(NGAS_BP.LE.0) THEN
        NGAS_BP = 1
      END IF
*
      NGAS_AL = NGAS - NGAS_AP
      NGAS_BL = NGAS - NGAS_BP
      IF(NTEST.GE.1000) THEN
        WRITE(6,'(A,4I4)') ' NGAS_AP, NGAS_BP, NGAS_AL, NGAS_BL = ',
     &                       NGAS_AP, NGAS_BP, NGAS_AL, NGAS_BL
      END IF
*. Dimension of the substrings
C          GET_DIM_GNSPGP(NGRPA,IGRPA,NSTPSM,IBSTPSM)
      IF(NTEST.GE.1000) WRITE(6,*) ' AP: '
      CALL GET_DIM_GNSPGP(NGAS_AP,ISPGPA,NST_AP,IDUM)
      IF(NTEST.GE.1000) WRITE(6,*) ' BP: '
      CALL GET_DIM_GNSPGP(NGAS_BP,ISPGPB,NST_BP,IDUM)
      IF(NGAS_AL.EQ.0) THEN
        CALL ISETVC(NST_AL,IZERO,NSMST)
        NST_AL(1) = 1
      ELSE
        IF(NTEST.GE.1000) WRITE(6,*) ' AL: '
        CALL GET_DIM_GNSPGP(NGAS_AL,ISPGPA(NGAS_AP+1),NST_AL,IDUM)
      END IF
      IF(NGAS_BL.EQ.0) THEN
        CALL ISETVC(NST_BL,IZERO,NSMST)
        NST_BL(1) = 1
      ELSE
        IF(NTEST.GE.1000) WRITE(6,*) ' BL: '
        CALL GET_DIM_GNSPGP(NGAS_BL,ISPGPB(NGAS_BP+1),NST_BL,IDUM)
      END IF
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) 'NST_AL, NST_BL '
       CALL IWRTMA(NST_AL, 1, NSMST,1,NSMST)
       CALL IWRTMA(NST_BL, 1, NSMST,1,NSMST)
      END IF
*. Max and sum
      MAX_AP = IMNMX(NST_AP,NSMST,2)
      MAX_BP = IMNMX(NST_BP,NSMST,2)
      MXNSTP = MAX(MAX_AP,NSMST,2)
      MAX_AL = IMNMX(NST_AL,NSMST,2)
      MAX_BL = IMNMX(NST_BL,NSMST,2)
      MXNSTL = MAX(MAX_AL,NSMST)
*
      NSTAP_AS = IELSUM(NST_AP, NSMST)
      NSTBP_AS = IELSUM(NST_BP, NSMST)
      NSTAL_AS = IELSUM(NST_AL, NSMST)
      NSTBL_AS = IELSUM(NST_BL, NSMST)
*
      MXNSTP_AS = MAX(NSTAP_AS,NSTBP_AS)
      MXNSTL_AS = MAX(NSTAL_AS,NSTBL_AS)
C              INT_PROD(I1,I2,NDIM)
      NSTSTP_AS = NSTAP_AS*NSTBP_AS
      NSTSTL_AS = NSTAL_AS*NSTBL_AS
      NELL = NLASTELA + NLASTELB
      NSTSTLOC_AS = NSTAL_AS*NSTBL_AS*NELL
*
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,4I6)') ' NSTAP_AS, NSTBP_AS, NSTAL_AS, NSTBL_AS =',
     &                       NSTAP_AS, NSTBP_AS, NSTAL_AS, NSTBL_AS 
        WRITE(6,'(A,2I6)') ' MXNSTP, MXNSTL = ', MXNSTP, MXNSTL
        WRITE(6,'(A,2I6)') 
     &  ' NSTSTP_AS, NSTSTL_AS = ', NSTSTP_AS, NSTSTL_AS
        WRITE(6,'(A,I6)') 
     &  ' NSTSTLOC_AS = ', NSTSTLOC_AS
      END IF
*
      RETURN
      END

      
     
      
      

      
C
C     NSTPTP_GAS_NEW(NGAS,ISPGRP,NSTSGP,NSMST,
C    &                      NSTSSPGP,IGRP,MXNSTR,
C    &                      NSMCLS,NSMCLSE,NSMCLSE1,NSTR_AS)

*
      FUNCTION NSTR_SPGP_ALLSYM(NGRP,NELGRP,NOBPGP)
*
* A superstring spanned by  NGRP groups with number of electrons
* and orbitals given by NELGRP, NOBGRP, respectively. Obtain
* number of strings, all sym
*
*. Jeppe Olsen, May 7, 2012
*
      INCLUDE 'implicit.inc'
*. Input:
      INTEGER NELGRP(NGRP), NOBGRP(NGRP)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from NSTR_SPGP_ALLSYM:'
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NOBGRP and NELGRP: '
        CALL IWRTMA(NOBGRP,1,NGRP,1,NGRP)
        CALL IWRTMA(NELGRP,1,NGRP,1,NGRP)
      END IF
*
      NSTR = 1
      DO IGRP = 1, NGRP
       NSTR = NSTR*IBION(NOBGRP(IGRP),NELGRP(IGRP))
      END DO
*
      NSTR_SPGP_ALLSYM = NSTR
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Number of strings obtained ', NSTR
      END IF
*
      RETURN
      END
      FUNCTION INT_PROD(I1,I2,NDIM)
*
*. Product of two integer arrays I1, I2
*
*. Jeppe Olsen, May 7, 2013
*
      INCLUDE 'implicit.inc'
*
      INTEGER I1(NDIM),I2(NDIM)
*
      IPROD = 0
      DO I = 1, NDIM
       IPROD = IPROD + I1(I)*I2(I)
      END DO
*
      INT_PROD = IPROD
*
      RETURN
      END 
      SUBROUTINE IWRTMATMAT(IAA,NBLR,NBLC,LBLR,LBLC,IB)
*
* A matrix of integer matrices is given in AA
*
*. Print!!
*
*. Printing dimension matrices running over
*. all symmetry blocks
*
*. Jeppe Olsen, May 14, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IAA(*)
      INTEGER LBLR(NBLR),LBLC(NBLC)
      INTEGER IB(NBLR,NBLC)
*
      WRITE(6,*) ' Blocked integer matrix: '
*
      DO IBLR = 1, NBLR
       DO IBLC = 1, NBLC
        NR = LBLR(IBLR)
        NC = LBLC(IBLC)
        IBB = IB(IBLR,IBLC)
*
        IF(NR*NC.NE.0) THEN
          WRITE(6,'(A,2I3)') 'Row and column block: ', IBLR,IBLC
          CALL IWRTMA(IAA(IBB),NR,NC,NR,NC)
        END IF
*
       END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GEN_NSDAB_FOR_ALL_OCCLS(
     &           NABSPGP_PER_OCCLS,IABSPGP_PER_OCCLS,
     &           N_SDAB_PER_OCCLS,N_CMAB_PER_OCCLS,
     &           N_SDAB_PER_OCCLS_MAX,N_CMAB_PER_OCCLS_MAX)
*
*. Determine:
*   Number of AB determinants and combinations per occupation class
*   Max of these
*
*. Note: AB determinants are the determinants of the string expansions
*        of CI-spaces. The number of these may differ from the number
*        of determinants obtained from the configurations, when MS < S.
*
*. Jeppe Olsen, May 17, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
*. Input
      INTEGER NABSPGP_PER_OCCLS(NOCCLS_MAX)
      INTEGER IABSPGP_PER_OCCLS(2,NOCCLS_MAX)
*. Output
      INTEGER N_SDAB_PER_OCCLS(NSMST,NOCCLS_MAX)
      INTEGER N_CMAB_PER_OCCLS(NSMST,NOCCLS_MAX)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Info from GEN_NSDAB_FOR_ALL_OCCLS '
        WRITE(6,*) ' ==================================='
        WRITE(6,*) 
      END IF
       
*
      IB_OCCLS = 1
      DO IOCCLS = 1, NOCCLS_MAX
        CALL NSDAB_FOR_OCCLS(NABSPGP_PER_OCCLS(IOCCLS),
     &                       IABSPGP_PER_OCCLS(1,IB_OCCLS),
     &                       N_SDAB_PER_OCCLS(1,IOCCLS),
     &                       N_CMAB_PER_OCCLS(1,IOCCLS) )
        IB_OCCLS = IB_OCCLS + NABSPGP_PER_OCCLS(IOCCLS)
      END DO
*
      N_SDAB_PER_OCCLS_MAX = IMNMX(N_SDAB_PER_OCCLS,
     &                       NSMST*NOCCLS_MAX,2)
      N_CMAB_PER_OCCLS_MAX = IMNMX(N_CMAB_PER_OCCLS,
     &                       NSMST*NOCCLS_MAX,2)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of SDABs per sym(row) and occlass(col):'
        WRITE(6,*)
        CALL IWRTMA(N_SDAB_PER_OCCLS,NSMST,NOCCLS_MAX,NSMST,NOCCLS_MAX)
        WRITE(6,*) ' Number of CMABs per sym(row) and occlass(col):'
*
        CALL IWRTMA(N_CMAB_PER_OCCLS,NSMST,NOCCLS_MAX,NSMST,NOCCLS_MAX)
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Largest number of AB dets per occlass ', N_SDAB_PER_OCCLS_MAX
        WRITE(6,*) 
     &  ' Largest number of AB cmbs per occlass ', N_CMAB_PER_OCCLS_MAX
      END IF
*
      RETURN
      END
      SUBROUTINE NSDAB_FOR_OCCLS(NABSPGP,IABSPGP,
     &           N_SDAB_FOR_OCCLS,N_CMAB_FOR_OCCLS)
*
*.  Obtain number of AB determinants for an occupation class
*
*. AB determinants: The determinants written in terms of 
*. strings. Does not depend on the number of min number of open shells, 
*. and differs therefore from the numbers obtained from configurations
*. (which are constructed to respect a min number of orbitals)
*.
*. Jeppe Olsen; May 16, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'wrkspc-static.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER IABSPGP(2,*)
*. Local scratch
      INTEGER IBLTP(MXPCSM)
*. Output
      INTEGER N_SDAB_FOR_OCCLS(NSMCI), N_CMAB_FOR_OCCLS(NSMCI)

*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NAB_CMB_FOR_OCCLS speaking'
        WRITE(6,*) ' =========================='
      END IF
*
      IATP = 1
      IBTP = 2
*
      NOCTPA =  NOCTYP(IATP)
      NOCTPB =  NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      CALL SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
*
      DO ISM = 1, NSMCI
*
        IF(NTEST.GE.1000) WRITE(6,*) ' ISM = ', ISM
        NSD = 0
        NCM = 0
*
        IDUM = 1
        CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,IBLTP,IDUM)
        IF(NTEST.GE.1000) WRITE(6,*) ' NABSPGP = ', NABSPGP
    
        DO JABSPGP = 1, NABSPGP
          JASPGP = IABSPGP(1,JABSPGP)
          JBSPGP = IABSPGP(2,JABSPGP)
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' JABSPGP, JASPGP, JBSPGP = ',
     &                 JABSPGP, JASPGP, JBSPGP
          CALL DIM_ABSPGP(JASPGP,JBSPGP,NSMST,
     &         int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &         IBLTP,ISMOST(1,ISM),NSDL,NCML)
C     DIM_ABSPGP(IASPGP,IBSPGP,NSMST,
C    &           NSSOA,NSSOB,IBLTP,ISMOST,NSDL,NCML)
          NSD = NSD + NSDL
          NCM = NCM + NCML
        END DO
        N_SDAB_FOR_OCCLS(ISM) = NSD
        N_CMAB_FOR_OCCLS(ISM) = NCM
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of SD and CM per symmetry for occ class '
        CALL IWRTMA(N_SDAB_FOR_OCCLS,1,NSMCI,1,NSMCI)
        WRITE(6,*)
        CALL IWRTMA(N_CMAB_FOR_OCCLS,1,NSMCI,1,NSMCI)
      END IF
*
      RETURN
      END
      SUBROUTINE DIM_ABSPGP(IASPGP,IBSPGP,NSMST,
     &           NSSOA,NSSOB,IBLTP,ISMOST,NSDL,NCML)
*
* Dimension of supergroup IASPGP, IBSPGP, ISM
*
      INCLUDE 'implicit.inc'
*. General input
      INTEGER NSSOA(NSMST,*), NSSOB(NSMST,*)
      INTEGER IBLTP(*), ISMOST(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from DIM_ABSPGP'
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' NSMST = ',  NSMST
        WRITE(6,*) ' ISMOST array '
        CALL IWRTMA3(ISMOST,1,NSMST,1,NSMST)
      END IF
      
      NSDL = 0
      NCML = 0
*
      ISYM = 0
      DO IASM = 1, NSMST
*
        IBSM = ISMOST(IASM)
        LASTR = NSSOA(IASM,IASPGP)
        LBSTR = NSSOB(IBSM,IBSPGP)
        NSDL = NSDL + LASTR*LBSTR
*
        IF(IBLTP(IASM).NE.0) THEN
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' IASM, IBSM = ', IASM, IBSM
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' IASPGP, IBSPGP,IBLTP(IASM) = ', 
     &                   IASPGP, IBSPGP,IBLTP(IASM)
            WRITE(6,*) ' LASTR, LBSTR = ', LASTR, LBSTR
          END IF
          IF(IBLTP(IASM).EQ.1) THEN
            ISYM = 0
          ELSE IF (IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          END IF
          IF(.NOT.(ISYM.EQ.1.AND.IASPGP.LT.IBSPGP)) THEN
            IF(ISYM.EQ.0) THEN
              NCML = NCML + LASTR*LBSTR
            ELSE
              IF(IASPGP.EQ.IBSPGP) THEN
               NCML = NCML + LASTR*(LASTR+1)/2
              ELSE
               NCML = NCML + LASTR*LBSTR
              END IF 
            END IF! ISYM = 0
          END IF! TTSS Block allowed
        END IF ! SS block allowed
      END DO! Loop over IASM
*
      IF(NTEST.GE.1000)
     & WRITE(6,*) ' NSDL, NCML = ', NSDL, NCML
*
      RETURN
      END
      SUBROUTINE EXP_CNFSPC(CIVECIN,CIVECUT,ICONF_OCC,NCONF_FOR_OPEN,
     &           MINOCC_IN,MAXOCC_IN,NOBCNF)
*
* An Output CI expansion is defined by ICONF_OCC, NCONF_FOR_OPEN
*
* An input expansion is defined by MINOCC_IN, MAXOCC_IN
*
* Expand input to output - at the moment it is assumed that output contains
* all input csf's
*
*. Jeppe Olsen, May 2013
*
      INCLUDE 'implicit.inc'
      REAL*8 INPROD
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input
      DIMENSION ICONF_OCC(*),NCONF_FOR_OPEN(*)
      DIMENSION CIVECIN(*)
*. Output
      DIMENSION CIVECUT(*)
*. Local scratch
      DIMENSION ICONFA(MXPORB), ICONFB(MXPORB)
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from EXP_CNFSPC '
        WRITE(6,*) ' ====================== '
      END IF
*
* The number of CSFs in OUT
*
      NCSF = 0
      DO IOPEN = 0, MAXOP
        NNCNF = NCONF_FOR_OPEN(IOPEN+1)
        NNCSF = NPCSCNF(IOPEN+1)
        NCSF = NCSF + NNCSF*NNCNF
      END DO
      IF(NTEST.GE.10) WRITE(6,*) ' NCSF = ', NCSF
*
      ZERO = 0.0D0
      CALL SETVEC(CIVECUT,ZERO,NCSF)
*
*. Loop over configurations for output expansion
      ICNF = 0
      ICSF = 0
      IB_IN = 1
      IB_UT = 1
*
      DO IOPEN = 0, MAXOP
        ITYP = IOPEN + 1
        ICL = (NACTEL - IOPEN) / 2
        IOCC = IOPEN + ICL
        IF( ITYP .EQ. 1 ) THEN
          ICNBS0 = 1
          IPBAS = 1
        ELSE
          ICNBS0 = ICNBS0 + NCONF_FOR_OPEN(IOPEN+1-1)*(NACTEL+IOPEN-1)/2
          IPBAS = IPBAS + NPCSCNF(IOPEN+1-1)*(IOPEN-1)
        END IF
*. Configurations of this type
        NNCNF = NCONF_FOR_OPEN(IOPEN+1)
        NNCSF = NPCSCNF(IOPEN+1)
        DO IC = 1, NNCNF
          ICNF = ICNF + 1
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
*. The configuration in action is ICONF_OCC(ICNBS)
*. Reorder first back to the 
C              REO_OB_CONFP(ICONFP_IN, ICONFP_UT,IREO,NOC)
          CALL REO_OB_CONFP(ICONF_OCC(ICNBS),ICONFB,IREO_MNMX_OB_ON,
     &                      IOCC)
*. Obtain accumulated form of 
C         PACK_TO_ACC_CONF(ICONFP,ICONFA,NOB,NOCOB,IWAY)
          CALL PACK_TO_ACC_CONF(ICONFB,ICONFA,NOBCNF,IOCC,1)
*. Check if configuration is in MINMAX space
C         IS_IACC_CONF_IN_MINMAX_SPC(IOCC,MIN_OCC,MAX_OCC,NORB)
          IM_IN = IS_IACC_CONF_IN_MINMAX_SPC(
     &            ICONFA,MINOCC_IN,MAXOCC_IN,
     &            NOBCNF)
          IF(IM_IN.EQ.1) THEN
            CALL COPVEC(CIVECIN(IB_IN),CIVECUT(IB_UT),NNCSF)
            IB_IN = IB_IN + NNCSF
            IB_UT = IB_UT + NNCSF
          ELSE
            IB_UT = IB_UT + NNCSF
          END IF
        END DO ! Loop over configurations in UT
      END DO ! Loop over NOPEN
*
*. Test: Norm of input and output vectors
*
      NCSF_IN = IB_IN - 1
      XNORM_IN = INPROD(CIVECIN,CIVECIN,NCSF_IN)
      XNORM_UT = INPROD(CIVECUT,CIVECUT,NCSF)
      IF(NTEST.GE.10) THEN
        WRITE(6,'(A,2E22.15)') 
     &  ' Norm of input and output vector ', XNORM_IN, XNORM_UT
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Obtained input and output vectors:'
        CALL WRTMAT(CIVECIN,1,NCSF_IN,1,NCSF_IN)
        WRITE(6,*)
        CALL WRTMAT(CIVECUT,1,NCSF,1,NCSF)
      END IF
C
      RETURN
      END
      SUBROUTINE PACK_TO_ACC_CONF(ICONFP,ICONFA,NOB,NOCOB,IWAY)
*
* Reform configuration between packed and accumulated form
* of configuration
*
*. IWAY = 1: Packed to accumulated 
*. IWAY = 2: Accumulated to packed 
*
*. Jeppe Olsen, May 2013
*
      INCLUDE 'implicit.inc'
*. Input or output
      INTEGER ICONFP(*), ICONFA(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from PACK_TO_ACC_CONF '
        WRITE(6,*) ' IWAY = ', IWAY
      END IF
*
      IF(IWAY.EQ.1) THEN
*. Packed => Accumulated
       NEL = 0
       IOC = 0
       DO JOB = 1, NOB
         IF(IOC.LE.NOCOB) THEN
           IF(ABS(ICONFP(IOC+1)).EQ.JOB) THEN
             IF(ICONFP(IOC+1).EQ.JOB) THEN
               NEL = NEL + 1
             ELSE
               NEL = NEL + 2
             END IF
             IOC = IOC + 1
           END IF
         END IF
         ICONFA(JOB) = NEL
       END DO
      ELSE
*. Accumulated => Packed
       STOP 'IWAY = 2 not Programmed in PACK_TO_ACC_CONF '
      END IF
*
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' ICONFP(input) and ICONFA(output) '
        ELSE IF (IWAY.EQ.2) THEN
          WRITE(6,*) ' ICONFP(output) and ICONFA(input) '
        END IF
        CALL IWRTMA(ICONFP,1,NOCOB,1,NOCOB)
        WRITE(6,*)
        CALL IWRTMA(ICONFA,1,NOB,1,NOB)
      END IF
*
      RETURN
      END
      SUBROUTINE REO_OB_CONFE(ICONFE_IN, ICONFE_UT,IREO_NO,NOB)
*
* Reorder occupations of configuration  according to orbital reordering IREO
*
* Input and output configurations are given in expanded form (occupation numbers)
*
*. Jeppe Olsen, May 30, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICONFE_IN(NOB), IREO_NO(NOB)
*. Output
      INTEGER ICONFE_UT(NOB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from REO_OB_CONFE '
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' IREO_NO array: '
        CALL IWRTMA(IREO_NO,1,NOB,1,NOB)
      END IF
*
      IZERO = 0
      CALL ISETVC(ICONFE_UT,IZERO,NOB)
*
       
      NOC_UT = 0
      DO IOB_IN = 1, NOB
        ICONFE_UT(IREO_NO(IOB_IN)) = ICONFE_IN(IOB_IN) 
      END DO
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Input and output configurations '
         CALL IWRTMA(ICONFE_IN,1,NOB,1,NOB)
         CALL IWRTMA(ICONFE_UT,1,NOB,1,NOB)
      END IF
*
      RETURN 
      END
      SUBROUTINE REO_OB_CONFP(ICONFP_IN, ICONFP_UT,IREO,NOC)
*
* Reorder occupations of configuration  according to orbital reordering IREO
* IREO gives new number for given old (input) number
*
* Input and output configurations are given in packed form (occupation numbers)
*. Dirty NOC*NOB scaling
*
*. Jeppe Olsen, May 30, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ICONFP_IN(NOC), IREO(*)
*. Output
      INTEGER ICONFP_UT(NOC)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from REO_OB_CONFP'
      END IF
*
*. Highest occupied orbital in in re-ordering
      JMAX = IREO(ABS(ICONFP_IN(1)))
      DO JOC = 2, NOC
       JOB_UT = IREO(ABS(ICONFP_IN(JOC)))
       IF(JOB_UT.GT.JMAX) JMAX = JOB_UT
      END DO
*
      JOC_UT = 0
      DO JOB_UT = 1, JMAX
      DO JOC_IN = 1, NOC
         IF(IREO(ABS(ICONFP_IN(JOC_IN))).EQ.JOB_UT) THEN
C         IF(IREO(JOC_IN).EQ.ABS(ICONFP_IN(JOB))) THEN
            JOC_UT = JOC_UT + 1
            ICONFP_UT(JOC_UT) = SIGN(JOB_UT,ICONFP_IN(JOC_IN))
          END IF
        END DO
      END DO
*  
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from REO_OB_CONFP '
        WRITE(6,*) ' Input and output configuration '
        CALL IWRTMA(ICONFP_IN,1,NOC,1,NOC)
        WRITE(6,*)
        CALL IWRTMA(ICONFP_UT,1,NOC,1,NOC)
      END IF
*
      RETURN 
      END
      SUBROUTINE GET_NSD_MINMAX_SPACE(MIN_OCC,MAX_OCC,MINMAX_ORB,ISYM,
     &           MS2X,MULTSX,
     &           NSD,NCM,NCSF,NCONF,LOCC)
* 
* Find Number of Determinants for a MINMAX expansion
* defined by MIN_OCC, MAX_OCC
*
* Jeppe Olsen, July 15, 2013
* Last modification; July 23, 2013; Jeppe Olsen; MINMAX_ORB added
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      DIMENSION MIN_OCC(*), MAX_OCC(*), MINMAX_ORB(*)
*. Local scratch
      INTEGER IREO_L(MXPORB)
*
      IDUM = 0  
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTNDMM')
*
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' GET_NSD_MINMAX in action '
       WRITE(6,*) ' ========================= '
      END IF
*
      CALL DIM_FOR_MAXMIN_OCC(MIN_OCC,MAX_OCC,MINMAX_ORB,NACOB,
     &     NINOB+1,ISYM,MINOP,MS2,IDC,MULTS,
     &     NCONF,NCSF,NSD,NCM,LOCC)
*
      IF(NTEST.GE.100) THEN
       WRITE(6,'(A,4(2X,4I9))') 
     & ' Number of CONFs, CSFs, CMs and SDs in MINMAX space ',
     &   NCONF,NCSF, NCM,NSD
       WRITE(6,'(A,I9)') ' Length of occupation array ', LOCC
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTNDMM')
      RETURN
      END
      SUBROUTINE DIM_FOR_MAXMIN_OCC(IOCC_MIN,IOCC_MAX,MINMAX_ORB,NACOB,
     &    IORB_OFF,ISYM,MINOP,MS2,IDC,MULTS,
     &    NCONF,NCSF,NSD,NCM,LOCC)
*
* Generate the number of determinants, combinations, csfs, confs for 
* a occupation space with given MIN and MAX
*
*
* Jeppe Olsen, July 2013
*
* Last modification; July 23, 2013; Jeppe Olsen; MINMAX_ORB added
*       
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
*
*.. Input
*
*. Min and max number of accumulated electrons
      INTEGER IOCC_MIN(NACOB),IOCC_MAX(NACOB), MINMAX_ORB(NACOB)
*. Local scratch
      INTEGER JOCC(2*MXPORB), JACOCC(MXPORB), JACOCC2(MXPORB)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==========================='
        WRITE(6,*) ' Entering DIM_FOR_MAXMIN_OCC'
        WRITE(6,*) ' ==========================='
        WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Min and Max accumulated occupations '
        WRITE(6,*)
        WRITE(6,*) ' Orbital Min. occ Max. occ '
        WRITE(6,*) ' =========================='
        DO IORB = 1, NACOB
          WRITE(6,'(3X,I4,2I3)') 
     &    IORB, IOCC_MIN(IORB), IOCC_MAX(IORB)
        END DO
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,'(A,2I3)') 
     &  ' Smallest number of singly occupied orbitals ',
     &  MINOP
        WRITE(6,*) ' Mults = ', MULTS
      END IF
*
*. Total number of electrons 
      NEL = IOCC_MAX(NACOB)
      IF(NTEST.GE.1000) WRITE(6,*) ' NEL = ', NEL
*. Initialize
      NCONF = 0
      NCSF = 0
      NSD = 0
      NCM = 0
      LOCC = 0
*. Loop over configurations in the form of accumulated occupations
      INI = 1
      NCONF = 0
 1000 CONTINUE
        CALL NEXT_CONF_FROM_MINMAX_OCC(JACOCC,
     &             IOCC_MIN,IOCC_MAX,INI,NONEW,NACOB)
        INI = 0
        IF(NONEW.EQ.0) THEN
*. Reform from accumulated to occupation form
          CALL REFORM_CONF_ACCOCC(JACOCC,JOCC,1,NACOB)
C              REFORM_CONF_ACCOCC(IACCOCC,IOCC,IWAY,NORB)
*. Reform to actual ordering of orbitals in JACOCC2
C              REO_OB_CONFE(ICONFP_IN, ICONFP_UT,IREO_NO,NOB)
          CALL REO_OB_CONFE(JOCC,JACOCC2,MINMAX_ORB,NACOB)
          CALL ICOPVE(JACOCC2,JOCC,NACOB)
C?        WRITE(6,*) ' Next conf in JOCC '
C?        CALL IWRTMA(JOCC,1,NACOB,1,NACOB)
*
*. Check symmetry and number of open orbitals for this space
          JSYM = ISYM_CONF(JOCC,NACOB,IORB_OFF)
C                ISYM_CONF(IOCC,NACOB,IORB_OFF)
          NOPEN = NOP_FOR_CONF_OCC(JOCC,NACOB)
          IF(NTEST.GE.10000)
     &    WRITE(6,*) ' Number of openorbitals ',  NOPEN
*
          IF(JSYM.EQ.ISYM.AND.NOPEN.GE.MINOP) THEN
*. A new configuration to be included, reform and save in packed form
            NCONF = NCONF + 1
            LOCC = LOCC + NOPEN + (NEL-NOPEN)/2
            NOPEN_ALPHA = (NOPEN + MS2)/2
            NSD = NSD + IBION(NOPEN,NOPEN_ALPHA)
            IF(IDC.EQ.1) THEN
              NCM = NSD
            ELSE 
              IF(NOPEN.NE.0) THEN
               NCM = NCM + IBION(NOPEN,NOPEN_ALPHA)/2
              ELSE
               NCM = NCM + 1
              END IF
            END IF !IDC switch
            NCSF = NCSF + IWEYLF(NOPEN,MULTS)
C?          WRITE(6,*) ' NOPEN, MULTS = ', NOPEN, MULTS
          END IF ! End if correct sym and number of open orbitals
      GOTO 1000
        END IF !End if nonew = 0
* 
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*)  ' ======== '
        WRITE(6,*)  ' Results: '
        WRITE(6,*)  ' ======== '
        WRITE(6,*)
        WRITE(6,*) ' Number of configurations of correct symmetry ',
     &       NCONF
        WRITE(6,*) ' Number of CSFs, SDs and CMs: ',
     &       NCSF, NSD, NCM
        WRITE(6,*) ' Length of occupation array ', LOCC
      END IF
*
      RETURN
      END
      SUBROUTINE ANACSF2(LUC,NOCCLS_SPC,IOCCLS_SPC,ISYM,
     &           CIVEC,ICONF_OCC,NCONF_FOR_OPEN,IPROCS,THRES,
     &           MAXTRM,IOUT)
*
* Analyze CI vector in CSF basis
* For CNFBAT = 2, i.e. form where information is not stored
*
*. Jeppe Olsen, July 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      REAL*8 INPROD
*. Input
      DIMENSION IOCCLS_SPC(NOCCLS_SPC)
      DIMENSION ICONF_OCC(*),NCONF_FOR_OPEN(*)
      INTEGER IPROCS(1)
      DIMENSION CIVEC(*)
*. Local scratch
      CHARACTER*3 PAT(MXPORB)
      INTEGER IOCCL(MXPORB)
      PARAMETER(NPOT = 10)
      INTEGER   N_IN_RANGE(NPOT+1)
      DIMENSION W_IN_RANGE(NPOT+1)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ANACSF')
      CALL LUCIAQENTER('ANACSF')
       

* =================================
* Printout of configuration weights
* =================================
*
      NTEST = 00
      IF(NTEST.GE.100) WRITE(6,*) ' ANACSF2, IOUT  = ', IOUT
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(1H ,A)') ' ====================== '
      WRITE(IOUT,'(1H ,A)') ' Info on configurations '
      WRITE(IOUT,'(1H ,A)') ' ====================== '
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(1H ,A)') 
     &' (Negative orbital index implies doubly occupied orbital)'
*
      ITRM = 0
      ILOOP = 0
*. Assume standard normalized matrix
      XNORM = 1.0D0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
*
3001  CONTINUE
       ILOOP = ILOOP + 1
       IF ( ILOOP  .EQ. 1 ) THEN
         XMAX = XNORM
         XMIN = XMAX/SQRT(10.0D0)
       ELSE
         XMAX = XMIN
         XMIN = XMIN/SQRT(10.0D0)
       END IF
       IF(XMIN .LT. THRES  ) XMIN =  THRES
C
       WRITE(IOUT,'(//A,E10.4,A,E10.4)')
     & '  Printout of contributions in interval  ',XMIN,' to ',XMAX
       WRITE(IOUT,'(A)')
     & '  ============================================================='
*
*. Loop over occupation classes
*
       CALL REWINO(LUC)
*
       NCIVAR = 0
       DO IIOCLS = 1, NOCCLS_SPC
        CALL FRMDSCN(CIVEC,1,-1,LUC)
COLD    IB_OCC = 1
COLD    IB_H0 = 1
COLD    IDO_REO = 1
COLD    IB_BLK = 1
COLD    NCSF_TOT = 0
        IOCLS = IOCCLS_SPC(IIOCLS)
C?      WRITE(6,*) ' IIOCLS, IOCLS = ', IIOCLS, IOCLS
C?      IF(NTEST.GE.10) WRITE(6,*) ' Output for IOCLS = ', IOCLS
*. Generate Conformation (only configurations are needed)
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCLS,0,ISYM)
        NCSF_OCCLS = IELSUM(NCS_FOR_OC_OP_ACT,MAXOP+1)
        NCIVAR = NCIVAR + NCSF_OCCLS
*
        IF(NTEST.GE.100) THEN
           WRITE(6,*) ' IIOCLS, IOCLS, NCSF_OCCLS = ',
     &                  IIOCLS, IOCLS, NCSF_OCCLS
        END IF
*. Loop over configurations and CSF's for given configuration
        ICNF = 0
        ICSF = 0
        INRANG = 0
        NOPRT = 0
        DO IOPEN = 0, MAXOP
          IF(NTEST.GE.100) WRITE(6,*) ' IOPEN = ', IOPEN
          ITYP = IOPEN + 1
          ICL = (NACTEL - IOPEN) / 2
          IOCC = IOPEN + ICL
*From previous
          ICLP = (NACTEL-(IOPEN-1))/2
          IOCCP = IOPEN-1+ICLP
*
          IF( ITYP .EQ. 1 ) THEN
            ICNBS0 = 1
            IPBAS = 1
          ELSE
            ICNBS0 = ICNBS0 + NCONF_FOR_OPEN(IOPEN+1-1)*(IOCCP)
            IPBAS = IPBAS + NPCSCNF(IOPEN+1-1)*(IOPEN-1)
          END IF
          IF(NTEST.GE.100) WRITE(6,*) ' ICNBS0, NCONF(Prev), IOCCP = ',
     &                 ICNBS0, NCONF_FOR_OPEN(IOPEN+1-1), IOCCP
*. Configurations of this type
          NNCNF = NCONF_FOR_OPEN(IOPEN+1)
          NNCSF = NPCSCNF(IOPEN+1)
          DO IC = 1, NNCNF
            IF(NTEST.GE.100) WRITE(6,*) ' IC = ', IC
            ICNF = ICNF + 1
            ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
*. Weight of CSF's in this configuration
            W = 0.0D0
            DO IICSF = 1, NNCSF
              ICSF = ICSF+1
              W = W +  CIVEC(ICSF)*CIVEC(ICSF)
            END DO
            SQ_W = SQRT(W)
            IF(XMAX.GE.SQ_W .AND. SQ_W.GT.XMIN ) THEN
              ITRM  = ITRM + 1
              INRANG = INRANG + 1
              IF( ITRM .LE. MAXTRM ) THEN
                WRITE(IOUT,'(A,E10.5)') 
     &          '  Square root of weight for conf ',SQ_W
                CALL ICOPVE(ICONF_OCC(ICNBS),IOCCL,IOCC)
                IF(NTEST.GE.100) WRITE(6,*) ' ICNBS0, ICNBS = ',
     &          ICNBS0,ICNBS
                WRITE(IOUT,*) ' Occupation of configuration: '
                IF(NINOB.NE.0) THEN
                  CALL IADD_SIGN_CONST(IOCCL,NINOB,IOCC)
C                      IADD_SIGN_CONST(IVEC,IADD,NELMNT)
                  WRITE(IOUT,'(4X,A,10(2X,I3),
     &                         (/,16X,10(2X,I3)))')
     &            ' (Inactive) ', (IOCCL(II),II = 1,IOCC)
                ELSE
                  WRITE(IOUT,'(4X,10(2X,I3))')
     &            (IOCCL(II),II = 1,IOCC)
                END IF !NINOB switch
              ELSE
                NOPRT = NOPRT + 1
              END IF! MAXTRM check
            END IF! XMAX/XMIN check
          END DO ! Loop over configs
        END DO ! Loop over IOPEN
       END DO ! Loop over occupation classes
      IF(XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 3001
*
* ==========================
* Printout of CSF-weights
* ==========================
*
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(1H ,A)') ' ============ '
      WRITE(IOUT,'(1H ,A)') ' Info on CSFs '
      WRITE(IOUT,'(1H ,A)') ' ============ '
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      ITRM = 0
      ILOOP = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
2001  CONTINUE
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = XNORM
        XMIN = XMAX/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      WRITE(IOUT,'(//A,E10.4,A,E10.4)')
     &'  Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(IOUT,'(A)')
     &'  =============================================================='
*. Loop over configurations and CSF's for given configuration
       CALL REWINO(LUC)
       DO IIOCLS = 1, NOCCLS_SPC
        CALL FRMDSCN(CIVEC,1,-1,LUC)
        IOCCLS = IOCCLS_SPC(IIOCLS)
        IF(NTEST.GE.100) WRITE(6,*) ' Output from IOCCLS = ', IOCCLS
*. Generate Conformation (only configurations are needed)
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCCLS,0,ISYM)
        NCSF_OCCLS = IELSUM(NCS_FOR_OC_OP_ACT,MAXOP+1)
        NCM_OCCLS = IELSUM(NCM_FOR_OC_OP_ACT,MAXOP+1)
C
        ICNF = 0
        ICSF = 0
        INRANG = 0
        NOPRT = 0
        DO 1000 IOPEN = 0, MAXOP
          ITYP = IOPEN + 1
          ICL = (NACTEL - IOPEN) / 2
          IOCC = IOPEN + ICL
*From previous
          ICLP = (NACTEL-(IOPEN-1))/2
          IOCCP = IOPEN-1+ICLP
*
          IF( ITYP .EQ. 1 ) THEN
            ICNBS0 = 1
            IPBAS = 1
          ELSE
            ICNBS0 = ICNBS0 + NCONF_FOR_OPEN(IOPEN+1-1)*IOCCP
            IPBAS = IPBAS + NPCSCNF(IOPEN+1-1)*(IOPEN-1)
          END IF
*. Configurations of this type
          NNCNF = NCONF_FOR_OPEN(IOPEN+1)
          NNCSF = NPCSCNF(IOPEN+1)
          DO 900  IC = 1, NNCNF
            ICNF = ICNF + 1
            ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
*. CSF's in this configuration
            DO 800 IICSF = 1, NNCSF
              ICSF = ICSF+1
              IF( XMAX .GE. ABS(CIVEC(ICSF)) .AND.
     &           ABS(CIVEC(ICSF)).GT. XMIN ) THEN
                ITRM  = ITRM + 1
                INRANG = INRANG + 1
                IF( ITRM .LE. MAXTRM ) THEN
C                      SPIN_COUPLING_PATTERN(IOCC,ISPIN,IPAT,NCL,NOP)
                  CALL SPIN_COUPLING_PATTERN(ICONF_OCC(ICNBS),
     &                 IPROCS(IPBAS+(IICSF-1)*IOPEN),PAT,ICL,IOPEN)
                  CNORM = CNORM + CIVEC(ICSF) ** 2
                  WRITE(IOUT,*) ' Coefficient of CSF ',ICSF,CIVEC(ICSF)
                  WRITE(IOUT,*) ' Occupation and spin coupling '
                  IF(NINOB.NE.0) THEN
                    WRITE(IOUT,'(4X,A,10(2X,I3),/,16X,10(2X,I3))') 
     &              ' (Inactive) ', 
     &              (ABS(ICONF_OCC(ICNBS-1+II))+NINOB,II = 1,IOCC)
                    WRITE(IOUT,'(5X,A,10(2X,A3),/,17X,10(2X,A3))') 
     &              '            ', (PAT(II),II=1, IOCC) 
                  ELSE
                    WRITE(IOUT,'(4X,A,10(2X,I3))') ' ',
     &              (ABS(ICONF_OCC(ICNBS-1+II)),II = 1,IOCC)
                    WRITE(IOUT,'(5X,A,10(2X,A3))') ' ',
     &              (PAT(II),II=1, IOCC) 
                  END IF! NINOB switch
                  NOPRT = NOPRT + 1
                END IF  !Itrm .le. maxtrm
              END IF ! In range
  800       CONTINUE
  900     CONTINUE
 1000   CONTINUE
       END DO !loop over occupation classes
       IF(INRANG .EQ. 0 ) WRITE(IOUT,*) '   ( no coefficients )'
      IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
      IF(NOPRT.NE.0) WRITE(IOUT,*)
     &' Number of coefficients not printed ', NOPRT
      WRITE(IOUT,'(//A,E15.8)')
     &'  Norm of printed CI vector .. ', CNORM
*
*
*
      WRITE(IOUT,'(/A)') '   Magnitude of CI coefficients '
      WRITE(IOUT,'(A/)') '  =============================='
*
      ZERO = 0.0D0 
      CALL SETVEC(W_IN_RANGE,ZERO,NPOT+1)
      IZERO = 0
      CALL ISETVC(N_IN_RANGE,IZERO,NPOT+1)
*
      ISUM = 0
*
      CALL REWINO(LUC)
      DO IIOCLS = 1, NOCCLS_SPC
        XMIN = 1.0D0
        CALL FRMDSCN(CIVEC,1,-1,LUC)
        IOCCLS = IOCCLS_SPC(IIOCLS)
        IF(NTEST.GE.100) WRITE(6,*) ' Output from IOCCLS = ', IOCCLS
*. Generate Conformation (only configurations are needed)
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCCLS,0,ISYM)
        NCSF_OCCLS = IELSUM(NCS_FOR_OC_OP_ACT,MAXOP+1)
        DO 200 IPOT = 0, 10
         XMAX = XMIN
         XMIN = XMIN * 0.1D0
C
         DO 180 IDET = 1, NCSF_OCCLS 
           IF( ABS(CIVEC(IDET)) .LE. XMAX  .AND.
     &         ABS(CIVEC(IDET)) .GT. XMIN ) THEN
                 N_IN_RANGE(IPOT+1) = N_IN_RANGE(IPOT+1) + 1
                 W_IN_RANGE(IPOT+1) = W_IN_RANGE(IPOT+1) 
     &                              + CIVEC(IDET) ** 2
                 ISUM = ISUM + 1
           END IF
  180    CONTINUE
  200   CONTINUE
      END DO ! Loop over occupation classes
C
      CLNORMT = 0.0D0
      DO IPOT = 0, 10
         INRANG = N_IN_RANGE(IPOT+1)
         CLNORM = W_IN_RANGE(IPOT+1)
         CLNORMT = CLNORMT + CLNORM
         IF (INRANG .GT. 0) THEN
           WRITE(IOUT,'(A,I2,A,I2,3X,I7,3X,E15.8,3X,E15.8)')
     &     '  10-',IPOT+1,' to 10-',IPOT,INRANG,CLNORM,CLNORMT
         END IF
      END DO
*
      WRITE(IOUT,*) ' Number of coefficients less than  10-11',
     &           ' is  ',NCIVAR - ISUM
C
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ANACSF')
      CALL LUCIAQEXIT('ANACSF')
      RETURN
      END
c $Id$
