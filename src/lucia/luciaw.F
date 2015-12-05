      SUBROUTINE GET_MINMAX_ADR_IN_CISPACE
     &           (MINAC,MAXAC,MINMAX_ORB,ISM,ISPC,IADR,NELMNT,ICNF)
*
* A space is given by a MINMAX distribution, MINMAX. 
* Obtain the addresses of components of this space in full space
* and the configurations (if CSFs are in action)
*
*. Jeppe Olsen, July 3, 2013
*
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER MINAC(*), MAXAC(*), MINMAX_ORB(*)
*. Output
      INTEGER IADR(*), ICNF(*)
*
      NTEST = 000
*
      IDUM = 0
      CALL LUCIAQENTER('GTMNAD')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTMNAD')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from GET_MINMAX_ADR_IN_FULL_SPACE'
        WRITE(6,*) ' ========================================'
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and space in action ', ISM, ISPC
      END IF
*
*. Standard def
*
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
*. Size of blocks (assumed in Z_BLKFO)
*
      IF(ISIMSYM.EQ.1.OR.ICISTR.EQ.2) THEN
        LBLOCK = MXSOOB_AS
      ELSE
        LBLOCK = MXSOOB
      END IF
      IF(NOCSF.EQ.0.OR.ICNFBAT.EQ.-2) THEN
CERR    LBLOCK  = MAX(NSD_FOR_OCCLS_MAX,LBLOCK)
        LBLOCK  = MAX(N_SDAB_PER_OCCLS_MAX,LBLOCK)
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      IF(NTEST.GE.100) WRITE(6,*) ' TEST: LBLOCK = ', LBLOCK
*
*
*
*. Information on blocks of CI-expansion
*
      ILTEST = 3006
      CALL Z_BLKFO_FOR_CISPACE(ISPC,ISM,LBLOCK,ICOMP,
     &     NTEST,NCBLOCK,NCBATCH,
     &   int_mb(KCIOIO),int_mb(KCBLTP),NCOCCLS_ACT,int_mb(KCIOCCLS_ACT),
     &     int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &     int_mb(KCIBT),
     &     int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
*. Space for strings
      IF(NOCSF.EQ.1) THEN
        CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR') !done
        CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR') !done
*. 
C       GET_MINMAX_ADR_IN_CISPACE_SD(
C    &           IADR,NDET_UT,MINAC,MAXAC,MINMAX_ORB,NSSOA,NSSOB,NOCTPA,NOCTPB,
C    &           IOCTPA,IOCTPB,NBLOCK,IBLOCK,
C    &           NAEL,NBEL,
C    &           IASTR,IBSTR,IBLTP,NSMST,
C    &           NGAS,NORB,NACOB,NINOB)
        CALL GET_MINMAX_ADR_IN_CISPACE_SD(
     &         IADR,NELMNT,MINAC,MAXAC,MINMAX_ORB,int_mb(KNSTSO(IATP)),
     &         int_mb(KNSTSO(IBTP)),
     &         NOCTPA,NOCTPB,IOCTPA,IOCTPB,NCBLOCK,int_mb(KCIBT),
     &         NAEL,NBEL,
     &         int_mb(KLASTR),int_mb(KLBSTR),int_mb(KCBLTP),NSMST,
     &         NGAS,NTOOB,NACOB,NINOB) 
      ELSE
C        GET_MINMAX_ADR_IN_CISPACE_CSF(IADR,NELMNT,MINAC,MAXAC,MINMAX_ORB,
C    &           NOCCLS_SPC,IOCCLS_SPC,ISYM,ICONF_OCC,NCONF_FOR_OPEN,
C    &           INCLUDE_CONFS,ICONF_OCC_SEL,NOP_CONF_SEL,NCONF_OCC_SEL)
C?       WRITE(6,*) ' NCONF_SUB(1) = ', NCONF_SUB
         CALL GET_MINMAX_ADR_IN_CISPACE_CSF(IADR,NELMNT,MINAC,MAXAC,
     &        MINMAX_ORB,NCOCCLS_ACT,int_mb(KCIOCCLS_ACT),ISM,
     &        int_mb(KICONF_OCC(ISM)),NCONF_PER_OPEN(1,ISM),1,
     &        int_mb(KSBCNFOCC),int_mb(KSBCNFOP),NCONF_SUB)

      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTMNAD')
      CALL LUCIAQEXIT('GTMNAD')
*
      RETURN
      END
      SUBROUTINE GET_MINMAX_ADR_IN_CISPACE_SD(
     &           IADR,NDET_UT,MINAC,MAXAC,MINMAX_ORB,NSSOA,NSSOB,
     &           NOCTPA,NOCTPB,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &           NAEL,NBEL,
     &           IASTR,IBSTR,IBLTP,NSMST,
     &           NGAS,NORB,NACOB,NINOB)
*
* Determine addresses of determinant in a CISPACE that 
* is in the MINMAX space defined by MINMAX
*
*
* Jeppe Olsen, July 2013  
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*· 
*. General input
     
*. Specific input
      INTEGER MINAC(NGAS),MAXAC(NGAS),MINMAX_ORB(*)
      INTEGER IBLOCK(8,NBLOCK)
      INTEGER NSSOA(NSMST,*), NSSOB(NSMST,*)  
      INTEGER IBLTP(*)
*. Scratch
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
*. Local scratch
      INTEGER IACC_CONF(MXPORB),IOCCX(MXPORB), IOCCX2(MXPORB)
      INTEGER IACC2_CONF(MXPORB)
*. Output
      DIMENSION IADR(*)
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' GET_MINMAX_ADR_IN_CISPACE_SD reporting:'
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        WRITE(6,*) ' NINOB, NACOB, NSEL = ', NINOB, NACOB, NSEL
      END IF
*
      IDET_IN = 0
      IDET_UT = 0
*
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
        IF(NTEST.GE.100) THEN
        WRITE(6,'(A,4I4)') 
     &  ' IATP, IBTP, IASM, IBSM = ', IATP, IBTP, IASM, IBSM
        END IF
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                         NORB,0,IDUM,IDUM)
*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        IBBAS = 1
        IABAS = 1
*
        DO  IB = IBBAS,IBBAS+NIB-1
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB - IBBAS + IABAS
          ELSE
            MINIA = IABAS
          END IF
          DO  IA = MINIA,IABAS+NIA-1
*
            IDET_IN = IDET_IN + 1
*. Is this determinent in MINMAX space
*. Accumulated form
C                IAIB_TO_ACCCONF(IA,IB,NAEL,NBEL,IACC,NACOB,NINOB)
            CALL IAIB_TO_ACCCONF
     &           (IASTR(1,IA),IBSTR(1,IB),NAEL, NBEL,
     &           IACC_CONF,NACOB,NINOB)
*. Return to standard (not accumulated )form in IOCCX
C                REFORM_CONF_ACCOCC(JACOCC,JOCC,1,NORBL)
            CALL REFORM_CONF_ACCOCC(IACC_CONF,IOCCX,1,NACOB)
*. Reorder orbitals to the order that is assumed in the min max arrays
C              REO_OB_CONFE(ICONFP_IN, ICONFP_UT,IREO_NO,NOB)
             CALL REO_OB_CONFE(IOCCX,IOCCX2,MINMAX_ORB,NACOB)
*. And put reordered configuration in accumulated form
            CALL REFORM_CONF_ACCOCC(IACC2_CONF,IOCCX2,2,NACOB)
            IM_IN = IS_IACC_CONF_IN_MINMAX_SPC 
     &              (IACC2_CONF,MINAC,MAXAC,NACOB)
C                   IS_IACC_CONF_IN_MINMAX_SPC(IOCC,MIN_OCC,MAX_OCC,NORB)
            IF(IM_IN.EQ.1) THEN
*. Enroll!
              IF(NTEST.GE.1000) 
     &        WRITE(6,*) ' Determinant in MINMAX space'
              IDET_UT = IDET_UT + 1
              IADR(IDET_UT) = IDET_IN 
            END IF
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
      END DO
*     ^ End of loop over blocks
      NDET_UT = IDET_UT
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Obtained number of dets in MINMAX-space ', IDET_UT
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Address of the obtained determinants'
        CALL IWRTMA(IADR,1,NDET_UT,1,NET_UT)
      END IF
*
      RETURN
      END
      SUBROUTINE IAIB_TO_ACCCONF(IA,IB,NAEL,NBEL,IACC,NACOB,NINOB)
*
* Alpha and beta strings are given, obtain corresponding 
* accumulated configuration
*
*. Note: configuration is only over active orbitals, whereas
*        strings has ninob + 1 as first active orbital
*
*. Jeppe Olsen, July 3, 2013
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IA(NAEL),IB(NBEL)
*. Output
      INTEGER IACC(NACOB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from IAIB_TO_ACCCONF '
        WRITE(6,*) ' =========================='
        WRITE(6,*)
        WRITE(6,*) ' NAEL, NBEL, NACOB = ', NAEL, NBEL,NACOB
      END IF
*
      IZERO = 0
      CALL ISETVC(IACC,IZERO,NACOB)
*
      DO JAEL = 1, NAEL
       IORB = IA(JAEL) - NINOB
       DO JORB = IORB, NACOB
         IACC(JORB) = IACC(JORB) + 1
       END DO
      END DO
*
      DO JBEL = 1, NAEL
       IORB = IB(JBEL) - NINOB
       DO JORB = IORB, NACOB
         IACC(JORB) = IACC(JORB) + 1
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Input: a- and b-strings: '
       CALL IWRTMA(IA,1,NAEL,1,NAEL)
       CALL IWRTMA(IB,1,NBEL,1,NBEL)
*
       WRITE(6,*) ' Output: accumulated occupation '
       CALL IWRTMA(IACC,1,NACOB,1,NACOB)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_IAIB_FOR_SEL_DETS(ISM,ISPC,ISEL,NSEL,IA,IB)
*
* Obtain the alpha- and  beta-strings for selected determinants in CI-space ISPC
*
*. Jeppe Olsen, July 4, 2013 (my old man would have turned 86 today...)
*
*
      INCLUDE 'implicit.inc' 
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Input
      INTEGER ISEL(NSEL)
*. Output
      INTEGER IA(*), IB(*)
*
      NTEST = 000
*
      IDUM = 0
      CALL LUCIAQENTER('GTSLDT')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTSLDT')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from GET_IAIB_FOR_SEL_DETS '
        WRITE(6,*) ' =================================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and space in action ', ISM, ISPC
        WRITE(6,*) ' Number of selected determinants ', NSEL
      END IF
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' Selected determinants: '
        CALL IWRTMA(ISEL,1,NSEL,1,NSEL)
      END IF
*
*. Standard def
*
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
*. Size of blocks (assumed in Z_BLKFO)
*
      IF(ISIMSYM.EQ.1.OR.ICISTR.EQ.2) THEN
        LBLOCK = MXSOOB_AS
      ELSE
        LBLOCK = MXSOOB
      END IF
      IF(NOCSF.EQ.0.OR.ICNFBAT.EQ.-2) THEN
CERR    LBLOCK  = MAX(NSD_FOR_OCCLS_MAX,LBLOCK)
        LBLOCK  = MAX(N_SDAB_PER_OCCLS_MAX,LBLOCK)
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      IF(NTEST.GE.100) WRITE(6,*) ' TEST: LBLOCK = ', LBLOCK
*
*. Information on blocks of CI-expansion
*
        ILTEST = 3006
        CALL Z_BLKFO_FOR_CISPACE(ISPC,ISM,LBLOCK,ICOMP,
     &       NTEST,NCBLOCK,NCBATCH,
     &       int_mb(KCIOIO),int_mb(KCBLTP),NCOCCLS_ACT,
     &       int_mb(KCIOCCLS_ACT),
     &       int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &       int_mb(KCIBT),
     &       int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
*. Space for strings
        CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR') !done
        CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR') !done
*. 
        CALL GET_IAIB_FOR_SEL_DETS_IN(
     &           ISEL,NSEL,IA,IB,int_mb(KNSTSO(IATP)),
     &           int_mb(KNSTSO(IBTP)),
     &           NOCTPA,NOCTPB,IOCTPA,IOCTPB,NCBLOCK,int_mb(KCIBT),
     &           NAEL,NBEL,
     &           int_mb(KLASTR),int_mb(KLBSTR),int_mb(KCBLTP),NSMST,
     &           NGAS,NTOOB,NACOB,NINOB) 
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GTSLDT')
      CALL LUCIAQEXIT('GTSLDT')
*
      RETURN
      END
      SUBROUTINE GET_IAIB_FOR_SEL_DETS_IN(
     &           ISEL,NSEL,IA_UT,IB_UT,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &           IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &           NAEL,NBEL,
     &           IASTR,IBSTR,IBLTP,NSMST,
     &           NGAS,NORB,NACOB,NINOB)
*
* Obtain alpha- and beta-strings for determinants with addresses given by ISEL
*
*
* Jeppe Olsen, July 2013  
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Specific input
      INTEGER ISEL(NSEL)
      INTEGER IBLOCK(8,NBLOCK)
      INTEGER NSSOA(NSMST,*), NSSOB(NSMST,*)  
      INTEGER IBLTP(*)
*. Scratch
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
*. Local scratch
      DIMENSION IACC_CONF(MXPORB)
*. Output
      DIMENSION IA_UT(NAEL,NSEL),IB_UT(NBEL,NSEL)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_IAIB_FOR_SEL_DETS_IN '
        WRITE(6,*) ' Requested dets: '
        CALL IWRTMA(ISEL,1,NSEL,1,NSEL)
      END IF
*
      IDET_IN = 0
      IDET_UT = 1
*
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
        IF(NTEST.GE.1000) THEN
        WRITE(6,'(A,4I4)') 
     &  ' IATP, IBTP, IASM, IBSM = ', IATP, IBTP, IASM, IBSM
        END IF
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                         NORB,0,IDUM,IDUM)
*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        IBBAS = 1
        IABAS = 1
*
        DO  IB = IBBAS,IBBAS+NIB-1
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB - IBBAS + IABAS
          ELSE
            MINIA = IABAS
          END IF
          DO  IA = MINIA,IABAS+NIA-1
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' IA, IB  = ', IA, IB
            END IF
*
            IDET_IN = IDET_IN + 1
            IF(IDET_IN.EQ.ISEL(IDET_UT)) THEN
*. Next det has been determined, enroll
              CALL ICOPVE(IASTR(1,IA),IA_UT(1,IDET_UT),NAEL)
              CALL ICOPVE(IBSTR(1,IB),IB_UT(1,IDET_UT),NBEL)
              IF(IDET_UT.EQ.NSEL) GOTO 1001
              IDET_UT = IDET_UT + 1
            END IF
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
      END DO
*     ^ End of loop over blocks
 1001 CONTINUE
*
*. Check that the required number of dets was obtained
*
      NDET_UT = IDET_UT
      IF(NDET_UT.NE.NSEL) THEN
        WRITE(6,*) ' Obtained number of dets differ from requested '
        WRITE(6,*) ' Obtained and requested dimensions: ', NDET_UT, NSEL
        STOP ' Obtained number of dets differ from requested '
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Obtained alpha- and beta-strings: '
        DO JSEL = 1, NSEL
          WRITE(6,*) ' Determinant ', ISEL(JSEL)
          WRITE(6,'(4X,10I4)') (IA_UT(IEL,JSEL),IEL = 1, NAEL )
          WRITE(6,'(4X,10I4)') (IB_UT(IEL,JSEL),IEL = 1, NBEL )
        END DO
      END IF
*
      RETURN
      END
      SUBROUTINE GET_SUBSPC_PRECOND_SPC(ISPC,ISM,ISEL,NSEL, 
     &           CBLK)
*
*
* Obtain the preconditioner subspace in the form of 
* a set of addresses of variables.
*
* It is assumed that diagonal has been calculated and stored on LUDIA
*
*. Jeppe Olsen, July 4, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'stinf.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Scratch holding a block of CI vector
      DIMENSION CBLK(*)
*. Output
      DIMENSION  ISEL(*) 
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBSPCN')
      CALL LUCIAQENTER('SBSPCN')
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_SUBSPC_PRECOND_SPC:'
        WRITE(6,*) ' ================================='
        WRITE(6,*)
        WRITE(6,*) ' ISBSPC_SEL = ', ISBSPC_SEL
        WRITE(6,*) ' MXP1, MXP2, MXQ = ', MXP1, MXP2, MXQ
      END IF 
*
*. Some general info
*
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
      MXDM = MXP1 + MXP2 + MXQ
*
*. Obtain the determinants to be included in the subspace
*
      IF(ISBSPC_SEL.EQ.1) THEN
*
*. Obtain subspace from lowest elements of CI diagonal
*
*. Local scratch
        CALL MEMMAN(KL1,3*(MXDM+1),'ADDL  ',1,'KL1   ')  !done
        CALL MEMMAN(KL2,2*(MXDM+1),'ADDL  ',2,'KL2   ')  !done
        CALL MEMMAN(KL3,2*(MXDM+1),'ADDL  ',2,'KL3   ')  !done
        WRITE(6,*) ' MXDM, KL1, KL2, KL3 = ', KL1, KL2, KL3
        LBLK = -1
*. And determine total subspace space 
        CALL FNDMND(LUDIA,LBLK,CBLK,MXDM,NPRDET,int_mb(KL1),
     &              dbl_mb(KL2),ISEL,dbl_mb(KL3),NTEST )

*
* Check for boundaries between P1, P2, and Q
*
* P1-P2
        IF(MXP1 .GT. 0 ) THEN
          IIDET = MXP1
 101      CONTINUE
          IF(ABS(dbl_mb(KL3-1+IIDET+1)-dbl_mb(KL3-1+IIDET))
     &         .LE. 0.000001D0) THEN
            IIDET = IIDET - 1
            GOTO 101
          END IF
          NP1 = IIDET
        ELSE
          NP1 = 0
        END IF
        IF(NTEST .GE. 2)
     &   WRITE(6,*) ' Actual dimension of P1 Space ', NP1
*. P2 - Q space
        IF(MXP2.GT.0) THEN
          IF(MXP1+MXP2.GE.NPRDET) THEN
            NP2 = NPRDET - NP1
          ELSE
            IIDET = MXP1 + MXP2
 102        CONTINUE
            IF( ABS(dbl_mb(KL3-1+IIDET+1)-dbl_mb(KL3-1+IIDET))
     &         .LE. 0.0000001) THEN
               IIDET = IIDET - 1
               GOTO 102
            END IF
            NP2 = IIDET - NP1
          END IF
        ELSE
          NP2 = 0
        END IF
        IF( NTEST .GE. 2 )
     &   WRITE(6,*) ' Actual dimension of P2 Space ', NP2
*. Q space
        IF(MXQ.NE.0) THEN
          NQ = MXP1 + MXP2 + MXQ - NP1 - NP2
        ELSE
          NQ = 0
        END IF
        IF( NTEST .GE. 2  )
     &   WRITE(6,*) ' Actual dimension of Q Space ', NQ
        NPVAR = NP1 + NP2
        NPRVAR = NP1 + NP2 + NQ
*. The determinants/CSFs should be delivered in ascending order, so sort
* ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
        CALL ORDINT(ISEL,int_mb(KL1),NP1,dbl_mb(KL2),0)
        CALL ICOPVE(int_mb(KL1),ISEL,NP1)
*. Should add for P2 and Q space when and if relevant
      ELSE IF (ISBSPC_SEL.EQ.2) THEN
*
*. Just choose the first elements
*
* No check that the dimensions are less or equal to dim of actual space..
        NP1 = MXP1
        NP2 = MXP2
        NQ = MXQ
        NPVAR = NP1 + NP2
        NPRVAR = NP1 + NP2 + NQ
C ISTVC2(IVEC,IBASE,IFACT,NDIM)
        CALL ISTVC2(ISEL,0,1,NPRVAR)
*
      ELSE IF (ISBSPC_SEL.EQ.3) THEN
*. A CI space is chosen as explicit preconditioner space
         WRITE(6,*) ' STOP: ISPSPC_SEL = 3 has not been programmed yet '
         STOP ' ISPSPC_SEL = 3 has not been programmed yet '
      ELSE IF (ISBSPC_SEL.EQ.4) THEN
*
* Obtain subspace from a MINMAX space
*
        IF(NOCSF.EQ.1) THEN
*. Define parameters connected with CSFs
          MULTS = MS2 + 1
          MINOP = 0
        END IF
*
        CALL GET_NSD_MINMAX_SPACE(ISBSPC_MINMAX(1,1),ISBSPC_MINMAX(1,2),
     &       ISBSPC_ORB,ISM,MS2,MULTS,NSD,NCM,NCSF,NCONF,LOCC)
C            GET_NSD_MINMAX_SPACE(MIN_OCC,MAX_OCC,ISYM,MS2X,MULTSX,
C    &           NSD,NCM,NCSF,NCONF,LOCC)
        IF(NOCSF.EQ.0) THEN
          NP1 = NCSF
        ELSE
          NP1 = NCM
        END IF
        NPRVAR = NP1
        NSEL = NP1
        NP2 = 0
        NQ = 0
*. And the space
C       GET_MINMAX_ADR_IN_CISPACE(MINAC,MAXAC,ISM,ISPC,IADR,NELMNT)
        CALL GET_MINMAX_ADR_IN_CISPACE(
     &       ISBSPC_MINMAX(1,1),ISBSPC_MINMAX(1,2),ISBSPC_ORB_INV,
     &       ISM,ISPC,ISEL,NSEL)
      END IF
      NSEL = NPRVAR
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Dimension of preconditioner subspace =', NP1
      END IF
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' And the addresses of the subspace variables'
        CALL IWRTMA(ISEL,1,NSEL,1,NSEL)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SBSPCN')
      CALL LUCIAQEXIT('SBSPCN')
*
      RETURN
      END
      SUBROUTINE GET_SUBSPC_PRECOND_MAT(ISPC,ISM,H0,ISEL,NSEL, 
     &           EIGVAL, EIGVEC)
*
* Obtain subspace preconditioner matrix for CI  
*
* The preconditioner space is assumed already determined and is 
* given by SEL, NSEL
*
*
* At the moment a single space preconditioner is assumed
*
* NP1, NP2, NQ transferred through common block
*
*. Jeppe Olsen, July 4, 2013, last change July 22, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Output
      DIMENSION H0(*), ISEL(*), EIGVAL(*), EIGVEC(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBSPMT')
      CALL LUCIAQENTER('SBSPMT')
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from GET_SUBSPC_PRECOND_MAT: '
        WRITE(6,*) ' =================================='
        WRITE(6,*)
        WRITE(6,*) ' ISPC, ISM = ', ISPC, ISM
        WRITE(6,*) ' Dimension of subspace = ', NSEL
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Addresses of subspace '
        CALL IWRTMA(ISEL,1,NSEL,1,NSEL)
      END IF 
      IF(NP2.NE.0.OR.NQ.NE.0) THEN
        WRITE(6,*) ' NP2 or NQ ne 0, ', NP2, NQ
        WRITE(6,*) ' Only P1 preconditioner in action '
        STOP ' NP2 or NQ ne 0 '
      END IF
*
*. Some general info
*
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
      MXDM = MXP1 + MXP2 + MXQ
*
      NPRVAR = NP1
      NSEL = NPRVAR
*
* Obtain the SD/CSFs defining the P-space
*
COLD  GET_IAIB_FOR_SEL_DETS(
COLD &           ISEL,NSEL,IA_UT,IB_UT,NSSOA,NSSOB,NOCTPA,NOCTPB,
COLD &           IOCTPA,IOCTPB,NBLOCK,IBLOCK,
COLD &           NAEL,NBEL,
COLD &           IASTR,IBSTR,IBLTP,NSMST,
COLD &           NGAS,NORB,NACOB,NINOB)
      IF(NOCSF.EQ.1) THEN
*. Obtain alpha and beta-strings for the selected determinants
        CALL MEMMAN(KLIASTR,NPRVAR*NAEL,'ADDL  ',1,'IASTR ') !done
        CALL MEMMAN(KLIBSTR,NPRVAR*NBEL,'ADDL  ',1,'IBSTR ') !done
C       GET_IAIB_FOR_SEL_DETS(ISM,ISPC,ISEL,NSEL,IA,IB)
        CALL GET_IAIB_FOR_SEL_DETS(ISM,ISPC,ISEL,NPRVAR,
     &        int_mb(KLIASTR),int_mb(KLIBSTR))
*
*. And obtain the corresponding Hamilton matrix
*
C       DIHDJ2_LUCIA_CONF
C    &  (IASTR,IBSTR,NIDET,JASTR,JBSTR,NJDET,NAEL,NBEL,IADOB,NORB,
C    &   IHORS,HAMIL,C,SIGMA,IWORK,ISYM,ECORE,ICOMBI,PSIGN,
C    &   NTERMS,NDIF0,NDIF1,NDIF2,I12OP,I_DO_ORBTRA,IORBTRA,
C    &   NTOOB,RJ,RK)
        XDUM = 0.0D0
        LSCR = 4*NTOOB + NSEL
        IF(PSSIGN.NE.0.0D0) THEN
          ICOMBI_L = 1
        ELSE
          ICOMBI_L = 0
        END IF
        CALL MEMMAN(KLSCR,LSCR,'ADDL  ',1,'LSCR  ')  !done
        XRJ = -1.0D0
        XRK = -1.0D0
*. In DIHDJ2 it is assumed that the I-and J-strings are in different 
*. arrays (I-strings are interchanged when using combinations). 
*. add extra copy if combinations are active
        IF(ICOMBI_L.EQ.0) THEN
          KLJASTR = KLIASTR
          KLJBSTR = KLIBSTR
        ELSE
          CALL MEMMAN(KLJASTR,NPRVAR*NAEL,'ADDL  ',1,'LJASTR') !done
          CALL MEMMAN(KLJBSTR,NPRVAR*NBEL,'ADDL  ',1,'LJBSTR') !done
          CALL ICOPVE(int_mb(KLIASTR),int_mb(KLJASTR),NPRVAR*NAEL)
          CALL ICOPVE(int_mb(KLIBSTR),int_mb(KLJBSTR),NPRVAR*NBEL)
        END IF
*
        CALL DIHDJ2_LUCIA_CONF(
     &       int_mb(KLIASTR),int_mb(KLIBSTR),NPRVAR,
     &       int_mb(KLJASTR),int_mb(KLJBSTR),NPRVAR,NAEL,NBEL,0,NTOOB,
     &       1, H0,XDUM,XDUM,int_mb(KLSCR),ISM,ECORE,ICOMBI_L,PSSIGN,
     &       NTERMS, NDIF0,NDIF1,NDIF2,2,0,IDUM,NTOOB,XRJ,XRK)
      ELSE
*
* CSF approach
*
C            CNHCN_FOR_CNLIST(ICNOCC,ICNOP,NCN,HCSF,ISCR,SCR,RJ,RK)
*
* Scratch
*
        NOP_MAX = IMNMX(int_mb(KSBCNFOP),NCONF_SUB,2)
        NPDT_MAX = NPDTCNF(NOP_MAX+1)
        WRITE(6,*) ' NOP_MAX, NPDT_MAX = ',  NOP_MAX, NPDT_MAX
        LISCR = 2*NPDT_MAX*NACTEL + NPDT_MAX + 6*NACOB
        LRSCR = 2*NPDT_MAX**2
        CALL MEMMAN(KLISCR,LISCR,'ADDL  ',1,'CNISCR') !done
        CALL MEMMAN(KLRSCR,LRSCR,'ADDL  ',2,'CNRSCR') !done
C?      WRITE(6,*) ' NCONF_SUB(2) = ', NCONF_SUB
        CALL CNHCN_FOR_CNLIST(int_mb(KSBCNFOCC),int_mb(KSBCNFOP),
     &       NCONF_SUB,H0,int_mb(KLISCR),dbl_mb(KLRSCR),XRJ,XRK)
C            CNHCN_FOR_CNLIST(ICNOCC,ICNOP,NCN,HCSF,ISCR,SCR,RJ,RK)
      END IF! Dets of CSFs are in use
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output subspace Hamilton matrix '
        CALL PRSYM(H0,NPRVAR)
      END IF
*
*. Diagonalize
*
*. Outpack matrix to complete form
      CALL TRIPAK(EIGVEC,H0,2,NPRVAR,NPRVAR)
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
*. and diagonalize
      CALL MEMMAN(KLSCRVEC,NSEL,'ADDL  ',2,'SCRVEC') !done
C          DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
      CALL DIAG_SYMMAT_EISPACK(EIGVEC,EIGVAL,dbl_mb(KLSCRVEC),NPRVAR,
     &     IRETURN)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Lowest subspace eigenvalues: '
        NPRINT = 10
      ELSE IF (NTEST.GE.1000) THEN
        WRITE(6,*) ' Subspace eigenvalues: '
        NPRINT = NPRVAL
      END IF
      IF(NTEST.GE.100) THEN
        CALL WRTMAT(EIGVAL,1,NPRINT,1,NPRINT)
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Subspace eigenvectors '
        CALL WRTMAT(EIGVEC,NPRVAR,NPRVAR,NPRVAR)
      END IF
*
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SBSPMT')
      CALL LUCIAQEXIT('SBSPMT')
*
      RETURN
      END
C KNCN_PER_OP_SM
      SUBROUTINE AVE_CSFDIA_CNF(LUDIA,LUDIA_AV,NOCCLS_SPC,IOCCLS_SPC,
     &           ISM,CIVEC,NCN_PER_OP_SM)
*
* A CSF diagonal of H is given on LUDIA
*. Average over elements belonging to the same configuration
*
*. Jeppe Olsen, July 19, 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'crun.inc'
*
      INTEGER NCN_PER_OP_SM(MAXOP+1,NIRREP,*)
*. Specific input
      INTEGER IOCCLS_SPC(NOCCLS_SPC)
*. Scratch
      DIMENSION CIVEC(*)
*
      IF(ICNFBAT.EQ.1) THEN
        WRITE(6,*) ' Average of CSF diag not programmed for CNFBAT = 1'
        STOP       ' Average of CSF diag not programmed for CNFBAT = 1'
      ELSE
*
       CALL REWINO(LUDIA)
       CALL REWINO(LUDIA_AV)
*
       DO IIOCCLS = 1, NOCCLS_SPC
         IOCCLS = IOCCLS_SPC(IIOCCLS)
         CALL FRMDSCN(CIVEC,1,-1,LUDIA)
         ICNBS = 1
         DO IOPEN = 0, MAXOP
           NNCNF = NCN_PER_OP_SM(IOPEN+1,ISM,IOCCLS)
           NPCSF = NPCSCNF(IOPEN+1)
           DO ICNF = 1, NNCNF
             DIASUM = ELSUM(CIVEC(ICNBS),NPCSF)
             AVE = DIASUM/FLOAT(NPCSF)
             CALL SETVEC(CIVEC(ICNBS),AVE,NPCSF)
             ICNBS = ICNBS + NPCSF
           END DO
         END DO
         LENGTH = ICNBS - 1
         CALL TODSCN(CIVEC,1,LENGTH,-1,LUDIA_AVE)
       END DO
      END IF !  CNFBAT switch
*
      RETURN
      END
      SUBROUTINE CNHCN_FOR_CNLIST(ICNOCC,ICNOP,NCN,HCSF,ISCR,SCR,RJ,RK)
*
* Calculate CI matrix for list of configurations specified by ICNOCC,ICNOP
*
*. Jeppe Olsen, July 2013
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cecore.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
*. Specific input
      INTEGER ICNOCC(*), ICNOP(*)
      DIMENSION RJ(*), RK(*)
*. Scratch through input
      DIMENSION ISCR(*), SCR(*)
*. Output
      DIMENSION HCSF(*)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CNHCN_FOR_CNLIST: '
        WRITE(6,*) ' =============================='
        WRITE(6,*)
        WRITE(6,*)  ' Number of configurations in action = ', NCN
      END IF
      
      CALL LUCIAQENTER('CNHCNL')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CNHCNL')
*. Largest number of open shells
      IOP_MAX = 0
      NCSF_T = 0
      DO ICN = 1, NCN
        IOP_MAX = MAX(IOP_MAX,ICNOP(ICN))
        NCSF_T = NCSF_T + NPCSCNF(ICNOP(ICN)+1)
      END DO
      NCSF_MAX = NPCSCNF(IOP_MAX+1)
*. Local memory for a H-matrix over a conf
      LHCNF = NCSF_MAX**2
      CALL MEMMAN(KLHCNF,LHCNF,'ADDL  ',2,'HCNF  ') !done
*
      IB_OCL = 1
      IB_CSL = 1
*
      DO ICNL = 1, NCN
        IOPL = ICNOP(ICNL)
        ICLL = (NACTEL - IOPL)/2
        IOCL = IOPL + ICLL
        NCSFL = NPCSCNF(IOPL+1)
        IB_CSR = 1
        IB_OCR = 1
        DO ICNR = 1, ICNL
          IOPR = ICNOP(ICNR)
          ICLR = (NACTEL - IOPR)/2
          IOCR = IOPR + ICLR
          NCSFR = NPCSCNF(IOPR+1)
          IF(ICNL.EQ.ICNR) THEN
            ISYMG = 1
          ELSE 
            ISYMG = 0
          END IF
*. For test
          ISYMG = 0

*. H-matrix over Confs in WORK(KLHCNF)
C             CNHCN_CSF_BLK(ICNL,IOPL,ICNR,IOPR,CNHCNM,IADOB,
C    &                     IPRODT,DTOC,I12OP,ISCR,SCR,ECORE,IONLY_DIAG,ISYMG,
C    &                     RJ, RK)
          CALL CNHCN_CSF_BLK(ICNOCC(IB_OCL),IOPL,ICNOCC(IB_OCR),IOPR,
     &                 dbl_mb(KLHCNF),NINOB,int_mb(KDFTP),int_mb(KDTOC),
     &                 2,ISCR,SCR,ECORE,0,ISYMG,RJ,RK)
*. Expand to complete matrix
C     EXTR_OR_CP_MAT(ABIG,LRBIG,LCBIG,ISYMBIG,
C    &                          ASMA,LRSMA,LCSMA,ISYMSMA,
C    &                          IOFFR,IOFFC,IEC)
          CALL EXTR_OR_CP_MAT(HCSF,NCSF_T,NCSF_T,1,
     &         dbl_mb(KLHCNF),NCSFL,NCSFR,ISYMG,IB_CSL,IB_CSR,2)

*. Update pointers
          IB_OCR = IB_OCR + IOCR
          IB_CSR = IB_CSR + NCSFR
        END DO ! Loop over ICNR
*. Update pointers
        IB_OCL = IB_OCL + IOCL
        IB_CSL = IB_CSL + NCSFL
      END DO ! Loop over ICNL


        
         
      CALL LUCIAQEXIT('CNHCNL')
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'CNHCNL')
*
      RETURN
      END
      SUBROUTINE EXTR_OR_CP_MAT(ABIG,LRBIG,LCBIG,ISYMBIG,
     &                          ASMA,LRSMA,LCSMA,ISYMSMA,
     &                          IOFFR,IOFFC,IEC)
*
*
* Copy or extract a smaller matrix, ASMA, from/to a larger matrix, ABIG
*
* IEC = 1 => Extract from big to small matrix
* IEC = 2 => Extract from small to big matrix
*
*. Jeppe Olsen, July 19, 2013
*
      INCLUDE 'implicit.inc'
*. Input or output
      DIMENSION ASMA(*),ABIG(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from EXTR_OR_CP_MAT'
        WRITE(6,*) ' ========================'
        WRITE(6,*) 
        WRITE(6,*) ' LRBIG, LCBIG = ', LRBIG, LCBIG
        WRITE(6,*) ' LRSAM, LCSMA = ', LRSMA, LCSMA 
        WRITE(6,*) ' IOFFR, IOFFC = ', IOFFR, IOFFC
        WRITE(6,*) ' ISYMBIG, ISYMSMA = ', ISYMBIG,ISYMSMA
        WRITE(6,*) ' IEC = ', IEC
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input small matrix '
        IF(ISYMSMA.EQ.0) THEN
          CALL WRTMAT(ASMA,LRSMA, LCSMA, LRSMA, LCSMA)
        ELSE
          CALL PRSYM(ASMA, LRSMA)
        END IF
      END IF
*
*
*. Extract small matrix from larger
*
      DO JC_SMA = 1, LCSMA
        IF(ISYMSMA.EQ.0) THEN
          JRMIN = 1
        ELSE
          JRMIN = JC
        END IF
        DO JR_SMA = JRMIN, LRSMA
          JC_BIG = JC_SMA + IOFFC - 1
          JR_BIG = JR_SMA + IOFFR - 1
          IF(ISYMBIG.EQ.0.OR.(ISYMBIG.EQ.1.AND.JR_BIG.GE.JC_BIG)) THEN
            IF(ISYMSMA.EQ.0) THEN
             IADR_SMA = (JC_SMA-1)*LRSMA+JR_SMA
            ELSE
             IADR_SMA = JR_SMA*(JR_SMA-1)/2 + JC_SMA
            END IF
            IF(ISYMBIG.EQ.0) THEN
             IADR_BIG = (JC_BIG-1)*LRBIG+JR_BIG
            ELSE
             IADR_BIG = JR_BIG*(JR_BIG-1)/2 + JC_BIG
            END IF
            IF(IEC.EQ.1) THEN
             ASMA(IADR_SMA) = ABIG(IADR_BIG)
            ELSE
             ABIG(IADR_BIG) =  ASMA(IADR_SMA)
            END IF
          END IF
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE GET_MINMAX_ADR_IN_CISPACE_CSF(IADR,NELMNT,MINAC,MAXAC,
     &           MINMAX_ORB,
     &           NOCCLS_SPC,IOCCLS_SPC,ISYM,ICONF_OCC,NCONF_FOR_OPEN,
     &           INCLUDE_CONFS,ICONF_OCC_SEL,NOP_CONF_SEL,NCONF_OCC_SEL)
*
* Address in CI space of CSF's belonging to a given MINMAX space
*
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
      DIMENSION MINAC(NACOB),MAXOC(NACOB),MINMAX_ORB(NACOB)
      DIMENSION IOCCLS_SPC(NOCCLS_SPC)
      DIMENSION ICONF_OCC(*),NCONF_FOR_OPEN(*)
*. Local scratch
      INTEGER IOCCL(MXPORB),IOCCL2(MXPORB)
*. Output 
      INTEGER IADR(*)
      INTEGER ICONF_OCC_SEL(*), NOP_CONF_SEL(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GTMNCS')
      CALL LUCIAQENTER('GTMNCS')
       

      IOUT = 6
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(IOUT,*)
        WRITE(IOUT,'(1H ,A)') ' ===================================== '
        WRITE(IOUT,'(1H ,A)') ' Info from GET_MINMAX_ADR_IN_CISPACE_CS'
        WRITE(IOUT,'(1H ,A)') ' ===================================== '
        WRITE(IOUT,*)
        WRITE(IOUT,*)
      END IF
*
*
*. Loop over occupation classes
*
*
      ISEL = 0
      NCIVAR = 0
      ICSF = 0
      NCONF_SEL = 0
      IBCONF_OCC_SEL = 1
      DO IIOCLS = 1, NOCCLS_SPC
        IOCLS = IOCCLS_SPC(IIOCLS)
*. Generate Conformation (only configurations are needed)
        CALL GEN_CNF_INFO_FOR_OCCLS(IOCLS,0,ISYM)
        NCSF_OCCLS = IELSUM(NCS_FOR_OC_OP_ACT,MAXOP+1)
        NCIVAR = NCIVAR + NCSF_OCCLS
*
        IF(NTEST.GE.200) THEN
           WRITE(6,*) ' IIOCLS, IOCLS, NCSF_OCCLS = ',
     &                  IIOCLS, IOCLS, NCSF_OCCLS
        END IF
*. Loop over configurations and CSF's for given configuration
        ICNBS0 = 1
        DO IOPEN = 0, MAXOP
          IF(NTEST.GE.200) WRITE(6,*) ' IOPEN = ', IOPEN
          ITYP = IOPEN + 1
          ICL = (NACTEL - IOPEN) / 2
          IOCC = IOPEN + ICL
*. Configurations of this type
          NNCNF = NCONF_FOR_OPEN(IOPEN+1)
          NNCSF = NPCSCNF(IOPEN+1)
          DO IC = 1, NNCNF
            IF(NTEST.GE.1000) WRITE(6,*) ' IC = ', IC
            ICNBS = ICNBS0 + (IC-1)*IOCC
            IF(NTEST.GE.1000) WRITE(6,*) ' IC, ICNBS = ', IC, ICNBS
*. Is this configuration in minmax?
*. Change first the configuration in the order required for the minmax check
*. packed to expanded:
C                REFORM_CONF_OCC2(ICONF_EXP,ICONF_PACK,NORBL,NOCOBL,IWAY)
*. Packed in ICONF_OCC => expanded in IOCCL
            CALL REFORM_CONF_OCC2(IOCCL,ICONF_OCC(ICNBS),NACOB,IOCC,2)
C                REO_OB_CONFE(ICONFE_IN, ICONFE_UT,IREO_NO,NOB)
*. Expanded in IOCCL => reordered expanded in IOCCL2
            CALL REO_OB_CONFE(IOCCL,IOCCL2,MINMAX_ORB,NACOB)
*. Reordered expanded in IOCCL2 to reordered packed in IOCCL
            CALL REFORM_CONF_OCC2(IOCCL2,IOCCL,NACOB,IOCC,1)
*. Reordered packed in IOCCL to reordered accumulated in IOCCL2
C                REFORM_PACK_TO_ACC_CONF(IP_CONF,IA_CONF,IWAY,NOCAB,NACOB)
            CALL REFORM_PACK_TO_ACC_CONF(IOCCL,IOCCL2,1,IOCC,NACOB)
            IM_IN = IS_IACC_CONF_IN_MINMAX_SPC
     &              (IOCCL2,MINAC,MAXAC,NACOB)
            IF(IM_IN.EQ.1) THEN
              CALL REFORM_PACK_TO_ACC_CONF(ICONF_OCC(ICNBS),IOCCL,1,
     &             IOCC,NACOB)
              IF(NTEST.GE.1000) 
     &        WRITE(6,*) ' Configuration in minmax space'
              DO JCSF = 1, NNCSF
               ICSF = ICSF + 1
               ISEL = ISEL + 1
               IADR(ISEL) = ICSF
              END DO
              IF(INCLUDE_CONFS.EQ.1) THEN
                CALL ICOPVE
     &          (ICONF_OCC(ICNBS),ICONF_OCC_SEL(IBCONF_OCC_SEL),IOCC)
                IBCONF_OCC_SEL = IBCONF_OCC_SEL + IOCC
                IF(NTEST.GE.1000) WRITE(6,*) ' IOCC, IBCONF_OCC_SEL = ', 
     &                                         IOCC, IBCONF_OCC_SEL
                NCONF_SEL = NCONF_SEL + 1
                IF(NTEST.GE.1000) WRITE(6,*)  ' NCONF_SEL = ', NCONF_SEL
                NOP_CONF_SEL(NCONF_SEL) = IOPEN
              END IF
            ELSE
              ICSF = ICSF + NNCSF
            END IF ! IM_IN
          END DO !loop over configurations
*. Update pointer
          ICNBS0 = ICNBS0 + NNCNF*IOCC
        END DO ! Loop over IOPEN
      END DO ! Loop over occupation classes
      NSEL = ISEL
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Number of selected CSFs and Confs', NSEL, NCONF_SEL
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Addresses of CSFs in space '
        CALL IWRTMA(IADR,1,NSEL,1,NSEL)
        WRITE(6,*) 
        WRITE(6,*) ' And the selected configurations:'
        IB = 1
        DO ICNF = 1, NCONF_SEL
          IOPEN = NOP_CONF_SEL(ICNF)
          IOCC = IOPEN + (NACTEL-IOPEN)/2
          CALL WRT_CONF(ICONF_OCC_SEL(IB),IOCC)
          IB = IB + IOCC
        END DO
      END IF
*
      IF(NSEL.NE.NELMNT.AND.NELMNT.GT.0) THEN
        WRITE(6,*) 
     &  ' Expected and actual number of CSFs differs: ', NELMNT,NSEL
          STOP 'Expected and actual number of CSFs differs'
      END IF

*
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'GTMNCS')
      CALL LUCIAQEXIT('GTMNCS')
      RETURN
      END
      SUBROUTINE REFORM_PACK_TO_ACC_CONF(IP_CONF,IA_CONF,IWAY,
     &           NOCOB,NACOB)
*
* Reform between packed and accumulated form of configurations
*
* IWAY = 1 PACK => Accumulated
*      = 2 Accumulated => Packed
*
*. Jeppe Olsen, July 22, 2013
*
      INCLUDE 'implicit.inc'
*. Input and output
      INTEGER IP_CONF(NOCOB),IA_CONF(NACOB)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Info from REFORM_PACK_TO_ACC_CONF'
        WRITE(6,*) ' ================================='
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' Packed to accumulated '
        ELSE IF (IWAY.EQ.2) THEN
         WRITE(6,*) ' Accumulated to packed '
        END IF
      END IF
*
      IF(IWAY.LT.1.OR.IWAY.GT.2) THEN
        WRITE(6,*) 
     &  ' REFORM_PACK_TO_ACC_CON: Illegal value of IWAY: ', IWAY
        STOP 'REFORM_PACK_TO_ACC_CON: Illegal value of IWAY '
      END IF
*
      IF(IWAY.EQ.1) THEN
*
*. Packed to accumulated
*
        IZERO = 0
        CALL ISETVC(IA_CONF,IZERO,NACOB)
        DO IOC = 1, NOCOB
          IOB = ABS(IP_CONF(IOC))
          IF(IP_CONF(IOC).GT.0) THEN
            IEL = 1
          ELSE
            IEL = 2
          END IF
          DO JOB = IOB, NACOB
            IA_CONF(JOB) = IA_CONF(JOB) + IEL
          END DO
        END DO
      ELSE
*
* Accumulated => packed
*
        IOC = 0
        DO IOB = 1, NACOB
          IF(IOB.EQ.1) THEN
            IEL = IA_CONF(1)
          ELSE
            IEL = IA_CONF(IOB)-IA_CONF(IOB-1)
          END IF
          IF(IEL.EQ.1) THEN
            IOC = IOC + 1
            IP_CONF(IOC) = IOB
          ELSE IF (IEL.EQ.2) THEN
            IOC = IOC + 1
            IP_CONF(IOC) = -IOB 
          END IF
        END DO
      NOCOB = IOC
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Packed configuration: '
        CALL IWRTMA(IP_CONF,1,NOCOB,1,NOCOB)
        WRITE(6,*)
        WRITE(6,*) ' Accumulated configuration '
        CALL IWRTMA(IA_CONF,1,NACOB,1,NACOB)
      END IF
*
      RETURN
      END

c $Id$
