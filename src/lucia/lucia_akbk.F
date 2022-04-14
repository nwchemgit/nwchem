      SUBROUTINE AKBKCI(JCMBSPC,IPRINT,
     &           EREF,ERROR_NORM_FINAL,CONV_F)
*
* Outer routine for performing AKBK calculation
*
*. Jeppe Olsen, Febr. 2012, Geneva
*
*                    Oct. 30, 2012; Jeppe Olsen; call to Z_BLKFO changed
* Last modification; July 2013; Jeppe Olsen; AKSIGMA2 path added
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'akbk.inc'
      INCLUDE 'cshift.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cstate.inc'  
      INCLUDE 'gasstr.inc'
      INCLUDE 'clunit.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'cprnt.inc'
      include 'cicisp.inc'
*
      LOGICAL CONV_F, CONV_LOCAL
      EXTERNAL AKBKSIGMA, AKBKSIGMA2
      real(8) KATI(4) !VOG
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'AKBKCT')
      CALL LUCIAQENTER('AKBKCT')
 
      if (ga_nodeid().eq.0) then
      WRITE(6,*) ' **************************************'
      WRITE(6,*) ' *                                    *'
      WRITE(6,*) ' * AKBK optimization control entered  *'
      WRITE(6,*) ' *                                    *'
      WRITE(6,*) ' * Version +0.1, Jeppe Olsen, July. 13*'
      WRITE(6,*) ' **************************************'
      end if
*
      NTEST = 1000 !VOG 000
      NTEST = MAX(NTEST,IPRINT)
*
* Test that CSFs are in use
      IF((NOCSF.EQ.1).and.(ga_nodeid().eq.0)) THEN
       WRITE(6,*) ' AKBK calculation with SD is under construction'
cVOG        WRITE(6,*) ' AKBK calculation with determinants is not enabled'
cVOG        WRITE(6,*) ' Please remove NOCSF keyword '
cVOG        STOP       ' AKBK calculation with determinants is not enabled'
      END IF
      IF((ICNFBAT.EQ.1.OR.ICISTR.EQ.1).and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' AKBK calculation with GAs is under construction'
cVOG        WRITE(6,*) ' AKBK calculation does not work with incore options'
cVOG        WRITE(6,*) ' Please make sure ICNFBAT .ge. 2 and ICISTR .ge. 2 '
cVOG        STOP       ' AKBK calculation does not work with incore options'
      END IF
*
cVOG  IAKBK_MET = 1
      if (ga_nodeid().eq.0) then
      WRITE(6,*) ' Allowed number of outer iterations ', MAXITM
      WRITE(6,*) ' Allowed number of inner iterations ', MAXIT 
*
cVOG  All are written in GAs
c     IF(IAKBK_MET.EQ.1) THEN
c       WRITE(6,*) ' Procedure with Q-vectors on DISC '
c     ELSE
c       WRITE(6,*) ' Procedure without Q-vectors on DISC '
c     END IF
*
      WRITE(6,*) 
     &' Convergence threshold for macroiterations: ', THRES_E
      end if
*
*. Memory for summarizer of information
*
      NITEM = 2
      LEN_SUMMARY = NITEM*(MAXITM+1)
cVOG  Important for AKBKMCSCF
cVOG  CALL MEMMAN(KL_SUMMARY,LEN_SUMMARY,'ADDL  ',2,'SUMMRY')
      IPSPC = JCMBSPC - 1
      IQSPC = JCMBSPC
      if (ga_nodeid().eq.0) then
      WRITE(6,*) ' P-space for expansion ', IPSPC
      WRITE(6,*) ' Q-space for expansion ', IQSPC
      end if
*. Number of CSFs in P and Q
      NPCSF = NCSF_PER_SYM_GN(IREFSM,IPSPC)
      NQCSF = NCSF_PER_SYM_GN(IREFSM,IQSPC)
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Initial energy = ', EREF
cVOG    WRITE(6,*) ' Number of CSFs in P-space = ', NPCSF
cVOG    WRITE(6,*) ' Number of CSFs in PQ-space = ', NQCSF
      END IF
      EREF_INI = EREF
*
* ==================================
* 0: Initial preparation of classes
* ==================================
* 
*
*. The occupation classes of the P-space are those with basespace IPSPC
*. identify and save thise
*
      CALL MEMMAN(KLOCCLS_P,NOCCLS_MAX,'ADDL  ',1,'POCCLS')
      CALL OCCLS_IN_CISPACE(NOCCLS_P,int_mb(KLOCCLS_P),
     &       NOCCLS_MAX,int_mb(KIOCCLS),NGAS,
     &       LCMBSPC(IPSPC),ICMBSPC(1,IPSPC),IGSOCCX,IPSPC)
      NCOCCLS = NOCCLS_P
*
*. And the initial identification of the classes in the Q-space
*
      CALL MEMMAN(KLOCCLS_Q,NOCCLS_MAX,'ADDL  ',1,'QOCCLS')
      CALL OCCLS_IN_CISPACE(NOCCLS_Q,int_mb(KLOCCLS_Q),
     &       NOCCLS_MAX,int_mb(KIOCCLS),NGAS,
     &       LCMBSPC(IQSPC),ICMBSPC(1,IQSPC),IGSOCCX,IQSPC)
      NSOCCLS = NOCCLS_Q
*. Well the Q-classes may contain terms that are not connected
*. to the P-classes by double excitations. Eliminate these
*
      CALL MEMMAN(KLRCONNECT,NOCCLS_MAX,'ADDL  ',1,'QOCCLS')
C     EXCIT_OCCLS_TO_OCCLS(NOCC_L,IOCC_L,NOCC_R,IOCC_R,
C    &           IOCCLS,NGAS,NEXCIT,NRCONNECT,IRCONNECT)
      MAX_EXCIT = 117
      if (ga_nodeid().eq.0)
     &WRITE(6,*) ' Max excitation level beween P and Q set to',
     &            MAX_EXCIT
      CALL EXCIT_OCCLS_TO_OCCLS(NOCCLS_P,int_mb(KLOCCLS_P),
     &    NOCCLS_Q,int_mb(KLOCCLS_Q),int_mb(KIOCCLS),NGAS,MAX_EXCIT,
     &    NRCONNECT,int_mb(KLRCONNECT))
      IF(NRCONNECT.LT.NOCCLS_Q) THEN
        if (ga_nodeid().eq.0) then
        WRITE(6,*) ' Q-space reduced to interacting occlasses'
        WRITE(6,*) ' Original and new number of Q occlasses ',
     &               NOCCLS_Q, NRCONNECT
        endif
        CALL ICOPVE(int_mb(KLRCONNECT),int_mb(KLOCCLS_P),NRCONNECT)
*. Should also be copied to Sigma occ classes?
      END IF
*
* Info on batches
*
      !
      ! VOG: My quick-and-dirty way to solve the problem
      !
cVOG  LBLOCK  = MXSOOB 
cVOG  LBLOCK  = NSD_FOR_OCCLS_MAX
cVOG  IF (NOCSF.EQ.1) CALL IAIBCM(1,LBLOCK)
      !
      ! VOG: Solution based on Jeppe's suggestions
      !
      CALL MAX_NOPEN_NOCCLS_FOR_CISPAC(IGSOCC(1,1),IGSOCC(1,2),
     &     MINOP,NOPEN_MAX,NOCCLS_MAX)
      write(*,*) 'Number of occupation classes in compound space',
     &            NOCCLS_MAX
      ! Generate the actual occupations of the various GASpaces that
      ! defines each occupation classes
      CALL GEN_OCCLS_FOR_CISPAC(IGSOCC(1,1), IGSOCC(1,2), -1,
     &     MINOP,int_mb(KIOCCLS),NOCCLS_MAX)
      ! Obtain the AB supergroups that constitute the various occupation
      ! classes
*. Obtain mapping from occupation classes to supergroups
      IATP = 1
      IBTP = 2
      NASPGP = NSPGPFTP(IATP)
      NBSPGP = NSPGPFTP(IBTP)
      IB_A = IBSPGPFTP(IATP)
      IB_B = IBSPGPFTP(IBTP)
*. First: offsets
      write(*,*) ' Before work -LBLOCK', LBLOCK
      CALL ABSPGP_TO_OCCLS(
     &     ISPGPFTP(1,IB_A),NASPGP,ISPGPFTP(1,IB_B),NBSPGP,
     &     NGAS,IGSOCC,NELFGP,
     &     NOCCLS_MAX,int_mb(KIOCCLS),N_ABSPGP_TOT2,
     &     int_mb(KIABSPGP_FOR_OCCLS),int_mb(KNABSPGP_FOR_OCCLS),
     &     int_mb(KIBABSPGP_FOR_OCCLS),1,ga_nnodes(),N_ABSPGP_MAX)
*. And then the actual mappings
      CALL ABSPGP_TO_OCCLS(
     &     ISPGPFTP(1,IB_A),NASPGP,ISPGPFTP(1,IB_B),NBSPGP,
     &     NGAS,IGSOCC,NELFGP,
     &     NOCCLS_MAX,int_mb(KIOCCLS),N_ABSPGP_TOT2,
     &     int_mb(KIABSPGP_FOR_OCCLS),int_mb(KNABSPGP_FOR_OCCLS),
     &     int_mb(KIBABSPGP_FOR_OCCLS),0,ga_nnodes(),N_ABSPGP_MAX)
      ! Order CI-vector, so AB supergroups so all supergroups belonging
      ! to the same occupation class are stored consecutively and one batch consists
      ! of one occupation class:
*
      ICOMP = 0
* Start be letting P be C-space and Q be S-space
      NCOCCLS = NOCCLS_P
      ICSPC = IPSPC
      ICSM = IREFSM
      ILTEST = 3006
      CALL ICOPVE(int_mb(KLOCCLS_P),int_mb(KCIOCCLS_ACT),NGAS*NOCCLS_P)
      CALL Z_BLKFO_FOR_CISPACE(ICSPC,ICSM,LBLOCK,ICOMP,
     &     IPRINT,NCBLOCK,NCBATCH,
     &     int_mb(KCIOIO),int_mb(KCBLTP),NCOCCLS_ACT,
     &     int_mb(KCIOCCLS_ACT),
     &     int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &     int_mb(KCIBT),
     &     int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
      NBLOCK_P = NCBLOCK
      NBATCH_P = NCBATCH
      NCBLK = NCBLOCK
      NCM_P = IELSUM(int_mb(KCLBLK),NBLOCK_P)
*
      NSOCCLS = NOCCLS_Q
      ISSPC = IQSPC
      ISSM = IREFSM
      ILTEST = 3006
      ! Trick code to use PART_CIV_OCC instead of PART_CIV for
      ! partitioning of vector
      NOCSF = 0
      CALL Z_BLKFO_FOR_CISPACE(ISSPC,ISSM,LBLOCK,ICOMP,
     &     IPRINT,NSBLOCK,NSBATCH,
     &     int_mb(KSIOIO),int_mb(KSBLTP),NSOCCLS_ACT,
     &     int_mb(KSIOCCLS_ACT),
     &     int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSLBLK),int_mb(KSI1BT),
     &     int_mb(KSIBT),
     &     int_mb(KSNOCCLS_BAT),int_mb(KSIBOCCLS_BAT),ILTEST)
      NOCSF = 1
      NBLOCK_Q = NSBLOCK
      NBATCH_Q = NSBATCH
      NSBLK = NSBLOCK
*      if (IPRINT.ge.10.and.ga_nodeid().eq.0) then
      write(*,*) ' After work -LBLOCK', LBLOCK
      NSD_FOR_OCCLS_MAX = LBLOCK
      if (ga_nodeid().eq.0) then
      WRITE(6,'(a,2I9)')  
     & ' Number of blocks and batches in P  ',NBLOCK_P, NBATCH_P
      WRITE(6,'(a,2I9)')  
     & ' Number of blocks and batches in PQ ',NBLOCK_Q, NBATCH_Q
      endif
*
      IF (NOCSF.EQ.1) THEN
*
* AKBK - NWChem version. No CSFs, just SDs
*
* Allocate three vectors: KCOMVEC1, KCOMVEC2, KVEC3
* Dimensions of KCOMVEC1, KCOMVEC2: P-space
*
       NSDP  = INT(XISPSM(IREFSM,1))
       NSDPQ = INT(XISPSM(IREFSM,2))
       NSDQ  = NSDPQ - NSDP
       NVAR  = NSDP
       IF ( .not. ga_create(mt_dbl, NSDP, 1, 'COMVEC1', 0, 0, KCOMVEC1))
     $    CALL ERRQUIT('AKBK: ga_create COMVEC1 failed', NSDP, GA_ERR)
       CALL ga_zero(KCOMVEC1)
!       IF (dra_read_section(.false., KVEC1P, 1, NSDP, 1, 1,
!     &                               LUC,    1, NSDP, 1, 1,ITASK).ne.0)
!     &    CALL ERRQUIT('dra error readin P vector',ITASK, 666)
!       IF (dra_wait(ITASK).ne.0) CALL ERRQUIT('dra read err',ITASK,666)
       IF ( .not. ga_create(mt_dbl, NSDP, 1, 'COMVEC2', 0, 0,KCOMVEC2))
     $    CALL ERRQUIT('AKBK: ga_create COMVEC2 failed', NSDP, GA_ERR)
       CALL ga_zero(KCOMVEC2)
*
* Allocate KVEC3
*
*       IF ( .not. ga_create( mt_dbl, NSDQ, 1, 'VEC3', 0, 0, KVEC3) )
*     $    CALL ERRQUIT('AKBK: ga_create VEC3 failed', NSDQ, GA_ERR)
*       CALL ga_zero(KVEC3)
*
*LOCAL
*
      CALL GET_3BLKS(KVEC1,KVEC2,KC2)
      KVEC3 = KC2
      KVEC1P = KVEC1
      KVEC2P = KVEC2
*      CALL MEMMAN(KVEC1P, NSDP,'ADDS  ',2,'VEC1P ')
*      CALL MEMMAN(KVEC2P, NSDP,'ADDS  ',2,'VEC2P ')
*      CALL MEMMAN(KC2,    NSDQ,'ADDS  ',2,'VEC3  ') !!!ACHTUNG!!!!
*
* Analyze P-space expansion (Do I need it?)
*
cVOG   CALL GASANA(dbl_mb(KVEC1P),NBLOCK,NSDP,NSDP,LUC,ICISTR)
*
* Allocate Q-space vector
*
!       IF ( .not. ga_create( mt_dbl, NSDQ, 1, 'VECQ', 0, 0, KVECQ) )
!     $    CALL ERRQUIT('AKBK: ga_create VECP failed', NSDQ, GA_ERR)
!       CALL ga_zero(KVECQ)
!       CALL ga_copy_patch('N', KVEC1P, 1, NSDP, 1, 1,
!     &                         KVECQ,  1, NSDP, 1, 1)
*
      IF(IAKBK_MET.EQ.2) THEN
*.Allocate extra block for incore
        write(*,*) 'NSD_FOR_OCCLS_MAX', NSD_FOR_OCCLS_MAX
        IF (.not.ga_create(mt_dbl,NSD_FOR_OCCLS_MAX,1,'VEC4',0,0,KVEC4))
     $    CALL ERRQUIT('AKBK: ga_create VEC4 failed',
     &                 NSD_FOR_OCCLS_MAX, GA_ERR)
        CALL ga_zero(KVEC4)
cVOG    CALL MEMMAN(KVEC4,NSD_FOR_OCCLS_MAX,'ADDL  ',2,'VEC4  ')
        LCSBLK = LENP
        IF(LENP.GT.LCSBLK) THEN
         if (ga_nodeid().eq.0) WRITE(6,*)
     &   ' Length of batch will be increased to accomodate P-expansion'
         LCSBLK = LENP
        END IF
      END IF
*
       IF (NTEST.GE.1000) THEN
        if (ga_nodeid().eq.0) then
        WRITE(6,'(a)')     ' ================================= '
        WRITE(6,'(a)')     '  Dimensions for AKBK Calculation  '
        WRITE(6,'(a)')     ' ================================= '
        WRITE(6,'(a,i12)') ' Number of CI Spaces    ', NCMBSPC
        WRITE(6,'(a,i12)') ' Dimension of P space   ', NSDP
        WRITE(6,'(a,i12)') ' Dimension of Q space   ', NSDQ
        WRITE(6,'(a,i12)') ' Dimension of PQ space  ', NSDPQ
        WRITE(6,'(a,i12)') ' Max Num of Occ. Classes', NOCCLS_MAX
        WRITE(6,'(a)')     ' Print initial vector:'
        endif
        if (ga_nodeid().eq.0) then
        CALL ga_print(KCOMVEC1)
        WRITE(6,'(a)')     ' ================================= '
        endif
       END IF
*
      ELSE
*
* AKBK - Olsen's version with CSFs
*
*
*. Allocate space: 
*
*. Using the low-disc version, the arrays should be able to hold the 
*. P-expansion (with one block expanded), so
*. The below is a bit sloppy and memory could be reduced
      IF(PSSIGN.EQ.0.0D0) THEN
       LENP = NCM_P
      ELSE
       LENP = 2*NCM_P
      END IF
*
      IF(IAKBK_MET.EQ.2.AND.LENP.GT.LCSBLK) THEN
        WRITE(6,*) 
     &  ' Length of batch will be increased to accomodate P-expansion'
        LCSBLK = LENP
      END IF
      IF(IAKBK_MET.EQ.2) THEN
*.Allocate extra block for incore
        CALL MEMMAN(KVEC4,NSD_FOR_OCCLS_MAX,'ADDL  ',2,'VEC4  ')
      END IF
      CALL GET_3BLKS(KVEC1,KVEC2,KC2)
      KVEC3 = KC2
      KVEC1P = KVEC1
      KVEC2P = KVEC2
*
      I_DO_COMHAM = 0
      IF(I_DO_COMHAM.EQ.1) THEN
         WRITE(6,*) ' Complete Hamiltonian matrix will be constructed '
         NVAR = NPCSF
         CALL MEMMAN(KLHMAT,NVAR*NVAR,'ADDL  ',2,'HMAT  ')
         ECOREL = 0.0D0
C             COMHAM(H,NVAR,NBLOCK,LBLOCK,VEC1,VEC2)
         CALL COMHAM(WORK(KLHMAT),NVAR,NCBLOCK,int_mb(KCLBLK),
     &               WORK(KVEC1P),WORK(KVEC2P),ECOREL)
         STOP ' Enforced stop after COMHAM'
      END IF
*
* Preconditioner in P space
*
      LBLK = - 1
*. Obtain Preconditioner in CSF-basis: form determined by IH0_CSF
c     WRITE(6,*) ' IH0_CSF = ', IH0_CSF
      IF(IH0_CSF.EQ.1) THEN
*
*  ================
*. Averaged Diagonal
*  ================
*
*. Obtain diagonal over determinants
cVOG        CALL REWINO(LUDIA)
        I12 = 2
        SHIFT = ECORE
        IUSE_EXP = 1
        CALL GASDIAT(WORK(KVEC2),LUDIA,SHIFT,ICISTR,I12,
     &               int_mb(KCBLTP),NCBLOCK,int_mb(KCIBT),IUSE_EXP)
        CALL FILEMAN_MINI(LUSC,'ASSIGN')
        CALL CSDIAG(WORK(KVEC1P),WORK(KVEC1P),
     &       NCONF_PER_OPEN(1,IREFSM),MAXOP,IREFSM,
     &       WORK(KSDREO_I(IREFSM)),NPCMCNF,NPCSCNF,IPRCSF,
     &       ICNFBAT,NCOCCLS_ACT,int_mb(KCIOCCLS_ACT),
     &       int_mb(KCLBT),LUDIA,LUSC)
        CALL COPVCD(LUSC,LUDIA,WORK(KVEC2P),1,-1)
        CALL FILEMAN_MINI(LUSC,'FREE  ')
      ELSE
*
* ==============================
*. Use exact diagonal or blocks
* ==============================
*
C?      WRITE(6,*) ' Ecore before GET_CSF.... ', ECORE
        NAEL = NELFTP(1) 
        NBEL = NELFTP(2)  
*
        CALL GET_CSF_H_PRECOND(NCONF_PER_OPEN(1,IREFSM),
     &        int_mb(KICONF_OCC(IREFSM)),KVEC1P,ECORE,
     &        LUDIA,NCOCCLS_ACT,int_mb(KCIOCCLS_ACT),IREFSM)
C      GET_CSF_H_PRECOND(NCONF_FOR_OPEN,ICONF_OCC,H0,
C     &           ECORE,LUDIA,NOCCLS_SPC,IOCCLS_SPC,ISYM)

C        GET_CSF_H_PRECOND(NEL,NCONF_FOR_OPEN,ICONF_OCC,H0,ECORE,
C    &   IFORM)
      END IF! IH0_CSF switch
      END IF! NOCSF - VOG
*
      MINST = 1
*. Since we are starting our with the solution from the previous space
      IRESTR = 1
      IF(IRESTR.EQ.0) THEN
        INICI = 0
      ELSE
        INICI = -1
      END IF
      LBLK = - 1
      NPRDET = 0
      IPERTOP = 0
      EADD = 0.0D0
*
      KSBEVC = 1
      KSBIDT = 1
      KH0 = 1
      KH0SCR = 1
*
      KLCLSC = 1
      KLCLSE = 1
      KLCLSCT = 1
      KLCLSET = 1
      KLCLSA = 1 
      KLCLSL = 1 
      KLCLSLR = 1 
      KLBLKA = 1
      KLCLSD = 1
      KLCLSDT = 1 
      KLCLSA2 = 1
      KLCBASC = 0
      KLEBASC = 0
      KLBLKCLS = 0
      KLCLSG = 0
      KLCLSGT = 0
*
*. Loop over macroiterations - iterations with fixed energy
*
      EREF_OLD = 3006.56
      CONV_F = .FALSE.
      DO IMAC = 1, MAXITM 
       if (ga_nodeid().eq.0) then
       WRITE(6,*) ' Info from macroiteration ', IMAC
       WRITE(6,*) ' ================================'
       endif
       IF(IMAC.GT.1) EREF_OLD = EREF
       SHIFT = - EREF
       IIUSEH0P = 0
       ILAST = -3006
C?     WRITE(6,*) ' IPRCOND = ', IPRECOND  
C?     WRITE(6,*) ' ICISTR, LBLK = ', ICISTR, LBLK
*. An additional file
       CALL FILEMAN_MINI(LU8,'ASSIGN')
       IF(IAKBK_MET.EQ.1) THEN
       CALL CIEIG5(AKBKSIGMA,INICI,EROOT,KCOMVEC1,KCOMVEC2,
     & MINST,LUDIA,LUC,LUHC,LUSC1,LUSC2,LUSC3,LUSC34,LUSC35,LU8,
     & NVAR,NBLK,NROOT,MXCIV,MAXIT,LUCIVI,IPRDIA,dbl_mb(KSBEVC),
     & NPRDET,dbl_mb(KH0),int_mb(KSBIDT),
     & MXP1,MXP2,MXQ,dbl_mb(KH0SCR),EADD,ICISTR,LBLK,
     & IDIAG,dbl_mb(KVEC3),THRES_E,
     & NBATCH,
     & int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     & int_mb(KCIBT),
     & int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSLBLK),int_mb(KSI1BT),
     & int_mb(KSIBT),
     & INIDEG,E_THRE,C_THRE,
     & E_CONV,C_CONV,ICLSSEL,int_mb(KLBLKCLS),NOCCLS,
     & dbl_mb(KLCLSC),dbl_mb(KLCLSE), dbl_mb(KLCLSCT),dbl_mb(KLCLSET),
     & int_mb(KLCLSA),int_mb(KLCLSL),dbl_mb(KLCLSLR),int_mb(KLBLKA),
     & dbl_mb(KLCLSD),dbl_mb(KLCLSDT),dbl_mb(KLCLSG),dbl_mb(KLCLSGT),
     & ISKIPEI,int_mb(KC2B),int_mb(KLCLSA2),
     & LBLOCK,IROOT_SEL,int_mb(KBASSPC),dbl_mb(KLEBASC),
     & dbl_mb(KLCBASC),NCMBSPC,MULSPCA,IPAT,LPAT,ISPC,NCNV_RT,
     & IPRECOND,IIUSEH0P,MPORENP_E,RNRM,CONV_LOCAL,ISBSPPR_ACT,
     & ILAST)
      ELSE
       CALL CIEIG5(AKBKSIGMA2,INICI,EROOT,KCOMVEC1,KCOMVEC2,
     & MINST,LUDIA,LUC,LUHC,LUSC1,LUSC2,LUSC3,LUSC34,LUSC35,LU8,
     & NVAR,NBLK,NROOT,MXCIV,MAXIT,LUCIVI,IPRDIA,dbl_mb(KSBEVC),
     & NPRDET,dbl_mb(KH0),int_mb(KSBIDT),
     & MXP1,MXP2,MXQ,dbl_mb(KH0SCR),EADD,ICISTR,LBLK,
     & IDIAG,dbl_mb(KVEC3),THRES_E,
     & NBATCH,
     & int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     & int_mb(KCIBT),
     & int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSLBLK),int_mb(KSI1BT),
     & int_mb(KSIBT),
     & INIDEG,E_THRE,C_THRE,
     & E_CONV,C_CONV,ICLSSEL,int_mb(KLBLKCLS),NOCCLS,
     & dbl_mb(KLCLSC),dbl_mb(KLCLSE), dbl_mb(KLCLSCT),dbl_mb(KLCLSET),
     & int_mb(KLCLSA),int_mb(KLCLSL),dbl_mb(KLCLSLR),int_mb(KLBLKA),
     & dbl_mb(KLCLSD),dbl_mb(KLCLSDT),dbl_mb(KLCLSG),dbl_mb(KLCLSGT),
     & ISKIPEI,int_mb(KC2B),int_mb(KLCLSA2),
     & LBLOCK,IROOT_SEL,int_mb(KBASSPC),dbl_mb(KLEBASC),
     & dbl_mb(KLCBASC),NCMBSPC,MULSPCA,IPAT,LPAT,ISPC,NCNV_RT,
     & IPRECOND,IIUSEH0P,MPORENP_E,RNRM,CONV_LOCAL,ISBSPPR_ACT,
     & ILAST)
      END IF
*
      CALL FILEMAN_MINI(LU8,'FREE  ')
       EREF = EROOT(NROOT)
       if (ga_nodeid().eq.0) WRITE(6,'(A,F22.12)') '   Energy: ', EREF
       ERROR_NORM_FINAL = RNRM(NROOT)
       dbl_mb(KL_SUMMARY-1+(IMAC-1)*NITEM+1) = ERROR_NORM_FINAL
       dbl_mb(KL_SUMMARY-1+(IMAC-1)*NITEM+2) = EREF
       if (ga_nodeid().eq.0) WRITE(6,'(A,E8.3)') '   Residual norm ',
     & ERROR_NORM_FINAL
*
       DELTA_E = EREF-EREF_OLD
       IF(ABS(DELTA_E).LE.THRES_E.AND.CONV_LOCAL) THEN
         CONV_F = .TRUE.
         GOTO 1001
       END IF
      END DO! loop over macroiterartions
 1001 CONTINUE
*
      IF(CONV_F) THEN
        if (ga_nodeid().eq.0) 
     &     WRITE(6,*) ' AKBK procedure converged '
        NIT = IMAC
      ELSE
        if (ga_nodeid().eq.0) 
     &     WRITE(6,*) ' AKBK procedure not converged '
        NIT = IMAC - 1
      END IF
*
      if (ga_nodeid().eq.0) then
      WRITE(6,*)
      WRITE(6,*)
     & ' Iter Residual norm        Energy '
      WRITE(6,*)
     & ' ==============================='
      DO IT = 1, NIT
        ENER =  dbl_mb(KL_SUMMARY-1+(IT-1)*NITEM+2) 
        RNRMI =  dbl_mb(KL_SUMMARY-1+(IT-1)*NITEM+1) 
        WRITE(6,'(1X,I3,4X,E8.3,F20.12)') IT, RNRMI, ENER
      END DO
      endif
*
*. Calculate density matrices
*
C?    IF(IDENSI.GE.1)  THEN
C?      WRITE(6,*) ' AKBK_DENSITY will be called '
C?      CALL AKBK_DENSITY(WORK(KRHO1),WORK(KRHO2),WORK(KVEC1),
C?   &       WORK(KVEC2),LUC,EREF,IDENSI)
C            AKBK_DENSITY(RHO1,RHO2,L,R,LUP,E,I12)
C?    END IF
*
      if (.not. ga_destroy(KCOMVEC2)) call
     &     errquit('GASCI: ga_destroy kcomvec2 failed', 0 , GA_ERR)
      if (.not. ga_destroy(KCOMVEC1)) call
     &     errquit('GASCI: ga_destroy kcomvec1 failed', 0 , GA_ERR)
      if (.not. ga_destroy(KVEC4)) call
     &     errquit('GASCI: ga_destroy kvec4 failed',    0 , GA_ERR)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'AKBKCT')
      CALL LUCIAQEXIT('AKBKCT')
*
      RETURN
      END
      SUBROUTINE EXCIT_OCCLS_TO_OCCLS(NOCC_L,IOCC_L,NOCC_R,IOCC_R,
     &           IOCCLS,NGAS,NEXCIT,NRCONNECT,IRCONNECT)
*
* Two sets, L and R, of occupation classes are given in the 
* form of pointers to the list of occupation classes. 
* Obtain those R classes that are connected to atleast one 
* of the L- classes by atleast a NEXCIT-fold excitation
*
*. Jeppe Olsen
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IOCC_L(NOCC_L),IOCC_R(NOCC_R)
      INTEGER IOCCLS(NGAS,*)
*. Output
      INTEGER IRCONNECT(*)
*
      NTEST = 000
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from EXCIT_OCCLS_TO_OCCLS '
        WRITE(6,*) ' ================================='
        WRITE(6,*)
        WRITE(6,*) ' Allowed excitation level = ', NEXCIT
        WRITE(6,*) ' NOCC_R, NOCC_L = ', NOCC_R, NOCC_L
      END IF
*
      NRCONNECT = 0
      DO IIR = 1, NOCC_R
       IR = IOCC_R(IIR)
       IM_CONNECTED = 0
       DO IIL = 1, NOCC_L
        IL = IOCC_L(IIL)
*. Differences in occupations between the classes
        IDIFF = 0
        DO IGAS = 1, NGAS
          IDIFF = IDIFF + ABS(IOCCLS(IGAS,IL)-IOCCLS(IGAS,IR))
        END DO
        IF(IDIFF.LE.2*NEXCIT) IM_CONNECTED = 1
       END DO!loop over IL
       IF(IM_CONNECTED.EQ.1) THEN
         NRCONNECT = NRCONNECT + 1
         IRCONNECT(NRCONNECT) = IR
       END IF
      END DO
*
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' Number of connected R-classes = ', NRCONNECT
      END IF
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' The connected classes '  
       CALL IWRTMA(IRCONNECT,1,NRCONNECT,1,NRCONNECT)
      END IF
*
      RETURN
      END
      SUBROUTINE AKBKSIGMA(CB,HCB,LUC,LUHC)
*. 
*. Obtain AKBK Hamiltonian times a vector on LUC
*
* A P and a Q space has been defined, so the eigenvalue problem reads:
*
* ( H(PP)    H(PQ) ) (CP)        (CP)
* (                ) (  )   = E  (  )
* ( H(QP)    H(QQ) ) (CQ)        (CQ)
* 
* Lowdin partition the matrix to obtain
*
* [(HPP - H(PQ) (H(QQ) - E )^(-1)H(QP)]CP = E CP
*
* HQQ is assumed so simple it may be diagonalized - here
* Diagonal or block diagonal over 
*
*. Input and output vectors are supposed to be in P-space(?)
*
*. It is assumed that ICNFBAT = 2, NOCSF = 0 has been checked outside
*
*. Jeppe Olsen, Geneva, Febr. 2012 (20 hours before take-off)
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cshift.inc'
      INCLUDE 'akbk.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc'
*
      integer CB,HCB, alo(2),ahi(2),blo(2),bhi(2)
      real(8) KATI(4)
*
      LBLK = -1
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',1,'AKBKSI')
      CALL LUCIAQENTER('AKBKSI')
      NTEST = 000
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Information from AKBKSIGMA'
        WRITE(6,*) ' =========================='
        WRITE(6,*)
      END IF 
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' LUC,  LUHC = ', LUC, LUHC
        WRITE(6,*) ' NSDP, NSDQ = ', NSDP, NSDQ
      END IF
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Input vector to AKBKSIGMA '
        CALL ga_print(CB)
!        KATI = 0.0d0
!        KATI(1) = -0.44231
!        KATI(2) =  0.89488
!        KATI(3) = -0.05952
!        KATI(4) =  0.00000
!        WRITE(*,*) (KATI(I), I = 1, 4)
!        WRITE(*,*) KATI
!        CALL ga_put(CB,1,4,1,1,KATI,1)
!!         if (dra_write_section(.false.,CB,1,4,1,1,
!!     &                         LUDIA,1,4,1,1,itask).ne.0)
!!     &        call errquit('dra error',itask,911)
!        WRITE(6,*) ' Input vector to AKBKSIGMA from MOLCAS'
!        CALL ga_print(CB)
cVOG    CALL WRTVCD(CB,LUC,1,LBLK)
      END IF
*
*
*. Scratch files for C and Sigma in determinant basis, HQP CP, (HQQ-E) HQP CP
* 
cVOG  CALL FILEMAN_MINI(LU_Q1,'ASSIGN')
cVOG  CALL FILEMAN_MINI(LU_Q2,'ASSIGN')
cVOG  CALL FILEMAN_MINI(LU_Q3,'ASSIGN')
*
* Allocate three scratch arrays: LU_Q1 : P-dim
*                                LU_Q2 : Q-dim
*                                LU_Q3 : Q-dim
*
      IF ( .not. ga_create( mt_dbl, NSDP, 1, 'LU_Q1', 0, 0, KLU_Q1) )
     $   CALL ERRQUIT('AKBK: ga_create LU_Q1 failed', NSDP, GA_ERR)
      IF ( .not. ga_create( mt_dbl, NSDPQ, 1, 'LU_Q2', 0, 0, KLU_Q2) )
     $   CALL ERRQUIT('AKBK: ga_create LU_Q2 failed', NSDPQ, GA_ERR)
      IF ( .not. ga_create( mt_dbl, NSDPQ, 1, 'LU_Q3', 0, 0, KLU_Q3) )
     $   CALL ERRQUIT('AKBK: ga_create LU_Q3 failed', NSDPQ, GA_ERR)
*
      CALL ga_zero(KLU_Q1)
      CALL ga_zero(KLU_Q2)
      CALL ga_zero(KLU_Q3)
*
*. Calculate HQP CP in PQ space and save on LU_Q2
*
      ICSPC = IPSPC
      ISSPC = IQSPC
cVOG  CALL REWINO(LUC)
cVOG  CALL REWINO(LU_Q2)
cVOG  CALL MV7(CB,HCB,LUC,LU_Q2,XDUM,XDUM)
      CALL MV7(CB,KLU_Q2,0,0,0,0)
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) 
     &  WRITE(6,*) ' H(QP)P CP accomplished'
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
       CALL ga_print(CB)
       CALL ga_print(KLU_Q2)
cVOG   CALL WRTVCD(CB,LU_Q2,1,LBLK)
      END IF
*
*. Well really, we wanted HQP CP only in Q space, remove P space terms, to LU_Q1
*
*. Array giving dimensions of occlass for given sym
*      CALL MEMMAN(KLLOCCLS_SM,NOCCLS_MAX,'ADDL  ',1,'OC_LSM')
*      IF(NTEST.GE.1000) WRITE(6,*) ' ISSM, NIRREP = ', ISSM, NIRREP
*C          EXTRROW(WORK(KLCIBT),8,8,NBLOCK,LBLOCK)
cVOG   CALL EXTRROW(int_mb(KNCS_FOR_OCCLS),ISSM,NIRREP,NOCCLS_MAX,
cVOG  &             int_mb(KLLOCCLS_SM))
cVOG  IF(NTEST.GE.1000) THEN
cVOG    WRITE(6,*) ' Number of CSF of right sym per occls'
cVOG    CALL IWRTMA(int_mb(KLLOCCLS_SM),1,NOCCLS_MAX,1,NOCCLS_MAX)
cVOG    write(*,*) 'NOCCLS_P, int_mb(KLOCCLS_P)',
cVOG &              NOCCLS_P, int_mb(KLOCCLS_P)
cVOG    write(*,*) 'NOCCLS_Q, int_mb(KLOCCLS_Q)',
cVOG &              NOCCLS_Q, int_mb(KLOCCLS_Q)
cVOG  END IF
* Be careful, it's working as SplitCAS at the moment      
      CALL ga_copy_patch('N', KLU_Q2, 1, NSDP, 1, 1,
     &                        KLU_Q1, 1, NSDP, 1, 1)
cVOG  CALL EXP_BLKVEC(LU_Q2,NOCCLS_P,int_mb(KLOCCLS_P),
cVOG &                LU_Q1,NOCCLS_Q,int_mb(KLOCCLS_Q),
cVOG &                int_mb(KLLOCCLS_SM),2,CB,XDUM,1,ICISTR,0)
!*
C     EXP_BLKVEC(LU_IN,NBLK_IN, IBLK_IN,
C    &                      LU_OUT,NBLK_OUT,IBLK_OUT,
C    &                      LBLK,ITASK,VEC,VEC_OUT,IREW,ICISTR,
C    &                      INCORE)
      IF(NTEST.GE.1000) WRITE(6,*) ' HQP CP  calculated'
*
*. Calculate diagonal(HQQ) and save on LU_Q3
*
cVOG CALL CSDIAGM(IREFSM,2,LU_Q3,CB,HCB)
      CALL CSDIAGM(IREFSM,2,0,CB,KLU_Q3)
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
       write(6,*) 'Diagonal HQQ is here'
       call ga_print(KLU_Q3)
      END IF
*
*. (HQQ - E)**-1 HQP CP on LU_Q2
*
      IPRECOND = 1
      LBLK = -1 
      IPNTR = 0
      H0 = 3006.56
      H0SCR = 3006.56
      IZERO = 0
!      NPQDM = NSDPQ
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) 
     &  WRITE(6,*) ' Shift :', SHIFT 
!      CALL H0M1TD(KLU_Q2,KLU_Q3,KLU_Q1,LBLK,IZERO,IPNTR,
!     &     H0,SHIFT,H0SCR,XH0IX,
!     &     IZERO,IZERO,IZERO,CB,HCB,NTEST,IPRECOND)
cVOG  CALL H0M1TD(LU_Q2,LU_Q3,LU_Q1,LBLK,IZERO,IPNTR,
cVOG &     H0,SHIFT,H0SCR,XH0IX,
cVOG &     IZERO,IZERO,IZERO,CB,HCB,NTEST,IPRECOND)
*      H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
      call ga_add_constant(KLU_Q3,SHIFT) 
      alo(1) = NSDP+1
      ahi(1) = NSDPQ
      alo(2) = 1
      ahi(2) = 1
      blo(1) = NSDP+1
      bhi(1) = NSDPQ
      blo(2) = 1
      bhi(2) = 1
      call ga_elem_divide_patch(KLU_Q2, alo, ahi,     
     &                          KLU_Q3, blo, bhi,
     &                          KLU_Q2, alo, ahi)
      alo(1) = 1
      ahi(1) = NSDP
      call nga_zero_patch(KLU_Q2, alo, ahi)
!      call ga_elem_divide(KLU_Q2,KLU_Q3,KLU_Q2)
       IF(NTEST.GE.1000) WRITE(6,*) ' Apres H0M1TD '
       CALL MEMCHK2('APH0M1')
       IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' (HQQ - E)**-1 HQP CP '
cVOG    CALL WRTVCD(CB,LU_Q2,1,LBLK)
        call ga_print(KLU_Q2)
       END IF
*
*. Multiply (HQQ - E)**-1 HQP CP (LU_Q2) with HQP and save on LU_Q3
*
      ISSPC = IPSPC
      ICSPC = IQSPC
      XDUM = 3006.56
cVOG  CALL REWINO(LU_Q3)
cVOG  CALL MV7(CB,HCB,LU_Q2,LU_Q3,XDUM,XDUM)
      CALL MV7(KLU_Q2,KLU_Q3,0,0,0,0)
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' HPQ (HQQ - E)**-1 HQP CP  '
cVOG    CALL WRTVCD(CB,LU_Q3,1,LBLK)
        call ga_print(KLU_Q2)
        call ga_print(KLU_Q3)
      END IF
*
* Obtain HPP CP and save on LU_Q1
*
      ISSPC = IPSPC
      ICSPC = IPSPC
cVOG  CALL REWINO(LU_Q1)
cVOG  CALL MV7(CB,HCB,LUC,LU_Q1,XDUM,XDUM)
      CALL MV7(CB,KLU_Q1,0,0,0,0)
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' HPP CP: '
cVOG    CALL WRTVCD(CB,LU_Q1,1,LBLK)
        call ga_print(KLU_Q1)
      END IF
*
* And the grand finale: 
*
* [(HPP - H(PQ) (H(QQ) - E )^(-1)H(QP)]CP
      ONE = 1.0D0
      ONEM = -1.0D0
C VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      LBLK = -1
cVOG  CALL VECSMD(CB,HCB,ONE,ONEM,LU_Q1,LU_Q3,LUHC,1,LBLK)
!      CALL VECSMD(CB,HCB,ONE,ONEM,LU_Q1,LU_Q3,LUHC,0,LBLK)
!PUT a ga operation gia KLU_Q1 - KLU_Q3 - NOT CHECKED!
!      call ga_add(ONE,KLU_Q1,ONEM,KLU_Q3,CB)
!      call ga_copy(CB,HCB)
!*** Be Sure about LU_Q3! It shoud be P-dim when you reach this point!
      call ga_add_patch(ONE,  KLU_Q1, 1, NSDP, 1, 1,
     &                  ONEM, KLU_Q3, 1, NSDP, 1, 1,
     &                           HCB, 1, NSDP, 1, 1)
*
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
       WRITE(6,*) ' Sigma vector from AKBKSIGMA:'
cVOG   CALL WRTVCD(CB,LUHC,1,LBLK)
       call ga_print(HCB)
      END IF
*
cVOG  CALL FILEMAN_MINI(LU_CDET,'FREE  ')
cVOG  CALL FILEMAN_MINI(LU_SDET,'FREE  ')
cVOG  CALL FILEMAN_MINI(LU_Q1,'FREE  ')
cVOG  CALL FILEMAN_MINI(LU_Q2,'FREE  ')
cVOG  CALL FILEMAN_MINI(LU_Q3,'FREE  ')
*
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) 
     &  WRITE(6,*) ' Leaving AKBKSIGMA'
C?    WRITE(6,*) ' Jeppe forced me to stop.... '
C?    STOP ' Jeppe forced me to stop.... '
*
      if (.not. ga_destroy(KLU_Q1)) call
     &     errquit('AKBK: ga_destroy KLU_Q1 failed', 0 , GA_ERR)
      if (.not. ga_destroy(KLU_Q2)) call
     &     errquit('AKBK: ga_destroy KLU_Q2 failed', 0 , GA_ERR)
      if (.not. ga_destroy(KLU_Q3)) call
     &     errquit('AKBK: ga_destroy KLU_Q3 failed', 0 , GA_ERR)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',1,'AKBKSI')
      CALL LUCIAQEXIT('AKBKSI')
*
      RETURN
      END
      SUBROUTINE CSDIAGM(ISM,ICS,LUDIA_A,VECUT,VECSCR)
*
* Set up the CSF preconditioner in CI space ISPC 
* and save on VECUT or LUDIA_A. It is assumed that the relevant block-and
* class info has been obtained and is stored in the KC*, KS*, or 
* KE arrays * (ICS=1,2,3, respectively)
*
*. Jeppe Olsen, Febr. 2012, Geneva, 8 hours before take off.
*  Last modification; Jeppe Olsen; July 29, 2013; ICISTR = 1 and ICS = 3 option added
*
* Incore requires: ICISTR = 1 + ICNFBAT = 1
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'ceblk.inc'
      INTEGER VECUT
*
      CALL LUCIAQENTER('CSDIAM')
*
      NTEST = 000
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN 
        WRITE(6,*) ' Entering CSDIAGM '
        WRITE(6,*) ' ================='
        WRITE(6,*) ' IH0_CSF, ICS = ', IH0_CSF, ICS
      END IF
*
      IF(ICS.NE.1.AND.ICS.NE.2.AND.ICS.NE.3) THEN
        WRITE(6,*) ' Illegal value of ICS in CSDIAGM = ', ICS
        STOP       ' Illegal value of ICS in CSDIAGM  '
      END IF
*
      IF (NOCSF.EQ.1) THEN
*
* AKBK - NWChem version.
*
       LUSC = -1
       I12 = 2
       SHIFT = ECORE
       IUSE_EXP = 1
       KXBLTP = KSBLTP
       KXIBT  = KSIBT
       NXBLOCK = NSBLK
       CALL GASDIAT(VECSCR,LUSC,SHIFT,ICISTR,I12,
     &              int_mb(KXBLTP),NXBLOCK,int_mb(KXIBT),IUSE_EXP)
*
      ELSE
*
* Original version
*
      IF(IH0_CSF.EQ.1) THEN
*
*  ================
*. Averaged Diagonal
*  ================
*
*. Obtain diagonal over determinants on a scratch file
*
        IF(ICISTR.NE.1) THEN
          CALL FILEMAN_MINI(LUSC,'ASSIGN')
          CALL REWINO(LUSC)
        ELSE
          LUSC = -1
        END IF
*
        I12 = 2
        SHIFT = ECORE
        IUSE_EXP = 1
*
        IF(ICS.EQ.1) THEN
          KXBLTP = KCBLTP
          KXIBT  = KCIBT
          NXBLOCK = NCBLK
        ELSE IF(ICS.EQ.2) THEN
          KXBLTP = KSBLTP
          KXIBT  = KSIBT
          NXBLOCK = NSBLK
        ELSE IF (ICS.EQ.3) THEN
*. use info from ceblk.inc
          KXBLTP = KEBLTP
          KXIBT = KEIBT
          NXBLOCK = NEBLK
        ELSE
          WRITE(6,*) ' Illegal value of ICS in CSDIAGM', ICS
          STOP       ' Illegal value of ICS in CSDIAGM'
        END IF
        IUSE_EXP = 1
        IF(NTEST.GE.1000)
     &  WRITE(6,*) ' NCBLK, NSBLK, NXBLOCK = ', NCBLK, NSBLK, NXBLOCK
        CALL GASDIAT(VECSCR,LUSC,SHIFT,ICISTR,I12,
     &               int_mb(KXBLTP),NXBLOCK,int_mb(KXIBT),IUSE_EXP)
*
* Average diagonal in each configuration
*
        IF(ICS.EQ.1) THEN
          NXIOCCLS_ACT = NCOCCLS
          KXIOCCLS_ACT = KCIOCCLS_ACT
          KXLBT = KCLBT
        ELSE IF(ICS .EQ. 2) THEN
          NXIOCCLS_ACT = NSOCCLS
          KXIOCCLS_ACT = KSIOCCLS_ACT
          KXLBT = KSLBT
        ELSE IF (ICS. EQ. 3) THEN
          NXIOCCLS_ACT = NEOCCLS
          KXIOCCLS_ACT = KEIOCCLS_ACT
          KXLBT = KELBT
        END IF
*
        CALL CSDIAG(VECUT,VECSCR,
     &       NCONF_PER_OPEN(1,ISM),MAXOP,ISM,
     &       dbl_mb(KSDREO_I(ISM)),NPCMCNF,NPCSCNF,IPRCSF,
     &       ICNFBAT,NXIOCCLS_ACT,dbl_mb(KXIOCCLS_ACT),
     &       int_mb(KXLBT),LUSC,LUDIA_A)
        IF(ICISTR.NE.1) CALL FILEMAN_MINI(LUSC,'FREE  ')
      ELSE
*
* ==============================
*. Use exact diagonal or blocks
* ==============================
*
        IF(ICS.EQ.1) THEN
          NXIOCCLS_ACT = NCOCCLS
          KXIOCCLS_ACT = KCIOCCLS_ACT
        ELSE 
          NXIOCCLS_ACT = NSOCCLS
          KXIOCCLS_ACT = KSIOCCLS_ACT
        END IF
*
C?      WRITE(6,*) ' Ecore before GET_CSF.... ', ECORE
        CALL GET_CSF_H_PRECOND(NCONF_PER_OPEN(1,ISM),
     &        int_mb(KICONF_OCC(ISM)),VECUT,ECORE,
     &        LUDIA_A,NXIOCCLS_ACT,dbl_mb(KXIOCCLS_ACT),ISM)
      END IF! IH0_CSF switch
      END IF! NOCSF - VOG
      CALL LUCIAQEXIT('CSDIAM')
*
      RETURN
      END
      SUBROUTINE GET_CQ_FROM_CP(LUCP,LUCQ,E,CB,HCB)
*
* Obtain CQ = - (HQQ - E)**(-1) CP
*
* With CP stored on LUCP and result stored at LUCQ
* ( Jeppe is in his systematic mode today....)
*
* Jeppe Olsen, Aug. 23, 2012, in Minneapolis
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cshift.inc'
      INCLUDE 'akbk.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc'
*. Scratch blocks required to hold largest blocks..
      DIMENSION CB(*), HCB(*)
*
      LBLK = -1
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',1,'GET_CQ')
      CALL LUCIAQENTER('GET_CQ')
      NTEST = 0000
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Information from GET_CQ'
        WRITE(6,*) ' =========================='
        WRITE(6,*)
        WRITE(6,*) ' LUCP, LUCQ = ', LUCP, LUCQ
        WRITE(6,*) ' IPSPC, IQSPC = ', IPSPC, IQSPC
      END IF
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Input vector to GET_CQ_FROM_CP '
        CALL WRTVCD(CB,LUCP,1,LBLK)
      END IF
*
*
*. Scratch files for C and Sigma in determinant basis, HQP CP, (HQQ-E) HQP CP
* 
      CALL FILEMAN_MINI(LU_Q1,'ASSIGN')
      CALL FILEMAN_MINI(LU_Q2,'ASSIGN')
      CALL FILEMAN_MINI(LU_Q3,'ASSIGN')
*
*. Calculate HQP CP in PQ space and save on LU_Q2
*
      ICSPC = IPSPC
      ISSPC = IQSPC
      CALL REWINO(LUCP)
      CALL REWINO(LU_Q2)
      CALL MV7(CB,HCB,LUCP,LU_Q2,XDUM,XDUM)
      IF(NTEST.GE.1000) WRITE(6,*) ' H(QP) CP accomplished'
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' H(QP) CP: '
       CALL WRTVCD(CB,LU_Q2,1,LBLK)
      END IF
*
*. Well really, we wanted HQP CP only in Q space, remove P space terms, 
*. save on LU_Q1
*
*. Array giving dimensions of occlass for given sym
      CALL MEMMAN(KLLOCCLS_SM,NOCCLS_MAX,'ADDL  ',1,'OC_LSM')
      IF(NTEST.GE.1000) WRITE(6,*) ' ISSM, NIRREP = ', ISSM, NIRREP
C          EXTRROW(WORK(KLCIBT),8,8,NBLOCK,LBLOCK)
      CALL EXTRROW(int_mb(KNCS_FOR_OCCLS),ISSM,NIRREP,NOCCLS_MAX,
     &             int_mb(KLLOCCLS_SM))
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Number of CSF of right sym per occls'
        CALL IWRTMA(int_mb(KLLOCCLS_SM),1,NOCCLS_MAX,1,NOCCLS_MAX)
      END IF
*
      CALL EXP_BLKVEC(LU_Q2,NOCCLS_P,int_mb(KLOCCLS_P),
     &                LU_Q1,NOCCLS_Q,int_mb(KLOCCLS_Q),
     &                int_mb(KLLOCCLS_SM),2,CB,XDUM,1,ICISTR,0)
*
C     EXP_BLKVEC(LU_IN,NBLK_IN, IBLK_IN,
C    &                      LU_OUT,NBLK_OUT,IBLK_OUT,
C    &                      LBLK,ITASK,VEC,VEC_OUT,IREW,ICISTR,
C    &                      INCORE)
      IF(NTEST.GE.1000) WRITE(6,*) ' HQP CP  calculated'
      IF(NTEST.GE.10000) THEN
       WRITE(6,*) ' H(QP) CP with P-terms removed: '
       CALL WRTVCD(CB,LU_Q1,1,LBLK)
      END IF
*
*. Calculate diagonal(HQQ) and save on LU_Q3
*
C          CSDIAGM(ISM,ICS,LUDIA_A)
C          CSDIAGM(ISM,ICS,LUDIA_A,VECUT,VECSCR)
      CALL CSDIAGM(IREFSM,2,LU_Q3,CB,HCB)
*
*. (HQQ - E)**-1 HQP CP on LU_Q2
*
      IPRECOND = 1
      LBLK = -1 
      IPNTR = 0
      H0 = 3006.56
      H0SCR = 3006.56
      IZERO = 0
      SHIFT = -E 
      IF(NTEST.GE.1000) WRITE(6,*) ' Shift :', SHIFT 
      CALL H0M1TD(LU_Q2,LU_Q3,LU_Q1,LBLK,IZERO,IPNTR,
     &     H0,SHIFT,H0SCR,XH0IX,
     &     IZERO,IZERO,IZERO,CB,HCB,NTEST,IPRECOND)
*     H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
      IF(NTEST.GE.1000) WRITE(6,*) ' Apres H0M1TD '
      CALL MEMCHK2('APH0M1')
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' (HQQ - E)**-1 HQP CP '
       CALL WRTVCD(CB,LU_Q2,1,LBLK)
      END IF
*. And the minus -1
      ONEM = -1.0D0
C SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
      CALL SCLVCD(LU_Q2,LUCQ,ONEM,CB,1,-1)
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The CQ vector '
        WRITE(6,*) ' ============= '
        CALL WRTVCD(CB,LUCQ,1,LBLK)
      END IF
*
      CALL FILEMAN_MINI(LU_Q1,'FREE')
      CALL FILEMAN_MINI(LU_Q2,'FREE')
      CALL FILEMAN_MINI(LU_Q3,'FREE')
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',1,'GET_CQ')
      CALL LUCIAQEXIT('GET_CQ')
*
      RETURN
      END
      SUBROUTINE GET_APPROX_DEN_AKBK_QQ(RHO1_QQ, RHO2_QQ,LUP,E,
     &                                 CB,HCB,I12)
*
* Obtain approximate densities for QQ interaction in AKBK approach
* for the Q vector corresponding of the P-vector given on LU_P
*
* The approximite densities are - those used to calculate energy, i.e.
*
*  RHO1_QQ(ij) = delta(ij) <0Q| E_ii |0Q>
*  RHO2_QQ(ijkl) = delta(ij)delta(kl) <0Q|E_ii E_kk |0Q>
*                 -delta(il)delta(kj) <0Q|a+is ais a+js ajs|0Q> (if I12 = 2)
*
*. Jeppe Olsen, Aug. 23, 2012 in Minneapolis
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'akbk.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cands.inc'
*. Output: Standard form of densities
      DIMENSION RHO1_QQ(*), RHO2_QQ(*)
*. Scratch for holding blocks of vectors
      DIMENSION CB(*), HCB(*)
*
      CALL MEMMAN(IDUM,IDUM,'MARK ',1,'AKQQDN')
      CALL LUCIAQENTER('AKQQCN')
* Some files
      CALL FILEMAN_MINI(LUQ,'ASSIGN')
      CALL FILEMAN_MINI(LUQ1,'ASSIGN')
*
*. 1: Generate the Q-vector from the given P-vector
*
C         GET_CQ_FROM_CP(LUCP,LUCQ,E,CB,HCB)
      CALL GET_CQ_FROM_CP(LUP,LUQ,E,CB,HCB)
*
*. Transform Q from CSF to SD basis and store on LUQ1
* 
*. requires some 
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
*. Block for storing complete or partial CI-vector
      IF(ISIMSYM.EQ.0) THEN
        LBLOCK = MXSOOB
      ELSE
        LBLOCK = MXSOOB_AS
      END IF
      IF(NOCSF.EQ.0) THEN
        LBLOCK  = MAX(NSD_FOR_OCCLS_MAX,MXSOOB)
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' TEST, LCSBLK, LBLOCK, MXSOOB  = ',
     &                     LCSBLK, LBLOCK, MXSOOB
      END IF
      ICOMP = 0
      
      CALL Z_BLKFO_FOR_CISPACE(IQSPC,ICSM,LBLOCK,ICOMP,
     &     NTEST,NQBLOCK,NQBATCH,
     &     int_mb(KCIOIO),int_mb(KCBLTP),NQOCCLS_ACT,
     &     int_mb(KCIOCCLS_ACT),
     &     int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &     int_mb(KCIBT),
     &     int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),0,ILTEST)

      CALL CSDTVCMN(CB,HCB,dbl_mb(KVEC3),
     &     1,0,ICSM,IQSPC,2,2,LUQ1,LUQ,NQOCCLS_ACT,
     &     int_mb(KCIOCCLS_ACT),int_mb(KCIBT),int_mb(KCLBT))
*
* 3: Generate diagonal contribitions to the density
*
      CALL DIAG_DET_TO_DENSI(CB,LUQ1,int_mb(KCBLTP),
     &     NQBLOCK,int_mb(KCIBT),RHO1_QQ, RHO2_QQ,I12)
C     DIAG_DET_TO_DENSI(CB,LUC,IBLTP,NBLOCK,IBLKFO,
C    &                             RHO1D,RHO2D,I12)
*
*. Clean up time
      CALL FILEMAN_MINI(LU_Q,'FREE')
      CALL FILEMAN_MINI(LU_Q1,'FREE')
      CALL MEMMAN(IDUM,IDUM,'FLUSM',1,'AKQQDN')
      CALL LUCIAQEXIT('AKQQCN')
*
      RETURN
      END
      SUBROUTINE DIAG_DET_TO_DENSI(CB,LUC,IBLTP,NBLOCK,IBLKFO,
     &                             RHO1D,RHO2D,I12)
*
*
* Diagonal contribution from determinants to one-and two-electron
* densities
* Jeppe Olsen, August 23, 2012
* Last revision, August 23 2012, Jeppe Olsen
*
c      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
*./ORBINP/: NACOB used
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc' 
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'glbbas.inc'
*
      DIMENSION IBLTP(*)
      DIMENSION IBLKFO(8,NBLOCK)
*
* ======
*.Output
* ======
      DIMENSION RHO1D(*), RHO2D(*)
*
      IDUM = 0
      CALL LUCIAQENTER('DIARHO')
      CALL MEMMAN(IDUM,  IDUM,    'MARK  ',IDUM,'DIARHO')
*
      NTEST = 1000
*
** Specifications of alpha - and beta types
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ==========================='
        WRITE(6,*) ' DIAG_DET_TO_DENSI speaking '
        WRITE(6,*) ' ==========================='
*
        WRITE(6,*) ' IATP IBTP NAEL NBEL = ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB = ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB = ', IOCTPA,IOCTPB
      END IF
*
**. Local memory
*
      IDUM = 0
*. Space for blocks of strings
      MAXA = IMNMX(int_mb(KNSTSO(IATP)),NSMST*NOCTPA,2)
      MAXB = IMNMX(int_mb(KNSTSO(IBTP)),NSMST*NOCTPB,2)
      CALL MEMMAN(KLASTR,MAXA*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MAXB*NBEL,'ADDL  ',1,'KLBSTR')
*
      CALL DIAG_DET_TO_DENSIS(NAEL,int_mb(KLASTR),NBEL,int_mb(KLBSTR),
     &       NACOB,CB,NSMST,int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)),
     &       LUC,IPRDIA,ICISTR,IBLTP,NBLOCK,IBLKFO,RHO1D,RHO2D,I12)

*.Flush local memory
      CALL MEMMAN(IDUM,  IDUM,    'FLUSM ',IDUM,'DIARHO')
      CALL LUCIAQEXIT('DIARHO')
*
      RETURN
      END
      SUBROUTINE DIAG_DET_TO_DENSIS(NAEL,IASTR,NBEL,IBSTR,
     &           NACOB,CB,NSMST,NSSOA,NSSOB,LUC,IPRNT,ICISTR,
     &           IBLTP,NBLOCK,IBLKFO,RHO1D,RHO2D,I12)
*
* Diagonal contributions to density matrix
*
* Jeppe Olsen, Aug. 23, Minneapolis
*
* Last revision, Aug 23 - 2012, Jeppe Olsen
*
      INCLUDE 'implicit.inc'
*.Input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
*
*. Scratch
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*), CB(*)
*. Output
      DIMENSION RHO1D(*), RHO2D(*) 
*
      NTEST = 1000
*
      IF( NTEST .GE. 100) THEN
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK               
          WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' IBLTP: '
        CALL IWRTMA(IBLTP,1,NSMST,1,NSMST)
        WRITE(6,*) ' I12 = ',I12
      END IF
*
      ITDET = 0
      LBLK = -1
      IF(LUC.NE.0) CALL REWINO(LUC)
*
      DO IBLK = 1, NBLOCK
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
        IF(NTEST.GE.1000) WRITE(6,*) ' IBLK,IATP,IBTP,IASM,IBSM=',
     &  IBLK,IATP,IBTP,IASM,IBSM
*. Read block of C in 
        IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(IDET,1,-1,LUC)
          NO_ZEROING = 0
          CALL FRMDSC2(CB,IDET,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' Number of elements readin ',IDET
          IDET = 0
        END IF
        IF(NTEST.GE.1000) WRITE(6,*) 'a: IBLK,IATP,IBTP,IASM,IBSM=',
     &  IBLK,IATP,IBTP,IASM,IBSM
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
        IF(NTEST.GE.1000) WRITE(6,*) 'b: IBLK,IATP,IBTP,IASM,IBSM=',
     &  IBLK,IATP,IBTP,IASM,IBSM
*. Read block of C in 
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF

        IBSTRT = 1                
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
          IBREL = IB - IBSTRT + 1
*
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
            IREST2 = 1
          ELSE
            IASTRT = 1                 
            IREST2 = 0
          END IF
          IASTOP = NSSOA(IASM,IATP) 
*
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            WRITE(6,*) ' IDET = ', IDET
            COEF = CB(IDET)
            WRITE(6,*) ' COEF = ', COEF
            IF(IREST2.EQ.0.OR.IA.EQ.IB) THEN
             WEIGHT = COEF**2
            ELSE
             WEIGHT = 2.0D0*COEF**2
            END IF
            IF(NTEST.GE.1000) WRITE(6,*) ' IA, IB, WEIGHT =',
     &      IA,IB, WEIGHT
*
*. Alpha contributions
*
            DO IAEL = 1, NAEL
              I = IASTR(IAEL,IA)
              II = (I-1)*NACOB + I
              RHO1D(II) = RHO1D(II) + WEIGHT
            END DO
*
*. beta contributions
*
            DO IBEL = 1, NBEL
              I = IBSTR(IBEL,IB)
              II = (I-1)*NACOB + I
              RHO1D(II) = RHO1D(II) + WEIGHT
            END DO
*
*alpha-alpha contributions
*
            IF(I12.EQ.2) THEN
             DO IAEL1 = 1, NAEL
             DO IAEL2 = 1, IAEL1
               I = IASTR(IAEL1,IA)
               K = IASTR(IAEL2,IA)
               II = (I-1)*NACOB + I
               KK = (K-1)*NACOB + I
               IIKK = II*(II-1)/2 + KK
               RHO2D(IIKK) =  RHO2D(IIKK) + WEIGHT
               IK = (K-1)*NACOB + I
               KI = (I-1)*NACOB + K
               KIIK = KI*(KI-1)/2 + IK
               RHO2D(KIIK) = RHO2D(KIIK) - WEIGHT
             END DO
             END DO
*
*beta-beta contributions
*
             DO IBEL1 = 1, NBEL
             DO IBEL2 = 1, IBEL1
               I = IBSTR(IBEL1,IB)
               K = IBSTR(IBEL2,IB)
               II = (I-1)*NACOB + I
               KK = (K-1)*NACOB + I
               IIKK = II*(II-1)/2 + KK
               RHO2D(IIKK) =  RHO2D(IIKK) + WEIGHT
               IK = (K-1)*NACOB + I
               KI = (I-1)*NACOB + K
               KIIK = KI*(KI-1)/2 + IK
               RHO2D(KIIK) = RHO2D(KIIK) - WEIGHT
             END DO
             END DO
*
* Alpha-beta contributions
*
             DO IAEL = 1, NAEL
             DO IBEL = 1, NBEL
               I = IASTR(IAEL,IA)
               K = IBSTR(IBEL,IB)
               II = (I-1)*NACOB + I
               KK = (K-1)*NACOB + I
               IF(II.GE.KK) THEN
                 IIKK = II*(II-1)/2 + KK
                 RHO2D(IIKK) =  RHO2D(IIKK) + WEIGHT
               END IF
               IK = (K-1)*NACOB + I
               KI = (I-1)*NACOB + K
               IF(KI.GE.IK) THEN
                 KIIK = KI*(KI-1)/2 + IK
                 RHO2D(KIIK) = RHO2D(KIIK) - WEIGHT
               END IF
             END DO
             END DO
            END IF ! I12 = 2
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over betastrings
      END DO
*     ^ End of loop over blocks
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Output from DIAG_DET_TO_DENSIS: '
        WRITE(6,*) ' --------------------------------' 
        WRITE(6,*)
        WRITE(6,*) ' Updated one-electron density '
        CALL WRTMAT(RHO1D,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*) ' Updated two-electron density '
        NAA = NACOB*(NACOB+1)/2
        CALL PRSYM(RHO2D,NAA)
      END IF
*
      RETURN
      END
      SUBROUTINE AKBK_DENSITY(RHO1,RHO2,L,R,LUP,E,I12)
*
* Obtain AKBK 1- and 2-body density(I12 = 2) matrices
*
*. Jeppe Olsen, Aug. 23 - 2012, Jeppe Olsen
* Last revision: Aug. 23 -2012, Jeppe Olsen
*
* The AKBK densities are given as 
*
* RHO1(ij) = <0^P!E_ij!0^P> + <0^Q!E_ij!0^P>  + <0^P!E_ij!0^Q> +
*               + Delta(ij) <0^Q!E_ii!0^Q> 
* 
* RHO2(ij,kl) = <0^P!E_ijE_kl - \delta(jk)E_il!0^P> 
*             + <0^P!E_ijE_kl - \delta(jk)E_il!0^Q> 
*             + <0^Q!E_ijE_kl - \delta(jk)E_il!0^P> 
*             + Delta(ij)delta(kl)  <0^Q!E_iiE_kk - \delta(jk)E_ik!0^Q>
*
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'akbk.inc'
      INCLUDE 'wrkspc-static.inc'
*. Output
      DIMENSION RHO1(*),RHO2(*)
*. Two scratch blocks
      REAL*8 L(*), R(*)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' AKBK_DENSITY at your service '
        WRITE(6,*) ' ============================='
      END IF
*
      LBLK = -1
      CALL LUCIAQENTER('AKBK_D')
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ', IDUMMY,'AKBK_D')
*. Local copies of Rho1 and Rho2
      LRHO1 = NACOB**2
      LRHO2 = (NACOB**2+1)*NACOB**2/2
      CALL MEMMAN(KLRHO1,LRHO1,'ADDL  ',2,'RHO1_L')
      CALL MEMMAN(KLRHO1B,LRHO1,'ADDL  ',2,'RHO1_B')
      IF(I12.EQ.2) THEN
        CALL MEMMAN(KLRHO2,LRHO2,'ADDL  ',2,'RHO2_L')
        CALL MEMMAN(KLRHO2B,LRHO2,'ADDL  ',2,'RHO2_L')
      END IF
*. And a few local files
      CALL FILEMAN_MINI(LUSC1,'ASSIGN')
      CALL FILEMAN_MINI(LUSC2,'ASSIGN')
      CALL FILEMAN_MINI(LUSC3,'ASSIGN')
*
*   Obtain P-P densities 
*  
C COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      CALL COPVCD(LUP,LUSC1,CB,1,LBLK)
C     DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR,EXPS2,IDOSRHO12,SRHO1,
C    &                  RHO2AA,RHO2AB,RHO2BB,IDONATORB)
      ICSPC = IPSPC
      ISSPC = IPSPC
      CALL DENSI2(I12,RHO1,RHO2,L,R,LUP,LUSC1,EXPS2,0,XDUM,
     &            XDUM,XDUM,XDUM,0)
*
*   Obtain P-Q densities 
*  
*. Q vector
      CALL GET_CQ_FROM_CP(LUP,LUSC1,E,L,R)
C     GET_CQ_FROM_CP(LUCP,LUCQ,E,CB,HCB)
      ISSPC = IPSPC
      ICSPC = IQSPC
      WRITE(6,*) ' DENSI2 will be called for PQ density '
      CALL DENSI2(I12,dbl_mb(KLRHO1),dbl_mb(KLRHO2),L,R,LUSC1,LUP,EXPS2,
     &     0,XDUM,XDUM,XDUM,XDUM,0)
      ONE = 1.0D0
      CALL VECSUM(RHO1,RHO1,dbl_mb(KLRHO1),ONE,ONE,LRHO1)
      IF(I12.EQ.2)
     &CALL VECSUM(RHO2,RHO2,dbl_mb(KLRHO2),ONE,ONE,LRHO2)
*
* Obtain Q-P densities - from  P-Q densities
*
      CALL TRPMAT(dbl_mb(KLRHO1),NACOB,NACOB,dbl_mb(KLRHO1B))
      CALL VECSUM(RHO1,RHO1,dbl_mb(KLRHO1B),ONE,ONE,LRHO1)
      IF(I12.EQ.2) THEN
        CALL CONJ_RHO2(NACOB,dbl_mb(KLRHO2),dbl_mb(KLRHO2B))
        CALL VECSUM(RHO2,RHO2,dbl_mb(KLRHO2B),ONE,ONE,LRHO2)
      END IF
*
* And then to the interesting: the Q-Q density
*
C     GET_APPROX_DEN_AKBK_QQ(RHO1_QQ, RHO2_QQ,LUP,E,CB,HCB,I12) 
      CALL GET_APPROX_DEN_AKBK_QQ(dbl_mb(KLRHO1),dbl_mb(KLRHO2),LUP,E,
     &    L,R,I12)
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' AKBK densities '
       WRITE(6,*) ' -------------- '
       WRITE(6,*) 
       WRITE(6,*) ' AKBK-rho1 '
       CALL WRTMAT(RHO1,NACOB,NACOB,NACOB,NACOB)
       IF(I12.EQ.2) THEN
         WRITE(6,*)
         WRITE(6,*) ' AKBK-rho2 '
         LAA = NACOB**2
         CALL PRSYM(RHO2,LAA)
       END IF
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ', IDUMMY,'AKBK_D')
      CALL LUCIAQEXIT('AKBK_D')
*
      RETURN
      END
      SUBROUTINE CONJ_RHO2(NACOB,RHO2_IN,RHO2_UT)
*
*. Conjugate two-electron density matrix RHO2_IN to obtain
*. RHO2_UT
* RHO2_UT(ji,lk) = RHO2_IN(ij,kl)
*
*. Jeppe Olsen, Aug. 23 2012, Minneapolis
*
      INCLUDE 'implicit.inc'  
*. Input
      DIMENSION RHO2_IN(*)
*. Output 
      DIMENSION  RHO2_UT(*)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from CONJ_RHO2'
        WRITE(6,*) ' ==================='
      END IF
*
      DO J = 1, NACOB
       DO I = 1, NACOB
        JI = (I-1)*NACOB + J
        IJ = (J-1)*NACOB + I
        DO L = 1, J
         IF(L.EQ.J) THEN
           KMAX = I
         ELSE
           KMAX = NACOB
         END IF
         DO K = 1, KMAX
           LK = (K-1)*NACOB + L
           KL = (L-1)*NACOB + K
           IF(IJ.GE.KL) THEN
             IJKL = IJ*(IJ-1)/2 + KL
           ELSE
             IJKL = KL*(KL-1)/2 + IJ
           END IF
           JILK = JI*(JI-1)/2 + KL
           RHO2_UT(JILK) = RHO2_IN(IJKL)
         END DO ! K
        END DO ! L
       END DO ! J
      END DO ! I
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' RHO2_IN: '
        LAA = NACOB**2
        CALL PRSYM(RHO2_IN,LAA)
        WRITE(6,*) ' RHO2_UT: '
        CALL PRSYM(RHO2_UT,LAA) 
      END IF
*
      RETURN
      END
      SUBROUTINE GET_FIA2(F2,FIA,IFIA2_OFFSS)
*
* Obtain First-order expansion matrix of FI + FA
*
* FIA(kappa)_mn = FIA(0)_mn + sum_(p>q) F2_mn,pq Kappa_pq, m geq n, p gt q
* m and n, as well as p and q, have identical symmetries
*
* F2_mn,pq  =  2 \delta(q,m) FIA_pn  + 2\delta(q,n) FIA_mp 
*           -  2 \delta(p,m) FIA_mq  - 2\delta(p,n) FIA_mq
*           +  (\delta(q,inactive) (4 (mn!pq) - (mp!qn)-(mq!pn))
*           -  (\delta(p,inactive) (4 (mn!qp) - (mq!pn)-(mp!qn))
*           +  (\delta(q,active) (Q(mnpq) + Q(nmpq))
*           -  (\delta(p,active) (Q(mnqp) + Q(nmqp))
*
* where Q(mnpq) =\sum_w (g(mnpw) - g(mqpn) D(qw) (q active)
*
*. Jeppe Olsen, Aug. 23, 2012, Minneapolis
*. Last revision: aug 24- 2012, Jeppe Olsen
*
*. Indeces corresponds to symmetry-type order, i.e. same order as integrals.
      INCLUDE 'implicit.inc'
      INTEGER P,Q,PM, PN,QM,QN, PQ, PT, QT
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'lucinp.inc'
*. Input
      DIMENSION FIA(*)
*. Output
      DIMENSION F2(*)
      INTEGER IFIA2_OFFSS(NSMOB,NSMOB)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' Info from  GET_FIA2 '
       WRITE(6,*) ' =================== '
       WRITE(6,*)
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'F2MAT ')
      CALL LUCIAQENTER('F2MAT')
*. Dimension of F2
      NR = 0
      NC = 0
      DO ISM = 1, NSMOB
        NR = NR + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
        NC = NC + NTOOBS(ISM)*(NTOOBS(ISM)-1)/2
      END DO
*. Offset to blocks of F with given symmetry of MN, PQ
      IMNPQ_OFF = 1
      DO IPQSM = 1, NSMOB
       DO IMNSM = 1, NSMOB
         NMN = NTOOBS(IMNSM)*(NTOOBS(IMNSM)+1)/2
         NPQ = NTOOBS(IPQSM)*(NTOOBS(IPQSM)-1)/2
         IFIA2_OFFSS(IMNSM,IPQSM) = IMNPQ_OFF
         IMNPQ_OFF = IMNPQ_OFF + NMN*NPQ
       END DO
      END DO
*
      ZERO = 0.0D0
      TWO = 2.0D0
      CALL SETVEC(F2,ZERO,NR*NC)
*
* =================
*. The delta terms
* =================
*
*. Contributes only to symmetry blocks with IMNSM = IPQSM - all four indeces have identical sym
      DO ISM = 1, NSMOB
       IF2_OFF =  IFIA2_OFFSS(ISM,ISM)
*. Offset to this symmetry block in FIA
       IF1_OFF = 1
       DO IISM = 1, ISM-1
        IF1_OFF = IF1_OFF + NTOOBS(IISM)*(NTOOBS(IISM)+1)/2
       END DO
*
       NO = NTOOBS(ISM)
       NOOR = NO*(NO+1)/2
       NOOC = NO*(NO-1)/2
*
*. Term1: 2 FIA_pn to F2(MN,PM)
*
       DO M = 1, NO
        DO N = 1, M
         DO P = M+1, NO
           MN = M*(M+1)/2 + N
           PM = (P-1)*(P-2)/2 + M
           MNPM = (PM-1)*NOOR + MN
           NP_F1 = (P-1)*(P-2)/2 + N
           F2(IF2_OFF-1+MNPM) = F2(IF2_OFF-1+MNPM) 
     &                        + TWO*FIA(IF1_OFF-1+NP_F1)
         END DO
        END DO
       END DO
*
* Term 2: 2FIA_mp to F2(MN,PN)
*
       DO M = 1, NO
        DO N = 1, M
         DO P = N+1, NO
          MN = M*(M-1)/2 + N
          PN = (P-1)*(P-2)/2 + N
          MNPN = (PN-1)*NOOR + MN
          MP_F1 = MAX(M,P)*(MAX(M,P)-1)/2 + MIN(M,P)
          F2(IF2_OFF-1+MNPN) =  F2(IF2_OFF-1+MNPN) 
     &                       + TWO*FIA(IF1_OFF-1+MP_F1)
         END DO
        END DO
       END DO
*
* Term 3: - 2FIA_qn to F2(MN,MQ)
*
       DO M = 1, NO
        DO N = 1, M
         DO Q = 1, M-1
          MN = M*(M-1)/2 + N
          MQ = (M-1)*(M-2)/2 + Q
          MNMQ = (MQ-1)*NOOR + MN
          NQ_F1 = MAX(N,Q)*(MAX(N,Q)-1)/2 + MIN(N,Q)
          F2(IF2_OFF-1+MNMQ) =  F2(IF2_OFF-1+MNMQ)  
     &                       - TWO*FIA(IF1_OFF+NQ_F1)
         END DO
        END DO
       END DO
*
* Term 4: -2FIA_mq to F2(MN,NQ)
*
       DO M = 1, NO  
        DO N = 1, NO
         DO Q = 1, N-1
           MN = M*(M-1)/2 + N
           NQ = (N-1)*(N-2)/2 + Q
           MNNQ = (NQ-1)*NOOR + MN
           MQ_F1 = M*(M-1)/2 + Q
           F2(IF2_OFF-1+MNNQ) =  F2(IF2_OFF-1+MNNQ)  
     &                        - TWO*FIA(IF1_OFF+MQ_F1)
         END DO
        END DO
       END DO
*
      END DO ! End of symmetry blocks of F2
*
* ==================================================================
* The terms
*           +  (\delta(q,inactive) (4 (mn!pq) - (mp!qn)-(mq!pn))
*           -  (\delta(p,inactive) (4 (mn!qp) - (mq!pn)-(mp!qn))
* ==================================================================
*
      DO IMNSM = 1, NSMOB
       DO IPQSM = 1, NSMOB
        NMN = NTOOBS(IMNSM)
        NPQ = NTOOBS(IPQSM)
        NOOR = NMN*(NMN+1)/2
        IF2_OFF =  IFIA2_OFFSS(IMNSM,IPQSM)
        DO M = 1, NMN
         DO N = 1, M
          DO P = 1, NPQ
           DO Q = 1, P-1
             IF(ITPFSO(Q).EQ.0) THEN
*. Q is Inactive 
               MN = M*(M-1)/2 + N
               PQ = (P-1)*(P-2)/2 + Q
               MNPQ = (PQ-1)*NOOR + MN
               MT = IREOST(M)
               NT = IREOST(N)
               PT = IREOST(P)
               QT = IREOST(Q)
C  4 (mn!pq) - (mp!qn)-(mq!pn))
               X = 4.0D0* GTIJKL_GN(MT,NT,PT,QT)-GTIJKL_GN(MT,PT,QT,NT)
     &           - GTIJKL_GN(MT,QT,PQ,NT)
               F2(MNPQ) = F2(MNPQ) + X
             END IF
             IF(ITPFSO(P).EQ.0) THEN
*. P is inactive
               MN = M*(M-1)/2 + N
               PQ = (P-1)*(P-2)/2 + Q
               MNPQ = (PQ-1)*NOOR + MN
               MT = IREOST(M)
               NT = IREOST(N)
               PT = IREOST(P)
               QT = IREOST(Q)
C   -(4 (mn!qp) - (mq!pn)-(mp!qn))
               X = 4.0D0*GTIJKL_GN(MT,NT,QT,PT) 
     &           - GTIJKL_GN(MT,QT,PT,NT)
     &           - GTIJKL_GN(MT,PT,QT,NT)
                F2(MNPQ) = F2(MNPQ) - X
             END IF
*
            END DO
           END DO
          END DO
         END DO ! loop over M,N,P,Q
        END DO
       END DO ! loop over IMNSM, IPQSM
*. The terms 
*     +  (\delta(q,active) (Q(mnpq) + Q(nmpq))
*     -  (\delta(p,active) (Q(mnqp) + Q(nmqp))
* with  Q(mnpq) =\sum_w (g(mnpw) - g(mqpn) D(qw) (q active)
*
*
      LQ = 0
      DO IMNSM = 1, NSMOB
       DO IPQSM = 1, NSMOB
        LQ = LQ + NTOOBS(IMNSM)**2*NTOOBS(IPQSM)*NACOBS(IPQSM)
       END DO
      END DO
      CALL MEMMAN(KLQ,LQ,'ADDL  ',2,'QMAT  ')
      CALL MEMMAN(KLQP,NSMOB**2,'ADDL  ',1,'IQ_OFF')
      CALL GET_QMAT(WORK(KLQ),WORK(KLQP))
      IOFFQ = 1
      IOFFF = 1
      DO IMNSM = 1, NSMOB
        DO IPQSM = 1, NSMOB
         NMN = NTOOBS(IMNSM)
         NOOR = NMN*(NMN+1)/2
         DO M = 1, NMN
          DO N = 1, NMN
           MN = M*(M-1)/2 + N
* +  (\delta(q,active) (Q(mnpq) + Q(nmpq))
           NQINA = NINOBS(IPQSM)
           DO P = 1, NTOOBS(IPQSM)
            DO Q = NQINA+1, MIN(NQINA + NACOBS(IMNSM),P)
              PQ = (P-1)*(P-2)/2 + Q
              MNPQ_F = (PQ-1)*NOOR + MN
              MNPQ_Q = (Q-NQINA-1)*NMN**2*NTOOBS(IPQSM)
     &               + (P-1)      *NMN**2
     &               + (N-1)      *NMN
     &               + M
              F2(IOFFF-1+MNPQ_F) = 
     &        F2(IOFFF-1+MNPQ_F)  + WORK(KLQ-1+IOFFQ-1+MNPQ_Q)
              NMPQ_Q = (Q-NQINA-1)*NMN**2*NTOOBS(IPQSM)
     &               + (P-1)      *NMN**2
     &               + (M-1)      *NMN
     &               + N
              F2(IOFFF-1+MNPQ_F) = 
     &        F2(IOFFF-1+MNPQ_F) + WORK(KLQ-1+IOFFQ-1+NMPQ_Q)
            END DO
           END DO ! Loops over pq
*     -  (\delta(p,active) (Q(mnqp) + Q(nmqp))
           DO P = NINOBS(IPQSM)+1, NINOBS(IPQSM)+NACOBS(IPQSM)
            NPINA = NINOBS(IPQSM)
            DO Q = 1, P
              PQ = (P-1)*(P-2)/2 + Q
              MNPQ_F = (PQ-1)*NOOR + MN
              MNQP_Q = (P-NQINA-1)*NMN**2*NTOOBS(IPQSM)
     &               + (Q-1)      *NMN**2
     &               + (N-1)      *NMN
     &               + M
              F2(IOFFF-1+MNPQ_F) = 
     &        F2(IOFFF-1+MNPQ_F)  + WORK(KLQ-1+IOFFQ-1+MNQP_Q)
              NMQP_Q = (P-NQINA-1)*NMN**2*NTOOBS(IPQSM)
     &               + (Q-1)      *NMN**2
     &               + (M-1)      *NMN
     &               + N
              F2(IOFFF-1+MNPQ_F) = 
     &        F2(IOFFF-1+MNPQ_F) + WORK(KLQ-1+IOFFQ-1+NMQP_Q)
            END DO
           END DO ! Loops over pq
*
          END DO
         END DO ! loops over MN
         IOFFQ = IOFFQ + NTOOBS(IMNSM)**2*NTOOBS(IPQSM)*NACOBS(IPQSM)
         IOFFF = IOFFF  
     &         + NTOOBS(IMNSM)*(NTOOBS(IMNSM)+1)/2
     &         * (NTOOBS(IPQSM)-1)*(NTOOBS(IPQSM)-2)/2
        END DO
      END DO ! Loops over IMNSM, IPQSM
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The F2-matrix as delivered '
        WRITE(6,*) ' ========================== '
        IOFF = 1
        DO IMNSM = 1, NSMOB 
         DO IPQSM = 1, NSMOB
          NMN = NTOOBS(IMNSM)*(NTOOBS(IMNSM)-1)/2
          NPQ = (NTOOBS(IPQSM)-1)*(NTOOBS(IPQSM)-2)/2
          WRITE(6,'(A,2I3)') ' Block with MNSM, PQSM = ',
     &    IMNSM,IPQSM
          CALL WRTMAT(F2,NMN,NPQ,NMN,NPQ)
          IOFF = IOFF + NMN*NPQ
         END DO
        END DO
      END IF
*

      CALL LUCIAQEXIT('F2MAT')
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'F2MAT ')
*
      RETURN
      END
      SUBROUTINE GET_QMAT(QM,IQ_OFF)
*
* Obtain Q-matrix: Q(m,n,p,q) = \sum_w ((mn!pq) - (mw!pn)) D(qw)
* First two indeces are required to have identical symmetry- as are the last two
* Note that the last index is active
*
* Jeppe Olsen, Aug. 23, 2012, Minneapolis
* Last revision, Aug. 23 - 2012, Jeppe Olsen
*
*. Initial version - celebrating the extreme power of current computers
*                    (positive spin on the expression: a very sloppy and inefficient initial code)
*
*
      INCLUDE 'implicit.inc'
      INTEGER P,Q,W,PT,QT,WT,PQ
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'glbbas.inc'
*. Output
      DIMENSION QM(*)
      INTEGER IQ_OFF(NSMOB,NSMOB)
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'QMAT ')
      NAA = NACOB**2
      CALL MEMMAN(KLRHO1,NAA,'ADDL  ',2,'RHO1_S')
      CALL MEMMAN(KLVEC,NACOB,'ADDL  ',2,'VEC_AC')
      CALL MEMMAN(KLVEC2,NACOB,'ADDL  ',2,'VEC_A2')
*
* =======
* Offsets
* =======
*
      IOFF = 1
      DO IMNSM = 1, NSMOB
       DO IPQSM = 1, NSMOB
         IQ_OFF(IMNSM,IPQSM) = IOFF
         NELMNT = 
     &   NTOOBS(IMNSM)*NTOOBS(IMNSM)*NTOOBS(IPQSM)*NACOBS(IPQSM)
         IOFF = IOFF + NELMNT
       END DO
      END DO
*
*======================================
* Reform Rho1 to Symmetry-blocked form
*======================================
*

C     REFORM_RHO1_TO_GNSM(RHO1_ST,RHO1_GNSM_ST,IWAY,IREO_GNSYM_TO_TS)
      CALL REFORM_RHO1_TO_GNSM(WORK(KRHO1),WORK(KLRHO1),1,
     &                         WORK(KIREO_GNSYM_TO_TS_ACOB))
*
* ===================
* And then the matrix 
* ===================
*
      DO IMNSM = 1, NSMOB
       DO IPQSM = 1, NSMOB
        IOFF = IQ_OFF(IMNSM,IPQSM)
        NTMN = NTOOBS(IMNSM)
        NTPQ = NTOOBS(IPQSM)
        NAPQ = NACOBS(IPQSM)
        DO M = 1, NTMN
        DO N = 1, NTMN
        DO P = 1, NTPQ
         MT = IREOST(M)
         NT = IREOST(N)
         PT = IREOST(P)
*. Collect integrals vec(w) = (mn!pw) - (mw!pn)
         DO W = 1, NAPQ
          WT = IREOST(W)
          WORK(KLVEC-1+W) = 
     &    GTIJKL_GN(MT,NT,PT,WT)-GTIJKL_GN(MT,QT,PT,NT)
         END DO
         IOFF_R1 = 1
         DO ISM = 1, IPQSM-1
           IOFF_R1 = IOFF_R1 + NACOBS(ISM)**2
         END DO
         CALL MATVCB(WORK(KLRHO1-1+IOFF_R1),WORK(KLVEC),WORK(KLVEC2),
     &        NAPQ,NAPQ,0)
         MN = (N-1)*NTMN + M
         DO Q = 1, NAPQ
           PQ = (Q-1)*NTPQ + P
           MNPQ = (PQ-1)*NTMN*NTMN+PQ
           QM(IOFF-1+MNPQ) = WORK(KLVEC2-1+Q)
         END DO! End of loop over Q
        END DO
        END DO
        END DO ! End of loops over M,N,P
       END DO ! Loop over IPQSM
      END DO ! Loop over IMNSM
*
      NTEST = 1000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============='
        WRITE(6,*) ' The Q-matrix '
        WRITE(6,*) ' ============='
        WRITE(6,*) 
        DO IMNSM = 1, NSMOB
          DO IPQSM = 1, NSMOB
            WRITE(6,'(A,2I3)') 
     &      ' Block with IMNSM, IPQSM = ', IMNSM,IPQSM
            IOFF = IQ_OFF(IMNSM,IPQSM)
            NMN = NTOOBS(IMNSM)**2
            NPQ = NTOOBS(IPQSM)*NACOBS(IPQSM)
            CALL  WRTMAT(QM(IOFF),NMN,NPQ,NMN,NPQ)
          END DO
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'QMAT ')
*
      RETURN
      END
      SUBROUTINE GET_CE_BLOCKS_OF_F2(NCEXC,ICEXC,NEEXC,IEEXC, 
     &           F2,F2CC,F2CE,F2EC,F2EE,IF2_OFF)
*
* Obtain C and E blocks of F from the complete F-matrix
*
*. Jeppe Olsen, Aug. 23 - 2012, Minneapolis
*  Last revision, Aug. 23 2012, Jeppe Olsen
*
      INCLUDE 'implicit.inc'
      INTEGER P,Q,PS,QS, PS_REL, QS_REL, PQ
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*. Input
      DIMENSION F2(*), IF2_OFF(NSMOB,NSMOB)
      INTEGER ICEXC(2,NCEXC),IEEXC(2,NEEXC)
*. Output
      DIMENSION F2CC(NCEXC*NCEXC),F2CE(NCEXC*NEEXC) 
      DIMENSION F2EE(NEEXC*NEEXC),F2EC(NEEXC*NCEXC)
*
      NTEST = 1000
*
*=====
*F2CC
*=====
*
      DO IMN_EXC = 1, NCEXC
      DO IPQ_EXC = 1, NCEXC
*
        M = ICEXC(1,IMN_EXC)
        N = ICEXC(2,IMN_EXC)
        P = ICEXC(1,IPQ_EXC)
        Q = ICEXC(2,IPQ_EXC)
*
        MS = IREOTS(M)
        NS = IREOTS(N)
        PS = IREOTS(P)
        QS = IREOTS(Q)
*
        IMNSM = ISMFSO(MS)
        IPQSM = ISMFSO(PS)
        NMN = NTOOBS(IMNSM)
        NMNMN = NMN*(NMN+1)/2
*
        MS_REL = MS - IBSO(IMNSM)+1
        NS_REL = NS - IBSO(IMNSM)+1
        PS_REL = PS - IBSO(IPQSM)+1
        QS_REL = QS - IBSO(IPQSM)+1
*
        MN = MS_REL*(MS_REL-1)/2 + NS_REL
        PQ = (PS_REL-1)*(PS_REL-2)/2 + QS_REL
        MNPQ = (PQ-1)*NMNMN + MN
*
        MNPQ_UT = (IMN_EXC-1)*NCEXC + IPQ_EXC
        F2CC(MNPQ_UT) = F2(IF2_OFF(IMNSM,IPQSM)-1+MNPQ)
      END DO
      END DO
*=====
*F2CE
*=====
*
      DO IMN_EXC = 1, NCEXC
      DO IPQ_EXC = 1, NEEXC
*
        M = ICEXC(1,IMN_EXC)
        N = ICEXC(2,IMN_EXC)
        P = IEEXC(1,IPQ_EXC)
        Q = IEEXC(2,IPQ_EXC)
*
        MS = IREOTS(M)
        NS = IREOTS(N)
        PS = IREOTS(P)
        QS = IREOTS(Q)
*
        IMNSM = ISMFSO(MS)
        IPQSM = ISMFSO(PS)
        NMN = NTOOBS(IMNSM)
        NMNMN = NMN*(NMN+1)/2
*
        MS_REL = MS - IBSO(IMNSM)+1
        NS_REL = NS - IBSO(IMNSM)+1
        PS_REL = PS - IBSO(IPQSM)+1
        QS_REL = QS - IBSO(IPQSM)+1
*
        MN = MS_REL*(MS_REL-1)/2 + NS_REL
        PQ = (PS_REL-1)*(PS_REL-2)/2 + QS_REL
        MNPQ = (PQ-1)*NMNMN + PQ
*
        MNPQ_UT = (IMN_EXC-1)*NCEXC + IPQ_EXC
        F2CE(MNPQ_UT) = F2(IF2_OFF(IMNSM,IPQSM)-1+MNPQ)
      END DO
      END DO
*
*=====
*F2EC
*=====
*
      DO IMN_EXC = 1, NEEXC
      DO IPQ_EXC = 1, NCEXC
*
        M = IEEXC(1,IMN_EXC)
        N = IEEXC(2,IMN_EXC)
        P = ICEXC(1,IPQ_EXC)
        Q = ICEXC(2,IPQ_EXC)
*
        MS = IREOTS(M)
        NS = IREOTS(N)
        PS = IREOTS(P)
        QS = IREOTS(Q)
*
        IMNSM = ISMFSO(MS)
        IPQSM = ISMFSO(PS)
        NMN = NTOOBS(IMNSM)
        NMNMN = NMN*(NMN+1)/2
*
        MS_REL = MS - IBSO(IMNSM)+1
        NS_REL = NS - IBSO(IMNSM)+1
        PS_REL = PS - IBSO(IPQSM)+1
        QS_REL = QS - IBSO(IPQSM)+1
*
        MN = MS_REL*(MS_REL-1)/2 + NS_REL
        PQ = (PS_REL-1)*(PS_REL-2)/2 + QS_REL
        MNPQ = (PQ-1)*NMNMN + PQ
*
        MNPQ_UT = (IMN_EXC-1)*NEEXC + IPQ_EXC
        F2EC(MNPQ_UT) = F2(IF2_OFF(IMNSM,IPQSM)-1+MNPQ)
      END DO
      END DO
*
*=====
*F2EE
*=====
*
      DO IMN_EXC = 1, NEEXC
      DO IPQ_EXC = 1, NEEXC
*
        M = IEEXC(1,IMN_EXC)
        N = IEEXC(2,IMN_EXC)
        P = IEEXC(1,IPQ_EXC)
        Q = IEEXC(2,IPQ_EXC)
*
        MS = IREOTS(M)
        NS = IREOTS(N)
        PS = IREOTS(P)
        QS = IREOTS(Q)
*
        IMNSM = ISMFSO(MS)
        IPQSM = ISMFSO(PS)
        NMN = NTOOBS(IMNSM)
        NMNMN = NMN*(NMN+1)/2
*
        MS_REL = MS - IBSO(IMNSM)+1
        NS_REL = NS - IBSO(IMNSM)+1
        PS_REL = PS - IBSO(IPQSM)+1
        QS_REL = QS - IBSO(IPQSM)+1
*
        MN = MS_REL*(MS_REL-1)/2 + NS_REL
        PQ = (PS_REL-1)*(PS_REL-2)/2 + QS_REL
        MNPQ = (PQ-1)*NMNMN + PQ
*
        MNPQ_UT = (IMN_EXC-1)*NEEXC + IPQ_EXC
        F2EE(MNPQ_UT) = F2(IF2_OFF(IMNSM,IPQSM)-1+MNPQ)
      END DO
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from GET_CE_BLOCKS_OF_F2 '
        WRITE(6,*) ' ================================'
        WRITE(6,*) 
        WRITE(6,*) ' F2CC: '
        CALL WRTMAT(F2CC,NCEXC,NCEXC,NCEXC,NCEXC)
        WRITE(6,*) ' F2CE: '
        CALL WRTMAT(F2CC,NCEXC,NEEXC,NCEXC,NEEXC)
        WRITE(6,*) ' F2EC '
        CALL WRTMAT(F2CC,NEEXC,NCEXC,NEEXC,NCEXC)
        WRITE(6,*) ' F2EE '
        CALL WRTMAT(F2CC,NEEXC,NEEXC,NEEXC,NEEXC)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_E1_E2_EFF(E1C,E1E,E2CC,E2EE,E2CE,E2EC,NC,NE,
     &           F2CC,F2CE,E1EFF, E2EFF)
*
* Generate the effective gradient and Hessian for the AKBK optimization
*
*. Jeppe Olsen, Aug. 24 - 2012, Minneapolis
*. Last revision, Aug. 24 2012, Jeppe Olsen
*
*. E1EFF = E1E - F2CE T * F2CC (T) ** -1 E1C
*. E2EFF = E2EE - F2CC (T) ** -1 E2CE
*.              - E2EC*F2CC ** -1 F2CE
*.              + F2CC (T) ** -1 E2CC F2CC ** -1 F2CE
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
*. Input
      DIMENSION E1C(NC),E1E(NE)
      DIMENSION E2CC(NC,NC),E2CE(NC,NE),E2EC(NE,NC),E2EE(NE,NE)
*
      NTEST = 1000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_E1_E2_EFF speaking '
        WRITE(6,*) ' ====================== '
      END IF
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT_EEF')
*
*. 1 Obtain F2CC ** -1 
*
      NMAX = MAX(NC,NE)
      CALL MEMMAN(KLF2CCI,NMAX**2,'ADDL  ',2,'F2CCI ')
      CALL MEMMAN(KLF2CP,NMAX**2,'ADDL  ',2,'F2CCCP')
      CALL MEMMAN(KLF2CPB,NMAX**2,'ADDL  ',2,'F2CCCB')
      CALL MEMMAN(KLCVEC,NMAX,'ADDL  ',2,'CVEC  ')
*
C     INVMAT(A,B,MATDIM,NDIM,ISING)
      CALL COPVEC(F2CC, WORK(KLF2CP),NC**2)
      CALL INVMAT(WORK(KLF2CP),WORK(KLF2CCI),NC,NC,ISING)
      IF(ISING.NE.0) THEN
        WRITE(6,*) ' Warning: problems with inversion of E2CC '
        WRITE(6,*) ' However, I am in a optimistic mood and continue'
      END IF
*
* 2: E1EFF = E1E - F2CE T * F2CC (T) ** -1 E1C
*
C MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
      CALL MATVCC(WORK(KLF2CCI),E1C,WORK(KLCVEC),NC,NC,1)
      CALL MATVCC(F2CE,WORK(KLCVEC),E1EFF,NC,NE,1)
      ONEM = -1.0D0
      ONE  = 1.0D0
      CALL VECSUM(E1EFF,E1EFF,E1E,ONEM,ONE,NC)
*
* 3 : E2EFF = E2EE -F2CE T F2CC (T) ** -1 E2CE
*.          - E2EC*F2CC ** -1 F2CE
*.          + F2CE(T) F2CC (T) ** -1 E2CC F2CC ** -1 F2CE
*
* a) E2EE - F2CE T F2CC (T) ** -1 E2CE
*
      CALL COPVEC(E2EE,E2EFF,NE**2)
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(WORK(KLF2CP), WORK(KLF2CCI),E2CE,
     &     NC,NE,NC,NC,NC,NE,FACTORC,FACTORAB,1)
      FACTORC = 1.0D0
      FACTORAB = -1.0D0
      CALL MATML7(E2EFF,F2CE,WORK(KLF2CP),
     &            NE,NE,NC,NE,NC,NE,
     &            FACTORC,FACTORAB,1)
*
* b) - E2EC*F2CC ** -1 F2CE
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(WORK(KLF2CP),E2EC, WORK(KLF2CCI),
     &     NC,NE,NC,NC,NC,NE,FACTORC,FACTORAB,0)
      FACTORC = 1.0D0
      FACTORAB = -1.0D0
      CALL MATML7(E2EFF,WORK(KLF2CP),F2CE,
     &     NE,NE,NE,NC,NC,NE,FACTORC,FACTORAB,0)
*
* c)  + F2CE(T) F2CC (T) ** -1 E2CC F2CC ** -1 F2CE
*
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(WORK(KLF2CP),WORK(KLF2CCI),F2CE,
     &     NC,NE,NC,NC,NC,NE,FACTORC,FACTORAB,0)
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(WORK(KLF2CPB),E2CC,WORK(KLF2CP),
     &     NC,NE,NC,NC,NC,NE,FACTORC,FACTORAB,0)
      FACTORC = 0.0D0
      FACTORAB = 1.0D0
      CALL MATML7(WORK(KLF2CP),WORK(KLF2CCI),WORK(KLF2CPB),
     &     NC,NE,NC,NC,NC,NE,FACTORC,FACTORAB,1)
      FACTORC = 1.0D0
      FACTORAB = 1.0D0
      CALL MATML7(E2EFF,F2CE,WORK(KLF2CP),
     &     NE,NE,NC,NE,NC,NE,FACTORC,FACTORAB,1)
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' Effective gradient'
       CALL WRTMAT(E1EFF,NE,1,NE,1)
       WRITE(6,*)
       WRITE(6,*) ' Effective Hessian '
       CALL WRTMAT(E2EFF,NE,NE,NE,NE)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GT_EEF')
      RETURN
      END
      SUBROUTINE LUCIA_AKMCSCF(IREFSM,IREFSPC_MCSCF,MAXMAC,MAXMIC,
     &                       EFINAL,CONVER,VNFINAL)
*
* Master routine for MCSCF optimization within the AK approximation
*
* Jeppe Olsen, Aug 24, 2012, Minneapolis
* Last revision, Aug. 23 2012, Jeppe Olsen

* Initial MO-INI transformation matrix is assumed set outside and is in MOMO
* Initial MO-AO transformation matrix is in MOAOIN
*
*. Output matrix is in
*   1) MOAOUT   - as it is the output matrix
*   2) MOAO_ACT - as it is the active matrix
*   3) MOAOIN   - as the integrals are in this basis ...
      INCLUDE 'wrkspc.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'intform.inc'
      INCLUDE 'cc_exc.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'akbk.inc'
      INCLUDE 'spinfo.inc'
*. Some indirect transfer
      COMMON/EXCTRNS/KLOOEXCC,KINT1_INI,KINT2_INI, IREFSPC_MCSCFL,
     &               IPRDIALL,IIUSEH0PL,MPORENP_EL,
     &               ERROR_NORM_FINALL,CONV_FL,
     &               I_DO_CI_IN_INNER_ACT
cVOG   COMMON/EXCTRNS/KLOOEXCC,KINT1_INI,KINT2_INI
* A bit of local scratch
      INTEGER ISCR(2), ISCR_NTS((7+MXPR4T)*MXPOBS)
*
      REAL*8
     &INPROD
      LOGICAL DISCH, CONV_INNER
*
      LOGICAL CONV_F,CONVER
      EXTERNAL EMCSCF_FROM_KAPPA
* Removing (incorrect) compiler warnings
      KINT2_FSAVE = 0
      IE2ARR_F = -1

      IDUMMY = 0
      CALL MEMMAN(IDUMMY, IDUMMY, 'MARK ', IDUMMY,'AKMCSC') 
      CALL LUCIAQENTER('AKMCS')
*
*. Local parameters defining optimization
*
*. reset kappa to zero in each inner or outer iteration
*
* IRESET_KAPPA_IN_OR_OUT = 1 => Reset kappa in each inner iteration
* IRESET_KAPPA_IN_OR_OUT = 2 => Reset kappa in each outer iteration
*
*. Use gradient or Brillouin vector (differs only when gradient is 
*  evaluated for Kappa ne. 0, ie. IRESET_KAPPA = 2
*
* I_USE_BR_OR_E1 = 1 => Use Brilloin vector
* I_USE_BR_OR_E2 = 2 => Use E1
      IRESET_KAPPA_IN_OR_OUT = 2
      I_USE_BR_OR_E1 = 2 
*. Largest allowed number of vectors in update
      NMAX_VEC_UPDATE = 50
*. Restrict orbital excitations in case of super-symmetry
      INCLUDE_ONLY_TOTSYM_SUPSYM = 1
*
      WRITE(6,*) ' *****************************************'
      WRITE(6,*) ' *                                       *'
      WRITE(6,*) ' * AK MCSCF optimization control entered *'
      WRITE(6,*) ' *                                       *'
      WRITE(6,*) ' * Version -1, Jeppe Olsen, August 12    *'
      WRITE(6,*) ' *****************************************'
      WRITE(6,*)
      WRITE(6,*) ' Occupation space: ', IREFSPC_MCSCF
      WRITE(6,*) ' Allowed number of outer iterations ', MAXMAC
      WRITE(6,*) ' Allowed number of inner iterations ', MAXMIC
*
* Some temporary(?) parameters to define transition from standard 
* to AK optimization
*
      I_DO_APPROX_DENS = 1
      I_DO_APPROX_CI = 1
      I_DO_CONSTR = 1
      WRITE(6,*) ' Temporary parameters: '
      WRITE(6,'(A,3I3)') 
     &'I_DO_APPROX_DENS, I_DO_APPROX_CI,I_DO_CONSTR ',
     & I_DO_APPROX_DENS, I_DO_APPROX_CI,I_DO_CONSTR   
*
      IF(I_USE_SUPSYM.EQ.1) THEN
        IF(INCLUDE_ONLY_TOTSYM_SUPSYM.EQ.1) THEN
          WRITE(6,*) 
     &   ' Excitations only between orbs with the same supersymmetry'
        ELSE
          WRITE(6,'(2X,A)') 
     &   'Excitations only between orbs with the same standard symmetry'
        END IF
      END IF
*
      WRITE(6,*)
      WRITE(6,*) ' MCSCF optimization method in action:'
      IF(IMCSCF_MET.EQ.1) THEN
        WRITE(6,*) '    One-step method NEWTON'
      ELSE  IF (IMCSCF_MET.EQ.2) THEN
        WRITE(6,*) '    Two-step method NEWTON'
      ELSE  IF (IMCSCF_MET.EQ.3) THEN
        WRITE(6,*) '    One-step method Update'
      ELSE  IF (IMCSCF_MET.EQ.4) THEN
        WRITE(6,*) '    Two-step method Update'
      END IF
*
      IF(IOOE2_APR.EQ.1) THEN
        WRITE(6,*) '    Orbital-Orbital Hessian constructed'
      ELSE IF(IOOE2_APR.EQ.2) THEN
        WRITE(6,*) 
     &  '    Diagonal blocks of Orbital-Orbital Hessian constructed'
      ELSE IF(IOOE2_APR.EQ.3) THEN
        WRITE(6,*) 
     &  '    Approx. diagonal of Orbital-Orbital Hessian constructed'
      END IF
*
*. Linesearch
*
      IF(IMCSCF_MET.LE.2) THEN
       IF(I_DO_LINSEA_MCSCF.EQ.1) THEN 
         WRITE(6,*) 
     &   '    Line search for Orbital optimization '
       ELSE IF(I_DO_LINSEA_MCSCF.EQ.0) THEN
         WRITE(6,*) 
     &   '    Line search when energy increases '
       ELSE IF(I_DO_LINSEA_MCSCF.EQ.-1) THEN
         WRITE(6,*) 
     &   '    Line search never carried out '
       END IF
      ELSE
*. Update method linesearch always used
        WRITE(6,*) 
     &  '    Line search for Orbital optimization '
      END IF
      IF(IMCSCF_MET.EQ.3.OR.IMCSCF_MET.EQ.4) THEN
        WRITE(6,'(A,I4)') 
     &  '     Max number of vectors in update space ', NMAX_VEC_UPDATE
      END IF
*
      IF(IRESET_KAPPA_IN_OR_OUT .EQ.1 ) THEN
        WRITE(6,*) 
     &  '       Kappa is reset to zero in each inner iteration '
      ELSE IF( IRESET_KAPPA_IN_OR_OUT .EQ.2 ) THEN
        WRITE(6,*) ' Kappa is reset to zero in each outer iteration '
      END IF
*
      IF(I_USE_BR_OR_E1.EQ.1) THEN
        WRITE(6,*) '    Brillouin vector in use'
      ELSE IF(I_USE_BR_OR_E1 .EQ.2) THEN
        WRITE(6,*) '    Gradient vector in use'
      END IF
*
      IF(NFRZ_ORB.NE.0) THEN
        WRITE(6,*) ' Orbitals frozen in MCSCF optimization: '
        CALL IWRTMA3(IFRZ_ORB,1,NFRZ_ORB,1,NFRZ_ORB)
      END IF
      
      I_MAY_DO_CI_IN_INNER_ITS = 0
      XKAPPA_THRES = 0.0D0
      MIN_OUT_IT_WITH_CI = 10
      IF(IMCSCF_MET.EQ.3) THEN
        I_MAY_DO_CI_IN_INNER_ITS = 1
        XKAPPA_THRES = 1.0D-1
        WRITE(6,'(A)') 
     &  '     CI - optimization in inner iterations starts when: '
        WRITE(6,'(A)')
     &  '       Hessian approximation is not shifted'
        WRITE(6,'(A,E8.2)') 
     &  '       Initial step is below ',  XKAPPA_THRES
        WRITE(6,'(A,I3)') 
     &  '     Outer iteration is atleast number ', MIN_OUT_IT_WITH_CI
      END IF
*
*. Initial allowed step length 
      STEP_MAX = 0.75D0
C     WRITE(6,*) ' Jeppe has reduced step to ', STEP_MAX
      TOLER = 1.1D0
      NTEST = 1000
      IPRNT= MAX(NTEST,IPRMCSCF)
*
      I_DO_NEWTON = 0
      I_DO_UPDATE = 0
      I_UPDATE_MET = 0
      IF(IMCSCF_MET.LE.2) THEN
        I_DO_NEWTON = 1
      ELSE IF (IMCSCF_MET.EQ.3.OR.IMCSCF_MET.EQ.4) THEN
        I_DO_UPDATE = 1
*. use BFGS update
        I_UPDATE_MET = 2
*. Update vectors will be kept in core
        DISCH = .FALSE.
        LUHFIL = -2810
      END IF
      WRITE(6,*) ' I_DO_NEWTON, I_DO_UPDATE = ', 
     &             I_DO_NEWTON, I_DO_UPDATE
*
*. Memory for information on convergence of iterative procedure
      NITEM = 4
      LEN_SUMMARY = NITEM*(MAXMAC+1)
      CALL MEMMAN(KL_SUMMARY,LEN_SUMMARY,'ADDL  ',2,'SUMMRY')
*. Memory for the initial set of MO integrals
      CALL MEMMAN(KINT1_INI,NINT1,'ADDL  ',2,'INT1_IN')
      CALL MEMMAN(KINT2_INI,NINT2,'ADDL  ',2,'INT2_IN')
*. And for two extra MO matrices 
      LEN_CMO =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,0)
      CALL MEMMAN(KLMO1,LEN_CMO,'ADDL  ',2,'MO1   ')
      CALL MEMMAN(KLMO2,LEN_CMO,'ADDL  ',2,'MO2   ')
      CALL MEMMAN(KLMO3,LEN_CMO,'ADDL  ',2,'MO3   ')
      CALL MEMMAN(KLMO4,LEN_CMO,'ADDL  ',2,'MO4   ')
*. And for storing MO coefficients from outer iteration
      CALL MEMMAN(KLMO_OUTER,LEN_CMO,'ADDL  ',2,'MOOUTE')
*. And initial set of MO's
      CALL MEMMAN(KLCMOAO_INI,LEN_CMO,'ADDL  ',2,'MOINI ')
*. Normal integrals accessed
      IH1FORM = 1
      I_RES_AB = 0
      IH2FORM = 1
*. CI not CC
      ICC_EXC = 0
*. Use of supersymmetry in generation of excitations?
      IF(I_USE_SUPSYM.EQ.1.AND.INCLUDE_ONLY_TOTSYM_SUPSYM.EQ.1) THEN
        I_RESTRICT_SUPSYM = 1
      ELSE 
        I_RESTRICT_SUPSYM = 0
      END IF
*
* ============================================
* Initialize info for Dual space (PQ ) set up
* ============================================
*
*. Code requires use of CSFS and ICNFBAT = 2, Check
*
* Test that CSFs are in use
      IF(NOCSF.EQ.1) THEN
        WRITE(6,*) ' AKBK calculation with determinants is not enabled'
        WRITE(6,*) ' Please remove NOCSF keyword '
        STOP       ' AKBK calculation with determinants is not enabled'
      END IF
      IF(ICNFBAT.EQ.1.OR.ICISTR.EQ.1) THEN
        WRITE(6,*) ' AKBK calculation does not work with incore options'
        WRITE(6,*) ' Please make sure ICNFBAT .ge. 2 and ICISTR .ge. 2 '
        STOP       ' AKBK calculation does not work with incore options'
      END IF
*
      IPSPC = IREFSPC_MCSCF-1
      IQSPC = IREFSPC_MCSCF
      WRITE(6,*) ' IPSPC, IQSPC = ', IPSPC, IQSPC
     
*. Number of CSFs in P and Q
      NPCSF = NCSF_PER_SYM_GN(IREFSM,IPSPC)
      NQCSF = NCSF_PER_SYM_GN(IREFSM,IQSPC)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Number of CSFs in P-space = ', NPCSF
        WRITE(6,*) ' Number of CSFs in PQ-space = ', NQCSF
      END IF
      EREF_INI = EREF
*
*=====================================
* Construction of orbital excitations 
*=====================================
*
* In the present context, we will have two sets of orbital excitations.
*   The standard orbital excitations (interspace), named e
*   The extra rotations connected with constraints (intraspace) named d
*
* We will need to work with thise together and alone. This is 
* realized by ordering the excitations so
*  1) The e-excitations come first
*  2) The c-excitations then show up
*
*
*. The number of e-excitations
      CALL MEMMAN(KLTTACT,(NGAS+2)**2,'ADDL  ',1,'TTACT ')
      CALL NONRED_TT_EXC(int_mb(KLTTACT),IREFSPC_MCSCF,1)
*. Number of Nonredundant orbital excitations
      KLOOEXC = 1
      KLOOEXCC= 1
      CALL NONRED_OO_EXC2(NOOEXC,int_mb(KLOOEXC),int_mb(KLOOEXCC),
     &     1,int_mb(KLTTACT),I_RESTRICT_SUPSYM,int_mb(KMO_SUPSYM),
     &     N_INTER_EXC,N_INTRA_EXC,1)
      NEOOEXC = N_INTER_EXC
      NCOOEXC = N_INTRA_EXC
*
      IF(NOOEXC.EQ.0) THEN
        WRITE(6,*) ' STOP: zero orbital excitations in MCSCF '
        STOP       ' STOP: zero orbital excitations in MCSCF '
      END IF
*.. And the excitations
      CALL MEMMAN(KLOOEXC,NTOOB*NTOOB,'ADDL  ',1,'OOEXC ')
      CALL MEMMAN(KLOOEXCC,2*NOOEXC,'ADDL  ',1,'OOEXCC')
*. Allow these parameters to be known outside
      KIOOEXC = KLOOEXC
      KIOOEXCC = KLOOEXCC
      CALL NONRED_OO_EXC2(NOOEXC,int_mb(KLOOEXC),int_mb(KLOOEXCC),
     &     1,int_mb(KLTTACT),I_RESTRICT_SUPSYM,int_mb(KMO_SUPSYM),
     &     N_INTER_EXC,N_INTRA_EXC,2)
*
* =========================================================================
* Memory allocation - why now: We needed to know the number of excitations
* =========================================================================
*
*. Memory for gradient 
      CALL MEMMAN(KLE1,NOOEXC,'ADDL  ',2,'E1_MC ')
*. And Brillouin matrix in complete form
      CALL MEMMAN(KLBR,LEN_CMO,'ADDL  ',2,'BR_MAT')
*. And an extra gradient
      CALL MEMMAN(KLE1B,NOOEXC,'ADDL  ',2,'E1B   ')
*. and a scratch vector for excitation
      CALL MEMMAN(KLEXCSCR,NOOEXC,'ADDL  ',2,'EX_SCR')
*. Memory for gradient and orbital-Hessian - if  required
      IF(IOOE2_APR.EQ.1) THEN
        LE2 = NOOEXC*(NOOEXC+1)/2
        CALL MEMMAN(KLE2,LE2,'ADDL  ',2,'E2_MC ')
*. For eigenvectors of orbhessian
        LE2F = NOOEXC**2
        CALL MEMMAN(KLE2F,LE2F,'ADDL  ',2,'E2_MCF')
*. and eigenvalues, scratch, kappa
        CALL MEMMAN(KLE2VL,NOOEXC,'ADDL  ',2,'EIGVAL')
      ELSE
        KLE2 = -1
        KLE2F = -1
        KLE2VL = -1
      END IF
      KLIBENV = -2810
      KCLKSCR = -2810
*
      IF(I_DO_UPDATE.EQ.1) THEN
*. Space for update procedure
*. Array defining envelope and a scratch vector
* and matrix
        CALL MEMMAN(KLIBENV,NOOEXC,'ADDL  ',2,'IBENV')
        CALL MEMMAN(KLCLKSCR,NOOEXC,'ADDL  ',2,'CLKSCR')
*. rank 2 matrices
        CALL MEMMAN(KLRANK2,4*NMAX_VEC_UPDATE,'ADDL  ',2,'RNK2MT')
* two vectors defining each rank two-space
        LENGTH_V = 2*NMAX_VEC_UPDATE*NOOEXC
        CALL MEMMAN(KLUPDVEC,LENGTH_V,'ADDL  ',2,'RNK2VC')
*. Vectors for saving previous kappa and gradient
        CALL MEMMAN(KLE1PREV,NOOEXC,'ADDL  ',2,'E1PREV')
        CALL MEMMAN(KLKPPREV,NOOEXC,'ADDL  ',2,'KPPREV')
      END IF
*. and scratch, kappa
      CALL MEMMAN(KLE2SC,NOOEXC,'ADDL  ',2,'EIGSCR')
      CALL MEMMAN(KLKAPPA,NOOEXC,'ADDL  ',2,'KAPPA ')
      CALL MEMMAN(KLSTEP, NOOEXC,'ADDL  ',2,'STEP  ')
*. Save the initial set of MO integrals 
      CALL COPVEC(WORK(KINT1O),WORK(KINT1_INI),NINT1)
      CALL COPVEC(WORK(KINT2) ,WORK(KINT2_INI),NINT2)
      IF(I_DO_CONSTR.EQ.1) THEN
*. Include space for treating orbital constraints: For the F2 matrix
*. and partitionings of F2, E2, and E1.
        LFE2 = (NOOEXC+NTOOB)*NOOEXC
        CALL MEMMAN(KLFE2,LFE2,'ADDL  ',2,'FE2   ')
        LCC = NCOOEXC*NCOOEXC
        LCE = NCOOEXC*NEOOEXC
        LEE = NEOOEXC*NEOOEXC
        CALL MEMMAN(KLFE2CC,LCC,'ADDL  ',2,'FE2CC ')
        CALL MEMMAN(KLFE2CE,LCE,'ADDL  ',2,'FE2CC ')
        CALL MEMMAN(KLFE2EC,LCE,'ADDL  ',2,'FE2CC ')
        CALL MEMMAN(KLFE2EE,LEE,'ADDL  ',2,'FE2CC ')
        CALL MEMMAN(KLE2CC,LCC,'ADDL  ',2,'E2CC  ')
        CALL MEMMAN(KLE2CE,LCE,'ADDL  ',2,'E2CC  ')
        CALL MEMMAN(KLE2EC,LCE,'ADDL  ',2,'E2CC  ')
        CALL MEMMAN(KLE2EE,LEE,'ADDL  ',2,'E2CC  ')
        CALL MEMMAN(KLE1C,NCOOEXC,'ADDL  ',2,'E1C   ')
        CALL MEMMAN(KLE1E,NEOOEXC,'ADDL  ',2,'E1E   ')
        LEN_H1 =  NDIM_1EL_MAT(1,NTOOBS,NTOOBS,NSMOB,1)
        CALL MEMMAN(KLFIA,LENH1,'ADDL  ',2,'FIA    ')
        CALL MEMMAN(KLF2OFF,NSMOB**2,'ADDL  ',1,'F2OFF ')
      END IF
*. Print will be reduced for densities
      IPRDEN_SAVE = IPRDEN
C     IPRDEN = 0
      IRESTR_SAVE = IRESTR
*
      IIUSEH0P = 0
      MPORENP_E = 0
      IPRDIAL = IPRMCSCF
*
      CONVER = .FALSE.
      CONV_F = .FALSE.
*. The various types of integral lists- should probably be made in
* def of lists
      IE2LIST_0F = 1
      IE2LIST_1F = 2
      IE2LIST_2F = 3
      IE2LIST_4F = 5
*. For integral transformation: location of MO coefs
      KKCMO_I = KMOMO
      KKCMO_J = KMOMO
      KKCMO_K = KMOMO
      KKCMO_L = KMOMO
*
      IF(I_DO_UPDATE.EQ.1) THEN
*. Define envelope for used orbital Hessian - pt complete
* is constructed so
        IONE = 1
        CALL ISETVC(WORK(KLIBENV),IONE,NOOEXC)
      END IF
*
*. Loop over outer iterations
*
* In summery
* 1: Norm of orbgradient
* 2: Norm of orbstep
* 3: Norm of CI after iterative procedure
* 4: Energy
*
*. Convergence is pt  energy change le THRES_E
*
      ZERO = 0.0D0
      NMAT_UPD = 0
*. Line search is not meaning full very close to convergence
      THRES_FOR_ENTER_LINSEA = 1.0D-8

      N_INNER_TOT = 0
      DO IOUT = 1, MAXMAC
*
        IF(IPRNT.GE.1) THEN
          WRITE(6,*)
          WRITE(6,*) ' ----------------------------------'
          WRITE(6,*) ' Output from outer iteration', IOUT
          WRITE(6,*) ' ----------------------------------'
          WRITE(6,*)
        END IF
        CALL MEMCHK2('ST_OUT')
        NOUTIT = IOUT
*
*. Transform integrals to current set of MO's
*
        IF(IPRNT.GE.10) WRITE(6,*) ' Integral transformation:' 
        KINT2 = KINT_2EMO
        CALL COPVEC(WORK(KINT1_INI),WORK(KINT1O),NINT1)
        CALL COPVEC(WORK(KINT2_INI),WORK(KINT2),NINT2)
*. Flag type of integral list to be obtained
C       IE2LIST_A, IOCOBTP_A,INTSM_A
*. Flag for integrals with Two  free index: energy + gradient+orb-Hessian
*. Check problem: raise!!
        IE2LIST_A = IE2LIST_2F
        IE2LIST_A = IE2LIST_4F
        IOCOBTP_A = 2
*. Check, end
        INTSM_A = 1
        CALL TRAINT
*
        CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
*. Calculate inactive Fockmatrix
*. Calculate inactive Fock matrix from integrals over initial orbitals
*
*. A problem with the modern integral structure: the code will look for 
*. a list of full two-electron integrals and will use this, rather than the 
*. above definition. Well, place pointer KINT2_INI at relevant place
        IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
        KINT2_FSAVE = KINT2_A(IE2ARR_F)
        KINT2_A(IE2ARR_F) = KINT2_INI
C            FI_FROM_INIINT(FI,CINI,H,EINAC,IHOLETP)
        CALL FI_FROM_INIINT(WORK(KFI),WORK(KMOMO),WORK(KH),
     &                      ECORE_HEX,3)
        ECORE = ECORE_ORIG + ECORE_HEX
        CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
        IF(NTEST.GE.10000) THEN
          WRITE(6,*) ' MCSCF: ECORE_ORIG, ECORE_HEX, ECORE(2) ',
     &                 ECORE_ORIG, ECORE_HEX, ECORE
        END IF
*. and   redirect integral fetcher back to actual integrals
        KINT2 = KINT_2EMO
        KINT2_A(IE2ARR_F) = KINT2_FSAVE
*. The diagonal will fetch J and K integrals using GTIJKL_GN,* 
*. prepare for this routine
        IE2ARRAY_A = IE2LIST_I(IE2LIST_IB(IE2LIST_A))
*
*. Perform CI - and calculate densities
*
        IF(IPRNT.GE.10) WRITE(6,*) ' CI: '
*. At most MAXMIC iterations
        IF(IOUT.NE.1) IRESTR = 1
     
        MAXIT_SAVE = MAXIT
        IF(I_DO_APPROX_CI.EQ.0) THEN
          CALL GASCI(IREFSM,IREFSPC_MCSCF,IPRDIAL,IIUSEH0P,
     &             MPORENP_E,EREF,ERROR_NORM_FINAL,CONV_F)  
        ELSE
           CALL AKBKCI(IREFSPC_MCSCF,IPRDIAL,
     &           EREF,ERROR_NORM_FINAL,CONV_F)
        END IF
        MAXIT = MAXIT_SAVE
        WRITE(6,*) ' Energy and residual from CI :', 
     &  EREF,ERROR_NORM_FINAL
        WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+3) = ERROR_NORM_FINAL
        WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+4) = EREF
        EOLD = EREF
        ENEW = EREF
*. (Sic)
*
        IF(IOUT.GT.1) THEN
*. Check for convergence
          DELTA_E = WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+4)-
     &              WORK(KL_SUMMARY-1+(IOUT-1-1)*NITEM+4)
          IF(IPRNT.GE.2) WRITE(6,'(A,E9.3)') 
     &    '  Change of energy between outer iterations = ', DELTA_E
          IF(ABS(DELTA_E).LE.THRES_E) CONVER = .TRUE.
        END IF
        IF(CONVER) THEN
          NOUTIT = NOUTIT-1
          IF(IPRNT.GE.1) THEN
            WRITE(6,*) ' MCSCF calculation has converged'
          END IF
          GOTO 1001
        END IF
*. A test
C       CALL EN_FROM_DENS(ENERGY,2,0)
        CALL EN_FROM_DENS(ENERGY2,2,0)
        WRITE(6,*) ' Energy from density matrices ', ENERGY2
*. The active Fock matrix
        IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
        KINT2_FSAVE = KINT2_A(IE2ARR_F)
        KINT2_A(IE2ARR_F) = KINT2_INI
        CALL FA_FROM_INIINT
     &       (WORK(KFA),WORK(KMOMO),WORK(KMOMO),dbl_mb(KRHO1),1)
*. Clean up
        KINT2 = KINT_2EMO
        KINT2_A(IE2ARR_F) = KINT2_FSAVE
*
*.======================================
*. Exact or approximate orbital Hessian 
*.======================================
*
*
*. Fock matrix in KF
          CALL FOCK_MAT_STANDARD(WORK(KF),2,WORK(KFI),WORK(KFA))
        IOOSM = 1
C            ORBHES(OOHES,IOOEXC,NOOEXC,IOOSM,ITTACT)
        IF(IOOE2_APR.EQ.1) THEN
          CALL ORBHES(WORK(KLE2),int_mb(KLOOEXC),NOOEXC,IOOSM,
     &         int_mb(KLTTACT))
          IF(NTEST.GE.1000) THEN
           WRITE(6,*) ' The orbital Hessian '
           CALL PRSYM(WORK(KLE2),NOOEXC)
          END IF
        END IF
*
        IF(I_DO_CONSTR.EQ.1) THEN
*. Calculate matrices to ensure that the constraints - FI + FA is 
*. block-diagonal
C     GET_FIA2(F2,FIA,IFIA2_OFFSS)
       END IF
   
*. 
        WRITE(6,*) ' Jeppe enforced me to to stop '
        STOP ' Jeppe enforced me to to stop '
*
*. Diagonalize to determine lowest eigenvalue
*
*. Outpack to complete form
        CALL TRIPAK(WORK(KLE2F),WORK(KLE2),2,NOOEXC,NOOEXC)
C            TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
*. Lowest eigenvalue
C            DIAG_SYMMAT_EISPACK(A,EIGVAL,SCRVEC,NDIM,IRETURN)
        CALL DIAG_SYMMAT_EISPACK(WORK(KLE2F),WORK(KLE2VL),
     &       WORK(KLE2SC),NOOEXC,IRETURN)
        IF(IRETURN.NE.0) THEN
           WRITE(6,*) 
     &     ' Problem with diagonalizing E2, IRETURN =  ', IRETURN
        END IF
        IF(IPRNT.GE.1000) THEN
          WRITE(6,*) ' Eigenvalues: '
          CALL WRTMAT(WORK(KLE2VL),1,NOOEXC,1,NOOEXC)
        END IF
*. Lowest eigenvalue
        E2VL_MN = XMNMX(WORK(KLE2VL),NOOEXC,1)
        IF(IPRNT.GE.2)  
     &  WRITE(6,*) ' Lowest eigenvalue of E2(orb) = ', E2VL_MN
*
*. Cholesky factorization orbital Hessian if required
*
        I_SHIFT_E2 = 0
        IF(I_DO_UPDATE.EQ.1) THEN
*. Cholesky factorization requires positive matrices.
*. add a constant to diagonal if needed
          XMINDIAG = 1.0D-4
          IF(E2VL_MN.LE.XMINDIAG) THEN
           ADD = XMINDIAG - E2VL_MN 
C               ADDDIA(A,FACTOR,NDIM,IPACK)
           CALL ADDDIA(WORK(KLE2),ADD,NOOEXC,1)
           I_SHIFT_E2 = 1
          END IF
C CLSKHE(AL,X,B,NDIM,IB,IALOFF,ITASK,INDEF)
C         WRITE(6,*) ' NOOEXC before CLSKHE = ', NOOEXC 
          CALL CLSKHE(WORK(KLE2),XDUM,XDUM,NOOEXC,WORK(KLIBENV),
     &         WORK(KLCLKSCR),1,INDEF)
          IF(INDEF.NE.0) THEN
            WRITE(6,*) ' Indefinite matrix in CKSLHE '
            STOP ' Indefinite matrix in CKSLHE '
          END IF
        END IF! Cholesky decomposition required
*
*
*. Finite difference check
*
        I_DO_FDCHECK = 0
        IF(I_DO_FDCHECK.EQ.1) THEN
*. First: Analytic gradient from Fock matrix - As kappa = 0, Brillouin vector
* = gradient
          CALL E1_FROM_F(WORK(KLE1),WORK(KF),1,int_mb(KLOOEXC),
     &                   int_mb(KLOOEXCC),
     &                   NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
*
          CALL MEMMAN(KLE1FD,NOOEXC,'ADDL  ',2,'E1_FD ')
          LE2 = NOOEXC*NOOEXC
          CALL MEMMAN(KLE2FD,LE2,   'ADDL  ',2,'E2_FD ')
          CALL SETVEC(WORK(KLE2VL),ZERO,NOOEXC)
          CALL GENERIC_GRA_HES_FD(E0,WORK(KLE1FD),WORK(KLE2FD),
     &         WORK(KLE2VL),NOOEXC,EMCSCF_FROM_KAPPA)
C              GENERIC_GRA_HES_FD(E0,E1,E2,X,NX,EFUNC)
*. Compare gradients
          ZERO = 0.0D0
          CALL CMP2VC(WORK(KLE1FD),WORK(KLE1),NOOEXC,ZERO)
*. transform Finite difference Hessian to packed form
          CALL TRIPAK(WORK(KLE2FD),WORK(KLE2F),1,NOOEXC,NOOEXC)
          LEN = NOOEXC*(NOOEXC+1)/2
          CALL CMP2VC(WORK(KLE2),WORK(KLE2F),LEN,ZERO)
              STOP ' Enforced stop after FD check'
        END IF
*       ^ End of finite difference check
*. Initialize sum of steps for outer iteration
        WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+2) = 0.0D0
*. Loop over Inner iterations, where orbitals are optimized
*. Initialize Kappa as zero
        IF(IRESET_KAPPA_IN_OR_OUT.EQ.2) THEN
          CALL SETVEC(WORK(KLKAPPA),ZERO,NOOEXC)
        END IF
*. Save MO's from start of each outer iteration
        CALL COPVEC(WORK(KMOMO),WORK(KMOREF),LEN_CMO)
*. Convergence Threshold for inner iterations
*. At the moment just chosen as the total convergence threshold
        THRES_E_INNER = THRES_E
        CONV_INNER = .FALSE.
        DO IINNER = 1, MAXMIC
          N_INNER_TOT = N_INNER_TOT + 1
*
          IF(IPRNT.GE.5) THEN
            WRITE(6,*)
            WRITE(6,*) ' Info from inner iteration = ', IINNER
            WRITE(6,*) ' ===================================='
            WRITE(6,*)
          END IF
*
          IF(IRESET_KAPPA_IN_OR_OUT.EQ.1) THEN
            CALL SETVEC(WORK(KLKAPPA),ZERO,NOOEXC)
          END IF
          E_INNER_OLD = EREF
          EOLD = ENEW
*
          IF(IINNER.NE.1) THEN
*
*. gradient integral transformation and Fock matrices
*
*. Flag type of integral list to be obtained:
*. Flag for integrals with one free index: energy + gradient
           IE2LIST_A = IE2LIST_1F
           IOCOBTP_A = 1
           INTSM_A = 1
           CALL TRAINT
           CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
*
*. Calculate inactive and active Fock matrix from integrals over 
*  initial orbitals
*. Redirect integral fetcher to initial integrals- for old storage mode
           KINT2 = KINT2_INI
*. A problem with the modern integral structure: the code will look for 
*. a list of full two-electron integrals and will use this, rather than the 
*. above definition. Well, place pointer KINT2_INI at relevant place
           IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
           KINT2_FSAVE = KINT2_A(IE2ARR_F)
           KINT2_A(IE2ARR_F) = KINT2_INI
C             FI_FROM_INIINT(FI,CINI,H,EINAC,IHOLETP)
           CALL FI_FROM_INIINT(WORK(KFI),WORK(KMOMO),WORK(KH),
     &                         ECORE_HEX,3)
           ECORE = ECORE_ORIG + ECORE_HEX
           CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' ECORE_ORIG, ECORE_HEX, ECORE(2) ',
     &                    ECORE_ORIG, ECORE_HEX, ECORE
           END IF
           CALL FA_FROM_INIINT
     &     (WORK(KFA),WORK(KMOMO),WORK(KMOMO),dbl_mb(KRHO1),1)
*. and   redirect integral fetcher back to actual integrals
           KINT2 = KINT_2EMO
           KINT2_A(IE2ARR_F) = KINT2_FSAVE
*. Fock matrix in KF
          CALL FOCK_MAT_STANDARD(WORK(KF),2,WORK(KFI),WORK(KFA))
          END IF ! IINNER .ne.1
*
*. Construct orbital gradient
*
          IF(IPRNT.GE.10) WRITE(6,*) ' Construction of E1: '
          XKAPPA_NORM = SQRT(INPROD(WORK(KLKAPPA),WORK(KLKAPPA),NOOEXC))
          IF(I_USE_BR_OR_E1.EQ.1.OR.XKAPPA_NORM.EQ.0.0D0) THEN
*. Brillouin vector from Fock matrix is used
           CALL E1_FROM_F(WORK(KLE1),WORK(KF),1,int_mb(KLOOEXC),
     &                   int_mb(KLOOEXCC),
     &                   NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
          ELSE
*. Calculate gradient at non-vanishing Kappa
*. Complete Brillouin matrix
C              GET_BRT_FROM_F(BRT,F)
          CALL GET_BRT_FROM_F(WORK(KLBR),WORK(KF))
C              E1_MCSCF_FOR_GENERAL_KAPPA(E1,F,KAPPA)
          CALL E1_MCSCF_FOR_GENERAL_KAPPA(WORK(KLE1),WORK(KLBR),
     &         WORK(KLKAPPA))
          END IF
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' E1, Gradient: '
            CALL WRTMAT(WORK(KLE1),1,NOOEXC,1,NOOEXC)
          END IF
*
          E1NRM = SQRT(INPROD(WORK(KLE1),WORK(KLE1),NOOEXC))
          IF(IPRNT.GE.2) WRITE(6,*) ' Norm of orbital gradient ', E1NRM
          WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+1) = E1NRM
*
* ==========================
* Two step Newton procedure
* ==========================
*
          IF(I_DO_NEWTON.EQ.1) THEN
*
*. Transform gradient to diagonal basis
*
*. (save original gradient)
            CALL COPVEC(WORK(KLE1),WORK(KLE1B),NOOEXC)
            CALL MATVCC(WORK(KLE2F),WORK(KLE1),WORK(KLE2SC),
     &           NOOEXC,NOOEXC,1)
            CALL COPVEC(WORK(KLE2SC),WORK(KLE1),NOOEXC)
*
*. Solve shifted NR equations with step control
*
*           SOLVE_SHFT_NR_IN_DIAG_BASIS(
*    &            E1,E2,NDIM,STEP_MAX,TOLERANCE,X,ALPHA)A
            CALL SOLVE_SHFT_NR_IN_DIAG_BASIS(WORK(KLE1),WORK(KLE2VL),
     &           NOOEXC,STEP_MAX,TOLER,WORK(KLSTEP),ALPHA,DELTA_E_PRED)
            XNORM_STEP = SQRT(INPROD(WORK(KLSTEP),WORK(KLSTEP),NOOEXC))
*. Is step close to max
            I_CLOSE_TO_MAX = 0 
            IF(0.8D0.LE.XNORM_STEP/STEP_MAX) I_CLOSE_TO_MAX  = 1
*
            WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+2) = 
     &      WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+2) + XNORM_STEP
            IF(IPRNT.GE.2) WRITE(6,'(A,2(2X,E12.5))')
     &      ' Norm of step and predicted energy change = ',
     &       XNORM_STEP, DELTA_E_PRED
*. transform step to original basis
            CALL MATVCC(WORK(KLE2F),WORK(KLSTEP),WORK(KLE2SC),
     &           NOOEXC,NOOEXC,0)
            CALL COPVEC(WORK(KLE2SC),WORK(KLSTEP),NOOEXC)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Step in original basis:'
              CALL WRTMAT(WORK(KLSTEP),1,NOOEXC,1,NOOEXC)
            END IF
*. Is direction down-hills
            E1STEP = INPROD(WORK(KLSTEP),WORK(KLE1B),NOOEXC)
            IF(IPRNT.GE.2) WRITE(6,'(A,E12.5)')
     &      ' < E1!Step> = ', E1STEP
            IF(E1STEP.GT.0.0D0) THEN
             WRITE(6,*) ' Warning: step is in uphill direction '
            END IF
*. Energy for rotated orbitals
*
            ONE = 1.0D0
            CALL VECSUM(WORK(KLKAPPA),WORK(KLKAPPA),WORK(KLSTEP),
     &      ONE,ONE,NOOEXC)
            XNORM2 = SQRT(INPROD(WORK(KLKAPPA),WORK(KLKAPPA),NOOEXC))
            WRITE(6,*) ' Norm of updated kappa step =', XNORM2
            ENERGY1 = EMCSCF_FROM_KAPPA(WORK(KLKAPPA))
            ENEW = ENERGY1
            WRITE(6,*) ' Energy for rotated orbitals', ENERGY1
*. Compare old and new energy to decide with to do
            DELTA_E_ACT = ENEW-EOLD
            E_RATIO = DELTA_E_ACT/DELTA_E_PRED  
            IF(IPRNT.GE.2) WRITE(6,'(A,3(2X,E12.5))') 
     &      ' Actual and predicted energy change, ratio ', 
     &      DELTA_E_ACT, DELTA_E_PRED,E_RATIO
*
            IF(E_RATIO.LT.0.0D0) THEN
             WRITE(6,*) ' Trustradius reduced '
             RED_FACTOR = 2.0D0
             STEP_MAX = STEP_MAX/RED_FACTOR
             WRITE(6,*) ' New trust-radius ', STEP_MAX
            END IF
            IF(E_RATIO.GT.0.8D0.AND.I_CLOSE_TO_MAX.EQ.1) THEN
             WRITE(6,*) ' Trustradius increased '
             XINC_FACTOR = 2.0D0
             STEP_MAX = STEP_MAX*XINC_FACTOR
             WRITE(6,*) ' New trust-radius ', STEP_MAX
            END IF
*
            IF((ABS(DELTA_E_ACT).GT.THRES_FOR_ENTER_LINSEA).AND.
     &         (I_DO_LINSEA_MCSCF.EQ.1.OR.
     &         I_DO_LINSEA_MCSCF.EQ.2.AND.EOLD.GT.ENEW)) THEN
*
*. line-search for orbital optimization
*
C                 LINES_SEARCH_BY_BISECTION(FUNC,REF,DIR,NVAR,XINI,
C    &            XFINAL,FFINAL,IKNOW,F0,FXINI)
*. Step was added to Kappa when calculating energy, get Kappa back
              ONE = 1.0D0
              ONEM = -1.0D0
              CALL VECSUM(WORK(KLKAPPA),WORK(KLKAPPA),WORK(KLSTEP),
     &        ONE,ONEM,NOOEXC)
              CALL LINES_SEARCH_BY_BISECTION(EMCSCF_FROM_KAPPA,
     &             WORK(KLKAPPA),WORK(KLSTEP),NOOEXC,ONE,XFINAL,FFINAL,
     &             2, EOLD, ENEW)
              ENEW = FFINAL
              IF(IPRNT.GE.2) WRITE(6,*) ' Line search value of X = ',
     &        XFINAL
              XKAPPA_NORM2 = 
     &        SQRT(INPROD(WORK(KLKAPPA),WORK(KLKAPPA),NOOEXC))
              CALL VECSUM(WORK(KLKAPPA),WORK(KLKAPPA),WORK(KLSTEP),
     &             ONE, XFINAL,NOOEXC)
            END IF! line search should be entered
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Updated total Kappa '
              CALL WRTMAT(WORK(KLKAPPA),1,NOOEXC,1,NOOEXC)
            END IF
          END IF! Newton method
          CALL MEMCHK2('AF_NEW')
          IF(I_DO_UPDATE.EQ.1) THEN
*
* ====================
*  Update procedure
* ====================
*
*. Update Hessian
            IF(IINNER.EQ.1) THEN
*. Just save current info
              CALL COPVEC(WORK(KLE1),WORK(KLE1PREV),NOOEXC)
              CALL COPVEC(WORK(KLKAPPA),WORK(KLKPPREV),NOOEXC)
              NMAT_UPD = 0
            ELSE
C             HESUPV (E2,AMAT,AVEC,
C    &                 X,E1,VEC2,
C    &                 VEC3,NVAR,IUPDAT,IINV,VEC1,NMAT,
C    &                 LUHFIL,DISCH,IHSAPR,IBARR,E2,VEC4)
C            HESUPV (HDIAG,A,AVEC,X,G,XPREV,GPREV,NVAR,
C    &                   IUPDAT,IINV,SCR,NMAT,LUHFIL,DISCH,
C    &                   IHSAPR,IB,E2,VEC4)

*. Update on inverse
              IINV = 1
*. Initial approximation is a cholesky factorized matrix
              IHSAPR = 3
              CALL HESUPV(WORK(KLE2),WORK(KLRANK2),WORK(KLUPDVEC),
     &             WORK(KLKAPPA),WORK(KLE1),WORK(KLKPPREV),
     &             WORK(KLE1PREV),NOOEXC,I_UPDATE_MET,IINV,
     &             WORK(KLCLKSCR),NMAT_UPD,LUHFIL,DISCH,IHSAPR,
     &             WORK(KLIBENV),WORK(KLE2),WORK(KLEXCSCR))
*. Forget the first(when starting out with exact Hessian)
              NMAT_UPD = NMAT_UPD + 1
            END IF! IINNER = 1
*
*. New search direction = step
*==============================
*
*. Inverse initial Hessian approximation times gradient
            IF(IHSAPR.EQ.1) THEN
*. Just inverse diagonal (in E2) times gradient
              CALL VVTOV(WORK(KLE2),WORK(KLE1),WORK(KLSTEP),NOOEXC)
            ELSE
              CALL COPVEC(WORK(KLE1),WORK(KLCLKSCR),NOOEXC)
C                  CLSKHE(AL,X,B,NDIM,IB,IALOFF,ITASK,INDEF)
              CALL CLSKHE(WORK(KLE2),WORK(KLSTEP),WORK(KLCLKSCR),
     &             NOOEXC,WORK(KLIBENV),WORK(KLEXCSCR),2,INDEF)
            END IF
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' Contribution from H(ini) to (-1) step:'
              CALL WRTMAT(WORK(KLSTEP),1,NOOEXC,1,NOOEXC)
            END IF
*. And the rank-two updates
            IF(NMAT_UPD.NE.0) THEN
C                SLRMTV(NMAT,NVAR,A,AVEC,NRANK,VECIN,VECOUT,IZERO,
C    &                  DISCH,LUHFIL)
              IZERO = 0
              CALL SLRMTV(NMAT_UPD,NOOEXC,WORK(KLRANK2),WORK(KLUPDVEC),
     &                    2,WORK(KLE1),WORK(KLSTEP),IZERO,DISCH,LUHFIL)
            END IF
*. And the proverbial minus 1
            ONEM = -1.0D0
            CALL SCALVE(WORK(KLSTEP),ONEM,NOOEXC)
*. Check norm and reduce to STEP_MAX if required
            STEP_NORM = SQRT(INPROD(WORK(KLSTEP),WORK(KLSTEP),NOOEXC))
            IF(STEP_NORM.GT.STEP_MAX) THEN
              FACTOR = STEP_MAX/STEP_NORM
              IF(IPRNT.GE.2) 
     &        WRITE(6,'(A,E8.2)') ' Step reduced by factor = ', FACTOR
              CALL SCALVE(WORK(KLSTEP),FACTOR,NOOEXC)
            END IF
*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' Step:'
              CALL WRTMAT(WORK(KLSTEP),1,NOOEXC,1,NOOEXC)
            END IF
*. Is direction down-hills
            E1STEP = INPROD(WORK(KLSTEP),WORK(KLE1),NOOEXC)
            IF(IPRNT.GE.2) WRITE(6,'(A,E12.5)')
     &      '  < E1!Step> = ', E1STEP
            IF(E1STEP.GT.0.0D0) THEN
             WRITE(6,*) ' Warning: step is in uphill direction '
             WRITE(6,*) ' Sign of step is changed '
             ONEM = -1.0D0
             CALL SCALVE(WORK(KLSTEP),ONEM,NOOEXC)
            END IF
            XNORM_STEP = SQRT(INPROD(WORK(KLSTEP),WORK(KLSTEP),NOOEXC))
            WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+2) = 
     &      WORK(KL_SUMMARY-1+(IOUT-1)*NITEM+2) + XNORM_STEP
            IF(IPRNT.GE.2) WRITE(6,'(A,E12.5)')
     &      '  Norm of step  = ', XNORM_STEP
*
*. Determine step length along direction
*. ======================================
*
*. Energy for rotated orbitals
*
            ONE = 1.0D0
            CALL VECSUM(WORK(KLKAPPA),WORK(KLKAPPA),WORK(KLSTEP),
     &      ONE,ONE,NOOEXC)
            XNORM2 = SQRT(INPROD(WORK(KLKAPPA),WORK(KLKAPPA),NOOEXC))
            WRITE(6,'(A,E12.5)') 
     &      '  Norm of total kappa = ', XNORM2
            ENERGY1 = EMCSCF_FROM_KAPPA(WORK(KLKAPPA))
            ENEW = ENERGY1
            WRITE(6,*) ' Energy for rotated orbitals', ENERGY1
*. Compare old and new energy to decide with to do
            DELTA_E_ACT = ENEW-EOLD
            IF(IPRNT.GE.2) WRITE(6,'(A,3(2X,E9.3))') 
     &      '  Actual energy change without linesearch ', DELTA_E_ACT
*
            IF((ABS(DELTA_E_ACT).GT.THRES_FOR_ENTER_LINSEA).AND.
     &         (I_DO_LINSEA_MCSCF.EQ.1.OR.
     &         I_DO_LINSEA_MCSCF.EQ.2.AND.EOLD.GT.ENEW)) THEN
*
*. line-search for orbital optimization
*
*. Step was added to Kappa when calculating energy, get Kappa back
              ONE = 1.0D0
              ONEM = -1.0D0
              CALL VECSUM(WORK(KLKAPPA),WORK(KLKAPPA),WORK(KLSTEP),
     &        ONE,ONEM,NOOEXC)
              CALL LINES_SEARCH_BY_BISECTION(EMCSCF_FROM_KAPPA,
     &             WORK(KLKAPPA),WORK(KLSTEP),NOOEXC,ONE,XFINAL,FFINAL,
     &             2, EOLD, ENEW)
              ENEW = FFINAL
              IF(IPRNT.GE.2) WRITE(6,'(A,E9.3)') 
     &        '  Step-scaling parameter from lineseach = ', XFINAL
              XKAPPA_NORM2 = 
     &        SQRT(INPROD(WORK(KLKAPPA),WORK(KLKAPPA),NOOEXC))
              CALL VECSUM(WORK(KLKAPPA),WORK(KLKAPPA),WORK(KLSTEP),
     &             ONE, XFINAL,NOOEXC)
              DELTA_E_ACT = ENEW-EOLD
              IF(IPRNT.GE.2) WRITE(6,'(A,3(2X,E9.3))') 
     &        '  Actual energy change with  linesearch ', DELTA_E_ACT
            END IF! line search should be entered
*    
            IF(ABS(DELTA_E_ACT).LT.THRES_E_INNER) THEN
             WRITE(6,*) ' Inner iterations converged '
             CONV_INNER = .TRUE.
            END IF
*
            IF(NTEST.GE.1000) THEN
               WRITE(6,*) ' Updated total Kappa '
               CALL WRTMAT(WORK(KLKAPPA),1,NOOEXC,1,NOOEXC)
            END IF
          END IF ! Update method
*
*=======================================
*. The new and improved MO-coefficients
*=======================================
*
*. Obtain exp(-kappa)
          CALL MEMCHK2('BE_NWM')
C              GET_EXP_MKAPPA(EXPMK,KAPPAP,IOOEXC,NOOEXC)
          CALL GET_EXP_MKAPPA(WORK(KLMO1),WORK(KLKAPPA),
     &                        int_mb(KLOOEXCC),NOOEXC)
          CALL MEMCHK2('AF_EMK')
*. And new MO-coefficients
          CALL MULT_BLOC_MAT(WORK(KLMO2),WORK(KMOREF),WORK(KLMO1),
     &         NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
          CALL COPVEC(WORK(KLMO2),WORK(KMOMO),LEN_CMO)
          CALL MEMCHK2('AF_ML1')
*. And the new MO-AO coefficients
C?        WRITE(6,*) '  KMOAO_ACT = ', KMOAO_ACT
          CALL MULT_BLOC_MAT(WORK(KMOAO_ACT),WORK(KMOAOIN),WORK(KMOMO),
     &       NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
          CALL MEMCHK2('AF_ML2')
          IF(IPRNT.GE.100) THEN
            WRITE(6,*) ' Updated MO-coefficients'
            CALL APRBLM2(WORK(KMOMO),NTOOBS,NTOOBS,NSMOB,0)
          END IF
          IF(IRESET_KAPPA_IN_OR_OUT.EQ.1) THEN
            CALL COPVEC(WORK(KMOMO),WORK(KMOREF),LEN_CMO)
          END  IF
          CALL MEMCHK2('AF_NWM')
*
*
*  ===========================================================
*. CI in inner its- should probably be moved (but not removed)
*  ===========================================================
*
          I_DO_CI_IN_INNER_ITS = 0
          IF(I_MAY_DO_CI_IN_INNER_ITS.EQ.1.AND.I_SHIFT_E2.EQ.0.AND.
     &       XNORM2.LT.XKAPPA_THRES.AND.IOUT.GE.MIN_OUT_IT_WITH_CI) THEN
            IF(IPRNT.GE.10) WRITE(6,*) ' CI in inner it '
            I_DO_CI_IN_INNER_ITS = 1
*
*. Transform integrals to current set of MO's
*
            IF(IPRNT.GE.10) WRITE(6,*) ' Integral transformation:' 
            KINT2 = KINT_2EMO
            CALL COPVEC(WORK(KINT1_INI),WORK(KINT1O),NINT1)
            CALL COPVEC(WORK(KINT2_INI),WORK(KINT2),NINT2)
*. Flag type of integral list to be obtained
C           IE2LIST_A, IOCOBTP_A,INTSM_A
*. Flag for integrals with zero free index: energy 
*. Problem: raise!!
            IE2LIST_A = IE2LIST_0F
            IE2LIST_A = IE2LIST_4F
            IOCOBTP_A = 2
            INTSM_A = 1
            CALL TRAINT
*
            CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
*. Calculate inactive Fock matrix from integrals over initial orbitals
            IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
            KINT2_FSAVE = KINT2_A(IE2ARR_F)
            KINT2_A(IE2ARR_F) = KINT2_INI
            CALL FI_FROM_INIINT(WORK(KFI),WORK(KMOMO),WORK(KH),
     &                          ECORE_HEX,3)
            ECORE = ECORE_ORIG + ECORE_HEX
            CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' ECORE_ORIG, ECORE_HEX, ECORE(2) ',
     &                     ECORE_ORIG, ECORE_HEX, ECORE
            END IF
            KINT2 = KINT_2EMO
            KINT2_A(IE2ARR_F) = KINT2_FSAVE
*. The diagonal will fetch J and K integrals using GTIJKL_GN,* 
*. prepare for this routine
            IE2ARRAY_A = IE2LIST_I(IE2LIST_IB(IE2LIST_A))
*
*. Perform CI - and calculate densities
*
            IF(IPRNT.GE.10) WRITE(6,*) ' CI: '
            IRESTR = 1
            MAXIT_SAVE = MAXIT
            CALL GASCI(IREFSM,IREFSPC_MCSCF,IPRDIAL,IIUSEH0P,
     &           MPORENP_E,EREF,ERROR_NORM_FINAL,CONV_F)  
            MAXIT = MAXIT_SAVE
            WRITE(6,*) ' Energy and residual from CI :', 
     &      EREF,ERROR_NORM_FINAL
            ENEW  = EREF
          END IF! CI in inner iterations
*
*. Obtain and block diagonalize FI+FA
*
          I_DIAG_FIFA = 0
          IF(I_DIAG_FIFA.EQ.1) THEN
            IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
            KINT2_FSAVE = KINT2_A(IE2ARR_F)
            KINT2_A(IE2ARR_F) = KINT2_INI
C                FI_FROM_INIINT(FI,CINI,H,EINAC,IHOLETP)
            CALL FI_FROM_INIINT(WORK(KFI),WORK(KMOMO),WORK(KH),
     &                      ECORE_HEX,3)
            ECORE = ECORE_ORIG + ECORE_HEX
            CALL FA_FROM_INIINT
     &      (WORK(KFA),WORK(KMOMO),WORK(KMOMO),dbl_mb(KRHO1),1)
*. Clean up
            KINT2_A(IE2ARR_F) = KINT2_FSAVE
*. Diagonalize FI+FA and save in KLMO2
            CALL VECSUM(WORK(KLMO1),WORK(KFI),WORK(KFA),ONE,ONE,NINT1)
            CALL DIAG_GASBLKS(WORK(KLMO1),WORK(KLMO2),
     &           IDUM,IDUM,IDUM,WORK(KLMO3),WORK(KLMO4),2)
*. And new MO-coefficients
            CALL MULT_BLOC_MAT(WORK(KLMO3),WORK(KMOMO),WORK(KLMO2),
     &           NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
            CALL COPVEC(WORK(KLMO3),WORK(KMOMO),LEN_CMO)
          END IF !FIFA should be diagonalized
*
          IF(CONV_INNER.AND.I_DO_CI_IN_INNER_ITS.EQ.1) THEN
            CONVER = .TRUE.
            GOTO 1001
          END IF
          IF(CONV_INNER) GOTO 901
        END DO !End of loop over inner iterations
 901    CONTINUE
        CALL MEMCHK2('EN_OUT')
      END DO
*     ^ End of loop over outer iterations
 1001 CONTINUE
      IF(CONVER) THEN
        WRITE(6,*) 
     &  ' Convergence of MCSCF was obtained in ', NOUTIT,' iterations'
      ELSE
        WRITE(6,*) 
     &  ' Convergence of MCSCF was not obtained in ', NOUTIT, 
     &  'iterations'
      END IF
      WRITE(6,'(A,I4)') 
     &'  Total number of inner iterations ', N_INNER_TOT
*
*
*. Finalize: Transform integrals to final MO's, obtain
*  norm of CI- and orbital gradient
*
*
*. Expansion of final orbitals in AO basis, pt in KLMO2
      CALL MULT_BLOC_MAT(WORK(KLMO2),WORK(KMOAOIN),WORK(KMOMO),
     &       NSMOB,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,0)
      CALL COPVEC(WORK(KLMO2),WORK(KMOAO_ACT),LEN_CMO)
      CALL COPVEC(WORK(KLMO2),WORK(KMOAOUT),LEN_CMO)
      WRITE(6,*) 
     &' Final MO-AO transformation stored in MOAOIN, MOAO_ACT, MOAOUT'
*. Integral transformation
      KINT2 = KINT_2EMO
      CALL COPVEC(WORK(KINT1_INI),WORK(KINT1O),NINT1)
      CALL COPVEC(WORK(KINT2_INI),WORK(KINT2),NINT2)
*. Flag for integrals with one free index: energy + gradient
      IE2LIST_A = IE2LIST_1F
      IE2LIST_A = IE2LIST_4F
      IOCOBTP_A = 1
      INTSM_A = 1
      CALL TRAINT
      CALL COPVEC(WORK(KINT1),WORK(KH),NINT1)
*. Calculate inactive Fockmatrix -
      KINT2 = KINT2_INI
      IF(ITRA_ROUTE.EQ.2) THEN
        IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
        KINT2_FSAVE = KINT2_A(IE2ARR_F)
        KINT2_A(IE2ARR_F) = KINT2_INI
      END IF
      CALL FI_FROM_INIINT(WORK(KFI),WORK(KMOMO),WORK(KH),
     &                    ECORE_HEX,3)
      IF(IPRNT.GE.100) WRITE(6,*) ' FI constructed '
      IF(ITRA_ROUTE.EQ.2) KINT2_A(IE2ARR_F) = KINT2_FSAVE
      CALL COPVEC(WORK(KFI),WORK(KINT1),NINT1)
      ECORE = ECORE_ORIG + ECORE_HEX
      KINT2 = KINT_2EMO
*. And 0 CI iterations with new integrals
      MAXIT_SAVE = MAXIT
      MAXIT = 1
      IRESTR = 1
*. and normal density print
      IPRDEN = IPRDEN_SAVE 
      CALL GASCI(IREFSM,IREFSPC_MCSCF,IPRDIA,IIUSEH0P,
     &            MPORENP_E,EREF,ERROR_NORM_FINAL_CI,CONV_F)
      EFINAL = EREF
      MAXIT = MAXIT_SAVE
*. Current orbital gradient
*. Active Fock matrix
      KINT2 = KINT2_INI
      IF(ITRA_ROUTE.EQ.2) THEN
        IE2ARR_F = IE2LIST_I(IE2LIST_IB(IE2LIST_FULL))
        KINT2_FSAVE = KINT2_A(IE2ARR_F)
        KINT2_A(IE2ARR_F) = KINT2_INI
      END IF
      CALL FA_FROM_INIINT
     &(WORK(KFA),WORK(KMOMO),WORK(KMOMO),dbl_mb(KRHO1),1)
      IF(IPRNT.GE.100) WRITE(6,*) ' FA constructed '
      KINT2 = KINT_2EMO
      IF(ITRA_ROUTE.EQ.2) KINT2_A(IE2ARR_F) = KINT2_FSAVE
*
      CALL FOCK_MAT_STANDARD(WORK(KF),2,WORK(KINT1),WORK(KFA))
      IF(IPRNT.GE.100) WRITE(6,*) ' F constructed '
      CALL E1_FROM_F(WORK(KLE1),WORK(KF),1,int_mb(KLOOEXC),
     &               int_mb(KLOOEXCC),
     &               NOOEXC,NTOOB,NTOOBS,NSMOB,IBSO,IREOST)
      E1NRM_ORB = SQRT(INPROD(WORK(KLE1),WORK(KLE1),NOOEXC))
      VNFINAL = E1NRM_ORB + ERROR_NORM_FINAL_CI
*
      IF(IPRORB.GE.2) THEN
        WRITE(6,*) 
     &  ' Final MOs in initial basis (not natural or canonical)'
        CALL APRBLM2(WORK(KMOMO),NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      IF(IPRORB.GE.1) THEN
        WRITE(6,*) 
     &  ' Final MOs in AO basis (not natural or canonical)'
        CALL PRINT_CMOAO(WORK(KLMO2))
      END IF
*
*. Projection of final occupied orbitals on initial set of occupied orbitals
*
*. Obtain initial and final occupied orbitals
      ISCR(1) = 0
      ISCR(2) = NGAS
      CALL MEMMAN(KLCOCC_INI,LEN_CMO,'ADDL  ',2,'COCC_IN')
      CALL MEMMAN(KLCOCC_FIN,LEN_CMO,'ADDL  ',2,'COCC_FI')
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,
      CALL CSUB_FROM_C(WORK(KMOAOIN),WORK(KLCOCC_INI),NOCOBS,ISCR_NTS,
     &                 2,ISCR,0)
      CALL CSUB_FROM_C(WORK(KLMO2),WORK(KLCOCC_FIN),NOCOBS,ISCR_NTS,
     &                 2,ISCR,0)
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,IONLY_DIM)
      WRITE(6,*) 
     &' Projecting final (MO2) on initial (MO1) occupied orbitals'
      CALL PROJ_ORBSPC_ON_ORBSPC(WORK(KLCOCC_INI),WORK(KLCOCC_FIN),
     &     NOCOBS,NOCOBS)
C     PROJ_ORBSPC_ON_ORBSPC(CMOAO1,CMOAO2,NMO1PSM,NMO2PSM)
*
*. Projection of final active orbitals on initial set of active orbitals
*
*. Obtain initial and final active orbitals
      ISCR(1) = NGAS
      CALL MEMMAN(KLCOCC_INI,LEN_CMO,'ADDL  ',2,'COCC_IN')
      CALL MEMMAN(KLCOCC_FIN,LEN_CMO,'ADDL  ',2,'COCC_FI')
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,
      CALL CSUB_FROM_C(WORK(KMOAOIN),WORK(KLCOCC_INI),NACOBS,ISCR_NTS,
     &                 1,ISCR,0)
      CALL CSUB_FROM_C(WORK(KLMO2),WORK(KLCOCC_FIN),NACOBS,ISCR_NTS,
     &                 1,ISCR,0)
C     CSUB_FROM_C(C,CSUB,LENSUBS,LENSUBTS,NSUBTP,ISUBTP,IONLY_DIM)
      WRITE(6,*) 
     &' Projecting final (MO2) on initial (MO1) active orbitals'
      CALL PROJ_ORBSPC_ON_ORBSPC(WORK(KLCOCC_INI),WORK(KLCOCC_FIN),
     &     NACOBS,NACOBS)
C     PROJ_ORBSPC_ON_ORBSPC(CMOAO1,CMOAO2,NMO1PSM,NMO2PSM)
*. Print summary
      CALL PRINT_MCSCF_CONV_SUMMARY(WORK(KL_SUMMARY),NOUTIT)
      WRITE(6,'(A,F20.12)') ' Final energy = ', EFINAL
      WRITE(6,'(A,F20.12)') ' Final norm of orbital gradient = ', 
     &                        E1NRM_ORB
*
C?    WRITE(6,*) ' E1NRM_ORB, ERROR_NORM_FINAL_CI = ',
C?   &             E1NRM_ORB, ERROR_NORM_FINAL_CI
C?    WRITE(6,*) ' Final energy = ', EFINAL

      CALL MEMMAN(IDUMMY, IDUMMY, 'FLUSM', IDUMMY,'AKMCSC') 
      CALL LUCIAQEXIT('AKMCS')
      RETURN
      END
      SUBROUTINE AKBKSIGMA2(CB,HCB,LUC,LUHC)
*. 
*. Obtain AKBK Hamiltonian times a vector on LUC
*
*  Version where storage of complete Q-vectors have been
*  eliminated
*
* A P and a Q space has been defined, so the eigenvalue problem reads:
*
* ( H(PP)    H(PQ) ) (CP)        (CP)
* (                ) (  )   = E  (  )
* ( H(QP)    H(QQ) ) (CQ)        (CQ)
* 
* Lowdin partition the matrix to obtain
*
* [(HPP - H(PQ) (H(QQ) - E )^(-1)H(QP)]CP = E CP
*
* HQQ is assumed so simple it may be diagonalized - here
* Diagonal or block diagonal over confs
*
*. Input and output vectors are in the P-space
*
*. It is assumed that ICNFBAT = 2 and NOCSF = 0- although pieces of code
*  easing the transition to in core are included
*
*. Jeppe Olsen, Aarhus, July 2013
*. Last revision; Aug. 2013; Jeppe Olsen; H0_CSF = 1 included
*
      INCLUDE 'implicit.inc'
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc-static.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'cshift.inc'
      INCLUDE 'akbk.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'cands.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cecore.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'ceblk.inc'
*
      CHARACTER*6 ITASK
*. Local scratch
      INTEGER ICONF_OCC(MXPORB)
cVOG  DIMENSION CB(*), HCB(*)
      integer CB, HCB
      LBLK = -1
      ITASK='SIGMA '
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',1,'AKBKSI')
      CALL LUCIAQENTER('AKBKSI')
*
      NTEST = 10
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Information from AKBKSIGMA2'
        WRITE(6,*) ' ==========================='
        WRITE(6,*)
        WRITE(6,*) ' IPSPC, IQSPC = ', IPSPC, IQSPC
        WRITE(6,*) ' NSDP , NSDQ  = ', NSDP , NSDQ
        WRITE(6,*) ' ISIMSYM      = ', ISIMSYM
        WRITE(6,*) ' IQINP        = ', IQINP
        WRITE(6,*) ' ICISTR       = ', ICISTR
        WRITE(6,*) ' NSBATCH      = ', NSBATCH
      END IF 
cVOG  IF(NTEST.GE.100) THEN
cVOG    WRITE(6,*) ' LUC, LUHC = ', LUC, LUHC
cVOG  END IF
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' Input vector to AKBKSIGMA2 '
        CALL ga_print(CB)
cVOG    CALL WRTVCD(CB,LUC,1,LBLK)
      END IF
*
* Allocate scratch files
*
*. For C and Sigma in det basis
cVOG  CALL FILEMAN_MINI(LU_CDET,'ASSIGN')
cVOG  CALL FILEMAN_MINI(LU_SDET,'ASSIGN')
      IF(ICISTR.EQ.1) THEN
       LU_CDET = 0
       LU_SDET = 0
      ELSE
       STOP ' Something is messy here...'
      END IF       
*
!      IF ( .not. ga_create( mt_dbl, NSDP, 1, 'LU_CDET', 0, 0, KLU_CDET))
!     $   CALL ERRQUIT('AKBK: ga_create LU_CDET failed', NSDP, GA_ERR)
!      IF ( .not. ga_create( mt_dbl, NSDQ, 1, 'LU_SDET', 0, 0, KLU_SDET))
!     $   CALL ERRQUIT('AKBK: ga_create LU_SDET failed', NSDQ, GA_ERR)
!      CALL ga_zero(KLU_CDET)
!      CALL ga_zero(KLU_SDET)
!      call ga_print(KLU_SDET)
*
*. Construct info about the block-structures in P and Q space
*
      ICSPC = IPSPC
      ISSPC = IQSPC
*
      ISSM = IREFSM
      ICSM = IREFSM
*
cVOG  ISFIRST = 1
      ISFIRST = 0
      ICFIRST = 1
*
      MAXK1_MX = 0
      LSCMAX_MX = 0
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
C Block for storing complete or partial CI-vector
      write(*,*) 'MXSOOB, MXSOOB_AS', MXSOOB, MXSOOB_AS
      IF(ISIMSYM.EQ.0) THEN
        LBLOCK = MXSOOB
      ELSE
        LBLOCK = MXSOOB_AS
      END IF
      ICOMP = 0
*. Allocate space for one occupation class
      IF(NOCSF.EQ.0.OR.ICNFBAT.EQ.-2) THEN
        LBLOCK  = MAX(NSD_FOR_OCCLS_MAX,MXSOOB)
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) ' TEST: LCSBLK, LBLOCK, MXSOOB  = ',
     &                     LCSBLK, LBLOCK, MXSOOB
      END IF

      NTESTL = 000
      IF(ISFIRST.EQ.1) THEN
        CALL Z_BLKFO_FOR_CISPACE(ISSPC,ISSM,LBLOCK,ICOMP,
     &       NTESTL,NSBLOCK,NSBATCH,
     &       int_mb(KSIOIO),int_mb(KSBLTP),NSOCCLS_ACT,
     &       int_mb(KSIOCCLS_ACT),
     &       int_mb(KSLBT),int_mb(KSLEBT),int_mb(KSLBLK),int_mb(KSI1BT),
     &       int_mb(KSIBT),
     &       int_mb(KSNOCCLS_BAT),int_mb(KSIBOCCLS_BAT),ILTEST)
        NSOCCLS = NSOCCLS_ACT
      END IF
      IF(ICFIRST.EQ.1) THEN
        CALL Z_BLKFO_FOR_CISPACE(ICSPC,ICSM,LBLOCK,ICOMP,
     &       NTESTL,NCBLOCK,NCBATCH,
     &       int_mb(KCIOIO),int_mb(KCBLTP),NCOCCLS_ACT,
     &       int_mb(KCIOCCLS_ACT),
     &       int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &       int_mb(KCIBT),
     &       int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
        NCOCCLS = NCOCCLS_ACT
      END IF
*
      NPBLK = NCBLOCK
      NPBATCH = NCBATCH
      IF(NTEST.GE.1000.and.ga_nodeid().eq.0)
     &WRITE(6,*) '  NCOCCLS, NPBLK = ', NCOCCLS, NPBLK
      NPSD = IELSUM(int_mb(KCLBLK),NPBLK)
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
* Obtain J and K matrices for later use
*
      CALL MEMMAN(KLJ,NTOOB**2,'ADDL  ',2,'IIJJ  ')
      CALL MEMMAN(KLK,NTOOB**2,'ADDL  ',2,'IJJI  ')
*
*. J(I,J) = (II!JJ), K(I,J) = (IJ!JI)
*
      CALL GTJK(dbl_mb(KLJ),dbl_mb(KLK),NTOOB,XDUM,IREOTS)
*
* Scratch space for diagonal

*. largest number of protype dets: that of maxop
cVOG  NPDT_MAX = NPCMCNF(MAXOP+1)
cVOG  NPCS_MAX = NPCSCNF(MAXOP+1)
cVOG  NPCMCNF and NPCSCNF are initiallized only for NOCSF = 0
cVOG  So, be careful!
*
      NEL = NELFTP(IATP) + NELFTP(2)
*
cVOG  IF(NTEST.GE.100.and.ga_nodeid().eq.0)
cVOG &WRITE(6,*) ' NPDT_NAX, NPCS_MAX = ', NPDT_MAX, NPCS_MAX
cVOG  LISCR = 2*NPDT_MAX*NEL + NPDT_MAX + 6*NACOB
cVOG  LRSCR = 2*NPDT_MAX**2
      liscr = 2*nocob + 6*nacob !Achtung!
      lrscr = 2*nocob**2        !Achtung!
      CALL MEMMAN(KLISCR,LISCR,'ADDL  ',1,'IS_CHC')
      CALL MEMMAN(KLRSCR,LRSCR,'ADDL  ',2,'RS_CHC')
cVOG  LCSFHCSF = NPCS_MAX**2
cVOG  CALL MEMMAN(KLCSFHCSF,LCSFHCSF,'ADDL  ',2,'CSHCS ')
*. For holding info on largest batch = occupation class
      NBLK_BAT_MAX = IMNMX(int_mb(KSLBT),NSBATCH,2)
      CALL MEMMAN(KEIBT,8*NBLK_BAT_MAX,'ADDL  ',2,'EBT   ')
      CALL MEMMAN(KELBT,1,'ADDL ',2,'ELBT  ')
      CALL MEMMAN(KEIOCCLS_ACT,1,'ADDL  ',2,'EIOCCL')
C     KEIBT,NEBLK,NEOCCLS,KEIOCCLS_ACT,KELBT
*
* 1: Transfer C from CSF to SD/CM basis 
*
cVOG  IF(ICNFBAT.EQ.1) THEN
*. In core
cVOG  CALL CSDTVCM(CB,WORK(KCOMVEC1_SD),WORK(KCOMVEC2_SD),
cVOG &               1,0,ICSM,ICSPC,2)
cVOG  ELSE
*. Not in core- write determinant expansion on LU_CDET
cVOG  CALL CSDTVCMN(CB,HCB,dbl_mb(KVEC3),
cVOG &      1,0,ICSM,ICSPC,2,2,LU_CDET,LUC,NCOCCLS_ACT,
cVOG &      int_mb(KCIOCCLS_ACT),int_mb(KCIBT),int_mb(KCLBT))
cVOG  END IF
!      call ga_copy(CB,KLU_CDET)
!      call ga_print(KLU_CDET)
* 2: HPP CP (in and out in SD/CM), save on LU_SDET
      ICSPC = IPSPC
      ISSPC = IQSPC
      IPRCIX = 0
cVOG  CALL RASSG3(CB,HCB,NCBATCH,int_mb(KCLBT),int_mb(KCLEBT),
cVOG &     int_mb(KCI1BT),int_mb(KCIBT),LU_CDET,LU_SDET,XDUM,XDUM,ECORE,
cVOG &     ITASK)
!~~~~~~~~~~~~~~~~~~~~~~COPIED FROM MV7~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!          CALL RASSG3(dbl_mb(KVEC1P),dbl_mb(KVEC2P),NSBATCH,
!     &         int_mb(KSLBT),int_mb(KSLEBT),
!     &         int_mb(KSI1BT),int_mb(KSIBT),LLUC,LLUHC,CB,HCB,ECORE,
!     &         CMV7TASK)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      CALL WRTMAT(dbl_mb(KVEC1P),1,NSDP, 1, 1)
      CALL WRTMAT(dbl_mb(KVEC2P),1,NSDP, 1, 1)
      CALL ga_print(CB)
      CALL ga_print(HCB)
      CALL RASSG3(dbl_mb(KVEC1P),dbl_mb(KVEC2P),NCBATCH,
     &     int_mb(KCLBT),int_mb(KCLEBT),
     &     int_mb(KCI1BT),int_mb(KCIBT),0,0,CB,HCB,ECORE,ITASK)
**** I AM AT THIS POINT!
      write(*,*) 'Out from the RASSG3(1) !!!'
      CALL WRTMAT(dbl_mb(KVEC1P),1,NSDP, 1, 1)
      CALL WRTMAT(dbl_mb(KVEC2P),1,NSDP, 1, 1)
      CALL ga_print(CB)
      CALL ga_print(HCB)
!      CALL ga_get(HCB, 1, NSDP, 1, 1, dbl_mb(LU_SDET),1)
!      CALL ga_copy(HCB,KLU_CDET)
!      CALL ga_zero(HCB)
!      CALL ga_print(KLU_CDET)
!      CALL ga_print(KLU_SDET)
*
* and next the Q-contributions
*
*. Notice: at the moment, a batch is by construction an occupation class
*.(see PART_CIV_OCC  where I_ONE_OCCLASS_PER_BATCH = 1)

      IBSBLOCK_BAT = 1
      DO IQBAT = 1, NSBATCH
C                  IFRMR(WORK,IROFF,IELMNT)
        IQ_OCCLS = int_mb(KSIOCCLS_ACT + IQBAT - 1)
        NSBLOCK_BAT = int_mb(KSLBT+IQ_OCCLS-1)
*. Is this occupation in P-space
C IS_IVEC_IN_LIST(IVEC,NELMNT,LIST,LLIST)
        IQINP = IS_IVEC_IN_LIST(IQ_OCCLS,1,int_mb(KCIOCCLS_ACT),
     &          NCOCCLS_ACT)
        IF(NTEST.GE.200) THEN
           WRITE(6,*)
           WRITE(6,*) ' ======================================'
           WRITE(6,*) ' Info for next Q-batch: '
           WRITE(6,*) ' IQBAT, IQ_OCCLS, IQINP = ',
     &                                IQBAT, IQ_OCCLS, IQINP
           WRITE(6,*) ' NSBLOCK_BAT = ', NSBLOCK_BAT
           WRITE(6,*) ' ======================================'
           WRITE(6,*)
        END IF
        IF(IQINP.EQ.0) THEN
*
* ===============================
* Calculate HQP CP for Q in batch
* ===============================
*
C         SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &           LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
C    &           ECORE,ITASK)
cVOG      CALL SBLOCK(NSBLOCK_BAT,int_mb(KSIBT),IBSBLOCK_BAT,CB,HCB,
cVOG &         LU_CDET,0,0,0,0,0,CV,ECORE,ITASK)
          CALL SBLOCK(NSBLOCK_BAT,int_mb(KSIBT),IBSBLOCK_BAT,
     &         dbl_mb(KVEC1P),dbl_mb(KVEC2P),
     &         -1,0,0,0,0,0,CB,ECORE,ITASK)
          write(*,*) 'Out from the SBLOCK(2)'
      CALL ga_print(CB)
      CALL ga_print(HCB)
      CALL WRTMAT(dbl_mb(KVEC1P), 1,NSDP, 1, 1)
      CALL WRTMAT(dbl_mb(KVEC2P), 1,NSDP, 1, 1)
!      CALL WRTMAT(dbl_mb(LU_CDET),1,NSDP, 1, 1)
!      CALL WRTMAT(dbl_mb(LU_SDET),1,NSDP, 1, 1)
*. The batch of HQP CP is returned in HCB
*
*
*. Transform HQP CP to CSF basis and save in CB
* 
*
*. Conf info
cVOG      CALL GEN_CNF_INFO_FOR_OCCLS(IQ_OCCLS,1,IREFSM)
cVOG      NCSF_OCCLS = IELSUM(NCS_FOR_OC_OP_ACT,MAXOP+1)
cVOG      NCM_OCCLS = IELSUM(NCM_FOR_OC_OP_ACT,MAXOP+1)
cVOG      NCM_OCCLS_AB =
cVOG &    int_mb(KNCMAB_FOR_OCCLS+(IQ_OCCLS-1)*NSMST+IREFSM-1)
*. HQP CP batch in CB
cVOG      NTESTL = 0
cVOG      CALL CSDTVCB(CB,HCB,2,dbl_mb(KDTOC),
cVOG &        WORK(KSDREO_I(IREFSM)),NCM_OCCLS,NCM_OCCLS_AB,
cVOG &        NCSF_OCCLS,NCN_FOR_OC_OP_ACT,
cVOG &        0,1,NTESTL)
*
* ===============================
*. (HQQ - E)**-1 HQP for Q batch
* ===============================
*
!      CALL ga_add_constant(KLU_Q3,SHIFT)
!      CALL ga_elem_divide(KLU_Q2,KLU_Q3,KLU_Q2)
c VOG VOG VOG PROBLEM, this part is written for CSFs!!! VOG VOG VOG
          CALL LUCIAQENTER('AKBKH0')
          IF(IH0_CSF.EQ.1) THEN
            KEBLTP = KCBLTP
            NEBLK = NSBLOCK_BAT
C                ICOPVE2(IIN,IOFF,NDIM,IOUT)
            CALL ICOPVE2(int_mb(KSIBT),(IBSBLOCK_BAT-1)*8+1,
     &           8*NSBLOCK_BAT,int_mb(KEIBT))
            CALL ICOPVE2(int_mb(KSLBT),IQBAT,1,dbl_mb(KELBT))
            CALL ICOPVE2(int_mb(KSIOCCLS_ACT),IQBAT,1,
     &                   int_mb(KEIOCCLS_ACT))
            ICISTR_SAVE = ICISTR
            ICNFBAT_SAVE = ICNFBAT
            ICISTR = 1
            ICNFBAT = 1
!            CALL CSDIAGM(IREFSM,3,-1,HCB,dbl_mb(KVEC3))
C                CSDIAGM(ISM,ICS,LUDIA_A,VECUT,VECSCR)
            ICISTR =ICISTR_SAVE
            ICNFBAT = ICNFBAT_SAVE
C     COMMON/EBLK/KEBLTP,KEIBT,NEBLK,NEOCCLS,KEIOCCLS_ACT,KELBT
          ELSE IF (IH0_CSF.EQ.2) THEN
            IONLY_DIAG = 1
            ISYMG = 1
            IB_OCC = 1
            IB_H0 = 1
            IDO_REO = 1
            IB_BLK = 1
            NCSF_TOT = 0
*. Construct diagonal in HCB
            DO IOPEN = 0, MAXOP
              ITYP =  IOPEN + 1
              IORB = (NEL+IOPEN)/2
              NCSF = NPCSCNF(IOPEN+1)
COLD          NCONF = NCONF_FOR_OPEN(IOPEN+1)
              NCONF = NCONF_PER_OPEN(IOPEN+1,IREFSM)
              DO  ICONF = 1, NCONF
*. Occupation
               CALL ICOPVE2(
     &              int_mb(KICONF_OCC(IREFSM)),IB_OCC,IORB,ICONF_OCC)
               IF(NTEST.GE.10000) THEN
                 WRITE(6,*) ' Occupation of next conf '
                 CALL IWRTMA(ICONF_OCC,1,IORB,1,IORB)
               END IF
*. Diagonal elements (delivered as diagonal elements of complete matrix even if ISYMG=1)
               CALL CNHCN_CSF_BLK(ICONF_OCC,IOPEN,
     &         ICONF_OCC,IOPEN,dbl_mb(KLCSFHCSF),NINOB,int_mb(KDFTP),
     &         dbl_mb(KDTOC),2,dbl_mb(KLISCR),dbl_mb(KLRSCR),ECORE,
     &         IONLY_DIAG,ISYMG,dbl_mb(KLJ),dbl_mb(KLK))
*. Extract diagonal elements 
!               CALL COPDIA(dbl_mb(KLCSFHCSF),HCB(IB_H0),NCSF,0)
C                   COPDIA(A,VEC,NDIM,IPACK)
               IF(NTEST.GE.10000) THEN
                 WRITE(6,*) 
     &           ' Added elements, IB_H0, NCSF = ', IB_H0, NCSF
cVOG              CALL WRTMAT(HCB(IB_H0),1,NSCF,1,NCSF)
               END IF
               IB_H0 = IB_H0 + NCSF
               IB_OCC = IB_OCC + IORB
               NCSF_TOT = NCSF_TOT + NCSF
              END DO ! Loop over configurations for given IOPEN
            END DO ! Loop over IOPEN
          END IF ! IH0_CSF switch
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Diagonal elements: '
!            CALL ga_print(HCB)
cVOG        CALL WRTMAT(HCB,1,NCSF_OCCLS,1,NCSF_OCCLS)
          END IF
C (HQQ-E)**-1 HQP CP in CB
C              DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
C?        WRITE(6,*)  ' Shift = ', Shift
!          CALL DIAVC2(CB,CB,HCB,SHIFT,NCSF_OCCLS)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' (Hqq-E)**-1 HQP CP '
            CALL ga_print(CB)
cVOG        CALL WRTMAT_EP(CB,1,NCSF_OCCLS,1,NCSF_OCCLS)
          END IF
*
* Transform (HQQ-E)**-1 HQP CP to SD basis and save in in KVEC4
*
          CALL LUCIAQEXIT('AKBKH0')
cVOG      NTESTL = 000
cVOG      CALL CSDTVCB(CB,WORK(KVEC4),1,dbl_mb(KDTOC),
cVOG &        WORK(KSDREO_I(IREFSM)),NCM_OCCLS,NCM_OCCLS_AB,
cVOG &        NCSF_OCCLS,NCN_FOR_OC_OP_ACT,
cVOG &        0,1,NTESTL)
* ===================================
*. HPQ (HQQ - E)**-1 HQP for Q batch
* ===================================
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*)
            WRITE(6,*) ' =========================================='
            WRITE(6,*) '  calc of HPQ (HQQ - E)**-1 HQP for Q batch'
            WRITE(6,*) ' =========================================='
            WRITE(6,*)
          END IF
* We have RHS in core, fool program to think we are doing in core
          ICISTR_SAVE = ICISTR
          ICISTR = 1
*
          ICSPC = IQSPC
          ISSPC = IPSPC
*
C              SWAP_CS_TTSSARRAYS
          CALL SWAP_CS_TTSSARRAYS
          IP_B_BLK = 1
          IP_B_CM  = 1
          DO IP_BATCH = 1, NPBATCH
           NP_BLK_BAT = int_mb(KSLBT+IP_BATCH-1)
           NP_CM_BAT  = int_mb(KSLEBT+IP_BATCH-1)
           IF(NTEST.GE.1000)
     &     WRITE(6,*) ' IP_BATCH, NP_BLK_BAT, NP_CM_BAT = ',
     &                  IP_BATCH, NP_BLK_BAT, NP_CM_BAT
!~~~~~~~~~~~~~~~~~~~~~~COPIED FROM RASSG3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!        CALL SBLOCK(LBATS(JBATS),IBATS(1,I1BATS(JBATS)),1,
!     &       CB,SB,LUC,0,0,0,0,0,CV,ECORE,ITASK,LEN)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
           ilen = NSD_FOR_OCCLS_MAX
           CALL SBLOCK(NP_BLK_BAT,int_mb(KSIBT),IP_B_BLK,
     &          dbl_mb(KVEC1P),dbl_mb(KVEC2P),
     &          -1,0,0,1,IQBAT,IQBAT,dbl_mb(KVEC4),ECORE,ITASK,ilen)
           CALL ga_print(KVEC4)
cVOG       CALL SBLOCK(NP_BLK_BAT,WORK(KSIBT),IP_B_BLK,CB,
cVOG &          HCB(IP_B_CM),
cVOG &          -1,0,0,1,IQBAT,IQBAT,WORK(KVEC4),ECORE,ITASK)
C               SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &          LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,CV,
C    &          ECORE,ITASK)
           IP_B_BLK = IP_B_BLK + NP_BLK_BAT
           IP_B_CM  = IP_B_CM  + NP_CM_BAT
          END DO
*
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) 
     &      ' Contribution to HPQ (HQQ - E)**-1 HQP for Q batch'
            CALL ga_print(HCB)
cVOG       CALL WRTMAT_EP(HCB,1,NPSD,1,NPSD)
          END IF
*
*. And restore standard order
          ICSPC = IPSPC 
          ISSPC = IQSPC
          CALL SWAP_CS_TTSSARRAYS
          ICISTR = ICISTR_SAVE
*. We have now another contribtion to HPP(eff) CP in CB, update sum
*. 
*. Update current sum
cVOG      CALL REWINO(LU_SDET)
C              FRMDSCN(VEC,NREC,LBLK,LU)
!          CALL FRMDSCN(CB,NPBLK,LBLK,LU_SDET)
          ONE = 1.0D0
          ONEM =  -1.0D0
!          CALL VECSUM(CB,CB,HCB,ONE,ONEM,NPSD)
cVOG      CALL REWINO(LU_SDET)
!          CALL TODSCN2(CB,NPBLK,int_mb(KCLBLK),LBLK,LU_SDET,1) 
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' Updated LU_SDET '
cVOG        CALL WRTVCD_EP(CB,LU_SDET,1,LBLK)
          END IF
        END IF ! Q-class differed from P-class
        IBSBLOCK_BAT = IBSBLOCK_BAT + NSBLOCK_BAT  
      END DO! Loop over occupation classes of Q
*
      IF(NTEST.GE.1000.and.NOCSF.eq.0) THEN
        WRITE(6,*) ' Final HPP(eff) CP in SD basis '
cVOG    CALL WRTMAT_EP(CB,1,NPSD,1,NPSD)
      END IF
*
*. We have HPP(eff) CP in SD-basis in LU_SDET, transform to CSF basis
*. and save in 
*
!      IF(ICNFBAT.EQ.1) THEN
*. In core
!        CALL CSDTVCM(CB,WORK(KCOMVEC1_SD),WORK(KCOMVEC2_SD),
!     &               2,0,ICSM,ICSPC,2)
!      ELSE
*. Not in core- write determinant expansion on LU_CDET
!       CALL CSDTVCMN(HCB,CB,dbl_mb(KVEC3),
!     &      2,0,ICSM,ICSPC,2,2,LU_SDET,LUHC,NCOCCLS_ACT,
!     &      int_mb(KCIOCCLS_ACT),int_mb(KCIBT),int_mb(KCLBT))
!      END IF
      IF(NTEST.GE.1000) WRITE(6,*) ' Home from CSDTVCMN '
*
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) '  Final HPP(eff) CP in CSF basis '
       WRITE(6,*) ' ================================'
       CALL ga_print(HCB)
cVOG   CALL WRTVCD(HCB, LUHC,1,LBLK)
C           WRTVCD(SEGMNT,LU,IREW,LBLK)
      END IF
*
cVOG  CALL FILEMAN_MINI(LU_CDET,'FREE  ')
cVOG  CALL FILEMAN_MINI(LU_SDET,'FREE  ')
*
!      if (.not. ga_destroy(KLU_CDET)) call
!     &     errquit('AKBK: ga_destroy KLU_CDET failed', 0 , GA_ERR)
!      if (.not. ga_destroy(KLU_SDET)) call
!     &     errquit('AKBK: ga_destroy KLU_SDET failed', 0 , GA_ERR)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',1,'AKBKSI')
      CALL LUCIAQEXIT('AKBKSI')
*
      RETURN
      END
      SUBROUTINE SWAP_CS_TTSSARRAYS
*
* The programmer has in his or her infinite wisdom decided to
* swap the C and S spaces. Swap the corresponding pointers to the TTSS arrays
*
*. Jeppe Olsen, July 25 2013, 
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cands.inc'
*
C KCLBT,KCLEBT,KCI1BT,KCIBT,KC2B,KCIOIO,KCBLTP,KCLBLK,
      CALL SWAP_TWO_INT(KCLBT,KSLBT)
      CALL SWAP_TWO_INT(KCLEBT,KSLEBT)
      CALL SWAP_TWO_INT(KCI1BT,KSI1BT)
      CALL SWAP_TWO_INT(KCIBT,KSIBT)
      CALL SWAP_TWO_INT(KC2B,KS2B)
      CALL SWAP_TWO_INT(KCIOIO,KSIOIO)
      CALL SWAP_TWO_INT(KCBLTP,KSBLTP)
      CALL SWAP_TWO_INT(KCLBLK,KSLBLK)
      CALL SWAP_TWO_INT(NCBATCH,NSBATCH)
*
      RETURN
      END
      SUBROUTINE SWAP_TWO_INT(I1,I2)
*
* Swap integers I1 and I2
*
*. Jeppe Olsen ( to avoid mistakes in swapping a number of integers...)
*
      INCLUDE 'implicit.inc'
*
      IX = I1
      I1 = I2
      I2 = IX
*
      RETURN
      END
*======================================================================*
      subroutine reduce_gas_spaces(trgas,densi1,densi1_bu,
     &                                   densi2,densi2_bu,only_dens)
*======================================================================*
*
*. Subroutine that defines which GAS spaces are part of the MCSCF.
*. It is used (1) for the MO transformation in MOROT_GS
*.        and (2) for the reordering of DENSI2 after MOROT_GS
*
*. The algorithm checks which GAS spaces should be transformed, and
*. changes MO and occupation pointers for the GAS spaces that will not
*. be transformed.
*. These GAS spaces will become either Inactive or Secondary
*
*  Konstantinos D. Vogiatzis, June 2015, Berkeley
*
*======================================================================*
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'glbbas.inc'
      include 'clunit.inc'
      include 'cgas.inc' 
      include 'orbinp.inc'
      include 'strinp.inc'
      include 'strbas.inc'
      include 'gasstr.inc'
      include 'lucinp.inc'
#include "global.fh"
#include "mafdecls.fh"
*
*. Input/Output
      integer   only_dens
      dimension densi1(*), densi1_bu(*)
      dimension densi2(*), densi2_bu(*)
*. Output
      integer trgas(ngas+2)
*. Scratch
      integer          gasocc(ngas,2), sumocc(2)
      integer          ioffden(3)
      integer          ntotmos(ngas)
      integer          maxocc(ngas)
      integer          ntest, ntrgas, nchk
      character(len=3) yesno
*
*. Initializations
*
      ! Get dimensions for RHO2
      call getnddens2(ndactorb,nddensi2,only_dens)
*
      NTEST = 00
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
       write(luout,*) ' Only_dens? ', only_dens
       write(luout,*) ' Initial Density '
       call  wrtmat(densi1,nacob,nacob,nacob,nacob)
       write(luout,*)
       if (only_dens.eq.1) call prsym(densi2,ndactorb*ndactorb)
      end if
*
*. Keep backups for pointers and vectors that you will mofify
*
      ntrgas   = ngas
      ngas_bu  = ngas
*
      nocob_bu = nocob
      nacob_bu = nacob
      ninob_bu = ninob
      nscob_bu = nscob
*
      ngas_gnsym_bu = 0
      nocobs_bu = 0
      nacobs_bu = 0
      ninobs_bu = 0
      nscobs_bu = 0
      do igensm = 1, ngensmob
       do igas= 1, ngas+2
        ngas_gnsym_bu(igensm,igas) = ngas_gnsym(igensm,igas-1)
       end do
       nocobs_bu(igensm) = nocobs(igensm)
       nacobs_bu(igensm) = nacobs(igensm)
       ninobs_bu(igensm) = ninobs(igensm)
       nscobs_bu(igensm) = nscobs(igensm)
      end do
*
      nobpts_gn_bu = 0
      do igensm = 1, ngensmob
       nobpts_gn_bu(0,igensm) = nobpts_gn(0,igensm)
       nobpts_gn_bu(ngas+1,igensm) = nobpts_gn(ngas+1,igensm)
       do igas = 1, ngas
        nobpts_gn_bu(igas,igensm) = nobpts_gn(igas,igensm)
       end do
      end do
*
      ! Symmetry => type reordering array
      do i = 1, ntoob
       ireost_bu(i) = ireost(i)
      end do
*
      ! For DENSI2
      if (only_dens.eq.1) then
       do igas= 1, ngas+2
        nobpt_bu(igas)   = nobpt(igas)
        ihpvgas_bu(igas) = ihpvgas(igas)
       end do
      end if
*
*. Store old density
*
      do i = 1, nacob**2
       densi1_bu(i) = densi1(i)
      end do
*
      if (only_dens.eq.1) then
       do i = 1, nddensi2
        densi2_bu(i) = densi2(i) !TODO PROSOXI!!!
       end do 
      end if
*
*. Calculate min and max occupation numbers of each GAS space from the
*. accumulated occupations of CI space 1
*
      icispc = 1
*
      sumocc(1) = 0
      sumocc(2) = 0
      gasocc(1,1) = igsoccx(1,1,icispc)
      gasocc(1,2) = igsoccx(1,2,icispc)
      if (ngas.gt.1) then
       do igas = 1, ngas-1
        sumocc(1) = sumocc(1) + gasocc(igas,1)
        sumocc(2) = sumocc(2) + gasocc(igas,2)
        gasocc(igas+1,1) = igsoccx(igas+1,1,icispc) - sumocc(2)
        gasocc(igas+1,2) = igsoccx(igas+1,2,icispc) - sumocc(1)
       end do
      end if
*
*. Identify non-filled GAS spaces
*
      ntotmos = 0
      do igas = 1, ngas
       do igensm = 1, ngensmob
        ntotmos(igas) = ntotmos(igas) + ngssh(igensm,igas)
       end do
       maxocc(igas) = ntotmos(igas) * 2
      end do
*
      trgas = 1
      ioffden(1) = 0
      ioffden(3) = 0
      ioffden(2) = 0
*
      do igas = 1, ngas
*
       !
       ! CASE 1: Empty GAS space belonging to Q-space
       !         Construct offsets for density resize ioffden(3)
       !
       if (gasocc(igas,2).eq.0) then
        trgas(igas) = 0
        ntrgas = ntrgas - 1
        ! Add MOs to secondary space
        do igensm = 1, ngensmob
         nocobs(igensm) = nocobs(igensm) - ngas_gnsym(igensm,igas)
         nacobs(igensm) = nacobs(igensm) - ngas_gnsym(igensm,igas)
         nscobs(igensm) = nscobs(igensm) + ngas_gnsym(igensm,igas)
         nocob = nocob - ngas_gnsym(igensm,igas)
         nacob = nacob - ngas_gnsym(igensm,igas)
         nscob = nscob + ngas_gnsym(igensm,igas)
         ioffden(3) = ioffden(3) + ngas_gnsym(igensm,igas)
*
         ngas_gnsym(igensm,ngas+1) = ngas_gnsym(igensm,ngas+1) +
     &                               ngas_gnsym(igensm,igas)
         ngas_gnsym(igensm,igas) = 0
        end do
*
        do igensm = 1, ngensmob
         nobpts_gn(ngas+1,igensm) = nobpts_gn(igas,igensm) +
     &                              nobpts_gn(ngas+1,igensm)
         nobpts_gn(igas,igensm) = 0
        end do
*
        if (only_dens.eq.0) then
        !Change MO pointer of itpfso to ngas+1 (Secondary)
         do i = 1, ntoob
          if (itpfso(i).eq.igas) itpfso(i) = ngas + 1
         end do
        end if
*
       end if
*
       !
       ! CASE 2: Fully occupied GAS space belonging to Q-space
       !         Construct offsets for density resize ioffden(1)
       !
       if (gasocc(igas,1).eq.maxocc(igas)) then
        trgas(igas) = 0
        ntrgas = ntrgas - 1
       if (ga_nodeid().eq.0)
     &   write(*,*) 'VOG, get rid of occup GAS of Q-Space', ntrgas
        ! Add MOs to the inactive space 
        do igensm = 1, ngensmob
         nacobs(igensm) = nacobs(igensm) - ngas_gnsym(igensm,igas)
         ninobs(igensm) = ninobs(igensm) + ngas_gnsym(igensm,igas)
         ninob  = ninob  + ngas_gnsym(igensm,igas)
         nacob  = nacob  - ngas_gnsym(igensm,igas)
         ioffden(1) = ioffden(1) + ngas_gnsym(igensm,igas)
*
        !TODO These should be checked for more GAS spaces in Q
         ngas_gnsym(igensm,0) = ngas_gnsym(igensm,igas) +
     &                          ngas_gnsym(igensm,0)
         ngas_gnsym(igensm,igas) = 0
        end do
*
        do igensm = 1, ngensmob
         nobpts_gn(0,igensm) = nobpts_gn(igas,igensm) +
     &                         nobpts_gn(0,igensm)
         nobpts_gn(igas,igensm) = 0
        end do
*
        if (only_dens.eq.0) then
        !Change MO pointer of itpfso to 0 (Inactive)
         do i = 1, ntoob
          if (itpfso(i).eq.igas) itpfso(i) = 0
         end do
        end if
*
       end if
*
       !
       ! CASE 3: GAS space that is treated as part of the active space
       !         for MCSCF and P-Space for AKBK
       !         The number of MOs are stored in ioffden(2)
       !
       if (gasocc(igas,2).ne.0.and.gasocc(igas,1).ne.maxocc(igas)) then
        do igensm = 1, ngensmob
         ioffden(2) = ioffden(2) + ngas_gnsym(igensm,igas)
        end do
       end if
      end do
*
*. Rearrange GAS spaces
*
      itrgas = 0
      do igas = 1, ngas
       !Check if the GAS is empty
       mosum = 0
       do igensm = 1, ngensmob
        mosum = mosum + ngas_gnsym(igensm,igas)
       end do
       if (mosum.ne.0) then
        itrgas = itrgas + 1
        do igensm = 1, ngensmob
         ngas_gnsym(igensm,itrgas) = ngas_gnsym(igensm,igas)
         nobpts_gn(itrgas,igensm)  = nobpts_gn(igas,igensm) 
        end do
       end if
      end do
*
      if (itrgas.ne.ntrgas.and.ga_nodeid().eq.0) then
       write(*,'(a,i3)') 'Inconsistent number of transformed GAS',
     &                    itrgas
       CALL ERRQUIT('reduce_gas_spaces: wrong number tranformed GAS',
     &              itrgas, 667)
      end if
*
*. Put the correct pointers for the secondary space
*
      do igensm = 1, ngensmob
       ngas_gnsym(igensm,itrgas+1) = ngas_gnsym(igensm,ngas+1)
       nobpts_gn(itrgas+1,igensm)  = nobpts_gn(ngas+1,igensm)
      end do
*
*. New number of GAS in reduced space
*
      ngas = ntrgas
*
*. Change MO pointer of itpfso for Active and Secondary 
*
      if (only_dens.eq.0) then
       ! GAS spaces that will be transformed (Active)
       itrgas = 0
       do igas = 1, ngas_bu
        if (trgas(igas).eq.1) then
         itrgas = itrgas + 1
         do i = 1, ntoob
          if (itpfso(i).eq.igas) itpfso(i) = itrgas
         end do
        end if
       end do
*
       ! Secondary space
       do i = 1, ntoob
        if (itpfso(i).eq.(ngas_bu + 1)) itpfso(i) = ngas + 1
       end do
*
      end if
*
*. Nulify array IOBPTS_SEL, it will be updated in DENSI2 for the proper
*. GAS spaces from IB_FOR_SEL_ORBSPC
*
      do igas = 1, ngas_bu
       do igensm = 1, ngensmob
        IOBPTS_SEL(igas,igensm) = 0
       end do
      end do
*
*. Reorder the "Symmetry => type reordering array"
*
      iorb = 0
      do igensm = 1, ngensmob
       do igas = 0, ngas+1
        ! Offset
        ioff = 0
        do ioffgas = 0, igas-1
         do ioffgensm = 1, ngensmob
          ioff = ioff + ngas_gnsym(ioffgensm,ioffgas)
         end do
        end do
        do ioffgensm = 1, igensm - 1
         ioff = ioff + ngas_gnsym(ioffgensm,igas)
        end do
        ! Add values
        do i = 1, ngas_gnsym(igensm,igas)
         iorb = iorb + 1
         ientry = i + ioff
         ireost(iorb) = ientry
        end do
       end do
      end do
*
*. Resize the density
*
      ioffden(2) = ioffden(2) + ioffden(1)
      ioffden(3) = ioffden(3) + ioffden(2)
      if (ioffden(3).ne.nacob_bu.and.ga_nodeid().eq.0) then
       write(*,'(a)') 'Inconsistent number of MOs and ioffden(3)'
       CALL ERRQUIT('reduce_gas_spaces: wrong number of MOs',ioffden(3),
     &              667)
      end if
*
      call resize_dens(densi1,densi1_bu,densi2,densi2_bu,ioffden,
     &                 only_dens)
*
*. Modify NOBPT and IHPVGAS, important for DENSI2
*
      if (only_dens.eq.1) then
       nobpt   = 0
       ihpvgas = 0
       itrgas  = 0
       do igas = 1, ngas_bu
        if (trgas(igas).eq.1) then
         itrgas = itrgas + 1
         nobpt(itrgas)   = nobpt_bu(igas)
         ihpvgas(itrgas) = ihpvgas_bu(igas)
        end if
       end do
      end if
*
*. Print some output
*
      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
       WRITE(LUOUT,'(A)')
       WRITE(LUOUT,'(A)')
       WRITE(LUOUT,'(A)') '         Modified GAS Spaces '
       WRITE(LUOUT,'(A)') '         =================== '
       WRITE(LUOUT,'(A)')
       WRITE(LUOUT,'(A)') '         Min. occ    Max. occ    Transform? '
       WRITE(LUOUT,'(A)') '         ================================== '
       DO igas = 1, ngas_bu
        yesno = 'YES'
        if (trgas(igas).eq.0) yesno = ' NO'
        WRITE(LUOUT,'(A,I2,3X,I3,9X,I3,10X,3A)')
     &   '   GAS',IGAS,gasocc(IGAS,1),gasocc(IGAS,2),yesno
       END DO
*
       WRITE(LUOUT,'(A)')
       WRITE(LUOUT,'(A)') '         Modified Orbital subspaces:'
       WRITE(LUOUT,'(A)') '         ========================== '
       WRITE(LUOUT,*)
       WRITE(LUOUT,'(1H ,A,10I4,A)')
     &  '        Irrep ',(I,I = 1,ngensmob)
       WRITE(LUOUT,'(1H ,A,2X,10A,A)')
     &  '        ===== ',('====',I = 1,ngensmob)
       WRITE(LUOUT,'(A,10I4)')
     &  '  Inactive     ',
     &   (ngas_gnsym(igensm,0),igensm = 1, ngensmob)
C      WRITE(LUOUT,'(A,10I4)')
C    &  '  Inactive     ',
C    &   (nobpts_gn(0,igensm),igensm = 1, ngensmob)
*
       DO IGAS = 1, ntrgas
        WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
     &  '  GAS',IGAS,'        ',
     &  (ngas_gnsym(igensm,IGAS), igensm = 1, ngensmob)
C       WRITE(LUOUT,'(A,I2,A,10I4,6X,2I6)')
C    &  '  GAS',IGAS,'        ',
C    &   (nobpts_gn(igas,igensm),igensm = 1, ngensmob)
       END DO
       WRITE(LUOUT,'(A,10I4)')
     &  '  Secondary    ',
     &   (ngas_gnsym(igensm,ntrgas+1),igensm = 1, ngensmob)
C      WRITE(LUOUT,'(A,10I4)')
C    &  '  Secondary    ',
C    &   (nobpts_gn(ntrgas+1,igensm),igensm = 1, ngensmob)
       WRITE(6,*)
       if(ntest.ge.100) then
*
        if (only_dens.eq.0) then
         write(luout,*) ' Modified ITPFSO string: '
         call iwrtma3(itpfso,1,ntoob,1,1)
         write(luout,*)
*
cVOG    else
         write(luout,*) ' Modified Density '
         call wrtmat(densi1,nacob,nacob,nacob,nacob)
         write(luout,*)
        else
         if (only_dens.eq.1) then
          ! Get reduced dimensions for RHO2
          call getnddens2(ndactorb,nddensi2,only_dens)
          call prsym(densi2,ndactorb*ndactorb)
         end if
         write(luout,*)
        end if
       write(*,'(a,10i4)') 'INACTIVE ',
     &          (ninobs(isym), isym = 1, ngensmob), ninob
       write(*,'(a,10i4)') 'ACTIVE   ',
     &          (nacobs(isym), isym = 1, ngensmob), nacob
       write(*,'(a,10i4)') 'OCCUPIED ',
     &          (nocobs(isym), isym = 1, ngensmob), nocob
       write(*,'(a,10i4)') 'SECONDARY',
     &          (nscobs(isym), isym = 1, ngensmob), nscob
       end if
      END IF
*
*. Final check
*
      do igensm = i, ngensmob
       nchk = 0
       do igas = 1, ntrgas
        nchk = nchk + ngas_gnsym(igensm,igas)
       end do
       nchk = nchk + ngas_gnsym(igensm,0) + ngas_gnsym(igensm,ntrgas+1)
       if (ntoobs(igensm).ne.nchk.and.ga_nodeid().eq.0) then
        write(*,'(a,i3)') 'Inconsistent number of MOs in IRREP', igensm
        CALL ERRQUIT('reduce_gas_spaces: wrong number of MOs', nchk,667)
       end if
      end do
*
      return
      end 
*======================================================================*
*======================================================================*
      subroutine resize_dens(densi1,densi1_bu,densi2,densi2_bu,ioffden,
     &                       only_dens)
*======================================================================*
*
*  Konstantinos D. Vogiatzis, June 2015, Berkeley
*
*======================================================================*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'glbbas.inc'
      include 'cgas.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'multd2h.inc'
*
*. Input
      integer   only_dens
      integer   nddensi2
      integer   ioffden(3)
      dimension densi1_bu(*)
      dimension densi2_bu(*)
*
*. Output
      dimension densi1(*)
      dimension densi2(*)
*
*. Local
      integer   ipqden(3)
*
*.
*
      ntest = 10 !VOG 0
*
*. Modify RHO1
*
      ! Nulify initial densi (RHO1) matrix
      do i = 1, nacob_bu**2
       densi1(i) = 0.0d0
      end do
      !
      ! Divide all MOs in the original GAS scheme into three groups:
      ! (1) MOs in GAS --> Now Inactive
      ! (2) MOs in GAS --> Still in GAS (transformed)
      ! (3) MOs in GAS --> Now Secondary
      !
      !TODO Check for more complicated cases (more GAS spaces)
      irow = 1
      icol = 1
      iadd = 0
*
      do i = 1, nacob_bu**2
       if (icol.gt.ioffden(1).and.icol.le.ioffden(2)  .and.
     &     irow.gt.ioffden(1).and.irow.le.ioffden(2)) then
        iadd = iadd + 1
        densi1(iadd) = densi1_bu(i)
       end if
       ! Change row
       if (icol.eq.nacob_bu) then
        icol = 0
        irow = irow + 1
       end if
       icol = icol + 1
      end do
*
*. Modify RHO2
*
      idorho2 = 0
      if (only_dens.eq.1.and.idorho2.eq.1) then
!      if (only_dens.eq.1) then
*
       ! Switch back to the original, large 2-e density matrix
       ngas_scr = ngas
       ngas = ngas_bu
       call getnddens2(ndactorb,nddensi2,only_dens)
       ngas = ngas_scr
*
       ! Nulify initial densi2 (RHO2) matrix
       do i = 1, nddensi2
        densi2(i) = 0.0d0
       end do
*
       do i = 1, 3
        ipqden(i) = ioffden(i)**2 + (ioffden(i)**2 +1)/2
       write(*,'(i1,a,i4,i5)') i,' ioffden, ipqden',ioffden(i),ipqden(i)
       end do
*      
*. Initialize pointers
*
       ioff = 0
       iadd = 0
       ij   = 1
*
*. Get diamensions
*
       nrow = ndactorb**2
*
*. Loop the low diagonal
*
      write(*,*) '   i,   j,imin,imax,jmin,jmax,istp,jstp,  ij'
       do i = 1, nrow
*
        write(*,*) i
        imin  = ioff + ioffden(1)
        imax  = ioff + ioffden(2)
        istop = ioff + ioffden(3)
*
        write(*,*) 'TI SKATA EXW', imin, imax
        if (i.le.imin) then
        write(*,'(9i5,a)') i, j, imin, imax, jmin, jmax,
     &                     istop,jstop,ij,' a'
         ij = ij + 1
         cycle
        else if (i.gt.imin.and.i.le.imax) then
*
         joff  = 0
         jstop = 0
*
         do j = 1, i
*
          jmin  = joff + ioffden(1)
          jmax  = joff + ioffden(2)
          jstop = joff + ioffden(3)
        write(*,'(9i5,a)') i, j, imin, imax, jmin, jmax,
     &                     istop,jstop,ij,' b'
*
          if (j.le.jmin) then
           ij = ij + 1
           cycle
          else if (j.gt.jmin.and.j.le.jmax) then
*
           ! Collect the RHO2 elements of the small matrix
           iadd = iadd + 1
           densi2(iadd) = densi2_bu(ij)
           if (ntest.ge.100) 
     &      write(*,'(a,2i4,e20.8)') 'Collect ',iadd,ij,densi2(iadd)
*
           if (j.eq.jstop) joff = joff + ioffden(3)
          else if (j.eq.jstop) then
*
           joff = joff + ioffden(3)
*
          end if ! Check colunm j
*
          if (i.ne.j) ij = ij +1
         end do ! End columns (j)
*
         if (i.eq.istop) ioff = ioff + ioffden(3)
        else if (i.le.istop) then
* 
         ioff = ioff + ioffden(3)
         if (i.eq.istop) ij = ij + ioff - 1 !ioffden(3) - 1
         write(*,'(9i5,a)') i, j, imin, imax, jmin, jmax,
     &                     istop,jstop,ij,' c'
*
        end if ! Check row i
*
        ij = ij + 1
       end do ! End rows (i)
*
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      subroutine restore_gas_spaces(densi1,densi1_bu,densi2,densi2_bu,
     &                              only_dens)
*======================================================================*
*
*  Konstantinos D. Vogiatzis, June 2015, Berkeley
*
*======================================================================*
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'cgas.inc'
      include 'glbbas.inc'
      include 'lucinp.inc'
#include "global.fh"
*
*. Input
      integer   only_dens
      dimension densi1_bu(*)
      dimension densi2_bu(*)
*
*. Output
      dimension densi1(*)
      dimension densi2(*)
*
*. Initializations
*
      ! Get dimensions for RHO2
      call getnddens2(ndactorb,nddensi2,only_dens)
*
      ngas  = ngas_bu
*
      nocob = nocob_bu
      nacob = nacob_bu
      ninob = ninob_bu
      nscob = nscob_bu
*
      do igensm = 1, ngensmob
       do igas= 1, ngas+2
        ngas_gnsym(igensm,igas-1) = ngas_gnsym_bu(igensm,igas)
       end do
       nocobs(igensm) = nocobs_bu(igensm)
       nacobs(igensm) = nacobs_bu(igensm)
       ninobs(igensm) = ninobs_bu(igensm)
       nscobs(igensm) = nscobs_bu(igensm)
      end do
*
      do igensm = 1, ngensmob
       nobpts_gn(0,igensm) = nobpts_gn_bu(0,igensm)
       nobpts_gn(ngas+1,igensm) = nobpts_gn_bu(ngas+1,igensm)
       do igas = 1, ngas
        nobpts_gn(igas,igensm) = nobpts_gn_bu(igas,igensm)
       end do
      end do
*
      ! Symmetry => type reordering array
      do i = 1, ntoob
       ireost(i) = ireost_bu(i)
      end do
*
*. RHO1
*
      do i = 1, nacob**2
       densi1(i) = densi1_bu(i)
      end do
*
*. RHO2
*
      if (only_dens.eq.1) then
       do i = 1, nddensi2
        densi2_bu(i) = densi2(i) !TODO PROSOXI!!!
       end do
*
      ! For DENSI2
       do igas= 1, ngas+2
        nobpt(igas)   = nobpt_bu(igas)
        ihpvgas(igas) = ihpvgas_bu(igas)
       end do
      end if
*
      NTEST = 0
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
       write(*,*) ' Initial Density (restored) '
       call wrtmat(densi1,nacob,nacob,nacob,nacob)
       write(*,*)
       if (only_dens.eq.1) call prsym(densi2,ndactorb*ndactorb)
       write(*,*)
       write(*,'(a,10i4)') 'INACTIVE ',
     &          (ninobs(isym), isym = 1, ngensmob), ninob
       write(*,'(a,10i4)') 'ACTIVE   ',
     &          (nacobs(isym), isym = 1, ngensmob), nacob
       write(*,'(a,10i4)') 'OCCUPIED ',
     &          (nocobs(isym), isym = 1, ngensmob), nocob
       write(*,'(a,10i4)') 'SECONDARY',
     &          (nscobs(isym), isym = 1, ngensmob), nscob
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      subroutine getnddens2(ndactorb,nddensi2,only_dens)
*======================================================================*
*
*  Konstantinos D. Vogiatzis, June 2015, back from Berkeley
*
*======================================================================*
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'cgas.inc' 
      include 'crun.inc'
#include "global.fh"
*
*. Input
      integer only_dens
*
*. Output
      integer ndactorb
      integer nddensi2
*
*. Local
      dimension idactspc(mxpngas)
*
      ntest = 10 !VOG 0
*
      ndactorb = 0
      if (only_dens.eq.1) then
       DO IGAS = 1, NGAS
        ITYP = IHPVGAS(IGAS)
        IDACTSPC(IGAS) = 0
        IF(ITYP.EQ.1.AND.IDENS_IN.EQ.1) IDACTSPC(IGAS) = 1
        IF(ITYP.EQ.2.AND.IDENS_SEC.EQ.1) IDACTSPC(IGAS) = 1
        IF(ITYP.EQ.3.AND.IDENS_AC.EQ.1) IDACTSPC(IGAS) = 1
       END DO
       if(ntest.ge.10.and.ga_nodeid().eq.0) then
        WRITE(6,*) ' Active orbital spaces in getnddens2 '
        CALL IWRTMA(IDACTSPC,1,NGAS,1,NGAS)
       end if
*
       do igas = 1, ngas
        if(idactspc(igas).eq.1) ndactorb = ndactorb + nobpt(igas)
       end do
      end if
*
      nddensi2 = ndactorb**2 *(ndactorb**2 + 1)/2
*
      if(ntest.ge.10.and.ga_nodeid().eq.0) then
       write(6,*) ' Number of active orbitals for densities',ndactorb
       write(6,*) ' Number of elements in RHO2             ',nddensi2
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      integer function ipq(i)
*======================================================================*
*
*. Input
      integer i
*
      if (i.le.0) then
       x = 0
      else
       x = i * (i + 1)/2
      end if
*
      ipq = int(x)
*
      return
      end
*======================================================================*
*======================================================================*
      subroutine fraggas(igasloop,ngasadd,imofromto)
*======================================================================*
*
*.  Fragment the full GASspace into smaller GAS spaces for increasign
*.  the number of CI blocks and reducing memory requirements 
*.  
*.  If I move one MO from a large GAS space, then I am using these
*.  pointers:
*.  imofromto(1): irrep of the MO that I will move
*.  imofromto(2): igas from which I will move the MO
*.  imofromto(3): igas to which I will move the MO
*.
*.  Konstantinos D. Vogiatzis, September 2015 
*.                             February  2016 GASCI update
*
*======================================================================*
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'cgas.inc'
      include 'cstate.inc'
      include 'crun.inc'
      include 'clunit.inc'
      include 'lucinp.inc'
#include "global.fh"
#include "mafdecls.fh"
*
*. Input
      integer igasloop
*
*. Input/Output
      integer ngasadd, imofromto(3)
*
*. Scratch
      integer igsoccupd(mxpngas,2,ncispc)
      integer ngsshupd(nirrep,mxpngas) 
      integer gasocc(ngas,2), sumocc(2)
      integer gathermos(nirrep)
      integer ntotmos(ngas)
      integer maxocc(ngas)
      integer totsds, ngastotmo, weyl, get_mo_occ_gas
      external weyl, get_mo_occ_gas
*
      ntest  = 100
*
*. Print Input
*
      if (ntest.eq.10.and.ga_nodeid().eq.0)
     & WRITE(*,'(a,i4)') ' Times that I have entered in fraggas:',
     &                   igasloop
      if(ntest.ge.100.and.ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ***********************************'
      WRITE(*,*) '   Input GAS Spaces for CI batches  '
      WRITE(*,*) ' ***********************************'
      WRITE(*,*)
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      WRITE(6,*)
      IF(IGSFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Gas space provided by courtesy of LUCIA:',  IGSFILL
      IF(ISECFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Secondary space provided by courtesy of LUCIA '
*
       WRITE(*,'(A,I3)')
     &  ' Number of Occupation spaces: ',NCISPC
       WRITE(*,*)
       DO ICISPC = 1, NCISPC
        WRITE(*,'(A,I3)')
     &  ' Bounds on accumulated occupations for space: ',ICISPC
        WRITE(*,'(A)')
     &   ' ====================================================== '
        WRITE(*,'(A)')
        WRITE(*,'(A)') '         Min. occ    Max. occ '
        WRITE(*,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(*,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
        END DO
       END DO
       WRITE(*,*)
      end if
*
      ntogas = 100
      nmoqs  = 0
*
*. Identify non-filled GAS spaces
*
       ntotmos = 0
      nntotmos = 0
      do igas = 1, ngas
       do igensm = 1, ngensmob
        ntotmos(igas) = ntotmos(igas) + ngssh(igensm,igas)
       end do
       maxocc(igas) = ntotmos(igas) * 2
       nntotmos = nntotmos + ntotmos(igas)
      end do
*
* Find which GAS spaces belong to P- and Q- spaces of AKBK
* by using the accumulated occupations of CI space 1
*
      lendp    = 0
      ngasp    = 0
      ngasq1   = 0
      ngasq2   = 0
      nmogasp  = 0
      nmogasq2 = 0
*
      do igas = 1, ngas
*
       !
       ! CASE 1: Fully occupied GAS space belonging to Q-space (Q1)
       !
       if ((igsoccx(igas,1,1).eq.maxocc(igas))
     &  .and.(igsoccx(igas,2,1).eq.maxocc(igas)).and.(lendp.eq.0)) then
        if ((ntest.ge.10).and.(ga_nodeid().eq.0))
     &     write(*,'(a,i4,a)') ' GAS', igas, ': Occupied GAS of Q-Space'
        ngasq1 = ngasq1 + 1
       !
       ! CASE 2: GAS space that is treated as part of the active space
       !         for MCSCF and P-Space for AKBK
       !
       else if ((igsoccx(igas,1,1).le.nactel)
     &     .and.(igsoccx(igas,2,1).le.nactel).and.(lendp.eq.0)) then
        if ((ntest.ge.10).and.(ga_nodeid().eq.0))
     &     write(*,'(a,i4,a)') ' GAS', igas, ': GAS of P-Space'
        ngasp = ngasp + 1
*
        do irrep = 1, nirrep
         nmogasp = nmogasp + ngssh(irrep,igas)
        end do
*
        if ((igsoccx(igas,1,1).eq.nactel).and.
     &     (igsoccx(igas,2,1).eq.nactel)) lendp = 1
       !
       ! CASE 3: Empty GAS space belonging to Q-space (Q2)
       !
       else if (lendp.eq.1) then
        if ((ntest.ge.10).and.(ga_nodeid().eq.0))
     &     write(*,'(a,i4,a)') ' GAS', igas, ': Empty GAS of Q-Space'
        ngasq2 = ngasq2 + 1
*
        do irrep = 1, nirrep
         nmogasq2 = nmogasq2 + ngssh(irrep,igas)
        end do
*
       end if
*
      end do
*
*. Check in how many subspaces I want to split the given GAS spaces
*. based on a rough estimate of the number of SDs (Weyl's formula)
*
!      totsds = weyl(nactel,nmogasp,ms2)
!      totsds = totsds / nirrep
!      idigits = int(log10(dble(totsds))+1)
!      if ((ntest.ge.10).and.(ga_nodeid().eq.0))
!     &  write(*,*) 'SDs from Weyl (per irrep)', totsds, idigits
*
*. Begin the fragmentation of the GAS spaces into smaller GAS
*
        ngasupd = 0
       ngsshupd = 0
      igsoccupd = 0
      gathermos = 0
      ishdirrep = 0
*
*. Two different types of calculations are considered:
*
*     ============
*       I. MCSCF
*     ============
*
      if (CSEQCI(1,2).ne.'AKBKCI') then
*
       if (ngasq1.ne.0.and.ngasq2.ne.0) then
        ngasq1 = ngasq1 + ngasq2
        write(*,'(a,i3)') 'MCSCF: I found MOs in Q-Space'
        CALL ERRQUIT('fraggas: wrong number of MOs', ngasq1,667)
       end if
*
       mxfraggas = mxpngas
       ngasold = ngas
       icispc = 1
*
       if (ngas.lt.mxfraggas) then
*
        if (ntest.ge.100.and.ga_nodeid().eq.0) write(*,*) 'imofromto', 
     &     (imofromto(i), i = 1, 3)
        ngasupd = ngas + 1
        !
        ! Copy MOs (irrep/GAS space) before the reduced GAS
        !
        do igas = 1, imofromto(3)-1
         do irrep = 1, nirrep
          ngsshupd(irrep,igas) = ngssh(irrep,igas)
         end do
        !
        ! Min/Max Occ
         igsoccupd(igas,1,icispc) = igsoccx(igas,1,icispc)
         igsoccupd(igas,2,icispc) = igsoccx(igas,2,icispc)
         if (igas.eq.imofromto(2)) then
          igsoccupd(igas,1,icispc) = igsoccupd(igas,1,icispc) - 2
!     &                               2*ngssh(imofromto(1),imofromto(2))
          if (((igas-1).gt.0).and.
     &       (igsoccupd(igas,1,icispc).lt.igsoccupd(igas-1,1,icispc)))
     &        igsoccupd(igas,1,icispc) = igsoccupd(igas-1,1,icispc)
          if (igsoccupd(igas,1,icispc).lt.0)
     &        igsoccupd(igas,1,icispc) = 0
         end if
        end do
        !
        ! Move target MOs in a new GAS space
        !
        ngsshupd(imofromto(1),imofromto(2)) = 
     &  ngsshupd(imofromto(1),imofromto(2)) - 1
        ngsshupd(imofromto(1),imofromto(3)) = 1
        !
        ! Min/Max Occ
        igsoccupd(imofromto(3),1,icispc) = 
     &    igsoccx(imofromto(3)-1,1,icispc)
        igsoccupd(imofromto(3),2,icispc) =  
     &    igsoccx(imofromto(3)-1,2,icispc)
        !
        ! Reduce the min occ of the intermediate GAS spaces,
        ! from imofromto(2) till imofromto(3)
        do igas = imofromto(2)+1, imofromto(3)-1
         igsoccupd(igas,1,icispc) = igsoccx(igas-1,1,icispc)
         igsoccupd(igas,2,icispc) = igsoccx(igas-1,2,icispc)
        end do
        !
        ! Copy MOs (irrep/GAS space) after the reduced GAS
        !
        if (imofromto(3).ne.ngasupd) then
         do igas = imofromto(3), ngas
          do irrep = 1, nirrep
           ngsshupd(irrep,igas+1) = ngssh(irrep,igas)
          end do
        !
        ! Min/Max Occ
          igsoccupd(igas+1,1,icispc) = igsoccx(igas,1,icispc)
          igsoccupd(igas+1,2,icispc) = igsoccx(igas,2,icispc)
         end do
        end if
*
*. Find from which GAS I will move the next orbital if I need to further
*. fragment the original active space
*
        MAXMOINGAS = 0
        DO IRREP = 1, NIRREP
         IF (IRREP.NE.imofromto(1)) THEN
          IF (ngssh(irrep,imofromto(2)).GT.MAXMOINGAS) THEN
           MAXMOINGAS = ngssh(irrep,imofromto(2))
           imofromto(1) = IRREP
          END IF
         END IF
        END DO
!        DO IRREP = 1, NIRREP
!         IF (IRREP.NE.imofromto(1)) THEN
!          IF (ngssh(irrep,imofromto(2)).GT.1) THEN
!           imofromto(1) = IRREP
!           exit
!          END IF
!         END IF
!        END DO
        !
        ! By moving only 1 MO in a new GAS, the imofromto(2) is always 1
        ! (i.e. we move always a MO from GAS 1)
        imofromto(2) = 1
        imofromto(3) = ngasupd + 1
*
*. Case that I have reached the max number of GAS spaces
*. add an MO to an existing one (not tested for GASCI yet!)
       else 
        write(*,*) ' I will add Orb in an existing GAS'
        if (ngasadd.eq.mxfraggas) ngasadd = 1
        ngasadd = ngasadd + 1
        ngasupd = ngasadd
       end if 
!*
!*. Copy...
!*
!       do irrep = 1, nirrep
!        do igas = 1, ngasold
!         ngsshupd(irrep,igas) = ngssh(irrep,igas) 
!        end do
!       end do
!*
!*. ...and move that orbital
!*
!!! igas=1 means that the user requested CASSCF, not GASSCF (fix it!)
!       ngsshupd(ifrgasirrep,1)       = ngsshupd(ifrgasirrep,1)       - 1
!       ngsshupd(ifrgasirrep,ngasupd) = ngsshupd(ifrgasirrep,ngasupd) + 1

*************************************************************************
*** OLD Implementation based on Weyl's formula
*************************************************************************
!       if (idigits.le.3) then
!        maxgas = 1
!       else if (idigits.eq.4) then
!        maxgas = 2
!       else if (idigits.eq.5) then
!        maxgas = 4
!       else if (idigits.eq.6) then
!        maxgas = 8
!       else if (idigits.gt.6) then
!        maxgas = mxpngas
!       end if
!!      maxgas = mxpngas
!*
!*. Start work
!*
!       navailgas = maxgas
!*
!       do igasp = 1, ngasp
!*
!*. Get max MO occupation for the first GAS
!*
!        ngastotmo = get_mo_occ_gas(nmogasp,navailgas)
!*
!        do irrep = 1, nirrep
!*
!         do ish = 1, ngssh(irrep,igasp)
!          gathermos(irrep) = gathermos(irrep) + 1
!          ishdirrep = ishdirrep + 1
!          if (ishdirrep.eq.ngastotmo) then
!           !
!           ! I have reached the max MO occupation of the new GAS
!           !
!           ngasupd = ngasupd + 1
!           do jrrep = 1, irrep
!            ngsshupd(jrrep,ngasupd) = gathermos(jrrep)
!           end do
!           !
!           ! Nullify and redefine MO occupation
!           !
!           gathermos = 0
!           ishdirrep = 0
!*
!           nmogasp = nmogasp - ngastotmo
!           navailgas = navailgas - 1
!           if (navailgas.ne.0)
!     &         ngastotmo = get_mo_occ_gas(nmogasp,navailgas)
!          end if
!         end do
!*
!*DISCARD?
!*
!*          ! nfrgsh = number of fragments for this GAS shell
!*          nfrgsh = int(ngssh(irrep,igasp) / (moratio + 1))
!*          do ish = 1, nfrgsh
!*           ! Move to the next GAS
!*           ngasupd = ngasupd + 1
!*           ngsshupd(irrep,ngasupd) = moratio + 1
!*           write(*,*) 'KATI'
!*          end do
!**
!*          ! Remaining MOs in separate GAS
!*          nremain = ngssh(irrep,igasp) - (nfrgsh*(moratio + 1))
!*          if (nremain.ne.0) then
!*           ! Move to the next GAS
!*           ngasupd = ngasupd + 1
!*           ngsshupd(irrep,ngasupd) = nremain
!*           write(*,*) 'KATI' 
!*          end if
!**
!*          ! Update remaining MOs to be distributed for next GAS/irrep
!*          nmogasp = nmogasp - ngssh(irrep,igasp)
!*          navailgas = maxgas - ngasupd
!*
!        end do ! End Irrep loop
!       end do ! End GAS P-space loop
*************************************************************************
*
*. Overwrite GAS information with updated values before exiting
*
       ngas = ngasupd
       ngssh = 0
       do igas = 1, ngas
        do irrep = 1, nirrep
         ngssh(irrep,igas) = ngsshupd(irrep,igas)
        end do
       end do
*
*. Overwrite Min/Max occ numbers
*
       igsoccx = 0
       mxoccgs = 0
       icispc = 1
*
       do igas = 1, ngas
        igsoccx(igas,1,icispc) = igsoccupd(igas,1,icispc)  
        igsoccx(igas,2,icispc) = igsoccupd(igas,2,icispc)
!        ! Min
!        if (igas.eq.ngas) then
!         igsoccx(igas,1,icispc) = nactel
!        else
!         igsoccx(igas,1,icispc) = 0
!        end if
!        ! Max
!        do irrep = 1, nirrep
!         mxoccgs = mxoccgs + 2*ngssh(irrep,igas)
!        end do
!        if (mxoccgs.ge.nactel) then
!         igsoccx(igas,2,icispc) = nactel
!        else
!         igsoccx(igas,2,icispc) = mxoccgs
!        end if
       end do
*
*     ============
*     II. SplitGAS
*     ============
*
      else
*
*. Start work for P-space
*
*. In how many GAS spaces shall I split the P-space?
*
       if (idigits.le.4) then
        maxgas = 1
       else if (idigits.eq.5) then
        maxgas = 2
       else if (idigits.eq.6) then
        maxgas = 4
       else if (idigits.gt.7) then
        maxgas = 8
       end if
*
       navailgas = maxgas
*
       do igasp = 1, ngasp
*
*. Get max MO occupation for the first GAS
*
        ngastotmo = get_mo_occ_gas(nmogasp,navailgas)
*
        do irrep = 1, nirrep
*
         do ish = 1, ngssh(irrep,igasp)
          gathermos(irrep) = gathermos(irrep) + 1
          ishdirrep = ishdirrep + 1
          if (ishdirrep.eq.ngastotmo) then
           !
           ! I have reached the max MO occupation of the new GAS
           !
           ngasupd = ngasupd + 1
           do jrrep = 1, irrep
            ngsshupd(jrrep,ngasupd) = gathermos(jrrep)
           end do
           !
           ! Nullify and redefine MO occupation
           !
           gathermos = 0
           ishdirrep = 0
*
           nmogasp = nmogasp - ngastotmo
           navailgas = navailgas - 1
           if (navailgas.ne.0)
     &         ngastotmo = get_mo_occ_gas(nmogasp,navailgas)
          end if
         end do
        end do ! End Irrep loop
       end do ! End GAS P-space loop
*
*. Store number of GAS for P-space
*
       ngas = ngasupd
*
*. Start work for Q-space
*
*. In how many GAS spaces shall I split the Q-space?
*
       if (nmogasq2.le.1000) then
        maxgasq = 8
       else
        maxgasq = mxpngas - maxgas
       end if
*
       navailgas = maxgasq
*
       gathermos = 0
       ishdirrep = 0
*
       do igasq2 = 1, ngasq2
*
*. Get max MO occupation for the first GAS
*
        ngastotmo = get_mo_occ_gas(nmogasq2,navailgas)
*
        do irrep = 1, nirrep
*
         do ish = 1, ngssh(irrep,ngasp+igasq2)
          gathermos(irrep) = gathermos(irrep) + 1
          ishdirrep = ishdirrep + 1
          if (ishdirrep.eq.ngastotmo) then
           !
           ! I have reached the max MO occupation of the new GAS
           !
           ngasupd = ngasupd + 1
           do jrrep = 1, irrep
            ngsshupd(jrrep,ngasupd) = gathermos(jrrep)
           end do
           !
           ! Nullify and redefine MO occupation
           !
           gathermos = 0
           ishdirrep = 0
*
           nmogasq2 = nmogasq2 - ngastotmo
           navailgas = navailgas - 1
           if (navailgas.ne.0)
     &         ngastotmo = get_mo_occ_gas(nmogasq2,navailgas)
          end if
         end do
        end do ! End Irrep loop
       end do ! End GAS Q-space loop
*
*. Overwrite GAS information with updated values
*
       ngssh = 0
       do igas =  1, ngas + ngasupd
        do irrep = 1, nirrep
         ngssh(irrep,igas) = ngsshupd(irrep,igas)
        end do
       end do
*
*. Store minimum/maximum ocupations
*
       igsoccx = 0
       mxoccgs = 0
*
*. CISPACE 1
*
       icispc = 1
*
       !
       ! Occupation numbers for CASCI step
       ! P-space
       !
       do igas = 1, ngas
        ! Min
        if (igas.eq.ngas) then
         igsoccx(igas,1,icispc) = nactel
        else
         igsoccx(igas,1,icispc) = 0
        end if
        ! Max
        do irrep = 1, nirrep
         mxoccgs = mxoccgs + 2*ngssh(irrep,igas)
        end do
        if (mxoccgs.ge.nactel) then
         igsoccx(igas,2,icispc) = nactel
        else
         igsoccx(igas,2,icispc) = mxoccgs
        end if
       end do
       !
       ! Occupation numbers for CASCI step
       ! Q-space has 0 occupation
       !
       do igas = ngas+1, ngasupd
        ! Min
        igsoccx(igas,1,icispc) = nactel
        ! Max
        igsoccx(igas,2,icispc) = nactel
       end do
*
*. CISPACE 2
*
       icispc = 2
*
       !
       ! Occupation numbers for SplitGAS
       ! P-space
       !
       do igas = 1, ngas
        ! Min
        if (igas.eq.ngas) then
         igsoccx(igas,1,icispc) = nactel - 2
        else
         igsoccx(igas,1,icispc) = 0
        end if
        ! Max
        do irrep = 1, nirrep
         mxoccgs = mxoccgs + 2*ngssh(irrep,igas)
        end do
        if (mxoccgs.ge.nactel) then
         igsoccx(igas,2,icispc) = nactel
        else
         igsoccx(igas,2,icispc) = mxoccgs
        end if
       end do
       !
       ! Occupation numbers for SplitGAS
       ! 2 particles in Q-space
       !
       do igas = ngas+1, ngasupd
        ! Min
        if (igas.eq.ngasupd) then
         igsoccx(igas,1,icispc) = nactel
        else
         igsoccx(igas,1,icispc) = nactel - 2
        end if
        ! Max
        igsoccx(igas,2,icispc) = nactel
       end do
*
*. Add number of GAS spaces in Q-space for the total number of GAS
*
       ngas = ngasupd
*
      end if
*
*. Print Output
*
      if(ntest.ge.10.and.ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ***********************************'
      WRITE(*,*) '  Updated GAS Spaces for CI batches '
      WRITE(*,*) ' ***********************************'
      WRITE(*,*)
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      WRITE(6,*)
      IF(IGSFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Gas space provided by courtesy of LUCIA:',  IGSFILL
      IF(ISECFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Secondary space provided by courtesy of LUCIA '
*
       WRITE(*,'(A,I3)')
     &  ' Number of Occupation spaces: ',NCISPC
       WRITE(*,*)
       DO ICISPC = 1, NCISPC
        WRITE(*,'(A,I3)')
     &  ' Bounds on accumulated occupations for space: ',ICISPC
        WRITE(*,'(A)')
     &   ' ====================================================== '
        WRITE(*,'(A)')
        WRITE(*,'(A)') '         Min. occ    Max. occ '
        WRITE(*,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(*,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
        END DO
       END DO
       WRITE(*,*)
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      integer function get_mo_occ_gas(nmogasp,navailgas)
      implicit none
*======================================================================*
*
*.  Find how many MOs each GAS should have (ngastotmo)
*.  ngastotmo can have MOs from different irreps
*
*.  Konstantinos D. Vogiatzis, October 2015
*
*======================================================================*
*
*. Input
      integer nmogasp, navailgas
*
*. Scratch
      integer moratio, moremo
*
*. Work
*
      moratio = int(nmogasp / navailgas)
      moremo = nmogasp - (moratio * navailgas)
*
      if (moremo.eq.0) then
       get_mo_occ_gas = moratio
      else
       get_mo_occ_gas = moratio + 1
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      integer function weyl(nele,nmo,spin)
      implicit none
#include "global.fh"
*======================================================================*
*
*.  Calculate the number of Slater Determinants with Weyl's formula
*.  without considering symmetry.
*.  n = number of electrons
*.  N = number of MOs in the active space
*.  S = Spin
*
*.  Konstantinos D. Vogiatzis, October 2015
*
*======================================================================*
*
*. Input
      integer nele, nmo, spin
*
*. Scratch
      integer inter, inter2, nom, fac1, fac2, denom
      double precision one, two, three, rweyl, rnele, rnmo, rspin
      integer ntest,inter1
      external vogfactorial
      integer vogfactorial
*
      ntest = 10
*
**. Work
*
      ! Transform integers to real
      rnele = dble(nele)
      rnmo  = dble(nmo)
      rspin = dble(spin)
*
*     First term = 2S + 1
*
      one = 2*rspin + 1.0d0
*
*     Second term = N! / ( (n/2 - S)! * (N + 1 - n/2 + S)! )
*
      nom = vogfactorial(nmo)
      inter1 = int(nele / 2 - spin)
      inter2 = nmo + 1 - inter1
      fac1 = vogfactorial(inter1)
      fac2 = vogfactorial(inter2)
      denom = fac1 * fac2
      two = dble(nom) / dble(denom)
*
*     Third term = (N+1)! / ( (n/2 + S + 1)! * (N - n/2 - S)! )
*
      nom = vogfactorial(nmo + 1)
      inter1 = int(nele / 2 + spin + 1)
      inter2 = nmo + 1 - inter1
      fac1 = vogfactorial(inter1)
      fac2 = vogfactorial(inter2)
      denom = fac1 * fac2
      three = dble(nom) / dble(denom)
*
*     Multiply the three terms for the number of SDs
*
      rweyl = one * two * three
      weyl = int(rweyl)
      if ((ntest.ge.10).and.(ga_nodeid().eq.0))
     &  write(*,*) 'Approximate number of SDs from Weyl:', weyl
*
      return
      end
*======================================================================*
*======================================================================*
      integer function vogfactorial(n)
      implicit none
*======================================================================*
*
*. Input
      integer n
*
*. Scratch
      integer i
*
      vogfactorial = 1
      do i = 1, n
       vogfactorial = vogfactorial * i
      end do 
*
      return
      end
*======================================================================*
*======================================================================*
      subroutine gather_mos_in_cas
*======================================================================*
*
*.  Gather all MOs of the active space distributed in various GAS 
*.  spaces into one space
*.
*.  Konstantinos D. Vogiatzis, October 2015
*
*======================================================================*
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'cgas.inc'
      include 'crun.inc'
      include 'clunit.inc'
      include 'lucinp.inc'
#include "global.fh"
#include "mafdecls.fh"
*
*. Scratch
      integer ngsshupd(nirrep), ngas_gnsymupd(nirrep,2)
*
      ntest = 1
*
*. Gather
*
      ngsshupd = 0
      ngas_gnsymupd = 0
      do irrep = 1, nirrep
       do igas =  1, ngas
        ngsshupd(irrep) = ngsshupd(irrep) + ngssh(irrep,igas)
        ngas_gnsymupd(irrep,1) =
     &  ngas_gnsymupd(irrep,1) + ngas_gnsym(irrep,igas)
       end do
      end do
      do irrep = 1, nirrep !Secondary MOs
       ngas_gnsymupd(irrep,2) = ngas_gnsym(irrep,ngas+1)
      end do
*
*. Backup
*
      ngas_bu = ngas
      ngas_gnsym_bu = 0
      ngssh_bu      = 0
      do igensm = 1, ngensmob
       do igas= 1, ngas+2
        ngas_gnsym_bu(igensm,igas) = ngas_gnsym(igensm,igas-1)
       end do
      end do
      do irrep = 1, nirrep
       do igas =  1, ngas
        ngssh_bu(irrep,igas) = ngssh(irrep,igas)
       end do
      end do
*
      nobpts_gn_bu = 0
      do igensm = 1, ngensmob
       nobpts_gn_bu(0,igensm) = nobpts_gn(0,igensm)
       nobpts_gn_bu(ngas+1,igensm) = nobpts_gn(ngas+1,igensm)
       do igas = 1, ngas
        nobpts_gn_bu(igas,igensm) = nobpts_gn(igas,igensm)
       end do
      end do
*
      itpfso_bu = 0
      do itoob = 1, ntoob
       itpfso_bu(itoob) = itpfso(itoob)
      end do
*
      icispc = 1
      igsoccx_bu = 0
      do igas =  1, ngas
       igsoccx_bu(igas,1,icispc) = igsoccx(igas,1,icispc)
       igsoccx_bu(igas,2,icispc) = igsoccx(igas,2,icispc)
      end do
*
*. Nullify
*
      ngssh = 0
      do irrep = 1, nirrep
       do igas =  1, ngas+1
        ngas_gnsym(irrep,ngas) = 0
        if (igas.le.ngas) nobpts_gn(igas,irrep) = 0
       end do
      end do
*
*. Copy!
*
      ! Change MO pointer of itpfso for Active and Secondary
      do itoob = 1, ntoob
       if ((itpfso(itoob).gt.0).and.(itpfso(itoob).ne.ngas+1)) then
        itpfso(itoob) = 1
       else if (itpfso(itoob).eq.ngas+1) then
        itpfso(itoob) = 2
       end if
      end do
*
      ! And the rest
      ngas = 1
      do irrep = 1, nirrep
       ngssh(irrep,ngas) = ngsshupd(irrep)
       ngas_gnsym(irrep,ngas)   = ngas_gnsymupd(irrep,1)
       ngas_gnsym(irrep,ngas+1) = ngas_gnsymupd(irrep,2)
       nobpts_gn(ngas,irrep) = ngsshupd(irrep)
       nobpts_gn(ngas+1,irrep) = nobpts_gn(ngas+2,irrep)
       nobpts_gn(ngas+2,irrep) = 0
      end do
*
      ! ...and do not forget the occupation numbers
      icispc = 1
      igsoccx = 0
      igsoccx(1,1,icispc) = nactel
      igsoccx(1,2,icispc) = nactel
*
*. Print Output
*
      if(ntest.ge.10.and.ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ***********************************************'
      WRITE(*,*) '  Gathered MOs in one GAS Space Before MOROT_GS '
      WRITE(*,*) ' ***********************************************'
      WRITE(*,*)
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      IF(ISECFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Secondary space provided by courtesy of LUCIA '
*
       WRITE(*,'(A,I3)')
     &  ' Number of Occupation spaces: ',NCISPC
       WRITE(*,*)
       DO ICISPC = 1, NCISPC
        WRITE(*,'(A,I3)')
     &  ' Bounds on accumulated occupations for space: ',ICISPC
        WRITE(*,'(A)')
     &   ' ====================================================== '
        WRITE(*,'(A)')
        WRITE(*,'(A)') '         Min. occ    Max. occ '
        WRITE(*,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(*,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
        END DO
       END DO
       WRITE(*,*)
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      subroutine return_values_after_morot
*======================================================================*
*
*.  After the successful transformation of the active MOs in MOROT_GS,
*.  return pointers and dimensions of variables corresponding to
*.  multiple GAS spaces
*.
*.  Konstantinos D. Vogiatzis, October 2015, after hitting the local
*.  bookstore with Jeppe
*
*======================================================================*
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'cgas.inc'
      include 'crun.inc'
      include 'clunit.inc'
      include 'lucinp.inc'
#include "global.fh"
#include "mafdecls.fh"
*
      ngas = ngas_bu
*
*. Backup
*
      ngssh = 0
      do igensm = 1, ngensmob
       do igas= 1, ngas+2
        ngas_gnsym(igensm,igas-1) = ngas_gnsym_bu(igensm,igas)
       end do
      end do
*
      do irrep = 1, nirrep
       do igas =  1, ngas
        ngssh(irrep,igas) = ngssh_bu(irrep,igas)
       end do
      end do
*
      nobpts_gn = 0
      do igensm = 1, ngensmob
       nobpts_gn(0,igensm)      = nobpts_gn_bu(0,igensm)
       nobpts_gn(ngas+1,igensm) = nobpts_gn_bu(ngas+1,igensm)
       do igas = 1, ngas
        nobpts_gn(igas,igensm)  = nobpts_gn_bu(igas,igensm)
       end do
      end do
*
      itpfso = 0
      do itoob = 1, ntoob
       itpfso(itoob) = itpfso_bu(itoob)
      end do
*
      icispc = 1
      igsoccx = 0
      do igas =  1, ngas
       igsoccx(igas,1,icispc) = igsoccx_bu(igas,1,icispc)
       igsoccx(igas,2,icispc) = igsoccx_bu(igas,2,icispc)
      end do
*
      ntest = 00
*
      if(ntest.ge.10.and.ga_nodeid().eq.0) then
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      WRITE(6,*)
      IF(IGSFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Gas space provided by courtesy of LUCIA:',  IGSFILL
      IF(ISECFILL.NE.0) WRITE(6,'(7X,A,I3)')
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      subroutine optimum_gas_for_ciblocks
*======================================================================*
*
*.  Write something
*
*.  Info on each TTSS block:
*
*.  KCIBT(1,*) : Alpha supergroup type
*.  KCIBT(2,*) : Beta supergroup type
*.  KCIBT(3,*) : Sym of alpha
*.  KCIBT(4,*) : Sym of beta
*.  KCIBT(5,*) : Offset of blk with respect to start of batch in expanded form
*.  KCIBT(6,*) : Offset of blk with respect to start of batch in packed form
*.  KCIBT(7,*) : Length of blk, expanded form
*.  KCIBT(8,*) : Length of blk, packed form
*
*.
*.  Konstantinos D. Vogiatzis, November 2015, Berkeley
*
*======================================================================*
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'cgas.inc'
      include 'crun.inc'
      include 'clunit.inc'
      include 'lucinp.inc'
      include 'cstate.inc'
      include 'cicisp.inc'
      include 'glbbas.inc' 
      include 'gasstr.inc'
      include 'strinp.inc'
#include "global.fh"
#include "mafdecls.fh"
*
*. Local
      integer imofromto(3) !imofromto(1) = irrep of the MO 
                           !imofromto(2) = igas FROM which we will move
                           !imofromto(3) = igas TO which we will move MO
*
      NTEST  = 100
      IPRCIX = 0
      IPRSTR = 0
      IPRORB = 0
*
* Original Values
*
      if(ntest.ge.10.and.ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*) '  Pre-calculating supergroups for Optimum CI batches '
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*)
      end if
*
      if(ntest.ge.100.and.ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ***********************************'
      WRITE(*,*) '   Input GAS Spaces for CI batches  '
      WRITE(*,*) ' ***********************************'
      WRITE(*,*)
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      WRITE(6,*)
      IF(IGSFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Gas space provided by courtesy of LUCIA:',  IGSFILL
      IF(ISECFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Secondary space provided by courtesy of LUCIA '
*
       WRITE(*,'(A,I3)')
     &  ' Number of Occupation spaces: ',NCISPC
       WRITE(*,*)
       DO ICISPC = 1, NCISPC
        WRITE(*,'(A,I3)')
     &  ' Bounds on accumulated occupations for space: ',ICISPC
        WRITE(*,'(A)')
     &   ' ====================================================== '
        WRITE(*,'(A)')
        WRITE(*,'(A)') '         Min. occ    Max. occ '
        WRITE(*,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(*,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
        END DO
       END DO
       WRITE(*,*)
      end if
*
*. Case that the user performs a CASSCF calculation and he/she has
*. requested only one GAS space:
*. LUCIA-at-NWCHEM now offers a completely free distribution of the MOs
*. belonging to the same irrep to individual GAS spaces!
*
*      IF (NGAS.EQ.1.AND.NIRREP.NE.1) CALL ONE_IRREP_PER_GAS
*      call ga_sync()
*
*. Cheat LUCIA for testing performance with large number of processors
*
      NNOD = ga_nnodes()
!      NNOD = 2048
*
*. Only for cases that are not with C1 symmetry
*
      IF (NGAS.LT.NIRREP) THEN
*
*. Obtain initial values for the MOs-to-be-moved pointers (imofromto)
*
       imofromto(1) = 0
       imofromto(2) = 0
       imofromto(3) = 0
       !
       ! GAS space from where we will move the MOs
       !
       nmoings = 0
       do irrep = 1, nirrep
        nmoings = nmoings + ngssh(irrep,1)
       end do
       if (nmoings.gt.1) then
        imofromto(2) = 1
       else
        STOP ' GAS 1 is empty!!! No MOs to move, I QUIT!!! '
       end if
       !
       ! GAS space in which we will move the MO
       !
       imofromto(3) = ngas + 1
       !
       ! MOs of the irrep that I will move to imofromto(3)
       !
       imaxmo = 0
       do irrep = 1, nirrep
        if (ngssh(irrep,imofromto(2)).gt.imaxmo) then
         imaxmo = ngssh(irrep,imofromto(2))
         imofromto(1) = irrep
        end if
       end do
*
      ELSE
       imofromto(2) = 1
      END IF
*
* ====================
* Begin the large Loop 
* ====================
*
*. Irrep from which one MO will be moved to another GAS
      ngasadd = 0
      nloop = 15 
      do igasloop = 1, nloop
*
*. Set memory pointer
*
       CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'OPTIMUM')
*
*. From shells to orbitals
       CALL ORBINF(LUOUT,IPRORB)
*
*. 1. Occupation classes
*. 2. Group of Strings
*. 3. Types (two are needed, not six)
*. 4. Supergroups
*
       CALL STRTYP_GAS(IPRSTR)
       CALL ALLO_ALLO_0
*. Divide orbital spaces into inactive/active/secondary
       CALL GASSPC
*. Symmetry information
       CALL SYMINF(IPRORB)
*
       CALL GEN_OCCLS_FOR_CISPAC(IGSOCC(1,1), IGSOCC(1,2), -1,
     &      int_mb(KIOCCLS),NOCCLS_MAX)
       if (ntest.ge.1000) call WRITE_OCCLS_AKBK(int_mb(KIOCCLS),
     &                                        NOCCLS_MAX,NGAS)
*
*. 5. Strings
*
       CALL STRINF_GAS(IPRSTR)
*
*. 6. TTSSs (from combinations of supergroups)
*
       CALL LCISPC(IPRCIX)
*
*. 7. Prepare allocations for Z_BLKFO_FOR_CISPACE
*
       LEN_IOIO= NSPGPFTP_MAX**2
       CALL MEMMAN(KCIOIO,LEN_IOIO,'ADDL  ',1,'CIOIO ')
*
       NTTS  = MXNTTS
       NSMST = NIRREP
       CALL MEMMAN(KCIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
       CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',1,'CBLTP ')
       CALL MEMMAN(KCLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
       CALL MEMMAN(KCLEBT,NTTS  ,'ADDL  ',1,'CLEBT ')
       CALL MEMMAN(KCLBLK,NTTS  ,'ADDL  ',1,'CLBLK ')
       CALL MEMMAN(KCI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
       CALL MEMMAN(KCNOCCLS_BAT,NOCCLS_MAX,'ADDL  ',1,'NOCBAT')
       CALL MEMMAN(KCIBOCCLS_BAT,NOCCLS_MAX,'ADDL  ',1,'NOCBAT')
*
*. 8. Prepare variables and call Z_BLKFO_FOR_CISPACE
*
       ICSPC  = 1
       ICSM   = IREFSM
       ICOMP  = 0
       ILTEST = -3006
*
*. Block for storing complete or partial CI-vector
       IF(ISIMSYM.EQ.1.OR.ICISTR.EQ.2) THEN
        LBLOCK = MXSOOB_AS
       ELSE
        LBLOCK = MXSOOB
       END IF
       LBLOCK = MAX(LBLOCK,LCSBLK)
       IF(NTEST.GE.10.and.ga_nodeid().eq.0)
     &         WRITE(6,*) ' TEST: LBLOCK = ', LBLOCK
*
*. First Z_FLOCK and CIV2 calls
*. Obtain TTSS batches and calculate those with very few blocks
* TODO: remove first Z_BLK is obsolete now
*
       CALL Z_BLKFO_FOR_CISPACE(ICSPC,ICSM,LBLOCK,ICOMP,
     &       NTEST,NCBLOCK,NCBATCH,int_mb(KCIOIO),int_mb(KCBLTP),
     &       NCOCCLS_ACT,int_mb(KCIOCCLS_ACT),
     &       int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &       int_mb(KCIBT),
     &       int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
*
       IBATCHCOUNT = NBATCHES_FEW_BLOCKS(NCBATCH,int_mb(KCLBT))
       NDISTBATCH = NCBATCH - IBATCHCOUNT
       IOPTBATCHPNODE = 1 !BE CAREFUL - might need to be revised or defined differently
       LBPNODE = INT(NCBATCH/NNOD)

!       if (ga_nodeid().eq.0) write(*,*) 'NCBATCH', NCBATCH
!       if (ga_nodeid().eq.0) write(*,*) 'IBATCHCOUNT', IBATCHCOUNT
!       if (ga_nodeid().eq.0) write(*,*) 'NCBLOCK', NCBLOCK
!       if (ga_nodeid().eq.0) write(*,*) 'NDISTBATCH', NDISTBATCH
!       if (ga_nodeid().eq.0) write(*,*) 'LBPNODE',LBPNODE,IOPTBATCHPNODE
*
*. Second Z_BLOCK and CIV2 calls with estimated number of blocks per batch
*
        CALL Z_BLKFO_FOR_CISPACE(ICSPC,ICSM,LBLOCK,ICOMP,
     &       NTEST,NCBLOCK,NCBATCH,int_mb(KCIOIO),int_mb(KCBLTP),
     &       NCOCCLS_ACT,int_mb(KCIOCCLS_ACT),
     &       int_mb(KCLBT),int_mb(KCLEBT),int_mb(KCLBLK),int_mb(KCI1BT),
     &       int_mb(KCIBT),
     &       int_mb(KCNOCCLS_BAT),int_mb(KCIBOCCLS_BAT),ILTEST)
*
        IBATCHCOUNT = NBATCHES_FEW_BLOCKS(NCBATCH,int_mb(KCLBT))
        NDISTBATCH = NCBATCH - IBATCHCOUNT
        LBPNODE = INT(NCBATCH/NNOD)
!        IF (LBPNODE.GE.IOPTBATCHPNODE) EXIT
c       END IF
*
* VOG Nov16: Check if you will split further the CAS after you have the final 
*            number of TTSS Blocks
*
       LBPNODE = INT(NCBATCH/NNOD)
       if (NTEST.GE.10.AND.ga_nodeid().eq.0)
     &     write(*,*) 'NCBATCH-final', NCBATCH
       if (NTEST.GE.10.AND.ga_nodeid().eq.0)
     &     write(*,*) 'NCBLOCK-final', NCBLOCK
       if (NTEST.GE.10.AND.ga_nodeid().eq.0)
     &      write(*,*) 'LBPNODE',LBPNODE,IOPTBATCHPNODE
       nfragagain = 0
        IF (LBPNODE.LT.IOPTBATCHPNODE) nfragagain = 1
*
*. Calculate blocks/nnodes ratio and print size of batches
*
c       call CI_batches_info(NCBATCH,int_mb(KCLBT),int_mb(KCIBT),NCBLOCK,
c     &                      nnod,ntest,nfragagain)
       call CI_batches_info_OLD(int_mb(KCIBT),NCBLOCK,nnod,ntest,
     &                          nfragagain,NGAS,NPARGAS)
       if (ntest.ge.100. and.ga_nodeid().eq.0)
     &    write(*,'(a,3i10)') 'nblocks, nnod, nfragagain',
     &                         NCBLOCK, NNOD, nfragagain
*
*. Flush memory
*
       call memman(idum,idum,'FLUSM ',idum,'OPTIMUM')
*
*. Condition if I should fragment further the GAS spaces and
*. recalculate the number of blocks
*
!      if (NCBLOCK.lt.NNOD) nfragagain = 1
!      if ((NCBLOCK.lt.NNOD).or.(nfragagain.eq.1)) then
!       nfragagain = 0 ! DANGER
      if (nfragagain.eq.1) then
       if (ntest.eq.10.and.ga_nodeid().eq.0)
     &   write(*,*) 'I will fragment the GAS space'
       !
       ! Check if you have a GAS with MOs in more than one IRREP
       !
       call fraggas(igasloop,ngasadd,imofromto)
      else
       exit
      end if
*
*. End Large Loop
*
      end do
*
*. Print Output
*
12345 continue
      if(ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ***********************************'
      WRITE(*,*) '  Updated GAS Spaces for CI batches '
      WRITE(*,*) ' ***********************************'
      WRITE(*,*)
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      WRITE(6,*)
      IF(IGSFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Gas space provided by courtesy of LUCIA:',  IGSFILL
      IF(ISECFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Secondary space provided by courtesy of LUCIA '
*
       WRITE(*,'(A,I3)')
     &  ' Number of Occupation spaces: ',NCISPC
       WRITE(*,*)
       DO ICISPC = 1, NCISPC
        WRITE(*,'(A,I3)')
     &  ' Bounds on accumulated occupations for space: ',ICISPC
        WRITE(*,'(A)')
     &   ' ====================================================== '
        WRITE(*,'(A)')
        WRITE(*,'(A)') '         Min. occ    Max. occ '
        WRITE(*,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(*,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
        END DO
       END DO
       WRITE(*,*)
       WRITE(*,*)  'Total number of TTSS BLOCKS :', NCBLOCK
       WRITE(*,*)  'Total number of TTSS BATCHES:', NCBATCH
       WRITE(*,*)
       WRITE(*,*) ' Exiting optimum_gas_for_blocks '
       WRITE(*,*)
      end if
*
      return
      end
*======================================================================*
*======================================================================*
      SUBROUTINE ONE_IRREP_PER_GAS
*======================================================================*
*
*.  Case that the user performs a CASSCF calculation and he/she has
*.  requested only one GAS space:
*.  LUCIA-at-NWCHEM now offers a completely free distribution of the MOs
*.  belonging to the same irrep to individual GAS spaces!
*
*.  Konstantinos D. Vogiatzis, May 2016, Minneapolis
*
*======================================================================*
*
      INCLUDE 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cgas.inc'
#include "global.fh"
*
*. Scratch
      integer NGSSHBU(NIRREP)
*
      NTEST = 100
*
*. Distribute MOs of the same IRREP to individual GAS spaces
*
      IF (NIRREP.EQ.8) THEN
       ILIMIT = 1
      ELSE
       ILIMIT = 0
      END IF
      ILIMIT = 0
*
      DO IRREP = 1, NIRREP-1
       IF (NGSSH(IRREP,NGAS).GT.ILIMIT) THEN 
        NGSSHBU = 0
        ISTRIRREP = IRREP+1
        DO JRREP = ISTRIRREP, NIRREP
         NGSSHBU(JRREP) = NGSSH(JRREP,NGAS)
         NGSSH(JRREP,NGAS) = 0
        END DO
        NGAS = NGAS + 1
        DO JRREP = 1, NIRREP
         IF (JRREP.LE.IRREP) THEN
          NGSSH(JRREP,NGAS) = 0
         ELSE
          NGSSH(JRREP,NGAS) = NGSSHBU(JRREP)
         END IF
        END DO
       END IF
      END DO
*
*. Modify Min/Max occ numbers
* 
      icispc = 1 ! Considering only the 1st CI Space for MCSCF
      IGSOCCX(NGAS,1,ICISPC) =  IGSOCCX(1,1,ICISPC) 
      IGSOCCX(NGAS,2,ICISPC) =  IGSOCCX(1,2,ICISPC)
      DO IGAS = 1, NGAS-1
       IGSOCCX(IGAS,1,ICISPC) = 0
       IGSOCCX(IGAS,2,ICISPC) = IGSOCCX(NGAS,2,ICISPC)
      END DO
*
*. Print Output
*
      if(ntest.ge.10.and.ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ****************************************'
      WRITE(*,*) '  Each GAS is occupied by only one IRREP '
      WRITE(*,*) ' ****************************************'
      WRITE(*,*)
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      WRITE(6,*)
*
      WRITE(*,'(A,I3)')
     &  ' Number of Occupation spaces: ',NCISPC
       WRITE(*,*)
       DO ICISPC = 1, NCISPC
        WRITE(*,'(A,I3)')
     &  ' Bounds on accumulated occupations for space: ',ICISPC
        WRITE(*,'(A)')
     &   ' ====================================================== '
        WRITE(*,'(A)')
        WRITE(*,'(A)') '         Min. occ    Max. occ '
        WRITE(*,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(*,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
        END DO
       END DO
       WRITE(*,'(A)')
      end if
*
      RETURN
      END
*======================================================================*
*======================================================================*
      SUBROUTINE BREAK_GAS()
*======================================================================*
*
*. Generate an additional GAS space by moving one MO from the GAS with
*. the most MOs for cases that I need more GAS spaces than the number
*. of IRREPs 
*. eg. when I need more than 8 GAS for a D2h case
*
*. Konstantinos D. Vogiatzis, Nov 2016, while sitting in a
*.                                      biochemistry lab at UCLA
* 
*======================================================================*
*
      INCLUDE 'implicit.inc'
      include 'mxpdim.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cgas.inc'
#include "global.fh"
*       
*. Scratch 
      integer NGSSHBU(NIRREP)
*
*. Initializations
*
      NTEST = 10
      IBREAKGAS = 0
*
*. Find GAS with the most MOs
*
      IBREAKGAS = 0
      MAXMO = 0
      DO IGAS = 1, NGAS
       ISUM = 0
       DO IRREP = 1, NIRREP
        ISUM = ISUM + NGSSH(IRREP,IGAS)
        IF (ISUM.GT.MAXMO) THEN
         IBREAKGAS = IGAS 
         MAXMO = ISUM
        END IF
       END DO
      END DO 
      IF (IBREAKGAS.EQ.0) 
     & STOP 'Empty rooms, where we learn to live without love...'
      IF (MAXMO.EQ.1)
     & STOP 'NOWHERE TO GO, NOWHERE TO HIDE !!!'
      if(ntest.ge.10.and.ga_nodeid().eq.0)
     & WRITE(*,*) ' *** I will break GAS', IBREAKGAS, 'in two GA spaces'
*
*. Move 1 MO from the GAS with the most MOs
*
      NGSSHBU = 0
      DO IRREP = 1, NIRREP
       IF (NGSSH(IRREP,IBREAKGAS).NE.0) IBREAKIRREP = IRREP
      END DO
      if(ntest.ge.10.and.ga_nodeid().eq.0)
     & WRITE(*,*) ' *** I will remove one MO from IRREP ', IBREAKIRREP
*
      NGAS = NGAS + 1
      NGSSH(IBREAKIRREP,IBREAKGAS) = NGSSH(IBREAKIRREP,IBREAKGAS) - 1
      NGSSH(IBREAKIRREP,NGAS) = 1
*
*. Modify Min/Max occ numbers
* 
      ICISPC = 1 ! Considering only the 1st CI Space for MCSCF
      IGSOCCX(NGAS,1,ICISPC)   = IGSOCCX(NGAS-1,1,ICISPC)
      IGSOCCX(NGAS,2,ICISPC)   = IGSOCCX(NGAS-1,1,ICISPC)
      IGSOCCX(NGAS-1,1,ICISPC) = 0
*
*. Print Output
*
      if(ntest.ge.10.and.ga_nodeid().eq.0) then
      WRITE(*,*)
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*) '  Updated GAS Spaces for CI batches from BREAK_GAS '
      WRITE(*,*) ' ****************************************************'
      WRITE(*,*)
      WRITE(*,'(A)') ' Orbital subspaces:'
      WRITE(*,'(A)') ' ================== '
      WRITE(*,*)
      WRITE(*,'(1H ,A,10I4,A)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(*,'(1H ,A,2X,10A,A)')
     &  '                ===== ',('====',I = 1,NIRREP)
      WRITE(*,'(A,10I4)')
     &  '        Inactive       ',
     &   ( NINASH(IRREP),IRREP = 1, NIRREP)

      DO IGAS = 1, NGAS
        WRITE(*,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      WRITE(*,'(A,10I4)')
     &  '        Secondary      ',
     &   ( NSECSH(IRREP),IRREP = 1, NIRREP)
      WRITE(6,*)
      IF(IGSFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Gas space provided by courtesy of LUCIA:',  IGSFILL
      IF(ISECFILL.NE.0) WRITE(6,'(7X,A,I3)')
     &' Secondary space provided by courtesy of LUCIA '
*
       WRITE(*,'(A,I3)')
     &  ' Number of Occupation spaces: ',NCISPC
       WRITE(*,*)
       DO ICISPC = 1, NCISPC
        WRITE(*,'(A,I3)')
     &  ' Bounds on accumulated occupations for space: ',ICISPC
        WRITE(*,'(A)')
     &   ' ====================================================== '
        WRITE(*,'(A)')
        WRITE(*,'(A)') '         Min. occ    Max. occ '
        WRITE(*,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(*,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
        END DO
       END DO
       WRITE(*,*)
      end if
*
      RETURN
      END
*======================================================================*
*======================================================================*
      SUBROUTINE WRITE_OCCLS_AKBK(IOCCLS,NOCCLS,NGAS)
*======================================================================*
*
*. Print the NOCCLS classes IOCCLS
*
*. Jeppe Olsen, Jan. 2012
*. K Vogiatzis, Nov. 2015, Berkeley
*
*======================================================================*
*
      INCLUDE 'implicit.inc'
#include "global.fh"
*
      INTEGER IOCCLS(NGAS,NOCCLS)
*
      if(ga_nodeid().eq.0) then
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
      end if
*
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE CI_batches_info_OLD(KCIBT,NCBLOCK,nnod,ntest,
     &                               nfragagain,NGAS,NPARGAS)
*======================================================================*
*
*. Criterion that was used successfully for the CAS(24,24) - Feb 2016
*. Re-introduced on Feb. 2017.
*
*. Konstantinos Vogiatzis, Feb. 2016, UMN
*.                         Feb. 2017, UTK
*
*======================================================================*
      INCLUDE 'implicit.inc'
#include "global.fh"
*
*. input
      integer NCBLOCK, nnod, ntest
      integer KCIBT(8,ncblock)
*
*. output
      integer nfragagain
*
*. Get information of CI blocks and print them (if requested)
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*) ' =================='
       write(*,*) '  CI_batches_info  '
      end if
      if (ga_nodeid().eq.0.and.ntest.ge.1000) then
       write(*,*)
       write(*,*) ' =================='
       write(*,*) '  BLK         SIZE '
       write(*,*) ' =================='
      end if
*
      isum = 0
      nmaxblk = 0
      do iblock = 1, ncblock
       if (ga_nodeid().eq.0.and.ntest.ge.1000)
     &     write(*,'(i8,i10)') iblock, KCIBT(8,iblock)
       isum = isum + KCIBT(8,iblock)
       if (KCIBT(8,iblock).ge.nmaxblk) nmaxblk =  KCIBT(8,iblock)
      end do
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*) ' =================='
       write(*,'(a5,i14)') '  SUM', isum
       write(*,'(a5,i14)') '  MAX', nmaxblk
       write(*,*) ' =================='
      end if
*
*. Calculate blocks/nnod ratio, criterion if I will fragment further
*. the GAS space
*
      nblkratio = int(isum/nnod)
      if (ga_nodeid().eq.0.and.ntest.ge.10)
     &    write(*,'(a,i20,2i9)') 'ISUM, NNOD, NBLKRATIO',
     &                            isum, nnod, nblkratio
      do iblock = 1, NCBLOCK
       if (KCIBT(8,iblock).gt.nblkratio) nfragagain = 1
      end do
      !
      ! Case that user has requested a MAX number of GAS spaces
      IF (NGAS.LT.NPARGAS) nfragagain = 1
*
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE CI_batches_info(NBATS,LBATS,KCIBT,NCBLOCK,
     &                           iconspa,iconspb,noctpa,noctpb,
     &                           nnod,ntest,
     &                           nfragagain)
*======================================================================*
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'stinf.inc'
#include "global.fh"
*
*. input
      integer NCBLOCK, NBATS, nnod, ntest
      integer LBATS(NBATS),KCIBT(8,ncblock)
      integer noctpa,noctpb
      integer iconspa(noctpa,noctpa),iconspb(noctpb,noctpb)
*
*. output
      integer nfragagain
*
*. Scratch
      integer NSDBATCH(NBATS)
*
*. Obtain number of occupied GAS spaces per supergroup
*
*. Number of supergroups
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
c      NGASSPGRPA = 0
c      NGASSPGRPB = 0
c      !
c      ! Loop over alpha supergroups
c      !
c      DO ISPGRP = 1, NOCTPA
c       DO IGAS = 1, NGAS
c        IF (NELFSPGP(IGAS,ISPGRP).NE.0)
c     &      NGASSPGRPA(ISPGRP) = NGASSPGRPA(ISPGRP) + 1
c       END DO
c      END DO
c      !
c      ! Loop over  beta supergroups
c      !
c      DO ISPGRP = 1, NOCTPB
c       DO IGAS = 1, NGAS
c        IF (NELFSPGP(IGAS,ISPGRP).NE.0)
c     &      NGASSPGRPB(ISPGRP) = NGASSPGRPB(ISPGRP) + 1
c       END DO
c      END DO
*
*. Get information of CI blocks and print them (if requested)
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*)
       write(*,'(a2,a)') '  ', repeat ('=', 56)
       write(*,*) 
     &  '  BLK     BATCH       SIZE    alpha SPGRP    beta SPGRP ',
     &  '  alpha occGAS    beta occGAS'
       write(*,'(a2,a)') '  ', repeat ('=', 56)
      end if
*
      isum = 0
      nmaxblk = 0
      IBA = 1
      ILCHK = LBATS(IBA)
      do iblock = 1, ncblock
       IA = KCIBT(1,iblock)
       IB = KCIBT(2,iblock)
       IAS = KCIBT(3,iblock)
       IBS = KCIBT(4,iblock)
       INTERACT = 0
       IWORK = 0
       do jblock = 1, ncblock 
         JA = KCIBT(1,jblock)
         JB = KCIBT(2,jblock)
         JAS = KCIBT(3,jblock)
         JBS = KCIBT(4,jblock)
         IF((ICONSPA(IA,JA).LE.1.AND.
     &        ICONSPB(IB,JB).LE.1     )   .OR.
     &       (ICONSPA(IA,JA).EQ.MXEXC.AND.
     &        IB.EQ.JB.AND.IBS.EQ.JBS) .OR.
     &       (ICONSPB(IB,JB).EQ.MXEXC.AND.
     &        IA.EQ.JA.AND.IAS.EQ.JAS)     ) then
                INTERACT = INTERACT+1
                IWORK = IWORK+KCIBT(8,iblock)*KCIBT(8,jblock)
         ENDIF
       enddo
       if (ga_nodeid().eq.0.and.ntest.ge.10) then
           write(*,'(i6,i10,i11,a,5i2,a,5i2,a)') iblock, IBA,
     &                                       KCIBT(8,iblock),
     &                        '  | ', (NELFSPGP(IGAS,IA), igas = 1, 5),
     &                        '  | ', (NELFSPGP(IGAS,IB), igas = 1, 5)
           write(*,'(a,i8,a, i20)') '     connections ',interact,
     &                              '  work i*j sum ', iwork
       endif
c    &     write(*,'(i8,2i10,i15)') iblock, KCIBT(8,iblock),
c    &                              KCIBT(1,iblock), KCIBT(2,iblock)
       isum = isum + KCIBT(8,iblock)
       if (KCIBT(8,iblock).ge.nmaxblk) nmaxblk =  KCIBT(8,iblock)
*
       IF (iblock.EQ.ILCHK) THEN
        IBA = IBA + 1
        ILCHK = ILCHK + LBATS(IBA)
       END IF
*
      end do
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,'(a2,a)') '  ', repeat ('=', 56)
       write(*,'(a5,i13)') '  SUM', isum
       write(*,'(a5,i13)') '  MAX', nmaxblk
       write(*,'(a2,a)') '  ', repeat ('=', 56)
      end if
*
*. New route based on the TTSS batches, not the TTSS blocks
*
      NSDBATCH = 0 ! Number of SDs in a given TTSS Batch
      if (ga_nodeid().eq.0.AND.ntest.ge.100) THEN
       write(*,*)
       write(*,*) 'Info for each TTSS batch'
       write(*,*)
       write(*,*) ' BATCH   BLOCKS       LENGTH'
      end if
*
*. Obtain size of each TTSS BATCH (number of SDs)
*
      II = 0
      DO I = 1, NBATS
       DO IBLK = 1, LBATS(I)
        II = II + 1
        NSDBATCH(I) = NSDBATCH(I) + KCIBT(8,II)
       END DO
       if (ga_nodeid().eq.0.and.ntest.ge.100)
     &     write(*,'(i6,i10,i20)') I,LBATS(I), NSDBATCH(I)
      END DO
*
*. Calculate blocks/nnod ratio, criterion if I will fragment further
*. the GAS space
*
      NFRAG = 1
      if (ga_nodeid().eq.0.and.ntest.ge.10) write(*,*) 'NFRAG', NFRAG
      nblkratio = int(isum/nnod)
      nblkratio = int(nblkratio/NFRAG) !TEST - best so far for small CAS: 8 or 16
                                       !       I keep it to 2 for testing
                                       !       large CAS > (20,20)
      nfragagain = 0
      if (ga_nodeid().eq.0.and.ntest.ge.10)
     & write(*,'(a,3i12)') 'isum, nnod, nblkratio', isum,nnod,nblkratio
*
      DO I = 1, NBATS
       IF (NSDBATCH(I).GT.nblkratio) nfragagain = 1
      END DO
*
cVOG  do iblock = 1, ncblock
cVOG   if (KCIBT(8,iblock).gt.nblkratio) nfragagain = 1
cVOG  end do
*
*. Check if number of TTSS batches smaller than number of nodes
*
      if (NBATS.lt.nnod) nfragagain = 1
*
*. Check if at least one GAS space has more than 1 MO
*
      lallsingle = 1
      do igas = 1, ngas
       nmoings = 0
       do irrep = 1, nirrep
        nmoings = nmoings + ngssh(irrep,igas)
       end do
       if (nmoings.gt.1) lallsingle = 0
       if (ga_nodeid().eq.0.and.ntest.ge.100)
     &   write(*,*) 'in the loop', igas, nmoings, lallsingle
      end do
*
      if (lallsingle.eq.1) nfragagain = 0
*
       if (ga_nodeid().eq.0.and.ntest.ge.10)
     &  write(*,*) 'nfragagain',nfragagain
*
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE distribute_CIblks(ILEN,NBATS,LBATS,IBATS,NBLK,nnod,
     &                             IMAP2)
*======================================================================*
*
*. The different levels of distribution of the CI blocks considered
*. (so far)
*
*. Zeroth-order: Based on the number of determinants
*.    1st-order: Based on the weight of each TTSS
*.               TTSS Weight = (# of connections) x (size of TTSS)
*.    2nd-order:
*
*
*. Variable that holds the batch/block sizes per node:
*  ibpnod(*,1): Number of TTSS blocks  per node
*  ibpnod(*,2): Number of TTSS batches per node
*
*. Konstantinos D. Vogiatzis, Nov. 2015, Minneapolis, after Bert
*.                                                    thretening my life
*
*======================================================================*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'strbas.inc'
#include "mafdecls.fh"
#include "global.fh"
*
*. input
      integer ILEN,NBATS,NBLK,nnod
      integer LBATS(NBATS),IBATS(8,NBLK)
*
*. input/output
c      integer IMAP2(nnod)
      integer IMAP2(NBATS)
*
*. Scratch
      integer ntotweight(nnod), nttssw(NBATS)
      integer iblksum(2), isize(nnod), ibpnod(nnod,2)
*
*. Initializations
*
      ntest = 00
      do i=1, NBATS
       imap2(i) = 0
      end do
      do inode = 1, nnod
c       imap2(inode) = 0
       isize(inode) = 0
      end do 
*
*. New route, based on the weight of each TTSS:
*
*. TTSS Weight = (# of connections) x (size of TTSS)
*
*
*. Initialize SPGRPCON subroutine
*
*
* Info for this internal space
*. type of alpha and beta strings - as H does not change
*. the number of electrons, I do not distinguish between spaces for C
*and S
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
*. connection matrices for supergroups
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DISTCI')
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*) ' ILEN    ', ILEN
       write(*,*) ' IOCTPA  ', IOCTPA
       write(*,*) ' NOCTPA  ', NOCTPA
       write(*,*) ' NOCTPB  ', NOCTPB
       write(*,*) ' NGAS    ', NGAS
       write(*,*) ' NBLK    ', NBLK
       write(*,*) ' NBATS   ', NBATS
      end if
*
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPA),ntest,ga_nodeid())
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              int_mb(KCONSPB),ntest,ga_nodeid())
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
         WRITE(6,*) ' Supergroup connection matrix '
         CALL IWRTMA(int_mb(KCONSPA),NOCTPA,NOCTPA,NOCTPA,NOCTPA)
         WRITE(6,*) ' Supergroup connection matrix '
         CALL IWRTMA(int_mb(KCONSPB),NOCTPB,NOCTPB,NOCTPB,NOCTPB)
      END IF
*
*. Only for PRINTING reasons
*
       call ga_sync() !DEC
       nnod = ga_nnodes()
       nfragagain = 0
       call CI_batches_info(NBATS,LBATS,IBATS,NBLK,
     &                      int_mb(KCONSPA),int_mb(KCONSPB),
     &                      NOCTPA,NOCTPB,
     &                      nnod,ntest,nfragagain)
*
*. For alpha beta strings
*
       IATP = 1
       IBTP = 2
*
*. Calculate the TTSS weights
*
      LFIRSTORDER = 0
      if (NTEST.GE.100.and.ga_nodeid().eq.0)
     &             write(*,*) ' Calculate the TTSS weights'
!      IF (LFIRSTORDER.EQ.0) THEN
       call ttssweight0(nttssw,nttsswsum,IBATS,NBATS,LBATS,NBLK)
!      ELSE
c       call ttssweight(nttssw,nttsswsum,
c     &                 int_mb(KCONSPA),int_mb(KCONSPB),
c     &                 NOCTPA,NOCTPB,IBATS,NBATS,LBATS,NBLK,NSMST,
c     &                 int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)))
!      END IF
cVOG  call ttssweight15(ttssw,int_mb(KCONSPA),int_mb(KCONSPB),
cVOG &                NOCTPA,NOCTPB,IBATS,NBLK)
cVOG  call ttssweight2(ttssw,NOCTPA,NOCTPB,IBATS,NBATS,LBATS,NBLK,
cVOG &                 int_mb(KNSTSO(IATP)),int_mb(KNSTSO(IBTP)))
cVOG  call ttssweight3(ttssw,IBATS,NBATS,LBATS,NBLK)
*
*.  Distribute TTSS Blocks in such a way that each proc has equal TTSS
*.  weight
*
*. Loop that varies the sum of the TTSS block weights per processor.
*. Important for cases that TTSS blocks are not distributed in the available nodes
*
      IF (LFIRSTORDER.EQ.0) THEN
      NPLUSONE = 0
      MINDIFF = nttsswsum
      ITHRESH = INT(2 * nttsswsum / NNOD)
110   CONTINUE
*
      !
      ! Initialize
      !
      iaccumsd = 0
      inodrem = nnod
      imap2(1) = 1
      ibpnod = 0
      isumsd = 0
      IIBLK = 0
      JBATS = 0
cUT
      ISUMW    = 0
      JBLK     = 0
      IACCUSD  = 0
      ISUMSD   = 0
      INODE    = 0
      NPLUSONE = NPLUSONE + 1
      IREMW    = nttsswsum
      INODREM  = NNOD
cMINMAX      LIMBLKCHK = limitblk
*
!      remw= 100.0d0 ! size of the remaining TTSS weights
*
*. NBPNMIN is the minimum number of TTSS batches that each processor will get
*
      NBPNMIN = INT(NBATS/NNOD)
      NREM    = MOD(NBATS,NNOD)
*
      limitblk = INT(nttsswsum / NNOD) !remw/NBATS * NBPNMIN
      limitblk = limitblk + NPLUSONE
*
      if (ga_nodeid().eq.0.and.ntest.ge.100) then
       write(*,*) 'NEW ALGORITHM FOR THE CI DISTRIBUTION'
       write(*,*) 'NBATS         ', NBATS
       write(*,*) 'limitblk      ', limitblk
       write(*,*) 'MINDIFF       ', MINDIFF
       write(*,*)
      end if
cUT      nwindow  = INT(limitblk/NFRAG)
cUT      limitblk = limitblk + nwindow
      !
cMINMAX      if (ga_nodeid().eq.0.and.ntest.ge.100) then
cMINMAX       write(*,*) 'NEW ALGORITHM FOR THE CI DISTRIBUTION'
cMINMAX       write(*,*) 'NBATS         ', NBATS
cMINMAX       write(*,*) 'limitblk      ', limitblk
cMINMAX       write(*,*) 'MINDIFF       ', MINDIFF
cMINMAX       write(*,*)
cMINMAX      end if
cMINMAX      !
cMINMAX      ! Begin work
cMINMAX      !
cMINMAX      DO JBATS = 1, NBATS
cMINMAX       !
cMINMAX       ! Add the TTSS weights
cMINMAX       ISUMW   = ISUMW   + NTTSSW(JBATS)
cMINMAX       ICHKBLK = ISUMW   + NTTSSW(JBATS + 1)
cMINMAX       !
cMINMAX       ! Count the number of CI batches distributed per processor
cMINMAX       IBPNOD(INODE+1,2) = IBPNOD(INODE+1,2) + 1
cMINMAX*
cMINMAX       DO IBLKS = 1, LBATS(JBATS)
cMINMAX        !
cMINMAX        ! Count the number of CI blocks distributed per processor
cMINMAX        IBPNOD(INODE+1,1) = IBPNOD(INODE+1,1) + 1
cMINMAX        !
cMINMAX        JBLK = JBLK + 1
cMINMAX        IACCUSD = IACCUSD + IBATS(8,JBLK)
cMINMAX        ISUMSD  = ISUMSD  + IBATS(8,JBLK)
cMINMAX       END DO ! End of the blocks per batch loop
cMINMAX*
cMINMAX       !
cMINMAX       ! Control
cMINMAX       IF (INODE.LT.(NNOD-1)) THEN
cMINMAX        IF (ICHKBLK.GT.LIMBLKCHK) THEN
cMINMAX         INODE = INODE + 1
cMINMAX         ISIZE(INODE) = ISUMSD
cMINMAX         IMAP2(INODE+1) = IACCUSD + 1
cMINMAX         !
cMINMAX         ! Dynamically modify the limit size of the block
cMINMAX         IREMW    = IREMW - ISUMW
cMINMAX         INODREM  = INODREM - 1
cMINMAX         LIMBLKCHK = INT(IREMW/INODREM)
cMINMAX         LIMBLKCHK = LIMBLKCHK + NPLUSONE
cMINMAX         !
cMINMAX         ! Nullify sums for next node
cMINMAX         ISUMW  = 0
cMINMAX         ISUMSD = 0
cMINMAX        END IF
cMINMAX       ELSE
cMINMAX        IF (JBLK.EQ.NBLK) ISIZE(INODE+1) = ISUMSD
cMINMAX       END IF
cMINMAX      END DO ! End of the distribution loop
      DO INODE = 1, NNOD
       !
       ! Sum TTSS weights
!       nchkttssw = 0
       nsumttssw = 0
       !
       DO I = 1, NBPNMIN
        JBATS = JBATS + 1
        nsumttssw = nsumttssw + nttssw(JBATS)
        !
        DO IBLK = 1, LBATS(JBATS)
         IIBLK = IIBLK + 1
         iaccumsd = iaccumsd + IBATS(8,IIBLK)
         isumsd   = isumsd   + IBATS(8,IIBLK)
        END DO
       END DO
       !
       ! Choose route
       IF (nsumttssw.lt.limitblk.and.NREM.NE.0) THEN
        ibpnod(INODE,2) = NBPNMIN+1
        NREM = NREM - 1
        JBATS = JBATS + 1
        nsumttssw = nsumttssw + nttssw(JBATS)
        DO IBLK = 1, LBATS(JBATS)
         IIBLK = IIBLK + 1
         iaccumsd = iaccumsd + IBATS(8,IIBLK)
         isumsd   = isumsd   + IBATS(8,IIBLK)
        END DO
       ELSE
        ibpnod(INODE,2) = NBPNMIN
       END IF
       !
       ! Add SD dimensions
       isize(INODE) = isumsd
       IF (INODE.NE.NNOD) imap2(inode+1) = iaccumsd + 1
*
       isumsd  = 0
*
      END DO
*
*. Sum the TTSS BLOCKS and weights per node
*
      j = 0
      ntotweight = 0
      do inode = 1, ga_nnodes()
       do i = 1, ibpnod(inode,2)
        j = j + 1
        ntotweight(inode) = ntotweight(inode) + nttssw(j)
        ibpnod(inode,1)   = ibpnod(inode,1) + LBATS(j)
       end do
      end do
*
*. Remaining IMAP2 enrties should hold the total number of SDs
*
      do i = nnod+1,NBATS
       imap2(i) = ilen
      end do
*
*. Minimization of the MIN/MAX difference for an optimum distribution
*
      !
      ! CI Block Weight statistics
      !
      minw = ntotweight(1)
      maxw = 0
      DO INODE = 1, ga_nnodes()
       if (ntotweight(inode).ge.maxw) maxw = ntotweight(inode)
       if (ntotweight(inode).lt.minw) minw = ntotweight(inode)
      END DO
      ndiffw = maxw - minw
      !
      ! Control
      !
      IF (ndiffw.LT.MINDIFF) THEN
       MINDIFF = ndiffw
       IOPTLIMIT = LIMITBLK
       NPLUSFINAL = NPLUSONE
       if (NTEST.GE.100.and.ga_nodeid().eq.0) write(*,*) 'NEW MIN DIFF',
     &                                   MINDIFF,LIMITBLK,NPLUSFINAL
      END IF
*     
*. Set a limit that allows to examine various LIMITBLK
*. and minimize the MIN/MAX difference
*
      IF (NPLUSONE.NE.ITHRESH) GOTO 110 !Arbitrary limit in order to avoid an infinite loop!
*
*. Once the big loop is over and we have the optimum LIMITBLK, repeat work
*
      LIMITBLK = IOPTLIMIT
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*) 'NEW ALGORITHM FOR THE CI DISTRIBUTION - FINAL STEP'
       write(*,*) 'NBATS         ', NBATS
       write(*,*) 'limitblk      ', limitblk
       write(*,*)
      end if   
      !
      ! Initializations
      !
      iaccumsd = 0
      inodrem = nnod
      imap2(1) = 1
      ibpnod = 0
      isumsd = 0
      IIBLK = 0
      JBATS = 0
      IACCUSD  = 0
      ISUMSD   = 0
      INODE    = 0
      IREMW    = nttsswsum
      INODREM  = NNOD
      NBPNMIN = INT(NBATS/NNOD)
      NREM    = MOD(NBATS,NNOD)
      !
      ! And start again
      !
      DO INODE = 1, NNOD
       !
       ! Sum TTSS weights
!       nchkttssw = 0
       nsumttssw = 0
       !
       DO I = 1, NBPNMIN
        JBATS = JBATS + 1
        nsumttssw = nsumttssw + nttssw(JBATS)
        !
        DO IBLK = 1, LBATS(JBATS)
         IIBLK = IIBLK + 1
         iaccumsd = iaccumsd + IBATS(8,IIBLK)
         isumsd   = isumsd   + IBATS(8,IIBLK)
        END DO
       END DO
       !
       ! Choose route
       IF (nsumttssw.lt.limitblk.and.NREM.NE.0) THEN
        ibpnod(INODE,2) = NBPNMIN+1
        NREM = NREM - 1
        JBATS = JBATS + 1
        nsumttssw = nsumttssw + nttssw(JBATS)
        DO IBLK = 1, LBATS(JBATS)
         IIBLK = IIBLK + 1
         iaccumsd = iaccumsd + IBATS(8,IIBLK)
         isumsd   = isumsd   + IBATS(8,IIBLK)
        END DO
       ELSE
        ibpnod(INODE,2) = NBPNMIN
       END IF
       !
       ! Add SD dimensions
       isize(INODE) = isumsd
       IF (INODE.NE.NNOD) imap2(inode+1) = iaccumsd + 1
*
       isumsd  = 0
*
      END DO
*
*. Sum the TTSS BLOCKS and weights per node
*
      j = 0
      ntotweight = 0
      do inode = 1, ga_nnodes()
       do i = 1, ibpnod(inode,2)
        j = j + 1
        ntotweight(inode) = ntotweight(inode) + nttssw(j)
        ibpnod(inode,1)   = ibpnod(inode,1) + LBATS(j)
       end do
      end do
*
*. Remaining IMAP2 enrties should hold the total number of SDs
*
      do i = nnod+1,NBATS
       imap2(i) = ilen
      end do
cMINMAX*     
cMINMAX*. Minimization of the MIN/MAX difference for an optimum distribution
cMINMAX*
cMINMAX      !
cMINMAX      ! CI Block Weight statistics
cMINMAX      !
cMINMAX      minw = ntotweight(1)
cMINMAX      maxw = 0
cMINMAX      DO INODE = 1, ga_nnodes()
cMINMAX       if (ntotweight(inode).ge.maxw) maxw = ntotweight(inode)
cMINMAX       if (ntotweight(inode).lt.minw) minw = ntotweight(inode)
cMINMAX      END DO
cMINMAX      ndiffw = maxw - minw
cMINMAX      !
cMINMAX      ! Control
cMINMAX      !
cMINMAX      IF (ndiffw.LT.MINDIFF) THEN
cMINMAX       MINDIFF = ndiffw
cMINMAX       IOPTLIMIT = LIMITBLK
cMINMAX       NPLUSFINAL = NPLUSONE
cMINMAX       if (NTEST.GE.100.and.ga_nodeid().eq.0) write(*,*) 'NEW MIN DIFF',
cMINMAX     &                                   MINDIFF,LIMITBLK,NPLUSFINAL
cMINMAX      END IF
cMINMAX*     
cMINMAX*. Set a limit that allows to examine various LIMITBLK
cMINMAX*. and minimize the MIN/MAX difference
cMINMAX*
cMINMAX      IF (NPLUSONE.NE.(INT(limitblk/2))) GOTO 110 !Arbitrary limit in order to avoid an infinite loop!
cMINMAX*
cMINMAX*. Once the big loop is over and we have the optimum LIMITBLK, repeat work
cMINMAX*
cMINMAX      LIMITBLK = IOPTLIMIT
cMINMAX      if (ga_nodeid().eq.0.and.ntest.ge.10) then
cMINMAX       write(*,*) 'NEW ALGORITHM FOR THE CI DISTRIBUTION - FINAL STEP'
cMINMAX       write(*,*) 'NBATS         ', NBATS
cMINMAX       write(*,*) 'limitblk      ', limitblk
cMINMAX       write(*,*)
cMINMAX      end if   
cMINMAX      !
cMINMAX      ! Initializations
cMINMAX      !
cMINMAX      iaccumsd = 0
cMINMAX      inodrem = nnod
cMINMAX      imap2(1) = 1
cMINMAX      ibpnod = 0
cMINMAX      isumsd = 0
cMINMAX      IIBLK = 0
cMINMAX      JBATS = 0
cMINMAX      ISUMW    = 0
cMINMAX      JBLK     = 0
cMINMAX      IACCUSD  = 0
cMINMAX      ISUMSD   = 0
cMINMAX      INODE    = 0
cMINMAX      IREMW    = nttsswsum
cMINMAX      INODREM  = NNOD
cMINMAX      LIMBLKCHK = limitblk
cMINMAX      !
cMINMAX      ! Begin work
cMINMAX      !
cMINMAX      DO JBATS = 1, NBATS
cMINMAX       !
cMINMAX       ! Add the TTSS weights
cMINMAX       ISUMW   = ISUMW   + NTTSSW(JBATS)
cMINMAX       ICHKBLK = ISUMW   + NTTSSW(JBATS + 1)
cMINMAX       !
cMINMAX       ! Count the number of CI batches distributed per processor
cMINMAX       IBPNOD(INODE+1,2) = IBPNOD(INODE+1,2) + 1
cMINMAX*
cMINMAX       DO IBLKS = 1, LBATS(JBATS)
cMINMAX        !
cMINMAX        ! Count the number of CI blocks distributed per processor
cMINMAX        IBPNOD(INODE+1,1) = IBPNOD(INODE+1,1) + 1
cMINMAX        !
cMINMAX        JBLK = JBLK + 1
cMINMAX        IACCUSD = IACCUSD + IBATS(8,JBLK)
cMINMAX        ISUMSD  = ISUMSD  + IBATS(8,JBLK)
cMINMAX       END DO ! End of the blocks per batch loop
cMINMAX*
cMINMAX       !
cMINMAX       ! Control
cMINMAX       IF (INODE.LT.(NNOD-1)) THEN
cMINMAX        IF (ICHKBLK.GT.LIMBLKCHK) THEN
cMINMAX         INODE = INODE + 1
cMINMAX         ISIZE(INODE) = ISUMSD
cMINMAX         IMAP2(INODE+1) = IACCUSD + 1
cMINMAX         !
cMINMAX         ! Dynamically modify the limit size of the block
cMINMAX         IREMW    = IREMW - ISUMW
cMINMAX         INODREM  = INODREM - 1
cMINMAX         LIMBLKCHK = INT(IREMW/INODREM)
cMINMAX         LIMBLKCHK = LIMBLKCHK + NPLUSFINAL
cMINMAX         !
cMINMAX         ! Nullify sums for next node
cMINMAX         ISUMW  = 0
cMINMAX         ISUMSD = 0
cMINMAX        END IF
cMINMAX       ELSE
cMINMAX        IF (JBLK.EQ.NBLK) ISIZE(INODE+1) = ISUMSD
cMINMAX       END IF
cMINMAX      END DO ! End of the distribution loop
cMINMAX*
cMINMAX*. Sum the TTSS weights per node
cMINMAX*
cMINMAX      j = 0
cMINMAX      ntotweight = 0
cMINMAX      do inode = 1, ga_nnodes()
cMINMAX       do i = 1, ibpnod(inode,2)
cMINMAX        j = j + 1
cMINMAX        ntotweight(inode) = ntotweight(inode) + nttssw(j)
cMINMAX       end do
cMINMAX      end do
cMINMAX*
cMINMAX*. Remaining IMAP2 enrties should hold the total number of SDs
cMINMAX*
cMINMAX      do i = nnod+1,NBATS
cMINMAX       imap2(i) = ilen
cMINMAX      end do
*
*. Print distribution output and statistics
*
      if (ga_nodeid().eq.0) then
       isum    = 0
       nsuw    = 0
       maxw    = 0
       nmaxblk = 0
       iblksum(2) = 0
*
       minw = ntotweight(1)
*
       write(*,*)
       write(*,'(a2,a)') '  ', repeat ('=', 81)
       write(*,'(a2,a81)') '  NODE           SIZE      ',
     &     'MAPPING   BATCH-PER-NODE  BLOCK-PER-NODE  TOTAL-WEIGHT'
       write(*,'(a2,a)') '  ', repeat ('=', 81)
       do inode = 1, ga_nnodes() !nblk
         write(*,'(i8,i14,i13,i9,i16,i20)') inode, isize(inode),
     &                                  IMAP2(inode), ibpnod(inode,2),
     &                                 ibpnod(inode,1),ntotweight(inode)
         !
         ! CI Block statistics
         !
         isum    = isum    + isize(inode)
         iblksum(1) = iblksum(1) + ibpnod(inode,1)
         iblksum(2) = iblksum(2) + ibpnod(inode,2)
         if (isize(inode).ge.nmaxblk) nmaxblk = isize(inode)
         if (inode.eq.1) then
          nminblk = isize(inode)
         else
          if (isize(inode).lt.nminblk) nminblk = isize(inode)
         end if
         !
         ! CI Block Weight statistics
         !
         nsuw  = nsuw    + ntotweight(inode)
         if (ntotweight(inode).ge.maxw) maxw = ntotweight(inode)
         if (ntotweight(inode).lt.minw) minw = ntotweight(inode)
       end do
       !
       ! Print Statistics
       !
       ndiffblk = nmaxblk - nminblk
       ndiffw   = maxw   - minw
       write(*,'(a2,a)') '  ', repeat ('=', 81)
       write(*,'(a5,i17,i22,i16,i20)') '  SUM',  isum, iblksum(2),
     &                                                 iblksum(1),nsuw
       write(*,'(a5,i17,i58)')         '  MIN',  nminblk,  minw
       write(*,'(a5,i17,i58)')         '  MAX',  nmaxblk,  maxw
       write(*,'(a6,i16,i58)')         '  DIFF', ndiffblk, ndiffw
       write(*,'(a2,a)') '  ', repeat ('=', 81)
       write(*,*)
      end if
*
*. Self-consistent distribution
*
      IF (ndiffw.LT.MINDIFF) MINDIFF = ndiffw
      ELSE
***********
       !
       ! New surboutine for 1st-order distribution
       ! TTSSweight is based on the number of connections
       ! TTSSweight is an integer
       !
       CALL DIST_FIRSTORDER(NTTSSW,IBATS,LBATS,NBLK,NBATS,
     &                      NNOD,IMAP2,ISIZE,ILEN) 
      END IF
***********
*
*. Release memeory
*
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'DISTCI')
*
      RETURN
      END
*======================================================================*
*======================================================================*
      SUBROUTINE DIST_FIRSTORDER(NTTSSW,IBATS,LBATS,NBLK,
     &                           NBATS,NNOD,IMAP2,ISIZE,ILEN)
*======================================================================*
*
*. Distribute the TTSS batches to the nodes by using the first-order
*. TTSS weights (integer)
*
*. Konstantinos D. Vogiatzis, September 2016, Knoxville
*
*TODO: Think to include (1) MAX of the weights 
*                       (2) loop that minimize te DIFF (max-min) 
*======================================================================*
*
      INCLUDE 'wrkspc.inc'
#include "mafdecls.fh"
#include "global.fh"
*
*. Input
      INTEGER NNOD, NBATS, NBLK, ILEN
      INTEGER NTTSSW(NBATS), LBATS(NBATS), IBATS(8,NBLK)
*
*. Output
      INTEGER IMAP2(NBATS), ISIZE(NNOD)
*
*. Scratch
      INTEGER IBLKSUM(2), IBPNOD(NNOD,2), NWEIGHTNOD(NNOD)
*
*. Print Input
*
      NTEST = 100
      IF (ga_nodeid().EQ.0.and.NTEST.GE.1000) then
       WRITE(*,*)
       WRITE(*,*) ' ***********************************************'
       WRITE(*,*) '  Distribution of the TTSS batches to the nodes '
       WRITE(*,*) ' ***********************************************'
       WRITE(*,*)
       WRITE(*,*) '            TTSS Weights per batch'
       WRITE(*,*)
       DO I = 1, NBATS
        WRITE(*,'(2I20)') I, NTTSSW(I)
       END DO 
       WRITE(*,*)
      END IF
*
*. Initializations
*
      !
      ! Dynamic variable that controls the distribution
      NFRAG = 100
      !
      ! Calculate NTOTWEIGHT
      NTOTWEIGHT = 0
!      TOTALF = 0
      DO I = 1, NBATS
       NTOTWEIGHT = NTOTWEIGHT + NTTSSW(I)
      END DO
*
*. Begin dynamic loop
*
120   CONTINUE
      !
      ! Threshold for adding TTSS Batches in a node
      LIMITBLK = INT(NTOTWEIGHT/NNOD)
      LIMITBLK = INT(LIMITBLK + LIMITBLK / NFRAG)
      !
      ! Nullify variables
      ISUMW   = 0
      ISUMSD  = 0
      IACCUSD = 0
      INODE   = 0
      IBPNOD  = 0
      ISIZE   = 0
      JBLK    = 0 
      !
      ! Except IMAP(1) !!!
      IMAP2(1) = 1
      !
      ! Control variables
      IREMW   = NTOTWEIGHT ! length of remaining weights
      INODREM = NNOD       ! number of remaining nodes
*
*. Print intermediate details
*
      IF (ga_nodeid().EQ.0.and.NTEST.GE.100) then
       WRITE(*,*) ' NEW ALGORITHM FOR THE CI DISTRIBUTION'
       WRITE(*,*) ' BASED ON THE 1ST-ORDER TTSS WEIGHT'
       WRITE(*,*)
       WRITE(*,*) ' NBatches     ', NBATS
       WRITE(*,*) ' NNodes       ', NNOD
       WRITE(*,*) ' NFRAG        ', NFRAG
       WRITE(*,*) ' TOTAL WEIGHT ', NTOTWEIGHT
       WRITE(*,*) ' NODE  LIMIT  ', LIMITBLK
       WRITE(*,*)
      END IF
*
*. Begin the distribution of the TTSS Batches
*
      DO JBATS = 1, NBATS
       !
       ! Add the TTSS weights
       ISUMW   = ISUMW   + NTTSSW(JBATS)
       ICHKBLK = ISUMW   + NTTSSW(JBATS + 1)
       !
       ! Count the number of CI batches distributed per processor
       IBPNOD(INODE+1,2) = IBPNOD(INODE+1,2) + 1 
*
       DO IBLKS = 1, LBATS(JBATS)
        !
        ! Count the number of CI blocks distributed per processor
        IBPNOD(INODE+1,1) = IBPNOD(INODE+1,1) + 1
        !
        JBLK = JBLK + 1
        IACCUSD = IACCUSD + IBATS(8,JBLK)
        ISUMSD  = ISUMSD  + IBATS(8,JBLK)
       END DO ! End of the blocks per batch loop
*
       !
       ! Control
       IF (INODE.LT.(NNOD-1)) THEN
        IF (ICHKBLK.GT.LIMITBLK) THEN
         INODE = INODE + 1
         ISIZE(INODE) = ISUMSD
         IMAP2(INODE+1) = IACCUSD + 1
         !
         ! Dynamically modify the limit size of the block
         IREMW    = IREMW - ISUMW
         INODREM  = INODREM - 1
         LIMITBLK = INT(IREMW/INODREM)
         LIMITBLK = INT(LIMITBLK + LIMITBLK / NFRAG)
         !
         ! Nullify sums for next node
         ISUMW  = 0
         ISUMSD = 0
        END IF
       ELSE
        IF (JBLK.EQ.NBLK) ISIZE(INODE+1) = ISUMSD
       END IF
      END DO ! End of the distribution loop
*
!*
!*. Begin dynamic loop
!*
!120   CONTINUE  
!      !
!      ! Nullify variables
!      IIBLK    = 0
!      JBATS    = 0
!      ISUMSD   = 0
!      IACCUMSD = 0
!      IBPNOD   = 0 
!      ! Except IMAP(1) !!!
!      IMAP2(1) = 1
!*
!      NBPNMIN = INT(NBATS/NNOD)
!      NREM    = MOD(NBATS,NNOD)
!*
!      ALIMITBLK = NTOTWEIGHT/NBATS * NBPNMIN
!      ALIMITBLK = ALIMITBLK + ALIMITBLK/NFRAG
!*
!*. Print intermediate details
!*
!      IF (ga_nodeid().EQ.0.and.NTEST.GE.100) then
!       WRITE(*,*) ' NEW ALGORITHM FOR THE CI DISTRIBUTION'
!       WRITE(*,*) ' BASED ON THE 1ST-ORDER TTSS WEIGHT'
!       WRITE(*,*)
!       WRITE(*,*) ' NBatches     ', NBATS
!       WRITE(*,*) ' NNodes       ', NNOD
!       WRITE(*,*) ' NBatchPerNode', NBPNMIN
!       WRITE(*,*) ' NREM(modulo) ', NREM
!       WRITE(*,*) ' NFRAG        ', NFRAG
!       WRITE(*,*) ' TOTAL WEIGHT ', NTOTWEIGHT
!       WRITE(*,*) ' BLOCK LIMIT  ', ALIMITBLK
!       WRITE(*,*)
!      END IF
!*
!*. Begin work
!*
!      DO INODE = 1, NNOD
!       !
!       ! Sum TTSS weights
!       !
!       NCHKTTSSW = 0
!       NSUMTTSSW = 0
!       !
!       DO I = 1, NBPNMIN
!        JBATS = JBATS + 1
!        NSUMTTSSW = NSUMTTSSW + NTTSSW(JBATS)
!        !
!        ! Add size of TTSS blocks for JBATS batch
!        DO IBLK = 1, LBATS(JBATS)
!         IIBLK = IIBLK + 1
!         IACCUMSD = IACCUMSD + IBATS(8,IIBLK)
!         ISUMSD   = ISUMSD   + IBATS(8,IIBLK)
!        END DO
!       END DO
!       !
!       ! Number of TTSS Batches per node
!       IBPNOD(INODE,2) = NBPNMIN
!       ! 
!       ! If the TTSS weight is lower than the threshold, 
!       ! add an additional TTSS Batch
!       !
!       IF (NSUMTTSSW.LT.ALIMITBLK.AND.NREM.NE.0) THEN
!        IBPNOD(INODE,2) = IBPNOD(INODE,2) + 1
!        NREM = NREM - 1
!        JBATS = JBATS + 1
!        NSUMTTSSW = NSUMTTSSW + NTTSSW(JBATS)
!        !
!        ! Add size of TTSS blocks for JBATS batch
!        DO IBLK = 1, LBATS(JBATS)
!         IIBLK = IIBLK + 1
!         IACCUMSD = IACCUMSD + IBATS(8,IIBLK)
!         ISUMSD   = ISUMSD   + IBATS(8,IIBLK)
!        END DO
!       END IF
!       !
!       ! Add SD dimensions
!       ISIZE(INODE) = ISUMSD
!       IF (INODE.NE.NNOD) IMAP2(inode+1) = IACCUMSD + 1
!*
!       ISUMSD = 0
!*
!      END DO
*
*. Check if the distribution has been done correclty.
*. If not, reduce NFRAG and repeat the process
*
      if (ga_nodeid().EQ.0) write(*,*) 'NREM,INODE,NNOD',NREM,INODE,NNOD
      IF ((NREM.NE.0) .OR. (INODE.NE.(NNOD-1))) THEN
       IF (ga_nodeid().EQ.0)
     &  WRITE(*,*) ' PROBLEM, TTSS BATCHES ARE MISSING!!!'
       NFRAG = NFRAG - 1
       IF (NFRAG.NE.0) THEN
        GOTO 120
       ELSE
        STOP       ' NFRAG = 0. Something BAD has happened here!!! '
       END IF
      END IF
*
*. *** END OF LOOP ***
*
*
*. Update variables
*
      NWEIGHTNOD = 0
      ICHECK = 0
      II = 1
      !
      ! Sum the TTSS weights per CI batch
      DO I = 1, NBATS
!      IF (ga_nodeid().EQ.0) write(*,*) 'NTTSSW', I, NTTSSW(I), NTTSSW(1)
       NWEIGHTNOD(II) = NWEIGHTNOD(II) + NTTSSW(I)
!      IF (ga_nodeid().EQ.0) write(*,*) 'NWEIGHTNOD',I, NWEIGHTNOD(II)
       ICHECK = ICHECK + 1
       IF (ICHECK.EQ.IBPNOD(II,2)) THEN
        II = II + 1
        ICHECK = 0
       END IF
      END DO
!      !
!      ! Calculate the number of TTSS blocks per node
!      J = 0
!      DO INODE = 1, NNOD
!       DO I = 1, IBPNOD(INODE,2)
!        J = J + 1
!        IBPNOD(INODE,1) = IBPNOD(INODE,1) + LBATS(J)
!       END DO
!      END DO
*
*. Remaining IMAP2 enrties should hold the total number of SDs
*
      DO I = NNOD+1,NBATS
       IMAP2(I) = ILEN
      END DO
*
*. PRINT OUTPUT
*
      IF (ga_nodeid().EQ.0) THEN
       !
       ! Prepare variables for statistics
       !
       ISUM    = 0
       NSUW    = 0
       MAXW    = 0
       NMAXBLK = 0
       IBLKSUM = 0
*
       MINW = NWEIGHTNOD(1)
*
       WRITE(*,*)
       WRITE(*,'(A2,A)') '  ', REPEAT ('=', 81)
       WRITE(*,'(A2,A81)') '  NODE           SIZE      ',
     &     'MAPPING   BATCH-PER-NODE  BLOCK-PER-NODE  TOTAL-WEIGHT'
       WRITE(*,'(A2,A)') '  ', REPEAT ('=', 81)
       DO INODE = 1, NNOD
         WRITE(*,'(I8,I14,I13,I9,I16,I20)') INODE, ISIZE(INODE),
     &                                 IMAP2(INODE), IBPNOD(INODE,2),
     &                                 IBPNOD(INODE,1),NWEIGHTNOD(INODE)
         !
         ! CI Block statistics
         !
         ISUM    = ISUM    + ISIZE(INODE)
         IBLKSUM(1) = IBLKSUM(1) + IBPNOD(INODE,1)
         IBLKSUM(2) = IBLKSUM(2) + ibpnod(inode,2)
         IF (ISIZE(INODE).GE.NMAXBLK) NMAXBLK = ISIZE(INODE)
         IF (INODE.EQ.1) THEN
          NMINBLK = ISIZE(INODE)
         ELSE
          IF (ISIZE(INODE).LT.NMINBLK) NMINBLK = ISIZE(INODE)
         END IF
         !
         ! CI Block Weight statistics
         !
         NSUW = NSUW + NWEIGHTNOD(INODE)
         IF (NWEIGHTNOD(INODE).GE.MAXW) MAXW = NWEIGHTNOD(INODE)
         IF (NWEIGHTNOD(INODE).LT.MINW) MINW = NWEIGHTNOD(INODE)
       END DO
       !
       ! Print Statistics
       !
       NDIFFBLK = NMAXBLK - NMINBLK
       NDIFFW   = MAXW    - MINW
       !
       WRITE(*,'(A2,A)') '  ', REPEAT ('=', 81)
       WRITE(*,'(A5,I17,I22,I16,I20)') '  SUM', ISUM, IBLKSUM(2),
     &                                                IBLKSUM(1), NSUW
       WRITE(*,'(A5,I17,I58)')        '  MIN',  NMINBLK,  MINW
       WRITE(*,'(A5,I17,I58)')        '  MAX',  NMAXBLK,  MAXW
       WRITE(*,'(A6,I16,I58)')        '  DIFF', NDIFFBLK, NDIFFW
       WRITE(*,'(A2,A)') '  ', REPEAT ('=', 81)
       WRITE(*,*)
      END IF
*
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE ttssweight0(ttssw,TOTWEIGHT,IBATS,NBATS,LBATS,NBLK)
*======================================================================*
*
*. 0th-order
*
*. TTSS weight based only on the number of SDs on each BATCH.
*. This new approach differs from the older on one major detail:
*. We don't consider the size of the TTSS blocks, but the size of the
*. TTSS batches.
*. The TTSS batches can contain more than one TTSS block!
*
*. Konstantinos D. Vogiatzis, May 2016, Minneapolis
*
*======================================================================*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "mafdecls.fh"
#include "global.fh"
*
*. input
      integer NBATS,NBLK
      integer LBATS(NBATS),IBATS(8,NBLK)
*
*. Output
      integer ttssw(NBATS)
*
*. Scratch
      integer I, II, IBLK
      integer ISDS(NBATS)
      integer TOTWEIGHT, MINBATCH
*
      integer ITP, IGAS, ISPGRP, jspgrp, ispgrpa, ispgrpb,
     &        idif, ncon,
     &        ioff, iblkspgrp, icon
c     real ASYM(NGAS)
*
*. Initializations
*
      ntest = 000
*
*. Obtain size of each TTSS BATCH (number of SDs)
*
      II = 0
      ISDS = 0
      DO I = 1, NBATS
       DO IBLK = 1, LBATS(I)
        II = II + 1
        ISDS(I) = ISDS(I) + IBATS(8,II)
       END DO
       !
       ! Diffenent 0th-order approximations driven by...
       ttssw(I) = ISDS(I)   ! Driven by size (number of SDs)
!      ttssw(I) = LBATS(I)  ! Driven by number of BLOCKS
!      ttssw(I) = 1         ! Driven by number of BATCHES
      END DO
*
*. Print Input
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*) '============================================ '
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        WRITE(6,*) '============================================ '
        DO ITP = 1, 2
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          write(*,*) ' IBSPGP', IBSPGP
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
        WRITE(6,*)
        !
        IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Info for each TTSS batch'
         DO I = 1, NBATS
          IA = ibats(1,I)
          IB = ibats(2,I)
          write(*,'(2i8,i15)'), I, LBATS(I), ttssw(I)
         END DO
         WRITE(6,*)
        END IF
      END IF
*
*. Total weight
*
      TOTWEIGHT = 0
      MINBATCH = ttssw(1)
      DO I = 1, NBATS
       TOTWEIGHT = TOTWEIGHT + ttssw(I)
       MINBATCH = min(ttssw(I),MINBATCH)
      END DO
*
*. Percentage
*
      DO I = 1, NBATS
       TTSSW(I) = (TOTWEIGHT/MINBATCH) * 100 * TTSSW(I) / TOTWEIGHT
      END DO
      TOTWEIGHT = 0.0d0
      DO I = 1, NBATS
       TOTWEIGHT = TOTWEIGHT + TTSSW(I)
      END DO
*
*. Print Output
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*)
       WRITE(6,*) '======================================= '
       write(*,*) ' TTSS Weights: Output from ttssweight0  '
       WRITE(6,*) '======================================= '
       write(*,*)
       do I = 1, NBATS
        write(*,'(i10,i10,a4)') I, ttssw(I), '   %'
       end do
       write(*,*) 'TOTAL:   ', TOTWEIGHT
       write(*,*)
      end if
*
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE ttssweight(nttssw,nttsswsum,iconspa,iconspb,
     &                      noctpa,noctpb,
     &                      ibats,nbats,lbats,nblk,NSMST,nssoa,nssob)
*======================================================================*
*
*. 1st-order
*
*. Calculation of the TTSS weight based on the
*
*. TTSS Weight = (# of connections) x (size of TTSS)
*  TTSS Weight = (sum(iconspa(iblock)) + sum(iconspb(iblock))) x
*  (ibats(8,iblock))
*
*
*. Konstantinos D. Vogiatzis, Mar. 2016, Minneapolis
*.                            Sep. 2016, Knoxville
*
*======================================================================*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cstate.inc'
#include "mafdecls.fh"
#include "global.fh"
*
*. input
      integer noctpa,noctpb,nbats,nblk, NSMST
      integer iconspa(noctpa,noctpa),iconspb(noctpb,noctpb)
      integer lbats(nbats),ibats(8,nblk)
!      DIMENSION NSSOA(NOCTPA,NSMST),NSSOB(NOCTPB,NSMST)
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
*
*. Output
      integer nttsswsum
      integer nttssw(nbats)
*
*. Scratch
      integer iblk, ispgrp, ncon, ioff, iblkspgrp,
     &        jatp, jbtp, jasm, jbsm, joff, jblk, IPERM
      DIMENSION LASM(4), LBSM(4), LATP(4), LBTP(4), LSGN(5), LTRP(5)
*
*. Initializations
*
      ntest = 100
      nttsswsum = 0
      !
      FACTOR00 = 1.000
      FACTOR11 = 1.000
      FACTOR20 = 1.000
*
*. Loop over all TTSS blocks
*
      iblk = 0
      DO I = 1, NBATS
       ncon = 0
       NTTSSW(I) = 0
       isize = 0
       DO II = 1, LBATS(I)
!        ncon = 0
        iblk = iblk + 1
        !
        JATP = ibats(1,iblk)
        JBTP = ibats(2,iblk)
        JASM = ibats(3,iblk)
        JBSM = ibats(4,iblk)
        JOFF = ibats(5,iblk)
        isize = isize + ibats(8,iblk)
        !
!        PS  = PSSIGN
!        PL  = 1
        ! PL  not needed for IDC = 1 or 2
!        ISTRFL = 1
!        CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
!     &             LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
!        DO IPERM = 1, NPERM
!         LLASM = LASM(IPERM)
!         LLBSM = LBSM(IPERM)
!         LLATP = LATP(IPERM)
!         LLBTP = LBTP(IPERM)
         !
         ! Loop over Sigma blocks
         !
         jblk = 0
         DO J = 1, NBATS
          DO JJ = 1, LBATS(J)
           jblk = jblk + 1
!           IF(ibats(1,jblk).GT.0) THEN
           INTERACT = 0
           IATP = ibats(1,jblk)
           IBTP = ibats(2,jblk)
           IASM = ibats(3,jblk)
           IBSM = ibats(4,jblk)
!           NIA = NSSOA(IASM,IATP)
!           NIB = NSSOB(IBSM,IBTP)
!           IF(NIA*NIB.NE.0) THEN
            !
            ! Are the two blocks connected by allowed excitation?
            MXEXC = 2
            IF((ICONSPA(IATP,JATP).LE.1.AND.
     &          ICONSPB(IBTP,JBTP).LE.1     )   .OR.
     &         (ICONSPA(IATP,JATP).EQ.MXEXC.AND.
     &          IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) .OR.
     &         (ICONSPB(IBTP,JBTP).EQ.MXEXC.AND.
     &          IATP.EQ.JATP.AND.IASM.EQ.JASM)     )
     &         INTERACT = 1
!            IF((ICONSPA(IATP,LLATP).LE.1.AND.
!     &          ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
!     &         (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
!     &          IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
!     &         (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
!     &          IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )
!     &         INTERACT = 1
            IF(INTERACT.NE.0) THEN
             ncon = ncon + 1
*
*. My way to calculate the TTSS weights
*
!             IF (ICONSPA(IATP,JATP).EQ.0.AND.
!     &           ICONSPB(IBTP,JBTP).EQ.0) 
!     &           NTTSSW(I) = NTTSSW(I) + INT(FACTOR00 * ibats(8,jblk))
!             IF (ICONSPA(IATP,JATP).EQ.1.AND.
!     &           ICONSPB(IBTP,JBTP).EQ.1)
!     &           NTTSSW(I) = NTTSSW(I) + INT(FACTOR11 * ibats(8,jblk))
!             IF (ICONSPA(IATP,JATP).EQ.2.AND.
!     &           ICONSPB(IBTP,JBTP).EQ.0)
!     &           NTTSSW(I) = NTTSSW(I) + INT(FACTOR20 * ibats(8,jblk))
!             IF (ICONSPA(IATP,JATP).EQ.0.AND.
!     &           ICONSPB(IBTP,JBTP).EQ.2)
!     &           NTTSSW(I) = NTTSSW(I) + INT(FACTOR20 * ibats(8,jblk))
*
*. Bert's suggestion for the TTSS weights: CI-size x sigma-size
*
!             NTTSSW(I) = NTTSSW(I) + ibats(8,iblk) * ibats(8,jblk)
c                   JJJ = JJJ + 1
c                   ICCTOS(JSBLOCK,JBLOCK) = 1
c                   ICWEIGHT(JJJ) = JSBLOCK
             if (ga_nodeid().eq.0.and.ntest.ge.1000)
     &        write(*,'(i4,i10,10i4)') jblk,
     &        ibats(8,jblk),
     &        ICONSPA(IATP,JATP),
     &        ICONSPB(IBTP,JBTP),IATP,JATP,IASM,JASM,IBTP,JBTP,IBSM,JBSM
            END IF
!           END IF
!          END IF
          !
          ! End of jbatch
          !
          END DO
         END DO
        ! 
        ! End of nperm
        !
!        END DO
C       !
C       ! Loop over alpha supergroups
C       !
C       iblkspgrp = ibats(1,iblk)
C       do ispgrp = 1, noctpa
C        ioff = noctpa * (iblkspgrp - 1)
C        ncon = ncon + conspa(ioff + ispgrp)
C        ntotcon = ntotcon + conspa(ioff + ispgrp)
C       end do
C       !
C       ! Loop over  beta supergroups
C       !
C       iblkspgrp = ibats(2,iblk)
C       do ispgrp = 1, noctpb
C        ioff = noctpb * (iblkspgrp - 1)
C        ncon = ncon + conspb(ispgrp)
C        ntotcon = ntotcon + conspa(ioff + ispgrp)
C       end do
*
*. Knecht's weights based on connectivity and the size of the CI vector
*
!       nttssw(I) = nttssw(I) + ncon * ibats(8,iblk)
!       if (ga_nodeid().eq.0.and.ntest.ge.100)
!     &  write(*,'(a,i4,4i12)') 'Knecht', I, iblk, ncon,
!     &                          ibats(8,iblk), nttssw(I)
       !
       ! End do loop for ibatch
       END DO ! End do loop for ibatch
*
*. My way to calculate the TTSS weights (part 2)
*
!       !
!       ! Multiply the TTSS weight with the number of TTSS batches
!       NTTSSW(I) = NTTSSW(I) * LBATS(I)
*
*. Number of Connections as the TTSS weight
       nttssw(I) = ncon * isize
       if (ga_nodeid().eq.0.and.ntest.ge.100)
     &     write(*,*) 'ncon', I, ncon, NTTSSW(I)
       nttsswsum = nttsswsum + nttssw(I) 
      END DO
!*
!*. NTTSSW = NTTSSW / NBLOCKS (test)
!*
!      DO I = 1, NBATS
!       nttssw(I) = INT(nttssw(I) /LBATS(I))
!      END DO
*
*. Print Output
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*)
       write(*,*)
       WRITE(6,*) '=========================== '
       write(*,*) ' Output from ttssweight1.0  '
       WRITE(6,*) '=========================== '
       WRITE(*,*)
       WRITE(*,*) '  Info for each TTSS batch  '
       WRITE(*,*)
       write(*,*) ' BATCH   BLOCKS              LENGTH       WEIGHT'
*
       II = 0
       DO I = 1, NBATS
        NSDB = 0 ! Number of SDs in a given TTSS Batch
        DO IBLK = 1, LBATS(I)
         II = II + 1
         NSDB = NSDB + ibats(8,II)
        END DO
        WRITE(*,'(i6,i10,i20,i20)') I, LBATS(I), NSDB, nttssw(i)
       END DO
       write(*,'(a,i55)') 'TOTAL              ', nttsswsum
       write(*,*)
      end if
*
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE ttssweight15(ttssw,conspa,conspb,noctpa,noctpb,
     &                       ibats,nblk)
*======================================================================*
*
*. 1.5-order
*
*  TTSS weight based on the alpha/beta connections
* 
*. Konstantinos D. Vogiatzis, Apr. 2016, Minneapolis
*
*======================================================================*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "mafdecls.fh"
#include "global.fh"
*
*. input
      integer noctpa,noctpb,ngas,nblk
      integer conspa(noctpa,noctpa),conspb(noctpb,noctpb)
      integer ibats(8,nblk)
*
*. Output
      real ttssw(nblk)
*
*. Scratch
      real totweight
      integer iblk, jblk, ispgrp, jspgrp, ispgrpa, ispgrpb, idif, ncon,
     &        ioff, iblkspgrp, icon, ntotweight
      integer ISPGRPCONAB(NOCTPA,NOCTPB)
*
*. Initializations
*
      ntest = 1000
      ttssw  = 0
      ntotweight = 0
*
*. Create alpha/beta connection matrix for string types
*
      !
      ! Loop over alpha supergroups
      ! 
      DO ISPGRP = 1, NOCTPA
       ISPGRPA = IOCTPA-1+ISPGRP
       !
       ! Loop over beta supergroups
       !
       DO JSPGRP = 1, NOCTPB
        ISPGRPB = IOCTPB-1+JSPGRP
        IDIF = 0
        !
        ! Calculate number of connections
        !
        DO IGAS = 1, NGAS
         IDIF = IDIF
     &   + ABS(NELFSPGP(IGAS,ISPGRPA)-NELFSPGP(IGAS,ISPGRPB))
        END DO
        NEXC = IDIF/2
        ISPGRPCONAB(ISPGRP,JSPGRP) = NEXC
       END DO
      END DO
*
*. Print Alpha/beta connection matrix
*
      IF(NTEST.GE.100.and.ga_nodeid().eq.0) THEN
        WRITE(6,*)
        WRITE(6,*) '======================================== '
        WRITE(6,*) ' Alpha/beta Supergroup connection matrix '
        WRITE(6,*) '======================================== '
        WRITE(6,*)
        WRITE(6,*) ' Supergroup connection matrix '
        CALL IWRTMA(ISPGRPCONAB,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
*. Obtain metric for each TTSS blocks
*
      !
      ! alpha/alpha connections
      !
      DO IBLK = 1, NBLK
       ISPGRPA = IBATS(1,IBLK)
       DO JBLK = 1, NBLK
        JSPGRPA = IBATS(1,JBLK)
*
        ICON = CONSPA(ISPGRPA,JSPGRPA)
        IF (ICON.EQ.0) ISUM = (NACOB*NACTEL)**2
        IF (ICON.EQ.1) ISUM =  NACOB*NACTEL
        IF (ICON.EQ.2) ISUM = 1
        IF (ICON.GT.2) ISUM = 0
        TTSSW(IBLK) = TTSSW(IBLK) + ISUM*IBATS(8,IBLK)
       END DO
      END DO
      !
      ! beta/beta connections
      !
      DO IBLK = 1, NBLK
       ISPGRPB = IBATS(2,IBLK)
       DO JBLK = 1, NBLK
        JSPGRPB = IBATS(2,JBLK)
*
        ICON = CONSPB(ISPGRPB,JSPGRPB)
        IF (ICON.EQ.0) ISUM = (NACOB*NACTEL)**2
        IF (ICON.EQ.1) ISUM =  NACOB*NACTEL
        IF (ICON.EQ.2) ISUM = 1
        IF (ICON.GT.2) ISUM = 0
        TTSSW(IBLK) = TTSSW(IBLK) + ISUM*IBATS(8,IBLK)
       END DO
      END DO
      !
      ! alpha/beta connections
      !
      DO IBLK = 1, NBLK
       ISPGRPA = IBATS(1,IBLK)
       DO JBLK = 1, NBLK
        ISPGRPB = IBATS(2,JBLK)
*
        ICON = ISPGRPCONAB(ISPGRPA,ISPGRPB)
!        write(*,*) 'ISPGRPA,ISPGRPB,ICON', ISPGRPA,ISPGRPB,ICON
        IF (ICON.EQ.0) ISUM = (NACOB*NACTEL)**2
        IF (ICON.EQ.1) ISUM =  NACOB*NACTEL
        IF (ICON.EQ.2) ISUM = 1
        IF (ICON.GT.2) ISUM = 0
        TTSSW(IBLK) = TTSSW(IBLK) + ISUM*IBATS(8,IBLK)
       END DO
      END DO
*
*. Total weight
*
      DO IBLK = 1, NBLK
       NTOTWEIGHT = NTOTWEIGHT + ttssw(IBLK)
      END DO
*
*. Percentage
*
      DO IBLK = 1, NBLK
       TTSSW(IBLK) = 100 * TTSSW(IBLK) / NTOTWEIGHT 
      END DO
      TOTWEIGHT = 0
      DO IBLK = 1, NBLK
       TOTWEIGHT = TOTWEIGHT + TTSSW(IBLK)
      END DO
*
*. Print Output
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*)
       WRITE(6,*) '========================== '
       write(*,*) ' Output from ttssweight1.5 ' 
       WRITE(6,*) '=========================== '
       write(*,*)
       do iblk = 1, nblk
        write(*,'(i10,f10.3,a4)') iblk, ttssw(iblk), '   %'
       end do
       write(*,*) 'TOTAL         ', TOTWEIGHT
       write(*,*)
      end if
**
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE ttssweight2(ttssw,noctpa,noctpb,
     &                       ibats,NBATS,LBATS,nblk,NSSOA,NSSOB)
*======================================================================*
*
*. 2nd-order
*
*. TTSS weight based on 1. the number of GAS spaces of each supergroup
*.                      2. the beta/alpha ratio of each supergroup
* 
* NSSOA : Number of strings per type and symmetry for alpha strings
* NSSOB : Number of strings per type and symmetry for beta strings
*
*. Konstantinos D. Vogiatzis, Apr. 2016, Minneapolis
*. after Jeppe visiting Minneapolis
*
*. Modified on May 2016 for using TTSS Batches
*
*======================================================================*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "mafdecls.fh"
#include "global.fh"
*
*. input
      integer noctpa,noctpb,nblk,NBATS
      integer LBATS(NBATS),ibats(8,nblk)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
*
*. Output
      real ttssw(NBATS)
*
*. Scratch
      real TOTWEIGHT
      integer ITP, IGAS, ISPGRP, jspgrp, ispgrpa, ispgrpb,
     &        idif, ncon,
     &        ioff, iblkspgrp, icon
      integer NGASSPGRPA(NOCTPA),NGASSPGRPB(NOCTPB)
      real ASYM(NGAS), AGASAB(nblk)
*
*. Initializations
*
      ntest = 10000
*
      ttssw = 1.0d0
      NGASSPGRPA = 0
      NGASSPGRPB = 0
*
      IF (ga_nodeid().eq.0) write(*,*)
     & 'TTSS weight = P E I R A M A T I Z O M A I '
*
*. Obtain number of occupied GAS spaces per supergroup
*
      !
      ! Loop over alpha supergroups
      !
      DO ISPGRP = 1, NOCTPA
       DO IGAS = 1, NGAS
        IF (NELFSPGP(IGAS,ISPGRP).NE.0)
     &      NGASSPGRPA(ISPGRP) = NGASSPGRPA(ISPGRP) + 1
       END DO
      END DO
      !
      ! Loop over  beta supergroups
      !
      DO ISPGRP = 1, NOCTPB
       DO IGAS = 1, NGAS
        IF (NELFSPGP(IGAS,ISPGRP).NE.0)
     &      NGASSPGRPB(ISPGRP) = NGASSPGRPB(ISPGRP) + 1
       END DO
      END DO
*
*. Calculate TTSS weight
*
      II = 0
      DO I = 1, NBATS
       ISIZE = 0
       DO IBLK = 1, LBATS(I)
        II = II + 1
        !
        ! Obtain number of alpha and beta strings per TTSS Batch
        !
        IATP = ibats(1,II)
        IBTP = ibats(2,II)
        IASM = ibats(3,II)
        IBSM = ibats(4,II)
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
        FACTOR = (REAL(NIB) / REAL(NIA))**2
        if (NTEST.GE.100.and.ga_nodeid().eq.0) write(*,'(a,2i10,f10.6)')
     &                        'FACTOR1', NIA, NIB, FACTOR
        !
        ! Obtain number of ooccupied alpha and beta GAS spaces
        !
        FACTORR = REAL(NGASSPGRPB(IBTP)) / REAL(NGASSPGRPA(IATP))
        if (NTEST.GE.100.and.ga_nodeid().eq.0) write(*,'(a,2i4,f10.6)')
     &                        'FACTOR2', IATP, IBTP, FACTORR
        !
        ISIZE = ibats(8,II) ! FTIAXTO!
        !
        ! Obtain weight for specific TTSS block
        weight = ISIZE * FACTOR * FACTORR
        if (NTEST.GE.100.and.ga_nodeid().eq.0) write(*,*)
     &    'Weight for Block', II, weight 
       END DO
       ttssw(I) = ttssw(I) + weight
      END DO
*      ! Update TTSS block weight
*      !
*      II = 0
*      DO I = 1, NBATS
*       DO IBLK = 1, LBATS(I)
*        II = II + 1
*        IA = ibats(1,II)
*        IB = ibats(2,II)
*        AX = 1.0d0
*        AGASAB(II) = (REAL(NGAS)/REAL(NGASSPGRPA(IA)) * 
*     &                REAL(NGAS)/REAL(NGASSPGRPB(IB)))**AX
*        ttssw(I) = ttssw(I) * AGASAB(II)
*       END DO
*      END DO
**
**. Calculate beta/alpha ratio for each TTSS block and GAS
**
*      II = 0
*      DO I = 1, NBATS
*       DO IBLK = 1, LBATS(I)
*        II = II + 1
*        IA = ibats(1,II)
*        IB = ibats(2,II)
*        ASYM = 0.0d0
*        DO IGAS = 1, NGAS
*         IF ((NELFSPGP(IGAS,IA).EQ.0).OR.(NELFSPGP(IGAS,IB).EQ.0)) THEN
*          ASYM(IGAS) = 1.0d0
*         ELSE
*          ASYM(IGAS) = REAL(NELFSPGP(IGAS,IB)) / REAL(NELFSPGP(IGAS,IA))
*         END IF
*         !
*         ! Update TTSS block weight
*         !
*         AX = 0.5d0
*         ttssw(I) = ttssw(I) * ASYM(IGAS)**AX
*        END DO
**
*        if (NTEST.GE.10000.and.ga_nodeid().eq.0)
*     &      write(*,'(i5,8f8.3)') II, (ASYM(IGAS), IGAS = 1, NGAS) 
*       END DO
*      END DO
**
**. Print Input
**
*      IF(NTEST.GE.10.and.ga_nodeid().eq.0) THEN
*        WRITE(6,*) '============================================ '
*        WRITE(6,*) ' Distribution of electrons in Active spaces '
*        WRITE(6,*) '============================================ '
*        DO ITP = 1, 2
*          WRITE(6,*) ' String type ', ITP
*          WRITE(6,*) ' Row : active space, Column: supergroup '
*          NSPGP = NSPGPFTP(ITP)
*          IBSPGP = IBSPGPFTP(ITP)
*          write(*,*) ' IBSPGP', IBSPGP
*          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
*        END DO
*        WRITE(6,*)
*        !
*        IF(NTEST.GE.1000) THEN
*         WRITE(6,*) ' Info for each TTSS block'
*         DO I = 1, NBLK
*          IA = ibats(1,I)
*          IB = ibats(2,I)
*          write(*,'(3i8,a,8i2,a,8i2,f10.6)')
*     &    I, NGASSPGRPA(IA), NGASSPGRPB(IB),
*     &    ' | ', (NELFSPGP(IGAS,IA), igas = 1, 8),
*     &    ' | ', (NELFSPGP(IGAS,IB), igas = 1, 8), AGASAB(I)
*         END DO
*         WRITE(6,*)
*        END IF
*      END IF
*c*
c*. Obtain metric for each TTSS blocks
c*
c      !
c      ! alpha/alpha connections
c      !
c      DO IBLK = 1, NBLK
c       ISPGRPA = IBATS(1,IBLK)
c       DO JBLK = 1, NBLK
c        JSPGRPA = IBATS(1,JBLK)
c*
c!        ICON = CONSPA(ISPGRPA,JSPGRPA)
c        IF (ICON.EQ.0) ISUM = (NACOB*NACTEL)**2
c        IF (ICON.EQ.1) ISUM =  NACOB*NACTEL
c        IF (ICON.EQ.2) ISUM = 1
c        IF (ICON.GT.2) ISUM = 0
c        TTSSW(IBLK) = TTSSW(IBLK) + ISUM*IBATS(8,IBLK)
c       END DO
c      END DO
c      !
c      ! beta/beta connections
c      !
c      DO IBLK = 1, NBLK
c       ISPGRPB = IBATS(2,IBLK)
c       DO JBLK = 1, NBLK
c        JSPGRPB = IBATS(2,JBLK)
c*
c!        ICON = CONSPB(ISPGRPB,JSPGRPB)
c        IF (ICON.EQ.0) ISUM = (NACOB*NACTEL)**2
c        IF (ICON.EQ.1) ISUM =  NACOB*NACTEL
c        IF (ICON.EQ.2) ISUM = 1
c        IF (ICON.GT.2) ISUM = 0
c        TTSSW(IBLK) = TTSSW(IBLK) + ISUM*IBATS(8,IBLK)
c       END DO
c      END DO
c      !
c      ! alpha/beta connections
c      !
c      DO IBLK = 1, NBLK
c       ISPGRPA = IBATS(1,IBLK)
c       DO JBLK = 1, NBLK
c        ISPGRPB = IBATS(2,JBLK)
c*
c!        ICON = ISPGRPCONAB(ISPGRPA,ISPGRPB)
c!        write(*,*) 'ISPGRPA,ISPGRPB,ICON', ISPGRPA,ISPGRPB,ICON
c        IF (ICON.EQ.0) ISUM = (NACOB*NACTEL)**2
c        IF (ICON.EQ.1) ISUM =  NACOB*NACTEL
c        IF (ICON.EQ.2) ISUM = 1
c        IF (ICON.GT.2) ISUM = 0
c        TTSSW(IBLK) = TTSSW(IBLK) + ISUM*IBATS(8,IBLK)
c       END DO
c      END DO
*
*. Total weight
*
      TOTWEIGHT = 0.0d0
      DO I = 1, NBATS
       TOTWEIGHT = TOTWEIGHT + ttssw(I)
      END DO
      if (ga_nodeid().eq.0) write(*,*) 'TOTWEIGHT', TOTWEIGHT
*
*. Percentage
*
      DO I = 1, NBATS
       TTSSW(I) = 100 * TTSSW(I) / TOTWEIGHT
      END DO
      TOTWEIGHT = 0.0d0
      DO I = 1, NBATS
       TOTWEIGHT = TOTWEIGHT + TTSSW(I)
      END DO
*
*. Print Output
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*)
       WRITE(6,*) '======================================= '
       write(*,*) ' TTSS Weights: Output from ttssweight2  '
       WRITE(6,*) '======================================= '
       write(*,*)
       do I = 1, NBATS
        write(*,'(i10,f10.3,a4)') I, ttssw(I), '   %'
       end do
       write(*,*) 'TOTAL:   ', TOTWEIGHT
       write(*,*)
      end if
*
      RETURN
      END
*======================================================================*
*======================================================================*
*
      SUBROUTINE ttssweight3(ttssw,ibats,NBATS,LBATS,nblk)
*======================================================================*
*
*. 3rd-order
*
*. TTSS weight based on precalculated timings from one CASCI iteration
*. with maximized number of connections
*
*. Konstantinos D. Vogiatzis, May 2016, Minneapolis
*
*======================================================================*
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
#include "mafdecls.fh"
#include "global.fh"
*
*. input
      integer nblk,NBATS
      integer LBATS(NBATS),ibats(8,nblk)
*
*. Output
      real ttssw(NBATS)
*
*. Scratch
      real TOTWEIGHT
      integer ITP, IGAS, ISPGRP, jspgrp, ispgrpa, ispgrpb,
     &        idif, ncon,
     &        ioff, iblkspgrp, icon
      real ASYM(NGAS), AGASAB(nblk)
*
*. Initializations
*
      ntest = 10000
*
      ttssw = 1.0d0
*
      IF (ga_nodeid().eq.0) write(*,*)
     & 'TTSS weight = precalculated timings'
*
*. Read precalculated timings from file
*
      iunit = 670
      open(iunit,file='timings',status='unknown')
      rewind(iunit)
      DO I = 1, NBATS
       READ(IUNIT,*) ttssw(I)
      END DO
      close(iunit)
*
*. Total weight
*
      TOTWEIGHT = 0.0d0
      DO I = 1, NBATS
       TOTWEIGHT = TOTWEIGHT + ttssw(I)
      END DO
      if (ga_nodeid().eq.0) write(*,*) 'TOTWEIGHT', TOTWEIGHT
*
*. Percentage
*
      DO I = 1, NBATS
       TTSSW(I) = 100 * TTSSW(I) / TOTWEIGHT
      END DO
      TOTWEIGHT = 0.0d0
      DO I = 1, NBATS
       TOTWEIGHT = TOTWEIGHT + TTSSW(I)
      END DO
*
*. Print Output
*
      if (ga_nodeid().eq.0.and.ntest.ge.10) then
       write(*,*)
       WRITE(6,*) '======================================= '
       write(*,*) ' TTSS Weights: Output from ttssweight3  '
       WRITE(6,*) '======================================= '
       write(*,*)
       do I = 1, NBATS
        write(*,'(i10,f10.3,a4)') I, ttssw(I), '   %'
       end do
       write(*,*) 'TOTAL:   ', TOTWEIGHT
       write(*,*)
      end if
*
      RETURN
      END
*======================================================================*
*======================================================================*
      integer function NBATCHES_FEW_BLOCKS(NCBATCH,LBATCH)
*======================================================================*
*
*. Input
      INTEGER NCBATCH
      INTEGER LBATCH(NCBATCH)
*
*. Scratch
      INTEGER I, ITHRS
*
      ITHRS = 5
*
      NBATCHES_FEW_BLOCKS = 0
      DO I = 1, NCBATCH
       IF (LBATCH(I).LE.ITHRS)
     &  NBATCHES_FEW_BLOCKS = NBATCHES_FEW_BLOCKS + 1
      END DO
*
      return
      end
*======================================================================*

c $Id$
