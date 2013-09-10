      SUBROUTINE DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR,EXPS2,IDOSRHO12,SRHO1,
     &                  RHO2AA,RHO2AB,RHO2BB)
*
* Density matrices between L and R
*
* I12 = 1 => only one-bodydensity
* I12 = 2 => one- and two-body-density matrices
*
* Jeppe Olsen,      Oct 94
* GAS modifications Aug 95
* Two body density added, '96
*
* Table-Block driven, June 97
* Spin density added, Jan. 99
* two-electron spindensities added, Sept. 2004 ( for singularities in IC...)
*
* Restriction of densities to active spaces allowed, Sept. 2005
* Prepared for explicit inactive/secondary orbital spaces, June 2010
*
* Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
* ijkl = ij*(ij-1)/2+kl, ij.ge.kl
*
* If the twobody density matrix is calculated, then also the
* expectation value of the spin is evalueated.
* The latter is realixed as
* S**2 
*      = S+S- + Sz(Sz-1)
*      = -Sum(ij) a+i alpha a+j beta a i beta a j alpha + Nalpha +
*        1/2(N alpha - N beta))(1/2(N alpha - Nbeta) - 1)
*
* If IDOSRHO12 = 1, spin density is also calculated
*
* if IDOSRHO12 = 2, then the spin-components of the 
* two-body density are also calculated
*
* RHO2AA(i,j,k,l) = <0!a+_ialpha a+_jalpha a_kalpha a_lalpha!0>
*                   i.ge.j, k.ge.l
* RHO2AA(i,j,k,l) = <0!a+_ibeta a+_jbeta a_kbeta a_lbeta!0>
*                   i.ge.j, k.ge.l
* RHO2AB(i,j,k,l) = <0!a+_ialpha a+_jbeta a_kbeta a_lalpha!0>
*                   No restrictions on i,j,k,l
*
* Call tree for densities : 
* =========================
*  DENSI2 - GASDN2 - GSDNBB2 - GSBBD1
*                            - GSBBD2A
*                            - GSBBD2B


      INCLUDE 'wrkspc.inc' 
      REAL*8 INPRDD
*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
      INCLUDE 'orbinp.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'glbbas.inc'
*
      INCLUDE 'csmprd.inc'
c      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
c      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
c     &              ADSXA(MXPOBS,2*MXPOBS),
c     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      INCLUDE 'lucinp.inc'
      INCLUDE 'clunit.inc'
*. Scratch for string information
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. Local scratch
      DIMENSION IDACTSPC(MXPNGAS)
* IDTFREORD Should give reordering full TS order to Density TS order 
      
*. Specific input 
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*),SRHO1(*)
*
      DIMENSION RHO2AA(*),RHO2AB(*),RHO2BB(*)
*
      NTEST = 0
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
      CALL QENTER('DENSI')
*
*. Divide orbital space into inactive(hole), active(valence) and secondary(particle)
*  based on space ISSPC. Result is stored in IHPVGAS
*. IF ISSPC and ICSPC differs in division into HPV, and some orbital spaces are
*. excluded one may run into problems
      CALL CC_AC_SPACES(ISSPC,IREFTYP)
      DO IGAS = 1, NGAS
        ITYP = IHPVGAS(IGAS)
        IDACTSPC(IGAS) = 0
        IF(ITYP.EQ.1.AND.IDENS_IN.EQ.1) IDACTSPC(IGAS) = 1
        IF(ITYP.EQ.2.AND.IDENS_SEC.EQ.1) IDACTSPC(IGAS) = 1
        IF(ITYP.EQ.3.AND.IDENS_AC.EQ.1) IDACTSPC(IGAS) = 1
      END DO
      IF(NTEST.GT.0) THEN
       WRITE(6,*) ' Active orbital spaces in DENSI '
       CALL IWRTMA(IDACTSPC,1,NGAS,1,NGAS)
      END IF
*. Number of orbitals for which densities will be constructed 
      NDACTORB = 0
      DO IGAS = 1, NGAS
        IF(IDACTSPC(IGAS).EQ.1) NDACTORB = NDACTORB + NOBPT(IGAS)
      END DO
      IF(NTEST.GT.0) 
     & WRITE(6,*) ' Number of active orbitals for densities',NDACTORB
*. Offsets to restricted set of orbital spaces
      CALL IB_FOR_SEL_ORBSPC(NOBPTS,NOBPS_SEL,IOBPTS_SEL,IDACTSPC,NGAS,
     &     MXPNGAS,NSMST)
C?    WRITE(6,*) ' NOBPS_SEL after call to IB_FOR... '
C?    CALL IWRTMA(NOBPS_SEL,1,NSMOB,1,NSMOB)
C     IB_FOR_SEL_ORBSPC(NOBPTS,NOBPS_SEL,IOBPTS_SEL,I_SEL,
C    &           NGAS,MXPNGAS,NSYM)
*. And reorder arrays for going betweeen  complete set of 
*. orbitals and restricted set of orbitals
      CALL IREO_DACT_TS(NOBPTS,IOBPTS_SEL,IDACTSPC,IREO_SELTF,
     &     IREO_FTSEL,NGAS,MXPNGAS,NSMOB)
C     IREO_DACT_TS(NOBPTS,IOBPTS_SEL,I_SEL,
C    &           IDTFREO,IFTDREO,NGAS,MXPNGAS,NSMOB)
*
      ZERO = 0.0D0
      CALL SETVEC(RHO1,ZERO ,NDACTORB ** 2 )
      IF(I12.EQ.2) 
     &CALL SETVEC(RHO2,ZERO ,NDACTORB ** 2 *(NDACTORB**2+1)/2)
*
      IF(IDOSRHO12.GE.1) THEN
        IDOSRHO1 = 1
      ELSE 
        IDOSRHO1 = 0
      END IF
*
      IF(IDOSRHO12.EQ.2) THEN
        IDOSRHO2 = 1
      ELSE 
        IDOSRHO2 = 0
      END IF
*.
      IF(IDOSRHO1.EQ.1) THEN  
        CALL SETVEC(SRHO1,ZERO,NDACTORB ** 2)
      END IF
*
      IF(IDOSRHO2.EQ.1) THEN 
        LEN_RHO2AA = (NDACTORB*(NDACTORB+1)/2)**2
        LEN_RHO2AB = NDACTORB ** 4
        CALL SETVEC(RHO2AA,ZERO,LEN_RHO2AA)
        CALL SETVEC(RHO2BB,ZERO,LEN_RHO2AA)
        CALL SETVEC(RHO2AB,ZERO,LEN_RHO2AB)
      END IF
*          
C?     WRITE(6,*) ' ISSPC ICSPC in DENSI2 ',ISSPC,ICSPC
*
* Info for this internal space
*. type of alpha and beta strings - as H does not change 
*. the number of electrons, I do not distinguish between spaces for C
*and S
      IF(ICSPC.LE.NCMBSPC) THEN
       IATP = 1
       IBTP = 2
      ELSE
       IATP = IALTP_FOR_GAS(ICSPC)
       IBTP = IBETP_FOR_GAS(ICSPC)
C?     WRITE(6,*) ' DENSI2 : IATP, IBTP = ', IATP, IBTP
      END IF
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. alpha and beta strings with an electron removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-1,IATPM1)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-1,IBTPM1)
*. alpha and beta strings with two electrons removed
      CALL  FIND_TYPSTR_WITH_TOTOCC(NAEL-2,IATPM2)
      CALL  FIND_TYPSTR_WITH_TOTOCC(NBEL-2,IBTPM2)

*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
      JATP = 1
      JBTP = 2
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offsets for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      ILSM = ISSM
      IRSM = ICSM
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. connection matrices for supergroups
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(WORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR          
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MAXA = MAX(MAXA,MAXA0)   
      MAXB = MAX(MAXB,MAXB0)   
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRDEN.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
*. replace with MXINKA !!!
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
C?    WRITE(6,*) ' DENSI2 : MAXI MAXK ', MAXI,MAXK
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(IPRDEN.GE.2) write(6,*) ' DENSI2 : MXSB MXTSOB MXSOOB ',
     &       MXSB,MXTSOB,MXSOOB 
      IF(ISIMSYM.NE.1) THEN
        LSCR1 = MXSOOB
      ELSE
        LSCR1 = MXSOOB_AS
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
      IF(ICISTR.EQ.1) THEN
        CALL MEMMAN(KCB,LSCR1,'ADDL  ',2,'KCB   ')
        CALL MEMMAN(KSB,LSCR1,'ADDL  ',2,'KSB   ')
      END IF
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
*
      CALL IAIBCM(ISSPC,WORK(KSIOIO))
      CALL IAIBCM(ISSPC,WORK(KCIOIO))
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESCPH(WORK(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSTFSMSPGP,MXPNSMST,
     &     NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &     NELFSPGP,
     &     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &     IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,MXADKBLK_AS,
     &     MX_NSPII)
      IF(IPRDEN.GE.2) THEN
        WRITE(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
      LSCR12 = MAX(LSCR1,2*LSCR2)
*. It is assumed that the third block already has been allocated, so
      KC2 = KVEC3
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL MEMMAN(KI1,  LSCR3       ,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KI2,  LSCR3       ,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KI3,  LSCR3       ,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KI4,  LSCR3       ,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
      CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
      CALL MEMMAN(KXI3S,LSCR3       ,'ADDL  ',2,'XI3S  ')
      CALL MEMMAN(KXI4S,LSCR3       ,'ADDL  ',2,'XI4S  ')
*. Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.0 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
* scratch space containing active one body
      CALL MEMMAN(KRHO1S,NDACTORB ** 2,'ADDL  ',2,'RHO1S ')
*. For natural orbitals
      CALL MEMMAN(KRHO1P,NDACTORB*(NDACTORB+1)/2,'ADDL  ',2,'RHO1P ')
      CALL MEMMAN(KXNATO,NDACTORB **2,'ADDL  ',2,'XNATO ')
*. Natural orbitals in symmetry blocks
      CALL MEMMAN(KRHO1SM,NDACTORB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KXNATSM,NDACTORB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KOCCSM,NDACTORB ,'ADDL  ',2,'RHO1S ')
*
*. Space for one block of string occupations and two arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
      DO K12 = 1, 1
        CALL MEMMAN(KLOCSTR(K12),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      END DO
      DO I1234 = 1, 2
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
*. Arrays for partitioning of Left vector = sigma 
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTL ,NTTS  ,'ADDL  ',1,'LBT_L  ')
      CALL MEMMAN(KLLEBTL,NTTS  ,'ADDL  ',1,'LEBT_L ')
      CALL MEMMAN(KLI1BTL,NTTS  ,'ADDL  ',1,'I1BT_L ')
      CALL MEMMAN(KLIBTL ,8*NTTS,'ADDL  ',1,'IBT_L  ')
      CALL MEMMAN(KLSCLFCL,NTTS, 'ADDL  ',2,'SCLF_L')
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),
     &     WORK(KLI1BTL),WORK(KLIBTL),0,ISIMSYM)
*. Number of BLOCKS
        NBLOCKL = IFRMR(WORK(KLI1BTL),1,NBATCHL)
     &         + IFRMR(WORK(KLLBTL),1,NBATCHL) - 1
*. Arrays for partitioning of Right  vector = C
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTR ,NTTS  ,'ADDL  ',1,'LBT_R  ')
      CALL MEMMAN(KLLEBTR,NTTS  ,'ADDL  ',1,'LEBT_R ')
      CALL MEMMAN(KLI1BTR,NTTS  ,'ADDL  ',1,'I1BT_R ')
      CALL MEMMAN(KLIBTR ,8*NTTS,'ADDL  ',1,'IBT_R  ')
      CALL MEMMAN(KLSCLFCR,NTTS, 'ADDL  ',2,'SCLF_R')
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),
     &     NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),
     &     WORK(KLI1BTR),WORK(KLIBTR),0,ISIMSYM)
*. Number of BLOCKS
        NBLOCKR = IFRMR(WORK(KLI1BTR),1,NBATCHR)
     &         + IFRMR(WORK(KLLBTR),1,NBATCHR) - 1
C?      WRITE(6,*) ' DENSI2T :NBLOCKR =',NBLOCKR

      IF(ICISTR.EQ.1) THEN
         WRITE(6,*) ' Sorry, ICISTR = 1 is out of fashion'
         WRITE(6,*) ' Switch to ICISTR = 2 - or reprogram '
         STOP' DENSI2T : ICISTR = 1 in use '
      ELSE IF(ICISTR.GE.2) THEN
        S2_TERM1 = 0.0D0
        CALL GASDN2(I12,RHO1,RHO2,L,R,L,R,WORK(KC2),
     &       WORK(KCIOIO),WORK(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),WORK(KCBLTP),WORK(KSBLTP),NACOB,
     &       WORK(KNSTSO(IATP)),WORK(KISTSO(IATP)),
     &       WORK(KNSTSO(IBTP)),WORK(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,      
     &       MAXK,MAXI,LSCR1,LSCR1,WORK(KCSCR),WORK(KSSCR),
     &       SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &       WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &       WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),WORK(KINSCR), 
     &       MXPOBS,IPRDEN,WORK(KRHO1S),LUL,LUR,
     &       PSSIGN,PSSIGN,WORK(KRHO1P),WORK(KXNATO),
     &       NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),WORK(KLI1BTL),
     &       WORK(KLIBTL),
     &       NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),WORK(KLI1BTR),
     &       WORK(KLIBTR),WORK(KCONSPA),WORK(KCONSPB),
     &       WORK(KLSCLFCL),WORK(KLSCLFCR),S2_TERM1,IUSE_PH,IPHGAS,
     &       IDOSRHO1,SRHO1,IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &       NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,IOBPTS_SEL,NINOB)
      END IF
C?    WRITE(6,*) ' Memcheck in densi2 after GASDN2'
C?    CALL MEMCHK
*
*
*. Add terms from hole-hole commutator
      IF(IUSE_PH.EQ.1) THEN 
*. Overlap between left and right vector
       XLR = INPRDD(L,R,LUR,LUL,1,-1)
       CALL RHO1_HH(RHO1,XLR)
      END IF
* Natural Orbitals 
C     NATORB2(RHO1,NSMOB,NOPSM_SEL,NOBPSM,NDACTORB,
C    &                  ISTREO,IFTDREO,XNAT,RHO1SM,OCCNUM,
C    &                  SCR,IPRDEN)
C     CALL NATORB2(RHO1,NSMOB,NOBPS_SEL,NTOOBS,NDACTORB,
C    &            IREOST,IREO_FTSEL,WORK(KXNATO),
C    &            WORK(KRHO1SM),WORK(KOCCSM),
C    &            WORK(KRHO1P),IPRDEN,NINOBS)
       CALL NATORB3(RHO1,NSMOB,NACOBS,NINOBS,NSCOBS,NINOB,NACOB,
     &              IREOST,WORK(KXNATO),WORK(KRHO1SM),WORK(KOCCSM),
     &              WORK(KRHO1P),IPRDEN)
C           NATORB3(RHO1,NSMOB,NACOBS,NINOBS,NSCOBS,
C    &                  NINOB,NACOB,ISTREO,XNAT,RHO1SM,OCCNUM,
C    &                  SCR,IPRDEN)


*
      IF(IPRDEN.GE.5) THEN
        WRITE(6,*) ' One-electron density matrix '
        WRITE(6,*) ' ============================'
        CALL WRTMAT(RHO1,NDACTORB,NDACTORB,NDACTORB,NDACTORB) 
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-electron density '
          CALL PRSYM(RHO2,NDACTORB**2)
        END IF
      END IF
*
      IF(I12.EQ.2) THEN
* <L!S**2|R>
        EXPS2 = S2_TERM1 + NAEL +
     &          0.5*(NAEL-NBEL)*(0.5*(NAEL-NBEL)-1)
        IF(IPRDEN.GT.0) THEN
          WRITE(6,*) ' Term 1 to S2 ', S2_TERM1
          WRITE(6,*) ' Expectation value of S2 ', EXPS2
        END IF
      ELSE
        EXPS2 = 0.0D0
      END IF
*
      IF(IDOSRHO1.EQ.1.AND.IPRDEN.GE.2) THEN
        WRITE(6,*) ' One-electron spindensity <0!E(aa) - E(bb)!0> '
        CALL WRTMAT(SRHO1,NDACTORB,NDACTORB,NDACTORB,NDACTORB)
      END IF
*
      IF(IPRDEN.GE.2.AND.IDOSRHO2.EQ.1) THEN
        WRITE(6,*) ' The RHO2AA(ij,kl) spin density '
        NDIM_AA = NDACTORB*(NDACTORB+1)/2
        CALL WRTMAT(RHO2AA,NDIM_AA,NDIM_AA,NDIM_AA,NDIM_AA)
        WRITE(6,*) ' The RHO2BB(ij,kl) spin density '
        CALL WRTMAT(RHO2BB,NDIM_AA,NDIM_AA,NDIM_AA,NDIM_AA)
        WRITE(6,*) ' The RHO2AB(ik,lj) spin density '
        NDIM_AB = NDACTORB*NDACTORB
        CALL WRTMAT(RHO2AB,NDIM_AB,NDIM_AB,NDIM_AB,NDIM_AB)
      END IF
*
      I_CHECK_SRHO2 = 1
      IF(IDOSRHO2.EQ.1.AND.I_CHECK_SRHO2.EQ.1) THEN
*. Obtain standard rho2 from rho2s and check
C             TEST_RHO2S(RHO2,RHO2AA,RHO2AB,RHO2BB,NORB)
         CALL TEST_RHO2S(RHO2,RHO2AA,RHO2AB,RHO2BB,NTOOB)
      END IF

*
*. Eliminate local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'DENSI ')
      CALL QEXIT('DENSI')
C     WRITE(6,*) ' Leaving DENSI '
      RETURN
      END
      SUBROUTINE GASDN2(I12,RHO1,RHO2,
     &           L,R,CB,SB,C2,ICOCOC,ISOCOC,ICSMOS,ISSMOS,
     &           ICBLTP,ISBLTP,NACOB,NSSOA,ISSOA,NSSOB,ISSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,
     &           MXPNGAS,NOBPTS,IOBPTS,MAXK,MAXI,LC,LS,
     &           CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &           MXPOBS,IPRNT,RHO1S,LUL,LUR,PSL,PSR,RHO1P,XNATO ,
     &           NBATCHL,LBATL,LEBATL,I1BATL,IBLOCKL,
     &           NBATCHR,LBATR,LEBATR,I1BATR,IBLOCKR,
     &           ICONSPA,ICONSPB,SCLFAC_L,SCLFAC_R,S2_TERM1,
     &           IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &           IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &           NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,IOBPTS_SEL,
     &           NINOB)
*
*
* Jeppe Olsen , Winter of 1991
* GAS modificatios, August 1995
*
* Table driven, June 97
*
* Last revision : Jan. 98 (IUSE_PH,IPHGAS added)
*                 Jan. 99 (IDOSRHO1,SRHO1 added)
*                 Sept.04 (IDOSRHO2,RHO2AA,RHO2AB,RHO2BB added)
*                 June 10 (NINOB added)
*
* =====
* Input
* =====
*
* I12    : = 1 => calculate one-electrondensity matrix
*          = 2 => calculate one-and two- electrondensity matrix
* RHO1   : Initial one-electron density matrix
* RHO2   : Initial two-electron density matrix
*
* ICOCOC : Allowed type combinations for C
* ISOCOC : Allowed type combinations for S(igma)
* ICSMOS : Symmetry array for C
* ISSMOS : Symmetry array for S
* ICBLTP : Block types for C
* ISBLTP : Block types for S
*
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* ISSOA : Offset for strings if given type and symmetry, alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* ISSOB : Offset for strings if given type and symmetry, beta strings
* NBEL  : Number of active beta electrons
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
*
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* RHO1S: Scratch array for one body
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
* The L and R vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB),ISOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST),ISSMOS(NSMST)
      INTEGER ICBLTP(*),ISBLTP(*)
      INTEGER NSSOA(NSMST,NOCTPA),ISSOA(NSMST,NOCTPA)
      INTEGER NSSOB(NSMST,NOCTPB),ISSOB(NSMST,NOCTPB)
      INTEGER SXSTSM(NSMSX,NSMST)
      INTEGER STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER ADSXA(MXPOBS,2*MXPOBS),ASXAD(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,NSMOB),IOBPTS(MXPNGAS,NSMOB)
      INTEGER NELFSPGP(MXPNGAS,*)
*. Info on batches and blocks
      INTEGER  LBATL(NBATCHL),LEBATL(NBATCHL),I1BATL(NBATCHL),
     &         IBLOCKL(8,*)
      INTEGER  LBATR(NBATCHR),LEBATR(NBATCHR),I1BATR(NBATCHR),
     &         IBLOCKR(8,*)
*. Interaction between supergroups
      INTEGER ICONSPA(NOCTPA,NOCTPA),ICONSPB(NOCTPB,NOCTPB)
*. Info on the orbital spaces that are active in density calculations
      INTEGER IDACTSPC(*),IDTFREORD(*),IFTDREORD(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),XI1S(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION RHO1S(*)
      DIMENSION SCLFAC_L(*),SCLFAC_R(*)
*.
      INTEGER LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      INTEGER RASM(4),RBSM(4),RATP(4),RBTP(4),RSGN(5),RTRP(5)
      REAL * 8 INPROD,L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*)
      DIMENSION RHO2AA(*),RHO2AB(*),RHO2BB(*),SRHO1(*)
      DIMENSION RHO1P(*),XNATO(*)
*
      CALL QENTER('GASDN')
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GASDN2 speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' NACOB,MAXK,NGAS,IDC,MXPOBS',
     &             NACOB,MAXK,NGAS,IDC,MXPOBS
        WRITE(6,*) ' LUL, LUR ', LUL,LUR
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Initial L vector '
        IF(LUL.EQ.0) THEN
          CALL WRTRS2(L,ISSMOS,ISBLTP,ISOCOC,NOCTPA,NOCTPB,
     &                NSSOA,NSSOB,NSMST)
        ELSE
          CALL WRTVCD(L,LUL,1,-1)
        END IF
        WRITE(6,*) ' Initial R vector '
        IF(LUR.EQ.0) THEN
          CALL WRTRS2(R,ICSMOS,ICBLTP,ICOCOC,NOCTPA,NOCTPB,
     &                NSSOA,NSSOB,NSMST)
        ELSE
          CALL WRTVCD(R,LUR,1,-1)
        END IF
      END IF
* Loop over batches over L blocks
      IF(LUL.NE.0) CALL REWINO(LUL)
      DO 10001 IBATCHL = 1, NBATCHL
*. Obtain L blocks
        NBLKL = LBATL(IBATCHL)
        IF(NTEST.GE.200)
     &    WRITE(6,*) ' Left batch, number of blocks',IBATCHL,NBLKL
        DO IIL  = 1,NBLKL                
          IL  = I1BATL(IBATCHL)-1+IIL
          IATP = IBLOCKL(1,IL)
          IBTP = IBLOCKL(2,IL)
          IASM = IBLOCKL(3,IL)
          IBSM = IBLOCKL(4,IL)
          IOFF = IBLOCKL(5,IL)
          IF(NTEST.GE.200)
     &    WRITE(6,*) 'IATP IBTP IASM IBSM',IATP,IBTP,IASM,IBSM
          ISCALE = 1
          IF(NTEST.GE.200)
     &    WRITE(6,*) 'IOFF ',IOFF
          CALL GSTTBL(L,SB(IOFF),IATP,IASM,IBTP,IBSM,ISOCOC,
     &                NOCTPA,NOCTPB,NSSOA,NSSOB,PSL,ISOOSC,IDC,
     &                PSL,LUL,C2,NSMST,ISCALE,SCLFAC_L(IL))
        END DO
*. Loop over batches  of R vector
        IF(LUR.NE.0) CALL REWINO(LUR)
        DO 9001 IBATCHR = 1, NBATCHR
*. Read R blocks into core
        NBLKR = LBATR(IBATCHR)
        IF(NTEST.GE.200)
     &    WRITE(6,*) ' Right batch, number of blocks',IBATCHR,NBLKR
        DO IIR  = 1,NBLKR                
          IR  = I1BATR(IBATCHR)-1+IIR       
          JATP = IBLOCKR(1,IR)
          JBTP = IBLOCKR(2,IR)
          JASM = IBLOCKR(3,IR)
          JBSM = IBLOCKR(4,IR)
          JOFF = IBLOCKR(5,IR)
          IF(NTEST.GE.200)
     &    WRITE(6,*) ' JATP JBTP JASM JBSM ',JATP,JBTP,JASM,JBSM
*. Read R blocks into core
*
*. Only blocks interacting with current batch of L are read in
*. Loop over L  blocks in batch
          DO IIL = 1, NBLKL
            IL  = I1BATL(IBATCHL)-1+IIL       
            IATP = IBLOCKL(1,IL)
            IBTP = IBLOCKL(2,IL)
            IASM = IBLOCKL(3,IL)
            IBSM = IBLOCKL(4,IL)
*. Well, permutations of L blocks
            ISTRFL = 1
            PL = 1.0D0
            CALL PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PL,
     &              LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
            DO IPERM = 1, NPERM
              IIASM = LASM(IPERM)
              IIBSM = LBSM(IPERM)
              IIATP = LATP(IPERM)
              IIBTP = LBTP(IPERM)

              IAEXC = ICONSPA(IIATP,JATP)
              IBEXC = ICONSPB(IIBTP,JBTP)
              IF(IAEXC.EQ.0.AND.IIASM.NE.JASM) IAEXC = 1
              IF(IBEXC.EQ.0.AND.IIBSM.NE.JBSM) IBEXC = 1
              IABEXC = IAEXC + IBEXC
              IF(IABEXC.LE.I12) THEN
                INTERACT = 1
              END IF
            END DO
          END DO
*.          ^ End of checking whether C-block is needed
          ISCALE = 1
          IF(INTERACT.EQ.1) THEN
            ISCALE = 1
            CALL GSTTBL(R,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                  NOCTPA,NOCTPB,NSSOA,NSSOB,PSR,ICOOSC,IDC,
     &                  PCL,LUR,C2,NSMST,ISCALE,SCLFAC_R(IR))
          ELSE
C             WRITE(6,*) ' TTSS for C block skipped  '
C             CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
            CALL IFRMDS(LBL,-1,1,LUR)
            CALL SKPRCD2(LBL,-1,LUR)
            SCLFAC_R(IR) = 0.0D0
          END IF
*
*
          IF(NTEST.GE.100) THEN
            IF(INTERACT.EQ.1) THEN
              WRITE(6,*) ' TTSS for C block read in  '
              CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
            ELSE
              WRITE(6,*) ' TTSS for C block skipped  '
              CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
            END IF
          END IF
        END DO
*. Loop over L and R blocks in batches and obtain  contribution from
* given L and R blocks
          DO 10000 IIL = 1, NBLKL
            IL  = I1BATL(IBATCHL)-1+IIL       
          IF(SCLFAC_L(IL).NE.0.0D0) THEN
            IATP = IBLOCKL(1,IL)
            IBTP = IBLOCKL(2,IL)
            IASM = IBLOCKL(3,IL)
            IBSM = IBLOCKL(4,IL)
            IOFF = IBLOCKL(5,IL)
*
            NIA = NSSOA(IASM,IATP)
            NIB = NSSOB(IBSM,IBTP)
*. Possible permutations of L blocks
            PL = 1.0D0
            CALL PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PL,
     &           LATP,LBTP,LASM,LBSM,LSGN,LTRP,NLPERM)
            DO 9999 ILPERM = 1, NLPERM
C             write(6,*) ' Loop 9999 ILPERM = ', ILPERM
              IIASM = LASM(ILPERM)
              IIBSM = LBSM(ILPERM)
              IIATP = LATP(ILPERM)
              IIBTP = LBTP(ILPERM)
              NIIA = NSSOA(IIASM,IIATP)
              NIIB = NSSOB(IIBSM,IIBTP)
*
              IF(LTRP(ILPERM).EQ.1) THEN
                LROW = NSSOA(LASM(ILPERM-1),LATP(ILPERM-1))
                LCOL = NSSOB(LBSM(ILPERM-1),LBTP(ILPERM-1))
                CALL TRPMT3(SB(IOFF),LROW,LCOL,C2)
                CALL COPVEC(C2,SB(IOFF),LROW*LCOL)
               END IF
              IF(LSGN(ILPERM).EQ.-1)
     &        CALL SCALVE(SB(IOFF),-1.0D0,NIA*NIB)

              DO 9000 IIR = 1, NBLKR
                IR  = I1BATR(IBATCHR)-1+IIR       
              IF(SCLFAC_R(IR).NE.0.0D0) THEN
                JATP = IBLOCKR(1,IR)
                JBTP = IBLOCKR(2,IR)
                JASM = IBLOCKR(3,IR)
                JBSM = IBLOCKR(4,IR)
                JOFF = IBLOCKR(5,IR)
*
                NJA = NSSOA(JASM,JATP)
                NJB = NSSOB(JBSM,JBTP)
*
                IAEXC = ICONSPA(JATP,IIATP)
                IBEXC = ICONSPB(JBTP,IIBTP)
*
                IF(IAEXC.EQ.0.AND.JASM.NE.IIASM) IAEXC = 1
                IF(IBEXC.EQ.0.AND.JBSM.NE.IIBSM) IBEXC = 1
                IABEXC = IAEXC + IBEXC
*
                IF(IABEXC.LE.I12) THEN
                  INTERACT = 1
                ELSE
                  INTERACT = 0
                END IF
*
                IF(INTERACT.EQ.1) THEN
*. Possible permutations of this block
                   CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,
     &                  PSR,PL,RATP,RBTP,RASM,RBSM,RSGN,RTRP,
     &                  NRPERM)
*. Well, spin permutations are simple to handle
* if there are two terms just calculate and and multiply with
* 1+PSL*PSR
                     IF(NRPERM.EQ.1) THEN
                       FACTOR = 1.0D0
                     ELSE
                       FACTOR = 1.0D0 +PSL*PSR
                     END IF
                     SCLFAC = FACTOR*SCLFAC_L(IL)*SCLFAC_R(IR)
                     IF(INTERACT.EQ.1.AND.SCLFAC.NE.0.0D0) THEN
                     IF(NTEST.GE.20) THEN
                       WRITE(6,*) ' RSDNBB will be called for '
                       WRITE(6,*) ' L block : '
                       WRITE(6,'(A,5I5)') 
     &                 ' IIASM IIBSM IIATP IIBTP',
     &                   IIASM,IIBSM,IIATP,IIBTP
                       WRITE(6,*) ' R  block : '
                       WRITE(6,'(A,5I5)') 
     &                 ' JASM JBSM JATP JBTP',
     &                   JASM,JBSM,JATP,JBTP
                       WRITE(6,*) ' IOFF,JOFF ', IOFF,JOFF
                       WRITE(6,*) ' SCLFAC = ', SCLFAC
                     END IF
                     CALL GSDNBB2(I12,RHO1,RHO2,
     &                    IIASM,IIATP,IIBSM,IIBTP,
     &                    JASM,JATP,JBSM,JBTP,NGAS,
     &                    NELFSPGP(1,IOCTPA-1+IIATP),
     &                    NELFSPGP(1,IOCTPB-1+IIBTP),
     &                    NELFSPGP(1,IOCTPA-1+JATP),
     &                    NELFSPGP(1,IOCTPB-1+JBTP),
     &                    NAEL,NBEL,IAGRP,IBGRP,
     &                    SB(IOFF),CB(JOFF),C2,
     &                    ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                    MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,
     &                    SSCR,CSCR,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                    X,NSMOB,NSMST,NSMSX,NSMDX,
     &                    NIIA,NIIB,NJA,NJB,MXPOBS,
     &                    IPRNT,NACOB,RHO1S,SCLFAC,
     &                    S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &                    IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &                    NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,
     &                    IOBPTS_SEL,NINOB)
                          IF(NTEST.GE.500) THEN
                            write(6,*) ' Updated rho1 '
                            call wrtmat(rho1,nacob,nacob,nacob,nacob)
                          END IF

*
                END IF
                END IF
                END IF
 9000         CONTINUE
*. End of loop over R blocks in Batch
 9999     CONTINUE
*. Transpose or scale L block to restore order ??
          IF(LTRP(NLPERM+1).EQ.1) THEN
            CALL TRPMT3(SB(IOFF),NIB,NIA,C2)
            CALL COPVEC(C2,SB(IOFF),NIA*NIB)
          END IF
          IF(LSGN(NLPERM+1).EQ.-1)
     &    CALL SCALVE(SB(IOFF),-1.0D0,NIA*NIB)
*
          END IF
10000     CONTINUE
*. End of loop over L blocks in batch
 9001   CONTINUE
*.      ^ End of loop over batches of R blocks
10001 CONTINUE
*.    ^ End of loop over batches of L blocks
      CALL QEXIT('GASDN')
      RETURN
      END
      SUBROUTINE GSDNBB2(I12,RHO1,RHO2,
     &                  IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                  NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  NAEL,NBEL,
     &                  IJAGRP,IJBGRP,
     &                  SB,CB,C2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  X,NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,RHO1S,
     &                  SCLFAC,S2_TERM1,IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,
     &                  IDOSRHO2,RHO2AA,RHO2AB,RHO2BB,
     &                  NDACTORB,IDACTSPC,IDTFREORD,IFTDREORD,
     &                  IOBPTS_SEL,NINOB)
*
* Contributions to density matrix from sigma block (iasm iatp, ibsm ibtp ) and
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS : Number of As'es
* IAOC : Occpation of each AS for alpha strings in L
* IBOC : Occpation of each AS for beta  strings in L
* JAOC : Occpation of each AS for alpha strings in R
* JBOC : Occpation of each AS for beta  strings in R
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* MXPNGAS : Largest number of As'es allowed by program
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* Rho1, RHo2 : Updated density blocks 
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*. Input  
      DIMENSION CB(*),SB(*)
      INTEGER IDACTSPC(*),IDTFREORD(*),IFTDREORD(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*. Output
      DIMENSION RHO1(*),RHO2(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION  I1(*),XI1S(*),I2(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION C2(*)
*
      CALL QENTER('GSDNB')
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      NTESTO= NTEST
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSDNBB2 :  R block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' GSDNBB2 :  L block '
        WRITE(6,*) ' ==================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
*
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in L '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in L '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in R '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in R '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB
* 
        WRITE(6,*) 'SCLFAC =',SCLFAC
      END IF
*
C?    WRITE(6,*) ' IOBPTS_SEL in GSDNBB2'
C?    CALL IWRTMA(IOBPTS_SEL,NGAS,NSMOB,MXPNGAS,NSMOB)
      IACTIVE = 0
*
      IF(IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
*  beta contribution to RHO1
* =============================
*
C?      WRITE(6,*) ' GSBBD1 will be called (beta)'
        IAB = 2
        CALL GSBBD1(RHO1,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &       NGAS,IBOC,JBOC,
     &       SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB,
     &       NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?      WRITE(6,*) ' GSBBD1 was called '
C?      WRITE(6,*) ' Memory check '
C?      CALL MEMCHK
*
* ================================
* beta-beta contribution to RHO2
* ================================
*
        IF(I12.EQ.2.AND.NBEL.GE.2) THEN
C?        WRITE(6,*) ' GSBBD2A will be called (beta)'
          CALL GSBBD2A(RHO2,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &         NGAS,IBOC,JBOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,IDOSRHO2,RHO2BB,
     &         NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?        WRITE(6,*) ' GSBBD2A was called '
*
C              GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
C    &         NGAS,ISEL,ICEL,SB,CB,
C    &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
C    &         NOBPTS,IOBPTS,MAXI,MAXK,
C    &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
C    &         NSMOB,NSMST,NSMSX,MXPOBS)
        END IF
      END IF
*
      IF(IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* =============================
*  alpha contribution to RHO1
* =============================
*
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
C?        WRITE(6,*) ' GSBBD1 will be called (alpha)'
        IAB = 1
        CALL GSBBD1(RHO1,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &       NGAS,IAOC,JAOC,SB,CB,
     &       ADSXA,SXSTST,STSTSX,MXPNGAS,
     &       NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &       SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &       NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &       IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB,
     &         NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?        WRITE(6,*) ' GSBBD1 was called '
        IF(I12.EQ.2.AND.NAEL.GE.2) THEN
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
C?        WRITE(6,*) ' GSBBD2A will be called (alpha)'
          CALL GSBBD2A(RHO2,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &         NGAS,IAOC,JAOC,SB,CB,
     &         ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,IDOSRHO2,RHO2AA,
     &         NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?        WRITE(6,*) ' GSBBD2A was called '
        END IF
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
* ===================================
*  alpha-beta contribution to RHO2
* ===================================
*
      IF(I12.EQ.2.AND.NAEL.GE.1.AND.NBEL.GE.1) THEN
*. Routine uses transposed blocks
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
C?      WRITE(6,*) ' GSBBD2B will be called '
        IUSEAB = 0
        CALL GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                    JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                    IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                    SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                    NOBPTS,IOBPTS,MAXK,
     &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                    NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                    SSCR,CSCR,NACOB,NTEST,SCLFAC,S2_TERM1,
     &                    IDOSRHO2,RHO2AB,
     &                    NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
C?      WRITE(6,*) ' GSBBD2B was called '
     &                    
C     GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
C    &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
C    &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
C    &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
C    &                  NOBPTS,IOBPTS,MAXK,
C    &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
C    &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
C    &                  CJRES,SIRES,NORB,NTEST)
        CALL TRPMT3(CB,NJB,NJA,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
        CALL TRPMT3(SB,NIB,NIA,C2)
        CALL COPVEC(C2,SB,NIB*NIA)
      END IF
*
      CALL QEXIT('GSDNB')
      RETURN
      END
      SUBROUTINE GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,SB,CB,
     &                  ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,SCLFAC,
     &                  IDOSRHO2,RHO2SS,
     &                  NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
*
* Contributions to two-electron density matrix from column excitations
*
* GAS version, '96 , Jeppe Olsen 
*              Sept. 04, 2-electron spin-densities added 
*
* =====
* Input
* =====
* RHO2  : two body density matrix to be updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces 
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO2 : Updated density block
* RHO2SS : Updated alpha-alpha ( or beta-beta) 2e- spindensity (if IDOSRHO2=1)
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S, I2,XI2S : For holding creations/annihilations
*              type and symmetry
*
* Jeppe Olsen, Fall of 96     
*              Calculating densities only over selected spaces, Sept. 05
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST), SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
      INTEGER IDACTSPC(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*.Output
      DIMENSION RHO2(*), RHO2SS(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
C     INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSBBD2A in action '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
      IFRST = 1
      JFRST = 1
*
* Type of single excitations that connects the two column strings
      CALL DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,NGAS,ISEL,ICEL)
*.Symmetry of Double excitation that connects IBSM and JBSM
*. For general use : STSTSX => STSTDX
      IDXSM = STSTSX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0) GOTO 2001
      IF(NTEST.GE.1000)    
     &WRITE(6,*) ' ISCSM,ICCSM IDXSM ', ISCSM,ICCSM,IDXSM
      DO 2000 IDXTP =  1, NDXTP
        ITYP = ITP(IDXTP)
        JTYP = JTP(IDXTP)
        KTYP = KTP(IDXTP)
        LTYP = LTP(IDXTP)
        IF(IDACTSPC(ITYP)+IDACTSPC(JTYP)+IDACTSPC(KTYP)+IDACTSPC(LTYP)
     &     .NE.4) GOTO 2000
        IF(NTEST.GE.1000)
     &  write(6,*) ' ITYP JTYP KTYP LTYP ', ITYP,JTYP,KTYP,LTYP
        DO 1950 IKOBSM = 1, NSMOB
          JLOBSM = SXDXSX(IKOBSM,IDXSM)
          IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
          KFRST = 1
*. Loop over of symmetry of i orbitals
          DO 1940 ISM = 1, NSMOB   
          KSM = ADSXA(ISM,IKOBSM)
          NI = NOBPTS(ITYP,ISM)
          NK = NOBPTS(KTYP,KSM)
          IF(NI.EQ.0.OR.NK.EQ.0) GOTO 1940
*. Loop over batches of j orbitals 
          DO 1930 JSM = 1, NSMOB   
          LSM = ADSXA(JSM,JLOBSM)
          NJ = NOBPTS(JTYP,JSM)
          NL = NOBPTS(LTYP,LSM)
          IF(NJ.EQ.0.OR.NL.EQ.0) GOTO 1930
*
          IOFF = IOBPTS_SEL(ITYP,ISM)
          JOFF = IOBPTS_SEL(JTYP,JSM)
          KOFF = IOBPTS_SEL(KTYP,KSM)
          LOFF = IOBPTS_SEL(LTYP,LSM)
*
          IF(IOFF.LT.KOFF) GOTO 1930
          IF(JOFF.LT.LOFF) GOTO 1930
*
*
* =========================================================================
*                    Use N-2 projection method
* =========================================================================
*
              IFIRST = 1
*. Loop over batches of I strings
              NPART = NROW/MAXI
              IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
              IF(NTEST.GE.2000)
     &        write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
              DO 1801 IIPART = 1, NPART
                IBOT = 1+(IIPART-1)*MAXI
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                KBOT = 1- MAXK
                KTOP = 0
 1800           CONTINUE
                  KBOT = KBOT + MAXK
                  KTOP = KTOP + MAXK
*
* =========================================================
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
* =========================================================
*
                  IONE = 1
                  JLBOFF = 1
                  IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                    NJL = NJ*(NJ+1)/2
                    JLSM = 1
                  ELSE
                    NJL = NJ * NL
                    JLSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
                  CALL QENTER('ADADS')
                  II12 = 1
                  K12 = 1
                  IONE = 1
                  CALL ADADST_GAS(IONE,JSM,JTYP,NJ,
     &                            IONE,LSM,LTYP,NL,
     &                        ICCTP,ICCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        JFRST,KFRST,II12,K12,SCLFAC)
                  JFRST = 0
                  KFRST = 0
*
                  CALL QEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
                  CALL QENTER('MATCG')
                  DO  IJL = 1, NJL
                    CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                    I1JL = (L-1)*NJ+J
*. JAN28
                    IF(JLSM.NE.0) THEN
                      IJLE = J*(J-1)/2+L
                    ELSE
                      IJLE = IJL
                    END IF
*. JAN28
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
C                   JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                    JLOFF = (JLBOFF-1+IJLE-1)*NKBTC*NIBTC+1
                    IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,I1JL),XI1S(1,I1JL))
                    END IF
                  END DO
                  CALL QEXIT ('MATCG')
*
*
* =========================================================
*
*. obtain sb(KB,IA,ik) = sum(IB)<KB!a kb a ib !IB>S(IA,IB)
*
* =========================================================
*
                  IONE = 1
                  IKBOFF = 1
                  IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                    NIK = NI*(NI+1)/2
                    IKSM = 1
                  ELSE
                    NIK = NI * NK
                    IKSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
CT                CALL QENTER('ADADS')
                  II12 = 2
                  K12 = 1
                  IONE = 1
                  IF(IFRST.EQ.1) KFRST = 1
                  ONE = 1.0D0
                  CALL ADADST_GAS(IONE,ISM,ITYP,NI,
     &                            IONE,KSM,KTYP,NK,
     &                        ISCTP,ISCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        IFRST,KFRST,II12,K12,ONE   )
                  IFRST = 0
                  KFRST = 0
*
CT                CALL QEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  I = 0
                  K = 1
*
CT                CALL QENTER('MATCG')
                  DO  IIK = 1, NIK
                    CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                    I1IK = (K-1)*NI+I
*. JAN28
                    IF(IKSM.NE.0) THEN
                      IIKE = I*(I-1)/2+K
                    ELSE
                      IIKE = IIK
                    END IF
*. JAN28
*.SB(IA,KB,ik) = +/-S(IA,a+ka+iIA)
C                   IKOFF = (IKBOFF-1+IIK-1)*NKBTC*NIBTC+1
                    IKOFF = (IKBOFF-1+IIKE-1)*NKBTC*NIBTC+1
                    IF(IKSM.EQ.1.AND.I.EQ.K) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(SSCR(IKOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,I1IK),XI1S(1,I1IK))
                    END IF
                  END DO
CT                CALL QEXIT ('MATCG')
*
*
* =================================================================
*
* RHO2C(ik,jl)  = RHO2C(ik,jl) - sum(Ia,Kb)SB(Ia,Kb,ik)*CB(Ia,Kb,jl)
*
* =================================================================
*
* The minus ??
*
* Well, the density matrices are constructed as 

* <I!a+i a+k aj al!> = -sum(K) <I!a+ia+k!K><J!aj al!K>, and
* the latter matrices are the ones we are constructing
*
              IOFF = IOBPTS_SEL(ITYP,ISM)
              JOFF = IOBPTS_SEL(JTYP,JSM)
              KOFF = IOBPTS_SEL(KTYP,KSM)
              LOFF = IOBPTS_SEL(LTYP,LSM)
              NTESTO = NTEST
C?            IF(IOFF.EQ.3.AND.JOFF.EQ.3.AND.KOFF.EQ.4.AND.LOFF.EQ.4)
C?   &            NTEST = 5000
                  LDUMMY = NKBTC*NIBTC
                  IF(NTEST.GE.2000) THEN
                    WRITE(6,*) ' CSCR matrix '
                    CALL WRTMAT(CSCR,LDUMMY,NJL,LDUMMY,NJL)
                    WRITE(6,*) ' SSCR matrix '
                    CALL WRTMAT(SSCR,LDUMMY,NIK,LDUMMY,NIK)
                  END IF
           
                  IF(IFIRST.EQ.1) THEN
                    FACTOR = 0.0D0
                  ELSE 
                    FACTOR = 1.0D0
                  END IF
C                 MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
                  LDUMMY = NKBTC*NIBTC
                  ONEM = -1.0D0
                  CALL MATML7(X,SSCR,CSCR,NIK,NJL,
     &                        LDUMMY,NIK,LDUMMY,NJL,
     &                        FACTOR,ONEM,1)
                  IFIRST = 0
                  IF(NTEST.GE.2000) THEN
                    WRITE(6,*) ' Updated X matrix'
                    CALL WRTMAT(X,NIK,NJL,NIK,NJL)
                  END IF

*
                IF(KEND.EQ.0) GOTO 1800
*. End of loop over partitionings of resolution strings
 1801         CONTINUE
*. Rho2(ik,jl) has been constructed for ik,jl belonging to
*. Scatter out to density matrix
              IOFF = IOBPTS_SEL(ITYP,ISM)
              JOFF = IOBPTS_SEL(JTYP,JSM)
              KOFF = IOBPTS_SEL(KTYP,KSM)
              LOFF = IOBPTS_SEL(LTYP,LSM)
C?            WRITE(6,*) ' ITYP, ISM, IOFF = ', ITYP, ISM, IOFF
C?            WRITE(6,*) ' JTYP, JSM, JOFF = ', JTYP, JSM, JOFF
C?            WRITE(6,*) ' KTYP, KSM, KOFF = ', KTYP, KSM, KOFF
C?            WRITE(6,*) ' LTYP, LSM, LOFF = ', LTYP, LSM, LOFF
              CALL ADTOR2(RHO2,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,
     &                    NDACTORB)
C                  ADTOR2(RHO2,RHO2T,ITYPE,
C    &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
               IF(IDOSRHO2.EQ.1) THEN
*. and add to spin-density 
              CALL ADTOR2S(RHO2SS,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,
     &                    NDACTORB)
               END IF

 1930       CONTINUE
 1940     CONTINUE
 1950   CONTINUE
 2000 CONTINUE
 2001 CONTINUE
*
      RETURN
      END
      SUBROUTINE GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXK,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  CJRES,SIRES,NORB,NTESTG,SCLFAC,S2_TERM1,
     &                  IDOSRHO2,RHO2AB,
     &                  NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
*
* alpha-beta contribution to two-particle density matrix 
* from given c-block and s-block.
*
* S2_TERM1 = - <L!a+i alpha a+jbeta a i beta a j alpha !R>
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
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* X : Space for block of two-electron integrals
*
* Jeppe Olsen, Fall of 1996
*
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
      INTEGER IOBPTS_SEL(MXPNGAS,*), IDACTSPC(NGAS)
*.Input
      DIMENSION CB(*),SB(*)
*. Output
      DIMENSION RHO2(*),RHO2AB(*)
*.Scratch
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION CJRES(*),SIRES(*)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
*
      CALL QENTER('GSD2B')
      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' GSBBD2B speaking '
        WRITE(6,*) ' ================== '
      END IF
C?    WRITE(6,*) ' NJAS NJB = ',NJA,NJB
C?    WRITE(6,*) ' IAGRP IBGRP = ', IAGRP,IBGRP
C?    WRITE(6,*) ' MXPNGAS = ', MXPNGAS
C?    WRITE(6,*) ' NSMOB = ', NSMOB
      IROUTE = 3
*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC)
      CALL SXTYP_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC)           
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        IF(IDACTSPC(ITYP)+IDACTSPC(JTYP).NE.2) GOTO 2001
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS_SEL(ITYP,ISM)
          JOFF = IOBPTS_SEL(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all Ka strings
*. a+j!ka> = +/-/0 * !Ja>
          CALL ADSTN_GAS(JSM,JTYP,JATP,JASM,IAGRP,
     &                   I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,
     &                   SCLFAC)
*. a+i!ka> = +/-/0 * !Ia>
          ONE    = 1.0D0
          CALL ADSTN_GAS(ISM,ITYP,IATP,IASM,IAGRP,
     &                   I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,
     &                   ONE   )
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
C             COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
              CALL COMPRS2LST(I1,XI1S,NJ,I3,XI3S,NI,NKASTR,NKAEFF)
          ELSE 
              NKAEFF = NKASTR
          END IF
            
*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK   
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
          DO 1801 IKABTC = 1, NKABTC
C?          write(6,*) ' Batch over kstrings ', IKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, NJ
              CALL GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*. Obtain S(ka,i,Ib) for Ka in batch
            DO II = 1, NI
              CALL GET_CKAJJB(SB,NI,NIA,SIRES,LKABTC,NIB,
     &             II,I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
              IF(IDACTSPC(KTYP)+IDACTSPC(LTYP).NE.2) GOTO 2000
*
              DO 1930 KSM = 1, NSMOB
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
C?              WRITE(6,*) ' Loop 1930, KSM LSM ',KSM,LSM
                KOFF = IOBPTS_SEL(KTYP,KSM)
                LOFF = IOBPTS_SEL(LTYP,LSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0  
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP) IKORD=1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
                ONE = 1.0D0
                CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE)
                IF(NKBSTR.EQ.0) GOTO 1930
*
*. Update two-electron density matrix
*  Rho2b(ij,kl) =  Sum(ka)S(Ka,i,Ib)<Ib!Eb(kl)!Jb>C(Ka,j,Jb)
*
                ZERO = 0.0D0
                CALL SETVEC(X,ZERO,NI*NJ*NK*NL)
*
C               WRITE(6,*) ' Before call to ABTOR2'
                CALL ABTOR2(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,X,NI,NJ,NK,NL,NKBSTR,
     &               I4,XI4S,I2,XI2S,IKORD)
*. contributions to Rho2(ij,kl) has been obtained, scatter out
C?              WRITE(6,*) ' Before call to ADTOR2'
C?              WRITE(6,*) ' RHO2B (X) matrix '
C?              call wrtmat(x,ni*nj,nk*nl,ni*nj,nk*nl)
*. Contribution to S2
                IF(KTYP.EQ.JTYP.AND.KSM.EQ.JSM.AND.
     &            ITYP.EQ.LTYP.AND.ISM.EQ.LSM) THEN
                  DO I = 1, NI
                    DO J = 1, NJ
                      IJ = (J-1)*NI+I
                      JI = (I-1)*NJ+J
                      NIJ = NI*NJ
                      S2_TERM1 = S2_TERM1-X((JI-1)*NIJ+IJ)
                    END DO
                  END DO
                END IF
         
     &             
C?            WRITE(6,*) ' ITYP, ISM, IOFF = ', ITYP, ISM, IOFF
C?            WRITE(6,*) ' JTYP, JSM, JOFF = ', JTYP, JSM, JOFF
C?            WRITE(6,*) ' KTYP, KSM, KOFF = ', KTYP, KSM, KOFF
C?            WRITE(6,*) ' LTYP, LSM, LOFF = ', LTYP, LSM, LOFF
                CALL ADTOR2(RHO2,X,2,
     &                NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NDACTORB)
                IF(IDOSRHO2.EQ.1) THEN
                  CALL ADTOR2S(RHO2AB,X,2,
     &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NDACTORB)
                END IF 
C?              write(6,*) ' updated density matrix '
C?              call prsym(rho2,NDACTORB*NDACTORB)

 1930         CONTINUE
 2000       CONTINUE
 1801     CONTINUE
*. End of loop over partitioning of alpha strings
 1940   CONTINUE
 2001 CONTINUE
*
 9999 CONTINUE
*
*
      CALL QEXIT('GSD2B')
      RETURN
      END
      SUBROUTINE GSBBD1(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,
     &                  IUSE_PH,IPHGAS,IDOSRHO1,SRHO1,IAB,
     &                  NDACTORB,IDACTSPC,IOBPTS_SEL,NINOB)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen 
* Particle-Hole version of Jan. 98
* Active orbital spaces added Sept. 05
*
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces 
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1 : Updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), ITSOB(*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
      INTEGER IDACTSPC(*)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*.Output
      DIMENSION RHO1(*), SRHO1(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(*),XI1S(*)
      DIMENSION I2(*),XI2S(*)
*.Local arrays ( assume MPNGAS = 16 ) !!! 
      DIMENSION ITP(16*16),JTP(16*16)
*
      DIMENSION IJ_REO(2),IJ_DIM(2),IJ_SM(2),IJ_TP(2),IJ_AC(2)
      DIMENSION IJ_OFF(2)
      DIMENSION ISCR(2)
      DIMENSION ICGRP(16),ISGRP(16)
*. Add or subtract for spindensity
      IF(IAB.EQ.1) THEN
        XAB = 1.0D0
      ELSE
        XAB = -1.0D0
      END IF
*.Local arrays
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GSBBD1 in action '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' ISCSM, ICCSM = ', ISCSM, ICCSM
*         
        WRITE(6,*) ' GSBBD1, sclfac ',SCLFAC
      END IF
*
      IFRST = 1
      JFRST = 1
*. Number of partitionings over column strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Groups defining supergroups
C          GET_SPGP_INF(ISPGP,ITP,IGRP)
      CALL GET_SPGP_INF(ICCTP,IGRP,ICGRP)
      CALL GET_SPGP_INF(ISCTP,IGRP,ISGRP)

* Type of single excitations that connects the two column strings
      CALL SXTYP2_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL,IPHGAS)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(NTEST.GE.1000)    
     &WRITE(6,*) ' ISCSM,ICCSM IJSM ', ISCSM,ICCSM,IJSM
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(IDACTSPC(ITYP)+IDACTSPC(JTYP).NE.2) GOTO 900
        IF(NTEST.GE.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
*. Hvilken vej skal vi valge, 
        NOP = 2
        IJ_AC(1) = 2
        IJ_AC(2) = 1
        IJ_TP(1) = ITYP
        IJ_TP(2) = JTYP
        IF(IUSE_PH.EQ.1) THEN
          CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TP,IJ_AC,IJ_REO,SIGNIJ)
        ELSE
          IJ_REO(1) = 1
          IJ_REO(2) = 2
          SIGNIJ = 1.0D0
        END IF
*
        ISCR(1) = IJ_AC(1)
        ISCR(2) = IJ_AC(2)
        IJ_AC(1) = ISCR(IJ_REO(1))
        IJ_AC(2) = ISCR(IJ_REO(2))
*
        ISCR(1) = ITYP
        ISCR(2) = JTYP
        IJ_TP(1) = ISCR(IJ_REO(1))
        IJ_TP(2) = ISCR(IJ_REO(2))

        DO 800 ISM = 1, NSMOB
*. new i and j so new intermediate strings
          KFRST = 1
*
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.1000) write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
          IBIORB = IOBPTS_SEL(ITYP,ISM)
          IBJORB = IOBPTS_SEL(JTYP,JSM)
*. Reorder 
* 
          ISCR(1) = ISM
          ISCR(2) = JSM
          IJ_SM(1) = ISCR(IJ_REO(1))
          IJ_SM(2) = ISCR(IJ_REO(2))
*
          ISCR(1) = NIORB
          ISCR(2) = NJORB
          IJ_DIM(1) = ISCR(IJ_REO(1))
          IJ_DIM(2) = ISCR(IJ_REO(2))
*
          ISCR(1) = IBIORB
          ISCR(2) = IBJORB
          IJ_OFF(1) = ISCR(IJ_REO(1))
          IJ_OFF(2) = ISCR(IJ_REO(2))
*

          IF(NTEST.GE.2000)
     &    WRITE(6,*) ' NIORB NJORB ', NIORB,NJORB
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          SCLFACS = SCLFAC*SIGNIJ
          IF(NTEST.GE.1000) 
     &    WRITE(6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
          CALL ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC(1))
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          ONE = 1.0D0
          CALL ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,
     &         I2,XI2S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC(1))
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I2,XI2S,IJ_DIM(1),
     &             NKASTR,NKAEFF)
          ELSE 
              NKAEFF = NKASTR
          END IF
C         WRITE(6,*) ' NKAEFF NKASTR', NKAEFF,NKASTR

*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = MIN(KTOP + MAXK,NKAEFF)
              IF(KTOP.EQ.NKAEFF) THEN
                KEND = 1
              ELSE
                KEND = 0
              END IF
              LKABTC = KTOP - KBOT +1

*. This is the place to start over partitioning of I strings
              DO 701 IIPART = 1, NIPART
                IBOT = (IIPART-1)*MAXI+1
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP - IBOT + 1
* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
                DO JJORB = 1,IJ_DIM(2)
                  ICGOFF = 1 + (JJORB-1)*LKABTC*NIBTC
                  CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                 LKABTC,I1(KBOT+(JJORB-1)*NKASTR),
     &                 XI1S(KBOT+(JJORB-1)*NKASTR) )
                END DO
* Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
                DO IIORB = 1,IJ_DIM(1)
*.Gather S Block
                  ISGOFF = 1 + (IIORB-1)*LKABTC*NIBTC
                  CALL MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,
     &                   LKABTC,I2(KBOT+(IIORB-1)*NKASTR),
     &                   XI2S(KBOT+(IIORB-1)*NKASTR) )
                END DO  
*
                NKI = LKABTC*NIBTC
                IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' CSCR and SSCR '
                 CALL WRTMAT(CSCR,IJ_DIM(2),NKI,IJ_DIM(2),NKI)
                 CALL WRTMAT(SSCR,IJ_DIM(1),NKI,IJ_DIM(1),NKI)
                END IF
*
*. And then the hard  work
                FACTORC = 0.0D0
                FACTORAB = 1.0D0
                CALL MATML7(RHO1S,SSCR,CSCR,IJ_DIM(1),IJ_DIM(2),NKI,
     &               IJ_DIM(1),NKI,IJ_DIM(2),FACTORC,FACTORAB,1)
*
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' Block to one-body density '
                  CALL WRTMAT(RHO1S,IJ_DIM(1),IJ_DIM(2),
     &                              IJ_DIM(1),IJ_DIM(2))
                END IF
*. Scatter out to complete matrix  
                DO JJORB = 1, IJ_DIM(2)
                  JORB = IJ_OFF(2)-1+JJORB
                  DO IIORB = 1, IJ_DIM(1)
                    IORB = IJ_OFF(1)-1+IIORB
                    RHO1((JORB-1)*NDACTORB+IORB) =
     &              RHO1((JORB-1)*NDACTORB+IORB) +
     &              RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)
                    IF(IDOSRHO1.EQ.1) THEN  
                      SRHO1((JORB-1)*NDACTORB+IORB) =
     &                SRHO1((JORB-1)*NDACTORB+IORB) +
     &                XAB*RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)   
                    END IF
                  END DO
                END DO
*               /\ End of hard work

  701     CONTINUE
*. /\ end of this I partitioning  
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
*. End of loop over I partitioninigs
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
C!    stop ' enforrced stop in RSBBD1 '
      RETURN
      END
      FUNCTION IB_H1(ISM,IHSM,NR,NC)
*
*. Offset to symmetryblock H(ISM,*)
*. Jeppe Olsen, Dec. 2000
*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
*
      INTEGER NR(*),NC(*)
*
      IB = 1
      DO IISM = 1, ISM - 1
        JJSM = MULTD2H(IISM,IHSM)
        IB = IB + NR(IISM)*NC(JJSM)
      END DO
*
      IB_H1 = IB
*
      RETURN
      END
      SUBROUTINE MULT_H1H2(H1,IH1SM,H2,IH2SM,H12,IH12SM)
*. Two set of one-electron integrals H1 and H2 are 
*. given as symmetrypacked complete quadratic matrices. 
*. Obtain product as symmetrypacked complete quadratic matrix
*
* Jeppe Olsen, Dec. 2000
*
      INCLUDE 'implicit.inc'
*.Input
      INCLUDE 'mxpdim.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'lucinp.inc'
*
      DIMENSION H1(*),H2(*)
*. Output
      DIMENSION H12(*)
*. 
      INCLUDE 'multd2h.inc'
*
C?    WRITE(6,*) ' Entering MULT_H1H2 '
C?    WRITE(6,*)' IH1SM, IH2SM = ', IH1SM, IH2SM
      IH12SM = MULTD2H(IH1SM,IH2SM)
*. Loop over symmetry blocks of H1H2
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IH12SM)
*. Connecting symmetry A(ISM,KSM)B(KSM,JSM)
        KSM = MULTD2H(ISM,IH1SM)
C?      WRITE(6,*) ' ISM, JSM, KSM =', ISM, JSM, KSM
*. Offsets to A(ISM,KSM) and B(KSM,JSM)
C              IB_H1(ISM,IHSM,NR,NC)
        IB_A = IB_H1(ISM,IH1SM,NTOOBS,NTOOBS)
        IB_B = IB_H1(KSM,IH2SM,NTOOBS,NTOOBS)
        IB_AB = IB_H1(ISM,IH12SM,NTOOBS,NTOOBS)
C?      WRITE(6,*) ' IB_A, IB_B, IB_AB =', IB_A,IB_B,IB_AB
        
*
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)
        NK = NTOOBS(KSM)
*
        FACTORC = 0.0D0
        FACTORAB = 1.0D0
        ZERO = 0.0D0
        CALL SETVEC(H12(IB_AB),ZERO,NI*NJ)
        CALL MATML7(H12(IB_AB),H1(IB_A),H2(IB_B),NI,NJ,NI,NK,NK,NJ,
     &              FACTORC,FACTORAB,0)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output matrix from MULTH1H2 '
        CALL PRHONE(H12,NTOOBS,IH12SM,NSMOB,0)
      END IF
*
      RETURN
      END
      SUBROUTINE ABEXP(A,IASM,B,IBSM,AB)
*
* Evaluate expectation value of product of 
* two one-electron operators 
*
* <0!AB!0> = sym(ijkl) A(ij)B(kl) d(ijkl) + sum(ij) rho1(ij) (AB)(ij)
*
* Jeppe Olsem Dec. 4, 2000 in Helsingfors 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*. Input, A and B are expexted to included in symmetrypacked 
*. complete  form
      DIMENSION A(*),B(*) 
*
      INCLUDE 'multd2h.inc'
*
C?    CALL CHK_D2
C?    CALL ETWO
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUMN,'MARK  ',IDUM,'ABEXP ')
*. Largest number of orbitals of given sym
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Allocate memory     
      LB_RHO2 = MXSOB**4
      CALL MEMMAN(KLRHO2B,LB_RHO2,'ADDL  ',2,'RHO2B ')
      L_VEC = MXSOB**2
      CALL MEMMAN(KLVEC,L_VEC,'ADDL  ',2,'VEC   ')
*
*. Two-electron contributions
*
      D2SUM = 0.0D0
      AB2 = 0.0D0
      IAB_SM = MULTD2H(IASM,IBSM)
      DO ISM = 1,  NSMOB
        JSM = MULTD2H(ISM,IASM)
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)  
        NIJ = NI*NJ
*. Offset for symmetryblock A(ISM,JSM)
        IA_OFF = 1
        DO IISM = 1, ISM-1 
          JJSM = MULTD2H(IASM,IISM)
          IA_OFF = IA_OFF + NTOOBS(IISM)*NTOOBS(JJSM)
        END DO
        DO KSM = 1, NSMOB
          LSM = MULTD2H(KSM,IBSM)
          NK = NTOOBS(KSM)
          NL = NTOOBS(LSM)
          NKL = NK*NL
*. Offset for symmetryblock B(KSM,LSM)
          IB_OFF = 1
          DO KKSM = 1, KSM - 1
            LLSM = MULTD2H(IBSM,KKSM)
            IB_OFF = IB_OFF + NTOOBS(KKSM)*NTOOBS(LLSM)
          END DO 
*. Fetch RHO2(ISM,JSM,KSM,LSM)
C         GETD2(RHO2B,ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,ICOUL)
          CALL GET_D2_SMBLK(WORK(KLRHO2B),ISM,JSM,KSM,LSM)
* sum(kl) RHO2(ij,kl) B(kl)
          FACTORC = 0.0D0
          FACTORAB = 1.0D0
          ZERO = 0.0D0
          CALL SETVEC(WORK(KLVEC),ZERO,NIJ)
          CALL MATML7(WORK(KLVEC),WORK(KLRHO2B),B(IB_OFF),NIJ,1,
     &                NIJ,NKL,NKL,1,FACTORC,FACTORAB,0)
*. sum(ij) (A(ij) (sum(kl) RHO2(ij,kl)B(kl) )
          AB2 = AB2 + INPROD(A(IA_OFF),WORK(KLVEC),NIJ)
        END DO
      END DO
C?    WRITE(6,*) ' Two-electron contribution ', AB2
*
*. One-electron part
*
* AB(IJ) RHO1(IJ)
C     MULT_H1H2(H1,IH1SM,H2,IH2SM,H12,IH12SM)
*. Yes there is space for a matrix in KLVEC
      CALL MULT_H1H2(A,IASM,B,IBSM,WORK(KLVEC),IABSM)
      AB1 = 0.0D0
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IABSM)
        NI = NTOOBS(ISM) 
        NJ = NTOOBS(JSM)
        NIJ = NI*NJ
*. Obtain density RHO1(ISM,JSM)
        CALL GET_D1_SMBLK(WORK(KLRHO2B),ISM,JSM)
        IH12_OF = IB_H1(ISM,IABSM,NTOOBS,NTOOBS)
C                 IB_H1(ISM,IHSM,NR,NC) 
        AB1 = AB1 + INPROD(WORK(KLRHO2B),
     &                     WORK(KLVEC-1+IH12_OF),NIJ)
      END DO
*
      AB = AB1 + AB2
*
      CALL MEMMAN(IDUM,IDUMN,'FLUSM ',IDUM,'ABEXP ')
      RETURN
      END
      SUBROUTINE CHK_D2
*
* Check trace of rho2 as delivered by SMBLK routine
*
* Jeppe Olsem Dec. 4, 2000 in Helsingfors 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      INCLUDE 'multd2h.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUMN,'MARK  ',IDUM,'ABEXP ')
*. Largest number of orbitals of given sym
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Allocate memory     
      LB_RHO2 = MXSOB**4
      CALL MEMMAN(KLRHO2B,LB_RHO2,'ADDL  ',2,'RHO2B ')
*
      D2SUM = 0.0D0
      DO ISM = 1, NSMOB
        DO KSM = 1, NSMOB
        NI = NTOOBS(ISM)
        NK = NTOOBS(KSM)
*. Fetch RHO2(ISM,JSM,KSM,LSM)
          CALL GET_D2_SMBLK(WORK(KLRHO2B),ISM,ISM,KSM,KSM)
*. Add to <0!E(ii)E(kk)-delta (ik)E (ik)!0>
            DO I = 1, NI
            DO K = 1, NK
             IIKK = (K-1)*NK*NI*NI
     &            + (K-1)*   NI*NI
     &            + (I-1)*      NI
     &            +  I
             D2SUM = D2SUM + WORK(KLRHO2B-1+IIKK)
           END DO
           END DO
        END DO
      END DO
      WRITE(6,*) ' Check of rho2 sum rule  = ', D2SUM
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ABEXP ')
      RETURN
      END
      SUBROUTINE ETWO
*
* Obtain 2-el contributions to energy from 
* SMBLK density 
*
* Jeppe Olsem Dec. 4, 2000 in Helsingfors 
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      REAL*8 INPROD
*
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
*
      INCLUDE 'multd2h.inc'
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUMN,'MARK  ',IDUM,'ABEXP ')
*. Largest number of orbitals of given sym
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Allocate memory     
      LB_RHO2 = MXSOB**4
      CALL MEMMAN(KLRHO2B,LB_RHO2,'ADDL  ',2,'RHO2B ')
*
      E2 = 0.0D0
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, NSMOB
        IJSM = MULTD2H(ISM,JSM)
        IJKSM = MULTD2H(IJSM,KSM)
        LSM = IJKSM
*
        NI = NTOOBS(ISM)
        NJ = NTOOBS(JSM)
        NK = NTOOBS(KSM)
        NL = NTOOBS(LSM)
*. Fetch RHO2(ISM,JSM,KSM,LSM)
        CALL GET_D2_SMBLK(WORK(KLRHO2B),ISM,JSM,KSM,LSM)
        IJKL = 0
        DO L = 1, NL
        DO K = 1, NK
        DO J = 1, NJ
        DO I = 1, NI
*
          IABS = IBSO(ISM)-1+I
          JABS = IBSO(JSM)-1+J
          KABS = IBSO(KSM)-1+K
          LABS = IBSO(LSM)-1+L
*
          IREO = IREOST(IABS)
          JREO = IREOST(JABS)
          KREO = IREOST(KABS)
          LREO = IREOST(LABS)
*
         XIJKL = GTIJKL(IREO,JREO,KREO,LREO) 
         IJKL = IJKL + 1
         E2= E2 + 0.5D0*WORK(KLRHO2B-1+IJKL)*XIJKL
         
        END DO
        END DO
        END DO 
        END DO
*       ^ End of loop over orbitals
      END DO
      END DO
      END DO
*     ^ End of loop over symmetries
      WRITE(6,*) ' Two-electron contribution to energy = ', E2
      RETURN
      END 
      SUBROUTINE PICO4(VEC1,VEC2,LU1,LU2,LU3,LU4,RNRM,EIG,FINEIG,MAXIT,
     &                 NBATCH,LLBATCHB,LLBATCHE,LBLOCK,IBLOCK,IPRTXX,
     &                 NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                 THRES_ET,THRES_EC,THRES_CC,
     &                 E_CONV,C_CONV,ICLSSEL,
     &                 IBLK_TO_CLS,NCLS,CLS_C,CLS_E,CLS_CT,CLS_ET,
     &                 ICLS_A,ICLS_L,RCLS_L,IBLKS_A,
     &                 CLS_DEL,CLS_DELT,ISKIPEI,
     &                 I2BLK,VEC3,ICLS_A2,MXLNG,IBASSPC,EBASC,CBASC,
     &                 NSPC,IMULTGR,IPAT,LPAT,IREFSPC,CONVER)
*
* Davidson algorithm , requires three blocks in core
*
* Only three vectors in on DISC
*
* Lu4 should only hold a batch of coefficients 
*
*
*
* Block processing version
*
* Jeppe Olsen Winter of 1996
*
* Revision of june 97, PICO2 => PICO3
*      modifications : Dynamic construction of of  batches
*                      Only relevant sigma blocks constructed
*              Oct, 98 : Info on base spaces added  
*              Aug  02 : RCLS_L Added 
*              June 03 : Microiterations in active subspace added

*
* Initial version - Only Diagonal preconditioner,
*
* Special version for NROOT = 1, MAXVEC = 2 !!
*
* Input :
* =======
*        LU1 : Initial  vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU2,LU3   : Scatch files
*        MAXIT     : Largest allowed number of iterations
*        NBATCH    : Number of batches of vector
*        LBATCHB   : Number of blocks in each batch
*        LBATCHE   : Number of elements  in each batch
*        IBLOCK    : Some additional informaition about the blocking 
*                    that this routine does not care about !!!!
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
*   THRES_ET   : Convergence criteria for eigenvalue
*
*   THRES_EC   : Threshold for second order energies for individual terms
*   THRES_CC   : Threshold for first  order wavefunction  for individual terms
*                
*
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL * 8   INPROD, INPRDD, INPRODB
*
      DIMENSION LLBATCHB(NBATCH),LLBATCHE(NBATCH),IBLOCK(*)
      DIMENSION LBLOCK(*)
*. I2BLK should atleast have length of number of
*. blocks 
      DIMENSION I2BLK(*)
*
      DIMENSION RNRM(MAXIT,1),EIG(MAXIT,1),FINEIG(1)
      DIMENSION H0(*),IPNTR(1)
      DIMENSION H0SCR(*)
*. Class and block information 
      DIMENSION IBLK_TO_CLS(*) 
      DIMENSION CLS_C(NCLS),CLS_E(NCLS),CLS_CT(NCLS),CLS_ET(NCLS)
      DIMENSION CLS_DEL(*), CLS_DELT(*)
*. Base CI spaces : CI space where a given class is introduced 
      DIMENSION IBASSPC(*),EBASC(*),CBASC(*)
      DIMENSION IPAT(*)
*.Initial VEC3      
      DIMENSION VEC3(*)
      INTEGER ICLS_A(NCLS), ICLS_L(NCLS),IBLKS_A(*),ICLS_A2(NCLS)
      DIMENSION RCLS_L(NCLS)
*
*     H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
      LOGICAL CONVER
*
C?    WRITE(6,*) ' Memchk at start of PICO3'
C?    CALL MEMCHK
      IPICO = 0
      IF(IPICO.NE.0) THEN
C?      WRITE(6,*) ' Perturbative solver '
C       MAXVEC = MIN(MAXVEC,2)
      ELSE IF(IPICO.EQ.0) THEN
C?      WRITE(6,*) ' Variational  solver '
      END IF
*
      WRITE(6,*) ' Number of spaces ', NSPC
C     WRITE(6,*) ' Map : Class => Base space '
C     CALL IWRTMA(IBASSPC,1,NCLS,1,NCLS)
      IF(ICLSSEL.EQ.1) THEN
        WRITE(6,*) ' Class selection will be performed '
        WRITE(6,*) ' Number of classes ', NCLS
C?      WRITE(6,*) ' Dimension of each class ( Integer )'
C?      CALL IWRTMA(ICLS_L,1,NCLS,1,NCLS)
C?      WRITE(6,*) ' Dimension of each class ( Real )'
C?      CALL WRTMAT(RCLS_L,1,NCLS,1,NCLS)
      END IF
*
      IF(IMULTGR.NE.0) THEN
        WRITE(6,*) ' Multispace method in use '
        WRITE(6,*)
        WRITE(6,*) ' Length of pattern ', LPAT
        WRITE(6,*) ' Pattern : '
        CALL IWRTMA(IPAT,1,LPAT,1,LPAT)
        WRITE(6,*) 
        WRITE(6,*) ' Reference space ', IREFSPC
      END IF
*
      IPRT = 10
      IOLSTM = 0
      IF(IPRT.GT.10.AND.IOLSTM.NE.0)
     &WRITE(6,*) ' Inverse iteration modified Davidson '
      IF(IPRT.GT.10.AND.IOLSTM.EQ.0)
     &WRITE(6,*) ' Normal Davidson method '
*
C?    WRITE(6,*) ' LU1 LU2 LU3 LU4 = ', LU1,LU2,LU3,LU4
      IF(IPRT.GE.20) THEN
        WRITE(6,*) ' Convergence threshold for eigenvalues', THRES_ET
        WRITE(6,*)
        WRITE(6,*) ' Elements of trial vectors discarded if '
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        WRITE(6,*) 
     &  '    Estimate of contribution to wavefunction is less than ',
     &  THRES_CC
        WRITE(6,*) 
     &  '    Estimate of contribution to Energy is less than ',
     &  THRES_EC
      END IF
      WRITE(6,*)
      IF(IPRT.GE.10)
     &WRITE(6,*) ' Max number of batches of vector ', NBATCH
*
*. Total number of blocks
      NBLOCKT = 0
      DO IBATCH = 1, NBATCH
        NBLOCKT = NBLOCKT + LLBATCHB(IBATCH)
      END DO
      IF(IPRT.GE.10)
     &WRITE(6,*) ' Max number of blocks ', NBLOCKT
      WRITE(6,*)
      TEST = 1.0D-6
      CONVER = .FALSE.
      IROOT = 1
      NROOT = 1
      ZERO = 0.0D0
      MAX_MICRO_IT = 3
*.    ^ Should be moved outside 
*
*. Play around with dynamic allocation of batches ...
      IDYNBATCH = 1
*
* ===================
*.Initial iteration
* ===================
*
      IF(MAXIT.EQ.0) THEN
        WRITE(6,*) ' Max number of iterations is zero'
        WRITE(6,*) ' I will just return from PICO3'
        RETURN
      END IF
*
      ITER = 1
      IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==============================='      
        WRITE(6,*) ' Info from iteration .... ', ITER
        WRITE(6,*) ' ==============================='      
        WRITE(6,*)
      END IF
      CALL QENTER('INIIT')
*. Obtain energy of initial vector
      IF(ISKIPEI.EQ.0) THEN
*. active classses of initial vector
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU1,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
*. Mark active blocks and find required number of batches
          CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
        ELSE
          NBATCHL = NBATCH
        END IF
*
        IF(IPRT.GE.10) WRITE(6,*) ' Number of batches for C ', NBATCHL
        E = 0.0D0
        IRESTRICT = 1
        IOFF = 1
        DO IBATCH = 1, NBATCHL
          LBATCHB = LLBATCHB(IBATCH)
          LBATCHE = LLBATCHE(IBATCH)
          IF(IPRT.GE.10) 
     &    WRITE(6,*) '  <Delta ! H ! Delta >, Batch : ',IBATCH
          CALL SETVEC(VEC1,ZERO,LBATCHE)
          CALL SBLOCK(LBATCHB,IBLOCK,IOFF,VEC2,VEC1,LU1,IRESTRICT,0,
     &                0,0,0,0)
          IF(IPRT.GE.200) THEN
            WRITE(6,*) ' Initial batch of S, number', IBATCH
            CALL WRTBLKN(VEC1,LBATCHB,LBLOCK(IOFF))
          END IF
*. Obtain corresponding C blocks
          CALL GET_BLOCKS_FROM_DISC
     &    (LU1,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC2,1)
          E = E + INPROD(VEC1,VEC2,LBATCHE)
          IF(IPRT.GE.200) THEN
            WRITE(6,*) ' Initial batch of C, number', IBATCH
            CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
          END IF
          IOFF = IOFF + LBATCHB
        END DO
        IF(IRESTRICT.EQ.1) E = 2*E
        EIG(1,IROOT) = E                                 
      ELSE IF(ISKIPEI.EQ.1) THEN
        E = FINEIG(IROOT) 
        WRITE(6,*) ' Initial energy obtained from previous calc as ',
     &             E+EIGSHF
        EAPR = E
        EIG(1,IROOT) = E                                 
      END IF
        
*
      IF(IPRT .GE. 3 ) THEN
        WRITE(6,'(A,I4)') ' Eigenvalues of initial iteration '
        WRITE(6,'(5F18.13)')
     &  ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
      END IF
      ITERX = 1
      CALL QEXIT('INIIT')
*
* ======================
*. Loop over iterations
* ======================
*
      DO 1000 ITER = 2, MAXIT
       IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==============================='      
        WRITE(6,*) ' Info from iteration .... ', ITER
        WRITE(6,*) ' ==============================='      
        WRITE(6,*)
       END IF
*. Allow loop over micro-iterations : 
*. In the first micro of a given iteration, the complete Sigma-vector
*. is calculated, in the following micro's, the Sigma vector is 
*. restricted to the space of the C-vectors
      DO 999, MICRO_IT = 1, MAX_MICRO_IT
*
       IF(IPRT  .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================================='      
        WRITE(6,*) ' Info from micro iteration .... ', MICRO_IT
        WRITE(6,*) ' ====================================='      
        WRITE(6,*)
       END IF
       CALL QENTER('PARTA')
*. Largest allowed basespace in this iteration in multispace method
       IF(IMULTGR.GT.0) THEN
         IBASSPC_MX = 1-IPAT(MOD(ITER-2,LPAT)+1)+IREFSPC
         IF(ITER.EQ.MAXIT) IBASSPC_MX = IREFSPC
         WRITE(6,*) ' Max allowed base space ', IBASSPC_MX
       ELSE
         IBASSPC_MX = 0         
       END IF
*
       ZERO = 0.0D0
       IF(ICLSSEL.EQ.1) THEN
         CALL SETVEC(CLS_CT,ZERO,NCLS)
         CALL SETVEC(CLS_ET,ZERO,NCLS)
         CALL SETVEC(CLS_DELT,ZERO,NCLS)
         CALL SETVEC(CLS_C ,ZERO,NCLS)
         CALL SETVEC(CLS_E ,ZERO,NCLS)
         CALL SETVEC(CLS_DEL,ZERO,NCLS)
*
         CALL SETVEC(EBASC,ZERO,NSPC)
         CALL SETVEC(CBASC,ZERO,NSPC)
*
       END IF
*. Active classes
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU1,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
*. Mark active blocks and find required number of batches
          CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
        ELSE
          NBATCHL = NBATCH
        END IF
       IF(IPRT.GE.10) WRITE(6,*) ' Number of batches for C ', NBATCHL
       EIGAPR = E
       CALL QEXIT('PARTA')
* ===============================
*. Obtain C(T) (H0-E)**-1 C
* ===============================
       CALL QENTER('PARTB')
       GAMMA = 0.0D0
       IOFF = 1
       CALL REWINO(LU1)
       DO IBATCH = 1, NBATCHL
         LBATCHB = LLBATCHB(IBATCH)
         LBATCHE = LLBATCHE(IBATCH)
*. Retrieve Batch of C
         NO_ZEROING = 0
         NO_ZEROING1= 1
         CALL FRMDSCN3(VEC1,LBATCHB,LBLK,LU1,NO_ZEROING1,I2BLK(IOFF),
     &                 LBLOCK(IOFF))
         CALL COPVEC(VEC1,VEC2,LBATCHE)
*. Multiply with (H0-E)** -1 
         FACTOR = -EIGAPR
         ITASK = 1
         CALL DIATERM_GAS(FACTOR,ITASK,VEC2,LBATCHB,IBLOCK,IOFF,0,
     &         NO_ZEROING1,I2BLK(IOFF))
         IF(NO_ZEROING1.EQ.0) THEN
           GAMMA = GAMMA + INPROD(VEC1,VEC2,LBATCHE)
         ELSE
            GAMMA = GAMMA + INPRODB(VEC1,VEC2,LBATCHB,LBLOCK(IOFF),
     &                      I2BLK(IOFF))
         END IF
         IOFF = IOFF + LBATCHB
       END DO
       IF(IPRT.GE.20)
     & WRITE(6,*) ' Gamma  calculated ',GAMMA
       CALL QEXIT('PARTB')

*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = (H0-E)-1 (H*C - EIGAPR*C) : Obtain in batches and save on DISC
*
C      EIGAPR = EIG(ITER-1,1)
       CALL QENTER('PARTC')
       XNELMNT = 0.0D0
       XNZERO = 0.0D0
*
       RNORM = 0.0D0
       DELTA = 0.0D0
       CHEDEL =0.0D0
       DELNORM = 0.0D0
*
       DELTAT = 0.0D0
       CHEDELT =0.0D0
       DELNORMT = 0.0D0
       ECC = 0.0D0
*
       CALL REWINO(LU2)
       IOFF = 1
*. Find partitioning of sigma
       IF(MICRO_IT.EQ.1) THEN
*. BLocks obttained by double excitations from classes in C
         NEXC  = 2
       ELSE  
*. Keep sigma-vector equal to C vector 
         NEXC = 0
       END IF
*
       CALL EXCCLS2(NCLS,ICLS_A,ICLS_A2,NEXC,IBASSPC_MX,IBASSPC)
*. Partitioning of sigma vector      
        CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A2,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
    
       IF(IPRT.GE.10) THEN
         WRITE(6,*) ' Number of batches for S ', NBATCHL
         WRITE(6,*)
       END IF
       CALL QEXIT('PARTC')
*
       DO IBATCH = 1, NBATCHL
         LBATCHB = LLBATCHB(IBATCH)
         LBATCHE = LLBATCHE(IBATCH)
         XNELMNT = XNELMNT + FLOAT(LBATCHE)
         IF(IPRT.GE.10) 
     &   WRITE(6,*) '  Delta, Batch : ',IBATCH
* Batch of HC in VEC1
         ZERO = 0.0D0
         CALL SETVEC(VEC1,ZERO,LBATCHE)
         CALL QENTER('PARTD')
         CALL SBLOCK(LBATCHB,IBLOCK,IOFF,VEC2,VEC1,LU1,0,0,
     &                0,0,0,0)
         CALL QEXIT('PARTD')
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' Batch of H C '
           CALL WRTBLKN(VEC1,LBATCHB,LBLOCK(IOFF)) 
         END IF
*. Retrieve Batch of C
         CALL QENTER('PARTE')
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC2,1)
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' C batch read in '
           CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
         END IF
*. Update energy
         ECC = ECC + INPROD(VEC1,VEC2,LBATCHE)
* Batch of (H-E)C in VEC1
         ONE = 1.0D0
         FACTOR = -EIGAPR
         CALL VECSUM(VEC1,VEC1,VEC2,ONE,FACTOR,LBATCHE)
         IF(IPRT.GE.500) THEN
           WRITE(6,*) ' Batch of (H - E ) C '
           CALL WRTBLKN(VEC1,LBATCHB,LBLOCK(IOFF))
         END IF
*. Norm of residual
         RNORM = RNORM + INPROD(VEC1,VEC1,LBATCHE)
*. Batch of (H0-E)-1(H-E)C  in VEC2
         CALL COPVEC(VEC1,VEC2,LBATCHE)
         FACTOR = -EIGAPR
         ITASK = 1
         I12 = 1
         CALL DIATERM_GAS(FACTOR,ITASK,VEC2,LBATCHB,IBLOCK,IOFF,0,0,0)
         DELNORMT = DELNORMT + INPROD(VEC2,VEC2,LBATCHE)
* C(H-E)(H0-E0)-1(H-E)C
         CHEDELT = CHEDELT + INPROD(VEC1,VEC2,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
*. Contributions divided into occupation classes, complete expansion
*. Wave function correction
           CALL CLASS_PROD3(VEC2,VEC2,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_CT)
*. Energy correction
           CALL CLASS_PROD3(VEC1,VEC2,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_ET)
         END IF
*.[(H0-E0)-1(H-E)C]_{truncated} and 
*.C(H-E) [(H0-E0)-1(H-E)C]_{truncated}
         ZERO = 0.0D0
         CALL SETVEC(VEC3,ZERO,LBATCHE)
         XNZERO = XNZERO + FLOAT(LBATCHE)
         DO I = 1, LBATCHE
           IF(ABS(VEC2(I)*VEC1(I)).GE.THRES_EC.OR.
     &        ABS(VEC2(I)).GE.THRES_CC           ) THEN
             CHEDEL = CHEDEL + VEC1(I)*VEC2(I)
             DELNORM = DELNORM + VEC2(I)*VEC2(I)
             VEC3(I)=VEC2(I)
             XNZERO = XNZERO - 1.0D0
           END IF
         END DO
         CALL QEXIT('PARTE')
*
         CALL QENTER('PARTF')
         IF(ICLSSEL.EQ.1) THEN
*. Contributions divided into occupation classes, truncated expansion
*. Wave function correction
           CALL CLASS_PROD3(VEC3,VEC3,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_C)
*. Energy correction
           CALL CLASS_PROD3(VEC1,VEC3,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_E)
         END IF
*. retrieve c batch from disc
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC1,1)
* C(H0-E0)-1(H-E)C
         DELTAT = DELTAT + INPROD(VEC1,VEC2,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
           CALL CLASS_PROD3(VEC1,VEC2,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_DELT)
         END IF
* C[(H0-E0)-1(H-E)C]{truncated}
         DELTA = DELTA + INPROD(VEC1,VEC3,LBATCHE)
         IF(ICLSSEL.EQ.1) THEN
           CALL CLASS_PROD3(VEC1,VEC3,IOFF,LBATCHB,IBLOCK,
     &                      IBLK_TO_CLS,NCLS,CLS_DEL)
         END IF
*
*. Write packed version to DISC
*. Pack out so zero blocks are given zero entries
         CALL TODSCNP(VEC3,LBATCHB,LBLOCK(IOFF),LBLK,LU2)
         IF(IPRT.GE.200) THEN
           WRITE(6,*) ' Batch of blocks of trial vector '
           CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
         END IF
         IOFF = IOFF + LBATCHB
         CALL QEXIT('PARTF')
       END DO
*      /\ End of loop over batches of correction vector
         CALL QENTER('PARTG')
         WRITE(6,'(A,E22.15)') 
     &   ' Number of zero elements in delta (before class trunc) ',
     &    XNZERO
         WRITE(6,'(A,E22.15)') 
     &   ' Number of nonzero terms in delta (before class trunc) ',
     &    XNELMNT-XNZERO
       WRITE(6,'(A,E25.16)') ' ECC+EIGSHF = ',ECC+EIGSHF
       CALL ITODS(-1,1,LBLK,LU2)
       IF(MICRO_IT.EQ.1 ) RNRM(ITER-1,1) = SQRT(RNORM)
*. (End of loop over batches of (H0-E)-1(H-E)C)
*. Predicted energy  
       IF(ICLSSEL.EQ.1) THEN
*. Energy and wave function per base space
         CALL QENTER('G1   ')
         CALL CLS_TO_BASE(CLS_E,EBASC,CLS_C,CBASC,NCLS,NSPC,
     &                    IBASSPC) 
         CALL QEXIT('G1   ')
*. decide which classes should be truncated
         CALL QENTER('PARTT')
         CALL CLASS_TRUNC(NCLS,ICLS_L,RCLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
     &                    E_CONV,ICLS_A,N_CLS_TRN,E_CLS_TRN,W_TRN)
         CALL QEXIT('PARTT')
*. Update corrections for class elimination
         N_ACT_CLS = 0
         CHEDEL = CHEDEL - (-1.0D0) * E_CLS_TRN
         DO JCLS = 1, NCLS
           IF(ICLS_A(JCLS).EQ.0) THEN
             DELTA = DELTA-CLS_DEL(JCLS) 
           ELSE
             N_ACT_CLS = N_ACT_CLS + 1
           END IF
         END DO
           
*. And do the truncation 
*. Truncation of classes => truncation of blocks
         IF(N_CLS_TRN.NE.0) THEN
           CALL QENTER('G2   ')
           CALL CLS_TO_BLK(NBLOCKT,IBLK_TO_CLS,ICLS_A,IBLKS_A)
           CALL QEXIT('G2   ')
*. from LU2 to LU3 and back to LU2
           CALL QENTER('G3   ')
           CALL ZAP_BLOCK_VEC(LU2,LBLK,IBLKS_A,VEC2,LU3)
           CALL QEXIT('G3   ')
         END IF
       ELSE
         N_CLS_TRN = 0
         E_CLS_TRN = 0.0D0
         W_CLS_TRN = 0.0D0
       END IF
*. Predicted energy  
       IF(GAMMA.NE.0.0D0) THEN
         E2PREDIT = - CHEDELT + DELTAT**2/GAMMA
         E2PREDI  = - CHEDEL  + DELTA * DELTAT /GAMMA
       ELSE
         E2PREDIT = - CHEDELT
         E2PREDI  = - CHEDEL 
         IF(ICLSSEL.EQ.1) THEN
           CALL ICOPVE(CLS_CT,CLS_C,NCLS)
           CALL ICOPVE(CLS_ET,CLS_E,NCLS)
         END IF
       END IF
*
       WRITE(6,*)
*. 
C?     WRITE(6,*) ' Information for untruncated expansion:'
C?     WRITE(6,*) ' ======================================'
C?     WRITE(6,*) 
C?   & ' CHEDELT DELTAT GAMMA ', CHEDELT,DELTAT,GAMMA
       IF(GAMMA.NE.0.0D0) THEN
         WRITE(6,*) 
     & ' Orthogonalization term to E2 (no trunc.)', DELTAT**2/GAMMA
       END IF
       WRITE(6,'(A,2E25.15)') 
     & ' Predicted energy(no truncation), change and total ', 
     &              E2PREDIT,EIGAPR+E2PREDIT+EIGSHF
       IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
C?       WRITE(6,*) 
C?       WRITE(6,*) ' Information for truncated expansion:'
C?       WRITE(6,*) ' ======================================'
C?       WRITE(6,*) 
C?   &   ' CHEDEL DELTA GAMMA ', CHEDEL,DELTA,GAMMA
         WRITE(6,*)
         IF(GAMMA.NE.0.0D0) THEN
           WRITE(6,*) 
     &   ' Orthogonalization term to E2 (trunc.)', DELTA**2/GAMMA
         END IF
         WRITE(6,'(A,2E25.15)') 
     &   ' Predicted energy (truncated), change and total ', 
     &                E2PREDI,EIGAPR+E2PREDI+EIGSHF
         WRITE(6,*)
         WRITE(6,*) ' Estimated square-norm of eliminated terms ',
     &   DELNORMT-DELNORM
         WRITE(6,'(A,E25.15)') 
     &   ' Estimated energy contributions of eliminated terms',
     &   E2PREDIT-E2PREDI
C        WRITE(6,'(A,E22.15)') 
C    &   ' Number of zero elements in delta (before class sel) ',
C    &    XNZERO
C        WRITE(6,'(A,E22.15)') 
C    &   ' Number of nonzero terms in delta (before class sel) ',
C    &    XNELMNT-XNZERO
       ELSE
       END IF
       CALL QEXIT('PARTG')
*
       IF(ICLSSEL.EQ.1) THEN
         IF(N_ACT_CLS .EQ. 0   ) THEN
           IF(IMULTGR.EQ.0.OR.IBASSPC_MX.EQ.IREFSPC) THEN
*. All classes were eliminated so we are home -and hopefully dry
           WRITE(6,*) ' No active classes  '    
           IF(MICRO_IT.EQ.1) THEN
             WRITE(6,*) ' I will therefore end the diagonalization'
             CONVER = .TRUE.
             GOTO 1001
           ELSE 
             WRITE(6,*) ' I will continue to the next macroiteration'
             GOTO 1000
           END IF
           ELSE
*. No active classes with this IBASSPC_MX, try next 
             EIG(ITER,1) =  EIG(ITER-1,1)
             ITERX = ITER
             WRITE(6,*) ' No active classes  '    
             WRITE(6,*) ' I go to the next iteration '
             GOTO 1000
           END IF
         END IF
       END IF
*
* ============================================
* 1.5 : Inverse Iteration Correction to Delta 
* ============================================
*
*
* Update delta to 
* -(H0-E0)-1(H-E)|0> + delta/gamma * (H0-E0)-1 |0>
*
* ( was +(H0-E0)-1(H-E)|0>)
* 
       CALL REWINO(LU3)
       CALL REWINO(LU1)
       CALL REWINO(LU2)
       IOFF = 1
       IF(IOLSTM.EQ.1.AND.ABS(GAMMA).GT.1.0D-6) THEN
         WRITE(6,*) ' Inverse iteration correction will be added '
         DO IBATCH = 1, NBATCHL
           LBATCHB = LLBATCHB(IBATCH)
           LBATCHE = LLBATCHE(IBATCH)
*. Retrieve Batch of C
           NO_ZEROING = 0
           CALL FRMDSCN3(VEC1,LBATCHB,LBLK,LU1,NO_ZEROING,I2BLK(IOFF),
     &                   LBLOCK(IOFF))
*. Multiply with (H0-E)** -1 
           FACTOR = -EIGAPR
           ITASK = 1
           CALL DIATERM_GAS(FACTOR,ITASK,VEC1,LBATCHB,IBLOCK,IOFF,0,
     &                       0,0)
*. Retrieve Batch of Delta
           NO_ZEROING = 0
           CALL FRMDSCN3(VEC2,LBATCHB,LBLK,LU2,NO_ZEROING,I2BLK(IOFF),
     &                   LBLOCK(IOFF))
*. And add
           FAC1 = -1.0D0
           FAC2 = DELTA/GAMMA
           CALL VECSUM(VEC2,VEC2,VEC1,FAC1,FAC2,LBATCHE)
*. Transfer to Disc
           CALL TODSCNP(VEC2,LBATCHB,LBLOCK(IOFF),LBLK,LU3)
           IOFF = IOFF + LBATCHB
         END DO
*        ^ End of loop over batches
         CALL ITODS(-1,1,LBLK,LU3)
*. Well, it is nice to have the correction vector on LU2 so 
         IREW = 1
         CALL COPVCD(LU3,LU2,VEC1,IREW,LBLK)
       END IF
*
* ===================================
* 2 : Calculate <Delta ! H ! Delta >
* ===================================
*
*. Active classes
        CALL QENTER('PARTH')
        IF(IDYNBATCH.EQ.1) THEN
          CALL FIND_ACTIVE_CLASSES(LU2,LBLK,IBLK_TO_CLS,
     &         ICLS_A,NCLS,VEC1)
          CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &         ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
        ELSE
          NBATCHL = NBATCH
        END IF
* Loop over batches of H !delta>
       IOFF = 1
       DELHDEL = 0.0D0
       IRESTRICT = 1
       WRITE(6,*) ' Number of batches for Delta = ', NBATCHL
       DO IBATCH = 1, NBATCHL
         IF(IPRT.GE.10) 
     &   WRITE(6,*) '  <Delta ! H ! Delta >, Batch : ',IBATCH
         LBATCHB = LLBATCHB(IBATCH)
         LBATCHE = LLBATCHE(IBATCH)
         CALL SETVEC(VEC1,ZERO,LBATCHE)
         CALL SBLOCK(LBATCHB,IBLOCK,IOFF,VEC2,VEC1,LU2,IRESTRICT,0,
     &                0,0,0,0)
*. Retrieve Batch of Delta
         CALL GET_BLOCKS_FROM_DISC
     &   (LU2,LBATCHB,IOFF,IBLOCK,NBLOCKT,VEC2,1)
         IF(IPRT.GE.1000) THEN
           WRITE(6,*) ' Blocks of delta retrieved '
           CALL WRTBLKN(VEC2,LBATCHB,LBLOCK(IOFF))
         END IF
         DELHDEL = DELHDEL + INPROD(VEC1,VEC2,LBATCHE)
         IOFF = IOFF + LBATCHB
       END DO
       CALL QEXIT('PARTH')
       IF(IRESTRICT.EQ.1) DELHDEL = 2.0D0*DELHDEL
*
* ===========================================
* 3 : Solve 2 by 2 problem : Nonorthogonal !!
* ===========================================) 
*
*      Norm of delta and overlap between 0 and delta
       S12 = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
       S22 = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
C      H11 = EIGAPR
       H11 = ECC
       IF(IOLSTM.EQ.1.AND.ABS(GAMMA).GT.1.0D-6) THEN
        H12 = -CHEDEL + DELTA*DELTAT/GAMMA
       ELSE 
        H12 = CHEDEL + EIGAPR*DELTA
       END IF
       H22 = DELHDEL
*
       S11 = 1.0D0
*
*.( H11  H12 ) (X1)       (S11    S12 )(X1)
* (          ) (  )   = E (           )( )
* ( H12  H22 ) (X2)       (S12    S22 )(X2)
*
* The eigenvalues
*
        A = S11*S22 -S12 **2
        B = 2*S12*H12-S11*H22-S22*H11
        C = H11*H22 -H12**2
*
        EA = -B/(2*A) - SQRT(B**2 - 4*A*C)/(2*A)
        EB = -B/(2*A) + SQRT(B**2 - 4*A*C)/(2*A)
*. And the lowest eigenvalue is 
        E = MIN(EA,EB)
*. The corresponding eigenvector
*. Intermediate normalization
        X1 = 1.0D0
COLD    X2 = -(H11-E*S11)/(H12-E*S12)
*. A stable form when H11 approx E
        X2 = -(H12-E*S12)/(H22-E*S22)
*. Normalized
        XNORM2 = S11*X1**2 + S22*X2**2 + 2.0*S12*X1*X2
        XNORM = SQRT(XNORM2)
        X1 = X1/XNORM
        X2 = X2/XNORM
*
        IF(IPRT.GE.10) THEN
          WRITE(6,*) 
          WRITE(6,*) ' 2 X 2 Generalized eigenvalue problem, H and S '
          WRITE(6,*) 
          WRITE(6,'(4X,E20.10)') H11
          WRITE(6,'(4X,2E20.10)') H12,H22
          WRITE(6,*) 
          WRITE(6,'(4X,E20.10)') S11
          WRITE(6,'(4X,2E20.10)') S12,S22
          WRITE(6,*)
*
          WRITE(6,'(A,E25.15)') 
     &    ' Lowest eigenvalue (with shift)', E + EIGSHF
          WRITE(6,*) ' Corresponding eigenvector ', X1,X2
        END IF

*. Save corresponding eigenvector on file LU1 ( first LU3, then COPY)
* VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
        IREW = 1
C       IF(THRES_EC.NE.0.0D0.OR.THRES_CC.NE.0.0D0) THEN
          CALL VECSMDP(VEC1,VEC2,X1,X2,LU1,LU2,LU3,IREW,LBLK)
          CALL COPVCD(LU3,LU1,VEC1,IREW,LBLK)
C       ELSE
C         CALL VECSMD(VEC1,VEC2,X1,X2,LU1,LU2,LU3,IREW,LBLK)
C         CALL COPVCD(LU3,LU1,VEC1,IREW,LBLK)
C       END IF
        IF(MICRO_IT.EQ.1 ) EIG(ITER,1) = E 
        ITERX = ITER
*. Convergence of complete iteration sequence ?
C    &                 NSPC,IMULTGR,IPAT,LPAT,IREFSPC)
      IF((IMULTGR.EQ.0.OR.IBASSPC_MX.EQ.IREFSPC).AND.MICRO_IT.EQ.1) THEN
        IF(ABS(EIG(ITER,1) - EIG(ITER-1,1)).LE.THRES_ET)
     &     CONVER = .TRUE.
      END IF
      IF(CONVER) GOTO 1001
  999 CONTINUE
*     ^ End of loop over micro-iterations 
 1000 CONTINUE
* ( End of loop over iterations )
 1001 CONTINUE
      ITER = ITERX
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*
      IF ( IPRT .GT. 1 ) THEN
        CALL REWINO(LU1)
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(6,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER-1
          IF(IMULTGR.EQ.0) THEN
            WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          ELSE 
            IDEL = 1-IPAT(MOD(I-1,LPAT)+1)
            IF(I.EQ.1.OR.IDEL.EQ.0.OR.I.EQ.MAXIT-1) THEN
              WRITE(6,1341) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1341         FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5,'  Full resid.')
            ELSE
              WRITE(6,1342) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1342         FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5,'  Partial resid.')
            END IF
          END IF
 1330     CONTINUE
          WRITE(6,1340) ITER,EIG(ITER,IROOT)+EIGSHF
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        DO 1601 IROOT = 1, NROOT
COLD       FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           FINEIG(IROOT) = EIG(ITER,IROOT)
 1601   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,E13.6,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER-1,IROOT)
 1607   CONTINUE
      END IF
C
*. Clean up : IBLOCK in original form
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      CALL REPART_CIV(IBLOCK,NBATCHL,LLBATCHB,LLBATCHE,MXLNG,
     &     ICLS_A,IBLK_TO_CLS,NCLS,NBLOCKT,LBLOCK)
*
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
      SUBROUTINE RSBB2BN2(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                   JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                   IAGRP,IBGRP,IOCTPA,IOCTPB, 
     &                   NGAS,IAOC,IBOC,JAOC,JBOC,
     &                   SB,CB,ADSXA,STSTSX,MXPNGASX,
     &                   NOBPTS,MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                   XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBSX,IUSEAB,
     &                   CJRES,SIRES,SCLFAC,NTESTG,
     &                   NSEL2E,ISEL2E,IUSE_PH,IPHGAS,XINT2,NSTFSMSPGP,
     &                   IASPGP,IBSPGP,JASPGP,JBSPGP)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
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
*
* Last change : Aug 2000
*. Some dull optimazation, July 2003
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INCLUDE 'mxpdim.inc'
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
      REAL*8 INPROD
      INTEGER NSTFSMSPGP(MXPNSMST,*)
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*), XINT2(*)
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
      INTEGER NKA_AS(MXPNSMST), NKB_AS(MXPNSMST)
*. Arrays for reorganization 
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C    &          SIGNREO,NADOP,NADDEL,IADDEL,ADSIGN)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'oper.inc'
      INCLUDE 'multd2h.inc'
*
      COMMON/XXTEST/ISETVECOPS(10)
*
      CALL QENTER('RS2B ')
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.500) THEN
*
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' RSBB2BN2 speaking '
        WRITE(6,*) ' ================ '
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
        WRITE(6,*) ' Memcheck at start of RSBB2BN '
        CALL MEMCHK
        WRITE(6,*) ' Memory check passed '
*
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Groups defining each supergroup
C     CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
C     CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
C     CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
C     CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)           
*
      ITEST = 0
      IF(ITEST.EQ.1.AND.(NIJTYP.EQ.0.AND.NKLTYP.EQ.0)) THEN
        WRITE(6,*) ' Unneccesary entrance to RSBB2BN2'
        WRITE(6,*) ' IAOC, IBOC,JAOC,JBOC : '
        WRITE(6,*) ' Occupation of IA '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of IB '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JA '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of JB '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' IATP, JATP, IBTP, JBTP = ',
     &               IATP, JATP, IBTP, JBTP
        STOP 
      END IF
*
      IAFRST = 1
      IBFRST = 1
      JAFRST = 1
      JBFRST = 1
*
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
*
        IJFIRST = 1
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Should N-1 or N+1 projection be used for alpha strings
          IJ_TYP(1) = ITYP
          IJ_TYP(2) = JTYP
          IJ_AC(1)  = 2
          IJ_AC(2) =  1
          NOP = 2
          IF(IUSE_PH.EQ.1) THEN
            CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,
     &           SIGNIJ)
          ELSE
*. Enforced a+ a
            IJ_REO(1) = 1
            IJ_REO(2) = 2
            SIGNIJ = 1.0D0
          END IF
*. Two choices here :
*  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> ( good old creation mapping)
*  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka>  + delta(i,j)                   
C?        WRITE(6,*) ' RSBB2BN : IOP_REO : ', (IOP_REO(II),II=1,2)
          IF(IJ_REO(1).EQ.1.AND.IJ_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
            IJAC = 2
            IOP2AC = 1
            SIGNIJ2 = SCLFAC
*
            IJ_DIM(1) = NI
            IJ_DIM(2) = NJ
            IJ_SYM(1) = ISM
            IJ_SYM(2) = JSM
            IJ_TYP(1) = ITYP
            IJ_TYP(2) = JTYP
*
            NOP1   = NI
            IOP1SM = ISM
            IOP1TP = ITYP
            NOP2   = NJ
            IOP2SM = JSM
            IOP2TP = JTYP
          ELSE
*. Terra Nova, annihilation map 
            IJAC = 1
            IOP2AC = 2
            SIGNIJ2 = -SCLFAC
*
            IJ_DIM(1) = NJ
            IJ_DIM(2) = NI
            IJ_SYM(1) = JSM
            IJ_SYM(2) = ISM
            IJ_TYP(1) = JTYP
            IJ_TYP(2) = ITYP
*
            NOP1   = NJ
            IOP1SM = JSM
            IOP1TP = JTYP
            NOP2   = NI
            IOP2SM = ISM
            IOP2TP = ITYP
          END IF
*
          IF(IJFIRST.EQ.1) THEN
*. Find supergroup type of Kstring
C  NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
            CALL NEWTYP(JATP+IOCTPA-1,IOP2AC,IOP2TP,1,KATP_ABS)
            CALL ICOPVE(NSTFSMSPGP(1,KATP_ABS),NKA_AS,NSMST)
            IJFIRST = 0
          END IF
*
*. Generate creation- or annihilation- mappings for all Ka strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          KASM = MULTD2H(JASM,IJ_SYM(2) )
          NKASTR = NKA_AS(KASM)
          CALL ADAST2_GAS(IJ_SYM(2),IJ_TYP(2),NGAS,JASPGP,JASM,
     &         I1,XI1S,NKASTR,IEND,JAFRST,KFRST,KACT,SIGNIJ2,IJAC,
     &         JAACT,1)
          IF(JAACT.EQ.1) JAFRST = 0
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          CALL ADAST2_GAS(IJ_SYM(1),IJ_TYP(1),NGAS,IASPGP,IASM,
     &         I3,XI3S,NKASTR,IEND,IAFRST,KFRST,KACT,ONE,IJAC,
     &         IAACT,2)
          IF(IAACT.EQ.1) IAFRST = 0
*. Compress list to common nonvanishing elements
          IDOCOMP = 0
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),
     &                        NKASTR,NKAEFF)
          ELSE 
              NKAEFF = NKASTR
          END IF
            
*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK   
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
*
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
CJTEST      WRITE(6,*) ' JTEST: Dimension of CJRES and SIRES ',
C    &      IJ_DIM(2)*LKABTC*NJB, IJ_DIM(1)*LKABTC*NIB
*. Obtain C(ka,J,JB) for Ka in batch
            I_SHOULD_CONSTRUCT_SICJ = 1
CM          IF(I_SHOULD_CONSTRUCT_SICJ.EQ.1) THEN
CM           DO JJ = 1, IJ_DIM(2)
CM             CALL GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,
CM   &              JJ,I1(KABOT+(JJ-1)*NKASTR),
CM   &              XI1S(KABOT+(JJ-1)*NKASTR))
CM           END DO
CM           IF(NTEST.GE.500) THEN
CM             WRITE(6,*) ' Updated CJRES as C(Kaj,Jb)'
CM             CALL WRTMAT(CJRES,NKASTR*NJ,NJB,NKASTR*NJ,NJB)
CM           END IF
*
CM           ISETVECOPS(2) = ISETVECOPS(2) + NIB*LKABTC*IJ_DIM(1)
CM           MXACJ=MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
CM           CALL SETVEC(SIRES,ZERO,NIB*LKABTC*IJ_DIM(1))
*
CM           I_SHOULD_CONSTRUCT_SICJ = 0
CM          END IF
*
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
              KLFIRST = 1
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
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = ADSXA(KSM,KLSM)
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' KSM, LSM', KSM, LSM
                END IF
                IF(LSM.EQ.0) GOTO 1930
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*
                IF(KL_REO(1).EQ.1.AND.KL_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
                  IOP4AC = 1
                  KLAC = 2
                  KL_DIM(1) = NK
                  KL_DIM(2) = NL
                  KL_SYM(1) = KSM
                  KL_SYM(2) = LSM
                  KL_TYP(1) = KTYP
                  KL_TYP(2) = LTYP
                ELSE
*. Terra Nova, annihilation map 
                  IOP4AC = 2
                  KLAC = 1
                  KL_DIM(1) = NL
                  KL_DIM(2) = NK
                  KL_SYM(1) = LSM
                  KL_SYM(2) = KSM
                  KL_TYP(1) = LTYP
                  KL_TYP(2) = KTYP
                END IF
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0  
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
                IF(KLFIRST.EQ.1) THEN
*. Find supergroup type of Kstring
C  NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
                  CALL NEWTYP(JBTP+IOCTPB-1,IOP4AC,KL_TYP(2),1,KBTP_ABS)
                  CALL ICOPVE(NSTFSMSPGP(1,KBTP_ABS),NKB_AS,NSMST)
                  KLFIRST = 0
                END IF
                KBSM = MULTD2H(JBSM,KL_SYM(2) )
                NKBSTR = NKB_AS(KBSM)
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
*. currently we are using creation mappings for kl
*. (Modify to use ADAST later )
                CALL ADAST2_GAS(KL_SYM(2),KL_TYP(2),NGAS,JBSPGP,JBSM,
     &               I2,XI2S,NKBSTR,IEND,JBFRST,KFRST,KACT,SIGNKL,KLAC,
     /               JBACT,3)
                IF(JBACT.EQ.1) JBFRST = 0
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADAST2_GAS(KL_SYM(1),KL_TYP(1),NGAS,IBSPGP,IBSM,
     &               I4,XI4S,NKBSTR,IEND,IBFRST,KFRST,KACT,ONE,KLAC,
     &               IBACT,4)
                IF(IBACT.EQ.1) IBFRST = 0
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (iop2 iop1 |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                ONE = 1.0D0
                IF(I_USE_SIMTRH .EQ.0 ) THEN
*. Normal integrals with conjugation symmetry
                  CALL GETINT(XINT,IJ_TYP(2),IJ_SYM(2),
     &                 IJ_TYP(1),IJ_SYM(1),
     &                 KL_TYP(1),KL_SYM(1),KL_TYP(2),KL_SYM(2),IXCHNG,
     &                 0,0,ICOUL,ONE,ONE)
                ELSE IF (I_USE_SIMTRH.EQ.1) THEN
C?              WRITE(6,*) ' I_USE_SIMTRH = ', I_USE_SIMTRH
*. Integrals does not have conjugation symmetry so be careful...
*. The following is not enough is particle hole symmetry is encountered
*. Obtain ( i j ! k l )
                  CALL GETINT(XINT,ITYP,ISM,JTYP,JSM,
     &                             KTYP,KSM,LTYP,LSM,
     &                        IXCHNG,0,0,ICOUL,ONE,ONE)
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
C               IJKL_DIM = IJ_DIM(1)*IJ_DIM(2)*KL_DIM(1)*KL_DIM(2)
C               IF(INPROD(XINT,XINT,IJKL_DIM).NE.0.0D0) THEN
            IF(I_SHOULD_CONSTRUCT_SICJ.EQ.1) THEN
             DO JJ = 1, IJ_DIM(2)
               CALL GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,
     &              JJ,I1(KABOT+(JJ-1)*NKASTR),
     &              XI1S(KABOT+(JJ-1)*NKASTR))
C              GET_CKAJJB(CB,NJ,NJA,CKAJJB,NKA,NJB,J,ISCA,SSCA)


             END DO
             IF(NTEST.GE.500) THEN
               WRITE(6,*) ' Updated CJRES as C(Kaj,Jb)'
               CALL WRTMAT(CJRES,NKASTR*NJ,NJB,NKASTR*NJ,NJB)
             END IF
*
             ISETVECOPS(2) = ISETVECOPS(2) + NIB*LKABTC*IJ_DIM(1)
             MXACJ=MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
             CALL SETVEC(SIRES,ZERO,NIB*LKABTC*IJ_DIM(1))
*
             I_SHOULD_CONSTRUCT_SICJ = 0
            END IF
                IROUTE = 3
                CALL SKICKJ(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),
     &               KL_DIM(1),KL_DIM(2),
     &               NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE )
C               END IF
*
                IF(NTEST.GE.500) THEN
                  WRITE(6,*) ' Updated Sires as S(Kai,Ib)'
                  CALL WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
                END IF
*
 1930         CONTINUE
*             ^ End of loop over KSM
 2000       CONTINUE
*           ^ End of loop over KLTYP
*
*. Scatter out from s(Ka,Ib,i)
*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' S(Ka,Ib,i) as S(Ka,Ibi)'
              CALL WRTMAT(SIRES,LKABTC,NIB*IJ_DIM(1),LKABTC,IJ_DIM(1))
            END IF
*. Was anything done ?
            IF(I_SHOULD_CONSTRUCT_SICJ.EQ.0) THEN
            DO II = 1, IJ_DIM(1)
              CALL ADD_SKAIIB(SB,IJ_DIM(1),NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
            END IF
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 1940   CONTINUE
*       ^ End of loop over ISM
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL QEXIT('RS2B ')
      RETURN
      END
      SUBROUTINE ADAST2_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,
     &                    I1,XI1S,NKSTR,IEND,IFRST,KFRST,KACT,SCLFAC,
     &                    IAC,IACT,IINDEX)
*
*
* Obtain creation or annihilation mapping
*
* IAC = 2 : Creation map
* a+IORB !KSTR> = +/-!ISTR> 
*
* IAC = 1 : Annihilation map
* a IORB !KSTR> = +/-!ISTR> 
*
* for orbitals of symmetry IOBSM and type IOBTP
* and Istrings defined by the NIGRP groups IGRP and symmetry ISPGPSM
* 
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR> 
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
* If info for the Istrings has been set up, IACT = 1
* IINDEX tells in which IOFFI info on Istrings should go
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version 
*               October 96   : Improved version
*               September 97 : annihilation mappings added
*                              I groups defined by IGRP
*               July 2003,   : Some dull optimization ...
*                              Notice that NKSTR is now input ...
*
*
* ======
*. Input
* ======
*
*./BIGGY
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'strinp.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'lucinp.inc'
*. Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER ISMFGS(MXPNGAS)
C     INTEGER MXVLI(MXPNGAS),MNVLI(MXPNGAS)
      INTEGER MXVLK(MXPNGAS),MNVLK(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      INTEGER KGRP(MXPNGAS)
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
      PARAMETER(LOFFI=8*8*8*8*8)
      COMMON/LOC_ADAST2/IOFFI(LOFFI,4),MXVLI(MXPNGAS,4),MNVLI(MXPNGAS,4)
*
      INCLUDE 'comjep.inc'
      INCLUDE 'multd2h.inc'
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension 
* (NKSTR,*), Where NKSTR is the number of K-strings of 
*  correct symmetry . Nk is provided by this routine.
*
C!    CALL QENTER('ADAST ')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADAST_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        WRITE(6,*) '  Symmetry of supergroup : ', ISPGPSM
        WRITE(6,*) ' SCLFAC = ', SCLFAC
*
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' Annihilation mapping '
        ELSE IF(IAC.EQ.2) THEN
          WRITE(6,*) ' Creation mapping '
        ELSE 
          WRITE(6,*) ' Unknown IAC parameter in ADAST ',IAC
          STOP       ' Unknown IAC parameter in ADAST '
        END IF
*
      END IF
*. A few preparations
      NORBTS= NOBPTS(IOBTP,IOBSM)
      NORBT= NOBPT(IOBTP)
      IACGAS = IOBTP
*
      IACT = 0
*. First orbital of given GASpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1 + NINOB
*. First orbital of given GASpace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM) 
*
*====================================================
*. K strings : Supergroup, symmetry and distributions
*====================================================
*
      IF(IAC.EQ.1) THEN
       IDELTA = +1
      ELSE
       IDELTA = -1
      END IF
*. Is required mapping contained within current set of maps?
*. a:) Is active GASpace included in IGRP - must be 
      IACGRP = 0
      DO JGRP = 1, NIGRP
       IF(IGSFGP(IGRP(JGRP)).EQ. IACGAS) IACGRP = JGRP
      END DO
*. Note : IACGRP is not the actual active group, it is the address of the
*         active group in IGRP
      IF(IACGRP.EQ.0) THEN
        WRITE(6,*) ' ADAST in problems '
        WRITE(6,*) ' Active GASpace not included in IGRP '
        WRITE(6,*) ' Active GASpace : ', IACGAS
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        STOP       ' ADAST : Active GASpace not included in IGRP '
      END IF
*. b:) active group in K strings
      NIEL = NELFGP(IGRP(IACGRP))
      NKEL = NIEL + IDELTA
      IF(NTEST.GE.1000) WRITE(6,*) ' NIEL and NKEL ',NIEL,NKEL
      IF(NKEL.EQ.-1.OR.NKEL.EQ.NOBPT(IACGAS)+1) THEN
*. No strings with this number of elecs - be happy : No work 
        NKSTR = 0
        KACT = 0
        KACGRP = 0
        GOTO 9999
      ELSE
*. Find group with NKEL electrons in IACGAS
        KACGRP = 0
        DO JGRP = IBGPSTR(IACGAS),IBGPSTR(IACGAS)+NGPSTR(IACGAS)-1
          IF(NELFGP(JGRP).EQ.NKEL) KACGRP = JGRP
        END DO
        IF(NTEST.GE.1000) WRITE(6,*) ' KACGRP = ',KACGRP
*. KACGRP is the Active group itself     
        IF(KACGRP.EQ.0) THEN
          WRITE(6,*)' ADAST : cul de sac, active K group not found'
          WRITE(6,*)' GAS space and number of electrons ',
     &               IACGAS,NKEL
          STOP      ' ADAST : cul de sac, active K group not found'
        END IF
      END IF
*. Okay active K group was found and is nontrivial
      KSM = MULTD2H(IOBSM,ISPGPSM)
*. The K supergroup
      CALL ICOPVE(IGRP,KGRP,NIGRP)
      KGRP(IACGRP) = KACGRP
*. Number of strings and symmetry distributions of K strings
      IF(NTEST.GE.1000) WRITE(6,*) 
     & ' KSM, NKSTR : ', KSM, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*. Last active space in K strings and number of strings per group and sym
      NGASL = 1
      DO JGRP = 1, NIGRP
       IF(NELFGP(KGRP(JGRP)).GT.0) NGASL = JGRP
      END DO
*. MIN/MAX for Kstrings
      DO JGRP = 1, NIGRP
        MNVLK(JGRP) =  MINMAX_SM_GP(1,KGRP(JGRP))
        MXVLK(JGRP) =  MINMAX_SM_GP(2,KGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) THEN
        write(6,*) 'MNVLK and MXVLK '
        CALL IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
        CALL IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
      END IF
*. (NKDIST_TOT is number of distributions, all symmetries )
* ==============
*. I Strings 
* ==============
*. Generate symmetry distributions of I strings with given symmetry
      NIGASL = 1
      DO JGRP = 1, NIGRP
        IF(NELFGP(IGRP(JGRP)).GT.0) NIGASL = JGRP
      END DO
*
      IF(IFRST.EQ.1) THEN
        DO IGAS = 1, NIGASL
          MNVLI(IGAS,IINDEX) = MINMAX_SM_GP(1,IGRP(IGAS))
          MXVLI(IGAS,IINDEX) = MINMAX_SM_GP(2,IGRP(IGAS))
        END DO
        CALL TS_SYM_PNT2(IGRP,NIGASL,
     &       MXVLI(1,IINDEX),MNVLI(1,IINDEX),ISPGPSM,
     &       IOFFI(1,IINDEX),LOFFI)
      END IF
      IACT = 1
*. Number of I strings per group and sym
*. Last entry in IGRP with a nonvanishing number of strings
*. Number of electrons before active space
      NELB = 0
      DO JGRP = 1, IACGRP-1
        NELB = NELB + NELFGP(IGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NELB = ', NELB
*
      ZERO =0.0D0
      IZERO = 0    
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
      DO IGAS = 1, NIGRP
        ISMFGS(IGAS) = 1
      END DO
 1000 CONTINUE
*. Next distribution
        CALL NEXT_SYM_DISTR(NGASL,MNVLK,MXVLK,ISMFGS,KSM,KFIRST,NONEW)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' Symmetry distribution ' 
          call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
        END IF
        IF(NONEW.EQ.1) GOTO 9999
        KFIRST = 0
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
C         NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
          NSTRIK = NSTRIK*NSTFSMGP(ISMFGS(IGAS),KGRP(IGAS))
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IACGRP)
        IACSM = MULTD2H(IOBSM,ISMFGS(IACGRP))
        ISMFGS(IACGRP) = IACSM
        IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,IOFFI(1,IINDEX),
     &             MXVLI(1,IINDEX),MNVLI(1,IINDEX))
        ISMFGS(IACGRP) = ISAVE
*. Number of strings before active GAS space
        NSTB = 1
        DO IGAS = 1, IACGRP-1
          NSTB = NSTB*NSTFSMGP(ISMFGS(IGAS),KGRP(IGAS))
        END DO
*. Number of strings After active GAS space
        NSTA = 1
        DO IGAS =  IACGRP+1, NIGRP
          NSTA = NSTA*NSTFSMGP(ISMFGS(IGAS),KGRP(IGAS))
        END DO
*. Number and offset for active group 
        NIAC  = NSTFSMGP(IACSM,IGRP(IACGRP))
        IIAC  = ISTFSMGP(IACSM,IGRP(IACGRP))
C       NIAC  = NACIST(IACSM)
C       IIAC =  IACIST(IACSM)
*
C       NKAC = NNSTSGP(ISMFGS(IACGRP),IACGRP)
C       IKAC = IISTSGP(ISMFGS(IACGRP),IACGRP)
        NKAC = NSTFSMGP(ISMFGS(IACGRP),KGRP(IACGRP))
        IKAC = ISTFSMGP(ISMFGS(IACGRP),KGRP(IACGRP))
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
        IF(NTEST.GE.1000) THEN
        write(6,*) ' nstb nsta niac nkac ',
     &               nstb,nsta,niac,nkac
        END IF
*. Obtain annihilation/creation mapping for all strings of this type
*. Are group mappings in expanded or compact form 
        IF(IAC.EQ.1.AND.ISTAC(KACGRP,2).EQ.0) THEN
          IEC = 2
          LROW_IN = NKEL
        ELSE 
          IEC = 1
          LROW_IN = NORBT
        END IF
        NKACT = NSTFGP(KACGRP)
*
        MXAADST = MAX(MXAADST,NKSTR*NORBTS)
        IF(NSTA*NSTB*NIAC*NKAC.NE.0)
     &  CALL ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,   
     &                 WORK(KSTSTM(KACGRP,1)),WORK(KSTSTM(KACGRP,2)),
     &                 IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                 NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC,IAC,
     &                 LROW_IN,IEC)
        KSTRBS = KSTRBS + NKSD     
        GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
*     JTEST = 1
      IF(JTEST.EQ.1) THEN
*. Chasing a bug, Jan 2008
        NNEG = 0
        DO IORB = 1, NORBTS
        DO KSTR = 1, NKSTR
          IF(I1((IORB-1)*NKSTR + KSTR).LT.0) THEN
           WRITE(6,*) ' Problem in ADAST2, negative adress'
           WRITE(6,*) ' IORB, KSTR, I1(KSTR,IORB) = ',
     &                  IORB, KSTR, I1((IORB-1)*NKSTR + KSTR)
          END IF
        END DO
        END DO
*
        IF(NNEG.NE.0) THEN
          CALL MEMCHK2('ADAST')
          STOP ' Problem in ADAST2' 
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADAST_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMAT(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
C!    CALL QEXIT('ADAST ')
      RETURN
      END
      SUBROUTINE TEST_RHO2S(RHO2,RHO2AA,RHO2AB,RHO2BB,NORB)
*
* Test two-body spin-density matrices by assembling 
* the standard two-body matrix from these and compare
*
* Jeppe Olsen, sitting in Oak Ridge, Sept. 2004
*
      INCLUDE 'implicit.inc'
*
      DIMENSION RHO2(*),RHO2AA(*),RHO2AB(*), RHO2BB(*)
*
*. Remember the form of the various densities : 
*
* RHO2(ij,kl) = sum_(s,s') <L!a+is a+ks' als' aj s!R> (ij.ge.kl)
* RHO2AA(iklj) = <L!a+iaa+ka aja ala!R> (i.ge.k, j.ge.l)
* RHO2BB(iklj) = <L!a+iba+kb ajb alb!R> (i.ge.k, j.ge.l)
* RHO2AB(iklj) = <L!a+iaa+kb alb aja!R> ( no restrictions)
*
      NTEST = 1000
      WRITE(6,*) ' Welcome to test of 2e-spindensities '
      WRITE(6,*) ' =================================== '
*
      NERROR = 0
      THRES = 1.0D-10
      DO J = 1, NORB
        DO L = 1, J
          DO I = 1, NORB
            IF(J.EQ.L) THEN
              MAXK = I
            ELSE
              MAXK = NORB
            END IF
            DO K = 1, MAXK
*. address in rho2
              IJ = (J-1)*NORB + I
              KL = (L-1)*NORB + K
              IJKL_RHO2 = IJ*(IJ-1)/2 + KL
*. addresses in rho2ab
              IKLJ_RHO2AB = (J-1)*NORB**3 + (L-1)*NORB**2
     &                    + (K-1)*NORB    +  I
              KIJL_RHO2AB = (L-1)*NORB**3 + (J-1)*NORB**2
     &                    + (I-1)*NORB    + K
*. adress in rho2ss
              IF(I.GT.K) THEN
                IK = I*(I-1)/2 + K
                SIGN = -1.0D0
              ELSE
                IK = K*(K-1)/2 + I
                SIGN = 1.0D0
              END IF
C             LJ = L*(L-1)/2 + J
              LJ = J*(J-1)/2 + L
              IKLJ_RHO2SS = (LJ-1)*NORB*(NORB+1)/2 + IK
*
              RHO2_1 = RHO2(IJKL_RHO2)
              RHO2_2 = RHO2AB(IKLJ_RHO2AB)
     &               + RHO2AB(KIJL_RHO2AB)
     &               + SIGN*RHO2AA(IKLJ_RHO2SS)
     &               + SIGN*RHO2BB(IKLJ_RHO2SS)
              IF(ABS(RHO2_1-RHO2_2).GT.THRES) THEN
                NERROR = NERROR + 1
                WRITE(6,*) ' Problem with spinden '
                WRITE(6,*) ' I,J,K,L= ', I,J,K,L
                WRITE(6,*) ' element from RHO2 and from RHO2S ',
     &          RHO2_1,RHO2_2
                WRITE(6,*) ' term from RHO2AB =',
     &          RHO2AB(IKLJ_RHO2AB)+ RHO2AB(KIJL_RHO2AB)
                WRITE(6,*) ' term from RHO2SS ',
     &          SIGN*(RHO2AA(IKLJ_RHO2SS)+RHO2BB(IKLJ_RHO2SS))
                WRITE(6,*) ' IKLJ_RHO2SS = ', IKLJ_RHO2SS
              END IF
            END DO
          END DO
        END DO
      END DO
*
      IF(NERROR.EQ.0) THEN
        WRITE(6,*) ' RHO2 and RHO2S are in agreement '
      ELSE
        WRITE(6,*) ' Number of differences between RHO2 and RHO2S',
     &               NERROR 
      END IF
*
      RETURN
      END
      SUBROUTINE COMHAM(H,NVAR,NBLOCK,LBLOCK,VEC1,VEC2,ECORE)
*
* Construct complete H matrix through a sequence 
* of direct CI iterations 
*
* BLocks assumed to be allocated outside
*
* Jeppe Olsen, April 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'clunit.inc'
*. Input
      INTEGER LBLOCK(NBLOCK)
*. Output
      DIMENSION H(NVAR,NVAR)
*. Scratch through argument list, should be able to hold largest SD TTS block
      DIMENSION VEC1(*), VEC2(*)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'COMHAM')
*.
      CALL MEMMAN(KLVEC1,NVAR,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,NVAR,'ADDL  ',2,'VEC2  ')
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from COMHAM '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' Number of parameters = ', NVAR
        WRITE(6,*) ' Number of blocks = ', NBLOCK
        WRITE(6,*) ' Length of each block : '
        CALL IWRTMA(LBLOCK,1,NBLOCK,1,NBLOCK)
      END IF
*
C      DO I = 1, NVAR
       DO I = 18, 18
       WRITE(6,*) ' Generating column = ', I
*. Create i'th unit vector
       ZERO = 0.0D0
       CALL SETVEC(WORK(KLVEC1),ZERO,NVAR)
       WORK(KLVEC1-1+I) = 1.0D0
*. Transfer to disc in TTS blocked form 
       CALL REWINO(LUSC1)
       CALL TODSCN(WORK(KLVEC1),NBLOCK,LBLOCK,-1,LUSC1)
C           TODSCN(VEC,NREC,LREC,LBLK,LU)
       CALL ITODS(-1,1,-1,LUSC1)
*. H * LUSC1 => LUHC
       CALL MV7(VEC1,VEC2,LUSC1,LUHC,0,0)
*. Read in He_i and save 
       CALL REWINO(LUHC)
       CALL FRMDSCN(H(1,I),NBLOCK,-1,LUHC)
       WRITE(6,*) ' Output sigma vector '
       CALL WRTVCD(VEC1,LUHC,1,-1)
       STOP ' Enforced stop in COMHAM '
C           FRMDSCN(VEC,NREC,LBLK,LU)
      END DO
*. Add core energy to Hamiltonian
      DO I = 1, NVAR
       H(I,I) = H(I,I) + ECORE
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' The complete H matrix '
        CALL WRTMAT(H,NVAR,NVAR,NVAR,NVAR)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ', IDUM,'COMHAM')
      RETURN
      END 
      SUBROUTINE IAIB_TO_ACCOCC(IATP,IAGRP,IBTP,IBGRP,IACCOCC)
*
* Numbers of IA and IB to accumulated occupation 
*
* Jeppe Olsen, April 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'gasstr.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'strbas.inc'
*. Output
      INTEGER IACCOCC(NGAS)
C?    WRITE(6,*) ' IATP, IBTP, IAGRP, IBGRP = ',
C?   &             IATP, IBTP, IAGRP, IBGRP 
      IATP_ABS = IATP + IBSPGPFTP(IAGRP) - 1
      IBTP_ABS = IBTP + IBSPGPFTP(IBGRP) - 1
C?    WRITE(6,*) ' IATP_ABS, IBTP_ABS ', IATP_ABS, IBTP_ABS
*
*. Occupation in each GAS space 
*
      IONE = 1
      CALL IVCSUM(IACCOCC,NELFSPGP(1,IATP_ABS),NELFSPGP(1,IBTP_ABS),
     &            IONE,IONE,NGAS)
C?    WRITE(6,*) ' Occupation before accumulation '
C?    CALL IWRTMA(IACCOCC,1,NGAS,1,NGAS)
*
*. And accumulated occupation 
*
      DO IGAS = 2, NGAS
        IACCOCC(IGAS) = IACCOCC(IGAS) + IACCOCC(IGAS-1)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Accumulated occupation from IAIB_TO_ACCOCC '
        CALL IWRTMA(IACCOCC,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END 
      FUNCTION IS_ACCOCC_IN_ACCOCC(IACCOCC,IACCOCC_TOT,NGAS,
     &                              NDIM_IACCOCC)
*
* Occupation of an occupation class is given in IACCOCC
* Is this occupation included in accumulated occupation IACCOCC_TOT ?
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER IACCOCC(NGAS),IACCOCC_TOT(NDIM_IACCOCC,2)
*
      INCLUDED = 1
      DO IGAS = 1, NGAS
        IF(IACCOCC(IGAS).LT.IACCOCC_TOT(IGAS,1).OR.
     &     IACCOCC(IGAS).GT.IACCOCC_TOT(IGAS,2)    ) INCLUDED = 0
      END DO
*
      IS_ACCOCC_IN_ACCOCC = INCLUDED 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  NDIM_IACCOCC = ',  NDIM_IACCOCC
        WRITE(6,*) ' Occupation to be tested : '
        CALL IWRTMA(IACCOCC,NGAS,1,NGAS,1)
        WRITE(6,*) ' Min and max of accumulated occupation '
        CALL IWRTMA(IACCOCC_TOT,NGAS,2,NDIM_IACCOCC,2)
        IF(INCLUDED.EQ.1) THEN
           WRITE(6,*) ' Occupation class is included '
        ELSE 
           WRITE(6,*) ' Occupation class is not included '
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE CMP_INFBLK_OCCSPC(INFBLK,NBLOCK,IOCCSPC,IBLK_IN_OCC,
     &                             NGAS,NVAR_IN_OCC,IFLAG)
*
* A CI space is defined by INFBLK. Check whether the blocks are
* included in accumulated occupation space and report back in
* IBLK_IN_OCC
*
* IF IFLAG = 1, only the number of parameters in space is obtained 
*
* Jeppe Olsen, April 2003
*
      INCLUDE 'implicit.inc'
      INCLUDE 'mxpdim.inc'
*. Input
      INTEGER INFBLK(8,NBLOCK), IOCCSPC(MXPNGAS,2)
*. Output 
      INTEGER IBLK_IN_OCC(NBLOCK)
*. Local scratch 
      INTEGER IOCCL(MXPNGAS)
*
      NVAR_IN_OCC = 0
      IAGRP = 1
      IBGRP = 2
      DO IBLOCK = 1, NBLOCK
*. Accumulated occupation 
        JATP = INFBLK(1,IBLOCK)
        JBTP = INFBLK(2,IBLOCK)
C?      WRITE(6,*) ' IBLOCK, JATP, JBTP = ', IBLOCK,JATP,JBTP
        CALL IAIB_TO_ACCOCC(JATP,IAGRP,JBTP,IBGRP,IOCCL)
C            IAIB_TO_ACCOCC(IAGRP,IATP,IBGRP,IBTP,IACCOCC)
*. Is accumulated occupation included ? 
        INCLUDED = IS_ACCOCC_IN_ACCOCC(IOCCL,IOCCSPC,NGAS,MXPNGAS)
C                  IS_ACCOCC_IN_ACCOCC(IACCOCC,IACCOCC_TOT,NGAS,
C    &                                              MXPNGAS)
        IF(INCLUDED.EQ.1) NVAR_IN_OCC = NVAR_IN_OCC + INFBLK(IBLOCK,8)
      
        IF(IFLAG.NE.1) IBLK_IN_OCC(IBLOCK) = INCLUDED 
      END DO
*
      NTEST = 00
      IF ( NTEST .GE. 100) THEN
        WRITE(6,*) ' Number of parameters in occspace ', NVAR_IN_OCC
        IF(IFLAG.NE.1) THEN
          WRITE(6,*) ' IBLK_IN_OCC array '
          CALL IWRTMA(IBLK_IN_OCC,1,NBLOCK,1,NBLOCK)
        END IF
      END IF
*
      RETURN
      END
      SUBROUTINE SDNUM_FOR_SELBLKS(INFBLK,NBLOCK,ISELBLK,ISELSD,NSELSD)
*
* A CI expansions is defined by INFBLK. Obtain the numbers of 
* the SD's in the blocks selected by nonvanishing numbers in ISELBLK
*
* Jeppe Olsen, April 2003
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER INFBLK(8,NBLOCK), ISELBLK(NBLOCK)
*. Output 
      INTEGER ISELSD(*)
*
      IB=1
      IB_SEL = 1
      DO IBLOCK = 1, NBLOCK
        L = INFBLK(8,IBLOCK)
        IF( ISELBLK(IBLOCK).EQ.1) THEN
          DO I = 1, L
            ISELSD(IB_SEL-1+I) = IB-1+I
          END DO
          IB_SEL = IB_SEL + L
        END IF
        IB = IB + L
      END DO
      NSELSD = IB_SEL - 1
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Addresses of SDs in selected blocks : '
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        CALL IWRTMA(ISELSD,1,NSELSD,1,NSELSD)
      END IF
*
      RETURN
      END
      SUBROUTINE DUMP_H_FOR_MRPT(H,NVAR,NBLOCK,INFBLK,LENBLK)
*
* Dump H matrices for Current MRPT program
* Dumped on file 95
*
* Jeppe Olsen, April 2003
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'cgas.inc'
*. Input
      INTEGER INFBLK(8,NBLOCK),LENBLK(NBLOCK)
      DIMENSION H(NVAR*NVAR)
*
      NTEST = 100
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DUMP_H')
*
*.  Obtain blocks in reference space 
*
      CALL MEMMAN(KLREFBLK,NBLOCK,'ADDL  ',1,'REFBLK')
      CALL CMP_INFBLK_OCCSPC(INFBLK,NBLOCK,IREFOCC_ACC,
     &     WORK(KLREFBLK),NGAS,NVAR_IN_OCC,0)
C     CMP_INFBLK_OCCSPC(INFBLK,NBLOCK,IOCCSPC,IBLK_IN_OCC,
C    &                             NGAS,NVAR_IN_OCC,IFLAG)
*. Obtain the corresponding determinants 
      CALL MEMMAN(KLREFSD,NVAR,'ADDL  ',1,'REFSD ')
C     SDNUM_FOR_SELBLKS(INFBLK,NBLOCK,ISELBLK,ISELSD,NSELSD)
      CALL SDNUM_FOR_SELBLKS(INFBLK,NBLOCK,WORK(KLREFBLK),WORK(KLREFSD),
     &                       NP)
      NQ = NVAR - NP
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' NP, NQ = ', NP,NQ 
      END IF
*. Augment list of P determinants with Q dets
C COMPL_LIST(ILIST,NIN,NTOT)
      CALL COMPL_LIST(WORK(KLREFSD),NP,NVAR)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' PQ list of determinants '
        CALL IWRTMA(WORK(KLREFSD),1,NVAR,1,NVAR)
      END IF
*. Reorder Hamiltonian to PQ order 
      CALL MEMMAN(KLHPQ,NVAR*NVAR,'ADDL  ',2,'HPQ   ')
      CALL GATMAT(WORK(KLHPQ),H,WORK(KLREFSD),NVAR,NVAR)
*. HPQ in packed form 
C   TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
*. Obtain 
*      (HPP   0 )
* H0 = (        )
*      ( 0   HQQ)
      CALL TRIPAK(WORK(KLHPQ),H,1,NVAR,NVAR)
      DO I = NP+1, NVAR
      DO J = 1, NP
        IJ = I*(I-1)/2+J
        H(IJ) = 0.0D0
      END DO
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' H0 matrix '
        CALL PRSYM(H,NVAR)
      END IF
      LU95 = 95
      CALL REWINO(LU95)
      DO I = 1, NVAR
        IJB = I*(I-1)/2
        WRITE(LU95,'(3E22.15)') (H(IJB+J),J=1,I)
      END DO
*     (0  HPQ)
* V = (      )
*     (HQP,0)
C     CALL TRIPAK(WORK(KLHPQ),H,1,NVAR,NVAR)
C     ZERO = 0.0D0
C     CALL ISETVC(H,ZERO,NP*(NP+1)/2)
C     DO I = NP+1, NVAR
C       DO J = NP+1,I
C         IJ = I*(I-1)/2 + J
C         H(IJ) = 0.0D0
C       END DO
C     END DO
*. and dump H0 to LU95
      CALL TRIPAK(WORK(KLHPQ),H,1,NVAR,NVAR)
      DO I = 1, NVAR
        IJB = I*(I-1)/2
        WRITE(LU95,'(3E22.15)') (H(IJB+J),J=1,I)
      END DO
* Mission completed 
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DUMP_H')
*
      RETURN
      END
      SUBROUTINE GATMAT(XMATO,XMATI,IGAT,NDIMO,NDIMI)
*
* MATO(I,J) = MATI(IREO(I),IREO(J))
*
      INCLUDE 'implicit.inc'
*. Input matrix
      DIMENSION XMATI(NDIMI,NDIMI)
*. Output
      DIMENSION XMATO(NDIMO,NDIMO)
*. Gather array
      INTEGER IGAT(NDIMO)
*
      DO J = 1, NDIMO
        JREO = IGAT(J)
        DO I = 1, NDIMO
          XMATO(I,J) = XMATI(IGAT(I),JREO)
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GATMAT : input and output matrices '
        CALL WRTMAT(XMATI,NDIMI,NDIMI,NDIMI,NDIMI)
        CALL WRTMAT(XMATO,NDIMO,NDIMO,NDIMO,NDIMO)
      END IF
*
      RETURN
      END
      SUBROUTINE ADIAJ_STR(IORB,JORB,ISTR_IN,NEL,LNUM,IZ,IREO,NTOOB,
     &                    ISTR_OUT,INUM,SIGN)
*
* Find string occupation (and address if LNUM.EQ.1) 
* for a+i aj |ISTR_IN> 
*
* It as assumed that the excitation is nonvanishing
*
* Jeppe Olsen, April 99
* 
      INCLUDE 'implicit.inc'
*. Input
      INTEGER ISTR_IN(NEL)
      INTEGER IZ(*),IREO(*)
*. Output
      INTEGER ISTR_OUT(NEL)
*
      NTEST = 000
C?    IF(NTEST.GE.100) THEN
C?      WRITE(6,*) 'ADIAJ_STR : Input string '
C?      CALL IWRTMA(ISTR_IN,1,NEL,1,NEL)
C?      WRITE(6,*) ' IORB, JORB ', IORB,JORB
C?    END IF
*
      IEL = 0
      IPLACED = 0
      DO IIEL = 1, NEL
        IF(ISTR_IN(IIEL).GT.IORB .AND. IPLACED.EQ.0 ) THEN
*. Add IORB
          IF(MOD(IIEL-1,2).EQ.0) THEN
            SIGNI =  1.0D0
          ELSE
            SIGNI = -1.0D0
          END IF
*
          IEL = IEL + 1
          ISTR_OUT(IEL) = IORB
          IPLACED = 1
        END IF
*
        IF( ISTR_IN(IIEL).NE.JORB) THEN
          IEL = IEL + 1
          ISTR_OUT(IEL) = ISTR_IN(IIEL)
        ELSE
          IF(MOD(IIEL-1,2).EQ.0) THEN
            SIGNJ = 1.0D0
          ELSE
            SIGNJ = -1.0D0
          END IF
        END IF
      END DO
*. Well, it could be that orbital i should be added as last elec
      IF(IPLACED.EQ.0) THEN
        ISTR_OUT(NEL) = IORB
        IF(MOD(NEL+1-1,2).EQ.0) THEN
          SIGNI = 1.0D0
        ELSE
          SIGNI = -1.0D0
        END IF
      END IF
      SIGN = SIGNI*SIGNJ
      IF(IORB.GT.JORB) SIGN = - SIGN
*
C?    WRITE(6,*) ' Output string '
C?    CALL IWRTMA(ISTR_OUT,1,NEL,1,NEL)
      IF(LNUM.EQ.1) THEN
       INUM = ISTRNM(ISTR_OUT,NTOOB,NEL,IZ,IREO,1)
      END IF
*. And address
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Output from ADIAJ_STR '
        WRITE(6,*) '======================'
        WRITE(6,*) ' Iorb and Jorb (a+(Iorb) a(Jorb)) :',IORB,JORB
        WRITE(6,*) ' Input and output strings '
        CALL IWRTMA(ISTR_IN ,1,NEL,1,NEL)
        CALL IWRTMA(ISTR_OUT,1,NEL,1,NEL)
        WRITE(6,*) ' Sign = ', sign
        IF(LNUM.EQ.1) WRITE(6,*) ' String number = ', INUM
      END IF
*
      RETURN
      END
      SUBROUTINE HCONFDIA_BBM(NAEL,NBEL,IJAGRP,IJBGRP,
     &           IASM,IATP,IAOC,NIA,IBSM,IBTP,IBOC,NIB,
     &           JASM,JATP,JAOC,NJA,JBSM,JBTP,JBOC,NJB,H,CB,SB)
*
* Outer routine for orbital conserving part of Ham times vector
*
* Jeppe Olsen, April 15, 1999 - Snowing in Aarhus
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cgas.inc'
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
      CALL QENTER('HCONF')
*. Fetch K matrix
C     CALL GTJK(XDUM1       ,H,NTOOB,XDUM3       ,IREOTS)
      CALL GTJK(WORK(KLZSCR),H,NTOOB,WORK(KLZSCR),IREOTS)
*
      CALL HCONFDIA_BBS(NAEL,NBEL,IJAGRP,IJBGRP,
     &           IASM,IATP,IAOC,NIA,IBSM,IBTP,IBOC,NIB,
     &           JASM,JATP,JAOC,NJA,JBSM,JBTP,JBOC,NJB,
     &           NTOOB,H,IOBPTS,NOBPT,IREOTS,WORK(KLZ(1)),
     &           WORK(KLZ(2)),WORK(KLZSCR),ISMFTO,NGAS,
     &           WORK(KLOCSTR(1)),WORK(KLOCSTR(2)),
     &           WORK(KLREO(1)),WORK(KLREO(2)),CB,SB) 
      CALL QEXIT('HCONF')
*
      RETURN
      END
      SUBROUTINE HCONFDIA_BBS(NAEL,NBEL,IJASPGP,IJBSPGP,
     &           IASM,IATP,IAOCC,NIA,IBSM,IBTP,IBOCC,NIB,
     &           JASM,JATP,JAOCC,NJA,JBSM,JBTP,JBOCC,NJB,
     &           NTOOB,RK,IOBPTS,NOBPT,IREOTS,
     &           IAZ,IBZ,IZSCR,ISMFTO,NGAS,
     &           JASTR_OC,JBSTR_OC,IAREO,IBREO,VECIN,VECOUT)
*
* Orbital occupation conserving part of Hamiltonian times vector.
* Part of this hamiltonian that does not conserve spin orbital 
* occupations
*

*
* The part of the Hamiltonian that conserves orbital occupations 
* is 
*
* H = sum_i h_ii + 0.5 sum_{ij}    E_{ii}E_{jj} (ii|jj) 
*                + 0.5 sum(i.ne.j) E_{ij}E_{ji} (ij|ji)
*
*
* and the part that conserves orbital occupatione but not 
* spin orbital occupations is
*
* sum(i.ne.j) a+ia aja a+jb aib (ij!ji) 
*  
* In search for a better preconditioner / H0 for EN pert, 
* Jeppe Olsen, April 99
*
*
* Notice : Present version works with supergroups from lists,
*          can therefore not work with passive/active division 

*
      INCLUDE 'implicit.inc'
      INCLUDE 'multd2h.inc'
      INCLUDE 'mxpdim.inc'
      INCLUDE 'gasstr.inc'
*.General input
      DIMENSION RK(NTOOB,NTOOB)
      INTEGER IOBPTS(MXPNGAS,*),NOBPT(*),IREOTS(*)
      INTEGER ISMFTO(*)
*. Specific input
      DIMENSION VECIN(NJA,NJB)
      DIMENSION IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
      DIMENSION IAREO(*), IBREO(*)
*. Scratch
      DIMENSION JASTR_OC(NAEL,*),JBSTR_OC(NBEL,*)
      DIMENSION IAZ(*),IBZ(*),IZSCR(*)
*. Local scratch
      INTEGER IASTR_OC(MXPORB),IBSTR_OC(MXPORB)
      INTEGER IXA(MXPORB),IXB(MXPORB),JXA(MXPORB),JXB(MXPORB)
*. Output
      DIMENSION VECOUT(NIA,NIB)
*
      NTEST =  000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from HCONFDIA_BBS'
        WRITE(6,*) ' IASM IATP ', IASM,IATP
        WRITE(6,*) ' JASM JATP ', JASM,JATP
        WRITE(6,*) ' IBSM IBTP ', IBSM,IBTP
        WRITE(6,*) ' JBSM JBTP ', JBSM,JBTP
        WRITE(6,*) ' IJASPGP, IJBSPGP', IJASPGP,IJBSPGP
        WRITE(6,*) ' NGAS = ', NGAS
        WRITE(6,*) ' Input block '
        CALL WRTMAT(VECIN,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' Initial value of sigma block '
        CALL WRTMAT(VECOUT,NIA,NIB,NIA,NIB)
*
        WRITE(6,*) 'IAOCC '
        CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
        WRITE(6,*) 'IBOCC '
        CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
        WRITE(6,*) 'JAOCC '
        CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
        WRITE(6,*) 'JBOCC '
        CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
*
      END IF
*
*
*. Obtain Reordering arrays for Ia,IB strings 
*
*. Arc weights for IA
      NTESTX = 0
      IATP_ABS = IBSPGPFTP(IJASPGP)-1+IATP
      JATP_ABS = IBSPGPFTP(IJASPGP)-1+JATP
      IBTP_ABS = IBSPGPFTP(IJBSPGP)-1+IBTP
      JBTP_ABS = IBSPGPFTP(IJBSPGP)-1+JBTP
      CALL WEIGHT_SPGP(IAZ,NGAS,NELFSPGP(1,IATP_ABS),NOBPT,
     &                 IZSCR,NTESTX)
*. Reorder array for IA strings
      CALL GETSTR_TOTSM_SPGP(IJASPGP,IATP,IASM,NAEL,NIASTR,
     &                       JASTR_OC,NTOOB,1,IAZ,IAREO)
*. Arc weight for IB
      CALL WEIGHT_SPGP(IBZ,NGAS,NELFSPGP(1,IBTP_ABS),NOBPT,
     &                 IZSCR,NTESTX)
*. Reorder array for IA strings
      CALL GETSTR_TOTSM_SPGP(IJBSPGP,IBTP,IBSM,NBEL,NIBSTR,
     &                       JBSTR_OC,NTOOB,1,IBZ,IBREO)
*
* String info for Ja, Jb : Actual string occ
*
*. Arc weight for JA
C     CALL WEIGHT_SPGP(JAZ,NGAS,JAOCC,NOBPT,ZSCR,NTESTX)
*. Occupation for JA strings
      IDUM  = 0
      CALL GETSTR_TOTSM_SPGP(IJASPGP,JATP,JASM,NAEL,NJASTR,
     &                       JASTR_OC,NTOOB,0,IDUM,IDUM)
*. Arc weight for JB
C     CALL WEIGHT_SPGP(JBZ,NGAS,JBOCC,NOBPT,ZSCR,NTEST)
*. Occupation for JB strings
      CALL GETSTR_TOTSM_SPGP(IJBSPGP,JBTP,JBSM,NBEL,NJBSTR,
     &                       JBSTR_OC,NTOOB,0,IDUM,IDUM)
*
      IJSM = MULTD2H(IASM,JASM)
*. Loop over orbital types of i and j
      DO ITP = 1, NGAS
      DO JTP = 1, NGAS
C?    WRITE(6,*) ' ITP JTP = ', ITP,JTP
*. Do a+ia a ja a+jb a ib connect string types
        IF(ITP.EQ.JTP) THEN
          IADEL = 0
          JADEL = 0
        ELSE
          IADEL = 1
          JADEL =-1
        END IF
        IBDEL = - IADEL
        JBDEL = - JADEL
*
        IAMOKAY = 1   
        DO KTP = 1, NGAS
          IF(KTP.NE.ITP.AND.KTP.NE.JTP) THEN
             IF(IAOCC(KTP).NE.JAOCC(KTP) .OR.
     &          IBOCC(KTP).NE.JBOCC(KTP)     ) IAMOKAY = 0
          END IF
        END DO
*
        IF(IAOCC(ITP).NE.JAOCC(ITP)+IADEL) IAMOKAY = 0 
        IF(IAOCC(JTP).NE.JAOCC(JTP)+JADEL) IAMOKAY = 0 
        IF(IBOCC(ITP).NE.JBOCC(ITP)+IBDEL) IAMOKAY = 0 
        IF(IBOCC(JTP).NE.JBOCC(JTP)+JBDEL) IAMOKAY = 0 
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' ITP JTP IAMOKAY : ', ITP,JTP,IAMOKAY
        END IF
        IF(IAMOKAY.EQ.1) THEN
*. Orbital range for I and J, 
          IOFF = IOBPTS(ITP,1)
          NIORB = NOBPT(ITP)
          JOFF = IOBPTS(JTP,1)
          NJORB  = NOBPT(JTP)
*
          IEND = IOFF + NIORB - 1
          JEND = JOFF + NJORB - 1
*
C?        WRITE(6,*) ' NJASTR NJBSTR ', NJASTR,NJBSTR
          DO JBSTR = 1, NJBSTR
            IZERO = 0
            CALL ISETVC(IXB,IZERO,NTOOB)
            CALL ISETVC(JXB,IZERO,NTOOB)
*. Expanded occupation in JB for i- and j- orbitals IXB, JXB
            DO IBEL = 1, NBEL
              IORB = JBSTR_OC(IBEL,JBSTR)
              IF(IOFF.LE.IORB.AND.IORB.LE.IEND) IXB(IORB)= 1
              IF(JOFF.LE.IORB.AND.IORB.LE.JEND) JXB(IORB)= 1
            END DO
*
            DO JASTR = 1, NJASTR
*. I orbitals not occupied in Ja and occupied in JB
             CALL ISETVC(IXA,IZERO,NTOOB)
             DO IEL = 1, NAEL
               IORB = JASTR_OC(IEL,JASTR)
               IF(IOFF.LE.IORB.AND.IORB.LE.IEND) IXA(IORB) = 1
             END DO
             NIACT = 0
             DO IORB = IOFF,IEND
              IF(IXA(IORB).EQ.0.AND.IXB(IORB).EQ.1) THEN
                NIACT = NIACT + 1
                IXA(NIACT) = IORB
              END IF
             END DO
*
*. Loop over J orbitals occupied on JA, Unoccupied in JB
             DO JEL = 1, NAEL
              JORB = JASTR_OC(JEL,JASTR)
              IF(JOFF.LE.JORB.AND.JORB.LE.JEND.AND.
     &           JXB(JORB).EQ.0) THEN
                JOBSM = ISMFTO(JORB)
                IOBSM = MULTD2H(IJSM,JOBSM)
                IF(MOD(JEL,2).EQ.0) THEN
                 SIGNJ = 1.0D0
                ELSE
                 SIGNJ = -1.0D0
                END IF
*. JORB is occupied in JA, unoccupied in JB 
*. Loop over Iorbital, check for sym
                DO IIORB = 1, NIACT
                 IORB = IXA(IIORB)
                 IF(ISMFTO(IORB).EQ.IOBSM.AND.IORB.NE.JORB) THEN
*. We have connection :  Find Ia and Ib strings
C                   ADIAJ_STR(IORB,JORB,ISTR_IN,NEL,LNUM,IZ,IREO,NTOOB,
C    &                    ISTR_OUT,INUM,SIGN)
                  CALL ADIAJ_STR(IORB,JORB,JASTR_OC(1,JASTR),NAEL,
     &                           1,IAZ,IAREO,NTOOB,IASTR_OC,IA,SIGNIA)
                  CALL ADIAJ_STR(JORB,IORB,JBSTR_OC(1,JBSTR),NBEL,
     &                           1,IBZ,IBREO,NTOOB,IBSTR_OC,IB,SIGNIB)
C                 XIJJI = RK(IREOTS(IORB),IREOTS(JORB))
                  XIJJI = RK(IORB,JORB)
                  VECOUT(IA,IB) = VECOUT(IA,IB)             
     &           + SIGNIA*SIGNIB*XIJJI * VECIN(JASTR,JBSTR) 
C?                WRITE(6,*) ' JASTR JBSTR IA IB', JASTR,JBSTR,IA,IB
C?                WRITE(6,*) ' SIGNIA,SIGNIB,XIJJI',SIGNIA,SIGNIB,XIJJI
C?                WRITE(6,*) ' VECOUT,VECIN', VECOUT(IA,IB),
C?   &                         VECIN(JASTR,JBSTR)
                 END IF
*                ^ End if IORB of interest
                END DO
*               ^ End of loop over IORB
              END IF
*             ^ End if JORB of interest
             END DO
*            ^ End of loop over JEL  
            END DO
*           ^ End of loop over JASTR
          END DO
*         ^ End of loop over JBSTR
        END IF
*       ^ End of combination ITP,JTP have connection
      END DO
      END DO
*     ^ End of loop over types of I and J
*
      RETURN
      END 
      SUBROUTINE HCONFINTV(LURHS,LUX,SHIFTG,SHIFT_DIAG,VECIN,VECOUT,
     &                LBLK,LUPROJ,LUPROJ2,LLUDIA)
*
* Solve  (HCONF+Shift)X = RHS
*
* Where HCONF is configuration conserving part of Hamiltonian
*
* If ICISTR.EQ.1 VECIN contains RHS, else RHS is assumed  on LURHS
* Output : solution is on LUX
*
*
* Jeppe Olsen, May 1999           
* 
*
      INCLUDE 'wrkspc.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cintfo.inc'
      INCLUDE 'clunit.inc'
      REAL*8  INPRDD
      LOGICAL CONVER
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
* SCRATCH files used : (LUSC3,LUSC34,LUSC35,LUSC37 ) <= Old
*                      LUSC36, LUSC37, LUSC38, LUSC39 <= New
*             
      INCLUDE 'oper.inc'
      INCLUDE 'crun.inc'
      INCLUDE 'cicisp.inc'
      INCLUDE 'strbas.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'stinf.inc'
      INCLUDE 'csm.inc'
      INCLUDE 'cecore.inc'
      COMMON/H_OCC_CONS/IH_OCC_CONS
      INCLUDE 'cshift.inc'
*
      EXTERNAL H0TVM
      DIMENSION ERROR(100)

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HCNINV')
*
      NTEST = 00
      IH_OCC_CONS = 1
*. Ensure that standard integrals are in WORK(KINT1)
CE    CALL SWAPVE(WORK(KINT1),WORK(KINT1O),NINT1)
*
* 2 : Solve linear set of equations
* ==================================
*
      ZERO = 0.0D0
      IF(LBLK.GT.0 ) THEN
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        STOP       ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
*
      ELSE IF(LBLK.LE.0)   THEN
*
*. Use path allowing segments of vectors
*
* corresponding eigenvalue info
        TEST = 
     &  SQRT(THRES_E) * SQRT(INPRDD(VECIN,VECOUT,LURHS,LURHS,1,-1))
        ILNPRT = NTEST
        CONVER = .FALSE.
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXIT_LOC ',MXIT_LOC
        SHIFT = SHIFTG 
        IPROJ = 0
C?      WRITE(6,*) ' LUX LURHS, LUSC38, LUSC38, LUSC39 LLUDIA',
C?   &               LUX,LURHS, LUSC38, LUSC38, LUSC39,LLUDIA 
C?      WRITE(6,*) ' LUPROJ2, LUPROJ ',  LUPROJ2, LUPROJ A
        WRITE(6,*) ' SHIFTG and SHIFT_DIAG', SHIFT,SHIFT_DIAG
C       SHIFTX = SHIFT_DIAG+ECORE_ORIG-ECORE
        SHIFTX = SHIFT_DIAG
        WRITE(6,*) ' SHIFTX before call to MICGCG ', SHIFTX
        CALL MICGCG(H0TVM,LUX,LURHS,LUSC37,LUSC38,LUSC39,LLUDIA,
     &              VECIN,VECOUT,MXIT_LOC,
     &              CONVER,TEST,SHIFTX,ERROR,NDIM,LUPROJ,LUPROJ2,
     &              VFINAL,ILNPRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Solution to linear set of Equations '
          CALL WRTVCD(VECIN,LUX,1,LBLK)
        END IF
*
      END IF
*
CE    CALL SWAPVE(WORK(KINT1),WORK(KINT1O),NINT1)
      IH_OCC_CONS = 0
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HCNINV')
      RETURN
      END 
      SUBROUTINE IB_FOR_SEL_ORBSPC(NOBPTS,NOBPS_SEL,IOBPTS_SEL,I_SEL,
     &           NGAS,MXPNGAS,NSYM)
*
* Obtain number of orbitals per symmetry and 
* offsets for orbitals in TS ordering when only 
* selected orbital spaces ( as defined by I_SEL) are included
*
*. Jeppe Olsen, Sept 2005
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NOBPTS(MXPNGAS,*),I_SEL(NGAS)
*. Output
      INTEGER IOBPTS_SEL(MXPNGAS,*), NOBPS_SEL(NSYM)
*
      IZERO = 0 
      CALL ISETVC(NOBPS_SEL,IZERO,NSYM)
      IOFF = 1
      DO IGAS = 1, NGAS
        DO ISYM = 1, NSYM
          IOBPTS_SEL(IGAS,ISYM) = IOFF
          IF(I_SEL(IGAS).EQ.1) THEN
             IOFF = IOFF + NOBPTS(IGAS,ISYM)
             NOBPS_SEL(ISYM) = NOBPS_SEL(ISYM) + NOBPTS(IGAS,ISYM)
          END IF
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Offsets for TS ordered selected orbspaces'
        CALL IWRTMA(IOBPTS_SEL,MXPNGAS,NSYM,MXPNGAS,NSYM)
        WRITE(6,*) ' Number of orbitals per sym, selected spaces'
        CALL IWRTMA(NOBPS_SEL,1,NSYM,1,NSYM)
      END IF
*
      RETURN
      END
      SUBROUTINE IREO_DACT_TS(NOBPTS,IOBPTS_SEL,I_SEL,
     &           IDTFREO,IFTDREO,NGAS,MXPNGAS,NSMOB)
*
* A set of active orbital spaces is given in I_SEL. Obtain 
* IDTFREL : Reorder array : Density order => Full TS order 
* IFTDREL : Reorder array : Full ST order => density  order 
*
      INCLUDE 'implicit.inc'
*. Input
      INTEGER NOBPTS(MXPNGAS,*),I_SEL(NGAS)
      INTEGER IOBPTS_SEL(MXPNGAS,*)
*. Output
      INTEGER IDTFREO(*),IFTDREO(*)
*
      IOB_TS = 0
      NDACTORB = 0
      DO IGAS = 1, NGAS
        DO ISYM = 1, NSMOB
          DO IOB = 1, NOBPTS(IGAS,ISYM)
            IOB_TS = IOB_TS + 1
            IF(I_SEL(IGAS).EQ.1) THEN  
              IOB_DACT = IOBPTS_SEL(IGAS,ISYM)-1+ IOB
              IFTDREO(IOB_TS) = IOB_DACT
              IDTFREO(IOB_DACT) = IOB_TS
              NDACTORB = NDACTORB + 1
            ELSE 
              IFTDREO(IOB_TS) = 0
            END IF
          END DO
        END DO
      END DO
      NOB_TOT = IOB_TS
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reorder array : FULL TS => Active '
        CALL IWRTMA(IFTDREO,NOB_TOT,1,NOB_TOT,1)
        WRITE(6,*) ' Reorder array : Active => FULL TS  '
        CALL IWRTMA(IDTFREO,NDACTORB,1,NDACTORB,1)
      END IF
*
      RETURN
      END
      SUBROUTINE EXTR_SYMBLK_ACTMAT(AIN,AOUT,IJSM)
*
* A matrix AIN is given in complete form over active orbitals, 
* symmetry ordered
*
* extract symmetry blocks with symmetry IJSM
*
* Jeppe Olsen, Feb. 2011 from REORHO1 
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      INCLUDE 'mxpdim.inc'
      INCLUDE 'lucinp.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'multd2h.inc'
*. Input
      DIMENSION AIN(NACOB,NACOB)
*. Output
      DIMENSION AOUT(*)
*
      IMOFF = 1
      DO ISM = 1, NSMOB
       JSM = MULTD2H(ISM,IJSM)
*. Offsets for active orbitals with symmetries ISM, JSM
       IOFF = 1
       DO IISM = 1, ISM -1
        IOFF = IOFF + NACOBS(IISM)
       END DO
       JOFF = 1
       DO JJSM = 1, JSM - 1
        JOFF = JOFF + NACOBS(JJSM)
       END DO
*
        NI  = NACOBS(ISM)
        NJ =  NACOBS(JSM)
        DO I = 1, NI
          DO J = 1, NJ
            AOUT(IMOFF-1+(J-1)*NI+I) = AIN(IOFF-1+I,JOFF-1+J)
          END DO
        END DO
        IMOFF = IMOFF + NI*NJ
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' EXTR_SYMBLK_ACTMAT in action '
        WRITE(6,*) ' Symmetry of blocks extracted ',IJSM
        WRITE(6,*) ' Input matrix '
        CALL WRTMAT(AIN,NACOB,NACOB,NACOB,NACOB)
        WRITE(6,*)
        WRITE(6,*) ' extracted blocks : '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
        CALL PRHONE(AOUT,NACOBS,IJSM,NSMOB,0)
      END IF
*
      RETURN
      END
