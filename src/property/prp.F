*
* $Id: prp.F,v 1.49 2000-05-05 18:55:02 mdupuis Exp $
*
      SUBROUTINE HND_PROPTY(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
      parameter  (mxefc=3*7*7*7)
      integer     rtdb, basis, geom
      logical     status
      LOGICAL     SOME
      LOGICAL     OUT 
      CHARACTER*8 WFNTYP,SCFTYP
      CHARACTER*8 SCF,RHF,UHF,MP2
      CHARACTER*8 RUNTYP
      CHARACTER*8 PROP
      CHARACTER*8 BLANK
      CHARACTER*4 IBLK 
      CHARACTER*4 IEFC 
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_RESTAR/NREST
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_RUNOPT/RUNTYP
      COMMON/HND_SYMTRY/INVT(48),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_MMTDIP/DIPOL,DMX,DMY,DMZ
      COMMON/HND_DATDPL/DPOL( 6)
      COMMON/HND_DATDHP/HPOL(10)
      COMMON/HND_DATFPL/FPOL( 6)
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      DATA SCF     /'SCF     '/
      DATA RHF     /'RHF     '/
      DATA UHF     /'UHF     '/
      DATA MP2     /'MP2     '/
      DATA ZERO    /0.0D+00/
      DATA PROP    /'PROPTY  '/
      DATA BLANK   /'        '/
      DATA IBLK    /'   '/
      DATA NODIP,NOQDP,NOOTP,NOPOP,NOSPIN,NOLOC,NODPL,NODHP,NOD2HP,
     1     NOFPL,NOFHP,NOF2HP,NOSOS,NOELP,NOELF,NOELFG,NODEN,
     2     NOGIAO,NOIGLO,NOSTON
C    3    /  0,    1,    1,    0,    1,    1,    1,    1,    1,
C    4       1,    1,    1,    1,    1,    1,    1,    1,
C    5       1,    1,    1/
     3    /  1,    1,    1,    1,    1,    1,    1,    1,    1,
     4       1,    1,    1,    1,    1,    1,    1,    1, 
     5       1,    1,    1/
C
      OUT   =.FALSE.
      SOME  =.TRUE.        
      SOME  =SOME.OR.OUT
      SOME  =SOME.and.ga_nodeid().eq.0
C
      ir=LuIN
      iw=LuOut
C
      status=rtdb_get(rtdb,'prop:dipole',MT_INT,1,nodip)
      status=rtdb_get(rtdb,'prop:quadrupole',MT_INT,1,noqdp)
      status=rtdb_get(rtdb,'prop:octupole',MT_INT,1,nootp)
      status=rtdb_get(rtdb,'prop:mulliken',MT_INT,1,nopop)
      status=rtdb_get(rtdb,'prop:spinpopulation',MT_INT,1,nospin)
      status=rtdb_get(rtdb,'prop:boyslocalization',MT_INT,1,noloc)
      status=rtdb_get(rtdb,'prop:esp',MT_INT,1,noelp)
      status=rtdb_get(rtdb,'prop:efield',MT_INT,1,noelf)
      status=rtdb_get(rtdb,'prop:efieldgrad',MT_INT,1,noelfg)
      status=rtdb_get(rtdb,'prop:electrondensity',MT_INT,1,noden)
      status=rtdb_get(rtdb,'prop:stoneanalysis',MT_INT,1,noston)
      status=rtdb_get(rtdb,'prop:iglo',MT_INT,1,noiglo)
      status=rtdb_get(rtdb,'prop:giao',MT_INT,1,nogiao)
C
      NEFC  =0
      IEFC  =IBLK
      RUNTYP=BLANK
      WFNTYP=SCF  
      SCFTYP=RHF       
C
      IF(NODHP.EQ.0.OR.NOD2HP.EQ.0) NODPL=0
      IF(NOFHP.EQ.0.OR.NOF2HP.EQ.0) NOFPL=0
      IF(SOME) WRITE(IW,9999) NODIP,NOQDP,NOOTP,NOPOP,NOSPIN,NOLOC,
     1                        NODPL,NODHP,NOD2HP,NOFPL,NOFHP,NOF2HP,
     2                        NOSOS,NOELP,NOELF,NOELFG,NODEN,NOGIAO,
     3                        NOIGLO,NOSTON
      IF(OUT ) WRITE(IW,9998) RUNTYP,SCFTYP
C
C     ----- INTERFACE OF NWCHEM AND HONDO -----
C
      CALL HND_PRP_INIT(GEOM,BASIS)
      CALL HND_PRP_NW_HND
      CALL HND_ATLABL(RTDB,GEOM,BASIS)
C
C     ----- DIPOLE MOMENT -----
C
      DIPOL=ZERO
      DMX  =ZERO
      DMY  =ZERO
      DMZ  =ZERO
      IF(NODIP.EQ.0.OR.NOLOC.EQ.0.OR.NODPL.EQ.0) 
     &  CALL HND_DIPOLE(RTDB,BASIS,GEOM)
C
C     ----- QUADRUPLOLE MOMENT -----
C
      IF(RUNTYP.NE.PROP.AND.NOQDP.EQ.0) 
     &  CALL HND_QDPOLE(RTDB,BASIS,GEOM)

C
C     ----- OCTUPOLE MOMENT -----
C
      IF(RUNTYP.NE.PROP.AND.NOOTP.EQ.0) 
     &  CALL HND_OTPOLE(RTDB,BASIS,GEOM)
C
C     ----- (HYPER)POLARIZABILITIES BY -SOS- -----
C
*     IF(NOSOS.EQ.0) CALL HND_SOSX
C
C     ----- DIPOLE POLARIZABILITY -----
C     ----- SYMMETRY TURNED OFF   -----
C
*     IF(NOSYM.NE.0) CALL HND_SYMOFF
C
*     DO I=1,6
*        DPOL(I)=ZERO
*     ENDDO
*     IF(NODPL.EQ.0.AND.NOFPL.NE.0) CALL HND_DPL
C
C     ----- DIPOLE HYPERPOLARIZABILITY -----
C
*     DO I=1,10
*        HPOL(I)=ZERO
*     ENDDO
*     IF(NODHP.EQ.0.AND.NOFHP.NE.0) CALL HND_DHP
C
C     ----- DIPOLE SECOND HYPERPOLARIZABILITY -----
C
*     IF(NOD2HP.EQ.0.AND.NOF2HP.NE.0) CALL HND_D2HP
C
C     ----- FREQUENCY DEPENDENT DIPOLE POLARIZABILITY -----
C
*     DO I=1,6
*        FPOL(I)=ZERO
*     ENDDO
*     IF(NOFPL.EQ.0) CALL HND_FPL
C
C     ----- FREQENCY DEPENDENT DIPOLE HYPERPOLARIZABILITY -----
C
*     IF(NOFHP.EQ.0) CALL HND_FHP
C
C     ----- FREQUENCY DEPENDENT DIPOLE SECOND HYPERPOLARIZABILITY -----
C
*     IF(NOF2HP.EQ.0) CALL HND_F2HP
C
C     ----- SYMMETRY TURNED ON AGAIN -----
C
*      IF(NOSYM.NE.0) CALL HND_SYMON
C
C     ----- MULLIKEN POPULATION ANALYSIS -----
C
      IF(NOPOP.EQ.0) CALL HND_MULKEN(RTDB,BASIS,GEOM)
C
C     ----- ATOMIC SPIN DENSITY -----
C
*     IF(RUNTYP.NE.PROP.AND.NOSPIN.EQ.0) 
*    1    CALL HND_SPIND(RTDB,BASIS,GEOM)
C
C     ----- ELECTROSTATIC POTENTIAL -----
C
      IF(RUNTYP.NE.PROP.AND.NOELP.EQ.0) CALL 
     1               HND_ELPMAP(RTDB,BASIS,GEOM)
C
C     ----- ELECTRIC FIELD -----
C
      IF(RUNTYP.NE.PROP.AND.NOELF.EQ.0) CALL 
     1               HND_ELFMAP(RTDB,BASIS,GEOM)
C
C     ----- ELECTRIC FIELD GRADIENT -----
C
      IF(RUNTYP.NE.PROP.AND.NOELFG.EQ.0) CALL 
     1               HND_EFGMAP(RTDB,BASIS,GEOM)
C
C     ----- ELECTRON AND SPIN DENSITIES -----
C
      IF(NODEN.EQ.0) CALL HND_ELDMAP(RTDB,BASIS,GEOM)
C
C     ----- GIAO -NMR- CHEMICAL SHIFTS -----
C
      IF(NOGIAO.EQ.0) CALL HND_GIAOX(RTDB,BASIS,GEOM)
C
C     ----- IGLO -NMR- CHEMICAL SHIFTS -----
C
      NOIGLO=1
      IF(NOIGLO.EQ.0) CALL HND_IGLOX(RTDB,BASIS,GEOM)
C
C     ----- STONE'S DISTRIBUTED MULTIPOLE ANALYSIS -----
C
      NOSTON=1
      IF(NOSTON.EQ.0) CALL HND_STONEX
C
C     ----- BOYS ORBITAL LOCALIZATION -----
C
      NOLOC=1
      IF(NOLOC.EQ.0) CALL HND_BOYLOC
      IF(NOLOC.EQ.0) CALL HND_MULKEN(RTDB,BASIS,GEOM)
C
      RETURN
 9999 FORMAT(/,' ----- property selection -----',/,
     1 ' nodip =',I2,' noqdp =',I2,' nootp =',I2,
     2 ' nopop =',I2,' nospin=',I2,
     3 ' noloc =',I2,' nodpl =',I2,/,' nodhp =',I2,' nod2hp=',I2,
     4 ' nofpl =',I2,' nofhp =',I2,' nof2hp=',I2,' nosos =',I2,
     5 ' noelp =',I2,/,' noelf =',I2,' noelfg=',I2,' noden =',I2,
     6 ' nogiao=',I2,' noiglo=',i2,' noston=',I2)
 9998 FORMAT(' runtyp = --',A8,'--',/,' scftyp = --',A8,'--')
      END
      SUBROUTINE HND_PRPDEN(DA,DB,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXIODA=255)
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_SQFILE/IJK,IPK
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      DIMENSION DA(1),DB(1)
      DATA SCF,UHF /'SCF     ','UHF     '/
      DATA MP2     /'MP2     '/
C
      CALL HND_DAREAD(IDAF,IODA,DA,L,16)
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) RETURN
      CALL HND_DAREAD(IDAF,IODA,DB,L,20)
C
      RETURN
      END
      FUNCTION hnd_prpDOT(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1)
      DATA ZERO,TWO /0.0D+00,2.0D+00/
      DUM=ZERO
      M=0
      DO 20 I=1,N
      DO 10 J=1,I
      M=M+1
   10 DUM=DUM+A(M)*B(M)*TWO
   20 DUM=DUM-A(M)*B(M)
      hnd_prpDOT=DUM
      RETURN
      END
c
      subroutine prp_cent(rtdb,geom,cm,am,c,nat)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb, geom, nat, iat, i, center
      double precision cm(3), am(nat), c(3,nat), sum
      double precision zero, tol, scale
c
      logical geom_get_user_scale, geom_check_handle
      external geom_get_user_scale, geom_check_handle
c
      DATA ZERO    /0.0D+00/
      DATA TOL     /1.0D-06/
c
      if (.not. rtdb_get(rtdb, 'prop:center', mt_int, 1, center))
     $  center = 1    ! default of center of charge
c
      if (center.eq.4) then   ! arbitrary point of expansion
        if (.not. rtdb_get(rtdb, 'prop:center_val',mt_dbl,3,cm))
     $  call errquit('prp_cent: rtdb_get center_val failed',555)
c
c get scale factor based on geometry UNITS keyword
c
        if (.not. geom_get_user_scale(geom,scale))
     $    call errquit('prp_cent: trouble getting units scale',555)
c
        do i = 1, 3             ! scale according to geometry units
          cm(i) = cm(i)*scale
        end do
        if(ga_nodeid().eq.0)
     &   write(6,9994) cm(1),cm(2),cm(3)
c
      elseif (center.eq.3) then   ! origin
        do i = 1, 3
          cm(i) = zero
        end do
        if(ga_nodeid().eq.0)
     &   write(6,9995) cm(1),cm(2),cm(3)
c
      elseif (center.eq.2) then   ! center of mass
        SUM=ZERO
        DO 20 I=1,3
   20   CM(I)=ZERO
        DO 30 IAT=1,NAT
          SUM=SUM+AM(IAT)
          DO 30 I=1,3
   30   CM(I)=CM(I)+AM(IAT)*C(I,IAT)
        IF( ABS(SUM).LT.TOL) GO TO 45
        DO 40 I=1,3
   40   CM(I)=CM(I)/SUM
   45   CONTINUE
        if(ga_nodeid().eq.0)
     &   write(6,9996) cm(1),cm(2),cm(3)
c
      elseif (center.eq.1) then   ! center of charge
        if(.not. geom_check_handle(geom, 'prp_cent'))
     $    call errquit('prp_cent:improper geometry',555)
        call geom_compute_values(geom)
        sum = zero
        do i = 1, ncenter(geom)
           sum = sum + charge(i,geom)
        end do
c
        if (sum .eq. zero) then ! System is charge neutral
          do i = 1, 3
            cm(i) = zero
          end do
        else
          do i = 1, 3
             cm(i) = ndipole(i,geom)/sum
          end do
        endif
        if(ga_nodeid().eq.0)
     &   write(6,9997) cm(1),cm(2),cm(3)
c
      else
        call errquit('prp_cent: improper center value',555)
      endif
c
      return
 9994 FORMAT(/,' Center of expansion (in au) is the arbitrary point',
     1 /,8X,' X = ',F15.7,' Y = ',F15.7,' Z = ',F15.7)
 9995 FORMAT(/,' Center of expansion (in au) is the origin',
     1 /,8X,' X = ',F15.7,' Y = ',F15.7,' Z = ',F15.7)
 9996 FORMAT(/,' Center of mass (in au) is the expansion point',
     1 /,8X,' X = ',F15.7,' Y = ',F15.7,' Z = ',F15.7)
 9997 FORMAT(/,' Center of charge (in au) is the expansion point',
     1 /,8X,' X = ',F15.7,' Y = ',F15.7,' Z = ',F15.7)
      end
c
      SUBROUTINE HND_DIPOLE(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      parameter  (mxefc=3*7*7*7)
      integer rtdb,basis,geom
      logical status
      PARAMETER   (MXIODA=255)
      PARAMETER   (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP, SCFTYP, RUNTYP
      CHARACTER*8 SCF, UHF, MP2, PROP
      LOGICAL     SOME,OUT
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/hnd_prpOPT/NOPROP
      COMMON/HND_RUNOPT/RUNTYP
      COMMON/HND_MMTDIP/DIPOL,DMX,DMY,DMZ
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION AM(MXATOM),CM(3),AXS(3,3)
      DIMENSION ERRMSG(3), dmtemp(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-DIPOLE-'/
      DATA SCF,UHF /'SCF     ','UHF     '/
      DATA MP2     /'MP2     '/
      DATA DEBYE   /2.54176568D+00/
      DATA ZERO    /0.0D+00/
      DATA TOL     /1.0D-06/
      DATA PROP    /'PROPTY  '/
      DATA KEYEFC  /4H EFC/
C
      OUT =.FALSE.
      SOME=.TRUE.
      SOME=SOME.OR.OUT
      IF(RUNTYP.EQ.PROP) SOME=.FALSE.
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM*NUM
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      NEED = L2+L2+L2+L2+L2+L3+L3
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED,'MEM INIT:prp:hnd_dipole:1',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_DIPOLE, MALLOC OF INIT  FAILED',911)
C
      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
C     ----- SET POINTERS FOR PARTITIONING OF CORE -----
C
      I10=     INIT
      I20=I10+L2
      I30=I20+L2
      I40=I30+L2
C     ----- MEMORY BLOCK FOR DA -----
      I50=I40+L2
C     ----- MOMORY BLOCK FOR DB -----
      I60=I50+L2
C     ----- MEMORY BLOCK FOR DA IN NWCHEM FORMAT -----
      I70=I60+L3
C     ----- MOMORY BLOCK FOR DA IN HONDO  FORMAT-----
      I80=I70+L3
C
      NCALL=1
   10 CONTINUE
      CALL HND_DIPAMS(AM,NCALL,NCODE,SOME)
      IF(NCODE.NE.0) GO TO 100
C
      call prp_cent(rtdb,geom,cm,am,c,nat)
C
C     ----- CALCULATE DIPOLE MOMENT INTEGRALS -----
C
      CALL HND_DIPINT(X(I10),X(I20),X(I30),CM,L1)
      CALL HND_DAWRIT(IDAF,IODA,X(I10),L2,33,NAV)
      CALL HND_DAWRIT(IDAF,IODA,X(I20),L2,34,NAV)
      CALL HND_DAWRIT(IDAF,IODA,X(I30),L2,35,NAV)
c
      if(out) then
         call hnd_prtr(x(i10),l1)
         call hnd_prtr(x(i20),l1)
         call hnd_prtr(x(i30),l1)
      endif
C
      IF(RUNTYP.EQ.PROP) GO TO 110
C
C     ----- ELECTRONIC CONTRIBUTION TO DIPOLE MOMENT -----
C
C     ----- DENSITY MATRIX FOR X(I40)=DA,X(I50)=DB -----
C           CALL hnd_prpDEN(X(I40),X(I50),L2)
C
c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i40),x(i50),x(i60),x(i70),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
c
      if(out) then
         call hnd_prtr(x(i40),l1)
         call hnd_prtr(x(i50),l1)
      endif
C
      DMX=-hnd_prpDOT(X(I40),X(I10),L1)
      DMY=-hnd_prpDOT(X(I40),X(I20),L1)
      DMZ=-hnd_prpDOT(X(I40),X(I30),L1)
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) GO TO 50
      DMX=-hnd_prpDOT(X(I50),X(I10),L1)+DMX
      DMY=-hnd_prpDOT(X(I50),X(I20),L1)+DMY
      DMZ=-hnd_prpDOT(X(I50),X(I30),L1)+DMZ
   50 CONTINUE
C
C     ----- NUCLEAR CONTRIBUTION
C
      DO 60 I=1,NAT
      DMX=DMX+ZAN(I)*(C(1,I)-CM(1))
      DMY=DMY+ZAN(I)*(C(2,I)-CM(2))
      DMZ=DMZ+ZAN(I)*(C(3,I)-CM(3))
   60 CONTINUE
C
C     ----- -EFC- CONTRIBUTION -----
C
      DMXEFC=ZERO
      DMYEFC=ZERO
      DMZEFC=ZERO
      IF(IEFC.NE.KEYEFC) GO TO 80
      DO 70 I=1,NEFC
      DMXEFC=DMXEFC+EFCZ(I)*(EFCC(1,I)-CM(1))
      DMYEFC=DMYEFC+EFCZ(I)*(EFCC(2,I)-CM(2))
      DMZEFC=DMZEFC+EFCZ(I)*(EFCC(3,I)-CM(3))
   70 CONTINUE
   80 CONTINUE
C
      DIPOL = SQRT(DMX*DMX+DMY*DMY+DMZ*DMZ)
      DIPEFC= SQRT(DMXEFC*DMXEFC+DMYEFC*DMYEFC+DMZEFC*DMZEFC)
      DMXTOT=DMX+DMXEFC
      DMYTOT=DMY+DMYEFC
      DMZTOT=DMZ+DMZEFC
      DIPTOT= SQRT(DMXTOT*DMXTOT+DMYTOT*DMYTOT+DMZTOT*DMZTOT)
      IF(SOME.and.ga_nodeid().eq.0)
     &   WRITE(IW,9995) DIPOL,DMX,DMXEFC,DMY,DMYEFC,DMZ,DMZEFC,
     1 DIPEFC,DIPTOT
c
c  write out dipole to the rtdb
c
      dmtemp(1) = dmx
      dmtemp(2) = dmy
      dmtemp(3) = dmz
      if (.not. rtdb_put(rtdb,'prop:dipval',mt_dbl,3,dmtemp))
     &  call errquit('hnd_dipole: rtdb failed for dipval',555)
c
      DMX=DMX*DEBYE
      DMY=DMY*DEBYE
      DMZ=DMZ*DEBYE
      DIPOL=DIPOL*DEBYE
      DMXEFC=DMXEFC*DEBYE
      DMYEFC=DMYEFC*DEBYE
      DMZEFC=DMZEFC*DEBYE
      DIPEFC=DIPEFC*DEBYE
      DIPTOT=DIPTOT*DEBYE
      IF(SOME.and.ga_nodeid().eq.0)
     &   WRITE(IW,9996) DIPOL,DMX,DMXEFC,DMY,DMYEFC,DMZ,DMZEFC,
     1 DIPEFC,DIPTOT
C
      NCALL=NCALL+1
      NCODE=1 
C     GO TO 10
C
  100 CONTINUE
C
C     ----- PROJECT DIPOLE COMPONENTS ONTO INERTIAL AXES -----
C
      IF(SOME.AND.NCODE.EQ.0.and.ga_nodeid().eq.0) 
     &   CALL HND_INRTIA(AM,CM,AXS,1)
      IF(SOME.AND.NCODE.EQ.0.and.ga_nodeid().eq.0)
     &   WRITE (IW,9999)
      IF(SOME.AND.NCODE.EQ.0.and.ga_nodeid().eq.0)
     &   CALL HND_DIPAXS(AXS,'A','B','C')
C
      IF(SOME.and.ga_nodeid().eq.0) WRITE(IW,9994)
C
  110 CONTINUE
C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_DIPOLE, MA_POP_STACK OF INIT FAILED',911)
C
      RETURN
C
 9999 FORMAT(/,2X,' DIPOLE MOMENTS ALONG INERTIAL AXES')
 9998 FORMAT(' NOT ENOUGH CORE IN -DIPOLE- . STOP ')
 9996 FORMAT(/,2X,' Dipole moment',F20.10,' Debye(s)',
     2 /,12X,' DMX',F20.10,' DMXEFC',F20.10,
     3 /,12X,' DMY',F20.10,' DMYEFC',F20.10,
     4 /,12X,' DMZ',F20.10,' DMZEFC',F20.10,
     5 /,2X,' -EFC- dipole ',F20.10,' DEBYE(S)',
     6 /,2X,' Total dipole ',F20.10,' DEBYE(S)')
 9995 FORMAT(/,2X,' Dipole moment',F20.10,' A.U.',
     2 /,12X,' DMX',F20.10,' DMXEFC',F20.10,
     3 /,12X,' DMY',F20.10,' DMYEFC',F20.10,
     4 /,12X,' DMZ',F20.10,' DMZEFC',F20.10,
     5 /, 2X,' -EFC- dipole ',F20.10,' A.U.',
     6 /, 2X,' Total dipole ',F20.10,' A.U.')
 9994 FORMAT(' 1 a.u. = 2.541766 Debyes ')
      END
      SUBROUTINE HND_DIPINT(XS,YS,ZS,CM,L1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDIP/XINT,YINT,ZINT,XINTX,YINTY,ZINTZ,T,XC,YC,ZC,
     1              X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      DIMENSION XS(1),YS(1),ZS(1)
      DIMENSION SX(225),SY(225),SZ(225)
      DIMENSION DIJ(225)
      DIMENSION  XIN(5,5), YIN(5,5), ZIN(5,5)
      DIMENSION XXIN(5,5),YYIN(5,5),ZZIN(5,5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION CM(3)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0E+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT=.FALSE.
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      XC=CM(1)
      YC=CM(2)
      ZC=CM(3)
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
C
      IJ=0
      DO 100 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 100 J=MINJ,JMAX
      IJ=IJ+1
      SX(IJ)=ZERO
      SY(IJ)=ZERO
      SZ(IJ)=ZERO
  100 CONTINUE
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 DIJ(IJ)=DUM2
C
C     ----- DIPOLE MOMENT INTEGRALS -----
C
      T = SQRT(AA)
      T1=ONE/T
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJT
      NJ=J
      DO 370 I=1,LIT
      NI=I
      CALL HND_DIPXYZ
       XIN(I,J)=XINT*T1
       YIN(I,J)=YINT*T1
       ZIN(I,J)=ZINT*T1
      XXIN(I,J)=XINTX*T1
      YYIN(I,J)=YINTY*T1
      ZZIN(I,J)=ZINTZ*T1
  370 CONTINUE
C
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 390 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      IJ=IJ+1
      SX(IJ)=SX(IJ)+DIJ(IJ)*(XXIN(IX,JX)* YIN(IY,JY)* ZIN(IZ,JZ))
      SY(IJ)=SY(IJ)+DIJ(IJ)*( XIN(IX,JX)*YYIN(IY,JY)* ZIN(IZ,JZ))
      SZ(IJ)=SZ(IJ)+DIJ(IJ)*( XIN(IX,JX)* YIN(IY,JY)*ZZIN(IZ,JZ))
  390 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
C     ----- SET UP DIPOLE MOMENT MATRICES -----
C
      IJ=0
      DO 7500 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 7500 J=MINJ,JMAX
      IJ=IJ+1
      LOCII=LOCI+I
      LOCJJ=LOCJ+J
      NN=(LOCII*(LOCII-1))/2+LOCJJ
      XS(NN)=SX(IJ)
      YS(NN)=SY(IJ)
      ZS(NN)=SZ(IJ)
 7500 CONTINUE
C
 8000 CONTINUE
 9000 CONTINUE
C
      IF(.NOT.OUT) GO TO 9100
      WRITE(IW,9999)
C     CALL HND_PRTRL(XS,L1)
      CALL HND_PRTR (XS,L1)
      WRITE(IW,9998)
C     CALL HND_PRTRL(YS,L1)
      CALL HND_PRTR (YS,L1)
      WRITE(IW,9997)
C     CALL HND_PRTRL(ZS,L1)
      CALL HND_PRTR (ZS,L1)
 9100 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,25(1H-),/,10X,'X-DIPOLE MOMENT INTEGRALS',/,
     1         10X,25(1H-))
 9998 FORMAT(/,10X,25(1H-),/,10X,'Y-DIPOLE MOMENT INTEGRALS',/,
     1         10X,25(1H-))
 9997 FORMAT(/,10X,25(1H-),/,10X,'Z-DIPOLE MOMENT INTEGRALS',/,
     1         10X,25(1H-))
      END
      SUBROUTINE HND_DIPXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
#include "hnd_whermt.fh"
c
      COMMON/HND_XYZDIP/XINT,YINT,ZINT,XINTX,YINTY,ZINTZ,T,XC,YC,ZC,
     1 X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      DIMENSION MIN(7),MAX(7)
      DATA MIN /1,2,4, 7,11,16,22/
      DATA MAX /1,3,6,10,15,21,28/
      DATA ZERO /0.0D+00/
C
      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      XINTX=ZERO
      YINTY=ZERO
      ZINTZ=ZERO
      NPTS=(NI+NJ-2+1)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 16 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)/T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (7,6,5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    6 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    7 GO TO (15,14,13,12,11,10,9,8),NJ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   12 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   13 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   14 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   15 CONTINUE
      XINT=XINT+PX
      YINT=YINT+PY
      ZINT=ZINT+PZ
      XINTX=XINTX+PX*(PTX-XC)
      YINTY=YINTY+PY*(PTY-YC)
      ZINTZ=ZINTZ+PZ*(PTZ-ZC)
   16 CONTINUE
      RETURN
      END
      SUBROUTINE HND_DIPAMS(BMASS,NCALL,NCODE,SOME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "nwc_const.fh"
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL SAME,SOME
      CHARACTER*8 ERRMSG
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      DIMENSION AMASS(MXATOM),BMASS(1)
      DIMENSION ERRMSG(3)
      DIMENSION AMS(106)
      DATA ERRMSG /'PROGRAM ',' STOP IN','-DIPAMS-'/
      DATA (AMS(I),I=1,54)  /
     1   1.007825D+00,  4.002600D+00,  7.016000D+00,  9.012180D+00,
     2  11.009310D+00, 12.000000D+00, 14.003070D+00, 15.994910D+00,
     3  18.998400D+00, 19.992440D+00, 22.989800D+00, 23.985040D+00,
     4  26.981530D+00, 27.976930D+00, 30.973760D+00, 31.972070D+00,
     5  34.968850D+00, 39.948000D+00, 38.963710D+00, 39.962590D+00,
     6  44.955920D+00, 47.900000D+00, 50.944000D+00, 51.940500D+00,
     7  54.938100D+00, 55.934900D+00, 58.933200D+00, 57.935300D+00,
     8  62.929800D+00, 63.929100D+00, 68.925700D+00, 73.921900D+00,
     9  74.921600D+00, 79.916500D+00, 78.918300D+00, 83.911500D+00,
     1  84.911700D+00, 87.905600D+00, 89.905400D+00, 89.904300D+00,
     2  92.906000D+00, 97.905500D+00, 97.000000D+00,101.903700D+00,
     3 102.904800D+00,105.903200D+00,106.904100D+00,113.903600D+00,
     4 114.904100D+00,119.902200D+00,120.903800D+00,129.906700D+00,
     5 126.904400D+00,131.904200D+00/
      DATA (AMS(I),I=55,106)  /
     1   132.9054D+00,137.9052D+00,138.9063D+00,139.9054D+00,
     2   140.9076D+00,141.9077D+00,144.9127D+00,151.9197D+00,
     3   152.9212D+00,157.9241D+00,158.9253D+00,163.9292D+00,
     4   164.9303D+00,165.9303D+00,168.9342D+00,173.9389D+00,
     5   174.9408D+00,179.9465D+00,180.9480D+00,183.9509D+00,
     6   186.9557D+00,191.9615D+00,192.9629D+00,194.9648D+00,
     7   196.9665D+00,201.9706D+00,204.9744D+00,207.9766D+00,
     8   208.9804D+00,208.9824D+00,209.9871D+00,222.0176D+00,
     9   223.0197D+00,226.0254D+00,227.0278D+00,232.0381D+00,
     1   231.0359D+00,238.0508D+00,237.0482D+00,244.0642D+00,
     2   243.0614D+00,247.0703D+00,247.0703D+00,251.0796D+00,
     3   252.0829D+00,257.0751D+00,258.0986D+00,259.1009D+00,
     4   260.1053D+00,261.1087D+00,  0.0000D+00,  0.0000D+00/
      DATA ZERO /0.0D+00/
C
      NAMELIST /MASS/ AMASS
C
      DO 10 IAT=1,NAT
         NUCZ=NUC(IAT)
         IF(NUCZ.GT.104) THEN
            WRITE(IW,9996)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(NUCZ.EQ.0) ZNUC=ZERO
         IF(NUCZ.GT.0) ZNUC=AMS(NUCZ)
   10    AMASS(IAT)=ZNUC
      IF(NCALL.EQ.1) GO TO 30
      IF(NCALL.GT.2)  GO TO 60
      IF(NCALL.EQ.2) REWIND IR
c%%
c%%   READ(IR,MASS,END=60,ERR=60)
c%%
C
      SAME=.TRUE.
      DO 20 IAT=1,NAT
      NUCZ=NUC(IAT)
      IF(NUCZ.EQ.0) ZNUC=ZERO
      IF(NUCZ.GT.0) ZNUC=AMS(NUCZ)
      IF(AMASS(IAT).EQ.ZNUC) GO TO 20
      SAME=.FALSE.
      GO TO 25
   20 CONTINUE
   25 CONTINUE
      IF(SAME) GO TO 60
C
   30 CONTINUE
      DO 40 IAT=1,NAT
   40 BMASS(IAT)= ABS(AMASS(IAT))
      if(some.and.ga_nodeid().eq.0) then
        WRITE(IW,9998)
        DO 50 IAT=1,NAT
   50   WRITE(IW,9997) IAT,ANAM(IAT),BNAM(IAT),BMASS(IAT)
      endif 
      NCODE=0
      RETURN
   60 NCODE=1
      RETURN
 9998 FORMAT(/,10X,21(1H-),/,10X,'Atomic weights (a.u.)',
     1       /,10X,21(1H-),/)
 9997 FORMAT(I5,5X,A8,A2,F15.5)
 9996 FORMAT(' NO ATOMIC MASS TABULATED FOR  -NUCZ.GT.104-  . STOP')
      END
      SUBROUTINE HND_INRTIA(AM,CM,PAXS,NOCOFM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
C
C     ----- PROGRAM TO GET INERTIAL AXES. CONSTRUCTS INERTIA -----
C           TENSOR, THEN CALLS DIAGONALIZATION ROUTINE TO GET
C           DIRECTION COSINES FOR PRINCIPAL INERTIA AXES AND
C           INERTIA MOMENTS. G. HURST, JUNE 1987.
C
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (ZERO = 0.0D+00)
      LOGICAL SOME
C     COMMON/HND_INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(MXATOM),
C    1                                         C(3,MXATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_OUTPUT/NPRINT
      DIMENSION AM(1),CM(3),PAXS(3,3)
      DIMENSION A(6),PMOM(3)
      CHARACTER*1 DXYZ(3),DABC(3)
      SAVE NRTFLG
      DATA NRTFLG /0/
      DATA DXYZ   /'X','Y','Z'/
      DATA DABC   /'A','B','C'/
C
      SOME=NPRINT.NE.-5
      IF(NRTFLG.GT.0) SOME=.FALSE.
      DO 10 I=1,6
         A(I)=ZERO
   10 CONTINUE
C
C     ----- CALCULATE INERTIA TENSOR (WRT TO CENTRE OF MASS) -----
C
      DO 20 IAT=1,NAT
         WT = AM(IAT)
         X = C(1,IAT) - CM(1)
         Y = C(2,IAT) - CM(2)
         Z = C(3,IAT) - CM(3)
         A(1) = A(1) + WT*( Y*Y + Z*Z )
         A(2) = A(2) - WT*X*Y
         A(3) = A(3) + WT*( X*X + Z*Z )
         A(4) = A(4) - WT*X*Z
         A(5) = A(5) - WT*Y*Z
         A(6) = A(6) + WT*( X*X + Y*Y )
   20 CONTINUE
C
C     ----- DIAGONALIZE TO FIND MOMENTS AND PRINCIPLE AXES -----
C
      NO=3
      CALL HND_DIAAXS(A,PAXS,PMOM,IA,NO,NO,NO)
      IF(SOME.AND.NOCOFM.EQ.0) WRITE(IW,9998) CM(1),CM(2),CM(3)
      IF(SOME) WRITE(IW,9999)
      IF(SOME) CALL HND_PRAXS(PAXS,PMOM,DXYZ,DABC,NO)
      NRTFLG = 1
      RETURN
 9999 FORMAT(/,8X,' PRINCIPLE MOMENTS OF INERTIA AND PRINCIPLE AXES',
     1       /,8X,'       (TRANSPOSE OF AXIS ROTATION MATRIX)')
 9998 FORMAT(/,2X,' CENTER OF EXPANSION',
     1 ' X = ',F15.7,' Y = ',F15.7,' Z = ',F15.7)
      END
      SUBROUTINE HND_DIAAXS(A,VEC,EIG,IA,NVEC,N,NDIM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     ----- ROUTINE TO SUBSTITUTE DIAGIV FOR DIAGONALIZATION -----
C           OF SYMMETRIC 3X3 MATRIX A  IN TRIANGULAR FORM
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION A(6),H(3,3),VEC(3,3),EIG(3),IA(1),BIG(6)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-DIAAXS-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA CONV   /1.0D-10/
      DATA MAXIT  /50/
C
      IF(N.EQ.3.AND.NDIM.EQ.3) GO TO 10
         WRITE(IW,9999)
         CALL HND_HNDERR(3,ERRMSG)
   10 CONTINUE
C
      DO 30 I=1,3
         DO 20 J=1,3
            VEC(J,I)=ZERO
   20    CONTINUE
         VEC(I,I)=ONE
   30 CONTINUE
      IJ=0
      DO 40 I=1,3
      DO 40 J=1,I
         IJ=IJ+1
         H(I,J)=A(IJ)
   40    H(J,I)=A(IJ)
      CALL HND_DIAJAC(H,VEC,EIG,N,N,BIG)
C
C     ----- CHECK FOR RIGHT HANDEDNESS, CORRECT IF NOT -----
C
      TEST =   VEC(1,3)*( VEC(2,1)*VEC(3,2) - VEC(3,1)*VEC(2,2) )
     1       + VEC(2,3)*( VEC(3,1)*VEC(1,2) - VEC(1,1)*VEC(3,2) )
     2       + VEC(3,3)*( VEC(1,1)*VEC(2,2) - VEC(2,1)*VEC(1,2) )
      IF(TEST.GT.ZERO) RETURN
      IF( ABS(EIG(1)-EIG(2)).GT.CONV) GO TO 60
         T = EIG(1)
         EIG(1) = EIG(2)
         EIG(2) = T
         DO 50 I=1,3
            T = VEC(I,1)
            VEC(I,1) = VEC(I,2)
            VEC(I,2) = T
   50    CONTINUE
         RETURN
   60 IF( ABS(EIG(2)-EIG(3)).GT.CONV) GO TO 80
         T = EIG(2)
         EIG(2) = EIG(3)
         EIG(3) = T
         DO 70 I=1,3
            T = VEC(I,2)
            VEC(I,2) = VEC(I,3)
            VEC(I,3) = T
   70    CONTINUE
         RETURN
   80 DO 90 I=1,3
         VEC(I,3) = - VEC(I,3)
   90 CONTINUE
      RETURN
 9999 FORMAT(/,' -DIAAXS- DIAGONALIZATION ONLY SET UP FOR 3X3 MATRIX')
      END
      SUBROUTINE HND_PRAXS(V,E,DIR1,DIR2,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT OUT AXIS ROTATION MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
      DIMENSION V(NDIM,1),E(1)
      CHARACTER*1 DIR1(3), DIR2(3)
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9068) (E(I),I = 1,3)
      IF(LIST.EQ.1) WRITE (IW,9168) (E(I),I = 1,3)
      IF(LIST.EQ.2) WRITE (IW,9268) (E(I),I = 1,3)
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (DIR2(I),I = 1,3)
      IF(LIST.EQ.1) WRITE (IW,9128) (DIR2(I),I = 1,3)
      IF(LIST.EQ.2) WRITE (IW,9228) (DIR2(I),I = 1,3)
      WRITE (IW,9008)
      DO 120 J = 1,3
      IF(LIST.EQ.0) WRITE (IW,9048) J,DIR1(J),(V(J,I),I = 1,3)
      IF(LIST.EQ.1) WRITE (IW,9148) J,DIR1(J),(V(J,I),I = 1,3)
      IF(LIST.EQ.2) WRITE (IW,9248) J,DIR1(J),(V(J,I),I = 1,3)
  120 CONTINUE
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(15X,10(7X,A1,3X))
 9048 FORMAT(I5,5X,A1,4X,10F11.5)
 9068 FORMAT(15X,10F11.3)
 9128 FORMAT(15X,7(9X,A1,5X))
 9148 FORMAT(I5,5X,A1,4X,7F15.10)
 9168 FORMAT(15X,7F15.10)
 9228 FORMAT(15X,7(9X,A1,5X))
 9248 FORMAT(I5,5X,A1,4X,7E15.8)
 9268 FORMAT(15X,7E15.8)
      END
      SUBROUTINE HND_DIPAXS(A,D1,D2,D3)
C
C     ----- PROGRAM TO ROTATE DIPOLE TENSOR COORDINATES -----
C           MATRIX A IS TRANSPONSE OF ROTATION MATRIX
C           OF DIRECTION COSINES. G. HURST, JUNE 1987.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (DEBYE = 2.54176568D+00)
      PARAMETER (ZERO = 0.0D+00)
      COMMON/HND_MMTDIP/DIPOL,DMX,DMY,DMZ
      COMMON/HND_IOFILE/IR,IW,IP
      CHARACTER*1 D1,D2,D3,D(3)
      DIMENSION A(3,3),DIP(3)
      D(1)=D1
      D(2)=D2
      D(3)=D3
      DIP(1)=DMX/DEBYE
      DIP(2)=DMY/DEBYE
      DIP(3)=DMZ/DEBYE
      WRITE(IW,9000)
      DO 20 NA=1,3
         ZZ=ZERO
         DO 10 I=1,3
            ZZ = ZZ + A(I,NA)*DIP(I)
   10    CONTINUE
         WRITE(IW,9001) D(NA),ZZ,ZZ*DEBYE
   20 CONTINUE
      WRITE(IW,9002)
      RETURN
 9000 FORMAT(/10X,'COMPONENT',9X,'BOHR*E',15X,'DEBYE')
 9001 FORMAT(13X,'MU ',A1,2F20.10)
 9002 FORMAT(1X)
      END
C
      SUBROUTINE HND_QDPOLE(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      parameter  (mxefc=3*7*7*7)
      integer rtdb,basis,geom
      logical status
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2
      LOGICAL OUT,QDPOUT
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION AM(MXATOM),CM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-QDPOLE'/
      DATA SCF,UHF /'SCF     ','UHF     '/
      DATA MP2     /'MP2     '/
      DATA BUCK    /1.344911D+00/
      DATA TOL     /1.0D-06/
      DATA ZERO    /0.0D+00/
      DATA TWO     /2.0D+00/
      DATA THREE   /3.0D+00/
      DATA KEYEFC  /4H EFC/
C
      QDPOUT=.TRUE.
      IF(.NOT.QDPOUT) RETURN
      OUT=.FALSE.      
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM*NUM
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      NEED = L2+L2+L2+L2+L2+L2+L2+L2+L3+L3
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED , 'MEM INIT:prp:hnd_qdpole:2',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_QDPOLE, MALLOC OF INIT  FAILED',911)

      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
C     ----- SET POINTERS FOR PARTITIONING OF CORE -----
C
      I10=  INIT
      I20=I10+L2
      I30=I20+L2
      I40=I30+L2
      I50=I40+L2
      I60=I50+L2
      I70=I60+L2
C     ----- MEMORY BLOCK FOR DENSITY MATRIX -----
      I80=I70+L2
      I90=I80+L2
C     ----- MEMORY BLOCK DA IN NWCHEM FMT -----
      I91=I90+L3
C     ----- MEMORY BLOCK DA IN NWCHEM FMT -----
      I92=I91+L3
C
      NCALL=1
   10 CONTINUE
      CALL HND_QDPAMS(AM,NCALL,NCODE)
      IF(NCODE.NE.0) GO TO 100
C
      call prp_cent(rtdb,geom,cm,am,c,nat)
C
C     ----- CALCULATE SECOND MOMENT INTEGRALS -----
C
      CALL HND_QDPINT(X(I10),X(I20),X(I30),
     1                X(I40),X(I50),X(I60),CM,L1)
C
C     ----- ELECTRONIC CONTRIBUTION TO SECOND MOMENTS -----
C
C*    CALL hnd_prpDEN(X(I70),X(I80),L2)
C     
C     ----- DENSITY MATRIX FOR X(I70)=DA,X(I80)=DB -----

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i70),x(i80),x(i90),x(i91),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
C
      QMXX=-hnd_prpDOT(X(I70),X(I10),L1)
      QMYY=-hnd_prpDOT(X(I70),X(I20),L1)
      QMZZ=-hnd_prpDOT(X(I70),X(I30),L1)
      QMXY=-hnd_prpDOT(X(I70),X(I40),L1)
      QMXZ=-hnd_prpDOT(X(I70),X(I50),L1)
      QMYZ=-hnd_prpDOT(X(I70),X(I60),L1)
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) GO TO 50
      QMXX=-hnd_prpDOT(X(I80),X(I10),L1)+QMXX
      QMYY=-hnd_prpDOT(X(I80),X(I20),L1)+QMYY
      QMZZ=-hnd_prpDOT(X(I80),X(I30),L1)+QMZZ
      QMXY=-hnd_prpDOT(X(I80),X(I40),L1)+QMXY
      QMXZ=-hnd_prpDOT(X(I80),X(I50),L1)+QMXZ
      QMYZ=-hnd_prpDOT(X(I80),X(I60),L1)+QMYZ
   50 CONTINUE
C
C     ----- < R**2 > = DIAMAGNETIC SUSCEPTIBILITY -----
C
      RSQUAR=-(QMXX+QMYY+QMZZ)
      if (ga_nodeid().eq.0) WRITE(IW,9992) RSQUAR
C
C     ----- NUCLEAR CONTRIBUTION
C
      DO 60 IAT=1,NAT
      QMXX=QMXX+ZAN(IAT)*(C(1,IAT)-CM(1))*(C(1,IAT)-CM(1))
      QMYY=QMYY+ZAN(IAT)*(C(2,IAT)-CM(2))*(C(2,IAT)-CM(2))
      QMZZ=QMZZ+ZAN(IAT)*(C(3,IAT)-CM(3))*(C(3,IAT)-CM(3))
      QMXY=QMXY+ZAN(IAT)*(C(1,IAT)-CM(1))*(C(2,IAT)-CM(2))
      QMXZ=QMXZ+ZAN(IAT)*(C(1,IAT)-CM(1))*(C(3,IAT)-CM(3))
      QMYZ=QMYZ+ZAN(IAT)*(C(2,IAT)-CM(2))*(C(3,IAT)-CM(3))
   60 CONTINUE
C
C     ----- -EFC- CONTRIBUTION -----
C
      QMXXEF=ZERO
      QMYYEF=ZERO
      QMZZEF=ZERO
      QMXYEF=ZERO
      QMXZEF=ZERO
      QMYZEF=ZERO
      IF(IEFC.NE.KEYEFC) GO TO 80
      DO 70 IAT=1,NEFC
      QMXXEF=QMXXEF+EFCZ(IAT)*(EFCC(1,IAT)-CM(1))*(EFCC(1,IAT)-CM(1))
      QMYYEF=QMYYEF+EFCZ(IAT)*(EFCC(2,IAT)-CM(2))*(EFCC(2,IAT)-CM(2))
      QMZZEF=QMZZEF+EFCZ(IAT)*(EFCC(3,IAT)-CM(3))*(EFCC(3,IAT)-CM(3))
      QMXYEF=QMXYEF+EFCZ(IAT)*(EFCC(1,IAT)-CM(1))*(EFCC(2,IAT)-CM(2))
      QMXZEF=QMXZEF+EFCZ(IAT)*(EFCC(1,IAT)-CM(1))*(EFCC(3,IAT)-CM(3))
      QMYZEF=QMYZEF+EFCZ(IAT)*(EFCC(2,IAT)-CM(2))*(EFCC(3,IAT)-CM(3))
   70 CONTINUE
   80 CONTINUE
C
      QMXXT =QMXX+QMXXEF
      QMYYT =QMYY+QMYYEF
      QMZZT =QMZZ+QMZZEF
      QMXYT =QMXY+QMXYEF
      QMXZT =QMXZ+QMXZEF
      QMYZT =QMYZ+QMYZEF
      if (ga_nodeid().eq.0) then
         WRITE(IW,9999) QMXX,QMXXEF,QMXXT,QMYY,QMYYEF,QMYYT,
     1               QMZZ,QMZZEF,QMZZT,QMXY,QMXYEF,QMXYT,
     2               QMXZ,QMXZEF,QMXZT,QMYZ,QMYZEF,QMYZT
         WRITE(IW,9995) QMXX*BUCK,QMXXEF*BUCK,QMXXT*BUCK,
     1               QMYY*BUCK,QMYYEF*BUCK,QMYYT*BUCK,
     2               QMZZ*BUCK,QMZZEF*BUCK,QMZZT*BUCK,
     3               QMXY*BUCK,QMXYEF*BUCK,QMXYT*BUCK,
     4               QMXZ*BUCK,QMXZEF*BUCK,QMXZT*BUCK,
     5               QMYZ*BUCK,QMYZEF*BUCK,QMYZT*BUCK
      endif
      DUM=QMXX+QMYY+QMZZ
      THXX=(THREE*QMXX-DUM)/TWO
      THYY=(THREE*QMYY-DUM)/TWO
      THZZ=(THREE*QMZZ-DUM)/TWO
      THXY= THREE*QMXY     /TWO
      THXZ= THREE*QMXZ     /TWO
      THYZ= THREE*QMYZ     /TWO
      DUM=QMXXEF+QMYYEF+QMZZEF
      THXXEF=(THREE*QMXXEF-DUM)/TWO
      THYYEF=(THREE*QMYYEF-DUM)/TWO
      THZZEF=(THREE*QMZZEF-DUM)/TWO
      THXYEF= THREE*QMXYEF     /TWO
      THXZEF= THREE*QMXZEF     /TWO
      THYZEF= THREE*QMYZEF     /TWO
      DUM=QMXXT +QMYYT +QMZZT
      THXXT =(THREE*QMXXT -DUM)/TWO
      THYYT =(THREE*QMYYT -DUM)/TWO
      THZZT =(THREE*QMZZT -DUM)/TWO
      THXYT = THREE*QMXYT      /TWO
      THXZT = THREE*QMXZT      /TWO
      THYZT = THREE*QMYZT      /TWO
      if (ga_nodeid().eq.0) then
         WRITE(IW,9997) THXX,THXXEF,THXXT,THYY,THYYEF,THYYT,
     1               THZZ,THZZEF,THZZT,THXY,THXYEF,THXYT,
     2               THXZ,THXZEF,THXZT,THYZ,THYZEF,THYZT
         WRITE(IW,9994) THXX*BUCK,THXXEF*BUCK,THXXT*BUCK,
     1               THYY*BUCK,THYYEF*BUCK,THYYT*BUCK,
     2               THZZ*BUCK,THZZEF*BUCK,THZZT*BUCK,
     3               THXY*BUCK,THXYEF*BUCK,THXYT*BUCK,
     4               THXZ*BUCK,THXZEF*BUCK,THXZT*BUCK,
     5               THYZ*BUCK,THYZEF*BUCK,THYZT*BUCK
      endif
C
      NCALL=NCALL+1
C     GO TO 10
C
  100 CONTINUE
C
      if (ga_nodeid().eq.0) WRITE(IW,9993)
C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_QDPOLE, MA_POP_STACK OF INIT
     & FAILED',911)
C
      RETURN
C
 9999 FORMAT(/,2X,' Second moments in atomic units',
     2 /,2X,' QMXX',F20.10,' QMXXEF',F20.10,' QMXXT',F20.10,
     2 /,2X,' QMYY',F20.10,' QMYYEF',F20.10,' QMYYT',F20.10,
     2 /,2X,' QMZZ',F20.10,' QMZZEF',F20.10,' QMZZT',F20.10,
     2 /,2X,' QMXY',F20.10,' QMXYEF',F20.10,' QMXYT',F20.10,
     2 /,2X,' QMXZ',F20.10,' QMXZEF',F20.10,' QMXZT',F20.10,
     2 /,2X,' QMYZ',F20.10,' QMYZEF',F20.10,' QMYZT',F20.10)
 9997 FORMAT(/,2X,' Quadrupole moments in atomic units',
     2 /,2X,' THXX',F20.10,' THXXEF',F20.10,' THXXT',F20.10,
     2 /,2X,' THYY',F20.10,' THYYEF',F20.10,' THYYT',F20.10,
     2 /,2X,' THZZ',F20.10,' THZZEF',F20.10,' THZZT',F20.10,
     2 /,2X,' THXY',F20.10,' THXYEF',F20.10,' THXYT',F20.10,
     2 /,2X,' THXZ',F20.10,' THXZEF',F20.10,' THXZT',F20.10,
     2 /,2X,' THYZ',F20.10,' THYZEF',F20.10,' THYZT',F20.10)
 9996 FORMAT(' Not enough core in -qdpole- . stop ')
 9995 FORMAT(/,2X,' Second moments in buckingham(s)',
     2 /,2X,' QMXX',F20.10,' QMXXEF',F20.10,' QMXXT',F20.10,
     2 /,2X,' QMYY',F20.10,' QMYYEF',F20.10,' QMYYT',F20.10,
     2 /,2X,' QMZZ',F20.10,' QMZZEF',F20.10,' QMZZT',F20.10,
     2 /,2X,' QMXY',F20.10,' QMXYEF',F20.10,' QMXYT',F20.10,
     2 /,2X,' QMXZ',F20.10,' QMXZEF',F20.10,' QMXZT',F20.10,
     2 /,2X,' QMYZ',F20.10,' QMYZEF',F20.10,' QMYZT',F20.10)
 9994 FORMAT(/,2X,' Quadrupole moments in buckingham(s)',
     2 /,2X,' THXX',F20.10,' THXXEF',F20.10,' THXXT',F20.10,
     2 /,2X,' THYY',F20.10,' THYYEF',F20.10,' THYYT',F20.10,
     2 /,2X,' THZZ',F20.10,' THZZEF',F20.10,' THZZT',F20.10,
     2 /,2X,' THXY',F20.10,' THXYEF',F20.10,' THXYT',F20.10,
     2 /,2X,' THXZ',F20.10,' THXZEF',F20.10,' THXZT',F20.10,
     2 /,2X,' THYZ',F20.10,' THYZEF',F20.10,' THYZT',F20.10)
 9993 FORMAT(/,' 1 a.u. = 1.344911 Buckinghams ',
     1                 '= 1.344911 10**(-26) esu*cm**2 ')
 9992 FORMAT(/,' < R**2 > = ',F10.6,' a.u. ',
     1 ' ( 1 a.u. = 0.280023 10**(-16) cm**2 ) ',/,
     2 ' ( also called diamagnetic susceptibility ) ')
      END
      SUBROUTINE HND_QDPINT(XXS,YYS,ZZS,XYS,XZS,YZS,CM,L1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZQDP/XINT,YINT,ZINT,XINTX,YINTY,ZINTZ,
     1                  XINTXX,YINTYY,ZINTZZ,T,XC,YC,ZC,
     2                  X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      DIMENSION XXS(1),YYS(1),ZZS(1),XYS(1),XZS(1),YZS(1)
      DIMENSION SXX(225),SYY(225),SZZ(225),SXY(225),SXZ(225),SYZ(225)
      DIMENSION DIJ(225)
      DIMENSION   XIN(5,5),  YIN(5,5),  ZIN(5,5)
      DIMENSION  XXIN(5,5), YYIN(5,5), ZZIN(5,5)
      DIMENSION XXXIN(5,5),YYYIN(5,5),ZZZIN(5,5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION CM(3)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT=.FALSE.        
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      XC=CM(1)
      YC=CM(2)
      ZC=CM(3)
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
C
      IJ=0
      DO 100 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 100 J=MINJ,JMAX
      IJ=IJ+1
      SXX(IJ)=ZERO
      SYY(IJ)=ZERO
      SZZ(IJ)=ZERO
      SXY(IJ)=ZERO
      SXZ(IJ)=ZERO
      SYZ(IJ)=ZERO
  100 CONTINUE
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 DIJ(IJ)=DUM2
C
C     ----- QUADRUPOLE MOMENT INTEGRALS -----
C
      T = SQRT(AA)
      T1=ONE/T
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJT
      NJ=J
      DO 370 I=1,LIT
      NI=I
      CALL HND_QDPXYZ
        XIN(I,J)=XINT*T1
        YIN(I,J)=YINT*T1
        ZIN(I,J)=ZINT*T1
       XXIN(I,J)=XINTX*T1
       YYIN(I,J)=YINTY*T1
       ZZIN(I,J)=ZINTZ*T1
      XXXIN(I,J)=XINTXX*T1
      YYYIN(I,J)=YINTYY*T1
      ZZZIN(I,J)=ZINTZZ*T1
  370 CONTINUE
C
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 390 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      IJ=IJ+1
      SXX(IJ)=SXX(IJ)+DIJ(IJ)*(XXXIN(IX,JX)*  YIN(IY,JY)*  ZIN(IZ,JZ))
      SYY(IJ)=SYY(IJ)+DIJ(IJ)*(  XIN(IX,JX)*YYYIN(IY,JY)*  ZIN(IZ,JZ))
      SZZ(IJ)=SZZ(IJ)+DIJ(IJ)*(  XIN(IX,JX)*  YIN(IY,JY)*ZZZIN(IZ,JZ))
      SXY(IJ)=SXY(IJ)+DIJ(IJ)*( XXIN(IX,JX)* YYIN(IY,JY)*  ZIN(IZ,JZ))
      SXZ(IJ)=SXZ(IJ)+DIJ(IJ)*( XXIN(IX,JX)*  YIN(IY,JY)* ZZIN(IZ,JZ))
      SYZ(IJ)=SYZ(IJ)+DIJ(IJ)*(  XIN(IX,JX)* YYIN(IY,JY)* ZZIN(IZ,JZ))
  390 CONTINUE
 6000 CONTINUE
 7000 CONTINUE
C
C     ----- SET UP QUADRUPOLE MOMENT MATRICES -----
C
      IJ=0
      DO 7500 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 7500 J=MINJ,JMAX
      IJ=IJ+1
      IDUM=LOCI+I
      JDUM=LOCJ+J
      NN=(IDUM*(IDUM-1))/2+JDUM
      XXS(NN)=SXX(IJ)
      YYS(NN)=SYY(IJ)
      ZZS(NN)=SZZ(IJ)
      XYS(NN)=SXY(IJ)
      XZS(NN)=SXZ(IJ)
      YZS(NN)=SYZ(IJ)
 7500 CONTINUE
C
 8000 CONTINUE
 9000 CONTINUE
C
      IF(.NOT.OUT) GO TO 9100
      WRITE(IW,9999)
      CALL HND_PRTR(XXS,L1)
      WRITE(IW,9998)
      CALL HND_PRTR(YYS,L1)
      WRITE(IW,9997)
      CALL HND_PRTR(ZZS,L1)
      WRITE(IW,9996)
      CALL HND_PRTR(XYS,L1)
      WRITE(IW,9995)
      CALL HND_PRTR(XZS,L1)
      WRITE(IW,9994)
      CALL HND_PRTR(YZS,L1)
 9100 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,26(1H-),/,10X,'XX-SECOND MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9998 FORMAT(/,10X,26(1H-),/,10X,'YY-SECOND MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9997 FORMAT(/,10X,26(1H-),/,10X,'ZZ-SECOND MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9996 FORMAT(/,10X,26(1H-),/,10X,'XY-SECOND MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9995 FORMAT(/,10X,26(1H-),/,10X,'XZ-SECOND MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9994 FORMAT(/,10X,26(1H-),/,10X,'YZ-SECOND MOMENT INTEGRALS',/,
     1         10X,26(1H-))
      END
      SUBROUTINE HND_QDPXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
#include "hnd_whermt.fh"
c
      COMMON/HND_XYZQDP/XINT,YINT,ZINT,XINTX,YINTY,ZINTZ,
     1 XINTXX,YINTYY,ZINTZZ,T,XC,YC,ZC,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      DIMENSION MIN(7),MAX(7)
      DATA MIN /1,2,4, 7,11,16,22/
      DATA MAX /1,3,6,10,15,21,28/
      DATA ZERO /0.0D+00/
C
      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      XINTX=ZERO
      YINTY=ZERO
      ZINTZ=ZERO
      XINTXX=ZERO
      YINTYY=ZERO
      ZINTZZ=ZERO
      NPTS=(NI+NJ-2+2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 16 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)/T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (7,6,5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    6 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    7 GO TO (15,14,13,12,11,10,9,8),NJ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   12 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   13 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   14 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   15 CONTINUE
      XINT=XINT+PX
      YINT=YINT+PY
      ZINT=ZINT+PZ
      XINTX=XINTX+PX*(PTX-XC)
      YINTY=YINTY+PY*(PTY-YC)
      ZINTZ=ZINTZ+PZ*(PTZ-ZC)
      XINTXX=XINTXX+PX*(PTX-XC)*(PTX-XC)
      YINTYY=YINTYY+PY*(PTY-YC)*(PTY-YC)
      ZINTZZ=ZINTZZ+PZ*(PTZ-ZC)*(PTZ-ZC)
   16 CONTINUE
      RETURN
      END
      SUBROUTINE HND_QDPAMS(BMASS,NCALL,NCODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "nwc_const.fh"
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL SAME
      CHARACTER*8 ERRMSG
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      DIMENSION AMASS(MXATOM),BMASS(1)
      DIMENSION ERRMSG(3)
      DIMENSION AMS(106)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-QDPAMS-'/
      DATA (AMS(I),I=1,54)  /
     1   1.007825D+00,  4.002600D+00,  7.016000D+00,  9.012180D+00,
     2  11.009310D+00, 12.000000D+00, 14.003070D+00, 15.994910D+00,
     3  18.998400D+00, 19.992440D+00, 22.989800D+00, 23.985040D+00,
     4  26.981530D+00, 27.976930D+00, 30.973760D+00, 31.972070D+00,
     5  34.968850D+00, 39.948000D+00, 38.963710D+00, 39.962590D+00,
     6  44.955920D+00, 47.900000D+00, 50.944000D+00, 51.940500D+00,
     7  54.938100D+00, 55.934900D+00, 58.933200D+00, 57.935300D+00,
     8  62.929800D+00, 63.929100D+00, 68.925700D+00, 73.921900D+00,
     9  74.921600D+00, 79.916500D+00, 78.918300D+00, 83.911500D+00,
     1  84.911700D+00, 87.905600D+00, 89.905400D+00, 89.904300D+00,
     2  92.906000D+00, 97.905500D+00, 97.000000D+00,101.903700D+00,
     3 102.904800D+00,105.903200D+00,106.904100D+00,113.903600D+00,
     4 114.904100D+00,119.902200D+00,120.903800D+00,129.906700D+00,
     5 126.904400D+00,131.904200D+00/
      DATA (AMS(I),I=55,106)  /
     1   132.9054D+00,137.9052D+00,138.9063D+00,139.9054D+00,
     2   140.9076D+00,141.9077D+00,144.9127D+00,151.9197D+00,
     3   152.9212D+00,157.9241D+00,158.9253D+00,163.9292D+00,
     4   164.9303D+00,165.9303D+00,168.9342D+00,173.9389D+00,
     5   174.9408D+00,179.9465D+00,180.9480D+00,183.9509D+00,
     6   186.9557D+00,191.9615D+00,192.9629D+00,194.9648D+00,
     7   196.9665D+00,201.9706D+00,204.9744D+00,207.9766D+00,
     8   208.9804D+00,208.9824D+00,209.9871D+00,222.0176D+00,
     9   223.0197D+00,226.0254D+00,227.0278D+00,232.0381D+00,
     1   231.0359D+00,238.0508D+00,237.0482D+00,244.0642D+00,
     2   243.0614D+00,247.0703D+00,247.0703D+00,251.0796D+00,
     3   252.0829D+00,257.0751D+00,258.0986D+00,259.1009D+00,
     4   260.1053D+00,261.1087D+00,  0.0000D+00,  0.0000D+00/
      DATA ZERO /0.0D+00/
C
      NAMELIST /MASS/ AMASS
C
      DO 10 IAT=1,NAT
         NUCZ=NUC(IAT)
         IF(NUCZ.GT.104) THEN
            WRITE(IW,9996)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(NUCZ.EQ.0) ZNUC=ZERO
         IF(NUCZ.GT.0) ZNUC=AMS(NUCZ)
   10    AMASS(IAT)=ZNUC
      IF(NCALL.EQ.1) GO TO 30
      IF(NCALL.GT.2)  GO TO 60
      IF(NCALL.EQ.2) REWIND IR
C%%   READ(IR,MASS,END=60,ERR=60)
C
      SAME=.TRUE.
      DO 20 IAT=1,NAT
      NUCZ=NUC(IAT)
      IF(NUCZ.EQ.0) ZNUC=ZERO
      IF(NUCZ.GT.0) ZNUC=AMS(NUCZ)
      IF(AMASS(IAT).EQ.ZNUC) GO TO 20
      SAME=.FALSE.
      GO TO 25
   20 CONTINUE
   25 CONTINUE
      IF(SAME) GO TO 60
C
   30 CONTINUE
      DO 40 IAT=1,NAT
   40 BMASS(IAT)= ABS(AMASS(IAT))
      if (ga_nodeid().eq.0) then
          WRITE(IW,9998)
          DO 50 IAT=1,NAT
   50     WRITE(IW,9997) IAT,ANAM(IAT),BNAM(IAT),BMASS(IAT)
      endif
      NCODE=0
      RETURN
   60 NCODE=1
      RETURN
 9998 FORMAT(/,10X,21(1H-),/,10X,'Atomic weights (a.u.)',
     1       /,10X,21(1H-),/)
 9997 FORMAT(I5,5X,A8,A2,F15.5)
 9996 FORMAT(' No atomic mass tabulated for  -nucz.gt.104-  . stop')
      END
C
      SUBROUTINE HND_OTPOLE(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      parameter  (mxefc=3*7*7*7)
      integer rtdb,basis,geom
      logical status
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2
      LOGICAL OUT,OTPOUT
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION AM(MXATOM),CM(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-OTPOLE'/
      DATA SCF,UHF /'SCF     ','UHF     '/
      DATA MP2     /'MP2     '/
      DATA OCTO    /0.711688D+00/
      DATA TOL     /1.0D-06/
      DATA ZERO    /0.0D+00/
      DATA TWO     /2.0D+00/
      DATA THREE   /3.0D+00/
      DATA FOUR    /4.0D+00/
      DATA FIVE    /5.0D+00/
      DATA KEYEFC  /4H EFC/
C
      OTPOUT=.TRUE. 
      IF(.NOT.OTPOUT) RETURN
      OUT=.FALSE.       
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM*NUM
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      NEED = L2+L2+L2+L2+L2+L2+L2+L2+L2+L2+L2+L2+L3+L3
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED , 'MEM INIT:prp:hnd_otpole:3',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_OTPOLE, MALLOC OF INIT  FAILED',911)

      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
C     ----- SET POINTERS FOR PARTITIONING OF CORE -----
C
      I10=  INIT
      I20=I10+L2
      I30=I20+L2
      I40=I30+L2
      I50=I40+L2
      I60=I50+L2
      I70=I60+L2
      I80=I70+L2
      I90=I80+L2
      I100=I90+L2
      I110=I100+L2
C     ----- MEMORY BLOCK FOR DENSITY MATRIX -----
      I120=I110+L2
      I130=I120+L2
C     ----- MEMORY BLOCK FOR DA IN NWCHEM FORMAT -----
      I140=I130+L3
C     ----- MOMORY BLOCK FOR DA IN HONDO  FORMAT-----
      I150=I140+L3
C
      NCALL=1
   10 CONTINUE
      CALL HND_OTPAMS(AM,NCALL,NCODE)
      IF(NCODE.NE.0) GO TO 100
C
      call prp_cent(rtdb,geom,cm,am,c,nat)
C
C     ----- CALCULATE THIRD MOMENT INTEGRALS -----
C
      CALL HND_OTPINT(X(I10),X(I20),X(I30),X(I40),X(I50),X(I60),
     1                X(I70),X(I80),X(I90),X(I100),
     2                CM,L1)
C
C     ----- ELECTRONIC CONTRIBUTION TO THIRD MOMENTS -----
C
C     ----- DENSITY MATRIX FOR X(I110)=DA,X(I120)=DB(0?) -----
C*    CALL hnd_prpDEN(X(I110),X(I120),L2)

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom, 
     1                      x(i110),x(i120),x(i130),x(i140),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
C
      RXXX=-hnd_prpDOT(X(I110),X(I10),L1)
      RYYY=-hnd_prpDOT(X(I110),X(I20),L1)
      RZZZ=-hnd_prpDOT(X(I110),X(I30),L1)
      RXXY=-hnd_prpDOT(X(I110),X(I40),L1)
      RXXZ=-hnd_prpDOT(X(I110),X(I50),L1)
      RYYX=-hnd_prpDOT(X(I110),X(I60),L1)
      RYYZ=-hnd_prpDOT(X(I110),X(I70),L1)
      RZZX=-hnd_prpDOT(X(I110),X(I80),L1)
      RZZY=-hnd_prpDOT(X(I110),X(I90),L1)
      RXYZ=-hnd_prpDOT(X(I110),X(I100),L1)
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) GO TO 50
      RXXX=-hnd_prpDOT(X(I120),X(I10),L1)+RXXX
      RYYY=-hnd_prpDOT(X(I120),X(I20),L1)+RYYY
      RZZZ=-hnd_prpDOT(X(I120),X(I30),L1)+RZZZ
      RXXY=-hnd_prpDOT(X(I120),X(I40),L1)+RXXY
      RXXZ=-hnd_prpDOT(X(I120),X(I50),L1)+RXXZ
      RYYX=-hnd_prpDOT(X(I120),X(I60),L1)+RYYX
      RYYZ=-hnd_prpDOT(X(I120),X(I70),L1)+RYYZ
      RZZX=-hnd_prpDOT(X(I120),X(I80),L1)+RZZX
      RZZY=-hnd_prpDOT(X(I120),X(I90),L1)+RZZY
      RXYZ=-hnd_prpDOT(X(I120),X(I100),L1)+RXYZ
   50 CONTINUE
C
C     ----- NUCLEAR CONTRIBUTION -----
C
      DO 60 IAT=1,NAT
      RXXX=RXXX+ZAN(IAT)*(C(1,IAT)-CM(1))*(C(1,IAT)-CM(1))
     1                                   *(C(1,IAT)-CM(1))
      RYYY=RYYY+ZAN(IAT)*(C(2,IAT)-CM(2))*(C(2,IAT)-CM(2))
     1                                   *(C(2,IAT)-CM(2))
      RZZZ=RZZZ+ZAN(IAT)*(C(3,IAT)-CM(3))*(C(3,IAT)-CM(3))
     1                                   *(C(3,IAT)-CM(3))
      RXXY=RXXY+ZAN(IAT)*(C(1,IAT)-CM(1))*(C(1,IAT)-CM(1))
     1                                   *(C(2,IAT)-CM(2))
      RXXZ=RXXZ+ZAN(IAT)*(C(1,IAT)-CM(1))*(C(1,IAT)-CM(1))
     1                                   *(C(3,IAT)-CM(3))
      RYYX=RYYX+ZAN(IAT)*(C(2,IAT)-CM(2))*(C(2,IAT)-CM(2))
     1                                   *(C(1,IAT)-CM(1))
      RYYZ=RYYZ+ZAN(IAT)*(C(2,IAT)-CM(2))*(C(2,IAT)-CM(2))
     1                                   *(C(3,IAT)-CM(3))
      RZZX=RZZX+ZAN(IAT)*(C(3,IAT)-CM(3))*(C(3,IAT)-CM(3))
     1                                   *(C(1,IAT)-CM(1))
      RZZY=RZZY+ZAN(IAT)*(C(3,IAT)-CM(3))*(C(3,IAT)-CM(3))
     1                                   *(C(2,IAT)-CM(2))
      RXYZ=RXYZ+ZAN(IAT)*(C(1,IAT)-CM(1))*(C(2,IAT)-CM(2))
     1                                   *(C(3,IAT)-CM(3))
   60 CONTINUE
C
C     ----- -EFC- CONTRIBUTION -----
C
      RXXXEF=ZERO
      RYYYEF=ZERO
      RZZZEF=ZERO
      RXXYEF=ZERO
      RXXZEF=ZERO
      RYYXEF=ZERO
      RYYZEF=ZERO
      RZZXEF=ZERO
      RZZYEF=ZERO
      RXYZEF=ZERO
      IF(IEFC.NE.KEYEFC) GO TO 80
      DO 70 IAT=1,NEFC
      RXXXEF=RXXXEF+EFCZ(IAT)*(EFCC(1,IAT)-CM(1))*(EFCC(1,IAT)-CM(1))
     1                                           *(EFCC(1,IAT)-CM(1))
      RYYYEF=RYYYEF+EFCZ(IAT)*(EFCC(2,IAT)-CM(2))*(EFCC(2,IAT)-CM(2))
     1                                           *(EFCC(2,IAT)-CM(2))
      RZZZEF=RZZZEF+EFCZ(IAT)*(EFCC(3,IAT)-CM(3))*(EFCC(3,IAT)-CM(3))
     1                                           *(EFCC(3,IAT)-CM(3))
      RXXYEF=RXXYEF+EFCZ(IAT)*(EFCC(1,IAT)-CM(1))*(EFCC(1,IAT)-CM(1))
     1                                           *(EFCC(2,IAT)-CM(2))
      RXXZEF=RXXZEF+EFCZ(IAT)*(EFCC(1,IAT)-CM(1))*(EFCC(1,IAT)-CM(1))
     1                                           *(EFCC(3,IAT)-CM(3))
      RYYXEF=RYYXEF+EFCZ(IAT)*(EFCC(2,IAT)-CM(2))*(EFCC(2,IAT)-CM(2))
     1                                           *(EFCC(1,IAT)-CM(1))
      RYYZEF=RYYZEF+EFCZ(IAT)*(EFCC(2,IAT)-CM(2))*(EFCC(2,IAT)-CM(2))
     1                                           *(EFCC(3,IAT)-CM(3))
      RZZXEF=RZZXEF+EFCZ(IAT)*(EFCC(3,IAT)-CM(3))*(EFCC(3,IAT)-CM(3))
     1                                           *(EFCC(1,IAT)-CM(1))
      RZZYEF=RZZYEF+EFCZ(IAT)*(EFCC(3,IAT)-CM(3))*(EFCC(3,IAT)-CM(3))
     1                                           *(EFCC(2,IAT)-CM(2))
      RXYZEF=RXYZEF+EFCZ(IAT)*(EFCC(1,IAT)-CM(1))*(EFCC(2,IAT)-CM(2))
     1                                           *(EFCC(3,IAT)-CM(3))
   70 CONTINUE
   80 CONTINUE
C
      RXXXT =RXXX+RXXXEF
      RYYYT =RYYY+RYYYEF
      RZZZT =RZZZ+RZZZEF
      RXXYT =RXXY+RXXYEF
      RXXZT =RXXZ+RXXZEF
      RYYXT =RYYX+RYYXEF
      RYYZT =RYYZ+RYYZEF
      RZZXT =RZZX+RZZXEF
      RZZYT =RZZY+RZZYEF
      RXYZT =RXYZ+RXYZEF
      if (ga_nodeid().eq.0) then
          WRITE(IW,9999) RXXX,RXXXEF,RXXXT,RYYY,RYYYEF,RYYYT,
     1               RZZZ,RZZZEF,RZZZT,RXXY,RXXYEF,RXXYT,
     2               RXXZ,RXXZEF,RXXZT,RYYX,RYYXEF,RYYXT,
     3               RYYZ,RYYZEF,RYYZT,RZZX,RZZXEF,RZZXT,
     4               RZZY,RZZYEF,RZZYT,RXYZ,RXYZEF,RXYZT
          WRITE(IW,9995) RXXX*OCTO,RXXXEF*OCTO,RXXXT*OCTO,
     1               RYYY*OCTO,RYYYEF*OCTO,RYYYT*OCTO,
     2               RZZZ*OCTO,RZZZEF*OCTO,RZZZT*OCTO,
     3               RXXY*OCTO,RXXYEF*OCTO,RXXYT*OCTO,
     4               RXXZ*OCTO,RXXZEF*OCTO,RXXZT*OCTO,
     5               RYYX*OCTO,RYYXEF*OCTO,RYYXT*OCTO,
     6               RYYZ*OCTO,RYYZEF*OCTO,RYYZT*OCTO,
     7               RZZX*OCTO,RZZXEF*OCTO,RZZXT*OCTO,
     8               RZZY*OCTO,RZZYEF*OCTO,RZZYT*OCTO,
     9               RXYZ*OCTO,RXYZEF*OCTO,RXYZT*OCTO
      endif
      OXXX=      RXXX-THREE*(RYYX+RZZX)/TWO
      OYYY=      RYYY-THREE*(RXXY+RZZY)/TWO
      OZZZ=      RZZZ-THREE*(RXXZ+RYYZ)/TWO
      OXXY=(FOUR*RXXY-RYYY-RZZY)/TWO
      OXXZ=(FOUR*RXXZ-RYYZ-RZZZ)/TWO
      OYYX=(FOUR*RYYX-RXXX-RZZX)/TWO
      OYYZ=(FOUR*RYYZ-RXXZ-RZZZ)/TWO
      OZZX=(FOUR*RZZX-RXXX-RYYX)/TWO
      OZZY=(FOUR*RZZY-RXXY-RYYY)/TWO
      OXYZ= FIVE*RXYZ/TWO
      OXXXEF=      RXXXEF-THREE*(RYYXEF+RZZXEF)/TWO
      OYYYEF=      RYYYEF-THREE*(RXXYEF+RZZYEF)/TWO
      OZZZEF=      RZZZEF-THREE*(RXXZEF+RYYZEF)/TWO
      OXXYEF=(FOUR*RXXYEF-RYYYEF-RZZYEF)/TWO
      OXXZEF=(FOUR*RXXZEF-RYYZEF-RZZZEF)/TWO
      OYYXEF=(FOUR*RYYXEF-RXXXEF-RZZXEF)/TWO
      OYYZEF=(FOUR*RYYZEF-RXXZEF-RZZZEF)/TWO
      OZZXEF=(FOUR*RZZXEF-RXXXEF-RYYXEF)/TWO
      OZZYEF=(FOUR*RZZYEF-RXXYEF-RYYYEF)/TWO
      OXYZEF= FIVE*RXYZEF/TWO
      OXXXT=      RXXXT-THREE*(RYYXT+RZZXT)/TWO
      OYYYT=      RYYYT-THREE*(RXXYT+RZZYT)/TWO
      OZZZT=      RZZZT-THREE*(RXXZT+RYYZT)/TWO
      OXXYT=(FOUR*RXXYT-RYYYT-RZZYT)/TWO
      OXXZT=(FOUR*RXXZT-RYYZT-RZZZT)/TWO
      OYYXT=(FOUR*RYYXT-RXXXT-RZZXT)/TWO
      OYYZT=(FOUR*RYYZT-RXXZT-RZZZT)/TWO
      OZZXT=(FOUR*RZZXT-RXXXT-RYYXT)/TWO
      OZZYT=(FOUR*RZZYT-RXXYT-RYYYT)/TWO
      OXYZT= FIVE*RXYZT/TWO
      if (ga_nodeid().eq.0) then
          WRITE(IW,9997) OXXX,OXXXEF,OXXXT,OYYY,OYYYEF,OYYYT,
     1               OZZZ,OZZZEF,OZZZT,OXXY,OXXYEF,OXXYT,
     2               OXXZ,OXXZEF,OXXZT,OYYX,OYYXEF,OYYXT,
     3               OYYZ,OYYZEF,OYYZT,OZZX,OZZXEF,OZZXT,
     4               OZZY,OZZYEF,OZZYT,OXYZ,OXYZEF,OXYZT
          WRITE(IW,9994) OXXX*OCTO,OXXXEF*OCTO,OXXXT*OCTO,
     1               OYYY*OCTO,OYYYEF*OCTO,OYYYT*OCTO,
     2               OZZZ*OCTO,OZZZEF*OCTO,OZZZT*OCTO,
     3               OXXY*OCTO,OXXYEF*OCTO,OXXYT*OCTO,
     4               OXXZ*OCTO,OXXZEF*OCTO,OXXZT*OCTO,
     5               OYYX*OCTO,OYYXEF*OCTO,OYYXT*OCTO,
     6               OYYZ*OCTO,OYYZEF*OCTO,OYYZT*OCTO,
     7               RZZX*OCTO,RZZXEF*OCTO,RZZXT*OCTO,
     8               RZZY*OCTO,RZZYEF*OCTO,RZZYT*OCTO,
     9               RXYZ*OCTO,RXYZEF*OCTO,RXYZT*OCTO
      endif
C
      NCALL=NCALL+1
C     GO TO 10
C
  100 CONTINUE
C
      if (ga_nodeid().eq.0) WRITE(IW,9993)
C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_OTPOLE, MA_POP_STACK OF INIT
     & FAILED',911)
C
      RETURN
C
 9999 FORMAT(/,2X,' Third moments in atomic units',
     2 /,2X,' RXXX',F20.10,' RXXXEF',F20.10,' RXXXT',F20.10,
     2 /,2X,' RYYY',F20.10,' RYYYEF',F20.10,' RYYYT',F20.10,
     2 /,2X,' RZZZ',F20.10,' RZZZEF',F20.10,' RZZZT',F20.10,
     2 /,2X,' RXXY',F20.10,' RXXYEF',F20.10,' RXXYT',F20.10,
     2 /,2X,' RXXZ',F20.10,' RXXZEF',F20.10,' RXXZT',F20.10,
     2 /,2X,' RYYX',F20.10,' RYYXEF',F20.10,' RYYXT',F20.10,
     2 /,2X,' RYYZ',F20.10,' RYYZEF',F20.10,' RYYZT',F20.10,
     2 /,2X,' RZZX',F20.10,' RZZXEF',F20.10,' RZZXT',F20.10,
     2 /,2X,' RZZY',F20.10,' RZZYEF',F20.10,' RZZYT',F20.10,
     2 /,2X,' RXYZ',F20.10,' RXYZEF',F20.10,' RXYZT',F20.10)
 9997 FORMAT(/,2X,' Octupole moments in atomic units',
     2 /,2X,' OXXX',F20.10,' OXXXEF',F20.10,' OXXXT',F20.10,
     2 /,2X,' OYYY',F20.10,' OYYYEF',F20.10,' OYYYT',F20.10,
     2 /,2X,' OZZZ',F20.10,' OZZZEF',F20.10,' OZZZT',F20.10,
     2 /,2X,' OXXY',F20.10,' OXXYEF',F20.10,' OXXYT',F20.10,
     2 /,2X,' OXXZ',F20.10,' OXXZEF',F20.10,' OXXZT',F20.10,
     2 /,2X,' OYYX',F20.10,' OYYXEF',F20.10,' OYYXT',F20.10,
     2 /,2X,' OYYZ',F20.10,' OYYZEF',F20.10,' OYYZT',F20.10,
     2 /,2X,' OZZX',F20.10,' OZZXEF',F20.10,' OZZXT',F20.10,
     2 /,2X,' OZZY',F20.10,' OZZYEF',F20.10,' OZZYT',F20.10,
     2 /,2X,' OXYZ',F20.10,' OXYZEF',F20.10,' OXYZT',F20.10)
 9996 FORMAT(' Not enough core in -otpole- . stop ')
 9995 FORMAT(/,2X,' Third moments in 10**(-34) esu*cm**3',
     2 /,2X,' RXXX',F20.10,' RXXXEF',F20.10,' RXXXT',F20.10,
     2 /,2X,' RYYY',F20.10,' RYYYEF',F20.10,' RYYYT',F20.10,
     2 /,2X,' RZZZ',F20.10,' RZZZEF',F20.10,' RZZZT',F20.10,
     2 /,2X,' RXXY',F20.10,' RXXYEF',F20.10,' RXXYT',F20.10,
     2 /,2X,' RXXZ',F20.10,' RXXZEF',F20.10,' RXXZT',F20.10,
     2 /,2X,' RYYX',F20.10,' RYYXEF',F20.10,' RYYXT',F20.10,
     2 /,2X,' RYYZ',F20.10,' RYYZEF',F20.10,' RYYZT',F20.10,
     2 /,2X,' RZZX',F20.10,' RZZXEF',F20.10,' RZZXT',F20.10,
     2 /,2X,' RZZY',F20.10,' RZZYEF',F20.10,' RZZYT',F20.10,
     2 /,2X,' RXYZ',F20.10,' RXYZEF',F20.10,' RXYZT',F20.10)
 9994 FORMAT(/,2X,' Octupole moments in 10**(-34) esu*cm**3',
     2 /,2X,' OXXX',F20.10,' OXXXEF',F20.10,' OXXXT',F20.10,
     2 /,2X,' OYYY',F20.10,' OYYYEF',F20.10,' OYYYT',F20.10,
     2 /,2X,' OZZZ',F20.10,' OZZZEF',F20.10,' OZZZT',F20.10,
     2 /,2X,' OXXY',F20.10,' OXXYEF',F20.10,' OXXYT',F20.10,
     2 /,2X,' OXXZ',F20.10,' OXXZEF',F20.10,' OXXZT',F20.10,
     2 /,2X,' OYYX',F20.10,' OYYXEF',F20.10,' OYYXT',F20.10,
     2 /,2X,' OYYZ',F20.10,' OYYZEF',F20.10,' OYYZT',F20.10,
     2 /,2X,' OZZX',F20.10,' OZZXEF',F20.10,' OZZXT',F20.10,
     2 /,2X,' OZZY',F20.10,' OZZYEF',F20.10,' OZZYT',F20.10,
     2 /,2X,' OXYZ',F20.10,' OXYZEF',F20.10,' OXYZT',F20.10)
 9993 FORMAT(/,' 1 a.u. = 0.711688 10**(-34) esu*cm**3 ')
      END
      SUBROUTINE HND_OTPINT(XXXS,YYYS,ZZZS,XXYS,XXZS,YYXS,YYZS,
     1                  ZZXS,ZZYS,XYZS,CM,L1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZOTP/XINT,YINT,ZINT,XINTX,YINTY,ZINTZ,
     1              XINTXX,YINTYY,ZINTZZ,XINXXX,YINYYY,ZINZZZ,
     2              T,XC,YC,ZC,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      DIMENSION XXXS(1),YYYS(1),ZZZS(1),XXYS(1),XXZS(1),YYXS(1),YYZS(1)
      DIMENSION ZZXS(1),ZZYS(1),XYZS(1)
      DIMENSION SXXX(225),SYYY(225),SZZZ(225)
      DIMENSION SXXY(225),SXXZ(225),SYYX(225),SYYZ(225)
      DIMENSION SZZX(225),SZZY(225),SXYZ(225)
      DIMENSION DIJ(225)
      DIMENSION    XIN(5,5),   YIN(5,5),   ZIN(5,5)
      DIMENSION   XXIN(5,5),  YYIN(5,5),  ZZIN(5,5)
      DIMENSION  XXXIN(5,5), YYYIN(5,5), ZZZIN(5,5)
      DIMENSION XXXXIN(5,5),YYYYIN(5,5),ZZZZIN(5,5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION CM(3)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      TOL=RLN10*ITOL
      OUT=NPRINT.EQ.6
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      XC=CM(1)
      YC=CM(2)
      ZC=CM(3)
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS -----
C
      IJ=0
      DO 100 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 100 J=MINJ,JMAX
      IJ=IJ+1
      SXXX(IJ)=ZERO
      SYYY(IJ)=ZERO
      SZZZ(IJ)=ZERO
      SXXY(IJ)=ZERO
      SXXZ(IJ)=ZERO
      SYYX(IJ)=ZERO
      SYYZ(IJ)=ZERO
      SZZX(IJ)=ZERO
      SZZY(IJ)=ZERO
      SXYZ(IJ)=ZERO
  100 CONTINUE
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 DIJ(IJ)=DUM2
C
C     ----- OCTUPOLE MOMENT INTEGRALS -----
C
      T = SQRT(AA)
      T1=ONE/T
      X0=AX
      Y0=AY
      Z0=AZ
      DO 370 J=1,LJT
      NJ=J
      DO 370 I=1,LIT
      NI=I
      CALL HND_OTPXYZ
         XIN(I,J)=XINT*T1
         YIN(I,J)=YINT*T1
         ZIN(I,J)=ZINT*T1
        XXIN(I,J)=XINTX*T1
        YYIN(I,J)=YINTY*T1
        ZZIN(I,J)=ZINTZ*T1
       XXXIN(I,J)=XINTXX*T1
       YYYIN(I,J)=YINTYY*T1
       ZZZIN(I,J)=ZINTZZ*T1
      XXXXIN(I,J)=XINXXX*T1
      YYYYIN(I,J)=YINYYY*T1
      ZZZZIN(I,J)=ZINZZZ*T1
  370 CONTINUE
C
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 390 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      IJ=IJ+1
      SXXX(IJ)=SXXX(IJ)+
     1         DIJ(IJ)*(XXXXIN(IX,JX)*   YIN(IY,JY)*   ZIN(IZ,JZ))
      SYYY(IJ)=SYYY(IJ)+
     1         DIJ(IJ)*(   XIN(IX,JX)*YYYYIN(IY,JY)*   ZIN(IZ,JZ))
      SZZZ(IJ)=SZZZ(IJ)+
     1         DIJ(IJ)*(   XIN(IX,JX)*   YIN(IY,JY)*ZZZZIN(IZ,JZ))
      SXXY(IJ)=SXXY(IJ)+
     1         DIJ(IJ)*( XXXIN(IX,JX)*  YYIN(IY,JY)*   ZIN(IZ,JZ))
      SXXZ(IJ)=SXXZ(IJ)+
     1         DIJ(IJ)*( XXXIN(IX,JX)*   YIN(IY,JY)*  ZZIN(IZ,JZ))
      SYYX(IJ)=SYYX(IJ)+
     1         DIJ(IJ)*(  XXIN(IX,JX)* YYYIN(IY,JY)*   ZIN(IZ,JZ))
      SYYZ(IJ)=SYYZ(IJ)+
     1         DIJ(IJ)*(   XIN(IX,JX)* YYYIN(IY,JY)*  ZZIN(IZ,JZ))
      SZZX(IJ)=SZZX(IJ)+
     1         DIJ(IJ)*(  XXIN(IX,JX)*   YIN(IY,JY)* ZZZIN(IZ,JZ))
      SZZY(IJ)=SZZY(IJ)+
     1         DIJ(IJ)*(   XIN(IX,JX)*  YYIN(IY,JY)* ZZZIN(IZ,JZ))
      SXYZ(IJ)=SXYZ(IJ)+
     1         DIJ(IJ)*(  XXIN(IX,JX)*  YYIN(IY,JY)*  ZZIN(IZ,JZ))
  390 CONTINUE
 6000 CONTINUE
 7000 CONTINUE
C
C     ----- SET UP OCTUPOLE MOMENT MATRICES -----
C
      IJ=0
      DO 7500 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 7500 J=MINJ,JMAX
      IJ=IJ+1
      LOCII=LOCI+I
      LOCJJ=LOCJ+J
      NN=(LOCII*(LOCII-1))/2+LOCJJ
      XXXS(NN)=SXXX(IJ)
      YYYS(NN)=SYYY(IJ)
      ZZZS(NN)=SZZZ(IJ)
      XXYS(NN)=SXXY(IJ)
      XXZS(NN)=SXXZ(IJ)
      YYXS(NN)=SYYX(IJ)
      YYZS(NN)=SYYZ(IJ)
      ZZXS(NN)=SZZX(IJ)
      ZZYS(NN)=SZZY(IJ)
      XYZS(NN)=SXYZ(IJ)
 7500 CONTINUE
C
 8000 CONTINUE
 9000 CONTINUE
C
      IF(.NOT.OUT) GO TO 9100
      WRITE(IW,9999)
      CALL HND_PRTRL(XXXS,L1)
      WRITE(IW,9998)
      CALL HND_PRTRL(YYYS,L1)
      WRITE(IW,9997)
      CALL HND_PRTRL(ZZZS,L1)
      WRITE(IW,9996)
      CALL HND_PRTRL(XXYS,L1)
      WRITE(IW,9995)
      CALL HND_PRTRL(XXZS,L1)
      WRITE(IW,9994)
      CALL HND_PRTRL(YYXS,L1)
      WRITE(IW,9993)
      CALL HND_PRTRL(YYZS,L1)
      WRITE(IW,9992)
      CALL HND_PRTRL(ZZXS,L1)
      WRITE(IW,9991)
      CALL HND_PRTRL(ZZYS,L1)
      WRITE(IW,9990)
      CALL HND_PRTRL(XYZS,L1)
 9100 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,26(1H-),/,10X,'XXX-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9998 FORMAT(/,10X,26(1H-),/,10X,'YYY-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9997 FORMAT(/,10X,26(1H-),/,10X,'ZZZ-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9996 FORMAT(/,10X,26(1H-),/,10X,'XXY-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9995 FORMAT(/,10X,26(1H-),/,10X,'XXZ-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9994 FORMAT(/,10X,26(1H-),/,10X,'YYX-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9993 FORMAT(/,10X,26(1H-),/,10X,'YYZ-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9992 FORMAT(/,10X,26(1H-),/,10X,'ZZX-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9991 FORMAT(/,10X,26(1H-),/,10X,'ZZY-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
 9990 FORMAT(/,10X,26(1H-),/,10X,'XYZ-THIRD MOMENT INTEGRALS',/,
     1         10X,26(1H-))
      END
      SUBROUTINE HND_OTPXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C
#include "hnd_whermt.fh"
c
      COMMON/HND_XYZOTP/XINT,YINT,ZINT,XINTX,YINTY,ZINTZ,
     1              XINTXX,YINTYY,ZINTZZ,XINXXX,YINYYY,ZINZZZ,
     2              T,XC,YC,ZC,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      DIMENSION MIN(7),MAX(7)
      DATA MIN /1,2,4, 7,11,16,22/
      DATA MAX /1,3,6,10,15,21,28/
      DATA ZERO /0.0D+00/
C
      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      XINTX=ZERO
      YINTY=ZERO
      ZINTZ=ZERO
      XINTXX=ZERO
      YINTYY=ZERO
      ZINTZZ=ZERO
      XINXXX=ZERO
      YINYYY=ZERO
      ZINZZZ=ZERO
      NPTS=(NI+NJ-2+3)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 16 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)/T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (7,6,5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    6 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    7 GO TO (15,14,13,12,11,10,9,8),NJ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   12 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   13 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   14 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   15 CONTINUE
      XINT=XINT+PX
      YINT=YINT+PY
      ZINT=ZINT+PZ
      XINTX=XINTX+PX*(PTX-XC)
      YINTY=YINTY+PY*(PTY-YC)
      ZINTZ=ZINTZ+PZ*(PTZ-ZC)
      XINTXX=XINTXX+PX*(PTX-XC)*(PTX-XC)
      YINTYY=YINTYY+PY*(PTY-YC)*(PTY-YC)
      ZINTZZ=ZINTZZ+PZ*(PTZ-ZC)*(PTZ-ZC)
      XINXXX=XINXXX+PX*(PTX-XC)*(PTX-XC)*(PTX-XC)
      YINYYY=YINYYY+PY*(PTY-YC)*(PTY-YC)*(PTY-YC)
      ZINZZZ=ZINZZZ+PZ*(PTZ-ZC)*(PTZ-ZC)*(PTZ-ZC)
   16 CONTINUE
      RETURN
      END
      SUBROUTINE HND_OTPAMS(BMASS,NCALL,NCODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "nwc_const.fh"
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL SAME
      CHARACTER*8 ERRMSG
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      DIMENSION AMASS(MXATOM),BMASS(1)
      DIMENSION ERRMSG(3)
      DIMENSION AMS(106)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-OTPAMS-'/
      DATA (AMS(I),I=1,54)  /
     1   1.007825D+00,  4.002600D+00,  7.016000D+00,  9.012180D+00,
     2  11.009310D+00, 12.000000D+00, 14.003070D+00, 15.994910D+00,
     3  18.998400D+00, 19.992440D+00, 22.989800D+00, 23.985040D+00,
     4  26.981530D+00, 27.976930D+00, 30.973760D+00, 31.972070D+00,
     5  34.968850D+00, 39.948000D+00, 38.963710D+00, 39.962590D+00,
     6  44.955920D+00, 47.900000D+00, 50.944000D+00, 51.940500D+00,
     7  54.938100D+00, 55.934900D+00, 58.933200D+00, 57.935300D+00,
     8  62.929800D+00, 63.929100D+00, 68.925700D+00, 73.921900D+00,
     9  74.921600D+00, 79.916500D+00, 78.918300D+00, 83.911500D+00,
     1  84.911700D+00, 87.905600D+00, 89.905400D+00, 89.904300D+00,
     2  92.906000D+00, 97.905500D+00, 97.000000D+00,101.903700D+00,
     3 102.904800D+00,105.903200D+00,106.904100D+00,113.903600D+00,
     4 114.904100D+00,119.902200D+00,120.903800D+00,129.906700D+00,
     5 126.904400D+00,131.904200D+00/
      DATA (AMS(I),I=55,106)  /
     1   132.9054D+00,137.9052D+00,138.9063D+00,139.9054D+00,
     2   140.9076D+00,141.9077D+00,144.9127D+00,151.9197D+00,
     3   152.9212D+00,157.9241D+00,158.9253D+00,163.9292D+00,
     4   164.9303D+00,165.9303D+00,168.9342D+00,173.9389D+00,
     5   174.9408D+00,179.9465D+00,180.9480D+00,183.9509D+00,
     6   186.9557D+00,191.9615D+00,192.9629D+00,194.9648D+00,
     7   196.9665D+00,201.9706D+00,204.9744D+00,207.9766D+00,
     8   208.9804D+00,208.9824D+00,209.9871D+00,222.0176D+00,
     9   223.0197D+00,226.0254D+00,227.0278D+00,232.0381D+00,
     1   231.0359D+00,238.0508D+00,237.0482D+00,244.0642D+00,
     2   243.0614D+00,247.0703D+00,247.0703D+00,251.0796D+00,
     3   252.0829D+00,257.0751D+00,258.0986D+00,259.1009D+00,
     4   260.1053D+00,261.1087D+00,  0.0000D+00,  0.0000D+00/
      DATA ZERO /0.0D+00/
C
      NAMELIST /MASS/ AMASS
C
      DO 10 IAT=1,NAT
         NUCZ=NUC(IAT)
         IF(NUCZ.GT.104) THEN
            WRITE(IW,9996)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(NUCZ.EQ.0) ZNUC=ZERO
         IF(NUCZ.GT.0) ZNUC=AMS(NUCZ)
   10    AMASS(IAT)=ZNUC
      IF(NCALL.EQ.1) GO TO 30
      IF(NCALL.GT.2)  GO TO 60
      IF(NCALL.EQ.2) REWIND IR
C%%   READ(IR,MASS,END=60,ERR=60)
C
      SAME=.TRUE.
      DO 20 IAT=1,NAT
      NUCZ=NUC(IAT)
      IF(NUCZ.EQ.0) ZNUC=ZERO
      IF(NUCZ.GT.0) ZNUC=AMS(NUCZ)
      IF(AMASS(IAT).EQ.ZNUC) GO TO 20
      SAME=.FALSE.
      GO TO 25
   20 CONTINUE
   25 CONTINUE
      IF(SAME) GO TO 60
C
   30 CONTINUE
      DO 40 IAT=1,NAT
   40 BMASS(IAT)= ABS(AMASS(IAT))
      if (ga_nodeid().eq.0) then
         WRITE(IW,9998)
         DO 50 IAT=1,NAT
   50    WRITE(IW,9997) IAT,ANAM(IAT),BNAM(IAT),BMASS(IAT)
      endif
      NCODE=0
      RETURN
   60 NCODE=1
      RETURN
 9998 FORMAT(/,10X,21(1H-),/,10X,'Atomic weights (a.u.)',
     1       /,10X,21(1H-),/)
 9997 FORMAT(I5,5X,A8,A2,F15.5)
 9996 FORMAT(' No atomic mass tabulated for  -nucz.gt.104-  . stop')
      END
      SUBROUTINE HND_MULKEN(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      integer   rtdb, basis, geom
      logical   status
      integer   g_s
      integer   g_dens(3)
      integer   ga_create_atom_blocked
      external  ga_create_atom_blocked
      PARAMETER (MXIODA=255)
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP, SCFTYP, RUNTYP
      CHARACTER*8 SCF, UHF, MP2, PROP
      CHARACTER*8 ELEC, ALPHA, BETA, ALL
      character*8 anam,bflab
      character*2 bnam
      LOGICAL LAST
      LOGICAL OUT,POPOUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/hnd_prpOPT/NOPROP
      COMMON/HND_RUNOPT/RUNTYP
      COMMON/HND_BASSPH/ISPHER
      COMMON/HND_MOLSPH/SPHLAB(MXBFN)
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MLKPOP/CHGMLK(MXATOM)
      COMMON/HND_BFNATM/LIMLOW(MXATOM),LIMSUP(MXATOM)
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION ZNUC(MXATOM)
      DIMENSION ERRMSG(3)
      DATA ERRMSG         /'PROGRAM ','STOP IN ','-MULKEN-'/
      DATA ALPHA,BETA,ALL /'   ALPHA','    BETA','     ALL'/
      DATA SCF,UHF        /'SCF     ','UHF     '/
      DATA MP2            /'MP2     '/
      DATA PROP           /'PROPTY  '/
      DATA ZERO           /0.0D+00/
      DATA NORB           /0/
      DATA MOFLG          /0/
C
      POPOUT=.TRUE.          
      IF(.NOT.POPOUT) RETURN
      OUT=.FALSE.
C
C     ----- READ NAMELIST /MLK/ -----
C
*     REWIND IR
*     READ(IR,MLK,END=10,ERR=10)
*  10 CONTINUE
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM*NUM
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      NEED = L3+L1+L2  + 
     &       (NAT*(NAT+1))/2+L2+L1+NAT  +
     &       L2+L1+NAT +
     &       L3+L1     +
     &       L3+NAT*NAT+NAT+NAT+L3+NAT+NAT  +
     &       L3+L3+L3+L3
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED , 'MEM INIT:prp:hnd_mulken:4',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_MULKEN, MALLOC OF INIT  FAILED',911)
C
      status=MA_verify_allocator_stuff()
C
      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
C     ----- SET POINTERS FOR PARTITIONING OF CORE -----
C
      I00=  INIT
      I01=I00+L3
      I10=I01+L1
      I20=I10+L2
C
      I30=I20
      I40=I30+(NAT*(NAT+1))/2
      I50=I40+L2
      I60=I50+L1
      I70=I60+NAT
C
      I41=I70
      I51=I41+L2
      I61=I51+L1
      I71=I61+NAT
C
      I42=I71+L3
      I52=I42+L1
C
      I53=I52+L3
      I63=I53+NAT*NAT
      I73=I63+NAT
      I83=I73+NAT
      I93=I83+L3
      I94=I93+NAT
      I95=I94+NAT
C     ----- MEMORY BLOCK FOR DA IN NWCHEM  FORMAT -----
      I96=I95+L3
C     ----- MEMORY BLOCK FOR DA IN HONDO   FORMAT-----
      I97=I96+L3
C     ----- MEMORY BLOCK FOR S  IN NWCHEM  FORMAT-----
      I98=I97+L3
C     ----- MEMORY BLOCK FOR S  IN HONDO   FORMAT-----
      I99=I98+L3
C
      DO IAT=1,NAT
         ZNUC(IAT)=NUC(IAT)
      ENDDO
C
C     ----- DETERMINE NUMBER OF BASIS FUNCTIONS PER ATOM -----
C
      CALL HND_ATMBFN
C
C     -----         READ IN OVERLAP MATRIX              -----
C
C     -  S- AT X(I10)
C
c
c     ----- get overlap matrix at x(i97) from nwchem -----
c
      g_s = ga_create_atom_blocked(geom, basis,'int1e: s')
      call ga_zero(g_s)
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)
      call schwarz_tidy()
      call int_terminate()
      call ga_get(g_s,1,num,1,num,x(i97),num)

c     ----- transform S from nwchem to hnd format -----
 
      call hnd_nwhnd_tran_s(x(i97),x(i98),num)

c     ----- destroy ga handle of overlap matrix -----
      status = ga_destroy ( g_s )
c    
c     ----- transform S from square to triangular format -----
c
      ij = 0
      do i = 1,num
         do j = 1,i
            ij = ij + 1
            x(i10 +ij-1) = x(i98 + (j-1)*(num)+(i-1))
         enddo
      enddo
C     CALL HND_DAREAD(IDAF,IODA,X(I10),L2,12)
      IF(OUT.and.ga_nodeid().eq.0) THEN
         WRITE(IW,9989)
         CALL HND_PRTR(X(I10),L1)
      ENDIF
C
C     ----- CALCULATE POPULATION FOR INDIVIDUAL ORBITALS -----
C
C            -S - AT X(I10)
C
c--   CALL HND_ORBPOP(rtdb,geom,basis,
c--  1                X(I10),X(I40),X(I42),
c--  2                IA,NORB,L1,L1,MOFLG,X(I00))

C
      IF(RUNTYP.EQ.PROP) GO TO 1000
C
C     -----            READ DENSITY MATRICES           -----
C
C               -DA- AT X(I40)
C               -DB- AT X(I41)
C
C*       CALL hnd_prpDEN(X(I40),X(I41),L2)

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i40),x(i41),x(i95),x(i96),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
      if(scftyp.eq.uhf) then
         wfntyp=scf
      endif
      if(out) then
         write(6,*) 'wfntyp = ',wfntyp
         write(6,*) 'scftyp = ',scftyp
      endif
C
      IPASS=1
      IDEN=I40
      IGOC=I50
      IGAC=I60
      LAST=.FALSE.
      ELEC=ALPHA
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) ELEC=All
  100 CONTINUE
C
C     ----- DO A MULLIKEN POPULATION ANALYSIS -----
C              CALCULATE OVERLAP POPULATION
C
      CALL HND_OVLPOP(X(IDEN),X(I10),L2)
  200 if (ga_nodeid().eq.0) WRITE(IW,9999) ELEC
      IF(OUT.and.ga_nodeid().eq.0) THEN
         WRITE(IW,9988)
         CALL HND_PRTR (X(IDEN),L1)
      ENDIF
C
C     ----- CALCULATE TOTAL GROSS POPULATION IN AO'S ----
C
      CALL HND_GROSSC(X(IDEN),X(IGOC),IA,L1)
      if (ga_nodeid().eq.0) then
         WRITE(IW,9998)
         WRITE(IW,9997) (I, BFLAB(I),X(I-1+IGOC),I=1,L1)
      endif
C
      if (ga_nodeid().eq.0) WRITE(IW,9996)
C
C     ----- COMPRESS FROM -AO- TO SHELLS -----
C
      CALL HND_SHLPOP(X(IGOC))
C
C     ----- COMPRESS FROM -AO- TO ATOMS -----
C
      CALL HND_ATMPOP(X(IDEN),IA,X(I30),NAT)
      IF(OUT) THEN
         CALL HND_PRTR(X(I30),NAT)
      ENDIF
C
C     ----- CALCULATE TOTAL GROSS POPULATION ON ATOMS -----
C
      CALL HND_GROSSC(X(I30),X(IGAC),IA,NAT)
      if (ga_nodeid().eq.0) then
         WRITE(IW,9995)
         WRITE(IW,9992) (I,ANAM(I),BNAM(I),ZNUC(I),ZAN(I),
     1                X(I-1+IGAC),I=1,NAT)
      endif
      DO IAT=1,NAT
         CHGMLK(IAT)=-(X(IAT+IGAC-1)-ZAN(IAT))
      ENDDO
C
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) GO TO 400
      IF(LAST) GO TO 400
      IF(IPASS.EQ.2) GO TO 300
      IPASS=2
      IDEN=I41
      IGOC=I51
      IGAC=I61
      ELEC=BETA
      GO TO 100
  300 CONTINUE
C
C     ----- CALCULATE ORBITAL AND ATOMIC SPIN DENSITIES -----
C
      DO I=1,L1
         X(I-1+I50)=X(I-1+I50)-X(I-1+I51)
      ENDDO
      DO I=1,NAT
         X(I-1+I60)=X(I-1+I60)-X(I-1+I61)
      ENDDO
      if (ga_nodeid().eq.0) then
         WRITE(IW,9994)
         WRITE(IW,9997) (I, BFLAB(I),X(I-1+I50),I=1,NUM)
      endif
      if (ga_nodeid().eq.0) then
          WRITE(IW,9993)
          WRITE(IW,9992) (I,ANAM(I),BNAM(I),ZNUC(I),ZAN(I),
     1                X(I-1+I60),I=1,NAT)
      endif
C
C     ----- DO ALL ELECTRONS -----
C
      IPASS=1
      LAST=.TRUE.
      ELEC=All
      DO I=1,L2
         X(I-1+I40)=X(I-1+I40)+X(I-1+I41)
      ENDDO
      IDEN=I40
      IGOC=I50
      IGAC=I60
      GO TO 200
  400 CONTINUE
C
C     ----- BOND INDEX AND VALENCY ANALYSIS -----
C
      DO I=1,L2
         X(I+I41-1)=ZERO
      ENDDO
C
C     CALL HND_DAREAD(IDAF,IODA,X(I10),L2,12)

c     ----- transform S from square to triangular format -----

      ij = 0
      do i = 1,num
         do j = 1,i
            ij = ij + 1
            x(i10 +ij-1) = x(i98 + (j-1)*(num)+(i-1))
         enddo
      enddo
C
C*    CALL hnd_prpDEN(X(I40),X(I41),L2)

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i40),x(i41),x(i95),x(i96),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
C
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
         DO I=1,L2
            DUMA=X(I+I40-1)
            DUMB=X(I+I41-1)
            X(I+I40-1)=DUMA+DUMB
            X(I+I41-1)=DUMA-DUMB
         ENDDO
      ENDIF
      CALL HND_APSPOP(X(I10),X(I40),X(I41),X(I51),X(I83),
     1            X(I53),X(I63),X(I73),X(I93),X(I94),IA,NAT,NUM)
C
 1000 CONTINUE
C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_MULKEN, MA_POP_STACK OF INIT
     & FAILED',911)
C
      RETURN
 9999 FORMAT(//,10X,28(1H-),/,10X,'Mulliken population analysis',
     1 10X,A8,' electrons',/,10X,28(1H-))
 9998 FORMAT(/,10X,'----- Total gross population in ao -----',/)
 9997 FORMAT(5(I5,1X,A8,F11.5))
 9996 FORMAT(/,10X,'----- Condensed to atoms -----',/)
 9995 FORMAT(/,10X,'----- Total gross population on atoms ----',/)
 9994 FORMAT(/,10X,'----- Spin population in ao -----',/)
 9993 FORMAT(/,10X,'----- Atomic spin population -----',/)
 9992 FORMAT(10X,I5,2X,A8,A2,F6.1,F6.1,F12.5)
 9991 FORMAT(' Not enough memory in -Mulken- . The program will bypass',
     1       ' This step and continue.')
 9990 FORMAT(/,10X,'----- TOTAL GROSS POPULATION IN SPH.HARM. -----',/)
 9989 FORMAT(/,10X,'----- OVERLAP MATRIX IN -MULKEN- -----')
 9988 FORMAT(/,10X,'----- OVERLAP POPULATION -----')
 9987 FORMAT(/,10X,'----- SPIN POPULATION IN SPH.HARM. -----',/)
      END
C
      SUBROUTINE HND_ATMBFN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_BFNATM/LIMLOW(MXATOM),LIMSUP(MXATOM)
C
      DO I=1,NAT
         LIMLOW(I)=NUM+1
         LIMSUP(I)=0
         DO J=1,NSHELL
            IF(KATOM(J).EQ.I) THEN
               LIMLOW(I)=MIN(LIMLOW(I),KLOC(J)                )
               LIMSUP(I)=MAX(LIMSUP(I),KLOC(J)+KMAX(J)-KMIN(J))
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_APSPOP(S,D,DS,P,PS,B,F,V,SB,FV,IA,NAT,NUM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "global.fh"
#include "nwc_const.fh"
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_BFNATM/LIMLOW(MXATOM),LIMSUP(MXATOM)
      DIMENSION S(1),D(1),P(NUM,1),B(NAT,1),F(1),V(1),IA(1)
      DIMENSION DS(1),PS(NUM,1),SB(1),FV(1)
      DIMENSION II(4),JJ(4)
      DATA ZERO /0.0D+00/
      DATA TWO  /2.0D+00/
      DATA BTOL /0.1D+00/
C
C     ----- CLOSED AND OPEN SHELL SYSTEMS POPULATION ANALYSIS -----
C
C     ORIGINAL REFERENCE FOR NDO SYSTEMS:
C     D.R. ARMSTRONG, P.G. PERKINS AND J.J.P. STEWART
C     J.C.S. DALTON TRANS. 838 (1973)
C     EXTENSION  FIRST PRINCIPLES: I. MAYER, CHEM. PHYS. LETT.
C     97,270(1983).
C     NUMBER OF FREE ELECTRONS: O.G.STRADELLA, H.O. VILLAR, E.A.CASTRO
C     THEORET. CHIM. ACTA(BERL.) 70, 67(1986)
C
C       D * S = P (NON - HERMITEAN)
C
      DO I=1,NUM
         DO J=1,NUM
            DUMS=ZERO
            DUM =ZERO
            DO K=1,NUM
               IKMAX=MAX0(I,K)
               IKMIN=MIN0(I,K)
               KJMAX=MAX0(K,J)
               KJMIN=MIN0(K,J)
               IK=(IKMAX*(IKMAX-1))/2+IKMIN
               KJ=(KJMAX*(KJMAX-1))/2+KJMIN
               DUMS=DUMS+DS(IK)*S(KJ)
               DUM =DUM + D(IK)*S(KJ)
            ENDDO
            PS(I,J)=DUMS
             P(I,J)=DUM
         ENDDO
      ENDDO
C
C     ----- -B-, VALENCY , AND FREE ELECTRON -----
C
      DO IAT=1,NAT
         LLI=LIMLOW(IAT)
         LUI=LIMSUP(IAT)
         DO JAT=1,NAT
            LLJ=LIMLOW(JAT)
            LUJ=LIMSUP(JAT)
            DUM=ZERO
            IF(LUI.NE.0.AND.LUJ.NE.0) THEN
               DO L = LLJ,LUJ
                  DO K = LLI,LUI
                     DUM=DUM+P(L,K)*P(K,L)+PS(L,K)*PS(K,L)
                  ENDDO
               ENDDO
            ENDIF
            B(IAT,JAT)=DUM
         ENDDO
         DUMS=ZERO
         DUM =ZERO
         DO L = LLI,LUI
            DUM=DUM+TWO*P(L,L)
            DO M = LLI,LUI
               DUMS=DUMS-PS(L,M)*PS(M,L)
               DUM =DUM - P(L,M)* P(M,L)
            ENDDO
         ENDDO
         V(IAT)=DUM
         DUM=ZERO
         DO JAT=1,NAT
            IF(IAT.NE.JAT) THEN
               DUM=DUM-B(IAT,JAT)
            ENDIF
         ENDDO
         F(IAT)=DUM+(V(IAT)+DUMS+B(IAT,IAT))/TWO
         B(IAT,IAT)=ZERO
      ENDDO
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
C
C     ----- PRINT ALL BOND INDICES -----
C
      MAX=0
   10 MIN=MAX+1
      MAX=MAX+8
      IF(MAX.GT.NAT) MAX=NAT
      if (ga_nodeid().eq.0) WRITE(IW,9995)
      DO I =1,NAT
       if (ga_nodeid().eq.0) WRITE(IW,9998) (I,J,B(I,J),J=MIN,MAX)
      ENDDO
      IF(MAX.LT.NAT) GO TO 10
C
C     ----- NOW PRINT ONLY THE LARGEST ONE -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9993)
      KK=0
      DO I=1,NAT
         DO J=I+1,NAT
            IF(B(I,J).GT.BTOL) THEN
               KK = KK+1
               II(KK) = I
               JJ(KK) = J
               IF(KK.EQ.4) THEN

                 if (ga_nodeid().eq.0)  WRITE(IW,9994)
     1              (II(K),ANAM(II(K)),JJ(K),ANAM(JJ(K)),
     1                            B(II(K),JJ(K)),K=1,4)

                  KK = 0
               ENDIF
            ENDIF
         ENDDO
         IF(KK.NE.0) THEN

            if (ga_nodeid().eq.0) WRITE(IW,9994)
     1          (II(K),ANAM(II(K)),JJ(K),ANAM(JJ(K)),
     1                      B(II(K),JJ(K)),K=1,KK)

         ENDIF
         KK=0
      ENDDO
C
C     ----- SUM BOND INDICES FOR A GIVEN ATOM -----
C
      DO JAT=1,NAT
         DUM=ZERO
         DO IAT=1,NAT
            DUM=DUM+B(IAT,JAT)
         ENDDO
         B(JAT,JAT)=DUM
      ENDDO
C
C     ----- V(A) = N(A) + SUM OF B(A,B) IS THE MULLIKEN CHARGE -----
C           SPIN(A) = SUM OF B(A,B)-V(A)
C
      DO IAT=1,NAT
         SB(IAT)=B(IAT,IAT)+F(IAT)
         FV(IAT)=V(IAT)-B(IAT,IAT)
      ENDDO
C
C     ----- PRINT ANALYSIS -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9997)
      DO I =1,NAT
         if (ga_nodeid().eq.0)
     & WRITE(IW,9996) I,ANAM(I),BNAM(I),V(I),F(I),B(I,I),SB(I),FV(I)
      ENDDO
C
      RETURN
 9999 FORMAT(/,10X,'----- Bond indices -----')
 9998 FORMAT(8(1X,I2,'-',I2,F10.5))
 9997 FORMAT(/,70X,'  Free electrons',5X,'        Valency',/,
     1         36X,'  Number of   ',3X,'   Sum of    ',
     2          6X,'+ Bond indices ',5X,'   - Bond indices    ',/,
     3         21X,'   Valency',5X,'Free electrons',4X,'Bond indices',
     4          5X,'=Mulliken charge',5X,'= Net spin population')
 9996 FORMAT(1X,I5,2X,A8,A2,3X,F10.5,5X,F10.5,7X,F10.5,8X,F10.5,
     1      14X,F10.5)
 9995 FORMAT(/)
 9994 FORMAT(5(I3,1X,A4,'-',I3,1X,A4,F10.5))
 9993 FORMAT(/,1X,'Large bond indices',/,1X,18(1H-))
      END
C
      SUBROUTINE HND_ORBPOP(rtdb,geom,basis,
     &                      S,V,T,IA,NPOP,L1,NDIM,MOFLG,SPI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXIODA=255)
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb,geom,basis
c*    logical  status
c*    integer  g_vecs
c*    character*255 movecs_in,movecs_out
c*    character*255 movecs
c*    logical  movecs_read
c*    external movecs_read
c*    integer  ga_create_atom_blocked
c*    external ga_create_atom_blocked
c*    double precision  occ(mxbfn)
      character*8 sphlab,dum
      character*8 anam,bflab
      character*2 bnam
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2, MP4, WFNMC, WFNCI
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_BASSPH/ISPHER
      COMMON/HND_MOLSPH/SPHLAB(MXBFN)
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_WFNOPT/WFNTYP
      dimension e(mxbfn)
      DIMENSION S(1),T(1),IA(1)
      DIMENSION   V(NDIM,1)
      DIMENSION SPI(NDIM,1)
      DATA ZERO    /0.0D+00/
      DATA SCF,UHF /'SCF     ','UHF     '/
      DATA MP2     /'MP2     '/
      DATA MP4     /'MP4     '/
      DATA WFNMC   /'MC      '/
      DATA WFNCI   /'CI      '/
C
      DBUG=.FALSE.
      if (ga_nodeid().eq.0) WRITE(IW,9999)
C
      IPASS=1
C
C     -----  FOR -MC- OR -CI- WAVEFUNCTIONS -----
C           SELECT ORBITAL SET FOR ANALYSIS
C
      IF(WFNTYP.EQ.WFNMC.OR.WFNTYP.EQ.WFNCI) THEN
         IF(MOFLG.EQ.0) THEN
            if (ga_nodeid().eq.0) WRITE(IW,9996)
            NDAF=19
         ELSE
            if (ga_nodeid().eq.0) WRITE(IW,9995)
            NDAF=15
         ENDIF
      ELSE
         NDAF=15
      ENDIF
C
C     -----  FOR -MP2- OR -MP4- CALCULATIONS   -----
C           ANALYSIS IS FOR THE -SCF- ORBITALS
C
      IF(WFNTYP.EQ.MP2.OR.WFNTYP.EQ.MP4) THEN
         if (ga_nodeid().eq.0) WRITE(IW,9993)
      ENDIF
C
      NA=0
      NORB=NA
      IF(NPOP.GT.0) NORB=NPOP
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF.and.ga_nodeid().eq.0).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF.and.ga_nodeid().eq.0))
     1    WRITE(IW,9998)
C
  100 CONTINUE
C
C     ----- READ ORBITALS -----
C
      L3=L1*L1
C     CALL HND_DAREAD(IDAF,IODA,V,L3,NDAF)

c     ----- figure output input/output MO vectors -----

c      ----- read eig-vec; eng-val in t for tmp. -----

       call hnd_get_eigens(rtdb,geom,basis,v,e,num)
c*     call hnd_prevl(v,e,norb,num,ndim)
C
C     ----- TRANSFORM TO SPHERICAL HARMONICS IF NEEDED -----
C     -----    MULTIPLY BY   SPH**(-1)   ON THE LEFT   -----
C
      IF(ISPHER.GE.0) THEN
         DO J=1,NORB
            DO I=1,NUM

C              DUM=ZERO
               DDUM=ZERO

               DO K=1,NUM

C                  DUM=DUM+SPI(I,K)*V(K,J)
                  DDUM=DDUM+SPI(I,K)*V(K,J)

               ENDDO

C              T(I)=DUM
               T(I)=DDUM

            ENDDO
            DO I=1,NUM
               V(I,J)=T(I)
            ENDDO
         ENDDO
C
C        ----- PRINT OUT ORBS. IN SPH.HARM. BASIS -----
C        ----- ( WITH CORRECT LABELS )            -----
C
         if (ga_nodeid().eq.0) WRITE(IW,9994)
c%% T = E ?
C        CALL HND_DAREAD(IDAF,IODA,T,L1,NDAF+2)
c%%
         DO I=1,NUM
            DUM      = BFLAB(I)
            BFLAB(I)=SPHLAB(I)
            SPHLAB(I)=DUM
         ENDDO
         CALL HND_PREVL(V,T,NORB,NUM,NDIM)
         DO I=1,NUM
            DUM      =SPHLAB(I)
            SPHLAB(I)= BFLAB(I)
            BFLAB(I)=DUM
         ENDDO
         if (ga_nodeid().eq.0) WRITE(IW,9992)

      ENDIF
C
C     ----- CALCULATE POPULATION FOR EACH ORBITAL  -----
C     ----- IF NEEDED, -S- IS ALREADY IN SPH.HARM. -----
C
      DO I=1,NORB
         DO J=1,NUM

            DDUM=ZERO
            DO K=1,NUM
               JJ=MAX0(J,K)
               KK=MIN0(J,K)
               JK=(JJ*(JJ-1))/2+KK
               DDUM=DDUM+V(K,I)*S(JK)*V(J,I)
            ENDDO
            T(J)=DDUM

         ENDDO
         DO J=1,NUM
            V(J,I)=T(J)
         ENDDO
      ENDDO
C
C     ----- PRINT OUT ( SWITCH LABELS IF SPH.HARM. ) -----
C
C     CALL HND_DAREAD(IDAF,IODA,T,L1,NDAF+2)
      IF(ISPHER.LT.0) THEN

C        CALL HND_PREVL(V,T,NORB,NUM,NDIM)
         call hnd_prevl(v,e,norb,num,ndim)

      ELSE

         DO I=1,NUM
            DUM      = BFLAB(I)
            BFLAB(I)=SPHLAB(I)
            SPHLAB(I)=DUM
         ENDDO
         CALL HND_PREVL(V,T,NORB,NUM,NDIM)
         DO I=1,NUM
            DUM      =SPHLAB(I)
            SPHLAB(I)= BFLAB(I)
            BFLAB(I)=DUM
         ENDDO

      ENDIF
C
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) GO TO 200
      IF(IPASS.EQ.2) GO TO 200
      IPASS=2
      NDAF=19
      NB=0
      NORB=NB
      IF(NPOP.GT.0) NORB=NPOP
      if (ga_nodeid().eq.0) WRITE(IW,9997)
      GO TO 100
  200 CONTINUE
C
      RETURN
 9999 FORMAT(/,10X,44(1H-),
     1       /,10X,'Mulliken populations over molecular orbitals',
     2       /,10X,44(1H-),/)
 9998 FORMAT(/,' ----- Alpha set ----- ')
 9997 FORMAT(/,' ----- Beta set ----- ')
 9996 FORMAT(  ' ( Core + natural orbitals of -MC- or -CI- wfn )',/)
 9995 FORMAT(  ' ( Core + valence orbitals of -MC- or -CI- wfn )',/)
 9994 FORMAT(/,10X,'----- Orbitals in sph.harm. basis -----')
 9993 FORMAT(  ' Analysis is for the -(CAS)SCF- wavefunction',
     1         ' In this -(MR)MP- calculation .'                )
 9992 FORMAT(/,10X,'----- Now, the analysis ...       -----')
      END
      SUBROUTINE HND_OVLPOP(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1)
      DO I=1,N
         A(I)=A(I)*B(I)
      ENDDO
      RETURN
      END
      SUBROUTINE HND_SHLPOP(AOPOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
#include "global.fh"
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (ZERO=0.0D+00)
C     LOGICAL DEBUG
      CHARACTER*8 ERRMSG
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_BASSPH/ISPHER
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION AOPOP(1)
      DIMENSION SHL(5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-SHLPOP-'/
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
      IF(ISPHER.LT.0) THEN
C
C     ----- IN -AO- BASIS -----
C
         DO IAT=1,NAT
            DO ISHL=1,5
               SHL(ISHL)=ZERO
            ENDDO
            DO ISH=1,NSHELL
               IF(KATOM(ISH).EQ.IAT) THEN
                  LIT=KTYPE(ISH)
                  MIN=KMIN(ISH)
                  MAX=KMAX(ISH)
                  LOC=KLOC(ISH)-MIN
                  IF(LIT.EQ.2) THEN
                     IF(MIN.EQ.1) THEN
                        SHL(1)=SHL(1)+AOPOP( 1+LOC)
                        MIN=MIN+1
                     ENDIF
                  ENDIF
                  DUM=ZERO
                  DO I=MIN,MAX
                     DUM=DUM+AOPOP(I+LOC)
                  ENDDO
                  SHL(LIT)=SHL(LIT)+DUM
               ENDIF
            ENDDO
            SUM=SHL(1)+SHL(2)+SHL(3)+SHL(4)+SHL(5)

            if (ga_nodeid().eq.0) WRITE(IW,9998)
     &         IAT,ANAM(IAT),BNAM(IAT),(SHL(I),I=1,5),SUM

         ENDDO
      ELSE
C
C     ----- IN -SPH.HARM.- BASIS -----
C
         DO IAT=1,NAT
            DO ISHL=1,5
               SHL(ISHL)=ZERO
            ENDDO
            DO ISH=1,NSHELL
               IF(KATOM(ISH).EQ.IAT) THEN
                  LIT=KTYPE(ISH)
                  MIN=KMIN(ISH)
                  MAX=KMAX(ISH)
                  LOC=KLOC(ISH)-MIN
                  IF(LIT.EQ.2) THEN
                     IF(MIN.EQ.1) THEN
                        SHL(1)=SHL(1)+AOPOP( 1+LOC)
                        MIN=MIN+1
                     ENDIF
                  ELSEIF(LIT.EQ.3) THEN
                     SHL(1)=SHL(1)+AOPOP(10+LOC)
                     MAX=MAX-1
                  ELSEIF(LIT.EQ.4) THEN
                     SHL(2)=SHL(2)+AOPOP(18+LOC)
     1                            +AOPOP(19+LOC)+AOPOP(20+LOC)
                     MAX=MAX-3
                  ELSEIF(LIT.EQ.5) THEN
                     SHL(1)=SHL(1)+AOPOP(30+LOC)
                     SHL(3)=SHL(3)+AOPOP(31+LOC)+AOPOP(32+LOC)
     1                            +AOPOP(33+LOC)+AOPOP(34+LOC)
     2                            +AOPOP(35+LOC)
                     MAX=MAX-6
                  ENDIF
                  DUM=ZERO
                  DO I=MIN,MAX
                     DUM=DUM+AOPOP(I+LOC)
                  ENDDO
                  SHL(LIT)=SHL(LIT)+DUM
               ENDIF
            ENDDO
            SUM=SHL(1)+SHL(2)+SHL(3)+SHL(4)+SHL(5)

            if (ga_nodeid().eq.0) WRITE(IW,9998) 
     &         IAT,ANAM(IAT),BNAM(IAT),(SHL(I),I=1,5),SUM

         ENDDO
      ENDIF
C
      RETURN
 9999 FORMAT(/,1X,'Total S,P,D,... shell population',/,1X,32(1H-),
     1       /,4X,'     Atom      ','    S     ','    P     ',
     2    '    D     ','    F     ','    G     ','  Total   ',
     3       /,1X,76(1H-))
 9998 FORMAT(1X,I4,1X,A8,A2,6F10.5)
      END
      SUBROUTINE HND_ATMPOP(A,IA,B,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      COMMON/HND_BFNATM/LIMLOW(MXATOM),LIMSUP(MXATOM)
      DIMENSION A(1),B(1),IA(1)
      DATA ZERO /0.0D+00/
C
      DO I=1,NAT
         DO J=1,I
            IJ=(I*(I-1))/2+J
            B(IJ)=ZERO
         ENDDO
      ENDDO
C
      DO I=1,NAT
         I1=LIMLOW(I)
         I2=LIMSUP(I)
         IF(I2.NE.0) THEN
            DO J=1,NAT
               J1=LIMLOW(J)
               J2=LIMSUP(J)
               IF(J2.NE.0) THEN
                  DUM=ZERO
                  DO K=I1,I2
                     DO L=J1,J2
                        KK=MAX0(K,L)
                        LL=MIN0(K,L)
                        KL=(KK*(KK-1))/2+LL
                        DUM=DUM+A(KL)
                     ENDDO
                  ENDDO
                  II=MAX0(I,J)
                  JJ=MIN0(I,J)
                  IJ=(II*(II-1))/2+JJ
                  B(IJ)=DUM
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_GROSSC(A,B,IA,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1),IA(1)
      DATA ZERO /0.0D+00/
C
      DO I=1,N
         DUM=ZERO
         DO J=1,N
            II=MAX0(I,J)
            JJ=MIN0(I,J)
            IJ=(II*(II-1))/2+JJ
            DUM=DUM+A(IJ)
         ENDDO
         B(I)=DUM
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_BOYLOC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     BOYLOC CALCULATES LOCALIZED MOLECULAR ORBITALS
C     BY THE METHOD OF BOYS. THIS IMPLEMENTATION IS
C     ESSENTIALLY QCPE PROGRAM 354, MODIFIED TO RUN
C     AS PART OF HONDO.
C
      PARAMETER (MXIODA=255)
      PARAMETER (MXBFN=3072)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
c%%
      CHARACTER*8 VECWRD, ENDWRD
c%%
      LOGICAL OUT
      LOGICAL NOBOYS
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_IJPAIR/IA(1)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(MXATOM),
     1                                         C(3,MXATOM)
      COMMON/HND_SCM/X(1)
      DIMENSION IFRZ(MXBFN),IPAIR(MXBFN)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-BOYLOC-'/
      DATA IDAR   /32/
      DATA VECWRD,ENDWRD /' $VEC   ',' $END   '/
C
      NAMELIST /BOYS/ NORB,IFRZ,IPAIR
C
C     ----- SET PARAMETERS -----
C
      OUT   =NPRINT.EQ. 6
      NOBOYS=NPRINT.EQ.-5
      IF(NOBOYS) RETURN
C
      WRITE(IW,9999)
      NORB=NA
      MAXORB=MXBFN
      DO 5 I=1,MAXORB
          IFRZ(I)=0
    5    IPAIR(I)=0

C
C     ----- READ NAMELIST - $BOYS - -----
C
      REWIND IR
      READ(IR,BOYS,END=10,ERR=10)
   10 CONTINUE
      NFRZ=0
      DO 30 I=1,MAXORB
      IF(IFRZ(I).EQ.0) GO TO 40
      NFRZ=NFRZ+1
   30 CONTINUE
   40 CONTINUE
C
C     ----- GET CORE MEMORY -----
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM* NUM
      M1= NORB
      M2=(NORB*(NORB+1))/2
C
      I10 = 1
      I20 = I10+MAX0(L1*M1,L2)
      I30 = I20+MAX0(L1*M1,L2)
      I40 = I30+L1
      I50 = I40
      I60 = I50+MAX0(3*M2,L2)
      I70 = I60+L1+L1+L1
      I80 = I70+L1
      I90 = I80+L1
      I100 = I90+L1
C
c%%
C     CALL HND_CMEM(LOADCM)
C     LAST = I100
C     LAST = LAST - 1
C     LOC10 = LOCCM(X(I10))
C     NEED = LOC10 + LAST
C     IF(NEED.LE.MAXCOR) GO TO 50
C     WRITE(IW,9995)
C     CALL HND_NEEDCM(LAST,NEED)
C     CALL HND_HNDERR(3,ERRMSG)
C  50 CONTINUE
c%%
C
      CALL HND_DAREAD(IDAF,IODA,X(I10),L3,15)
      LENGTH = M1*L1
      CALL HND_DAWRIT(IDAF,IODA,X(I10),LENGTH,IDAR,NAV)
C
C     ----- CALL HND_LOCALIZATION PROGRAM -----
C
      CALL HND_BOYORB(NORB,X(I10),X(I10),X(I20),X(I20),X(I20),X(I30),
     1     IA,X(I50),X(I60),X(I80),X(I90),L1,L2,M1,M2,
     1     IFRZ,NFRZ)
      NDUM=L1*M1
      DO 110 IDUM=1,NDUM
  110 X(IDUM+I10-1)=X(IDUM+I20-1)
C
C     ----- GENERATE ANTI-BONDING ORBITALS -----
C
      NPAIR=0
      DO 120 I=1,512
      IF(IPAIR(I).EQ.0) GO TO 130
      NPAIR=NPAIR+1
  120 CONTINUE
  130 CONTINUE
      IF(NPAIR.EQ.0) GO TO 140
C
      CALL HND_BOYBFN
      CALL HND_BOYGVB(NORB,NPAIR,NUM,IPAIR,NAT,X(I10))
C
  140 CONTINUE
C
C     ----- PRINT AND PUNCH -----
C
      NORB=NORB+NPAIR
      CALL HND_PRSQL(X(I10),NORB,L1,L1)
c--   WRITE(IP,9997)
c--   WRITE(IP,9996) VECWRD
c--   CALL HND_PUSQL(X(I10),NORB,L1,L1)
c--   WRITE(IP,9996) ENDWRD
C
      NA=NORB
      NB=NORB
      CALL HND_DAWRIT(IDAF,IODA,X(I10),L3,15,NAV)
C
C     ----- RETURN CORE MEMORY -----
C
c%%
C     CALL HND_CMEM(NGOTCM)
C     IF(NGOTCM.NE.LOADCM) CALL HND_SETC(LOADCM)
c%%
      RETURN
 9999 FORMAT(/,10X,23(1H-),/,10X,'BOYS LOCALIZED ORBITALS',/,10X,
     1 23(1H-) )
 9997 FORMAT(' ----- BOYS LOCALIZED ORBITALS ----- ')
 9996 FORMAT(A8)
 9995 FORMAT(' NOT ENOUGH CORE IN -BOYLOC- . STOP ')
      END
      SUBROUTINE HND_BOYORB(NORB,CC,CL,B,DMAO,IORD,IIR,IA,RIJ,RI,
     1 QPIX,QPJX,L1,L2,M1,M2,IFRZ,NFRZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(MXATOM),
     1                                       CRD(3,MXATOM)
      COMMON/HND_OUTPUT/NPRINT
      DIMENSION CC(L1,M1),CL(L1,M1),DMAO(L2),IIR(L1)
      DIMENSION QPIX(L1),QPJX(L1)
      DIMENSION IA(1),B(L1,M1),IORD(L1),RI(L1,3),RIJ(M2,3)
      DIMENSION IFRZ(1)
      DATA IDAR  /32/
      DATA FAC   /2.54176568D+00/
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA TWO   /2.0D+00/
      DATA FOUR  /4.0D+00/
      DATA TOL1  /1.0D-10/
      DATA TOL2  /1.0D-08/
      DATA TOL3  /1.0D-06/
      DATA MAXIT /500/
      DATA NPASS /3/
C
      OUT=NPRINT.EQ.6
C
C     ----- READ MOLECULAR ORBITALS                             -----
C     ----- SET UP FROZEN ORBITALS                              -----
C     ----- CALCULATE (I/MX/J) FOR MO'S I AND J AND MX=MX,MY,MZ -----
C     ----- RI ARE DIAGONAL ELEMENTS RIJ ARE THE OFF DIAGONAL   -----
C
      NOCC=NORB*NUM
      CALL HND_DAREAD(IDAF,IODA,CC,NOCC,IDAR)
      IF(OUT) CALL HND_PRSQL(CC,NORB,NUM,L1)
C
      IF (NFRZ .EQ. 0) GO TO 240
      DO 220 I = 1,NFRZ
      K = IFRZ(I)
      DO 200 J = 1,NUM
  200 CC(J,K) = ZERO
  220 CONTINUE
C
  240 CONTINUE
      DO 300 IXYZ=1,3
      CALL HND_DAREAD(IDAF,IODA,DMAO,L2,IDAR+IXYZ)
      MM = 0
      DO 280 I = 1,NORB
      DO 280 J = 1,I
      IF (I .NE. J) MM = MM+1
      SUM = ZERO
      DO 260 K = 1,NUM
      DO 260 L = 1,NUM
      KK=MAX0(K,L)
      LL=MIN0(K,L)
      KL=(KK*(KK-1))/2+LL
  260 SUM = SUM+CC(K,I)*CC(L,J)*DMAO(KL)
      IF (I .EQ. J) RI(I,IXYZ) = SUM
      IF (I .NE. J) RIJ(MM,IXYZ) = SUM
  280 CONTINUE
  300 CONTINUE
C
C     ----- INITIALIZE RANDOM NUMBER GENERATOR -----
C           INITIALIZE ROTATION MATRIX
C
      XX=HND_BOYRND(ONE,CC,L1)
      DO 340 I = 1,NORB
      DO 320 J = 1,I
      CL(I,J) = ZERO
  320 CL(J,I) = ZERO
  340 CL(I,I) = ONE
C
C     ----- LOOP OVER PASSES -----
C
      NREDO = 0
  350 CONTINUE
C
C     ----- LOOP OVER 2*2 ROTATIONS FOR THIS PASS -----
C
      ITER = 0
      SHIFT =  DATAN(ONE)
  360 CHANGE = ZERO
      ITER = ITER+1
      DO 380 I = 1,NORB
  380 IIR(I) = I
      NNN = NORB
      DO 400 I = 1,NORB
      XX = HND_BOYRND(ZERO,CC,L1)
      III = (XX* DBLE(NNN)+ONE)
      IORD(I) = IIR(III)
      IIR(III) = IIR(NNN)
      NNN = NNN-1
  400 CONTINUE
C
C     ----- FOR EACH PAIR OF ORBITALS A TWO DIMENSIONAL -----
C           UNITARY TRANSFORMATION IS PERFORMED.
C           THE TRANSFORMATION IS
C           PSI PRIME(I) =  COS(T)*PSI(I)+SIN(T)*PSI(J)
C           PSI PRIME(J) = -SIN(T)*PSI(I)+COS(T)*PSI(J)
C           THE BOYS METHOD REQUIRES THAT -T- BE SUCH AS
C           TO MAXIMIZE THE SUM OF THE SQUARES OF THE
C           ONE CENTER MOLECULAR ORBITAL DIPOLE MOMENTS
C
      DO 920 III = 1,NORB
      I = IORD(III)
      IM = I-1
      JM = 1
      IJM = IM*(IM-1)/2+1
      RM = ZERO
      TM = ZERO
      SM = ZERO
      CM = ONE
      DO 580 J = 1,NORB
      IF (I-J) 440,580,420
  420 IJ = IM*(IM-1)/2+J
      GO TO 460
  440 IJ = (J-1)*(J-2)/2+I
  460 T = ZERO
      TX = ZERO
      DO 480 IXYZ = 1,3
      T = T+FOUR*RIJ(IJ,IXYZ)**2-RI(I,IXYZ)**2-RI(J,IXYZ)**2+
     1      TWO*RI(I,IXYZ)*RI(J,IXYZ)
  480 TX = TX+RIJ(IJ,IXYZ)*(RI(J,IXYZ)-RI(I,IXYZ))
      IF(( ABS(T).LE.TOL1).AND.( ABS(TX).LE.TOL1)) GO TO 580
      TX = FOUR*TX
      T =  DATAN2(TX,T)
      T = T/FOUR
      ASIGN = ONE
      IF (T .GT. ZERO) ASIGN = -ONE
      T = T+ASIGN*SHIFT
      ITIM = 0
  500 ITIM = ITIM + 1
      C = DCOS(T)
      S = DSIN(T)
      RIN = ZERO
      DO 520 IXYZ = 1,3
      QPI = C*C*RI(I,IXYZ)+S*S*RI(J,IXYZ)+TWO*C*S*RIJ(IJ,IXYZ)
      QPJ = C*C*RI(J,IXYZ)+S*S*RI(I,IXYZ)-TWO*C*S*RIJ(IJ,IXYZ)
      RIN = RIN+QPI*QPI+QPJ*QPJ-RI(I,IXYZ)**2-RI(J,IXYZ)**2
  520 CONTINUE
      TTEST =  ABS(T)-SHIFT
      IF(( ABS(T).LE.TOL2).OR.( ABS(TTEST).LE.TOL2)) GO TO 560
      IF (RIN.GE.(-TOL2)) GO TO 560
      IF (ITIM.LE.1) GO TO 540
      WRITE (IW,9008) I,J
      WRITE (IW,9028) T,S,C,RIN
      GO TO 940
  540 CONTINUE
      ASIGN = ONE
      IF (T.GT.ZERO) ASIGN = -ONE
      T = T+ASIGN*SHIFT
      GO TO 500
  560 CONTINUE
      IF (RIN.LE.RM) GO TO 580
      IJM = IJ
      RM = RIN
      JM = J
      SM = S
      CM = C
      TM = T
  580 CONTINUE
      T = TM
      RIN = RM
      S = SM
      C = CM
      J = JM
      IJ = IJM
      CHANGE = CHANGE+T*T
      DO 880 IXYZ = 1,3
      QPI = C*C*RI(I,IXYZ)+S*S*RI(J,IXYZ)+TWO*C*S*RIJ(IJ,IXYZ)
      QPJ = C*C*RI(J,IXYZ)+S*S*RI(I,IXYZ)-TWO*C*S*RIJ(IJ,IXYZ)
      QPIJ = (C*C-S*S)*RIJ(IJ,IXYZ)+C*S*(RI(J,IXYZ)-RI(I,IXYZ))
      DO 720 K = 1,NORB
      IF (I-K) 600,720,620
  600 IK = (K-1)*(K-2)/2+I
      GO TO 640
  620 IK = (I-1)*(I-2)/2+K
  640 IF (J-K) 660,720,680
  660 JK = (K-1)*(K-2)/2+J
      GO TO 700
  680 JK = (J-1)*(J-2)/2+K
  700 QPIX(K) = C*RIJ(IK,IXYZ)+S*RIJ(JK,IXYZ)
      QPJX(K) = C*RIJ(JK,IXYZ)-S*RIJ(IK,IXYZ)
  720 CONTINUE
      DO 860 K = 1,NORB
      IF (I-K) 740,860,760
  740 IK = (K-1)*(K-2)/2+I
      GO TO 780
  760 IK = (I-1)*(I-2)/2+K
  780 IF (J-K) 800,860,820
  800 JK = (K-1)*(K-2)/2+J
      GO TO 840
  820 JK = (J-1)*(J-2)/2+K
  840 RIJ(IK,IXYZ) = QPIX(K)
      RIJ(JK,IXYZ) = QPJX(K)
  860 CONTINUE
      RIN = RIN+QPI+QPJ-RI(I,IXYZ)-RI(J,IXYZ)
      RI(I,IXYZ) = QPI
      RI(J,IXYZ) = QPJ
  880 RIJ(IJ,IXYZ) = QPIJ
C
C     ----- INCLUDE 2*2 ROTATION IN ROTATION MATRIX -----
C
      DO 900 K = 1,NORB
      C1 =  C*CL(K,I)+S*CL(K,J)
      C2 = -S*CL(K,I)+C*CL(K,J)
      CL(K,I) = C1
  900 CL(K,J) = C2
  920 CONTINUE
C
C     ----- IF CONVERGENCE HAS NOT BEEN REACHED          -----
C           START ANOTHER SERIES OF TWO CENTER ROTATIONS
C
      IF(NORB.EQ.1) CHANGE =  SQRT(TWO*CHANGE/(NORB         ))
      IF(NORB.GT.1) CHANGE =  SQRT(TWO*CHANGE/(NORB*(NORB-1)))
      IF (ITER .GT. MAXIT) GO TO 940
      IF (ITER .GE. 1 .AND. OUT) WRITE (IW,9048) ITER, CHANGE
      IF (CHANGE .GE. TOL1) GO TO 360
  940 CONTINUE
      IF(OUT) WRITE (IW,9068) ITER
      IF (ITER .LT. MAXIT .AND. CHANGE .LE. TOL3) GO TO 950
      NREDO = NREDO+1
      IF (NREDO .LT. NPASS) GO TO 350
      WRITE (IW,9088)
  950 CONTINUE
C
C     ----- PRINT CENTROIDS OF CHARGE -----
C
      IF(OUT) THEN
         WRITE (IW,9048) ITER, CHANGE
         WRITE (IW,9188)
         DO 960 I=1,NORB
  960       WRITE(IW,9168) I,(RI(I,K),K=1,3)
      ENDIF
C
C     ----- NORMALIZE ROTATION MATRIX -----
C
      DO 970 I = 1,NORB
      DO 970 J = 1,NORB
      B(J,I) = CL(J,I)
  970 CONTINUE
      DO 980 J=1,NORB
      DUM=ZERO
      DO 975 I=1,NORB
  975 DUM=DUM+B(I,J)**2
      DUM=ONE/ SQRT(DUM)
      DO 980 I=1,NORB
      B(I,J)=B(I,J)*DUM
  980 CONTINUE
C
C     ----- FORCE SYMMETRY -----
C
      DO 985 J=1,NORB
      DO 985 I=1,NORB
      IF( ABS(B(I,J)).LT.TOL3) B(I,J)=ZERO
  985 CONTINUE
      DO 995 J=1,NORB
      DO 995 JJ=J,NORB
      IF(JJ.EQ.J) GO TO 995
      DO 990 I=1,NORB
      DUM= ABS(B(I,J))- ABS(B(I,JJ))
      IF( ABS(DUM).GT.TOL3) GO TO 990
      B(I,JJ)= DSIGN(B(I,J),B(I,JJ))
  990 CONTINUE
  995 CONTINUE
C
C     ----- APPLY ROTATION MATRIX -----
C
      CALL HND_DAREAD(IDAF,IODA,CC,NOCC,IDAR)
      IF(OUT) CALL HND_PRSQ(B,NORB,NORB,L1)
      DO 1060 I = 1,NORB
      DO 1000 J = 1,NUM
      RI(J,1) = ZERO
 1000 CONTINUE
      DO 1020 J = 1,NORB
      DO 1020 K = 1,NUM
      RI(K,1) = RI(K,1)+CC(K,J)*B(J,I)
 1020 CONTINUE
      DO 1040 J = 1,NUM
      B(J,I) = RI(J,1)
 1040 CONTINUE
 1060 CONTINUE
      IF (NFRZ .EQ. 0) GO TO 1120
      CALL HND_DAREAD(IDAF,IODA,CC,NOCC,IDAR)
      DO 1100 I = 1,NFRZ
      K = IFRZ(I)
      DO 1080 J = 1,NUM
 1080 B(J,K) = CC(J,K)
 1100 CONTINUE
 1120 CONTINUE
C
C     ----- GENERATE INDIVIDUAL MO CONTRIBUTIONS -----
C           TO THE DIPOLE MOMENT
C
      DO 1240 IXYZ = 1,3
      CALL HND_DAREAD(IDAF,IODA,CL,L2,IDAR+IXYZ)
      DO 1200 I = 1,NORB
      SUM = ZERO
      DO 1180 K = 1,NUM
      DO 1180 L = 1,NUM
      KK=MAX0(K,L)
      LL=MIN0(K,L)
      KL=(KK*(KK-1))/2+LL
      SUM = SUM + B(K,I)*B(L,I)*CL(KL,1)
 1180 CONTINUE
      RI(I,IXYZ) = SUM
 1200 CONTINUE
      FACT = -FAC
      DO 1220 I = 1,NORB
      RI(I,IXYZ) = RI(I,IXYZ)*FACT
 1220 CONTINUE
 1240 CONTINUE
      IF(OUT) THEN
         WRITE (IW,9148)
         DO 1260 I = 1,NORB
            WRITE (IW,9168) I,(RI(I,IXYZ),IXYZ = 1,3)
 1260    CONTINUE
      ENDIF
C
      RETURN
 9008 FORMAT(
     1     53H NO ROTATION INCREASES INTEGRALS  --- PROGRAM ABORTED,/,
     1     10X,3HI= ,I3,5X,3HJ= ,I3,///)
 9028 FORMAT(5X,8HTHETA = ,G20.10/5X,12HSIN(THETA)= ,F10.7,
     1     15H   COS(THETA)= ,F10.7/5X,29HTOTAL CHANGE TO THIS POINT = ,
     1     G20.10)
 9048 FORMAT(10X,17HITERATION NUMBER ,I3,12H   CHANGE = ,G20.10)
 9068 FORMAT (25X,21HPROGRAM STOPPED AFTER ,I3,11H ITERATIONS/)
 9088 FORMAT(/,' ----- NOTE, LOCALIZATION HAS BEEN UNSUCCESSFUL -----',
     1       /,' ----- EVEN AFTER THREE PASSES. PRINTED ORBITALS-----',
     2       /,' ----- ARE NOT CONVERGED LOCALIZED ORBITALS.    -----')
 9148 FORMAT(
     1 /,5X,'MO CONTRIBUTIONS TO THE ELECTRIC DIPOLE MOMENT (DEBYE)',
     1 /,5X,54(1H-),
     1 /,4X,3HNO.,5X,7X,1HX,14X,1HY,14X,1HZ,/)
 9168 FORMAT(I4,10X,F10.5,5X,F10.5,5X,F10.5)
 9188 FORMAT(
     1 /,5X,'CENTROIDS OF LOCALIZED MO (DIAGONAL ELEMENTS) (AU)',
     1 /,5X,50(1H-),
     1 /,4X,3HNO.,5X,7X,1HX,14X,1HY,14X,1HZ,/)
 9228 FORMAT(2X,I2,1H=,10(2X,F8.5)/6(5X,10(2X,F8.5)/))
      END
      FUNCTION HND_BOYRND(XX,D,L1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      COMMON/HND_INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(MXATOM),
     1                                         C(3,MXATOM)
      DIMENSION D(L1,L1)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
      DATA PI   /3.1416D+00/
C
      U = ZERO
      IF (XX .EQ. ZERO) GO TO 100
      N = IABS(NA-NUM)+1
      M = MIN0(N+5,NA)
      XY = D(N,M)* DATAN(ONE)
      U = (PI+XY)**5
      IU = U
      XY =  DBLE( IU)
      U = U-XY
      HND_BOYRND = U
      RETURN
  100 U = (PI+U)**5
      IU = U
      XY =  DBLE( IU)
      U = U-XY
      HND_BOYRND = U
      RETURN
      END
      SUBROUTINE HND_BOYBFN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)      
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_BFNBOY/LIMLOW(MXATOM),LIMSUP(MXATOM)
C
      DO I=1,NAT
         LIMLOW(I)=NUM+1
         LIMSUP(I)=0
         DO J=1,NSHELL
            IF(KATOM(J).EQ.I) THEN
               LIMLOW(I)=MIN(LIMLOW(I),KLOC(J)                )
               LIMSUP(I)=MAX(LIMSUP(I),KLOC(J)+KMAX(J)-KMIN(J))
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_BOYGVB(NORB,NPAIR,NUM,IPAIR,NAT,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      COMMON/HND_BFNBOY/LIMLOW(MXATOM),LIMSUP(MXATOM)
      DIMENSION V(NUM,1),IPAIR(1)
      DATA ZERO /0.0D+00/
C
      DO M=1,NPAIR
         IORB=IPAIR(M)
C
C     ----- FIND ATOM WITH LARGEST MO COEFFICIENT -----
C
         DMAX=ZERO
         DO JAT=1,NAT
            MIN=LIMLOW(JAT)
            MAX=LIMSUP(JAT)
            DUM=ZERO
            IF(MAX.NE.0) THEN
               DO K=MIN,MAX
                  DUM=DUM+V(K,IORB)**2
               ENDDO
            ENDIF
            IF(DUM.GT.DMAX) THEN
               IAT=JAT
               DMAX=DUM
            ENDIF
         ENDDO
C
C     ----- FLIP SIGN OF COEFFICIENTS ON THAT ATOM -----
C           TO CREATE THE GVB CORRELATING ORBITAL
C
         DO K=1,NUM
            V(K,NORB+M)=V(K,IORB)
         ENDDO
         MIN=LIMLOW(IAT)
         MAX=LIMSUP(IAT)
         DO K=MIN,MAX
            V(K,NORB+M)=-V(K,IORB)
         ENDDO
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE HND_ELPMAP(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      parameter  (mxefc=3*7*7*7)
      integer rtdb,basis,geom
      logical status
      PARAMETER (MXIODA=255)
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2
      character*8 anam,bflab
      character*2 bnam
      LOGICAL SOME,OUT
C     INTEGER SPIND
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_GRDPAR/NXYZG,NXYG,NXG,IGUNIT,UX(3),UY(3),UZ(3),
     1                  XYZG0(3)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_WFNOPT/WFNTYP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION AM(MXATOM),CM(3)
      DIMENSION MO(100),MOA(100),MOB(100)
      DIMENSION XYZMAP(3,100)
      DIMENSION XRANGE(3),YRANGE(3),ZRANGE(3),XVEC(3),YVEC(3),ZVEC(3)
      DIMENSION ORIGIN(3)
      DIMENSION IATOM(4)
      DIMENSION ERRMSG(3)
      EQUIVALENCE (MO(1),MOA(1))
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ELPMAP-'/
      DATA KEYEFC /4H EFC/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TENM3  /1.0D-03/
      DATA TEN12  /1.0D+12/
      DATA UHF    /'UHF     '/
      DATA MP2    /'MP2     '/
      DATA SCF    /'SCF     '/
      DATA MOA    /100*0/
      DATA MOB    /100*0/
      DATA XYZMAP /300*1.0D+12/
      DATA NPTXYZ /100/
      DATA IUNIT  /-1/
      DATA BOHR   /0.52917715D+00/
      DATA XRANGE /0.0D+00,0.0D+00,1.0D+00/
      DATA YRANGE /0.0D+00,0.0D+00,1.0D+00/
      DATA ZRANGE /0.0D+00,0.0D+00,1.0D+00/
      DATA XVEC   /1.0D+00,0.0D+00,0.0D+00/
      DATA YVEC   /0.0D+00,1.0D+00,0.0D+00/
      DATA ZVEC   /0.0D+00,0.0D+00,1.0D+00/
      DATA ORIGIN /0.0D+00,0.0D+00,0.0D+00/
      DATA IATOM  /0,0,0,0/
c%%%
C
*     NAMELIST /MAP/ IUNIT,XYZMAP,XRANGE,YRANGE,ZRANGE,XVEC,YVEC,ZVEC,
*    1               ORIGIN,IATOM,MO,MOA,MOB
C
*     REWIND IR
*     READ(IR,MAP,END=10,ERR=10)
* 10  CONTINUE
c%%%
C
C     ----- CALCULATE ELECTROSTATIC POTENTIAL -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
      if (ga_nodeid().eq.0) WRITE(IW,9994)
      SOME=.TRUE.
      OUT =.FALSE.     
C
C     ----- SET POINTERS FOR PARTITIONING OF MEMORY -----
C
      L1 = NUM
      L2 = (L1*(L1+1))/2
      L3 = L1*L1
C note: NPTMAX calculated based on zero grid points 
      NPTMAX = NAT +1
C
      I00 =          1
      I10 = I00    + NPTMAX*3
      I20 = I10    + NPTMAX
      I30 = I20    + (NUM*(NUM+1))/2
      I40 = I30    + (NUM*(NUM+1))/2
C     ----- MEMORY BLOCK FOR DA IN NWCHEM FORMAT -----
      I50 = I40+L3
C     ----- MOMORY BLOCK FOR DA IN HONDO  FORMAT-----
      I60 = I50+L3
      LAST= I60
C     NEED= LAST   - 1
      NEED= LAST
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED , 'MEM INIT:prp:hnd_elpmap:5',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_DIPOLE, MALLOC OF INIT  FAILED',911)

      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
C     ----- SET POINTERS FOR PARTITIONING OF MEMORY -----
C
      I00 =          INIT
      I10 = I00    + NPTMAX*3
      I20 = I10    + NPTMAX
      I30 = I20    + (NUM*(NUM+1))/2
      I40 = I30    + (NUM*(NUM+1))/2
C     ----- MEMORY BLOCK FOR DA IN NWCHEM FORMAT -----
      I50 = I40+L3
C     ----- MOMORY BLOCK FOR DA IN HONDO  FORMAT-----
      I60 = I50+L3
C
C     ----- DENSITY MATRIX FOR X(I20)=DA,X(I30)=DB -----
C

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i20),x(i30),x(i40),x(i50),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
C
C     ----- DEFINE POINTS FOR CALCULATION -----
C           1. GRID POINTS
C           2. SPECIFIC POINTS -XYZMAP-
C           3. NUCLEI
C           4. CENTER OF MASS
C
c%%
      NPTMAX=0
      NXGRD =0
      NYGRD =0
      NZGRD =0
*     CALL HND_GRDMAP(X(INIT+3*NPTMAX),NPTGRD,NXGRD,NYGRD,NZGRD,
*    1            XRANGE,YRANGE,ZRANGE,XVEC,YVEC,ZVEC,
*    2            ORIGIN,IATOM,IUNIT)

      NPTGRD=0
      NPTMAX=NPTMAX+NPTGRD
c%%
C
      IF(NPTGRD.GT.0.and.ga_nodeid().eq.0) THEN
         WRITE(IW,9985) NXGRD,NYGRD,NZGRD,
     1    IATOM,ORIGIN,XVEC,XRANGE(1),XRANGE(2),XRANGE(3),
     2                 YVEC,YRANGE(1),YRANGE(2),YRANGE(3),
     3                 ZVEC,ZRANGE(1),ZRANGE(2),ZRANGE(3)
      ELSE
         if (ga_nodeid().eq.0) WRITE(IW,9986)
      ENDIF
C
      FAC=ONE
      IF(IUNIT.GT.0) FAC=BOHR
      NPT=0
      DO 20 IPT=1,NPTXYZ
      IF(XYZMAP(1,IPT).EQ.TEN12.OR.XYZMAP(2,IPT).EQ.TEN12.OR.
     1   XYZMAP(3,IPT).EQ.TEN12) GO TO 20
      NPT=NPT+1
      X(INIT  +3*(NPT-1+NPTMAX))=XYZMAP(1,IPT)/FAC
      X(INIT+1+3*(NPT-1+NPTMAX))=XYZMAP(2,IPT)/FAC
      X(INIT+2+3*(NPT-1+NPTMAX))=XYZMAP(3,IPT)/FAC

   20 CONTINUE

      NPTMAX=NPTMAX+NPT

C
      DO 30 IAT=1,NAT
      X(INIT  +3*(IAT-1+NPTMAX))=C(1,IAT)
      X(INIT+1+3*(IAT-1+NPTMAX))=C(2,IAT)
      X(INIT+2+3*(IAT-1+NPTMAX))=C(3,IAT)

   30 CONTINUE
      NPTMAX=NPTMAX+NAT
C
      NCALL=1
      CALL HND_DIPAMS(AM,NCALL,NCODE,OUT)
      XCM=ZERO
      YCM=ZERO
      ZCM=ZERO
      TCM=ZERO
      DO 40 IAT=1,NAT
      XCM=XCM+AM(IAT)*C(1,IAT)
      YCM=YCM+AM(IAT)*C(2,IAT)
      ZCM=ZCM+AM(IAT)*C(3,IAT)
   40 TCM=TCM+AM(IAT)
      CM(1)=XCM/TCM
      CM(2)=YCM/TCM
      CM(3)=ZCM/TCM
      X(INIT  +3*(1-1+NPTMAX))=CM(1)
      X(INIT+1+3*(1-1+NPTMAX))=CM(2)
      X(INIT+2+3*(1-1+NPTMAX))=CM(3)
      NPTMAX=NPTMAX+1
C
C
C     ----- GET TOTAL DENSITY MATRIX -----
C
C     CALL HND_DAREAD(IDAF,IODA,X(I20),L2,16)
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
C        CALL HND_DAREAD(IDAF,IODA,X(I30),L2,20)
         DO 120 I=1,L2
            X(I+I20-1)=X(I+I20-1)+X(I+I30-1)
  120    CONTINUE
      ENDIF

C
C     ----- CALCULATE ELECTRONIC CONTRIBUTION AT ALL POINTS -----
C
      CALL HND_ELPOTC(NPTMAX,X(I00),X(I10),X(I20))
C
C     ----- OUTPUT RESULTS OF ALL POINTS -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9997)
C
*     WRITE(IP,9991)
*     WRITE(IP,9990) XRANGE(1),XRANGE(2),XRANGE(3),NXGRD
*     WRITE(IP,9990) YRANGE(1),YRANGE(2),YRANGE(3),NYGRD
*     WRITE(IP,9990) ZRANGE(1),ZRANGE(2),ZRANGE(3),NZGRD
C
      DO 230  IPT=1,NPTMAX
         XP = X(1+3*(IPT-1)+I00-1)
         YP = X(2+3*(IPT-1)+I00-1)
         ZP = X(3+3*(IPT-1)+I00-1)
C
C     ----- FORM NUCLEAR CONTRIBUTION -----
C
         ELPOTN = ZERO
         DO 210 I = 1,NAT
            XN = C(1,I) - XP
            YN = C(2,I) - YP
            ZN = C(3,I) - ZP
            RR =  SQRT(XN*XN + YN*YN + ZN*ZN)
            IF(ga_nodeid().eq.0.and.
     1         (RR.LT.TENM3).AND.(IPT.LT.(NPTMAX-1-NAT).OR.
     1                          IPT.GT.(NPTMAX-1    )    ) )
     2      WRITE(IW,9996) XP,YP,ZP,I
            IF(RR.LT.TENM3) GO TO 210
               ELPOTN = ELPOTN + ZAN(I)/RR
  210    CONTINUE
C
C     ----- FORM -EFC- CONTRIBUTION -----
C
      IF(IEFC.EQ.KEYEFC) THEN
         DO 220 I = 1,NEFC
            XN = EFCC(1,I) - XP
            YN = EFCC(2,I) - YP
            ZN = EFCC(3,I) - ZP
            RR =  SQRT(XN*XN + YN*YN + ZN*ZN)
            IF(RR.LT.TENM3.and.ga_nodeid().eq.0)
     &         WRITE(IW,9993) XP,YP,ZP,I
            IF(RR.LT.TENM3) GO TO 220
               ELPOTN = ELPOTN + EFCZ(I)/RR
  220    CONTINUE
      ENDIF
C
         ELPOTE = X(IPT+I10-1)
         ELPOTT = - ELPOTE + ELPOTN
         X(IPT+I10-1)=ELPOTT
         IF(ga_nodeid().eq.0.and.IPT.GT.NPTGRD.OR.OUT)
     1   WRITE(IW,9995) IPT,XP,YP,ZP,ELPOTT,ELPOTE
  230    CONTINUE
c--   IF(NPTGRD.GT.0.and.ga_nodeid().eq.0)
c--  &   WRITE(IP,9989) (X(IPT+I10-1),IPT=1,NPTGRD)
c--   if(ga_nodeid().eq.0) WRITE(IP,9988)
C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_DIPOLE, MA_POP_STACK OF INIT
     & FAILED',911)
C
      RETURN
 9999 FORMAT(/,10X,45(1H-),
     1       /,10X,'Electrostatic potential/diamagnetic shielding',
     2       /,10X,45(1H-),/)
 9998 FORMAT(' Not enough core in -elpmap-')
 9997 FORMAT(3X,'Point',6X,'X',9X,'Y',9X,'Z',5X,'Potential(a.u.)',
     1       3X,'Diamagnetic shielding(a.u.)')
 9996 FORMAT(' --- Warning - electrostatic potential at ',
     1 3F10.5,' . contribution from nucleus ',I3,' ignored')
 9995 FORMAT(1X,I5,2X,3F10.5,F15.6,6X,F15.6)
 9994 FORMAT(' 1 a.u. = 9.07618 esu/cm ( or statvolts ) ')
 9993 FORMAT(' --- Warning - electrostatic potential at ',
     1 3f10.5,' . contribution from  -efc-  ',i3,' ignored')
 9991 FORMAT(' $elpmap')
 9990 FORMAT(3F13.6,I5)
 9989 FORMAT(6F13.6,2X)
 9988 FORMAT(' $end   ')
 9986 FORMAT(/,' No grid defined. ')
 9985 FORMAT(/,1X,'( ',I4,' BY ',I4,' BY ',I4,' ) Grid defined.',
     1 /,' Iatom  = ',4I4,
     2 /,' Origin = ',3F10.5,
     3 /,' XVEC   = ',3F10.5,' XMIN=',F8.2,' XMAX=',F8.2,' XSTEP=',F8.2,
     4 /,' YVEC   = ',3F10.5,' YMIN=',F8.2,' YMAX=',F8.2,' YSTEP=',F8.2,
     5 /,' ZVEC   = ',3F10.5,' ZMIN=',F8.2,' ZMAX=',F8.2,' ZSTEP=',F8.2,
     6 /)
      END
      SUBROUTINE HND_ELPOTC(NPT,XYZPT,ELPOT,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     &                  NI,NJ,CCX,CCY,CCZ
      DIMENSION ELPOT(1),XYZPT(3,1),DAB(1)
      DIMENSION DIJ(225)
      DIMENSION XV(5,5,5),YV(5,5,5),ZV(5,5,5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA PI212 /1.1283791670955D+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT=NPRINT.EQ.6
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      DO 10 IPT=1,NPT
   10 ELPOT(IPT)=ZERO
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT-2)/2+1
C
      IJ=0
      DO 100 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 100 J=MINJ,JMAX
      IJ=IJ+1
  100 CONTINUE
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 DIJ(IJ)=DUM2
C
C     ----- ELECTROSTATIC POTENTIAL -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 IPT=1,NPT
      ZNUC=ONE
      CX=XYZPT(1,IPT)
      CY=XYZPT(2,IPT)
      CZ=XYZPT(3,IPT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      WW=W(IROOT)*ZNUC
      TT=ONE/(AA+UU)
      TX= SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LIT
      NI=I
      CALL HND_SXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT*WW
  410 CONTINUE
  420 CONTINUE
C
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 440 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUM=ZERO
      DO 430 IROOT=1,NROOTS
  430 DUM=DUM+XV(IX,JX,IROOT)*YV(IY,JY,IROOT)*ZV(IZ,JZ,IROOT)
      DUM=DUM*(AA1*PI212)
      IJ=IJ+1
      DUM=DUM*DIJ(IJ)
      LOCII=LOCI+I
      LOCJJ=LOCJ+J
      IJDEN=(LOCII*(LOCII-1))/2+LOCJJ
      DEN=DAB(IJDEN)
      IF(.NOT.IANDJ.OR.J.NE.I) DEN=DEN+DEN
      ELPOT(IPT)=ELPOT(IPT)+DUM*DEN
  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
 8000 CONTINUE
 9000 CONTINUE
C
      RETURN
      END
      SUBROUTINE HND_ELFMAP(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      parameter  (mxefc=3*7*7*7)
      integer rtdb,basis,geom
      logical status
      PARAMETER (MXIODA=255)
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL SOME,OUT
C     INTEGER SPIND
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_WFNOPT/WFNTYP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DATA KEYEFC /4H EFC/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TENM3  /1.0D-03/
      DATA UHF    /'UHF     '/
      DATA MP2    /'MP2     '/
      DATA SCF    /'SCF     '/
C
C     ----- CALCULATE ELECTRIC FIELD -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
      if (ga_nodeid().eq.0) WRITE(IW,9994)
      SOME=.TRUE.
      OUT =.FALSE.       
C
C     ----- SET BLOCK FOR PARTITIONING OF MEMORY -----
C
      NPTMAX = NAT
      L1     = NUM
      L2     = (L1*(L1+1))/2
      L3     = L1*L1
      I00    =          1
      I10    = I00    + NPTMAX*3
      I20    = I10    + NPTMAX*3
      I30    = I20    + (NUM*(NUM+1))/2
      I40    = I30    + (NUM*(NUM+1))/2
      I50    = I40    + L3
      I60    = I50    + L3
      LAST   = I60
      NEED   = LAST   - 1
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED , 'MEM INIT:prp:hnd_elfmap:6',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_DIPOLE, MALLOC OF INIT  FAILED',911)

      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
C     ----- SET POINTERS FOR PARTITIONING OF MEMORY -----
C
      I00    =          INIT
      I10    = I00    + NPTMAX*3
      I20    = I10    + NPTMAX*3
      I30    = I20    + (NUM*(NUM+1))/2
      I40    = I30    + (NUM*(NUM+1))/2
      I50    = I40    + L3
      I60    = I50    + L3
C
C     ----- DEFINE POINTS FOR CALCULATION -----
C           1. NUCLEI
C
      NPTMAX=0
C
      NPTMAX=NPTMAX+NAT
      DO 30 IAT=1,NAT
         X(INIT  +3*(IAT-1))=C(1,IAT)
         X(INIT+1+3*(IAT-1))=C(2,IAT)
         X(INIT+2+3*(IAT-1))=C(3,IAT)
   30 CONTINUE
C
C     ----- DENSITY MATRIX FOR X(I20)=DA,X(I30)=DB -----
C

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i20),x(i30),x(i40),x(i50),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)

C
C     ----- GET TOTAL DENSITY MATRIX -----
C
C     CALL HND_DAREAD(IDAF,IODA,X(I20),L2,16)

      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
C        CALL HND_DAREAD(IDAF,IODA,X(I30),L2,20)
         DO 120 I=1,L2
            X(I+I20-1)=X(I+I20-1)+X(I+I30-1)
  120    CONTINUE
      ENDIF
C
C     ----- CALCULATE ELECTRIC FIELD AT ALL POINTS -----
C
      CALL HND_ELFLDC(NPTMAX,X(I00),X(I10),X(I20))
C
C     ----- OUTPUT RESULTS OF ALL POINTS -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9997)
C
      DO 230  IPT=1,NPTMAX
         XP = X(1+3*(IPT-1)+I00-1)
         YP = X(2+3*(IPT-1)+I00-1)
         ZP = X(3+3*(IPT-1)+I00-1)
C
C     ----- FORM NUCLEAR CONTRIBUTION -----
C
         ELFLDX = ZERO
         ELFLDY = ZERO
         ELFLDZ = ZERO
         DO 210 I = 1,NAT
            XN = C(1,I) - XP
            YN = C(2,I) - YP
            ZN = C(3,I) - ZP
            RR =  SQRT(XN*XN + YN*YN + ZN*ZN)
            IF(RR.LT.TENM3) GO TO 210
               RR3=RR*RR*RR
               ELFLDX = ELFLDX - ZAN(I)*XN/RR3
               ELFLDY = ELFLDY - ZAN(I)*YN/RR3
               ELFLDZ = ELFLDZ - ZAN(I)*ZN/RR3
  210    CONTINUE
C
C     ----- FORM -EFC- CONTRIBUTION -----
C
      IF(IEFC.EQ.KEYEFC) THEN
         DO 220 I = 1,NEFC
            XN = EFCC(1,I) - XP
            YN = EFCC(2,I) - YP
            ZN = EFCC(3,I) - ZP
            RR =  SQRT(XN*XN + YN*YN + ZN*ZN)
            IF(RR.LT.TENM3.and.ga_nodeid().eq.0)
     &         WRITE(IW,9993) XP,YP,ZP,I
            IF(RR.LT.TENM3) GO TO 220
               RR3=RR*RR*RR
               ELFLDX = ELFLDX - EFCZ(I)*XN/RR3
               ELFLDY = ELFLDY - EFCZ(I)*YN/RR3
               ELFLDZ = ELFLDZ - EFCZ(I)*ZN/RR3
  220    CONTINUE
      ENDIF
C
         ELFLDX  = X(1+3*(IPT-1)+I10-1) + ELFLDX
         ELFLDY  = X(2+3*(IPT-1)+I10-1) + ELFLDY
         ELFLDZ  = X(3+3*(IPT-1)+I10-1) + ELFLDZ
         ELFLD   = SQRT(ELFLDX*ELFLDX + ELFLDY*ELFLDY + ELFLDZ*ELFLDZ)
         if (ga_nodeid().eq.0) WRITE(IW,9995) 
     &      IPT,XP,YP,ZP,ELFLDX,ELFLDY,ELFLDZ,ELFLD
  230    CONTINUE

C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_DIPOLE, MA_POP_STACK OF INIT
     & FAILED',911)
C
      RETURN
 9999 FORMAT(/,10X,14(1H-),/,10X,'Electric field',
     1       /,10X,14(1H-),/)
 9998 FORMAT(' Not enough core in -ELFMAP-')
 9997 FORMAT(3X,' Atom',6X,'X',9X,'Y',9X,'Z',8X,1X,
     1 15X,'Electric field (a.u.)',/,
     2 39X,7X,'X',14X,'Y',14X,'Z',11X,'Field',/,
     3 2X,96(1H-))
 9996 FORMAT(' --- Warning - electric field at ',
     1 3F10.5,' . contribution from nucleus ',I3,' ignored')
 9995 FORMAT(1X,I5,2X,3F10.5,4F15.6)
 9994 FORMAT(' 1 a.u. = 0.171524 10**(-8) dyn/esu ',/)
 9993 FORMAT(' --- warning - electric field at ',
     1 3F10.5,' . contribution from  -EFC-  ',I3,' ignored')
      END
      SUBROUTINE HND_ELFLDC(NPT,XYZPT,ELFLD,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_rys.fh"
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL SOME,OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ
     1                              ,CX,CY,CZ
      DIMENSION ELFLD(3,1),XYZPT(3,1),DAB(1)
      DIMENSION DIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION  XV(5,5,5), YV(5,5,5), ZV(5,5,5)
      DIMENSION DXV(5,5,5),DYV(5,5,5),DZV(5,5,5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ELFLDC-'/
      DATA MAXRYS /5/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      TOL=RLN10*ITOL
      OUT =NPRINT.EQ.6
      SOME=NPRINT.NE.-5
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      DO 10 IPT=1,NPT
      DO 10 I=1,3
   10 ELFLD(I,IPT)=ZERO
C
      NDER=1
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      LOCII=LOCI+I
      LOCJJ=LOCJ+J
      NN=(LOCII*(LOCII-1))/2+LOCJJ
      DEN=DAB(NN)
      IF(.NOT.IANDJ.OR.J.NE.I) DEN=DEN+DEN
      IJ=IJ+1
  360 DIJ(IJ)=DUM2*DEN
C
C     ----- ELECTRIC FIELD TERM -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 IPT=1,NPT
      ZNUC=ONE
      CX=XYZPT(1,IPT)
      CY=XYZPT(2,IPT)
      CZ=XYZPT(3,IPT)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
      UU=U(IROOT)*AA
      WW=W(IROOT)*ZNUC
      WW=WW*(UU+UU)
      TT=ONE/(AA+UU)
      TX= SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LIT
      NI=I
      CALL HND_DSXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT*WW
      CALL HND_DVXYZ
      DXV(I,J,IROOT)=XINT
      DYV(I,J,IROOT)=YINT
      DZV(I,J,IROOT)=ZINT*WW
  410 CONTINUE
  420 CONTINUE
C
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 440 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMX=DUMX+DXV(IX,JX,IROOT)* YV(IY,JY,IROOT)* ZV(IZ,JZ,IROOT)
      DUMY=DUMY+ XV(IX,JX,IROOT)*DYV(IY,JY,IROOT)* ZV(IZ,JZ,IROOT)
      DUMZ=DUMZ+ XV(IX,JX,IROOT)* YV(IY,JY,IROOT)*DZV(IZ,JZ,IROOT)
  430 CONTINUE
      IJ=IJ+1
      DUM=PI212*AA1*DIJ(IJ)
      ELFLD(1,IPT)=ELFLD(1,IPT)+DUMX*DUM
      ELFLD(2,IPT)=ELFLD(2,IPT)+DUMY*DUM
      ELFLD(3,IPT)=ELFLD(3,IPT)+DUMZ*DUM
  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
 8000 CONTINUE
 9000 CONTINUE
C
      RETURN
 9997 FORMAT(' IN -ELFLDC- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
      END
      SUBROUTINE HND_EFGMAP(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      parameter  (mxefc=3*7*7*7)
      integer rtdb,basis,geom
      logical status
      PARAMETER (MXIODA=255)
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL SOME,OUT
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_WFNOPT/WFNTYP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION A(6),VEC(3,3),EIG(3)
      DATA KEYEFC /4H EFC/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA THREE  /3.0D+00/
      DATA TENM3  /1.0D-03/
      DATA ANGLE  /180.0D+00/
      DATA UHF    /'UHF     '/
      DATA MP2    /'MP2     '/
      DATA SCF    /'SCF     '/
C
C     ----- CALCULATE ELECTRIC FIELD GRADIENT -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
      if (ga_nodeid().eq.0) WRITE(IW,9994)
      SOME=.TRUE.
      OUT =.FALSE.
C
      PI  = ACOS(-ONE)
      DEG = ANGLE/PI
C
C     ----- SET BLOCK FOR PARTITIONING OF MEMORY -----
C
      NPTMAX = NAT
      L1     = NUM
      L2     = (L1*(L1+1))/2
      L3     = L1*L1
      I00    =          1
      I10    = I00    + NPTMAX*3
      I20    = I10    + NPTMAX*6
      I30    = I20    + (NUM*(NUM+1))/2
      I40    = I30    + (NUM*(NUM+1))/2
      I50    = I40    + L3
      I60    = I50    + L3
      LAST   = I60
      NEED   = LAST   - 1
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED , 'MEM INIT:prp:hnd_efgmap:7',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_DIPOLE, MALLOC OF INIT  FAILED',911)

      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
C     ----- SET POINTERS FOR PARTITIONING OF MEMORY -----
C
      I00    =          INIT
      I10    = I00    + NPTMAX*3
      I20    = I10    + NPTMAX*6
      I30    = I20    + (NUM*(NUM+1))/2
      I40    = I30    + (NUM*(NUM+1))/2
      I50    = I40    + L3
      I60    = I50    + L3
C
C     ----- DEFINE POINTS FOR CALCULATION -----
C           1. NUCLEI
C
      NPTMAX=0
      NPTMAX=NPTMAX+NAT
      DO 30 IAT=1,NAT
         X(INIT  +3*(IAT-1))=C(1,IAT)
         X(INIT+1+3*(IAT-1))=C(2,IAT)
         X(INIT+2+3*(IAT-1))=C(3,IAT)
   30 CONTINUE
C
C     ----- DENSITY MATRIX FOR X(I20)=DA,X(I30)=DB(0?) -----
C

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i20),x(i30),x(i40),x(i50),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)

C
C     ----- GET TOTAL DENSITY MATRIX -----
C
C     CALL HND_DAREAD(IDAF,IODA,X(I20),L2,16)

      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
         CALL HND_DAREAD(IDAF,IODA,X(I30),L2,20)
         DO 120 I=1,L2
            X(I+I20-1)=X(I+I20-1)+X(I+I30-1)
  120    CONTINUE
      ENDIF
C
C     ----- CALCULATE ELECTRIC FIELD GRADIENT AT ALL POINTS -----
C
      CALL HND_ELFGRD(NPTMAX,X(I00),X(I10),X(I20))
C
C     ----- OUTPUT RESULTS OF ALL POINTS -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9997)
C
      DO 230  IPT=1,NPTMAX
         XP = X(1+3*(IPT-1)+I00-1)
         YP = X(2+3*(IPT-1)+I00-1)
         ZP = X(3+3*(IPT-1)+I00-1)
C
C     ----- FORM NUCLEAR CONTRIBUTION -----
C
         EFGRXX = ZERO
         EFGRYY = ZERO
         EFGRZZ = ZERO
         EFGRXY = ZERO
         EFGRXZ = ZERO
         EFGRYZ = ZERO
         DO 210 I = 1,NAT
            XN = C(1,I) - XP
            YN = C(2,I) - YP
            ZN = C(3,I) - ZP
            RR =  SQRT(XN*XN + YN*YN + ZN*ZN)
            IF(RR.LT.TENM3) GO TO 210
               RR5=RR*RR*RR*RR*RR
               EFGRXX = EFGRXX - ZAN(I)*XN*XN/RR5
               EFGRYY = EFGRYY - ZAN(I)*YN*YN/RR5
               EFGRZZ = EFGRZZ - ZAN(I)*ZN*ZN/RR5
               EFGRXY = EFGRXY - ZAN(I)*XN*YN/RR5
               EFGRXZ = EFGRXZ - ZAN(I)*XN*ZN/RR5
               EFGRYZ = EFGRYZ - ZAN(I)*YN*ZN/RR5
  210    CONTINUE
C
C     ----- FORM -EFC- CONTRIBUTION -----
C
      IF(IEFC.EQ.KEYEFC) THEN
         DO 220 I = 1,NEFC
            XN = EFCC(1,I) - XP
            YN = EFCC(2,I) - YP
            ZN = EFCC(3,I) - ZP
            RR =  SQRT(XN*XN + YN*YN + ZN*ZN)
            IF(RR.LT.TENM3.and.ga_nodeid().eq.0) 
     &        WRITE(IW,9993) XP,YP,ZP,I
            IF(RR.LT.TENM3) GO TO 220
               RR5=RR*RR*RR*RR*RR
               EFGRXX = EFGRXX - EFCZ(I)*XN*XN/RR5
               EFGRYY = EFGRYY - EFCZ(I)*YN*YN/RR5
               EFGRZZ = EFGRZZ - EFCZ(I)*ZN*ZN/RR5
               EFGRXY = EFGRXY - EFCZ(I)*XN*YN/RR5
               EFGRXZ = EFGRXZ - EFCZ(I)*XN*ZN/RR5
               EFGRYZ = EFGRYZ - EFCZ(I)*YN*ZN/RR5
  220    CONTINUE
      ENDIF
C
       EFGRXX  = X(1+6*(IPT-1)+I10-1)/THREE + EFGRXX
       EFGRYY  = X(2+6*(IPT-1)+I10-1)/THREE + EFGRYY
       EFGRZZ  = X(3+6*(IPT-1)+I10-1)/THREE + EFGRZZ
       EFGRXY  = X(4+6*(IPT-1)+I10-1)/THREE + EFGRXY
       EFGRXZ  = X(5+6*(IPT-1)+I10-1)/THREE + EFGRXZ
       EFGRYZ  = X(6+6*(IPT-1)+I10-1)/THREE + EFGRYZ
C
       EFGXX  =   TWO*EFGRXX - EFGRYY - EFGRZZ
       EFGYY  =   TWO*EFGRYY - EFGRXX - EFGRZZ
       EFGZZ  =   TWO*EFGRZZ - EFGRXX - EFGRYY
       EFGXY  = THREE*EFGRXY
       EFGXZ  = THREE*EFGRXZ
       EFGYZ  = THREE*EFGRYZ
C
       if (ga_nodeid().eq.0) WRITE(IW,9995)
     &    IPT,XP,YP,ZP,EFGXX,EFGYY,EFGZZ,EFGXY,EFGXZ,EFGYZ
C
       X(1+6*(IPT-1)+I10-1) = EFGXX
       X(2+6*(IPT-1)+I10-1) = EFGYY
       X(3+6*(IPT-1)+I10-1) = EFGZZ
       X(4+6*(IPT-1)+I10-1) = EFGXY
       X(5+6*(IPT-1)+I10-1) = EFGXZ
       X(6+6*(IPT-1)+I10-1) = EFGYZ
C
  230  CONTINUE
C
      if (ga_nodeid().eq.0) WRITE(IW,9992)
C
      DO 240 IPT=1,NPTMAX
      A(1) = X(1+6*(IPT-1)+I10-1)
      A(2) = X(4+6*(IPT-1)+I10-1)
      A(3) = X(2+6*(IPT-1)+I10-1)
      A(4) = X(5+6*(IPT-1)+I10-1)
      A(5) = X(6+6*(IPT-1)+I10-1)
      A(6) = X(3+6*(IPT-1)+I10-1)
      CALL HND_DIAEFG(A,VEC,EIG,3)
      ETA  = ABS( (EIG(1)-EIG(2)) / EIG(3) )
      if (ga_nodeid().eq.0) WRITE(IW,9991)
     &   IPT,EIG(1),EIG(2),EIG(3),ETA
  240 CONTINUE
C
      if (ga_nodeid().eq.0) WRITE(IW,9990)
C
      DO 260 IPT=1,NPTMAX
      A(1) = X(1+6*(IPT-1)+I10-1)
      A(2) = X(4+6*(IPT-1)+I10-1)
      A(3) = X(2+6*(IPT-1)+I10-1)
      A(4) = X(5+6*(IPT-1)+I10-1)
      A(5) = X(6+6*(IPT-1)+I10-1)
      A(6) = X(3+6*(IPT-1)+I10-1)
      CALL HND_DIAEFG(A,VEC,EIG,3)
      DO 250 J=1,3
      DO 250 I=1,3
  250 VEC(I,J)= ACOS(VEC(I,J))*DEG
      if (ga_nodeid().eq.0) WRITE(IW,9989) IPT
      if (ga_nodeid().eq.0) WRITE(IW,9988) ((VEC(I,J),J=1,3),I=1,3)
  260 CONTINUE
C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_DIPOLE, MA_POP_STACK OF INIT
     & FAILED',911)
C
      RETURN
 9999 FORMAT(/,10X,23(1H-),/,10X,'Electric field gradient',
     1       /,10X,23(1H-),/)
 9998 FORMAT(' Not enough core in -efgmap-')
 9997 FORMAT(1X,' Atom',6X,'X',9X,'Y',9X,'Z',8X,1X,
     1 19X,'Electric field gradient in molecular frame (a.u.)',/,
     2 39X,6X,'XX',13X,'YY',13X,'ZZ',13X,'XY',13X,'XZ',13X,'YZ',/,
     3 2X,126(1H-))
 9996 FORMAT(' --- Warning - electric field gradient at ',
     1 3F10.5,' . contribution from nucleus ',I3,' ignored')
 9995 FORMAT(1X,I5,2X,3F10.5,6F15.6)
 9994 FORMAT(' 1 a.u. = 0.324123 10**(16) esu/cm**3 ',
     1       ' ( or statvolts/cm**2 )',' = 0.97174 10**(22) v/m**2 ',/)
 9993 FORMAT(' --- Warning - electric field gradient at ',
     1 3F10.5,' . contribution from  -efc-  ',I3,' ignored')
 9992 FORMAT(/,1X,' Atom',3x,'Principal components of the -EFG- ',
     1       'tensor (a.u.)',9x,'Asymmetry parameter eta',/,
     2       2X,86(1H-))
 9991 FORMAT(1X,I5,3F15.6,14X,F15.6)
 9990 FORMAT(/,' Orientation ( in degrees ) of the principal axis',
     1         ' of -EFG- tensor w.r.t absolute (molecular) frame')
 9989 FORMAT(/,' Atom = ',I5,/,' ------ ')
 9988 FORMAT(1X,3F15.6)
      END
      SUBROUTINE HND_DIAEFG(A,VEC,EIG,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(6),H(3,3),VEC(3,3),EIG(3),BIG(6)
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA CONV   /1.0D-05/
      DATA MAXIT  /50/
C
      DO 20 I=1,3
         DO 10 J=1,3
   10       VEC(J,I)=ZERO
   20    VEC(I,I)=ONE
      IJ=0
      DO 25 I=1,3
      DO 25 J=1,I
         IJ=IJ+1
         H(I,J)=A(IJ)
   25    H(J,I)=A(IJ)
      CALL HND_DIAJAC(H,VEC,EIG,NDIM,NDIM,BIG)
C
C     ----- PUT THE PRINCIPAL COMPONENTS OF THE TENSOR AND -----
C           THE CORRESPONDING EIGENVECTORS IN INCREASING
C           ORDER ACCORDING TO THE ABSOLUTE VALUE OF THE
C           EIGENVALUE ( THE PRINCIPAL COMPONENT ) .
C
      DO 50 I=1,3
         JJ=I
         DO 30 J=I,3
            IF( ABS(EIG(J)).LT. ABS(EIG(JJ))) JJ=J
   30       CONTINUE
         IF(JJ.EQ.I) GO TO 50
C
         XX=EIG(JJ)
         EIG(JJ)=EIG(I)
         EIG(I)=XX
         DO 40 J=1,3
            XX=VEC(J,JJ)
            VEC(J,JJ)=VEC(J,I)
            VEC(J,I)=XX
   40       CONTINUE
C
   50    CONTINUE
      RETURN
      END
      SUBROUTINE HND_ELFGRD(NPT,XYZPT,ELFGR,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
#include "hnd_rys.fh"
#include "hnd_tol.fh"
      CHARACTER*8 ERRMSG
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL SOME,OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDER/XINT,YINT,ZINT,TX,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ
     1                              ,CX,CY,CZ
      DIMENSION W2(6),W4(6)
      DIMENSION XYZPT(3,1),ELFGR(6,1),DAB(1)
      DIMENSION DIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION   XV(5,5,6),  YV(5,5,6),  ZV(5,5,6)
      DIMENSION  DXV(5,5,6), DYV(5,5,6), DZV(5,5,6)
      DIMENSION DDXV(5,5,6),DDYV(5,5,6),DDZV(5,5,6)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ELFGRD-'/
      DATA MAXRYS /6/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA FOUR   /4.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      TOL=RLN10*ITOL
      OUT =NPRINT.EQ.6
      SOME=NPRINT.NE.-5
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      DO 10 IPT=1,NPT
      DO 10 I=1,6
   10 ELFGR(I,IPT)=ZERO
C
      NDER=2
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      MROOTS=(LIT+LJT+NDER-2)/2+1
      IF(MROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,MROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      LOCII=LOCI+I
      LOCJJ=LOCJ+J
      NN=(LOCII*(LOCII-1))/2+LOCJJ
      DEN=DAB(NN)
      IF(.NOT.IANDJ.OR.J.NE.I) DEN=DEN+DEN
      IJ=IJ+1
  360 DIJ(IJ)=DUM2*DEN
C
C     ----- ELECTRIC FIELD GRADIENT TERM -----
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 IPT=1,NPT
      ZNUC=ONE
      CX=XYZPT(1,IPT)
      CY=XYZPT(2,IPT)
      CZ=XYZPT(3,IPT)
      YY=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      CALL HND_DROOT
      DO 420 IROOT=1,MROOTS
      UU=U9(IROOT)*AA
      U2=UU
      U4=UU*UU
      WW=W9(IROOT)*ZNUC
      W2(IROOT)=WW*U2*TWO
      W4(IROOT)=WW*U4*FOUR
      TT=ONE/(AA+UU)
      TX= SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LIT
      NI=I
      CALL HND_DSXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT
      CALL HND_DVXYZ
      DXV(I,J,IROOT)=XINT
      DYV(I,J,IROOT)=YINT
      DZV(I,J,IROOT)=ZINT
      CALL HND_DDVXYZ
      DDXV(I,J,IROOT)=XINT
      DDYV(I,J,IROOT)=YINT
      DDZV(I,J,IROOT)=ZINT
  410 CONTINUE
  420 CONTINUE
C
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 440 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DUMXX=ZERO
      DUMYY=ZERO
      DUMZZ=ZERO
      DUMXY=ZERO
      DUMXZ=ZERO
      DUMYZ=ZERO
      DO 430 IROOT=1,MROOTS
      DUMX =DUMX+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMY =DUMY+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMZ =DUMZ+
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W2(IROOT)
      DUM  =
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMXX=DUMXX-DUM+
     1   DDXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMYY=DUMYY-DUM+
     1     XV(IX,JX,IROOT)*DDYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMZZ=DUMZZ-DUM+
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*DDZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMXY=DUMXY+
     1    DXV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMXZ=DUMXZ+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
      DUMYZ=DUMYZ+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
  430 CONTINUE
      IJ=IJ+1
      DUM=PI212*AA1*DIJ(IJ)
      ELFGR(1,IPT)=ELFGR(1,IPT)+DUMXX*DUM
      ELFGR(2,IPT)=ELFGR(2,IPT)+DUMYY*DUM
      ELFGR(3,IPT)=ELFGR(3,IPT)+DUMZZ*DUM
      ELFGR(4,IPT)=ELFGR(4,IPT)+DUMXY*DUM
      ELFGR(5,IPT)=ELFGR(5,IPT)+DUMXZ*DUM
      ELFGR(6,IPT)=ELFGR(6,IPT)+DUMYZ*DUM
  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
 8000 CONTINUE
 9000 CONTINUE
C
      RETURN
 9997 FORMAT(' IN -ELFGRD- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS = ',3I3)
      END
      SUBROUTINE HND_ELDMAP(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXBFN=3072)
      integer  rtdb,geom,basis
      logical  status
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP, SCFTYP, RUNTYP
      CHARACTER*8 RHF, UHF, MP2, GVB, SCF, PROP
      LOGICAL SOME,OUT,DBUG
      INTEGER SPIND
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_GRDPAR/NXYZG,NXYG,NXG,IGUNIT,UX(3),UY(3),UZ(3),
     1       XYZG0(3)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_SCFWFN/CICOEF(2,12),F(25),ALPHA(325),BETA(325),
     1 NO(10),NCO,NSETO,NPAIR
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_RUNOPT/RUNTYP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION AM(MXATOM),CM(3)
      DIMENSION MO(100),MOA(100),MOB(100)
      DIMENSION XYZMAP(3,100)
      DIMENSION XRANGE(3),YRANGE(3),ZRANGE(3),XVEC(3),YVEC(3),ZVEC(3)
      DIMENSION ORIGIN(3)
      DIMENSION IATOM(4)
      DIMENSION ERRMSG(3)
      EQUIVALENCE (MO(1),MOA(1))
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ELDMAP-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA TENM3  /1.0D-03/
      DATA TEN12  /1.0D+12/
      DATA PROP   /'PROPTY  '/
      DATA RHF    /'RHF     '/
      DATA UHF    /'UHF     '/
      DATA MP2    /'MP2     '/
      DATA GVB    /'GVB     '/
      DATA SCF    /'SCF     '/
      DATA MOA    /100*0/
      DATA MOB    /100*0/
      DATA MAXMO  /100/
      DATA XYZMAP /300*1.0D+12/
      DATA NPTXYZ /100/
      DATA IUNIT  /-1/
      DATA BOHR   /0.52917715D+00/
      DATA XRANGE /0.0D+00,0.0D+00,1.0D+00/
      DATA YRANGE /0.0D+00,0.0D+00,1.0D+00/
      DATA ZRANGE /0.0D+00,0.0D+00,1.0D+00/
      DATA XVEC   /1.0D+00,0.0D+00,0.0D+00/
      DATA YVEC   /0.0D+00,1.0D+00,0.0D+00/
      DATA ZVEC   /0.0D+00,0.0D+00,1.0D+00/
      DATA ORIGIN /0.0D+00,0.0D+00,0.0D+00/
      DATA IATOM  /0,0,0,0/
C
*     NAMELIST /MAP/ IUNIT,XYZMAP,XRANGE,YRANGE,ZRANGE,XVEC,YVEC,ZVEC,
*    1               ORIGIN,IATOM,MO,MOA,MOB
C
*     REWIND IR
*     READ(IR,MAP,END=10,ERR=10)
* 10  CONTINUE
C
C     ----- CALCULATE ELECTRON OR SPIN DENSITY -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
      if (ga_nodeid().eq.0) WRITE(IW,9996)
      DBUG=.FALSE.
      OUT =.FALSE.       
      OUT =OUT.OR.DBUG
      SOME=.TRUE.
      SOME=SOME.OR.OUT
      IF(RUNTYP.EQ.PROP) THEN
         NA =0
         NCO=NA
      ENDIF
C
C     ----- DEFINE POINTS FOR CALCULATION -----
C           1. GRID POINTS
C           2. SPECIFIC POINTS -XYZMAP-
C           3. NUCLEI
C           4. CENTER OF MASS
C
      NPTMAX=0
      NPTGRD=0
      NXGRD =0
      NYGRD =0
      NZGRD =0
*     CALL HND_GRDMAP(X(1+3*NPTMAX),NPTGRD,NXGRD,NYGRD,NZGRD,
*    1            XRANGE,YRANGE,ZRANGE,XVEC,YVEC,ZVEC,
*    2            ORIGIN,IATOM,IUNIT)
      NPTMAX=NPTMAX+NPTGRD
C
      IF(NPTGRD.GT.0.and.ga_nodeid().eq.0) THEN
         WRITE(IW,9985) NXGRD,NYGRD,NZGRD,
     1    IATOM,ORIGIN,XVEC,XRANGE(1),XRANGE(2),XRANGE(3),
     2                 YVEC,YRANGE(1),YRANGE(2),YRANGE(3),
     3                 ZVEC,ZRANGE(1),ZRANGE(2),ZRANGE(3)
      ELSE
         if (ga_nodeid().eq.0) WRITE(IW,9986)
      ENDIF
C
      FAC=ONE
      IF(IUNIT.GT.0) FAC=BOHR
      NPT=0
      DO 20 IPT=1,NPTXYZ
      IF(XYZMAP(1,IPT).EQ.TEN12.OR.XYZMAP(2,IPT).EQ.TEN12.OR.
     1   XYZMAP(3,IPT).EQ.TEN12) GO TO 20
      NPT=NPT+1
*     X(INIT  +3*(NPT-1+NPTMAX))=XYZMAP(1,IPT)/FAC
*     X(INIT+1+3*(NPT-1+NPTMAX))=XYZMAP(2,IPT)/FAC
*     X(INIT+2+3*(NPT-1+NPTMAX))=XYZMAP(3,IPT)/FAC
   20 CONTINUE
      NPTMAX=NPTMAX+NPT
C
      NPTMAX=NPTMAX+NAT
C
      NCALL=1
      CALL HND_DIPAMS(AM,NCALL,NCODE,OUT)
      XCM=ZERO
      YCM=ZERO
      ZCM=ZERO
      TCM=ZERO
      DO 40 IAT=1,NAT
      XCM=XCM+AM(IAT)*C(1,IAT)
      YCM=YCM+AM(IAT)*C(2,IAT)
      ZCM=ZCM+AM(IAT)*C(3,IAT)
   40 TCM=TCM+AM(IAT)
      CM(1)=XCM/TCM
      CM(2)=YCM/TCM
      CM(3)=ZCM/TCM
C
      NPTMAX=NPTMAX+1
C
C     ----- we are going to use the general case -----
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM* NUM
      NEED= L2 + L2 + L3 + L3 + 4*NPTMAX
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED , 'MEM INIT:prp:hnd_eldmap:8',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_DIPOLE, MALLOC OF INIT  FAILED',911)

      DO I=1,NEED
         X(I+INIT-1)=ZERO
      ENDDO
C
      NPTMAX=0           
      DO 30 IAT=1,NAT
      X(INIT  +3*(IAT-1+NPTMAX))=C(1,IAT)
      X(INIT+1+3*(IAT-1+NPTMAX))=C(2,IAT)
      X(INIT+2+3*(IAT-1+NPTMAX))=C(3,IAT)
   30 CONTINUE
      NPTMAX=NPTMAX+NAT
      X(INIT  +3*(1-1+NPTMAX))=CM(1)
      X(INIT+1+3*(1-1+NPTMAX))=CM(2)
      X(INIT+2+3*(1-1+NPTMAX))=CM(3)
      NPTMAX=NPTMAX+1
C
      IF(DBUG) THEN
         write(6,*) 'need   = ',need
         write(6,*) 'runtyp = ',runtyp
         write(6,*) 'wfntyp = ',wfntyp
         write(6,*) 'scftyp = ',scftyp
      ENDIF
C
C--   IF(WFNTYP.EQ.SCF.AND.SCFTYP.NE.GVB) GO TO 200
C
C     ----- GENERAL CASE = ONE PARTICLE REDUCED DENSITY MATRIX -----
C           IS NOT DIAGONAL.
C
      IF(RUNTYP.EQ.PROP) THEN
         if (ga_nodeid().eq.0) WRITE(IW,9983)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(DBUG) THEN
         write(6,*) 'general case ...'
      ENDIF
C
C     ----- SET POINTERS FOR PARTITIONING OF MEMORY -----
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM* NUM
      I00 =          INIT
      I10 = I00    + NPTMAX*3
      I20 = I10    + NPTMAX
      I30 = I20    + (NUM*(NUM+1))/2
      I40 = I30    + (NUM*(NUM+1))/2
C     ----- MEMORY BLOCK FOR DA IN NWCHEM FORMAT -----
      I50 = I40    + L3
C     ----- MOMORY BLOCK FOR DA IN HONDO  FORMAT-----
      I60 = I50    + L3
C
C     ----- GET TOTAL DENSITY MATRIX -----
C
C     CALL HND_DAREAD(IDAF,IODA,X(I20),L2,16)

c     ----- calculate density matrix -----

      call hnd_prp_get_dens(rtdb,basis,geom,
     1                      x(i20),x(i30),x(i40),x(i50),
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)

      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
         CALL HND_DAREAD(IDAF,IODA,X(I30),L2,20)
         DO 120 I=1,L2
            DUMA      =X(I+I20-1)
            DUMB      =X(I+I30-1)
            X(I+I20-1)=DUMA+DUMB             
            X(I+I30-1)=DUMA-DUMB             
  120    CONTINUE
      ENDIF
C
C     ----- CALCULATE ELECTRONIC CONTRIBUTION AT ALL POINTS -----
C
      CALL HND_ELDENC(NPTMAX,X(I00),X(I10),X(I20))
C
C     ----- OUTPUT RESULTS FOR ALL POINTS -----
C
      if (ga_nodeid().eq.0) then
         WRITE(IW,9993)
         WRITE(IW,9997)
C
*        WRITE(IP,9991)
*        WRITE(IP,9990) XRANGE(1),XRANGE(2),XRANGE(3),NXGRD
*        WRITE(IP,9990) YRANGE(1),YRANGE(2),YRANGE(3),NYGRD
*        WRITE(IP,9990) ZRANGE(1),ZRANGE(2),ZRANGE(3),NZGRD
      endif
C
      DO 130  IPT=1,NPTMAX
         XP = X(1+3*(IPT-1)+I00-1)
         YP = X(2+3*(IPT-1)+I00-1)
         ZP = X(3+3*(IPT-1)+I00-1)
         ELDENE = X(IPT+I10-1)
         IF(ga_nodeid().eq.0.and.IPT.GT.NPTGRD.OR.OUT)
     &      WRITE(IW,9995) IPT,XP,YP,ZP,ELDENE
  130    CONTINUE
c--   IF(ga_nodeid().eq.0.and.NPTGRD.GT.0) 
c--  &     WRITE(IP,9989) (X(IPT+I10-1),IPT=1,NPTGRD)
c--   if (ga_nodeid().eq.0) WRITE(IP,9988)
C
C     ----- SPIN DENSITY NOW -----
C
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
         CALL HND_ELDENC(NPTMAX,X(I00),X(I10),X(I30))
         if (ga_nodeid().eq.0) then
            WRITE(IW,9992)
            WRITE(IW,9997)
         endif
         DO 140  IPT=1,NPTMAX
            XP = X(1+3*(IPT-1)+I00-1)
            YP = X(2+3*(IPT-1)+I00-1)
            ZP = X(3+3*(IPT-1)+I00-1)
            ELDENE = X(IPT+I10-1)
            IF(ga_nodeid().eq.0.and.IPT.GT.NPTGRD.OR.OUT)
     &         WRITE(IW,9995) IPT,XP,YP,ZP,ELDENE
  140    CONTINUE
      ENDIF
      GO TO 500
  200 CONTINUE
C
C     ----- SPECIAL CASE = ONE PARTICLE REDUCED DENSITY MATRIX -----
C           IS DIAGONAL.
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM* NUM
      I00 =    INIT
      I10 =I00+NPTMAX*3
      I11 =I10+NPTMAX
      I12 =I11+NPTMAX
      I20 =I12+NPTMAX
      I30 =I20+NUM*NUM
      I40  = I30    + NUM*NUM
C     ----- MEMORY BLOCK FOR C IN NWCHEM FORMAT -----
      I50  = I40+L3
C     ----- MOMORY BLOCK FOR E IN NECHEM FORMAT-----
      I60  = I50+L1
      NEED=LAST-1
C
C     ----- READ MOLECULAR ORBITALS -----
C
C*    CALL HND_DAREAD(IDAF,IODA,X(I20),L3,15)

c     ----- get mo vectors from file -----

c*      g_vecs=ga_create_atom_blocked (geom, basis, 'mo vectors')

c     ----- figure output input/output MO vectors -----

c*     if (.not.rtdb_cget(rtdb,'scf:input  vectors',1,movecs_in))
c*   $      movecs_in = 'atomic'
c*     if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
c*   $      movecs_out = ' '
c*     if (movecs.eq.' ') then
c*       if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore') then
c           movecs_out = 'movecs'
c*       else 
c*          movecs_out = movecs_in
c*       endif
c*     endif
c*     movecs = movecs_out

c      ----- read eig-vec x(i40); eng-val in x(i50) for tmp. -----

c*     status = movecs_read (movecs, 1, noc, x(i50), g_vecs)
c*     call ga_get(g_vecs, 1, num, 1, num, x(i40), num)
c*     status = ga_destroy ( g_vecs )
c%%%
c*     call hnd_nwhnd_tran(x(i40),x(i20),num)
       call hnd_get_eigens(rtdb,geom,basis,x(i20),x(i50),num)
c%%%
c
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
      CALL HND_DAREAD(IDAF,IODA,X(I30),L3,19)
      ENDIF
C
C     ----- INDIVIDUAL (ALPHA) ORBITAL DENSITY -----
C
      DO 230 IMO=1,MAXMO
      IF(MO(IMO).EQ.0) GO TO 230
      if (ga_nodeid().eq.0) then
         IF(MO(IMO).GT.0) WRITE(IW,9994) ABS(MO(IMO))
         IF(MO(IMO).GT.0) WRITE(IW,9997)
         IF(MO(IMO).LT.0) WRITE(IW,9982) ABS(MO(IMO))
         IF(MO(IMO).LT.0) WRITE(IW,9980)
C
*        WRITE(IP,9991)
*        WRITE(IP,9990) XRANGE(1),XRANGE(2),XRANGE(3),NXGRD
*        WRITE(IP,9990) YRANGE(1),YRANGE(2),YRANGE(3),NYGRD
*        WRITE(IP,9990) ZRANGE(1),ZRANGE(2),ZRANGE(3),NZGRD
      endif
C
      CALL HND_ELWFNC(NPTMAX,X(I00),X(I12),X(I20+
     1                L1*(ABS(MO(IMO))-1)))
C
      DO 220 IPT=1,NPTMAX
      IF(MO(IMO).GT.0) X(IPT+I12-1)=X(IPT+I12-1)**2
      XP=X(1+3*(IPT-1)+I00-1)
      YP=X(2+3*(IPT-1)+I00-1)
      ZP=X(3+3*(IPT-1)+I00-1)
      ELDENE=X(IPT+I12-1)
      IF(ga_nodeid().eq.0.and.IPT.GT.NPTGRD.OR.OUT)
     &   WRITE(IW,9995) IPT,XP,YP,ZP,ELDENE
  220 CONTINUE
c--   IF(ga_nodeid().eq.0.and.NPTGRD.GT.0) 
c--  &   WRITE(IP,9989) (X(IPT+I12-1),IPT=1,NPTGRD)
c--   if (ga_nodeid().eq.0) WRITE(IP,9988)
  230 CONTINUE
C
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) GO TO 300
C
C     ----- INDIVIDUAL (BETA) ORBITAL DENSITY -----
C
      DO 250 IMO=1,MAXMO
      if (ga_nodeid().eq.0) then
         IF(MOB(IMO).EQ.0) GO TO 250
         IF(MOB(IMO).GT.0) WRITE(IW,9984) ABS(MOB(IMO))
         IF(MOB(IMO).GT.0) WRITE(IW,9997)
         IF(MOB(IMO).LT.0) WRITE(IW,9981) ABS(MOB(IMO))
         IF(MOB(IMO).LT.0) WRITE(IW,9980)
C
*        WRITE(IP,9991)
*        WRITE(IP,9990) XRANGE(1),XRANGE(2),XRANGE(3),NXGRD
*        WRITE(IP,9990) YRANGE(1),YRANGE(2),YRANGE(3),NYGRD
*        WRITE(IP,9990) ZRANGE(1),ZRANGE(2),ZRANGE(3),NZGRD
      endif
C
      CALL HND_ELWFNC(NPTMAX,X(I00),X(I12),X(I30+
     1                L1*(ABS(MOB(IMO))-1)))
C
      DO 240 IPT=1,NPTMAX
      IF(MOB(IMO).GT.0) X(IPT+I12-1)=X(IPT+I12-1)**2
      XP=X(1+3*(IPT-1)+I00-1)
      YP=X(2+3*(IPT-1)+I00-1)
      ZP=X(3+3*(IPT-1)+I00-1)
      ELDENE=X(IPT+I12-1)
      IF(ga_nodeid().eq.0.and.IPT.GT.NPTGRD.OR.OUT) 
     &  WRITE(IW,9995) IPT,XP,YP,ZP,ELDENE
  240 CONTINUE
c--   IF(ga_nodeid().eq.0.and.NPTGRD.GT.0) 
c--  &  WRITE(IP,9989) (X(IPT+I12-1),IPT=1,NPTGRD)
c--   if (ga_nodeid().eq.0) WRITE(IP,9988)
  250 CONTINUE
C
  300 CONTINUE
C
C     ----- TOTAL ELECTRON DENSITY AND SPIN DENSITY -----
C
      DO 310 IPT=1,NPTMAX
      X(IPT+I10-1)=ZERO
  310 X(IPT+I11-1)=ZERO
C
      DO 330 IMO=1,NA
      CALL HND_ELWFNC(NPTMAX,X(I00),X(I12),X(I20+L1*(IMO-1)))
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.RHF).AND.IMO.LE.NCO) OCC=TWO
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.RHF).AND.IMO.GT.NCO) OCC=ONE
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)               ) OCC=ONE
      DO 320 IPT=1,NPTMAX
      X(IPT+I10-1)=X(IPT+I10-1)+OCC*X(IPT+I12-1)**2
  320 X(IPT+I11-1)=X(IPT+I11-1)+OCC*X(IPT+I12-1)**2
  330 CONTINUE
C
      IF((WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1   (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)     ) GO TO 360
      DO 350 IMO=1,NB
      CALL HND_ELWFNC(NPTMAX,X(I00),X(I12),X(I30+L1*(IMO-1)))
      OCC=ONE
      DO 340 IPT=1,NPTMAX
      X(IPT+I10-1)=X(IPT+I10-1)+OCC*X(IPT+I12-1)**2
  340 X(IPT+I11-1)=X(IPT+I11-1)-OCC*X(IPT+I12-1)**2
  350 CONTINUE
  360 CONTINUE
C
C     ----- PRINT TOTAL ELECTRON AND SPIN DENSITY -----
C
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
       if (ga_nodeid().eq.0) then
         WRITE(IW,9992)
         WRITE(IW,9997)
C
*        WRITE(IP,9987)
*        WRITE(IP,9990) XRANGE(1),XRANGE(2),XRANGE(3),NXGRD
*        WRITE(IP,9990) YRANGE(1),YRANGE(2),YRANGE(3),NYGRD
*        WRITE(IP,9990) ZRANGE(1),ZRANGE(2),ZRANGE(3),NZGRD
       endif
C
      DO 370 IPT=1,NPTMAX
      XP=X(1+3*(IPT-1)+I00-1)
      YP=X(2+3*(IPT-1)+I00-1)
      ZP=X(3+3*(IPT-1)+I00-1)
      SPNDEN=X(IPT+I11-1)
      IF(ga_nodeid().eq.0.and.IPT.GT.NPTGRD.OR.OUT)
     &   WRITE(IW,9995) IPT,XP,YP,ZP,SPNDEN
  370 CONTINUE
c--   IF(ga_nodeid().eq.0.and.NPTGRD.GT.0) 
c--  &   WRITE(IP,9989) (X(IPT+I11-1),IPT=1,NPTGRD)
c--   if (ga_nodeid().eq.0) WRITE(IP,9988)
      ENDIF
C
      if (ga_nodeid().eq.0) then
         WRITE(IW,9993)
         WRITE(IW,9997)
C
*        WRITE(IP,9991)
*        WRITE(IP,9990) XRANGE(1),XRANGE(2),XRANGE(3),NXGRD
*        WRITE(IP,9990) YRANGE(1),YRANGE(2),YRANGE(3),NYGRD
*        WRITE(IP,9990) ZRANGE(1),ZRANGE(2),ZRANGE(3),NZGRD
      endif
C
      DO 380 IPT=1,NPTMAX
      XP=X(1+3*(IPT-1)+I00-1)
      YP=X(2+3*(IPT-1)+I00-1)
      ZP=X(3+3*(IPT-1)+I00-1)
      ELDENE=X(IPT+I10-1)
      IF(ga_nodeid().eq.0.and.IPT.GT.NPTGRD.OR.OUT)
     &   WRITE(IW,9995) IPT,XP,YP,ZP,ELDENE
  380 CONTINUE
c--   IF(ga_nodeid().eq.0.and.NPTGRD.GT.0) 
c--  &  WRITE(IP,9989) (X(IPT+I10-1),IPT=1,NPTGRD)
c--   if(ga_nodeid().eq.0) WRITE(IP,9988)
C
  500 CONTINUE

C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_DIPOLE, MA_POP_STACK OF INIT
     & FAILED',911)
C
      RETURN
 9999 FORMAT(/,10X,21(1H-),/,10X,'Electron/spin density',
     1       /,10X,21(1H-),/)
 9998 FORMAT(' Not enough core in -ELDMAP-')
 9997 FORMAT(3X,'Point',6X,'X',9X,'Y',9X,'Z',6X,'Density (a.u.)')
 9996 FORMAT(' 1 a.u. = 4.80286 10**(-10) esu ')
 9995 FORMAT(1X,I5,2X,3F10.5,F15.6)
 9994 FORMAT(/,' Density for (alpha) orbital = ',i5,/,
     1         ' ----------------------------------- ')
 9993 FORMAT(/,' Total electron density ',/,' ---------------------- ')
 9992 FORMAT(/,' Total spin density     ',/,' ------------------     ')
 9991 FORMAT(' $ELDMAP')
 9990 FORMAT(3F13.6,I5)
 9989 FORMAT(6F13.6,2X)
 9988 FORMAT(' $end   ')
 9987 FORMAT(' $SPNMAP')
 9986 FORMAT(/,' No grid defined. ')
 9985 FORMAT(/,1X,'( ',I4,' BY ',I4,' BY ',I4,' ) Grid defined.',
     1 /,' Iatom  = ',4I4,
     2 /,' Origin = ',3F10.5,
     3 /,' Xvec   = ',3F10.5,' XMIN=',F8.2,' XMAX=',F8.2,' XSTEP=',F8.2,
     4 /,' Yvec   = ',3F10.5,' YMIN=',F8.2,' YMAX=',F8.2,' YSTEP=',F8.2,
     5 /,' Zvec   = ',3F10.5,' ZMIN=',F8.2,' ZMAX=',F8.2,' ZSTEP=',F8.2,
     6 /)
 9984 FORMAT(/,' DENSITY FOR ( BETA) ORBITAL = ',I5,/,
     1         ' ----------------------------------- ')
 9983 FORMAT(/,' FOR -RUNTYP- .EQ. -PROPTY  - THERE IS NO TOTAL',
     1         ' DENSITY. STOP')
 9982 FORMAT(/,' WAVEFUNCTION FOR (ALPHA) ORBITAL = ',I5,/,
     1         ' ---------------------------------------- ')
 9981 FORMAT(/,' WAVEFUNCTION FOR ( BETA) ORBITAL = ',I5,/,
     1         ' ---------------------------------------- ')
 9980 FORMAT(3X,'POINT',6X,'X',9X,'Y',9X,'Z',6X,'  WFN   (A.U.)')
      END
      SUBROUTINE HND_ELWFNC(NPT,XYZPT,ELDEN,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
C     LOGICAL IANDJ
C     LOGICAL NORM,DOUBLE
      LOGICAL NORM
      LOGICAL OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZWFN/XWFN,YWFN,ZWFN,X0,Y0,Z0,XI,YI,ZI,NI
      DIMENSION ELDEN(1),XYZPT(3,1),V(1)
      DIMENSION DI(15)
      DIMENSION XG(5),YG(5),ZG(5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA TWO   /2.0D+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT=NPRINT.EQ.6
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      DO 10 IPT=1,NPT
   10 ELDEN(IPT)=ZERO
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- DENSITY FACTOR -----
C
      INUM=0
      DO 230 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI
      GO TO 220
  120 DUM1=CPI
      GO TO 220
  130 DUM1=CDI
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      INUM=INUM+1
  230 DI(INUM)=DUM1*V(LOCI+I)
C
C     ----- LOOP OVER POINTS TO BE EVALUATED -----
C
      DO 500 IPT=1,NPT
      X0 = XYZPT(1,IPT)
      Y0 = XYZPT(2,IPT)
      Z0 = XYZPT(3,IPT)
C
      DUM = AI*((X0-XI)**2+(Y0-YI)**2+(Z0-ZI)**2)
      IF(DUM.GT.TOL) GO TO 500
      FAC = EXP(-DUM)
C
C     ----- WAVEFUNCTION VALUES -----
C
      DO 370 I=1,LIT
      NI=I
      CALL HND_WFNXYZ
      XG(I)=XWFN
      YG(I)=YWFN
      ZG(I)=ZWFN
  370 CONTINUE
C
      INUM=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      INUM=INUM+1
      WFN=FAC*DI(INUM)*XG(IX)*YG(IY)*ZG(IZ)
      ELDEN(IPT)=ELDEN(IPT) + WFN
  390 CONTINUE
C
  500 CONTINUE
C
 7000 CONTINUE
C
 9000 CONTINUE
C
      RETURN
      END
      SUBROUTINE HND_WFNXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_XYZWFN/XWFN,YWFN,ZWFN,X0,Y0,Z0,XI,YI,ZI,NI
      DATA ONE /1.0D+00/
C
      XWFN=ONE
      YWFN=ONE
      ZWFN=ONE
      PTXI=X0-XI
      PTYI=Y0-YI
      PTZI=Z0-ZI
      GO TO (7,6,5,4,3,2,1),NI
    1 XWFN=XWFN*PTXI
      YWFN=YWFN*PTYI
      ZWFN=ZWFN*PTZI
    2 XWFN=XWFN*PTXI
      YWFN=YWFN*PTYI
      ZWFN=ZWFN*PTZI
    3 XWFN=XWFN*PTXI
      YWFN=YWFN*PTYI
      ZWFN=ZWFN*PTZI
    4 XWFN=XWFN*PTXI
      YWFN=YWFN*PTYI
      ZWFN=ZWFN*PTZI
    5 XWFN=XWFN*PTXI
      YWFN=YWFN*PTYI
      ZWFN=ZWFN*PTZI
    6 XWFN=XWFN*PTXI
      YWFN=YWFN*PTYI
      ZWFN=ZWFN*PTZI
    7 CONTINUE
      RETURN
      END
      SUBROUTINE HND_ELDENC(NPT,XYZPT,ELDEN,DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "hnd_tol.fh"
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      LOGICAL IANDJ
      LOGICAL NORM,DOUBLE
      LOGICAL OUT
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_XYZDEN/XINT,YINT,ZINT,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ
      DIMENSION ELDEN(1),XYZPT(3,NPT),DAB(1)
      DIMENSION DIJ(225)
      DIMENSION XIN(5,5),YIN(5,5),ZIN(5,5)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA TWO   /2.0D+00/
      DATA RLN10 /2.30258D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
C
      OUT=NPRINT.EQ.6
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      DO 10 IPT=1,NPT
   10 ELDEN(IPT)=ZERO
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
      IJ=0
      DO 100 I=MINI,MAXI
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 100 J=MINJ,JMAX
      IJ=IJ+1
  100 CONTINUE
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      IJ=IJ+1
  360 DIJ(IJ)=DUM2
C
C     ----- ELECTRON DENSITY -----
C
      DO 500 IPT=1,NPT
      X0 = XYZPT(1,IPT)
      Y0 = XYZPT(2,IPT)
      Z0 = XYZPT(3,IPT)
C
      DUM = AA*((X0-AX)**2+(Y0-AY)**2+(Z0-AZ)**2)
      IF(DUM.GT.TOL) GO TO 500
      FAC = EXP(-DUM)
C
C     ----- DENSITY VALUES -----
C
      DO 370 J=1,LJT
      NJ=J
      DO 370 I=1,LIT
      NI=I
      CALL HND_DENXYZ
      XIN(I,J)=XINT
      YIN(I,J)=YINT
      ZIN(I,J)=ZINT
  370 CONTINUE
C
      IJ=0
      DO 390 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 390 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      IJ=IJ+1
      ELD=FAC*DIJ(IJ)*XIN(IX,JX)*YIN(IY,JY)*ZIN(IZ,JZ)
      LOCII=LOCI+I
      LOCJJ=LOCJ+J
      NN=(LOCII*(LOCII-1))/2+LOCJJ
      DEN=DAB(NN)
      IF(.NOT.IANDJ.OR.J.NE.I) DEN=DEN+DEN
      ELDEN(IPT)=ELDEN(IPT)+ELD*DEN
  390 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
 8000 CONTINUE
 9000 CONTINUE
C
      RETURN
      END
      SUBROUTINE HND_DENXYZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HND_XYZDEN/XINT,YINT,ZINT,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ
      DATA ONE /1.0D+00/
C
      XINT=ONE
      YINT=ONE
      ZINT=ONE
      PTXI=X0-XI
      PTYI=Y0-YI
      PTZI=Z0-ZI
      PTXJ=X0-XJ
      PTYJ=Y0-YJ
      PTZJ=Z0-ZJ
      GO TO (7,6,5,4,3,2,1),NI
    1 XINT=XINT*PTXI
      YINT=YINT*PTYI
      ZINT=ZINT*PTZI
    2 XINT=XINT*PTXI
      YINT=YINT*PTYI
      ZINT=ZINT*PTZI
    3 XINT=XINT*PTXI
      YINT=YINT*PTYI
      ZINT=ZINT*PTZI
    4 XINT=XINT*PTXI
      YINT=YINT*PTYI
      ZINT=ZINT*PTZI
    5 XINT=XINT*PTXI
      YINT=YINT*PTYI
      ZINT=ZINT*PTZI
    6 XINT=XINT*PTXI
      YINT=YINT*PTYI
      ZINT=ZINT*PTZI
    7 GO TO (14,13,12,11,10,9,8),NJ
    8 XINT=XINT*PTXJ
      YINT=YINT*PTYJ
      ZINT=ZINT*PTZJ
    9 XINT=XINT*PTXJ
      YINT=YINT*PTYJ
      ZINT=ZINT*PTZJ
   10 XINT=XINT*PTXJ
      YINT=YINT*PTYJ
      ZINT=ZINT*PTZJ
   11 XINT=XINT*PTXJ
      YINT=YINT*PTYJ
      ZINT=ZINT*PTZJ
   12 XINT=XINT*PTXJ
      YINT=YINT*PTYJ
      ZINT=ZINT*PTZJ
   13 XINT=XINT*PTXJ
      YINT=YINT*PTYJ
      ZINT=ZINT*PTZJ
   14 CONTINUE
      RETURN
      END
      SUBROUTINE HND_SPIND(RTDB,BASIS,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      integer rtdb,basis,geom
      logical status
      PARAMETER (MXIODA=255)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
      CHARACTER*8 WFNTYP, SCFTYP
      CHARACTER*8 SCF, UHF, MP2
      LOGICAL SOME,OUT
      LOGICAL NOSPIN
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCM/X(1)
      DIMENSION A(6),VEC(3,3),EIG(3)
      DIMENSION ISOTOP(18),GN(18)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- SPIND-'/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA THREE  /3.0D+00/
      DATA FOUR   /4.0D+00/
      DATA UHF    /'UHF     '/
      DATA MP2    /'MP2     '/
      DATA SCF    /'SCF     '/
      DATA ANGLE  /180.D+00/
      DATA GE     /2.002319304386D+00/
      DATA EMF    /1836.152701D+00/
      DATA VL     /137.0359895D+00/
      DATA AUEV   /27.2113961D+00/
      DATA EVMHZ  /2.41798836D+08/
      DATA GMHZ   /2.8025D+00/
      DATA GN /5.58556D+00  , -4.255D+00 ,2.1709333D+00,-0.7800666D+00,
     1         1.7923333D+00,  1.4048D+00,0.4036D+00   ,-0.75748D+00  ,
     2         5.2576D+00   , -0.4412D+00,1.4783333D+00,-0.34212D+00  ,
     3         1.45656D+00  ,-1.11106D+00,2.2634D+00   , 0.4288666D+00,
     4         0.5478866D+00,-0.3714285D+00/
      DATA ISOTOP / 1, 3, 7, 9,11,13,14,17,19,
     1             21,23,25,27,29,31,33,35,39/
C
      NOSPIN=.FALSE.
      IF(NOSPIN) RETURN
      NOSPIN=(WFNTYP.NE.SCF.OR.SCFTYP.NE.UHF).AND.
     1       (WFNTYP.NE.MP2.OR.SCFTYP.NE.UHF)
      IF(NOSPIN) WRITE(IW,9993)
      IF(NOSPIN) RETURN
C
      SOME=.TRUE.
      OUT =NPRINT.EQ.6
      WRITE(IW,9999)
      WRITE(IW,9985)
C
C     --- CALCULATE CONSTANTS AND CONVERSION TERMS ---
C
      PI   = ACOS(-ONE)
      FAC  =(FOUR*PI/THREE)
      BETAE=ONE/(TWO*VL)
      BETAN=BETAE/EMF
      CONVF=AUEV*EVMHZ
      CON  =GE*BETAE*BETAN*CONVF
      DEG  =ANGLE/PI
C
C     ----- DEFINE POINTS FOR CALCULATION -----
C
      DO 10 IAT=1,NAT
      X(1+3*(IAT-1))=C(1,IAT)
      X(2+3*(IAT-1))=C(2,IAT)
      X(3+3*(IAT-1))=C(3,IAT)
   10 CONTINUE
      NPTMAX=NAT
C
C     ----- CALCULATE SPIN DENSITY AT THE NUCLEI    -----
C           THIS IS FOR THE ISOTROPIC INTERACTION
C           ( FERMI CONTACT TERM ) OF THE HYPERFINE
C           COUPLING TENSOR.
C
C     ----- SPECIAL CASE = ONE PARTICLE REDUCED DENSITY MATRIX -----
C           IS DIAGONAL.
C
      L1= NUM
      L2=(NUM*(NUM+1))/2
      L3= NUM* NUM
      I00 =    1
      I10 =I00+NPTMAX*3
      I11 =I10+NPTMAX
      I20 =I11+NPTMAX
      I30 =I20+NUM*NUM
      LAST=I30+NUM*NUM
      NEED=LAST-1
C
C     ----- GET CORE MEMORY -----
C
c%%
C     CALL HND_CMEM(LOADCM)
C     IF(NEED.LE.MAXCOR) GO TO 20
C     WRITE(IW,9998)
C     CALL HND_NEEDCM(LAST,NEED)
C     CALL HND_HNDERR(3,ERRMSG)
C  20 CONTINUE
c%%
C
C     ----- READ MOLECULAR ORBITALS -----
C
      CALL HND_DAREAD(IDAF,IODA,X(I20),L3,15)
      CALL HND_DAREAD(IDAF,IODA,X(I30),L3,19)
      DO 310 IPT=1,NPTMAX
  310 X(IPT+I10-1)=ZERO
C
c--
      NA=0
c--
      DO 330 IMO=1,NA
      CALL HND_ELWFNC(NPTMAX,X(I00),X(I11),X(I20+L1*(IMO-1)))
      DO 320 IPT=1,NPTMAX
  320 X(IPT+I10-1)=X(IPT+I10-1)+X(IPT+I11-1)**2
  330 CONTINUE
C
c--
      NB=0
c--
      DO 350 IMO=1,NB
      CALL HND_ELWFNC(NPTMAX,X(I00),X(I11),X(I30+L1*(IMO-1)))
      DO 340 IPT=1,NPTMAX
  340 X(IPT+I10-1)=X(IPT+I10-1)-X(IPT+I11-1)**2
  350 CONTINUE
C
C     ----- PRINT SPIN DENSITY AND -AISO- AT THE NUCLEI -----
C
      WRITE(IW,9994)
      WRITE(IW,9997)
C
      DO 370 IPT=1,NPTMAX
      IAZ=ZAN(IPT)
      IF(IAZ.LE.18) THEN
      GNU=GN(IAZ)
      ELSE
      GNU=ONE
      ENDIF
C
      XP=X(1+3*(IPT-1)+I00-1)
      YP=X(2+3*(IPT-1)+I00-1)
      ZP=X(3+3*(IPT-1)+I00-1)
      SPNDEN=X(IPT+I10-1)
      AISO  =TWO*FAC*CON*GNU*SPNDEN
      WRITE(IW,9995) IPT,IAZ,XP,YP,ZP,SPNDEN,AISO,AISO/GMHZ
  370 CONTINUE
C
C     ----- CALCULATE SPIN-DIPOLAR TERM AT THE NUCLEI -----
C           THIS IS FOR THE ANISOTROPIC INTERACTION
C           OF THE HYPERFINE COUPLING TENSOR.
C
      I00 =    1
      I10 =I00+NPTMAX*3
      I11 =I10+NPTMAX
      I20 =I11+NPTMAX*6
      I30 =I20+(NUM*(NUM+1))/2
      LAST=I30+(NUM*(NUM+1))/2
      NEED=LAST-1
C
C     ----- GET CORE MEMORY -----
C
c%%
C     CALL HND_CMEM(LOADCM)
C     IF(NEED.LE.MAXCOR) GO TO 410
C     WRITE(IW,9998)
C     CALL HND_NEEDCM(LAST,NEED)
C     CALL HND_HNDERR(3,ERRMSG)
C 410 CONTINUE
c%%
C
C     ----- GET SPIN DENSITY MATRIX -----
C
      CALL HND_DAREAD(IDAF,IODA,X(I20),L2,16)
      IF((WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF).OR.
     1   (WFNTYP.EQ.MP2.AND.SCFTYP.EQ.UHF)    ) THEN
         CALL HND_DAREAD(IDAF,IODA,X(I30),L2,20)
         DO 420 I=1,L2
            X(I+I20-1)=X(I+I20-1)-X(I+I30-1)
  420    CONTINUE
      ENDIF
C
C     ----- CALCULATE SPIN-DIPOLAR TERM AT ALL POINTS -----
C
      CALL HND_ELFGRD(NPTMAX,X(I00),X(I11),X(I20))
C
      WRITE(IW,9996)
      WRITE(IW,9992)
C
      DO 430 IPT=1,NPTMAX
      IAZ=ZAN(IPT)
      XP=X(1+3*(IPT-1)+I00-1)
      YP=X(2+3*(IPT-1)+I00-1)
      ZP=X(3+3*(IPT-1)+I00-1)
      SPNDXX  = X(1+6*(IPT-1)+I11-1) + FAC*X(IPT+I10-1)
      SPNDYY  = X(2+6*(IPT-1)+I11-1) + FAC*X(IPT+I10-1)
      SPNDZZ  = X(3+6*(IPT-1)+I11-1) + FAC*X(IPT+I10-1)
      SPNDXY  = X(4+6*(IPT-1)+I11-1)
      SPNDXZ  = X(5+6*(IPT-1)+I11-1)
      SPNDYZ  = X(6+6*(IPT-1)+I11-1)
      WRITE(IW,9991) IPT,IAZ,XP,YP,ZP,SPNDXX,SPNDYY,SPNDZZ,
     1                                SPNDXY,SPNDXZ,SPNDYZ
      X(1+6*(IPT-1)+I11-1) = SPNDXX
      X(2+6*(IPT-1)+I11-1) = SPNDYY
      X(3+6*(IPT-1)+I11-1) = SPNDZZ
      X(4+6*(IPT-1)+I11-1) = SPNDXY
      X(5+6*(IPT-1)+I11-1) = SPNDXZ
      X(6+6*(IPT-1)+I11-1) = SPNDYZ
  430 CONTINUE
C
      WRITE(IW,9990)
C
      DO 440 IPT=1,NPTMAX
      IAZ=ZAN(IPT)
      IF(IAZ.LE.18) THEN
      GNU=GN(IAZ)
      ELSE
      GNU=ONE
      ENDIF
      A(1) = X(1+6*(IPT-1)+I11-1)
      A(2) = X(4+6*(IPT-1)+I11-1)
      A(3) = X(2+6*(IPT-1)+I11-1)
      A(4) = X(5+6*(IPT-1)+I11-1)
      A(5) = X(6+6*(IPT-1)+I11-1)
      A(6) = X(3+6*(IPT-1)+I11-1)
      CALL HND_DIASPN(A,VEC,EIG,3)
      FAC1 = CON*GNU
      FAC2 = CON*GNU/GMHZ
      WRITE(IW,9989) IPT,IAZ,EIG(1),EIG(2),EIG(3),
     1                       EIG(1)*FAC1,EIG(2)*FAC1,EIG(3)*FAC1,
     2                       EIG(1)*FAC2,EIG(2)*FAC2,EIG(3)*FAC2
  440 CONTINUE
C
      WRITE(IW,9988)
C
      DO 460 IPT=1,NPTMAX
      A(1) = X(1+6*(IPT-1)+I11-1)
      A(2) = X(4+6*(IPT-1)+I11-1)
      A(3) = X(2+6*(IPT-1)+I11-1)
      A(4) = X(5+6*(IPT-1)+I11-1)
      A(5) = X(6+6*(IPT-1)+I11-1)
      A(6) = X(3+6*(IPT-1)+I11-1)
      CALL HND_DIASPN(A,VEC,EIG,3)
      DO 450 J=1,3
      DO 450 I=1,3
  450 VEC(I,J)= ACOS(VEC(I,J))*DEG
      WRITE(IW,9987) IPT
      WRITE(IW,9986) ((VEC(I,J),J=1,3),I=1,3)
  460 CONTINUE
C
C     -----  RESET CORE MEMORY -----
C
c%%
C     CALL HND_CMEM(NGOTCM)
C     IF(NGOTCM.NE.LOADCM) CALL HND_SETC(LOADCM)
c%%
      RETURN
 9999 FORMAT(/,10X,22(1H-),/,10X,'HYPERFINE INTERACTIONS',
     1       /,10X,22(1H-),/)
 9998 FORMAT(' NOT ENOUGH CORE IN - SPIND-')
 9997 FORMAT(3X,'ATOM ','ZNUC ',4X,'X',9X,'Y',9X,'Z',9X,
     1          'DENSITY (A.U.)',4X,'AISO(MHZ)',5X,'AISO(GAUSS)')
 9996 FORMAT(/,' ANISOTROPIC INTERACTION (SPIN-DIPOLAR TERM)',/,
     1         ' ----------------------- -------------------')
 9995 FORMAT(1X,I5,I4,2X,3F10.5,2X,F13.6,5X,F13.6,2X,F13.6)
 9994 FORMAT(/,' TOTAL SPIN DENSITY (FERMI CONTACT TERM)',/,
     1         ' ------------------ --------------------')
 9993 FORMAT(/,' HYPERFINE INTERACTIONS AVAILABLE FOR -UHF- ONLY . ')
 9992 FORMAT(3X,'ATOM ','ZNUC ',4X,'X',9X,'Y',9X,'Z',8X,1X,
     1 32X,'SPIN-DIPOLAR TERM (A.U.)',/,
     2 43X,6X,'XX',13X,'YY',13X,'ZZ',13X,'XY',13X,'XZ',13X,'YZ')
 9991 FORMAT(1X,I5,I4,2X,3F10.5,6F15.6)
 9990 FORMAT(/,12X,'PRINCIPAL COMP. OF -HF- TENSOR (AU) ',3X,2X,
     1          8X,'ANISOTROPIC -HFCC- (MHZ)',3X,2X,
     2          8X,'ANISOTROPIC -HFCC- (GAUSS)')
 9989 FORMAT(1X,I5,I4,3F13.6,2X,3F13.6,2X,3F13.6)
 9988 FORMAT(/,' ORIENTATION ( IN DEGREES ) OF THE PRINCIPAL AXIS',
     1         ' OF HYPERFINE TENSOR W.R.T ABSOLUTE (MOLECULAR) FRAME')
 9987 FORMAT(/,' ATOM = ',I5,/,' ------ ')
 9986 FORMAT(1X,3F15.6)
 9985 FORMAT(' ..... FOR ATOMS WITH Z.GT.18 , NEITHER THE -AISO- ',
     1       ' NOR THE ANISOTROPIC -HFCC- HAVE BEEN MULTIPLIED',
     2       ' BY -G(N)- ')
      END
      SUBROUTINE HND_DIASPN(A,VEC,EIG,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(6),H(3,3),VEC(3,3),EIG(3),BIG(6)
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA CONV   /1.0D-05/
      DATA MAXIT  /50/
C
      DO 20 I=1,3
         DO 10 J=1,3
   10       VEC(J,I)=ZERO
   20    VEC(I,I)=ONE
      IJ=0
      DO 25 I=1,3
      DO 25 J=1,I
         IJ=IJ+1
         H(I,J)=A(IJ)
   25    H(J,I)=A(IJ)
      CALL HND_DIAJAC(H,VEC,EIG,NDIM,NDIM,BIG)
      RETURN
      END
C
      SUBROUTINE HND_GRDMAP(XYZ,NXYZG,NXG,NYG,NZG,
     1                  XRANGE,YRANGE,ZRANGE,XVEC,YVEC,ZVEC,
     2                  ORIGIN,IATOM,IUNIT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 ERRMSG
      character*8 anam,bflab
      character*2 bnam
      LOGICAL PUNCH
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      DIMENSION XYZ(3,1)
      DIMENSION UX(3),UY(3),UZ(3),XYZG0(3)
      DIMENSION XRANGE(3),YRANGE(3),ZRANGE(3)
      DIMENSION ORIGIN(3),XVEC(3),YVEC(3),ZVEC(3)
      DIMENSION IATOM(4)
      DIMENSION COOR(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-GRDMAP-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA TM8    /1.0D-08/
      DATA TM2    /1.0D-02/
      DATA BOHR   /0.52917715D+00/
C
C       DETERMINE PARAMETERS TO GENERATE A 3-DIMENSIONAL RECTANGULAR
C       GRID OF POINTS.
C       PARAMETERS INCLUDE -
C
C         UX(I),UY(I),UZ(I) : UNIT VECTORS DEFINING GRID ORIENTATION
C         XYZG0             : COORDINATES OF POINT (XMIN,YMIN,ZMIN)
C         NXG               : NUMBER OF INCREMENTS ALONG X-AXIS
C         NXYG              : NXG  * NUMBER OF INCREMENTS ALONG Y-AXIS
C         NXYZG             : NXYG * NUMBER OF INCREMENTS ALONG Z-AXIS
C
C       INPUT PARAMETERS INCLUDE -
C
C         XRANGE(I) : XMIN, XMAX, XSTEP: POINTS ALONG X-AXIS OF THE GRID
C         YRANGE(I) : YMIN, YMAX, YSTEP: POINTS ALONG Y-AXIS OF THE GRID
C         ZRANGE(I) : ZMIN, ZMAX, ZSTEP: POINTS ALONG Z-AXIS OF THE GRID
C         IUNIT     : UNITS (BOHR OR ANGS)
C         IATOM(I)  : SERIAL NUMBER OF ATOMS DEFINING GRID PLANES
C         ORIGIN(I) : COORDINATES OF POINT DEFINING CUBE ORIGIN
C         XVEC(I)   : COORDINATES OF POINT DEFINING CUBE X-AXIS
C         YVEC(I)   : COORDINATES OF POINT DEFINING CUBE +Y DIRECTION
C
C       THE CUBE Z-DIRECTION IS DEFINED AS THE NORMAL ON THE PLANE
C       DEFINED BY X- AND Y-AXIS.
C
C       TO DEFINE THE CUBE OF THE GRID, EITHER IATOM(I=1,4) OR ORIGIN,
C       XVEC, YVEC AND ZVEC, OR A COMBINATION OF THE TWO ARE INPUT.
C       FOR EVERY NON-ZERO IATOM(I), THE FOLLOWING SUBSTITUTION IS MADE
C
C         IATOM(1) -> ORIGIN    IATOM(2) -> XVEC
C                               IATOM(3) -> YVEC
C                               IATOM(4) -> ZVEC
C
C       ----- FOR EXAMPLE, FOR THE INPUT
C
C         $MAP
C          IATOM=0,3,0,0, ORIGIN=0.,0.,1.,
C          YVEC=0.,2.,0., ZVEC=0.,0.,5., IUNIT=-1,
C         $END
C
C       XVEC(I=1,3) ARE SET EQUAL TO THE COORDINATES OF ATOM 3.
C
C       THE GRID CENTER IS AT -ORIGIN- . THE UNITS ARE A.U.
C
      PUNCH=.FALSE.
C
      IF(IUNIT.GT.0) THEN
         DO 10 I=1,3
            XRANGE(I)=XRANGE(I)/BOHR
            YRANGE(I)=YRANGE(I)/BOHR
            ZRANGE(I)=ZRANGE(I)/BOHR
            XVEC(I)=XVEC(I)/BOHR
            YVEC(I)=YVEC(I)/BOHR
            ZVEC(I)=ZVEC(I)/BOHR
            ORIGIN(I)=ORIGIN(I)/BOHR
   10    CONTINUE
      ENDIF
C
      XMIN =XRANGE(1)
      XMAX =XRANGE(2)
      XSTEP=XRANGE(3)
      YMIN =YRANGE(1)
      YMAX =YRANGE(2)
      YSTEP=YRANGE(3)
      ZMIN =ZRANGE(1)
      ZMAX =ZRANGE(2)
      ZSTEP=ZRANGE(3)
C
      XSTEP= ABS(XSTEP)
      YSTEP= ABS(YSTEP)
      ZSTEP= ABS(ZSTEP)
      IF (XSTEP.LT.TM2) XSTEP = TM2
      IF (YSTEP.LT.TM2) YSTEP = TM2
      IF (ZSTEP.LT.TM2) ZSTEP = TM2
C
C     --- IF PLANE DEFINED BY ATOMS, GET ATOMIC COORDINATES ---
C
      IATOM1=IATOM(1)
      IATOM2=IATOM(2)
      IATOM3=IATOM(3)
      DO 20 I=1,3
         IF(IATOM1 .NE. 0) ORIGIN(I) = C(I,IATOM1)
         IF(IATOM2 .NE. 0) XVEC(I)   = C(I,IATOM2)
         IF(IATOM3 .NE. 0) YVEC(I)   = C(I,IATOM3)
   20 CONTINUE
C
C     --- DETERMINE VECTORS DEFINING GRID X-AXIS AND POSITIVE-Y DIRECTIO
C
      DO 30 I=1,3
         UX(I) = XVEC(I) - ORIGIN(I)
         UY(I) = YVEC(I) - ORIGIN(I)
   30 CONTINUE
C
C     --- DETERMINE VECTOR DEFINING GRID Z-AXIS
C
      UZ(1) = UX(2)*UY(3) - UX(3)*UY(2)
      UZ(2) = UX(3)*UY(1) - UX(1)*UY(3)
      UZ(3) = UX(1)*UY(2) - UX(2)*UY(1)
C
C     --- CHECK FOR COLINEARITY OF DEFINING VECTORS
C
      UMOD  = SQRT(UZ(1)*UZ(1) + UZ(2)*UZ(2) + UZ(3)*UZ(3))
      IF (UMOD .LE. TM8) THEN
         WRITE(IW,9999)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     --- NORMALIZE VECTOR DEFINING GRID Z-AXIS
C
      UZ(1) = UZ(1)/UMOD
      UZ(2) = UZ(2)/UMOD
      UZ(3) = UZ(3)/UMOD
C
C     --- DETERMINE VECTOR DEFINING GRID Y-AXIS
C
      UY(1) = UZ(2)*UX(3) - UZ(3)*UX(2)
      UY(2) = UZ(3)*UX(1) - UZ(1)*UX(3)
      UY(3) = UZ(1)*UX(2) - UZ(2)*UX(1)
C
C     --- NORMALIZE VECTOR DEFINING GRID Y-AXIS
C
      UMOD  = SQRT(UY(1)*UY(1) + UY(2)*UY(2) + UY(3)*UY(3))
      UY(1) = UY(1)/UMOD
      UY(2) = UY(2)/UMOD
      UY(3) = UY(3)/UMOD
C
C     --- NORMALIZE VECTOR DEFINING GRID X-AXIS
C
      UMOD  = SQRT(UX(1)*UX(1) + UX(2)*UX(2) + UX(3)*UX(3))
      UX(1) = UX(1)/UMOD
      UX(2) = UX(2)/UMOD
      UX(3) = UX(3)/UMOD
C
C     --- DETERMINE NUMBER OF INCREMENTS IN X,Y,Z DIRECTIONS
C
      NXG  = INT((XMAX-XMIN)/XSTEP) + 1
      NYG  = INT((YMAX-YMIN)/YSTEP) + 1
      NZG  = INT((ZMAX-ZMIN)/ZSTEP) + 1
      XMAXP = XMIN + XSTEP*DBLE(NXG-1)
      YMAXP = YMIN + YSTEP*DBLE(NYG-1)
      ZMAXP = ZMIN + ZSTEP*DBLE(NZG-1)
      IF ((XMAX-XMAXP).GT.(XSTEP/TWO)) NXG=NXG+1
      IF ((YMAX-YMAXP).GT.(YSTEP/TWO)) NYG=NYG+1
      IF ((ZMAX-ZMAXP).GT.(ZSTEP/TWO)) NZG=NZG+1
      IF (NXG.LT.0 .OR. NYG.LT.0 .OR. NZG.LT.0) THEN
         IF (NXG.LT.0) WRITE(IW,9998)
         IF (NYG.LT.0) WRITE(IW,9997)
         IF (NZG.LT.0) WRITE(IW,9996)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      XMAX = XMIN + XSTEP*DBLE(NXG-1)
      YMAX = YMIN + YSTEP*DBLE(NYG-1)
      ZMAX = ZMIN + ZSTEP*DBLE(NZG-1)
      XRANGE(2) = XMAX
      YRANGE(2) = YMAX
      ZRANGE(2) = ZMAX
C
C     --- COMPUTE COORDINATES OF STARTING POINT (XMIN,YMIN,ZMIN)
C
      DO 40 I=1,3
         XYZG0(I) = ORIGIN(I) + XMIN*UX(I) + YMIN*UY(I) + ZMIN*UZ(I)
   40 CONTINUE
C
C     --- CALCULATE COORDINATES OF ATOMS IN NEW COORDINATES SYSTEM ---
C
      DO 50 IAT=1,NAT
         COOR(1) = (C(1,IAT)-ORIGIN(1))*UX(1) +
     1             (C(2,IAT)-ORIGIN(2))*UX(2) +
     2             (C(3,IAT)-ORIGIN(3))*UX(3)
         COOR(2) = (C(1,IAT)-ORIGIN(1))*UY(1) +
     1             (C(2,IAT)-ORIGIN(2))*UY(2) +
     2             (C(3,IAT)-ORIGIN(3))*UY(3)
         COOR(3) = (C(1,IAT)-ORIGIN(1))*UZ(1) +
     1             (C(2,IAT)-ORIGIN(2))*UZ(2) +
     2             (C(3,IAT)-ORIGIN(3))*UZ(3)
         NUCZ = ZAN(IAT)
   50 CONTINUE
C
C     --- RENORMALIZE UNIT VECTORS TO STEP SIZE ---
C
      DO 60 I=1,3
         UX(I) = UX(I)*XSTEP
         UY(I) = UY(I)*YSTEP
         UZ(I) = UZ(I)*ZSTEP
   60 CONTINUE
C
      NXYZG=NXG*NYG*NZG
      IF(NXYZG.EQ.1) NXYZG=0
      IF(NXYZG.EQ.0) GO TO 100
C
      IF(PUNCH) WRITE(IP,9995)
      IF(PUNCH) WRITE(IP,9994)
C
      IXYZG=0
      IZG=0
   70 IZG=IZG+1
      IYG=0
   80 IYG=IYG+1
      IXG=0
   90 IXG=IXG+1
C
      IXYZG=IXYZG+1
      XG = XYZG0(1) + (IXG-1)*UX(1) + (IYG-1)*UY(1) + (IZG-1)*UZ(1)
      YG = XYZG0(2) + (IXG-1)*UX(2) + (IYG-1)*UY(2) + (IZG-1)*UZ(2)
      ZG = XYZG0(3) + (IXG-1)*UX(3) + (IYG-1)*UY(3) + (IZG-1)*UZ(3)
      XYZ(1,IXYZG)=XG
      XYZ(2,IXYZG)=YG
      XYZ(3,IXYZG)=ZG
C
      IF(PUNCH) WRITE(IP,9993) XG,YG,ZG
C
      IF(IXG.LT.NXG) GO TO 90
      IF(IYG.LT.NYG) GO TO 80
      IF(IZG.LT.NZG) GO TO 70
C
      IF(PUNCH) WRITE(IP,9992)
C
  100 CONTINUE
C
C     ----- PUNCH ATOMIC COORDINATES -----
C
C     WRITE(IP,9991)
C     DO 110 IAT=1,NAT
C     ZNUC=NUC(IAT)
C 110 WRITE(IP,9990) ANAM(IAT),BNAM(IAT),ZNUC,(C(I,IAT),I=1,3)
C     WRITE(IP,9989)
C
      RETURN
 9999 FORMAT(1H ,'--- $GRID - DEFINING VECTORS ARE COLINEAR, STOP')
 9998 FORMAT(' --- IN $MAP , INCONSISTENT -XRANGE- DEFINITION, STOP')
 9997 FORMAT(' --- IN $MAP , INCONSISTENT -YRANGE- DEFINITION, STOP')
 9996 FORMAT(' --- IN $MAP , INCONSISTENT -ZRANGE- DEFINITION, STOP')
 9995 FORMAT(' $MAP   ')
 9994 FORMAT(' XYZMAP=')
 9993 FORMAT(8X,F10.5,1H,,F10.5,1H,,F10.5,1H,)
 9992 FORMAT(' $END   ')
 9991 FORMAT(' $XYZ   ')
 9990 FORMAT(A8,A2,F5.0,3F20.10)
 9989 FORMAT(' $END   ')
      END
      SUBROUTINE HND_STONEX
      implicit DOUBLE PRECISION (a-h,o-z)
#include "nwc_const.fh"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     data $stone
c     -----------
c
c     function
c     ........ the data group defines the expansion points for the
c              distributed multipole decomposition of the electrostatic
c              potential.  points may be defined at atom centers,
c              bond midpoints, bond overlap centers of charge, molecular
c              center of mass, or any other user defined point.
c              the several types of input lines (format-free) are
c              described below. enter as many lines as you wish,
c              terminated by a line with $end.
c
c              the distributed multipole method takes the multipolar
c              expansion of each overlap charge density defined by two
c              gaussian primitives, and translates it from the
c              center of charge of the overlap density to the nearest
c              expansion point.
c
c     atom  i  <name>
c
c              'atom'    a keyword indicating that an atom is selected
c                        as an expansion center.
c
c               i        the number of the atom, corresponding to an
c                        atom in the $basis input.
c
c              <name>    an optional name for the atom.  if not entered
c                        the name will be set to the name used in the
c                        $basis input.
c
c     bond  i  j  <name>
c
c              'bond'    keyword indicating that a bond midpoint is
c                        selected as an expansion center.
c
c               i,j      the indices of the atoms defining the bond,
c                        corresponding to two atoms in the $basis input.
c
c              <name>    an optional name for the bond midpoint.
c                        if not entered the name will be set to 'bond'.
c
c     olap  i  j  <name>
c
c              'olap'    keyword indicating that the center of charge
c                        of a bond overlap distibution is selected as
c                        an expansion center.
c
c               i,j      the indices of the atoms defining the bond,
c                        corresponding to two atom in the $basis input.
c
c              <name>    an optional name for the overlap expansion
c                        point.  if omitted, it will be set to 'olap'.
c
c     point  x  y  z  <name>
c
c              'point'   keyword indicating that a user defined point
c                        is selected as an expansion point.
c
c               x,y,z    the coordinates of the point (in a.u.)
c
c              <name>    an optional name for the point.  if omitted,
c                        it will be set to 'point'
c
c     cmass              a line with 'cmass' as the keyword will
c                        select the center of mass of the molecule
c                        as an expansion point.  no other input on
c                        the line is necessary.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter (mxioda=255)
      parameter (nummx=3072)
      parameter (MXATOM=NW_MAX_ATOM)
      character*8   errmsg
      character*255 input
      character*8   word1,word2,word3,word4,word5
      character*8   atom,bond,point,olap,cmass
      character*8   wrdkey,wrdstn,wrdend,endwrd
      character*8   blank,debug
      character*8   anam,bflab
      character*2   bnam
      character*8   pnttyp
      character*8   wrd1,wrd2,wrd3,wrd4,wrd5
c%%
      logical dbug
      logical some
c%%
      COMMON/HND_infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
c%%
      COMMON/HND_mollab/anam(MXATOM),bnam(MXATOM),bflab(nummx)
      COMMON/HND_runlab/title(10)
      COMMON/HND_iofile/ir,iw,ip
      COMMON/HND_dafile/idaf,nav,ioda(mxioda)
      COMMON/HND_stnpnt/pnt(4,256),pnttyp(256),npts
      COMMON/HND_output/nprint
      COMMON/HND_scm/x(1)
      dimension cm(3),amass(MXATOM)
      dimension errmsg(3)
      equivalence (wrd1,word1),(wrd2,word2),(wrd3,word3),
     1            (wrd4,word4),(wrd5,word5)
c
      data errmsg          /'program ','stop in ','-stonex-'/
      data wrdstn          /' $STONE '/
      data wrdend          /' $END   '/
      data cmass           /'CMASS'/
      data atom,bond,point /'ATOM','BOND','POINT'/
      data olap            /'OLAP'/
      data debug           /'debug'/
      data blank           /'        '/
      data zero            /0.0d+00/
      data pt5             /0.5d+00/
      data maxpnt          /256/
c
      dbug=.false.
      write(iw,9999)
c
      l1=num
      l2=(num*(num+1))/2
c
c     ----- set pointers for partitioning of core -----
c
      i10=1
      i20=i10+l2
      i30=i20+l2
      i40=i30+l2
      i50=i40+l2
      i60=i50+l2
c
c     ----- set up defaults : use all the atoms -----
c           as expansion centers.
c
      word1=atom
      ipnt =nat
      do 1 i=1,ipnt
c%%
c        pnt(1,i)=anam(i)
c%%
         pnt(2,i)=c(1,i)
         pnt(3,i)=c(2,i)
         pnt(4,i)=c(3,i)
    1    pnttyp(i)=wrd1
c
c     ----- get density and overlap matrices -----
c           -d- at x(i10)
c           -s- at x(i20)
c           -x- at x(i30)
c           -y- at x(i40)
c           -z- at x(i50)
c
      CALL HND_daread(idaf,ioda,x(i10),l2,16)
      CALL HND_daread(idaf,ioda,x(i20),l2,12)
c
c     ----- compute dipole matrices -----
c           -x- at x(i30)
c           -y- at x(i40)
c           -z- at x(i50)
c
      do 2 i=1,3
    2    cm(i)=zero
      CALL HND_dipint(x(i30),x(i40),x(i50),cm,l1)
c
c     ----- compute coordinates of center of mass -----
c
      ncall=1
      CALL HND_dipams(amass,ncall,ncode,some)
      sum=zero
      do 3 i=1,3
    3    cm(i)=zero
      do 4 iat=1,nat
         sum=sum+amass(iat)
         do 4 i=1,3
    4       cm(i)=cm(i)+amass(iat)*c(i,iat)
      do 5 i=1,3
    5    cm(i)=cm(i)/sum
c
c     ----- read input data for distributed multipole analysis -----
c
      rewind ir
   10 read(ir,9994,end=20,err=20) wrdkey
      go to 30
   20 write(iw,9997)
      go to 200
   30 if(wrdkey.ne.wrdstn) go to 10
c
      ipnt=0
   40 continue
      CALL HND_rdfree(ir,input,iend)
      if(iend.ne.0) then
         CALL HND_hnderr(3,errmsg)
      endif
      word1=blank
      word2=blank
      word3=blank
      word4=blank
      word5=blank
      read(input,*,end=50,err=50) word1,word2,word3,word4,word5
   50 continue
c
      endwrd=word1
c%%% need to be ported
c     CALL HND_wrdnam(endwrd)
c%%%
      if(endwrd.eq.wrdend) go to 200
c
      if(word1.eq.debug) then
         dbug=.true.
c
c     ----- -atom- -----
c
      elseif(word1.eq.atom ) then
         ipnt=ipnt+1
         if(ipnt.gt.maxpnt) then
            write(iw,9991) ipnt,maxpnt
            CALL HND_hnderr(3,errmsg)
         endif
         word1=blank
         iatm =0
         word3=blank
         read(input,*,end=110,err=110) word1,iatm,word3
  110    if(iatm.le.0) then
            write(iw,9993) word1,iatm
            write(iw,9992)
            return
         endif

         if(word3.eq.blank) wrd3=anam(iatm)
c%%
c        pnt(1,ipnt)=wrd3
c%%
         do 120 i=1,3
  120       pnt(i+1,ipnt)=c(i,iatm)
         pnttyp(ipnt)=wrd1
         go to 40
c
c     ----- -bond- -----
c
      elseif(word1.eq.bond ) then
         ipnt=ipnt+1
         if(ipnt.gt.maxpnt) then
            write(iw,9991) ipnt,maxpnt
            CALL HND_hnderr(3,errmsg)
         endif
         word1=blank
         iatm =0
         jatm =0
         word4=blank
         read(input,*,end=130,err=130) word1,iatm,jatm,word4
  130    if(iatm.le.0.or.jatm.le.0) then
            write(iw,9993) word1,iatm,jatm
            write(iw,9992)
            return
         endif
         if(word4.eq.blank) word4=bond
c%%
c        pnt(1,ipnt)=wrd4
c%%
         do 140 i=1,3
  140       pnt(i+1,ipnt)=pt5*(c(i,iatm)+c(i,jatm))
         pnttyp(ipnt)=wrd1
         go to 40
c
c     ----- -olap- -----
c
      elseif(word1.eq.olap ) then
         ipnt=ipnt+1
         if(ipnt.gt.maxpnt) then
            write(iw,9991) ipnt,maxpnt
            CALL HND_hnderr(3,errmsg)
         endif
         word1=blank
         iatm =0
         jatm =0
         word4=blank
         read(input,*,end=150,err=150) word1,iatm,jatm,word4
  150    if(iatm.le.0) then
            write(iw,9993) word1,iatm,jatm
            write(iw,9992)
            return
         endif
         if(word4.eq.blank) word4=olap
c%%
c        pnt(1,ipnt)=wrd4
c%%
         CALL HND_stnolp(x(i10),x(i20),x(i30),x(i40),x(i50),
     1               pnt(2,ipnt),iatm,jatm)
         pnttyp(ipnt)=wrd1
         go to 40
c
c     ----- -point- -----
c
      elseif(word1.eq.point) then
         ipnt=ipnt+1
         if(ipnt.gt.maxpnt) then
            write(iw,9991) ipnt,maxpnt
            CALL HND_hnderr(3,errmsg)
         endif
         xpnt=zero
         ypnt=zero
         zpnt=zero
         read(input,*,end=160,err=160) word1,xpnt,ypnt,zpnt,word5
  160    pnt(1+1,ipnt)=xpnt
         pnt(2+1,ipnt)=ypnt
         pnt(3+1,ipnt)=zpnt
         if(word5.eq.blank) word5=point
c%%
c        pnt(  1,ipnt)=wrd5
c%%
         pnttyp(ipnt)=wrd1
         go to 40
c
c     ----- -cmass- -----
c
      elseif(word1.eq.cmass) then
         ipnt=ipnt+1
         if(ipnt.gt.maxpnt) then
            write(iw,9991) ipnt,maxpnt
            CALL HND_hnderr(3,errmsg)
         endif
c%%
c        pnt(  1,ipnt)=wrd1
c%%
         do 170 i=1,3
  170       pnt(i+1,ipnt)=cm(i)
         pnttyp(ipnt)=wrd1
         go to 40
      else
         write(iw,9996) word1
         CALL HND_hnderr(3,errmsg)
      endif
c
c     ----- print distributed points -----
c
  200 continue
      npts=ipnt
      write(iw,9989)

      do 210 i=1,nat
  210    write(iw,9995) i,anam(i),(c(j,i),j=1,3)

      write(iw,9998)
      do 220 i=1,npts
  220    write(iw,9995) i,(pnt(j,i),j=1,4),pnttyp(i)
c
c     ----- perform Stone's multipole analysis -----
c
      CALL HND_stnint(x(i10),dbug)
c
c     ----- calculate electrostatic potential -----
c              from distributed multipoles
c                     on a 3D-grid.
c
      CALL HND_stnmap(dbug)
c
      return
 9999 format(/,10x,30(1h-),
     1       /,10x,'DISTRIBUTED MULTIPOLE ANALYSIS',
     2       /,10x,30(1h-))
 9998 format(/,'     EXPANSION POINTS',
     1       /,'     ----------------',
     2       /,'   #  name    ',12x,'x',17x,'y',17x,'z',7x,'type',/)
 9997 format(' ----- no -$stone- input found.',
     1       ' atomic centers used as expansion points. -----')
 9996 format(' ----- error in -stonex-. illegal expansion point ,',
     1       ' type = ',a8)
 9995 format(1x,i3,2x,a8,3(2x,f16.7),2x,a8)
 9994 format(a8)
 9993 format(' input error in -$stone- . type,iatm,jatm = ',a8,2i5)
 9992 format(' bypass -distributed multipole analysis- .')
 9991 format(' too many expansions centers for',
     1       ' -distributed multipole analysis- . stop.',/,
     2       ' ipnt , maxpnt = ',2i5)
 9989 format(/,'     ATOMIC CENTERS',
     1       /,'     --------------',
     2       /,'   #  name    ',12x,'x',17x,'y',17x,'z',/)
      END
      SUBROUTINE HND_STNOLP(D,S,X,Y,Z,C,IAT,JAT)
      implicit DOUBLE PRECISION (a-h,o-z)
c
c     ----- determine center of overlap distribution -----
c
      logical iandj
      COMMON/HND_nshel/
     1ex(2048),cs(2048),cp(2048),cd(2048),cf(2048),cg(2048),
     1kstart(512),katom(512),ktype(512),kng(512),kloc(512),kmin(512),
     2kmax(512),nshell
      dimension x(*),y(*),z(*),d(*),s(*),c(3)
      data zero,one,two /0.0d+00,1.0d+00,2.0d+00/
c
      cx=zero
      cy=zero
      cz=zero
      ss=zero
      ix=max0(iat,jat)
      jx=min0(iat,jat)
c
c     ----- ishell -----
c
      do 9000 ii=1,nshell
      i=katom(ii)
      if(i.ne.ix)go to 9000
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell -----
c
      do 8000 jj=1,ii
      j=katom(jj)
      if(j.ne.jx)go to 8000
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
      iandj=ii.eq.jj
      max=maxj
      do 7500 i=mini,maxi
      li=loci+i
      in=(li*(li-1))/2
      if(iandj) max=i
      do 7500 j=minj,max
      lj=locj+j
      jn=lj+in
      ss=ss+s(jn)*d(jn)
      cx=cx+x(jn)*d(jn)
      cy=cy+y(jn)*d(jn)
      cz=cz+z(jn)*d(jn)
 7500 continue
 8000 continue
 9000 continue
      c(1)=cx/ss
      c(2)=cy/ss
      c(3)=cz/ss
      return
      END
      SUBROUTINE HND_STNINT(DMAT,DBUG)
      implicit DOUBLE PRECISION (a-h,o-z)
#include "hnd_tol.fh"
#include "nwc_const.fh"
      parameter (nummx=3072)
      parameter (MXATOM=NW_MAX_ATOM)
      logical iandj,out,norm,double,dbug
      character*8 anam,bflab
      character*2 bnam
      COMMON/HND_output/nprint
      COMMON/HND_iofile/ir,iw,ip
      COMMON/HND_runlab/title(10)
      COMMON/HND_mollab/anam(mxatom),bnam(mxatom),bflab(nummx)
      COMMON/HND_infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      COMMON/HND_nshel/
     1ex(2048),cs(2048),cp(2048),cd(2048),cf(2048),cg(2048),
     1kstart(512),katom(512),ktype(512),kng(512),kloc(512),kmin(512),
     2kmax(512),nshell
      COMMON/HND_stnpnt/pnt(4,256),pnttyp(256),npts
      COMMON/HND_stndma/chg(256),x(256),y(256),z(256),
     1      xx(256), yy(256), zz(256), xy(256), xz(256), yz(256),
     2     xxx(256),yyy(256),zzz(256),xxy(256),xxz(256),xyy(256),
     3     yyz(256),xzz(256),yzz(256),xyz(256)
      COMMON/HND_xyzstn/tx0,ty0,tz0,tx1,ty1,tz1,tx2,ty2,tz2,tx3,ty3,tz3,
     1 t,px0,py0,pz0,pxi,pyi,pzi,pxj,pyj,pzj,ni,nj
      dimension dij(225),ijx(225),ijy(225),ijz(225)
      dimension x0(25),y0(25),z0(25)
      dimension x1(25),y1(25),z1(25)
      dimension x2(25),y2(25),z2(25)
      dimension x3(25),y3(25),z3(25)
      dimension jx(35),jy(35),jz(35),ix(35),iy(35),iz(35)
      dimension keq(256)
      dimension dmat(*)
      data thresh    /1.0d-06/
      data zero,one  /0.0d+00,1.0d+00/
      data two,three /2.0d+00,3.0d+00/
      data five      /5.0d+00/
      data sqrt3     /1.73205080756888d+00/
      data sqrt5     /2.23606797749979d+00/
      data rln10     /2.30258d+00/
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     2          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     3          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     1         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     2         21, 1, 1,16,16, 6, 1, 6, 1,11,
     3         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     2          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     3          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     1          1,16, 1, 6, 1,11,11, 1, 6, 6,
     2          1,21, 1, 6, 1,16,16, 1, 6,11,
     3          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     2          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     3          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     1          1, 1,16, 1, 6, 1, 6,11,11, 6,
     2          1, 1,21, 1, 6, 1, 6,16,16, 1,
     3         11,11, 6, 6,11/
c
      tol=rln10*itol
      out=nprint.eq.6
      norm=normf.ne.1.or.normp.ne.1
c
c     ----- zero moment arrays -----
c
      do 10 i=1,npts
         chg(i)=zero
         x(i)=zero
         y(i)=zero
         z(i)=zero
         xx(i)=zero
         yy(i)=zero
         zz(i)=zero
         xy(i)=zero
         xz(i)=zero
         yz(i)=zero
         xxx(i)=zero
         yyy(i)=zero
         zzz(i)=zero
         xxy(i)=zero
         xxz(i)=zero
         xyy(i)=zero
         yyz(i)=zero
         xzz(i)=zero
         yzz(i)=zero
         xyz(i)=zero
   10    continue
c
c     ----- read in density matrix -----
c
      nn=0
      do 30 i=1,num
         do 20 j=1,i
            nn=nn+1
   20       dmat(nn)=two*dmat(nn)
   30    dmat(nn)=dmat(nn)/two
c
c     ----- distributed multipole analysis -----
c
c     ----- ishell -----
c
      do 9000 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell -----
c
      do 8000 jj=1,ii
      if(dbug) then
         write(iw,9277) ii,jj
      endif
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      iandj=ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ij=0
      max=maxj
      do 40 i=mini,maxi
      nx=ix(i)
      ny=iy(i)
      nz=iz(i)
      if(iandj) max=i
      do 40 j=minj,max
      ij=ij+1
      ijx(ij)=nx+jx(j)
      ijy(ij)=ny+jy(j)
      ijz(ij)=nz+jz(j)
   40 continue
c
c     ----- i primitive -----
c
      jgmax=j2
      do 7000 ig=i1,i2
      igp=ig-i1+1
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
c
c     ----- j primitive -----
c
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      jgp=jg-j1+1
      aj=ex(jg)
      aa=ai+aj
      dum=aj*arri/aa
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      ax=(axi+aj*xj)/aa
      ay=(ayi+aj*yj)/aa
      az=(azi+aj*zj)/aa
c
c     ----- data for translation to chosen centers -----
c
      ieq=1
      keq(1)=1
      ceq=one
      ipt=1
      if(npts.ge.2) then
         dist=(ax-pnt(2,1))**2+(ay-pnt(3,1))**2+(az-pnt(4,1))**2
         do 45 npt=2,npts
         test=(ax-pnt(2,npt))**2+(ay-pnt(3,npt))**2+(az-pnt(4,npt))**2
         if(test.gt.dist)go to 45
         dist=test
         ipt=npt
   45    continue
c
c     ----- equivalent centers for translation -----
c
         keq(1)=ipt
         ieq=1
         do 50 npt=1,npts
         if(npt.eq.ipt)go to 50
         test=(ax-pnt(2,npt))**2+(ay-pnt(3,npt))**2+(az-pnt(4,npt))**2
         if( abs(test-dist).gt.thresh)go to 50
         ieq=ieq+1
         keq(ieq)=npt
   50    continue
         ceq=one/dble(ieq)
      endif
c
c     ----- density factor -----
c
      double=iandj.and.ig.ne.jg
      max=maxj
      nn=0
      do 310 i=mini,maxi
      go to ( 70, 80,180,180, 90,180,180,100,180,180,
     1       110,180,180,120,180,180,180,180,180,130),i
   70 dum1=csi*fac
      go to 180
   80 dum1=cpi*fac
      go to 180
   90 dum1=cdi*fac
      go to 180
  100 if(norm) dum1=dum1*sqrt3
      go to 180
  110 dum1=cfi*fac
      go to 180
  120 if(norm) dum1=dum1*sqrt5
      go to 180
  130 if(norm) dum1=dum1*sqrt3
  180 if(iandj) max=i
      do 310 j=minj,max
      go to (190,200,300,300,210,300,300,220,300,300,
     1       230,300,300,240,300,300,300,300,300,250),j
  190 dum2=dum1*csj
      if(.not.double) go to 300
      if(i.gt.1) go to 195
      dum2=dum2+dum2
      go to 300
  195 dum2=dum2+csi*cpj*fac
      go to 300
  200 dum2=dum1*cpj
      if(double) dum2=dum2+dum2
      go to 300
  210 dum2=dum1*cdj
      if(double) dum2=dum2+dum2
      go to 300
  220 if(norm) dum2=dum2*sqrt3
      go to 300
  230 dum2=dum1*cfj
      if(double) dum2=dum2+dum2
      go to 300
  240 if(norm) dum2=dum2*sqrt5
      go to 300
  250 if(norm) dum2=dum2*sqrt3
  300 nn=nn+1
  310 dij(nn)=dum2
c
c     ----- moment integrals -----
c
      t= sqrt(aa)
c
c     ----- translate everything to center of charge distribution -----
c
      px0=zero
      py0=zero
      pz0=zero
      pxi=xi-ax
      pyi=yi-ay
      pzi=zi-az
      pxj=xj-ax
      pyj=yj-ay
      pzj=zj-az
      in=-5
      do 340 i=1,lit
      in=in+5
      ni=i
      do 340 j=1,ljt
      jn=in+j
      nj=j
      CALL HND_stnxyz
      x0(jn)=tx0/t
      y0(jn)=ty0/t
      z0(jn)=tz0/t
      x1(jn)=tx1/t
      y1(jn)=ty1/t
      z1(jn)=tz1/t
      x2(jn)=tx2/t
      y2(jn)=ty2/t
      z2(jn)=tz2/t
      x3(jn)=tx3/t
      y3(jn)=ty3/t
      z3(jn)=tz3/t
  340 continue
c
      max=maxj
      nn=0
      do 500 i=mini,maxi
      li=loci+i
      in=(li*(li-1))/2
      if(iandj)max=i
      do 500 j=minj,max
      lj=locj+j
      if(dbug) then
         write(iw,9278) li,igp,lj,jgp
      endif
      jn=lj+in
      nn=nn+1
      nx=ijx(nn)
      ny=ijy(nn)
      nz=ijz(nn)
      dm=-dij(nn)*dmat(jn)
      ts=dm*x0(nx)*y0(ny)*z0(nz)
      tx=dm*x1(nx)*y0(ny)*z0(nz)
      ty=dm*x0(nx)*y1(ny)*z0(nz)
      tz=dm*x0(nx)*y0(ny)*z1(nz)
      txx=dm*x2(nx)*y0(ny)*z0(nz)
      tyy=dm*x0(nx)*y2(ny)*z0(nz)
      tzz=dm*x0(nx)*y0(ny)*z2(nz)
      txy=dm*x1(nx)*y1(ny)*z0(nz)
      txz=dm*x1(nx)*y0(ny)*z1(nz)
      tyz=dm*x0(nx)*y1(ny)*z1(nz)
      txxx=dm*x3(nx)*y0(ny)*z0(nz)
      tyyy=dm*x0(nx)*y3(ny)*z0(nz)
      tzzz=dm*x0(nx)*y0(ny)*z3(nz)
      txxy=dm*x2(nx)*y1(ny)*z0(nz)
      txxz=dm*x2(nx)*y0(ny)*z1(nz)
      txyy=dm*x1(nx)*y2(ny)*z0(nz)
      tyyz=dm*x0(nx)*y2(ny)*z1(nz)
      txzz=dm*x1(nx)*y0(ny)*z2(nz)
      tyzz=dm*x0(nx)*y1(ny)*z2(nz)
      txyz=dm*x1(nx)*y1(ny)*z1(nz)
      do 490 neq=1,ieq
      iq=keq(neq)
      cx=ax-pnt(2,iq)
      cy=ay-pnt(3,iq)
      cz=az-pnt(4,iq)
      fchg=ceq*ts
      chg(iq)=chg(iq)+fchg
      if(dbug) then
         write(iw,9279) iq,dij(nn),dmat(jn),fchg
      endif
      x(iq)=x(iq)+ceq*(tx+cx*ts)
      y(iq)=y(iq)+ceq*(ty+cy*ts)
      z(iq)=z(iq)+ceq*(tz+cz*ts)
      xx(iq)=xx(iq)+ceq*(txx+two*cx*tx+cx*cx*ts)
      yy(iq)=yy(iq)+ceq*(tyy+two*cy*ty+cy*cy*ts)
      zz(iq)=zz(iq)+ceq*(tzz+two*cz*tz+cz*cz*ts)
      xy(iq)=xy(iq)+ceq*(txy+cx*ty+tx*cy+cx*cy*ts)
      xz(iq)=xz(iq)+ceq*(txz+cx*tz+tx*cz+cx*cz*ts)
      yz(iq)=yz(iq)+ceq*(tyz+cy*tz+ty*cz+cy*cz*ts)
      xxx(iq)=xxx(iq)+ceq*(txxx+three*cx*txx+three*cx*cx*tx+cx*cx*cx*ts)
      yyy(iq)=yyy(iq)+ceq*(tyyy+three*cy*tyy+three*cy*cy*ty+cy*cy*cy*ts)
      zzz(iq)=zzz(iq)+ceq*(tzzz+three*cz*tzz+three*cz*cz*tz+cz*cz*cz*ts)
      xxy(iq)=xxy(iq)+ceq*(txxy+two*cx*txy+cx*cx*ty+
     1              txx*cy+two*tx*cx*cy+cx*cx*cy*ts)
      xxz(iq)=xxz(iq)+ceq*(txxz+two*cx*txz+cx*cx*tz+
     1              txx*cz+two*tx*cx*cz+cx*cx*cz*ts)
      xyy(iq)=xyy(iq)+ceq*(txyy+two*txy*cy+tx*cy*cy+
     1              cx*tyy+two*cx*cy*ty+cx*cy*cy*ts)
      yyz(iq)=yyz(iq)+ceq*(tyyz+two*cy*tyz+cy*cy*tz+
     1              tyy*cz+two*ty*cy*cz+cy*cy*cz*ts)
      xzz(iq)=xzz(iq)+ceq*(txzz+two*txz*cz+tx*cz*cz+
     1              cx*tzz+two*cx*cz*tz+cx*cz*cz*ts)
      yzz(iq)=yzz(iq)+ceq*(tyzz+two*tyz*cz+ty*cz*cz+
     1              cy*tzz+two*cy*cz*tz+cy*cz*cz*ts)
      xyz(iq)=xyz(iq)+ceq*(txyz+txy*cz+txz*cy+tx*cy*cz+
     1              cx*tyz+cx*ty*cz+cx*cy*tz+cx*cy*cz*ts)
  490 continue
  500 continue
c
 6000 continue
 7000 continue
c
 8000 continue
 9000 continue
c
c     ----- print and punch the results -----
c
      write(iw,9200)
      write(ip,9201)
c%%
C     do 9202 i=1,nat
C     write(iw,9211)i,anam(i),zan(i)
C     write(ip,9216)anam(i),zan(i),(c(j,i),j=1,3)
C9202 continue
c%%
      write(iw,9203)
      do 9205 i=1,npts
      write(iw,9210)i,pnt(1,i),chg(i)
      write(ip,9215)pnt(1,i),chg(i),(pnt(j,i),j=2,4)
 9205 continue
      write(ip,9206)
c
      write(iw,9220)
      write(ip,9221)
      do 9225 i=1,npts
      write(iw,9230)i,pnt(1,i),x(i),y(i),z(i)
      write(ip,9215)(pnt(j,i),j=1,4)
      write(ip,9236)x(i),y(i),z(i)
 9225 continue
      write(ip,9206)
c
      write(iw,9240)
      write(ip,9241)
      do 9245 i=1,npts
      write(iw,9250)i,pnt(1,i),xx(i),yy(i),zz(i),xy(i),xz(i),yz(i)
      write(ip,9215)(pnt(j,i),j=1,4)
      write(ip,9256)xx(i),yy(i),zz(i),xy(i),xz(i),yz(i)
 9245 continue
      write(ip,9206)
      write(iw,9242)
      do 9246 i=1,npts
         xxi=xx(i)
         yyi=yy(i)
         zzi=zz(i)
         xyi=xy(i)
         xzi=xz(i)
         yzi=yz(i)
         dum=xxi+yyi+zzi
         xxi=three*xxi-dum
         yyi=three*yyi-dum
         zzi=three*zzi-dum
         xyi=three*xyi
         xzi=three*xzi
         yzi=three*yzi
      write(iw,9250)i,pnt(1,i),xxi,yyi,zzi,xyi,xzi,yzi
 9246 continue
c
      write(iw,9260)
      write(ip,9261)
      do 9265 i=1,npts
      write(iw,9270)i,pnt(1,i),xxx(i),yyy(i),zzz(i),xxy(i),xxz(i),
     1                         xyy(i),yyz(i),xzz(i),yzz(i),xyz(i)
      write(ip,9215)(pnt(j,i),j=1,4)
      write(ip,9276)xxx(i),yyy(i),zzz(i),xxy(i),xxz(i),xyy(i),
     1              yyz(i),xzz(i),yzz(i),xyz(i)
 9265 continue
      write(ip,9206)
      write(iw,9262)
      do 9266 i=1,npts
         xxxi=xxx(i)
         yyyi=yyy(i)
         zzzi=zzz(i)
         xxyi=xxy(i)
         xxzi=xxz(i)
         xyyi=xyy(i)
         yyzi=yyz(i)
         xzzi=xzz(i)
         yzzi=yzz(i)
         xyzi=xyz(i)
         dumx=xxxi+xyyi+xzzi
         dumy=xxyi+yyyi+yzzi
         dumz=xxzi+yyzi+zzzi
         xxxi=five*xxxi-three*dumx
         yyyi=five*yyyi-three*dumy
         zzzi=five*zzzi-three*dumz
         xxyi=five*xxyi-three*dumy
         xxzi=five*xxzi-three*dumz
         xyyi=five*xyyi-three*dumx
         yyzi=five*yyzi-three*dumz
         xzzi=five*xzzi-three*dumx
         yzzi=five*yzzi-three*dumy
         xyzi=five*xyzi
      write(iw,9270)i,pnt(1,i),xxxi,yyyi,zzzi,xxyi,xxzi,
     1                         xyyi,yyzi,xzzi,yzzi,xyzi
 9266 continue
c
      return
c
 9200 format(/'  NET CHARGES AT  ATOMS',
     1       /'  ---------------------',
     1       /'      name       charge',/)
 9201 format(' $EFC          stone analysis')
 9203 format(/'  NET CHARGES AT POINTS',
     1       /'  ---------------------',
     1       /'      name       charge',/)
 9206 format(' $END   ')
 9210 format(1x,i3,1x,a8,2x,f9.5)
 9211 format(1x,i3,1x,'z',a8,1x,f9.5)
 9215 format(a8,4(2x,f10.6))
 9216 format('z',a8,1x,4(f10.6,2x))
 9220 format(/'  FIRST MOMENTS AT POINTS',
     1       /'  -----------------------',
     1       /'      name        x           y           z',/)
 9221 format(' $EFD          stone analysis')
 9230 format(1x,i3,1x,a8,2x,f9.5,2x,f9.5,2x,f9.5)
 9236 format(3(2x,f10.6))
 9240 format(/'  SECOND MOMENTS AT POINTS',
     1       /'  ------------------------',
c%%
cd   2       /'      name        xx         yy         zz         xy'
     2       /'      name        xx         yy         zz         xy',
c%%
     3        '         xz         yz',/)
 9241 format(' $EFQ          stone analysis')
 9242 format(/'  QUADRUPOLE MOMENTS AT POINTS',
     1       /'  ----------------------------',
c%%
cd   2       /'      name        xx         yy         zz         xy'
     2       /'      name        xx         yy         zz         xy',
c%%
     3        '         xz         yz',/)
 9250 format(1x,i3,1x,a8,2x,f9.5,2x,f9.5,2x,f9.5,
     1                   2x,f9.5,2x,f9.5,2x,f9.5)
 9256 format(6(2x,f10.6))
 9260 format(/'  THIRD MOMENTS AT POINTS',
     1       /'  -----------------------',
     1       /'      name        xxx        yyy        zzz        xxy',
     1        '        xxz        xyy',
     1        '        yyz        xzz        yzz        xyz',/)
 9261 format(' $EFO          stone analysis')
 9262 format(/'  OCTUPOLE MOMENTS AT POINTS',
     1       /'  --------------------------',
     1       /'      name        xxx        yyy        zzz        xxy',
     1        '        xxz        xyy',
     1        '        yyz        xzz        yzz        xyz',/)
 9270 format(1x,i3,1x,a8,2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5,
     1                   2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)
 9276 format(6(2x,f10.6))
 9277 format(' ishell ',i3,'  jshell ',i3)
 9278 format('    func ',i2,'  prim ',i2,'   func ',i2,'  prim ',i2)
 9279 format('    cent,dij,dmat,chg = ',i2,3f12.5)
      END
      SUBROUTINE HND_STNXYZ
      implicit DOUBLE PRECISION (a-h,o-z)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
#include "hnd_whermt.fh"
c
      COMMON/HND_xyzstn/tx0,ty0,tz0,tx1,ty1,tz1,tx2,ty2,tz2,tx3,ty3,tz3,
     1              t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension min(6),max(6)
      data min /1,2,4,7,11,16/
      data max /1,3,6,10,15,21/
      data zero /0.0d+00/
      tx0=zero
      ty0=zero
      tz0=zero
      tx1=zero
      ty1=zero
      tz1=zero
      tx2=zero
      ty2=zero
      tz2=zero
      tx3=zero
      ty3=zero
      tz3=zero
      npts=(ni+nj-2+3)/2+1
      imin=min(npts)
      imax=max(npts)
      do 13 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)/t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 go to (12,11,10,9,8,7,6),nj
    6 px=px*bx
      py=py*by
      pz=pz*bz
    7 px=px*bx
      py=py*by
      pz=pz*bz
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 continue
      tx0=tx0+px
      ty0=ty0+py
      tz0=tz0+pz
      tx1=tx1+px*ptx
      ty1=ty1+py*pty
      tz1=tz1+pz*ptz
      tx2=tx2+px*ptx*ptx
      ty2=ty2+py*pty*pty
      tz2=tz2+pz*ptz*ptz
      tx3=tx3+px*ptx*ptx*ptx
      ty3=ty3+py*pty*pty*pty
      tz3=tz3+pz*ptz*ptz*ptz
   13 continue
      return
      END
      SUBROUTINE HND_STNMAP(DBUG)
      implicit DOUBLE PRECISION (a-h,o-z)
      character*8 errmsg
      logical some,dbug
      COMMON/HND_memory/maxcor,maxlcm
      COMMON/HND_output/nprint
      COMMON/HND_iofile/ir,iw,ip
      COMMON/HND_grdpar/nxyzg,nxyg,nxg,igunit,ux(3),uy(3),uz(3),xyzg0(3)
      COMMON/HND_scm/x(1)
      dimension mo(100),moa(100),mob(100)
      dimension xyzmap(3,100)
      dimension xrange(3),yrange(3),zrange(3),xvec(3),yvec(3),zvec(3)
      dimension origin(3)
      dimension iatom(4)
      dimension errmsg(3)
      equivalence (mo(1),moa(1))
      data errmsg /'program ','stop in ','-STNMAP-'/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data tenm3  /1.0d-03/
      data ten12  /1.0d+12/
      data moa    /100*0/
      data mob    /100*0/
      data xyzmap /300*1.0d+12/
      data nptxyz /100/
      data iunit  /-1/
      data bohr   /0.52917715d+00/
      data xrange /0.0d+00,0.0d+00,1.0d+00/
      data yrange /0.0d+00,0.0d+00,1.0d+00/
      data zrange /0.0d+00,0.0d+00,1.0d+00/
      data xvec   /1.0d+00,0.0d+00,0.0d+00/
      data yvec   /0.0d+00,1.0d+00,0.0d+00/
      data zvec   /0.0d+00,0.0d+00,1.0d+00/
      data origin /0.0d+00,0.0d+00,0.0d+00/
      data iatom  /0,0,0,0/
c
      namelist /map/ iunit,xyzmap,xrange,yrange,zrange,xvec,yvec,zvec,
     1               origin,iatom,mo,moa,mob
c
      rewind ir
      read(ir,map,end=10,err=10)
  10  continue
c
c     ----- calculate electrostatic potential from -----
c           distributed multipole analysis.
c
      write(iw,9999)
      write(iw,9994)
      some=nprint.ne.-5
c
c     ----- define points for calculation -----
c           1. grid points
c           2. specific points -xyzmap-
c
      nptmax=0
      CALL HND_grdmap(x(1+3*nptmax),nptgrd,nxgrd,nygrd,nzgrd,
     1            xrange,yrange,zrange,xvec,yvec,zvec,
     2            origin,iatom,iunit)
      nptmax=nptmax+nptgrd
c
      if(nptgrd.gt.0) then
         write(iw,9985) nxgrd,nygrd,nzgrd,
     1    iatom,origin,xvec,xrange(1),xrange(2),xrange(3),
     2                 yvec,yrange(1),yrange(2),yrange(3),
     3                 zvec,zrange(1),zrange(2),zrange(3)
      else
         write(iw,9986)
         return
      endif
c
      fac=one
      if(iunit.gt.0) fac=bohr
      npt=0
      do 20 ipt=1,nptxyz
      if(xyzmap(1,ipt).eq.ten12.or.xyzmap(2,ipt).eq.ten12.or.
     1   xyzmap(3,ipt).eq.ten12) go to 20
      npt=npt+1
      x(1+3*(npt-1+nptmax))=xyzmap(1,ipt)/fac
      x(2+3*(npt-1+nptmax))=xyzmap(2,ipt)/fac
      x(3+3*(npt-1+nptmax))=xyzmap(3,ipt)/fac
   20 continue
      nptmax=nptmax+npt
c
c     ----- set pointers for partitioning of memory -----
c
      i00 =          1
      i10 = i00    + nptmax*3
      last= i10    + nptmax
      need= last   - 1
c
c     ------ get core memory ------
c
c%%
C     CALL HND_cmem(loadcm)
C     if(need.le.maxcor) go to 110
C     write(iw,9998)
C     CALL HND_needcm(last,need)
C     CALL HND_hnderr(3,errmsg)
C 110 continue
c%%
c
c     ----- electrostatic potential from distributed multipoles -----
c
      write(iw,9997)
c
      write(ip,9991)
      write(ip,9990) xrange(1),xrange(2),xrange(3),nxgrd
      write(ip,9990) yrange(1),yrange(2),yrange(3),nygrd
      write(ip,9990) zrange(1),zrange(2),zrange(3),nzgrd
c
      do 120  ipt=1,nptmax
         xp = x(1+3*(ipt-1)+i00-1)
         yp = x(2+3*(ipt-1)+i00-1)
         zp = x(3+3*(ipt-1)+i00-1)
         elpott = hnd_stnesp(xp,yp,zp)
         x(ipt+i10-1)=elpott
         if(ipt.gt.nptgrd.or.dbug)
     1   write(iw,9995) ipt,xp,yp,zp,elpott
  120    continue
      if(nptgrd.gt.0) write(ip,9989) (x(ipt+i10-1),ipt=1,nptgrd)
      write(ip,9988)
c
c     -----  reset core memory -----
c
c%%
C     CALL HND_cmem(ngotcm)
C     if(ngotcm.ne.loadcm) CALL HND_setc(loadcm)
c%%
      return
 9999 format(/,10x,55(1h-),
     1       /,10x,'ELECTROSTATIC POTENTIAL MAP ',
     2             'FROM DISTRIBUTED MULTIPOLES',
     3       /,10x,55(1h-),/)
 9998 format(' not enough core in -stnmap-')
 9997 format(3x,'point',6x,'x',9x,'y',9x,'z',5x,'potential(a.u.)')
 9995 format(1x,i5,2x,3f10.5,f15.6)
 9994 format(' 1 a.u. = 9.07618 esu/cm ( or statvolts ) ')
 9991 format(' $STNMAP')
 9990 format(3f13.6,i5)
 9989 format(6f13.6,2x)
 9988 format(' $END   ')
 9986 format(/,' no grid defined. return. ')
 9985 format(/,1x,'( ',i4,' by ',i4,' by ',i4,' ) grid defined.',
     1 /,' iatom  = ',4i4,
     2 /,' origin = ',3f10.5,
     3 /,' xvec   = ',3f10.5,' xmin=',f8.2,' xmax=',f8.2,' xstep=',f8.2,
     4 /,' yvec   = ',3f10.5,' ymin=',f8.2,' ymax=',f8.2,' ystep=',f8.2,
     5 /,' zvec   = ',3f10.5,' zmin=',f8.2,' zmax=',f8.2,' zstep=',f8.2,
     6 /)
      END
      FUNCTION HND_STNESP(CX,CY,CZ)
      implicit DOUBLE PRECISION (a-h,o-z)
#include "nwc_const.fh"
      parameter (MXATOM=NW_MAX_ATOM)
      COMMON/HND_iofile/ir,iw,ip
      COMMON/HND_infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      COMMON/HND_stnpnt/pnt(4,256),pntyp(256),npts
      COMMON/HND_stndma/chg(256),efd(3,256),efq(6,256),efo(10,256)
      data zero,pt5   /0.0D+00,0.5D+00/
      data one,six    /1.0D+00,6.0D+00/
      data three,five /3.0D+00,5.0D+00/
      data two        /2.0D+00/
      data cscren     /1.0d+00/
      data ascren     /0.5d+00/
      data tenm3      /1.0d-03/
c
      epot=zero
c
c     -- distributed multipole contribution ( note damping factor ) --
c
      do 10 ipt=1,npts
      xi=-pnt(2,ipt)+cx
      yi=-pnt(3,ipt)+cy
      zi=-pnt(4,ipt)+cz
      rr=xi*xi+yi*yi+zi*zi
      r1= sqrt(rr)
      if(r1.le.tenm3) then
         write(iw,9999) cx,cy,cz,ipt
         go to 10
      endif
      r3=rr*r1
      r5=rr*r3
      r7=rr*r5
      damp=cscren* exp(-ascren*rr)
c
c     ----- contribution from -charge- -----
c
      epot=epot+(one-damp)*chg(ipt)/r1
c
c     ----- contribution from -dipole- -----
c
      term=xi*efd(1,ipt)+yi*efd(2,ipt)+zi*efd(3,ipt)
      epot=epot+(one-damp)*term/r3
c
c     ----- contribution from -quadrupole- -----
c
      xx=efq(1,ipt)
      yy=efq(2,ipt)
      zz=efq(3,ipt)
      xy=efq(4,ipt)
      xz=efq(5,ipt)
      yz=efq(6,ipt)
      xxi=xi*xi
      yyi=yi*yi
      zzi=zi*zi
      xyi=xi*yi
      xzi=xi*zi
      yzi=yi*zi
      term=zero
      term=term+xx*(three*xxi-rr)
      term=term+yy*(three*yyi-rr)
      term=term+zz*(three*zzi-rr)
      term=term+xy*three*xyi*two
      term=term+xz*three*xzi*two
      term=term+yz*three*yzi*two
      epot=epot+(one-damp)*pt5*term/(three*r5)
c
c     ----- contribution from -octupole- -----
c
      xxx=efo(1,ipt)
      yyy=efo(2,ipt)
      zzz=efo(3,ipt)
      xxy=efo(4,ipt)
      xxz=efo(5,ipt)
      xyy=efo(6,ipt)
      yyz=efo(7,ipt)
      xzz=efo(8,ipt)
      yzz=efo(9,ipt)
      xyz=efo(10,ipt)
      xxxi=xxi*xi
      yyyi=yyi*yi
      zzzi=zzi*zi
      xxyi=xxi*yi
      xxzi=xxi*zi
      xyyi=xi*yyi
      yyzi=yyi*zi
      xzzi=xi*zzi
      yzzi=yi*zzi
      xyzi=xi*yi*zi
      rrxi=rr*xi
      rryi=rr*yi
      rrzi=rr*zi
      term=zero
      term=term+xxx*(five*xxxi-three*rrxi)
      term=term+yyy*(five*yyyi-three*rryi)
      term=term+zzz*(five*zzzi-three*rrzi)
      term=term+xxy*(five*xxyi-rryi)*three
      term=term+xxz*(five*xxzi-rrzi)*three
      term=term+xyy*(five*xyyi-rrxi)*three
      term=term+yyz*(five*yyzi-rrzi)*three
      term=term+xzz*(five*xzzi-rrxi)*three
      term=term+yzz*(five*yzzi-rryi)*three
      term=term+xyz*(five*xyzi)*six
      epot=epot+(one-damp)*pt5*term/(five*r7)
   10 continue
c
c     ----- nuclear contribution ( note damping factor ) -----
c
      do 20 iat=1,nat
      xi=-c(1,iat)+cx
      yi=-c(2,iat)+cy
      zi=-c(3,iat)+cz
      rr=xi*xi+yi*yi+zi*zi
      r1= sqrt(rr)
      if(r1.le.tenm3) then
         write(iw,9998) cx,cy,cz,iat
         go to 20
      endif
      damp=cscren* exp(-ascren*rr)
      epot=epot+(one-damp)*zan(iat)/r1
   20 continue
c
      hnd_stnesp=epot
      return
 9999 format(' --- warning - electrostatic potential at ',
     1 3F10.5,' . contribution from multipole ',I3,' ignored')
 9998 format(' --- warning - electrostatic potential at ',
     1 3F10.5,' . contribution from nucleus   ',I3,' ignored')
      END
C
C     ----- SUBROUTINES BELOW TAKEN FROM CTL.F IN HONDO -----
C
      SUBROUTINE HND_PREVL(V,E,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "global.fh"
#include "nwc_const.fh"
C
C     ----- PRINT OUT E AND V-MATRICES
C
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXBFN=3072)
      character*8  anam,bflab,cbf
      character*2  bnam
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_MOLXYZ/CC(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      DIMENSION V(NDIM,1),E(1)
      DIMENSION IC(5),CBF(5),C(5)
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      if (ga_nodeid().eq.0) then 
         WRITE (IW,9008)
         IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
         IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
         IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
         WRITE (IW,9008)
         IF(LIST.EQ.0) WRITE (IW,9068) (E(I),I = IMIN,IMAX)
         IF(LIST.EQ.1) WRITE (IW,9168) (E(I),I = IMIN,IMAX)
         IF(LIST.EQ.2) WRITE (IW,9268) (E(I),I = IMIN,IMAX)
         WRITE (IW,9008)
         DO 120 J = 1,N
         IF(LIST.EQ.0) WRITE (IW,9048)
     &      J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
         IF(LIST.EQ.1) WRITE (IW,9148)
     &      J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
         IF(LIST.EQ.2) WRITE (IW,9248) 
     &     J,BFLAB(J),(V(J,I),I = IMIN,IMAX)

  120    CONTINUE
      endif
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      if (ga_nodeid().eq.0) WRITE(IW,9328) J,E(J)
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I

      CBF(II)=BFLAB(I)

      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      if (ga_nodeid().eq.0) WRITE(IW,9348)
     &   (IC(IDUM),CBF(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      if (ga_nodeid().eq.0) WRITE(IW,9348)
     &   (IC(IDUM),CBF(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(15X,10(3X,I4,4X))
 9048 FORMAT(I5,2X,A8,10F11.5)
 9068 FORMAT(15X,10F11.5)
 9128 FORMAT(15X,7(5X,I4,6X))
 9148 FORMAT(I5,2X,A8,7F15.10)
 9168 FORMAT(15X,7F15.10)
 9228 FORMAT(15X,7(5X,I4,6X))
 9248 FORMAT(I5,2X,A8,7E15.8)
 9268 FORMAT(15X,7E15.8)
 9328 FORMAT(/,' Vec =',I5,' Eig/occ =',F11.5,/,1X,31(1H-))
 9348 FORMAT(5(I5,1X,A8,F11.5))
      END
      SUBROUTINE HND_DIAJAC(A,V,D,N,NDIM,T)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NDIM,NDIM),V(NDIM,NDIM),D(NDIM),T(NDIM,2)

      CALL HND_JACOBI(A,V,D,N,NDIM,NROT,T(1,1),T(1,2))
      CALL HND_JACSRT(V,D,N,NDIM)

      RETURN
      END
      SUBROUTINE HND_SYMOFF
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*5 IRPA
      PARAMETER (MXSYM =48)
      COMMON/HND_SYMREP/IRPNAM(14),IPA(14),LAMBDA(14),LAMBD0(14),
     1              IADDR1(14),IADDR2(14),IADDR3(14)
      COMMON/HND_SYMBLK/NIRRED,NSALC
      COMMON/HND_SYMTRY/INVT(MXSYM),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMSAV/NIRP1,LAB1,NSALC1,NDIM1,NDEG1
      DATA IRPA /'A   '/
C
      IF(NIRRED.EQ.1) RETURN
C
      NIRP1    =NIRRED
      LAB1     =IRPNAM(1)
      NSALC1   =IPA(1)
      NDIM1    =LAMBDA(1)
      NDEG1    =LAMBD0(1)
      NT       =1
      NIRRED   =1
c%%
C     IRPNAM(1)=IRPA
c%%
      IPA(1)   =NSALC
      LAMBDA(1)=1
      LAMBD0(1)=1
C
      RETURN
      END
      SUBROUTINE HND_SYMON
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MXSYM =48)
      COMMON/HND_SYMREP/IRPNAM(14),IPA(14),LAMBDA(14),LAMBD0(14),
     1              IADDR1(14),IADDR2(14),IADDR3(14)
      COMMON/HND_SYMBLK/NIRRED,NSALC
      COMMON/HND_SYMTRY/INVT(MXSYM),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMSAV/NIRP1,LAB1,NSALC1,NDIM1,NDEG1
C
      IF(NIRP1.EQ.1) RETURN
C
      NT       =NTMAX
      NIRRED   =NIRP1
      IRPNAM(1)=LAB1
      IPA(1)   =NSALC1
      LAMBDA(1)=NDIM1
      LAMBD0(1)=NDEG1
C
      RETURN
      END
      SUBROUTINE HND_PUSQL(V,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- PUNCH OUT A SQUARE MATRIX WITH ORDERING LABELS -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION V(NDIM,1)
      DO 120 J = 1,M
      IC = 0
      MAX = 0
  100 MIN = MAX+1
      MAX = MAX+5
      IC = IC+1
      IF (MAX .GT. N) MAX = N
      IF (N   .LE.99) WRITE (IP,9008) J,IC,(V(I,J),I = MIN,MAX)
      IF (N   .GT.99) WRITE (IP,9028) J,IC,(V(I,J),I = MIN,MAX)
      IF (MAX .LT. N) GO TO 100
  120 CONTINUE
      RETURN
 9008 FORMAT(I2,I3,5E15.8)
 9028 FORMAT(I3,I2,5E15.8)
      END
C
C     ----- SUBROUTINES ABOVE TAKEN FROM CTL.F IN HONDO -----
C


C
C     ----- SUBROUTINES BELOW TAKEN FROM SYM.F IN HONDO -----
C
      SUBROUTINE HND_SPHMAT(SP,NDIM,IFINV,IFTSP)
      IMPLICIT REAL*8(A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
      CHARACTER*8 ERRMSG
      LOGICAL DBUG
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_SPHERI/PSHELL(3,3),DSHELL(6,6),FSHELL(10,10),
     1              GSHELL(15,15),
     1              PIHELL(3,3),DIHELL(6,6),FIHELL(10,10),GIHELL(15,15)
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION SP(NDIM,1)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-SPHMAT-'/
C
      DBUG=.FALSE.
C
      DO J=1,NUM
         DO I=1,NUM
            SP(I,J)=ZERO
         ENDDO
         SP(J,J)=ONE
      ENDDO
C
C     ----- BRANCH ON -IFINV- AND -IFTSP- -----
C
      IF(IFINV.EQ.-1.AND.IFTSP.EQ.0) THEN
C
C     ----- INVERSE -----
C
         DO II=1,NSHELL
            LIT =KTYPE(II)
            MIN =KMIN(II)
            MAX =KMAX(II)
            LOC =KLOC(II)
            IF(LIT.EQ.1) THEN
            ELSEIF(LIT.EQ.2) THEN
               IF(MIN.EQ.1) THEN
                  DO J=1,3
                     DO I=1,3
                        SP(I+LOC  ,J+LOC  )=PIHELL(I,J)
                     ENDDO
                  ENDDO
               ELSE
                  DO J=1,3
                     DO I=1,3
                        SP(I+LOC-1,J+LOC-1)=PIHELL(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ELSEIF(LIT.EQ.3) THEN
               DO J=1,6
                  DO I=1,6
                     SP(I+LOC-1,J+LOC-1)=DIHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.4) THEN
               DO J=1,10
                  DO I=1,10
                     SP(I+LOC-1,J+LOC-1)=FIHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.5) THEN
               DO J=1,15
                  DO I=1,15
                     SP(I+LOC-1,J+LOC-1)=GIHELL(I,J)
                  ENDDO
               ENDDO
            ELSE
               WRITE(IW,9999)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ENDDO
C
      ELSEIF(IFINV.EQ.-1.AND.IFTSP.NE.0) THEN
C
C     ----- INVERSE TRANSPOSE -----
C
         DO II=1,NSHELL
            LIT =KTYPE(II)
            MIN =KMIN(II)
            MAX =KMAX(II)
            LOC =KLOC(II)
            IF(LIT.EQ.1) THEN
            ELSEIF(LIT.EQ.2) THEN
               IF(MIN.EQ.1) THEN
                  DO J=1,3
                     DO I=1,3
                        SP(J+LOC  ,I+LOC  )=PIHELL(I,J)
                     ENDDO
                  ENDDO
               ELSE
                  DO J=1,3
                     DO I=1,3
                        SP(J+LOC-1,I+LOC-1)=PIHELL(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ELSEIF(LIT.EQ.3) THEN
               DO J=1,6
                  DO I=1,6
                     SP(J+LOC-1,I+LOC-1)=DIHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.4) THEN
               DO J=1,10
                  DO I=1,10
                     SP(J+LOC-1,I+LOC-1)=FIHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.5) THEN
               DO J=1,15
                  DO I=1,15
                     SP(J+LOC-1,I+LOC-1)=GIHELL(I,J)
                  ENDDO
               ENDDO
            ELSE
               WRITE(IW,9999)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ENDDO
C
      ELSEIF(IFINV.NE.-1.AND.IFTSP.EQ.0) THEN
C
C     -----                   -----
C
         DO II=1,NSHELL
            LIT =KTYPE(II)
            MIN =KMIN(II)
            MAX =KMAX(II)
            LOC =KLOC(II)
            IF(LIT.EQ.1) THEN
            ELSEIF(LIT.EQ.2) THEN
               IF(MIN.EQ.1) THEN
                  DO J=1,3
                     DO I=1,3
                        SP(I+LOC  ,J+LOC  )=PSHELL(I,J)
                     ENDDO
                  ENDDO
               ELSE
                  DO J=1,3
                     DO I=1,3
                        SP(I+LOC-1,J+LOC-1)=PSHELL(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ELSEIF(LIT.EQ.3) THEN
               DO J=1,6
                  DO I=1,6
                     SP(I+LOC-1,J+LOC-1)=DSHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.4) THEN
               DO J=1,10
                  DO I=1,10
                     SP(I+LOC-1,J+LOC-1)=FSHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.5) THEN
               DO J=1,15
                  DO I=1,15
                     SP(I+LOC-1,J+LOC-1)=GSHELL(I,J)
                  ENDDO
               ENDDO
            ELSE
               WRITE(IW,9999)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ENDDO
C
      ELSEIF(IFINV.NE.-1.AND.IFTSP.NE.0) THEN
C
C     -----         TRANSPOSE -----
C
         DO II=1,NSHELL
            LIT =KTYPE(II)
            MIN =KMIN(II)
            MAX =KMAX(II)
            LOC =KLOC(II)
            IF(LIT.EQ.1) THEN
            ELSEIF(LIT.EQ.2) THEN
               IF(MIN.EQ.1) THEN
                  DO J=1,3
                     DO I=1,3
                        SP(J+LOC  ,I+LOC  )=PSHELL(I,J)
                     ENDDO
                  ENDDO
               ELSE
                  DO J=1,3
                     DO I=1,3
                        SP(J+LOC-1,I+LOC-1)=PSHELL(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ELSEIF(LIT.EQ.3) THEN
               DO J=1,6
                  DO I=1,6
                     SP(J+LOC-1,I+LOC-1)=DSHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.4) THEN
               DO J=1,10
                  DO I=1,10
                     SP(J+LOC-1,I+LOC-1)=FSHELL(I,J)
                  ENDDO
               ENDDO
            ELSEIF(LIT.EQ.5) THEN
               DO J=1,15
                  DO I=1,15
                     SP(J+LOC-1,I+LOC-1)=GSHELL(I,J)
                  ENDDO
               ENDDO
            ELSE
               WRITE(IW,9999)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ENDDO
C
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998) IFINV,IFTSP
         CALL HND_PRSQ(SP,NUM,NUM,NDIM)
      ENDIF
C
      RETURN
 9999 FORMAT(' -SPHMAT- LIMITED TO -G- SHELLS . STOP .')
 9998 FORMAT(/,1X,'SPHERICAL HARMONICS TRANSFORMATION MATRIX',
     1            ' -IFINV- = ',I5,' -IFTSP- = ',I5,
     2       /,1X,41(1H-))
      END
      SUBROUTINE HND_SPHPRT
      IMPLICIT REAL*8 (A-H,O-Z)
#include "global.fh"
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      CHARACTER*4 NAMXYZ,NAMSPH,NAMLAB
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_SPHERI/PSHELL(3,3),DSHELL(6,6),FSHELL(10,10),
     1              GSHELL(15,15),
     1              PIHELL(3,3),DIHELL(6,6),FIHELL(10,10),GIHELL(15,15)
      DIMENSION NF(5),NAMXYZ(35),NAMSPH(35),NAMLAB(35)
      DATA NF     /1,3,6,10,15/
      DATA NAMXYZ /'   S','   X','   Y','   Z','  XX','  YY','  ZZ',
     1             '  XY','  XZ','  YZ',' XXX',' YYY',' ZZZ',' XXY',
     2             ' XXZ',' YYX',' YYZ',' ZZX',' ZZY',' XYZ','XXXX',
     3             'YYYY','ZZZZ','XXXY','XXXZ','YYYX','YYYZ','ZZZX',
     4             'ZZZY','XXYY','XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/
      DATA NAMSPH /'  S ','P X ','P Y ','P Z ','D Z2','D X2','D XY',
     1             'D XZ','D YZ','D S ','  F1','  F2','  F3','  F4',
     2             '  F5','  F6','  F7','F P1','F P2','F P3','  G1',
     3             '  G2','  G3','  G4','  G5','  G6','  G7','  G8',
     4             '  G9','G S ','G D1','G D2','G D3','G D4','G D5'/
      DATA NAMLAB /'  S ','  P1','  P2','  P3','  D1','  D2','  D3',
     1             '  D4','  D5','  DS','  F1','  F2','  F3','  F4',
     2             '  F5','  F6','  F7','F P1','F P2','F P3','  G1',
     3             '  G2','  G3','  G4','  G5','  G6','  G7','  G8',
     4             '  G9','G S ','G D1','G D2','G D3','G D4','G D5'/
      DATA ONE    /1.0D+00/
C
C     ----  PRINT SPHERICAL HARMONICS DEFINITION -----
C
      if (ga_nodeid().eq.0) WRITE(IW,9999)
C
      MAXTYP=0
      DO ISH=1,NSHELL
         IF(KTYPE(ISH).GT.MAXTYP) MAXTYP=KTYPE(ISH)
      ENDDO
C
      DO 100 ISPH=1,MAXTYP
         N=NF(ISPH)
         GO TO (10,20,30,40,50),ISPH
C
C     ---- S ----
C
  10     J=1
            if (ga_nodeid().eq.0) WRITE(IW,9998) 
     1          J,NAMSPH(J   ),J,NAMLAB(J   ),
     1                     ONE,NAMXYZ(1)
         GO TO 100
C
C     ---- P ----
C
  20     DO J=1,N
            MAX=0
  21        MIN=MAX+1
            MAX=MAX+6
            IF(MAX.GT.N) MAX=N
            IF(MIN.LT.6) THEN
               if (ga_nodeid().eq.0) WRITE(IW,9998) 
     &             J,NAMSPH(J+ 1),J,NAMLAB(J+ 1),
     1                        (PSHELL(I,J),NAMXYZ(I+ 1),I=MIN,MAX)
            ELSE
               if (ga_nodeid().eq.0) WRITE(IW,9997)
     &            (PSHELL(I,J),NAMXYZ(I+ 1),I=MIN,MAX)
            ENDIF
            IF(MAX.LT.N) GO TO 21
         ENDDO
         GO TO 100
C
C     ---- D ----
C
  30     DO J=1,N
            MAX=0
  31        MIN=MAX+1
            MAX=MAX+6
            IF(MAX.GT.N) MAX=N
            IF(MIN.LT.6) THEN
               if (ga_nodeid().eq.0) WRITE(IW,9998)
     &            J,NAMSPH(J+4),J,NAMLAB(J+4),
     1                        (DSHELL(I,J),NAMXYZ(I+ 4),I=MIN,MAX)
            ELSE
               if (ga_nodeid().eq.0) WRITE(IW,9997)
     &            (DSHELL(I,J),NAMXYZ(I+ 4),I=MIN,MAX)
            ENDIF
            IF(MAX.LT.N) GO TO 31
         ENDDO
         GO TO 100
C
C     ---- F ----
C
  40     DO J=1,N
            MAX=0
  41        MIN=MAX+1
            MAX=MAX+6
            IF(MAX.GT.N) MAX=N
            IF(MIN.LT.6) THEN
               if (ga_nodeid().eq.0) WRITE(IW,9998)
     &            J,NAMSPH(J+10),J,NAMLAB(J+10),
     1                        (FSHELL(I,J),NAMXYZ(I+10),I=MIN,MAX)
            ELSE
               if (ga_nodeid().eq.0) WRITE(IW,9997) 
     &            (FSHELL(I,J),NAMXYZ(I+10),I=MIN,MAX)
            ENDIF
            IF(MAX.LT.N) GO TO 41
         ENDDO
         GO TO 100
C
C     ---- G ----
C
  50     DO J=1,N
            MAX=0
  51        MIN=MAX+1
            MAX=MAX+6
            IF(MAX.GT.N) MAX=N
            IF(MIN.LT.6) THEN
               if (ga_nodeid().eq.0) WRITE(IW,9998)
     &            J,NAMSPH(J+20),J,NAMLAB(J+20),
     1                        (GSHELL(I,J),NAMXYZ(I+20),I=MIN,MAX)
            ELSE
               if (ga_nodeid().eq.0) WRITE(IW,9997)
     &            (GSHELL(I,J),NAMXYZ(I+20),I=MIN,MAX)
            ENDIF
            IF(MAX.LT.N) GO TO 51
         ENDDO
         GO TO 100
C
  100 CONTINUE
C
      RETURN
 9999 FORMAT(/,1X,'Spherical harmonics',/,1X,19(1H-))
 9998 FORMAT(5H CHI(,I4,2X,A4,3H) =,1X,5H CHI(,I4,2X,A4,3H) =,
     1           1X,6(F7.3,3H * ,A4))
 9997 FORMAT(37X,1X,6(F7.3,3H * ,A4))
      END
C
C     ----- SUBROUTINES ABOVE TAKEN FROM SYM.F IN HONDO -----
C
      SUBROUTINE HND_PRSQL(V,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "global.fh"
#include "nwc_const.fh"
C
C     ----- PRINT OUT A SQUARE MATRIX WITH LABELS
C
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXBFN=3072)
      character*8 anam,bflab,cbf
      character*2 bnam
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM),BFLAB(MXBFN)
      DIMENSION V(NDIM,1)
      DIMENSION IC(5),CBF(5),C(5)
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
C
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      if (ga_nodeid().eq.0) then
         WRITE (IW,9008)
         IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
         IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
         IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
         WRITE (IW,9008)
         DO 120 J = 1,N
         IF(LIST.EQ.0) WRITE (IW,9048) J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
         IF(LIST.EQ.1) WRITE (IW,9148) J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
         IF(LIST.EQ.2) WRITE (IW,9248) J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
  120    CONTINUE
      endif
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      if (ga_nodeid().eq.0) WRITE(IW,9328) J
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      CBF(II)=BFLAB(I)
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      if (ga_nodeid().eq.0) WRITE(IW,9348) 
     &   (IC(IDUM),CBF(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      if (ga_nodeid().eq.0) WRITE(IW,9348) 
     &   (IC(IDUM),CBF(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(15X,10(3X,I4,4X))
 9048 FORMAT(I5,2X,A8,10F11.5)
 9128 FORMAT(15X,7(5X,I4,6X))
 9148 FORMAT(I5,2X,A8,7F15.10)
 9228 FORMAT(15X,7(5X,I4,6X))
 9248 FORMAT(I5,2X,A8,7E15.8)
 9328 FORMAT(/,' Vec =',I5,/,1X,10(1H-))
 9348 FORMAT(5(I5,1X,A8,F11.5))
      END
      SUBROUTINE HND_JACOBI(A,V,D,N,NP,NROT,B,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 ERRMSG
      COMMON/IOFILE/IR,IW
      DIMENSION A(NP,NP),V(NP,NP),D(NP),B(NP),Z(NP)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-JACOBI-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA PT2    /0.2D+00/
      DATA PT5    /0.5D+00/
      DATA HUNDRD /1.0D+02/
C
      DO 12 IP=1,N
         DO 11 IQ=1,N
            V(IQ,IP)=ZERO
   11    CONTINUE
         V(IP,IP)=ONE
   12 CONTINUE
      DO 13 IP=1,N
         B(IP)=A(IP,IP)
         D(IP)=B(IP)
         Z(IP)=ZERO
   13 CONTINUE
      NROT=0
      IF(N.EQ.1) RETURN
      DO 24 I=1,50
         SM=ZERO
         DO 15 IP=1,N-1
            DO 14 IQ=IP+1,N
               SM=SM+ABS(A(IP,IQ))
   14       CONTINUE
   15    CONTINUE
         IF(SM.EQ.ZERO) RETURN
         IF(I.LT.4) THEN
            THRES=PT2*SM/N**2
         ELSE
            THRES=ZERO
         ENDIF
         DO 22 IP=1,N-1
            DO 21 IQ=IP+1,N
               G=HUNDRD*ABS(A(IP,IQ))
               IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP))).AND.
     1                         (ABS(D(IQ))+G.EQ.ABS(D(IQ)))) THEN
                  A(IP,IQ)=ZERO
               ELSE IF(ABS(A(IP,IQ)).GT.THRES) THEN
                  H=D(IQ)-D(IP)
                  IF(ABS(H)+G.EQ.ABS(H)) THEN
                     T=A(IP,IQ)/H
                  ELSE
                     THETA=PT5*H/A(IP,IQ)
                     T=ONE/(ABS(THETA)+SQRT(ONE+THETA**2))
                     IF(THETA.LT.ZERO) T=-T
                  ENDIF
                  C=ONE/SQRT(ONE+T**2)
                  S=T*C
                  TAU=S/(ONE+C)
                  H=T*A(IP,IQ)
                  Z(IP)=Z(IP)-H
                  Z(IQ)=Z(IQ)+H
                  D(IP)=D(IP)-H
                  D(IQ)=D(IQ)+H
                  A(IP,IQ)=ZERO
                  DO 16 J=1,IP-1
                     G=A(J,IP)
                     H=A(J,IQ)
                     A(J,IP)=G-S*(H+G*TAU)
                     A(J,IQ)=H+S*(G-H*TAU)
   16             CONTINUE
                  DO 17 J=IP+1,IQ-1
                     G=A(IP,J)
                     H=A(J,IQ)
                     A(IP,J)=G-S*(H+G*TAU)
                     A(J,IQ)=H+S*(G-H*TAU)
   17             CONTINUE
                  DO 18 J=IQ+1,N
                     G=A(IP,J)
                     H=A(IQ,J)
                     A(IP,J)=G-S*(H+G*TAU)
                     A(IQ,J)=H+S*(G-H*TAU)
   18             CONTINUE
                  DO 19 J=1,N
                     G=V(J,IP)
                     H=V(J,IQ)
                     V(J,IP)=G-S*(H+G*TAU)
                     V(J,IQ)=H+S*(G-H*TAU)
   19             CONTINUE
                  NROT=NROT+1
               ENDIF
   21       CONTINUE
   22    CONTINUE
         DO 23 IP=1,N
            B(IP)=B(IP)+Z(IP)
            D(IP)=B(IP)
            Z(IP)=ZERO
   23    CONTINUE
   24 CONTINUE
      WRITE(IW,*) '50 ITERATIONS SHOULD NEVER HAPPEN'
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
      END
      SUBROUTINE HND_JACSRT(V,D,N,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(NP,NP),D(NP)
      IF(N.EQ.1) RETURN
      DO 13 I=1,N-1
         K=I
         P=D(I)
         DO 11 J=I+1,N
            IF(D(J).GE.P) THEN
               K=J
               P=D(J)
            ENDIF
   11    CONTINUE
         IF(K.NE.I) THEN
            D(K)=D(I)
            D(I)=P
            DO 12 J=1,N
               P=V(J,I)
               V(J,I)=V(J,K)
               V(J,K)=P
   12       CONTINUE
         ENDIF
   13 CONTINUE
      RETURN
      END
      SUBROUTINE HND_SPHSET
      IMPLICIT REAL*8 (A-H,O-Z)
#include "hnd_tol.fh"
      LOGICAL DBUG
C
C     -----  THIS ROUTINE DEFINES THE SPHERICAL HARMONICS -----
C     -----  AND THE INVERSE TRANSFORMATION               -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_SPHERI/PSHELL(3,3),DSHELL(6,6),FSHELL(10,10),
     1                  GSHELL(15,15),
     1              PIHELL(3,3),DIHELL(6,6),FIHELL(10,10),
     1                  GIHELL(15,15)
      COMMON/HND_FRAME/U1,U2,U3,V1,V2,V3,W1,W2,W3
      DIMENSION  PSPH(3,3), DSPH(6,6), FSPH(10,10), GSPH(15,15)
      DIMENSION IPSPH(3,3),IDSPH(6,6),IFSPH(10,10),IGSPH(15,15)
      DIMENSION BUFF(15,15),INDX(15),VV(15)
      DIMENSION ROT(15,15)
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA SQRT7 /2.64575131106459D+00/
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA IPSPH/
     1  1,  0,  0,
     2  0,  1,  0,
     3  0,  0,  1/
C
      DATA IDSPH/
     1 -1, -1,  2,  0,  0,  0,
     2  1, -1,  0,  0,  0,  0,
     3  0,  0,  0,  1,  0,  0,
     4  0,  0,  0,  0,  1,  0,
     5  0,  0,  0,  0,  0,  1,
     6  1,  1,  1,  0,  0,  0/
C
      DATA IFSPH/
     1   1,  0,  0,  0,  0, -3,  0,  0,  0,  0,
     2   0,  0,  0,  0,  1,  0, -1,  0,  0,  0,
     3  -1,  0,  0,  0,  0, -1,  0,  4,  0,  0,
     4   0,  0,  2,  0, -3,  0, -3,  0,  0,  0,
     5   0, -1,  0, -1,  0,  0,  0,  0,  4,  0,
     6   0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
     7   0, -1,  0,  3,  0,  0,  0,  0,  0,  0,
     8   1,  0,  0,  0,  0,  1,  0,  1,  0,  0,
     9   0,  1,  0,  1,  0,  0,  0,  0,  1,  0,
     1   0,  0,  1,  0,  1,  0,  1,  0,  0,  0/
C
      DATA IGSPH/
     1  1,  1,  0,  0,  0,  0,  0,  0,  0, -6,  0,  0,  0,  0,  0,
     2  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,
     3 -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,
     4  0,  0,  0,  0, -3,  0,  0,  4,  0,  0,  0,  0,  0, -3,  0,
     5  3,  3,  8,  0,  0,  0,  0,  0,  0,  6,-24,-24,  0,  0,  0,
     6  0,  0,  0,  0,  0,  0, -3,  0,  4,  0,  0,  0, -3,  0,  0,
     7  0,  0,  0, -1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  6,
     8  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -3,  0,  0,
     9  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     1  1,  1,  1,  0,  0,  0,  0,  0,  0,  2,  2,  2,  0,  0,  0,
     2 -1, -1,  2,  0,  0,  0,  0,  0,  0, -2,  1,  1,  0,  0,  0,
     3  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,
     4  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,
     5  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,
     6  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,  0,  0/
C
      DBUG=.FALSE.
      IF(DBUG) WRITE(IW,9999)
C
C     ----- INITIALIZATIONS -----
C
      DO 10 I=1,3
      DO  5 J=1,3
        PSPH(I,J)=ZERO
      PSHELL(I,J)=ZERO
  5   PIHELL(I,J)=ZERO
 10   PIHELL(I,I)=ONE
      DO 20 I=1,6
      DO 15 J=1,6
        DSPH(I,J)=ZERO
      DSHELL(I,J)=ZERO
 15   DIHELL(I,J)=ZERO
 20   DIHELL(I,I)=ONE
      DO 30 I=1,10
      DO 25 J=1,10
        FSPH(I,J)=ZERO
      FSHELL(I,J)=ZERO
 25   FIHELL(I,J)=ZERO
 30   FIHELL(I,I)=ONE
      DO 40 I=1,15
      DO 35 J=1,15
        GSPH(I,J)=ZERO
      GSHELL(I,J)=ZERO
 35   GIHELL(I,J)=ZERO
 40   GIHELL(I,I)=ONE
C
C     ----- P SPHERICAL HARMONICS -----
C
      DO 60 I=1,3
      ANORM=ZERO
      DO 50 J=1,3
      AA=DBLE(IPSPH(J,I))
      ANORM=ANORM+AA**2
  50  CONTINUE
      ANORM=ONE/ SQRT(ANORM)
      DO 60 J=1,3
      AA=DBLE(IPSPH(J,I))
      PSPH(J,I)=ANORM*AA
  60  CONTINUE
      IF(DBUG) CALL HND_SYMPRP(PSPH, 3, 3, 3)
C
C     ----- LOAD TRANSFORMATION OF FRAMES IN -ROT- -----
C
C     NEW X
      ROT(1,1)=U1
      ROT(2,1)=U2
      ROT(3,1)=U3
C     NEW Y
      ROT(1,2)=V1
      ROT(2,2)=V2
      ROT(3,2)=V3
C     NEW Z
      ROT(1,3)=W1
      ROT(2,3)=W2
      ROT(3,3)=W3
C
C     ----- EXPRESS SPHERICAL HARMONICS IN ABSOLUTE FRAME -----
C
      DO 80 I=1,3
      DO 80 J=1,3
      A=ZERO
      DO 70 K=1,3
      A=A+ROT(I,K)*PSPH(K,J)
  70  CONTINUE
        BUFF(I,J)=A
  80  PSHELL(I,J)=A
C
C     ----- DEFINE INVERSE TRANSFORMATION -----
C
      CALL HND_LUDCMP(BUFF,15, 3,INDX,VV,D,IERR)
      DO 90 I=1,3
      CALL HND_LUBKSB(BUFF,15, 3,INDX,PIHELL(1,I))
  90  CONTINUE
C
C     PRINT
C
      IF(DBUG) CALL HND_SYMPRP(PSHELL, 3, 3, 3)
      IF(DBUG) CALL HND_SYMPRP(PIHELL, 3, 3, 3)
C
C     ----- D SPHERICAL HARMONICS -----
C
      DO 120 I=1,6
      ANORM=ZERO
      DO 100 J=1,6
      AA=DBLE(IDSPH(J,I))
      ANORM=ANORM+AA**2
 100  CONTINUE
      ANORM=ONE/ SQRT(ANORM)
      DO 120 J=1,6
      AA=DBLE(IDSPH(J,I))
      DSPH(J,I)=ANORM*AA
 120  CONTINUE
      IF(DBUG) CALL HND_SYMPRP(DSPH, 6, 6, 6)
C
C     ----- LOAD TRANSFORMATION OF FRAMES IN -ROT- -----
C
C     NEW X2
      ROT(1,1)=U1*U1
      ROT(2,1)=U2*U2
      ROT(3,1)=U3*U3
      ROT(4,1)=U1*U2+U1*U2
      ROT(5,1)=U1*U3+U1*U3
      ROT(6,1)=U2*U3+U2*U3
C     NEW Y2
      ROT(1,2)=V1*V1
      ROT(2,2)=V2*V2
      ROT(3,2)=V3*V3
      ROT(4,2)=V1*V2+V1*V2
      ROT(5,2)=V1*V3+V1*V3
      ROT(6,2)=V2*V3+V2*V3
C     NEW Z2
      ROT(1,3)=W1*W1
      ROT(2,3)=W2*W2
      ROT(3,3)=W3*W3
      ROT(4,3)=W1*W2+W1*W2
      ROT(5,3)=W1*W3+W1*W3
      ROT(6,3)=W2*W3+W2*W3
C     NEW XY
      ROT(1,4)=U1*V1
      ROT(2,4)=U2*V2
      ROT(3,4)=U3*V3
      ROT(4,4)=U1*V2+U2*V1
      ROT(5,4)=U1*V3+U3*V1
      ROT(6,4)=U2*V3+U3*V2
C     NEW XZ
      ROT(1,5)=U1*W1
      ROT(2,5)=U2*W2
      ROT(3,5)=U3*W3
      ROT(4,5)=U1*W2+U2*W1
      ROT(5,5)=U1*W3+U3*W1
      ROT(6,5)=U2*W3+U3*W2
C     NEW YZ
      ROT(1,6)=V1*W1
      ROT(2,6)=V2*W2
      ROT(3,6)=V3*W3
      ROT(4,6)=V1*W2+V2*W1
      ROT(5,6)=V1*W3+V3*W1
      ROT(6,6)=V2*W3+V3*W2
C
C     ----- EXPRESS SPHERICAL HARMONICS IN ABSOLUTE FRAME -----
C
      DO 160 I=1,6
      DO 160 J=1,6
      A=ZERO
      DO 150 K=1,6
      A=A+ROT(I,K)*DSPH(K,J)
 150  CONTINUE
        BUFF(I,J)=A
 160  DSHELL(I,J)=A
C
C     ----- ACCOUNT FOR NORMALIZATION OF CARTESIAN FUNCTIONS -----
C
      IF(NORMF.EQ.0.OR.NORMP.EQ.0) THEN
      DO 170 I=1,6
      DSHELL(4,I)=DSHELL(4,I)/SQRT3
      DSHELL(5,I)=DSHELL(5,I)/SQRT3
      DSHELL(6,I)=DSHELL(6,I)/SQRT3
        BUFF(4,I)=DSHELL(4,I)
        BUFF(5,I)=DSHELL(5,I)
        BUFF(6,I)=DSHELL(6,I)
  170 CONTINUE
      END IF
C
C     ----- DEFINE INVERSE TRANSFORMATION -----
C
      CALL HND_LUDCMP(BUFF,15, 6,INDX,VV,D,IERR)
      DO 180 I=1,6
      CALL HND_LUBKSB(BUFF,15, 6,INDX,DIHELL(1,I))
 180  CONTINUE
C
C     PRINT
C
      IF(DBUG) CALL HND_SYMPRP(DSHELL,6,6,6)
      IF(DBUG) CALL HND_SYMPRP(DIHELL,6,6,6)
C
C     ----- F SPHERICAL HARMONICS -----
C
      DO 220 I=1,10
      ANORM=ZERO
      DO 200 J=1,10
      AA=DBLE(IFSPH(J,I))
      ANORM=ANORM+AA**2
 200  CONTINUE
      ANORM=ONE/ SQRT(ANORM)
      DO 220 J=1,10
      AA=DBLE(IFSPH(J,I))
      FSPH(J,I)=ANORM*AA
 220  CONTINUE
      IF(DBUG) CALL HND_SYMPRP(FSPH,10,10,10)
C
C     ----- LOAD TRANSFORMATION OF FRAMES IN -ROT- -----
C
C     NEW XXX
      ROT(1,1) =U1*U1*U1
      ROT(2,1) =U2*U2*U2
      ROT(3,1) =U3*U3*U3
      ROT(4,1) =U2*U1*U1+ U1*U2*U1+ U1*U1*U2
      ROT(5,1) =U3*U1*U1+ U1*U3*U1+ U1*U1*U3
      ROT(6,1) =U1*U2*U2+ U2*U1*U2+ U2*U2*U1
      ROT(7,1) =U3*U2*U2+ U2*U3*U2+ U2*U2*U3
      ROT(8,1) =U1*U3*U3+ U3*U1*U3+ U3*U3*U1
      ROT(9,1) =U2*U3*U3+ U3*U2*U3+ U3*U3*U2
      ROT(10,1)=U1*U2*U3+ U1*U3*U2+ U2*U1*U3+ U2*U3*U1+
     1          U3*U1*U2+ U3*U2*U1
C     NEW YYY
      ROT(1,2) =V1*V1*V1
      ROT(2,2) =V2*V2*V2
      ROT(3,2) =V3*V3*V3
      ROT(4,2) =V2*V1*V1+ V1*V2*V1+ V1*V1*V2
      ROT(5,2) =V3*V1*V1+ V1*V3*V1+ V1*V1*V3
      ROT(6,2) =V1*V2*V2+ V2*V1*V2+ V2*V2*V1
      ROT(7,2) =V3*V2*V2+ V2*V3*V2+ V2*V2*V3
      ROT(8,2) =V1*V3*V3+ V3*V1*V3+ V3*V3*V1
      ROT(9,2) =V2*V3*V3+ V3*V2*V3+ V3*V3*V2
      ROT(10,2)=V1*V2*V3+ V1*V3*V2+ V2*V1*V3+ V2*V3*V1+
     1          V3*V1*V2+ V3*V2*V1
C     NEW ZZZ
      ROT(1,3) =W1*W1*W1
      ROT(2,3) =W2*W2*W2
      ROT(3,3) =W3*W3*W3
      ROT(4,3) =W2*W1*W1+ W1*W2*W1+ W1*W1*W2
      ROT(5,3) =W3*W1*W1+ W1*W3*W1+ W1*W1*W3
      ROT(6,3) =W1*W2*W2+ W2*W1*W2+ W2*W2*W1
      ROT(7,3) =W3*W2*W2+ W2*W3*W2+ W2*W2*W3
      ROT(8,3) =W1*W3*W3+ W3*W1*W3+ W3*W3*W1
      ROT(9,3) =W2*W3*W3+ W3*W2*W3+ W3*W3*W2
      ROT(10,3)=W1*W2*W3+ W1*W3*W2+ W2*W1*W3+ W2*W3*W1+
     1          W3*W1*W2+ W3*W2*W1
C     NEW XXY
      ROT(1,4) =U1*U1*V1
      ROT(2,4) =U2*U2*V2
      ROT(3,4) =U3*U3*V3
      ROT(4,4) =U2*U1*V1+ U1*U2*V1+ U1*U1*V2
      ROT(5,4) =U3*U1*V1+ U1*U3*V1+ U1*U1*V3
      ROT(6,4) =U1*U2*V2+ U2*U1*V2+ U2*U2*V1
      ROT(7,4) =U3*U2*V2+ U2*U3*V2+ U2*U2*V3
      ROT(8,4) =U1*U3*V3+ U3*U1*V3+ U3*U3*V1
      ROT(9,4) =U2*U3*V3+ U3*U2*V3+ U3*U3*V2
      ROT(10,4)=U1*U2*V3+ U1*U3*V2+ U2*U1*V3+ U2*U3*V1+
     1          U3*U1*V2+ U3*U2*V1
C     NEW XXZ
      ROT(1,5) =U1*U1*W1
      ROT(2,5) =U2*U2*W2
      ROT(3,5) =U3*U3*W3
      ROT(4,5) =U2*U1*W1+ U1*U2*W1+ U1*U1*W2
      ROT(5,5) =U3*U1*W1+ U1*U3*W1+ U1*U1*W3
      ROT(6,5) =U1*U2*W2+ U2*U1*W2+ U2*U2*W1
      ROT(7,5) =U3*U2*W2+ U2*U3*W2+ U2*U2*W3
      ROT(8,5) =U1*U3*W3+ U3*U1*W3+ U3*U3*W1
      ROT(9,5) =U2*U3*W3+ U3*U2*W3+ U3*U3*W2
      ROT(10,5)=U1*U2*W3+ U1*U3*W2+ U2*U1*W3+ U2*U3*W1+
     1          U3*U1*W2+ U3*U2*W1
C     NEW YYX
      ROT(1,6) =V1*V1*U1
      ROT(2,6) =V2*V2*U2
      ROT(3,6) =V3*V3*U3
      ROT(4,6) =V2*V1*U1+ V1*V2*U1+ V1*V1*U2
      ROT(5,6) =V3*V1*U1+ V1*V3*U1+ V1*V1*U3
      ROT(6,6) =V1*V2*U2+ V2*V1*U2+ V2*V2*U1
      ROT(7,6) =V3*V2*U2+ V2*V3*U2+ V2*V2*U3
      ROT(8,6) =V1*V3*U3+ V3*V1*U3+ V3*V3*U1
      ROT(9,6) =V2*V3*U3+ V3*V2*U3+ V3*V3*U2
      ROT(10,6)=V1*V2*U3+ V1*V3*U2+ V2*V1*U3+ V2*V3*U1+
     1          V3*V1*U2+ V3*V2*U1
C     NEW YYZ
      ROT(1,7) =V1*V1*W1
      ROT(2,7) =V2*V2*W2
      ROT(3,7) =V3*V3*W3
      ROT(4,7) =V2*V1*W1+ V1*V2*W1+ V1*V1*W2
      ROT(5,7) =V3*V1*W1+ V1*V3*W1+ V1*V1*W3
      ROT(6,7) =V1*V2*W2+ V2*V1*W2+ V2*V2*W1
      ROT(7,7) =V3*V2*W2+ V2*V3*W2+ V2*V2*W3
      ROT(8,7) =V1*V3*W3+ V3*V1*W3+ V3*V3*W1
      ROT(9,7) =V2*V3*W3+ V3*V2*W3+ V3*V3*W2
      ROT(10,7)=V1*V2*W3+ V1*V3*W2+ V2*V1*W3+ V2*V3*W1+
     1          V3*V1*W2+ V3*V2*W1
C     NEW ZZX
      ROT(1,8) =W1*W1*U1
      ROT(2,8) =W2*W2*U2
      ROT(3,8) =W3*W3*U3
      ROT(4,8) =W2*W1*U1+ W1*W2*U1+ W1*W1*U2
      ROT(5,8) =W3*W1*U1+ W1*W3*U1+ W1*W1*U3
      ROT(6,8) =W1*W2*U2+ W2*W1*U2+ W2*W2*U1
      ROT(7,8) =W3*W2*U2+ W2*W3*U2+ W2*W2*U3
      ROT(8,8) =W1*W3*U3+ W3*W1*U3+ W3*W3*U1
      ROT(9,8) =W2*W3*U3+ W3*W2*U3+ W3*W3*U2
      ROT(10,8)=W1*W2*U3+ W1*W3*U2+ W2*W1*U3+ W2*W3*U1+
     1          W3*W1*U2+ W3*W2*U1
C     NEW ZZY
      ROT(1,9) =W1*W1*V1
      ROT(2,9) =W2*W2*V2
      ROT(3,9) =W3*W3*V3
      ROT(4,9) =W2*W1*V1+ W1*W2*V1+ W1*W1*V2
      ROT(5,9) =W3*W1*V1+ W1*W3*V1+ W1*W1*V3
      ROT(6,9) =W1*W2*V2+ W2*W1*V2+ W2*W2*V1
      ROT(7,9) =W3*W2*V2+ W2*W3*V2+ W2*W2*V3
      ROT(8,9) =W1*W3*V3+ W3*W1*V3+ W3*W3*V1
      ROT(9,9) =W2*W3*V3+ W3*W2*V3+ W3*W3*V2
      ROT(10,9)=W1*W2*V3+ W1*W3*V2+ W2*W1*V3+ W2*W3*V1+
     1          W3*W1*V2+ W3*W2*V1
C     NEW XYZ
      ROT(1,10) =U1*V1*W1
      ROT(2,10) =U2*V2*W2
      ROT(3,10) =U3*V3*W3
      ROT(4,10) =U2*V1*W1+ U1*V2*W1+ U1*V1*W2
      ROT(5,10) =U3*V1*W1+ U1*V3*W1+ U1*V1*W3
      ROT(6,10) =U1*V2*W2+ U2*V1*W2+ U2*V2*W1
      ROT(7,10) =U3*V2*W2+ U2*V3*W2+ U2*V2*W3
      ROT(8,10) =U1*V3*W3+ U3*V1*W3+ U3*V3*W1
      ROT(9,10) =U2*V3*W3+ U3*V2*W3+ U3*V3*W2
      ROT(10,10)=U1*V2*W3+ U1*V3*W2+ U2*V1*W3+ U2*V3*W1+
     1           U3*V1*W2+ U3*V2*W1
C
C     ----- EXPRESS SPHERICAL HARMONICS IN ABSOLUTE FRAME -----
C
      DO 260 I=1,10
      DO 260 J=1,10
      A=ZERO
      DO 250 K=1,10
      A=A+ROT(I,K)*FSPH(K,J)
 250  CONTINUE
        BUFF(I,J)=A
 260  FSHELL(I,J)=A
C
C     ----- ACCOUNT FOR NORMALIZATION OF CARTESIAN FUNCTIONS -----
C
      IF(NORMF.EQ.0.OR.NORMP.EQ.0) THEN
      DO 270 I=1,10
      FSHELL( 4,I)=FSHELL( 4,I)/SQRT5
      FSHELL( 5,I)=FSHELL( 5,I)/SQRT5
      FSHELL( 6,I)=FSHELL( 6,I)/SQRT5
      FSHELL( 7,I)=FSHELL( 7,I)/SQRT5
      FSHELL( 8,I)=FSHELL( 8,I)/SQRT5
      FSHELL( 9,I)=FSHELL( 9,I)/SQRT5
      FSHELL(10,I)=FSHELL(10,I)/(SQRT5*SQRT3)
        BUFF( 4,I)=FSHELL( 4,I)
        BUFF( 5,I)=FSHELL( 5,I)
        BUFF( 6,I)=FSHELL( 6,I)
        BUFF( 7,I)=FSHELL( 7,I)
        BUFF( 8,I)=FSHELL( 8,I)
        BUFF( 9,I)=FSHELL( 9,I)
        BUFF(10,I)=FSHELL(10,I)
  270 CONTINUE
      ENDIF
C
C     ----- DEFINE INVERSE TRANSFORMATION -----
C
      CALL HND_LUDCMP(BUFF,15,10,INDX,VV,D,IERR)
      DO 280 I=1,10
      CALL HND_LUBKSB(BUFF,15,10,INDX,FIHELL(1,I))
 280  CONTINUE
C
C     PRINT
C
      IF(DBUG) CALL HND_SYMPRP(FSHELL,10,10,10)
      IF(DBUG) CALL HND_SYMPRP(FIHELL,10,10,10)
C
C     ----- G SPHERICAL HARMONICS -----
C
      DO 320 I=1,15
      ANORM=ZERO
      DO 300 J=1,15
      AA=DBLE(IGSPH(J,I))
      ANORM=ANORM+AA**2
 300  CONTINUE
      ANORM=ONE/ SQRT(ANORM)
      DO 320 J=1,15
      AA=DBLE(IGSPH(J,I))
      GSPH(J,I)=ANORM*AA
 320  CONTINUE
      IF(DBUG) CALL HND_SYMPRP(GSPH,15,15,15)
C
C     ----- LOAD TRANSFORMATION OF FRAMES IN -ROT- -----
C
C     NEW XXXX
      ROT(1,1)  =U1*U1*U1*U1
      ROT(2,1)  =U2*U2*U2*U2
      ROT(3,1)  =U3*U3*U3*U3
      ROT(4,1)  =U2*U1*U1*U1+ U1*U2*U1*U1+ U1*U1*U2*U1+ U1*U1*U1*U2
      ROT(5,1)  =U3*U1*U1*U1+ U1*U3*U1*U1+ U1*U1*U3*U1+ U1*U1*U1*U3
      ROT(6,1)  =U1*U2*U2*U2+ U2*U1*U2*U2+ U2*U2*U1*U2+ U2*U2*U2*U1
      ROT(7,1)  =U3*U2*U2*U2+ U2*U3*U2*U2+ U2*U2*U3*U2+ U2*U2*U2*U3
      ROT(8,1)  =U1*U3*U3*U3+ U3*U1*U3*U3+ U3*U3*U1*U3+ U3*U3*U3*U1
      ROT(9,1)  =U2*U3*U3*U3+ U3*U2*U3*U3+ U3*U3*U2*U3+ U3*U3*U3*U2
      ROT(10,1) =U1*U1*U2*U2+ U2*U2*U1*U1+ U1*U2*U1*U2+ U1*U2*U2*U1
     1         + U2*U1*U1*U2+ U2*U1*U2*U1
      ROT(11,1) =U1*U1*U3*U3+ U3*U3*U1*U1+ U1*U3*U1*U3+ U1*U3*U3*U1
     1         + U3*U1*U1*U3+ U3*U1*U3*U1
      ROT(12,1) =U2*U2*U3*U3+ U3*U3*U2*U2+ U2*U3*U2*U3+ U2*U3*U3*U2
     1         + U3*U2*U2*U3+ U3*U2*U3*U2
      ROT(13,1) =U1*U1*U2*U3+ U1*U1*U3*U2+ U1*U2*U1*U3+ U1*U3*U1*U2
     1         + U1*U2*U3*U1+ U1*U3*U2*U1+ U2*U1*U1*U3+ U3*U1*U1*U2
     2         + U2*U1*U3*U1+ U3*U1*U2*U1+ U2*U3*U1*U1+ U3*U2*U1*U1
      ROT(14,1) =U2*U2*U1*U3+ U2*U2*U3*U1+ U2*U1*U2*U3+ U2*U3*U2*U1
     1         + U2*U1*U3*U2+ U2*U3*U1*U2+ U1*U2*U2*U3+ U3*U2*U2*U1
     2         + U1*U2*U3*U2+ U3*U2*U1*U2+ U1*U3*U2*U2+ U3*U1*U2*U2
      ROT(15,1) =U3*U3*U1*U2+ U3*U3*U2*U1+ U3*U1*U3*U2+ U3*U2*U3*U1
     1         + U3*U1*U2*U3+ U3*U2*U1*U3+ U1*U3*U3*U2+ U2*U3*U3*U1
     2         + U1*U3*U2*U3+ U2*U3*U1*U3+ U1*U2*U3*U3+ U2*U1*U3*U3
C     NEW YYYY
      ROT(1,2)  =V1*V1*V1*V1
      ROT(2,2)  =V2*V2*V2*V2
      ROT(3,2)  =V3*V3*V3*V3
      ROT(4,2)  =V2*V1*V1*V1+ V1*V2*V1*V1+ V1*V1*V2*V1+ V1*V1*V1*V2
      ROT(5,2)  =V3*V1*V1*V1+ V1*V3*V1*V1+ V1*V1*V3*V1+ V1*V1*V1*V3
      ROT(6,2)  =V1*V2*V2*V2+ V2*V1*V2*V2+ V2*V2*V1*V2+ V2*V2*V2*V1
      ROT(7,2)  =V3*V2*V2*V2+ V2*V3*V2*V2+ V2*V2*V3*V2+ V2*V2*V2*V3
      ROT(8,2)  =V1*V3*V3*V3+ V3*V1*V3*V3+ V3*V3*V1*V3+ V3*V3*V3*V1
      ROT(9,2)  =V2*V3*V3*V3+ V3*V2*V3*V3+ V3*V3*V2*V3+ V3*V3*V3*V2
      ROT(10,2) =V1*V1*V2*V2+ V2*V2*V1*V1+ V1*V2*V1*V2+ V1*V2*V2*V1
     1         + V2*V1*V1*V2+ V2*V1*V2*V1
      ROT(11,2) =V1*V1*V3*V3+ V3*V3*V1*V1+ V1*V3*V1*V3+ V1*V3*V3*V1
     1         + V3*V1*V1*V3+ V3*V1*V3*V1
      ROT(12,2) =V2*V2*V3*V3+ V3*V3*V2*V2+ V2*V3*V2*V3+ V2*V3*V3*V2
     1         + V3*V2*V2*V3+ V3*V2*V3*V2
      ROT(13,2) =V1*V1*V2*V3+ V1*V1*V3*V2+ V1*V2*V1*V3+ V1*V3*V1*V2
     1         + V1*V2*V3*V1+ V1*V3*V2*V1+ V2*V1*V1*V3+ V3*V1*V1*V2
     2         + V2*V1*V3*V1+ V3*V1*V2*V1+ V2*V3*V1*V1+ V3*V2*V1*V1
      ROT(14,2) =V2*V2*V1*V3+ V2*V2*V3*V1+ V2*V1*V2*V3+ V2*V3*V2*V1
     1         + V2*V1*V3*V2+ V2*V3*V1*V2+ V1*V2*V2*V3+ V3*V2*V2*V1
     2         + V1*V2*V3*V2+ V3*V2*V1*V2+ V1*V3*V2*V2+ V3*V1*V2*V2
      ROT(15,2) =V3*V3*V1*V2+ V3*V3*V2*V1+ V3*V1*V3*V2+ V3*V2*V3*V1
     1         + V3*V1*V2*V3+ V3*V2*V1*V3+ V1*V3*V3*V2+ V2*V3*V3*V1
     2         + V1*V3*V2*V3+ V2*V3*V1*V3+ V1*V2*V3*V3+ V2*V1*V3*V3
C     NEW ZZZZ
      ROT(1,3)  =W1*W1*W1*W1
      ROT(2,3)  =W2*W2*W2*W2
      ROT(3,3)  =W3*W3*W3*W3
      ROT(4,3)  =W2*W1*W1*W1+ W1*W2*W1*W1+ W1*W1*W2*W1+ W1*W1*W1*W2
      ROT(5,3)  =W3*W1*W1*W1+ W1*W3*W1*W1+ W1*W1*W3*W1+ W1*W1*W1*W3
      ROT(6,3)  =W1*W2*W2*W2+ W2*W1*W2*W2+ W2*W2*W1*W2+ W2*W2*W2*W1
      ROT(7,3)  =W3*W2*W2*W2+ W2*W3*W2*W2+ W2*W2*W3*W2+ W2*W2*W2*W3
      ROT(8,3)  =W1*W3*W3*W3+ W3*W1*W3*W3+ W3*W3*W1*W3+ W3*W3*W3*W1
      ROT(9,3)  =W2*W3*W3*W3+ W3*W2*W3*W3+ W3*W3*W2*W3+ W3*W3*W3*W2
      ROT(10,3) =W1*W1*W2*W2+ W2*W2*W1*W1+ W1*W2*W1*W2+ W1*W2*W2*W1
     1         + W2*W1*W1*W2+ W2*W1*W2*W1
      ROT(11,3) =W1*W1*W3*W3+ W3*W3*W1*W1+ W1*W3*W1*W3+ W1*W3*W3*W1
     1         + W3*W1*W1*W3+ W3*W1*W3*W1
      ROT(12,3) =W2*W2*W3*W3+ W3*W3*W2*W2+ W2*W3*W2*W3+ W2*W3*W3*W2
     1         + W3*W2*W2*W3+ W3*W2*W3*W2
      ROT(13,3) =W1*W1*W2*W3+ W1*W1*W3*W2+ W1*W2*W1*W3+ W1*W3*W1*W2
     1         + W1*W2*W3*W1+ W1*W3*W2*W1+ W2*W1*W1*W3+ W3*W1*W1*W2
     2         + W2*W1*W3*W1+ W3*W1*W2*W1+ W2*W3*W1*W1+ W3*W2*W1*W1
      ROT(14,3) =W2*W2*W1*W3+ W2*W2*W3*W1+ W2*W1*W2*W3+ W2*W3*W2*W1
     1         + W2*W1*W3*W2+ W2*W3*W1*W2+ W1*W2*W2*W3+ W3*W2*W2*W1
     2         + W1*W2*W3*W2+ W3*W2*W1*W2+ W1*W3*W2*W2+ W3*W1*W2*W2
      ROT(15,3) =W3*W3*W1*W2+ W3*W3*W2*W1+ W3*W1*W3*W2+ W3*W2*W3*W1
     1         + W3*W1*W2*W3+ W3*W2*W1*W3+ W1*W3*W3*W2+ W2*W3*W3*W1
     2         + W1*W3*W2*W3+ W2*W3*W1*W3+ W1*W2*W3*W3+ W2*W1*W3*W3
C     NEW XXXY
      ROT(1,4)  =U1*U1*U1*V1
      ROT(2,4)  =U2*U2*U2*V2
      ROT(3,4)  =U3*U3*U3*V3
      ROT(4,4)  =U2*U1*U1*V1+ U1*U2*U1*V1+ U1*U1*U2*V1+ U1*U1*U1*V2
      ROT(5,4)  =U3*U1*U1*V1+ U1*U3*U1*V1+ U1*U1*U3*V1+ U1*U1*U1*V3
      ROT(6,4)  =U1*U2*U2*V2+ U2*U1*U2*V2+ U2*U2*U1*V2+ U2*U2*U2*V1
      ROT(7,4)  =U3*U2*U2*V2+ U2*U3*U2*V2+ U2*U2*U3*V2+ U2*U2*U2*V3
      ROT(8,4)  =U1*U3*U3*V3+ U3*U1*U3*V3+ U3*U3*U1*V3+ U3*U3*U3*V1
      ROT(9,4)  =U2*U3*U3*V3+ U3*U2*U3*V3+ U3*U3*U2*V3+ U3*U3*U3*V2
      ROT(10,4) =U1*U1*U2*V2+ U2*U2*U1*V1+ U1*U2*U1*V2+ U1*U2*U2*V1
     1         + U2*U1*U1*V2+ U2*U1*U2*V1
      ROT(11,4) =U1*U1*U3*V3+ U3*U3*U1*V1+ U1*U3*U1*V3+ U1*U3*U3*V1
     1         + U3*U1*U1*V3+ U3*U1*U3*V1
      ROT(12,4) =U2*U2*U3*V3+ U3*U3*U2*V2+ U2*U3*U2*V3+ U2*U3*U3*V2
     1         + U3*U2*U2*V3+ U3*U2*U3*V2
      ROT(13,4) =U1*U1*U2*V3+ U1*U1*U3*V2+ U1*U2*U1*V3+ U1*U3*U1*V2
     1         + U1*U2*U3*V1+ U1*U3*U2*V1+ U2*U1*U1*V3+ U3*U1*U1*V2
     2         + U2*U1*U3*V1+ U3*U1*U2*V1+ U2*U3*U1*V1+ U3*U2*U1*V1
      ROT(14,4) =U2*U2*U1*V3+ U2*U2*U3*V1+ U2*U1*U2*V3+ U2*U3*U2*V1
     1         + U2*U1*U3*V2+ U2*U3*U1*V2+ U1*U2*U2*V3+ U3*U2*U2*V1
     2         + U1*U2*U3*V2+ U3*U2*U1*V2+ U1*U3*U2*V2+ U3*U1*U2*V2
      ROT(15,4) =U3*U3*U1*V2+ U3*U3*U2*V1+ U3*U1*U3*V2+ U3*U2*U3*V1
     1         + U3*U1*U2*V3+ U3*U2*U1*V3+ U1*U3*U3*V2+ U2*U3*U3*V1
     2         + U1*U3*U2*V3+ U2*U3*U1*V3+ U1*U2*U3*V3+ U2*U1*U3*V3
C     NEW XXXZ
      ROT(1,5)  =U1*U1*U1*W1
      ROT(2,5)  =U2*U2*U2*W2
      ROT(3,5)  =U3*U3*U3*W3
      ROT(4,5)  =U2*U1*U1*W1+ U1*U2*U1*W1+ U1*U1*U2*W1+ U1*U1*U1*W2
      ROT(5,5)  =U3*U1*U1*W1+ U1*U3*U1*W1+ U1*U1*U3*W1+ U1*U1*U1*W3
      ROT(6,5)  =U1*U2*U2*W2+ U2*U1*U2*W2+ U2*U2*U1*W2+ U2*U2*U2*W1
      ROT(7,5)  =U3*U2*U2*W2+ U2*U3*U2*W2+ U2*U2*U3*W2+ U2*U2*U2*W3
      ROT(8,5)  =U1*U3*U3*W3+ U3*U1*U3*W3+ U3*U3*U1*W3+ U3*U3*U3*W1
      ROT(9,5)  =U2*U3*U3*W3+ U3*U2*U3*W3+ U3*U3*U2*W3+ U3*U3*U3*W2
      ROT(10,5) =U1*U1*U2*W2+ U2*U2*U1*W1+ U1*U2*U1*W2+ U1*U2*U2*W1
     1         + U2*U1*U1*W2+ U2*U1*U2*W1
      ROT(11,5) =U1*U1*U3*W3+ U3*U3*U1*W1+ U1*U3*U1*W3+ U1*U3*U3*W1
     1         + U3*U1*U1*W3+ U3*U1*U3*W1
      ROT(12,5) =U2*U2*U3*W3+ U3*U3*U2*W2+ U2*U3*U2*W3+ U2*U3*U3*W2
     1         + U3*U2*U2*W3+ U3*U2*U3*W2
      ROT(13,5) =U1*U1*U2*W3+ U1*U1*U3*W2+ U1*U2*U1*W3+ U1*U3*U1*W2
     1         + U1*U2*U3*W1+ U1*U3*U2*W1+ U2*U1*U1*W3+ U3*U1*U1*W2
     2         + U2*U1*U3*W1+ U3*U1*U2*W1+ U2*U3*U1*W1+ U3*U2*U1*W1
      ROT(14,5) =U2*U2*U1*W3+ U2*U2*U3*W1+ U2*U1*U2*W3+ U2*U3*U2*W1
     1         + U2*U1*U3*W2+ U2*U3*U1*W2+ U1*U2*U2*W3+ U3*U2*U2*W1
     2         + U1*U2*U3*W2+ U3*U2*U1*W2+ U1*U3*U2*W2+ U3*U1*U2*W2
      ROT(15,5) =U3*U3*U1*W2+ U3*U3*U2*W1+ U3*U1*U3*W2+ U3*U2*U3*W1
     1         + U3*U1*U2*W3+ U3*U2*U1*W3+ U1*U3*U3*W2+ U2*U3*U3*W1
     2         + U1*U3*U2*W3+ U2*U3*U1*W3+ U1*U2*U3*W3+ U2*U1*U3*W3
C     NEW YYYX
      ROT(1,6)  =V1*V1*V1*U1
      ROT(2,6)  =V2*V2*V2*U2
      ROT(3,6)  =V3*V3*V3*U3
      ROT(4,6)  =V2*V1*V1*U1+ V1*V2*V1*U1+ V1*V1*V2*U1+ V1*V1*V1*U2
      ROT(5,6)  =V3*V1*V1*U1+ V1*V3*V1*U1+ V1*V1*V3*U1+ V1*V1*V1*U3
      ROT(6,6)  =V1*V2*V2*U2+ V2*V1*V2*U2+ V2*V2*V1*U2+ V2*V2*V2*U1
      ROT(7,6)  =V3*V2*V2*U2+ V2*V3*V2*U2+ V2*V2*V3*U2+ V2*V2*V2*U3
      ROT(8,6)  =V1*V3*V3*U3+ V3*V1*V3*U3+ V3*V3*V1*U3+ V3*V3*V3*U1
      ROT(9,6)  =V2*V3*V3*U3+ V3*V2*V3*U3+ V3*V3*V2*U3+ V3*V3*V3*U2
      ROT(10,6) =V1*V1*V2*U2+ V2*V2*V1*U1+ V1*V2*V1*U2+ V1*V2*V2*U1
     1         + V2*V1*V1*U2+ V2*V1*V2*U1
      ROT(11,6) =V1*V1*V3*U3+ V3*V3*V1*U1+ V1*V3*V1*U3+ V1*V3*V3*U1
     1         + V3*V1*V1*U3+ V3*V1*V3*U1
      ROT(12,6) =V2*V2*V3*U3+ V3*V3*V2*U2+ V2*V3*V2*U3+ V2*V3*V3*U2
     1         + V3*V2*V2*U3+ V3*V2*V3*U2
      ROT(13,6) =V1*V1*V2*U3+ V1*V1*V3*U2+ V1*V2*V1*U3+ V1*V3*V1*U2
     1         + V1*V2*V3*U1+ V1*V3*V2*U1+ V2*V1*V1*U3+ V3*V1*V1*U2
     2         + V2*V1*V3*U1+ V3*V1*V2*U1+ V2*V3*V1*U1+ V3*V2*V1*U1
      ROT(14,6) =V2*V2*V1*U3+ V2*V2*V3*U1+ V2*V1*V2*U3+ V2*V3*V2*U1
     1         + V2*V1*V3*U2+ V2*V3*V1*U1+ V1*V2*V2*U3+ V3*V2*V2*U1
     2         + V1*V2*V3*U2+ V3*V2*V1*U2+ V1*V3*V2*U2+ V3*V1*V2*U2
      ROT(15,6) =V3*V3*V1*U2+ V3*V3*V2*U1+ V3*V1*V3*U2+ V3*V2*V3*U1
     1         + V3*V1*V2*U3+ V3*V2*V1*U3+ V1*V3*V3*U2+ V2*V3*V3*U1
     2         + V1*V3*V2*U3+ V2*V3*V1*U3+ V1*V2*V3*U3+ V2*V1*V3*U3
C     NEW YYYZ
      ROT(1,7)  =V1*V1*V1*W1
      ROT(2,7)  =V2*V2*V2*W2
      ROT(3,7)  =V3*V3*V3*W3
      ROT(4,7)  =V2*V1*V1*W1+ V1*V2*V1*W1+ V1*V1*V2*W1+ V1*V1*V1*W2
      ROT(5,7)  =V3*V1*V1*W1+ V1*V3*V1*W1+ V1*V1*V3*W1+ V1*V1*V1*W3
      ROT(6,7)  =V1*V2*V2*W2+ V2*V1*V2*W2+ V2*V2*V1*W2+ V2*V2*V2*W1
      ROT(7,7)  =V3*V2*V2*W2+ V2*V3*V2*W2+ V2*V2*V3*W2+ V2*V2*V2*W3
      ROT(8,7)  =V1*V3*V3*W3+ V3*V1*V3*W3+ V3*V3*V1*W3+ V3*V3*V3*W1
      ROT(9,7)  =V2*V3*V3*W3+ V3*V2*V3*W3+ V3*V3*V2*W3+ V3*V3*V3*W2
      ROT(10,7) =V1*V1*V2*W2+ V2*V2*V1*W1+ V1*V2*V1*W2+ V1*V2*V2*W1
     1         + V2*V1*V1*W2+ V2*V1*V2*W1
      ROT(11,7) =V1*V1*V3*W3+ V3*V3*V1*W1+ V1*V3*V1*W3+ V1*V3*V3*W1
     1         + V3*V1*V1*W3+ V3*V1*V3*W1
      ROT(12,7) =V2*V2*V3*W3+ V3*V3*V2*W2+ V2*V3*V2*W3+ V2*V3*V3*W2
     1         + V3*V2*V2*W3+ V3*V2*V3*W2
      ROT(13,7) =V1*V1*V2*W3+ V1*V1*V3*W2+ V1*V2*V1*W3+ V1*V3*V1*W2
     1         + V1*V2*V3*W1+ V1*V3*V2*W1+ V2*V1*V1*W3+ V3*V1*V1*W2
     2         + V2*V1*V3*W1+ V3*V1*V2*W1+ V2*V3*V1*W1+ V3*V2*V1*W1
      ROT(14,7) =V2*V2*V1*W3+ V2*V2*V3*W1+ V2*V1*V2*W3+ V2*V3*V2*W1
     1         + V2*V1*V3*W2+ V2*V3*V1*W1+ V1*V2*V2*W3+ V3*V2*V2*W1
     2         + V1*V2*V3*W2+ V3*V2*V1*W2+ V1*V3*V2*W2+ V3*V1*V2*W2
      ROT(15,7) =V3*V3*V1*W2+ V3*V3*V2*W1+ V3*V1*V3*W2+ V3*V2*V3*W1
     1         + V3*V1*V2*W3+ V3*V2*V1*W3+ V1*V3*V3*W2+ V2*V3*V3*W1
     2         + V1*V3*V2*W3+ V2*V3*V1*W3+ V1*V2*V3*W3+ V2*V1*V3*W3
C     NEW ZZZX
      ROT(1,8)  =W1*W1*W1*U1
      ROT(2,8)  =W2*W2*W2*U2
      ROT(3,8)  =W3*W3*W3*U3
      ROT(4,8)  =W2*W1*W1*U1+ W1*W2*W1*U1+ W1*W1*W2*U1+ W1*W1*W1*U2
      ROT(5,8)  =W3*W1*W1*U1+ W1*W3*W1*U1+ W1*W1*W3*U1+ W1*W1*W1*U3
      ROT(6,8)  =W1*W2*W2*U2+ W2*W1*W2*U2+ W2*W2*W1*U2+ W2*W2*W2*U1
      ROT(7,8)  =W3*W2*W2*U2+ W2*W3*W2*U2+ W2*W2*W3*U2+ W2*W2*W2*U3
      ROT(8,8)  =W1*W3*W3*U3+ W3*W1*W3*U3+ W3*W3*W1*U3+ W3*W3*W3*U1
      ROT(9,8)  =W2*W3*W3*U3+ W3*W2*W3*U3+ W3*W3*W2*U3+ W3*W3*W3*U2
      ROT(10,8) =W1*W1*W2*U2+ W2*W2*W1*U1+ W1*W2*W1*U2+ W1*W2*W2*U1
     1         + W2*W1*W1*U2+ W2*W1*W2*U1
      ROT(11,8) =W1*W1*W3*U3+ W3*W3*W1*U1+ W1*W3*W1*U3+ W1*W3*W3*U1
     1         + W3*W1*W1*U3+ W3*W1*W3*U1
      ROT(12,8) =W2*W2*W3*U3+ W3*W3*W2*U2+ W2*W3*W2*U3+ W2*W3*W3*U2
     1         + W3*W2*W2*U3+ W3*W2*W3*U2
      ROT(13,8) =W1*W1*W2*U3+ W1*W1*W3*U2+ W1*W2*W1*U3+ W1*W3*W1*U2
     1         + W1*W2*W3*U1+ W1*W3*W2*U1+ W2*W1*W1*U3+ W3*W1*W1*U2
     2         + W2*W1*W3*U1+ W3*W1*W2*U1+ W2*W3*W1*U1+ W3*W2*W1*U1
      ROT(14,8) =W2*W2*W1*U3+ W2*W2*W3*U1+ W2*W1*W2*U3+ W2*W3*W2*U1
     1         + W2*W1*W3*U2+ W2*W3*W1*U2+ W1*W2*W2*U3+ W3*W2*W2*U1
     2         + W1*W2*W3*U2+ W3*W2*W1*U2+ W1*W3*W2*U2+ W3*W1*W2*U2
      ROT(15,8) =W3*W3*W1*U2+ W3*W3*W2*U1+ W3*W1*W3*U2+ W3*W2*W3*U1
     1         + W3*W1*W2*U3+ W3*W2*W1*U3+ W1*W3*W3*U2+ W2*W3*W3*U1
     2         + W1*W3*W2*U3+ W2*W3*W1*U3+ W1*W2*W3*U3+ W2*W1*W3*U3
C     NEW ZZZY
      ROT(1,9)  =W1*W1*W1*V1
      ROT(2,9)  =W2*W2*W2*V2
      ROT(3,9)  =W3*W3*W3*V3
      ROT(4,9)  =W2*W1*W1*V1+ W1*W2*W1*V1+ W1*W1*W2*V1+ W1*W1*W1*V2
      ROT(5,9)  =W3*W1*W1*V1+ W1*W3*W1*V1+ W1*W1*W3*V1+ W1*W1*W1*V3
      ROT(6,9)  =W1*W2*W2*V2+ W2*W1*W2*V2+ W2*W2*W1*V2+ W2*W2*W2*V1
      ROT(7,9)  =W3*W2*W2*V2+ W2*W3*W2*V2+ W2*W2*W3*V2+ W2*W2*W2*V3
      ROT(8,9)  =W1*W3*W3*V3+ W3*W1*W3*V3+ W3*W3*W1*V3+ W3*W3*W3*V1
      ROT(9,9)  =W2*W3*W3*V3+ W3*W2*W3*V3+ W3*W3*W2*V3+ W3*W3*W3*V2
      ROT(10,9) =W1*W1*W2*V2+ W2*W2*W1*V1+ W1*W2*W1*V2+ W1*W2*W2*V1
     1         + W2*W1*W1*V2+ W2*W1*W2*V1
      ROT(11,9) =W1*W1*W3*V3+ W3*W3*W1*V1+ W1*W3*W1*V3+ W1*W3*W3*V1
     1         + W3*W1*W1*V3+ W3*W1*W3*V1
      ROT(12,9) =W2*W2*W3*V3+ W3*W3*W2*V2+ W2*W3*W2*V3+ W2*W3*W3*V2
     1         + W3*W2*W2*V3+ W3*W2*W3*V2
      ROT(13,9) =W1*W1*W2*V3+ W1*W1*W3*V2+ W1*W2*W1*V3+ W1*W3*W1*V2
     1         + W1*W2*W3*V1+ W1*W3*W2*V1+ W2*W1*W1*V3+ W3*W1*W1*V2
     2         + W2*W1*W3*V1+ W3*W1*W2*V1+ W2*W3*W1*V1+ W3*W2*W1*V1
      ROT(14,9) =W2*W2*W1*V3+ W2*W2*W3*V1+ W2*W1*W2*V3+ W2*W3*W2*V1
     1         + W2*W1*W3*V2+ W2*W3*W1*V2+ W1*W2*W2*V3+ W3*W2*W2*V1
     2         + W1*W2*W3*V2+ W3*W2*W1*V2+ W1*W3*W2*V2+ W3*W1*W2*V2
      ROT(15,9) =W3*W3*W1*V2+ W3*W3*W2*V1+ W3*W1*W3*V2+ W3*W2*W3*V1
     1         + W3*W1*W2*V3+ W3*W2*W1*V3+ W1*W3*W3*V2+ W2*W3*W3*V1
     2         + W1*W3*W2*V3+ W2*W3*W1*V3+ W1*W2*W3*V3+ W2*W1*W3*V3
C     NEW XXYY
      ROT(1,10) =U1*U1*V1*V1
      ROT(2,10) =U2*U2*V2*V2
      ROT(3,10) =U3*U3*V3*V3
      ROT(4,10) =U2*U1*V1*V1+ U1*U2*V1*V1+ U1*U1*V2*V1+ U1*U1*V1*V2
      ROT(5,10) =U3*U1*V1*V1+ U1*U3*V1*V1+ U1*U1*V3*V1+ U1*U1*V1*V3
      ROT(6,10) =U1*U2*V2*V2+ U2*U1*V2*V2+ U2*U2*V1*V2+ U2*U2*V2*V1
      ROT(7,10) =U3*U2*V2*V2+ U2*U3*V2*V2+ U2*U2*V3*V2+ U2*U2*V2*V3
      ROT(8,10) =U1*U3*V3*V3+ U3*U1*V3*V3+ U3*U3*V1*V3+ U3*U3*V3*V1
      ROT(9,10) =U2*U3*V3*V3+ U3*U2*V3*V3+ U3*U3*V2*V3+ U3*U3*V3*V2
      ROT(10,10)=U1*U1*V2*V2+ U2*U2*V1*V1+ U1*U2*V1*V2+ U1*U2*V2*V1
     1         + U2*U1*V1*V2+ U2*U1*V2*V1
      ROT(11,10)=U1*U1*V3*V3+ U3*U3*V1*V1+ U1*U3*V1*V3+ U1*U3*V3*V1
     1         + U3*U1*V1*V3+ U3*U1*V3*V1
      ROT(12,10)=U2*U2*V3*V3+ U3*U3*V2*V2+ U2*U3*V2*V3+ U2*U3*V3*V2
     1         + U3*U2*V2*V3+ U3*U2*V3*V2
      ROT(13,10)=U1*U1*V2*V3+ U1*U1*V3*V2+ U1*U2*V1*V3+ U1*U3*V1*V2
     1         + U1*U2*V3*V1+ U1*U3*V2*V1+ U2*U1*V1*V3+ U3*U1*V1*V2
     2         + U2*U1*V3*V1+ U3*U1*V2*V1+ U2*U3*V1*V1+ U3*U2*V1*V1
      ROT(14,10)=U2*U2*V1*V3+ U2*U2*V3*V1+ U2*U1*V2*V3+ U2*U3*V2*V1
     1         + U2*U1*V3*V2+ U2*U3*V1*V2+ U1*U2*V2*V3+ U3*U2*V2*V1
     2         + U1*U2*V3*V2+ U3*U2*V1*V2+ U1*U3*V2*V2+ U3*U1*V2*V2
      ROT(15,10)=U3*U3*V1*V2+ U3*U3*V2*V1+ U3*U1*V3*V2+ U3*U2*V3*V1
     1         + U3*U1*V2*V3+ U3*U2*V1*V3+ U1*U3*V3*V2+ U2*U3*V3*V1
     2         + U1*U3*V2*V3+ U2*U3*V1*V3+ U1*U2*V3*V3+ U2*U1*V3*V3
C     NEW XXZZ
      ROT(1,11) =U1*U1*W1*W1
      ROT(2,11) =U2*U2*W2*W2
      ROT(3,11) =U3*U3*W3*W3
      ROT(4,11) =U2*U1*W1*W1+ U1*U2*W1*W1+ U1*U1*W2*W1+ U1*U1*W1*W2
      ROT(5,11) =U3*U1*W1*W1+ U1*U3*W1*W1+ U1*U1*W3*W1+ U1*U1*W1*W3
      ROT(6,11) =U1*U2*W2*W2+ U2*U1*W2*W2+ U2*U2*W1*W2+ U2*U2*W2*W1
      ROT(7,11) =U3*U2*W2*W2+ U2*U3*W2*W2+ U2*U2*W3*W2+ U2*U2*W2*W3
      ROT(8,11) =U1*U3*W3*W3+ U3*U1*W3*W3+ U3*U3*W1*W3+ U3*U3*W3*W1
      ROT(9,11) =U2*U3*W3*W3+ U3*U2*W3*W3+ U3*U3*W2*W3+ U3*U3*W3*W2
      ROT(10,11)=U1*U1*W2*W2+ U2*U2*W1*W1+ U1*U2*W1*W2+ U1*U2*W2*W1
     1         + U2*U1*W1*W2+ U2*U1*W2*W1
      ROT(11,11)=U1*U1*W3*W3+ U3*U3*W1*W1+ U1*U3*W1*W3+ U1*U3*W3*W1
     1         + U3*U1*W1*W3+ U3*U1*W3*W1
      ROT(12,11)=U2*U2*W3*W3+ U3*U3*W2*W2+ U2*U3*W2*W3+ U2*U3*W3*W2
     1         + U3*U2*W2*W3+ U3*U2*W3*W2
      ROT(13,11)=U1*U1*W2*W3+ U1*U1*W3*W2+ U1*U2*W1*W3+ U1*U3*W1*W2
     1         + U1*U2*W3*W1+ U1*U3*W2*W1+ U2*U1*W1*W3+ U3*U1*W1*W2
     2         + U2*U1*W3*W1+ U3*U1*W2*W1+ U2*U3*W1*W1+ U3*U2*W1*W1
      ROT(14,11)=U2*U2*W1*W3+ U2*U2*W3*W1+ U2*U1*W2*W3+ U2*U3*W2*W1
     1         + U2*U1*W3*W2+ U2*U3*W1*W2+ U1*U2*W2*W3+ U3*U2*W2*W1
     2         + U1*U2*W3*W2+ U3*U2*W1*W2+ U1*U3*W2*W2+ U3*U1*W2*W2
      ROT(15,11)=U3*U3*W1*W2+ U3*U3*W2*W1+ U3*U1*W3*W2+ U3*U2*W3*W1
     1         + U3*U1*W2*W3+ U3*U2*W1*W3+ U1*U3*W3*W2+ U2*U3*W3*W1
     2         + U1*U3*W2*W3+ U2*U3*W1*W3+ U1*U2*W3*W3+ U2*U1*W3*W3
C     NEW YYZZ
      ROT(1,12) =V1*V1*W1*W1
      ROT(2,12) =V2*V2*W2*W2
      ROT(3,12) =V3*V3*W3*W3
      ROT(4,12) =V2*V1*W1*W1+ V1*V2*W1*W1+ V1*V1*W2*W1+ V1*V1*W1*W2
      ROT(5,12) =V3*V1*W1*W1+ V1*V3*W1*W1+ V1*V1*W3*W1+ V1*V1*W1*W3
      ROT(6,12) =V1*V2*W2*W2+ V2*V1*W2*W2+ V2*V2*W1*W2+ V2*V2*W2*W1
      ROT(7,12) =V3*V2*W2*W2+ V2*V3*W2*W2+ V2*V2*W3*W2+ V2*V2*W2*W3
      ROT(8,12) =V1*V3*W3*W3+ V3*V1*W3*W3+ V3*V3*W1*W3+ V3*V3*W3*W1
      ROT(9,12) =V2*V3*W3*W3+ V3*V2*W3*W3+ V3*V3*W2*W3+ V3*V3*W3*W2
      ROT(10,12)=V1*V1*W2*W2+ V2*V2*W1*W1+ V1*V2*W1*W2+ V1*V2*W2*W1
     1         + V2*V1*W1*W2+ V2*V1*W2*W1
      ROT(11,12)=V1*V1*W3*W3+ V3*V3*W1*W1+ V1*V3*W1*W3+ V1*V3*W3*W1
     1         + V3*V1*W1*W3+ V3*V1*W3*W1
      ROT(12,12)=V2*V2*W3*W3+ V3*V3*W2*W2+ V2*V3*W2*W3+ V2*V3*W3*W2
     1         + V3*V2*W2*W3+ V3*V2*W3*W2
      ROT(13,12)=V1*V1*W2*W3+ V1*V1*W3*W2+ V1*V2*W1*W3+ V1*V3*W1*W2
     1         + V1*V2*W3*W1+ V1*V3*W2*W1+ V2*V1*W1*W3+ V3*V1*W1*W2
     2         + V2*V1*W3*W1+ V3*V1*W2*W1+ V2*V3*W1*W1+ V3*V2*W1*W1
      ROT(14,12)=V2*V2*W1*W3+ V2*V2*W3*W1+ V2*V1*W2*W3+ V2*V3*W2*W1
     1         + V2*V1*W3*W2+ V2*V3*W1*W2+ V1*V2*W2*W3+ V3*V2*W2*W1
     2         + V1*V2*W3*W2+ V3*V2*W1*W2+ V1*V3*W2*W2+ V3*V1*W2*W2
      ROT(15,12)=V3*V3*W1*W2+ V3*V3*W2*W1+ V3*V1*W3*W2+ V3*V2*W3*W1
     1         + V3*V1*W2*W3+ V3*V2*W1*W3+ V1*V3*W3*W2+ V2*V3*W3*W1
     2         + V1*V3*W2*W3+ V2*V3*W1*W3+ V1*V2*W3*W3+ V2*V1*W3*W3
C     NEW XXYZ
      ROT(1,13) =U1*U1*V1*W1
      ROT(2,13) =U2*U2*V2*W2
      ROT(3,13) =U3*U3*V3*W3
      ROT(4,13) =U2*U1*V1*W1+ U1*U2*V1*W1+ U1*U1*V2*W1+ U1*U1*V1*W2
      ROT(5,13) =U3*U1*V1*W1+ U1*U3*V1*W1+ U1*U1*V3*W1+ U1*U1*V1*W3
      ROT(6,13) =U1*U2*V2*W2+ U2*U1*V2*W2+ U2*U2*V1*W2+ U2*U2*V2*W1
      ROT(7,13) =U3*U2*V2*W2+ U2*U3*V2*W2+ U2*U2*V3*W2+ U2*U2*V2*W3
      ROT(8,13) =U1*U3*V3*W3+ U3*U1*V3*W3+ U3*U3*V1*W3+ U3*U3*V3*W1
      ROT(9,13) =U2*U3*V3*W3+ U3*U2*V3*W3+ U3*U3*V2*W3+ U3*U3*V3*W2
      ROT(10,13)=U1*U1*V2*W2+ U2*U2*V1*W1+ U1*U2*V1*W2+ U1*U2*V2*W1
     1         + U2*U1*V1*W2+ U2*U1*V2*W1
      ROT(11,13)=U1*U1*V3*W3+ U3*U3*V1*W1+ U1*U3*V1*W3+ U1*U3*V3*W1
     1         + U3*U1*V1*W3+ U3*U1*V3*W1
      ROT(12,13)=U2*U2*V3*W3+ U3*U3*V2*W2+ U2*U3*V2*W3+ U2*U3*V3*W2
     1         + U3*U2*V2*W3+ U3*U2*V3*W2
      ROT(13,13)=U1*U1*V2*W3+ U1*U1*V3*W2+ U1*U2*V1*W3+ U1*U3*V1*W2
     1         + U1*U2*V3*W1+ U1*U3*V2*W1+ U2*U1*V1*W3+ U3*U1*V1*W2
     2         + U2*U1*V3*W1+ U3*U1*V2*W1+ U2*U3*V1*W1+ U3*U2*V1*W1
      ROT(14,13)=U2*U2*V1*W3+ U2*U2*V3*W1+ U2*U1*V2*W3+ U2*U3*V2*W1
     1         + U2*U1*V3*W2+ U2*U3*V1*W2+ U1*U2*V2*W3+ U3*U2*V2*W1
     2         + U1*U2*V3*W2+ U3*U2*V1*W2+ U1*U3*V2*W2+ U3*U1*V2*W2
      ROT(15,13)=U3*U3*V1*W2+ U3*U3*V2*W1+ U3*U1*V3*W2+ U3*U2*V3*W1
     1         + U3*U1*V2*W3+ U3*U2*V1*W3+ U1*U3*V3*W2+ U2*U3*V3*W1
     2         + U1*U3*V2*W3+ U2*U3*V1*W3+ U1*U2*V3*W3+ U2*U1*V3*W3
C     NEW YYXZ
      ROT(1,14) =V1*V1*U1*W1
      ROT(2,14) =V2*V2*U2*W2
      ROT(3,14) =V3*V3*U3*W3
      ROT(4,14) =V2*V1*U1*W1+ V1*V2*U1*W1+ V1*V1*U2*W1+ V1*V1*U1*W2
      ROT(5,14) =V3*V1*U1*W1+ V1*V3*U1*W1+ V1*V1*U3*W1+ V1*V1*U1*W3
      ROT(6,14) =V1*V2*U2*W2+ V2*V1*U2*W2+ V2*V2*U1*W2+ V2*V2*U2*W1
      ROT(7,14) =V3*V2*U2*W2+ V2*V3*U2*W2+ V2*V2*U3*W2+ V2*V2*U2*W3
      ROT(8,14) =V1*V3*U3*W3+ V3*V1*U3*W3+ V3*V3*U1*W3+ V3*V3*U3*W1
      ROT(9,14) =V2*V3*U3*W3+ V3*V2*U3*W3+ V3*V3*U2*W3+ V3*V3*U3*W2
      ROT(10,14)=V1*V1*U2*W2+ V2*V2*U1*W1+ V1*V2*U1*W2+ V1*V2*U2*W1
     1         + V2*V1*U1*W2+ V2*V1*U2*W1
      ROT(11,14)=V1*V1*U3*W3+ V3*V3*U1*W1+ V1*V3*U1*W3+ V1*V3*U3*W1
     1         + V3*V1*U1*W3+ V3*V1*U3*W1
      ROT(12,14)=V2*V2*U3*W3+ V3*V3*U2*W2+ V2*V3*U2*W3+ V2*V3*U3*W2
     1         + V3*V2*U2*W3+ V3*V2*U3*W2
      ROT(13,14)=V1*V1*U2*W3+ V1*V1*U3*W2+ V1*V2*U1*W3+ V1*V3*U1*W2
     1         + V1*V2*U3*W1+ V1*V3*U2*W1+ V2*V1*U1*W3+ V3*V1*U1*W2
     2         + V2*V1*U3*W1+ V3*V1*U2*W1+ V2*V3*U1*W1+ V3*V2*U1*W1
      ROT(14,14)=V2*V2*U1*W3+ V2*V2*U3*W1+ V2*V1*U2*W3+ V2*V3*U2*W1
     1         + V2*V1*U3*W2+ V2*V3*U1*W2+ V1*V2*U2*W3+ V3*V2*U2*W1
     2         + V1*V2*U3*W2+ V3*V2*U1*W2+ V1*V3*U2*W2+ V3*V1*U2*W2
      ROT(15,14)=V3*V3*U1*W2+ V3*V3*U2*W1+ V3*V1*U3*W2+ V3*V2*U3*W1
     1         + V3*V1*U2*W3+ V3*V2*U1*W3+ V1*V3*U3*W2+ V2*V3*U3*W1
     2         + V1*V3*U2*W3+ V2*V3*U1*W3+ V1*V2*U3*W3+ V2*V1*U3*W3
C     NEW ZZXY
      ROT(1,15) =W1*W1*U1*V1
      ROT(2,15) =W2*W2*U2*V2
      ROT(3,15) =W3*W3*U3*V3
      ROT(4,15) =W2*W1*U1*V1+ W1*W2*U1*V1+ W1*W1*U2*V1+ W1*W1*U1*V2
      ROT(5,15) =W3*W1*U1*V1+ W1*W3*U1*V1+ W1*W1*U3*V1+ W1*W1*U1*V3
      ROT(6,15) =W1*W2*U2*V2+ W2*W1*U2*V2+ W2*W2*U1*V2+ W2*W2*U2*V1
      ROT(7,15) =W3*W2*U2*V2+ W2*W3*U2*V2+ W2*W2*U3*V2+ W2*W2*U2*V3
      ROT(8,15) =W1*W3*U3*V3+ W3*W1*U3*V3+ W3*W3*U1*V3+ W3*W3*U3*V1
      ROT(9,15) =W2*W3*U3*V3+ W3*W2*U3*V3+ W3*W3*U2*V3+ W3*W3*U3*V2
      ROT(10,15)=W1*W1*U2*V2+ W2*W2*U1*V1+ W1*W2*U1*V2+ W1*W2*U2*V1
     1         + W2*W1*U1*V2+ W2*W1*U2*V1
      ROT(11,15)=W1*W1*U3*V3+ W3*W3*U1*V1+ W1*W3*U1*V3+ W1*W3*U3*V1
     1         + W3*W1*U1*V3+ W3*W1*U3*V1
      ROT(12,15)=W2*W2*U3*V3+ W3*W3*U2*V2+ W2*W3*U2*V3+ W2*W3*U3*V2
     1         + W3*W2*U2*V3+ W3*W2*U3*V2
      ROT(13,15)=W1*W1*U2*V3+ W1*W1*U3*V2+ W1*W2*U1*V3+ W1*W3*U1*V2
     1         + W1*W2*U3*V1+ W1*W3*U2*V1+ W2*W1*U1*V3+ W3*W1*U1*V2
     2         + W2*W1*U3*V1+ W3*W1*U2*V1+ W2*W3*U1*V1+ W3*W2*U1*V1
      ROT(14,15)=W2*W2*U1*V3+ W2*W2*U3*V1+ W2*W1*U2*V3+ W2*W3*U2*V1
     1         + W2*W1*U3*V2+ W2*W3*U1*V2+ W1*W2*U2*V3+ W3*W2*U2*V1
     2         + W1*W2*U3*V2+ W3*W2*U1*V2+ W1*W3*U2*V2+ W3*W1*U2*V2
      ROT(15,15)=W3*W3*U1*V2+ W3*W3*U2*V1+ W3*W1*U3*V2+ W3*W2*U3*V1
     1         + W3*W1*U2*V3+ W3*W2*U1*V3+ W1*W3*U3*V2+ W2*W3*U3*V1
     2         + W1*W3*U2*V3+ W2*W3*U1*V3+ W1*W2*U3*V3+ W2*W1*U3*V3
C
C     ----- EXPRESS SPHERICAL HARMONICS IN ABSOLUTE FRAME -----
C
      DO 360 I=1,15
      DO 360 J=1,15
      A=ZERO
      DO 350 K=1,15
      A=A+ROT(I,K)*GSPH(K,J)
 350  CONTINUE
        BUFF(I,J)=A
 360  GSHELL(I,J)=A
C
C     ----- ACCOUNT FOR NORMALIZATION OF CARTESIAN FUNCTIONS -----
C
      IF(NORMF.EQ.0.OR.NORMP.EQ.0) THEN
      DO 370 I=1,15
      GSHELL( 4,I)=GSHELL( 4,I)/SQRT7
      GSHELL( 5,I)=GSHELL( 5,I)/SQRT7
      GSHELL( 6,I)=GSHELL( 6,I)/SQRT7
      GSHELL( 7,I)=GSHELL( 7,I)/SQRT7
      GSHELL( 8,I)=GSHELL( 8,I)/SQRT7
      GSHELL( 9,I)=GSHELL( 9,I)/SQRT7
      GSHELL(10,I)=GSHELL(10,I)*SQRT3/(SQRT7*SQRT5)
      GSHELL(11,I)=GSHELL(11,I)*SQRT3/(SQRT7*SQRT5)
      GSHELL(12,I)=GSHELL(12,I)*SQRT3/(SQRT7*SQRT5)
      GSHELL(13,I)=GSHELL(13,I)/(SQRT7*SQRT5)
      GSHELL(14,I)=GSHELL(14,I)/(SQRT7*SQRT5)
      GSHELL(15,I)=GSHELL(15,I)/(SQRT7*SQRT5)
        BUFF( 4,I)=GSHELL( 4,I)
        BUFF( 5,I)=GSHELL( 5,I)
        BUFF( 6,I)=GSHELL( 6,I)
        BUFF( 7,I)=GSHELL( 7,I)
        BUFF( 8,I)=GSHELL( 8,I)
        BUFF( 9,I)=GSHELL( 9,I)
        BUFF(10,I)=GSHELL(10,I)
        BUFF(11,I)=GSHELL(11,I)
        BUFF(12,I)=GSHELL(12,I)
        BUFF(13,I)=GSHELL(13,I)
        BUFF(14,I)=GSHELL(14,I)
        BUFF(15,I)=GSHELL(15,I)
  370 CONTINUE
      ENDIF
C
C     ----- DEFINE INVERSE TRANSFORMATION -----
C
      CALL HND_LUDCMP(BUFF,15,15,INDX,VV,D,IERR)
      DO 380 I=1,15
      CALL HND_LUBKSB(BUFF,15,15,INDX,GIHELL(1,I))
  380 CONTINUE
C
C     PRINT
C
      IF(DBUG) CALL HND_SYMPRP(GSHELL,15,15,15)
      IF(DBUG) CALL HND_SYMPRP(GIHELL,15,15,15)
      RETURN
 9999 FORMAT(/,' DEFINITION OF SPHERICAL HARMONICS',/,
     1         ' ---------------------------------')
      END
c
      SUBROUTINE HND_SYMPRP(V,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LISTNG/LIST
      DIMENSION V(NDIM,1)
      DIMENSION IC(5),C(5)
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
C
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.LT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      IF(LIST.EQ.0) WRITE (IW,9048) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      WRITE(IW,9328) J
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
 9328 FORMAT(/,' VEC =',I5,/,1X,10(1H-))
 9348 FORMAT(5(I5,F11.5))
      END
      SUBROUTINE HND_LUDCMP(A,LDA,N,IB,VV,D,IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*8 ERRMSG
      DIMENSION A(LDA,1),IB(N),VV(N)
      DIMENSION ERRMSG(3)
      DATA ERRMSG   /'PROGRAM ','STOP IN ','-LUDCMP-'/
      DATA TINY     /1.0D-20/
      DATA ZERO,ONE /0.0D+00,1.0D+00/
C
      IERR=0
      D=ONE
      DO 12 I=1,N
         AAMAX=ZERO
         DO 11 J=1,N
            IF( ABS(A(I,J)).GT.AAMAX) AAMAX= ABS(A(I,J))
   11    CONTINUE
         IF(AAMAX.EQ.ZERO) THEN
            IERR=1
            RETURN
         ENDIF
         VV(I)=ONE/AAMAX
   12 CONTINUE
      DO 19 J=1,N
         IF(J.GT.1) THEN
            DO 14 I=1,J-1
               SUM=A(I,J)
               IF(I.GT.1) THEN
                  DO 13 K=1,I-1
                     SUM=SUM-A(I,K)*A(K,J)
   13             CONTINUE
                  A(I,J)=SUM
               ENDIF
   14       CONTINUE
         ENDIF
         AAMAX=ZERO
         DO 16 I=J,N
            SUM=A(I,J)
            IF(J.GT.1) THEN
               DO 15 K=1,J-1
                  SUM=SUM-A(I,K)*A(K,J)
   15          CONTINUE
               A(I,J)=SUM
            ENDIF
            DUM=VV(I)* ABS(SUM)
            IF(DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
   16    CONTINUE
         IF(J.NE.IMAX) THEN
            DO 17 K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
   17       CONTINUE
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
         IB(J)=IMAX
         IF(J.NE.N) THEN
            IF(A(J,J).EQ.ZERO) A(J,J)=TINY
            DUM=ONE/A(J,J)
            DO 18 I=J+1,N
               A(I,J)=A(I,J)*DUM
   18       CONTINUE
         ENDIF
   19 CONTINUE
      IF(A(N,N).EQ.ZERO) A(N,N)=TINY
      RETURN
      END
      SUBROUTINE HND_LUBKSB(A,LDA,N,IB,B)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(LDA,1),IB(N),B(N)
      DATA ZERO /0.0D+00/
C
      II=0
      DO 12 I=1,N
         LL=IB(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF(II.NE.0) THEN
            DO 11 J=II,I-1
               SUM=SUM-A(I,J)*B(J)
   11       CONTINUE
         ELSE IF (SUM.NE.ZERO) THEN
            II=I
         ENDIF
         B(I)=SUM
   12 CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         IF(I.LT.N) THEN
            DO 13 J=I+1,N
               SUM=SUM-A(I,J)*B(J)
   13       CONTINUE
         ENDIF
         B(I)=SUM/A(I,I)
   14 CONTINUE
      RETURN
      END
      SUBROUTINE HND_ATLABL(rtdb,geom,basis)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "nwc_const.fh"
      integer rtdb,geom,basis
      parameter (max_center = 500)
      integer          ic 
      integer          nat !counts no. of centers as input
      character*16     tag(max_center)
      character*2      symbol
      character*16     element
      integer          atn
      logical          status
      logical          geom_tag_to_element
      external         geom_tag_to_element
      logical          geom_cent_tag
      external         geom_cent_tag
      logical          debug
      character*8 errmsg
      CHARACTER*4 BFNAM(35)
      CHARACTER*4 ATMLAB(106),BONDF,LABEL
      CHARACTER*8 BFL
      CHARACTER*2 SLAB
      CHARACTER*8 ANAM
      character*8 bflab
      character*2 blank
      CHARACTER*2 BNAM
      PARAMETER (MXPRIM=2048)
      PARAMETER (MXSHEL=512)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      PARAMETER (MXBFN =3072)
      PARAMETER (SLAB=' S')
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_basnum/num
      common/hnd_nshel/ex(mxprim),cs(mxprim),cp(mxprim),
     1                 cd(mxprim),cf(mxprim),cg(mxprim),
     2                 kstart(mxshel),katom(mxshel),ktype(mxshel),
     3                 kng(mxshel),kloc(mxshel),kmin(mxshel),
     4                 kmax(mxshel),nshell
      common/hnd_mollab/anam(mxatom),bnam(mxatom),bflab(mxbfn)
C     COMMON/hnd_MOLNUC/NUC(MXATOM)
      common/hnd_iofile/ir,iw,ip
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ATLABL-'/
C
      data bfnam/'   S','   X','   Y','   Z',
     1           '  XX','  YY','  ZZ','  XY','  XZ','  YZ',
     2           ' XXX',' YYY',' ZZZ',' XXY',' XXZ',
     3           ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',
     4           'XXXX','YYYY','ZZZZ','XXXY','XXXZ',
     5           'YYYX','YYYZ','ZZZX','ZZZY','XXYY',
     6           'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/
      DATA BONDF/' BF '/
C
c                           104 truel elements, and 2 dummies
c
      data atmlab/'H   ','HE  ','LI  ','BE  ','B   ','C   ',
     1            'N   ','O   ','F   ','NE  ','NA  ','MG  ',
     2            'AL  ','SI  ','P   ','S   ','CL  ','AR  ',
     3            'K   ','CA  ','SC  ','TI  ','V   ','CR  ',
     4            'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ',
     5            'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ',
     6            'RB  ','RU  ','Y   ','ZR  ','NB  ','MO  ',
     7            'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ',
     8            'IN  ','SN  ','SB  ','TE  ','I   ','XE  ',
     9            'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ',
     1            'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ',
     2            'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ',
     3            'TA  ','W   ','RE  ','OS  ','IR  ','PT  ',
     4            'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ',
     5            'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ',
     6            'PA  ','U   ','NP  ','PU  ','AM  ','CM  ',
     7            'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ',
     8            'LR  ','RF  ','X   ','BQ  '/
C
      data blank/'  '/
c
c         ----- fill basis function symbol table -----
c
      debug = .false.
c
      n = 0
      lastat=0
      do 120 ii = 1,nshell
         iat = katom(ii)
         j   = zan(iat)
         if (j .gt. 106) j = 106
         if (j .gt. 0) label=atmlab(j)
         if (j .le. 0) label=bondf
         mini = kmin(ii)
         maxi = kmax(ii)
         do 110 i = mini,maxi
            n = n+1
            if ((lastat.eq.iat).or.(ktype(ii).ne.1)) then
               write(unit=bfl,fmt='(a4,a4)') label,bfnam(i)
            else
               write(unit=bfl,fmt='(a3,i3,a2)') label,iat,slab
            endif
            read (unit=bfl,fmt='(a8)') bflab(n)
  110       continue
         lastat=iat
  120    continue

         if(debug.and.ga_nodeid().eq.0) then
            do i=1,num
               print*,'bflab(i)=',bflab(i)
            enddo
         endif
c
c        ----- fill atom symbol table -----
c

      do ic=1,nat
         status = geom_cent_tag(geom,ic,tag(ic))
         status = geom_tag_to_element(tag(ic),symbol,element,atn)


         anam(ic)=tag(ic)
         bnam(ic)=blank

         if(debug.and.ga_nodeid().eq.0) then
            print*,'ic=',ic
            print*,'tag(ic)=',tag(ic)
            print*,'symbol=',symbol
            print*,'element=',element
            print*,'atn=',atn
            print*,'anam(ic)=',anam(ic)
         endif
      enddo
C
      RETURN
      END

c*    SUBROUTINE HND_WRDNAM(WRDKEY)
c*    IMPLICIT REAL*8 (A-H,O-Z)
c*    CHARACTER*9 WORD
c*    WRITE(WORD,FMT='(1X,A8)') WRDKEY
c*     READ(WORD,FMT='(   A8)') WRDKEY
c*     RETURN
c*    END
      subroutine hnd_get_eigens(rtdb,geom,basis,vecs,e,num)
c
c     only RHF(close shell) is considered 
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer           mxnbf
      parameter         (mxnbf=3072)
      integer           num
      integer           rtdb,geom,basis
      logical           status
      integer           g_vecs
      character*255     movecs_in,movecs_out
      character*255     movecs
      logical           movecs_read
      external          movecs_read
      integer           ga_create_atom_blocked
      external          ga_create_atom_blocked
      double precision  vecs(num*num) 
      double precision  e(num)
      double precision  occ(mxnbf)
c
      integer          i_init,init,need 
      integer          num3
      integer          i
c*    integer          i,j,ij
      integer          i00,i10
c
      double precision x
      double precision zero
c
      dimension    x(1)
      equivalence (x(1),dbl_mb(1))
c
      logical          debug
      data zero   /0.0d+00/
c
      debug=.true. 
c
      num3=num*num
c     
      need=num3+num3
c
c     allocate memory block
c
      if(.not.ma_push_get(mt_dbl,need ,
     &    'mem init:prp:hnd_get_eigens:9',
     &    i_init,init))
     & call errquit('hnd_get_eigens, malloc of init  failed',911)
c
      i00=init
      i10=i00+ num3
c
      do 10 i=1,need
   10 x(i+init-1)=zero
c
c     get mo vectors from file 
c
      g_vecs=ga_create_atom_blocked (geom, basis, 'mo vectors')
c
c     check output input/output MO vectors 
c
       if (.not.rtdb_cget(rtdb,'scf:input  vectors',1,movecs_in))
     $      movecs_in = 'atomic'
       if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_out))
     $      movecs_out = ' '
       movecs = ' '
       if (movecs.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore') then
            movecs_out = 'movecs'
         else 
            movecs_out = movecs_in
         endif
       endif
       movecs = movecs_out
c
c      read eig-vec=>x(i00); eng-val=>x(i10) 
c
       status = movecs_read (movecs, 1, occ, e, g_vecs)
       call ga_get(g_vecs, 1, num, 1, num,x(i00), num)
       status = ga_destroy ( g_vecs )
c
c      tranform eigen-V matrix from nwchem to hondo shell order
c
      if (debug) print*,'eigen vectors in nwchem format'
      if (debug) call hnd_sprtr(x(i00),num)
      call hnd_nwhnd_tran_vecs(x(i00),vecs,num)
      if (debug) print*,'eigen vectors in hondo  format'
      if (debug) call hnd_sprtr(vecs,num)
c
c     release memory block         
c
      if (.not.ma_pop_stack(i_init))
     & call errquit('hnd_get_eigens, ma_pop_stack of init
     & failed',911)
c
      return
      end
c
      subroutine hnd_nwhnd_tran_vecs(dnwc,dhnd,ndim)
      implicit double precision (a-h,o-z)
#include "global.fh"
      parameter (mxnbf =2048)
      common/hnd_nwtohnd/inw_to_hnd(mxnbf)
      common/hnd_iofile/ir,iw,ip
      common/hnd_facntoh/fac_nwthnd(mxnbf)
      dimension dnwc(ndim,ndim),dhnd(ndim,ndim)
c
      logical out
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... nwhnd_tran_vecs '
      endif
c
c     ----- eigen-vectors from nwchem to hondo shell order-----
c
      do j=1,ndim
        do i=1,ndim

           dhnd(inw_to_hnd(i),j)=
     &     dnwc(i,j)/(fac_nwthnd(i))

        enddo
      enddo
c
      return 
      end
      subroutine hnd_nwhnd_tran_s(dnwc,dhnd,ndim)
      implicit double precision (a-h,o-z)
#include "global.fh"
      parameter (mxnbf =2048)
      common/hnd_nwtohnd/inw_to_hnd(mxnbf)
      common/hnd_iofile/ir,iw,ip
      common/hnd_facntoh/fac_nwthnd(mxnbf)
      dimension dnwc(ndim,ndim),dhnd(ndim,ndim)
c
      logical out
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... nwhnd_tran ... '
      endif
c
c     ----- matrices tran. from nwchem to hondo shell order-----
c
      do j=1,ndim
        do i=1,ndim

           dhnd(inw_to_hnd(i),inw_to_hnd(j))=
     &     dnwc(i,j)*fac_nwthnd(i)*fac_nwthnd(j)

        enddo
      enddo
c
      return 
      end
      subroutine hnd_prp_get_dens(rtdb,basis,geom,
     1                            da,db,d1nw,d2nw,
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "schwarz.fh"
#include "util.fh"
c
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      double precision da,db,d1nw,d2nw
      integer               nclosed ,nopen ,nvirt 
      integer               nclosedu,nopenu,nvirtu
c
      logical  int_normalize
      external int_normalize

      integer       rtdb              ! [input] database handle
      integer       geom, basis       ! handles
      logical       status
      integer       g_dens
      character*255 title
      character*255 movecs_in, movecs_out
      character*8   scftyp 
      dimension     g_dens(3)         ! max # of density matrices
      dimension     da(*),db(*),d1nw(*),d2nw(*)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title)) title = ' '
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'prop:input vectors',  1, 
     $     movecs_in )) then
           call util_file_name('movecs', .false.,.false.,movecs_in )
      endif
      call util_file_name_resolve(movecs_in , .false.)
      if (.not. rtdb_cget(rtdb, 'prop:output vectors',  1,
     $     movecs_out)) then
           call util_file_name('movecs', .false.,.false.,movecs_out)
      endif
      call util_file_name_resolve(movecs_out, .false.)
c
      call hnd_prp_dens_read(rtdb,basis,geom,movecs_in,
     1                       da,db,d1nw,d2nw,
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
c
      return
      end
      subroutine hnd_prp_dens_read(rtdb,basis,geom,movecs,
     1                             da,db,d1nw,d2nw,
     2                      scftyp,
     3                      nclosed ,nopen ,nvirt ,
     4                      nclosedu,nopenu,nvirtu)
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"

      logical       status
      integer       rtdb, basis, geom
      character*(*) movecs
      character*4   scftype
      character*8   scftyp 
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      logical  movecs_read_header, movecs_read
      external movecs_read_header, movecs_read

      integer g_dens,     ! density matrices (up to 3)
     $        g_molec 
  
      integer max_sh_bf, max_at_bf 
      integer nat, nbf, nopen , nclosed , nvirt 
      integer           nopenu, nclosedu, nvirtu  
      integer ndens, nbf2, nocc, i, iset
      character*255 title_molec, basis_name_molec
      character*20  scftype_molec
      integer       nbf_molec, nsets_molec, nmo_molec(2)
      integer                               nmo
c
      integer   l_eval, k_eval, l_occ, k_occ
      dimension g_dens(3) ! max # of density matrices
c
      double precision da(*),db(*),d1nw(*),d2nw(*)
c
      logical dbug
      dbug=.false.
c
c     get information about basis set
c
      status = geom_ncent(geom,nat)
      if ( .not.status )
     & call errquit('dens_read: could not get natoms',0)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) 
     & call errquit('dens_read: could not get nbf',0)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status )
     & call errquit('dens_read: could not get max_sh_bf',0)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status )
     & call errquit('dens_read: could not get max_at_bf',0)
c
c     read vector file header
c
      if (.not. movecs_read_header(movecs, title_molec,
     &    basis_name_molec, scftype_molec, nbf_molec,
     &    nsets_molec, nmo_molec, 2))
     &    call errquit('dens_read: failed to read header', 911)
      if(dbug) then
         write(6,9999)
     1     movecs
     2   ,title_molec
     3   ,basis_name_molec
     4   ,scftype_molec
     5   ,nbf_molec
     6   ,nsets_molec
     7   ,nmo_molec
      endif
c
      if (nsets_molec.eq.1) then
         scftype = 'RHF'
         scftyp  = 'RHF     '
      else
         scftype = 'UHF'
         scftyp  = 'UHF     '
      endif
      call inp_ucase(scftype)

C     # of eigenvalues and density matrices
      if ( scftype .eq. 'UHF' ) then
        nbf2  = 2 * nbf
        ndens = 3
        nmo   = max(nmo_molec(1),nmo_molec(2))
      else if ( scftype .eq. 'RHF' ) then
        nbf2  = nbf
        ndens = 1
        nmo   =     nmo_molec(1)
      endif

C     allocate and initialize global and local memory

C     occupation numbers
      status = MA_push_get(MT_DBL, nbf2,'occ. numbers', l_occ, 
     $                                                  k_occ)
      if (.not.status) 
     $     call errquit('property: could not allocate l_occ',911)
C     eigenvalues
      status = MA_push_get(MT_DBL, nbf2,'MO eval', l_eval, k_eval)
      if (.not.status) 
     $     call errquit('property: could not allocate l_eval',911)
      
      status = ga_create(MT_DBL, nbf_molec, nmo,' molec MOs',
     $                           nbf_molec, nmo,g_molec)
      if (.not.status)
     $     call errquit('property: GA for molecular MOs',911)

      if(dbug) then
         write(6,*) 'in hnd_prop_dens_read'
         write(6,*) 'scftype',scftype
         write(6,*) 'nbf,nbf2',nbf,nbf2
         write(6,*) 'nmo,ndens',nmo,ndens
      endif

c
c     read ( alpha ) orbitals
c
      iset=1
      status = movecs_read(movecs,iset,
     $         dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
      if (.not.status)
     $  call errquit('property: could not read molecular MOs',911)
       
      nclosed = 0
      nopen   = 0
      nvirt   = 0
      do i=1,nmo_molec(1)
         if(dbl_mb(k_occ+i-1).eq.2.0d0) then
            nclosed=nclosed+1
         elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
            nopen=nopen+1
         else
            nvirt=nvirt+1
         endif
      enddo
      if(dbug) then
         write(6,9998) nclosed , nopen , nvirt
      endif
c
c     read ( beta ) orbitals
c
      if( nsets_molec.gt.1 ) then
            iset=2
         status = movecs_read(movecs,iset,
     $            dbl_mb(k_occ),dbl_mb(k_eval),g_molec)
         if (.not.status)
     $     call errquit('property: could not read molecular MOs',911)

         nclosedu= 0
         nopenu  = 0
         nvirtu  = 0
         do i=1,nmo_molec(2)
            if(dbl_mb(k_occ+i-1).eq.2.0d0) then
               nclosedu=nclosedu+1
            elseif(dbl_mb(k_occ+i-1).ne.0.0d0) then
               nopenu=nopenu+1
            else
               nvirtu=nvirtu+1
            endif
         enddo
         if(dbug) then
            write(6,9997) nclosedu, nopenu, nvirtu
         endif
      endif

      status = ga_destroy ( g_molec )
c     global density
      do i=1, ndens
        g_dens(i) = ga_create_atom_blocked ( geom, basis, 
     $              'density matrix' )
      enddo

      if(dbug) then
         write(6,*) 'in hnd_prop_dens_read'
         write(6,*) 'scftype',scftype
         write(6,*) 'nbf,nbf2',nbf,nbf2
         write(6,*) 'nclosed,nopen',nclosed,nopen
         write(6,*) 'nclosedu,nopenu',nclosedu,nopenu
      endif

c
      call hnd_prop_dens_make (geom, basis, g_dens, 
     $     dbl_mb(k_eval), dbl_mb(k_occ), ndens, nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     da, db, d1nw, d2nw)

c
c     free temporary arrays
c
      status = MA_pop_stack ( l_eval)
      status = MA_pop_stack ( l_occ )

      return
 9999 format(
     1   /,'movecs             = ',a24,
     2   /,'title_molec        = ',a24,
     3   /,'basis_name_molec   = ',a24,
     4   /,'scftype_molec      = ',a24,
     5   /,'nbf_molec          = ',i12,
     6   /,'nsets_molec        = ',i12,
     7   /,'nmo_molec          = ',i12,i12
     8      )
 9998 format(' nclosed , nopen , nvirt  = ',3i10)
 9997 format(' nclosedu, nopenu, nvirtu = ',3i10)
      end
      subroutine hnd_prop_dens_make (geom, basis, g_dens, 
     $     evals, occ, ndens, nbf, nbf2, 
     $     nclosed, nopen, nclosedu, nopenu, scftype, movecs,
     $     da, db, d1nw, d2nw )

c     calculate the density matrices from the movecs 

      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer geom, basis, g_dens, ndens
      integer nbf, nbf2, nclosed, nopen, nclosedu, nopenu
      double precision evals, occ
      double precision da, db, d1nw, d2nw
      character*4 scftype
      character*(*) movecs

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical  movecs_read, file_read_ga
      external movecs_read, file_read_ga

      dimension evals(nbf2), occ(nbf2), g_dens(ndens)
      dimension da(*), db(*), d1nw(*), d2nw(*)

C     local variables
      integer nocc, l_ind, k_ind, l_ind2, k_ind2, g_vecs, g_vecs2
      integer i, j, ij
      integer iopen
      logical status
c
      logical dbug
      dbug=.false.

      if(dbug) then
         write(6,*) 'in hnd_prop_dens_make'
         write(6,*) 'scftype',scftype
         write(6,*) 'nbf,nbf2',nbf,nbf2
         write(6,*) 'nclosed,nopen',nclosed,nopen
         write(6,*) 'nclosedu,nopenu',nclosedu,nopenu
      endif

C  get MO vectors from file 

      g_vecs = ga_create_atom_blocked ( geom, basis, 'mo vectors' )
      if ( scftype.eq. 'UHF' ) then
        g_vecs2 = ga_create_atom_blocked (geom, basis, 'uhf mo vectors')
      endif

      status = movecs_read (movecs, 1, occ, evals, g_vecs)
      if (.not.status) then
        call errquit('prop_dens_make: could not read mo vectors', 110)
      endif

      if (scftype .eq. 'UHF') then
        if (.not. movecs_read(movecs, 2, occ(nbf+1), evals(nbf+1),
     $       g_vecs2)) then
        call errquit('prop_dens_make: could not read uhf vectors', 110)
        endif
      endif

      if (scftype .eq. 'RHF') then
c       density matrix

        if(nclosed.eq.0) then
           call ga_zero(g_dens(1))
        else
           call ga_matmul_patch('n', 't', 2.0D0, 0.0d0  ,
     $          g_vecs,    1, nbf,            1, nclosed,
     $          g_vecs,    1, nclosed,        1, nbf    ,
     $          g_dens(1), 1, nbf,            1, nbf    )
        endif
        if (nopen .gt. 0) then
           do iopen=1,nopen
             call ga_matmul_patch('n', 't', occ(nclosed+iopen), 1.0d0,
     $       g_vecs,1, nbf, nclosed+iopen, nclosed+iopen,
     $       g_vecs,nclosed+iopen, nclosed+iopen,1,nbf  ,
     $       g_dens(1), 1, nbf,              1,nbf  )
           enddo
        endif
        call ga_symmetrize(g_dens(1))

        call ga_get(g_dens(1),1,nbf,1,nbf,d1nw,nbf)
 
c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              da(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
           call hnd_prtr(da,nbf)
        endif

      else  ! UHF
c       density matrix  D_up and D_down

        call ga_zero(g_dens(1))
        call ga_zero(g_dens(2))
        if(nopen .gt.0) then
           do iopen=1,nopen
              call ga_matmul_patch('n', 't', occ(iopen), 1.0d0,
     $           g_vecs,    1, nbf,        iopen , iopen          ,
     $           g_vecs,iopen, iopen          , 1, nbf,
     $           g_dens(1), 1, nbf,             1, nbf)
           enddo
        endif
        if(nopenu.gt.0) then
           do iopen=1,nopenu
              call ga_matmul_patch('n', 't', occ(nbf+iopen), 1.0d0,
     $           g_vecs2,    1, nbf,       iopen , iopen          ,
     $           g_vecs2,iopen,iopen          , 1, nbf,
     $           g_dens(2),  1, nbf,            1, nbf)
           enddo
        endif

C       density total

        call ga_dadd ( 1.D0, g_dens(1), 1.D0, g_dens(2), g_dens(3) )
        call ga_symmetrize(g_dens(1))
        call ga_symmetrize(g_dens(2))
        call ga_symmetrize(g_dens(3))

        call ga_get(g_dens(1),1,nbf,1,nbf,d1nw,nbf)
c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              da(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
           call hnd_prtr(da,nbf)
        endif

        call ga_get(g_dens(2),1,nbf,1,nbf,d1nw,nbf)
c       hondo to nwchem mapping

        call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              db(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
           call hnd_prtr(db,nbf)
        endif

c       
      endif

c

      status = ga_destroy ( g_vecs )
      if ( scftype .eq. 'UHF' ) then
        status = ga_destroy ( g_vecs2 )
      endif
      do i=1,3
         status = ga_destroy ( g_dens(i) )
      enddo
      return
      end
      subroutine hnd_prp_init(geom,basis)
      implicit double precision (a-h,o-z)
c 
#include "nwc_const.fh"
#include "hnd_tol.fh"
      parameter (zero=0.0d+00)
      parameter (mxatom=nw_max_atom)
      parameter (mxprim=2048)
      parameter (mxshel=512)
      parameter (mxsym=48)
      parameter (mxnbf =2048)
      parameter (mxtype=5)
      logical  status
      integer  geom 
      integer  basis
      integer  num
      logical  geom_ncent
      external geom_ncent
      logical  geom_create
      external geom_create
      logical  geom_rtdb_load
      external geom_rtdb_load
      logical  geom_cart_get
      external geom_cart_get
      logical  bas_numcont
      external bas_numcont
      logical  bas_numbf
      external bas_numbf
      logical  bas_ce2cnr
      external bas_ce2cnr
      logical  bas_cn2bfr
      external bas_cn2bfr
      character*16 atom_tags
c
      logical dbug
      logical out
      dimension atom_tags(mxatom)
      dimension mini(mxtype),maxi(mxtype)
c
c     ----- routine initializes HONDO common blocks -----
c
      common/inp_file/iread,iwrite
      common/hnd_iofile/ir,iw,ip
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_basnum/num
      common/hnd_nshel/ex(mxprim),cs(mxprim),cp(mxprim),
     1                 cd(mxprim),cf(mxprim),cg(mxprim),
     2                 kstart(mxshel),katom(mxshel),ktype(mxshel),
     3                 kng(mxshel),kloc(mxshel),kmin(mxshel),
     4                 kmax(mxshel),nshell
      common/hnd_shlord/modshl(mxshel),invshl(mxshel)
      common/hnd_lcapid/nap,iap
      common/hnd_symtry/invt(48),nt,ntmax,ntwd,nosym
      common/hnd_symshl/mapshl(mxshel,mxsym),mapctr(mxatom,mxsym)
c
#include "global.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"

      data mini / 1, 2, 5,11,21/
      data maxi / 1, 4,10,20,35/

#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
      dbug=.false.
      out =.false.
      ir=iread
      iw=iwrite
      ip=7
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... hnd_prp_init ... '
      endif
      if(dbug) then
         write(iw,*) ' ... geom  = ',geom
         write(iw,*) ' ... basis = ',basis
      endif
c
c     ----- get basis set info + parameters -----
c
      if (.not. geom_ncent (geom , natoms)) call errquit
     $     ('hnd_prp_init: geom_ncent  failed for geom', geom)
      if(dbug) then
         write(iw,*) ' ... natoms = ', natoms
      endif
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('hnd_prp_init: bas_numcont failed for basis', basis)

      if (.not. bas_numbf(basis, num)) call errquit
     $     ('hnd_prp_init: bas_numbf failed for basis', basis)
      if(dbug) then
         write(iw,*) ' ... nshell = ', nshell
         write(iw,*) ' ... num    = ', num
      endif

c
c     ----- loop over atoms -----
c
      istart=1
      do iat=1,natoms
         status=bas_ce2cnr(basis,iat,ishlo,ishhi)
         if(.not.status) call errquit
     $       ('hnd_prp_init: bas_ce to cnr/bfr failed', iat)
         if(dbug) then
            write(iw,*) ' ... iat, ishlo, ishhi = ',iat,ishlo,ishhi
         endif
c 
c     ----- loop over shells on atoms -----
c
         do ish=ishlo,ishhi
            status=bas_cn2bfr(basis,ish,ishbflo,ishbfhi)
            if(.not.status) call errquit
     $          ('hnd_prp_init: bas_cn to bfr failed', ish)
            ioff=basis+BASIS_HANDLE_OFFSET
            ishu=sf_ibs_cn2ucn(ish,ioff)
            ityp=infbs_cont(CONT_TYPE ,ishu,ioff)
            ngi =infbs_cont(CONT_NPRIM,ishu,ioff)
            icon=infbs_cont(CONT_ICFP ,ishu,ioff)
            iexp=infbs_cont(CONT_IEXP ,ishu,ioff)
c
            katom(ish) =iat
            kloc(ish)  =ishbflo
            ktype(ish) =ityp+1
            kmin(ish)  =mini(ityp+1)
            kmax(ish)  =maxi(ityp+1)
            kng(ish)   =ngi
            kstart(ish)=istart
c
c     ----- loop over primitives in shells -----
c
            do ig=1,ngi
               ex(ig+istart-1)=zero
               cs(ig+istart-1)=zero
               cp(ig+istart-1)=zero
               cd(ig+istart-1)=zero
               cf(ig+istart-1)=zero
               cg(ig+istart-1)=zero
               cont_expo     =dbl_mb(mb_exndcf(ig+iexp-1,ioff))
               cont_coeff    =dbl_mb(mb_exndcf(ig+icon-1,ioff))
               ex(ig+istart-1)=cont_expo
               if(ityp.eq.0) then
                  cs(ig+istart-1)=cont_coeff
               elseif(ityp.eq.1) then
                  cp(ig+istart-1)=cont_coeff
               elseif(ityp.eq.2) then
                  cd(ig+istart-1)=cont_coeff
               elseif(ityp.eq.3) then
                  cf(ig+istart-1)=cont_coeff
               elseif(ityp.eq.4) then
                  cg(ig+istart-1)=cont_coeff
               else
                  call errquit('hnd_prp_init: illegal -ityp- ', ityp)
               endif
               if(dbug) then
                  write(iw,*) ' ish,ig,ityp,exp,con= ',ish,ig,ityp,
     $                         cont_expo,cont_coeff
               endif
            enddo
            istart=istart+ngi
c
         enddo
c
      enddo
c
c     ----- print /nshel/ -----
c
      if(dbug) then
         write(iw,*) ' ... kstart= ',(kstart(ish),ish=1,nshell)
         write(iw,*) ' ... katom = ',( katom(ish),ish=1,nshell)
         write(iw,*) ' ... ktype = ',( ktype(ish),ish=1,nshell)
         write(iw,*) ' ... kng   = ',(   kng(ish),ish=1,nshell)
         write(iw,*) ' ... kloc  = ',(  kloc(ish),ish=1,nshell)
         write(iw,*) ' ... kmin  = ',(  kmin(ish),ish=1,nshell)
         write(iw,*) ' ... kmax  = ',(  kmax(ish),ish=1,nshell)
      endif
c
c     ----- /baspar/ -----
c
      normf=0
      normp=0
      itol=15
      if(dbug) then
         write(iw,*) ' ... normf,normp,itol = ',normf,normp,itol
      endif
c
c     ----- /molxyz/ -----
c
      status=geom_ncent(geom,nat)
      status=geom_cart_get(geom,nat,atom_tags,c,zan)
      if(.not.status) call errquit
     $    (' int_1e_ga_hnd_init: geom_cart_get failed.',911)
      do iat=1,nat
         nuc(iat)=nint(zan(iat))
      enddo
      if(dbug) then
         do iat=1,nat
            write(iw,*) ' ... i,Z,c = ',iat,zan(iat),c(1,iat),c(2,iat),
     $                                               c(3,iat)
         enddo
      endif
c
c     ----- initialization for 2e integral derivatives -----
c
      do i=1,mxshel
         modshl(i)=i
         invshl(i)=i
      enddo
c
c     nap = 1
c     iap = 1
      nap = ga_nnodes()
      iap = ga_nodeid() + 1
c
      nt=1
      do i=1,mxshel
         mapshl(i,1)=i
      enddo
c
      do i=1,mxatom
         mapctr(i,1)=i
      enddo
c
c     ----- load Hermite's roots and weights -----
c
      call hnd_whermt
c
      return
      end
      subroutine hnd_prp_nw_hnd
      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
c
c     ----- order of components in shells in HONDO -----
c
c            x    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
c    1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
c    2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
c    3              3, 1, 3, 2, 2/
c            y    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
c    1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
c    2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
c    3              1, 3, 2, 3, 2/
c            z    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
c    1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
c    2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
c    3              3, 3, 2, 2, 3/
c
c     ----- order of components in shells in NWCHEM -----
c
c            x    / 1, 2, 1, 1, 3, 2, 2, 1, 1, 1,
c    1              4, 3, 3, 2, 2, 2, 1, 1, 1, 1,
c    2              5, 4, 4, 3, 3, 3, 2, 2, 2, 2,
c    3              1, 1, 1, 1, 1/
c            y    / 1, 1, 2, 1, 1, 2, 1, 3, 2, 1,
c    1              1, 2, 1, 3, 2, 1, 4, 3, 2, 1,
c    2              1, 2, 1, 3, 2, 1, 4, 3, 2, 1,
c    3              5, 4, 3, 2, 1/
c            z    / 1, 1, 1, 2, 1, 1, 2, 1, 2, 3,
c    1              1, 1, 2, 1, 2, 3, 1, 2, 3, 4,
c    2              1, 1, 2, 1, 2, 3, 1, 2, 3, 4,
c    3              1, 2, 3, 4, 5/
c
c     ----- order of components in shells in HONDO -----
c
c                    1      2      3      4
c                '   S','   X','   Y','   Z',
c                    5      6      7      8      9     10
c    1           '  XX','  YY','  ZZ','  XY','  XZ','  YZ',
c                   11     12     13     14     15  
c    2           ' XXX',' YYY',' ZZZ',' XXY',' XXZ',
c                   16     17     18     19     20
c    3           ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',
c                   21     22     23     24     25
c    4           'XXXX','YYYY','ZZZZ','XXXY','XXXZ',
c                   26     27     28     29     30
c    5           'YYYX','YYYZ','ZZZX','ZZZY','XXYY',
c                   31     32     33     34     35
c    6           'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/
c
c     ----- order of components in shells in NWCHEM -----
c
c                    1      2      3      4
c                '   S','   X','   Y','   Z',
c                    5      6      7      8      9     10
c    1           '  XX','  XY','  XZ','  YY','  YZ','  ZZ',
c                   11     12     13     14     15     
c    2           ' XXX',' XXY',' XXZ',' XYY',' XYZ',
c                   16     17     18     19     20
c    3           ' XZZ',' YYY',' YYZ',' YZZ',' ZZZ',
c                   21     22     23     24     25
c    4           'XXXX','XXXY','XXXZ','XXYY','XXYZ',
c                   26     27     28     29     30
c    5           'XXZZ','XYYY','XYYZ','XYZZ','XZZZ',
c                   31     32     33     34     35
c    6           'YYYY','YYYZ','YYZZ','YZZZ','ZZZZ'/
c
      parameter (mxatom=nw_max_atom)
      parameter (mxprim=2048)
      parameter (mxshel=512)
      parameter (mxnbf =2048)
      parameter (sqrt3 = 1.73205080756888d+00)
      parameter (sqrt5 = 2.23606797749979d+00)
      parameter (sqrt7 = 2.64575131106459d+00)
      parameter (sq5d3 = sqrt5/sqrt3)
      common/hnd_iofile/ir,iw,ip
      common/hnd_nshel/ex(mxprim),cs(mxprim),cp(mxprim),
     1                 cd(mxprim),cf(mxprim),cg(mxprim),
     2                 kstart(mxshel),katom(mxshel),ktype(mxshel),
     3                 kng(mxshel),kloc(mxshel),kmin(mxshel),
     4                 kmax(mxshel),nshell
      common/hnd_nwtohnd/inw_to_hnd(mxnbf)
      common/hnd_hndtonw/ihnd_to_nw(mxnbf)
      common/hnd_facntoh/fac_nwthnd(mxnbf)
c
#include "global.fh"
c
      logical dbug
      logical out 
c
      integer map_ntoh(35)
      double precision fac_ntoh(35)
c
c     ----- mapping of shell components from nwchem to hondo -----
c
      data map_ntoh /1,     2,     3,     4,
     1               5,     8,     9,     6,    10,    7,
     2              11,    14,    15,    16,    20,
     3              18,    12,    17,    19,    13,  
     4              21,    24,    25,    30,    33,
     5              31,    26,    34,    35,    28,
     6              22,    27,    32,    29,    23 /

      data fac_ntoh /1,     1,     1,     1,
     1               1, sqrt3, sqrt3,     1, sqrt3,    1,
     2               1, sqrt5, sqrt5, sqrt5, sqrt3,
     3           sqrt5,     1, sqrt5, sqrt5,     1,  
     4               1, sqrt7, sqrt7, sq5d3, sqrt3,
     5           sq5d3, sqrt7, sqrt3, sqrt3, sqrt7,
     6               1, sqrt7, sq5d3, sqrt7,     1/
c
      dbug=.false.
      out =.false.
c
      if(out.and.ga_nodeid().eq.0) then
         write(iw,*) ' in routine ... nw_hnd ... '
      endif
c
c     ----- loop over shells -----
c
     
      do ish = 1,nshell
         mini =  kmin(ish)
         maxi =  kmax(ish)
         loci =  kloc(ish)
         do im = mini,maxi

c     ----- pointer of component of shell in nwchem ----

            inw  = loci + im - mini

c     ----- pointer of component of shell in hondo -----

            ihnd = loci + map_ntoh(im) -mini

c     ----- switch from nwchem to hondo -----

            inw_to_hnd(inw ) = ihnd
            fac_nwthnd(inw)  = fac_ntoh(im)

c     ----- switch from hondo to nwchem -----

            ihnd_to_nw(ihnd) = inw 

         enddo
      enddo
c
      return
      end
