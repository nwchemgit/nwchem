C> \ingroup stpr_priv
C> @{
      subroutine stpr_select_grtdb(nat,rtdb,selbuf)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c::passed
      integer rtdb, nat
      integer selbuf(3,nat)
      external stpr_data ! For T3D linker
c
      if (.not.rtdb_get(rtdb,'stepper:select',MT_INT,(3*nat),selbuf))
     &    call errquit
     &    ('stpr_select_grtdb: rtdb_get failed',911, RTDB_ERR)
      end
      SUBROUTINE stpr_face (
     &     HESS, GRAD, COORD, CHARGE, SELECT, ATMASS,
     &     AMS, G, C, ETOT, TRUST, 
     &     ENERGY, CONVGE, CONVGG, CONVGGM, 
     &     NATD, NMODE,
     &     TRACK, BCKSTP, RPATH , rtdb, step_number)
c     
c     $Id$
c     
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
#include "util.fh"
      integer rtdb, him, kim, step_number, maxiter
      LOGICAL BCKSTP, TRACK, RPATH
      LOGICAL SECO, SWSEC, SWFIR
      LOGICAL SETM
      logical oselect, osetmass
c     CHARACTER LINE*7
      COMMON / STAPES / ITAP05, ITAP06, ITAP41,
     &     ITAP44, ITAP45
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      COMMON / SWITCH / SECO, SWSEC, SWFIR
      COMMON / THRESH / RCONVG,RSTD
      COMMON / MASS / TOTM, NUMAS
#include "mafdecls.fh"
#include "chow.fh"
      integer select
      logical oprint_info
      DIMENSION HESS(NAT3,NAT3), GRAD(NAT3), COORD(3,NATOM),
     &     CHARGE(NATOM), SELECT(NAT3), C(3,NATOM),
     &     ATMASS(NATOM), AMS(NUMAS), G(NATD*3)
C     
C     OMIN - minimization walk
C     TS - transition state walk
C     
C     STD - steepest descent walk
C     The step length (norm of the x vector) can not overpass RSTD
C     TROT and TRTR can enlarge the step length to the acceptable
C     RSTD limit
C     If they are false the STD step will be scaled down to RSTD
C     (if the step is too long) but it can not be scaled up
C     (if the step is shorter than RSTD)
C     TROT will allow scaling up in the region of a negative
C     NMODE eigenvalue
C     TRTR will allow scaling up along the whole walk
C     Using TRTR is dangerous since you can miss the minimum
C     
C     SECO - analytic second derivatives in the first step only
C     SWSEC - switch to second derivatives; then bomb
C     SWFIR - switch to first derivatives; then bomb
C     
C     Set default values.
C     
*     write(6,*)' nat3 = ',nat3
      oprint_info = util_print('information',print_low)
      oselect = .false.
      osetmass = .false.
      call stpr_set_defaults(OMIN,TS,STD,TROT,TRTR,MWC,SECO,SWFIR,
     &     SWSEC,RSTD,TRUSTC,TRUSTM,TRUST,NMODE,TRACK,RPATH,RCONVG,
     &     SETM,convggm,convgg,convge,maxiter)
      call stpr_load_rtdb(rtdb,omin,ts,std,trot,trtr,mwc,seco, swfir,
     &     swsec,track,rpath,setm,oselect,osetmass,nmode,rconvg,rstd,
     &     trustc, trustm, trust, convggm, convgg, convge, maxiter)
c     
      if (oselect) then
         call stpr_select_grtdb(natom,rtdb,select)
      endif
c     
      if(convggm.eq.0.0d0) CONVGGM=8.0D-4
      if(convgg .eq.0.0d0) CONVGG =1.0D-2
      if(convge .eq.0.0d0) CONVGE =1.0D-4
      AMUAU = (1.6605656D0/9.109534D0)*1.0D+04
      if (osetmass) then
         write(6,*)'  atomic number and masses read '
         call stpr_setmass_grtdb(natom,rtdb,him,kim,atmass)
         do 00003 i=1,natom
            if (atmass(i).eq.0.0d00) then
               if(int_mb((kim-1+i)).gt.numas) then
                  WRITE(6,*)'  Atomic number out of range. '
                  WRITE(6,*)'  Masses for atomic numbers .LE. ',
     &                 NUMAS,' only. '
                  WRITE(6,*)'  You must input the mass. '
                  call errquit('stepper:face fatal error',911,
     &       INPUT_ERR)
               endif
               atmass(i) = ams(int_mb((kim-1+i)))
            else
               atmass(i) = atmass(i)*amuau
            endif
00003    continue
         setm = .false.
         write(6,*)'  Atomic masses: '
         write(6,*) (atmass(i),i=1,natom)
         if(.not.MA_free_heap(him))
     &        call errquit
     &        ('stpr_face: failed to free dummy read array',911, MA_ERR)
         
      endif
      if (TS)  omin = .false.
      if (std) omin = .false.
c     
c     OMIN = .TRUE.
c     TS = .FALSE.
c     STD = .FALSE.
c     TROT = .FALSE.
c     TRTR = .FALSE.
c     MWC = .FALSE.
c     SECO = .FALSE.
c     SWFIR = .FALSE.
c     SWSEC = .FALSE.
c     RSTD = 5.D-02
c     TRUSTC = 0.1D0
c     TRUSTM = 5.D0
c     TRUST = 0.D0
c     NMODE = 0
c     TRACK = .FALSE.
c     RPATH = .FALSE.
c     RCONVG = 1.D-02
c     SETM = .TRUE.
c     
c     C
c     C     Read input.
c     C
c     CALL stpr_locate (ITAP05, '&STEPPER  ', IERROR)
c     IF (IERROR .EQ. 0) THEN
c     1   READ (ITAP05,'(A)',ERR=2,END=2) LINE
c     IF (LINE(1:4) .EQ. '&END') THEN
c     GO TO 2
c     ELSE IF (LINE(1:3) .EQ. 'MIN') THEN
c     OMIN = .TRUE.
c     WRITE(ITAP06,*)'  Keyword MIN found. '
c     WRITE(ITAP06,*)'  MIN = ',OMIN
c     ELSE IF (LINE(1:2) .EQ. 'TS') THEN
c     TS = .TRUE.
c     OMIN = .FALSE.
c     WRITE(ITAP06,*)'  Keyword TS found. '
c     WRITE(ITAP06,*)'  TS = ',TS
c     ELSE IF (LINE(1:3) .EQ. 'STD') THEN
c     STD = .TRUE.
c     OMIN = .FALSE.
c     WRITE(ITAP06,*)'  Keyword STD found. '
c     WRITE(ITAP06,*)'  STD = ',STD
c     ELSE IF (LINE(1:4) .EQ. 'RSTD') THEN
c     READ (ITAP05,*,END=2,ERR=2) RSTD
c     WRITE(ITAP06,*)'  Keyword RSTD found. '
c     WRITE(ITAP06,*)'  RSTD = ',RSTD
c     ELSE IF (LINE(1:4) .EQ. 'TROT') THEN
c     TROT = .TRUE.
c     WRITE(ITAP06,*)'  Keyword TROT found. '
c     WRITE(ITAP06,*)'  TROT = ',TROT
c     ELSE IF (LINE(1:4) .EQ. 'TRTR') THEN
c     TRTR = .TRUE.
c     WRITE(ITAP06,*)'  Keyword TRTR found. '
c     WRITE(ITAP06,*)'  TRTR = ',TRTR
c     ELSE IF (LINE(1:3) .EQ. 'MWC') THEN
c     MWC = .TRUE.
c     WRITE(ITAP06,*)'  Keyword MWC found. '
c     WRITE(ITAP06,*)'  MWC = ',MWC
c     ELSE IF (LINE(1:4) .EQ. 'SECO') THEN
c     SECO = .TRUE.
c     WRITE(ITAP06,*)'  Keyword SECO found. '
c     WRITE(ITAP06,*)'  SECO = ',SECO
c     ELSE IF (LINE(1:5) .EQ. 'SWSEC') THEN
c     SWSEC = .TRUE.
c     WRITE(ITAP06,*)'  Keyword SWSEC found. '
c     WRITE(ITAP06,*)'  SWSEC = ',SWSEC
c     ELSE IF (LINE(1:5) .EQ. 'SWFIR') THEN
c     SWFIR = .TRUE.
c     WRITE(ITAP06,*)'  Keyword SWFIR found. '
c     WRITE(ITAP06,*)'  SWFIR = ',SWFIR
c     ELSE IF (LINE(1:5) .EQ. 'TRUST') THEN
c     READ (ITAP05,*,END=2,ERR=2) TRUST
c     WRITE(ITAP06,*)'  Keyword TRUST found. '
c     WRITE(ITAP06,*)'  TRUST = ',TRUST
c     ELSE IF (LINE(1:5) .EQ. 'NMODE') THEN
c     READ (ITAP05,*,END=2,ERR=2) NMODE
c     WRITE(ITAP06,*)'  Keyword NMODE found. '
c     WRITE(ITAP06,*)'  NMODE = ',NMODE
c     ELSE IF (LINE(1:5) .EQ. 'TRACK') THEN
c     TRACK = .TRUE.
c     WRITE(ITAP06,*)'  Keyword TRACK found. '
c     WRITE(ITAP06,*)'  TRACK = ',TRACK
c     ELSE IF (LINE(1:6) .EQ. 'SELECT') THEN
c     READ (ITAP05,*,END=2,ERR=2) SELECT
c     WRITE(ITAP06,*)'  Keyword SELECT found. '
c     NMODE = 1
c     WRITE(ITAP06,*)'  SELECT: '
c     WRITE(6,1001)
c     DO 30 I = 1,NATOM
c     JSTART = 3*(I-1) + 1
c     JEND = JSTART + 2
c     WRITE(6,1002)I,(SELECT(J),J=JSTART,JEND)
c     30     CONTINUE
c     1001     FORMAT('  ATOM:',12x,'X',14x,'Y',14x,'Z')
c     1002     FORMAT('  ',I3,6X,3(1X,F14.10))
c     c          CALL stpr_matout(SELECT,3,NATOM,3)
c     ELSE IF (LINE(1:5) .EQ. 'RPATH') THEN
c     RPATH = .TRUE.
c     WRITE(ITAP06,*)'  Keyword RPATH found. '
c     ELSE IF (LINE(1:6) .EQ. 'RCONVG') THEN
c     READ (ITAP05,*,END=2,ERR=2) RCONVG
c     WRITE(ITAP06,*)'  Keyword RCONVG found. '
c     WRITE(ITAP06,*)'  RCONVG = ',RCONVG
c     ELSE IF (LINE(1:7) .EQ. 'SETMASS') THEN
c     WRITE(6,*)'  atomic number and masses read from',
c     &              ' standard in. '
c     DO 3 I = 1,NATOM
c     READ(ITAP05,*,ERR=2,END=2) IREAD,ATMASS(I)
c     write(6,*)'  atom ',i,' mass ',atmass(i)
c     IF(ATMASS(I).EQ.0.0D00) THEN
c     IF(IREAD.GT.NUMAS)THEN
c     WRITE(6,*)'  Atomic number out of range. '
c     WRITE(6,*)'  Masses for atomic numbers .LE. ',
c     &                            NUMAS,' only. '
c     WRITE(6,*)'  You must input the mass. '
c     call errquit('stepper:face fatal error',911)
c     ENDIF
c     ATMASS(I) = AMS(IREAD)
c     ELSE
c     ATMASS(I)=ATMASS(I)*AMUAU
c     ENDIF
c     3     CONTINUE
c     SETM = .FALSE.
c     write(6,*)'  Atomic masses: '
c     write(6,*) (atmass(i),i=1,natom)
c     END IF
c     GO TO 1
c     ELSE
c     c       WRITE (ITAP06,*)'  Label &STEPPER not found. '
c     c       WRITE (ITAP06,*)'  Default input used. '
c     ICNVG = -99
c     END IF
c     2 CONTINUE
      IF (MWC.AND.SETM) THEN
C     
C     Get charges
C     
c     DO 4 I = 1, NATOM
c     
c     no current need for mass weighted hess ... temp fix - jan.
c     
c     CHARGE(I) = ZAN(I)
c     4     CONTINUE
         WRITE(6,'(/A/)')
     +        ' Default masses from atomic number (charge) used. '
         DO 5 I = 1,NATOM
            IDUM = INT(CHARGE(I))
            ATMASS(I) = AMS(IDUM)
            WRITE(6,*)'  Atom ',I,' Mass ',ATMASS(I)
 5       CONTINUE
      ENDIF
      IF (MWC) THEN
C     
C     Calculate total mass.
C     
         TOTM = 0.D0
         DO 6 I = 1,NATOM
            TOTM = TOTM + ATMASS(I)
 6       CONTINUE
         WRITE(6,*)'  Total mass is ',TOTM,' a.u.'
      ENDIF
      IF (TRUST.EQ.0.D0) THEN
         IF (MWC) THEN
            TRUST=TRUSTM
         ELSE
            TRUST=TRUSTC
         ENDIF
      ENDIF
      if (oprint_info) then
         IF(TS)THEN
            WRITE(ITAP06,*)'  Step sought for transition state. '
         ELSE IF (OMIN) THEN
            WRITE(ITAP06,*)'  Step sought for minimum. '
         ELSE IF (STD) THEN
            WRITE(ITAP06,*)'  Finding steepest descent step. '
         ENDIF
      endif
      IF ((NMODE.NE.0).AND..NOT.(TS.OR.STD)) THEN
         WRITE(ITAP06,*)'  keyword TS or STD required. '
         call errquit('stepper:face fatal error',911, INPUT_ERR)
      ENDIF
      IF ((TS.OR.STD).AND.(NMODE.EQ.0)) THEN
         NMODE = 1
      ENDIF
      IF (STD) THEN
         TRACK = .TRUE.
      ENDIF
      IF (RPATH.AND.OMIN) THEN
         RPATH = .FALSE.
         WRITE(ITAP06,*)'  Keyword RPATH ignored. '
      ENDIF
C     
C     Get constants from SUPERMOLECULE or set here.
C     
C     
C     IWCTR = step counter.
C     
      IWCTR = IWCTR + 1
      step_number = iwctr
      if (step_number.ge.maxiter) step_number = -step_number
C     
C     ICALC = order of derivatives.
C     
      ICALC = 1
      IF(RPATH.AND.(ICALC.NE.2)) THEN
         RPATH = .FALSE.
         WRITE(ITAP06,*)'  Keyword RPATH ignored. '
      ENDIF
      IF (SECO.AND.(IWCTR.GT.1)) THEN
         SECO = .FALSE.
         WRITE(6,*)'  KEYWORD SECO IGNORED'
         IF ((ICALC.EQ.1).AND.SECO) THEN
            SECO = .FALSE.
            WRITE(6,*)'  KEYWORD SECO IGNORED'
         ENDIF
      ENDIF
      IF ((ICALC.EQ.1).AND.SWFIR) THEN
         SWFIR = .FALSE.
         WRITE(6,*)'  KEYWORD SWFIR IGNORED'
      ENDIF
      IF ((ICALC.EQ.2).AND.SWSEC) THEN
         SWSEC = .FALSE.
         WRITE(6,*)'  KEYWORD SWSEC IGNORED'
      ENDIF
      IF((SECO.AND.SWSEC).OR.(SECO.AND.SWFIR).OR.
     &     (SWSEC.AND.SWFIR)) THEN
         WRITE(6,*)'  TOO MANY SWITCHING KEYWORDS'
         call errquit('stepper:face fatal error',911, INPUT_ERR)
      ENDIF
C     
C     Get energy.
C     
      ENERGY = ETOT
      CALL stpr_wlkchk(ETOT,CONVGE,BCKSTP,NMODE)
      IF(BCKSTP)RETURN
C     
C     Get coordinates.
C     
      DO 8 I = 1, 3
         DO 7 N = 1, NATOM
            COORD(I,N) = C(I,N)
 7       CONTINUE
    8 CONTINUE
C     
C     Get gradients and if calculated, the analytic hessian.
C     
      DO 9 I = 1, NAT3
         GRAD(I) = G(I)
    9 CONTINUE
*      grad_max = dabsmax(nat3,grad)
*      write(6,*)' grad max in stpr_face ',grad_max
*      write(6,*)' grad ',(grad(i), i=1,nat3)
      IF (ICALC .EQ. 2) THEN
         WRITE(6,*)'  Second derivatives NOT available. '
         call errquit('stepper:face fatal error',911, INPUT_ERR)
C     
C     Check to make sure it is symmetric
C     
         CALL stpr_chkhess(HESS)
C     
C     Construct triangular hessian.
C     
         CALL stpr_mktri (HESS, HESS)
      ELSEIF (ICALC .EQ. 1) THEN
         if (oprint_info) then
            WRITE(6,*)'  Update Hessian used. '
         endif
      ELSE
         WRITE(6,*)'  Unknown calculation type. Program aborts. '
         call errquit('stepper:face fatal error',911, INPUT_ERR)
      END IF
      RETURN
      END
      subroutine stpr_set_defaults(OMIN,TS,STD,TROT,TRTR,MWC,SECO,
     &    SWFIR,SWSEC,RSTD,TRUSTC,TRUSTM,TRUST,NMODE,TRACK,RPATH,
     &    RCONVG,SETM,CONVGGM,CONVGG,CONVGE,maxiter)
c
      implicit none
      logical omin, ts, std, trot, trtr, mwc, seco, swfir, swsec
      logical track, rpath, setm
      integer nmode, maxiter
      double precision rconvg, rstd, trustc, trustm, trust
      double precision convggm, convgg, convge
c
      OMIN = .TRUE.
      TS = .FALSE.
      STD = .FALSE.
      TROT = .FALSE.
      TRTR = .FALSE.
      MWC = .FALSE.
      SECO = .FALSE.
      SWFIR = .FALSE.
      SWSEC = .FALSE.
      TRACK = .FALSE.
      RPATH = .FALSE.
      SETM = .TRUE.
c
      NMODE = 0
      maxiter = 20
c
      RCONVG = 1.0D-02
      RSTD   = 5.0D-02
      TRUSTC = 0.1D+00
      TRUSTM = 5.0D+00
      TRUST  = 0.0D+00
      convggm= 0.0d+00
      convgg = 0.0d+00
      convge = 0.0d+00
c
      end
      subroutine stpr_select_srtdb(nat,rtdb,selbuf)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c::passed
      integer rtdb, nat
      integer selbuf(3,nat)
c
      if (.not.rtdb_put(rtdb,'stepper:select',MT_INT,(3*nat),selbuf))
     &    call errquit
     &    ('stpr_select_srtdb: rtdb_put failed',911, RTDB_ERR)
      end
      subroutine stpr_setmass_srtdb(nat,rtdb,imass,xmass)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c::passed
      integer rtdb, nat
      integer imass(nat)
      double precision xmass(nat)
c
      if (.not.rtdb_put(rtdb,'stepper:setmass:int',MT_INT,nat,imass))
     &    call errquit
     &    ('stpr_setmass_srtdb: rtdb_put failed',911, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'stepper:setmass:dbl',MT_DBL,nat,xmass))
     &    call errquit
     &    ('stpr_setmass_srtdb: rtdb_put failed',911, RTDB_ERR)
      end
      subroutine stpr_setmass_grtdb(nat,rtdb,himass,kimass,xmass)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c::passed
      integer rtdb, nat
      integer himass, kimass
      double precision xmass(nat)
c::local
      integer my_ma_type
      integer my_nat
c
      if (.not.rtdb_ma_get
     &    (rtdb,'stepper:setmass:int',my_ma_type,my_nat,himass))
     &    call errquit
     &    ('stpr_setmass_grtdb: rtdb_ma_get failed',911, RTDB_ERR)
      if(.not.((my_ma_type.eq.MT_INT).and.(my_nat.eq.nat))) then
        write(6,*)' my_ma_type .ne. MT_INT',my_ma_type,MT_INT
        write(6,*)' my_nat .ne. nat       ',my_nat, nat
        call errquit('stpr_setmass_grtdb: rtdb_ma_get failure',911,
     &       RTDB_ERR)
      endif
      if (.not.MA_get_index(himass,kimass))
     &    call errquit('stpr_setmass_grtdb: MA_get_index error',911,
     &       MA_ERR)
      if (.not.rtdb_get(rtdb,'stepper:setmass:dbl',MT_DBL,nat,xmass))
     &    call errquit
     &    ('stpr_setmass_grtdb: rtdb_get failed',911, RTDB_ERR)
      end
      subroutine stpr_store_rtdb(rtdb,omin,ts,std,trot,trtr,mwc,seco,
     &    swfir,
     &    swsec,track,rpath,setm,oselect,osetmass,nmode,rconvg, rstd,
     &    trustc, trustm, trust, convggm, convgg, convge, maxiter)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stpr_inP.fh"
c::passed
      integer rtdb
      logical omin,  ts,   std,     trot,  trtr
      logical mwc,   seco, swfir,   swsec, track
      logical rpath, setm, oselect, osetmass
      integer nmode, maxiter
      double precision rconvg, rstd, trustc, trustm, trust
      double precision convggm, convgg, convge
c::local
      integer hlog,klog,    hdbl,kdbl,  hint, kint
c
c allocate memory
c.. logs array 
      if(.not.MA_Push_Get(MT_LOG,ns_log,'stepper logs',hlog,klog))
     &    call errquit
     &    ('stpr_store_rtdb: push for logs failed?',911, MA_ERR)
c.. one one int no allocate now
      if(.not.MA_Push_get(MT_INT,ns_int,'stepper ints',hint,kint))
     &    call errquit
     &    ('stpr_store_rtdb: push for ints failed?',911, MA_ERR)
      if(.not.MA_Push_get(MT_DBL,ns_dbl,'stepper dbls',hdbl,kdbl))
     &    call errquit
     &    ('stpr_store_rtdb: push for dbls failed?',911, MA_ERR)
c
      log_mb( klog)       = omin ! OVERRIDDEN BELOW (see stpr_load_rtdb)
      log_mb((klog +  1)) = ts   ! OVERRIDDEN BELOW (see stpr_load_rtdb)
      log_mb((klog +  2)) = std
      log_mb((klog +  3)) = trot
      log_mb((klog +  4)) = trtr
      log_mb((klog +  5)) = mwc
      log_mb((klog +  6)) = seco
      log_mb((klog +  7)) = swfir
      log_mb((klog +  8)) = swsec
      log_mb((klog +  9)) = track
      log_mb((klog + 10)) = rpath
      log_mb((klog + 11)) = setm
      log_mb((klog + 12)) = oselect
      log_mb((klog + 13)) = osetmass
c
      int_mb(kint)    = nmode
      int_mb(kint+1)   = maxiter
c
      dbl_mb( kdbl) = rconvg
      dbl_mb((kdbl +  1)) = rstd
      dbl_mb((kdbl +  2)) = trustc
      dbl_mb((kdbl +  3)) = trustm
      dbl_mb((kdbl +  4)) = trust
      dbl_mb((kdbl +  5)) = convgg
      dbl_mb((kdbl +  6)) = convge
      dbl_mb((kdbl +  7)) = convggm
c
c     The values of stepper:ts and stepper:omin are also overriden by
c     task_optimize and task_saddle.
c
      if (.not.rtdb_put(rtdb,'stepper:ts',MT_LOG,1,ts))
     &    call errquit('stpr_store_rtdb: rtdb_put failed',911, RTDB_ERR)
c
      if (.not.rtdb_put(rtdb,'stepper:omin',MT_LOG,1,omin))
     &    call errquit('stpr_store_rtdb: rtdb_put failed',911, RTDB_ERR)
c
      if (.not.rtdb_put
     &    (rtdb,'stepper:logicalvals',MT_LOG,ns_log,log_mb(klog)))
     &    call errquit
     &    ('stpr_store_rtdb: rtdb_put failed',911, RTDB_ERR)
c
      if (.not.rtdb_put
     &    (rtdb,'stepper:integervals',MT_INT,ns_int,int_mb(kint)))
     &    call errquit
     &    ('stpr_store_rtdb: rtdb_put failed',911, RTDB_ERR)
c
      if (.not.rtdb_put
     &    (rtdb,'stepper:doublevals',MT_DBL,ns_dbl,dbl_mb(kdbl)))
     &    call errquit
     &    ('stpr_store_rtdb: rtdb_put failed',911, RTDB_ERR)
c
      if(.not.MA_Pop_Stack(hdbl))
     &    call errquit
     &    ('stpr_store_rtdb: error freeing hdbl',911, RTDB_ERR)
c
      if(.not.MA_Pop_Stack(hint))
     &    call errquit
     &    ('stpr_store_rtdb: error freeing hint',911, MA_ERR)
c
      if(.not.MA_Pop_Stack(hlog))
     &    call errquit
     &    ('stpr_store_rtdb: error freeing hlog',911, MA_ERR)
c
*      if(rtdb_print(rtdb,.true.))write(6,*)'store'
c
      end
      subroutine stpr_load_rtdb(rtdb,omin,ts,std,trot,trtr,mwc,seco,
     &    swfir,
     &    swsec,track,rpath,setm,oselect,osetmass,nmode,rconvg, rstd,
     &    trustc, trustm, trust, convggm, convgg, convge, maxiter)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stpr_inP.fh"
c::passed
      integer rtdb
      logical omin,  ts,   std,     trot,  trtr
      logical mwc,   seco, swfir,   swsec, track
      logical rpath, setm, oselect, osetmass
      integer nmode, maxiter
      double precision rconvg, rstd, trustc, trustm, trust
      double precision convggm, convgg, convge
c::local
      logical status
      integer hlog,klog,    hdbl,kdbl
      integer hint, kint
      integer my_ma_type, my_len
c
*      if(rtdb_print(rtdb,.true.))write(6,*)'load 1'
      status = rtdb_ma_get
     &    (rtdb,'stepper:logicalvals',my_ma_type,my_len,hlog)
c
c.. status=.false. ==> either no input or rtdb is corrupt
      if(.not.status) goto 00001
c
      if(.not.((my_ma_type.eq.MT_LOG).and.(my_len.eq.ns_log))) then
        write(6,*)' read type .ne. expected type ',my_ma_type,MT_LOG
        write(6,*)' read length .ne. expected    ',my_len,ns_log
        call errquit('stpr_load_rtdb: rtdb_ma_get fail',911,
     &       RTDB_ERR)
      endif
      if(.not.MA_get_index(hlog,klog))
     &    call errquit('stpr_load_rtdb: ma_get_index fail',911, MA_ERR)
*      if(rtdb_print(rtdb,.true.))write(6,*)'load 2'
c
c... currently only two integer input nmode/maxiter
      if(.not.rtdb_ma_get
     &    (rtdb,'stepper:integervals',my_ma_type,my_len,hint))
     &    call errquit
     &    ('stpr_load_rtdb:rtdb_ma_get failed',911, RTDB_ERR)
      if (.not.MA_get_index(hint,kint))
     &    call errquit('stpr_load_rtdb: ma_get_index fail',911, MA_ERR)
c
      if(.not.((my_ma_type.eq.MT_INT).and.(my_len.eq.ns_int)))then
        write(6,*)' read type .ne. expected type ',my_ma_type,MT_INT
        write(6,*)' read length .ne. expected    ',my_len,ns_int
        call errquit('stpr_load_rtdb: rtdb_ma_get fail',911, RTDB_ERR)
      endif
c
*      if(rtdb_print(rtdb,.true.))write(6,*)'load 3'
      if (.not.rtdb_ma_get
     &    (rtdb,'stepper:doublevals',my_ma_type,my_len,hdbl))
     &    call errquit
     &    ('stpr_load_rtdb: rtdb_ma_get failed',911, RTDB_ERR)
c
      if(.not.((my_ma_type.eq.MT_DBL).and.(my_len.eq.ns_dbl)))then
        write(6,*)' read type .ne. expected type ',my_ma_type,MT_DBL
        write(6,*)' read length .ne. expected    ',my_len,ns_dbl
        call errquit('stpr_load_rtdb: rtdb_ma_get fail',911, RTDB_ERR)
      endif
      if(.not.MA_get_index(hdbl,kdbl))
     &    call errquit('stpr_load_rtdb: ma_get_index fail',911,
     &       RTDB_ERR)
c
c
      if (.not.rtdb_get(rtdb,'stepper:ts',MT_LOG,1,ts))
     &    call errquit('stpr_load_rtdb: rtdb_put failed',911,
     &       RTDB_ERR)
c
      if (.not.rtdb_get(rtdb,'stepper:omin',MT_LOG,1,omin))
     &    call errquit('stpr_load_rtdb: rtdb_put failed',911,
     &       RTDB_ERR)
c
***   omin     = log_mb( klog)          ! NOW DONE INDIVIDUALLY
***   ts       = log_mb((klog +  1))    ! NOW DONE INDIVIDUALLY
      std      = log_mb((klog +  2)) 
      trot     = log_mb((klog +  3)) 
      trtr     = log_mb((klog +  4)) 
      mwc      = log_mb((klog +  5)) 
      seco     = log_mb((klog +  6)) 
      swfir    = log_mb((klog +  7)) 
      swsec    = log_mb((klog +  8)) 
      track    = log_mb((klog +  9)) 
      rpath    = log_mb((klog + 10)) 
      setm     = log_mb((klog + 11)) 
      oselect  = log_mb((klog + 12)) 
      osetmass = log_mb((klog + 13)) 
c
      nmode   = int_mb(kint)
      maxiter = int_mb(kint+1)
c
      rconvg = dbl_mb( kdbl) 
      rstd   = dbl_mb((kdbl +  1)) 
      trustc = dbl_mb((kdbl +  2)) 
      trustm = dbl_mb((kdbl +  3)) 
      trust  = dbl_mb((kdbl +  4)) 
      convgg = dbl_mb((kdbl +  5)) 
      convge = dbl_mb((kdbl +  6)) 
      convggm= dbl_mb((kdbl +  7)) 
c
c.... free heap memory
        if(.not.MA_free_heap(hdbl))
     &      call errquit
     &      ('stpr_load_rtdb: error freeing hdbl',911, MA_ERR)
c
        if(.not.MA_free_heap(hint))
     &      call errquit
     &      ('stpr_load_rtdb: error freeing hint',911, MA_ERR)
c
        if(.not.MA_free_heap(hlog))
     &      call errquit
     &      ('stpr_load_rtdb: error freeing hlog',911, MA_ERR)
c
00001 continue
c
      call util_flush(6)
c
      end
      block data stpr_data
      implicit none
      COMMON / CFACE / IWCTR,NATOM,ICALC
      integer iwctr, natom, icalc
c
      data iwctr /0/
      data natom /0/
      data icalc /0/
c
      end
C> @}
