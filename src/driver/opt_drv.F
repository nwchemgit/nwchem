      logical function drv_opt(rtdb)
C$Id: opt_drv.F,v 1.71 1999-06-10 02:05:06 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer  rtdb
      integer  idriver
      logical  drv_optx
      external drv_optx
c
      call util_print_push
      call util_print_rtdb_load(rtdb, 'driver')
      call ecce_print_module_entry('driver')
c
      idriver=0
c
      drv_opt = drv_optx(rtdb,idriver)
c
      if (drv_opt) then
         call ecce_print_module_exit('driver', 'ok')
      else
         call ecce_print_module_exit('driver', 'failed')
      endif
c
      call movecs_ecce_print_on() ! Restore MO printing
      call util_print_pop
c
      end
      logical function drv_optx(rtdb,idriver)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "coptopt.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb, idriver
c     
      integer geom
      integer istep
      logical converged
      double precision start    ! Tracks time used in last step
      integer required          ! Time required
      logical  driver_converged, task_gradient
      external driver_converged, task_gradient
c
c     Read input, load /coptopt/, get geometry
c
      call driver_initialize(rtdb, geom)
c     
c     ----- initialize hessian matrix -----
c     
      call driver_hss_init(rtdb,geom)
c
c     Iterative loop
c
      start = util_wallsec()
      do istep = 1, nptopt
         if (oprint) write(6,1) istep
 1       format(/,10x,'--------',/,10x,'Step',i4,/,10x,'--------')
         if ((ga_nodeid() .eq. 0) .and.
     $        util_print('geometry',print_default)) then
            if (.not. geom_print(geom)) call errquit('driver: geom?',0)
         endif
c     
c     Save old energy and gradient, and compute new energy and gradient
c
         penergy = energy
         call dcopy(nvar, g, 1, gp, 1)
c
         if (.not. task_gradient(rtdb))
     $        call errquit('driver: task_gradient failed',0)
         if (.not. rtdb_get(rtdb, 'task:gradient', mt_dbl, ncart,
     $        gx)) call errquit('driver: could not get gradient',0)
         if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $        call errquit('driver: could not get energy',0)
c
c     Disable printing to ecce of movecs after first point
c
      call movecs_ecce_print_off()
c     
c     If necessary transform gradient to internal coordinates
c     
         if(zcoord) then
            call geom_bandbi(geom)
            call driver_tfgx(geom) ! g() will be set to constrained gq()
         else
            call driver_cart_pmat(geom)
            call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
         endif
c     
c     Update the Hessian
c     
         if (istep .ge. 2) call driver_hessian_update()
c
c     Compute step/gradient info and print for user
c     (for current energy & gradient, and the previous alpha*step).
c
         call driver_compute_info()
         call driver_print(geom,istep)
c     
c     Check for convergence and time remaining
c     
         if (driver_converged()) then
            converged = .true.
            goto 100
         endif
         required = int(1.2d0*(util_wallsec() - start)) + 1
         if (.not. util_test_time_remaining(rtdb,required)) goto 200
         start = util_wallsec()
c
c     Generate a new search direction
c
         call driver_search_dir(geom) ! Fills in ds(), recommends alpha
c     
c     Either do line search or take recommended step.
c     
         if (linopt .gt. 0) then
            call driver_line_search(rtdb, geom) ! Updates geom,alpha,ds,energy
         else
            call driver_take_step(rtdb, geom) ! Updates geom
         endif
c
      enddo                     ! End of iterative loop
 200  write(6,201)
 201  format(/,1x,63('-')/,1x,'Failed to converge in maximum number',
     $     ' of steps or available time'/,1x,63('-')/)
      converged = .false.
c
 100  if (converged) then
         write(6,101) 
 101     format(/,6x,22('-'),/,6x,'Optimization converged',/,
     $        6x,22('-'),/)
      endif
c
      if (oprint) then
         if (zcoord) call util_print_enable('geomzmat')
         if (.not. geom_print(geom)) call errquit
     $        ('hnd_opt_drv: geom_print?',0)
         if (util_print('bonds',print_default)) then
            if (.not.geom_print_distances(geom)) call errquit(
     &           'hnd_opt_drv: geom_print_distances failed',911)
         endif
         if (util_print('angles',print_default)) then
            if (.not.geom_print_angles(geom)) call errquit(
     &           'hnd_opt_drv: geom_print_angles failed',911)
         endif
      endif
c
      if (.not.geom_destroy(geom))
     &     call errquit('hnd_opt: geom_destroy?', 911)
c
      call ga_sync()
c
      drv_optx=converged
c     
      end
      subroutine driver_tfgx(geom)
      implicit none
#include "coptopt.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
      integer geom
c
c     Transform the cartesian gradient in gx() into (redundant)
c     internal coordinates in gq().  Apply constraints etc.
c     putting result into g().
c     
c     Assumes that bandbi has already been called.
c
c     Only gq() and g() are modified.
c
c     gq() = gradient in internals with no constants/constraints applied
c     g()  = gradient in internals with    constants/constraints applied
c
      integer l_bi, k_bi, i, j, num
      double precision sum
      character*8 test
      logical geom_zmt_get_varinfo
c
      if (.not. ma_push_get(mt_dbl, ncart*nzvar,'binv',l_bi, k_bi))
     $     call errquit('driver: memory for binv', ncart*nzvar)
      call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nzvar)
      if (odebug) then
         write(6,*) ' B^-1 '
         call output(dbl_mb(k_bi), 1, ncart, 1, nvar, nvar, nvar, 1)
         write(6,*) ' Cartesian gradient'
         call doutput(gx, 1, ncart, 1, 1, ncart, 1, 1)
      endif
      call dgemv('t', ncart, nzvar, 1.0d0, dbl_mb(k_bi), ncart,
     $     gx, 1, 0.0d0, gq, 1)
      if (.not. ma_pop_stack(l_bi)) call errquit('driver: ma?',99)
      if (odebug) then
         write(6,*) ' Internal coordinate gradient '
         call doutput(gq, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     Apply constants and constraints to the gradient before
c     it is used for upating the Hessian.  
c
c     Currently constants and constraints only happen in the
c     user input non-redundant coordinates ... when we
c     eventually do this for the redundant coordinates
c     then it will have to transform the gradient to & from
c     the redundant set.  DON'T DO THIS YET.
c
      if (nzfrz .gt. 0) then
         do i = 1, nzfrz
            gq(izfrz(i)) = 0d0
         enddo
      endif
c
      if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $     call errquit('hnd_opt_tfgx: varinfo?',0)
      do i = 1, nzvar
         if (zvarname(i).ne.' ') then
            sum = 0d0
            num = 0
            test = zvarname(i)
            do j = 1, nzvar
               if (zvarname(j).eq.test) then 
                  sum = sum + g(j)*zvarsign(j)
                  num = num + 1
               endif
            enddo
            if (num .gt. 1) then
               sum = sum / dble(num)
               do j = 1, nzvar
                  if (zvarname(j).eq.test) then 
*     write(6,77) j, test, zvarsign(j), g(j),sum
*     77                  format(' Constraining : ',i3,2x,a,2x,f3.1,2f10.6)
                     g(j) = sum*zvarsign(j)
                     zvarname(j) = ' '
                  endif
               enddo
            endif
         endif
      enddo
c
c     Just stomped on zvarname ... restore it
c
      if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $     call errquit('hnd_opt_tfgx: varinfo?',0)
c
c     Finally copy gq() into g()
c
      call dcopy(nzvar, gq, 1, g, 1)
c
      end
      subroutine hnd_opt_hss_init(rtdb,geom)
      call errquit('hnd_opt_hss_init: not yet ported to new code!',0)
      end
      subroutine driver_hss_init(rtdb,geom)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "coptopt.fh"
c
      integer rtdb, geom

c
      logical geom_zmt_get_izmat, geom_zmt_get_nizmat
      external geom_zmt_get_izmat, geom_zmt_get_nizmat
c
      double precision zero
      parameter (zero=0.0d+00)
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)     
      parameter (mxcoor=1500)
c
c     These commons are used in the internal coordinate guess
c
      integer nuc
      COMMON/HND_MOLNUC/NUC(MXATOM)
      double precision c, zan
      integer natom
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),natom
      integer nnzmat, nnzvar, nnvar
      common/hnd_zmtpar/nnzmat,nnzvar,nnvar
c
      integer l_hess, k_hess, l_zmat, k_zmat, i
      logical old_hessian
      character*16 atom_tags(mxatom)
c
      nnzmat = nzmat
      nnzvar = nzvar
      nnvar  = nzvar
      if (.not. geom_ncent(geom,natom))
     1       call errquit('hnd_opt: geom_ncent?',911)
c
      old_hessian=.false.
      if(inhess.eq.0) call driver_opt_check_hess(nvar, old_hessian)
      if (old_hessian) then
         if (oprint) write(6,*) ' Using old Hessian'
         return
      else
         if (oprint) write(6,*) ' Using diagonal Hessian '
      endif
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_hess, k_hess)) call errquit
     $     ('driver_init_hess: failed allocating hessian',nvar**2)
c
      if (.not. zcoord) then
c
c     Cartesians are easy
c
         call dfill(nvar**2, 0.0d0, dbl_mb(k_hess), 1)
         call dfill(nvar, 0.5d0, dbl_mb(k_hess), nvar+1)
      else
c
c     Internals.
c
         call geom_bandbi(geom) ! SHOULD USE COMPUTE_ZMATRIX
         if (.not. ma_push_get(mt_dbl, nzvar, 'zmat',
     $        l_zmat, k_zmat)) call errquit
     $        ('driver_init_hess: failed allocating zmat ', nzvar)
         call geom_hnd_get_data('zmat',dbl_mb(k_zmat),nzvar)
c
         if (.not. ma_verify_allocator_stuff()) call errquit('h',0)
c
C     ----- INTERNALLY GENERATED (DIAGONAL) INITIAL HESSIAN -----
C
c     this routine uses /molnuc/, /zmtpar/, /molxyz/ which 
c     must be initialized before calling
c
         if (.not. geom_cart_get(geom, nat, atom_tags, 
     $        c, zan)) call errquit('driver_init_hess: geom?',geom)
         do i = 1, nat
            nuc(i) = nint(zan(i)) ! Won't work if an ECP is loaded
         enddo
         nvar = nzvar
         call hnd_opt_hss_init_zmt(dbl_mb(k_hess),
     $        nzvar, dbl_mb(k_zmat), izmat, .false.)
         if (odebug) then
            write(6,*) ' The initial hessian '
            call output(dbl_mb(k_hess),1,nzvar,1,nzvar,nzvar,nzvar,1)
         endif
      endif
c
      call geom_hnd_put_data('drv.hess', dbl_mb(k_hess), nvar*nvar)
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('driver_init_hess ma corrupt',0)
c
      END
      SUBROUTINE HND_OPT_HSS_INIT_ZMT(HSSQ,NINT,ZMAT,IZMAT,SOME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (ZERO=0.0D+00,      ONE=1.0D+00,
     1     QUARTR=0.25D+00,   PT5=0.5D+00,
     2     THIRD=ONE/3.0D+00, TOHB2=6.423D-02,
     3     TOANGS=0.52917724924D+00)
      PARAMETER (MXATOM=500)
      LOGICAL OUT 
      LOGICAL SOME
      LOGICAL DBUG
      LOGICAL BADZ
      CHARACTER*8 HGUESS
      CHARACTER*8 GUESS
      CHARACTER*8 GUESS1
      CHARACTER*8 GUESS2
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION HSSQ(NINT,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION IROW(86),DPARM(6,6),APARM(6,6),RCOV(86)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-HSSGSS-'/
C     
C     --- DATA TABLES FOR BADGER'S RULES (SEE J.C.P. 2, 128(1934)) ---
C     A GENERALIZED BADGER'S RULE, AND THESE PARAMETERS ARE FROM
C     D.R.HERSCHBACH, V.W.LAURIE, J.CHEM.PHYS. 35, 458-463(1961).
C     
      DATA IROW  /2*1,8*2,8*3,18*4,18*5,32*6/
      DATA APARM /1.26D+00,1.66D+00,1.84D+00,1.98D+00,2.03D+00,2.03D+00,
     1     1.66D+00,1.91D+00,2.28D+00,2.35D+00,2.33D+00,2.50D+00,
     2     1.84D+00,2.28D+00,2.41D+00,2.52D+00,2.61D+00,2.60D+00,
     3     1.98D+00,2.35D+00,2.52D+00,2.58D+00,2.66D+00,2.75D+00,
     4     2.03D+00,2.33D+00,2.61D+00,2.66D+00,2.85D+00,2.76D+00,
     5     2.03D+00,2.50D+00,2.60D+00,2.75D+00,2.75D+00,3.00D+00/
      DATA DPARM /2.50D-02,0.30D+00,0.38D+00,0.49D+00,0.51D+00,0.81D+00,
     1     0.30D+00,0.68D+00,0.74D+00,0.85D+00,0.68D+00,0.97D+00,
     2     0.38D+00,0.74D+00,1.18D+00,1.02D+00,1.28D+00,0.84D+00,
     3     0.49D+00,0.85D+00,1.02D+00,1.41D+00,0.86D+00,1.14D+00,
     4     0.51D+00,0.68D+00,1.28D+00,0.86D+00,1.62D+00,1.25D+00,
     5     0.81D+00,0.97D+00,0.84D+00,1.14D+00,1.25D+00,1.83D+00/
C     
C     COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
C     EXCEPT VAN DER WAALS RADII FOR HE,NE,AR,KR (SAME SOURCE),
C     AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT,RN
C     
      DATA (RCOV(NUCZ),NUCZ=1,2)/0.30D+00,1.22D+00/
      DATA (RCOV(NUCZ),NUCZ=3,10)
     1     /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     2     0.70D+00,0.66D+00,0.58D+00,1.60D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     1     /1.66D+00,1.36D+00,1.25D+00,1.17D+00,
     2     1.10D+00,1.04D+00,0.99D+00,1.91D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     1     /2.03D+00,1.74D+00,
     2     1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     3     1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     4     1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.98D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     1     /2.22D+00,1.92D+00,
     2     1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     3     1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     4     1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,2.09D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     1     /2.35D+00,1.98D+00,
     2     1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D0,
     3     1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D0,
     4     1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     5     1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     6     1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,2.20D+00/
C     
      DATA GUESS,GUESS1,GUESS2 /'GUESS   ','GUESS1  ','GUESS2  '/
C     
      DIST(IAT,JAT)=SQRT((C(1,IAT)-C(1,JAT))**2+
     1     (C(2,IAT)-C(2,JAT))**2+
     1     (C(3,IAT)-C(3,JAT))**2) * TOANGS
      BADZ(NUCZ) = NUCZ.LE.0.OR.NUCZ.GT.86
C     
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      IF(OUT) THEN
         WRITE(IW,9168)
         WRITE(IW,9188) NAT,NVAR,NZVAR,NINT,NZMAT
      ENDIF
C     
C     ----- THIS ROUTINE MAKES A DIAGONAL GUESS OF THE HESSIAN -----
C     IN INTERNAL COORDINATE SPACE.
C     
      DO J=1,NZVAR
         DO I=1,NZVAR
            HSSQ(I,J)=ZERO
         ENDDO
      ENDDO
C     
C     -----          INTERNAL COORDINATES                      -----
C     
C     FOR MORE INFORMATION ON -GUESS1- OR -GUESS2- SEE
C     H.B.SCHLEGEL, THEORET.CHIM.ACTA, 66, 333-340 (1984).
C     T.H.FISCHER, J.ALMLOF, J.PHYS.CHEM. 96, 9786-9774 (1992).
C     -GUESS1- AND -GUESS2- ARE UNDOCUMENTED WAYS TO GET PURE
C     BADGER'S RULES, OR PURE FISCHER/ALMLOF RULES.
C     THE DEFAULT -GUESS- USES F/A FOR STRETCHES AND BENDS,
C     BUT NOT FOR TORSIONS OR OUT-OF-PLANE BENDS.
C     
      IF(DBUG) THEN
         WRITE(IW,9208) (IZMAT(IZ),IZ=1,NZMAT)
         WRITE(IW,9128) ( ZMAT(IZ),IZ=1,NZVAR)
      ENDIF
      HGUESS=GUESS2
      IF(OUT) THEN
         WRITE(IW,9008) HGUESS
      ENDIF
C     
      IDEGEN=0
      IZ=1
      DO 190 I=1,NZVAR
         ITYPE=IZMAT(IZ)
         IZ1 = 0
         IZ2 = 0
         IZ3 = 0
         IZ4 = 0
         IZA = 0
         IZB = 0
         IZC = 0
         IZD = 0
C     
C     BOND STRETCH, B-A
C     
         if (itype .eq. 1) then
            IZ1 = IZMAT(IZ+1)
            IZ2 = IZMAT(IZ+2)
            IZA = NUC(IZ1)
            IZB = NUC(IZ2)
            IF(BADZ(IZA)  .OR.  BADZ(IZB)) THEN
               HSSQ(I,I) = THIRD
               GO TO 180
            ENDIF
C     
            IROWA = IROW(IZA)
            IROWB = IROW(IZB)
            AIJ = APARM(IROWA,IROWB)
            DIJ = DPARM(IROWA,IROWB)
            RRR = TOANGS*ZMAT(I)
            HII = (AIJ-DIJ)/(RRR-DIJ)
            HESS1 = TOHB2 * HII * HII * HII
C     
            RAB = DIST(IZ1,IZ2)
            CAB = RCOV(IZA) + RCOV(IZB)
            HESS2 = 0.3601 * EXP(-1.944*(RAB - CAB))
C     
            IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS2
            IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
            IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
C     
C     ANGLE BEND, B-A-C
C     
         else if (itype.eq.2 .or. itype.eq.7) then
            IZ1 = IZMAT(IZ+1)
            IZ2 = IZMAT(IZ+2)
            IZ3 = IZMAT(IZ+3)
            IZB = NUC(IZ1)
            IZA = NUC(IZ2)
            IZC = NUC(IZ3)
            IF(BADZ(IZA)  .OR.  BADZ(IZB)  .OR.  BADZ(IZC)) THEN
               HSSQ(I,I) = QUARTR
               GO TO 180
            ENDIF
C     
            HESS1 = QUARTR
            IF(IZB.EQ.1  .OR.  IZC.EQ.1) HESS1 = 0.16D+00
C     
            RAB = DIST(IZ2,IZ1)
            RAC = DIST(IZ2,IZ3)
            CAB = RCOV(IZA) + RCOV(IZB)
            CAC = RCOV(IZA) + RCOV(IZC)
            HESS2 = 0.089 + 0.11 * (CAB*CAC)**0.42
     1           * EXP(-0.44*(RAB + RAC - CAB - CAC))
C     
            IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS2
            IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
            IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
C     
C     TORSION, D-A-B-C
C     
         else if (itype .eq. 3) then
            IZ1 = IZMAT(IZ+1)
            IZ2 = IZMAT(IZ+2)
            IZ3 = IZMAT(IZ+3)
            IZ4 = IZMAT(IZ+4)
            IZD = NUC(IZ1)
            IZA = NUC(IZ2)
            IZB = NUC(IZ3)
            IZC = NUC(IZ4)
            IF(BADZ(IZA).OR.BADZ(IZB).OR.BADZ(IZC).OR.BADZ(IZD))THEN
               HSSQ(I,I) = QUARTR
               GO TO 180
            ENDIF
C     
            HESS1 = QUARTR
C     
C     COUNT BONDS TO ATOMS A AND B, EXCEPT FOR A-B ITSELF
C     
            NBA = 0
            NBB = 0
            DO 135 J=1,NAT
               IZJ = NUC(J)
               IF(IZJ.LE.0) GO TO 135
               IF(IZJ.LE.86) THEN
                  CAJ = RCOV(IZA) + RCOV(IZJ)
                  CBJ = RCOV(IZB) + RCOV(IZJ)
               ELSE
                  CAJ = RCOV(IZA) + 1.70D+00
                  CBJ = RCOV(IZB) + 1.70D+00
               ENDIF
               IF(DIST(IZ2,J)  .LT.  1.2D+00*CAJ  .AND.
     1              J.NE.IZ2  .AND.  J.NE.IZ3) NBA = NBA+1
               IF(DIST(IZ3,J)  .LT.  1.2D+00*CBJ  .AND.
     1              J.NE.IZ3  .AND.  J.NE.IZ2) NBB = NBB+1
 135        CONTINUE
C     
C     THE FACTOR 14.0 IN PREEXP SEEMS LARGE, SO BRACKET THE RESULT
C     
            RAB = DIST(IZ2,IZ3)
            CAB = RCOV(IZA) + RCOV(IZB)
            PREEXP = 14.0 * (NBA+NBB)**0.57 / (RAB*CAB)**4
            HESS2 = 0.0015 + PREEXP * EXP(-2.85*(RAB-CAB))
            IF(HESS2.GT.0.30D+00) HESS2 = 0.30D+00
            IF(HESS2.LT.0.05D+00) HESS2 = 0.05D+00
C     
            IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS1
            IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
            IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
C     
C     OUT OF PLANE BEND, X-A-B-C, X IS HOOKED TO A
C     
         else if (itype .eq. 4) then
            IZ1 = IZMAT(IZ+1)
            IZ2 = IZMAT(IZ+2)
            IZ3 = IZMAT(IZ+3)
            IZ4 = IZMAT(IZ+4)
            IZX = NUC(IZ1)
            IZA = NUC(IZ2)
            IZB = NUC(IZ3)
            IZC = NUC(IZ4)
            IF(BADZ(IZA) .OR. BADZ(IZB) .OR. BADZ(IZC) .OR. 
     $           BADZ(IZX)) THEN
               HSSQ(I,I) = QUARTR
               GO TO 180
            ENDIF
C     
            HESS1 = QUARTR
C     
            CAX = RCOV(IZA) + RCOV(IZX)
            CAB = RCOV(IZA) + RCOV(IZB)
            CAC = RCOV(IZA) + RCOV(IZC)
            RAX = DIST(IZ1,IZ2)
            PREEXP = 0.0061 * (CAB*CAC)**0.80 * (COS(ZMAT(I)))**4.00
            HESS2 = 0.0025 + PREEXP * EXP(-3.00*(RAX-CAX))
C     
            IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS1
            IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
            IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
C     
C     ALL OTHER COORDINATE TYPES...
C     
         else 
            HSSQ(I,I)=QUARTR
         endif
C     
 180     CONTINUE
         IF(ITYPE.EQ.1) IZ=IZ+3
         IF(ITYPE.EQ.2) IZ=IZ+4
         IF(ITYPE.EQ.3) IZ=IZ+5
         IF(ITYPE.EQ.4) IZ=IZ+5
         IF(ITYPE.EQ.5) THEN
            IDEGEN=IDEGEN+1
            IF(IDEGEN.EQ.2) IZ=IZ+4
            IF(IDEGEN.EQ.2) IDEGEN=0
         ENDIF
         IF(ITYPE.EQ.6) IZ=IZ+6
         IF(ITYPE.EQ.7) IZ=IZ+5 ! RJH WAS 7 ... only place 7 appeared????
         IF(OUT) THEN
            WRITE(IW,9148) I,ITYPE,IZ1,IZ2,IZ3,IZ4,
     1           IZA,IZB,IZC,IZD,HSSQ(I,I)
         ENDIF
 190  CONTINUE
C     
c--   
c     DO J=1,NZVAR
c     DO I=1,NZVAR
c     HSSQ(I,J)=ZERO
c     ENDDO
c     HSSQ(J,J)=0.5D+00
c     ENDDO
c--   
      IF(OUT) THEN
         WRITE(IW,9048)
         CALL HND_PRSQ(HSSQ,NZVAR,NZVAR,NINT)
      ENDIF
C     
      RETURN
 9208 FORMAT(' IZMAT = ',12I4)
 9188 FORMAT(' IN  OPT_HSS-INIT_ZMT ,', 
     1     ' NAT,NVAR,NZVAR,NINT,NZMAT = ',5I4)
 9168 FORMAT(' IN  OPT_HSS_INIT_ZMT , GUESS INTERNAL HESSIAN ')
 9148 FORMAT(' I,ITYPE,HSSQ(I,I) = ',2I5,8I3,F12.8)
 9128 FORMAT('  ZMAT = ',5F7.3)
 9048 FORMAT(/,1X,38(1H-),
     1     /,1X,'INITIAL -INTERN.COORD.- HESSIAN MATRIX',
     2     /,1X,38(1H-))
 9008 FORMAT(' IN  OPT_HSS_INIT_ZMT , HGUESS = ',A8) 
      END
      subroutine driver_opt_check_hess(nvar, old_hessian)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
c
      integer nvar
      logical old_hessian
      character*255 filename
c
      integer m
c
c     Look at an existing hessian file and verify it
c
      call util_file_name('drv.hess',
     1     .false.,.false.,filename)
c
      if (ga_nodeid() .eq. 0) then
         open(32,file=filename,form='unformatted',status='old',err=10)
         read(32,err=11) m
         if (m.ne.nvar*nvar) goto 11
         close(32)
         old_hessian = .true.
         goto 20
c
 11      close(32)
 10      old_hessian = .false.
      endif
c
 20   call ga_brdcst(323, old_hessian, mitob(1), 0)
c
      end
      subroutine driver_opt_del_hess()
      implicit none
#include "util.fh"
c
c     Delete the Hessian information restart file.
c     
      character*255 opt_hess_fil
c
      call util_file_name('drv.hess',
     1     .false.,.false.,opt_hess_fil)
      call util_file_unlink(opt_hess_fil)
c
      if (util_print('information',print_low)) then
         write(6,*)
         write(6,*) ' Deleted DRIVER restart files '
         write(6,*)
      endif
c
      end
      subroutine driver_initialize(rtdb, geom)
      implicit none
#include "coptopt.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb
      integer geom              ! [output] 
c     
c     This routine initializes the common /coptopt/ and
c     also creates and returns the geometry handle
c     
      integer i, j, num
      logical status
      character*8 source, test
      logical geom_zmt_get_source, geom_zmtmak, geom_zmt_get_nizmat,
     $     geom_zmt_get_izmat, geom_zmt_get_nzvar, geom_zmt_get_nzfrz,
     $     geom_zmt_get_izfrz, geom_zmt_get_varinfo
c     
      oprint = util_print('information', print_low)
     $     .and. (ga_nodeid() .eq. 0)
      odebug = util_print('debug', print_debug)
c     
c     Before loading info regenerate the autoz coordinates if
c     the user has requested this
c     
      if (rtdb_get(rtdb,'driver:redoautoz',mt_log,1,status)) then
         if (ga_nodeid() .eq. 0) then
            status = rtdb_parallel(.false.)
            status = rtdb_delete(rtdb,'driver:redoautoz')
            call driver_opt_del_hess() ! Clear hessian info
            if (.not. geom_create(geom, 'geometry'))
     &           call errquit('hnd_opt: geom_create?', 911)
            if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &           call errquit('hnd_opt: no geometry ', 911)
            if (.not.geom_zmt_get_source(geom,source))
     $           call errquit('hnd_optx: zmat source ',geom)
c     
            if (source .eq. 'autoz') then
               if (oprint) write(6,777)
 777           format(/' Regenerating redundant internal coordinates'/)
               if (.not. geom_zmtmak(geom,oprint))
     $              call errquit('driver: regeneration of autoz failed',
     $              0)
               if (oprint) write(6,*)
               if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &              call errquit('hnd_opt: failed storing geometry',0)
            endif
            if (.not. geom_destroy(geom))
     $           call errquit('driver: geom_destroy?',0)
            status = rtdb_parallel(.true.)
         endif
         call ga_sync()
      endif
c     
c     ----- parameters for optimization driver -----
c     
      if (.not. rtdb_get(rtdb,'driver:trust',mt_dbl,1,trust))
     $     trust=0.3d0
      if (.not. rtdb_get(rtdb,'driver:cvgopt',mt_dbl,1,cvgopt))
     $     cvgopt=8.0d-04
      if (.not. rtdb_get(rtdb,'driver:nptopt',mt_int,1,nptopt))
     $     nptopt=20
      if (.not. rtdb_get(rtdb,'driver:inhess',mt_int,1,inhess))
     $     inhess=0
      if (.not. rtdb_get(rtdb,'driver:linopt',mt_int,1,linopt)) 
     $     linopt=0
      if (.not. rtdb_get(rtdb,'driver:modupd',mt_int,1,modupd))
     $     modupd=1
      if (.not. rtdb_get(rtdb,'driver:moddir',mt_int,1,moddir)) 
     $     moddir=1
      if (.not. rtdb_get(rtdb,'driver:modsad',mt_int,1,modsad)) 
     $     modsad=0
      if (.not. rtdb_get(rtdb,'driver:sadstp',mt_dbl,1,sadstp)) 
     $     sadstp=0.35d+00
      if (.not. rtdb_get(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $     oqstep = .true.
c     
c     Save a  copy of the initial geometry so we can analyze what 
c     happened during the optimization
c     
      if (ga_nodeid() .eq. 0) then
         status = rtdb_parallel(.false.)
         if (.not. geom_create(geom, 'geometry'))
     &        call errquit('hnd_opt: geom_create?', 911)
         if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &        call errquit('hnd_opt: no geometry ', 911)
         if (.not. geom_rtdb_store(rtdb, geom, 'driverinitial'))
     &        call errquit('hnd_opt: no geometry ', 911)
         if (.not. geom_destroy(geom))
     $        call errquit('driver: geom_destroy?',0)
         status = rtdb_parallel(.true.)
      endif
      call ga_sync()
c     
c     Load the geometry info
c     
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('hnd_opt: geom_create?', 911)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('hnd_opt: no geometry ', 911)
      if (.not. geom_ncent(geom,nat))
     $     call errquit('hnd_opt: natoms?',nat)
c     
c     Zmatrix info
c     
      if (.not.geom_zmt_get_source(geom,source))
     $     call errquit('hnd_optx: zmat source ',geom)
      zcoord  = (source .ne. ' ') ! Any kind of zmatrix
      zmatrix = (source .eq. 'user') ! zmatrix from user
      nzvar = 0
      nzfrz = 0
      if (zcoord) then
         if(.not.geom_zmt_get_nizmat(geom,nzmat))
     &        call errquit('geom_input: geom_zmt_get_nizmat failed',0)
         if(.not.geom_zmt_get_izmat(geom, izmat,nzmat))
     1        call errquit('hnd_opt_hss_init:zmt_get_izmat failed?',911)
         if(.not.geom_zmt_get_nzvar(geom,nzvar))
     &        call errquit('geom_input: geom_zmt_get_nzvar failed',0)
         call geom_bandbi(geom)
         if (.not. zmatrix) oqstep = .false. ! Only qstep in user z-matrix
         if(.not.geom_zmt_get_nzfrz(geom,nzfrz))
     &        call errquit('geom_input: geom_zmt_get_nzfrz failed',0)
         if(.not.geom_zmt_get_izfrz(geom,izfrz,nzfrz))
     $        call errquit('geom_input: geom_zmt_get_izfrz  failed',0)
      endif
      if (ga_nodeid().eq.0)  then
         write(6,9996) cvgopt,nptopt,linopt,inhess,modupd,
     $        modsad,moddir, sadstp, source, oqstep
 9996    format(' optimization convergence threshold    (cvgopt) = ',
     $        f10.6,
     1        /,' maximum number of energy evaluations  (nptopt) = ',
     $        i4,/,
     2        ' max.num. of points in linear search   (linopt) = ',i4,/,
     3        ' initial guess option for opt.hessian  (inhess) = ',i4,/,
     4        ' hessian update selection option       (modupd) = ',i4,/,
     5        ' optimization to saddle point          (modsad) = ',i4,/,
     6        ' initial eigen-mode followed to saddle (moddir) = ',i4,/,
     7        ' maximum step size to saddle           (sadstp) = ',f6.2/
     7        ' source of zmatrix                              =   ',a8/
     8        ' line search in internals              (oqstep) =    ',
     $        l1)
      endif
c     
      if ((ga_nodeid().eq.0) .and.
     $     util_print('information', print_low))  then
         write(6,9994)
 9994    format(/,10x,19(1h-),
     1        /,10x,'Energy Minimization',
     2        /,10x,19(1h-),/)
      endif
c     
c     Nvar is the no. of variables in the optimization
c     
      ncart = 3*nat
      if (zcoord) then
         nvar = nzvar
      else
         nvar = 3*nat
      endif
c     
c     Nindep is the no. of independent variables ... not needed?
c     
      if (zmatrix) then
         nindep = nzvar
      else
         if (nat .eq. 1) then
            nindep = 0
         else if (nat .eq. 2) then
            nindep = 1
         else
            nindep = 3*nat - 6  ! Linear molecules ????
         endif
      endif
c     
c     ----- frozen coordinates ? if so .... -----
c     
      if (oprint) then
         if(nzfrz.gt.0) then
            write(6,3) nzfrz, (izfrz(i),i=1,nzfrz) 
 3          format(/' Number of frozen variables =',i5/
     $           20(1x,15i5,/,:))
            if (.not. oqstep) write(6,*) 'freezing of variables',
     $           ' will not be exact since qstep is not used'
         endif
         if (zmatrix) then
            if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $           call errquit('hnd_optx: varinfo?',0)
            write(6,*)
            write(6,*) '  names of Z-matrix variables '
            do i = 1, nzvar, 5
               write(6,5432) (j,zvarname(j),j=i,min(nzvar,i+4))
 5432          format(1x,5(i4,1x,a8,2x))
            enddo
            write(6,*)
            write(6,*) 'variables with the same non-blank name',
     $           ' are constrained to be equal'
            if (.not. oqstep) write(6,*) 'constraints will be not be',
     $           ' exact since qstep is not used'
            write(6,*)
            do i = 1, nzvar
               if (zvarname(i).ne.' ') then
                  num = 0
                  test = zvarname(i)
                  do j = 1, nzvar
                     if (zvarname(j).eq.test) then 
                        num = num + 1
                     endif
                  enddo
                  if (num .gt. 1) then
                     write(6,5433) test(1:inp_strlen(test))
 5433                format(1x,'    variables named ',a8,' = ',$)
                     do j = 1, nzvar
                        if (zvarname(j).eq.test) then 
                           if (zvarsign(j).gt.0d0) then
                              write(6,5434) j
                           else
                              write(6,5434) -j
                           endif
 5434                      format(i5,$)
                           zvarname(j) = ' '
                        endif
                     enddo
                     write(6,*)
                  endif
               endif
            enddo
         endif
      endif
c     
c     We'll be printing the zmatrix and gradient so no
c     need for geom print to do this
c
      if (zcoord) call util_print_disable('geomzmat')
c
      energy = 0d0
      penergy= 0d0
      call dfill(max_nvar, 0d0, ds, 1)
      call dfill(max_nvar, 0d0, gx, 1)
      call dfill(max_nvar, 0d0, gq, 1)
      call dfill(max_nvar, 0d0, g, 1)
      call dfill(max_nvar, 0d0, gp, 1)
c     
      end
      subroutine driver_hessian_update()
      implicit none
#include "coptopt.fh"
#include "util.fh"
#include "mafdecls.fh"
c     
c     Update the current Hessian in the optimization variables using 
c     .   gp() - the gradient at the previous point
c     .    g() - the gradient at the current point
c     .   ds() - the search direction
c     .  alpha - the step in the search direction
c     
c     Only the Hessian is modified.
c     
      double precision hds(max_nvar)
      double precision dsds, dshds, dsdg
      integer l_hess, k_hess, i, j, iz
      integer ind
      ind(i,j) = k_hess + i + (j-1)*nvar - 1
c     
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('driver_hessian_update: memory for hessian',nvar**2)
      call geom_hnd_get_data('drv.hess',dbl_mb(k_hess), nvar**2)
c     
c     Form bits and pieces that are needed
c     
      call dgemv('n',nvar,nvar,1d0,dbl_mb(k_hess),nvar,
     $     ds,1,0d0,hds,1)
      dshds = ddot(nvar, ds, 1, hds, 1)
      dsds  = ddot(nvar, ds, 1,  ds, 1)
      dsdg  = 0d0
      do i = 1, nvar
         dsdg = dsdg + ds(i)*(g(i) - gp(i))
      enddo
c     
      if(modupd.le.1) then       
c     
c     ----- -bfgs- update -----
c     
         if(abs(dsdg).gt.1d-8) then
            do i=1,nvar
               do j=1,nvar
                  dbl_mb(ind(i,j))=dbl_mb(ind(i,j))
     $                 + (g(i)-gp(i))*(g(j)-gp(j))/dsdg
     1                 - hds(i)* hds(j)/dshds
               enddo   
            enddo   
         endif
      else      
c     
c     ----- -psb- update -----
c     
         do i=1,nvar
            do j=1,nvar
               dbl_mb(ind(i,j))=dbl_mb(ind(i,j))
     $              + ((g(i)-gp(i))-hds(i))*ds(j)/dsds
     1              + ((g(j)-gp(j))-hds(j))*ds(i)/dsds
     2              - ds(i)*ds(j)*(dsdg-dshds)/(dsds*dsds)
            enddo    
         enddo    
      endif    
c     
c     ----- damp out frozen coordinates -----
c     Since the gradient is already zeroed in the frozen coordinates
c     then the update should not have introduced any junk, but
c     just in case do this.  
c
c     Actually need to do this to the projected & shifted matrix.
c     
      do iz=1,nzfrz
         i=izfrz(iz)     
         do j=1,nvar
            dbl_mb(ind(i,j)) = 0d0
            dbl_mb(ind(j,i)) = 0d0
         enddo
      enddo    
      do iz=1,nzfrz
         i=izfrz(iz)     
         dbl_mb(ind(i,i)) = 0d0
      enddo
c     
      call geom_hnd_put_data('drv.hess',dbl_mb(k_hess), nvar**2)
      if (.not. ma_pop_stack(l_hess)) call errquit
     $     ('driver_hessian_update: ma?',0)
c     
      end
      subroutine driver_search_dir(geom)
      implicit none
#include "coptopt.fh"
#include "mafdecls.fh"
      integer geom
c
c     Put into ds() a search direction in the optimization
c     variables (internal or cartesian) based upon the
c     current gradient, g(), and Hessian.  Apply constraints.
c
c     Apply step restrictions by recommending an initial
c     value for the line search parameter alpha.
c
c     Only ds() is modified.
c
      integer i, info
      integer l_hess, k_hess, l_work, k_work, lenwork
      double precision  e(max_nvar) ! Hessian eigenvalues
      double precision pg(max_nvar) ! P.g
      double precision gv(max_nvar) ! Gradient along eigenvectors
      double precision dv(max_nvar) ! Step along eigenvectors
      double precision dsmax    ! Max. value of current step (smax is prev.)
c
      double precision trustds  ! Restriction of step in opt. variable
      parameter (trustds=0.2d0)
c
c     Get the Hessian and gradient with appropriate projects
c     applied following Peng, Ayala, Schlegel and Frisch so that
c     redundant internal modes are shifted to high eigenvalues.
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('driver_search_dir: memory for hessian',nvar**2)
      call driver_project_hess_grad(dbl_mb(k_hess), pg)
c
c     Diagonalize the Hessian.  Should really do the generalized
c     eigenvalue problem since the underlying basis is not independent
c     (if we are using autoz). NOT YET BEING DONE.  
c
      lenwork = max(nvar**2,100)
      if (.not. ma_push_get(mt_dbl, lenwork, 'work',
     $     l_work, k_work)) call errquit
     $     ('driver_search_dir: memory for hessian', lenwork)
      call dsyev('v', 'u', nvar, dbl_mb(k_hess), 
     $     nvar, e, dbl_mb(k_work), lenwork, info)
      if (info .ne. 0) call errquit
     $     ('driver_search_dir: dsyev failed',info)
      if (odebug) then
         write(6,*) ' Eigenvalues of the Hessian '
         call doutput(e, 1, nvar, 1, 1, nvar, 1, 1)
         write(6,*) ' Eigenvectors of the Hessian '
         call output(dbl_mb(k_hess), 1, nvar, 1, nvar, nvar, nvar, 1)
      endif
c
c     Now have eigenvalues in e, eigenvectors in dbl_mb(k_hess).
c
c     Project the gradient onto the Hessian eigenvectors
c
      call dgemv('t', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $     pg, 1, 0d0, gv, 1)
      if (odebug) then
         write(6,*) ' Gradient in spectral form '
         call doutput(gv, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     Form the step applying trust radius on each mode.
c
      do i = 1, nvar
         if (abs(e(i)) .gt. 1d-8) then
            dv(i) = -gv(i) / abs(e(i))
         else
c
c     Linearly dependent modes, and rotations and translations
c     have been shifted up.  Zero modes here are either accidental,
c     or dependencies that we don't know about.  As long as there is
c     not much gradient in the mode zero the step.  If there is
c     some gradient, then steepest descent with step restriction.
c
c     This would be more sensible in normalized coordinates.
c
            if (abs(gv(i)) .gt. 1d-2*cvgopt) then
               dv(i) = sign(trust,-gv(i))
            else
               dv(i) = 0d0
            endif
            if (oprint) write(6,1) i, e(i), dv(i)
 1          format(' Restricting step in mode',i5,
     $           ' due to tiny eigenvalue=', 1p,d8.1,
     $           ' step=',d8.1)
         endif
         if (abs(dv(i)) .gt. trust) then 
            if (oprint) write(6,2) i, e(i), dv(i), sign(trust,dv(i))
 2          format(' Restricting large step in mode',i5, 
     $           ' eval=', 1p,d8.1,' step=',d8.1,' new=',d8.1)
            dv(i) = sign(trust,dv(i))
         endif
      enddo
      if (odebug) then
         write(6,*) ' Step in spectral form '
         call doutput(dv, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     Transform back to optimization space
c
      call dgemv('n', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $     dv, 1, 0d0, ds, 1)
      if (odebug) then
         write(6,*) ' Step in optimization variables'
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     Restriction of the overall step based on
c     
c     1.  Did we go uphill on the last step?
c     2.  Is a component of the step large in the optimization space.
c
      dsmax = 0d0
      do i = 1, nvar
         dsmax = max(dsmax, abs(ds(i)))
      enddo
c
      alpha = 1d0
      if (dsmax .gt. trustds) then
         alpha = trustds/dsmax
         if (oprint) write(6,7) alpha
 7       format('  Restricting overall step due to large component.',
     $        ' alpha=', f6.2)
      endif
      if ((energy-penergy) .gt. eprec) then
         alpha = 0.5d0*alpha
         if (oprint) write(6,8) alpha
 8       format('  Restricting overall step due to uphill motion.  ',
     $        ' alpha=', f6.2)
      endif
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('driver_search_dir: ma?',0)
c
c     Edo seems to have encountered a case where difference processors
c     generated different steps.  To prevent this, broadcast the 
c     critical info to everyone.
c
      call ga_brdcst(1,ds,8*nvar,0)
      call ga_brdcst(2,alpha,8,0)
c
      end
      subroutine driver_line_search(rtdb,geom)
      implicit none
#include "coptopt.fh"
#include "geom.fh"
#include "util.fh"
      integer rtdb, geom
c
c     Starting with a value of alpha recommended by whatever
c     generated the search direction, perform an approximate line 
c     search using no more than linopt steps.  Accept the first 
c     point that is lower in energy than the starting point.
c
c     Enter with the supposition that any preconditioning of
c     the gradient is positive-definite so that the minimum
c     is found by a positive step along the search direction.
c
c     Return with alpha, energy and geom as the accepted point.
c
c     Only alpha, energy and geom are modified.
c
      double precision enew, anew, a, dsgrad
      integer lstep
      double precision driver_energy_step
c
      dsgrad = ddot(nvar, ds, 1, g, 1)
      if (dsgrad*alpha .ge. 0d0) then
         write(6,*) ' Warning ... line search gradient +ve ',
     $        alpha, dsgrad
         call errquit('fred',0)
      endif
c
      do lstep = 1, linopt
c
c     Compute energy at alpha*ds
c
         enew = driver_energy_step(rtdb, geom)
         write(6,*) ' ALPHA ', lstep, alpha, enew
         if ((enew - energy) .lt. eprec ) goto 1000 ! Downhill.  Accept.
c
c     Must have overshot.  Downhill somewhere in [0,alpha].
c
c     Fit to a parabola using e(0), g(0), e(alpha)
c
         a = (enew-energy-alpha*dsgrad)/alpha**2
         if (a .le. 0d0) then
c
c     Not positive curvature. Bisect.
c
            write(6,*) ' Negative curvature ... bisecting'
            alpha = 0.5d0*alpha
         else
c
c     Parabolic extrapolation
c
            anew = -dsgrad/(2d0*a)
            if (anew .ge. alpha) then 
               write(6,*) ' Overshooting ... bisecting'
               alpha = 0.5d0*alpha ! Downhill somewhere in [0,alpha]
            else
               alpha = anew
            endif
         endif
      enddo
c
c     Mmmm ... if we fell out of the loop we did not find a downhill 
c     step.  Not much to do about this ... just take the last point 
c     which will have been a small step
c
c     Jump to 1000 if went downhill.  Set energy and geometry
c
 1000 energy = enew
      call driver_take_step(rtdb, geom)
c
      end
      double precision function driver_energy_step(rtdb, geom)
      implicit none
#include "geom.fh"
#include "rtdb.fh"
#include "coptopt.fh"
#include "mafdecls.fh"
      integer rtdb, geom
c
c     Compute the energy by displacing the cartesian coordinates
c     in geom by a step of alpha.  After computing the energy
c     restore the original geometry.
c
c     Nothing is modified.  The energy is returned.
c
      double precision xold(max_cart)
      double precision aenergy
      logical task_energy
      external task_energy
c
c     Get original coordinates so can restore them
c
      if (.not. geom_cart_coords_get(geom, xold))
     $     call errquit('driver_energy_step: coordinates?',geom)
c
c     Take the step
c
      call driver_take_step(rtdb, geom)
c
c     Energy at displaced point
c
      if (.not. task_energy(rtdb)) call errquit
     $     ('driver_energy_step: energy failed',0)
      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, aenergy))
     $     call errquit('driver_energy_step: rtdb?',0)
c      
c     Restore the starting geometry
c
      if (.not. geom_cart_coords_set(geom, xold))
     $     call errquit('driver_energy_step: coordinates?',geom)
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $     call errquit('driver_energy_step: grs?',geom)
c
      driver_energy_step = aenergy
c
      end
      subroutine driver_take_step(rtdb, geom)
      implicit none
#include "coptopt.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb, geom
c
c     Update the geometry in geom and in the database
c     'geometry' by taking the step
c     alpha*ds() in the optimization variables
c
c     The geom is modified, and xmax/xrms are computed.
c
      double precision xold(max_cart), xnew(max_cart)
      integer i, l_bi, k_bi
c
c     Get original coordinates
c
      if (.not. geom_cart_coords_get(geom, xold))
     $     call errquit('driver_energy_step: coordinates?',geom)
c
c     Take the step
c
      if (zcoord) then
c
c     First order step only for now.
c
         if (.not. ma_push_get(mt_dbl, ncart*nzvar,'binv',l_bi, k_bi))
     $        call errquit('driver_step: memory for binv', ncart*nzvar)
         call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nzvar)
         call dgemv('n', ncart, nzvar, alpha, dbl_mb(k_bi), ncart,
     $        ds, 1, 0.0d0, xnew, 1)
         if (.not. ma_pop_stack(l_bi)) call errquit('driver: ma?',99)
         if (odebug) then
            write(6,*) ' First-order internal coord step in cartesians'
            call doutput(xnew, 1, ncart, 1, 1, ncart, 1, 1)
         endif
         do i = 1, ncart
            xnew(i) = xold(i) + xnew(i)
         enddo
         if (.not. geom_cart_coords_set(geom, xnew))
     $        call errquit('driver_energy_step: coordinates?',geom)
      else
         do i = 1, ncart
            xnew(i) = xold(i) + alpha*ds(i)
         enddo
         if (.not. geom_cart_coords_set(geom, xnew))
     $        call errquit('driver_energy_step: coordinates?',geom)
      endif
c
c     Must ensure the geometry has the required symmetry.  Should use
c     an error criterion consistent with the step size.
c
      call sym_geom_project(geom, 0.1d0)
c
c     Compute the maximum and RMS cartesian displacements
c     
      xmax = 0d0
      xrms = 0d0
      do i = 1, ncart
         xmax = max(xmax, abs(xold(i)-xnew(i)))
         xrms = xrms + (xold(i)-xnew(i))**2
      enddo
      xrms = sqrt(xrms/dble(ncart))
c
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $     call errquit('driver_energy_step: grs?',geom)
c
      end
      subroutine driver_print(geom, istep)
      implicit none
#include "coptopt.fh"
#include "global.fh"
#include "util.fh"
      integer geom, istep
c
c     Print out stuff
c
      integer i
      double precision diff, de
      logical geom_print_zmatrix, oprint_grad
c
      de = 0d0
      if (istep .gt. 1) de = energy-penergy
      if (oprint) write(6,1) istep, energy, de,
     $     gmax, grms, xrms, smax, util_wallsec()
 1    format(
     $     /,1x,' Step   New Energy   Delta E   Gmax',
     $     '     Grms     Xrms     Smax   Walltime',
     $     /,1x,' ---- ------------- -------- --------',
     $     ' -------- -------- -------- --------',
     $     /,1x,i5,f14.6,1p,d9.1,0p,4f9.5,f9.1/)
c
c     Optionally print out the gradient in internal coordinates.
c
      oprint_grad = zcoord .and. (ga_nodeid().eq.0) .and. 
     $     util_print('internal gradient', print_default)
      if (oprint_grad) then
         diff = 0d0             ! See if anything has been constrained
         do i = 1, nvar
            diff = diff + abs(g(i)-gq(i))
         enddo
         if (diff .gt. 1d-6)
     $        write(6,*) ' Internal coordinate gradient before',
     $        ' constraints'
         if (.not. geom_print_zmatrix(geom, gq, 'Gradient',.true.))
     $        call errquit('tfgx: print internal gradient?',geom)
         if (diff .gt. 1d-6) then
            write(6,*) ' Internal coordinate gradient after',
     $           ' constraints'
            if (.not. geom_print_zmatrix(geom, gq, 'Gradient',.true.))
     $           call errquit('tfgx: print internal gradient?',geom)
         endif
      endif
c
      end
      logical function driver_converged()
      implicit none
#include "coptopt.fh"
c
c     Return true if we have converged
c
c     Nothing is modified.  Assumes driver_compute_info()
c     has been called.
c
      driver_converged = (gmax .lt. cvgopt)
c
      end
      subroutine driver_cart_pmat(geom)
      implicit none
#include "coptopt.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer geom
c
c     Compute the cartesian equivalent of the P = G.G^-1 matrix
c     which projects to and from the linearly independent
c     set of coordinates.  In the cartesian case P is the complement
c     of the projector onto the rotations and translations
c
c     Only the P matrix is generated.  Nothing is modified.
c
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Minor little catch is that if some atoms are being frozen
c     we are no longer invariant to translations or rotations.
c     NOT CHECKING THIS YET.
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     For ease of use we also write out a unit Binv matrix.
c
      double precision centroid(3), x, y, z, xx, yy, zz, fx
      double precision coords(3,max_cent)
      double precision work(max_cart,6)
      integer i, j, k, l_pmat, k_pmat, i3
      integer ind
      ind(i,j) = k_pmat + i-1 + (j-1)*ncart
c
      if (.not. geom_cart_coords_get(geom, coords))
     $     call errquit('driver_cart_pmat: geom?',geom)
c
c     Construct normalized vectors in work in the direction
c     of the rotations and translations.
c
      call dfill(3, 0.0d0, centroid, 1)
      do i = 1, nat
         do k = 1, 3
            centroid(k) = centroid(k) + coords(k,i)/nat
         enddo
      enddo
c     
      do k = 1, 3               ! x, y, z translations
         call dfill(ncart, 0.0d0, work(1,k), 1)
         call dfill(nat, sqrt(1.0d0/nat), work(k,k), 3)
      enddo
      do k = 4, 6               ! x, y, z rotations
         do i = 1, nat
            x = coords(1,i) - centroid(1)
            y = coords(2,i) - centroid(2)
            z = coords(3,i) - centroid(3)
            if (k .eq. 4) then
               xx = 0.0d0
               yy = -z
               zz =  y
            else if (k .eq. 5) then
               xx =  z
               yy =  0.0d0
               zz = -x
            else if (k .eq. 6) then
               xx = -y
               yy =  x
               zz =  0.0d0
            endif
            i3 = (i-1)*3
            work(i3+1,k) = xx
            work(i3+2,k) = yy
            work(i3+3,k) = zz
         enddo
         do j = 1, k-1
            fx = ddot(ncart, work(1,j), 1, work(1,k), 1)
            call daxpy(ncart, -fx, work(1,j), 1, work(1,k), 1)
         enddo
         fx = sqrt(ddot(ncart, work(1,k), 1, work(1,k), 1))
         if (fx . gt. 1d-6) then
            call dscal(ncart, 1.0d0/fx, work(1,k), 1)
         else
            call dfill(ncart, 0.0d0, work(1,k), 1)
         endif
      enddo
c     
c     The project is then 1 - V.VT where V is in work
c
      if (.not. ma_push_get(mt_dbl, ncart**2, 'pmat',
     $     l_pmat, k_pmat)) call errquit
     $     ('driver_cart_pmat: memory for pmat',ncart**2)
c
c     Form unit matrix
c
      call dfill(ncart**2, 0d0, dbl_mb(k_pmat), 1)
      call dfill(ncart, 1d0, dbl_mb(k_pmat), ncart+1)
c
c     Store dummy unit matrix for Binv ... the cartesian 
c     gradient should already be invariant to rotations and translations.
c
      call geom_hnd_put_data('b^-1', dbl_mb(k_pmat), ncart**2)
c
c     Finish P
c
      do i = 1, ncart
         do j = 1, ncart
            do k = 1, 6
               dbl_mb(ind(j,i)) = dbl_mb(ind(j,i)) - 
     $              work(j,k)*work(i,k)
            enddo
         enddo
      enddo
c
      call geom_hnd_put_data('p', dbl_mb(k_pmat), ncart**2)
c
      if (odebug) then
         write(6,*) ' Cartesian P matrix'
         call output(dbl_mb(k_pmat),1,ncart,1,ncart,ncart,ncart,1)
      endif
c
      if (.not. ma_chop_stack(l_pmat)) call errquit
     $     ('driver_cart_bmat: ma?',0)
c
      end
      subroutine driver_project_hess_grad(hess, pg)
      implicit none
#include "coptopt.fh"
#include "mafdecls.fh"
      double precision 
     $     hess(nvar,nvar),     ! returns projected & shifted Hessian
     $     pg(nvar)             ! returns projected gradient
c
c     Project and shift the Hessian and gradient following Peng et al.
c
      integer l_pmat, k_pmat, l_work, k_work, i
      double precision big
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'work',
     $     l_work, k_work)) call errquit
     $     ('driver_proj_h_g: memory for pmat',nvar**2)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'pmat',
     $     l_pmat, k_pmat)) call errquit
     $     ('driver_proj_h_g: memory for work',nvar**2)
c
      call geom_hnd_get_data('drv.hess',hess, nvar**2)
      if (odebug) then
         write(6,*) ' Hessian before projection'
         call output(hess, 1, nvar, 1, nvar, nvar, nvar, 1)
         write(6,*) ' Gradient before projection'
         call doutput(g, 1, nvar, 1, 1, nvar, 1, 1)
      endif
      call geom_hnd_get_data('p',dbl_mb(k_pmat), nvar**2)
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('freddy',0)
c
c     PG
c
      call dgemv('n',nvar, nvar, 1d0, dbl_mb(k_pmat), nvar,
     $     g, 1, 0d0, pg, 1)
      if (odebug) then
         write(6,*) ' Gradient after projection'
         call doutput(g, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     PHP + 1000*(1-P)
c
      call dgemm('n', 'n', nvar, nvar, nvar, 1d0, dbl_mb(k_pmat), nvar,
     $     hess, nvar, 0d0, dbl_mb(k_work), nvar)
      call dgemm('n', 'n', nvar, nvar, nvar, 1d0, dbl_mb(k_work), nvar,
     $     dbl_mb(k_pmat), nvar, 0d0, hess, nvar)
      if (odebug) then
         write(6,*) ' Hessian after projection before shift'
         call output(hess, 1, nvar, 1, nvar, nvar, nvar, 1)
      endif
c     
      big = 1000d0
      call daxpy(nvar*nvar, -big, dbl_mb(k_pmat), 1, hess, 1)
      do i = 1, nvar
         hess(i,i) = hess(i,i) + big
      enddo
      if (odebug) then
         write(6,*) ' Hessian after projection & shift'
         call output(hess, 1, nvar, 1, nvar, nvar, nvar, 1)
      endif
c
      if (.not. ma_chop_stack(l_work)) call errquit
     $     ('driver_p_h_g:ma?',0)
c
      end
      subroutine driver_compute_info()
      implicit none
#include "coptopt.fh"
c
c     Compute stuff used for printing and convergence tests
c
c     gmax = maxmimum gradient element in optimization variables
c     grms = rms grad
c     smax = maximum step in opt. var
c     srms = rms step
c     
c     xrms and xmax are computed by driver_take_step from the
c     first-order step.
c
      integer i
c
      grms = 0d0
      srms = 0d0
      gmax  = 0d0
      smax  = 0d0
c
      do i = 1, nvar
         grms = grms + g(i)*g(i)
         srms = srms + ds(i)*ds(i)*alpha*alpha
         gmax  = max(gmax, abs(g(i)))
         smax  = max(smax, abs(ds(i)*alpha))
      enddo
      grms = sqrt(grms/dble(nvar))
      srms = sqrt(srms/dble(nvar))
c
      end

         
