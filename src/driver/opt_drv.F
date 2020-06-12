C>
C> \ingroup driver
C> @{
C>
C> \file opt_drv.F
C> The top-level Driver code
C>
C> @}
C>
C> \ingroup driver
C> @{
C>
      logical function driver(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb
c     
      integer geom, geomi, i, iat
      integer istep
      logical converged, status
      double precision start    ! Tracks time used in last step
      integer required          ! Time required
      logical ophigh
      logical  driver_converged, task_gradient, 
     $     geom_compute_zmatrix, geom_print_zmatrix
      external driver_converged, task_gradient,
     $     geom_compute_zmatrix, geom_print_zmatrix
c
      ophigh = util_print('high', print_high)
c
c     Read input, load /coptopt/, get geometry
c
      call driver_initialize(rtdb, geom)

c
c     Optionally check the gradient by finite difference
c
      if (ocheckgrad) call driver_check_grad(rtdb,geom)
c     
      call driver_hss_init(rtdb,geom)
c
c     Iterative loop
c
      do istep = 1, nptopt+1    ! +1 since first pass thru loop is not a step
         start = util_wallsec()
         if (oprint) write(6,1) istep-1
 1       format(/,10x,'--------',/,10x,'Step',i4,/,10x,'--------')
         if ((ga_nodeid() .eq. 0) .and.
     $        util_print('geometry',print_default)) then
            if (.not. geom_print(geom)) call errquit('driver: geom?',0,
     &       GEOM_ERR)
         endif
c     
c     Save old energy, gradient, searchdir and alpha, and compute 
c     new energy and gradient.
c
         alphap  = alpha        ! Used in conjugacy and starting LS
         energyp3= energyp2
         energyp2= energyp
         energyp = energy       ! Used for convergence and step restriction
         call dcopy(nvar, g, 1, gp, 1) ! Used for Hessian update
         call dcopy(nvar, ds, 1, dsp, 1) ! Used for conjugacy
c
         if (.not. task_gradient(rtdb))
     $        call errquit('driver: task_gradient failed',0, GEOM_ERR)
         call opt_get_grad(rtdb,geom) ! Into gx
         if (ostress2) then
           call sym_lattice_stress2(geom,gx(nat_real*3+1))
         endif
         if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $        call errquit('driver: could not get energy',0, RTDB_ERR)
c
c     Disable printing to ecce of movecs after first point
c
         call movecs_ecce_print_off()
c     
c     Zero the gradient associated with atoms frozen in cartesians
c     
       if (ga_nodeid().eq.0.and.ophigh) 
     %        write(luout,*) 'Zeroing constrained gradient'

         if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
            do iat = 1, nat_real
               if (.not. oactive(iat)) then
                  do i = 1, 3
c                     write(6,*) ' Zeroing constrained gradient', 
c     $                    iat, i, gx((iat-1)*3+i)
                     gx((iat-1)*3+i) = 0.0
                  end do
               end if
            end do
         end if
c     
c     If necessary transform gradient to internal coordinates
c     
         if(zcoord) then
            call geom_bandbi(geom)
            call driver_tfgx(geom) ! g() will be set to constrained gq()
         else
            call driver_cart_pmat(rtdb, geom)
            call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
         endif
c     
c     Update the Hessian
c     
         if (istep .gt. 1) call driver_hessian_update(geom,rtdb)
c
c     Compute step/gradient info and print for user
c     (for current energy & gradient, and the previous alpha*step).
c
         call driver_compute_info()
         call driver_print(geom,istep)
c     
c     Check for convergence
c     
         if (istep.gt.1 .and. driver_converged()) then
            converged = .true.
            goto 100
         endif

        if(istep.gt.nptopt) goto 200
c
c     Generate a new search direction
c
         if (modsad .eq. 0) then
            call driver_opt_search_dir(geom) ! Fills in ds(), recommends alpha
         else
            call driver_sad_search_dir(rtdb,geom,istep) ! Fills in ds(), alpha
         endif
c     
c     Either do line search or take recommended step.
c     
         if (linopt.le.0) then
            call driver_take_step(rtdb, geom) ! Updates geom using alpha*ds
         else if (linopt .eq. 1) then
            call driver_line_search1(rtdb, geom) ! Updates geom,alpha
         else if (linopt .eq. 2) then
            call driver_line_search2(rtdb, geom) ! Updates geom,alpha
         else
            if (istep.gt.1 .and. (
     >             (energy-energyp) .gt. eprec .or.
     >             (energyp-energyp2) .gt. eprec)) then
               call driver_line_search2(rtdb, geom) ! More aggressive, if uphill
            else
               call driver_line_search1(rtdb, geom)
            endif
         endif
c
c     We have now taken a step.  Replace the approximate step taken
c     by the exact step in case update of internals was not exact.
c     This is only significant if big steps are being taken
c     (e.g., saves 30 iterations in h2o5 test case).
c     
         call driver_compute_actual_step(geom)
c
c     Check time before next iteration
c
         required = int(1.2d0*(util_wallsec() - start)) + 1
         if (.not. util_test_time_remaining(rtdb,required)) goto 200
c
      enddo                     ! End of iterative loop
c      istep = istep - 1         ! Since we fell out
 200  if (oprint) write(6,201)
 201  format(/,1x,63('-')/,1x,'Failed to converge in maximum number',
     $     ' of steps or available time'/,1x,63('-')/)
      converged = .false.
c
 100  if (converged) then
         if (oprint) write(6,101) 
 101     format(/,6x,22('-'),/,6x,'Optimization converged',/,
     $        6x,22('-'),/)
      endif
c
      if (ga_nodeid().eq.0 .and. util_print('finish',print_low)) then
c
c     Print out final info and geometry
c
         call driver_print(geom,istep)
         if (.not. geom_print(geom)) call errquit
     $        ('hnd_opt_drv: geom_print?',0, GEOM_ERR)
c
c     Analyze the change in internal coordinates
c     (overwrites ds and g)
c
         if (zcoord) then
            status = rtdb_parallel(.false.)
            if (.not. geom_create(geomi, 'driverinitial'))
     &           call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
            if (.not. geom_rtdb_load(rtdb, geomi, 'driverinitial'))
     &           call errquit('hnd_opt: no initial geometry ', 911,
     &       RTDB_ERR)
            if (.not. geom_compute_zmatrix(geomi, ds)) ! Original coords
     $           call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
            if (.not. geom_compute_zmatrix(geom, g)) ! Final coords
     $           call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
            do i = 1, nvar
               g(i) = g(i) - ds(i)
               if (g(i) .gt. 180d0) then
                  g(i) = g(i) - 360d0
               else if (g(i) .lt. -180d0) then
                  g(i) = g(i) + 360d0
               endif
            enddo
            write(6,*)
            call util_print_centered(6,
     $           'Final and change from initial internal coordinates ',
     $           40, .true.)
            write(6,*)
            if (.not. geom_print_zmatrix(geom, g, 'Change ',.true.))
     $           call errquit('driver:print geom change?',geom,
     &       GEOM_ERR)
            if (.not. geom_destroy(geomi))
     $           call errquit('driver:geomi corrupt',0, GEOM_ERR)
            status = rtdb_parallel(.true.)
         endif
c
         if (util_print('bonds',print_default)) then
            if (.not.geom_print_distances(geom)) call errquit(
     &           'hnd_opt_drv: geom_print_distances failed',911,
     &       GEOM_ERR)
         endif
         if (util_print('angles',print_default)) then
            if (.not.geom_print_angles(geom)) call errquit(
     &           'hnd_opt_drv: geom_print_angles failed',911,
     &       GEOM_ERR)
         endif
      endif
c
c     Clean up and go home
c
      if (.not.geom_destroy(geom))
     &     call errquit('hnd_opt: geom_destroy?', 911, GEOM_ERR)
c
      driver=converged
      if (converged) then
         call ecce_print_module_exit('driver', 'ok')
      else
         call ecce_print_module_exit('driver', 'failed')
      endif
c
      call movecs_ecce_print_on() ! Restore MO printing
      call util_print_pop
c
      call ga_sync()
c
      end
      subroutine driver_tfgx(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
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
      integer l_bi, k_bi, num
c
c     gq = gx*B^-1
c
      num = max(ncart*nvar,nvar**2)
      if (.not. ma_push_get(mt_dbl, num,'binv',l_bi, k_bi))
     $     call errquit('driver: memory for binv', num, MA_ERR)
      call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nvar)
      if (odebug) then
         write(6,*) ' B^-1 '
         call output(dbl_mb(k_bi), 1, ncart, 1, nvar, nvar, nvar, 1)
         write(6,*) ' Cartesian gradient'
         call doutput(gx, 1, ncart, 1, 1, ncart, 1, 1)
      endif
      call dgemv('t', ncart, nvar, 1.0d0, dbl_mb(k_bi), ncart,
     $     gx, 1, 0.0d0, gq, 1)
      if (odebug) then
         write(6,*) ' Internal coordinate gradient '
         call doutput(gq, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     g = C*gq 
c
c     Projecting with C only works for constraints in non-redundant
c     variables ... must project with P in redundant case.
c
      call geom_hnd_get_data('p', dbl_mb(k_bi), nvar**2)
      call dgemv('n', nvar, nvar, 1.0d0, dbl_mb(k_bi), nvar,
     $     gq, 1, 0.0d0, g, 1)
      if (odebug) then
         write(6,*) ' Internal coordinate gradient after constraints'
         call doutput(gq, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
      if (.not. ma_pop_stack(l_bi)) call errquit('driver: ma?',99,
     &       MA_ERR)
c
      end
      subroutine driver_hss_init(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
c
      integer rtdb, geom

c
      logical geom_zmt_get_izmat, geom_zmt_get_nizmat, 
     $     geom_compute_zmatrix
      external geom_zmt_get_izmat, geom_zmt_get_nizmat,
     $     geom_compute_zmatrix
c
      double precision zero
      parameter (zero=0.0d+00)
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=nw_max_atom)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=nw_max_zmat)     
      parameter (mxcoor=nw_max_coor)
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
      double precision hscale, ascale, bscale, tscale, amat(3,3)
c
      integer l_hess, k_hess, l_zmat, k_zmat, l_izmat, k_izmat, i, j, k
      integer l_c, k_c, l_t, k_t, iat, jat, kat, lat, l
      integer l_shrhess, k_shrhess
      logical old_hessian
      character*16 atom_tags(mxatom)
c
      nnzmat = nzmat
      nnzvar = nzvar
      nnvar  = nzvar
      if (.not. geom_ncent(geom,natom))
     1       call errquit('hnd_opt: geom_ncent?',911, GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_hess, k_hess)) call errquit
     $     ('driver_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      old_hessian=.false.
      if(inhess.ne.1) call driver_opt_check_hess(nvar, old_hessian)
      if (oprint) write(6,*)
      if (old_hessian) then
         if (oprint) write(6,*) 'Using old Hessian from',
     $        ' previous optimization'
         goto 999
      else if (inhess.eq.2) then
         call driver_hess_cart_guess()
         if (oprint) write(6,*)
     $        ' Using Cartesian Hessian from previous frequency',
     $        ' calculation'
         goto 999
      else
         if (oprint) write(6,*) 'Using diagonal initial Hessian '
      endif
c
      if (.not. zcoord) then
c
c     Cartesians are easy
c
         call dfill(nvar**2, 0.0d0, dbl_mb(k_hess), 1)
         call hnd_opt_hss_init_cart(dbl_mb(k_hess),
     &                     odebug,isystype,geom,nvar)
!    Now the lattice hessian, if any
         if (isystype .ne. 0 .and. (ostress .or. ostress2)) then
            call HND_OPT_HSS_INIT_LATTICE(dbl_mb(k_hess),nvar,isystype,
     $           ostress,ostress2,geom)
            if (oprint) then 
               write(6,*) ' The initial hessian '
               call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
            end if
         end if
!        Symmetry
         if (nvar .eq. 3*natom) then
           call sym_hess_symmetrize(geom,dbl_mb(k_hess))
         else
           if (.not. ma_push_get(mt_dbl, (natom*3)**2, 'shrhess',
     $        l_shrhess, k_shrhess)) call errquit
     $        ('driver_init_hess: failed allocating shrunk hessian ',
     $        (natom*3)**2, MA_ERR)
           do i = 1, 3*natom
             do j = 1, 3*natom
               dbl_mb(k_shrhess+j-1 + (i-1)*(3*natom)) = 
     $         dbl_mb(k_hess   +j-1 + (i-1)*(nvar    ))
             end do
           end do
           call sym_hess_symmetrize(geom,dbl_mb(k_shrhess))
           do i = 1, 3*natom
             do j = 1, 3*natom
               dbl_mb(k_hess   +j-1 + (i-1)*(nvar    )) =
     $         dbl_mb(k_shrhess+j-1 + (i-1)*(3*natom))
             end do
           end do
         endif
      else
c
c     Internals.
c
         if (.not. ma_push_get(mt_dbl, nzvar, 'zmat',
     $        l_zmat, k_zmat)) call errquit
     $        ('driver_init_hess: failed allocating zmat ', nzvar,
     &       MA_ERR)
         if (.not. ma_push_get(mt_int, nzmat, 'izmat',
     $        l_izmat, k_izmat)) call errquit
     $        ('driver_init_hess: failed allocating izmat ', nzmat,
     &       MA_ERR)
         if (.not. geom_compute_zmatrix(geom, dbl_mb(k_zmat)))
     $        call errquit('driver_init_hess: zmatrix?',0, GEOM_ERR)
         if (.not. geom_zmt_get_izmat(geom,int_mb(k_izmat),nzmat))
     $        call errquit('driver_init_hess: izmat?',0, GEOM_ERR)
C  The code assume Cartesians in many places:
         if (ostress .or. ostress2) call errquit(
     $     'driver: lattice optimization does not work with a zmatrix',
     $     0,GEOM_ERR)
c
         if (.not. ma_verify_allocator_stuff()) call errquit('h',0,
     &       MA_ERR)
c
C     ----- INTERNALLY GENERATED (DIAGONAL) INITIAL HESSIAN -----
C
c     this routine uses /molnuc/, /zmtpar/, /molxyz/ which 
c     must be initialized before calling
c
         if (.not. geom_cart_get(geom, nat, atom_tags, 
     $        c, zan)) call errquit('driver_init_hess: geom?',geom,
     &       GEOM_ERR)
         do i = 1, nat
            nuc(i) = nint(zan(i))
         enddo
         nvar = nzvar
c
c     The angle scaling factors are responsible for a factor of two
c     speedup in most organics.  Clearly we need a better force field.
c
c     The chosen values have been tested on a large no. of 3-21g* SCF
c     optimizations, and a smaller no. of 6-31g* SCF calculations.
c
         if (.not. rtdb_get(rtdb,'driver:bscale',mt_dbl,1,bscale))
     $        bscale = 1d0
         if (.not. rtdb_get(rtdb,'driver:ascale',mt_dbl,1,ascale))
     $        ascale = 0.25d0
         if (.not. rtdb_get(rtdb,'driver:tscale',mt_dbl,1,tscale))
     $        tscale = 0.1d0
         if (oprint) then
            write(6,789) bscale, ascale, tscale
 789        format(' Scaling for Hessian diagonals: bonds =',f5.2,
     $           '  angles =',f5.2,'  torsions =',f5.2)
         endif
         call hnd_opt_hss_init_zmt(dbl_mb(k_hess),
     $        nzvar, dbl_mb(k_zmat), int_mb(k_izmat), .false.,
     $        bscale, ascale, tscale,odebug)
         if (odebug) then
            write(6,*) ' The initial hessian '
            call output(dbl_mb(k_hess),1,nzvar,1,nzvar,nzvar,nzvar,1)
         endif
      endif
c
c     Artificially break degeneracies so that accidentally degenerate
c     modes are split and therefore step restriction along modes
c     is well defined
c
      do i = 1, nvar
         dbl_mb(k_hess+i-1 + (i-1)*nvar) = 
     $        dbl_mb(k_hess+i-1 + (i-1)*nvar) + dble(i-1)*1e-7
      end do
c
      call geom_hnd_put_data('drv.hess', dbl_mb(k_hess), nvar*nvar)
c
c     Apply constants, constraints and overall scaling
c
 999  if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_c, k_c)) call errquit
     $     ('driver_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_t, k_t)) call errquit
     $     ('driver_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      call geom_hnd_get_data('drv.hess', dbl_mb(k_hess), nvar*nvar)
c
c     Used to use c here ... now use p
c
      call geom_hnd_get_data('p',dbl_mb(k_c), nvar*nvar)
c
      if (odebug) then
         write(6,*) ' Initial Hessian before P'
         call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
      endif
      call dgemm('n','n',nvar,nvar,nvar,1d0,dbl_mb(k_c),nvar,
     $     dbl_mb(k_hess),nvar,0d0,dbl_mb(k_t),nvar)
      call dgemm('n','t',nvar,nvar,nvar,1d0,dbl_mb(k_t),nvar,
     $     dbl_mb(k_c),nvar,0d0,dbl_mb(k_hess),nvar)
      if (odebug) then
         write(6,*) ' Initial Hessian after P'
         call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
      endif
c
      if (nactive .ne. nat_real) then
c     
c     We are in cartesian coordinates and some have been frozen.
c     Since there is no redundancy or coupling we just need
c     to make sure that the initial Hessian does not couple
c     frozen with unfrozen variables and we are OK.
c
         do iat = 1, nat_real
            if (.not. oactive(iat)) then
               do i = 1+(iat-1)*3, iat*3
                  do j = 1, nvar
                     dbl_mb(k_hess+j-1+(i-1)*nvar) = 0d0
                     dbl_mb(k_hess+i-1+(j-1)*nvar) = 0d0
                  enddo
                  dbl_mb(k_hess+i-1+(i-1)*nvar) = 1d0
               enddo
            endif
         enddo
      endif
c
      if (.not. rtdb_get(rtdb,'driver:hscale',mt_dbl,1,hscale))
     $     hscale = 1d0
      call dscal(nvar*nvar, hscale, dbl_mb(k_hess), 1)
      if (oprint .and. hscale.ne.1d0) write(6,78) hscale
 78   format(' Scaling initial hessian by ',f6.2)
c
      call geom_hnd_put_data('drv.hess',dbl_mb(k_hess), nvar*nvar)
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('driver_init_hess ma corrupt',0, MA_ERR)
c
      END
      SUBROUTINE HND_OPT_HSS_INIT_ZMT(HSSQ,NINT,ZMAT,IZMAT,SOME,
     $     bscale, ascale, tscale,odebug)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER (ZERO=0.0D+00,      ONE=1.0D+00,
     1     QUARTR=0.25D+00,   PT5=0.5D+00,
     2     THIRD=ONE/3.0D+00, TOHB2=6.423D-02,
     3     TOANGS=0.52917724924D+00)
      PARAMETER (MXATOM=nw_max_atom)
      LOGICAL odebug
      LOGICAL SOME
      LOGICAL BADZ1,BADZ2
      CHARACTER*8 HGUESS
      CHARACTER*8 GUESS
      CHARACTER*8 GUESS1
      CHARACTER*8 GUESS2
      double precision C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION HSSQ(NINT,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION IROW(86),DPARM(6,6),APARM(6,6),RCOV(118)
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
     5     2.03D+00,2.50D+00,2.60D+00,2.75D+00,2.76D+00,3.00D+00/
      DATA DPARM /2.50D-02,0.30D+00,0.38D+00,0.49D+00,0.51D+00,0.81D+00,
     1     0.30D+00,0.68D+00,0.74D+00,0.85D+00,0.68D+00,0.97D+00,
     2     0.38D+00,0.74D+00,1.18D+00,1.02D+00,1.28D+00,0.84D+00,
     3     0.49D+00,0.85D+00,1.02D+00,1.41D+00,0.86D+00,1.14D+00,
     4     0.51D+00,0.68D+00,1.28D+00,0.86D+00,1.62D+00,1.25D+00,
     5     0.81D+00,0.97D+00,0.84D+00,1.14D+00,1.25D+00,1.83D+00/
C     
C     COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
C     AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT
C     HE,NE,AR,KR,RN are from:
C     Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés,
C     Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez
C     (2008). "Covalent radii revisited". Dalton Trans.: 2832–2838.
C     doi:10.1039/b801115j.
C     
      DATA (RCOV(NUCZ),NUCZ=1,2)/0.30D+00,0.46D+00/
      DATA (RCOV(NUCZ),NUCZ=3,10)
     1     /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     2     0.70D+00,0.66D+00,0.58D+00,0.67D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     1     /1.66D+00,1.36D+00,1.25D+00,1.17D+00,
     2     1.10D+00,1.04D+00,0.99D+00,0.96D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     1     /2.03D+00,1.74D+00,
     2     1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     3     1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     4     1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.17D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     1     /2.22D+00,1.92D+00,
     2     1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     3     1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     4     1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,1.31D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     1     /2.35D+00,1.98D+00,
     2     1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D0,
     3     1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D0,
     4     1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     5     1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     6     1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,1.42D+00/
!     These are from: P. Pyykkö, M. Atsumi (2009). "Molecular Single-Bond
!     Covalent Radii for Elements 1-118". Chemistry: A European Journal 15:
!     186–197. doi:10.1002/chem.200800987.
      DATA (RCOV(NUCZ),NUCZ=87,99)
     1    / 2.18D+00, 1.73D+00, 1.53D+00, 1.43D+00, 1.38D+00, 1.34D+00,
     2      1.36D+00, 1.35D+00, 1.35D+00, 1.36D+00, 1.39D+00, 1.40D+00,
     3      1.40D+00/
!     These are from: Beatriz Cordero, Verónica Gómez, Ana E.
!     Platero-Prats, Marc Revés, Jorge Echeverría, Eduard Cremades,
!     Flavia Barragán and Santiago Alvarez (2008). "Covalent radii
!     revisited". Dalton Trans.: 2832–2838. doi:10.1039/b801115j.
      DATA (RCOV(NUCZ),NUCZ=100,118)
     1    /1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00,
     2     1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, 1.29D+00, 1.28D+00,
     3     1.21D+00, 1.22D+00, 1.36D+00, 1.43D+00, 1.62D+00, 1.75D+00,
     4     1.65D+00, 1.57D+00 /
C     
      DATA GUESS,GUESS1,GUESS2 /'GUESS   ','GUESS1  ','GUESS2  '/
C     
C    These are fortran statement functions
C
      DIST(IAT,JAT)=SQRT((C(1,IAT)-C(1,JAT))**2+
     1     (C(2,IAT)-C(2,JAT))**2+
     1     (C(3,IAT)-C(3,JAT))**2) * TOANGS
      BADZ1(NUCZ) = NUCZ.LE.0.OR.NUCZ.GT.86
      BADZ2(NUCZ) = NUCZ.LE.0.OR.NUCZ.GT.118
C     
      IF(odebug) THEN
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
      IF(odebug) THEN
         WRITE(IW,9208) (IZMAT(IZ),IZ=1,NZMAT)
         WRITE(IW,9128) ( ZMAT(IZ),IZ=1,NZVAR)
      ENDIF
      HGUESS=GUESS2
      IF(odebug) THEN
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
            IF(BADZ1(IZA)  .OR.  BADZ1(IZB)) THEN
               hess1 = THIRD
            ELSE
C     
              IROWA = IROW(IZA)
              IROWB = IROW(IZB)
              AIJ = APARM(IROWA,IROWB)
              DIJ = DPARM(IROWA,IROWB)
              RRR = TOANGS*ZMAT(I)
              HII = (AIJ-DIJ)/(RRR-DIJ)
              HESS1 = TOHB2 * HII * HII * HII
            ENDIF
C     
            IF(BADZ2(IZA)  .OR.  BADZ2(IZB)) THEN
              hess2 = THIRD
            ELSE
              RAB = DIST(IZ1,IZ2)
              CAB = RCOV(IZA) + RCOV(IZB)
              HESS2 = 0.3601 * EXP(-1.944*(RAB - CAB))
            ENDIF
*
*     rjh
*
            hess2 = hess2 * bscale
            hess1 = hess1 * bscale
*
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
C     
            HESS1 = QUARTR
            IF(IZB.EQ.1  .OR.  IZC.EQ.1) HESS1 = 0.16D+00
C     
            IF(BADZ2(IZA)  .OR.  BADZ2(IZB)  .OR.  BADZ2(IZC)) THEN
               HESS2 = QUARTR
            ELSE
               RAB = DIST(IZ2,IZ1)
               RAC = DIST(IZ2,IZ3)
               CAB = RCOV(IZA) + RCOV(IZB)
               CAC = RCOV(IZA) + RCOV(IZC)
               HESS2 = 0.089 + 0.11 * (CAB*CAC)**0.42
     1           * EXP(-0.44*(RAB + RAC - CAB - CAC))
            ENDIF
c
*
*     rjh ... trial scaling based on observation that the
*     diagonal seems too large for both angles and torsions
*     Could it be that the original units were wrong?
*     THIS SCALING GIVES A 2X SPEEDUP FOR MANY SYSTEMS!
*
*     generalized to user input
*
            hess1 = hess1 * ascale
            hess2 = hess2 * ascale
*
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
C     
            HESS1 = QUARTR
C     
C     COUNT BONDS TO ATOMS A AND B, EXCEPT FOR A-B ITSELF
C     
            IF(BADZ2(IZA).OR.BADZ2(IZB).OR.BADZ2(IZC).OR.BADZ2(IZD))THEN
               HESS2= QUARTR
            ELSE
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
 135          CONTINUE
C     
C     THE FACTOR 14.0 IN PREEXP SEEMS LARGE, SO BRACKET THE RESULT
C     
              RAB = DIST(IZ2,IZ3)
              CAB = RCOV(IZA) + RCOV(IZB)
              PREEXP = 14.0 * (NBA+NBB)**0.57 / (RAB*CAB)**4
              HESS2 = 0.0015 + PREEXP * EXP(-2.85*(RAB-CAB))
              IF(HESS2.GT.0.30D+00) HESS2 = 0.30D+00
              IF(HESS2.LT.0.05D+00) HESS2 = 0.05D+00
            END IF
*
*     rjh ... trial scaling based on observation that the
*     diagonal seems too large for both angles and torsions
*     Could it be that the original units were wrong?
*     THIS SCALING GIVES A 2X SPEEDUP FOR MANY SYSTEMS!
*
            hess1 = hess1 * tscale
            hess2 = hess2 * tscale
*
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
C     
            HESS1 = QUARTR
C     
            IF(BADZ2(IZA) .OR. BADZ2(IZB) .OR. BADZ2(IZC) .OR.
     $           BADZ2(IZX)) THEN
              HESS2 = QUARTR
            ELSE
              CAX = RCOV(IZA) + RCOV(IZX)
              CAB = RCOV(IZA) + RCOV(IZB)
              CAC = RCOV(IZA) + RCOV(IZC)
              RAX = DIST(IZ1,IZ2)
              PREEXP = 0.0061 * (CAB*CAC)**0.80 * (COS(ZMAT(I)))**4.00
              HESS2 = 0.0025 + PREEXP * EXP(-3.00*(RAX-CAX))
            ENDIF
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
         IF(ITYPE.EQ.7) IZ=IZ+5 ! RJH WAS 7 ... only place 7 appeared?
         IF(odebug) THEN
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
      IF(odebug) THEN
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
#include "mafdecls.fh"
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
 20   call ga_brdcst(323, old_hessian, ma_sizeof(mt_log,1,mt_byte), 0)
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
#include "errquit.fh"
#include "nwc_const.fh"
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
      integer i, j, num, ma_type, nactive_atoms, l_actlist
      logical ignore
      character*80 title
      character*8 source, test
      character*32 theory
      logical geom_zmt_get_source, geom_zmtmak, geom_zmt_get_nizmat,
     $     geom_zmt_get_nzvar, geom_zmt_get_nzfrz, geom_compute_zmatrix,
     $     geom_zmt_get_izfrz, geom_zmt_get_varinfo
      logical opt_geom_cart_coords_get
c     
      call util_print_push
      call util_print_rtdb_load(rtdb, 'driver')
      call ecce_print_module_entry('driver')
      oprint = util_print('information', print_low)
     $     .and. (ga_nodeid() .eq. 0)
      odebug = util_print('debug', print_debug)
c
      if (oprint) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,
     $        'NWChem Geometry Optimization',
     $        40,.true.)
         write(6,*)
         write(6,*)
      endif
c
      if (rtdb_cget(rtdb,'title',1,title)) then
         if (oprint) then
            write(6,*)
            write(6,*)
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
      endif
c     
c     Before loading info regenerate the autoz coordinates if
c     the user has requested this
c     
      if (rtdb_get(rtdb,'driver:redoautoz',mt_log,1,ignore)) then
         if (ga_nodeid() .eq. 0) then
            ignore = rtdb_parallel(.false.)
            ignore = rtdb_delete(rtdb,'driver:redoautoz')
            call driver_opt_del_hess() ! Clear hessian info
            if (.not. geom_create(geom, 'geometry'))
     &           call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
            if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &           call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
            if (.not.geom_zmt_get_source(geom,source))
     $           call errquit('hnd_optx: zmat source ',geom, GEOM_ERR)
c     
            if (source .eq. 'autoz') then
               if (oprint) write(6,777)
 777           format(/' Regenerating redundant internal coordinates'/)
               if (.not. geom_zmtmak(rtdb, geom, oprint))
     $              call errquit('driver: regeneration of autoz failed',
     $              0, GEOM_ERR)
               if (oprint) write(6,*)
               if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &              call errquit('hnd_opt: failed storing geometry',0,
     &       RTDB_ERR)
            endif
            if (.not. geom_destroy(geom))
     $           call errquit('driver: geom_destroy?',0, GEOM_ERR)
            ignore = rtdb_parallel(.true.)
         endif
         call ga_sync()
      endif
c     
c     ----- parameters for optimization driver -----
c     
      if (.not. rtdb_get(rtdb,'driver:modsad',mt_int,1,modsad)) 
     $     modsad=0
      if (.not. rtdb_get(rtdb,'driver:trust',mt_dbl,1,trust)) then
         if (modsad .eq. 0) then
            trust = 0.3d0
         else
            trust = 0.1d0
         endif
      endif
      if (.not. rtdb_cget(rtdb,'driver:xyz',1,xyz))
     $     xyz = ' '
      if (.not. rtdb_get(rtdb,'driver:eprec',mt_dbl,1,eprec)) then
         if (.not. rtdb_cget(rtdb,'task:theory',1,theory))
     $        theory = ' '
         if (inp_compare(.false.,theory,'dft')) then
            eprec = 5e-6
         else
            eprec = 1e-7
         endif
      endif
      if (.not. rtdb_get(rtdb,'driver:gmax_tol',mt_dbl,1,gmax_tol))
     $     gmax_tol = 0.00045d0
      if (.not. rtdb_get(rtdb,'driver:grms_tol',mt_dbl,1,grms_tol))
     $     grms_tol = 0.0003d0
      if (.not. rtdb_get(rtdb,'driver:xmax_tol',mt_dbl,1,xmax_tol))
     $     xmax_tol = 0.0018d0
      if (.not. rtdb_get(rtdb,'driver:xrms_tol',mt_dbl,1,xrms_tol))
     $     xrms_tol = 0.0012d0
      if (.not. rtdb_get(rtdb,'driver:nptopt',mt_int,1,nptopt))
     $     nptopt=20
      if (.not. rtdb_get(rtdb,'driver:inhess',mt_int,1,inhess))
     $     inhess=0
      if (.not. rtdb_get(rtdb,'driver:linopt',mt_int,1,linopt)) 
     $     linopt=1
      if (.not. rtdb_get(rtdb,'driver:moddir',mt_int,1,moddir)) 
     $     moddir=0
      if (.not. rtdb_get(rtdb,'driver:vardir',mt_int,1,vardir)) 
     $     vardir = 0
      ignore = rtdb_delete(rtdb,'driver:vardir') ! Forget for restart
      if (.not. rtdb_get(rtdb,'driver:modsad',mt_int,1,modsad)) 
     $     modsad=0
      if (.not. rtdb_get(rtdb,'driver:sadstp',mt_dbl,1,sadstp)) 
     $     sadstp=0.1d0
      if (.not. rtdb_get(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $     oqstep = .true.
      if (.not. rtdb_get(rtdb,'driver:ofirstneg',mt_log,1,ofirstneg))
     $     ofirstneg = .true.
      if (.not. rtdb_get(rtdb,'driver:opcg',mt_log,1,opcg))
     $     opcg = .false.
      if (.not. rtdb_get(rtdb,'driver:modupd',mt_int,1,modupd)) then
         if (modsad .eq. 0) then
            modupd = 1          ! BFGS update for minimization
         else
            modupd = 2          ! PSB update for saddle point
         endif
      endif
      if (.not. rtdb_get(rtdb,'driver:ocheckgrad',mt_log,1,ocheckgrad))
     $     ocheckgrad = .false.

      if (.not. rtdb_get(rtdb,'includestress',mt_log,1,ostress)) then
         ostress = .false.
      end if
      if (.not. rtdb_get(rtdb,'includelattice',mt_log,1,ostress2)) then
         ostress2 = .false.
      end if
      if ((.not.ostress).and.(ostress2)) ostress2 = .false.
      if ((ostress)     .and.(ostress2)) ostress  = .false.

      if (.not.rtdb_get(rtdb,'driver:diagonal_hessian',
     $                   mt_log,1,diagonal_hessian))
     $   diagonal_hessian = .false.
      
c
c     Force sensible options
c
      if (modsad .eq. 0) then
         modupd = 1             ! BFGS update for minimization
      else
         linopt = 0             ! No line search for saddle
      endif
c     
c     Save a  copy of the initial geometry so we can analyze what 
c     happened during the optimization
c     
      if (ga_nodeid() .eq. 0) then
         ignore = rtdb_parallel(.false.)
         if (.not. geom_create(geom, 'geometry'))
     &        call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &        call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
         if (.not. geom_rtdb_store(rtdb, geom, 'driverinitial'))
     &        call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
         if (.not. geom_destroy(geom))
     $        call errquit('driver: geom_destroy?',0, GEOM_ERR)
         ignore = rtdb_parallel(.true.)
      endif
      call ga_sync()
c     
c     Load the geometry info
c     
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
      if (.not. geom_ncent(geom,nat))
     $     call errquit('hnd_opt: natoms?',nat, GEOM_ERR)
      call grad_active_atoms(rtdb, nat, oactive, nactive)
      if (.not. geom_systype_get(geom, isystype))
     $     call errquit('driver: systype?',0, GEOM_ERR)
c     
c     Zmatrix info
c     
      if (.not.geom_zmt_get_source(geom,source))
     $     call errquit('hnd_optx: zmat source ',geom, GEOM_ERR)
      zcoord  = (source .ne. ' ') ! Any kind of zmatrix
      zmatrix = (source .eq. 'user') ! zmatrix from user
      nzvar = 0
      nzfrz = 0
      if (zcoord) then
         if(.not.geom_zmt_get_nizmat(geom,nzmat))
     $        call errquit('geom_input: geom_zmt_get_nizmat failed',0,
     &       GEOM_ERR)
         if(.not.geom_zmt_get_nzvar(geom,nzvar))
     $        call errquit('geom_input: geom_zmt_get_nzvar failed',0,
     &       GEOM_ERR)
         call geom_bandbi(geom)
         if (.not. zmatrix) oqstep = .false. ! Only qstep in user z-matrix
         if(.not.geom_zmt_get_nzfrz(geom,nzfrz))
     $        call errquit('geom_input: geom_zmt_get_nzfrz failed',0,
     &       GEOM_ERR)
         if(.not.geom_zmt_get_izfrz(geom,izfrz))
     $        call errquit('geom_input: geom_zmt_get_izfrz  failed',0,
     &       GEOM_ERR)
      endif
      if (oprint) then
         write(6,1) gmax_tol, grms_tol, xmax_tol, xrms_tol, trust,
     $        sadstp, eprec,
     $        nptopt, inhess, linopt, modupd, modsad, moddir, 
     $        vardir, ofirstneg, opcg, source
 1       format(
     $        ' maximum gradient threshold         (gmax) = ', f10.6,/,
     $        ' rms gradient threshold             (grms) = ', f10.6,/,
     $        ' maximum cartesian step threshold   (xmax) = ', f10.6,/,
     $        ' rms cartesian step threshold       (xrms) = ', f10.6,/,
     $        ' fixed trust radius                (trust) = ', f10.6,/,
     $        ' maximum step size to saddle      (sadstp) = ', f10.6,/,
     $        ' energy precision                  (eprec) = ', 1p,d9.1,
     $        0p,/,
     $        ' maximum number of steps          (nptopt) = ', i4,/,
     $        ' initial hessian option           (inhess) = ', i4,/,
     $        ' line search option               (linopt) = ', i4,/,
     $        ' hessian update option            (modupd) = ', i4,/,
     $        ' saddle point option              (modsad) = ', i4,/
     $        ' initial eigen-mode to follow     (moddir) = ', i4,/,
     $        ' initial variable to follow       (vardir) = ', i4,/,
     $        ' follow first negative mode     (firstneg) = ', l4,/,
     $        ' apply conjugacy                    (opcg) = ', l4,/,
     $        ' source of zmatrix                         =   ', a8,/)
         if (modsad .eq. 0) then
            write(6,9994)
 9994       format(/,10x,19('-'),
     1           /,10x,'Energy Minimization',
     2           /,10x,19('-'),/)
         else
            write(6,9995)
 9995       format(/,10x,23('-'),
     1           /,10x,'Transition State Search',
     2           /,10x,23('-'),/)
         endif
         if (ostress) then
            write(6,*) ' INCLUDING STRESS !!!!!!!!!!!!!!!!'
            if (isystype.eq.0) call errquit('NOT A PERIODIC SYSTEM',0,
     &       GEOM_ERR)
         endif
         if (ostress2) then
            write(6,*) ' INCLUDING LATTICE GRADIENTS !!!!!'
            if (isystype.eq.0) call errquit('NOT A PERIODIC SYSTEM',0,
     &       GEOM_ERR)
         endif

         if (diagonal_hessian) then
            write(6,*) ' USING A DIAGONAL HESSIAN !!!!!!!!!!!!!!!!'
         end if

         call util_flush(6)
      endif
c     
c     Nvar is the no. of variables in the optimization
c
c     If we are optimizing the unit cell parameters then we pretend
c     there there are 3 more atoms which will parameterize the
c     unit cell.
c     
      nat_real = nat
      if (ostress)  nat = nat + 3
      if (ostress2) nat = nat + 2
      ncart = 3*nat
      if (zcoord) then
         nvar = nzvar
      else
         nvar = ncart
      endif
c     
c     ----- frozen coordinates ? if so .... -----
c     
      if (oprint) then
         if (zcoord) then
            if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $           call errquit('hnd_optx: varinfo?',0, GEOM_ERR)
            write(6,*)
            write(6,*) 'Names of Z-matrix variables '
            do i = 1, nzvar, 5
               write(6,5432) (j,zvarname(j),j=i,min(nzvar,i+4))
 5432          format(1x,5(i4,1x,a8,2x))
            enddo
            write(6,*)
            write(6,*) 'Variables with the same non-blank name',
     $           ' are constrained to be equal'
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
         if(nzfrz.gt.0) then
            write(6,3) nzfrz
 3          format(/' Number of frozen variables =',i5)
            write(6,33) (izfrz(i),i=1,nzfrz) 
 33         format( '           Frozen variables =',8i5)
         endif
      endif
c
      if (zcoord) then
c     have already called bandbi but must ensure no cartesian
c     constraints have been imposed since then cannot use internals
         if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type, 
     $        nactive_atoms, l_actlist) .or.
     $        rtdb_ma_get(rtdb, 'geometry:inactlist', ma_type,
     $        nactive_atoms, l_actlist)
     $       ) then
            if (ga_nodeid() .eq. 0) then
               write(6,666)
 666           format(/,/,
     $              '!!! You are using internal coordinates for the ',
     $              'optimization but applying',/,
     $              '!!! constraints on the Cartesian coordinates by',
     $              ' only computing gradients',/,
     $              '!!! for a subset of the atoms (via geometry:act',
     $              'list).  This does not ',/,
     $              '!!! currently work.  Either,',/, '!!!',/,
     $              '!!!   1) Retain the Cartesian constraints and p',
     $              'erform the optimization ',/,
     $              '!!!      entirely in Cartesians by specifying t',
     $              'he coordinates in ',/,
     $              '!!!      cartesians and using the NOAUTOZ keywo',
     $              'rd on the GEOMETRY directive',/,
     $              '!!!',/,'!!! or',/, '!!!',/,
     $              '!!!   2) Remove the Cartesian constraints and p',
     $              'erform the optimization',/,
     $              '!!!      in redundant internal coordinates.  To',
     $              ' remove the Cartesian',/,
     $              '!!!      constraints remove the "set geometry:a',
     $              'ctlist ..." directive from',/,
     $              '!!!      startup jobs, or use "unset geometry:a',
     $              'ctlist" in restart jobs.',/)
               call util_flush(6)
               call errquit('driver: cartesian constraints+internals',0,
     &       GEOM_ERR)
            endif
            call ga_sync()
         endif
      else
         call driver_cart_pmat(rtdb, geom)
      endif
c
c     We'll be printing the zmatrix and gradient so no
c     need for geom print to do this
c
      if (zcoord) call util_print_disable('geomzmat')
c
      lsmode = ' '
      energy = 0d0
      energyp= 0d0
      energyp2=0d0
      energyp3=0d0
      alphap = 1d0
      alpha  = 1d0
      gmax   = 0d0
      grms   = 0d0
      smax   = 0d0
      srms   = 0d0
      xmax   = 0d0
      xrms   = 0d0
      call dfill(max_nvar, 0d0, ds, 1)
      call dfill(max_nvar, 0d0,dsp, 1)
      call dfill(max_nvar, 0d0, gx, 1)
      call dfill(max_nvar, 0d0, gq, 1)
      call dfill(max_nvar, 0d0,  g, 1)
      call dfill(max_nvar, 0d0, gp, 1)
c
      if (zcoord) then
         if (.not. geom_compute_zmatrix(geom, sp))
     $        call errquit('driver: zmatrix?',0, GEOM_ERR)
      else
c     FRACTIONAL?
         if (.not. opt_geom_cart_coords_get(geom, sp))
     $        call errquit('driver: geom?',0, GEOM_ERR)
      endif
c     
      end
      subroutine driver_hessian_update(geom,rtdb)
      implicit none

#include "errquit.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "util.fh"
#include "mafdecls.fh"
c     
c     Update the current Hessian in the optimization variables using 
c     .   gp() - the gradient at the previous point
c     .    g() - the gradient at the current point
c     .   ds() - the previous search direction
c     .  alpha - the step in the previous search direction
c     
c     Only the Hessian is modified.
c     
      integer geom, rtdb
      double precision hds(max_nvar)
      double precision dsds, dshds, dsdg
      integer l_hess, k_hess, i, j, l_dstmp, k_dstmp
      integer ind, iend, iat
      integer istrss
      double precision strss1,strss2,dum1,dum2,dum3,dum4
      logical redo_hessian

! This is a statement function
      ind(i,j) = k_hess + i + (j-1)*nvar - 1
c save ds
      if (.not. ma_push_get(mt_dbl, nvar, 'dstmp',
     $     l_dstmp, k_dstmp)) call errquit
     $     ('driver_hessian_update: memory for ds',nvar,
     &       GEOM_ERR)
      call dcopy(nvar, ds, 1, dbl_mb(k_dstmp), 1)
c
      if (alpha .eq. 0d0) call errquit
     $     ('driver_hessian_update: zero step?',0, GEOM_ERR)
      call dscal(nvar, alpha, ds, 1)
C
C    For fractional coordinates, this code assumes that
C    atoms do not wrap around. So we unwrap step
C    This is why we save ds, and then restore it
C
      if (isystype .ne. 0 .and. .not. zcoord) then
        iend = nvar/3
        if (ostress) iend = (nvar-9)/3
        if (ostress2)iend = (nvar-6)/3
        do i = 1, iend
          if(isystype .ge. 1) then
            if ((ds((3*i)-2)) .gt. 0.5d0) ds((3*i)-2)=ds((3*i)-2)-1.0d0 
            if ((ds((3*i)-2)) .lt.-0.5d0) ds((3*i)-2)=ds((3*i)-2)+1.0d0
          endif
          if(isystype .ge. 2) then
            if ((ds((3*i)-1)) .gt. 0.5d0) ds((3*i)-1)=ds((3*i)-1)-1.0d0
            if ((ds((3*i)-1)) .lt.-0.5d0) ds((3*i)-1)=ds((3*i)-1)+1.0d0
          endif
          if(isystype .ge. 3) then
            if ((ds((3*i)-0)) .gt. 0.5d0) ds((3*i)-0)=ds((3*i)-0)-1.0d0
            if ((ds((3*i)-0)) .lt.-0.5d0) ds((3*i)-0)=ds((3*i)-0)+1.0d0
          endif
        enddo
      endif
c     
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('driver_hessian_update: memory for hessian',nvar**2,
     &       GEOM_ERR)
      call geom_hnd_get_data('drv.hess',dbl_mb(k_hess), nvar**2)
c
c     Form bits and pieces that are needed
c     
      call dgemv('n',nvar,nvar,1d0,dbl_mb(k_hess),nvar,
     $     ds,1,0d0,hds,1)
c
      dshds = ddot(nvar, ds, 1, hds, 1)
      dsds  = ddot(nvar, ds, 1,  ds, 1)
      dsdg  = 0d0
      do i = 1, nvar
         dsdg = dsdg + ds(i)*(g(i) - gp(i))
      enddo
c     
c        
c     ---- Limit Lattice-Atom cross-terms ----
c              
      if (isystype .ne. 0 .and. ostress) then
        istrss = nvar-9
        strss1 = 0.25d0 ! Lattice: cross-terms
        strss2 = 0.95d0 ! Forgetfulness
      elseif(isystype .ne. 0 .and. ostress2) then
        istrss = nvar-6
        strss1 = 0.25d0
        strss2 = 0.95d0
      else
        istrss = nvar+1 ! ie. never
        strss1 = 1.0d0
        strss2 = 1.0d0
      endif  

      if(abs(dsdg).gt.1d-8) then
         do i=1,nvar
            if (i.gt.istrss) then
              dum1 = strss1
              dum2 = strss2
            else
              dum1 = 1.0d0
              dum2 = 1.0d0
            endif
            do j=1,nvar
               if (j.gt.istrss) then
                 dum3 = strss1*dum1
                 dum4 = strss2*dum2
               else 
                 dum3 = dum1
                 dum4 = dum2
               endif
               if(modupd.le.1) then
c     ----- -bfgs- update -----
                  dbl_mb(ind(i,j))=dbl_mb(ind(i,j))*dum4
     $                 + ((g(i)-gp(i))*(g(j)-gp(j))/dsdg
     1                 - hds(i)* hds(j)/dshds)*dum3
               else
c     ----- -psb- update -----
                  dbl_mb(ind(i,j))=(dbl_mb(ind(i,j)))*dum4
     $                 + (((g(i)-gp(i))-hds(i))*ds(j)/dsds
     1                 + ((g(j)-gp(j))-hds(j))*ds(i)/dsds
     2                 - ds(i)*ds(j)*(dsdg-dshds)/(dsds*dsds))*dum3
               endif
            enddo   
         enddo   
      else
c     Last step very very tiny - clean up lattice terms, since we are probably stuck
         do i=1,istrss
           do j=istrss+1,nvar
              dbl_mb(ind(i,j))=0
              dbl_mb(ind(j,i))=0
           enddo
         enddo
         do i=istrss+1,nvar
           do j=i+1,nvar
              dbl_mb(ind(i,j))=0
              dbl_mb(ind(j,i))=0
           enddo
         enddo
         do i=istrss+1,nvar
           dbl_mb(ind(i,i))=max(dbl_mb(ind(i,i)),5.0d0)
         enddo
      endif
c
      if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
c     
c     We are in cartesian coordinates and some have been frozen.
c     Since there is no redundancy or coupling we just need
c     to make sure that the Hessian does not couple
c     frozen with unfrozen variables and we are OK.
c
         do iat = 1, nat_real
            if (.not. oactive(iat)) then
               do i = 1+(iat-1)*3, iat*3
                  do j = 1, nvar
                     dbl_mb(k_hess+j-1+(i-1)*nvar) = 0d0
                     dbl_mb(k_hess+i-1+(j-1)*nvar) = 0d0
                  enddo
                  dbl_mb(k_hess+i-1+(i-1)*nvar) = 1d0
               enddo
            endif
         enddo
      endif


*     ************************************************
*     *** reinitialize the hessian all the time    ***
*     *** This option is used to try to avoid      ***
*     *** hessians that update poorly and thus     ***
*     *** avoid (higher-order) saddle points       ***
*     ************************************************
      redo_hessian = .false.
      if (redo_hessian) then
        call driver_hss_reinit(rtdb,geom,l_hess, k_hess)
      endif

*     ************************************************
*     *** making the hessian diagonal all the time ***
*     *** This option is used to try to avoid      ***
*     *** (higher-order) saddle points             ***
*     ************************************************
      if (diagonal_hessian) then
         do i=1,nvar
            do j=1,nvar
               if (i.ne.j) then
                  dbl_mb(ind(i,j))  = 0.0d0  ! Diagonal only
               else
                  dbl_mb(ind(i,j)) = abs(dbl_mb(ind(i,j)))
               end if
            end do
         end do
      end if

c
      call geom_hnd_put_data('drv.hess',dbl_mb(k_hess), nvar**2)
      if (.not. ma_pop_stack(l_hess)) call errquit
     $     ('driver_hessian_update: ma?',0, MA_ERR)
c     restore ds
      call dcopy(nvar,dbl_mb(k_dstmp), 1, ds, 1)
      if (.not. ma_pop_stack(l_dstmp)) call errquit
     $     ('driver_hessian_update: ma?',0, MA_ERR)
c     
      end
      subroutine driver_opt_search_dir(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
      integer geom
c
c     This routine for minimization
c
c     Put into ds() a search direction in the optimization
c     variables (internal or cartesian) based upon the
c     current gradient, g(), and Hessian.  Apply constraints.
c
c     Apply step restrictions by recommending an initial
c     value for the line search parameter alpha.
c
c     Only alpha and ds() are modified.
c
      integer i, iat,j, top
*     integer info
      integer l_hess, k_hess, l_work, k_work, lenwork
      double precision  e(max_nvar) ! Hessian eigenvalues
      double precision pg(max_nvar) ! P.g
      double precision gv(max_nvar) ! Gradient along eigenvectors
      double precision dv(max_nvar) ! Step along eigenvectors
      double precision dsmax    ! Max. value of current step (smax is prev.)
c
      double precision beta, s0g0, s0g1, s1g0, s1g1, numerator,
     $     denominator
      double precision bohr, deg ! For printing purposes
      double precision trustds  ! Restriction of step in opt. variable
      double precision lattice(6), scaler(3) ! periodic scaling
      double precision dum1,dum2,dum3
      double precision smalleig
      parameter (smalleig = 1.0d-8)
      logical geom_print_zmatrix
      external geom_print_zmatrix
      logical ophigh
      logical geom_lattice_get
      external geom_lattice_get
c
      ophigh = util_print('high', print_high)
c     
c     Enforce frozen atoms in cartesians
c        
      if (ga_nodeid().eq.0.and.ophigh)
     $     write(luout,*) 'Zeroing constrained gradient'

      if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
         do iat = 1, nat_real
            if (.not. oactive(iat)) then
               do i = 1, 3
c                  write(6,*) ' Zeroing constrained cartesian', 
c     $                 iat, i, dv((iat-1)*3+i)
                  g((iat-1)*3+i) = 0.0
               end do
            end if
         end do
      end if
c
c     Enforce symmetry on lattices - this is only done here just to make sure
c     that the fixed angles have zero gradients, because otherwise the
c     geometry optimization will never converge.
c
      if (ostress2) call sym_lattice_stress2(geom,g(3*nat_real+1))
c
c     Get the Hessian and gradient with appropriate projectors
c     applied following Peng, Ayala, Schlegel and Frisch so that
c     redundant internal modes are shifted to high eigenvalues.
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('driver_search_dir: memory for hessian',nvar**2, MA_ERR)
      call driver_project_hess_grad(dbl_mb(k_hess), pg)
      if (ostress2) call sym_lattice_stress2(geom,pg(3*nat_real+1))
c
c     Diagonalize the Hessian.  Should really do the generalized
c     eigenvalue problem since the underlying basis is not independent
c     (if we are using autoz). NOT YET BEING DONE.  
c
c     To cause degenerate eigenvalues to be resolved into symmetry
c     adapted combinations use Jacobi not dsyev and screen out junk
c
      lenwork = max(nvar**2,100)
      if (.not. ma_push_get(mt_dbl, lenwork, 'work',
     $     l_work, k_work)) call errquit
     $     ('driver_search_dir: memory for hessian', lenwork, MA_ERR)
      do i = 0, nvar**2-1
         if (abs(dbl_mb(k_hess+i)).lt.1d-8) dbl_mb(k_hess+i) = 0d0
      enddo
*      call dsyev('v', 'u', nvar, dbl_mb(k_hess), 
*     $     nvar, e, dbl_mb(k_work), lenwork, info)
*      if (info .ne. 0) call errquit
*     $     ('driver_search_dir: dsyev failed',info, 0)
      call util_jacobi(nvar, dbl_mb(k_hess), nvar, e)
      if (odebug .or. (util_print('hvecs',print_never) 
     $     .and. ga_nodeid().eq.0)) then
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
         if (e(i) .gt. smalleig) then
            dv(i) = -gv(i) / abs(e(i)) ! +ve curvature
         else if (e(i) .lt. -smalleig) then
            dum1 = -gv(i) / abs(e(i)) 
            dum2 = sign(0.3d0*trust,-gv(i)) ! upper limit
            dum3 = sign(0.03d0*trust,-gv(i)) ! lower limit
            if (abs(dum1) .gt. abs(dum2)) then
              dv(i) = dum2
              if (oprint) write(6,664) i, e(i), gv(i), dv(i)
 664             format(' Limiting step in negative mode',i5,
     $           '  eval=', 1p,d8.1,' grad=',d8.1,' step=',d8.1)
            else if (abs(dum1) .lt. abs(dum3)) then
              dv(i) = dum3
              if (oprint) write(6,665) i, e(i), gv(i), dv(i)
 665             format(' Forcing step in negative mode',i5,
     $           '  eval=', 1p,d8.1,' grad=',d8.1,' step=',d8.1)
            else
              dv(i) = dum1
              if (oprint) write(6,666) i, e(i), gv(i), dv(i)
 666             format(' Taking step in negative mode',i5,
     $           '  eval=', 1p,d8.1,' grad=',d8.1,' step=',d8.1)
            endif
         else
c
c     Linearly dependent modes, and rotations and translations
c     have been shifted up.  Zero modes here are either accidental,
c     or dependencies that we don't know about.  In this case just
c     use steepest descent with step restriction.
c
c     This would be more sensible in normalized coordinates.
c
            dv(i) = sign(min(trust,abs(gv(i))),-gv(i))
            if (oprint) write(6,1) i, e(i), dv(i)
 1          format(' Using steepest descent step in mode',i5,
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
c     Uncomment following for conjugacy.  Dubious benefit with quasi-Newton.
c
c     Conjugacy ... sometimes helps the end game.  However, if the
c     last line search was very incomplete, or we are not doing line
c     searches then we should not do this.
c
      beta = 0d0
      if (linopt.ne.0  .and.  lsmode.eq.'ok' .and. opcg) then
         s0g1 = ddot(nvar, dsp, 1, g, 1)
         s0g0 = ddot(nvar, dsp, 1,gp, 1)
         s1g1 = ddot(nvar, ds, 1, g, 1)
         s1g0 = ddot(nvar, ds, 1,gp, 1)
         numerator = s1g1 - s1g0
         denominator = s0g1 - s0g0
         if (abs(denominator) .gt. 1d-8) then
            beta = -numerator/denominator
            if (beta .lt. 0d0) then 
               beta = 0d0
            endif
         endif
         if (oprint .and. beta.ne.0d0) write(6,88) beta
 88      format(' Applying conjugacy beta = ',1p,d8.1)
      endif
      call daxpy(nvar, beta, dsp, 1, ds, 1)
c
c     Enforce symmetry
c
      call driver_symmetrize_step(geom)
c
c     Enforce frozen atoms in cartesians
c

      if (ga_nodeid().eq.0.and.ophigh) 
     $     write(luout,*) 'Zeroing constrained gradient'

      if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
         do iat = 1, nat_real
            if (.not. oactive(iat)) then
               do i = 1, 3
c                  write(6,*) ' Zeroing constrained cartesian', 
c     $                 iat, i, ds((iat-1)*3+i)
                  ds((iat-1)*3+i) = 0.0
               end do
            end if
         end do
      end if
c
c     Restriction of the overall step based on
c     
c     1.  Did we go uphill on the last step?
c     2.  Is a component of the step large in the optimization space.
c
      dsmax = 0.0d0
      if(isystype.gt.0 .and. (.not. zcoord)) then ! Fractional coordinates
        if (.not. geom_lattice_get(geom, lattice))
     $     call errquit('driver: failed to get lattice',0,0)

        do i = 1, isystype
           scaler(i) = lattice(i)
        enddo
        do i = isystype+1,3
           scaler(i) = 1.0d0
        enddo
        j=0
        top = nvar
        if(ostress) top = nvar - 9
        if(ostress2)top = nvar - 6
        do i = 1, top
           j=j+1
           if (j.gt.3) j=1
           dsmax = max(dsmax, abs(ds(i))*scaler(j))
        enddo
!       Note that the linear optimizer can increase the step size by 4X
        if (ostress) then
           do i = top+1,nvar
              j = ((i-top)+2)/(3) ! 1,1,1,2,2,2,3,3,3
              dsmax = max(dsmax, abs(ds(i)/scaler(j))*24)!1.25% max, if trust=.3
              dsmax = max(dsmax, abs(ds(i)))
           enddo
           do i = top+1,nvar,3
              j = ((i-top)+2)/(3) ! 1,1,1,2,2,2,3,3,3
              dum1 = abs(ds(i))+abs(ds(i+1))+abs(ds(i+2)) ! Total change in a vector
              dsmax = max(dsmax, abs(dum1/scaler(j))*18)!1.75% max, if trust=.3
              dsmax = max(dsmax, abs(dum1))
           enddo
        endif
        if (ostress2) then
           do i = top+1,top+3! Lengths
              j=i-top
              dsmax = max(dsmax, abs(ds(i)/scaler(j))*12)!2.5% max, if trust=.3
           enddo
           do i = top+4,nvar ! Angles
              dsmax = max(dsmax, abs(ds(i))*6.87548d0)!2.5 degrees max if trust=.3
           enddo
        endif
      else
        do i = 1, nvar
          dsmax = max(dsmax, abs(ds(i)))
        enddo
      endif
c
      trustds =trust
c
c     If we are not applying conjugacy and only doing approximate LS
c     then the LS could be exploring a similar mode to the last step, 
c     thus use alphap as a guide.  sqrt(alphap) converges to 1.0 when iterated.
c     If we are using conjugacy then our seach step will be orthogonal
c     to the last one so we cannot use alphap.  Exact NR is automatically
c     conjugate (and will also have unit step).
c
c     But ... after many tests conclude that obtain slightly better speed
c     and overall stability by taking initial NR step of alpha=1
c
c$$$      if (beta .ne. 0d0) then
c$$$         alpha = 1d0
c$$$      else
c$$$         alpha = sqrt(abs(alphap))
c$$$         alpha = min(2.5d0,alpha)
c$$$         alpha = max(0.4d0,alpha)
c$$$      endif
c
      alpha = 1d0               !!!!!!!!!!!!!!!!!!1
c
      if (dsmax.gt.trustds .or. (energy-energyp).gt.eprec)
     $     alpha = min(1d0,alpha)
      if (dsmax .gt. trustds) then
         alpha = trustds/dsmax
         if (oprint) write(6,7) alpha
 7       format(' Restricting overall step due to large component.',
     $        ' alpha=', f6.2)
      endif
      if ((energy-energyp) .gt. eprec) then
         alpha = min(0.5d0,0.5d0*alpha)
         if ((energyp-energyp2) .gt. eprec) then
c          two bad steps in a row, scale alpha
           alpha = min(0.3d0,0.5d0*alpha)
           if ((energyp2-energyp3) .gt. eprec) then
c            3 bad steps in a row, scale alpha
             alpha = min(0.2d0,0.6d0*alpha)
           endif
         endif
         if (oprint) write(6,8) alpha
      else if ((energyp-energyp2) .gt. eprec) then
         alpha = min(0.9d0,0.9d0*alpha) ! Did ok this step, but not last step
         if ((energyp2-energyp3) .gt. eprec) then
c          Did OK this step, but not last two steps
           alpha = min(0.7d0,0.8d0*alpha)
         endif
         if (oprint) write(6,8) alpha
      else if ((energy-energyp3) .gt. eprec) then
         alpha = min(0.9d0,0.9d0*alpha) ! Higher than 3 steps ago, but this and last step are OK
         if (oprint) write(6,8) alpha
      endif
 8       format(' Restricting overall step due to uphill motion.  ',
     $        ' alpha=', f6.2)
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('driver_search_dir: ma?',0, MA_ERR)
c
c     Edo seems to have encountered a case where different processors
c     generated different steps.  To prevent this, broadcast the 
c     critical info to everyone.
c
      call ga_brdcst(1,ds,8*nvar,0)
      call ga_brdcst(2,alpha,8,0)
c
      if (util_print('searchdir',print_high) .and. 
     $     ga_nodeid().eq.0) then
         write(6,*)
         write(6,*) '       The search direction'
         if (zcoord) then
            bohr = 0.52917715d0
            deg  = 0.52917715d0*180d0/(4d0*atan(1d0))
            call geom_zmat_ico_scale(geom, ds, bohr, deg)
            if (.not. geom_print_zmatrix(geom, ds, 'Search  ',.true.))
     $           call errquit('driver:print search_dir?',geom, GEOM_ERR)
            call geom_zmat_ico_scale(geom, ds, 1d0/bohr, 1d0/deg)
         else
            call output(ds,1,3,1,nat,3,nat,1)
         endif
         write(6,*)
         call util_flush(6)
      endif
c
      end
      subroutine driver_line_search1(rtdb,geom)
      implicit none
#include "nwc_const.fh"
#include "coptopt.fh"
#include "geom.fh"
#include "util.fh"
      integer rtdb, geom
c
c     Starting with a value of alpha recommended by whatever
c     generated the search direction, perform an approximate line 
c     search.
c
c     Return with alpha, geom as the accepted point.
c
c     Only alpha and geom are modified.
c
      double precision e0, e1, e2p, dsgrad
      double precision hess, a0, a1, a2
      double precision driver_energy_step
c
      dsgrad = ddot(nvar, ds, 1, g, 1)
      if (dsgrad*alpha .ge. 0d0) then
         if (oprint) write(6,*) 
     $        ' Warning ... line search gradient +ve ', alpha, dsgrad
         alpha = -alpha
      endif
c
      e0 = energy
      a0 = 0d0
c
c     Compute energy at alpha*ds
c
      e1 = driver_energy_step(rtdb, geom)
      a1 = alpha
c
c     Fit to a parabola using e0, g0, e1
c
      hess = (e1-e0-alpha*dsgrad)/alpha**2
c
      if (abs(e1-e0) .lt. eprec) then
c
c     Insufficient precision in the energy to do anything else
c
c     If this condition persists we should switch off the LS.
c
         a2 = alpha 
         lsmode = 'fixed'
      else if (hess .gt. 0d0) then
c
c     Positive curvature.  Compute NR step, and ...
c     1) If we went downhill take the step with some restriction
c     2) If we went uphill then take the step if it is in [0,alpha]
c     .  otherwise take alpha/2.
c
         a2 = -dsgrad/(2d0*hess)
         if ((e1 - e0) .le. -eprec) then ! Downhill
            lsmode = 'downhill'
            if (abs(a2) .gt. 4d0*abs(alpha)) then 
               lsmode = 'restrict'
               a2 = sign(4d0*abs(alpha),a2)
            endif
         else                   ! Uphill
            lsmode = 'bracket'
            if (abs(a2) .gt. abs(alpha)) then
               lsmode = 'confused'
               a2 = alpha/2d0
            endif
         endif
      else
c
c     Negative curvature.
c     1.  If we went downhill, step cautiously ahead.
c     2.  If we went uphill then something is really confused.
c     .   Just accept the current point.
c
         if ((e1 - e0) .lt. eprec ) then ! Downhill
            lsmode = 'negative'
            a2 = 2d0*alpha
         else
            lsmode = 'punt'
            a2 = alpha
         endif
      endif
c
c     Now have a2.  If it is close to alpha or does
c     not give us much energy lowering then just
c     return alpha so that any SCF/MCSCF/DFT does not have to
c     be repeated before computing the next gradient.
c
c     Be sure to have reduced line search gradient by 10x
c     and that we are working within available precision.
c
c     Added accept of alpha if predicted energy lowering is
c     less than energy precision.  Should diminish thrashing
c     around on flat surfaces.
c
      e2p = e0 + dsgrad*a2 + hess*a2**2
***     $     ((e2p-e1) .gt. 0.01*(e1-e0)) .or.
***     $     (abs(e2p-e1) .lt. eprec) )then
      if ( (abs(a2-alpha) .lt. abs(0.1d0*alpha)) ) then
         lsmode = 'accept'
         a2 = alpha
      endif
      e2p = e0 + dsgrad*a2 + hess*a2**2
c
      if (oprint) write(6,1) alpha, dsgrad, hess, e1, lsmode, a2, e2p
 1    format(' Line search: ',/, '     step=',f5.2, ' grad=',1p,d8.1,0p,
     $     ' hess=',1p,d8.1,0p,' energy=',f14.6,' mode=',a,/,
     $     ' new step=',f5.2, 19x,'predicted energy=',f14.6)
      alpha = a2
c
c     Mode will be used by search_dir on the next iteration to
c     figure out if the LS was complete or not.  To make this decision
c     easier reassign mode to be either ok or notok
c
      if (lsmode.eq.'bracket' .or. lsmode.eq.'downhill' .or. 
     $     lsmode.eq.'accept') then
         lsmode = 'ok'          ! Seemingly good situation
      else
         lsmode = 'notok'
      endif
c
      call driver_take_step(rtdb, geom)
c
c     There's code in the graveyard to do a more thorough search
c     for steps that have been restricted, but it actually seems
c     slightly counterproductive for most tests.  Dunno why.
c
      end
      subroutine driver_line_search2(rtdb,geom)
      implicit none
#include "nwc_const.fh"
#include "coptopt.fh"
#include "geom.fh"
#include "util.fh"
      integer rtdb, geom
c
c     Starting with a value of alpha recommended by whatever
c     generated the search direction, perform a more accurate 
c     line search than that done by search1.
c
c     Return with alpha, geom as the accepted point.
c
c     Only alpha and geom are modified.
c
c     This code does the same as driver_line_search1, but then
c     verifies that the new step is better (CPU used).  If it is, then
c     then return.  If the new step is not better, but the
c     original alpha value was good, then return that.
c
c     If both attempts are worse than the previous step, then the 'fun'
c     begins. Try in order, and return if step is good:
c     1.  reduce alpha to 1/4 its smallest value
c     2.  try minus that alpha value. If that works then:
c         a.  try 3 times that, and if that works then:
c         b.  try 3 times that.
c     3.  if ds*g is positive, then cut that ds in half, then try -half, then 0
c         this step is for situtations were delta E != dE for deltaX > 0
c     4.  lastly do a brute force line search
c
c     #2 & #3 can actually work because of limiting motion in some 
c     values of ds because of the hessian, trust, etc. can actually
c     result in a step that is going in the wrong direction.
c     Note that the #3 attempt will be symmetry corrected, etc.
c
      double precision e0, e1, e2, e2p, dsgrad
      double precision hess, a0, a1, a2
      double precision driver_energy_step
      double precision dssave(max_cart)
      integer i, j
      double precision walka(5), walke(5),dum
      logical success
c
      dsgrad = ddot(nvar, ds, 1, g, 1)
      if (dsgrad*alpha .ge. 0d0) then
         if (oprint) write(6,*) 
     $        ' Warning ... line search gradient +ve ', alpha, dsgrad
         alpha = -alpha
      endif
c
      e0 = energy
      a0 = 0d0
      walke(1) = e0
      walka(1) = a0
c
c     Compute energy at alpha*ds
c
      e1 = driver_energy_step(rtdb, geom)
      a1 = alpha
      walke(2) = e1
      walka(2) = a1
c
c     Fit to a parabola using e0, g0, e1
c
      hess = (e1-e0-alpha*dsgrad)/alpha**2
c
      if (abs(e1-e0) .lt. eprec) then
c
c     Insufficient precision in the energy to do anything else
c
c     If this condition persists we should switch off the LS.
c
         a2 = alpha 
         lsmode = 'fixed'
      else if (hess .gt. 0d0) then
c
c     Positive curvature.  Compute NR step, and ...
c     1) If we went downhill take the step with some restriction
c     2) If we went uphill then take the step if it is in [0,alpha]
c     .  otherwise take alpha/2.
c
         a2 = -dsgrad/(2d0*hess)
         if ((e1 - e0) .lt. eprec ) then ! Downhill
            lsmode = 'downhill'
            if (abs(a2) .gt. 4d0*abs(alpha)) then 
               lsmode = 'restrict'
               a2 = sign(4d0*abs(alpha),a2)
            endif
         else                   ! Uphill
            lsmode = 'bracket'
            if (abs(a2) .gt. abs(alpha)) then
               lsmode = 'confused'
               a2 = alpha/2d0
            endif
         endif
      else
c
c     Negative curvature.
c     1.  If we went downhill, step cautiously ahead.
c     2.  If we went uphill then something is really confused.
c     .   Just accept the current point.
c
         if ((e1 - e0) .lt. eprec ) then ! Downhill
            lsmode = 'negative'
            a2 = 2d0*alpha
         else
            lsmode = 'punt'
            a2 = alpha
         endif
      endif
c
c     Now have a2.  If it is close to alpha or does
c     not give us much energy lowering then just
c     return alpha so that any SCF/MCSCF/DFT does not have to
c     be repeated before computing the next gradient.
c
c     75% of energy lowering reduced to 99% so that we
c     are sure to have reduced line search gradient by 10x.
c
c     Added accept of alpha if predicted energy lowering is
c     less than energy precision.  Should diminish thrashing
c     around on flat surfaces.
c
      e2p = e0 + dsgrad*a2 + hess*a2**2
      if ( (abs(a2-alpha) .lt. abs(0.1d0*alpha)) .or.
     $     ((e2p-e1) .gt. 0.01*(e1-e0)) .or.
     $     (abs(e2p-e1) .lt. eprec) )then
         lsmode = 'accept'
         a2 = alpha
      endif
      e2p = e0 + dsgrad*a2 + hess*a2**2
c
      if (oprint) write(6,1) alpha, dsgrad, hess, e1, lsmode, a2, e2p
 1    format(' Line search: ',/, '     step=',f5.2, ' grad=',1p,d8.1,0p,
     $     ' hess=',1p,d8.1,0p,' energy=',f14.6,' mode=',a,/,
     $     ' new step=',f5.2, 19x,'predicted energy=',f14.6)
      alpha = a2
c
c     Mode will be used by search_dir on the next iteration to
c     figure out if the LS was complete or not.  To make this decision
c     easier reassign mode to be either ok or notok
c
      if (lsmode.eq.'bracket' .or. lsmode.eq.'downhill' .or. 
     $     lsmode.eq.'accept') then
         lsmode = 'ok'          ! Seemingly good situation
      else
         lsmode = 'notok'
      endif
c
      if (alpha .eq. a1) then
         e2 = e1
         walke(3) = e2 + 1d10 ! push duplicate to end
         walka(3) = alpha + 1d10
      else
         e2 = driver_energy_step(rtdb, geom)
         walke(3) = e2
         walka(3) = alpha
      endif
c
c     There are six possible energy orders
c     012,021,102,120,210,201
c
      if (e2 .le. e1 .and. e2 .le. e0) then ! 210&201
        ! Seconds step was best
      else if (e1 .le. e0) then ! 120&102
        alpha = a1 ! First step was best
        if (oprint) write(6,2) e2p, e2,alpha
 2         format( 'warning: line search expected E=',f14.6,' found E=',
     &             f14.6,' Backtracking to alpha=',f5.2)
      else ! 012&021
        ! This is not good, the previous geometry is lower than both
        lsmode = 'notok'
        if (abs(a2) .gt. abs(a1)) then
           alpha = 0.25d0*a1 ! Energy going up as you get smaller alpha!?
        else
           alpha = 0.25d0*a2
        endif
        e1 = driver_energy_step(rtdb, geom)
        walke(4) = e1
        walka(4) = alpha
        if (e1 .le. e0) then
           ! Got something good
           if (oprint) write(6,3) e2p, e2,alpha
 3           format( 'warning: line search wanted E=',f14.6,' found E=',
     &             f14.6,' Backtracked to alpha=',f5.2)
        else
           alpha = -alpha/2.0d0 ! Desperate - is possible because of gradient limiting, etc.
           e1 = driver_energy_step(rtdb, geom)
           walke(5) = e1
           walka(5) = alpha
           if (e1 .le. e0) then 
             e0 = e1
             alpha = alpha*3.0d0 ! Try a bigger step
             e1 = driver_energy_step(rtdb, geom)
             if (e1 .ge. e0) then
                alpha = alpha/3.0d0 ! Didn't work
             else
                e0 = e1
                alpha = alpha*3.0d0 ! Try a bigger step AGAIN!
                e1 = driver_energy_step(rtdb, geom)
                if (e1 .ge. e0) then
                   alpha = alpha/3.0d0 ! Didn't work
                endif
             endif
             if (oprint) write(6,4) e2p, e2,alpha
 4           format( 'warning: line search wanted E=',f14.6,' found E=',
     &             f14.6,' Flipped sign to alpha=',f5.2)
           else
              alpha=-alpha ! negative failed
              j = 0
              do i = 1, nvar
                dssave(i) = ds(i)
                dsgrad = ds(i)*g(i)
                if (dsgrad*alpha .gt. 0d0) then
                   j = j + 1
                   ds(i)=ds(i)*0.5d0
                endif
              enddo
              if ( j .ne. 0) then ! Try with new ds
                e1 = driver_energy_step(rtdb, geom)
                if (e1 .ge. e0) then ! Change ds again
                  do i = 1, nvar
                    dsgrad = ds(i)*g(i)
                    if (dsgrad*alpha .gt. 0d0) then
                       ds(i)=-ds(i)
                    endif
                  enddo
                  e1 = driver_energy_step(rtdb, geom)
                endif
                if (e1 .ge. e0) then ! Change ds again
                  do i = 1, nvar
                    dsgrad = dssave(i)*g(i)
                    if (dsgrad*alpha .gt. 0d0) then
                       ds(i)=0
                    endif
                  enddo
                  e1 = driver_energy_step(rtdb, geom)
                endif
              else
                e1 = e0 + 100 ! Fake it
              endif
              if (e1 .lt. e0) then
                  if (oprint) write(6,6) e2p, e2
 6          format( 'warning: line search wanted E=',f14.6,' found E=',
     &             f14.6,' changed step and alpha')
              else
!               Now we just linear search
                do i = 1, nvar
                   ds(i) = dssave(i)
                enddo
                do i = 1, 5
                  do j = i+1, 5 ! Evil bubble-ish sort
                    if (walka(j) .lt. walka(i)) then
                      dum = walka(j)
                      walka(j) = walka(i)
                      walka(i) = dum
                      dum = walke(j)
                      walke(j) = walke(i)
                      walke(i) = dum
                    endif
                  enddo
                enddo
!               We know that the alpha=0 point is lowest energy, and we
!               know that is ia walka(2), since we only try one negative point
!               So, walka(4&5) are not relevent
!               On steps 1,4,7,10,... do right side
!               On steps 2,5,8,11,... do left side
!               On steps 3,6,9,12,... do side which is lower in energy
!               We set step to 10% of distance, until we get a lower point,
!               then we use 40%
!               We keep looking even after we find a good point, because this
!               is clearly a bad step, so we work it out so that the next step
!               is not bad.
                success = .false.
 777            continue
                do i = 1, 12
                 if((mod(i,3) .eq. 1) .or.
     >            (mod(i,3) .eq. 0 .and. walke(3) .lt. walke(1))) then ! Right
                  if (.not. success) then
                    alpha = (walka(3)-walka(2))*0.1d0 + walka(2)
                  else
                    alpha = (walka(3)-walka(2))*0.4d0 + walka(2)
                  endif
                  e1 = driver_energy_step(rtdb, geom)
                  if (e1 .gt. walke(2)) then ! failed
                    walke(3) = e1
                    walka(3) = alpha
                  else
                    walke(1) = walke(2)
                    walka(1) = walka(2)
                    walke(2) = e1
                    walka(2) = alpha
                    success = .true.
                  endif
                 else ! Left
                  if (.not. success) then
                    alpha = (walka(2)-walka(1))*0.9d0 + walka(1)
                  else
                    alpha = (walka(2)-walka(1))*0.6d0 + walka(1)
                  endif
                  e1 = driver_energy_step(rtdb, geom)
                  if (e1 .gt. walke(2)) then ! failed
                    walke(1) = e1
                    walka(1) = alpha
                  else
                    walke(3) = walke(2)
                    walka(3) = walka(2)
                    walke(2) = e1
                    walka(2) = alpha
                    success = .true.
                  endif
                 endif
                enddo
                if (success) then
                  alpha = walka(2)
                  if (oprint) write(6,9) e2p, e2,alpha
 9           format( 'warning: line search wanted E=',f14.6,' found E=',
     &             f14.6,' Changed alpha to alpha=',f5.2)
                else
                   alpha = walka(3)*0.5d0
                   if (oprint) write(6,10) e2p, e2,alpha
 10          format( 'warning: line search wanted E=',f14.6,' found E=',
     &             f14.6,' Giving up, set alpha=',f5.2)
                endif
             endif
           endif
        endif
      endif
c
      call driver_take_step(rtdb, geom)
c
c     There's code in the graveyard to do a more thorough search
c     for steps that have been restricted, but it actually seems
c     slightly counterproductive for most tests.  Dunno why.
c
      end
      double precision function driver_energy_step(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
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
      logical opt_geom_cart_coords_get
      logical opt_geom_cart_coords_set
c
c     Get original coordinates so can restore them
c
c     FRACTIONAL?
      if (.not. opt_geom_cart_coords_get(geom, xold))
     $     call errquit('driver_energy_step: coordinates?',geom,
     &       GEOM_ERR)
c
c     Take the step
c
      call driver_take_step(rtdb, geom)
c
c     Energy at displaced point
c
      if (.not. task_energy(rtdb)) call errquit
     $     ('driver_energy_step: energy failed',0, GEOM_ERR)
      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, aenergy))
     $     call errquit('driver_energy_step: rtdb?',0, RTDB_ERR)
c      
c     Restore the starting geometry
c
c     FRACTIONAL?
      if (.not. opt_geom_cart_coords_set(geom, xold))
     $     call errquit('driver_energy_step: coordinates?',geom,
     &       GEOM_ERR)
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $     call errquit('driver_energy_step: grs?',geom, RTDB_ERR)
c
      driver_energy_step = aenergy
c
      end
      subroutine driver_take_step(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
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
c     The geom is modified, and xmax/xrms are computed from the
c     first-order step.
c
      double precision xold(max_cart), xnew(max_cart), err
      double precision aaa(3,3)
      integer i, l_bi, k_bi
      logical opt_geom_cart_coords_get
      logical opt_geom_cart_coords_set
      external opt_geom_cart_coords_get,opt_geom_cart_coords_set
      integer sym_number_ops
      external sym_number_ops
c
c     Get original coordinates
c
c     FRACTIONAL?
      if (.not. opt_geom_cart_coords_get(geom, xold))
     $     call errquit('driver_energy_step: coordinates?',geom,
     &       GEOM_ERR)
c
c     Take the step
c
      if (zcoord) then
c
c     Compute first order step only so we can compute xmax/xrms
c
         if (.not. ma_push_get(mt_dbl, ncart*nzvar,'binv',l_bi, k_bi))
     $        call errquit('driver_step: memory for binv', ncart*nzvar,
     &       MA_ERR)
         call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nzvar)
         call dgemv('n', ncart, nzvar, alpha, dbl_mb(k_bi), ncart,
     $        ds, 1, 0.0d0, xnew, 1)
         if (.not. ma_pop_stack(l_bi)) call errquit('driver: ma?',99,
     &       MA_ERR)
         if (odebug) then
            write(6,*) ' First-order internal coord step in cartesians'
            call doutput(xnew, 1, ncart, 1, 1, ncart, 1, 1)
         endif
         call sym_grad_symmetrize(geom, xnew) ! This symmetrizes the step.
         do i = 1, ncart
            xnew(i) = xold(i) + xnew(i)
         enddo
c
c     Now take the step for real in internals
c
c     FRACTIONAL?
         call geom_update_cart_from_int(geom, ds, alpha, err)
c
*         if (.not. opt_geom_cart_coords_set(geom, xnew))
*     $        call errquit('driver_energy_step: coordinates?',geom,
*     &       GEOM_ERR)
      else
         do i = 1, ncart
            xnew(i) = alpha*ds(i)
         enddo
         if( sym_number_ops(geom) .ne. 0) then
          if (isystype.ne.0 .and. (.not. zcoord)) then
            if (.not. geom_frac_to_cart(geom, xnew))
     $           call errquit('driver: frac_to_cart?',0,0)
          endif
          ! This assumes cartesian - internally converts to fractional!
          call sym_grad_symmetrize(geom, xnew)
          if (isystype.ne.0 .and. (.not. zcoord)) then
              if (.not. geom_cart_to_frac(geom, xnew))
     $         call errquit('driver: cart_to_frac?',0,0)
          endif
         endif
         do i = 1, ncart
            xnew(i) = xold(i) + xnew(i)
         enddo
c     FRACTIONAL?
         if (.not. opt_geom_cart_coords_set(geom, xnew))
     $        call errquit('driver_energy_step: coordinates?',geom,
     &       GEOM_ERR)
      endif
c
c     Must ensure the geometry has the required symmetry even after
c     enforcing it on the step.  Should use
c     an error criterion consistent with the step size.
c     For periodic systems this maps back to central cell
c
      call sym_geom_project(geom, trust)
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
     $     call errquit('driver_energy_step: grs?',geom, RTDB_ERR)
c
      end
      subroutine driver_print(geom, istep)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
      integer geom, istep
c
c     Print out stuff
c
      integer i
      double precision diff, de
      logical geom_print_zmatrix, oprint_grad
      character*9 cvg1, cvg2, cvg3, cvg4
      character*1 mark
      character*255 filename, dir
      logical geom_print_xyz
c
      de = 0d0
      if (istep .gt. 1) de = energy-energyp
      cvg1 = ' '
      cvg2 = ' '
      cvg3 = ' '
      cvg4 = ' '
      if (gmax .lt. gmax_tol) cvg1 = '     ok  '
      if (grms .lt. grms_tol) cvg2 = '     ok  '
      if (xrms .lt. xrms_tol) cvg3 = '     ok  '
      if (xmax .lt. xmax_tol) cvg4 = '     ok  '
C     Don't print OK on first step, if value is not real
      if (istep.eq.1 .and.(xrms.gt.xrms_tol .or. gmax.gt.gmax_tol)) then
        if (xrms .eq. 0) cvg3 = ' '
        if (xmax .eq. 0) cvg4 = ' '
      endif
c
      if (oprint) then 
         mark = '@'
         if (istep .gt. 1) mark = ' '
         write(luout,1) mark, mark
         mark = '@'
         write(luout,2) mark, istep-1, energy, de,
     $     gmax, grms, xrms, xmax, util_wallsec(),
     $     cvg1, cvg2, cvg3, cvg4
 1       format(
     $        /,a1,' Step       Energy      Delta E   Gmax',
     $        '     Grms     Xrms     Xmax   Walltime',
     $        /,a1,' ---- ---------------- -------- --------',
     $        ' -------- -------- -------- --------')
 2       format(
     $        a1,i5,f17.8,1p,d9.1,0p,4f9.5,f9.1,/,
     $        1x,5x,17x,9x,4a9,/)
      endif
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
     $        write(luout,*) ' Internal coordinate gradient before',
     $        ' constraints'
         if (.not. geom_print_zmatrix(geom, gq, 'Gradient',.true.))
     $        call errquit('tfgx: print internal gradient?',geom,
     &       GEOM_ERR)
         if (diff .gt. 1d-6) then
            write(luout,*) ' Internal coordinate gradient after',
     $           ' constraints'
            if (.not. geom_print_zmatrix(geom, g, 'Gradient',.true.))
     $           call errquit('tfgx: print internal gradient?',geom,
     &       GEOM_ERR)
         endif
      endif
c
      if (ga_nodeid().eq.0 .and. xyz.ne.' ') then
         dir      = ' '
         filename = ' ' 
         call util_directory_name(dir, .false., 0)
         write(filename,13) dir(1:inp_strlen(dir)), 
     $        xyz(1:inp_strlen(xyz)),istep-1
 13      format(a,'/',a,'-',i3.3,'.xyz')
         open(88,file=filename,form='formatted',status='unknown',
     $        err=133)
         if (.not. geom_print_xyz(geom, 88))
     $        call errquit('driver: geom_print_xyz?',0, GEOM_ERR)
         close(88,status='keep',err=133)
      endif
c
      return
 133  call errquit('driver: error opening/closing xyz file',0, GEOM_ERR)
c
      end
      logical function driver_converged()
      implicit none
#include "nwc_const.fh"
#include "coptopt.fh"
c
c     Return true if we have converged
c
c     Nothing is modified.  Assumes driver_compute_info()
c     has been called.
c     gmax_tol,            ! [user] tolerance for max internal gradient
c     grms_tol,            ! [user] tolerance for rms internal gradient
c     xrms_tol,            ! [user] tolerance for rms cartesian step
c     xmax_tol,            ! [user] tolerance for max cartesian step
c
      driver_converged = 
     $     ((gmax .lt. gmax_tol) .and. (grms .lt. grms_tol) .and.
     $     (xrms .lt. xrms_tol)  .and. (xmax .lt. xmax_tol))
     $     .or. 
     $     ((gmax.lt.0.01d0*gmax_tol) .and. (grms.lt.0.01d0*grms_tol))

c
      end
      subroutine driver_cart_pmat(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer rtdb, geom
c     
c     Compute the cartesian equivalent of the P = G.G^-1 matrix
c     which projects to and from the linearly independent
c     set of coordinates.  In the cartesian case P is the complement
c     of the projector onto the rotations and translations
c     For ease of use we also write out a unit Binv matrix.
c     
c     Only the P/Binv matrices are generated.  Nothing is modified.
c
c     Minor little catch is that if some atoms are being frozen
c     we are no longer invariant to translations or rotations.
c     
c     RTDB is used to look for frozen atoms.
c     
      double precision centroid(3), x, y, z, xx, yy, zz, fx
      double precision coords(3,max_cent)
      double precision work(max_cart,6)
      integer i, j, k, l_pmat, k_pmat, i3, ma_type, nelem
      character*26 date
      integer ind
      logical omm
      logical task_qmmm
      logical opt_geom_cart_coords_get
      ind(i,j) = k_pmat + i-1 + (j-1)*ncart
c     
c     FRACTIONAL?
      if (.not. opt_geom_cart_coords_get(geom, coords))
     $     call errquit('driver_cart_pmat: geom?',geom, GEOM_ERR)
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
     $     ('driver_cart_pmat: memory for pmat',ncart**2, GEOM_ERR)
c     
c     Form unit matrix
c     
      call dfill(ncart**2, 0d0, dbl_mb(k_pmat), 1)
      call dfill(ncart, 1d0, dbl_mb(k_pmat), ncart+1)
c     
c     Store dummy unit matrix for B, Binv ... the cartesian 
c     gradient should already be invariant to rotations and translations.
c     Also store dummy unit matrix for cmat (constraints)
c     
      call geom_hnd_put_data('b', dbl_mb(k_pmat), ncart**2)
      call geom_hnd_put_data('b^-1', dbl_mb(k_pmat), ncart**2)
      call geom_hnd_put_data('c', dbl_mb(k_pmat), ncart**2)
c     
      if (.not.rtdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     &    task_qmmm = .false.

      if( .not. rtdb_get(rtdb,'task:mm',mt_log,1,omm))
     &    omm = .false.

      if ( rtdb_get_info(rtdb, 'geometry:actlist', ma_type, 
     $     nelem, date) .or.
     $     rtdb_get_info(rtdb, 'geometry:inactlist', ma_type, 
     $     nelem, date) .or.
     $     isystype .ne. 0 .or.
     $     task_qmmm .or. 
     $     omm .or. 
     $     geom_extbq_on() ) then
c     
c     Some atoms are frozen or we have a periodic system so don't have
c     invariance ...  also store unit matrix for P.
c     or we have QMMM calculation here
c     
         call geom_hnd_put_data('p', dbl_mb(k_pmat), ncart**2)
c     
      else
c     
c     Finish P
c     
         do i = 1, ncart
            do j = 1, ncart
               do k = 1, 6
                  dbl_mb(ind(j,i)) = dbl_mb(ind(j,i)) - 
     $                 work(j,k)*work(i,k)
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
      endif
c     
      if (.not. ma_chop_stack(l_pmat)) call errquit
     $     ('driver_cart_bmat: ma?',0, MA_ERR)
c     
      end
      subroutine driver_project_hess_grad(hess, pg)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "mafdecls.fh"
      double precision 
     $     hess(nvar,nvar),     ! returns projected & shifted Hessian
     $     pg(nvar)             ! returns projected gradient
c
c     Project and shift the Hessian and gradient following Peng et al.
c
c     Nothing else is changed.
c
      integer l_pmat, k_pmat, l_work, k_work, i
      double precision big
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'work',
     $     l_work, k_work)) call errquit
     $     ('driver_proj_h_g: memory for pmat',nvar**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'pmat',
     $     l_pmat, k_pmat)) call errquit
     $     ('driver_proj_h_g: memory for work',nvar**2, MA_ERR)
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
     $     call errquit('freddy',0, MA_ERR)
c
c     PG
c
      call dgemv('n',nvar, nvar, 1d0, dbl_mb(k_pmat), nvar,
     $     g, 1, 0d0, pg, 1)
      if (odebug) then
         write(6,*) ' Gradient after projection'
         call doutput(pg, 1, nvar, 1, 1, nvar, 1, 1)
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
     $     ('driver_p_h_g:ma?',0, MA_ERR)
c
      end
      subroutine driver_compute_info()
      implicit none
#include "nwc_const.fh"
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
      subroutine driver_hess_cart_guess()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "inp.fh"
c     
c     Read in cartesian Hessian and transform it as necessary
c     to internal coordinates (neglecting the component due to
c     the derivative) and writing the result to the hessian file.
c     
c     Reads file in vib_vib format using vib_vib filename default
c     Note the default filename is set in task_freq
c     filenames must be made identical.
c     
c     Format of vib file is ascii lower triangular elements only.
c     
      integer h_unit
      parameter (h_unit=47)
      character*255 fname
      double precision x
      integer i,j
      integer l_bi, k_bi, l_hc, k_hc, l_hq, k_hq
c     
      if (.not. ma_push_get(mt_dbl, ncart*nvar, 'binv',
     $     l_bi, k_bi)) call errquit
     $     ('driver_hess_cart_guess: ma?', ncart*nvar, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart',
     $     l_hc, k_hc)) call errquit
     $     ('driver_hess_cart_guess: ma?', ncart**2, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart2',
     $     l_hq, k_hq)) call errquit
     $     ('driver_hess_cart_guess: ma?', nvar**2, MA_ERR)
c     
      if (ga_nodeid().eq.0) then
         call util_file_name('hess',.false.,.false.,fname)
         open(unit=h_unit,file=fname,form='formatted',status='unknown',
     $        err=99990,access='sequential')
         rewind h_unit
         do i = 1,ncart
            do j = 1,i
               read(h_unit,10000,err=99992,end=99992) x
               dbl_mb(k_hc+(i-1)*ncart+(j-1)) = x
               dbl_mb(k_hc+(j-1)*ncart+(i-1)) = x
            enddo
         enddo
         close(unit=h_unit,status='keep')
      endif
      call ga_brdcst(1,dbl_mb(k_hc),8*ncart**2,0)
c     
      call geom_hnd_get_data('b^-1', dbl_mb(k_bi), nvar*ncart)
      call dgemm('n', 'n', ncart, nvar, ncart, 1d0, dbl_mb(k_hc), ncart,
     $     dbl_mb(k_bi), ncart, 0d0, dbl_mb(k_hq), ncart)
      call dgemm('t', 'n', nvar, nvar, ncart, 1d0, dbl_mb(k_bi), ncart,
     $     dbl_mb(k_hq), ncart, 0d0, dbl_mb(k_hc), nvar)
c
      do i = 1,nvar
         do j = 1,i
            x = (dbl_mb(k_hc+(i-1)*nvar+(j-1)) + 
     $           dbl_mb(k_hc+(j-1)*nvar+(i-1))) * 0.5d0
            dbl_mb(k_hc+(i-1)*nvar+(j-1)) = x
            dbl_mb(k_hc+(j-1)*nvar+(i-1)) = x
         enddo
      enddo
c
      call geom_hnd_put_data('drv.hess',dbl_mb(k_hc), nvar**2)
c
      if (.not. ma_chop_stack(l_bi))
     $     call errquit('driver_hess_cart_guess: ma corrupt?',0, MA_ERR)
c
      return
10000 format(f30.15)
99990 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as unknown file'
      call errquit('driver_hess_cart: fatal error', 911, GEOM_ERR)
99991 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as new file'
      call errquit('driver_hess_cart: fatal error', 911, GEOM_ERR)
99992 write(6,*)' error in reading <',fname(1:inp_strlen(fname)),
     $     '> as hessian file'
      call errquit('driver_hess_cart: fatal error', 911, GEOM_ERR)
      end
      subroutine driver_sad_search_dir(rtdb,geom,istep)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
      integer rtdb
      integer geom
      integer istep
c     
c     This routine for saddle point search
c     
c     Put into ds() a search direction in the optimization
c     variables (internal or cartesian) based upon the
c     current gradient, g(), and Hessian.  Apply constraints.
c     
c     Only alpha and ds() are modified.  Saves moddir to RTDB
c     
c     On first geometry optimization step need to identify
c     the search direction.
c     
      integer i
*     integer info
      integer l_hess, k_hess, l_work, k_work, lenwork
      double precision  e(max_nvar) ! Hessian eigenvalues
      double precision pg(max_nvar) ! P.g
      double precision gv(max_nvar) ! Gradient along eigenvectors
      double precision dv(max_nvar) ! Step along eigenvectors
      double precision saddir(max_nvar) ! Tracks mode for saddle ... saved
      double precision evalp    ! Remembers previous search eval ... saved
      double precision dsmax    ! Max. value of current step (smax is prev.)
c     
      integer nzero, negeig, idamax, mod2, ntest
      double precision bohr, deg ! For printing purposes
      double precision trustds  ! Restriction of step in opt. variable
      double precision gvrms, dir, tmp1, tmp2, ovtol
c
      double precision smalleig
      parameter (smalleig = 1d-4)
c
      logical geom_print_zmatrix, omatchneg
      save saddir, evalp
c     
      omatchneg = .true.
      ovtol     = 0.7d0
c     
c     Get the Hessian and gradient with appropriate projects
c     applied following Peng, Ayala, Schlegel and Frisch so that
c     redundant internal modes are shifted to high eigenvalues.
c     
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('driver_search_dir: memory for hessian',nvar**2, GEOM_ERR)
      call driver_project_hess_grad(dbl_mb(k_hess), pg)
c     
c     Diagonalize the Hessian.  Should really do the generalized
c     eigenvalue problem since the underlying basis is not independent
c     (if we are using autoz). NOT YET BEING DONE.  
c     
c     To cause degenerate eigenvalues to be resolved into symmetry
c     adapted combinations use Jacobi not dsyev and screen out junk
c     
      lenwork = max(nvar**2,100)
      if (.not. ma_push_get(mt_dbl, lenwork, 'work',
     $     l_work, k_work)) call errquit
     $     ('driver_search_dir: memory for hessian', lenwork, GEOM_ERR)
      do i = 0, nvar**2-1
         if (abs(dbl_mb(k_hess+i)).lt.1d-8) dbl_mb(k_hess+i) = 0d0
      enddo
c$$$  call dsyev('v', 'u', nvar, dbl_mb(k_hess), 
c$$$  $     nvar, e, dbl_mb(k_work), lenwork, info)
c$$$  if (info .ne. 0) call errquit
c$$$  $     ('driver_search_dir: dsyev failed',info, GEOM_ERR)
      call util_jacobi(nvar, dbl_mb(k_hess), nvar, e)
      if (odebug .or. (util_print('hvecs',print_never) 
     $     .and. ga_nodeid().eq.0)) then
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
c     
      gvrms = ddot(nvar, gv, 1, gv, 1)
      gvrms = sqrt(gvrms/dble(nvar))
      if (odebug) then
         write(6,*) ' Gradient in spectral form ', gvrms
         call doutput(gv, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c     
c     Analyze the structure of the Hessian.  Don't count
c     modes with significant gradients as having zero eigenvalues.
c     
      nzero =0
      negeig=0
      do i=1,nvar
         if(abs(e(i)).lt.smalleig.and.abs(gv(i)).lt.gmax_tol) then
            nzero=nzero+1
         else if(e(i) .lt. -smalleig) then
            negeig=negeig+1
         endif
      enddo     
      if (oprint) then
         write(6,11) nvar-negeig-nzero, negeig, nzero
 11      format(' Hessian eigenvalues: positive=',i4,
     $        ' negative=',i3,' zero=',i3)
      endif
c     
c     Translations and rotations have already been shifted to
c     the right so zero modes are rather weird.  However, they
c     can arise if dummy centers have been specified in user
c     input z-matrices.  They can also arise by accident, or
c     if a cartesian coordinate has been frozen but partial rot/trans
c     invariance remains. Shift remaining zero modes to the right.
c     
      if(nzero.gt.0) then
         do i=negeig+nzero+1,nvar ! COPY LOGIC INCOMPLETE WITH GV MOD ABOVE !!
            e(i-nzero)  = e(i)
            gv(i-nzero) = gv(i)
            call dcopy(nvar, dbl_mb(k_hess+(i-1)*nvar), 1,
     $           dbl_mb(k_hess+(i-nzero-1)*nvar), 1)
         enddo
         do i=nvar-nzero+1,nvar
            e(i)  = 1000d0
            gv(i) = 0d0
            call dfill(nvar, 0d0, dbl_mb(k_hess+(i-1)*nvar), 1)
         enddo
      endif
c     
c     Look for the internal coordinate or mode to be followed.  On the
c     first iteration step along the mode/variable specified by the user.
c     Otherwise do maximum overlap to track the mode ... must track the
c     mode not the overall step since the step in the modes being minimized
c     might be bigger than the uphill step.  
c     
c     Also, might not want to take the first negative mode that appears 
c     since it might be a transverse mode rather than the one of interest.
c     However, since we might not be tracking the correct mode the default
c     option is to just follow the first -ve mode that appears.  This
c     can be disabled with nofirstneg.
c     
      dir = 1d0                 ! Sign of very first step ... +ve default
      if (istep .eq. 1) then
         if (vardir .ne. 0) then
            evalp = 99d0        ! Any +ve value will do
            if (vardir .lt. 0) dir = -1d0
            vardir = abs(vardir)
            call dfill(nvar, 0d0, saddir, 1)
            call dfill(nvar, 0d0, ds, 1)
            saddir(vardir) = 1d0
            ds(vardir) = sign(sadstp,pg(vardir))*dir
            if (oprint) then
               write(6,776) vardir, pg(vardir), dir, ds(vardir)
 776           format(' Initial step taken in variable ', 
     $              i5,' grad=', f9.4,' dir=', f4.1, ' step=',f9.4)
            endif
         else if (moddir .ne. 0) then
            if (moddir .lt. 0) dir = -1d0
            moddir = abs(moddir)
            evalp = e(moddir)
            call dcopy(nvar,dbl_mb(k_hess+(moddir-1)*nvar),1,saddir,1)
         else
c
c     Default first step is uphill ... makes most sense if we have
c     just done a constrained optimization.  Screen out just the
c     largest few components.
c
c     Screening with just /3d0 fails for Claisen rearrangement 
c     case in peng & schelgel using some coord systems.
c
            tmp1 = abs(pg(idamax(nvar,pg,1))) / 5d0  ! was 10d0
            do i = 1, nvar
               saddir(i) = 0d0
*
*     Tried scaling gradient by 0.1 but it breaks things!!!
*     A big step seems to be good.
*
               if (abs(pg(i)) .gt. tmp1) saddir(i) = pg(i)
            enddo
            call dcopy(nvar, saddir, 1, ds, 1)
            if (oprint) write(6,*) 'Initial step taken uphill'
c     Normalize saddir for future max overlap
            tmp1 = 1d0/sqrt(ddot(nvar,saddir,1,saddir,1))
            call dscal(nvar, tmp1, saddir, 1)
         endif
      else if (ofirstneg .and. negeig.eq.1) then
         if (oprint) write(6,34) 
     $        ddot(nvar, saddir, 1, dbl_mb(k_hess), 1)
 34      format(' Hessian has correct structure - ',
     $        'following negative mode, overlap=',1p,d8.1)
         moddir = 1
      else
c     
c     Maximum overlap for tracking the mode.  
c     
c     Modified to 
c     a) If we were following a -ve mode just match against the 
c     .  negative modes, unless omatchneg is false.
c     .  (omatchneg=.false. recovers previous behaviour)
c     b) match the lowest mode that has an overlap which is 
c     .  greater than ovtol(=0.5)*max_overlap.
c     .  (ovtol=1.0 recovers previous behaviour)
c     
         call dgemv('t', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $        saddir, 1, 0d0, dv, 1)
         do i = 1, nvar         ! This to avoid matching rot/trans/constraint
            if (e(i) .gt. 1d2) dv(i) = 0d0
         enddo
         if (util_print('modeoverlap',print_debug) .and. 
     $        ga_nodeid().eq.0) then
            write(6,*) ' Overlap of modes with searchdir'
            call output(dv, 1, nvar, 1, 1, nvar, 1, 1)
         endif
c     
         ntest = nvar
         if (omatchneg .and. (negeig.gt.1) .and. (evalp.lt.0d0)) then
            if (oprint) write(6,*) 
     $           'Maximum overlap only with negative modes'
            ntest = negeig
         endif
c     
         moddir = idamax(ntest, dv, 1)
         do i = 1, ntest
            if (abs(dv(i)) .gt. abs(dv(moddir))*ovtol) then
               moddir = i
               goto 304
            endif
         enddo
 304     continue
c     
c     Print out first two maximum overlaps to warn of ambiguity
c     
         tmp1 = abs(dv(moddir))
         dv(moddir) = 0d0
         mod2 = idamax(ntest, dv, 1)
         tmp2 = abs(dv(mod2))
         if (oprint) write(6,98) moddir, tmp1, mod2, tmp2
 98      format(' Maximum overlap of modes: mode=', i3, 1p, 
     $        ' overlap=',d8.1, 0p, '  mode=', i3, 1p,' overlap=',d8.1)
         if (tmp1 .eq. 0d0) call errquit('driver_sad_dir: no mode?',0,
     &       GEOM_ERR)
c
      endif
      if (moddir .ne. 0) then
         call dcopy(nvar, dbl_mb(k_hess+(moddir-1)*nvar), 1, saddir, 1)
         evalp = e(moddir)
      endif
      if (util_print('sadmode',print_default).and.ga_nodeid().eq.0)then
         write(6,*) 
         write(6,*) 'The mode being followed to the saddle point'
         if (zcoord) then
            if (.not.geom_print_zmatrix(geom, saddir,'Sadmode',.true.))
     $           call errquit('driver: printing sadmode?',0, GEOM_ERR)
         else
            write(6,*)
            write(6,88) (saddir(i),i=1,nvar)
 88         format(1x,8f9.4)
         endif
      endif
c     
c     Save moddir so that it is picked up for a restart
c     
      if (.not. rtdb_put(rtdb,'driver:moddir',mt_int,1,moddir)) 
     $     call errquit('driver_sad_search_dir: rtdb?',0, RTDB_ERR)
c     
c     First step might be uphill or along a variable.  Otherwise
c     we take the approx. augmented newton step forcing uphill
c     motion along the identified mode.
c
      if (moddir .gt. 0) then
c     
c     Form the step applying trust radius on each mode.
c     Augmented Hessian step is similar to replacing
c     eigenvalue e(i) with (e(i) + gv(i)**2/e(i))
c     
         call dfill(nvar, 0d0, dv, 1)
         do i = 1, nvar
            if (i .eq. moddir) then
               if (e(i) .ge. 0d0) then ! +ve curvature ... fixed step uphill
c     
c     +ve curvature and first step ... force in the direction
c     requested by the user.  Otherwise go uphill.
c     
                  if (istep.eq.1) then
                     dv(i) = sadstp*dir
                  else
                     dv(i) = sign(sadstp,gv(i))
                  endif
c     
                  if (oprint) write(6,667) i, e(i), gv(i), dv(i)
 667              format(' Forcing    uphill step in mode',i5,
     $                 ' eval=', 1p,d8.1,' grad=',d8.1,
     $                 ' new=',d8.1)
               else
*     dv(i) = -gv(i) / e(i)
                  dv(i) = -gv(i)*dir/(e(i)+max(gvrms**2,gv(i)**2)/e(i))
                  if (abs(dv(i)) .gt. sadstp) then 
                     if (oprint) write(6,2) i, e(i), dv(i), 
     $                    sign(sadstp,dv(i))
                     dv(i) = sign(sadstp,dv(i))
                  endif
               endif
            else
               if (e(i) .gt. smalleig) then ! +ve curvature
*     dv(i) = -gv(i) / e(i)
                  dv(i) = -gv(i) / (e(i) + max(gvrms**2,gv(i)**2)/e(i))
               else if (e(i) .lt. -smalleig) then ! -ve curvature 
                  dv(i) = sign(trust,-gv(i)) ! Used to be 0.3*trust
                  if (oprint) write(6,666) i, e(i), gv(i), dv(i)
 666              format(' Forcing  downhill step in mode',i5,
     $                 ' eval=', 1p,d8.1,' grad=',d8.1,' new=',
     $                 d8.1)
               else
                  dv(i) = sign(min(trust,abs(gv(i))),-gv(i))
                  if (oprint) write(6,1) i, e(i), dv(i)
 1                format(' Using steepest descent step in mode',i5,
     $                 ' due to tiny eigenvalue=', 1p,d8.1,
     $                 ' step=',d8.1)
               endif
               if (abs(dv(i)) .gt. trust) then 
                  if (oprint) write(6,2) i, e(i), dv(i), 
     $                 sign(trust,dv(i))
 2                format(' Restricting large step in mode',i5, 
     $                 ' eval=', 1p,d8.1,' step=',d8.1,' new=',d8.1)
                  dv(i) = sign(trust,dv(i))
               endif
            endif
         enddo
         if (oprint) then
            write(6,12) moddir,e(moddir),gv(moddir),dv(moddir)
 12         format(' Search mode=',i2,' eval=',1p,d8.1,' grad=',d8.1,
     $           ' step=',d8.1)
         endif
         if (odebug) then
            write(6,*) ' Step in spectral form '
            call doutput(dv, 1, nvar, 1, 1, nvar, 1, 1)
         endif
c     
c     Transform back to optimization space
c     
         call dgemv('n', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $        dv, 1, 0d0, ds, 1)
      endif
c     
      call driver_symmetrize_step(geom)
c     
c     Restriction of the overall step based on
c     
c     1.  Is a component of the step large in the optimization space.
c     
      dsmax = 0d0
      do i = 1, nvar
         dsmax = max(dsmax, abs(ds(i)))
      enddo
c     
      trustds = min(trust,sadstp)
c     
      alpha = 1d0
c     
      if (dsmax .gt. trustds) then
         alpha = trustds/dsmax
         if (oprint) write(6,7) alpha
 7       format(' Restricting overall step due to large component.',
     $        ' alpha=', f6.2)
      endif
c     
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('driver_search_dir: ma?',0, MA_ERR)
c     
c     Edo seems to have encountered a case where different processors
c     generated different steps.  To prevent this, broadcast the 
c     critical info to everyone.
c     
      call ga_brdcst(1,ds,8*nvar,0)
      call ga_brdcst(2,alpha,8,0)
c     
      if (util_print('searchdir',print_high) .and. 
     $     ga_nodeid().eq.0) then
         write(6,*)
         write(6,*) '       The search direction'
         if (zcoord) then
            bohr = 0.52917715d0
            deg  = 0.52917715d0*180d0/(4d0*atan(1d0))
            call geom_zmat_ico_scale(geom, ds, bohr, deg)
            if (.not. geom_print_zmatrix(geom, ds, 'Search  ',.true.))
     $           call errquit('driver:print search_dir?',geom, GEOM_ERR)
            call geom_zmat_ico_scale(geom, ds, 1d0/bohr, 1d0/deg)
         else
            call output(ds,1,ncart,1,1,ncart,1,1)
         endif
         write(6,*)
         call util_flush(6)
      endif
c     
      end
      subroutine driver_symmetrize_step(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "mafdecls.fh"
      integer geom
      logical geom_frac_to_cart,geom_cart_to_frac
      integer sym_number_ops
      external geom_frac_to_cart,geom_cart_to_frac,sym_number_ops
c
c     Force symmetry upon the current search direction by projecting
c     ds() onto symmetric component in cartesians
c
c     Updates ds().
c
      double precision dx(max_nvar)
      integer k_bi, l_bi
c
c     1) dx = dq*B-1
c     2) symmetrize dx
c     3) dq = dx*B
c
      if (sym_number_ops(geom) .eq. 0) return
      if (ostress) then
        call sym_lattice_stress(geom,ds(nvar-8))
      endif
      if (ostress2) then
        call sym_lattice_stress2(geom,ds(nvar-5))
      endif
! Thankfully, the code below can just ignore stress terms at end of ds
      if (.not. zcoord) then
        if (isystype.ne.0 .and. (.not. zcoord)) then
           if (.not. geom_frac_to_cart(geom, ds))
     $          call errquit('driver: frac_to_cart?',0,0)
        endif
        ! This assumes cartesian - internally converts to fractional!
        call sym_grad_symmetrize(geom, ds)
        if (isystype.ne.0 .and. (.not. zcoord)) then
            if (.not. geom_cart_to_frac(geom, ds))
     $       call errquit('driver: cart_to_drac?',0,0)
        endif
        return
      endif

      if (.not. ma_push_get(mt_dbl, ncart*nvar,'binv',l_bi, k_bi))
     $     call errquit('driver_sym: memory for binv', ncart*nvar,
     &       MA_ERR)
      call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nvar)
      if (odebug) then
         write(6,*) ' Symmetrize step - initial q '
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
      call dgemv('n', ncart, nvar, 1d0, dbl_mb(k_bi), ncart,
     $     ds, 1, 0.0d0, dx, 1)
      call sym_grad_symmetrize(geom, dx)
      call geom_hnd_get_data('b', dbl_mb(k_bi), ncart*nvar)
      call dgemv('t', ncart, nvar, 1d0, dbl_mb(k_bi), ncart,
     $     dx, 1, 0.0d0, ds, 1)
      if (odebug) then
         write(6,*) ' Symmetrize step - final q '
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
      if (.not. ma_pop_stack(l_bi)) call errquit('driver_sym:ma',0,
     &       MA_ERR)
c
      end
      subroutine driver_compute_actual_step(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "geom.fh"
      integer geom
c
c     We have now taken a step.  Since the non-linear transformations
c     involved in taking a step may not have been done exactly replace
c     ds() with the actual step taken so that the Hessian may be precisely
c     updated.  Updates ds(), sp().  Divides by alpha so the step is still
c     alpha*ds().
c
c     This has little effect on most calculations but for (h2o)5 it
c     reduces the number of iterations from 99 to 69.
c
      logical geom_compute_zmatrix
      double precision bohr, deg
      integer i
      logical opt_geom_cart_coords_get
c
      if (odebug) then
         write(6,*) ' Expected ds '
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
      call dcopy(nvar, sp, 1, ds, 1) ! Old coordinates into ds()
      if (zcoord) then
         if (.not. geom_compute_zmatrix(geom, sp)) ! New coordinates into sp()
     $        call errquit('driver: zmatrix?',0, GEOM_ERR)
      else
c     FRACTIONAL?
         if (.not. opt_geom_cart_coords_get(geom, sp))
     $        call errquit('driver: geom?',0, GEOM_ERR)
      endif
      do i = 1, nvar 
         ds(i) = sp(i) - ds(i)
      enddo
      if (zcoord) then
         bohr = 0.52917715d0
         deg  = 0.52917715d0*180d0/(4d0*atan(1d0))
         call geom_zmat_sane_step(geom, ds)
         call geom_zmat_ico_scale(geom, ds, 1d0/bohr, 1d0/deg)
      endif
      call dscal(nvar, 1d0/alpha, ds, 1)
      if (odebug) then
         write(6,*) ' Actual ds '
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
      end
      logical function opt_geom_cart_coords_get(geom, coords)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
      integer geom
      double precision coords(*)
c
c     If we are doing a periodic system and not using internals
c     then we want the fractional coordinates.  Otherwise cartesian.
c
c     If we are including stress append the amatrix
c
      if (.not. geom_cart_coords_get(geom, coords))
     $     call errquit('driver: geom cart?',0, GEOM_ERR)
c
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_cart_to_frac(geom, coords))
     $           call errquit('driver: cart_to_frac?',0, GEOM_ERR)
      endif
c
      if (ostress) then
         if (.not. geom_amatrix_get(geom, coords(3*nat_real+1)))
     $        call errquit('driver: failed to get amatrix',0,0)
      endif
      if (ostress2) then
         if (.not. geom_lattice_get(geom, coords(3*nat_real+1)))
     $        call errquit('driver: failed to get lattice',0,0)
      endif
c
      opt_geom_cart_coords_get = .true.
c
      end
      logical function opt_geom_cart_coords_set(geom, coords)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
      integer geom
      double precision coords(*)
c
c     If we are doing a periodic system and not using internals
c     then we want the fractional coordinates.  Otherwise cartesian.
c
      logical geom_amatrix_set
      external geom_amatrix_set
c
      if (ostress) then
         if (.not. geom_amatrix_set(geom, coords(3*nat_real+1)))
     $        call errquit('driver: failed to set amatrix',0,0)
      endif
      if (ostress2) then
         if (.not. geom_lattice_set(geom, coords(3*nat_real+1)))
     $        call errquit('driver: failed to set lattice',0,0)
      endif

      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_frac_to_cart(geom, coords))
     $           call errquit('driver: frac_to_cart?',0,0)
      endif
      if (.not. geom_cart_coords_set(geom, coords))
     $     call errquit('driver: geom cart?',0,0)
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_cart_to_frac(geom, coords))
     $        call errquit('driver: cart_to_frac?',0,0)
      endif
c
      opt_geom_cart_coords_set = .true.
c
      end
      subroutine opt_get_grad(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb, geom
      character*32 theory
c
c     Get the gradient.
c     
c     If the optimization is supposed to be happening in fractional
c     coordinates convert the gradients from cartesians.
c
c     If we are including stress append the cell param gradients
c     
      logical geom_grad_cart_to_frac,value
c
      if (.not. rtdb_get(rtdb, 'task:gradient', mt_dbl, ncart,
     $     gx)) call errquit('driver: could not get gradient',0,0)
      if (isystype .ne. 0) then
         if (.not. geom_grad_cart_to_frac(geom, gx))
     $        call errquit('driver: cart_to_frac?',0,0)
      end if
      if (ostress) then
         if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $   call errquit('driver: stress theory not specified',0,RTDB_ERR)
         if (theory.eq.'pspw') then
          if (.not. rtdb_get(rtdb, 'pspw:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('driver: could not get stress',0,0)
         else if (theory.eq.'band') then
          if (.not. rtdb_get(rtdb, 'band:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('driver: could not get stress',0,0)
         else if (theory.eq.'paw') then
          if (.not. rtdb_get(rtdb, 'paw:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('driver: could not get stress',0,0)
         else
           call errquit('driver: no stress in theory',0,RTDB_ERR)
         end if
         if (isystype.eq.2) then
            !gx(3*nat_real+1) = 0.0
            !gx(3*nat_real+2) = 0.0
            gx(3*nat_real+3) = 0.0
            !gx(3*nat_real+4) = 0.0
            !gx(3*nat_real+5) = 0.0
            gx(3*nat_real+6) = 0.0
            gx(3*nat_real+7) = 0.0
            gx(3*nat_real+8) = 0.0
            gx(3*nat_real+9) = 0.0
         end if
         if (isystype.eq.1) then
            !gx(3*nat_real+1) = 0.0
            gx(3*nat_real+2) = 0.0
            gx(3*nat_real+3) = 0.0
            gx(3*nat_real+4) = 0.0
            gx(3*nat_real+5) = 0.0
            gx(3*nat_real+6) = 0.0
            gx(3*nat_real+7) = 0.0
            gx(3*nat_real+8) = 0.0
            gx(3*nat_real+9) = 0.0
         end if
      endif

      if (ostress2) then
         if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $   call errquit('driver: stress theory not specified',0,RTDB_ERR)
         if (theory.eq.'pspw') then
          if (.not. rtdb_get(rtdb, 'pspw:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('driver: could not get stress',0,0)
         else if (theory.eq.'band') then
          if (.not. rtdb_get(rtdb, 'band:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('driver: could not get stress',0,0)
         else if (theory.eq.'paw') then
          if (.not. rtdb_get(rtdb, 'paw:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('driver: could not get stress',0,0)
         else
           call errquit('driver: no stress in theory',0,RTDB_ERR)
         end if
         if (isystype.eq.2) then
            !gx(3*nat_real+1) = 0.0d0
            !gx(3*nat_real+2) = 0.0d0
            gx(3*nat_real+3) = 0.0d0
            gx(3*nat_real+4) = 0.0d0
            gx(3*nat_real+5) = 0.0d0
            !gx(3*nat_real+6) = 0.0d0
         end if
         if (isystype.eq.1) then
            !gx(3*nat_real+1) = 0.0d0
            gx(3*nat_real+2) = 0.0d0
            gx(3*nat_real+3) = 0.0d0
            gx(3*nat_real+4) = 0.0d0
            gx(3*nat_real+5) = 0.0d0
            gx(3*nat_real+6) = 0.0d0
         end if
      endif

c
      end
      subroutine driver_check_grad(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb,geom
c
      integer i,ia,k
      double precision e0, ep, em, grad, hess, step
      double precision step0
c
      logical task_gradient, task_energy, opt_geom_cart_coords_get,
     $     opt_geom_cart_coords_set
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) ' Checking the gradient '
         write(6,*)
         call util_flush(6)
      end if
c
      if (.not. geom_rtdb_store(rtdb,geom,'checkgrad')) 
     $     call errquit('driver: failed saving geom for checkgrad',0,0)
c
c     Check if fract_cart is OK
c
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. opt_geom_cart_coords_get(geom,sp))
     $      call errquit('driver: failed get coords for checkgrad',0,0)
         write(6,*) ' These are the fractional coords'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
         if (.not. geom_frac_to_cart(geom, sp))
     $       call errquit('driver: frac_to_cart?',0,0)
         write(6,*) ' These are the cartesian coords'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
         if (.not. geom_cart_to_frac(geom, sp))
     $           call errquit('driver: cart_to_frac?',0,0)
         write(6,*) ' These are the fractional coords regenerated'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
c
         if (.not. opt_geom_cart_coords_set(geom,sp))
     $      call errquit('driver: failed set coords for checkgrad',0,0)
         if (.not. opt_geom_cart_coords_get(geom,sp))
     $      call errquit('driver: failed get coords for checkgrad',0,0)
         write(6,*) ' These are the fractional coords after set'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
c
         if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &      call errquit('hnd_opt: failed storing geometry',0,0)
         if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &      call errquit('hnd_opt: no geometry ', 911,0)
         if (.not. opt_geom_cart_coords_get(geom,sp))
     $      call errquit('driver: failed get coords for checkgrad',0,0)
         write(6,*) ' These are the fractional coords after st/lo'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
      end if
c
      if (.not. task_gradient(rtdb))
     $     call errquit('driver: task_gradient failed',0,0)
      call opt_get_grad(rtdb,geom) ! Into gx

      if (ostress2) then
        call sym_lattice_stress2(geom,gx(nat_real*3+1))
      endif

      if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, e0))
     $     call errquit('driver: could not get energy',0,0)

      if (.not.rtdb_get(rtdb,"driver:checkgradstep",mt_dbl,1,step0)) 
     >  step0 = 0.01
c
      if (ga_nodeid() .eq. 0) 
     >   write(6,8) "comp","anal-g","num-g","error","hess","step","de" 
 8       format(1x,'checkgrad ',A5,4A16,2x,2A12)
c
c      do i = 1, ncart
      i = 0
      do ia = 1, ncart/3
      do k = 1, 3
         i = i +1
         if(oactive(ia)) then
c
         step = step0
c
 10      continue
         if (.not. geom_rtdb_load(rtdb,geom,'checkgrad')) 
     $        call errquit('driver: failed load geom for checkgrad',0,0)
         if (.not. opt_geom_cart_coords_get(geom,sp))
     $    call errquit('driver: failed get coords for checkgrad',0,0)
         sp(i) = sp(i) + step
         if (.not. opt_geom_cart_coords_set(geom,sp))
     $    call errquit('driver: failed set coords for checkgrad',0,0)
         if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &    call errquit('hnd_opt: failed storing geometry',0,0)
         if (.not. task_energy(rtdb)) call errquit
     $        ('driver_energy_step: energy failed',0,0)
         if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, ep))
     $     call errquit('driver: could not get energy',0,0)
         if (abs(ep-e0) .lt. 1e-6) then
            write(6,*) ' Increasing the step ', ep-e0, step
            step = step*10.0d0
            goto 10
         else if (abs(ep-e0) .gt. 1e-2) then
            write(6,*) ' Decreasing the step ', ep-e0, step
            step = step/3.0
            goto 10
         end if
c
         if (.not. geom_rtdb_load(rtdb,geom,'checkgrad')) 
     $        call errquit('driver: failed load geom for checkgrad',0,0)
         if (.not. opt_geom_cart_coords_get(geom,sp))
     $     call errquit('driver: failed get coords for checkgrad',0,0)
         sp(i) = sp(i) - step
         if (.not. opt_geom_cart_coords_set(geom,sp))
     $     call errquit('driver: failed set coords for checkgrad',0,0)
         if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &        call errquit('hnd_opt: failed storing geometry',0,0)
         if (.not. task_energy(rtdb)) call errquit
     $        ('driver_energy_step: energy failed',0,0)
         if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, em))
     $        call errquit('driver: could not get energy',0,0)
c
         grad = (ep - em) / (2.d0*step)
         hess = (ep+em-2.0d0*e0) / (2.d0*step*step)
c
         if (ga_nodeid() .eq. 0) 
     >     write(6,7) i, gx(i), grad,abs(gx(i)-grad), hess, step, ep-e0
 7       format(1x,'checkgrad ',i5,2f16.8,e16.8,f16.8,2x,2f12.6)
         end if
      enddo
      enddo
c
      if (.not. geom_rtdb_load(rtdb,geom,'checkgrad')) 
     $     call errquit('driver: failed load geom for checkgrad',0,0)
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &     call errquit('hnd_opt: failed storing geometry',0,0)
c
      end


      subroutine HND_OPT_HSS_INIT_LATTICE(HSSQ,nvar,isystype,
     $           ostress,ostress2,geom)
      implicit none
C     passed in arguments
      integer nvar
      double precision HSSQ(nvar,nvar)
      integer isystype
      logical ostress,ostress2
      integer geom
c     includes and common blocks
#include "errquit.fh"
c     local variables
      double precision volume, bulk, lattice(3,2)
      double precision amat(3,3)
      integer i,j,k,l,kat,lat,iat,jat
      integer nat_real
      integer grp_num
      double precision dum

C     functions that get called
      logical geom_volume_get,geom_lattice_get,geom_amatrix_get
      external geom_volume_get,geom_lattice_get,geom_amatrix_get
      integer geom_get_group_number
      external geom_get_group_number
c

      if (.not. geom_amatrix_get(geom, amat))
     $           call errquit('geom_amatrix_get failed:', 0, GEOM_ERR)
      grp_num=geom_get_group_number(geom)
      if (ostress) then
            bulk=5
            nat_real = (nvar/3)-3
            if (.not. geom_volume_get(geom,volume))
     $         call errquit('geom_volume_get failed',0, GEOM_ERR)
            do iat = 1, 3
             do jat = 1, 3
              do i = 1, 3
               do j = 1, 3
                 if((iat .eq. jat) .and. (i .eq. j))then
!  The next lines are correct if Energy is porportional to PV, and we
!  set second-derivative of Pressure term to 10.0, and pick a value for
!  the bulk modulus (note that diamond has a bulk modulus of 442 GPa,
!  which is converts to Hartrees/bohr^3 with 0.000034, thus giving
!  0.015.  We use a much larger value, because we are not near equilibrium.
                  if (i.eq.1) Then
                     k = 2
                     l = 3
                  else if (i.eq.2) Then
                     k = 1
                     l = 3
                  else
                     k = 1
                     l = 2
                  endif
                  if (iat.eq.1) Then
                     kat = 2
                     lat = 3
                  else if (iat.eq.2) Then
                     kat = 1
                     lat = 3
                  else
                     kat = 1
                     lat = 2
                  endif
                  HSSQ((iat+nat_real-1)*3 + i, (jat+nat_real-1)*3+j) =
     $                1.0 +
     $                (2*bulk/volume)*abs(amat(k,kat)*amat(l,lat)-
     $                                    amat(k,lat)*amat(l,kat))
!  The next lines are correct if Energy is porportional to Volume
                 else if ((iat .eq. jat) .or. (i .eq. j)) then
                   HSSQ((iat+nat_real-1)*3+i,(jat+nat_real-1)*3+j)=0.0
                 else
!  The next lines are correct if Energy is porportional to Volume
                  k = (1+2+3)-(i+j) ! the other x,y,z coordinate
                  kat = (1+2+3) - (iat+jat) ! the other 'atom'
                  HSSQ((iat+nat_real-1)*3 + i, (jat+nat_real-1)*3+j) = 
     $                   (amat(k,kat))*0.001
!  I think the above term should be symmetry correct to start with. KRG
                 end if
               end do
              end do
             end do
           end do
           if (isystype.eq.2) then
              do i=1,9
                 do j=1,9
                    if (i.ne.j) then
                       HSSQ((nat_real)*3+i,(nat_real)*3+j) = 0.0d0
                    end if
                 end do
              end do
           end if
           if (isystype.eq.1) then
              do i=1,9
                 do j=1,9
                    if (i.ne.j) then
                       HSSQ((nat_real)*3+i,(nat_real)*3+j) = 0.0d0
                    end if
                 end do
              end do
           end if

      else if (ostress2) then
           nat_real = (nvar/3)-2
           if (.not. geom_lattice_get(geom, lattice))
     $       call errquit('driver: failed to get lattice',0,0)
           if (.not. geom_volume_get(geom,volume))
     $         call errquit('geom_volume_get failed',0, GEOM_ERR)
           do i = 1, 3
            do j = 1, 3
              ! this is right for E=PV at P=0, basically to give cross-terms
              bulk=0.0001
              HSSQ((nat_real)*3 + i, (nat_real)*3+j) = 
     $               -volume*bulk/abs(lattice(i,1)*lattice(j,1))
              if (i .ne. j) then
                   HSSQ((nat_real)*3 + i, (nat_real)*3+j) = 
     $               2*HSSQ((nat_real)*3 + i, (nat_real)*3+j)
              endif
              ! The above is all negative, and the off-diagonal is largest
              ! So, mostly use term that assumes the systems is harmonic
              ! bonds basically, so it scales inversly with system size.  Add 1 
              ! to avoid terrible things with small values
              if (i.eq.j) then
                 HSSQ((nat_real)*3 + i, (nat_real)*3+j) =
     $                  1 + 0.1d0*volume/abs(lattice(i,1)**2)
              endif
            end do
           end do
           do i = 1, 3
             j = i ! No cross terms-Angles are only the same when fixed
             HSSQ((nat_real+1)*3 + i, (nat_real+1)*3+j) = 300.0 ! Angles are stiff
           end do
!          Symmetry of lengths
           if (grp_num.lt.75) then
!            Triclinic Monoclinic Orthorhombic
!            no lengths are the same
           else if (grp_num.ge.75 .and. grp_num .lt. 195) then
!            Tetragonal Trigonal Hexagonal
!            First two lengths are same
             dum =(HSSQ((nat_real)*3 + 1, (nat_real)*3+1) +
     $             HSSQ((nat_real)*3 + 2, (nat_real)*3+2))/2
             HSSQ((nat_real)*3 + 1, (nat_real)*3+1) = dum
             HSSQ((nat_real)*3 + 2, (nat_real)*3+2) = dum 
             dum =(HSSQ((nat_real)*3 + 1, (nat_real)*3+2) +
     $             HSSQ((nat_real)*3 + 2, (nat_real)*3+1))/2
             dum = max(dum,0.05d0) ! at least some coupling
             HSSQ((nat_real)*3 + 1, (nat_real)*3+2) = dum
             HSSQ((nat_real)*3 + 2, (nat_real)*3+1) = dum
           else if ((grp_num.ge.195).and.(grp_num.lt.231)) then
!            Cubic
!            All three lengths are the same
             dum = 0
             do i = 1, 3
              dum=dum+HSSQ((nat_real)*3+i,(nat_real)*3+i)
             enddo
             dum = dum/3
             do i = 1, 3
              HSSQ((nat_real)*3+i,(nat_real)*3+i) = dum
             enddo
             dum = 0
             do i = 1, 3
              do j = 1, 3
               if(i.ne.j)then
                dum=dum+HSSQ((nat_real)*3+i,(nat_real)*3+j)
               endif
              enddo
             enddo
             dum = dum/6
             dum = max(dum,0.05d0) ! at least some coupling
             do i = 1, 3
              do j = 1, 3
               if (i.ne.j) then
                HSSQ((nat_real)*3+i,(nat_real)*3+j) = dum
               endif
              enddo
             enddo
           endif
           if (isystype.eq.2) then
              do i=1,3
                 do j=1,3
                 if ((i.ne.j).and.((i.eq.3).or.(j.eq.3))) then
                    HSSQ((nat_real)*3+i,(nat_real)*3+j) = 0.0d0
                 end if
                 end do
              end do
           end if
           if (isystype.eq.1) then
              do i=1,3
                 do j=1,3
                 if ((i.ne.1).or.(j.ne.1)) then
                    HSSQ((nat_real)*3+i,(nat_real)*3+j) = 0.0d0
                 end if
                 end do
              end do
           end if
      endif
      return
      end

      subroutine hnd_opt_hss_init_cart(HSSQ,odebug,isystype,geom,nvar)
C
C     This code is a copy of Kurt Glaesemann's Cartesian Hessian
C     code that was added to GAMESS file statpt.src in:
C 29 MAR 06 - KRG - HESGES: REPLACE CARTESIAN GUESS
C     Code copied in February 2011 by Kurt into NWChem.
C
C     This guess just connects all the atoms in the systems
C     with bonds, like a gigantic over-determined Z-matrix.
C     This has proven to work very well in GAMESS.
C
C     This code differs from GAMESS code in that we add the
C     xscale term to reduce inter-atom terms, and it supports
C     periodic systems, and it add distance cutoff.
C
C This routine assumes that the Hessian is zero on entry
C
      implicit none
C     Commons and paramters
#include "nwc_const.fh"
#include "errquit.fh"
      double precision TOANGS
      integer MXATOM
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (MXATOM=nw_max_atom)
      integer ir, iw
      COMMON/HND_IOFILE/IR,IW
      integer nuc
      COMMON/HND_MOLNUC/NUC(MXATOM)
      double precision c,zan
      integer nat
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
C     Passed in variables
      integer isystype,geom,nvar
      double precision HSSQ
      DIMENSION HSSQ(nvar,nvar)
      LOGICAL odebug
C     Data structures
      integer nucz
      double precision rcov(118)
      DATA (RCOV(NUCZ),NUCZ=1,2)/0.30D+00,0.46D+00/
      DATA (RCOV(NUCZ),NUCZ=3,10)
     1     /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     2     0.70D+00,0.66D+00,0.58D+00,0.67D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     1     /1.66D+00,1.36D+00,1.25D+00,1.17D+00,
     2     1.10D+00,1.04D+00,0.99D+00,0.96D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     1     /2.03D+00,1.74D+00,
     2     1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     3     1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     4     1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.17D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     1     /2.22D+00,1.92D+00,
     2     1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     3     1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     4     1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,1.31D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     1     /2.35D+00,1.98D+00,
     2     1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D0,
     3     1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D0,
     4     1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     5     1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     6     1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,1.42D+00/
      DATA (RCOV(NUCZ),NUCZ=87,99)
     1    / 2.18D+00, 1.73D+00, 1.53D+00, 1.43D+00, 1.38D+00, 1.34D+00,
     2      1.36D+00, 1.35D+00, 1.35D+00, 1.36D+00, 1.39D+00, 1.40D+00,
     3      1.40D+00/
      DATA (RCOV(NUCZ),NUCZ=100,118)
     1    /1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00,
     2     1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, 1.29D+00, 1.28D+00,
     3     1.21D+00, 1.22D+00, 1.36D+00, 1.43D+00, 1.62D+00, 1.75D+00,
     4     1.65D+00, 1.57D+00 /
C     Variables
      double precision LONG, TOOLONG, TOOLONG2
      double precision rx,ry,rz,rr,rab,rxi,ryi,rzi,cab
      double precision BONDH(3,3),HESS2
      double precision amat(3,3),tmp(3,3),amatangs(3,3)
      integer lowx,highx,lowy,highy,lowz,highz
      integer IPOST, JPOST, II, JJ, KK, I, J, ix, iy, iz
      integer iza,izb
      logical badz
      character*16 atom_tags(mxatom)
      double precision xscale
C     External function calls
      logical geom_amatrix_get,geom_cart_get
      external geom_amatrix_get,geom_cart_get
C     
C    This is a fortran statement functions
C
      BADZ(NUCZ) = NUCZ.LE.0.OR.NUCZ.GT.118
C     
      IF(odebug) THEN
         WRITE(6,*) "Generating Cartesian Hessian Guess"
      ENDIF
      if (.not. geom_cart_get(geom, nat, atom_tags, c, zan))
     &      call errquit('Call to geom_cart_get failed',geom,GEOM_ERR)
      do i = 1, nat
         nuc(i) = nint(zan(i))
      enddo
C
C     Inter-atom terms can mess up periodic systems.
C     Scale non-periodic just because
C     Should note that this type of matrix guess without scaling
C       will have at least one zero eigenvalue, and in general should
C       probably have 3 for translations and 3 for rotations (the rotations
C       are not zero with periodic boundaries).  The use of xscale shifts
C       the off-diagonal down, which removes the zeros.  Since we use
C       an inverse Hessian update scheme, that might explain why this
C       scaling is needed - the inverse of zero is not easy to calculate.
C
      IF (isystype .ne. 0) then
        xscale = 0.8d0
      else
        xscale = 0.9d0
      endif
!     Any periodic images to look for?
!     We only go +/- one unit cell
      IF (isystype .eq. 0) then      ! Not periodic
        lowx = 0
        highx= 0
        lowy = 0
        highy= 0
        lowz = 0
        highz= 0
      else if (isystype .eq. 1) then ! Periodic in X
        lowx =-1
        highx= 1
        lowy = 0
        highy= 0
        lowz = 0
        highz= 0
      else if(isystype .eq. 2) then ! Periodic in XY
        lowx =-1
        highx= 1
        lowy =-1
        highy= 1
        lowz = 0
        highz= 0
      else if(isystype .eq. 3) then ! Periodic in XYZ
        lowx =-1
        highx= 1
        lowy =-1
        highy= 1
        lowz =-1
        highz= 1
      end if
      if (isystype .ne. 0) then
       if (.not. geom_amatrix_get(geom, amat))
     $           call errquit('geom_amatrix_get: a', 0, GEOM_ERR)
        do i=1,3
         do j=1,3
          amatangs(i,j) = amat(i,j)*toangs
         end do
        end do
      else
        do i=1,3
         do j=1,3
          amat(i,j) = 0
          amatangs(i,j) = 0
         end do
        end do
      endif
      LONG = 2.5 ! 0.0028 = 0.3601D+00 * EXP(-1.944D+00*( LONG ))
      TOOLONG = LONG + 2*2.2d0 ! The 2.2 is from BADZ below
      TOOLONG2=TOOLONG**2
      DO I=1,NAT
       DO J=1,I-1
         IZA = NUC(I)
         IZB = NUC(J)
         IF(BADZ(IZA)) THEN
           CAB = 2.2D0
         ELSE
           CAB = RCOV(IZA)
         END IF
         IF(BADZ(IZB)) THEN
           CAB = CAB + 2.2D0
         ELSE
           CAB = CAB + RCOV(IZB)
         END IF
!        Vector from I to J inside the box - in Angstroms
         RXI = (C(1,I)-C(1,J))*TOANGS
         RYI = (C(2,I)-C(2,J))*TOANGS
         RZI = (C(3,I)-C(3,J))*TOANGS
         DO ix = lowx,highx
          do iy = lowy,highy
           do iz = lowz,highz
         IF(ix.eq.0 .and. iy.eq.0 .and. iz.eq.0 .and. i.eq.j) goto 999
         RX = (RXI + ix*amatangs(1,1)+iy*amatangs(1,2)+iz*amatangs(1,3))
         RY = (RYI + ix*amatangs(2,1)+iy*amatangs(2,2)+iz*amatangs(2,3))
         RZ = (RZI + ix*amatangs(3,1)+iy*amatangs(3,2)+iz*amatangs(3,3))
         RR = RX*RX + RY*RY + RZ*RZ
         IF(RR .GT. TOOLONG2) GOTO 999
         RAB = SQRT(RR)
         RAB = RAB - CAB
         IF (RAB .GT. LONG) GOTO 999
C        HESS2 IS HESSIAN WITH RESPECT TO RAB (SEE FISCHER/ALMLOF)
         HESS2 = 0.3601D+00 * EXP(-1.944D+00*RAB)
C        NOW CONVERT TO CARTESIANS BY CHAIN RULE
         BONDH(1,1) = (RX*RX*HESS2)/RR
         BONDH(1,2) = (RX*RY*HESS2)/RR
         BONDH(1,3) = (RX*RZ*HESS2)/RR
         BONDH(2,2) = (RY*RY*HESS2)/RR
         BONDH(2,3) = (RY*RZ*HESS2)/RR
         BONDH(3,3) = (RZ*RZ*HESS2)/RR
         BONDH(2,1) = BONDH(1,2)
         BONDH(3,1) = BONDH(1,3)
         BONDH(3,2) = BONDH(2,3)
C
C        Convert to fractional coordinates
C
         IF (isystype .ne. 0) then
          do kk = 1, 3
           do ii = 1, 3
            tmp(ii,kk) = 0.0d0
           end do
          end do
          do kk = 1, 3
           do ii = 1, 3
            do jj = 1, 3
             tmp(ii,kk) = tmp(ii,kk) + amat(jj,ii)*bondh(jj,kk)
            end do
           end do
          end do
          do kk = 1, 3
           do ii = 1, 3
            bondh(ii,kk) = 0.0d0
           end do
          end do
          do kk = 1, 3
           do ii = 1, 3
            do jj = 1, 3
             bondh(ii,kk) = bondh(ii,kk) + amat(jj,ii)*tmp(kk,jj)
            end do
           end do
          end do
         end if
         IPOST = (I - 1)*3
         JPOST = (J - 1)*3
         DO II = 1,3
          DO JJ = 1,3
           HSSQ(IPOST+II,IPOST+JJ)=HSSQ(IPOST+II,IPOST+JJ)+BONDH(II,JJ)
           HSSQ(JPOST+II,JPOST+JJ)=HSSQ(JPOST+II,JPOST+JJ)+BONDH(II,JJ)
           HSSQ(IPOST+II,JPOST+JJ)=HSSQ(IPOST+II,JPOST+JJ)-BONDH(II,JJ)
     &          *xscale
           HSSQ(JPOST+II,IPOST+JJ)=HSSQ(JPOST+II,IPOST+JJ)-BONDH(II,JJ)
     &          *xscale
          END DO
         END DO
  999      end do
          end do
         end do
       END DO
      END DO
C     Add some to diagonal for stability, just in case of lone floating atom
      IF (isystype .ne. 0) then
       do i = 1, nat
          ipost = (i - 1)*3
          do j = 1, 3
             HSSQ(ipost+j,ipost+j) = max(HSSQ(ipost+j,ipost+j),
     &               0.4d0*(amat(1,j)**2+amat(2,j)**2+amat(3,j)**2))
          end do
        end do
      ELSE
        DO I=1,3*NAT
          HSSQ(I,I)=MAX(HSSQ(I,I),0.2D0)
        END DO
      ENDIF

      IF(odebug) THEN
         WRITE(6,*) "Cartesian Hessian Guess",NVAR
         CALL HND_PRSQ(HSSQ,NVAR,NVAR,NVAR)
      ENDIF
      END

!************************************************************************

      subroutine sym_lattice_stress2(geom,lstress)
      implicit none
#include "geom.fh"
!     include files
!     calling args
      integer geom
      double precision lstress(0:5)
!     local variables
      integer i
      integer grp_num
      integer crystal
      double precision lstress_sym(0:5),radtodeg,lat(6)
      integer Triclinic,Monoclinic,Orthorhombic
      integer Tetragonal,Trigonal,Hexagonal,Cubic
      parameter (Triclinic=1,Monoclinic=2,Orthorhombic=3)
      parameter (Tetragonal=4,Trigonal=5,Hexagonal=6,Cubic=7)
      logical value,is_conventional
c      logical is90(3),g120,isab,isac,isbc
      logical is_A,is_C,is_R,is_F,is_I
      integer isystype
      character*1 ctr_type

!     Called functions
      integer geom_get_group_number
      external geom_get_group_number
      logical  geom_is_conventional
      external geom_is_conventional
      character*1 geom_primitive_center
      external    geom_primitive_center
!
      if (.not. geom_lattice_get(geom, lat))
     $        call errquit('driver: failed to get lat',0,0)


      radtodeg = 180.0d0/(4.0d0*datan(1.0d0))
      value = geom_systype_get(geom,isystype)
      grp_num=geom_get_group_number(geom)
      is_conventional = geom_is_conventional(geom)
      ctr_type = geom_primitive_center(geom)
      is_conventional = ctr_type.eq.'x'
      is_A =  ctr_type.eq.'A'
      is_C =  ctr_type.eq.'C'
      is_F =  ctr_type.eq.'F'
      is_I =  ctr_type.eq.'I'
      is_R =  ctr_type.eq.'R'
c      write(*,*) "is_conventional,ctr_type=",is_conventional,ctr_type,
c     >  is_A,is_C,is_F,is_I,is_R

c      value = geom_lattice_get(geom, lat)
c      is90(1) = (dabs((lat(4)*radtodeg)-90.0d0).lt.1.0d-3)
c      is90(2) = (dabs((lat(5)*radtodeg)-90.0d0).lt.1.0d-3)
c      is90(3) = (dabs((lat(6)*radtodeg)-90.0d0).lt.1.0d-3)
c      g120    = (dabs((lat(6)*radtodeg)-120.0d0).lt.1.0d-3)
c      isab    = (dabs(lat(1)-lat(2)).lt.1.0d-3)
c      isac    = (dabs(lat(1)-lat(3)).lt.1.0d-3)
c      isbc    = (dabs(lat(2)-lat(3)).lt.1.0d-3)

      if (grp_num.lt.3)                         crystal = Triclinic
      if ((grp_num.ge.3  ).and.(grp_num.lt.16 ))crystal = Monoclinic
      if ((grp_num.ge.16 ).and.(grp_num.lt.75 ))crystal = Orthorhombic
      if ((grp_num.ge.75 ).and.(grp_num.lt.143))crystal = Tetragonal
      if ((grp_num.ge.143).and.(grp_num.lt.168))crystal = Trigonal
      if ((grp_num.ge.168).and.(grp_num.lt.195))crystal = Hexagonal
      if ((grp_num.ge.195).and.(grp_num.lt.231))crystal = Cubic
      if (crystal.eq.Triclinic) then
        lstress_sym(0) = lstress(0)
        lstress_sym(1) = lstress(1) 
        lstress_sym(2) = lstress(2)
        lstress_sym(3) = lstress(3)
        lstress_sym(4) = lstress(4)
        lstress_sym(5) = lstress(5)

      else if (crystal.eq.Monoclinic) then
        if (is_conventional) then
           lstress_sym(0) = lstress(0)
           lstress_sym(1) = lstress(1)
           lstress_sym(2) = lstress(2)
           lstress_sym(3) = 0.0D0
           lstress_sym(4) = lstress(4)
           lstress_sym(5) = 0.0D0
        else
           !** C-centered primitive **
           lstress_sym(0) = (lstress(0)+lstress(1))/2.0d0
           lstress_sym(1) = (lstress(0)+lstress(1))/2.0d0
           lstress_sym(2) = lstress(2)
           lstress_sym(3) = lstress(3)
           lstress_sym(4) = lstress(4)
           lstress_sym(5) = lstress(5)
        end if

      else if (crystal.eq.Orthorhombic) then
        if (is_conventional) then
           lstress_sym(0) = lstress(0)
           lstress_sym(1) = lstress(1)
           lstress_sym(2) = lstress(2)
           lstress_sym(3) = 0.0D0
           lstress_sym(4) = 0.0D0
           lstress_sym(5) = 0.0D0
        else
           !** C-centered primitive **
c           if (isab.and.is90(1).and.is90(2).and.(.not.is90(3))) then
           if (is_C) then
              lstress_sym(0) = (lstress(0)+lstress(1))/2.0d0
              lstress_sym(1) = (lstress(0)+lstress(1))/2.0d0
              lstress_sym(2) = lstress(2)
              lstress_sym(3) = 0.0d0
              lstress_sym(4) = 0.0d0
              lstress_sym(5) = lstress(5)
           !** A-centered primitive ***
c           else if (isbc.and.
c     >              (.not.is90(1)).and.
c     >              is90(2).and.
c     >              is90(3).and.isbc) then
           else if (is_A) then
              lstress_sym(0) = lstress(0)
              lstress_sym(1) = (lstress(1)+lstress(2))/2.0d0
              lstress_sym(2) = (lstress(1)+lstress(2))/2.0d0
              lstress_sym(3) = lstress(3)
              lstress_sym(4) = 0.0d0
              lstress_sym(5) = 0.0d0
           !** I-centered primitive **
           !else if (isab.and.isac.and.isbc) then
           else if (is_I) then
              lstress_sym(0) = (lstress(0)+lstress(1)+lstress(2))/3.0d0
              lstress_sym(1) = (lstress(0)+lstress(1)+lstress(2))/3.0d0
              lstress_sym(2) = (lstress(0)+lstress(1)+lstress(2))/3.0d0
              lstress_sym(3) = lstress(3)
              lstress_sym(4) = lstress(4)
              lstress_sym(5) = lstress(5)
           !** F-centered primitive - everything changes**
           else
              lstress_sym(0) = lstress(0)
              lstress_sym(1) = lstress(1)
              lstress_sym(2) = lstress(2)
              lstress_sym(3) = lstress(3)
              lstress_sym(4) = lstress(4)
              lstress_sym(5) = lstress(5)
           end if
        end if

      else if (crystal.eq.Tetragonal) then
        if (is_conventional) then
           lstress_sym(0) = (lstress(0)+lstress(1))/2.0D0
           lstress_sym(1) = (lstress(0)+lstress(1))/2.0D0
           lstress_sym(2) = lstress(2)
           lstress_sym(3) = 0.0D0
           lstress_sym(4) = 0.0D0
           lstress_sym(5) = 0.0D0
        else
           !** I-centered primitive **
           lstress_sym(0) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
           lstress_sym(1) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
           lstress_sym(2) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
           lstress_sym(3) = (lstress(3)+lstress(4))/2.0d0
           lstress_sym(4) = (lstress(3)+lstress(4))/2.0d0
           lstress_sym(5) = lstress(5)
        end if

      else if (crystal.eq.Trigonal) then
        if (is_conventional) then
           lstress_sym(0) = (lstress(0)+lstress(1))/2.0D0
           lstress_sym(1) = (lstress(0)+lstress(1))/2.0D0
           lstress_sym(2) = lstress(2)
           lstress_sym(3) = 0.0D0
           lstress_sym(4) = 0.0D0
           lstress_sym(5) = 0.0D0
        else
           !** R-centered primitive **
           lstress_sym(0) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
           lstress_sym(1) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
           lstress_sym(2) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
           lstress_sym(3) = (lstress(3)+lstress(4)+lstress(5))/3.0D0
           lstress_sym(4) = (lstress(3)+lstress(4)+lstress(5))/3.0D0
           lstress_sym(5) = (lstress(3)+lstress(4)+lstress(5))/3.0D0
        end if
      else if (crystal.eq.Hexagonal) then
        if (is_conventional) then
           lstress_sym(0) = (lstress(0)+lstress(1))/2.0D0
           lstress_sym(1) = (lstress(0)+lstress(1))/2.0D0
           lstress_sym(2) = lstress(2)
           lstress_sym(3) = 0.0D0
           lstress_sym(4) = 0.0D0
           lstress_sym(5) = 0.0D0
        else
           !** Should never happen **
           lstress_sym(0) = lstress(0)
           lstress_sym(1) = lstress(1)
           lstress_sym(2) = lstress(2)
           lstress_sym(3) = lstress(3)
           lstress_sym(4) = lstress(4)
           lstress_sym(5) = lstress(5)
        end if
      else if (crystal.eq.Cubic) then
        lstress_sym(0) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
        lstress_sym(1) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
        lstress_sym(2) = (lstress(0)+lstress(1)+lstress(2))/3.0D0
        lstress_sym(3) = 0.0D0
        lstress_sym(4) = 0.0D0
        lstress_sym(5) = 0.0D0
        !** I-centered primitive 0 equal dist, all angs=109.471**
        !** F-centered primitive - equal dist, all angs=60**
      endif
      do i = 0, 5
        lstress(i) = lstress_sym(i)
      end do
      if (isystype.eq.2) then
         !lstress(0) = 0.0d0
         !lstress(1) = 0.0d0
         lstress(2) = 0.0d0
         lstress(3) = 0.0d0
         lstress(4) = 0.0d0
         !lstress(5) = 0.0d0
      end if
      if (isystype.eq.1) then
         !lstress(0) = 0.0d0
         lstress(1) = 0.0d0
         lstress(2) = 0.0d0
         lstress(3) = 0.0d0
         lstress(4) = 0.0d0
         lstress(5) = 0.0d0
      end if
      return
      end

!************************************************************************

      subroutine sym_lattice_stress(geom,astress)
      implicit none
#include "geom.fh"
!     include files
!     calling args
      integer geom
      double precision astress(3,3)
!     local variables
      integer i, j
      double precision lstress(6),lattnew(6)
      double precision amat_old(3,3),amat_new(3,3)
      integer Triclinic,Monoclinic,Orthorhombic
      integer Tetragonal,Trigonal,Hexagonal,Cubic
      parameter (Triclinic=1,Monoclinic=2,Orthorhombic=3)
      parameter (Tetragonal=4,Trigonal=5,Hexagonal=6,Cubic=7)
      integer grp_num, crystal
      double precision dum
      double precision deg90,deg120
      logical value
      integer isystype
!     Called functions 
      integer geom_get_group_number
      external geom_get_group_number
!       
      value = geom_systype_get(geom,isystype)
      grp_num=geom_get_group_number(geom)
!        
      if (grp_num.lt.3)                         crystal = Triclinic
      if ((grp_num.ge.3  ).and.(grp_num.lt.16 ))crystal = Monoclinic
      if ((grp_num.ge.16 ).and.(grp_num.lt.75 ))crystal = Orthorhombic
      if ((grp_num.ge.75 ).and.(grp_num.lt.143))crystal = Tetragonal
      if ((grp_num.ge.143).and.(grp_num.lt.168))crystal = Trigonal
      if ((grp_num.ge.168).and.(grp_num.lt.195))crystal = Hexagonal
      if ((grp_num.ge.195).and.(grp_num.lt.231))crystal = Cubic

      if (crystal .eq. Triclinic) return
! Get old amatrix
      if (.not. geom_amatrix_get(geom, amat_old))
     $        call errquit('driver: failed to get amatrix',0,0)
! Generate proposed amatrix
      do i = 1, 3
       do j = 1, 3
         amat_new(j,i) = amat_old(j,i) + astress(j,i)
       enddo
      enddo
! Generate proposed lattice vectors
      call xlattice_abc_abg(lattnew(1),lattnew(2),lattnew(3),lattnew(4),
     &              lattnew(5),lattnew(6),amat_new)
! Symmetrize proposed lattice
      deg90 = 2d0*atan(1d0)
      deg120= (8d0/3d0)*atan(1d0)
      if (crystal.eq.Triclinic) then
      else if (crystal.eq.Monoclinic) then
        lattnew(4) = deg90
        lattnew(6) = deg90
      else if (crystal.eq.Orthorhombic) then
        lattnew(4) = deg90
        lattnew(5) = deg90
        lattnew(6) = deg90
      else if (crystal.eq.Tetragonal) then
        dum = (lattnew(1)+lattnew(2))/2.0D0
        lattnew(1) = dum
        lattnew(2) = dum
        lattnew(4) = deg90
        lattnew(5) = deg90
        lattnew(6) = deg90
      else if (crystal.eq.Trigonal) then
        dum = (lattnew(1)+lattnew(2))/2.0D0
        lattnew(1) = dum
        lattnew(2) = dum
        lattnew(4) = deg90
        lattnew(5) = deg90
        lattnew(6) = deg120
      else if (crystal.eq.Hexagonal) then
        dum = (lattnew(1)+lattnew(2))/2.0D0
        lattnew(1) = dum
        lattnew(2) = dum
        lattnew(4) = deg90
        lattnew(5) = deg90
        lattnew(6) = deg120
      else if (crystal.eq.Cubic) then
        dum = (lattnew(1)+lattnew(2)+lattnew(3))/3.0D0
        lattnew(1) = dum
        lattnew(2) = dum
        lattnew(3) = dum
        lattnew(4) = deg90
        lattnew(5) = deg90
        lattnew(6) = deg90
      endif
! Generate symmetrized proposed amatrix
!   Apply lengths
      do i = 1, 3
        dum = sqrt(amat_new(1,i)**2 +amat_new(2,i)**2 +amat_new(3,i)**2)
        dum = lattnew(i)/dum
        if (abs(dum - 1.0d0) .gt. 1.0d-20) then
         do j = 1, 3
          amat_new(j,i) = amat_new(j,i) * dum
         end do
        endif
      end do
!   Apply Angles - Only Monoclinic has both variable angles and Symmetry
!     It also only has one angle with non-zero stress (others fixed at 90 deg)
!     Note that groups with no lattice symmetry bailed earlier
!  Enforce alpha
! TODO
!  Enforce beta
! TODO
!  Enforce gamma
! TODO
! Generate symmetrized amatrix stress
      do i = 1, 3
       do j = 1, 3
         astress(j,i) = amat_new(j,i) - amat_old(j,i)
         if (abs(astress(j,i)) .lt. 1.0d-10) astress(j,i)=0.0d0
       enddo
      enddo
      if (isystype.eq.2) then
         !astress(1,1) = 0.0d0
         !astress(2,1) = 0.0d0
         astress(3,1) = 0.0d0
         !astress(1,2) = 0.0d0
         !astress(2,2) = 0.0d0
         astress(3,2) = 0.0d0
         astress(1,3) = 0.0d0
         astress(2,3) = 0.0d0
         astress(3,3) = 0.0d0
      end if
      if (isystype.eq.1) then
         !astress(1,1) = 0.0d0
         astress(2,1) = 0.0d0
         astress(3,1) = 0.0d0
         astress(1,2) = 0.0d0
         astress(2,2) = 0.0d0
         astress(3,2) = 0.0d0
         astress(1,3) = 0.0d0
         astress(2,3) = 0.0d0
         astress(3,3) = 0.0d0
      end if
      return
      end
C>
C> @}
      subroutine driver_hss_reinit(rtdb,geom,l_hess, k_hess)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
c
      integer rtdb, geom

c
      logical geom_zmt_get_izmat, geom_zmt_get_nizmat, 
     $     geom_compute_zmatrix
      external geom_zmt_get_izmat, geom_zmt_get_nizmat,
     $     geom_compute_zmatrix
c
      double precision zero
      parameter (zero=0.0d+00)
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=nw_max_atom)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=nw_max_zmat)     
      parameter (mxcoor=nw_max_coor)
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
      double precision hscale, ascale, bscale, tscale, amat(3,3)
c
      integer l_hess, k_hess, l_zmat, k_zmat, l_izmat, k_izmat, i, j, k
      integer l_c, k_c, l_t, k_t, iat, jat, kat, lat, l, l_kurt, k_kurt
      integer l_shrhess, k_shrhess
      character*16 atom_tags(mxatom)
c
      nnzmat = nzmat
      nnzvar = nzvar
      nnvar  = nzvar
      if (.not. geom_ncent(geom,natom))
     1       call errquit('hnd_opt: geom_ncent?',911, GEOM_ERR)
c     create a chop point
      if (.not. ma_push_get(mt_dbl, 1, 'kurt',
     $     l_kurt, k_kurt)) call errquit
     $     ('driver_hss_reinit: failed allocating ',1,
     &       MA_ERR)

c
      if (.not. zcoord) then
c
c     Cartesians are easy
c
         call dfill(nvar**2, 0.0d0, dbl_mb(k_hess), 1)
         call hnd_opt_hss_init_cart(dbl_mb(k_hess),
     &                     odebug,isystype,geom,nvar)
!    Now the lattice hessian, if any
         if (isystype .ne. 0 .and. (ostress .or. ostress2)) then
            call HND_OPT_HSS_INIT_LATTICE(dbl_mb(k_hess),nvar,isystype,
     $           ostress,ostress2,geom)
            if (oprint) then 
               write(6,*) ' The initial hessian '
               call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
            end if
         end if
!        Symmetry
         if (nvar .eq. 3*natom) then
           call sym_hess_symmetrize(geom,dbl_mb(k_hess))
         else
           if (.not. ma_push_get(mt_dbl, (natom*3)**2, 'shrhess',
     $        l_shrhess, k_shrhess)) call errquit
     $        ('driver_hss_reinit: failed allocating shrunk hessian ',
     $        (natom*3)**2, MA_ERR)
           do i = 1, 3*natom
             do j = 1, 3*natom
               dbl_mb(k_shrhess+j-1 + (i-1)*(3*natom)) = 
     $         dbl_mb(k_hess   +j-1 + (i-1)*(nvar    ))
             end do
           end do
           call sym_hess_symmetrize(geom,dbl_mb(k_shrhess))
           do i = 1, 3*natom
             do j = 1, 3*natom
               dbl_mb(k_hess   +j-1 + (i-1)*(nvar    )) =
     $         dbl_mb(k_shrhess+j-1 + (i-1)*(3*natom))
             end do
           end do
         endif
      else
c
c     Internals.
c
         if (.not. ma_push_get(mt_dbl, nzvar, 'zmat',
     $        l_zmat, k_zmat)) call errquit
     $        ('driver_hss_reinit: failed allocating zmat ', nzvar,
     &       MA_ERR)
         if (.not. ma_push_get(mt_int, nzmat, 'izmat',
     $        l_izmat, k_izmat)) call errquit
     $        ('driver_hss_reinit: failed allocating izmat ', nzmat,
     &       MA_ERR)
         if (.not. geom_compute_zmatrix(geom, dbl_mb(k_zmat)))
     $        call errquit('driver_hss_reinit: zmatrix?',0, GEOM_ERR)
         if (.not. geom_zmt_get_izmat(geom,int_mb(k_izmat),nzmat))
     $        call errquit('driver_hss_reinit: izmat?',0, GEOM_ERR)
C  The code assume Cartesians in many places:
         if (ostress .or. ostress2) call errquit(
     $     'driver: lattice optimization does not work with a zmatrix',
     $     0,GEOM_ERR)
c
         if (.not. ma_verify_allocator_stuff()) call errquit('h',0,
     &       MA_ERR)
c
C     ----- INTERNALLY GENERATED (DIAGONAL) INITIAL HESSIAN -----
C
c     this routine uses /molnuc/, /zmtpar/, /molxyz/ which 
c     must be initialized before calling
c
         if (.not. geom_cart_get(geom, nat, atom_tags, 
     $        c, zan)) call errquit('driver_hss_reinit: geom?',geom,
     &       GEOM_ERR)
         do i = 1, nat
            nuc(i) = nint(zan(i))
         enddo
         nvar = nzvar
c
c     The angle scaling factors are responsible for a factor of two
c     speedup in most organics.  Clearly we need a better force field.
c
c     The chosen values have been tested on a large no. of 3-21g* SCF
c     optimizations, and a smaller no. of 6-31g* SCF calculations.
c
         if (.not. rtdb_get(rtdb,'driver:bscale',mt_dbl,1,bscale))
     $        bscale = 1d0
         if (.not. rtdb_get(rtdb,'driver:ascale',mt_dbl,1,ascale))
     $        ascale = 0.25d0
         if (.not. rtdb_get(rtdb,'driver:tscale',mt_dbl,1,tscale))
     $        tscale = 0.1d0
         if (oprint) then
            write(6,789) bscale, ascale, tscale
 789        format(' Scaling for Hessian diagonals: bonds =',f5.2,
     $           '  angles =',f5.2,'  torsions =',f5.2)
         endif
         call hnd_opt_hss_init_zmt(dbl_mb(k_hess),
     $        nzvar, dbl_mb(k_zmat), int_mb(k_izmat), .false.,
     $        bscale, ascale, tscale,odebug)
         if (odebug) then
            write(6,*) ' The initial hessian '
            call output(dbl_mb(k_hess),1,nzvar,1,nzvar,nzvar,nzvar,1)
         endif
      endif
c
c     Artificially break degeneracies so that accidentally degenerate
c     modes are split and therefore step restriction along modes
c     is well defined
c
      do i = 1, nvar
         dbl_mb(k_hess+i-1 + (i-1)*nvar) = 
     $        dbl_mb(k_hess+i-1 + (i-1)*nvar) + dble(i-1)*1e-7
      end do
c
      call geom_hnd_put_data('drv.hess', dbl_mb(k_hess), nvar*nvar)
c
c     Apply constants, constraints and overall scaling
c
 999  if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_c, k_c)) call errquit
     $     ('driver_hss_reinit: failed allocating hessian',nvar**2,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_t, k_t)) call errquit
     $     ('driver_hss_reinit: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      call geom_hnd_get_data('drv.hess', dbl_mb(k_hess), nvar*nvar)
c
c     Used to use c here ... now use p
c
      call geom_hnd_get_data('p',dbl_mb(k_c), nvar*nvar)
c
      if (odebug) then
         write(6,*) ' Initial Hessian before P'
         call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
      endif
      call dgemm('n','n',nvar,nvar,nvar,1d0,dbl_mb(k_c),nvar,
     $     dbl_mb(k_hess),nvar,0d0,dbl_mb(k_t),nvar)
      call dgemm('n','t',nvar,nvar,nvar,1d0,dbl_mb(k_t),nvar,
     $     dbl_mb(k_c),nvar,0d0,dbl_mb(k_hess),nvar)
      if (odebug) then
         write(6,*) ' Initial Hessian after P'
         call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
      endif
c
      if (nactive .ne. nat_real) then
c     
c     We are in cartesian coordinates and some have been frozen.
c     Since there is no redundancy or coupling we just need
c     to make sure that the initial Hessian does not couple
c     frozen with unfrozen variables and we are OK.
c
         do iat = 1, nat_real
            if (.not. oactive(iat)) then
               do i = 1+(iat-1)*3, iat*3
                  do j = 1, nvar
                     dbl_mb(k_hess+j-1+(i-1)*nvar) = 0d0
                     dbl_mb(k_hess+i-1+(j-1)*nvar) = 0d0
                  enddo
                  dbl_mb(k_hess+i-1+(i-1)*nvar) = 1d0
               enddo
            endif
         enddo
      endif
c
      if (.not. rtdb_get(rtdb,'driver:hscale',mt_dbl,1,hscale))
     $     hscale = 1d0
      call dscal(nvar*nvar, hscale, dbl_mb(k_hess), 1)
      if (oprint .and. hscale.ne.1d0) write(6,78) hscale
 78   format(' Scaling initial hessian by ',f6.2)
c
      call geom_hnd_put_data('drv.hess',dbl_mb(k_hess), nvar*nvar)
c
      if (.not. ma_chop_stack(l_kurt)) call errquit
     $     ('driver_hss_reinit ma corrupt',0, MA_ERR)
c
      END
