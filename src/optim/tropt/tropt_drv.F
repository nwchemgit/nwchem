      logical function tropt(rtdb)
*
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb
c     
      integer geom, geomi, i, iat
      integer istep
      logical converged, status
      double precision start    ! Tracks time used in last step
cjmc
      double precision oldx(max_cart)
      logical tropt_geom_cart_coords_get
      logical tropt_geom_cart_coords_set
cjmc
      integer required          ! Time required
      logical ophigh
      logical  tropt_converged, task_gradient, 
     $     geom_compute_zmatrix, geom_print_zmatrix
      external tropt_converged, task_gradient,
     $     geom_compute_zmatrix, geom_print_zmatrix
c
      ophigh = util_print('high', print_high)
c
c     Read input, load /ctropt/, get geometry
c
      call tropt_initialize(rtdb, geom)
c
c     Optionally check the gradient by finite difference
c
      if (ocheckgrad) call tropt_check_grad(rtdb,geom)
c
c
c     Energy and Gradient
c
      if (.not. task_gradient(rtdb))
     $    call errquit('tropt: task_gradient failed',0, GEOM_ERR)
      call tropt_get_grad(rtdb,geom) ! Into gx
      if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $    call errquit('tropt: could not get energy',0, RTDB_ERR)
c     
c     Zero the gradient associated with atoms frozen in cartesians
c     
      if (ga_nodeid().eq.0.and.ophigh) 
     $    write(6,*) 'Zeroing constrained gradient'

      if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
        do iat = 1, nat_real
          if (.not. oactive(iat)) then
            do i = 1, 3
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
        call tropt_tfgx(geom) ! g() will be set to constrained gq()
      else
        call tropt_cart_pmat(rtdb, geom)
        call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
      endif
c     
c     Initial hessian
c
      call tropt_hss_init(rtdb,geom)
c
c     Iterative loop
c
      do istep = 1, nptopt+1    ! +1 since first pass thru loop is not a step
c
c     iteration exceded
c
c
        if(istep.gt.nptopt) goto 200
c
        start = util_wallsec()
        if (oprint) write(6,1) istep-1
 1      format(/,10x,'--------',/,10x,'Step',i4,/,10x,'--------')
        if ((ga_nodeid() .eq. 0) .and.
     $       util_print('geometry',print_default)) then
           if (.not. geom_print(geom)) call errquit('tropt: geom?',0,
     $       GEOM_ERR)
        endif
c
c     Disable printing to ecce of movecs after first point
c
        call movecs_ecce_print_off()
c
c     Compute step/gradient info and print for user
c     (for current energy & gradient, and the previous alpha*step).
c
        call tropt_compute_info()
        call tropt_print(geom,istep,.false.)
c     
c     Check for convergence
c     
        if (tropt_converged(istep)) then
           converged = .true.
           goto 100
        endif
c
        redostep = .true.
c    
c     Save old energy, gradient and coordinates
C
        energyp = energy       ! Used for convergence and step restriction
        call dcopy(nvar, g, 1, gp, 1) ! Used for Hessian update
        if (.not. tropt_geom_cart_coords_get(geom, oldx))
     $      call errquit('tropt: geom?',0, geom_err) ! obtain previous cart
c
c     *** Perform the step until it is good enough ***
c
        do while (redostep)
c
c     generate a new search direction
c
          call tropt_pickstp(rtdb,geom,istep) ! fills in ds()
c     
c     take recommended step.
c     
          call tropt_take_step(rtdb, geom) ! Updates geom using ds
c
c     We have now taken a step.  Replace the approximate step taken
c     by the exact step in case update of internals was not exact.
c     
          call tropt_compute_actual_step(geom)
c
c     Energy and Gradient
c
          if (.not. task_gradient(rtdb))
     $        call errquit('tropt: task_gradient failed',0, GEOM_ERR)
          call tropt_get_grad(rtdb,geom) ! Into gx
          if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $        call errquit('tropt: could not get energy',0, RTDB_ERR)
c     
c     Zero the gradient associated with atoms frozen in cartesians
c     
          if (ga_nodeid().eq.0.and.ophigh) 
     $        write(6,*) 'Zeroing constrained gradient'
          if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
            do iat = 1, nat_real
              if (.not. oactive(iat)) then
                do i = 1, 3
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
            call tropt_tfgx(geom) ! g() will be set to constrained gq()
          else
            call tropt_cart_pmat(rtdb, geom)
            call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
          endif
c
c     update the trust radius
c
          call tropt_radius()
c
c     update the hessian
c
          call tropt_hessian_update()
c
c       reload previous fiels if necessary
c
          if (redostep) then
            call dcopy(nvar, gp, 1, g, 1) ! used for hessian update
            if (.not. tropt_geom_cart_coords_set(geom, oldx))
     $         call errquit('tropt: geom?',0, geom_err) ! reload previous cart
            if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $         call errquit('tropt: grs?',geom, RTDB_ERR)
            write(6,*) "Non good step, calculation step iteration",istep
c
          end if
c
c       End of while step calculation 
c
        enddo                     ! End of while loop
c
c     Check time before next iteration
c
         required = int(1.2d0*(util_wallsec() - start)) + 1
         if (.not. util_test_time_remaining(rtdb,required)) goto 200
c
c
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
         call tropt_print(geom,istep,.true.)
         if (.not. geom_print(geom)) call errquit
     $        ('hnd_opt_drv: geom_print?',0, GEOM_ERR)
c
c     Analyze the change in internal coordinates
c     (overwrites ds and g)
c
         if (zcoord) then
            status = rtdb_parallel(.false.)
            if (.not. geom_create(geomi, 'troptinitial'))
     &           call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
            if (.not. geom_rtdb_load(rtdb, geomi, 'troptinitial'))
     &           call errquit('hnd_opt: no initial geometry ', 911,
     &       RTDB_ERR)
            if (.not. geom_compute_zmatrix(geomi, ds)) ! Original coords
     $           call errquit('tropt_u_c_f_i: geom?',0, GEOM_ERR)
            if (.not. geom_compute_zmatrix(geom, g)) ! Final coords
     $           call errquit('tropt_u_c_f_i: geom?',0, GEOM_ERR)
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
     $           call errquit('tropt:print geom change?',geom,
     &       GEOM_ERR)
            if (.not. geom_destroy(geomi))
     $           call errquit('tropt:geomi corrupt',0, GEOM_ERR)
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
      tropt=converged
      if (converged) then
         call ecce_print_module_exit('tropt', 'ok')
      else
         call ecce_print_module_exit('tropt', 'failed')
      endif
c
      call movecs_ecce_print_on() ! Restore MO printing
      call util_print_pop
c
      call ga_sync()
c
      end
      subroutine tropt_tfgx(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
     $     call errquit('tropt: memory for binv', num, MA_ERR)
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
      if (.not. ma_pop_stack(l_bi)) call errquit('tropt: ma?',99,
     &       MA_ERR)
c
      end
      subroutine tropt_hss_init(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
      integer l_hess, k_hess, l_zmat, k_zmat, l_izmat, k_izmat, i, j
      integer l_c, k_c, l_t, k_t, iat
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
     $     ('tropt_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      old_hessian=.false.
      if(inhess.ne.1) call tropt_check_hess(nvar, old_hessian)
      if (oprint) write(6,*)
      if (old_hessian) then
        if (inhess.eq.2) then
           call tropt_hess_cart_guess()
           if (oprint) write(6,*)
     $        ' Using Cartesian Hessian from previous frequency',
     $        ' calculation'
        else 
          if (oprint) write(6,*) 'Using old Hessian from',
     $                           ' previous optimization'
        end if
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
         call dfill(nvar, 0.5d0, dbl_mb(k_hess), nvar+1)
         if (isystype .ne. 0) then
            if (.not. geom_amatrix_get(geom, amat))
     $           call errquit('geom_frac_to_cart: a', 0, GEOM_ERR)
            do iat = 1, nat_real
               do i = 1, 3
                  dbl_mb(k_hess + (iat-1)*3 + (i-1) +
     $                 ((iat-1)*3 + (i-1))*nat*3) = 0.5*amat(i,i)**2
               end do
            end do
            write(6,*) ' The initial hessian '
            call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
         end if
      else
c
c     Internals.
c
         if (.not. ma_push_get(mt_dbl, nzvar, 'zmat',
     $        l_zmat, k_zmat)) call errquit
     $        ('tropt_init_hess: failed allocating zmat ', nzvar,
     &       MA_ERR)
         if (.not. ma_push_get(mt_int, nzmat, 'izmat',
     $        l_izmat, k_izmat)) call errquit
     $        ('tropt_init_hess: failed allocating izmat ', nzmat,
     &       MA_ERR)
         if (.not. geom_compute_zmatrix(geom, dbl_mb(k_zmat)))
     $        call errquit('tropt_init_hess: zmatrix?',0, GEOM_ERR)
         if (.not. geom_zmt_get_izmat(geom,int_mb(k_izmat),nzmat))
     $        call errquit('tropt_init_hess: izmat?',0, GEOM_ERR)
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
     $        c, zan)) call errquit('tropt_init_hess: geom?',geom,
     &       GEOM_ERR)
         do i = 1, nat
            nuc(i) = nint(zan(i)) ! Won't work if an ECP is loaded
         enddo
         nvar = nzvar
c
c     The angle scaling factors are responsible for a factor of two
c     speedup in most organics.  Clearly we need a better force field.
c
c     The chosen values have been tested on a large no. of 3-21g* SCF
c     optimizations, and a smaller no. of 6-31g* SCF calculations.
c
         if (.not. rtdb_get(rtdb,'tropt:bscale',mt_dbl,1,bscale))
     $        bscale = 1d0
         if (.not. rtdb_get(rtdb,'tropt:ascale',mt_dbl,1,ascale))
     $        ascale = 0.25d0
         if (.not. rtdb_get(rtdb,'tropt:tscale',mt_dbl,1,tscale))
     $        tscale = 0.1d0
         if (oprint) then
            write(6,789) bscale, ascale, tscale
 789        format(' Scaling for Hessian diagonals: bonds =',f5.2,
     $           '  angles =',f5.2,'  torsions =',f5.2)
         endif
         call tropt_hss_init_zmt(dbl_mb(k_hess),
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
      call geom_hnd_put_data('tropt.hess', dbl_mb(k_hess), nvar*nvar)
c
c     Apply constants, constraints and overall scaling
c
 999  if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_c, k_c)) call errquit
     $     ('tropt_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_t, k_t)) call errquit
     $     ('tropt_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      call geom_hnd_get_data('tropt.hess', dbl_mb(k_hess), nvar*nvar)
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
         do iat = 1, nat
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
      if (.not. rtdb_get(rtdb,'tropt:hscale',mt_dbl,1,hscale))
     $     hscale = 1d0
      call dscal(nvar*nvar, hscale, dbl_mb(k_hess), 1)
      if (oprint .and. hscale.ne.1d0) write(6,78) hscale
 78   format(' Scaling initial hessian by ',f6.2)
c
      call geom_hnd_put_data('tropt.hess',dbl_mb(k_hess), nvar*nvar)
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('tropt_init_hess ma corrupt',0, MA_ERR)
c
      END
      SUBROUTINE TROPT_HSS_INIT_ZMT(HSSQ,NINT,ZMAT,IZMAT,SOME,
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
      LOGICAL BADZ
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
      DIMENSION IROW(86),DPARM(6,6),APARM(6,6),RCOV(86)
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
      subroutine tropt_check_hess(nvar, old_hessian)
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
      double precision big,x
c
      big = 1d6
c
c     Look at an existing hessian file and verify it
c
      call util_file_name('tropt.hess',.false.,.false.,filename)
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
CJMC  *** Correct for hessian reading from freq ***
 
 20   call util_file_name('hess',.false.,.false.,filename)
c
      x = big
      if (ga_nodeid() .eq. 0) then
         open(32,file=filename,form='formatted',status='unknown',
     $        err=30,access='sequential')
         read(32,500,err=31,end=31) x
         if (x.eq.big) goto 31
         close(32)
         old_hessian = .true.
         goto 40
c
 31      close(32)
 30      old_hessian = .false.
      endif
CJMC  *** Correct for hessian reading from freq ***
c
 40   call ga_brdcst(323, old_hessian, ma_sizeof(mt_log,1,mt_byte), 0)
c
      return
 500  format(f30.15)
      end
      subroutine tropt_del_hess()
      implicit none
#include "util.fh"
c
c     Delete the Hessian information restart file.
c     
      character*255 opt_hess_fil
c
      call util_file_name('tropt.hess',
     1     .false.,.false.,opt_hess_fil)
      call util_file_unlink(opt_hess_fil)
c
      if (util_print('information',print_low)) then
         write(6,*)
         write(6,*) ' Deleted TROPT restart files '
         write(6,*)
      endif
c
      end
      subroutine tropt_initialize(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb
      integer geom              ! [output] 
c     
c     This routine initializes the common /ctropt/ and
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
      logical tropt_geom_cart_coords_get
c     
      call util_print_push
      call util_print_rtdb_load(rtdb, 'tropt')
      call ecce_print_module_entry('tropt')
      oprint = util_print('information', print_low)
     $     .and. (ga_nodeid() .eq. 0)
      odebug = util_print('debug', print_debug)
     $     .and. (ga_nodeid() .eq. 0)
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
      if (rtdb_get(rtdb,'tropt:redoautoz',mt_log,1,ignore)) then
         if (ga_nodeid() .eq. 0) then
            ignore = rtdb_parallel(.false.)
            ignore = rtdb_delete(rtdb,'tropt:redoautoz')
            call tropt_del_hess() ! Clear hessian info
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
               if (.not. geom_zmtmak(geom,oprint))
     $              call errquit('tropt: regeneration of autoz failed',
     $              0, GEOM_ERR)
               if (oprint) write(6,*)
               if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &              call errquit('hnd_opt: failed storing geometry',0,
     &       RTDB_ERR)
            endif
            if (.not. geom_destroy(geom))
     $           call errquit('tropt: geom_destroy?',0, GEOM_ERR)
            ignore = rtdb_parallel(.true.)
         endif
         call ga_sync()
      endif
c     
c     ----- parameters for optimization tropt -----
c     
      if (.not. rtdb_get(rtdb,'tropt:modsad',mt_int,1,modsad)) 
     $     modsad=0
cjmc
      if (modsad .eq. 0) then
         trust = 0.3d0
      else
         trust = 0.1d0
      endif
      if (.not. rtdb_get(rtdb,'tropt:maxdr',mt_dbl,1,maxdr))
     $     maxdr = 0.3d0
      trust = min(maxdr,trust)
cjmc
      if (.not. rtdb_cget(rtdb,'tropt:xyz',1,xyz))
     $     xyz = ' '
      if (.not. rtdb_get(rtdb,'tropt:eprec',mt_dbl,1,eprec)) then
         if (.not. rtdb_cget(rtdb,'task:theory',1,theory))
     $        theory = ' '
         if (inp_compare(.false.,theory,'dft')) then
            eprec = 5e-6
         else
            eprec = 1e-7
         endif
      endif
      if (.not. rtdb_get(rtdb,'tropt:gmax_tol',mt_dbl,1,gmax_tol))
     $     gmax_tol = 0.00045d0
      if (.not. rtdb_get(rtdb,'tropt:grms_tol',mt_dbl,1,grms_tol))
     $     grms_tol = 0.0003d0
      if (.not. rtdb_get(rtdb,'tropt:xmax_tol',mt_dbl,1,xmax_tol))
     $     xmax_tol = 0.0018d0
      if (.not. rtdb_get(rtdb,'tropt:xrms_tol',mt_dbl,1,xrms_tol))
     $     xrms_tol = 0.0012d0
      if (.not. rtdb_get(rtdb,'tropt:nptopt',mt_int,1,nptopt))
     $     nptopt=20
      if (.not. rtdb_get(rtdb,'tropt:inhess',mt_int,1,inhess))
     $     inhess=0
      if (.not. rtdb_get(rtdb,'tropt:linopt',mt_int,1,linopt)) 
     $     linopt=1
      if (.not. rtdb_get(rtdb,'tropt:moddir',mt_int,1,moddir)) 
     $     moddir=0
      if (.not. rtdb_get(rtdb,'tropt:modsad',mt_int,1,modsad)) 
     $     modsad=0
      if (.not. rtdb_get(rtdb,'tropt:sadstp',mt_dbl,1,sadstp)) 
     $     sadstp=0.1d0
      if (.not. rtdb_get(rtdb,'tropt:oqstep',mt_log,1,oqstep))
     $     oqstep = .true.
      if (.not. rtdb_get(rtdb,'tropt:modupd',mt_int,1,modupd)) then
         if (modsad .eq. 0) then
            modupd = 1          ! BFGS update for minimization
         else
            modupd = 2          ! PSB update for saddle point
         endif
      endif
      if (.not. rtdb_get(rtdb,'tropt:ocheckgrad',mt_log,1,ocheckgrad))
     $     ocheckgrad = .false.

      if (.not. rtdb_get(rtdb,'includestress',mt_log,1,ostress)) then
         ostress = .false.
      end if
      if (.not. rtdb_get(rtdb,'includelattice',mt_log,1,ostress2)) then
         ostress2 = .false.
      end if
      if ((.not.ostress).and.(ostress2)) ostress2 = .false.
      if ((ostress)     .and.(ostress2)) ostress  = .false.
      
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
         if (.not. geom_rtdb_store(rtdb, geom, 'troptinitial'))
     &        call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
         if (.not. geom_destroy(geom))
     $        call errquit('tropt: geom_destroy?',0, GEOM_ERR)
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
     $     call errquit('tropt: systype?',0, GEOM_ERR)
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
         write(6,1) gmax_tol, grms_tol, xmax_tol, xrms_tol, maxdr,
     $        eprec,
     $        nptopt, inhess, modupd, modsad, moddir, 
     $        source
 1       format(
     $        ' maximum gradient threshold         (gmax) = ', f10.6,/,
     $        ' rms gradient threshold             (grms) = ', f10.6,/,
     $        ' maximum cartesian step threshold   (xmax) = ', f10.6,/,
     $        ' rms cartesian step threshold       (xrms) = ', f10.6,/,
     $        ' maximum step size                 (maxdr) = ', f10.6,/,
     $        ' energy precision                  (eprec) = ', 1p,d9.1,
     $        0p,/,
     $        ' maximum number of steps          (nptopt) = ', i4,/,
     $        ' initial hessian option           (inhess) = ', i4,/,
     $        ' hessian update option            (modupd) = ', i4,/,
     $        ' saddle point option              (modsad) = ', i4,/
     $        ' initial eigen-mode to follow     (moddir) = ', i4,/,
     $        ' source of zmatrix                (source) = ', a8,/)
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
               call errquit('tropt: cartesian constraints+internals',0,
     &       GEOM_ERR)
            endif
            call ga_sync()
         endif
      else
         call tropt_cart_pmat(rtdb, geom)
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
     $        call errquit('tropt: zmatrix?',0, GEOM_ERR)
      else
c     FRACTIONAL?
         if (.not. tropt_geom_cart_coords_get(geom, sp))
     $        call errquit('tropt: geom?',0, GEOM_ERR)
      endif
c     
      end
      subroutine tropt_hessian_update()
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
      double precision hds(max_nvar)
      double precision dsds, dshds, dsdg
      integer l_hess, k_hess, i, j
      integer ind
      double precision ddot
      external ddot
      ind(i,j) = k_hess + i + (j-1)*nvar - 1
c
      if (alpha .eq. 0d0) call errquit
     $     ('tropt_hessian_update: zero step?',0, GEOM_ERR)
      call dscal(nvar, alpha, ds, 1)
c     
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('tropt_hessian_update: memory for hessian',nvar**2,
     &       GEOM_ERR)
      call geom_hnd_get_data('tropt.hess',dbl_mb(k_hess), nvar**2)
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
      if(modupd.le.1) then       
c     
c     ----- -bfgs- update -----
c     
         if(dsdg.gt.1d-14) then
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
         if (abs(dsdg).gt.1d-8) then
            do i=1,nvar
               do j=1,nvar
                  dbl_mb(ind(i,j))=dbl_mb(ind(i,j))
     $                 + ((g(i)-gp(i))-hds(i))*ds(j)/dsds
     1                 + ((g(j)-gp(j))-hds(j))*ds(i)/dsds
     2                 - ds(i)*ds(j)*(dsdg-dshds)/(dsds*dsds)
               enddo    
            enddo    
         endif
      endif    
c
      call geom_hnd_put_data('tropt.hess',dbl_mb(k_hess), nvar**2)
      if (.not. ma_pop_stack(l_hess)) call errquit
     $     ('tropt_hessian_update: ma?',0, MA_ERR)
c
      call dscal(nvar, 1d0/alpha, ds, 1)
c     
      end
c
      double precision function tropt_energy_step(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
      logical tropt_geom_cart_coords_get
      logical tropt_geom_cart_coords_set
c
c     Get original coordinates so can restore them
c
c     FRACTIONAL?
      if (.not. tropt_geom_cart_coords_get(geom, xold))
     $     call errquit('tropt_energy_step: coordinates?',geom,
     &       GEOM_ERR)
c
c     Take the step
c
      call tropt_take_step(rtdb, geom)
c
c     Energy at displaced point
c
      if (.not. task_energy(rtdb)) call errquit
     $     ('tropt_energy_step: energy failed',0, GEOM_ERR)
      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, aenergy))
     $     call errquit('tropt_energy_step: rtdb?',0, RTDB_ERR)
c      
c     Restore the starting geometry
c
c     FRACTIONAL?
      if (.not. tropt_geom_cart_coords_set(geom, xold))
     $     call errquit('tropt_energy_step: coordinates?',geom,
     &       GEOM_ERR)
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $     call errquit('tropt_energy_step: grs?',geom, RTDB_ERR)
c
      tropt_energy_step = aenergy
c
      end
      subroutine tropt_take_step(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
      logical tropt_geom_cart_coords_get
      logical tropt_geom_cart_coords_set
c
c     Get original coordinates
c
c     FRACTIONAL?
      if (.not. tropt_geom_cart_coords_get(geom, xold))
     $     call errquit('tropt_energy_step: coordinates?',geom,
     &       GEOM_ERR)
c
c     Take the step
c
      if (zcoord) then
c
c     Compute first order step only so we can compute xmax/xrms
c
         if (.not. ma_push_get(mt_dbl, ncart*nzvar,'binv',l_bi, k_bi))
     $      call errquit('tropt_take: memory for binv', ncart*nzvar,
     &      MA_ERR)
         call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nzvar)
         call dgemv('n', ncart, nzvar, alpha, dbl_mb(k_bi), ncart,
     $        ds, 1, 0.0d0, xnew, 1)
         if (.not. ma_pop_stack(l_bi)) call errquit('tropt: ma?',99,
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
*         if (.not. tropt_geom_cart_coords_set(geom, xnew))
*     $        call errquit('tropt_energy_step: coordinates?',geom,
*     &       GEOM_ERR)
      else
         do i = 1, ncart
            xnew(i) = alpha*ds(i)
         enddo
         call sym_grad_symmetrize(geom, xnew)
         do i = 1, ncart
            xnew(i) = xold(i) + xnew(i)
         enddo
c     FRACTIONAL?
         if (.not. tropt_geom_cart_coords_set(geom, xnew))
     $        call errquit('tropt_energy_step: coordinates?',geom,
     &       GEOM_ERR)
      endif
c
c     Must ensure the geometry has the required symmetry even after
c     enforcing it on the step.  Should use
c     an error criterion consistent with the step size.
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
     $     call errquit('tropt_energy_step: grs?',geom, RTDB_ERR)
c
      end
      subroutine tropt_print(geom, istep, closefile)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"

      integer geom, istep
      logical closefile
c
c     Print out stuff
c
      integer i
      double precision diff, de
      logical geom_print_zmatrix, oprint_grad
      character*9 cvg1, cvg2, cvg3, cvg4
      character*1 mark
      character*255 filename, dir
      logical mol_geom_print_xyz
      external mol_geom_print_xyz
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
c
      if (oprint) then 
         mark = '@'
         if (istep .gt. 1) mark = ' '
         write(6,1) mark, mark
         mark = '@'
         write(6,2) mark, istep-1, energy, de,
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
     $        write(6,*) ' Internal coordinate gradient before',
     $        ' constraints'
         if (.not. geom_print_zmatrix(geom, gq, 'Gradient',.true.))
     $        call errquit('tfgx: print internal gradient?',geom,
     &       GEOM_ERR)
         if (diff .gt. 1d-6) then
            write(6,*) ' Internal coordinate gradient after',
     $           ' constraints'
            if (.not. geom_print_zmatrix(geom, g, 'Gradient',.true.))
     $           call errquit('tfgx: print internal gradient?',geom,
     &       GEOM_ERR)
         endif
      endif
c
c     Print a trajectory file
c 
       if (ga_nodeid().eq.0 .and. xyz.ne.' ') then
         dir      = ' '
         filename = ' ' 
         call util_directory_name(dir, .false., 0)
         write(filename,13) dir(1:inp_strlen(dir)), 
     $        xyz(1:inp_strlen(xyz))
 13      format(a,'/',a,'.xyz')
         if (istep.eq.1) then
           open(88,file=filename,form='formatted',status='unknown',
     $          access='sequential',err=133)
           rewind(88)
         end if
c
         if (.not. mol_geom_print_xyz(geom, 88, energy))
     $       call errquit('tropt: mol_geom_print_xyz?',0, GEOM_ERR)
         call util_flush(88)
c
         if (closefile) close(88,status='keep',err=133)
      end if
c
      return
 133  call errquit('tropt: error opening/closing xyz file',0, GEOM_ERR)
c
      end
      logical function tropt_converged(istep)
      implicit none
#include "nwc_const.fh"
#include "ctropt.fh"
#include "global.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
      integer istep
      double precision de
c
c     Return true if we have converged
c
c     Nothing is modified.  Assumes tropt_compute_info()
c     has been called.
c     gmax_tol,            ! [user] tolerance for max internal gradient
c     grms_tol,            ! [user] tolerance for rms internal gradient
c     xrms_tol,            ! [user] tolerance for rms cartesian step
c     xmax_tol,            ! [user] tolerance for max cartesian step
c
      de = abs(energy-energyp)
      tropt_converged = 
     $     ((gmax .lt. gmax_tol) .and. (grms .lt. grms_tol) .and.
     $     (xrms .lt. xrms_tol)  .and. (xmax .lt. xmax_tol))
     $     .or. 
     $     ((gmax.lt.0.01d0*gmax_tol) .and. (grms.lt.0.01d0*grms_tol))
CJMC
c 
c    *** Convergence over energy criterion
c
     $     .or.
     $      ((gmax .lt. gmax_tol) .and. (grms .lt. grms_tol).and.
     $      (de .lt. eprec))
     $     .or.
c 
c    *** Initial structure already a critical point
c
     $      ((gmax .lt. gmax_tol) .and. (grms .lt. grms_tol).and.
     $      (istep .eq. 1))
CJMC
      if ((gmax .lt. gmax_tol) .and. (grms .lt. grms_tol).and.
     $      (de .lt. eprec)) then
        if (ga_nodeid().eq.0) then
          write(luout,*) "Convergence reached over gradient and energy"
          write(luout,*) "Energy change = " , de
          write(luout,*) "Energy precision = " , eprec
        endif
      end if

      call ga_brdcst(1001,tropt_converged,ma_sizeof(mt_log,1,mt_byte),0)
c
      end
      subroutine tropt_cart_pmat(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
      logical task_qmmm
      logical tropt_geom_cart_coords_get
      double precision ddot
      external ddot
      ind(i,j) = k_pmat + i-1 + (j-1)*ncart
c     
c     FRACTIONAL?
      if (.not. tropt_geom_cart_coords_get(geom, coords))
     $     call errquit('tropt_cart_pmat: geom?',geom, GEOM_ERR)
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
     $     ('tropt_cart_pmat: memory for pmat',ncart**2, GEOM_ERR)
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

      if ( rtdb_get_info(rtdb, 'geometry:actlist', ma_type, 
     $     nelem, date) .or.
     $     rtdb_get_info(rtdb, 'geometry:inactlist', ma_type, 
     $     nelem, date) .or.
     $     isystype .ne. 0 .or.
     $     task_qmmm .or. 
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
     $     ('tropt_cart_bmat: ma?',0, MA_ERR)
c     
      end
      subroutine tropt_project_hess_grad(hess, pg)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
     $     ('tropt_proj_h_g: memory for pmat',nvar**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'pmat',
     $     l_pmat, k_pmat)) call errquit
     $     ('tropt_proj_h_g: memory for work',nvar**2, MA_ERR)
c
      call geom_hnd_get_data('tropt.hess',hess, nvar**2)
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
     $     ('tropt_p_h_g:ma?',0, MA_ERR)
c
      end
      subroutine tropt_compute_info()
      implicit none
#include "nwc_const.fh"
#include "ctropt.fh"
c
c     Compute stuff used for printing and convergence tests
c
c     gmax = maxmimum gradient element in optimization variables
c     grms = rms grad
c     smax = maximum step in opt. var
c     srms = rms step
c     
c     xrms and xmax are computed by tropt_take_step from the
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
      subroutine tropt_hess_cart_guess()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
     $     ('tropt_hess_cart_guess: ma?', ncart*nvar, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart',
     $     l_hc, k_hc)) call errquit
     $     ('tropt_hess_cart_guess: ma?', ncart**2, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart2',
     $     l_hq, k_hq)) call errquit
     $     ('tropt_hess_cart_guess: ma?', nvar**2, MA_ERR)
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
      call geom_hnd_put_data('tropt.hess',dbl_mb(k_hc), nvar**2)
c
      if (.not. ma_chop_stack(l_bi))
     $     call errquit('tropt_hess_cart_guess: ma corrupt?',0, MA_ERR)
c
      return
10000 format(f30.15)
99990 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as unknown file'
      call errquit('tropt_hess_cart: fatal error', 911, GEOM_ERR)
99991 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as new file'
      call errquit('tropt_hess_cart: fatal error', 911, GEOM_ERR)
99992 write(6,*)' error in reading <',fname(1:inp_strlen(fname)),
     $     '> as hessian file'
      call errquit('tropt_hess_cart: fatal error', 911, GEOM_ERR)
      end
      subroutine tropt_symmetrize_step(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "mafdecls.fh"
      integer geom
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
      if (ostress) return       ! Yikes!
      if (ostress2) return       ! Yikes!
c
      if (.not. ma_push_get(mt_dbl, ncart*nvar,'binv',l_bi, k_bi))
     $     call errquit('tropt_sym: memory for binv', ncart*nvar,
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
      if (.not. ma_pop_stack(l_bi)) call errquit('tropt_sym:ma',0,
     &       MA_ERR)
c
      end
      subroutine tropt_compute_actual_step(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
      logical tropt_geom_cart_coords_get
c
      if (odebug) then
         write(6,*) ' Expected ds '
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
      call dcopy(nvar, sp, 1, ds, 1) ! Old coordinates into ds()
      if (zcoord) then
         if (.not. geom_compute_zmatrix(geom, sp)) ! New coordinates into sp()
     $        call errquit('tropt: zmatrix?',0, GEOM_ERR)
      else
c     FRACTIONAL?
         if (.not. tropt_geom_cart_coords_get(geom, sp))
     $        call errquit('tropt: geom?',0, GEOM_ERR)
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
      logical function tropt_geom_cart_coords_get(geom, coords)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
      integer geom
      double precision coords(*)
c
c     If we are doing a periodic system and not using internals
c     then we want the fractional coordinates.  Otherwise cartesian.
c
c     If we are including stress append the amatrix
c
      if (.not. geom_cart_coords_get(geom, coords))
     $     call errquit('tropt: geom cart?',0, GEOM_ERR)
c
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_cart_to_frac(geom, coords))
     $           call errquit('tropt: frac_to_cart?',0, GEOM_ERR)
      endif
c
      if (ostress) then
         if (.not. geom_amatrix_get(geom, coords(3*nat_real+1)))
     $        call errquit('tropt: failed to get amatrix',0,0)
      endif
      if (ostress2) then
         if (.not. geom_lattice_get(geom, coords(3*nat_real+1)))
     $        call errquit('tropt: failed to get lattice',0,0)
      endif
c
      tropt_geom_cart_coords_get = .true.
c
      end
      logical function tropt_geom_cart_coords_set(geom, coords)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
     $        call errquit('tropt: failed to set amatrix',0,0)
      endif
      if (ostress2) then
         if (.not. geom_lattice_set(geom, coords(3*nat_real+1)))
     $        call errquit('tropt: failed to set lattice',0,0)
      endif

      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_frac_to_cart(geom, coords))
     $           call errquit('tropt: frac_to_cart?',0,0)
      endif
      if (.not. geom_cart_coords_set(geom, coords))
     $     call errquit('tropt: geom cart?',0,0)
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_cart_to_frac(geom, coords))
     $        call errquit('tropt: frac_to_cart?',0,0)
      endif
c
      tropt_geom_cart_coords_set = .true.
c
      end
      subroutine tropt_get_grad(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
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
      logical geom_grad_cart_to_frac
c
      if (.not. rtdb_get(rtdb, 'task:gradient', mt_dbl, ncart,
     $     gx)) call errquit('tropt: could not get gradient',0,0)
      if (isystype .ne. 0) then
         if (.not. geom_grad_cart_to_frac(geom, gx))
     $        call errquit('tropt: frac_to_cart?',0,0)
      end if
      if (ostress) then
         if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $   call errquit('tropt: stress theory not specified',0,RTDB_ERR)
         if (theory.eq.'pspw') then
          if (.not. rtdb_get(rtdb, 'pspw:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('tropt: could not get stress',0,0)
         else if (theory.eq.'band') then
          if (.not. rtdb_get(rtdb, 'band:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('tropt: could not get stress',0,0)
         else if (theory.eq.'paw') then
          if (.not. rtdb_get(rtdb, 'paw:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('tropt: could not get stress',0,0)
         else
           call errquit('tropt: no stress in theory',0,RTDB_ERR)
         end if
      endif

      if (ostress2) then
         if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $   call errquit('tropt: stress theory not specified',0,RTDB_ERR)
         if (theory.eq.'pspw') then
          if (.not. rtdb_get(rtdb, 'pspw:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('tropt: could not get stress',0,0)
         else if (theory.eq.'band') then
          if (.not. rtdb_get(rtdb, 'band:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('tropt: could not get stress',0,0)
         else if (theory.eq.'paw') then
          if (.not. rtdb_get(rtdb, 'paw:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('tropt: could not get stress',0,0)
         else
           call errquit('tropt: no stress in theory',0,RTDB_ERR)
         end if
      endif

c
      end
      subroutine tropt_check_grad(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb,geom
c
      integer i,ia,k
      double precision e0, ep, em, grad, hess, step
      double precision step0
c
      logical task_gradient, task_energy, tropt_geom_cart_coords_get,
     $     tropt_geom_cart_coords_set
c
      if (ga_nodeid() .eq. 0) then
         write(6,*) 
         write(6,*) ' Checking the gradient '
         write(6,*)
         call util_flush(6)
      end if
c
      if (.not. geom_rtdb_store(rtdb,geom,'checkgrad')) 
     $     call errquit('tropt: failed saving geom for checkgrad',0,0)
c
c     Check if fract_cart is OK
c
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. tropt_geom_cart_coords_get(geom,sp))
     $      call errquit('tropt: failed get coords for checkgrad',0,0)
         write(6,*) ' These are the fractional coords'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
         if (.not. geom_frac_to_cart(geom, sp))
     $       call errquit('tropt: frac_to_cart?',0,0)
         write(6,*) ' These are the cartesian coords'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
         if (.not. geom_cart_to_frac(geom, sp))
     $           call errquit('tropt: frac_to_cart?',0,0)
         write(6,*) ' These are the fractional coords regenerated'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
c
         if (.not. tropt_geom_cart_coords_set(geom,sp))
     $      call errquit('tropt: failed set coords for checkgrad',0,0)
         if (.not. tropt_geom_cart_coords_get(geom,sp))
     $      call errquit('tropt: failed get coords for checkgrad',0,0)
         write(6,*) ' These are the fractional coords after set'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
c
         if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &      call errquit('hnd_opt: failed storing geometry',0,0)
         if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &      call errquit('hnd_opt: no geometry ', 911,0)
         if (.not. tropt_geom_cart_coords_get(geom,sp))
     $      call errquit('tropt: failed get coords for checkgrad',0,0)
         write(6,*) ' These are the fractional coords after st/lo'
         call output(sp,1,nat_real*3,1,1,nat_real*3,1,1)
      end if
c
      if (.not. task_gradient(rtdb))
     $     call errquit('tropt: task_gradient failed',0,0)
      call tropt_get_grad(rtdb,geom) ! Into gx
      if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, e0))
     $     call errquit('tropt: could not get energy',0,0)

      if (.not.rtdb_get(rtdb,"tropt:checkgradstep",mt_dbl,1,step0)) 
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
     $        call errquit('tropt: failed load geom for checkgrad',0,0)
         if (.not. tropt_geom_cart_coords_get(geom,sp))
     $    call errquit('tropt: failed get coords for checkgrad',0,0)
         sp(i) = sp(i) + step
         if (.not. tropt_geom_cart_coords_set(geom,sp))
     $    call errquit('tropt: failed set coords for checkgrad',0,0)
         if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &    call errquit('hnd_opt: failed storing geometry',0,0)
         if (.not. task_energy(rtdb)) call errquit
     $        ('tropt_energy_step: energy failed',0,0)
         if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, ep))
     $     call errquit('tropt: could not get energy',0,0)
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
     $        call errquit('tropt: failed load geom for checkgrad',0,0)
         if (.not. tropt_geom_cart_coords_get(geom,sp))
     $     call errquit('tropt: failed get coords for checkgrad',0,0)
         sp(i) = sp(i) - step
         if (.not. tropt_geom_cart_coords_set(geom,sp))
     $     call errquit('tropt: failed set coords for checkgrad',0,0)
         if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &        call errquit('hnd_opt: failed storing geometry',0,0)
         if (.not. task_energy(rtdb)) call errquit
     $        ('tropt_energy_step: energy failed',0,0)
         if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, em))
     $        call errquit('tropt: could not get energy',0,0)
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
     $     call errquit('tropt: failed load geom for checkgrad',0,0)
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &     call errquit('hnd_opt: failed storing geometry',0,0)
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tropt_pickstp(rtdb,geom,istep)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
      integer rtdb
      integer geom
      integer istep
c
c     this routine for minimization
c
c     put into ds() a search direction in the optimization
c     variables (internal or cartesian) based upon the
c     current gradient, g(), and hessian.  apply constraints.
c
c     apply step restrictions by recommending an initial
c     value for the line search parameter alpha.
c
c     only alpha and ds() are modified.
c
      integer i, iat
cjmc
      integer iact, bmin, lowest, nneg
      double precision lambda, tau, term
      double precision tiniest
      double precision big
      parameter (tiniest = 1.0e-14) 
      parameter (big = 1d6)
      logical restricted
      double precision tmp1, tmp2, ovtol
      integer idamax, mod2
cjmc
*     integer info
      integer l_hess, k_hess, l_work, k_work, lenwork
      double precision eigval(max_nvar) ! hessian eigenvalues
      double precision pg(max_nvar) ! p.g
      double precision gv(max_nvar) ! gradient along eigenvectors
      double precision dv(max_nvar) ! step along eigenvectors
      double precision dsmax    ! max. value of current step (smax is prev.)
c
      double precision beta, s0g0, s0g1, s1g0, s1g1, numerator,
     $     denominator
      double precision bohr, deg ! for printing purposes
      double precision trustds  ! restriction of step in opt. variable
      logical geom_print_zmatrix
      logical ophigh
      double precision ddot
      external ddot
c
c     get the hessian and gradient with appropriate projectors
c     applied following peng, ayala, schlegel and frisch so that
c     redundant internal modes are shifted to high eigenvalues.
c
      ophigh = util_print('high', print_high)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('tropt_pickstp: memory for hessian',nvar**2, ma_err)
      call tropt_project_hess_grad(dbl_mb(k_hess), pg)
c
c     diagonalize the hessian.  should really do the generalized
c     eigenvalue problem since the underlying basis is not independent
c     (if we are using autoz). not yet being done.  
c
c     to cause degenerate eigenvalues to be resolved into symmetry
c     adapted combinations use jacobi not dsyev and screen out junk
c
      lenwork = max(nvar**2,100)
      if (.not. ma_push_get(mt_dbl, lenwork, 'work',
     $     l_work, k_work)) call errquit
     $     ('tropt_pickstp: memory for hessian', lenwork, ma_err)
      do i = 0, nvar**2-1
         if (abs(dbl_mb(k_hess+i)).lt.1d-8) dbl_mb(k_hess+i) = 0d0
      enddo
c
c     have eigenvalues in eigval, eigenvectors in dbl_mb(k_hess).
c
      call util_jacobi(nvar, dbl_mb(k_hess), nvar, eigval)
      if (odebug .or. (util_print('hvecs',print_never) 
     $     .and. ga_nodeid().eq.0)) then
         write(6,*) ' Eigenvalues of the hessian '
         call doutput(eigval, 1, nvar, 1, 1, nvar, 1, 1)
         write(6,*) ' Eigenvectors of the hessian '
         call output(dbl_mb(k_hess), 1, nvar, 1, nvar, nvar, nvar, 1)
      endif
c
c     project the gradient onto the hessian eigenvectors
c
      call dgemv('t', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $     pg, 1, 0d0, gv, 1)
      if (odebug) then
         write(6,*) ' Gradient in spectral form '
         call doutput(gv, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     *** count negative eigenvalues of hessian matrix ***
c
      nneg = 0
      do iact=1,nvar
        if (eigval(iact).lt.0.0) then
          if (abs(eigval(iact)).gt.1.0d-8) nneg = nneg + 1
        endif
      end do
c
c     *** eigenvector following for transition state search ***
c
      if (modsad.eq.1) then
c
        ovtol   = 0.7d0
        moddir  = abs(moddir)
        hfollow = moddir
        if ((istep.eq.1).and.(moddir.ne.0))
     $    call dcopy(nvar,dbl_mb(k_hess+(moddir-1)*nvar),1,saddir,1)
c     
c     Maximum overlap for tracking the mode.  
c     
        if (moddir.ne.0) then
          call dgemv('t', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $         saddir, 1, 0d0, dv, 1)
          do i = 1, nvar         ! This to avoid matching rot/trans/constraint
            if (eigval(i) .gt. 1d2) dv(i) = 0d0
          enddo
          if (util_print('modeoverlap',print_debug) .and. 
     $        ga_nodeid().eq.0) then
            write(6,*) ' Overlap of modes with searchdir'
            call output(dv, 1, nvar, 1, 1, nvar, 1, 1)
          endif
c     
          moddir = idamax(nvar, dv, 1)
          do i = 1, nvar
            if (abs(dv(i)) .gt. abs(dv(moddir))*ovtol) then
              moddir = i
              goto 304
            endif
          enddo
 304      continue
c     
c     Print out first two maximum overlaps to warn of ambiguity
c     
          tmp1 = abs(dv(moddir))
          dv(moddir) = 0d0
          mod2 = idamax(nvar, dv, 1)
          tmp2 = abs(dv(mod2))
          if (oprint) write(6,98) moddir, tmp1, mod2, tmp2
 98       format(' Maximum overlap of modes: mode=', i3, 1p, 
     $          ' overlap=',d8.1, 0p, ' mode=', i3, 1p,' overlap=',d8.1)
          if (tmp1 .eq. 0d0) call errquit('tropt_sad_dir: no mode?',0,
     &        GEOM_ERR)
c
          call dcopy(nvar, dbl_mb(k_hess+(moddir-1)*nvar), 1, saddir, 1)
          hfollow = moddir
c
        else
c
          moddir  = 1
          hfollow = moddir
          call dcopy(nvar, dbl_mb(k_hess), 1, saddir, 1)
c
        endif
c
        if (util_print('sadmode',print_default)
     $     .and.ga_nodeid().eq.0)then
          write(6,*) 
          write(6,*) 'The mode being followed to the saddle point'
          if (zcoord) then
            if (.not.geom_print_zmatrix(geom, saddir,'Sadmode',.true.))
     $          call errquit('tropt: printing sadmode?',0, GEOM_ERR)
          else
            write(6,*)
            write(6,88) (saddir(i),i=1,nvar)
 88         format(1x,8f9.4)
          endif
        endif
c     
c     Save moddir so that it is picked up for a restart
c     
        if (.not. rtdb_put(rtdb,'tropt:moddir',mt_int,1,moddir)) 
     $       call errquit('tropt_pickstp: rtdb?',0, RTDB_ERR)
c     
      end if
!
!     *** find lowest (countable) eigenvalue position ***
!     ***    assuming a non-shifted hessian matrix    ***
!
      do iact=1,nvar
        if (abs(eigval(iact)).gt.1.0d-8) then
          lowest = iact
          goto 100
        end if
      end do
c
  100 continue
c
      if (modsad.eq.1) then
        do iact=1,nvar
          if (iact.ne.hfollow) then
            if (abs(eigval(iact)).gt.1.0d-8) then
              bmin = iact
              lowest = hfollow
              goto 200
            end if
          end if
        end do
      end if
c
  200 continue
c
!
!     *** test for the "hard case", i.e. f1 = 0 ***
!
      if (abs(gv(lowest)).lt.1.e-8) then
!
!     *** set lambda to lowest eigenvalue ***
!
        lambda = -eigval(lowest)
        if (modsad.eq.1) then
          if (eigval(hfollow).gt.-eigval(bmin)) lambda =  eigval(lowest)
          if (eigval(hfollow).lt.-eigval(bmin)) lambda = -eigval(lowest)
        end if
!
!     *** calculate reduced constrained step ***
!
        call dfill(nvar, 0d0, dv, 1)
        do iact=1,nvar
          if (iact.ne.lowest) then
            dv(iact) = -gv(iact)/(eigval(iact) + lambda + tiniest)
          end if
        end do
!
!     *** check step length and modify, if necessary ***
!
        term =  ddot(nvar, dv, 1,  dv, 1)
        if (term.le.trust**2) then
          tau = sqrt(trust**2 - term)
          dv(lowest) = tau
          if (oprint) write(6,*) "hard case encountered"
          slength = sqrt(ddot(nvar, dv, 1,  dv, 1))
          goto 300
        end if
      end if
!
!     *** quasi-newton step type ***
!
      restricted = .false.
      lambda = 0.0
!
!     *** calculate unconstrained quasi-newton step ***
!
      call tropt_raphson(gv,dv,eigval,lambda)
!
!     *** calculate step length ***
!
      slength = sqrt(ddot(nvar, dv, 1,  dv, 1))
!
!     *** check step length and hessian eigenvalue spectrum ***
!
      if (slength.gt.trust) then
        restricted = .true.
        if (oprint) write (6,5000) slength,trust
 5000   format ('newton step ',g10.3,' exceeded trust radius ',g10.3)
      else if ((modsad.eq.1).and.(nneg.ne.1)) then
        restricted = .true.
        if (oprint) write (6,5000) slength,trust
        write  (6,5010)
 5010   format ('the hessian spectrum is not correct for ts search')
      else if ((modsad.ne.1).and.(nneg.ne.0)) then
        restricted = .true.
        if (oprint) write (6,5020)
 5020   format ('the hessian spectrum is not positive definite')
      end if
!
!     *** restrict step size, if necessary ***
!
      if (restricted) then
        call tropt_lambda(gv,eigval,lambda)
        call tropt_raphson(gv,dv,eigval,lambda)
        slength = sqrt(ddot(nvar, dv, 1,  dv, 1))
      end if
!
!     *** print shift factor(s) ***
!
 300  continue

      if (oprint) then
        write (6,5030) lambda
 5030   format ('lambda for step: ',g10.3)
        write (6,5040) slength
 5040   format ('step magnitude : ',g10.3)
      end if

      if (odebug) then
         write(6,*) ' Step in spectral form '
         call doutput(dv, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     transform back to optimization space
c
      call dgemv('n', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $     dv, 1, 0d0, ds, 1)
      if (odebug) then
         write(6,*) ' Step in optimization variables'
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
c     enforce symmetry
c
      call tropt_symmetrize_step(geom)
c
c     enforce frozen atoms in cartesians
c
      if (ga_nodeid().eq.0.and.ophigh) 
     $     write(6,*) 'Zeroing constrained gradient'

      if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
         do iat = 1, nat
            if (.not. oactive(iat)) then
               do i = 1, 3
                  ds((iat-1)*3+i) = 0.0
               end do
            end if
         end do
      end if
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('tropt_pickstp: ma?',0, MA_ERR)
c
c     edo seems to have encountered a case where different processors
c     generated different steps.  to prevent this, broadcast the 
c     critical info to everyone.
c
      call ga_brdcst(1,ds,8*nvar,0)
      call ga_brdcst(2,alpha,8,0)
c
      if (util_print('searchdir',print_high) .and. 
     $     ga_nodeid().eq.0) then
         write(6,*)
         write(6,*) '       the search direction'
         if (zcoord) then
            bohr = 0.52917715d0
            deg  = 0.52917715d0*180d0/(4d0*atan(1d0))
            call geom_zmat_ico_scale(geom, ds, bohr, deg)
            if (.not. geom_print_zmatrix(geom, ds, 'search  ',.true.))
     $         call errquit('tropt:print search_dir?',geom, geom_err)
            call geom_zmat_ico_scale(geom, ds, 1d0/bohr, 1d0/deg)
         else
            call output(ds,1,3,1,nat,3,nat,1)
         endif
         write(6,*)
         call util_flush(6)
      endif
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tropt_radius()
!  update trust radius for restricted step optimization.
      implicit none
#include "nwc_const.fh"
#include "ctropt.fh"

       double precision eps,tscale
       double precision tiniest
       double precision mindr
       double precision rmlow,rmup

       parameter (eps = 1d-4)
       parameter (tscale = 2.0)
       parameter (tiniest = 1.0d-14) 
       parameter (mindr = 0.001)
       parameter (rmlow = 0.25, rmup = 0.75)

       double precision echange,ratio,upper,epredict
!
!    *** assume the step is good enough
!
        redostep = .false.
!
!     *** calculate energy change ***
!
        echange = energyp - energy
!
!     *** quick return for small uphill steps ***
!
        if ((modsad.eq.0).and.
     $     ((echange.lt.0.0).and.(abs(echange).lt.eps))) return
!
!     *** calculate energy ratio ***
!
        call tropt_etaylor(epredict)
        ratio = echange/max(abs(epredict),tiniest)
!
!     *** decrease/increase trust radius ***
!
        if (modsad.eq.1) then
          ratio = abs(ratio)
          upper = 1.0 + rmup
          if ((ratio.lt.rmlow).or.(ratio.gt.upper)) then
            trust = trust/tscale
          else
            if (slength.le.(0.8*trust)) then
              trust = trust
            else
              trust = tscale*trust
            end if
          end if
        else 
!
!    *** the step sucks keep trying
!
          if (ratio.lt.0.0) redostep = .true.
          if (ratio.gt.rmup) then
            if (slength.le.(0.8*trust)) then
              trust = trust
            else
              trust = tscale*trust
            end if 
          else if ((ratio.ge.rmlow).and.(ratio.le.rmup)) then
            trust = trust    
          else
            trust = trust/tscale    
          end if
        end if
!
!     *** avoid infinite rejecttion ***
!
        if ((trust.lt.mindr).and.(redostep)) then
          redostep = .false.
          if (oprint) then
            write(6,*) "Can not reduce the step anymore"
            write(6,*) "Forcing to accept the step, hope it improves"
          end if
        end if
!
!     *** bound trust radius ***
!
        trust = max(trust,mindr)
        trust = min(trust,maxdr)

        if (oprint) then
          ratio = max(-10.d0,min(ratio,10.d0))
          write (6,5000) -echange,ratio,trust
 5000     format (/,t2,'energy change = ',f10.6,   
     $        /,t2,'energy ratio  = ',f10.6,  
     $        /,t2,'trust radius  = ',f10.6)
        end if

      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tropt_etaylor(epredict)
!  calculate predicted change by the quadratic taylor series expansion.

      implicit none

#include "errquit.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "util.fh"
#include "mafdecls.fh"

      integer l_hess, k_hess
      double precision epredict
      double precision gamma(nvar)
      double precision ddot
      external ddot
!
!     *** calculate predicted energy change ***
!
      call dfill(nvar, 0d0, gamma, 1)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $   l_hess, k_hess)) call errquit
     $   ('tropt_etaylor: memory for hessian',nvar**2,
     $      geom_err)
      call geom_hnd_get_data('tropt.hess',dbl_mb(k_hess), nvar**2)
c
c     form bits and pieces that are needed
c     
      call dgemv('n',nvar,nvar,1d0,dbl_mb(k_hess),nvar,
     $     ds,1,0d0,gamma,1)
c
      epredict = -ddot(nvar, gp, 1, ds, 1)
     $           - 0.5*ddot(nvar, ds, 1, gamma,1 )
c
      if (.not. ma_pop_stack(l_hess)) call errquit
     $     ('tropt_etaylor: ma?',0, MA_ERR)
c
      if (oprint) write (6,5000) epredict
5000    format ('predicted energy change: ',f20.9)
c
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tropt_raphson(de1,dr,eigval,lambda)
! do newton-raphson step with hessian eigenvalue shift.

      implicit none
#include "nwc_const.fh"
#include "ctropt.fh"

      double precision eps
      parameter (eps = 1d-5)
      double precision tiniest
      parameter (tiniest = 1.0d-14) 

      double precision de1(nvar),dr(nvar),eigval(nvar)

      integer iact
      double precision lambda,term
!
!     *** calculate displacement vector ***
!
      call dfill(nvar, 0d0, dr, 1)

      do iact=1,nvar
        if ((modsad.eq.1).and.(iact.eq.hfollow)) then
          dr(iact) = -de1(iact)/(eigval(iact) - lambda + tiniest)
        else
          dr(iact) = -de1(iact)/(eigval(iact) + lambda + tiniest)
        end if

        if (abs(eigval(iact)).lt.eps) then
          if (nvar.gt.1) dr(iact) = 0.0
        end if

      end do
      
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tropt_lambda(f,eigval,lambda)
! calculate hessian eigenvalue shift factor from 
! trust-region quadratic approximation.

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer maxiter
      parameter (maxiter = 100)

      double precision big,eps
      parameter (big = 1d6, eps = 1d-10)

      double precision f(nvar),eigval(nvar)

      logical error
      integer iact,iter,loop
      double precision check,dfn,dx,fn,itlamb,lambda,lower,upper
      double precision tmp1,tmp2
!
!     *** initialization ***
!
      error = .false.
      upper = big
      if (modsad.eq.1) then
        tmp1 = eigval(hfollow)
        eigval(hfollow) = big
        tmp2 = big
        do iact = 1, nvar
          tmp2 = min(tmp2,eigval(iact))
        enddo
        tmp2 = -tmp2
        eigval(hfollow) = tmp1
        lower = max(tmp1,tmp2)
        itlamb = lower + 0.05
      else
        if (eigval(1).ge.0.0) then
          itlamb = 0.05
          lower = 0.0
        else
          itlamb = -eigval(1) + 0.05
          lower = -eigval(1)
        end if
      end if
!
!     *** iteration loop ***
!
      iter = 0

      do loop=1,maxiter
!
!     *** calculate fn(lambda) function ***
!
        fn = 0.0
        do iact=1,nvar
          if ((modsad.eq.1).and.(iact.eq.hfollow)) then
            fn = fn + (f(iact)/(eigval(iact) - itlamb))**2
          else
            fn = fn + (f(iact)/(eigval(iact) + itlamb))**2
          end if
        end do
c
c     *** derivative of fn(lambda) with respect to lambda ***
c
        dfn = 0.0
        do iact=1,nvar
          if ((modsad.eq.1).and.(iact.eq.hfollow)) then
            dfn = dfn - f(iact)**2/(eigval(iact) - itlamb)**3
          else
            dfn = dfn + f(iact)**2/(eigval(iact) + itlamb)**3
          end if
        end do
!
!     *** calculate next estimate of lambda ***
!
        dx = fn/dfn
        lambda = itlamb + (sqrt(fn)/trust - 1.0)*dx
!
!     *** select new lambda value ***
!
        if ((lambda.le.upper).and.(lambda.gt.lower)) then

          if (lambda.lt.itlamb) then
            upper = itlamb
          else
            lower = itlamb
          end if

          itlamb = lambda

        else

          if (lambda.lt.itlamb) then
            upper = itlamb
          else
            lower = itlamb
          end if

          if (upper.eq.big) then
            itlamb = lambda + 0.05
          else
            itlamb = 0.5*(upper + lower)
          end if

        end if

        iter = iter + 1
        check = abs(sqrt(fn) - trust)
        if (check.le.eps) go to 100

      end do
!
!     *** iteration failed ***
!
      error = .true. 
      
      if (oprint) write (6,5000) 
5000  format ('lambda iteration did not converge')

 100  continue
!      
!     *** check lambda value ***
!
      if (modsad.eq.1) then
        lower = max(tmp1,tmp2)
        if (lambda.le.lower) then
          error = . true.
          if (oprint) write (6,5010) 
 5010     format ('lambda is not appropriate for ts search')
        end if
      else
        if ((eigval(1).lt.0.0).and.(lambda.lt.-eigval(1))) then
          error = .true.
          if (oprint) write (6,5020)
5020      format ('lambda < smallest hessian eigenvalue')
        else if ((eigval(1).gt.0.0).and.(lambda.lt.0.0)) then
          error = .true.
          if (oprint) write (6,5030) 
5030      format ('lambda < zero and h is positive definite')
        end if
      end if
!
!     *** lambda value is acceptable ***
!
      if (.not.error) then 
        if (oprint) write (6,5040) iter
5040    format ('lambda iteration converged in ',i3,' cycles')
      end if
!
!     *** print output ***
!
!
!     *** failed lambda calculation ***
!      
      if (error) then
        write (6,'(a,f16.8)') 'lambda', lambda
        write (6,'(a,f16.8)') 'trust radii', trust
        call errquit('tropt_lambda: fatal error', 911, geom_err)
      end if

      end 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccc
c $Id$
