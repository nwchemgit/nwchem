C> \ingroup task
C> @{
      logical function task_mepgs(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb
c     
      integer geom
      logical ignore
      logical firstpass
      double precision energyts
      logical status
      double precision start    ! Tracks time used in last step
      logical mepgs_freq, mepgs_opt
      external mepgs_freq, mepgs_opt
      logical  task_gradient, task_energy
      external task_gradient, task_energy
c
c     Disable printing to ecce of movecs after this point
c
      call movecs_ecce_print_off()
      firstpass = .true.
c
 1000 continue
c
      stotal = 0d0
c
c     Read input, load /cmepgs/, get geometry
c
      call mepgs_initialize(rtdb, geom, firstpass)
c
c     **** Obtain initial energy ****
c
      if (firstpass) then
        if (.not. task_energy(rtdb))
     $      call errquit('mepgs: task_energy failed',0, GEOM_ERR)
        if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $      call errquit('mepgs: could not get energy',0, RTDB_ERR)
        energyts = energy
      end if
c
c     **** Print trayectory file ****
c
      energy = energyts
      call mepgs_path(geom, .true., .false.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Displace TS along selected mode ! has to be negative          c
c     Assumes the user is sure on the TS nature and that a          c
c     preferently a frequency analysis has already been performed   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   
c     **** Deallocate geom ****
c
      if (.not.geom_destroy(geom))
     &  call errquit('hnd_opt: geom_destroy?', 911, GEOM_ERR)
c
c     **** Move away from TS *****
c
      if (firstpass) then
        if (.not. mepgs_freq(rtdb))
     $    call errquit('mepgs: mepgs_freq failed',0, GEOM_ERR)
      end if
c
c     *** Reallocate geom info ****
c
      if (.not. geom_create(geom, 'geometry'))
     &  call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
c
c     **** Select side to traverse **** 
c
      if (forward) then
        if (.not. geom_rtdb_load(rtdb, geom, 'ircforward'))
     &    call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
        if (ga_nodeid().eq.0) write(6,5000) 
5000    format(/,10x,24('-'),/,10x,'Forward IRC optimization',
     $         /,10x,24('-'))
      else if (backward) then 
        if (.not. geom_rtdb_load(rtdb, geom, 'ircbackward'))
     &    call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
        if (ga_nodeid().eq.0) write(6,5100) 
5100    format(/,10x,25('-'),/,10x,'Backward IRC optimization',
     $         /,10x,25('-'))
      end if
c
c     **** Store selected side ****
c
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $    call errquit('gsopt_energy_step: grs?',geom, RTDB_ERR)
      if (.not. geom_ncent(geom,nat))
     $  call errquit('hnd_opt: natoms?',nat, GEOM_ERR)
        call grad_active_atoms(rtdb, nat, oactive, nactive)
      if (.not. geom_systype_get(geom, isystype))
     $  call errquit('mepgs: systype?',0, GEOM_ERR)
c
c     **** Energy and Gradient ****
c
      call mepgs_gra(rtdb, geom)
      if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $    call errquit('mepgs: could not get energy',0, RTDB_ERR)
c     
c     **** Construct projector   ****
c     
      call gsopt_cart_pmat(rtdb, geom)
      call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
c
c     **** Print trayectory file ****
c
      call mepgs_path(geom, .false., .false.)
c
c     **** Check energy decrease agreement ****
c
      if (ga_nodeid().eq.0) write(6,6000) -evib
6000  format(/,10x,'Expected change in energy', 10x, f12.6)
      if (ga_nodeid().eq.0) write(6,6100) energy - energyts
6100  format(/,10x,'Obtained change in energy', 10x, f12.6)
      if ((energyts - energy) .lt. 0.0d0) then
        if (ga_nodeid().eq.0) write(6,6200) 
6200    format(/,1x,25('-')/,1x,'The energy has increased'/,1x,25('-') )
        ircdone = .false.
        goto 2000
      end if
c
c     **** Obtain displacement step ****
c
      call gsopt_compute_actual_step(geom)
c     
c     **** Read Initial hessian for IRC and OPT loops ****
c     ****  to be able to handle separate updates     ****
c
      call mepgs_hss_init(rtdb,geom)
c
c     **** initialization ****
c
      flip = .false.
      ircdone = .false.
c   
c     **** Deallocate geom ****
c
      if (.not.geom_destroy(geom))
     &   call errquit('hnd_opt: geom_destroy?', 911, GEOM_ERR)
c
c     ****  Gonzalez & Schlegel Iterative loop ****
c
      if (.not. mepgs_opt(rtdb))
     $   call errquit('mepgs: could not optimize',0, RTDB_ERR)
c
c     **** Obtain second side ****
c
      if (ircboth) then
        forward   = .false.
        ircboth   = .false.
        backward  = .true.
        firstpass = .false.
        goto 1000
      end if
c
 2000 continue
c
      task_mepgs=ircdone
c
      call ga_sync()
c
      end
C> @}
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      subroutine mepgs_initialize(rtdb, geom, start)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb
      integer geom              ! [output] 
      logical start              ! [output] 
c     
c     This routine initializes the common /cmepgs/ and
c     also creates and returns the geometry handle
c     
      integer i, j, num, ma_type, nactive_atoms, l_actlist
      logical ignore
      character*80 title
      character*8 source, test
      character*32 theory
      logical gsopt_geom_cart_coords_get
c     
      call util_print_push
      call util_print_rtdb_load(rtdb, 'mepgs')
      call ecce_print_module_entry('mepgs')
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
c     ----- parameters for optimization mepgs -----
c     
      if (.not. rtdb_get(rtdb,'mepgs:evib',mt_dbl,1,evib))
     $     evib = 0.0001d0
      if (.not. rtdb_get(rtdb,'mepgs:stride',mt_dbl,1,stride))
     $     stride = 0.1d0
      if (.not. rtdb_cget(rtdb,'mepgs:xyz',1,xyz))
     $     xyz = ' '
      if (.not. rtdb_get(rtdb,'mepgs:inhess',mt_int,1,inhess))
     $     inhess=0
      if (.not. rtdb_get(rtdb,'mepgs:nircopt',mt_int,1,nircopt)) 
     $     nircopt=250
      if (.not. rtdb_get(rtdb,'ircgs:mswg',mt_log,1,mswg))
     $     mswg = .false.
c
c     **** Select side to traverse ****
c
      if (start) then
        ircboth  = .true.
        forward  = .true.
        if (rtdb_get(rtdb,'mepgs:backward',mt_log,1,backward)) then 
             backward = .true.
             forward  = .false.
             ircboth  = .false.
        end if
        if (rtdb_get(rtdb,'mepgs:forward',mt_log,1,forward)) then 
             backward = .false.
             forward  = .true.
             ircboth  = .false.
        end if
c     
c     Save a  copy of the TS geometry 
c     
        if (ga_nodeid() .eq. 0) then
           ignore = rtdb_parallel(.false.)
           if (.not. geom_create(geom, 'geometry'))
     &          call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
           if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &          call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
           if (.not. geom_rtdb_store(rtdb, geom, 'tsreference'))
     &          call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
           if (.not. geom_destroy(geom))
     $          call errquit('mepgs: geom_destroy?',0, GEOM_ERR)
           ignore = rtdb_parallel(.true.)
        endif
      end if
      call ga_sync()
c     
c     Load the geometry info
c     
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
      if (start) then
        if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &      call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
      else if (.not. start) then
        if (.not. geom_rtdb_load(rtdb, geom, 'tsreference'))
     &      call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
      end if

      if (.not. geom_ncent(geom,nat))
     $     call errquit('hnd_opt: natoms?',nat, GEOM_ERR)
CJMC mass
      if (start .and. mswg) then
        if (.not. geom_masses_get(geom, nat, atmass))
     &      call errquit('ircgs: geom_masses_get failed',911, GEOM_ERR)
      end if
CJMC mass
c
      call grad_active_atoms(rtdb, nat, oactive, nactive)
      if (.not. geom_systype_get(geom, isystype))
     $     call errquit('mepgs: systype?',0, GEOM_ERR)
c
      if (oprint) then
        if (ga_nodeid().eq.0) then
         write(6,1) evib, stride, nircopt, inhess, mswg
 1       format(
     $      ' energy decrease                  (evib) = ', 1p,d9.1,0p,/,
     $      ' initial stride                 (stride) = ', 1p,d9.1,0p,/,
     $      ' maximum number of steps       (nircopt) = ', i4,/,
     $      ' initial hessian option         (inhess) = ', i4,/,
     $      ' mass weight coordinates          (mswg) = ', l4)
         write(6,9994)
 9994    format(/,10x,36('-'),
     1          /,10x,'Gonzalez & Schlegel IRC Optimization',
     2          /,10x,36('-'),/)
c
         call util_flush(6)
        end if
      endif
c     
c     Nvar is the no. of variables in the optimization
c
c     If we are optimizing the unit cell parameters then we pretend
c     there there are 3 more atoms which will parameterize the
c     unit cell.
c     
      nat_real = nat
      ncart = 3*nat
      nvar = ncart
      call gsopt_cart_pmat(rtdb, geom)
c
      energy    = 0d0
      energyref = 0d0
      alpha     = 1d0
      gmax      = 0d0
      grms      = 0d0
      smax      = 0d0
      srms      = 0d0
      xmax      = 0d0
      xrms      = 0d0
      call dfill(max_nvar, 0d0, ds, 1)
      call dfill(max_nvar, 0d0,dsp, 1)
      call dfill(max_nvar, 0d0, gx, 1)
      call dfill(max_nvar, 0d0, gq, 1)
      call dfill(max_nvar, 0d0,  g, 1)
      call dfill(max_nvar, 0d0, oldgra, 1)
      call dfill(max_nvar, 0d0, sp, 1)
c
      if (.not. gsopt_geom_cart_coords_get(geom, sp))
     $        call errquit('mepgs: geom?',0, GEOM_ERR)
c     
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccc
      subroutine mepgs_proj_grad(pgref)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "mafdecls.fh"
      double precision pgref(nvar) ! returns projected gradient
c
c     Nothing else is changed.
c
      integer l_pmat, k_pmat
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'pmat',
     $     l_pmat, k_pmat)) call errquit
     $     ('mepgs_proj_h_g: memory for work',nvar**2, MA_ERR)
c
      call geom_hnd_get_data('p',dbl_mb(k_pmat), nvar**2)
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('freddy',0, MA_ERR)
c
c     PG
c
      call dgemv('n',nvar, nvar, 1d0, dbl_mb(k_pmat), nvar,
     $     g, 1, 0d0, pgref, 1)
c
      if (.not. ma_chop_stack(l_pmat)) call errquit
     $     ('mepgs_p_h_g:ma?',0, MA_ERR)
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
      subroutine mepgs_cent(rtdb, geom, geoma, pgref, sfactor, string)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb, geom, geoma
      double precision pgref(nvar)
      double precision sfactor
      character*(*) string 
c
c     Update the geometry in cent and in the database
c     'center' by taking the step
c     Update the geometry in geom and in the database
c     'geometry' by taking the step
c
      double precision pgnorm
      double precision xold(max_cart), xnew(max_cart)
      integer i, iat, l_bi, k_bi
      logical gsopt_geom_cart_coords_get
      logical gsopt_geom_cart_coords_set
      logical ophigh
      double precision ddot
      external ddot
c
      ophigh = util_print('high', print_high)
c
c
      pgnorm = sqrt(ddot(nvar, pgref, 1,  pgref, 1))
      call dcopy(nvar, pgref, 1, ds, 1) 
      call dscal(nvar, sfactor*stride/pgnorm, ds, 1)
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
      if (mswg) then
        call mwcoord(ds, nvar, .false.) 
      end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
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
      call ga_brdcst(1,ds,ma_sizeof(mt_dbl,nvar,mt_byte),0)
c
c     Get original coordinates
c
      if (.not. gsopt_geom_cart_coords_get(geom, xold))
     $     call errquit('mepgs_energy_step: coordinates?',geom,
     &       GEOM_ERR)
c
      call dcopy(ncart, ds, 1, xnew, 1) 
c
      call daxpy(ncart, 1.0d00, xold, 1, xnew, 1)
c     FRACTIONAL?
      if (.not. gsopt_geom_cart_coords_set(geoma, xnew))
     $    call errquit('mepgs_energy_step: coordinates?', string,
     &    GEOM_ERR)
c
      if (.not. geom_rtdb_store(rtdb, geoma, string))
     $     call errquit('mepgs_energy_step: grs?',geom, RTDB_ERR)
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function mepgs_cosang(avec,bvec,angle)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "util.fh"
#include "mafdecls.fh"
c
      logical angle
      double precision avec(nvar), bvec(nvar)
c
      double precision ctheta, factor(2)
      double precision ddot
      external ddot
c
      mepgs_cosang = 0.0d0
c
      factor(1) = ddot(nvar, avec, 1, bvec, 1)
c
      factor(2) = ddot(nvar, avec, 1, avec, 1)*
     $            ddot(nvar, bvec, 1, bvec, 1)
c
      factor(2) = sqrt(factor(2))
c
      ctheta = factor(1)/factor(2)
c
      if (abs(ctheta).gt.1.0d0) ctheta = dsign(1.d0,ctheta) 
c
      if (angle) then
        mepgs_cosang = acos(ctheta)
      else
        mepgs_cosang = ctheta
      end if
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mepgs_hss_init(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
c
      integer rtdb, geom

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
     $     ('mepgs_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      old_hessian=.false.
      call gsopt_check_hess(nvar, old_hessian)
      if (oprint) write(6,*)
      if (old_hessian) then
        call mepgs_hess_cart_guess()
        if (oprint) write(6,*)
     $     ' Using Cartesian Hessian from previous frequency',
     $     ' calculation'
      else
        if (oprint) write(6,*) 'Not restart Hessian? '
      endif
c
c     Apply constants, constraints and overall scaling
c
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_c, k_c)) call errquit
     $     ('mepgs_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian',
     $     l_t, k_t)) call errquit
     $     ('mepgs_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      call geom_hnd_get_data('irc.hess', dbl_mb(k_hess), nvar*nvar)
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
      if (.not. rtdb_get(rtdb,'gsopt:hscale',mt_dbl,1,hscale))
     $     hscale = 1d0
      call dscal(nvar*nvar, hscale, dbl_mb(k_hess), 1)
      if (oprint .and. hscale.ne.1d0) then
         if (ga_nodeid().eq.0)  write(6,78) hscale
      end if
 78   format(' Scaling initial hessian by ',f6.2)
c
      call geom_hnd_put_data('irc.hess',dbl_mb(k_hess), nvar*nvar)
CJMC prepare OPT start hessian
      call geom_hnd_put_data('gsopt.hess',dbl_mb(k_hess), nvar*nvar)
CJMC prepare OPT start hessian
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('mepgs_init_hess ma corrupt',0, MA_ERR)
c
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mepgs_hess_cart_guess()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
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
      double precision x,factor
      integer i,j,iat,jat
      integer l_bi, k_bi, l_hc, k_hc, l_hq, k_hq
c     
      if (.not. ma_push_get(mt_dbl, ncart*nvar, 'binv',
     $     l_bi, k_bi)) call errquit
     $     ('mepgs_hess_cart_guess: ma?', ncart*nvar, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart',
     $     l_hc, k_hc)) call errquit
     $     ('mepgs_hess_cart_guess: ma?', ncart**2, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart2',
     $     l_hq, k_hq)) call errquit
     $     ('mepgs_hess_cart_guess: ma?', nvar**2, MA_ERR)
c     
      if (ga_nodeid().eq.0) then
         call util_file_name('hess',.false.,.false.,fname)
         open(unit=h_unit,file=fname,form='formatted',status='unknown',
     $        err=99990,access='sequential')
         rewind h_unit
         do i = 1,ncart
            iat = (i+2)/3
            do j = 1,i
               jat = (j+2)/3
               read(h_unit,10000,err=99992,end=99992) x

cDMR mass-weight Hessian
               if (mswg) then
                 factor = sqrt(atmass(iat))*sqrt(atmass(jat))
                 x = x/factor
               endif
cDMR

               dbl_mb(k_hc+(i-1)*ncart+(j-1)) = x
               dbl_mb(k_hc+(j-1)*ncart+(i-1)) = x
            enddo
         enddo
         close(unit=h_unit,status='keep')
      endif

      call ga_brdcst(1,dbl_mb(k_hc),
     $     ma_sizeof(mt_dbl,ncart**2,mt_byte),0)
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
      call geom_hnd_put_data('irc.hess',dbl_mb(k_hc), nvar**2)
c
      if (.not. ma_chop_stack(l_bi))
     $     call errquit('mepgs_hess_cart_guess: ma corrupt?',0, MA_ERR)
c
      return
10000 format(f30.15)
99990 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as unknown file'
      call errquit('mepgs_hess_cart: fatal error', 911, GEOM_ERR)
99991 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as new file'
      call errquit('mepgs_hess_cart: fatal error', 911, GEOM_ERR)
99992 write(6,*)' error in reading <',fname(1:inp_strlen(fname)),
     $     '> as hessian file'
      call errquit('mepgs_hess_cart: fatal error', 911, GEOM_ERR)
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mepgs_path(geom, openfile, closefile)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
      integer geom
      logical openfile, closefile
      character*255 filename, dir
      logical mol_geom_print_xyz
      external mol_geom_print_xyz
c
c     Print a trajectory file
c 
      if (ga_nodeid().eq.0 .and. xyz.ne.' ') then
         dir      = ' '
         filename = ' '
         call util_directory_name(dir, .false., 0)
CJMC
         if (forward) then
           write(filename,12) dir(1:inp_strlen(dir)),
     $          xyz(1:inp_strlen(xyz))
 12        format(a,'/',a,'.fxyz')
         else if (backward) then
           write(filename,13) dir(1:inp_strlen(dir)),
     $          xyz(1:inp_strlen(xyz))
 13        format(a,'/',a,'.bxyz')
         end if
CJMC
         if (openfile) then
           open(88,file=filename,form='formatted',status='unknown',
     $          access='sequential',err=133)
           rewind(88)
         end if
c
         if (.not. mol_geom_print_xyz(geom, 88, energy))
     $       call errquit('mepgs_path: mol_geom_print_xyz?',0, GEOM_ERR)
         call util_flush(88)
c
         if (closefile) close(88,status='keep',err=133)
      end if
c
      return
 133  call errquit('mepgs_path: error open/close xyz file',0, GEOM_ERR)
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
      subroutine mepgs_gra(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "rtdb.fh"
#include "geom.fh"
c
      integer rtdb, geom
c
      integer iat, ixyz
      logical ophigh
      logical  task_gradient
      external task_gradient
c
      ophigh = util_print('high', print_high)
c
      if (.not. task_gradient(rtdb))
     $    call errquit('mepgs: task_gradient failed',0, GEOM_ERR)
      call gsopt_get_grad(rtdb, geom) ! Into gx
c     
c     Zero the gradient associated with atoms frozen in cartesians
c     
      if (ga_nodeid().eq.0.and.ophigh) 
     $    write(6,*) 'Zeroing constrained gradient'
      if ((.not. zcoord) .and. (nactive .ne. nat_real)) then
        do iat = 1, nat_real
          if (.not. oactive(iat)) then
            do ixyz = 1, 3
               gx((iat-1)*3+ixyz) = 0.0
            end do
          end if
        end do
      end if
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mepgs_chk(geom, irccyc)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "geom.fh"
c     
      integer geom, irccyc
c
      integer ivar
      double precision abvec(max_nvar) 
      double precision bcvec(max_nvar) 
      double precision newx(max_nvar) 
      double precision ogrms
      double precision graang, cosang, stpang, rmsdif
      double precision eps
      parameter (eps = 1d-3)
      double precision echange
      double precision stemp
c
      character*1 mark
      character*72 string
c
      double precision mepgs_cosang
      logical gsopt_geom_cart_coords_get
c

      if (.not. gsopt_geom_cart_coords_get(geom, newx))
     $    call errquit('mepgs: geom?',0, geom_err) 
c
      call dcopy(nvar, oldgeo, 1, abvec, 1) 
      call daxpy(nvar, -1.0d0, center, 1, abvec, 1)
c
      call dcopy(nvar, newx, 1, bcvec, 1) 
      call daxpy(nvar, -1.0d0, center, 1, bcvec, 1)
c
      cosang = mepgs_cosang(oldgra, g, .false.)
      stpang = mepgs_cosang(abvec, bcvec, .true.)
      svalue = abs(stpang*stride/2d0)
      stemp = stotal + svalue
c
      grms = 0d0
      ogrms = 0d0
      srms = 0d0
      gmax  = 0d0
      smax  = 0d0
c
      do ivar = 1, nvar
        ogrms = ogrms + oldgra(ivar)*oldgra(ivar)
        grms  =  grms +  g(ivar)*g(ivar)
        gmax  = max(gmax, abs(g(ivar)))
      enddo
      ogrms = sqrt(ogrms/dble(nvar))
      grms  = sqrt(grms/dble(nvar))
c
      rmsdif = grms - ogrms
      if ((rmsdif.lt.0d0).and.(stemp.gt.1d0)) flip = .true.
      if (abs(rmsdif).lt.eps) rmsdif = 0d0
c
      echange = energy - energyref
c
cDMR check if we are already flat
c
      if ((echange.gt.0.01d0*eps).and.(stemp.gt.1d0))
     $  nbounded = nbounded + 1 
      if ((stride.le.0.002d0).and.(stemp.gt.1d0)) nbounded = nbounded+1
c
      if (flip .and. (grms.le.grms_tol) .and. (irccyc.gt.5)) then
        if ((gmax.le.gmax_tol).and.(abs(echange).lt.eprec)) then
          ircdone = .true.
        elseif (abs(echange).lt.0.1d0*eprec) then
          ircdone = .true.
        elseif (stride.le.0.002d0) then
          ircdone = .true.
        endif
      endif

      if(nbounded.gt.5) then
        ircdone=.true.
      endif
c
c       Assume step was good
c
      redogs = .false.
c
      if(echange.gt.eps) then
        redogs= .true.
      elseif(stpang.le.eps) then
        redogs =.true.
      end if
c
      if (.not. redogs) then
        stotal = stotal + svalue
        mark = '&'
        if (ga_nodeid().eq.0) then
          write(luout,1)
          write(luout,2) mark,irccyc-1,energy,svalue,stotal,grms,gmax,
     $                   echange              
          string = repeat(" ",72)
          if (grms.lt.grms_tol) write(string(47:48),"(A2)") "ok"
          if (gmax.lt.gmax_tol) write(string(55:56),"(A2)") "ok"
          if (abs(echange).lt.eprec) write(string(66:67),"(A2)") "ok"
          write(luout,3) string
        endif
 1      format(
     $   /,'   Point      Energy     svalue   stotal    grms  ',
     $        '   gmax      Delta E  ',
     $   /,'  -------  ------------ -------- -------- --------',
     $        ' -------- ------------')
 2      format(
     $        a1,i5,f17.8,4f9.5,2x,1p,d10.2)
 3      format(a72,/)       
      end if
c
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      subroutine mepgs_dist(geom)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "geom.fh"
c     
      integer geom
c
      integer ivar
      double precision newx(max_nvar) 
      logical gsopt_geom_cart_coords_get
c
      if (.not. gsopt_geom_cart_coords_get(geom, newx))
     $    call errquit('mepgs: geom?',0, geom_err) 
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
      if (mswg) then
        call mwcoord(center, nvar, .true.)
        call mwcoord(  newx, nvar, .true.)
      end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
      distcx = 0d0
      do ivar=1, nvar
        distcx = distcx + (center(ivar) - newx(ivar))**2
      end do
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
      if (mswg) then
        call mwcoord(center, nvar, .false.)
        call mwcoord(  newx, nvar, .false.)
      end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function mepgs_opt(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb
c     
      integer geom
      integer reference
      integer imepgs
      logical ignore
      double precision mepgs_cosang
      double precision start    ! Tracks time used in last step
cjmc
      double precision minst    ! Tracks time used in last step
      parameter(minst = 0.002d0)
      double precision pgref(max_nvar) 
      logical gsopt_geom_cart_coords_get
      logical gsopt_geom_cart_coords_set
      logical gsopt
cjmc
      integer required          ! Time required
      logical ophigh
      logical  task_gradient, task_energy
      external task_gradient, task_energy
c
c     **** Initialization ****
c
      if (.not. geom_create(geom, 'geometry'))
     &  call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &  call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
      if (.not. geom_ncent(geom,nat))
     $  call errquit('hnd_opt: natoms?',nat, GEOM_ERR)
        call grad_active_atoms(rtdb, nat, oactive, nactive)
      if (.not. geom_systype_get(geom, isystype))
     $  call errquit('mepgs: systype?',0, GEOM_ERR)
      call gsopt_get_grad(rtdb, geom) ! Into gx
      nbounded = 0
c
c     ****  Gonzalez & Schlegel Iterative loop ****
c
      do imepgs = 1, nircopt+1   
c
c     **** Iteration exceded ****
c
        if (imepgs.gt.nircopt) goto 200
c
        start = util_wallsec()
        if (oprint.and.ga_nodeid().eq.0) write(6,1) imepgs-1
 1      format(/,10x,11('-'),/,10x,'GS Step',i4,/,10x,11('-'))
        if ((ga_nodeid() .eq. 0) .and.
     $       util_print('geometry',print_default)) then
           if (.not. geom_print(geom)) call errquit('mepgs: geom?',0,
     $       GEOM_ERR)
        endif
c    
c     Save old energy, gradient and coordinates
c
        energyref = energy    
        call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
        call dcopy(nvar, g, 1, oldgra, 1) 
        call dcopy(nvar, g, 1, gp, 1) 
        if (.not. gsopt_geom_cart_coords_get(geom, oldgeo))
     $      call errquit('mepgs: geom?',0, geom_err) 
c
c       Project gradient
c
        if (mswg) call mwgrad(g, nvar, .true.)
        call mepgs_proj_grad(pgref)
        if (mswg) call mwgrad(g, nvar, .false.)
c
 1000   continue
c
c     **** Allocate and prepare field for reference (center)  ****
c
          if (.not. geom_create(reference, 'center'))
     &        call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
          if (.not. geom_rtdb_load(rtdb, reference, 'geometry'))
     &        call errquit('hnd_opt: no initial geometry ',911,RTDB_ERR)
c
c     **** Calculate center, half stride ****
c
          call mepgs_cent(rtdb, geom, reference,pgref, -0.5d0, 'center')
c
c     **** Obtain center coordinates ****
c
          if (.not. gsopt_geom_cart_coords_get(reference, center))
     $        call errquit('mepgs: geom?',0, geom_err) 
c
c     **** Deallocate reference (center) ****
c
          if (.not. geom_destroy(reference))
     $        call errquit('mepgs:reference corrupt',0, GEOM_ERR)
c
c     **** Calculate x^(l), complete stride ****
c
          call mepgs_cent(rtdb, geom, geom, pgref, -1.0d0, 'geometry')
c
c     **** Energy and Gradient ****
c
          call mepgs_gra(rtdb, geom)
          if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $      call errquit('gsopt: could not get energy',0, RTDB_ERR)
          call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
c
c     **** Update the OPT hessian ****
c
          call dcopy(nvar, oldgeo, 1, sp, 1)
          call gsopt_compute_actual_step(geom)
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
          if (mswg) then
            call mwcoord( ds, nvar, .true.)
            call mwgrad(   g, nvar, .true.)
            call mwgrad(  gp, nvar, .true.)
          end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
          call gsopt_hessian_update()
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
          if (mswg) then
            call mwcoord( ds, nvar, .false.)
            call mwgrad(   g, nvar, .false.)
            call mwgrad(  gp, nvar, .false.)
          end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
c     **** Compute initial distance ****
c
          call mepgs_dist(geom)
          ctrust2 = distcx
c   
c     **** Deallocate geom ****
c
         if (.not.geom_destroy(geom))
     &     call errquit('hnd_opt: geom_destroy?', 911, GEOM_ERR)
c
c     **** Perform constrained optimization ****
c
          if (.not. gsopt(rtdb))
     $       call errquit('mepgs: gsopt failed',0, GEOM_ERR)
c
c    *** Reallocate geom info ****
c
          if (.not. geom_create(geom, 'geometry'))
     &      call errquit('hnd_opt: geom_create?', 911, GEOM_ERR)
          if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &      call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
          if (.not. geom_ncent(geom,nat))
     $      call errquit('hnd_opt: natoms?',nat, GEOM_ERR)
            call grad_active_atoms(rtdb, nat, oactive, nactive)
          if (.not. geom_systype_get(geom, isystype))
     $      call errquit('mepgs: systype?',0, GEOM_ERR)
c
c     **** Compute final distance ****
c
          call mepgs_dist(geom)
          if (abs(distcx - ctrust2) .gt. 1d-4) 
     $      call errquit('mepgs: norm not preserved', 911, GEOM_ERR)
c
c     **** obtain current gradient ****
c
          call gsopt_get_grad(rtdb, geom) ! Into gx
          call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
c
c     **** Compute convergence info ****
c
          call mepgs_chk(geom, imepgs)
c
c     **** reload previous fiels if necessary ****
c
          if (redogs) then
            if (oprint) write(6,*) "Rejecting the optimized Point"
            if (stride.le.minst) then 
              ircdone = .true. 
              goto 100
            end if
            stride = stride/2d0
            if (.not. gsopt_geom_cart_coords_set(geom, oldgeo))
     $         call errquit('gsopt: geom?',0, geom_err) ! reload previous cart
            if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $         call errquit('gsopt: grs?',geom, RTDB_ERR)
            goto 1000
          else
            call updstride(geom)
          end if
c
c    Check for convergence
c     
        if (ircdone) goto 100
c
c     Update the IRC hessian
c
        call dcopy(nvar, oldgeo, 1, sp, 1)
        call gsopt_compute_actual_step(geom)
c
c     Print a trajectory file
c
        call mepgs_path(geom, .false., .false.)
c
c     Check time before next iteration
c
        required = int(1.2d0*(util_wallsec() - start)) + 1
        if (.not. util_test_time_remaining(rtdb,required)) goto 200
c
c
      enddo                     ! End of iterative loop
c
c     **** Failed to converge ****
c
 200  if (oprint.and.ga_nodeid().eq.0) write(6,201)
 201  format(/,1x,63('-')/,1x,'Failed to converge in maximum number',
     $     ' of steps or available time'/,1x,63('-')/)
      ircdone = .false.
c
c     **** Procedure finished ****
c
 100  if (ircdone) then
         if (oprint.and.ga_nodeid().eq.0) write(6,101)
 101     format(/,6x,22('-'),/,6x,'IRC Optimization converged',/,
     $        6x,22('-'),/)
      endif
c
c     Print out final info and geometry
c
      if (ga_nodeid().eq.0 .and. util_print('finish',print_low)) then
        call mepgs_path(geom, .false., .true.)
        if (.not. geom_print(geom)) call errquit
     $     ('hnd_opt_drv: geom_print?',0, GEOM_ERR)
c
      end if
c
c    **** Save geometry ****
c
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $   call errquit('gsopt: grs?',geom, RTDB_ERR)
c
c     Clean up and go home
c
      if (.not.geom_destroy(geom))
     &     call errquit('hnd_opt: geom_destroy?', 911, GEOM_ERR)
c
      mepgs_opt=ircdone
      if (ircdone) then
         call ecce_print_module_exit('mepgs', 'ok')
      else
         call ecce_print_module_exit('mepgs', 'failed')
      endif
c
      call movecs_ecce_print_on() ! Restore MO printing
      call util_print_pop
c
      call ga_sync()
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC                GS optimization                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function gsopt(rtdb)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "rtdb.fh"
#include "geom.fh"
c     
      integer rtdb
c     
      integer geom, i, iat
      integer istep
      logical converged
      double precision start    ! Tracks time used in last step
cjmc
      double precision oldx(max_cart)
      logical gsopt_geom_cart_coords_get
      logical gsopt_geom_cart_coords_set
cjmc
      integer required          ! Time required
      logical ophigh
      logical  gsopt_converged, task_gradient
      external gsopt_converged, task_gradient
c
      ophigh = util_print('high', print_high)
c
c     Read input, load /cgsopt/, get geometry
c
      call gsopt_initialize(rtdb, geom)
c
c     Energy and Gradient
c
      call gsopt_get_grad(rtdb, geom) ! Into gx
      if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $    call errquit('gsopt: could not get energy',0, RTDB_ERR)
c     
c     Construct proyector
c
      call gsopt_cart_pmat(rtdb, geom)
      call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
c     
c     Initial hessian
c
      call gsopt_hss_init(rtdb,geom)
c
c     Iterative loop
c
      do istep = 1, nptopt+1    ! +1 since first pass thru loop is not a step
c
c     iteration exceded
c
        if(istep.gt.nptopt) goto 200
c
        start = util_wallsec()
        if (oprint.and. ga_nodeid().eq.0) write(6,1) istep-1
 1      format(/,10x,'--------',/,10x,'Step',i4,/,10x,'--------')
        if ((ga_nodeid() .eq. 0) .and.
     $       util_print('geometry',print_default)) then
           if (.not. geom_print(geom)) call errquit('gsopt: geom?',0,
     $       GEOM_ERR)
        endif
c
c     Compute step/gradient info and print for user
c     (for current energy & gradient, and the previous alpha*step).
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
        if (mswg) then
          call mwcoord(radius, nvar, .true.)
          call mwcoord(    ds, nvar, .true.)
          call  mwgrad(     g, nvar, .true.)
        end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
        call gsopt_compute_info()
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
        if (mswg) then
          call mwcoord(radius, nvar, .false.)
          call mwcoord(    ds, nvar, .false.)
          call  mwgrad(     g, nvar, .false.)
        end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
        call gsopt_print(geom, istep)
c     
c     Check for convergence
c     
        if (gsopt_converged(istep)) then
           converged = .true.
           goto 100
        endif
c
c     Save old energy, gradient and coordinates
C
        energyp = energy       ! Used for convergence and step restriction
        call dcopy(nvar, g, 1, gp, 1) ! Used for Hessian update
        if (.not. gsopt_geom_cart_coords_get(geom, oldx))
     $      call errquit('gsopt: geom?',0, geom_err) 
c
c     generate a new search direction
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
        if (mswg) call mwgrad(g, nvar, .true.)
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
        call gsopt_pickstp(rtdb, geom, istep) ! fills in ds()
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
        if (mswg) call mwgrad(g, nvar, .false.)
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c     
c     take recommended step.
c     
        call gsopt_take_step(rtdb, geom) ! Updates geom using ds
c
c     We have now taken a step.  Replace the approximate step taken
c     by the exact step in case update of internals was not exact.
c     
        call gsopt_compute_actual_step(geom)
c
c     Energy and Gradient
c
        call mepgs_gra(rtdb, geom)
        if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $      call errquit('gsopt: could not get energy',0, RTDB_ERR)
c     
c     Construct proyector
c     
        call gsopt_cart_pmat(rtdb, geom)
        call dcopy(ncart, gx, 1, g, 1) ! g() set to gx()
c
c     update the hessian
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
        if (mswg) then
          call mwcoord(ds, nvar, .true.)
          call mwgrad(  g, nvar, .true.)
          call mwgrad( gp, nvar, .true.)
        end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
        call gsopt_hessian_update()
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
        if (mswg) then
          call mwcoord(ds, nvar, .false.)
          call mwgrad(  g, nvar, .false.)
          call mwgrad( gp, nvar, .false.)
        end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
c     Check time before next iteration
c
        required = int(1.2d0*(util_wallsec() - start)) + 1
        if (.not. util_test_time_remaining(rtdb,required)) goto 200
c
c
      enddo                     ! End of iterative loop
c      istep = istep - 1         ! Since we fell out
 200  if (oprint.and.ga_nodeid().eq.0) write(6,201)
 201  format(/,1x,63('-')/,1x,'Failed to converge in maximum number',
     $     ' of steps or available time'/,1x,63('-')/)
      converged = .false.
c
 100  if (converged) then
         if (oprint.and.ga_nodeid().eq.0) write(6,101) 
 101     format(/,6x,22('-'),/,6x,'Optimization converged',/,
     $        6x,22('-'),/)
      endif
c
      if (ga_nodeid().eq.0 .and. util_print('finish',print_low)) then
c
c     Print out final info and geometry
c
         call gsopt_print(geom, istep)
         if (.not. geom_print(geom)) call errquit
     $        ('hnd_opt_drv: geom_print?',0, GEOM_ERR)
c
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
      gsopt=converged
      if (converged) then
         call ecce_print_module_exit('gsopt', 'ok')
      else
         call ecce_print_module_exit('gsopt', 'failed')
      endif
c
      call movecs_ecce_print_on() ! Restore MO printing
      call util_print_pop
c
      call ga_sync()
c
      end
c
      subroutine gsopt_hss_init(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
c
      integer rtdb, geom

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
     $     ('gsopt_init_hess: failed allocating hessian',nvar**2,
     &       MA_ERR)
c
      old_hessian=.false.
      if(inhess.ne.1) call gsopt_check_hess(nvar, old_hessian)
      if (oprint) write(6,*)
      if (old_hessian) then
        if (inhess.eq.2) then
           call gsopt_hess_cart_guess()
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
c     Cartesians are easy
c
      call dfill(nvar**2, 0.0d0, dbl_mb(k_hess), 1)
      call dfill(nvar, 0.5d0, dbl_mb(k_hess), nvar+1)
      if (isystype .ne. 0) then
         if (.not. geom_amatrix_get(geom, amat))
     $        call errquit('geom_frac_to_cart: a', 0, GEOM_ERR)
         do iat = 1, nat_real
            do i = 1, 3
               dbl_mb(k_hess + (iat-1)*3 + (i-1) +
     $              ((iat-1)*3 + (i-1))*nat*3) = 0.5*amat(i,i)**2
            end do
         end do
         if (odebug) then
           write(6,*) ' The initial hessian '
           call output(dbl_mb(k_hess),1,nvar,1,nvar,nvar,nvar,1)
         end if
      end if
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
      call geom_hnd_put_data('gsopt.hess', dbl_mb(k_hess), nvar*nvar)
c
c     Apply constants, constraints and overall scaling
c
 999  if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian', l_c, k_c)) 
     $   call errquit('gsopt_init_hess: failed allocating hessian',
     $                nvar**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hessian', l_t, k_t)) 
     $   call errquit('gsopt_init_hess: failed allocating hessian',
     $                nvar**2, MA_ERR)
c
      call geom_hnd_get_data('gsopt.hess', dbl_mb(k_hess), nvar*nvar)
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
      if (.not. rtdb_get(rtdb,'gsopt:hscale',mt_dbl,1,hscale))
     $     hscale = 1d0
      call dscal(nvar*nvar, hscale, dbl_mb(k_hess), 1)
      if (oprint .and. hscale.ne.1d0) then
         if (ga_nodeid().eq.0) write(6,78) hscale
      end if
 78   format(' Scaling initial hessian by ',f6.2)
c
      call geom_hnd_put_data('gsopt.hess',dbl_mb(k_hess), nvar*nvar)
c
      if (.not. ma_chop_stack(l_hess)) call errquit
     $     ('gsopt_init_hess ma corrupt',0, MA_ERR)
c
      end
      subroutine gsopt_check_hess(nvar, old_hessian)
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
      call util_file_name('gsopt.hess',.false.,.false.,filename)
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
      subroutine gsopt_del_hess()
      implicit none
#include "util.fh"
c
c     Delete the Hessian information restart file.
c     
      character*255 opt_hess_fil
c
      call util_file_name('gsopt.hess',
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
      subroutine gsopt_initialize(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb
      integer geom              ! [output] 
c     
c     This routine initializes the common /cgsopt/ and
c     also creates and returns the geometry handle
c     
      integer i, j, num, ma_type, nactive_atoms, l_actlist
      integer ivar
      logical ignore
      character*80 title
      character*8 source, test
      character*32 theory
      logical gsopt_geom_cart_coords_get
c     
      call util_print_push
      call util_print_rtdb_load(rtdb, 'gsopt')
      call ecce_print_module_entry('gsopt')
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
c     ----- parameters for optimization gsopt -----
c     
      if (.not. rtdb_get(rtdb,'gsopt:modsad',mt_int,1,modsad)) 
     $     modsad=0
cjmc
      trust = sqrt(ctrust2)
cjmc
      if (.not. rtdb_cget(rtdb,'mepgs:xyz',1,xyz))
     $     xyz = ' '
      if (.not. rtdb_get(rtdb,'gsopt:eprec',mt_dbl,1,eprec)) then
         if (.not. rtdb_cget(rtdb,'task:theory',1,theory))
     $        theory = ' '
         if (inp_compare(.false.,theory,'dft')) then
            eprec = 5.0d-6
         else
            eprec = 1.0d-6
         endif
      endif
      if (.not. rtdb_get(rtdb,'gsopt:gmax_tol',mt_dbl,1,gmax_tol))
     $     gmax_tol = 0.00045d0
      if (.not. rtdb_get(rtdb,'gsopt:grms_tol',mt_dbl,1,grms_tol))
     $     grms_tol = 0.0003d0
      if (.not. rtdb_get(rtdb,'gsopt:xmax_tol',mt_dbl,1,xmax_tol))
     $     xmax_tol = 0.0018d0
      if (.not. rtdb_get(rtdb,'gsopt:xrms_tol',mt_dbl,1,xrms_tol))
     $     xrms_tol = 0.0012d0
      if (.not. rtdb_get(rtdb,'gsopt:nptopt',mt_int,1,nptopt))
     $     nptopt=50
      if (.not. rtdb_get(rtdb,'gsopt:inhess',mt_int,1,inhess))
     $     inhess=0
      if (.not. rtdb_get(rtdb,'gsopt:linopt',mt_int,1,linopt)) 
     $     linopt=1
      if (.not. rtdb_get(rtdb,'gsopt:moddir',mt_int,1,moddir)) 
     $     moddir=0
      if (.not. rtdb_get(rtdb,'gsopt:modsad',mt_int,1,modsad)) 
     $     modsad=0
      if (.not. rtdb_get(rtdb,'gsopt:sadstp',mt_dbl,1,sadstp)) 
     $     sadstp=0.1d0
      if (.not. rtdb_get(rtdb,'gsopt:oqstep',mt_log,1,oqstep))
     $     oqstep = .true.
      if (.not. rtdb_get(rtdb,'gsopt:modupd',mt_int,1,modupd)) then
         if (modsad .eq. 0) then
            modupd = 1          ! BFGS update for minimization
         else
            modupd = 2          ! PSB update for saddle point
         endif
      endif
      if (.not. rtdb_get(rtdb,'gsopt:ocheckgrad',mt_log,1,ocheckgrad))
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
         if (.not. geom_rtdb_store(rtdb, geom, 'gsoptinitial'))
     &        call errquit('hnd_opt: no geometry ', 911, RTDB_ERR)
         if (.not. geom_destroy(geom))
     $        call errquit('gsopt: geom_destroy?',0, GEOM_ERR)
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
     $     call errquit('gsopt: systype?',0, GEOM_ERR)
c
      if (oprint.and.ga_nodeid().eq.0) then
         write(6,1) gmax_tol, grms_tol, xmax_tol, xrms_tol, 
     $        eprec,
     $        nptopt, inhess, modupd
 1       format(
     $        ' maximum gradient threshold         (gmax) = ', f10.6,/,
     $        ' rms gradient threshold             (grms) = ', f10.6,/,
     $        ' maximum cartesian step threshold   (xmax) = ', f10.6,/,
     $        ' rms cartesian step threshold       (xrms) = ', f10.6,/,
     $        ' energy precision                  (eprec) = ', 1p,d9.1,
     $        0p,/,
     $        ' maximum number of steps          (nptopt) = ', i4,/,
     $        ' initial hessian option           (inhess) = ', i4,/,
     $        ' hessian update option            (modupd) = ', i4,/)
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
      nvar = ncart
c     
c
      call gsopt_cart_pmat(rtdb, geom)
c
      energy = 0d0
      energyp= 0d0
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
      call dfill(max_nvar, 0d0, radius, 1)
c
      if (.not. gsopt_geom_cart_coords_get(geom, sp))
     $    call errquit('gsopt: geom?',0, GEOM_ERR)
      call dcopy(nvar, sp, 1, radius, 1)
      call daxpy(nvar, -1.0d0, center, 1, radius, 1)
c     
      end
      subroutine gsopt_hessian_update()
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
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
     $     ('gsopt_hessian_update: zero step?',0, GEOM_ERR)
      call dscal(nvar, alpha, ds, 1)
c     
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('gsopt_hessian_update: memory for hessian',nvar**2,
     &       GEOM_ERR)
      call geom_hnd_get_data('gsopt.hess',dbl_mb(k_hess), nvar**2)
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
         if(dsdg.gt.1d-8) then
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
      call geom_hnd_put_data('gsopt.hess',dbl_mb(k_hess), nvar**2)
      call geom_hnd_put_data('irc.hess',dbl_mb(k_hess), nvar**2)
      if (.not. ma_pop_stack(l_hess)) call errquit
     $     ('gsopt_hessian_update: ma?',0, MA_ERR)
c
c     
      end
c
      subroutine gsopt_take_step(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
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
      logical gsopt_geom_cart_coords_get
      logical gsopt_geom_cart_coords_set
c
c     Get original coordinates
c
c     FRACTIONAL?
      if (.not. gsopt_geom_cart_coords_get(geom, xold))
     $     call errquit('gsopt_energy_step: coordinates?',geom,
     &       GEOM_ERR)
c
c     Take the step
c  
       call dcopy(ncart, ds, 1, xnew, 1)
       call dscal(nvar, alpha, xnew, 1)
c       call sym_grad_symmetrize(geom, xnew)
       call daxpy(ncart, 1.0d0, xold, 1, xnew, 1)
c     FRACTIONAL?
       if (.not. gsopt_geom_cart_coords_set(geom, xnew))
     $     call errquit('gsopt_energy_step: coordinates?',geom,
     &                  GEOM_ERR)
c
c     Must ensure the geometry has the required symmetry even after
c     enforcing it on the step.  Should use
c     an error criterion consistent with the step size.
c
c      call sym_geom_project(geom, trust)
c
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     $     call errquit('gsopt_energy_step: grs?',geom, RTDB_ERR)
c
c     Compute the maximum and RMS cartesian displacements
c    
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
      if (mswg) then
        call mwcoord(xold, nvar, .true.)
        call mwcoord(xnew, nvar, .true.)
      end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
      xmax = 0d0
      xrms = 0d0
      do i = 1, ncart
         xmax = max(xmax, abs(xold(i)-xnew(i)))
         xrms = xrms + (xold(i)-xnew(i))**2
      enddo
      xrms = sqrt(xrms/dble(ncart))
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
      if (mswg) then
        call mwcoord(xold, nvar, .false.)
        call mwcoord(xnew, nvar, .false.)
      end if
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
      end
      subroutine gsopt_print(geom, istep)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"

      integer geom, istep
c
c     Print out stuff
c
      integer i
      double precision de
      character*9 cvg1, cvg2, cvg3, cvg4
      character*1 mark
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
         if (ga_nodeid().eq.0) write(6,1) mark, mark
         mark = '@'
         if (ga_nodeid().eq.0) write(6,2) mark, istep-1, energy, de,
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
      end
      logical function gsopt_converged(optcyc)
      implicit none
#include "nwc_const.fh"
#include "cgsopt.fh"
c
      integer optcyc
      double precision de
c
c     Return true if we have converged
c
c     Nothing is modified.  Assumes gsopt_compute_info()
c     has been called.
c     gmax_tol,            ! [user] tolerance for max internal gradient
c     grms_tol,            ! [user] tolerance for rms internal gradient
c     xrms_tol,            ! [user] tolerance for rms cartesian step
c     xmax_tol,            ! [user] tolerance for max cartesian step
c
      de = abs(energy-energyp)
      gsopt_converged = 
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
     $      (de .lt. eprec) .and. (optcyc.gt.1))
     $     .or.
     $     ((xmax.lt.xmax_tol).and.(xrms.lt.xrms_tol).and.
     $      (de .lt. eprec) .and. (optcyc.gt.1))
     $      .or. 
c 
c    *** Initial structure already a critical point
c
     $      ((gmax.lt.0.1d0*gmax_tol).and.(grms.lt.0.1d0*grms_tol).and.
     $      (optcyc .eq. 1))
CJMC
      if ((gmax .lt. gmax_tol) .and. (grms .lt. grms_tol).and.
     $      (de .lt. eprec)) then
        write(6,*) "Convergence reached over gradient and energy"
        write(6,*) "Energy change = " , de
        write(6,*) "Energy precision = " , eprec
      end if
c
      end
      subroutine gsopt_cart_pmat(rtdb, geom)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
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
      integer iat,jat
      character*26 date
      integer ind
      logical task_qmmm
      logical gsopt_geom_cart_coords_get
      double precision mass,norm
      double precision ddot
      external ddot
      ind(i,j) = k_pmat + i-1 + (j-1)*ncart
c     
c     FRACTIONAL?
      if (.not. gsopt_geom_cart_coords_get(geom, coords))
     $     call errquit('gsopt_cart_pmat: geom?',geom, GEOM_ERR)
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
cDMR
cDMR  mass-weight projector, if needed
cDMR
      do k = 1, 3               ! x, y, z translations
         call dfill(ncart, 0.0d0, work(1,k), 1)
         if (mswg) then
           norm = sum(atmass(1:nat))
           work(k:3*nat:3,k) = sqrt(atmass(1:nat)/norm)
         else
           call dfill(nat, sqrt(1.0d0/nat), work(k,k), 3)
         endif
      enddo

      do k = 4, 6               ! x, y, z rotations
         do i = 1, nat
            mass = 1d0
            if (mswg) mass = sqrt(atmass(i))
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
            work(i3+1,k) = mass*xx
            work(i3+2,k) = mass*yy
            work(i3+3,k) = mass*zz
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
     $     ('gsopt_cart_pmat: memory for pmat',ncart**2, GEOM_ERR)
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
     $     ('gsopt_cart_bmat: ma?',0, MA_ERR)
c     
      end
      subroutine gsopt_project_hess_grad(hess, pg)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
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
     $     ('gsopt_proj_h_g: memory for pmat',nvar**2, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvar**2, 'pmat',
     $     l_pmat, k_pmat)) call errquit
     $     ('gsopt_proj_h_g: memory for work',nvar**2, MA_ERR)
c
      call geom_hnd_get_data('gsopt.hess',hess, nvar**2)
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
         write(6,*) ' Hessian after projection '
         call output(hess, 1, nvar, 1, nvar, nvar, nvar, 1)
      endif
c
      if (.not. ma_chop_stack(l_work)) call errquit
     $     ('gsopt_p_h_g:ma?',0, MA_ERR)
c
      end
      subroutine gsopt_compute_info()
      implicit none
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
#include "util.fh"
      double precision desphere(max_nvar)
      double precision zeta(max_nvar)
      double precision norm
      double precision ddot
      external ddot
c
c     Compute stuff used for printing and convergence tests
c
c     gmax = maxmimum gradient element in optimization variables
c     grms = rms grad
c     smax = maximum step in opt. var
c     srms = rms step
c     
c     xrms and xmax are computed by gsopt_take_step from the
c     first-order step.
c
      integer i
c
      grms  = 0d0
      srms  = 0d0
      gmax  = 0d0
      smax  = 0d0
c
      call dfill(nvar, 0d0, zeta, 1)
      call dcopy(nvar, radius, 1, zeta, 1)
      call daxpy(nvar, 1.0d00, ds, 1, zeta, 1)
      norm =  ddot(nvar, zeta, 1, zeta, 1)
      norm =  ddot(nvar, g, 1, zeta, 1)/norm
      call dfill(nvar, 0d0, desphere, 1)
      call dcopy(nvar, g, 1, desphere, 1)
      call daxpy(nvar, -norm, zeta, 1, desphere, 1)
c
      do i = 1, nvar
         grms = grms + desphere(i)*desphere(i)
         srms = srms + ds(i)*ds(i)
         gmax  = max(gmax, abs(desphere(i)))
         smax  = max(smax, abs(ds(i)))
      enddo
      grms = sqrt(grms/dble(nvar))
      srms = sqrt(srms/dble(nvar))
c
      end
      subroutine gsopt_hess_cart_guess()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
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
      double precision factor
      integer i,j,iat,jat
      integer l_bi, k_bi, l_hc, k_hc, l_hq, k_hq
c     
      if (.not. ma_push_get(mt_dbl, ncart*nvar, 'binv',
     $     l_bi, k_bi)) call errquit
     $     ('gsopt_hess_cart_guess: ma?', ncart*nvar, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart',
     $     l_hc, k_hc)) call errquit
     $     ('gsopt_hess_cart_guess: ma?', ncart**2, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, max(ncart**2,nvar**2), 'hcart2',
     $     l_hq, k_hq)) call errquit
     $     ('gsopt_hess_cart_guess: ma?', nvar**2, MA_ERR)
c     
      if (ga_nodeid().eq.0) then
         call util_file_name('hess',.false.,.false.,fname)
         open(unit=h_unit,file=fname,form='formatted',status='unknown',
     $        err=99990,access='sequential')
         rewind h_unit
         do i = 1,ncart
            iat = (i+2)/3
            do j = 1,i
               jat = (j+2)/3
               read(h_unit,10000,err=99992,end=99992) x
cDMR mass-weight Hessian
               if (mswg) then
                 factor = sqrt(atmass(iat))*sqrt(atmass(jat))
                 x = x/factor
               endif
cDMR
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
      call geom_hnd_put_data('gsopt.hess',dbl_mb(k_hc), nvar**2)
c
      if (.not. ma_chop_stack(l_bi))
     $     call errquit('gsopt_hess_cart_guess: ma corrupt?',0, MA_ERR)
c
      return
10000 format(f30.15)
99990 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as unknown file'
      call errquit('gsopt_hess_cart: fatal error', 911, GEOM_ERR)
99991 write(6,*)' could not open <',fname(1:inp_strlen(fname)),
     $     '> as new file'
      call errquit('gsopt_hess_cart: fatal error', 911, GEOM_ERR)
99992 write(6,*)' error in reading <',fname(1:inp_strlen(fname)),
     $     '> as hessian file'
      call errquit('gsopt_hess_cart: fatal error', 911, GEOM_ERR)
      end


      subroutine gsopt_compute_actual_step(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
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
      integer i
      logical gsopt_geom_cart_coords_get
c
      if (odebug) then
         write(6,*) ' Expected ds '
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
      call dcopy(nvar, sp, 1, ds, 1) ! Old coordinates into ds()
      if (.not. gsopt_geom_cart_coords_get(geom, sp))
     $    call errquit('gsopt: geom?',0, GEOM_ERR)
c
      do i = 1, nvar
         ds(i) = sp(i) - ds(i)
      enddo
      if (odebug) then
         write(6,*) ' Actual ds '
         call doutput(ds, 1, nvar, 1, 1, nvar, 1, 1)
      endif
c
      end
      logical function gsopt_geom_cart_coords_get(geom, coords)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
      integer geom
      double precision coords(*)
c
c     If we are doing a periodic system and not using internals
c     then we want the fractional coordinates.  Otherwise cartesian.
c
c     If we are including stress append the amatrix
c
      if (.not. geom_cart_coords_get(geom, coords))
     $     call errquit('gsopt: geom cart?',0, GEOM_ERR)
c
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_cart_to_frac(geom, coords))
     $           call errquit('gsopt: frac_to_cart?',0, GEOM_ERR) 
       endif
c
      if (ostress) then
         if (.not. geom_amatrix_get(geom, coords(3*nat_real+1)))
     $        call errquit('gsopt: failed to get amatrix',0,0)
      endif
      if (ostress2) then
         if (.not. geom_lattice_get(geom, coords(3*nat_real+1)))
     $        call errquit('gsopt: failed to get lattice',0,0)
      endif
c
      gsopt_geom_cart_coords_get = .true.
c
      end
      logical function gsopt_geom_cart_coords_set(geom, coords)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
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
     $        call errquit('gsopt: failed to set amatrix',0,0)
      endif
      if (ostress2) then
         if (.not. geom_lattice_set(geom, coords(3*nat_real+1)))
     $        call errquit('gsopt: failed to set lattice',0,0)
      endif

      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_frac_to_cart(geom, coords))
     $           call errquit('gsopt: frac_to_cart?',0,0)
      endif
      if (.not. geom_cart_coords_set(geom, coords))
     $     call errquit('gsopt: geom cart?',0,0)
      if (isystype.ne.0 .and. (.not. zcoord)) then
         if (.not. geom_cart_to_frac(geom, coords))
     $        call errquit('gsopt: frac_to_cart?',0,0)
      endif
c
      gsopt_geom_cart_coords_set = .true.
c
      end


      subroutine gsopt_get_grad(rtdb,geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
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
     $     gx)) call errquit('gsopt: could not get gradient',0,0)
      if (isystype .ne. 0) then
         if (.not. geom_grad_cart_to_frac(geom, gx))
     $        call errquit('gsopt: frac_to_cart?',0,0)
      end if
      if (ostress) then
         if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $   call errquit('gsopt: stress theory not specified',0,RTDB_ERR)
         if (theory.eq.'pspw') then
          if (.not. rtdb_get(rtdb, 'pspw:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('gsopt: could not get stress',0,0)
         else if (theory.eq.'band') then
          if (.not. rtdb_get(rtdb, 'band:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('gsopt: could not get stress',0,0)
         else if (theory.eq.'paw') then
          if (.not. rtdb_get(rtdb, 'paw:stress', mt_dbl, 9,
     $        gx(3*nat_real+1))) call errquit
     $        ('gsopt: could not get stress',0,0)
         else
           call errquit('gsopt: no stress in theory',0,RTDB_ERR)
         end if
      endif

      if (ostress2) then
         if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $   call errquit('gsopt: stress theory not specified',0,RTDB_ERR)
         if (theory.eq.'pspw') then
          if (.not. rtdb_get(rtdb, 'pspw:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('gsopt: could not get stress',0,0)
         else if (theory.eq.'band') then
          if (.not. rtdb_get(rtdb, 'band:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('gsopt: could not get stress',0,0)
         else if (theory.eq.'paw') then
          if (.not. rtdb_get(rtdb, 'paw:lstress', mt_dbl, 6,
     $        gx(3*nat_real+1))) call errquit
     $        ('gsopt: could not get stress',0,0)
         else
           call errquit('gsopt: no stress in theory',0,RTDB_ERR)
         end if
      endif
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gsopt_pickstp(rtdb, geom, istep)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
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
      integer iact, lowest
      double precision lambda, tau, term
      double precision tiniest
      double precision big
      parameter (tiniest = 1.0d-14)
      parameter (big = 1d6)
cjmc
*     integer info
      integer l_hess, k_hess, l_work, k_work, lenwork
      double precision eigval(max_nvar) ! hessian eigenvalues
      double precision pg(max_nvar) ! p.g
      double precision pr(max_nvar) ! p.radius
      double precision gv(max_nvar) ! gradient along eigenvectors
      double precision dv(max_nvar) ! step along eigenvectors
      double precision gc(max_nvar) ! step along eigenvectors
      double precision coords(max_nvar) ! step along eigenvectors
      double precision dsmax    ! max. value of current step (smax is prev.)
c
      double precision beta, s0g0, s0g1, s1g0, s1g1, numerator,
     $     denominator
      double precision bohr, deg ! for printing purposes
      double precision trustds  ! restriction of step in opt. variable
      logical ophigh
      logical gsopt_geom_cart_coords_get
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
     $     ('gsopt_pickstp: memory for hessian',nvar**2, ma_err)
      call gsopt_project_hess_grad(dbl_mb(k_hess), pg)
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
     $     ('gsopt_pickstp: memory for hessian', lenwork, ma_err)
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
c
c     *** calculate constrained saddle gradient ***
c
      call dfill(nvar, 0d0,     pr, 1)
      call dfill(nvar, 0d0, radius, 1)
c
      if (.not. gsopt_geom_cart_coords_get(geom, coords))
     $   call errquit('tropt: geom?',0, GEOM_ERR)
c
      call dcopy(nvar, coords, 1, radius, 1)
      call daxpy(nvar, -1.0d0, center, 1, radius, 1)
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
      if (mswg) call mwcoord(radius, nvar, .true.)
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
      call dgemv('t', nvar, nvar, 1d0, dbl_mb(k_hess), nvar,
     $     radius, 1, 0d0, pr, 1)
c
c     *** calculate combined projected cartesian gradient ***
c
      do iact=1,nvar
        gc(iact) = eigval(iact)*pr(iact) - gv(iact)
      end do
!
!     *** find lowest (countable) eigenvalue position ***
!     ***    assuming a non-shifted hessian matrix    ***
!
      do iact=1,nvar
        if (abs(eigval(iact)).gt.1.0d-10) then
          lowest = iact
          goto 100
        end if
      end do
c
  100 continue
!
!     *** test for the "hard case", i.e. f1 + b1*d1 = 0 ***
!
      if (abs(gc(lowest)).lt.1.d-10) then
!
!     *** set lambda to lowest eigenvalue ***
!
        lambda = -eigval(lowest) + 1.0d-4
!
!     *** calculate reduced constrained step ***
!
        call dfill(nvar, 0d0, dv, 1)
        do iact=1,nvar
          if (iact.ne.lowest) then
            dv(iact) = gv(iact) + lambda*pr(iact)
            dv(iact) = -dv(iact)/(eigval(iact) + lambda)
          end if
        end do
!
!     *** check step length and modify, if necessary ***
!
        term =  ddot(nvar, dv, 1,  dv, 1)
        if (term.le.0.25d0*stride**2) then
          tau = sqrt(0.25d0*stride**2 - term)
          dv(lowest) = tau
          if (oprint) write(6,*) "hard case encountered"
          slength = sqrt(ddot(nvar, dv, 1,  dv, 1))
          goto 200
        end if
      end if
!
      lambda = 0.0d0
!
!     *** restrict step size ***
!
      call gsopt_lambda(gc, eigval, lambda, stride)
      call dfill(nvar, 0d0, dv, 1)
      do iact=1,nvar
        dv(iact) = gv(iact) + lambda*pr(iact)
        dv(iact) = -dv(iact)/(eigval(iact) + lambda)
      end do
      slength = sqrt(ddot(nvar, dv, 1,  dv, 1))
!
!     *** print shift factor(s) ***
!
 200  continue

      if (oprint) then
        if (ga_nodeid().eq.0) write (6,5030) lambda
 5030   format ('lambda for step: ',g10.3)
        if (ga_nodeid().eq.0) write (6,5040) slength
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
c      call gsopt_symmetrize_step(geom)
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
     $     ('gsopt_pickstp: ma?',0, MA_ERR)
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
         call output(ds,1,3,1,nat,3,nat,1)
         write(6,*)
         call util_flush(6)
      endif
c
c
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
      if (mswg) then
        call mwcoord(radius, nvar, .false.)
        call mwcoord(    ds, nvar, .false.)
      end if 
CCCCCCCCCCCCCCCCCCCCCC
CCCC     MASS     CCCC
CCCCCCCCCCCCCCCCCCCCCC
c
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gsopt_raphson(de1,dr,eigval,lambda)
! do newton-raphson step with hessian eigenvalue shift.

      implicit none
#include "nwc_const.fh"
#include "cgsopt.fh"

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
        dr(iact) = -de1(iact)/(eigval(iact) + lambda + tiniest)
      end do
      
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gsopt_lambda(f,eigval,lambda,stride)
! calculate hessian eigenvalue shift factor from 
! trust-region quadratic approximation.

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer maxiter
      parameter (maxiter = 100)

      double precision big,eps,stride
      parameter (big = 1d6, eps = 1d-6)

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
c
c     *** search interval definition ***
c
      do iact=1,nvar
        if (abs(eigval(iact)).gt.1.0d-10) then
          lower = -eigval(iact)
          goto 100
        end if
      end do

  100 continue
c
      itlamb = lower + 1.0d-5
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
          fn = fn + (f(iact)/(eigval(iact) + itlamb))**2
        end do
c
c     *** derivative of fn(lambda) with respect to lambda ***
c
        dfn = 0.0
        do iact=1,nvar
          dfn = dfn + f(iact)**2/(eigval(iact) + itlamb)**3
        end do
!
!     *** calculate next estimate of lambda ***
!
        dx = fn/dfn
        lambda = itlamb + (sqrt(fn)/(0.5d0*stride) - 1.0d0)*dx
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
            itlamb = lambda + 0.05d0
          else
            itlamb = 0.5d0*(upper + lower)
          end if

        end if

        iter = iter + 1
        check = abs(sqrt(fn) - 0.5d0*stride)
        if (check.le.eps) go to 200

      end do
!
!     *** iteration failed ***
!
      error = .true. 
      
      if (oprint.and.ga_nodeid().eq.0) write (6,5000) 
5000  format ('lambda iteration did not converge')

 200  continue
!      
!     *** failed lambda calculation ***
!      
      if (error) then
        write (6,'(a,f16.8)') 'lambda', lambda
        write (6,'(a,f16.8)') 'trust radii', 0.5d0*stride
        call errquit('gsopt_lambda: fatal error', 911, geom_err)
      end if

      end 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccc
      logical function mepgs_freq(rtdb)
c     
c     **** A copy of task_freq ****
c  
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
c
      logical task_hessian
      external task_hessian
c     
      integer rtdb
c     
      integer nat, geom
      logical status
      character*(nw_max_path_len) filehess
c     
      double precision cpu, wall
c 
      logical      ignore
      logical o_reuse
c    
      call ecce_print_module_entry('task frequencies')
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_freq: failed to invalidate status',0,
     &       RTDB_ERR)
      if (ga_nodeid().eq.0 .and.
     $    util_print('task_freq', print_low)) then
        write(LuOut,*)
        write(LuOut,*)
        call util_print_centered(6,
     $      'NWChem Nuclear Hessian and Frequency Analysis',
     $      40,.true.)
        write(LuOut,*)
      endif
*
      if (rtdb_get(rtdb,'vib:reuse',mt_log,1,ignore)) then
        o_reuse = ignore
      else
        o_reuse = .false.
      endif
*
      if (.not.(o_reuse)) then
        status = task_hessian(rtdb)
        if (.not.status) call errquit
     &      ('task_freq: task_hessian failed',911, CALC_ERR)
      else
        if (ga_nodeid().eq.0)
     &    call util_print_centered(LuOut,
     &        'reusing previously generated Hessian',
     &        40,.true.)
        status = .true.
      endif
*
      ignore = rtdb_parallel(.false.)
      if ((ga_nodeid()).eq.0) then
        if (o_reuse) then
          if (rtdb_cget(rtdb,'vib:reuse_hessian_file',1,filehess)) then
            write(LuOut,*)' re-using hessian in file ',
     &          filehess(1:inp_strlen(filehess))
          else
* in case of manual restart and renaming of hess file
            call util_file_name('hess',  .false., .false.,filehess)
          endif
        else
          if (.not. rtdb_cget(rtdb, 'task:hessian file name', 1,
     $        filehess)) call errquit
     $        ('task_freq: failed reading hessian filename from rtdb',0,
     &       RTDB_ERR)
        endif
c
c     create/load reference geometry
c
        if (.not.geom_create(geom,'geometry')) call errquit
     $      ('task_freq:geom_create failed?',1, GEOM_ERR)
        if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $      call errquit
     $      ('task_freq:geom_rtdb_load failed?',2, GEOM_ERR)
        if (.not. geom_ncent(geom,nat)) call errquit
     $      ('task_freq:geom_ncent failed?',3, GEOM_ERR)
        if (.not. geom_destroy(geom)) call errquit
     $      ('task_freq:geom_destroy failed?',911, GEOM_ERR)
        call mepgs_vib(rtdb,filehess,.true.,
     $      0,.false.,0,.false.,nat)
      endif
      call ga_sync()
      ignore = rtdb_parallel(.true.)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - cpu
c
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_freq: failed storing cputime',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_freq: failed storing walltime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_freq: failed to set status',0,
     &       RTDB_ERR)
c
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      mepgs_freq = status
c
c
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
      subroutine mepgs_vib(rtdb,hess_file,in_file,hess_ma,in_ma,
     &    hess_ga,in_ga,natomin)
c
c     **** A copy of vib_vib routine *****
c 
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
#include "errquit.fh"
      LOGICAL PROJEC,ZEROPE,HESOUT,INTERN
*      CHARACTER*7 INPFIL
      INTEGER NATOM, NAT3, NHESS, NHESST
      COMMON /cvib_HESS/ NATOM,NAT3,NHESS,NHESST    ! hessian information
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c:: passed
      integer rtdb             ! [input] rtdb handle
      character*(*) hess_file  ! [input] name of file storing lower triangular packed hessian
      integer hess_ma          ! [input] MA handle to square hessian
      integer hess_ga          ! [input] GA handle to square hessian
      logical in_file          ! [input] hessian is in file get it there
      logical in_ma            ! [input] hessian is in MA array
      logical in_ga            ! [input] hessian is in GA array
      integer natomin          ! [input] number of atoms
c
      logical status
c
      integer i_core, h_core, iii, ioldlabs, ivc, itot
      integer nels, npri, ihess, icoord, ihesst, ihesstcp,
     &    ihessp, iegval, iegvec, iddpol, iddpolq, intense
      integer imass, iscr, i10, i20, i30, i40,i_w1,l_w1,i_w2,l_w2
      double precision dbl_tmp
      logical first_pass
      character*255 dipole_file
      logical dipole_file_exists
      logical animation_on
c
      first_pass = .true.
*... check input logic 
      status =           (in_file.and.(in_ma.or.in_ga))
      status = status.or.(in_ma.and.(in_file.or.in_ga))
      status = status.or.(in_ga.and.(in_file.or.in_ma))
      if (status) then
        write(luout,*)' ERROR: more than one source for hessian '
        write(luout,*)' in_file :',in_file
        write(luout,*)' in_ma   :',in_ma
        write(luout,*)' in_ga   :',in_ga
        call errquit(' vib_vib: error ',911, UNKNOWN_ERR)
      endif
      if (in_ga)
     &    call errquit
     &    ('vib_vib: ga access to hessian not implemented yet',911,
     &       CAPMIS_ERR)
C
C Zero core
C
      call vib_setup ! subroutine to set up some constants
      NATOM  =  natomin ! number of atoms in species.
      IF (NATOM.LE.1) THEN      ! check for incorrect number of atoms
          WRITE(6,*)' You want to calculate the vibrational ',
     +              'frequencies for ',NATOM,' atoms?'
          WRITE(6,*)' Unfortunately this is not possible '
          CALL errquit('vib_vib: bomb',911, INPUT_ERR)
      ENDIF
      NAT3   =  NATOM*3         ! 3-N (as in degrees of freedom)
      NHESS  =  NAT3*NAT3       ! dimension of hessian
      NHESST =  NAT3*(NAT3+1)/2 ! dimension of lower triangular hessian
      NELS   =  7*MAX(3*NATOM-6,1)
      NPRI = 0
C
C Calculate pointers
C
      IHESS    =  1                ! square hessian
      IHESST   =  IHESS  + NHESS   ! lower-tri Hessian
      ihesstcp =  IHESST + NHESST  ! copy of lower-tri hessian
      ICOORD   =  IHESSTcp + NHESST  ! geometrical coordinates
      IMASS    =  ICOORD + NAT3    ! mass of each atom
      IEGVAL   =  IMASS  + NATOM   ! eigenvalues from Hessian matrix
      IEGVEC   =  IEGVAL + NAT3    ! eigenvectors from Hessian matrix
      ISCR     =  IEGVEC + NHESS   ! dynamic bottom of core array
      IHESSP   =  ISCR   + 8*NAT3  ! addition of scratch space needed
C--------The following are pointers for GAMESS internal coordinate subroutines.
C
      I10      =  IHESSP + NAT3*NAT3 ! space for zmat
      I20      =  I10    + NELS
      I30      =  I20    + NAT3*NAT3 ! Space to represent internal coord. Hessian
      I40      =  I30    + NAT3*NAT3 !
      Iddpol   =  I40    + 8*NAT3    ! derivative dipole in cartesians
      Iddpolq  =  Iddpol + 3*NAT3    ! derivative dipole in normal modes
      Intense  =  Iddpolq+ 3*NAT3    ! intensities
      ITOT     =  Intense + 3*natom*4
      itot = itot + 2*natom+1 + 6*nat3 ! extra for call to rdinp
c
      if (.not.ma_push_get
     &    (MT_DBL,itot,' core for vib ',h_core, i_core))
     &    call errquit('vib_vib: ma_push_get failed ',911, MA_ERR)
C
C Reset pointers for MA array
C
      IHESS    =  i_core           ! square hessian
      IHESST   =  IHESS  + NHESS   ! lower-tri Hessian
      ihesstcp =  IHESST + NHESST  ! copy lower-tri Hessian
      ICOORD   =  IHESSTcp + NHESST  ! geometrical coordinates
      IMASS    =  ICOORD + NAT3    ! mass of each atom
      IEGVAL   =  IMASS  + NATOM   ! eigenvalues from Hessian matrix
      IEGVEC   =  IEGVAL + NAT3    ! eigenvectors from Hessian matrix
      ISCR     =  IEGVEC + NHESS   ! dynamic bottom of core array
      IHESSP   =  ISCR   + 8*NAT3  ! addition of scratch space needed
C--------The following are pointers for GAMESS internal coordinate subroutines.
C
      I10      =  IHESSP + NAT3*NAT3 ! space for zmat
      I20      =  I10    + NELS
      I30      =  I20    + NAT3*NAT3 ! Space to represent internal coord. Hessian
      I40      =  I30    + NAT3*NAT3 !
      ioldlabs =  I40    + 8*NAT3    !
      ivc      =  ioldlabs + 2*natom + 1
      iddpol   =  ivc + 6*nat3
      iddpolq  =  iddpol + 3*nat3
      intense  =  iddpolq + 3*nat3
      Itot     =  intense + 3*natom*4
c
c read/load hessian and form triangle/square as needed
c
      if (in_ma) then
        ihess = hess_ma   ! simply reset ptr to dbl_mb
*        form triangle
        call vib_dtrngl(dbl_mb(ihess),dbl_mb(ihesst),nat3,nat3)
      endif
      if (in_file) then
        open(unit=69,file=hess_file,form='formatted',status='old',
     &      err=99900,access='sequential')
        do iii = 0,(nhesst-1)
          read(69,*,err=99901,end=99902)dbl_tmp
          dbl_mb(ihesst+iii) = dbl_tmp
        enddo
        close(unit=69,status='keep')
        call vib_dsquar(dbl_mb(ihesst),dbl_mb(ihess),nat3,nat3)
      endif
      call util_file_name('fd_ddipole',.false., .false.,dipole_file)
      dipole_file_exists = .false.
      inquire(file=dipole_file,exist=dipole_file_exists)
      if (dipole_file_exists) then
        open(unit=70,file=dipole_file,form='formatted',status='old',
     &      err=89900,access='sequential')
        do iii = 0,((3*nat3)-1)
          read(70,*,err=89901,end=89902) dbl_tmp
          dbl_mb(iddpol+iii) = dbl_tmp
        enddo
        close(unit=70,status='keep')
      endif
00001 continue
      write(luout,*)
      write(luout,*)
      if (.not. first_pass) then
        WRITE(luout,*)
     &      '       Vibrational analysis via the FX method '
        write(luout,*)
     &      ' --- with translations and rotations projected out ---'
        write(luout,*)
     &      ' --- via the Eckart algorithm                      ---'
      endif
      if (first_pass) then
c
c save a copy of hesst
c
        call dcopy(nhesst,dbl_mb(ihesst),1,dbl_mb(ihesstcp),1)
      else
c
c restore copy of hesst and hess
c
        call dcopy(nhesst,dbl_mb(ihesstcp),1,dbl_mb(ihesst),1)
        call vib_dsquar(dbl_mb(ihesst),dbl_mb(ihess),nat3,nat3)
      endif
C
C Read in user input and tape10 arrays.  NO INPUT REQUIRED NOW
C      Note: ! scratch pointer for atom charges (real NAT words)
C              and atom lables (real 2*nat words)
      call vib_rdinp(
     &    dbl_mb(ihess),dbl_mb(ihesst),dbl_mb(icoord),
     &    dbl_mb(imass),dbl_mb(iscr),  dbl_mb(ioldlabs),
     &    dbl_mb(i10),nels,projec,zerope,hesout,intern,
     &    rtdb,first_pass)
      if (projec) then
      if (.not.ma_push_get
     &    (MT_DBL,nat3*nat3,' w1 ',l_w1, i_w1))
     &    call errquit('vib_vib: ma_push_get failed ',911, MA_ERR)
      if (.not.ma_push_get
     &    (MT_DBL,nat3*nat3,' w2 ',l_w2, i_w2))
     &    call errquit('vib_vib: ma_push_get failed ',911, MA_ERR)
        call vib_eckart( dbl_mb(ihess), dbl_mb(ihessp), dbl_mb(ihesst),
     &      dbl_mb(icoord),  dbl_mb(ivc), dbl_mb(i_w1),dbl_mb(i_w2))
        if(.not.ma_chop_stack(l_w1)) call errquit(
     '       ' vib_vib: machopstack failed',1, MA_ERR)
      end if
* rak dfill
      CALL Dfill(NAT3,0.0d00,DBL_MB(ISCR),1) ! zero scratch used
C
      CALL mepgs_vib_hmass(DBL_MB(IHESST),DBL_MB(IMASS)) ! mass weight and scale hessian
C
C Diagonalize mass-weighted, scaled hessian matrix
C     Note: ! scratch pointer for givens (real 5*NAT3 words)
c use hessp as scratch now calling rsg
C
      call dscal(nhesst,1.0d3,dbl_mb(ihesst),1)
      CALL vib_CALLG(DBL_MB(IHESSt),nhesst,DBL_MB(IHESSP),
     &    dbl_mb(iscr),dbl_mb(iscr+nat3),DBL_MB(IEGVAL),
     &    DBL_MB(IEGVEC), NAT3,NAT3)
      call dscal(nat3,1.0d-3,dbl_mb(iegval),1)
C
C
C
      CALL vib_NMASS(DBL_MB(IEGVEC),DBL_MB(IMASS)) ! "unmass" weight the normal modes.
C
C *** Note: DBL_MB(IHESST) now destroyed if needed reinitialize from DBL_MB(IHESS)
C
* rak dfill
      call dfill(5*nat3,0.0d00,dbl_mb(iscr),1)    ! zero scratch used
C
CJMC      CALL vib_WRTFREQ(rtdb,DBL_MB(IEGVAL),NAT3,ZEROPE,NPRI) ! Write out the zero-point energy
C
      CALL vib_CLEAN(DBL_MB(IEGVEC),NAT3*NAT3,1.0D-27) ! CLEAN eigenvectors
      if (.not. first_pass .and .dipole_file_exists) then
        call vib_intense(rtdb,dbl_mb(iegvec),dbl_mb(iegval),natom,
     &      dbl_mb(iddpol),dbl_mb(iddpolq),dbl_mb(intense),
     &      first_pass)
      endif
c
      if (.not.first_pass) then
* if any negative eigenvalues print out steps in their direction
        call mepgs_vib_istep(
     &        rtdb,nat3,natom,
     &        dbl_mb(iegvec),dbl_mb(iegval),
     &        dbl_mb(icoord),dbl_mb(iscr),
     &        dbl_mb(iscr+nat3),dbl_mb(iscr+(2*nat3)))

      endif
c
      if (first_pass) then
        first_pass = .false.
        goto 00001
      endif
C
      if (.not.ma_pop_stack(h_core)) call errquit
     &    ('vib_rdinp ma_pop failed',911, MA_ERR)
      return
89900 continue
      write(luout,*)'dipole_file => ',dipole_file
      call errquit('vib_vib: error opening file: "dipole_file"',811,
     &       DISK_ERR)
89901 continue
      write(luout,*)'dipole_file => ',dipole_file
      call errquit('vib_vib: error reading file: "dipole_file"',811,
     &       DISK_ERR)
89902 continue
      write(luout,*)'dipole_file => ',dipole_file
      call errquit
     & ('vib_vib: unexpected EOF when reading file: "dipole_file"',811,
     &       DISK_ERR)
99900 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit('vib_vib: error opening file: "hess_file"',911,
     &       DISK_ERR)
99901 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit('vib_vib: error reading file: "hess_file"',911,
     &       DISK_ERR)
99902 continue
      write(luout,*)'hess_file => ',hess_file
      call errquit
     & ('vib_vib: unexpected EOF when reading file: "hess_file"',911,
     &       DISK_ERR)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      subroutine mepgs_vib_istep(rtdb,nat3,natom,
     &    eigenvecs,eigenvals,coords,steps,rawstep,master)
c
c     **** Copy of vib_istep
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
      double precision ddot
      external ddot
*
      integer rtdb    ! [input] rtdb handle
      integer natom   ! [input] number of atoms
      integer nat3    ! [input] 3*number of atoms
      double precision eigenvecs(nat3,nat3) ! [input](xyz&atom,mode)
      double precision eigenvals(nat3)      ! [input] (mode)
      double precision master(3,natom)    ! [scratch] original coordintates
      double precision coords(3,natom)    ! [scratch] coords after step
      double precision rawstep(3,natom)  ! [scratch] step generated by vector 
      double precision steps(3,natom)  ! [scratch] step generated by vector 
c
      integer imode,ivec,iatom,ixyz
      integer geom
      double precision scale
      double precision xyz(3),charge
      double precision length_of_step
      character*16 tag
      character*10 units
      intrinsic sqrt
c
CJMC
      double precision delta
      double precision atmass(natom)
CJMC
      double precision thresh
      parameter (thresh=1.0d-2)
c::-statement function
      logical is_it_close_to  
      double precision value,test
      intrinsic abs
*---          is value close to test?
      is_it_close_to(value,test) = (abs(value-test).lt.thresh)
c
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('vib_istep: geom create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('vib_istep: geom_rtdb_load failed',911, RTDB_ERR)
      if (.not.geom_cart_coords_get(geom,master)) call errquit
     &    ('vib_istep: geom_get_cart_coords failed',911, GEOM_ERR)
      if (.not.geom_get_user_scale(geom,scale)) call errquit
     &    ('vib_istep: geom_get_user_scale failed',911, GEOM_ERR)
      if (.not.geom_get_user_units(geom,units)) call errquit
     &    ('vib_istep: geom_get_user_units failed',911, GEOM_ERR)
c
      imode = 1
      ivec = 0
      write(luout,10000)imode,eigenvals(imode)
      call dfill(nat3,0.0d00,rawstep,1)
      do iatom = 1, natom
        do ixyz = 1,3
          ivec = ivec+1
          rawstep(ixyz,iatom) = eigenvecs(ivec,imode)
        enddo
      enddo
CJMC
      delta = sqrt(2.0d0*evib/abs(eigenvals(imode)))
CJMC 
      call dcopy(nat3,rawstep,1,steps,1)
      call dscal(nat3, delta, steps, 1)
      length_of_step = sqrt(ddot(nat3,steps,1,steps,1))
      write(luout,10001)length_of_step,units
      do iatom=1,natom
        if(.not.geom_cent_get(geom,iatom,tag,xyz,charge))
     &    call errquit('vib_istep: geom_cent_get failed',911, GEOM_ERR)
        write(luout,10002)iatom,tag,charge,
     &      (steps(ixyz,iatom),ixyz=1,3)
      enddo
      write(luout,10003)
c
c      **** Store "forward" direction ****
c
      call dcopy(nat3,master,1,coords,1)
      call daxpy(nat3,1.0d00,steps,1,coords,1)
      if (.not.geom_cart_coords_set(geom,coords)) call errquit
     &   ('vib_istep: geom_cart_coords_set failed',911, GEOM_ERR)
      if (.not. geom_rtdb_store(rtdb, geom, 'ircforward'))
     $     call errquit('mepgs_vib_istep: grs?',geom, RTDB_ERR)
c
c      **** Store "backward" direction ****
c
      call dcopy(nat3,master,1,coords,1)
      call daxpy(nat3,-1.0d00,steps,1,coords,1)
      if (.not.geom_cart_coords_set(geom,coords)) call errquit
     &   ('vib_istep: geom_cart_coords_set failed',911, GEOM_ERR)
      if (.not. geom_rtdb_store(rtdb, geom, 'ircbackward'))
     $    call errquit('mepgs_vib_istep: grs?',geom, RTDB_ERR)
c
c     **** Deallocate geom ****
c
      if (.not.geom_destroy(geom)) call errquit
     &    ('vib_istep: geom_destroy failed',911, GEOM_ERR)
c
c
10000 format(/,/,/,1x,78('='),/,6x,'Negative Nuclear Hessian Mode',
     &      i5,2x,'Eigenvalue = ',f9.2,' a.u.',/,1x,78('-'))
10001 format(2x,' Raw step length:',f7.3,1x,a10,';',
     &    2x, 'The Raw step for this mode is:')
10002 format(' ',i4,' ',a16,' ',f10.4,3f15.8)
10003 format(78('-'))
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mepgs_vib_HMASS(HESST,ATMASS)
* $Id$
C
C  This routine mass weights and scales the Hessian matrix
C       The scaling is done to avoid numerical problems in the
C       diagonalization routine
C
      IMPLICIT NONE ! REAL*8 (A-H,O-Z)
#include "cmepgs.fh"
      INTEGER NAT, NAT3, NHESS, NHESST
      COMMON /CVIB_HESS/ NAT,NAT3,NHESS,NHESST   ! HESSIAN INFORMATION
c
      double precision HESST(NHESST) ! lower triangular Hessian
      double precision ATMASS(NAT) ! mass of the atoms
c
      double precision fact, scale
      integer ii, jj, jjend, iatii, iatjj, idum
      double precision mass_ii, mass_jj
C
C      set up function for locating i,j elements packed canonically as ij
C
      integer i, j, isym2, iatom
      ISYM2(I,J)=MAX(I,J)*((MAX(I,J))-1)/2 + MIN(I,J)
      IATOM(I)  = (I+2)/3   ! function for coordinate I is on atom IATOM
C
      DO 00100 II = 1,NAT3 ! loop over coordinates
        JJEND = II
        IATII = IATOM(II) ! coordinate II is for atom IATII
        DO 00200 JJ = 1,JJEND ! loop over coordinates
          IDUM = ISYM2(II,JJ) ! get canonical index
          IATJJ = IATOM(JJ) ! coordinate JJ is for atom IATJJ
          mass_ii = atmass(iatii)
          mass_jj = atmass(iatjj)
          if (abs(mass_ii).lt.1.0d-01) mass_ii  = 1.0d05
          if (abs(mass_jj).lt.1.0d-01) mass_jj  = 1.0d05
          FACT = SQRT(mass_ii)*SQRT(mass_jj) ! mass weight
          HESST(IDUM) = HESST(IDUM)/FACT ! weight Hessian
00200   CONTINUE
00100 CONTINUE

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mwcoord(vector, ndeg, put_mass)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
      integer ndeg
      double precision vector(ndeg)
      logical put_mass
c
      integer ipos, iatom, ixyz
c
      if (put_mass) then
        ipos = 0
        do iatom=1, nat
          do ixyz=1,3
            ipos = ipos + 1
            vector(ipos) = vector(ipos)*sqrt(atmass(iatom))
          end do
        end do
      else if (.not. put_mass) then
        ipos = 0
        do iatom=1, nat
          do ixyz=1,3
            ipos = ipos + 1
            vector(ipos) = vector(ipos)/sqrt(atmass(iatom))
          end do
        end do
      end if
c
      end
      subroutine mwgrad(vector, ndeg, put_mass)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cmepgs.fh"
#include "cgsopt.fh"
      integer ndeg
      double precision vector(ndeg)
      logical put_mass
c
      integer ipos, iatom, ixyz
c
      if (put_mass) then
        ipos = 0
        do iatom=1, nat
          do ixyz=1,3
            ipos = ipos + 1
            vector(ipos) = vector(ipos)/sqrt(atmass(iatom))
          end do
        end do
      else if (.not. put_mass) then
        ipos = 0
        do iatom=1, nat
          do ixyz=1,3
            ipos = ipos + 1
            vector(ipos) = vector(ipos)*sqrt(atmass(iatom))
          end do
        end do
      end if
c
      end

      subroutine updstride(geom)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "geom.fh"
      integer geom, l_hess, k_hess
      double precision gamma(nvar),delta(nvar),newx(max_nvar)
      double precision epredict,echange,ratio
      logical,external :: gsopt_geom_cart_coords_get
      if (.not. ma_push_get(mt_dbl, nvar**2, 'hess',
     $     l_hess, k_hess)) call errquit
     $     ('gsopt_pickstp: memory for hessian',nvar**2, ma_err)
      call geom_hnd_get_data('gsopt.hess',dbl_mb(k_hess), nvar**2)
      if (.not. gsopt_geom_cart_coords_get(geom, newx))
     $    call errquit('mepgs: geom?',0, geom_err) 

      if (mswg) then
        call mwcoord(newx, nvar, .true.)
        call mwcoord(oldgeo, nvar, .true.)
        call mwgrad(g, nvar, .true.)
      endif

      delta(1:nvar) = newx(1:nvar) - oldgeo(1:nvar)
      call dgemv('n',nvar,nvar,1d0,dbl_mb(k_hess),nvar,delta,1,0d0,
     $            gamma,1)
      epredict = -dot_product(g(1:nvar),delta(1:nvar))-
     $            0.5d0*dot_product(delta(1:nvar),gamma)
      echange = energyref - energy
      ratio = echange/epredict


      if (ratio.gt.0.75d0) then
        stride = 2d0*stride
      elseif(ratio.lt.0.25d0) then
        stride = stride/2d0
      endif

      stride = min(stride,0.3d0)

      if (mswg) then
        call mwcoord(newx, nvar, .false.)
        call mwcoord(oldgeo, nvar, .false.)
        call mwgrad(g, nvar, .false.)
      endif

      if (.not.ma_chop_stack(l_hess)) then
        call errquit('updstride: could not chop stack',0,MA_ERR)
      endif

      end
