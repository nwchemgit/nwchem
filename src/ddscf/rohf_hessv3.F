c ... jochen: this comes from rohf_hessv2.F but has everything
c     "v2" replaced by "v3", for use with frequency-dependent response
c     and is called from the cphf_solve3 part of the cphf code
c
c ... jochen: further mods were made to accomodate the situation that
c     we might have damping in the response. That causes all quantities to
c     have an imaginary part, too

      subroutine rohf_hessv3(acc, g_x, g_ax, g_x_im, g_Ax_im,
     &   omega, limag, lifetime, gamwidth, ncomp)
      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "stdio.fh"
#include "util.fh"
#include "global.fh"
c     
c     $Id$
c
c ... jochen: these two arrays now have two components:
      integer g_x(2)  ! [input]  A-matrix elements for density matrix
      integer g_ax(2) ! [output] Perturbed Fock operator
c ... jochen: also, we might have imaginary components:
      integer g_x_im(2)  ! [input]  A-matrix elements, Im
      integer g_ax_im(2) ! [output] Perturbed Fock operator, Im

      double precision acc, omega, gamwidth
      logical limag, lifetime
c     
      integer gtype,grow,gcol,growp,gcolp, ipm, ncomp
      logical oprint, debug
c
c     ================================================================
c
      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

      if (debug) write (6,*)
     &   'rohf_hessv3: limag, omega, lifetime, gamwidth',
     &   limag, omega, lifetime, gamwidth
c
c     Check for debug 
c     
      oprint= util_print('rohf_hessv2',print_debug)
      if (crohf_init_flag.ne.1)
     $     call errquit('rohf_hessv3: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
c
c ... jochen: use first component for the dimension checks.
c     the second component MUST have the same dimensions
c     otherwise there will be problems
      call ga_inquire(g_x(1),gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv3: invalid vector length',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_ax(1),gtype,growp,gcolp)
      if (growp.ne.crohf_vlen)
     $     call errquit('rohf_hessv3: invalid vector length',0,
     &       UNKNOWN_ERR)
      if (gcol.ne.gcolp)
     $     call errquit('rohf_hessv3: invalid no. of vectors',0,
     &       UNKNOWN_ERR)
c     
c     Call internal routine
c  
      if (debug) write (6,*) 'calling rohf_hessv_xx3'
      call rohf_hessv_xx3( basis, geom, nbf, nmo,
     $     nclosed, nopen,
     $     pflg, g_movecs, oskel, noskew,
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $     acc, lshift, g_x, g_ax, g_x_im, g_Ax_im, omega, limag,
     &   lifetime, gamwidth, ncomp)
      if (debug) write (6,*) 'back from rohf_hessv_xx3'
c   DIM/QM JEM
c      if(ldimqm .and. lrsp .and. dimqm_on) then
c        call dimqm_addop(g_x, g_x_im, g_ax, g_ax_im, ncomp, 
c     $                   limag, lifetime)
c      end if
c
c     Zap numbers much smaller than acc to ensure hard zeroes 
c     remain unpolluted ... cannot use a threshold larger than the
c     integral accuracy since can break symmetry in non-abelian groups
c     Also must ensure that the threshold tends to zero to permit
c     tight convergence.
c
c ... jochen: screen components
      do ipm = 1,ncomp
        call ga_screen(g_ax(ipm),
     &       max(min(acc*acc,acc*0.01d0,1d-12),1d-16))
        if (lifetime) call ga_screen(g_ax_im(ipm), 
     &       max(min(acc*acc,acc*0.01d0,1d-12),1d-16))
      enddo
c
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rohf_hessv_xx3( basis, geom, nbf, nmo, nclosed, nopen, 
     $     pflg,
     $     g_movecs, oskel, noskew, g_fcv, g_fpv, g_fcp,
     $     acc, lshift, g_x, g_ax, g_x_im, g_Ax_im, omega, limag,
     &     lifetime, gamwidth, ncomp)
C     $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
c     
      integer basis, geom
      integer nbf, nmo, nclosed, nopen
      integer pflg
      integer g_movecs
      logical oskel, noskew
      integer g_fcv, g_fpv, g_fcp
      double precision acc
      double precision lshift
c ... jochen: input arrays g_x and g_Ax have two components here
      integer g_x(2), g_ax(2), vlen, nvec, g_tmp, gtype, ipm
      integer g_x_im(2)  ! [input]  A-matrix elements, Im
      integer g_ax_im(2) ! [output] Perturbed Fock operator, Im

      double precision omega, gamwidth, wls, wlsim
      logical limag, lifetime
      integer ncomp
      double precision omg(ncomp),gam(ncomp)
      logical debug
c     DIM/QM JEM
      logical ldimqm
c
c     =================================================================
c
      debug = (.false. .and. ga_nodeid().eq.0) ! for code development
c
      if (debug) write (6,*) 'hessv3: omega =',omega
      if (debug) write (6,*) 'hessv3: limag =',limag
      if (debug) write (6,*)
     &   'hessv3: lifetime, gamwidth =',lifetime, gamwidth
c
      do ipm = 1,ncomp
        call ga_zero(g_Ax(ipm))
        if (lifetime) call ga_zero(g_Ax_im(ipm))
      end do
      if (pflg.gt.2 .or. pflg.le.0) then
         call errquit('rohf_hessv_xx: pflg invalid ', pflg,
     &       UNKNOWN_ERR)
      endif
c
c ... jochen: to be consistent with the preconditioner, where
c     the level shift is added, we need to do the same thing here
c     and also add and subtract the frequency times 4 (it is times
c     4 because of the factors of 4 in rohf_hessv_1e and in the
c     preconditioner)
c     During a response calculation, pflg is equal to 2
c
c     what do we do here? Compare Gauss' paper Eqs. (32) and (135):
c     The lhs of the CPHF equations contain a term
c     (e_a - e_i -/+ omega) U_ai. First, we initialize g_Ax with
c     the term proportional to omega, then we add the delta-e term
c     (the e's are the orbital energies, calculated in hessv_1e as
c     the diagonal of the Fock matrix transformed to the MO basis)
c
      if (pflg .gt. 0) then
c
         omg(1)=-omega
         omg(2)= omega
         gam(1)=-gamwidth
         gam(2)= gamwidth
        if (.not.lifetime) then
c         no damping: initialize Ax with terms proportional omega
         do ipm=1,ncomp
          wls   = lshift + 4d0 * omg(ipm)
c
c         HvD: g_x_im does not exist if (.not.lifetime) (see
c         cphf_solve3.F line 406) so the code on the next line fails. I
c         am assuming that a non-existing imaginary part equates to it
c         being zero. Hence the copy and scale operation on the
c         subsequent two lines should do what is required.
c         call ga_dadd(wls,g_x(ipm),wlsim,g_x_im(ipm),g_ax(ipm))
          call ga_copy(g_x(ipm),g_ax(ipm))
          call ga_scale(g_ax(ipm),wls)
         enddo
        else                    ! lifetime
c         take care of damping here: Re and Im are coupled by gamwidth
         do ipm=1,ncomp
          wls   = lshift + 4d0 * omg(ipm)
          wlsim =  -4d0 * gam(ipm)
          call ga_dadd(wls,g_x(ipm),wlsim,g_x_im(ipm),g_ax(ipm))
          wls   =   4d0 * gam(ipm)
          wlsim = lshift + 4d0 * omg(ipm)
          call ga_dadd(wls,g_x(ipm),wlsim,g_x_im(ipm),g_ax_im(ipm))
         enddo
        endif                   ! .not.lifetime

        call ga_sync()
        if (debug) write (6,*) 'calling rohf_hessv_1e'
c 
c       next: add (e_a - e_i) times A (also called U) matrix to Ax
        call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $     g_fcv, g_fpv, g_fcp,
     $     g_x(1), g_ax(1))
        if (ncomp.gt.1)
     &     call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $     g_fcv, g_fpv, g_fcp,
     $     g_x(2), g_ax(2))
        if (lifetime) then
          call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $       g_fcv, g_fpv, g_fcp,
     $       g_x_im(1), g_ax_im(1))
          if (ncomp.gt.1)
     &       call rohf_hessv_1e( basis, geom, nmo, nclosed, nopen,
     $       g_fcv, g_fpv, g_fcp,
     $       g_x_im(2), g_ax_im(2))
        endif                   ! lifetime
        if (debug) write (6,*) 'back from rohf_hessv_1e'

      endif                     ! pflg.gt.0

      if (pflg .gt. 1) then
c
c       the next call basically uses the current guess for the solution
c       vector x (in g_x, which is the perturbed density matrix in the
c       MO basis) and calculates the perturbed Fock operator in the MO basis.
c       real and imaginary part of that Fock operator can be handled
c       separately here
c
        if (ncomp.gt.1) then    ! call 2e code for dynamic case
          if (debug) write (6,*) 'calling rohf_hessv_2e3'
          call rohf_hessv_2e3( basis, geom, nbf, nmo, nclosed, nopen,
     $       g_movecs, oskel, noskew, g_x, acc, ! was min(1d-6,acc)
     $       g_ax, limag)
          if (debug) write (6,*) 'back from rohf_hessv_2e3'
          if (lifetime) then
            if (debug) write (6,*) 'calling rohf_hessv_2e3 for Im part'
            call rohf_hessv_2e3( basis, geom, nbf, nmo, nclosed, nopen,
     $         g_movecs, oskel, noskew, g_x_im, acc, ! was min(1d-6,acc)
     $         g_ax_im, limag)
            if (debug) write (6,*) 'back from rohf_hessv_2e3 Im part'
          endif                 ! lifetime
        else                    ! call static 2e code
          call rohf_hessv_2e2( basis, geom, nbf, nmo, nclosed, nopen,
     $       g_movecs, oskel, noskew, g_x(1), acc, ! was min(1d-6,acc)
     $       g_ax(1))
          if (lifetime) then
            call rohf_hessv_2e2( basis, geom, nbf, nmo, nclosed, nopen,
     $         g_movecs, oskel, noskew, g_x_im(1), acc, ! was min(1d-6,acc)
     $         g_ax_im(1))
          endif ! lifetime
        endif                   ! ncomp
        
      endif                     ! pflg.gt.1
c   DIM/QM JEM
c   This call is placed differently than in rohf_hessv2 because we
c   require knowledge of the real and imaginary simultaneously.
      if (.not.rtdb_get(bgj_get_rtdb_handle(),'dimqm:lrsp',MT_LOG, 1,
     $                  ldimqm)) ldimqm = .false.
      if(ldimqm) then
        call dimqm_addop(g_x, g_x_im, g_ax, g_ax_im, ncomp, 
     $                   limag, lifetime)
      end if
c     
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rohf_hessv_2e3( basis, geom, nbf, nmo, nclosed, nopen,
     $     g_movec, oskel, noskew, g_x, acc, g_ax, limag)
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "cscfps.fh"
#include "rtdb.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "case.fh"
c     
c     Return the ROHF orbital 2e-Hessian vector product, g_ax = A * g_x
c
c ... jochen: modified version of rohf_hessv_2e2 which keeps track
c     of two sets of input vectors that couple via the density matrix.
c     one could likely save some memory here by re-using temp arrays
c ... jochen: Also made modifications to calculate imaginary terms due
c     to finite lifetime damping
c     
c     
ccccccccccccccc This code does NOT work for open shell!!!!!ccccccccccccccccc
c
c
      integer basis, geom       ! basis & geom handle
      integer nbf, nclosed, nopen ! Basis size and occupation
      integer nmo               ! No. of linearly dependent MOs
      integer g_movec           ! MO coefficients
      logical oskel
      integer g_x(2)            ! Argument
      double precision acc      ! Accuracy of "Fock" construction
      integer g_Ax(2)           ! Hessian product
      logical limag             ! imaginary perturbation?
c     
      integer nvir, voff, xoff, xend
      integer ooff, oend
c     
      integer nfock, ifock
      integer g_tmp1, g_tmp2, g_tmp3
c
      integer ivec, nvec, gtype, vlen
      integer k_jfac(2),l_jfac(2),k_kfac(2),l_kfac(2)
c
      integer g_dcv,  g_dcp,  g_dpv, g_dens(2), g_fock(2)
      double precision tol2e
      logical odebug
      logical oprint
      logical noskew
c
      integer dims(3), chunk(3), alo(3), ahi(3), blo(2), bhi(2)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision one, zero, mone, four, half, mhalf, two, mtwo
      double precision itol_floor, itol_ceil
      double precision xc_xfac1
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)
c
c ... jochen: additional declarations added:
      integer ipm ! counter for density matrix components
      character*(255) cstemp
      integer g_pmats(2), g_pmata(2), g_h1mat(2)
      double precision pre_factor, tenm6
      parameter (tenm6 = 1d-6)
      logical debug
c
      integer g_work
c
c     =================================================================
c
      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

      if (debug) write (6,*) 'hessv3: limag',limag
c
      if (debug) call ga_print(g_Ax(1))
      if (debug) call ga_print(g_Ax(2))
c
      oprint= util_print('rohf_hessv2',print_debug)
      xc_xfac1=bgj_kfac() ! amount of exact exchange in the functional
      if (debug) write (6,*) 'xc_xfac1: ', xc_xfac1
c     
c     This get's cleaned up a lot when the Fock build accepts
c     a single multiple dimension GA for its input.
c
c     Lots more to do on memory management.
c
      if (nopen.ne.0) call errquit
     $     ('rohf_h2e3: does not work for open shells',nopen,
     &       UNKNOWN_ERR)
c
c ... jochen: do the inquire only on the first component of g_x.
c     if the second component is not the same we are in trouble
      call ga_inquire(g_x(1),gtype,vlen,nvec)
      odebug = util_print('rohf_hessv', print_debug)
      if (debug) write (6,*) 'rohf_h2e3 step 1'
c     
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1
      ooff = nclosed + 1
      oend = nclosed + nopen
c
c     allocate Mas for jfacs
c
      do ipm = 1,2
        if (.not.MA_Push_Get(mt_dbl,nvec*2,'hessv jfacs',
     &     l_jfac(ipm),k_jfac(ipm)))
     &     call errquit('hessv2: cannot allocate jfac',nvec*2,
     &     MA_ERR)
        if (.not.MA_Push_Get(mt_dbl,nvec*2,'hessv kfacs',
     &     l_kfac(ipm),k_kfac(ipm)))
     &     call errquit('hessv2: cannot allocate kfac',nvec*2,
     &     MA_ERR)
      end do
      if (debug) write (6,*) 'rohf_h2e3 step 2'

c     Allocate global for work space
c     
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp2')
      if (nopen .gt. 0) then
         g_tmp3 = ga_create_atom_blocked(geom, basis, 'rohf_h2e: tmp3')
      endif
      if (debug) write (6,*) 'rohf_h2e3 step 3'
c
c     Allocate arrays and generate the AO densities
c ... jochen: the _first_ dimension of the density matrix is the
c     number of components for the perturbation (usually 3)
c
      dims(1) = nbf
      dims(2) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
c
c ... jochen: create pmats and pmata arrays as well as g_h1mat
      do ipm = 1,2
        write(cstemp,'(a,i1)') 'pmats_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmats(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmats(ipm))
        write(cstemp,'(a,i1)') 'pmata_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmata(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmata(ipm))
        write(cstemp,'(a,i1)') 'h1mat_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_h1mat(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_h1mat(ipm))
      enddo
      if (debug) write (6,*) 'rohf_h2e3 step 4'
c
      dims(1) = nvec
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1     
c
c ... jochen: allocate two components for g_dens here
c
      do ipm = 1,2
        if (.not. nga_create (MT_DBL, 3, dims, 'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('rohf_h2e: could not allocate g_dens',555,
     &     GA_ERR)
        call ga_zero(g_dens(ipm))
      end do                    ! ipm = 1,2
      if (debug) write (6,*) 'rohf_h2e3 step 5' 
c
      if (.not. nga_create (MT_DBL, 3, dims, 'CPKS dens',chunk, g_work))
     & call errquit('rohf_h2e: could not allocate g_work',555,GA_ERR)
      call ga_zero(g_work)
c
      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcv')
      if (nopen .gt. 0) then
         g_dcp = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dcp')
         g_dpv = ga_create_atom_blocked(geom, basis, 'rohf_h2e: dpv')
      endif
      if (debug) write (6,*) 'rohf_h2e3 step 6'
c
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      nfock = 0
      do ivec = 1, nvec
        nfock          = nfock + 1
c
c ... jochen: note that the arrays jfac and kfac will be modified
c       in subroutine shell_fock_build: jfac * 4 and kfac * 2.
c       here, we need to assign a set of factors for a symmetric
c       (ipm = 1) and an antisymmetric density matrix (ipm = 2)
c
        dbl_mb(k_jfac(1)+nfock-1)    =  1.0d0 ! CV
        dbl_mb(k_jfac(2)+nfock-1)    =  0.0d0 ! CV
        
        dbl_mb(k_kfac(1)+nfock-1)    = -0.5d0*xc_xfac1 ! CV
        dbl_mb(k_kfac(2)+nfock-1)    = -0.5d0*xc_xfac1 ! CV

c       jfac(nfock)    =  1.0d0   ! CV
c       kfac(nfock)    = 1.0d0   ! CV
c       jfac(nfock)    =  4.0d0   ! CV
c       kfac(nfock)    = -1.0d0   ! CV
        if (nopen .gt. 0) then
c ... jochen: the next part does not yet work. I left it as in hessv2 but
c         added the array indices for jfac and kfac
          nfock          = nfock + 1
          if (noskew) then
            dbl_mb(k_jfac(1)+nfock-1)    =  1.0d0 ! CP
          else
            dbl_mb(k_jfac(2)+nfock-1)    =  0.0d0 ! CP
          endif
          dbl_mb(k_kfac(1)+nfock-1)    = -0.5d0*xc_xfac1 ! CP
          nfock          = nfock + 1
          dbl_mb(k_jfac(1)+nfock-1) =  0.0d0 ! PV
          dbl_mb(k_kfac(1)+nfock-1) =  0.5d0*xc_xfac1 ! PV
        endif
c       
c       Compute CV, PV & CP "densities" from argument vector
c       
c ... jochen: skip this part and place a subroutine call instead.
c       it calculates the perturbed density matrix in the AO basis.
c       I keep this source code here for reference; it is left
c       unmodified from the version of rohf_hessv2 that this
c       subroutine was created from.
c       
c       note: for closed shell system nfock = nvec
        if (debug) write (6,*) 'rohf_h2e3 step 7 ivec=',ivec
        goto 666
c ... jochen end
c       
        call ga_zero(g_tmp1)
        call ga_zero(g_tmp2)
        call ga_zero(g_dcv)
        call ga_vec_to_mat(g_tmp1, 1, nvir, 1, nclosed,
     $     g_x, 1, ivec)
        call ga_matmul_patch('n','t',one,zero,
     $     g_tmp1, 1, nvir, 1, nclosed,
     $     g_movec, 1, nclosed, 1, nbf,
     $     g_tmp2, 1, nvir, 1, nbf )
        call ga_matmul_patch('n','n',one,zero,
     $     g_movec, 1, nbf, voff, nmo,
     $     g_tmp2, 1, nvir, 1, nbf,
     $     g_dcv, 1, nbf, 1, nbf)
c       
        if (nopen .gt. 0) then
c         
c         PV density
c         
          xoff = nclosed*nvir + 1
          xend = nclosed*nvir + nopen*nvir
          call ga_zero(g_tmp1)
          call ga_zero(g_tmp2)
          call ga_zero(g_dpv)
          call ga_copy_patch('n',
     $       g_x, xoff, xend, ivec, ivec,
     $       g_tmp1, 1, nvir, 1, nopen )
          call ga_matmul_patch('n', 't', one, zero,
     $       g_tmp1, 1, nvir, 1, nopen,
     $       g_movec, ooff, oend, 1, nbf,
     $       g_tmp2, 1, nvir, 1, nbf )
          call ga_matmul_patch('n', 'n', one, zero,
     $       g_movec, 1, nbf, voff, nmo,
     $       g_tmp2, 1, nvir, 1, nbf,
     $       g_dpv, 1, nbf, 1, nbf )
c         
c         CP density
c         
          call ga_zero(g_tmp1)
          call ga_zero(g_tmp2)
          call ga_zero(g_dcp)
          xoff = (nclosed+nopen)*nvir + 1
          xend = (nclosed+nopen)*nvir + nclosed*nopen
          call ga_copy_patch('n', g_x, xoff, xend, ivec, ivec,
     $       g_tmp1, 1, nopen, 1, nclosed )
          call ga_matmul_patch('n', 't', one, zero,
     $       g_tmp1, 1, nopen, 1, nclosed,
     $       g_movec, 1, nclosed, 1, nbf,
     $       g_tmp2, 1, nopen, 1, nbf )
          call ga_matmul_patch('n', 'n', one, zero,
     $       g_movec, 1, nbf, ooff, oend,
     $       g_tmp2, 1, nopen, 1, nbf,
     $       g_dcp, 1, nbf, 1, nbf )
c         
        endif
c        
c       We actually use dcp+dpv and dcp-dpv to build fock matrix
c       
        if (nopen .gt. 0) then
          call ga_dadd(one, g_dcp, one, g_dpv, g_tmp1)
          call ga_dadd(one, g_dcp,mone, g_dpv, g_tmp2)
          call ga_copy(g_tmp1, g_dcp)
          call ga_copy(g_tmp2, g_dpv)
        endif
c       
        
        if (noskew) then
          call ga_symmetrize(g_dcv)
        else
          call ga_antisymmetrize(g_dcv)
        endif
        if (nopen .gt. 0) then
          if (noskew) then
            call ga_symmetrize(g_dpv)
            call ga_symmetrize(g_dcp)
          else
            call ga_antisymmetrize(g_dpv)
            call ga_antisymmetrize(g_dcp)
          endif
        endif
c       
        alo(1) = ivec
        ahi(1) = ivec
        call nga_copy_patch('N',g_dcv,blo,bhi,g_dens,alo,ahi)
        if (nopen .gt. 0) then
c         Note that this is NOT correct!
          call nga_copy_patch('N',g_dpv,blo,bhi,g_dens,alo,ahi)
          call nga_copy_patch('N',g_dcp,blo,bhi,g_dens,alo,ahi)
        endif
c       
c ... jochen: if the preceeding code were exectuted then
c       g_dens would now hold the perturbed density matrix in the.
c       AO basis. Instead, we place
c       a call to one of the aoresponse subroutines:
        
  666   continue
        
        do ipm = 1,2
          call ga_zero(g_h1mat(ipm))
          call ga_vec_to_mat(g_h1mat(ipm), 1, nvir, 1, nclosed,
     $       g_x(ipm), 1, ivec) ! g_h1mat now holds A-matrix for
        end do                  ! ipm = 1,2
        call ga_sync()
        if (debug) write (6,*) 'calling CalcPerturbedTDPmat1'
c       note: the last argument tells it not to use an occ-occ
c             block to build the density marix.
        call CalcPerturbedTDPmat1
     &     (2, g_pmats, g_pmata, g_h1mat, g_movec, nbf, nclosed,
     &     nvir, nmo, .false., .false.,
     &     limag, .false.)  ! density matrix -> pmats
        if (debug) write (6,*) 'back from CalcPerturbedTDPmat1'

c next 2 lines for debugging only, to force uncoupled CPKS
c
ccc        call ga_zero(g_pmats(1))
ccc        call ga_zero(g_pmats(2))


        call ga_zero(g_pmata(1))
        call ga_zero(g_pmata(2))
c
c       
c       Instead of P(+) and P(-) which are both non-symmetric for
c       non-zero frequency
c       we will work with a symmetrized (S) and an antisymmetrized (A)
c       component, calculate F(S) and F(A), respectively, and construct
c       the Fock operators F(+/-) afterwards from F(S) +/- F(A).
c       If it works for the skew-symmetric density matrix of NMR then
c       it should work for this problem here, too

c       note: here is one of those scalings by 1/4
c       needed to get the correct final results
        call ga_sync()
        call ga_scale(g_pmats(1),0.25d0)
        call ga_scale(g_pmats(2),0.25d0)
c       
        alo(1) = ivec
        ahi(1) = ivec
c
c       we need to take care here of the symmetry of the density
c       matrices depending on whether the perturbation is real
c       or purely imaginary.
c
c       this works for real, symmetric, perturbations
c
        if (.not.limag) then
c       calculate P(S) = [P(+) + P(-)]/2
        call nga_add_patch (0.5d0, g_pmats(1), blo, bhi,
     &     0.5d0, g_pmats(2), blo, bhi,
     &     g_dens(1), alo, ahi)
c       caluclate P(A) = [-P(+) + P(-)]/2  (wrong results
c                                          with opposite sign ...)
        call nga_add_patch (-0.5d0, g_pmats(1), blo, bhi,
     &     0.5d0, g_pmats(2), blo, bhi,
     &     g_dens(2), alo, ahi)
c
        else
c
c       this here is for imaginary, antisymmetric, perturbations
c       calculate P(S) = [P(+) - P(-)]/2
        call nga_add_patch (0.5d0, g_pmats(1), blo, bhi,
     &   -0.5d0, g_pmats(2), blo, bhi,
     &     g_dens(1), alo, ahi)
c       caluclate P(A) = -[P(+) + P(-)]/2  ! sign needs to be determined
c                                          
        call nga_add_patch (-0.5d0, g_pmats(1), blo, bhi,
     &     -0.5d0, g_pmats(2), blo, bhi,
     &     g_dens(2), alo, ahi)
c
      endif                     ! limag
      call ga_sync()

c       
        if (debug) then 
          write (6,*) 'rohf_h2e3 step 8 ivec=',ivec
          write (6,*) 'P(S) and P(A) after CalcPerturbedTDPmat1'
          call nga_print_patch(g_dens(1),alo,ahi,1)
          call nga_print_patch(g_dens(2),alo,ahi,1)
        endif
c
      enddo                     ! ivec = 1,nvec

c      stop  ! FA-12-02-13

c     
      if(oprint) then
        call ga_print(g_dens(1))
        call ga_print(g_dens(2))
      endif
c     
      if (.not. ga_destroy(g_tmp1)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_dcv)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      if (nopen .gt. 0) then
         if (.not. ga_destroy(g_tmp3)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_dcp)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
         if (.not. ga_destroy(g_dpv)) call errquit('rohf_h2e: GA?',0,
     &       GA_ERR)
      endif
c
      if (debug) write (6,*) 'rohf_h2e3 step 9'
c
c     Fock build. dims are the same as for g_pmat
c
c ... jochen: fock matrix has now two components. allocate both
      do ipm = 1,2
        if (.not. nga_create (MT_DBL, 3, dims, 'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_fock',555,
     &     GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo                     ! ipm = 1,2
      if (debug) write (6,*) 'rohf_h2e3 step 10'
c
c ... jochen: next line is a comment from a previous version:
c     noskew is true, except for nmr perturbations (skew symmetric)
c
c ... jochen: we call the fock matrix builder twice here, once with
c     a symmetric, and once with an antisymmetric density matrix
c
c     last argument is "noskew", true if P-matrix is symmetric,
c     false if P-matrix is antisymmetric. The output Fock matrices
c     will be symmetrized and antisymmetrized, respectively.
c
c     case I: symmetric density matrix
      if (.not.cam_exch) then
        call shell_fock_build(geom, basis, 0, nfock, 
     $   dbl_mb(k_jfac(1)), dbl_mb(k_kfac(1)),
     $   tol2e, g_dens(1), g_fock(1), .true.)
         if (debug) write (6,*) 'rohf_h2e3 step 11'
      else   ! for attenuated calculations
        call shell_fock_build_cam(geom, basis, 0, nfock, 
     $   dbl_mb(k_jfac(1)), dbl_mb(k_kfac(1)),
     $   tol2e, g_dens(1), g_fock(1), .true.)
      end if
c
c     case II: antisymmetric density matrix. note that jfac for this
c     case is zero, i.e. there is no Coulomb term because there is
c     no density. However, there are still exchange contributions from
c     HF exchange. For pure DFT we just skip the call to save CPU time.
c
      if (abs(xc_xfac1).gt.tenm6) then
       if (.not.cam_exch) then
        call shell_fock_build(geom, basis, 0, nfock, 
     $   dbl_mb(k_jfac(2)), dbl_mb(k_kfac(2)),
     $   tol2e, g_dens(2), g_fock(2), .false.)
         if (debug) write (6,*) 'shell_fock_build called with P(A)'
       else  ! for attenuated calculations
        call shell_fock_build_cam(geom, basis, 0, nfock, 
     $   dbl_mb(k_jfac(2)), dbl_mb(k_kfac(2)),
     $   tol2e, g_dens(2), g_fock(2), .false.)
       end if
      else 
        call ga_zero (g_fock(2))
        if (debug) write (6,*) 'shell_fock_build skipped for P(A)'
      end if
      if (debug) write (6,*) 'rohf_h2e3 step 12'
c
c ... jochen: next, we need to undo the symmetrization in order to
c     recover F(+) and F(-) separately. THIS HAS NOT YET BEEN ADAPTED
c     TO THE ROHF CASE!!!
c     apparently the density matrices are not needed any further

      ifock = 0
      do ivec = 1,nvec
        ifock = ifock + 1
        alo(1) = ivec
        ahi(1) = ivec 
c       use g_pmats for temp storage of the fock matrices
        call nga_copy_patch ('N',g_fock(1),alo,ahi,g_pmats(1),blo,bhi)
        call nga_copy_patch ('N',g_fock(2),alo,ahi,g_pmats(2),blo,bhi)

        if (debug) then
          write (6,*) 'F(S) and F(A), ivec = ',ivec         
          call ga_print(g_pmats(1))
          call ga_print(g_pmats(2))
        end if
c
c ... jochen: the following should not make a diff. for PURE DFT
c        call ga_zero(g_pmats(2)) ! for testing ONLY !
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        pre_factor = 1.0d0
        call ga_sync()
        if (.not.limag) then
c         real perturbation: 
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_fock(1), alo, ahi)
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       -pre_factor, g_pmats(2), blo, bhi,
     &       g_fock(2), alo, ahi)
        else
c         imaginary perturbation:
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_fock(1), alo, ahi)
          call nga_add_patch (-pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_fock(2), alo, ahi)  
        end if

c
        if (debug) write (6,*) 'rohf_h2e3 step 13'

        if (debug) call nga_print_patch(g_fock(1),alo,ahi,1)
        if (debug) call nga_print_patch(g_fock(2),alo,ahi,1)
c
c       open-shell parts are yet to be implemented here ....
c       
      enddo                     ! ivec = 1,nvec
        
c     
c ... jochen: the next statement has a wonderfully self-explanatory
c     exit message that must ke kept this way by all means
      if (oskel) 
     $   call errquit('rohf_h2e3: no way',0, UNKNOWN_ERR)
c     
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'rohf_h2e3: tmp1')
c
c ... jochen: I think there is a bug in hessv2: g_dcv is used below
c     but it has already been deallocated. So we'll allocate it agai
      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e3: dcv')
c
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'rohf_h2e3: tmp2')
      if (nopen .gt. 0) then
        g_tmp3 = ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp3')
      endif
      if (debug) write (6,*) 'rohf_h2e3 step 14'
c     
      ifock = 0
c     
c     start loop over components of perturbing field
c     
      call ga_sync()
      do ivec = 1, nvec
        ifock = ifock + 1       ! History ... copy into tmp arrays
        alo(1) = ivec
        ahi(1) = ivec
c       
c ... jochen: add loop over Fock matrix components +/- here
c       
        do ipm = 1,2
c         
          call nga_copy_patch('N',g_fock(ipm),alo,ahi,g_tmp1,blo,bhi)
          call nga_copy_patch('N',g_dens(ipm),alo,ahi,g_dcv,blo,bhi)
c         
c ... jochen: Fock matrix is now in g_tmp1,
c         density matrix in g_dcv.
c         For the open shell code below, I assume that the ifock index
c         counts additional alpha and beta components of the
c         ROHF Fock matrix
c         
          if (nopen .gt. 0) then 
c ... jochen: this should not be executed.
c           I haven't adapted this
c           except for assigning (ipm) components to prevent compiling
c           problems
            if (debug) write (6,*) 'WAAAAHHHHH this code must not be
     &         executed: hessv3'
            if (ipm.eq.1) ifock = ifock + 1
            alo(1) = ifock
            ahi(1) = ifock
            call nga_copy_patch('N',g_fock(ipm),alo,ahi,
     &         g_tmp2,blo,bhi)
            call nga_copy_patch('N',g_dens(ipm),alo,ahi,g_dcp,blo,bhi)
            if (ipm.eq.1) ifock = ifock + 1
            alo(1) = ifock
            ahi(1) = ifock
            call nga_copy_patch('N',g_fock(ipm),alo,ahi,
     &         g_tmp3,blo,bhi)
            call nga_copy_patch('N',g_dens(ipm),alo,ahi,
     &         g_dpv,blo,bhi)
          endif                 ! nopen.gt.0

          if (debug) write (6,*) 'rohf_h2e3 step 15 ivec,ipm',ivec, ipm
c         
c          P      =  4(ij|kl) - (ik|jl) - (il|kj)
c           ij,kl
c     
c          K      =  (ik|jl) + (il|kj)
c           ij,kl
c     
c          cv         cv          pv   cp
c          Z   =  2P.[D  ]  +  P.[D  + D  ]
c     
c          pv          cv           cp   pv
c          Z   =  0.5d0*Z   + 0.5*K.[D  - D  ]
c          
c          cp          cv           cp   pv
c          Z   =  0.5d0*Z   - 0.5*K.[D  - D  ]
c          
c         Add the Fock matrices together overwriting the density
c         matrices to form the results above
c         
          if (debug) write (6,*) 'g_tmp1,g_dcv',g_tmp1,g_dcv
          call ga_copy(g_tmp1, g_dcv)
          if (debug) write (6,*) 'rohf_h2e3 step 16 ivec,ipm',ivec, ipm
          call ga_sync()
          call ga_dscal(g_dcv, two)
          call ga_sync()
          if (nopen .gt. 0) then
            call ga_dadd(one,  g_dcv,  one, g_tmp2, g_dcv)
            call ga_dadd(half, g_dcv, half, g_tmp3, g_dpv)
            call ga_dadd(half, g_dcv,mhalf, g_tmp3, g_dcp)
          endif
          if (debug) write (6,*) 'rohf_h2e3 step 17 ivec,ipm',ivec, ipm
c         
c         Closed-Virtual bit
c         
          call ga_matmul_patch('n','n',two,zero,
     $       g_dcv, 1, nbf, 1, nbf, ! Fock matrix
     $       g_movec, 1, nbf, 1, nclosed, ! MO coefficients
     $       g_tmp1, 1, nbf, 1, nclosed )    
          call ga_sync()
          call ga_matmul_patch('t','n',one,zero,
     $       g_movec, voff, nmo, 1, nbf, ! MO coefficients
     $       g_tmp1, 1, nbf, 1, nclosed, ! result from step 1 
     $       g_tmp2, 1, nvir, 1, nclosed ) ! vir-occ Fock matrix
          call ga_sync()
          xoff = 1
          xend = nvir*nclosed
c         
c ... jochen: the next line has one of those factors of 4
c         again ... here we add 4 times F_ai to the
c         vector g_ax here
          
          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed,
     $       g_ax(ipm), xoff, ivec, four, '+')
          if (debug) write (6,*) 'rohf_h2e3 step 18 ivec,ipm',ivec, ipm
          
*         call ga_dadd_patch( four, g_tmp2, 1, nvir, 1, nclosed,
*         $        one, g_ax, xoff, xend, ivec, ivec,
*         $        g_ax, xoff, xend, ivec, ivec )
c         
c ... jochen: next block only for non-closed-shell systems
c         for which this routine is not supposed to work anyway,
c         according to all those coments left by other developers.
c         I haven't made any changes in the nopen.gt.0 part
c... jochen: by now, g_Ax(ipm) should hold the desired result
c         
          if (nopen .gt. 0) then
            if (debug) write (6,*) 'WAAAAHHHHH this code must not be
     &         executed: hessv3'
c           
c           Open-Virtual bit
c           
            call ga_matmul_patch('n','n',two,zero,
     $         g_dpv, 1, nbf, 1, nbf,
     $         g_movec, 1, nbf, ooff, oend,
     $         g_tmp1, 1, nbf, 1, nopen )
            call ga_matmul_patch('t','n',one,zero,
     $         g_movec, voff, nmo, 1, nbf,
     $         g_tmp1, 1, nbf, 1, nopen,
     $         g_tmp2, 1, nvir, 1, nopen )
            xoff = xend + 1
            xend = xend + nvir*nopen
            call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nopen,
     $         g_ax(ipm), xoff, ivec, four, '+')
*           call ga_dadd_patch( four, g_tmp2, 1, nvir, 1, nopen,
*           $           one, g_ax, xoff, xend, ivec, ivec,
*           $           g_ax, xoff, xend, ivec, ivec)
c           
c           Closed-Open bit
c           
            call ga_matmul_patch('n','n',two,zero,
     $         g_dcp, 1, nbf, 1, nbf,
     $         g_movec, 1, nbf, 1, nclosed,
     $         g_tmp1, 1, nbf, 1, nclosed )
            call ga_zero(g_tmp2)
            call ga_matmul_patch('t','n',one,zero,
     $         g_movec, ooff, oend, 1, nbf,
     $         g_tmp1, 1, nbf, 1, nclosed,
     $         g_tmp2, 1, nopen, 1, nclosed )
            xoff = xend + 1
            xend = xend + nopen*nclosed
            call ga_mat_to_vec(g_tmp2, 1, nopen, 1, nclosed,
     $         g_ax(ipm), xoff, ivec, four, '+')
*           call ga_dadd_patch( four, g_tmp2, 1, nopen, 1, nclosed,
*           $           one, g_ax, xoff, xend, ivec, ivec,
*           $           g_ax, xoff, xend, ivec, ivec)
            
          endif                 ! nopen.gt.0
          if (debug) write (6,*) 'rohf_h2e3 step 19 ivec,ipm',ivec, ipm
          
        enddo                   ! ipm = 1,2 loop over +/- components
      enddo                     ! ivec = 1,nvec loop over field dirs.
c
c     dallocate memory:
c
      if (.not.ma_chop_stack(l_jfac(1)))
     &   call errquit('rohf_hessv3: cannot chop stack',0, MA_ERR)

      if (debug) call ga_print(g_Ax(1))
      if (debug) call ga_print(g_Ax(2))
      
c...  jochen: deallocate temp arrays added earlier
      do ipm = 1,2
        if (.not.ga_destroy(g_pmats(ipm)))
     &     call 
     &     errquit('rohf_hessv3: ga_destroy failed g_pmats',
     &     0,GA_ERR)
        
        if (.not.ga_destroy(g_pmata(ipm)))
     &     call 
     &     errquit('rohf_hessv3: ga_destroy failed g_pmata',
     &     0,GA_ERR)

       if (.not.ga_destroy(g_h1mat(ipm)))
     &     call 
     &     errquit('rohf_hessv3: ga_destroy failed g_h1mat',
     &     0,GA_ERR)

      enddo                     ! ipm = 1,2
c
c ... jochen: also added a loop over ipm here:     
      do ipm = 1,2
        if (.not. ga_destroy(g_dens(ipm))) call errquit
     &     ('rohf_h2e: GA?',0,
     &     GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_h2e: GA?',0,
     &     GA_ERR)
      enddo                     ! ipm = 1,2
c     
      if (.not. ga_destroy(g_tmp1)) call errquit('rohf:destroy?',0,
     &   GA_ERR)
c ... jochen: 
      if (.not. ga_destroy(g_dcv)) call errquit('rohf:destroy?',0,
     &   GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('rohf:destroy?',0,
     &   GA_ERR)
      if (nopen .gt. 0) then
        if (.not. ga_destroy(g_tmp3)) call errquit('rohf:destroy?',0,
     &     GA_ERR)
      endif
      if (.not. ga_destroy(g_work)) call errquit('rohf:destroy?',0,
     &   GA_ERR)
c     
      end
