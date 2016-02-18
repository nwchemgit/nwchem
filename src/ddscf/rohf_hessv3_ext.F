c ... jochen: this comes from rohf_hessv2.F but has everything
c     "v2" replaced by "v3", for use with frequency-dependent response
c     and is called from the cphf_solve3 part of the cphf code
c
c ... jochen: further mods were made to accomodate the situation that
c     we might have damping in the response. That causes all quantities to
c     have an imaginary part, too

      subroutine rohf_hessv3_ext(
     &                       acc, 
     &                       g_x, 
     &                       g_ax, 
     &                       g_x_im, 
     &                       g_Ax_im,
     &                       omega, 
     &                       limag, 
     &                       lifetime, 
     &                       gamwidth, 
     &                       ncomp)
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c     Note.- Equivalent to rohf_hessv3() (routine written by 
c            J. Autschbach) but using optimizing routines.
 
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
      external rohf_hessv_xx3_ext
c
c     ================================================================

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development
c      debug = (.true. .and. ga_nodeid().eq.0) ! for code development

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

      call rohf_hessv_xx3_ext( basis, geom, nbf, nmo,
     $     nclosed, nopen,
     $     pflg, 
     &     g_movecs, 
     &     oskel, noskew,
     $     crohf_g_fcv, crohf_g_fpv, crohf_g_fcp,
     $     acc, lshift, 
     &     g_x, g_ax, 
     &     g_x_im, g_Ax_im, 
     &     omega, limag,
     &     lifetime, gamwidth, ncomp)

      if (debug) write (6,*) 'back from rohf_hessv_xx3'
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

      subroutine rohf_hessv_xx3_ext( 
     &     basis, geom, 
     &     nbf, nmo, nclosed, nopen, 
     $     pflg,
     $     g_movecs, 
     &     oskel, noskew, 
     &     g_fcv, g_fpv, g_fcp,
     $     acc, 
     &     lshift, 
     &     g_x, g_ax, 
     &     g_x_im, g_Ax_im, 
     &     omega, 
     &     limag,
     &     lifetime, 
     &     gamwidth, 
     &     ncomp)
c     Note.- Improvements to this subroutine done by
c            Fredy W. Aquino, Northwestern University (Oct 2012)
c            Using rohf_hessv_2e2_opt() instead of old rohf_hessv_2e2()
c                  rohf_hessv_2e3_opt() instead of old rohf_hessv_2e3() 

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

      double precision omega, gamwidth, 
     &                 wls, wlsim
      logical limag, lifetime
      integer ncomp
      double precision omg(ncomp),gam(ncomp),coeffw
      logical debug
      external rohf_hessv_2e3_opt 

c
c     =================================================================

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

c      debug = (.true. .and. ga_nodeid().eq.0) ! for code development
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

      if (pflg .gt. 0) then
        coeffw=4.0d0 ! r-dft
        omg(1)=-omega
        omg(2)= omega
        gam(1)=-gamwidth
        gam(2)= gamwidth
        if (.not.lifetime) then
c        no damping: initialize Ax with terms proportional omega
         do ipm=1,ncomp
          wls = lshift + coeffw * omg(ipm)
          call ga_dadd(wls,g_x(ipm),0.d0,g_ax(ipm),g_ax(ipm))
         enddo ! end-loop-ipm
        else                    ! lifetime
c        take care of damping here: Re and Im are coupled by gamwidth
         do ipm=1,ncomp
          wls   = lshift + coeffw * omg(ipm)
          wlsim = -coeffw * gam(ipm) 
          call ga_dadd(wls  ,g_x(ipm),
     &                 wlsim,g_x_im(ipm),
     &                       g_ax(ipm))
          wls   =  coeffw * gam(ipm) 
          wlsim = lshift + coeffw * omg(ipm)
          call ga_dadd(wls  ,g_x(ipm),
     &                 wlsim,g_x_im(ipm),
     &                       g_ax_im(ipm))  
         enddo ! end-lopp-ipm
        endif                   ! .not.lifetime
        call ga_sync()
        if (debug) write (6,*) 'calling rohf_hessv_1e'

c ============== debug g_ax ==================== START
        if (debug) then
          do ipm=1,ncomp
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_ax_re-0-c(',ipm,')------ START' 
            call ga_print(g_ax(ipm))
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_ax_re-0-c(',ipm,')------ END'
          enddo ! end-loop-ipm
         if (lifetime) then
          do ipm=1,2
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_ax_im-0-c(',ipm,')------ START' 
           call ga_print(g_ax_im(ipm))
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_ax_im-0-c(',ipm,')------ END'
          enddo ! end-loop-ipm
         endif ! end-if-lifetime
        endif ! end-if-debug
c ============== debug g_ax ==================== END
c 
c       next: add (e_a - e_i) times A (also called U) matrix to Ax
        do ipm=1,ncomp
         call rohf_hessv_1e(basis,geom,        ! in : handles
     &                      nmo,nclosed,nopen, ! in : (nmo,nocc) nopen=0 for closed shell
     $                      g_fcv,g_fpv,g_fcp, ! in : densities
     $                      g_x(ipm),          ! in : g_x
     &                      g_ax(ipm))         ! out: 1e contrib to Ax
        enddo ! end-loop-ipm
        if (lifetime) then
        do ipm=1,ncomp
         call rohf_hessv_1e(basis,geom,        ! in : handles
     &                      nmo,nclosed,nopen, ! in : (nmo,nocc) nopen=0 for closed shell
     $                      g_fcv,g_fpv,g_fcp, ! in : densities
     $                      g_x_im(ipm),       ! in : g_x_im
     &                      g_ax_im(ipm))      ! out: 1e contrib to Ax_im
        enddo ! end-loop-ipm
        endif ! end-if-lifetime
      endif                     ! pflg.gt.0
c ============== debug g_ax ==================== START
      if (debug) then
       do ipm=1,ncomp
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-0-d(',ipm,')------ START' 
         call ga_print(g_ax(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-0-d(',ipm,')------ END'
       enddo ! end-loop-ipm
       if (lifetime) then
       do ipm=1,2
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-0-d(',ipm,')------ START' 
         call ga_print(g_ax_im(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-0-d(',ipm,')------ END'
       enddo ! end-loop-ipm
       endif ! end-if-lifetime
      endif ! end-if-debug
c ============== debug g_ax ==================== END
      if (pflg .gt. 1) then
c
c       the next call basically uses the current guess for the solution
c       vector x (in g_x, which is the perturbed density matrix in the
c       MO basis) and calculates the perturbed Fock operator in the MO basis.
c       real and imaginary part of that Fock operator can be handled
c       separately here

        if (ncomp.gt.1) then    ! call 2e code for dynamic case

           call rohf_hessv_2e3_opt(
     &                  basis,   ! in: basis handle
     &                  geom,    ! in: geom  handle
     &                  nbf,     ! in: nr. basis functions
     &                  nmo,     ! in: nr. MOs
     &                  nclosed, ! in: nr. double occupied MOs
     &                  nopen,   ! in: nr. single occupied MOs
     $                  g_movecs,! in: MO coefficients
     &                  oskel,   ! in: =.true. ->
     &                  noskew,  ! in: =.true. -> symmetric density matrix
     &                  g_x,     ! in: 
     &                  g_x_im,  ! in: 
     &                  g_ax,    ! in/out: Hessian product
     &                  g_ax_im, ! in/out: Hessian product
     &                  acc,     ! in: accuracy Fock construction
     &                  limag,   ! in: =.true. -> imaginary component allowed
     &                  lifetime)! in: =.true. -> RE-IM =.false -> RE

        else                    ! call static 2e code

             call rohf_hessv_2e2_opt(
     &                              basis,     ! in : basis handle
     &                              geom,      ! in : geom  handle
     &                              nbf,       ! in : nr. basis functions
     &                              nmo,       ! in : nr. MOs vecs
     &                              nclosed,   ! in : nr. occupied MOs 
     &                              nopen,     ! in : nr. open shells (unpaired e's)
     $                              g_movecs,  ! in : MO vec coeffs 
     &                              oskel,     ! in :
     &                              noskew,    ! in : symm density ?
     &                              acc,       ! in : accuracy Fock construction   
     &                              g_x(1),    ! in : RE g_x 
     &                              g_x_im(1), ! in : IM g_x
     &                              g_ax(1),   ! in/ou : RE g_ax Hessian product
     &                              g_ax_im(1),! in/ou : IM g_ax Hessian product
     &                              lifetime)

        endif                   ! ncomp
        
      endif                     ! pflg.gt.1
c     
      end

      subroutine rohf_hessv3_cmplx(
     &                  acc,       ! in : accuracy
     &                  g_z,       ! in : z
     &                  g_Az1,     ! in : Az1
     &                  nsub,      ! in
     &                  omega,     ! in :
     &                  limag,     ! in :
     &                  lifetime,  ! in : 
     &                  gamwidth,  ! in :
     &                  ncomp,     ! in : nr. components (+/-)
     &                  iter_cphf) ! in: iteration nr. in cphf
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c     Note.-  Modifying rohf_hessv3 to handle complex vars

      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "stdio.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "global.fh"
c     
c     $Id$
      integer ncomp
      integer g_z(ncomp), ! [input]  A-matrix elements for density matrix
     &        g_Az(ncomp) ! [output] Perturbed Fock operator
      integer g_Az1,
     &        nsub
      double precision acc, omega, gamwidth
      logical limag, lifetime    
      integer gtype,grow,gcol,
     &        growp,gcolp,ipm,iter_cphf
      logical oprint, debug
      external rohf_hessv_xx3_cmplx
      

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
      call ga_inquire(g_z(1),gtype,grow,gcol) 

      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv3_cmplx: invalid vector length-1',
     &                  0,UNKNOWN_ERR)

      goto 177 ! skip-for-the moment
      if (growp.ne.crohf_vlen)
     $     call errquit('rohf_hessv3_cmplx: invalid vector length-2',
     &                  0,UNKNOWN_ERR)
      if (gcol.ne.gcolp)
     $     call errquit('rohf_hessv3_cmplx: invalid no. of vectors',
     &                  0,UNKNOWN_ERR)
 177  continue
c     
c     Call internal routine
c  
      if (debug) write (6,*) 'calling rohf_hessv_xx3_cmplx'

      call rohf_hessv_xx3_cmplx(
     &     g_z,         ! in :
     &     g_Az1,
     &     nsub,
     &     ncomp,       ! in : nr. components
     &     basis,       ! in : basis handle
     &     geom,        ! in : geom  hande
     &     nbf,         ! in : nr. basis functions
     &     nmo,         ! in : nr. MOs
     $     nclosed,     ! in : nr. occupied MOs
     &     nopen,       ! in : nr. 1/2 occupied MOs (=0 for closed shells)
     $     pflg, 
     &     g_movecs,    ! in : MO coefficients 
     &     oskel, 
     &     noskew,
     $     crohf_g_fcv, ! in : closed-virtual density
     &     crohf_g_fpv, ! in : open-virtual   density (not used)
     &     crohf_g_fcp, ! in : closed-open    density (not used)
     $     acc, 
     &     lshift, 
     &     omega,
     &     limag,
     &     lifetime,
     &     gamwidth,
     &     iter_cphf)  ! in : iteration nr. in cphf cycle

      if (debug) write (6,*) 'back from rohf_hessv_xx3_cmplx'

c
c     Zap numbers much smaller than acc to ensure hard zeroes 
c     remain unpolluted ... cannot use a threshold larger than the
c     integral accuracy since can break symmetry in non-abelian groups
c     Also must ensure that the threshold tends to zero to permit
c     tight convergence.
c
c ... jochen: screen components
c ++++++++++++++++++++++++++++++++++++++++++
c =====> WARNING-UNFINISHED ROUTINE <=======
c ++++++++++++++++++++++++++++++++++++++++++
c     I need to adapt it for complex (later) --> zapping routine
c      do ipm = 1,ncomp
c        call ga_screen(g_ax(ipm),
c     &                 max(min(acc*acc,acc*0.01d0,1d-12),1d-16))
c        if (lifetime) call ga_screen(g_ax_im(ipm), 
c     &       max(min(acc*acc,acc*0.01d0,1d-12),1d-16))
c      enddo

      end

      subroutine rohf_hessv_xx3_cmplx(
     &               g_z,     ! in :
     &               g_Az1,   ! ou : product: A x z
     &               nsub,    ! in : nr. subspace
     &               ncomp,   ! in : nr. components
     &               basis,   ! in : basis handle
     &               geom,    ! in : geom  hande
     &               nbf,     ! in : nr. basis functions
     &               nmo,     ! in : nr. MOs
     &               nclosed, ! in : nr. occupied MOs
     &               nopen,   ! in : nr. 1/2 occupied MOs (=0 for closed shells)
     $               pflg,
     $               g_movecs,! in : MO coefficients 
     &               oskel, 
     &               noskew, 
     &               g_fcv,   ! in : closed-virtual density
     &               g_fpv,   ! in : open-virtual   density (not used)
     &               g_fcp,   ! in : closed-open    density (not used)
     $               acc, 
     &               lshift, 
     &               omega, 
     &               limag,
     &               lifetime, 
     &               gamwidth,
     &               iter_cphf)
c
c     Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c     Date   : 03-24-12
c     -> modified from rohf_hessv_xx3()
c     $Id$

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bgj.fh"
      integer ipm,ncomp,iter_cphf
      integer ncompmx
      parameter(ncompmx=2)
      integer g_z(ncompmx), ! [input]  A-matrix elements
     &        g_Az(ncompmx) ! [output] Perturbed Fock operator 
      integer g_Az1,nsub,
     &        alo(3),ahi(3),
     &        n1,maxsub,
     &        m1,m2,p1,p2,
     &        pretty
      integer g_movecs   
      integer g_fcv,
     &        g_fpv,
     &        g_fcp
      integer g_xreim,g_Axreim ! scratch GA arrays
      integer basis,geom
      integer nbf,nmo,nvir,
     &        nclosed,nopen,nreim,n
      integer pflg
      double precision acc,lshift
      integer vlen, nvec,g_tmp,gtype
      double precision omega,gamwidth,
     &                 omg(ncompmx), 
     &                 gam(ncompmx),
     &                 wls,wlsim
      double complex wls_cmplx
      logical oskel, noskew
      logical limag,lifetime
      logical debug
c     DIM/QM JEM
      integer g_x(ncomp), g_x_im(ncomp)
      integer g_Ax(ncomp), g_Ax_im(ncomp)
      logical ldimqm
c
      external rohf_hessv_2e3_opt_cmplx,
     &         rohf_hessv_2e2_opt_cmplx,
     &         getreorim,getreorim1,
     &         conv2complex3,conv2complex4
c
c     =================================================================

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

      if (debug) write (6,*) 'hessv3: omega =',omega
      if (debug) write (6,*) 'hessv3: limag =',limag
      if (debug) write (6,*)
     &   'hessv3: lifetime, gamwidth =',lifetime, gamwidth

c ----- Clear sub-block ---- START
       call ga_inquire(g_Az1,gtype,n1,maxsub)
       call ga_inquire(g_z(1),gtype,n,nvec) ! get (n1,nvec)
c
        alo(1)=1
        ahi(1)=n1
        alo(2)=nsub+1
        ahi(2)=nsub+nvec
        call nga_zero_patch(g_Az1,alo,ahi)
c ----- Clear sub-block ---- END
c
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
      if (pflg .gt. 0) then
         omg(1)=-omega
         omg(2)= omega
         gam(1)=-gamwidth
         gam(2)= gamwidth
c        take care of damping here: Re and Im are coupled by gamwidth
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
          wls   = lshift + 4.0d0 * omg(ipm)
          wlsim = -4d0 * gam(ipm) 
          wls_cmplx=dcmplx(wls,-wlsim)
          call ga_copy_patch('n',g_z(ipm),1 ,n ,1 ,nvec,
     &                           g_Az1   ,m1,m2,p1,p2)
          call ga_scale_patch(g_Az1,m1,m2,p1,p2,wls_cmplx)
          m1=m1+n
          m2=m2+n
         enddo ! end-lopp-ipm

        if (debug) then
          m1=1
          m2=n
          p1=nsub+1
          p2=nsub+nvec
          do ipm=1,ncomp
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_Az1-c(',ipm,')------ START' 
            call ga_print(g_Az1)
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_Az1-c(',ipm,')------ END'
           m1=m1+n
           m2=m2+n
          enddo ! end-loop-ipm
        endif ! end-if-debug
c 
c       next: add (e_a - e_i) times A (also called U) matrix to Ax
       call ga_inquire(g_z(1),gtype,n,nvec) ! get (n,nvec)
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx3_cmplx: g_xreim',0,0,g_xreim))
     $   call errquit('hessv_xx3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)     
        if (.not. ga_create(MT_DBL,n,nvec,
     &     'hessv_xx3_cmplx: g_xreim',0,0,g_Axreim))
     $   call errquit('hessv_xx3_cmplx: failed allocating g_xreim', 
     &                n,GA_ERR)  

       nvir = nmo - nclosed - nopen   
       do nreim=1,2 ! loop in RE,IM
        do ipm=1,ncomp
         call getreorim(g_xreim,  ! out : real or im arr
     &                  g_z(ipm), ! in  : = complx(g_xre,g_xim)
     &                  nvir,     ! in  : nr. virtual  MOs
     &                  nclosed,  ! in  : nr. occupied MOs
     &                  nreim)    ! in  : =1 -> re =2 -> im
         call getreorim1(g_Axreim,! out : real or im arr
     &                   g_Az1,   ! in  : = complx(g_xre,g_xim)
     &                   nsub,    ! in  : subblock index
     &                   ipm,     ! in  : = 1,2 to access slctd component
     &                   nvir,    ! in  : nr. virtual  MOs
     &                   nclosed, ! in  : nr. occupied MOs
     &                   nreim)   ! in  : =1 -> re =2 -> im
         call rohf_hessv_1e(
     &                  basis,geom,        ! in : handles
     &                  nmo,nclosed,nopen, ! in : (nmo,nocc) nopen=0 for closed shell
     $                  g_fcv,g_fpv,g_fcp, ! in : densities
     &                  g_xreim,
     &                  g_Axreim)
c ++++++++ update g_Az +++++++++ START
         call conv2complex4(g_Az1,   ! out: = history matrix complex
     &                      g_Axreim,! in : real      arr
     &                      nsub,    ! in  : subblock index
     &                      ipm,     ! in  : = 1,2 to access slctd component
     &                      nvir,    ! in  : nr. virtual  MOs
     &                      nclosed, ! in  : nr. occupied MOs
     &                      nreim)   ! in  : =1 -> re =2 -> im
c ++++++++ update g_Az +++++++++ END
        enddo ! end-loop-ipm
       enddo ! end-loop-nreim
        if (.not. ga_destroy(g_xreim))  call errquit
     &     ('hessv_xx3_cmplx: g_xreim',0, GA_ERR)
        if (.not. ga_destroy(g_Axreim)) call errquit
     &     ('hessv_xx3_cmplx: g_xreim',0, GA_ERR)
      endif                     ! pflg.gt.0
c ============== debug g_ax ==================== START
      if (debug) then
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_Az1-d------ START' 
            call ga_print(g_Az1)
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_Az1-d------ END'
      endif ! end-if-debug
c ============== debug g_ax ==================== END
      if (pflg .gt. 1) then
c
c       the next call basically uses the current guess for the solution
c       vector x (in g_x, which is the perturbed density matrix in the
c       MO basis) and calculates the perturbed Fock operator in the MO basis.
c       real and imaginary part of that Fock operator can be handled
c       separately here

        if (ncomp.gt.1) then    ! call 2e code for dynamic case
           call rohf_hessv_2e3_opt_cmplx(
     &                  g_z,      ! in :
     &                  g_Az1,    ! out: history of g_Az
     &                  nsub,     ! in :
     &                  basis,    ! in : basis handle
     &                  geom,     ! in : geom  handle
     &                  nbf,      ! in : nr. basis functions
     &                  nmo,      ! in : nr. MOs
     &                  nclosed,  ! in : nr. double occupied MOs
     &                  nopen,    ! in : nr. single occupied MOs
     $                  g_movecs, ! in : MO coefficients
     &                  oskel,    ! in : =.true. ->
     &                  noskew,   ! in : =.true. -> symmetric density matrix
     &                  acc,      ! in : accuracy Fock construction
     &                  limag,    ! in : =.true. -> imaginary component allowed
     &                  lifetime) ! in : =.true. -> RE-IM =.false -> RE
        else                      ! call static 2e code
           call rohf_hessv_2e2_opt_cmplx(
     &                          g_z(1),  ! in :
     &                          g_Az1,   ! out: history of g_Az
     &                          nsub,    ! in :
     &                          basis,   ! in : basis handle
     &                          geom,    ! in : geom  handle
     &                          nbf,     ! in : nr. basis functions
     &                          nmo,     ! in : nr. MOs vecs
     &                          nclosed, ! in : nr. occupied MOs 
     &                          nopen,   ! in : nr. open shells (unpaired e's)
     $                          g_movecs,! in : MO vec coeffs 
     &                          oskel,   ! in :
     &                          noskew,  ! in : symm density ?
     &                          acc,     ! in : accuracy Fock construction   
     &                          lifetime)
        endif                   ! ncomp
        
      endif                     ! pflg.gt.1
c
c   DIM/QM JEM
      if (.not.rtdb_get(bgj_get_rtdb_handle(), 'dimqm:lrsp', MT_LOG,
     $                  1, ldimqm)) ldimqm = .false.
      if(ldimqm) then
c       Transforming complex arrays to 2 real arrays to fit DIM structure
c       This could be cleaned up with a substantial rework of dimqm_addop
        do ipm = 1, ncomp
          if (.not. ga_create(MT_DBL,n,nvec,'dimqm: g_x',
     $                        0,0,g_x(ipm)))
     $      call errquit('dimqm: failed allocating g_x',n,GA_ERR)
          if (.not. ga_create(MT_DBL,n,nvec,'dimqm: g_x',
     $                        0,0,g_x_im(ipm)))
     $      call errquit('dimqm: failed allocating g_x_im',n,GA_ERR)
          if (.not. ga_create(MT_DBL,n,nvec,'dimqm: g_Ax',
     $                        0,0,g_Ax(ipm)))
     $      call errquit('dimqm: failed allocating g_Ax',n,GA_ERR)
          if (.not. ga_create(MT_DBL,n,nvec,'dimqm: g_Ax_im',
     $                        0,0,g_Ax_im(ipm)))
     $      call errquit('dimqm: failed allocating g_Ax_im',n,GA_ERR)
          call getreorim(g_x(ipm), ! out : real or im arr
     &                   g_z(ipm), ! in  : = complx(g_xre,g_xim)
     &                   nvir,     ! in  : nr. virtual  MOs
     &                   nclosed,  ! in  : nr. occupied MOs
     &                   1)        ! in  : =1 -> re =2 -> im
          call getreorim(g_x_im(ipm), ! out : real or im arr
     &                   g_z(ipm),    ! in  : = complx(g_xre,g_xim)
     &                   nvir,        ! in  : nr. virtual  MOs
     &                   nclosed,     ! in  : nr. occupied MOs
     &                   2)           ! in  : =1 -> re =2 -> im
          call getreorim1(g_Ax(ipm),  ! out : real or im arr
     &                    g_Az1,      ! in  : = complx(g_xre,g_xim)
     &                    nsub,       ! in  : subblock index
     &                    ipm,        ! in  : = 1,2 to access slctd component
     &                    nvir,       ! in  : nr. virtual  MOs
     &                    nclosed,    ! in  : nr. occupied MOs
     &                    1)          ! in  : =1 -> re =2 -> im
          call getreorim1(g_Ax_im(ipm), ! out : real or im arr
     &                    g_Az1,        ! in  : = complx(g_xre,g_xim)
     &                    nsub,         ! in  : subblock index
     &                    ipm,          ! in  : = 1,2 to access slctd component
     &                    nvir,         ! in  : nr. virtual  MOs
     &                    nclosed,      ! in  : nr. occupied MOs
     &                    2)   ! in  : =1 -> re =2 -> im
        end do
        call dimqm_addop(g_x, g_x_im, g_Ax, g_Ax_im,
     &                   ncomp, limag, lifetime)
        do ipm = 1,ncomp
          call conv2complex4(g_Az1,     ! out: = history matrix complex
     &                       g_Ax(ipm), ! in : real      arr
     &                       nsub,      ! in  : subblock index
     &                       ipm,       ! in  : = 1,2 to access slctd component
     &                       nvir,      ! in  : nr. virtual  MOs
     &                       nclosed,   ! in  : nr. occupied MOs
     &                       1)         ! in  : =1 -> re =2 -> im
          call conv2complex4(g_Az1,        ! out: = history matrix complex
     &                       g_Ax_im(ipm), ! in : real      arr
     &                       nsub,         ! in  : subblock index
     &                       ipm,          ! in  : = 1,2 to access slctd component
     &                       nvir,         ! in  : nr. virtual  MOs
     &                       nclosed,      ! in  : nr. occupied MOs
     &                       2)            ! in  : =1 -> re =2 -> im
          if (.not. ga_destroy(g_x(ipm)))  call errquit
     &      ('dimqm: failed destroying g_x',ipm, GA_ERR)
          if (.not. ga_destroy(g_x_im(ipm)))  call errquit
     &      ('dimqm: failed destroying g_x_im',ipm, GA_ERR)
          if (.not. ga_destroy(g_Ax(ipm)))  call errquit
     &      ('dimqm: failed destroying g_Ax',ipm, GA_ERR)
          if (.not. ga_destroy(g_Ax_im(ipm)))  call errquit
     &      ('dimqm: failed destroying g_Ax_im',ipm, GA_ERR)

        end do
        
      end if
c     
      end

      subroutine rohf_hessv_2e3_opt_cmplx(
     &                  g_z,
     &                  g_Az1,   ! in: (n1,maxsub) history of Az matrix (large matrix)
     &                  nsub,    ! in: point to (n1,nvec) block to be updated in g_Az1
     &                  basis,   ! in: basis handle
     &                  geom,    ! in: geom  handle
     &                  nbf,     ! in: nr. basis functions
     &                  nmo,     ! in: nr. MOs
     &                  nclosed, ! in: nr. double occupied MOs
     &                  nopen,   ! in: nr. single occupied MOs
     $                  g_movec, ! in: MO coefficients
     &                  oskel,   ! in: =.true. ->
     &                  noskew,  ! in: =.true. -> symmetric density matrix
     &                  acc,     ! in: accuracy Fock construction
     &                  limag,   ! in: =.true. -> imaginary component allowed
     &                  lifetime)! in: =.true. -> RE-IM =.false -> RE
c
c   Purpose: Optimization of rohf_hessv_2e3()
c   Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c   Date   : 03-28-12
c   Note1.- Modifying rohf_hessv_2e3, reducing computation of
c          two-electron integrals by putting together 
c          symmetric+antisymmetric perturbed densities: g_dens(ipm) ipm=1,2
c          and doing a single call to shell_fock_build() when using
c          the routine shell_fock_build2().
c   Note2.- size(g_Az1)=(n1,maxsub)  n1=ncomp*n  maxsumb=maxiter*nvec
c          ncomp=2 (+/-) n=nvir*nocc maxiter=10 (usually) nvec=3 (x,y,z)
c          nsub, will point to next (n1,nvec) block to be updated 
c          --> The purpose of using g_Az1 is to reduce usage of memory.
c              by doing it I skip using g_Az(ipm) ipm=2 (n1,nvec) complex matrix.

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
ccccccccccccccc This code does NOT work for open shell!!!!!ccccccccccccccccc
      integer g_z(2),g_Az(2)
      integer g_Az1,nsub
      integer m1,m2,p1,p2,pretty
      integer basis, geom       ! basis & geom handle
      integer nclosed,nvir,nopen! Basis size and occupation
      integer nbf,nmo           ! No. of linearly dependent MOs
      integer g_movec           ! MO coefficients
      logical oskel

      double precision acc      ! Accuracy of "Fock" construction
      logical limag             ! imaginary perturbation?    
      integer voff,xoff
      integer nmul,gtype,
     &        ivec,vlen,nvec
      integer g_tmp,g_tmp1,
     &        g_dens(4),g_fock(4),
     &        g_xreim(2)
      double precision tol2e
      logical odebug,oprint,noskew
      integer dims(3),chunk(3), 
     &        alo(3),ahi(3), 
     &        blo(2),bhi(2)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         get_undosymm_fock,update_ax_fock,
     &         shell_fock_build2,getreorim,
     &         update_gz_reorim,
     &         update_gz_reorim1 
      double precision one,zero,mone,four,half,mhalf,two,mtwo
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)
      integer ipm ! counter for density matrix components
      character*(255) cstemp
      integer g_pmats(2),g_pmata(2),g_h1mat(2),
     &        cnt,ind,indx(2,2),
     &        npol,nset,nblock,ncomp
      double precision tenm6,coef(2,2)
      parameter (tenm6 = 1d-6)
      logical lifetime,debug
      data npol /1/ ! for restricted calculations
      data indx /1,2, ! indx(1,1),indx(1,2)
     &           3,4/ ! indx(2,1),indx(2,2)
      ncomp=2 ! using two components
      call ga_inquire(g_z(1),gtype,vlen,nvec) ! out: nvec,vlen
      do ipm = 1,ncomp
       if (.not. ga_create(MT_DBL,vlen,nvec, 
     &      'hessv_2e3_opt_cmplx: g_xreim',0,0,g_xreim(ipm)))
     $   call errquit('rhessv_2e3_opt_cmplx: failed alloc g_xreim',
     &                nvec,GA_ERR)
      enddo ! end-loop-ipm
      nmul=1
      if (npol.eq. 2) nmul=2
      nset  =1 ! for RE          
      nblock=2 ! for RE
      if (lifetime) then
      nset  =2 ! for RE-IM
      nblock=4 ! for RE-IM
      endif

      if (oskel) 
     $   call errquit('rohf_h2e3: no way',0, UNKNOWN_ERR)

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,2001) limag,lifetime,nset,nblock
 2001   format('(limag,lifetime,nset,nblock)=(',
     &           L1,',',L1,',',i3,',',i3,')')
       endif
      endif ! end-if-debug

      oprint= util_print('rohf_hessv2',print_debug)

      if (nopen.ne.0) call errquit
     $     ('rohf_h2e3: does not work for open shells',nopen,
     &       UNKNOWN_ERR)
      odebug = util_print('rohf_hessv', print_debug)    
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1

      dims(1)  = nbf
      dims(2)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1

c ------ create scratch GA arrays (g_pmats,g_mata,g_h1mat---START
      do ipm = 1,ncomp
        write(cstemp,'(a,i1)') 'pmats_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmats(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmats(ipm))
        write(cstemp,'(a,i1)') 'pmata_',ipm
        write(cstemp,'(a,i1)') 'h1mat_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_h1mat(ipm))) call 
     &     errquit('rohf_h2e3: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_h1mat(ipm))
      enddo ! end-loop-ipm
c ------ create scratch GA arrays (g_pmats,g_mata,g_h1mat---END
      dims(1)  = nvec
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1   
      call ga_sync()  
      do ipm = 1,nblock ! =2 or 4
c ... allocate g_dens=[g_dens_re g_dens_im]
        if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_dens',555,
     &     GA_ERR)
        call ga_zero(g_dens(ipm))
c ... allocate g_fock=[g_fock_re g_fock_im]
        if (.not. nga_create (MT_DBL,3,dims,'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_fock',555,
     &     GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo ! end-loop-ipm
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF get_dens_reorim-RE'
        endif ! end-if-debug
c ---- Copy g_z --> g_x_reim ------ START
        do ipm=1,ncomp
         call ga_zero(g_xreim(ipm))
         call getreorim(g_xreim(ipm),! out : real or im arr
     &                  g_z(ipm),    ! in  : = complx(g_xre,g_xim)
     &                  nvir,        ! in  : nr. virtual  MOs
     &                  nclosed,     ! in  : nr. occupied MOs
     &                  1)           ! in  : =1 -> re =2 -> im
        enddo ! end-loop-ipm
c ---- Copy g_z --> g_x_reim ------ END
      call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    1,      ! in   : =1 1st block RE
     &                    g_xreim,! in   : RE
     &                    g_movec,! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : for ipol=1 -> istart=1
     &                    nclosed,! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    npol,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array - NOT USED 
     &                    g_h1mat)! in   : scratch GA array   

      if (lifetime) then
c ---- Copy g_z --> g_x_reim ------ START
        do ipm=1,ncomp
         call ga_zero(g_xreim(ipm))
         call getreorim(g_xreim(ipm),! out : real or im arr
     &                  g_z(ipm),    ! in  : = complx(g_xre,g_xim)
     &                  nvir,        ! in  : nr. virtual  MOs
     &                  nclosed,     ! in  : nr. occupied MOs
     &                  2)           ! in  : =1 -> re =2 -> im
        enddo ! end-loop-ipm
c ---- Copy g_z --> g_x_reim ------ END

       call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    2,      ! in   : =2 2nd block IM
     &                    g_xreim,! in   : IM
     &                    g_movec,! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : for ipol=1 -> istart=1
     &                    nclosed,! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    npol,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array - NOT USED   
     &                    g_h1mat)! in   : scratch GA array  
 
      endif ! end-if-lifetime
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_xreim(ipm))) call errquit(
     &     'rohf_hessv3: ga_destroy failed g_xreim',0,GA_ERR)  
       if (.not.ga_destroy(g_h1mat(ipm)))  call errquit(
     &      'rohf_hessv3: ga_destroy failed g_h1mat',0,GA_ERR)     
      enddo ! end-loop-ipm

      if (debug) then
        do ipm=1,nblock
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------START'
         call ga_print(g_dens(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------END'
        enddo ! end-loop-ipm
      endif ! end-if-debug

      call shell_fock_build2(g_fock,   ! out: Fock    matrices
     &                       g_dens,   ! in : density matrices
     &                       geom,     ! in : geom  handle
     &                       basis,    ! in : basis handle
     &                       nbf,      ! in : nr. basis functions
     &                       nvec,     ! in : nr. vecs (x,y,z)
     &                       npol,     ! in : npol=1 for R-DFT =2 for U-DFT
     &                       ncomp,    ! in : nr. components
     &                       nblock,   ! in : nr. of g_dens,g_fock blocks
     &                       .true.,   ! in : =.true. for symm dens
     &                       tol2e,    ! in :
     &                       debug)    ! in : = .true. -> debugging printouts

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-0(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-0(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug

      call get_undosymm_fock(
     &            g_fock,  ! in/ou: fock matrix
     &            nset,    ! in   : =1 g_x is real, =2 g_x is complex (g_x_re,g_x_im)
     &            nvec,    ! in   : nr. directions (x,y,z)
     &            nbf,     ! in   : nr. basis functions
     &            npol,    ! in   : nr. polarizations
     &            nmul,    ! in   : =1 npol=1 =2 npol=2 (acc. JK terms)
     &            g_pmats, ! in   : scratch GA array
     &            limag)   ! in   : =.true. imaginary comp. exists

c ------- Remove GA arrays:
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_pmats(ipm))) call errquit(
     &     'rohf_hessv3: ga_destroy failed g_pmats',0,GA_ERR)       
      enddo ! end-loop-ipm

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug
        
      if (debug) then
          if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_Az1-0------ START' 
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_Az1-0------ END'
      endif ! end-if-debug
c     
c     start loop over components of perturbing field
c   
      g_tmp = ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp')
      g_tmp1= ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp1')
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nclosed  
      do cnt=1,nset
       do ivec = 1, nvec
        alo(1) = ivec
        ahi(1) = ivec     
        do ipm = 1,ncomp ! loop over Fock matrix components +/- here     
          ind=indx(ipm,cnt)

          if (debug) then
            if (ga_nodeid().eq.0) then
             write(*,117) cnt,ivec,ipm,ind
 117         format('XX:(cnt,ivec,ipm,ind)=(',
     &              i3,',',i3,',',i3,',',i3,')')
            endif
          endif ! end-if-debug
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
c         Closed-Virtual bit   
          if (debug) then
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_fck -------- START'
           call ga_print(g_fock(ind)) 
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_fck -------- END'
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_vecs -------- START'
           call ga_print(g_movec) 
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_vecs -------- END'
          endif ! end-if-debug

          call ga_zero(g_tmp) 
          call nga_matmul_patch('n','n',four,zero,
     &                          g_fock(ind),alo,ahi,
     &                          g_movec    ,blo,bhi,
     &                          g_tmp      ,blo,bhi)

          if (debug) then
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- FnnCno -------- START'
           call ga_print(g_tmp) 
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- FnnCno -------- END'
          endif ! end-if-debug

          call ga_zero(g_tmp1)
          call ga_matmul_patch('t','n',one,zero,
     $                         g_movec,voff,nmo ,1,nbf,     ! MO coefficients
     $                         g_tmp  ,1   ,nbf ,1,nclosed, ! result from step 1 
     $                         g_tmp1 ,1   ,nvir,1,nclosed) ! vir-occ Fock matrix

         if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,3701) cnt,ivec,ipm
3701     format('----- CvnFnnCno(',i3,',',i3,',',i3,')------START')
          endif
           call ga_print(g_tmp1)
          if (ga_nodeid().eq.0) then
           write(*,3702) cnt,ivec,ipm
3702     format('----- CvnFnnCno(',i3,',',i3,',',i3,')------END')
          endif
         endif ! end-if-debug

          if      (cnt.eq.1) then

           if (debug) then
            if (ga_nodeid().eq.0)        
     &      write(*,*) '------- g_Az1-re-BEF(',ipm,')------ START' 
            call ga_print(g_Az1)
            if (ga_nodeid().eq.0)
     &      write(*,*) '------- g_Az1-re-BEF(',ipm,')------ END'
           endif ! end-if-debug

c Note.- The operation below does:
c        g_ax_re= g_ax_re + 4 [4 C^T F C]  --> I am not sure if this is right.

          call update_gz_reorim1(g_Az1,  ! out: = complx(g_xre,g_xim)
     &                           g_tmp1, ! in : real      arr 
     &                           1,      ! in  : =1 -> re =2 -> im
     &                           nsub,   ! in : index to sub-block in g_z
     &                           ipm,    ! in : = 1 or 2 index for component
     &                           vlen,   ! in : = nocc*nvir
     &                           four,   ! in  : scaling factor
     &                           nvir,
     &                           nclosed,
     &                           ivec)

          if (debug) then
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az1-re-AFT(',ipm,')------ START' 
           call ga_print(g_Az1)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az1-re-AFT(',ipm,')------ END'
          endif ! end-if-debug

          else if (cnt.eq.2) then

          if (debug) then
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az1-im-BEF(',ipm,')------ START' 
           call ga_print(g_Az1)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az1-im-BEF(',ipm,')------ END'
          endif ! end-if-debug

          call update_gz_reorim1(g_Az1,  ! out: = complx(g_xre,g_xim)
     &                           g_tmp1, ! in : real      arr 
     &                           2,      ! in  : =1 -> re =2 -> im
     &                           nsub,   ! in : index to sub-block in g_z
     &                           ipm,    ! in : = 1 or 2 index for component
     &                           vlen,   ! in : = nocc*nvir
     &                           four,   ! in  : scaling factor
     &                           nvir,
     &                           nclosed,
     &                           ivec)

          if (debug) then
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az1-im-AFT(',ipm,')------ START' 
           call ga_print(g_Az1)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az1-im-AFT(',ipm,')------ END'
          endif ! end-if-debug

          endif ! end-if-cnt   
        enddo ! end-loop-ipm lopp in +/- components
       enddo ! end-loop-ivec-loop in field directions
      enddo ! end-loop-cnt
      if (debug) then
       do ipm=1,ncomp
           if (ga_nodeid().eq.0)        
     &     write(*,*) '------- g_Az1-1(',ipm,')------ START' 
           call ga_print(g_Az1)
           if (ga_nodeid().eq.0)
     &     write(*,*) '------- g_Az1-1(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug

      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_hessv3: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm    
        if (.not.ga_destroy(g_tmp)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
        if (.not.ga_destroy(g_tmp1)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
      return     
      end

      subroutine rohf_hessv_2e3_opt(
     &                  basis,   ! in: basis handle
     &                  geom,    ! in: geom  handle
     &                  nbf,     ! in: nr. basis functions
     &                  nmo,     ! in: nr. MOs
     &                  nclosed, ! in: nr. double occupied MOs
     &                  nopen,   ! in: nr. single occupied MOs
     $                  g_movec, ! in: MO coefficients
     &                  oskel,   ! in: =.true. ->
     &                  noskew,  ! in: =.true. -> symmetric density matrix
     &                  g_x_re,  ! in: 
     &                  g_x_im,  ! in: 
     &                  g_ax_re, ! in/out: Hessian product
     &                  g_ax_im, ! in/out: Hessian product
     &                  acc,     ! in: accuracy Fock construction
     &                  limag,   ! in: =.true. -> imaginary component allowed
     &                  lifetime)! in: =.true. -> RE-IM =.false -> RE
c
c   Purpose: Optimization of rohf_hessv_2e3()
c   Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c   Date   : 03-28-12
c   Note.- Modifying rohf_hessv_2e3, reducing computation of
c          two-electron integrals by putting together 
c          symmetric+antisymmetric perturbed densities: g_dens(ipm) ipm=1,2
c          and doing a single call to shell_fock_build() when using
c          the routine shell_fock_build2().

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
ccccccccccccccc This code does NOT work for open shell!!!!!ccccccccccccccccc
      integer basis, geom       ! basis & geom handle
      integer nclosed,nvir,nopen! Basis size and occupation
      integer nbf,nmo           ! No. of linearly dependent MOs
      integer g_movec           ! MO coefficients
      logical oskel
      integer g_x_re(2),        !
     &        g_x_im(2)         ! Argument
      integer g_ax_re(2),       ! Hessian product
     &        g_ax_im(2)        ! Hessian product
      double precision acc      ! Accuracy of "Fock" construction
      logical limag             ! imaginary perturbation?    
      integer voff,xoff
      integer ivec,nvec,nmul,gtype,vlen
      integer g_tmp,g_tmp1,g_dens(4),g_fock(4)
      double precision tol2e
      logical odebug,oprint,noskew
      integer dims(3),chunk(3), 
     &        alo(3),ahi(3), 
     &        blo(2),bhi(2)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked,
     &         get_undosymm_fock,update_ax_fock,
     &         shell_fock_build2    
      double precision one,zero,mone,four,half,mhalf,two,mtwo
      double precision itol_floor, itol_ceil
      parameter(itol_floor=1.d-15, itol_ceil=1.d-3)
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)
      integer ipm ! counter for density matrix components
      character*(255) cstemp
      integer g_pmats(2),g_pmata(2),g_h1mat(2),
     &        cnt,ind,indx(2,2),
     &        npol,nset,nblock,ncomp
      double precision tenm6,coef(2,2)
      parameter (tenm6 = 1d-6)
      logical lifetime,debug
      data npol /1/ ! for restricted calculations
      data indx /1,2, ! indx(1,1),indx(1,2)
     &           3,4/ ! indx(2,1),indx(2,2)
      call ga_inquire(g_x_re(1),gtype,vlen,nvec) ! out: nvec,vlen

      ncomp=2 ! using two components
      nmul=1
      if (npol.eq. 2) nmul=2
      nset  =1 ! for RE          
      nblock=2 ! for RE
      if (lifetime) then
      nset  =2 ! for RE-IM
      nblock=4 ! for RE-IM
      endif

      if (oskel) 
     $   call errquit('rohf_h2e3: no way',0, UNKNOWN_ERR)

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,2001) limag,lifetime,nset,nblock
 2001   format('(limag,lifetime,nset,nblock)=(',
     &           L1,',',L1,',',i3,',',i3,')')
       endif
      endif ! end-if-debug

      oprint= util_print('rohf_hessv2',print_debug)

      if (nopen.ne.0) call errquit
     $     ('rohf_h2e3: does not work for open shells',nopen,
     &       UNKNOWN_ERR)
      odebug = util_print('rohf_hessv', print_debug)    
      tol2e = min(max(acc,itol_floor),itol_ceil)
      nvir = nmo - nclosed - nopen
      voff = nclosed + nopen + 1

      dims(1)  = nbf
      dims(2)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
c ------ create scratch GA arrays (g_pmats,g_mata,g_h1mat---START
      do ipm = 1,ncomp
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
      enddo ! end-loop-ipm
c ------ create scratch GA arrays (g_pmats,g_mata,g_h1mat---END
      dims(1)  = nvec
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1     
      do ipm = 1,nblock ! =2 or 4
c ... allocate g_dens=[g_dens_re g_dens_im]
        if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_dens',555,
     &     GA_ERR)
        call ga_zero(g_dens(ipm))
c ... allocate g_fock=[g_fock_re g_fock_im]
        if (.not. nga_create (MT_DBL,3,dims,'Fockv',chunk,
     &     g_fock(ipm)))
     &     call errquit('rohf_h2e3: could not allocate g_fock',555,
     &     GA_ERR)
        call ga_zero(g_fock(ipm))
      enddo ! end-loop-ipm

        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF get_dens_reorim-RE'
        endif ! end-if-debug

      call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    1,      ! in   : =1 1st block RE
     &                    g_x_re, ! in   : RE
     &                    g_movec,! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : for ipol=1 -> istart=1
     &                    nclosed,! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    npol,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array   
     &                    g_h1mat)! in   : scratch GA array   

         if (debug) then
          if (ga_nodeid().eq.0)
     &     write(*,*) 'BEF get_dens_reorim-IM'
         endif ! end-if-debug

      if (lifetime) then

       call get_dens_reorim(
     &                    g_dens, ! in/ou: perturbed density matrix
     &                    2,      ! in   : =2 2nd block IM
     &                    g_x_im, ! in   : IM
     &                    g_movec,! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    1,      ! in   : for ipol=1 -> istart=1
     &                    nclosed,! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    npol,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array   
     &                    g_h1mat)! in   : scratch GA array   

      endif ! end-if-lifetime

      if (debug) then
        do ipm=1,nblock
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------START'
         call ga_print(g_dens(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------g_dens(',ipm,')-------END'
        enddo ! end-loop-ipm
      endif ! end-if-debug

      call shell_fock_build2(g_fock,  ! out: Fock    matrices
     &                       g_dens,  ! in : density matrices
     &                       geom,    ! in : geom  handle
     &                       basis,   ! in : basis handle
     &                       nbf,     ! in : nr. basis functions
     &                       nvec,    ! in : nr. vecs (x,y,z)
     &                       npol,    ! in : npol=1 for R-DFT =2 for U-DFT
     &                       ncomp,   ! in : nr. components
     &                       nblock,  ! in : nr. of g_dens,g_fock blocks
     &                       .true.,  ! in : =.true. for symm dens
     &                       tol2e,   ! in :
     &                       debug)   ! in : = .true. -> debugging printouts

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-0(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-0(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug

      call get_undosymm_fock(
     &            g_fock,  ! in/ou: fock matrix
     &            nset,    ! in   : =1 g_x is real, =2 g_x is complex (g_x_re,g_x_im)
     &            nvec,    ! in   : nr. directions (x,y,z)
     &            nbf,     ! in   : nr. basis functions
     &            npol,    ! in   : nr. polarizations
     &            nmul,    ! in   : =1 npol=1 =2 npol=2 (acc. JK terms)
     &            g_pmats, ! in   : scratch GA array
     &            limag)   ! in   : =.true. imaginary comp. exists

c ------- Remove GA arrays:
      do ipm = 1,ncomp
       if (.not.ga_destroy(g_pmats(ipm))) call errquit(
     &     'rohf_hessv3: ga_destroy failed g_pmats',0,GA_ERR)       
       if (.not.ga_destroy(g_pmata(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_pmata',0,GA_ERR)
       if (.not.ga_destroy(g_h1mat(ipm)))  call errquit(
     &      'rohf_hessv3: ga_destroy failed g_h1mat',0,GA_ERR)
      enddo ! end-loop-ipm

      if (debug) then
       do ipm=1,nblock
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ START' 
         call ga_print(g_fock(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_fock-unsym(',ipm,')------ END'
       enddo ! end-loop-ipm
      endif ! end-if-debug
        
      if (debug) then
       do ipm=1,ncomp
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-0(',ipm,')------ START' 
         call ga_print(g_ax_re(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-0(',ipm,')------ END'
       enddo ! end-loop-ipm
       if (lifetime) then
       do ipm=1,ncomp
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-0(',ipm,')------ START' 
         call ga_print(g_ax_im(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-0(',ipm,')------ END'
       enddo ! end-loop-ipm
       endif ! end-if-lifetime
      endif ! end-if-debug
c     
c     start loop over components of perturbing field
c   
      g_tmp = ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp')
      g_tmp1= ga_create_atom_blocked(geom, basis,'rohf_h2e3: tmp1')
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nclosed  
      do cnt=1,nset
       do ivec = 1, nvec
        alo(1) = ivec
        ahi(1) = ivec      
        do ipm = 1,ncomp ! loop over Fock matrix components +/- here     
          ind=indx(ipm,cnt)

           if (debug) then
            if (ga_nodeid().eq.0) then
             write(*,117) cnt,ivec,ipm,ind
 117         format('XX:(cnt,ivec,ipm,ind)=(',
     &              i3,',',i3,',',i3,',',i3,')')
            endif
           endif ! end-if-debug
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
c         Closed-Virtual bit   
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- START'
         call ga_print(g_fock(ind)) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_fck -------- END'
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- START'
         call ga_print(g_movec) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- g_vecs -------- END'
        endif ! end-if-debug
          call ga_zero(g_tmp)    
          call nga_matmul_patch('n','n',four,zero,
     &                          g_fock(ind),alo,ahi,
     &                          g_movec    ,blo,bhi,
     &                          g_tmp      ,blo,bhi)
        if (debug) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- FnnCno -------- START'
         call ga_print(g_tmp) 
         if (ga_nodeid().eq.0)
     &    write(*,*) '--------- FnnCno -------- END'
        endif ! end-if-debug

          call ga_zero(g_tmp1)
          call ga_matmul_patch('t','n',one,zero,
     $                         g_movec,voff,nmo ,1,nbf,     ! MO coefficients
     $                         g_tmp  ,1   ,nbf ,1,nclosed, ! result from step 1 
     $                         g_tmp1 ,1   ,nvir,1,nclosed) ! vir-occ Fock matrix

         if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,3701) cnt,ivec,ipm
3701     format('----- CvnFnnCno(',i3,',',i3,',',i3,')------START')
          endif
           call ga_print(g_tmp1)
          if (ga_nodeid().eq.0) then
           write(*,3702) cnt,ivec,ipm
3702     format('----- CvnFnnCno(',i3,',',i3,',',i3,')------END')
          endif
         endif ! end-if-debug

          if      (cnt.eq.1) then

          if (debug) then
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_ax-re-BEF-------- START'
           call ga_print(g_ax_re(ipm)) 
           if (ga_nodeid().eq.0)
     &     write(*,*) '--------- g_ax-re-BEF-------- END'
          endif ! end-if-debug

c Note.- The operation below does:
c        g_ax_re= g_ax_re + 4 [4 C^T F C]  --> I am not sure if this is right.
          call ga_mat_to_vec(g_tmp1,1,nvir,1,nclosed,
     $                       g_ax_re(ipm),1,ivec,four,'+')

          if (debug) then
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-re-AFT-------- START'
           call ga_print(g_ax_re(ipm)) 
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-re-AFT-------- END'
          endif ! end-if-debug

          else if (cnt.eq.2) then

          if (debug) then
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-im-BEF-------- START'
           call ga_print(g_ax_im(ipm)) 
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-im-BEF-------- END'
          endif ! end-if-debug

          call ga_mat_to_vec(g_tmp1,1,nvir,1,nclosed,
     $                       g_ax_im(ipm),1,ivec,four,'+')

          if (debug) then
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-im-AFT-------- START'
           call ga_print(g_ax_im(ipm)) 
           if (ga_nodeid().eq.0)
     &      write(*,*) '--------- g_ax-im-AFT-------- END'
          endif ! end-if-debug

          endif ! end-if-cnt   
        enddo ! end-loop-ipm lopp in +/- components
       enddo ! end-loop-ivec-loop in field directions
      enddo ! end-loop-cnt

      if (debug) then
       do ipm=1,ncomp
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_re-1(',ipm,')------ START' 
         call ga_print(g_ax_re(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_re-1(',ipm,')------ END'
       enddo ! end-loop-ipm
       if (lifetime) then
       do ipm=1,ncomp
         if (ga_nodeid().eq.0)        
     &    write(*,*) '------- g_ax_im-1(',ipm,')------ START' 
         call ga_print(g_ax_im(ipm))
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_ax_im-1(',ipm,')------ END'
       enddo ! end-loop-ipm
       endif ! end-if-lifetime
      endif ! end-if-debug

      do ipm = 1,nblock
        if (.not. ga_destroy(g_dens(ipm))) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_dens',0,GA_ERR)
        if (.not. ga_destroy(g_fock(ipm))) call errquit
     &     ('rohf_hessv3: ga_destroy failed g_fock',0,GA_ERR)
      enddo ! end-loop-ipm    
        if (.not.ga_destroy(g_tmp)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
        if (.not.ga_destroy(g_tmp1)) call errquit(
     &      'rohf_hessv3: ga_destroy failed g_tmp',0,GA_ERR)
      return     
      end

      subroutine get_dens_reorim(
     &                    g_dens, ! out  : perturbed density matrix
     &                    cnt,    ! in/ou: counter of g_dens, =1 or 2
     &                    g_x,    ! in   : 
     &                    g_movec,! in   : MO coefficients
     &                    nbf,    ! in   : nr. basis functions
     &                    nmo,    ! in   : nr. MOs
     &                    istart, ! in   : shift nocc-nvirt block
     &                    nocc,   ! in   : nr. occupied MOs
     &                    nvir,   ! in   : nr. virtual  MOs 
     &                    nvec,   ! in   : nr. directions (x,y,z)
     &                    ipol,   ! in   : nr. polarizations
     &                    limag,  ! in   : = .true. imaginary allowed   
     &                    g_pmats,! in   : scratch GA array
     &                    g_pmata,! in   : scratch GA array - NOT USED
     &                    g_h1mat)! in   : scratch GA array   
c
c Purpose: Calculate perturbed density matrix
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c Date   : 03-15-12

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
      integer g_x(2)            ! Argument: g_x_re or g_x_im
      integer g_dens(*)         ! size= 2 RE only or 4 RE+IM
      integer nbf,nmo,nocc,nvir
      integer g_movec           ! MO coefficients
      integer dims(3),chunk(3), 
     &        alo(3),ahi(3), 
     &        blo(2),bhi(2)
      integer ivec,nvec,ipm,ncomp,
     &        ipol, ! =1 for Alpha =2 for Beta
     &        shift,cnt,ind,indx(2,2)
      character*(255) cstemp
      integer g_pmats(2),g_pmata(2),g_h1mat(2)
      integer istart,iend
      double precision tenm6,coef(2,2)
      parameter (tenm6 = 1d-6)
      logical limag,debug
      double precision one, zero, mone, 
     &                 four, half, mhalf, two, mtwo
      data indx /1,2, ! indx(1,1),indx(2,1)
     &           3,4/ ! indx(1,2),indx(2,2) 
      data ncomp/2/
      parameter (one=1.0d0, mone=-1.0d0, zero=0.0d0, four=4.0d0)
      parameter (half=0.5d0, mhalf=-0.5d0, two=2.0d0, mtwo=-2.0d0)
      external ga_vec_to_mat,
     &         CalcPerturbedTDPmat1_opt

      debug=.false. ! allow debugging printouts

c ----- Construct coeffs for P(S),P(A) ------- START
      coef(1,1)= 0.5d0
      coef(1,2)= 0.5d0
      coef(2,1)=-0.5d0
      coef(2,2)= 0.5d0
      if (limag) then
      coef(1,1)= 0.5d0
      coef(1,2)=-0.5d0
      coef(2,1)=-0.5d0
      coef(2,2)=-0.5d0
      endif ! end-if-limag
      if (debug) then
       if (ga_nodeid().eq.0) then
        write(*,10) limag,
     &              coef(1,1),coef(1,2),
     &              coef(2,1),coef(2,2)
 10     format('(limag,coeff)=(',L1,',',f15.8,',',
     &         f15.8,',',f15.8,',',f15.8,')')
       endif 
      endif
c ----- Construct coeffs for P(S),P(A) ------- END
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      shift=(ipol-1)*nvec
      iend = istart + nocc*nvir - 1     
      do ivec = 1, nvec
c       
c       Compute CV, PV & CP "densities" from argument vector
c       
c ... jochen: skip this part and place a subroutine call instead.
c       it calculates the perturbed density matrix in the AO basis.
c       I keep this source code here for reference; it is left
c       unmodified from the version of rohf_hessv2 that this
c       subroutine was created from. 
        do ipm = 1,ncomp
          call ga_zero(g_h1mat(ipm))
          call ga_copy_patch('n', ! Reshape vector into matrix 
     $                       g_x(ipm)    ,istart,iend,ivec,ivec,
     $                       g_h1mat(ipm),1     ,nvir,1   ,nocc)

        enddo ! end-loop-ipm

       if (debug) then
       do ipm=1,2
        if (ga_nodeid().eq.0)
     &   write(*,*) '----------g_h1mat(',ipm,')-----START'
        call ga_print(g_h1mat(ipm))
        if (ga_nodeid().eq.0)
     &   write(*,*) '----------g_h1mat(',ipm,')-----END'
       enddo ! end-loop-ip
       endif ! end-if-debug

        if (debug) then
         if (ga_nodeid().eq.0)
     &   write(*,*) 'In rohf_hessv_2e3: BEF CalcPerturbedTDPmat1'
         if (ga_nodeid().eq.0) then
          write(*,667) 2,nbf,nocc,nvir,nmo,limag
 667      format('(ncomp,nbf,nclosed,nvir,nmo,limag)=(',
     &          i3,',',i3,',',i3,',',i3,',',i3,',',L1,')')
         endif
        endif ! end-if -debug

              call CalcPerturbedTDPmat1_opt(
     &                 2,        ! in : nr. components to calculate
     &                 g_pmats,  ! out: density matrix      symmetrized
     &                 g_pmata,  ! out: density matrix  antisymmetrized
     &                 g_h1mat,  ! in :   perturbed MO coefficients
     &                 g_movec,  ! in : unperturbed MO coefficients
     &                 nbf,      ! in : nr. AOs
     &                 nocc,     ! in : nr. occupied MOs
     &                 nvir,     ! in : nr. virtual  MOs
     &                 nmo,      ! in : nr. MOs
     &                 .false.,  ! in : = .true. calc. (symm,antisymm)=(pmats,pmata)
     &                 .false.,  ! in : = .true. static response, dynamic otherwise
     &                  limag,   ! in : = .true. if amat is imaginary instead of real
     &                 .false.)  ! in : = .true. if amat contains occ-occ

         if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-1-------- START'
             call ga_print(g_pmats(1))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-1-------- END'      
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-2-------- START'
             call ga_print(g_pmats(2))
            if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_pmats-2-------- END'    
         endif ! end-if-debug

c next 2 lines for debugging only, to force uncoupled CPKS
c        call ga_zero(g_pmata(1))
c        call ga_zero(g_pmata(2))
c       
c       Instead of P(+) and P(-) which are both non-symmetric for
c       non-zero frequency
c       we will work with a symmetrized (S) and an antisymmetrized (A)
c       component, calculate F(S) and F(A), respectively, and construct
c       the Fock operators F(+/-) afterwards from F(S) +/- F(A).
c       If it works for the skew-symmetric density matrix of NMR then
c       it should work for this problem here, too
c       note: here is one of those ominous scalings by 1/4
c       needed to get the correct final results
        call ga_scale(g_pmats(1),0.25d0)
        call ga_scale(g_pmats(2),0.25d0)     
        alo(1) = shift+ivec
        ahi(1) = alo(1)
c       we need to take care here of the symmetry of the density
c       matrices depending on whether the perturbation is real
c       or purely imaginary.
c
c       this works for real, symmetric, perturbations
c       calculate P(S) = [ P(+) + P(-)]/2
c       calculate P(A) = [-P(+) + P(-)]/2  (wrong results
c                                          with opposite sign ...)
         do ipm=1,ncomp
          ind=indx(ipm,cnt)

          if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,11) cnt,ipm,ind,ivec,coef(ipm,1),coef(ipm,2)
 11        format('check-ind: (cnt,ipm,ind,ivec)=(',
     &            i3,',',i3,',',i3,',',i3,')',
     &            'coeff12=(',f15.8,',',f15.8,')')
          endif
          endif ! end-if-debug

          call nga_add_patch(coef(ipm,1),g_pmats(1) ,blo,bhi,
     &                       coef(ipm,2),g_pmats(2) ,blo,bhi,
     &                                   g_dens(ind),alo,ahi)
          if (debug) then
            if (ga_nodeid().eq.0) then
             write(*,*) '---g_dens-acc(',ind,',',ivec,')-------START'
            endif
            call ga_print(g_dens(ind))
            if (ga_nodeid().eq.0)
     &      write(*,*) '----g_dens-acc(',ind,',',ivec,')-------END'
          endif ! end-if-debug

         enddo ! end-loop-ipm   
      enddo ! end-loop-ivec 
      return
      end
