      subroutine cphf_solve3(rtdb, omega, lifetime, gamwidth)
c
c  Written by J. Autschbach, SUNY Buffalo
c  Extension to spin-unrestricted case 
c          by F. Aquino,     Northwestern University (Oct 2012)
c          03-15-12
c 
c jochen: modified version of cphf_solve2 to deal with 
c         frequency dependent properties. for this, array g_grad
c         will now be of dimension 2. g_work holds a copy of
c         that array and will also have dimension 2
c jochen: Further extended this subroutine to deal with
c         finite-lifetime damping. Arrays dealing with imaginary
c         parts have been added, indicated by suffix "_im"
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
#include "crohf.fh"
#include "cuhf.fh"
#include "stdio.fh"
#include "cosmo.fh"
c#include "dimqm.fh"
c!!! BGJ
#include "bgj.fh"
c!!! BGJ
      integer rtdb              ! [input] database handle
      double precision omega    ! [input] frequency
      logical lifetime          ! [input] damping switch
      double precision gamwidth ! [input] damping constant
c
      double precision eone,etwo,enrep,energy
      integer g_grad(2), g_work(2)
      integer g_grad_im(2), g_work_im(2)
      integer g_grad_cmplx(2), ! for complex calc.
     &        g_work_cmplx(2)  ! for complex calc.
      integer max_precond_iter,max_precond_sub
      double precision precond_acc
      integer precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
      logical oprint_solve
      
c ... jochen: modified the definition of some of the externals
c     here
      external rohf_hessv_precond3,rohf_hessv,rohf_hessv3
      external rohf_hessv_precond3_damp
      external uhf_precond,uhf_hessv,uhf_hessv2,
     &         uhf_hessv_precond3,
     &         uhf_hessv_precond3_damp,
     &         uhf_hessv3,
     &         uhf_hessv3_cmplx,
     &         uhf_hessv3_cmplx1
       logical oprint, oprint_conv, ocosmo
c     
      character*255 rhs_title
      character*(nw_max_path_len) cphf_rhs(2)
      character*(nw_max_path_len) cphf_sol(2)

      character*(nw_max_path_len) cphf_rhs_im(2)
      character*(nw_max_path_len) cphf_sol_im(2)

      logical file_write_ga,file_read_ga, file_read_ga_info
      external file_write_ga,file_read_ga, file_read_ga_info
      integer vlen, nvirt, nmos(2), nbf_vecs, nsets
      integer bases(2), nbases, nvec, ncol
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      double precision residual
      double precision ecosmo
      logical movecs_read_header, cphf_poliz
      external movecs_read_header, cphf_poliz
      logical int_normalize
      external int_normalize

      external ga_lkain_2cpl3_damp_cmplx_redmem,
     &         ga_lkain_2cpl3_damp_cmplx_redmem1,
     &         ga_lkain_2cpl3_ext,
     &         conv2complex,rohf_hessv3_cmplx,
     &         ga_lkain_2cpl3_damp

      integer nocc1(2),nvir1(2)
      external get_cconjugate
      double precision omg(2),  ! (-omega,omega)
     &                 gam(2),  ! (-gamwidth,gamwidth)
     &                 solveshift
      logical hf_job, debug, limag, lstatic
      integer ipm,ncomp,ndata,m1,m2,shift
      double precision small
      logical use_dimqm,lrsp,ldimqm,lfirst,dimqm_go
      logical xc_gotxc,xc_gothfx
      external xc_gotxc,xc_gothfx
      parameter (small = 1d-8)
c
c     ==================================================================

      debug = (.false. .and. ga_nodeid().eq.0)

c     debug = .true. ! to allow debugging printouts

      lstatic = (abs(omega).lt.small)
      use_dimqm=util_module_avail("dimqm")
      hf_job = (.not. xc_gotxc()).and.(.not.xc_gothfx())

c
c     Print levels assume that CPHF is being invoked from within
c     another application.  By default be informative.
c     At low print, just indicate convergence or failure.

      oprint = util_print('information',print_default) .and.
     $     ga_nodeid() .eq. 0
      oprint_conv = util_print('information',print_low) .and.
     $     ga_nodeid() .eq. 0
      oprint_solve = util_print('solve',print_default)
c
c     determine how many components we should treat:

      if (.not. rtdb_get(rtdb, 'cphf:ncomp', mt_int, 1, ncomp)) call
     &   errquit('cphf_solve3: failed to read ncomp ', 0, RTDB_ERR)

      if (ncomp.gt.1 .and. lstatic) then ! something's fishy here
        call errquit
     &     ('cphf_solve3: static response but more than one component?',
     &     ncomp, RTDB_ERR)
      elseif (ncomp.eq.1 .and. .not.lstatic) then ! fishy fishy
        call errquit
     &     ('cphf_solve3: dynamic response with one component?',
     &     ncomp, RTDB_ERR)
      endif
c
c     Must use the same vectors as the converged SCF

      if (hf_job) then
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $        call errquit('cphf: no SCF MO vectors',0, RTDB_ERR)
      else
         if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $        call errquit('cphf_solve3: no DFT MO vectors',0, RTDB_ERR)
#if 0
c !!! Hack to set up scftype properly for later code !!! ? need this?
         if (bgj_restricted()) then
            scftype = 'RHF'
         else
            scftype = 'UHF'
         endif
#endif
      endif
c      write(*,*)'*** cphf_solve3: scftype ',scftype
c     
      call scf_init(rtdb)
      if (.not.int_normalize(rtdb,basis))
     &    call errquit('cphf_solve3: int_normalize failed',911,
     &       UNKNOWN_ERR)
      nbases = 1
      bases(1) = basis
      if (bgj_have_j_fit()) then
        nbases = 2
        bases(2) = bgj_CD_bas_han()
        write(*,*) 'cphf_solve3: bases',bases
      endif
      call int_init(rtdb,nbases,bases)
      call schwarz_init(geom, basis)
      oskel=sym_number_ops(geom).gt.0
      call int_1e_cache_ga(basis,oskel)
c
c     Get the MO vectors
c
      if(.not.movecs_read_header(movecs_in,title_vecs,basis_vecs,
     $     scftype_vecs, nbf_vecs,nsets,nmos,2))
     $     call errquit('cphf: problem with movecs_read_header',1,
     &       UNKNOWN_ERR)
cnew hack replacing bgj-isms
      if (nsets.eq.1) then
         scftype = 'RHF'
      else
         scftype = 'UHF'
      endif
      if(nbf.ne.nbf_vecs)call errquit('cphf: MO vectors mismatch?',0,
     &       UNKNOWN_ERR)
      nmo = nmos(1)

      if (nmo .ne. nbf) then
         if (.not. ga_destroy(g_movecs(1))) call errquit('cphf: ga',1,
     &       GA_ERR)
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: tmp MOs a',
     $        0, 0, g_movecs(1))) call errquit('scf_v_g: tmp MOs', 0,
     &       GA_ERR)
         if (scftype .eq. 'UHF') then
            if (.not. ga_destroy(g_movecs(2)))call errquit('cphf: ga',2,
     &       GA_ERR)
*ga:1:0
            if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: tmp MOs b',
     $           0, 0, g_movecs(2))) call errquit('scf_v_g: tmp MOs', 0,
     &       GA_ERR)
         endif
      endif
c         
      call scf_movecs_read

      if (scftype.eq.'UHF') then 
         vlen = nalpha*(nmo-nalpha) + nbeta*(nmo-nbeta)
      else if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         nvirt = nmo - nclosed - nopen
         vlen = nclosed*(nopen+nvirt) + nopen*nvirt
      else
         call errquit('cphf: unknown SCF type',0, INPUT_ERR)
      endif

      lshift=0.0d0
c
c note: the original default value for "tol" was 1d-4
c
      if (.not. rtdb_get(rtdb, 'cphf:thresh', mt_dbl, 1, 
     $     precond_tol)) precond_tol=1d-4
      if (.not. rtdb_get(rtdb, 'cphf:maxiter', mt_int, 1,
     $     max_precond_iter)) max_precond_iter = 50
      maxiter = max_precond_iter
c
c ... jochen: open second rhs file first and check if it has data
      if (ncomp.gt.1) then
        call cphf_fname('cphf_rhs_2',cphf_rhs(2)) 
        if (.not. file_read_ga_info(cphf_rhs(2),rhs_title,ncol,nvec))
     $     call errquit('cphf_solve3: reading rhs 2 info failed?',0,
     &     UNKNOWN_ERR)
      endif
      call cphf_fname('cphf_rhs',cphf_rhs(1))
      if (.not. file_read_ga_info(cphf_rhs(1), rhs_title, ncol, nvec))
     $     call errquit('cphf_solve3: reading rhs info failed?',0,
     &       UNKNOWN_ERR)
      if (vlen .ne. ncol) 
     $     call errquit('cphf_solve3: vlen != ncol', ncol,
     &       UNKNOWN_ERR)
c
c ... jochen: open data file(s) with "imaginary" data
      if (lifetime) then
        if (ncomp.gt.1) then
          call cphf_fname('cphf_rhs_im_2',cphf_rhs_im(2)) 
          if (.not. file_read_ga_info(cphf_rhs_im(2), rhs_title,
     +       ncol, nvec))
     $       call errquit
     &       ('cphf_solve3: reading rhs_im 2 info failed?',0,
     &       UNKNOWN_ERR)
        endif
        call cphf_fname('cphf_rhs_im',cphf_rhs_im(1)) 
        if (.not. file_read_ga_info(cphf_rhs_im(1), rhs_title,
     +     ncol, nvec))
     $     call errquit('cphf_solve3: reading rhs_im 1 info failed?',0,
     &     UNKNOWN_ERR)        
      endif                     ! lifetime
c
c      control size of the subspace
c
       max_precond_sub = 10*nvec
      if (.not. rtdb_get(rtdb, 'cphf:maxsub', mt_int, 1,
     1  max_precond_sub)) max_precond_sub = 10*nvec
c
c     allocate g_grad: rhs vector for CPKS
c
      do ipm=1,ncomp
       if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: rhs',
     $     0,nvec,g_grad(ipm)))
     $     call errquit('cphf_solve3: ga_cre grad',vlen, GA_ERR)
      enddo ! end-loop-ipm
      if (lifetime) then
       do ipm=1,ncomp
       if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: rhs',
     $     0,nvec,g_grad_im(ipm)))
     $     call errquit('cphf_solve3: ga_cre grad',vlen, GA_ERR)
       enddo ! end-loop-ipm
      endif 
c
      pflg = 2 ! this variable is used in rohf_hessv_xx
      if (oprint) then
         write(luout,170) scftype
 170     format('    scftype = ',a)
         if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
            write(luout,171) nclosed, nopen
 171        format('    nclosed = ', i4/
     $             '      nopen = ', i4)
         else
            write(luout,172) nalpha, nbeta
 172        format('     nalpha = ', i4/
     $             '      nbeta = ', i4)
         endif
          
c ... jochen: no. of variables is 2*vlen in dynamic case,
c        and 4*vlen if dampling is applied
         if (lifetime) then
           write(luout,173) 2*ncomp*vlen, nvec, precond_tol,
     &        max_precond_iter
         else
           write(luout,173) ncomp*vlen, nvec, precond_tol,
     &        max_precond_iter
         endif
 173     format('  variables = ', i8/
     $          '  num. vecs = ', i8/
     $          '  tolerance = ', 1p,d9.1,0p/
     $          '  max. iter = ', i4/)
      endif
c
*ga:1:0
c
c ... jochen: here we read the CPKS rhs vector generated elsewhere
c     (A- or U-matrices, not yet divided by orbital energy diffs)

      do ipm=1,ncomp
       if(.not.file_read_ga(cphf_rhs(ipm),g_grad(ipm)))
     $     call errquit('cphf_solve3: could not read cphf_rhs 1',1,
     &       UNKNOWN_ERR)
      enddo ! end-loop-ipm
      if (lifetime) then
      do ipm=1,ncomp
       if(.not.file_read_ga(cphf_rhs_im(ipm),g_grad_im(ipm)))
     $     call errquit('cphf_solve3: could not read cphf_rhs 1',1,
     &       UNKNOWN_ERR)
      enddo ! end-loop-ipm
      endif ! end-if-lifetime

      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve3: work1',
     $   0,0,g_work(1)))
     &   call errquit('cphf_solve3: ga_cre work1',0, GA_ERR)
c     
c     Compute necessary Fock matrices
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve3: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 1))
     $     call errquit('cphf_solve3: rtdb_put of calc_type failed',0,
     &   RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $   call errquit('cphf_solve3: rtdb_put of j_derfit failed',0,
     &   RTDB_ERR)
c !!! BGJ
c
c     Get the cosmo contribution

      if (.not. rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo))
     $     ocosmo = .false.
      if (ocosmo) then
         cosmo_on = .true.
         cosmo_phase = 2
      endif
c     DIM/QM JEM
      if (use_dimqm) then
         call dimqm_used(ldimqm)
         call dimqm_getlrsp(lrsp)
      else
         ldimqm=.false.
         lrsp=.false.
      endif
      if (ldimqm .and. lrsp) then
        lfirst = .true.
        call dimqm_setlfirst(lfirst)
        dimqm_go = .false.
        call dimqm_setgo(dimqm_go)
      end if

      if(scftype .eq. 'UHF')then
        call uhf_init(rtdb)
        call uhf_energy(rtdb,g_movecs,eone,etwo,enrep,ecosmo,energy,
     $                  g_work(1))
      else
        call rohf_init(rtdb)
        call rohf_energy(g_movecs,eone,etwo,enrep,ecosmo,energy,
     $                   g_work(1),rtdb)
      endif
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $   call errquit('cphf_solve3: rtdb_put of xc_active failed',1,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('cphf_solve3: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
c !!! BGJ
c
c     == check the residual ==
      residual = sqrt(ga_ddot(g_work(1), g_work(1)))
      if(.not.rtdb_put(rtdb,'cphf_solve:residual',mt_dbl,1,residual))
     $   call errquit('cphf_solve: rtdb_put of residual failed',0,
     &     RTDB_ERR)
c
      if (ga_nodeid().eq.0) then
       write(luout,*) 'SCF residual: ', residual
      endif
      if(residual .gt.1d-2) then
         call ga_sync()
         call errquit('cphf_solve3:SCF residual greater than 1d-2',
     &              0,UNKNOWN_ERR)
      endif
c
c     destroy temp work array used for Fock matrix calculation
      if (.not. ga_destroy(g_work(1))) call errquit('cphf_solve3:ga?',0,
     &   GA_ERR)
c     
c ... jochen: allocate ncomp components for array g_work here and take
c             care of arrays needed in a damping calculation, too

      do ipm = 1,ncomp
        if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: work',
     $     0, 0, g_work(ipm)))
     &     call errquit('cphf_solve3: ga_cre work',ipm,
     &     GA_ERR)
        if (lifetime) then
          if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: work_im',
     $       0, 0, g_work_im(ipm)))
     +       call errquit('cphf_solve3: ga_cre work im',ipm,
     &       GA_ERR)
        endif
      enddo ! end-loop-ipm
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve3: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('cphf_solve3: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve3: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
c !!! BGJ

          do ipm = 1,ncomp
           call ga_copy_patch('n',g_grad(ipm),1,vlen,1,nvec, 
     $                            g_work(ipm),1,vlen,1,nvec)
           if (lifetime) then
           call ga_copy_patch('n',g_grad_im(ipm),1,vlen,1,nvec, 
     $                            g_work_im(ipm),1,vlen,1,nvec)
           endif
          enddo ! end-loop-ipm
c 

      if (scftype .eq. 'UHF') then

c        For property related CPHF. NMR has imaginary perturbation
c        which makes matrices skew symmetric. Setting noskew to .false.
c        enforces that.
         if (.not. rtdb_get(rtdb, 'cphf:skew', mt_log, 1, 
     $        noskew_uhf)) noskew_uhf = .true.

c        the dynamic density matrix builder needs to know if the perturbation
c        is imaginary (or antisymmetric). This should have been written
c        to the rtdb at some time (e.g. in aoresponse)

         if (.not. rtdb_get(rtdb, 'cphf:imag', mt_log, 1,
     &      limag)) then
           limag = .false.
           if (ga_nodeid().eq.0) then
             write (luout,*)
     &        'WARNING from cphf_solve3: could not read limag from RTDB'
             write (luout,*)
     &        'Assuming REAL Symmetric perturbation.'
           endif
         endif

         
          solveshift= 0.0d0
          omg(1)=-omega
          omg(2)= omega
          gam(1)=-gamwidth
          gam(2)= gamwidth
          if (.not.lifetime) then
           do ipm=1,ncomp
            call uhf_hessv_precond3(
     &                     g_work(ipm),
     &                     omg(ipm))         
           enddo ! end-loop-ipm
          else
           do ipm=1,ncomp
            call uhf_hessv_precond3_damp(
     &                     g_work(ipm), 
     &                     g_work_im(ipm),
     &                     omg(ipm),
     &                     gam(ipm))
           enddo ! end-loop-ipm
          endif
c  jbecca START -- Need to create the int. mat. here for UHF
         if (ldimqm .and. lrsp) call dimqm_prep(rtdb, omega)
c  jbecca END


         if (.not.lifetime) then
           if (debug) write (6,*) 'calling ga_lkain_2cpl3_ext'

           call ga_lkain_2cpl3_ext(
     &                      rtdb,
     &                      g_work   , g_grad, 
     &                      g_work_im, g_grad_im, 
     &                      uhf_hessv3,
     &                      uhf_hessv_precond3,
     $                      precond_tol,
     &                      max_precond_sub,
     &                      maxiter,
     &                      .false., 
     &                      oprint_solve, 
     &                      omega,
     &                      limag, 
     &                      lifetime, 
     &                      gamwidth, 
     &                      ncomp)


           if (debug) write (6,*) 'exiting ga_lkain_2cpl3_ext'

         else

           if (debug) write (6,*) 'calling ga_lkain_2cpl3_damp'

      do ipm = 1,ncomp
        if(.not.ga_create(mt_dcpl,vlen,nvec,'cphf_solve3: work_cmplx',
     $     0, 0, g_work_cmplx(ipm)))
     &     call errquit('cphf_solve3: g_work_cmplx',ipm,
     &     GA_ERR)
        if(.not.ga_create(mt_dcpl,vlen,nvec,'cphf_solve3: work_cmplx',
     $     0, 0, g_grad_cmplx(ipm)))
     &     call errquit('cphf_solve3: g_work_cmplx',ipm,
     &     GA_ERR)
      enddo ! end-loop-ipm
c ----- Create complex GA arrays ----------- END
           call conv2complex(g_work_cmplx,! out: = complx(g_xre,g_xim)
     &                       g_work,      ! in : real      arr
     &                       g_work_im,   ! in : imaginary arr
     &                       vlen,        ! in : n    rows
     &                       nvec,        ! in : nvec columns
     &                       ncomp)       ! in : nr. components
           call conv2complex(g_grad_cmplx,! out: = complx(g_xre,g_xim)
     &                       g_grad,      ! in : real      arr
     &                       g_grad_im,   ! in : imaginary arr
     &                       vlen,        ! in : n    rows
     &                       nvec,        ! in : nvec columns
     &                       ncomp)       ! in : nr. components

           nocc1(1)=nalpha
           nocc1(2)=nbeta
           nvir1(1)=nmo-nocc1(1)
           nvir1(2)=nmo-nocc1(2)

           call ga_lkain_2cpl3_damp_cmplx_redmem1(
     &              rtdb,
     &              g_work_cmplx,             ! in/out: solution
     &              g_grad_cmplx,             ! in    : b (of Ax=b)  
     &              uhf_hessv3_cmplx1,         ! in: product routine: gen Az
     &              uhf_hessv_precond3_damp,  ! in: routines
     &              precond_tol,
     &              max_precond_sub, ! jochen changed this from 11*nvec
     &              maxiter, 
     &              .false., 
     &              oprint_solve,
     &              omega, 
     &              limag,
     &              lifetime,  ! damp means complex, it is redundant
     &              gamwidth, 
     &              ncomp,    ! ncomp=2 (+/-)
     &              2,
     &              nvir1, 
     &              nocc1)


           call conv2reim(g_work,      ! out : real      arr
     &                    g_work_im,   ! out : imaginary arr
     &                    g_work_cmplx,! in  : = complx(g_xre,g_xim)
     &                    vlen,        ! in  : n    rows
     &                    nvec,        ! in  : nvec columns
     &                    ncomp)       ! in  : nr. components

           do ipm=1,ncomp
            if(.not.ga_destroy(g_work_cmplx(ipm))) call errquit(
     $     'cphf_solve3: could not destroy work_cmplx',ipm, GA_ERR)
            if(.not.ga_destroy(g_grad_cmplx(ipm))) call errquit(
     $     'cphf_solve3: could not destroy grad_cmplx',ipm, GA_ERR)
           enddo ! end-loop-ipm
 401       continue
           if (debug) write (6,*) 'exiting ga_lkain_2cpl3_damp'
         endif
         noskew_uhf = .true.

       else 

c        -----------------------
c        RHF section starts here
c        -----------------------
c
c        For property related CPHF. NMR has imaginary perturbation
c        which makes matrices skew symmetric. Setting noskew to .false.
c        enforces that.
         
c        ... jochen: our density matrix is neither symmetric nor
c        antisymmetric for omega.ne.0. the whole thing
c        will probably work with RHF only anyway so I didn't put a similar
c        comment up in the UHF section
c        
         if (.not. rtdb_get(rtdb, 'cphf:skew', mt_log, 1, 
     $      noskew)) noskew = .true.
c        the dynamic density matrix builder needs to know if the perturbation
c        is imaginary (or antisymmetric). This should have been written
c        to the rtdb at some time (e.g. in aoresponse)

         if (.not. rtdb_get(rtdb, 'cphf:imag', mt_log, 1,
     &      limag)) then
           limag = .false.
           if (ga_nodeid().eq.0) then
             write (luout,*)
     &        'WARNING from cphf_solve3: could not read limag from RTDB'
             write (luout,*)
     &        'Assuming REAL Symmetric perturbation.'
           endif
         endif

         if (debug) then
           write (luout,*) 'cphf_solve3: omega,imag,lifetime,gamwidth',
     &        omega, limag,lifetime,gamwidth,ncomp
         endif

c ... jochen: call precond3 instead of precond next line. Originally,
c        the second argument was set to zero. Now we call the pre-
c        conditioner twice, once with a shift of -omega, once with
c        +omega. That should yield the correct A(+/-) matrices
c ... jochen: added preconditioner that takes care of damping


          omg(1)=-omega
          omg(2)= omega
          gam(1) = -gamwidth
          gam(2) = gamwidth
          
          if (.not.lifetime) then
           do ipm=1,ncomp
            call rohf_hessv_precond3(g_work(ipm),omg(ipm))
           enddo ! end-loop-ip
          else
            do ipm=1,ncomp
              if (debug) write (luout,*) 'cphf_solve3: calling'//
     &          'precond3_damp with ipm, omega, gamma',
     &          ipm, omg(ipm), gam(ipm)
              call rohf_hessv_precond3_damp(
     &          g_work(ipm),g_work_im(ipm),omg(ipm),gam(ipm))
            enddo ! end-loop-ipm
          endif

        
c ... jochen: note: ga_lkain* is the equation solver. the
c        Fock matrix builder is in the rohf_hessv? piece, the
c        precoditioner is rohf_hessv_precon?. We should use
c        versions "3" here that goes along with cphf_solve3
c
c        in the following call, g_work is the initial guess, that
c        is why we have already divided by the energy denominators
c        in the preconditioner.
c        g_grad is the constant rhs vector
c
c       DIM/QM JEM
         if (ldimqm .and. lrsp) call dimqm_prep(rtdb, omega)

         if (.not.lifetime) then
           if (debug) write (6,*) 'calling ga_lkain_2cpl3'


c           call ga_lkain_2cpl3(rtdb,
           call ga_lkain_2cpl3_ext(rtdb,
     &                         g_work   , g_grad, 
     &                         g_work_im, g_grad_im, 
     &                         rohf_hessv3,rohf_hessv_precond3,
     $                         precond_tol,max_precond_sub,maxiter,
     &                         .false., oprint_solve, 
     &                         omega,
     &                         limag, lifetime, 
     &                         gamwidth, ncomp)

           if (debug) write (6,*) 'exiting ga_lkain_2cpl3'

         else

c --- Transfer (g_work,g_work_im) --> g_work_cmplx
c              (g_grad,g_grad_im) --> g_grad_cmplx
c ----- Create complex GA arrays ----------- START
           do ipm = 1,ncomp
            if(.not.ga_create(mt_dcpl,vlen,nvec,
     &         'cphf_solve3: work_cmplx',0,0,g_work_cmplx(ipm)))
     &      call errquit('cphf_solve3: g_work_cmplx',ipm,
     &                   GA_ERR)
            if(.not.ga_create(mt_dcpl,vlen,nvec,
     &          'cphf_solve3: work_cmplx',0,0,g_grad_cmplx(ipm)))
     &      call errquit('cphf_solve3: g_work_cmplx',ipm,
     &                   GA_ERR)
           enddo ! end-loop-ipm
c ----- Create complex GA arrays ----------- END
           call conv2complex(g_work_cmplx,! out: = complx(g_xre,g_xim)
     &                       g_work,      ! in : real      arr
     &                       g_work_im,   ! in : imaginary arr
     &                       vlen,        ! in : n    rows
     &                       nvec,       ! in : nvec columns
     &                       ncomp)       ! in : nr. components
           call conv2complex(g_grad_cmplx,! out: = complx(g_xre,g_xim)
     &                       g_grad,      ! in : real      arr
     &                       g_grad_im,   ! in : imaginary arr
     &                       vlen,        ! in : n    rows
     &                       nvec,       ! in : nvec columns
     &                       ncomp)       ! in : nr. components


           nocc1(1)=nclosed
           nvir1(1)=nmo-nocc1(1)

           call ga_lkain_2cpl3_damp_cmplx_redmem(
     &              rtdb,
     &              g_work_cmplx,             ! in/out: solution
     &              g_grad_cmplx,             ! in    : b (of Ax=b)  
     &              rohf_hessv3_cmplx,        ! in: product routine: gen Az
     &              rohf_hessv_precond3_damp, ! in: routines
     &              precond_tol,
     &              max_precond_sub, ! jochen: changed this from 11*nvec,
     &              maxiter, 
     &              .false., oprint_solve,
     &              omega, 
     &              limag,
     &              lifetime,  ! damp means complex, it is redundant
     &              gamwidth, 
     &              ncomp,    ! ncomp=2 (+/-)
     &              1,
     &              nvir1, 
     &              nocc1)

           call conv2reim(g_work,      ! out : real      arr
     &                    g_work_im,   ! out : imaginary arr
     &                    g_work_cmplx,! in  : = complx(g_xre,g_xim)
     &                    vlen,        ! in  : n    rows
     &                    nvec,        ! in  : nvec columns
     &                    ncomp)       ! in  : nr. components

           do ipm=1,ncomp
            if(.not.ga_destroy(g_work_cmplx(ipm))) call errquit(
     $     'cphf_solve3: could not destroy work_cmplx',ipm, GA_ERR)
            if(.not.ga_destroy(g_grad_cmplx(ipm))) call errquit(
     $     'cphf_solve3: could not destroy grad_cmplx',ipm, GA_ERR)
           enddo ! end-loop-ipm

           if (debug) write (6,*) 'exiting ga_lkain_2cpl3_damp'
         endif
         noskew = .true.
c
c ... jochen: note that on output the rhs is already divided by
c        the orbital energy diffs from preconditioner
c        calls inside the solver

      endif
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $   call errquit('cphf_solve3: rtdb_put of xc_active failed',0,
     &   RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $   call errquit('cphf_solve3: rtdb_put of calc_type failed',0,
     &   RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $   call errquit('cphf_solve3: rtdb_put of j_derfit failed',0,
     &   RTDB_ERR)
c ======== writing (g_b,g_rhs_sol) ======== START
      call cphf_fname('cphf_sol',cphf_sol(1))
      if(.not.file_write_ga(cphf_sol(1),g_work(1)))
     $     call errquit('cphf_solve3: could not write cphf_sol',1,
     &       UNKNOWN_ERR)
c ... jochen: also write second component of solution vector to file
      if (ncomp.gt.1) then
        call cphf_fname('cphf_sol_2',cphf_sol(2))
        if(.not.file_write_ga(cphf_sol(2),g_work(2)))
     $     call errquit('cphf_solve3: could not write cphf_sol',2,
     &     UNKNOWN_ERR)
      endif
c ... jochen: and write imaginary parts to file as well
      if (lifetime) then
        call cphf_fname('cphf_sol_im',cphf_sol_im(1))
        if(.not.file_write_ga(cphf_sol_im(1),g_work_im(1)))
     $     call errquit('cphf_solve3: could not write cphf_sol_im',1,
     &       UNKNOWN_ERR)
        if (ncomp.gt.1) then
          call cphf_fname('cphf_sol_im_2',cphf_sol_im(2))
          if(.not.file_write_ga(cphf_sol_im(2),g_work_im(2)))
     $       call errquit
     &       ('cphf_solve3: could not write cphf_sol_im_2',1,
     &       UNKNOWN_ERR)
        endif
      endif ! lifetime
c ======== writing (g_b,g_rhs_sol) ======== END
      do ipm = 1,ncomp
        if(.not.ga_destroy(g_work(ipm)))call errquit(
     $     'cphf_solve3: could not destroy work array',ipm, GA_ERR)
        if(.not.ga_destroy(g_grad(ipm)))call errquit(
     $     'cphf_solve3: could not destroy gradient handle',ipm, GA_ERR)
        if (lifetime) then
          if(.not.ga_destroy(g_work_im(ipm)))call errquit(
     $       'cphf_solve3: could not destroy Im work array',ipm, GA_ERR)
          if(.not.ga_destroy(g_grad_im(ipm)))call errquit(
     $       'cphf_solve3: could not destroy Im gradient handle',
     &       ipm, GA_ERR)
        endif                   ! lifetime
      enddo                     ! ipm

      call int_1e_uncache_ga()
      call schwarz_tidy()
      call int_terminate()
      if(scftype .eq. 'UHF') then
         call uhf_tidy(rtdb)
      else if(scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         call rohf_tidy
      endif
      call scf_tidy(rtdb)
      if (ldimqm .and. lrsp) then
        call dimqm_tidy()
      endif
      return
      end

      subroutine prepare_guess(
     &       g_work_re,   ! out: RE guess
     &       g_work_im,   ! out: IM guess
     &       g_grad_re,   ! in : RE (g_b,g_z)
     &       g_grad_im,   ! in : RE (g_b,g_z)
     &       vlen,        ! in : nrows guess
     &       nvec,        ! in : ncols guess
     &       ncomp,       ! in : nr. components
     &       lifetime,    ! in : = T complex   = F real
     &       use_my_guess,! in : = T guess-file =F guess=en-scl*g_b 
     &       rtdb)        ! in : rtdb handle
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: using a guess or not
c Date   : 03-15-12

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"   
#include "stdio.fh"
#include "rtdb.fh"
      integer rtdb,
     &        vlen,nvec,ncomp,ipm,
     &        m1,m2,shift,
     &        g_work_re(ncomp),
     &        g_work_im(ncomp),  
     &        g_grad_re(ncomp),
     &        g_grad_im(ncomp)
      logical use_my_guess,lifetime
c Note.- use_my_guess, stored in rtdb and used
c        to indicate if we are using a vector solution
c        taken from a file or not
c        = T, then it reads vector from a file in location
c        that calls cphf3 (e.g. fiao_f1_movecs)
c        after the vector is read it will be stored in
c        2nd subspace (n1 x nvec) of g_grad GA array which
c        is being read in cphf_solve3.F
c        This idea could be generalized to:
c        cphfX, cphf_solvX.F.
c        The use of reading a guess from file is to
c        get a vector guess obtained from a previous calculation
c        a) In order to finish the convergence, or
c        b) In order to start a new calc. on same molecular
c           system but with parameters (e.g. damping) slightly
c           changed.
c        = F, it does the standard way to guess:
c         ==> g_z = (energy scaling) g_b 
c             (for iteratively linear solving A z = b)
c        FA-05-6-12        
c        Note.- if using guess stored in a file
c               contained in g_grad 2nd subspace
c               subspace  (n1,i1) i1=4,5,6=x,y,z
c               n1= [sum_i nocc(i)*nvir(i)] * ncomp
         if (.not.(use_my_guess)) then ! if-use_my_guess--START
          m1=1
          m2=nvec
         else   ! use guess stored in g_grad 2nd subspace
           if (ga_nodeid().eq.0)
     &      write(*,*) 'cphf_solve3: Using guess from file ...'
          shift=3 ! skip first subspace that corresponds to g_b
          m1=shift+1
          m2=shift+nvec
         endif ! ------------------------if-use_my_guess--END
          do ipm = 1,ncomp
           call ga_copy_patch('n',g_grad_re(ipm),1,vlen,m1,m2, 
     $                            g_work_re(ipm),1,vlen,1 ,nvec)
           if (lifetime) then
           call ga_copy_patch('n',g_grad_im(ipm),1,vlen,m1,m2, 
     $                            g_work_im(ipm),1,vlen,1 ,nvec)
           endif
          enddo ! end-loop-ipm
      return
      end
