      subroutine cphf_solve3(rtdb, omega, lifetime, gamwidth)
c 
c jochen: modified version of cphf_solve2 to deal with 
c         frequency dependent properties. for this, array g_grad
c         will now be of dimension 2. g_work holds a copy of
c         that array and will also have dimension 2
c jochen: Further extended this subroutine to deal with
c         finite-lifetime damping. Arrays dealing with imaginary
c         parts have been added, indicated by suffix "_im"
*
* $Id: cphf_solve3.F,v 1.4 2009-02-19 22:13:05 niri Exp $
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
      integer max_precond_iter,max_precond_sub
      double precision precond_acc
      integer precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
      logical oprint_solve
c$$$  logical ga_it_ls_mrhs
c$$$  external ga_it_ls_mrhs
c ... jochen: modified the definition of some of the externals
c     here
      external rohf_hessv_precond3,rohf_hessv, rohf_hessv3
      external rohf_hessv_precond3_damp
      external uhf_precond,uhf_hessv, uhf_hessv2
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
      logical movecs_read_header, cphf_poliz
      external movecs_read_header, cphf_poliz
      logical int_normalize
      external int_normalize
c
      logical hf_job, debug, limag, lstatic
      integer ipm, ncomp
      double precision small
      parameter (small = 1d-8)
c
c     ==================================================================
c
      debug = (.false. .and. ga_nodeid().eq.0)
c
      lstatic = (abs(omega).lt.small)
c
      hf_job = .not. bgj_have_xc()      
c
c     Print levels assume that CPHF is being invoked from within
c     another application.  By default be informative.
c     At low print, just indicate convergence or failure.
c
      oprint = util_print('information',print_default) .and.
     $     ga_nodeid() .eq. 0
      oprint_conv = util_print('information',print_low) .and.
     $     ga_nodeid() .eq. 0
      oprint_solve = util_print('solve',print_default)
c
c     determine how many components we should treat:
c
      if (.not. rtdb_get(rtdb, 'cphf:ncomp', mt_int, 1, ncomp)) call
     &   errquit('cphf_solve3: failed to read ncomp ', 0, RTDB_ERR)
c
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
c
      if (hf_job) then
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $        call errquit('cphf: no SCF MO vectors',0, RTDB_ERR)
      else
         if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $        call errquit('cphf_solve3: no DFT MO vectors',0, RTDB_ERR)
c !!! Hack to set up scftype properly for later code !!! ? need this?
         if (bgj_restricted()) then
            scftype = 'RHF'
         else
            scftype = 'UHF'
         endif
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
        write(*,*) __FILE__,': bases',bases
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
      if(nbf.ne.nbf_vecs)call errquit('cphf: MO vectors mismatch?',0,
     &       UNKNOWN_ERR)
      nmo = nmos(1)
c
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
c
      if (scftype.eq.'UHF') then 
         vlen = nalpha*(nmo-nalpha) + nbeta*(nmo-nbeta)
      else if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         nvirt = nmo - nclosed - nopen
         vlen = nclosed*(nopen+nvirt) + nopen*nvirt
      else
         call errquit('cphf: unknown SCF type',0, INPUT_ERR)
      endif
c
      lshift=0.0d0
c
c note: the original default value for "tol" was 1d-4
      if (.not. rtdb_get(rtdb, 'cphf:thresh', mt_dbl, 1, 
     $     precond_tol)) precond_tol=1d-4
      if (.not. rtdb_get(rtdb, 'cphf:maxiter', mt_int, 1,
     $     max_precond_iter)) max_precond_iter = 20
c
c ... jochen: open second rhs file first and check if it has data
      if (ncomp.gt.1) then
        call util_file_name('cphf_rhs_2',.true.,.true.,cphf_rhs(2)) 
        if (.not. file_read_ga_info(cphf_rhs(2), rhs_title, ncol, nvec))
     $     call errquit('cphf_solve3: reading rhs 2 info failed?',0,
     &     UNKNOWN_ERR)
      endif
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs(1))         
      if (.not. file_read_ga_info(cphf_rhs(1), rhs_title, ncol, nvec))
     $     call errquit('cphf_solve3: reading rhs info failed?',0,
     &       UNKNOWN_ERR)
      if (vlen .ne. ncol) 
     $     call errquit('cphf_solve3: vlen != ncol', ncol,
     &       UNKNOWN_ERR)

c ... jochen: open data file(s) with "imaginary" data
      if (lifetime) then
        if (ncomp.gt.1) then
          call util_file_name('cphf_rhs_im_2',.true.,.true.,
     +       cphf_rhs_im(2)) 
          if (.not. file_read_ga_info(cphf_rhs_im(2), rhs_title,
     +       ncol, nvec))
     $       call errquit
     &       ('cphf_solve3: reading rhs_im 2 info failed?',0,
     &       UNKNOWN_ERR)
        endif
        call util_file_name('cphf_rhs_im',.true.,.true.,
     +     cphf_rhs_im(1)) 
        if (.not. file_read_ga_info(cphf_rhs_im(1), rhs_title,
     +     ncol, nvec))
     $     call errquit('cphf_solve3: reading rhs_im 1 info failed?',0,
     &     UNKNOWN_ERR)        
      endif                     ! lifetime
c
c     allocate g_grad: rhs vector for CPKS
c
      if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: rhs',
     $     0,nvec,g_grad(1)))
     $     call errquit('cphf_solve3: ga_cre grad',vlen, GA_ERR)
c
c ... jochen: allocate second component of g_grad, and imaginary
c     components, if necessary for a computation with damping

      if (ncomp.gt.1) then
        if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: rhs2',
     $   0,nvec,g_grad(2)))
     $   call errquit('cphf_solve3: ga_cre grad 2',vlen, GA_ERR)    
      endif
      if (lifetime) then
        if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: rhs_im1',
     $     0,nvec,g_grad_im(1)))
     $     call errquit('cphf_solve3: ga_cre grad_im 1',vlen, GA_ERR)  
        if (ncomp.gt.1) then
          if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve3: rhs_im2',
     $       0,nvec,g_grad_im(2)))
     $       call errquit('cphf_solve3: ga_cre grad_im 2',vlen, GA_ERR) 
        endif
      endif                     ! lifetime
c
c
      pflg=2
      max_precond_sub = 20
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
c ... jochen here we read the CPKS rhs vector generated elsewhere
c     (A- or U-matrices, not yet divided by orbital energy diffs)
c
      if(.not.file_read_ga(cphf_rhs(1), g_grad(1)))
     $     call errquit('cphf_solve3: could not read cphf_rhs 1',1,
     &       UNKNOWN_ERR)
      if (ncomp.gt.1) then
        if(.not.file_read_ga(cphf_rhs(2), g_grad(2)))
     $     call errquit('cphf_solve3: could not read cphf_rhs 2',1,
     &     UNKNOWN_ERR)
      endif

c ... jochen: also read imaginary components if necessary
c             (should be zero before the preconditioner call)
      if (lifetime) then
        if(.not.file_read_ga(cphf_rhs_im(1), g_grad_im(1)))
     $     call errquit('cphf_solve3: could not read cphf_rhs 1',1,
     &     UNKNOWN_ERR)       
        if (ncomp.gt.1) then
          if(.not.file_read_ga(cphf_rhs_im(2), g_grad_im(2)))
     $       call errquit('cphf_solve3: could not read cphf_rhs 2',1,
     &       UNKNOWN_ERR)
        endif
      endif                     ! lifetime
c     
*     ga:1:0
      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve3: work1',
     $   0,0,g_work(1)))
     &   call errquit('cphf_solve3: ga_cre work1',0, GA_ERR)
c
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
c
      if (.not. rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo))
     $     ocosmo = .false.
      if (ocosmo) then
         cosmo_on = .true.
         cosmo_phase = 2
      endif
c
      if(scftype .eq. 'UHF')then
        call uhf_init(rtdb)
        call uhf_energy(g_movecs,eone,etwo,enrep,energy,g_work(1))
      else
        call rohf_init(rtdb)
        call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_work(1),
     ,     rtdb)
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
      residual = sqrt(ga_ddot(g_work(1), g_work(1)))
      if (residual .gt. 1d-4) then
        if (ga_nodeid() .eq. 0) then
          write(luout,*) ' WARNING: SCF residual = ', residual
          call util_flush(6)
        endif
        if(residual .gt.1d-2) then
          call ga_sync()
          call errquit('cphf_solve3:the available MOs do not satisfy'//
     $       ' the SCF equations',0, UNKNOWN_ERR)
        endif
      endif
c     
c     destroy temp work array used for Fock matrix calculation
      if (.not. ga_destroy(g_work(1))) call errquit('cphf_solve3:ga?',0,
     &   GA_ERR)
c     
c ... jochen: allocate ncomp components for array g_work here and take
c             care of arrays needed in a damping calculation, too
c
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
      end do
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
      if (scftype .eq. 'UHF') then
c
c       *******************************
c       UHF section, NOT used presently
c       (rhf section is below)
c       *******************************

         call errquit ('cphf_solv3 open shell not implemented',
     &     0, UNKNOWN_ERR)
c        THIS PIECE OF CODE SHOULD NOT BE EXECUTED IN RESPONSE CALC
C        BECAUSE THIS HAS NOT BEEN ADAPTED TO FREQUENCY AND DAMPING
c
c        For property related CPHF. NMR has imaginary perturbation
c        which makes matrices skew symmetric. Setting noskew to .false.
c        enforces that.
c
         if (.not. rtdb_get(rtdb, 'cphf:skew', mt_log, 1, 
     $        noskew_uhf)) noskew_uhf = .true.
c
c        if (nvec .ne. 1) call errquit('not yet',nvec)
         call ga_zero(g_work(1))   ! Essential to get default initial guess
c         if(.not.ga_iter_lsolve(vlen, max_precond_iter,
c     $        max_precond_sub, precond_tol, 
c     $        uhf_precond, uhf_hessv, .false., oprint_solve, 
c     $        g_grad, g_work, precond_acc, precond_iter))
c     $        call errquit('cphf did not converge',precond_iter)
         call ga_copy(g_grad(1), g_work(1))
         call uhf_precond(g_work(1),0.0d0)
         call ga_lkain(rtdb,g_work(1),g_grad(1),uhf_hessv2,uhf_precond,
     $      precond_tol, 10*nvec, 20, .false., oprint_solve)
         noskew_uhf = .true.
       else 
c
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
         do ipm = 1,ncomp
           call ga_copy(g_grad(ipm), g_work(ipm))
           if (lifetime) then
             call ga_copy(g_grad_im(ipm), g_work_im(ipm))
           endif
         enddo

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

         if (debug) then
           write (6,*) 'cphf_solve3: before preconditioner:'
           do ipm = 1,ncomp
             call ga_print(g_work(ipm))
             if (lifetime) then
               call ga_print(g_work_im(ipm))
             endif
           enddo
         endif

         if (.not.lifetime) then
           call rohf_hessv_precond3(g_work(1),-omega)
           if (ncomp.gt.1) call rohf_hessv_precond3(g_work(2), omega)
           if (debug) then
             write (6,*) 'cphf_solve3: after preconditioner:'
             call ga_print(g_work(1))
             if (ncomp.gt.1) call ga_print(g_work(2))
           endif
         else
           call rohf_hessv_precond3_damp(
     &        g_work(1),g_work_im(1),-omega, gamwidth)
           if (ncomp.gt.1) call rohf_hessv_precond3_damp(
     &        g_work(2),g_work_im(2), omega, gamwidth)
           if (debug) then
             write (6,*) 'cphf_solve3: after preconditioner:'
             call ga_print(g_work(1))
             if (ncomp.gt.1) call ga_print(g_work(2))
             call ga_print(g_work_im(1))
             if (ncomp.gt.1) call ga_print(g_work_im(2))
           endif
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

         if (.not.lifetime) then
           if (debug) write (6,*) 'calling ga_lkain_2cpl3'
           call ga_lkain_2cpl3(rtdb,
     &        g_work, g_grad, g_work_im, g_grad_im, rohf_hessv3,
     &        rohf_hessv_precond3,
     $        precond_tol,10*nvec,maxiter,.false., oprint_solve, omega,
     &        limag, lifetime, gamwidth, ncomp)
           if (debug) write (6,*) 'exiting ga_lkain_2cpl3'
         else
           if (debug) write (6,*) 'calling ga_lkain_2cpl3_damp'
           call ga_lkain_2cpl3_damp(rtdb,
     &        g_work, g_grad, g_work_im, g_grad_im, rohf_hessv3,
     &        rohf_hessv_precond3_damp,
     $        precond_tol,10*nvec,maxiter,.false., oprint_solve, omega,
     &        limag, lifetime, gamwidth, ncomp)
           if (debug) write (6,*) 'exiting ga_lkain_2cpl3_damp'
         endif
         noskew = .true.
c
c ... jochen: note that on output the rhs is already divided by
c        the orbital energy diffs from preconditioner
c        calls inside the solver
c
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
c !!! BGJ
c
c     if (oprint_conv) then
c        write(luout,174) precond_acc, precond_iter
c174     format(/' CPHF converged to ',1p,d9.1,0p,' in ',i4,' iters'/)
c        call util_flush(6)
c     endif
c
c$$$  if(.not.ga_it_ls_mrhs(cuhf_vlen,1,max_iter,
c$$$  $        max_sub,accuracy, 
c$$$  $        uhf_precond,uhf_hessv,oprint, 
c$$$  $        g_grad,g_work,precond_acc,precond_iter))
c$$$  $        call errquit('cphf did not converge',1)
c$$$  pol=ga_ddot(g_work,g_grad)
c$$$  if(ga_nodeid().eq.0)then
c$$$  write(luout,911) pol
c$$$  911        format(//,15x,40('-'),
c$$$  $           /,15x,'UHF Polarizability is:',f16.6,
c$$$  $           /,15x,40('-'),/)
c$$$  endif
c
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol(1))
      if(.not.file_write_ga(cphf_sol(1),g_work(1)))
     $     call errquit('cphf_solve3: could not write cphf_sol',1,
     &       UNKNOWN_ERR)
c ... jochen: also write second component of solution vector to file
      if (ncomp.gt.1) then
        call util_file_name('cphf_sol_2',.true.,.true.,cphf_sol(2))
        if(.not.file_write_ga(cphf_sol(2),g_work(2)))
     $     call errquit('cphf_solve3: could not write cphf_sol',2,
     &     UNKNOWN_ERR)
      endif
c ... jochen: and write imaginary parts to file as well
      if (lifetime) then
        call util_file_name('cphf_sol_im',.true.,.true.,
     &     cphf_sol_im(1))
        if(.not.file_write_ga(cphf_sol_im(1),g_work_im(1)))
     $     call errquit('cphf_solve3: could not write cphf_sol_im',1,
     &       UNKNOWN_ERR)
        if (ncomp.gt.1) then
          call util_file_name('cphf_sol_im_2',.true.,.true.,
     &       cphf_sol_im(2))
          if(.not.file_write_ga(cphf_sol_im(2),g_work_im(2)))
     $       call errquit
     &       ('cphf_solve3: could not write cphf_sol_im_2',1,
     &       UNKNOWN_ERR)
        endif
      endif ! lifetime

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
c
      call int_1e_uncache_ga()
      call schwarz_tidy()
      call int_terminate()
      if(scftype .eq. 'UHF') then
         call uhf_tidy(rtdb)
      else if(scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
         call rohf_tidy
      endif
      call scf_tidy(rtdb)
c
      return
c
c     if (.not. cphf_poliz(rtdb))     ! Never executed.
c    $     call errquit(' cphf: not possible!',0)
c
      end
c
