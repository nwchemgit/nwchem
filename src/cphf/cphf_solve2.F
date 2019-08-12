      subroutine cphf_solve2(rtdb)
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
      double precision eone,etwo,enrep,energy
      integer g_grad, g_work,
     &        g_work_im,g_grad_im ! dummy
      integer max_precond_iter,max_precond_sub
      double precision precond_acc
      integer precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
      logical oprint_solve
c$$$  logical ga_it_ls_mrhs
c$$$  external ga_it_ls_mrhs
      logical use_my_guess
      integer m1,m2,shift,nvec1
      external rohf_hessv_precond,rohf_hessv, rohf_hessv2
      external uhf_precond,uhf_hessv, uhf_hessv2,
     &         prepare_guess
      logical oprint, oprint_conv, ocosmo
c     
      character*255 rhs_title
      character*(nw_max_path_len) cphf_rhs
      character*(nw_max_path_len) cphf_sol
      logical file_write_ga,file_read_ga, file_read_ga_info
      external file_write_ga,file_read_ga, file_read_ga_info
      integer vlen, nvirt, nmos(2), nbf_vecs, nsets
      integer bases(2), nbases, nvec, ncol,ndata
      character*255 title_vecs, basis_vecs
      character*20 scftype_vecs
      double precision residual
      double precision ecosmo
      logical movecs_read_header, cphf_poliz
      external movecs_read_header, cphf_poliz
      logical int_normalize
      external int_normalize
      external ga_lkain_ext
     
c
      logical use_dimqm,lrsp,ldimqm,lfirst,dimqm_go
      logical do_cphf
      logical hf_job
      logical xc_gotxc,xc_gothfx
      external xc_gotxc,xc_gothfx
c
      hf_job = (.not. xc_gotxc()).and.(.not.xc_gothfx())
      use_dimqm=util_module_avail("dimqm")
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
c     Must use the same vectors as the converged SCF
c
      if (hf_job) then
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $        call errquit('cphf: no SCF MO vectors',0, RTDB_ERR)
      else
         if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $        call errquit('cphf_solve: no DFT MO vectors',0, RTDB_ERR)
      endif
c      write(*,*)'*** cphf_solve: scftype ',scftype
c     
      call scf_init(rtdb)
      if (.not.int_normalize(rtdb,basis))
     &    call errquit('cphf_solve: int_normalize failed',911,
     &       UNKNOWN_ERR)
      nbases = 1
      bases(1) = basis
#if 0
      if (bgj_have_j_fit()) then
        nbases = 2
        bases(2) = bgj_CD_bas_han()
        if (ga_nodeid().eq.0) then
          write(LuOut,*) 'cphf_solve2: bases',bases
        endif
      endif
#endif
      call int_init(rtdb,nbases,bases)
      call schwarz_init(geom, basis)
      oskel=sym_number_ops(geom).gt.0
      call int_1e_cache_ga(basis,oskel)
c
c     Get the MO vectors

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
      if (.not. rtdb_get(rtdb, 'cphf:lshift', mt_dbl, 1,
     $     lshift)) lshift=0.0
c jochen: previously it tried to read cphf:precond_tol
c         I changed that so that we can set this by input
      if (.not. rtdb_get(rtdb, 'cphf:thresh', mt_dbl, 1, 
     $     precond_tol)) precond_tol=1d-4

      if (.not. rtdb_get(rtdb,'cphf:maxiter', mt_int, 1,
     $     maxiter)) maxiter = 50
      call cphf_fname('cphf_rhs',cphf_rhs)        
      if (.not. file_read_ga_info(cphf_rhs, rhs_title, ncol, nvec1))
     $     call errquit('cphf_solve: reading rhs info failed?',0,
     &       UNKNOWN_ERR)
      if (vlen .ne. ncol) 
     $     call errquit('cphf_solve: vlen != ncol', ncol,
     &       UNKNOWN_ERR)

      if (.not. rtdb_get(rtdb, 'cphf2-aores:ndata',mt_int,1, 
     $     ndata)) ndata=1
c      if (ga_nodeid().eq.0) then
c       write(*,10) ndata
c  10   format('Check ndata=',i3)
c      endif
      if (ndata.eq.2) then
       nvec=nvec1/2 ! real nvec=3
       if(.not.ga_create(mt_dbl,vlen,nvec1,'cphf_solve: rhs',
     $      0,nvec1,g_grad))
     $      call errquit('cphf_solve: ga_cre grad',vlen, GA_ERR)
      else
       nvec=nvec1
       if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve: rhs',
     $      0,nvec,g_grad))
     $      call errquit('cphf_solve: ga_cre grad',vlen, GA_ERR)
      endif

      max_precond_sub = 10*nvec
      if (.not. rtdb_get(rtdb, 'cphf:maxsub', mt_int, 1,
     $     max_precond_sub)) max_precond_sub = 10*nvec
c
      if (.not. rtdb_get(rtdb, 'cphf:do_cphf', mt_log, 1,
     $     do_cphf)) do_cphf=.true.
c
      pflg=2
      if (oprint) then
            write(luout,170) scftype
 170        format('  scftype          = ',4x,a)
         if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
            write(luout,171) nclosed, nopen
 171        format('  nclosed          = ', i8/
     $             '  nopen            = ', i8)
         else
            write(luout,172) nalpha, nbeta
 172        format('  nalpha           = ', i8/
     $             '  nbeta            = ', i8)
         endif
            write(luout,173) vlen,nvec,precond_tol,lshift,maxiter,
     &                        max_precond_sub
 173        format('  variables        = ', i8/
     $             '  # of vectors     = ', i8/
     $             '  tolerance        = ', d8.2/
     $             '  level shift      = ', d8.2/
     $             '  max iterations   = ', i8/
     $             '  max subspace     = ', i8/)
      endif
c
*ga:1:0
      if(.not.file_read_ga(cphf_rhs, g_grad))
     $     call errquit('cphf_solve: could not read cphf_rhs',1,
     &       UNKNOWN_ERR)

*ga:1:0
      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve: grad',
     $     0,0,g_work))call errquit('cphf_solve: ga_cre grad',0, GA_ERR)
c
c     Compute necessary Fock matrices
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 1))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
c !!! BGJ
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
     $                   g_work)
      else
         call rohf_init(rtdb)
         call rohf_energy(g_movecs,eone,etwo,enrep,ecosmo,energy,g_work,
     ,        rtdb)
      endif
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',1,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
c !!! BGJ
c
c     == check the residual ==
      residual = sqrt(ga_ddot(g_work, g_work))
      if(.not.rtdb_put(rtdb,'cphf_solve:residual',mt_dbl,1,residual))
     $   call errquit('cphf_solve: rtdb_put of residual failed',0,
     &     RTDB_ERR)
c
      if (ga_nodeid().eq.0) then  
       write(luout,*) 'SCF residual: ', residual
      endif
      if(residual .gt.1d-2) then
         call ga_sync()
         call errquit('cphf_solve2:SCF residual greater than 1d-2',
     &              0,UNKNOWN_ERR)
      endif

      if (.not. ga_destroy(g_work)) call errquit('cphf_solve:ga?',0,
     &       GA_ERR)

      if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve: work',
     $     0, 0, g_work))call errquit('cphf_solve: ga_cre work',0,
     &       GA_ERR)
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
c !!! BGJ
c 000000000000000 prepare guess 000000000000 START
c HvD: prepare_guess must always be called as it initializes g_work
      if (.not. rtdb_get(rtdb, 'cphf2-aores:guess', 
     &    mt_log,1,use_my_guess)) use_my_guess = .false.
      call prepare_guess(
     &               g_work,      ! out: RE guess
     &               g_work_im,   ! out: dummy
     &               g_grad,      ! in : RE (g_b,g_z)
     &               g_grad_im,   ! in : dummy
     &               vlen,        ! in : nrows guess
     &               nvec,        ! in : ncols guess
     &               1,           ! in : nr. components
     &               .false.,     ! in : lifetime= T complex   = F real
     &               use_my_guess,! in : = T guess-file =F guess=(energy-scaling)*g_b 
     &               rtdb)        ! in : rtdb handle
c Note.- en-scl*g_b is unpropertly called "preconditioning"
c 000000000000000 prepare guess 000000000000 END

c     DIM/QM JEM
c     Prep DIM/QM for response      
      if (ldimqm .and. lrsp) then
        call dimqm_prep(rtdb, -1d0)
      end if

      if (scftype .eq. 'UHF') then
c
c        For property related CPHF. NMR has imaginary perturbation
c        which makes matrices skew symmetric. Setting noskew to .false.
c        enforces that.
         if (.not. rtdb_get(rtdb, 'cphf:skew', mt_log, 1, 
     $        noskew_uhf)) noskew_uhf = .true.
         if (.not.(use_my_guess)) then ! doing default guess
          call uhf_precond(g_work,0.0d0)
         endif 
c        Call the solver
         if (do_cphf) then
         if (ndata.eq.2) then

          call ga_lkain_ext(rtdb,g_work,g_grad,uhf_hessv2,uhf_precond,
     $        precond_tol, max_precond_sub, maxiter, .false., 
     &        oprint_solve)
         else 

          call ga_lkain(rtdb,g_work, g_grad, uhf_hessv2, uhf_precond,
     $        precond_tol, max_precond_sub, maxiter, .false., 
     &        oprint_solve)
         endif
      else
         if(ga_nodeid().eq.0) write(luout,*) ' !!! Skipping CPHF !!! '
      endif
         noskew_uhf = .true.
      else 
c        For property related CPHF. NMR has imaginary perturbation
c        which makes matrices skew symmetric. Setting noskew to .false.
c        enforces that.
         if (.not. rtdb_get(rtdb, 'cphf:skew', mt_log, 1, 
     $        noskew)) noskew = .true.
         if (.not.(use_my_guess)) then ! doing default guess
          call rohf_hessv_precond(g_work,0.0d0)
         endif
c        Call the solver
         if (do_cphf) then
         if (ndata.eq.2) then

          call ga_lkain_ext(rtdb,g_work, g_grad, rohf_hessv2,
     R         rohf_hessv_precond,
     $         precond_tol, max_precond_sub, maxiter, .false., 
     &         oprint_solve)
         else

          call ga_lkain(rtdb,g_work, g_grad, rohf_hessv2,
     R         rohf_hessv_precond,
     $         precond_tol, max_precond_sub, maxiter, .false., 
     &         oprint_solve)
         endif
         else
         if(ga_nodeid().eq.0) write(luout,*) ' !!! Skipping CPHF !!! '
         endif
         noskew = .true.
      endif
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
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
      call cphf_fname('cphf_sol',cphf_sol)
       if (ndata.eq.2) then
c == copy solution (g_work) --> g_grad 2nd subspace===START
          shift=3 ! skip first subspace that corresponds to g_b
          m1=shift+1
          m2=shift+nvec
          call ga_copy_patch('n',g_work,1,vlen,1 ,nvec,
     &                           g_grad,1,vlen,m1,m2)
c == copy solution (g_work) --> g_grad 2nd subspace===END
        if(.not.file_write_ga(cphf_sol,g_grad))
     $     call errquit('cphf_solve: could not write cphf_sol',1,
     &       UNKNOWN_ERR)
       else
        if(.not.file_write_ga(cphf_sol,g_work))
     $     call errquit('cphf_solve: could not write cphf_sol',1,
     &       UNKNOWN_ERR)
       endif
      if(.not.ga_destroy(g_work))call errquit(
     $     'cphf_solve: could not destroy work array',1, GA_ERR)
      if(.not.ga_destroy(g_grad))call errquit(
     $     'cphf_solve: could not destroy gradient handle',1, GA_ERR)

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
