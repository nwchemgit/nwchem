      subroutine cphf_solve2(rtdb)
*
* $Id: cphf_solve2.F,v 1.4 2002-01-16 17:23:18 windus Exp $
*
      implicit none
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
c!!! BGJ
#include "bgj.fh"
c!!! BGJ
      integer rtdb              ! [input] database handle
      double precision eone,etwo,enrep,energy
      integer g_grad, g_work
      integer max_precond_iter,max_precond_sub
      double precision precond_acc
      integer precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
      logical oprint_solve
c$$$  logical ga_it_ls_mrhs
c$$$  external ga_it_ls_mrhs
      external rohf_hessv_precond,rohf_hessv, rohf_hessv2
      external uhf_precond,uhf_hessv, uhf_hessv2
      logical oprint, oprint_conv
c     
      character*255 rhs_title
      character*(nw_max_path_len) cphf_rhs
      character*(nw_max_path_len) cphf_sol
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
      logical hf_job
      integer nroot, leval, ieval, iroot
      logical cphf_tddft, cphf_triplet
      cphf_tddft = .false. 
      cphf_triplet = .false. ! true for triplet, false for singlet
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
c     Must use the same vectors as the converged SCF
c
      if (hf_job) then
         if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $        call errquit('cphf: no SCF MO vectors',0)
      else
         if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $        call errquit('cphf_solve: no DFT MO vectors',0)
c !!! Hack to set up scftype properly for later code !!! ? need this?
         if (bgj_restricted()) then
            scftype = 'RHF'
         else
            scftype = 'UHF'
         endif
      endif
c      write(*,*)'*** cphf_solve: scftype ',scftype
c     
      call scf_init(rtdb)
      if (.not.int_normalize(rtdb,basis))
     &    call errquit('cphf_solve: int_normalize failed',911)
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
     $     call errquit('cphf: problem with movecs_read_header',1)
      if(nbf.ne.nbf_vecs)call errquit('cphf: MO vectors mismatch?',0)
      nmo = nmos(1)
c
      if (nmo .ne. nbf) then
         if (.not. ga_destroy(g_movecs(1))) call errquit('cphf: ga',1)
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: tmp MOs a',
     $        0, 0, g_movecs(1))) call errquit('scf_v_g: tmp MOs', 0)
         if (scftype .eq. 'UHF') then
            if (.not. ga_destroy(g_movecs(2)))call errquit('cphf: ga',2)
*ga:1:0
            if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: tmp MOs b',
     $           0, 0, g_movecs(2))) call errquit('scf_v_g: tmp MOs', 0)
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
         call errquit('cphf: unknown SCF type',0)
      endif
c
      lshift=0.0d0
      if (.not. rtdb_get(rtdb, 'cphf:thresh', mt_dbl, 1, 
     $     precond_tol)) precond_tol=1d-4
      if (.not. rtdb_get(rtdb, 'cphf:maxiter', mt_int, 1,
     $     max_precond_iter)) max_precond_iter = 100
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)         
      if (.not. file_read_ga_info(cphf_rhs, rhs_title, ncol, nvec))
     $     call errquit('cphf_solve: reading rhs info failed?',0)
      if (vlen .ne. ncol) 
     $     call errquit('cphf_solve: vlen != ncol', ncol)
      if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve: rhs',
     $     0,nvec,g_grad))
     $     call errquit('cphf_solve: ga_cre grad',vlen)
c
      pflg=2
      max_precond_sub = 20
      if (oprint) then
         write(6,170) scftype
 170     format('    scftype = ',a)
         if (scftype .eq. 'RHF' .or. scftype.eq.'ROHF') then
            write(6,171) nclosed, nopen
 171        format('    nclosed = ', i4/
     $             '      nopen = ', i4)
         else
            write(6,172) nalpha, nbeta
 172        format('     nalpha = ', i4/
     $             '      nbeta = ', i4)
         endif
         write(6,173) vlen, nvec, precond_tol, max_precond_iter
 173     format('  variables = ', i8/
     $          '  num. vecs = ', i8/
     $          '  tolerance = ', 1p,d9.1,0p/
     $          '  max. iter = ', i4/)
      endif
c
*ga:1:0
c
      if(.not.file_read_ga(cphf_rhs, g_grad))
     $     call errquit('cphf_solve: could not read cphf_rhs',1)
c
*ga:1:0
      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve: grad',
     $     0,0,g_work))call errquit('cphf_solve: ga_cre grad',0)
c
c     Compute necessary Fock matrices
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 1))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0)
c !!! BGJ
      if(scftype .eq. 'UHF')then
         call uhf_init(rtdb)
         call uhf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
      else
         call rohf_init(rtdb)
         call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_work,
     ,        rtdb)
      endif
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',1)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
c !!! BGJ
c
      residual = sqrt(ga_ddot(g_work, g_work))
      if (residual .gt. 1d-4) then
         if (ga_nodeid() .eq. 0) then
            write(6,*) ' SCF residual = ', residual
            call util_flush(6)
         endif
         call ga_sync()
         call errquit('cphf_solve:the available MOs do not satisfy'//
     $        ' the SCF equations',0)
      endif
c
      if (.not. ga_destroy(g_work)) call errquit('cphf_solve:ga?',0)
      if(.not.ga_create(mt_dbl,vlen,nvec,'cphf_solve: work',
     $     0, 0, g_work))call errquit('cphf_solve: ga_cre work',0)
c
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .true.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0)
c !!! BGJ
      if (scftype .eq. 'UHF') then
c        if (nvec .ne. 1) call errquit('not yet',nvec)
         call ga_zero(g_work)   ! Essential to get default initial guess
c         if(.not.ga_iter_lsolve(vlen, max_precond_iter,
c     $        max_precond_sub, precond_tol, 
c     $        uhf_precond, uhf_hessv, .false., oprint_solve, 
c     $        g_grad, g_work, precond_acc, precond_iter))
c     $        call errquit('cphf did not converge',precond_iter)
         call ga_copy(g_grad, g_work)
         call uhf_precond(g_work,0.0d0)
         call ga_lkain(g_work, g_grad, uhf_hessv2, uhf_precond,
     $        precond_tol, 10*nvec, 20, .false., oprint_solve)
      else 
c$$$         if(.not.ga_iter_lsolve(vlen, max_precond_iter,
c$$$     $        max_precond_sub, precond_tol, 
c$$$     $        rohf_hessv_precond, rohf_hessv, .false., oprint_solve, 
c$$$     $        g_grad, g_work, precond_acc,precond_iter))
c$$$     $        call errquit('cphf did not converge', precond_iter)
         call ga_copy(g_grad, g_work)
         call rohf_hessv_precond(g_work,0.0d0)
         call ga_lkain(g_work, g_grad, rohf_hessv2, rohf_hessv_precond,
     $        precond_tol, 10*nvec, 20, .false., oprint_solve)
      endif
c !!! BGJ
      if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of xc_active failed',0)
      if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 0))
     $     call errquit('cphf_solve: rtdb_put of calc_type failed',0)
      if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('cphf_solve: rtdb_put of j_derfit failed',0)
c !!! BGJ
c
      if (oprint_conv) then
         write(6,174) precond_acc, precond_iter
 174     format(/' CPHF converged to ',1p,d9.1,0p,' in ', i4, ' iters'/)
         call util_flush(6)
      endif
c
c$$$  if(.not.ga_it_ls_mrhs(cuhf_vlen,1,max_iter,
c$$$  $        max_sub,accuracy, 
c$$$  $        uhf_precond,uhf_hessv,oprint, 
c$$$  $        g_grad,g_work,precond_acc,precond_iter))
c$$$  $        call errquit('cphf did not converge',1)
c$$$  pol=ga_ddot(g_work,g_grad)
c$$$  if(ga_nodeid().eq.0)then
c$$$  write(6,911) pol
c$$$  911        format(//,15x,40('-'),
c$$$  $           /,15x,'UHF Polarizability is:',f16.6,
c$$$  $           /,15x,40('-'),/)
c$$$  endif
c
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      if(.not.file_write_ga(cphf_sol,g_work))
     $     call errquit('cphf_solve: could not write cphf_sol',1)
      if(.not.ga_destroy(g_work))call errquit(
     $     'cphf_solve: could not destroy work array',1)
      if(.not.ga_destroy(g_grad))call errquit(
     $     'cphf_solve: could not destroy gradient handle',1)
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
