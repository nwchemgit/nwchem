      subroutine cphf_solve(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
#include "crohf.fh"
#include "cuhf.fh"
      integer rtdb              ! [input] database handle
      double precision eone,etwo,enrep,energy,pol
      integer g_grad,g_work,g_gradrhf,g_workrhf
      integer g_tmp,g_scratch
      integer max_precond_iter,max_precond_sub
      double precision precond_acc
      integer precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
      logical oprint_solve
c$$$  logical ga_it_ls_mrhs
c$$$  external ga_it_ls_mrhs
      external rohf_hessv_precond,rohf_hessv
      external uhf_precond,uhf_hessv
      logical oprint,do_rohf
c     
      character*256 cphf_rhs
      character*256 cphf_sol
      integer cphf_sol_u
      parameter(cphf_sol_u=18)
      logical file_write_ga,file_read_ga
      external file_write_ga,file_read_ga
      integer ierr
      integer k_buf,l_buf
c     
c     
c$$$  if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype))
c$$$  $     call errquit('scf:scftype wrong',1)
c$$$  if(scftype.eq.'UHF')do_rohf=.false.
c$$$  if(scftype.eq.'ROHF')do_rohf=.true.
c$$$  if(scftype.eq.'RHF')do_rohf=.true.
c$$$      call scf_init(rtdb)
c$$$  if(do_rohf)then
c$$$  call scf_pstat_init(rtdb)
c$$$  call int_normalize(basis)
c$$$  call int_init(rtdb,1,basis)
c$$$  call schwarz_init(geom, basis)
c$$$  call int_1e_cache_ga(basis,oskel)
c$$$  call scf_vectors_guess(rtdb)
c$$$  call rohf_init(rtdb)
c$$$  if(.not.ga_create(mt_dbl,crohf_vlen,1,'rohf_nr_solve: grad',
c$$$  $        1,1,g_grad)) call errquit('rohf_nr_solve: ga_cre grad',0)
c$$$  if(nopen.gt.0)then
c$$$  call ga_copy(crohf_g_fcv,crohf_g_fcp)
c$$$  call ga_copy(crohf_g_fcv,crohf_g_fpv)
c$$$  endif
c$$$  call rohf_get_grad(g_grad)
c$$$  lshift=0.0d0
c$$$  accuracy=1d-3
c$$$  pflg=2
c$$$  max_iter=20
c$$$  max_sub=10
c$$$  oprint=.false.
c$$$  if(.not.ga_create(mt_dbl,crohf_vlen,1,'rohf_nr_solve: grad',
c$$$  $        1,1,g_work)) call errquit('rohf_nr_solve: ga_cre grad',0)
c$$$  call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
c$$$  oskel=.false.
c$$$  call ga_zero(g_work)
c$$$  if(.not.ga_iter_lsolve(crohf_vlen,max_iter,
c$$$  $        max_sub,accuracy, 
c$$$  $        rohf_hessv_precond,rohf_hessv,oprint, 
c$$$  $        g_grad,g_work,precond_acc,precond_iter))
c$$$  $        call errquit('cphf did not converge',1)
c$$$  c$$$         if(.not.ga_it_ls_mrhs(crohf_vlen,1,max_iter,
c$$$  c$$$     $        max_sub,accuracy, 
c$$$  c$$$     $        rohf_hessv_precond,rohf_hessv,oprint, 
c$$$  c$$$     $        g_grad,g_work,precond_acc,precond_iter))
c$$$  c$$$     $        call errquit('cphf did not converge',1)
c$$$  c$$$         pol=ga_ddot(g_work,g_grad)
c$$$  c$$$         if(ga_nodeid().eq.0)then
c$$$  c$$$            write(6,999) pol
c$$$  c$$$ 999        format(//,15x,40('-'),
c$$$  c$$$     $           /,15x,'RHF/ROHF Polarizability is:',f16.6,
c$$$  c$$$     $           /,15x,40('-'),/)
c$$$  c$$$         endif
c$$$  call int_1e_uncache_ga()
c$$$  call schwarz_tidy()
c$$$  call int_terminate()
c$$$  call rohf_tidy
c$$$  call scf_tidy (rtdb)
c$$$  if(.not.ga_destroy(g_grad))call errquit(
c$$$  $        'could not destroy gradient handle',1)
c$$$  if(.not.ga_destroy(g_work))call errquit(
c$$$  $        'could not destroy work array',1)
c$$$  else
      call scf_init(rtdb)
      call int_normalize(basis)
      call int_init(rtdb,1,basis)
      call schwarz_init(geom, basis)
      call int_1e_cache_ga(basis,oskel)
      call scf_vectors_guess(rtdb)
c      if(scftype.eq.'UHF')then
      call uhf_init(rtdb)
      if(.not.cuhf_init_flag)
     $     call errquit(
     $     'uhf_hessv: uhf internal block invalid',0)
      if(.not.ga_create(mt_dbl,cuhf_vlen,1,'uhf_nr_solve: grad',
     $     1,1,g_grad))call errquit(
     $     'uhf_nr_solve: ga_cre grad',0)
c      endif
c      if(scftype.eq.'RHF')then
      call rohf_init(rtdb)
c$$$         if(.not.ga_create(mt_dbl,crohf_vlen*2,1,
c$$$     $        'rohf_nr_solve:grad',1,1,g_grad))
c$$$     $        call errquit(
c$$$     $        'could not alloc g_grad for rhf in cphf',0)
      if(.not.ga_create(mt_dbl,crohf_vlen,1,
     $     'rohf_nr_solve:grad',1,1,g_gradrhf))
     $     call errquit(
     $     'could not alloc g_gradrhf in cphf',0)
c      endif
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)         
      if(.not.file_read_ga(cphf_rhs,g_grad))
     $     call errquit('could not read cphf_rhs',1)
      if(scftype.eq.'RHF')then
         if(.not.ma_push_get(mt_dbl,crohf_vlen,'buf',l_buf,k_buf))
     $        call errquit('could not alloc buf',1)
         call ga_get(g_grad,1,crohf_vlen,1,1,dbl_mb(k_buf),1)
         call ga_put(g_gradrhf,1,crohf_vlen,1,1,dbl_mb(k_buf),1)
c$$$         if(.not.ma_pop_stack(l_buf))call errquit(
c$$$     $        'could not pop stack for l_buf',1)
      endif
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      call waopen(cphf_sol_u,cphf_sol,1,0,ierr)
      if(ierr.ne.0)call errquit(
     $     'could not open cphf_sol file',1)
      lshift=0.0d0
      precond_tol=1d-6
      pflg=2
      max_precond_iter=20
      max_precond_sub=10
      oprint=.true.
      if(.not.ga_create(mt_dbl,cuhf_vlen,1,'uhf_nr_solve: grad',
     $     1,1,g_work))call errquit('uhf_nr_solve: ga_cre grad',0)
      if(scftype.eq.'UHF')then
         call uhf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
         oskel=.false.
         call ga_zero(g_work)
         if(.not.ga_iter_lsolve(cuhf_vlen,max_precond_iter,
     $        max_precond_sub,precond_tol, 
     $        uhf_precond,uhf_hessv,oprint_solve, 
     $        g_grad,g_work,precond_acc,precond_iter))
     $        call errquit('cphf did not converge',1)
      endif
      if(scftype.eq.'RHF')then
         if(.not.ga_create(mt_dbl,crohf_vlen,1,
     $        'rohf_nr_solve: grad',1,1,g_workrhf))
     $        call errquit('uhf_nr_solve: ga_cre grad',0)
         call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_workrhf)
         oskel=.false.
         call ga_zero(g_work)
         if(.not.ga_iter_lsolve(crohf_vlen,max_precond_iter,
     $        max_precond_sub,precond_tol, 
     $        rohf_hessv_precond,rohf_hessv,oprint_solve, 
     $        g_gradrhf,g_workrhf,precond_acc,precond_iter))
     $        call errquit('cphf did not converge',1)
      endif
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
      if(scftype.eq.'RHF')then
         call ga_dscal(g_workrhf,2.0d0)
         call ga_get(g_workrhf,1,crohf_vlen,1,1,dbl_mb(k_buf),1)
         call ga_put(g_work,1,crohf_vlen,1,1,dbl_mb(k_buf),1)
         call ga_put(g_work,crohf_vlen+1,cuhf_vlen,1,1,dbl_mb(k_buf),1)
      endif
      if(.not.file_write_ga(cphf_sol,g_work))
     $     call errquit('could not write cphf_sol',1)
      if(.not.ga_destroy(g_work))call errquit(
     $     'could not destroy work array',1)
c$$$      if(.not.ga_destroy(g_workrhf))call errquit(
c$$$     $     'could not destroy work array',1)
      call int_1e_uncache_ga()
      call schwarz_tidy()
      call int_terminate()
c      if(scftype.eq.'UHF')
c      call uhf_tidy(rtdb)
c      if(scftype.eq.'RHF')
c      call rohf_tidy
c      call scf_tidy(rtdb)
      if(.not.ga_destroy(g_grad))call errquit(
     $     'could not destroy gradient handle',1)
      if(.not.ga_destroy(g_gradrhf))call errquit(
     $     'could not destroy gradient handle',1)
c$$$  endif
      end
