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
      integer g_grad,g_work
      integer g_tmp,g_scratch
      integer max_iter,max_sub
      double precision accuracy,precond_acc
      integer precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
c$$$      logical ga_it_ls_mrhs
c$$$      external ga_it_ls_mrhs
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
c
c
      if(.not.rtdb_cget(rtdb,'scf:scftype',1,scftype))
     $     call errquit('scf:scftype wrong',1)
      if(scftype.eq.'UHF')do_rohf=.false.
      if(scftype.eq.'ROHF')do_rohf=.true.
      if(scftype.eq.'RHF')do_rohf=.true.
      call scf_init(rtdb)
      if(do_rohf)then
         call scf_pstat_init(rtdb)
         call int_normalize(basis)
         call int_init(rtdb,1,basis)
         call schwarz_init(geom, basis)
         call int_1e_cache_ga(basis,oskel)
         call scf_vectors_guess(rtdb)
         call rohf_init(rtdb)
         if(.not.ga_create(mt_dbl,crohf_vlen,1,'rohf_nr_solve: grad',
     $        1,1,g_grad)) call errquit('rohf_nr_solve: ga_cre grad',0)
         if(nopen.gt.0)then
            call ga_copy(crohf_g_fcv,crohf_g_fcp)
            call ga_copy(crohf_g_fcv,crohf_g_fpv)
         endif
         call rohf_get_grad(g_grad)
         lshift=0.0d0
         accuracy=1d-3
         pflg=2
         max_iter=20
         max_sub=10
         oprint=.false.
         if(.not.ga_create(mt_dbl,crohf_vlen,1,'rohf_nr_solve: grad',
     $        1,1,g_work)) call errquit('rohf_nr_solve: ga_cre grad',0)
         call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
         oskel=.false.
         call ga_zero(g_work)
         if(.not.ga_iter_lsolve(crohf_vlen,max_iter,
     $        max_sub,accuracy, 
     $        rohf_hessv_precond,rohf_hessv,oprint, 
     $        g_grad,g_work,precond_acc,precond_iter))
     $        call errquit('cphf did not converge',1)
c$$$         if(.not.ga_it_ls_mrhs(crohf_vlen,1,max_iter,
c$$$     $        max_sub,accuracy, 
c$$$     $        rohf_hessv_precond,rohf_hessv,oprint, 
c$$$     $        g_grad,g_work,precond_acc,precond_iter))
c$$$     $        call errquit('cphf did not converge',1)
c$$$         pol=ga_ddot(g_work,g_grad)
c$$$         if(ga_nodeid().eq.0)then
c$$$            write(6,999) pol
c$$$ 999        format(//,15x,40('-'),
c$$$     $           /,15x,'RHF/ROHF Polarizability is:',f16.6,
c$$$     $           /,15x,40('-'),/)
c$$$         endif
         call int_1e_uncache_ga()
         call schwarz_tidy()
         call int_terminate()
         call rohf_tidy
         call scf_tidy (rtdb)
         if(.not.ga_destroy(g_grad))call errquit(
     $        'could not destroy gradient handle',1)
         if(.not.ga_destroy(g_work))call errquit(
     $        'could not destroy work array',1)
      else
         call int_normalize(basis)
         call int_init(rtdb,1,basis)
         call schwarz_init(geom, basis)
         call int_1e_cache_ga(basis,oskel)
         call scf_vectors_guess(rtdb)
         call uhf_init(rtdb)
         if(.not.cuhf_init_flag)
     $        call errquit('uhf_hessv: uhf internal block invalid',0)
         if(.not.ga_create(mt_dbl,cuhf_vlen,1,'uhf_nr_solve: grad',
     $        1,1,g_grad)) call errquit('uhf_nr_solve: ga_cre grad',0)
         call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)         
         if(.not.file_read_ga(cphf_rhs,g_grad))
     $        call errquit('could not read cphf_rhs',1)
         call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
         call waopen(cphf_sol_u,cphf_sol,1,0,ierr)
         if(ierr.ne.0)call errquit(
     $        'could not open cphf_sol file',1)
         lshift=0.0d0
         accuracy=1d-6
         pflg=2
         max_iter=20
         max_sub=10
         oprint=.true.
         if(.not.ga_create(mt_dbl,cuhf_vlen,1,'uhf_nr_solve: grad',
     $        1,1,g_work))call errquit('uhf_nr_solve: ga_cre grad',0)
         call uhf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
         oskel=.false.
         call ga_zero(g_work)
         if(.not.ga_iter_lsolve(cuhf_vlen,max_iter,
     $        max_sub,accuracy, 
     $        uhf_precond,uhf_hessv,oprint, 
     $        g_grad,g_work,precond_acc,precond_iter))
     $        call errquit('cphf did not converge',1)
c$$$         if(.not.ga_it_ls_mrhs(cuhf_vlen,1,max_iter,
c$$$     $        max_sub,accuracy, 
c$$$     $        uhf_precond,uhf_hessv,oprint, 
c$$$     $        g_grad,g_work,precond_acc,precond_iter))
c$$$     $        call errquit('cphf did not converge',1)
c$$$         pol=ga_ddot(g_work,g_grad)
c$$$         if(ga_nodeid().eq.0)then
c$$$            write(6,911) pol
c$$$ 911        format(//,15x,40('-'),
c$$$     $           /,15x,'UHF Polarizability is:',f16.6,
c$$$     $           /,15x,40('-'),/)
c$$$         endif
         call int_1e_uncache_ga()
         call schwarz_tidy()
         call int_terminate()
         call uhf_tidy(rtdb)
         call scf_tidy(rtdb)
         if(.not.ga_destroy(g_grad))call errquit(
     $        'could not destroy gradient handle',1)
         if(.not.file_write_ga(cphf_sol,g_work))
     $        call errquit('could not write cphf_sol',1)
         if(.not.ga_destroy(g_work))call errquit(
     $        'could not destroy work array',1)
      endif
      end
