      subroutine cphf_solve(rtdb)
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
      external rohf_hessv_precond,rohf_hessv
      external uhf_precond,uhf_hessv
      logical oprint
c     
      character*(nw_max_path_len) cphf_rhs
      character*(nw_max_path_len) cphf_sol
      logical file_write_ga,file_read_ga
      external file_write_ga,file_read_ga
      integer vlen, nvirt
c     
      call scf_init(rtdb)
      call int_normalize(basis)
      call int_init(rtdb,1,basis)
      call schwarz_init(geom, basis)
      call int_1e_cache_ga(basis,oskel)
      call scf_vectors_guess(rtdb)
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
      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve: rhs',1,1,g_grad))
     $     call errquit('cphf_solve: ga_cre grad',vlen)
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)         
      if(.not.file_read_ga(cphf_rhs, g_grad))
     $     call errquit('cphf_solve: could not read cphf_rhs',1)
c
      lshift=0.0d0
      precond_tol=1d-6
      pflg=2
      max_precond_iter=20
      max_precond_sub=10
      oprint=.true.
c
      if(.not.ga_create(mt_dbl,vlen,1,'cphf_solve: grad',
     $     1,1,g_work))call errquit('cphf_solve: ga_cre grad',0)
c
c     Compute necessary Fock matrices
c
      oskel=.false.
      if(scftype .eq. 'UHF')then
         call uhf_init(rtdb)
         call uhf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
      else
         call rohf_init(rtdb)
         call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
      endif
c
      if (sqrt(ga_ddot(g_work, g_work)) .gt. 1d-4) 
     $     call errquit('cphf_solve:the available MOs do not satisfy'//
     $     ' the SCF equations',0)
c
      call ga_zero(g_work)      ! Essential
c
      if (scftype .eq. 'UHF') then
         if(.not.ga_iter_lsolve(vlen, max_precond_iter,
     $        max_precond_sub, precond_tol, 
     $        uhf_precond, uhf_hessv, oprint_solve, 
     $        g_grad, g_work, precond_acc, precond_iter))
     $        call errquit('cphf did not converge',precond_iter)
      else 
         if(.not.ga_iter_lsolve(vlen, max_precond_iter,
     $        max_precond_sub, precond_tol, 
     $        rohf_hessv_precond, rohf_hessv, oprint_solve, 
     $        g_grad, g_work, precond_acc,precond_iter))
     $        call errquit('cphf did not converge', precond_iter)
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
c
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      if(.not.file_write_ga(cphf_sol,g_work))
     $     call errquit('cphf_solve: could not write cphf_sol',1)
      if(.not.ga_destroy(g_work))call errquit(
     $     'cphf_solve: could not destroy work array',1)
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
      if(.not.ga_destroy(g_grad))call errquit(
     $     'cphf_solve: could not destroy gradient handle',1)
c
      end
