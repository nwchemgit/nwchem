      subroutine cphf(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "cscf.fh"
#include "crohf.fh"
      integer rtdb              ! [input] database handle
      logical status
      double precision eone,etwo,enrep,energy,pol
      integer g_grad,g_x,g_y,g_z,g_tmp,g_work
      integer max_iter,max_sub
      double precision accuracy,precond_acc,precond_iter
      logical ga_iter_lsolve
      external ga_iter_lsolve
      external rohf_hessv_precond, rohf_hessv
      logical oprint
      if(.not.context_push('cphf'))
     $     call errquit('cphf: context_push failed',0)
c     
      status=rtdb_parallel(.true.) ! Broadcast reads to all processes
c     
      if(ga_nodeid().eq.0)then
        call util_print_centered(6,'NWCHEM CPHF Module',40,.true.)
       endif
c     
      call scf_init(rtdb)
      call int_normalize(basis)
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
      call scf_vectors_guess(rtdb)
      call rohf_init(rtdb)
      if(.not.ga_create(mt_dbl,crohf_vlen,1,'rohf_nr_solve: grad',
     $     1,1,g_grad)) call errquit('rohf_nr_solve: ga_cre grad',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_x)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_y)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_z)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_tmp)) call errquit('ga',0)
      call int_dip_ga(basis,basis,g_x,g_y,g_z)
c      call ga_print(g_x)
      call two_index_transf(g_x,g_movecs,
     $     g_movecs,g_tmp,crohf_g_fcv)
      if(nopen.gt.0)then
         call ga_copy(crohf_g_fcv,crohf_g_fcp)
         call ga_copy(crohf_g_fcv,crohf_g_fpv)
      endif
c      call ga_summarize(1)
      call rohf_get_grad(g_grad)
      lshift=0.0d0
      accuracy=1d-3
      pflg=2
      max_iter=20
      max_sub=10
      oprint=.false.
      if(.not.ga_create(mt_dbl,crohf_vlen,1,'rohf_nr_solve: grad',
     $     1,1,g_work)) call errquit('rohf_nr_solve: ga_cre grad',0)
      call rohf_energy(g_movecs,eone,etwo,enrep,energy,g_work)
c      call ga_print(g_work)
      oskel=.false.
      call ga_zero(g_work)
      if(.not.ga_iter_lsolve(crohf_vlen,max_iter,
     $     max_sub,accuracy, 
     $     rohf_hessv_precond,rohf_hessv,oprint, 
     $     g_grad,g_work,precond_acc,precond_iter))
     $     call errquit('cphf did not converge',1)
      pol=ga_ddot(g_work,g_grad)
      if(ga_nodeid().eq.0)then
         write(6,999) pol
 999     format(//,15x,31('-'),
     $          /,15x,'Polarizability is:',f16.10,
     $          /,15x,31('-'),/)
      endif
      call schwarz_tidy()
      call int_terminate()
      call rohf_tidy
      call scf_tidy 
      if(.not.ga_destroy(g_grad))call errquit(
     $     'could not destroy gradient handle',1)
      if(.not.ga_destroy(g_x))call errquit(
     $     'could not destroy < |x| > handle',1)
      if(.not.ga_destroy(g_y))call errquit(
     $     'could not destroy < |y| > handle',1)
      if(.not.ga_destroy(g_z))call errquit(
     $     'could not destroy < |z| > handle',1)
      if(.not.ga_destroy(g_tmp))call errquit(
     $     'could not destroy tmp handle',1)
      if(.not.ga_destroy(g_work))call errquit(
     $     'could not destroy work array',1)
      if(.not.context_pop('cphf'))
     $     call errquit('context_pop failed',0)
      end
