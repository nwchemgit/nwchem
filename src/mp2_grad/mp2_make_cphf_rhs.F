      subroutine mp2_make_cphf_rhs(transform,basis,nmo,nbf,
     $     noa,nva,nob,nvb,irs_a,irs_b,work,map_a,map_b,
     $     num_oa,num_ob,num_va,num_vb,
     $     sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,
     $     sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $     g_vecs_a,g_vecs_b,g_lai_a,g_lai_b,g_rhs)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      logical transform
      integer irs_a(*),irs_b(*),work(*),map_a(*),map_b(*)
      integer sym_lo_oa(0:7),sym_hi_oa(0:7)
      integer sym_lo_va(0:7),sym_hi_va(0:7)
      integer sym_lo_ob(0:7),sym_hi_ob(0:7)
      integer sym_lo_vb(0:7),sym_hi_vb(0:7)
      integer num_oa(0:*),num_va(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer nmo,nbf,noa,nva,nob,nvb
      integer basis
      character*256 cphf_rhs
      integer cphf_rhs_u
      parameter(cphf_rhs_u=17)
      integer g_vecs_a,g_vecs_b,g_lai_a,g_lai_b
      integer g_rhs
      integer ierr
      logical file_write_ga
      external file_write_ga
      integer g_lai_a_tr,g_lai_b_tr
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)
      call waopen(cphf_rhs_u,cphf_rhs,1,0,ierr)
      if(ierr.ne.0)call errquit(
     $     'could not open cphf_rhs file',1)
      if(.not.ga_create(mt_dbl,nva,noa,'temp',1,1,g_lai_a_tr))
     $     call errquit('could not alloc tran array',1)
      if(.not.ga_create(mt_dbl,nvb,nob,'temp',1,1,g_lai_b_tr))
     $     call errquit('could not alloc tran array',1)
      call ga_transpose(g_lai_a,g_lai_a_tr)
      call ga_transpose(g_lai_b,g_lai_b_tr)
      call cphf_rhs_reshape(g_rhs,g_lai_a_tr,g_lai_b_tr,noa,nva,
     $     nob,nvb,nmo)
      if(.not.file_write_ga(cphf_rhs,g_rhs))
     $     call errquit('could not write cphf_rhs',1)
      if(.not.ga_destroy(g_lai_a_tr))call errquit(
     $     'could not destroy lai_a_tr handle',1)
      if(.not.ga_destroy(g_lai_b_tr))call errquit(
     $     'could not destroy lai_b_tr handle',1)
      end
c
c
c
      subroutine cphf_rhs_reshape(g_rhs,g_lai_a,g_lai_b,noa,
     $     nva,nob,nvb,nmo)
      implicit none
      integer g_rhs,g_lai_a,g_lai_b
      integer noa,nva,nob,nvb,nmo
      integer ioff
c
      call ga_copy_patch('n',g_lai_a,1,nva,1,noa,
     $     g_rhs,1,nva*noa,1,1)
      ioff=nva*noa
      call ga_copy_patch('n',g_lai_b,1,nvb,1,nob,
     $     g_rhs,ioff+1,ioff+(nvb*nob),1,1)
      end

