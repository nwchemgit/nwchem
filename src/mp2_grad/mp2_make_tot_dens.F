      subroutine mp2_make_tot_dens(g_p_ij_a,g_p_ij_b,g_p_ab_a,
     $     g_p_ab_b,g_rhs,g_p_tot_a,g_p_tot_b,noa,nva,nob,nvb,nbf,
     $     nmo,
     $     g_vecs_a,g_vecs_b,basis,scftype)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      character*(*) scftype
      integer nbf,nmo,noa,nva,nob,nvb
      integer g_p_ij_a,g_p_ij_b
      integer g_p_ab_a,g_p_ab_b
      integer g_rhs
      integer g_p_tot_a,g_p_tot_b
      integer g_scratch,g_vecs_a,g_vecs_b,g_tmp
      integer i,a,indx
      double precision temp
      integer basis
      double precision one,zero
      logical file_write_ga
      external file_write_ga
      character*256 p_a_ao
      character*256 p_b_ao
      integer base
c
      call util_file_name('p_a_ao',.true.,.true.,p_a_ao)
      call util_file_name('p_b_ao',.true.,.true.,p_b_ao)
      one=1.0d0
      zero=0.0d0
c
      call ga_dscal(g_rhs,2.0d0) ! Since our hessian is twice Frisch's
      call ga_zero(g_p_tot_a)
      call ga_zero(g_p_tot_b)
      call ga_copy_patch('n',g_rhs,1,nva*noa,1,1,g_p_tot_a,
     $     noa+1,nmo,1,noa)
c
      call ga_sync()
      do a=1+ga_nodeid(),nva,ga_nnodes()
         do i=1,noa
            indx=(i-1)*nva+a
            call ga_get(g_rhs,indx,indx,1,1,temp,1)
            call ga_put(g_p_tot_a,i,i,a+noa,a+noa,temp,1)
         enddo
      enddo
      call ga_copy_patch('n',g_p_ij_a,1,noa,1,noa,g_p_tot_a,
     $     1,noa,1,noa)
      call ga_copy_patch('n',g_p_ab_a,1,nva,1,nva,g_p_tot_a,
     $     noa+1,nmo,noa+1,nmo)
c
      if (scftype .eq. 'UHF') then
         base = (nva*noa)
      else
         base = 0               ! Alpha is the same as beta
      endif
      call ga_copy_patch('n',g_rhs,base+1,base+(nvb*nob),
     $     1,1,g_p_tot_b,nob+1,nmo,1,nob)
      do a=1+ga_nodeid(),nvb,ga_nnodes()
         do i=1,nob
            indx=(i-1)*nvb+a+base
            call ga_get(g_rhs,indx,indx,1,1,temp,1)
            call ga_put(g_p_tot_b,i,i,a+nob,a+nob,temp,1)
         enddo
      enddo
c
      call ga_copy_patch('n',g_p_ij_b,1,nob,1,nob,g_p_tot_b,
     $     1,nob,1,nob)
      call ga_copy_patch('n',g_p_ab_b,1,nvb,1,nvb,g_p_tot_b,
     $     nob+1,nmo,nob+1,nmo)
c
      if(.not.ga_create(mt_dbl,nmo,nbf,'scratch',1,1,g_scratch))
     $     call errquit('mp2_grad: scratch',1)
      if(.not.ga_create(mt_dbl,nbf,nbf,'scratch',1,1,g_tmp))
     $     call errquit('mp2_grad: scratch',2)
c
      call ga_dgemm('n','t',nmo,nbf,nmo,one,g_p_tot_a,g_vecs_a,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,one,g_vecs_a,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
c
      if(.not.file_write_ga(p_a_ao,g_tmp))
     $     call errquit('could not write p_tot_a',1)
c
      call ga_dgemm('n','t',nmo,nbf,nmo,one,g_p_tot_b,g_vecs_b,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,one,g_vecs_b,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
c
      if(.not.file_write_ga(p_b_ao,g_tmp))
     $     call errquit('could not write p_tot_b',1)
c
      if(.not.ga_destroy(g_scratch))call errquit(
     $     'could not destroy g_scratch',1)
      if(.not.ga_destroy(g_tmp))call errquit(
     $     'could not destroy g_scratch',1)
c
      end


