      subroutine mp2_make_tot_dens(g_p_ij_a,g_p_ij_b,g_p_ab_a,
     $     g_p_ab_b,g_rhs,g_p_tot_a,g_p_tot_b,noa,nva,nob,nvb,nbf,
     $     g_vecs_a,g_vecs_b,basis,rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf,noa,nva,nob,nvb
      integer g_p_ij_a,g_p_ij_b
      integer g_p_ab_a,g_p_ab_b
      integer g_rhs
      integer g_p_tot_a,g_p_tot_b
      integer g_scratch,g_vecs_a,g_vecs_b,g_tmp
      integer i,a,indx
      double precision temp
      integer basis,rtdb
      double precision one,zero
      integer ierr
      logical file_write_ga
      external file_write_ga
      character*256 p_a_ao
      character*256 p_b_ao
      integer p_a_ao_u,p_b_ao_u
      parameter(p_a_ao_u=19)
      parameter(p_b_ao_u=20)
c
c$$$      write(6,*) ' PIJ_A'
c$$$      call ga_print(g_p_ij_a)
c$$$      write(6,*) ' PAB_A'
c$$$      call ga_print(g_p_ab_a)
c$$$      write(6,*) ' PAI_A'
c$$$      call ga_print(g_rhs)
      call util_file_name('p_a_ao',.true.,.true.,p_a_ao)
      call waopen(p_a_ao_u,p_a_ao,1,0,ierr)
      if(ierr.ne.0)call errquit(
     $     'could not open p_a_ao file',1)
      call util_file_name('p_b_ao',.true.,.true.,p_b_ao)
      call waopen(p_b_ao_u,p_b_ao,1,0,ierr)
      if(ierr.ne.0)call errquit(
     $     'could not open p_b_ao file',1)
      one=1.0d0
      zero=0.0d0
      call int_init(rtdb,1,basis)
      call ga_dscal(g_rhs,2.0d0)
      call ga_zero(g_p_tot_a)
      call ga_zero(g_p_tot_b)
      call ga_copy_patch('n',g_rhs,1,nva*noa,1,1,g_p_tot_a,
     $     noa+1,nbf,1,noa)
c$$$      call ga_copy_patch('t',g_rhs,1,nva*noa,1,1,g_p_tot_a,
c$$$     $     1,noa,noa+1,nbf)
      do i=1,noa
         do a=1,nva
            indx=(i-1)*nva+a
            call ga_get(g_rhs,indx,indx,1,1,temp,1)
            call ga_put(g_p_tot_a,i,i,a+noa,a+noa,temp,1)
         enddo
      enddo
      call ga_copy_patch('n',g_p_ij_a,1,noa,1,noa,g_p_tot_a,
     $     1,noa,1,noa)
      call ga_copy_patch('n',g_p_ab_a,1,nva,1,nva,g_p_tot_a,
     $     noa+1,nbf,noa+1,nbf)
c
*      write(6,*) ' TOTAL P_A'
*      call ga_print(g_p_tot_a)
c
      call ga_copy_patch('n',g_rhs,(nva*noa)+1,(nva*noa)+(nvb*nob)
     $     ,1,1,g_p_tot_b,nob+1,nbf,1,nob)
      do i=1,nob
         do a=1,nvb
            indx=(i-1)*nvb+a+(nva*noa)
            call ga_get(g_rhs,indx,indx,1,1,temp,1)
            call ga_put(g_p_tot_b,i,i,a+nob,a+nob,temp,1)
         enddo
      enddo
c$$$      call ga_copy_patch('t',g_rhs,(nva*noa)+1,(nva*noa)+(nvb*nob)
c$$$     $     ,1,1,g_p_tot_b,1,nob,nob+1,nbf)
      call ga_copy_patch('n',g_p_ij_b,1,nob,1,nob,g_p_tot_b,
     $     1,nob,1,nob)
      call ga_copy_patch('n',g_p_ab_b,1,nvb,1,nvb,g_p_tot_b,
     $     nob+1,nbf,nob+1,nbf)
c$$$c...make the overlap matrix
c$$$      if(.not.ga_create(mt_dbl,nbf,nbf,'overlap',1,1,g_over))
c$$$     $     call errquit('could not create overlap',1)
c$$$      call int_1e_ga(basis,basis,g_over,'overlap',.false.)
c$$$c...make the inverse of the MO coefficients
c$$$      if(.not.ga_create(mt_dbl,nbf,nbf,'inv_mo',1,1,g_vecs_a_inv))
c$$$     $     call errquit('inv_mo',1)
c$$$      if(.not.ga_create(mt_dbl,nbf,nbf,'inv_mo',1,1,g_vecs_b_inv))
c$$$     $     call errquit('inv_mo',1)
      if(.not.ga_create(mt_dbl,nbf,nbf,'scratch',1,1,g_scratch))
     $     call errquit('scratch',1)
      if(.not.ga_create(mt_dbl,nbf,nbf,'scratch',1,1,g_tmp))
     $     call errquit('scratch',1)
c$$$      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_over,g_vecs_a,
c$$$     $     zero,g_vecs_a_inv)
c$$$      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_over,g_vecs_b,
c$$$     $     zero,g_vecs_b_inv)
c$$$c...transform to AO basis and write
c$$$      call two_index_transf(g_p_tot_a,g_vecs_a_inv,g_vecs_a_inv,
c$$$     $     g_scratch,g_p_tot_a)
c$$$      call two_index_transf(g_p_tot_b,g_vecs_b_inv,g_vecs_b_inv,
c$$$     $     g_scratch,g_p_tot_b)
c$$$      call int_terminate()
c     
c$$$      write(6,*) 'total density in MO'
c$$$      call ga_print(g_p_tot_a)
      call ga_dgemm('n','t',nbf,nbf,nbf,one,g_p_tot_a,g_vecs_a,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_vecs_a,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
c      call ga_zero(g_tmp)
c$$$      write(6,*) 'total density in AO'
c$$$      call ga_print(g_tmp)
      if(.not.file_write_ga(p_a_ao,g_tmp))
     $     call errquit('could not write p_tot_a',1)
c$$$c
      call ga_dgemm('n','t',nbf,nbf,nbf,one,g_p_tot_b,g_vecs_b,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_vecs_b,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
c      call ga_zero(g_tmp)
      if(.not.file_write_ga(p_b_ao,g_tmp))
     $     call errquit('could not write p_tot_b',1)
c
      if(.not.ga_destroy(g_scratch))call errquit(
     $     'could not destroy g_scratch',1)
      if(.not.ga_destroy(g_tmp))call errquit(
     $     'could not destroy g_scratch',1)
c
c$$$      if(.not.ga_destroy(g_over))call errquit(
c$$$     $     'could not destroy g_over',1)
c$$$      if(.not.ga_destroy(g_vecs_a_inv))call errquit(
c$$$     $     'could not destroy g_vecs_a_inv',1)
c$$$      if(.not.ga_destroy(g_vecs_b_inv))call errquit(
c$$$     $     'could not destroy g_vecs_b_inv',1)
      end


