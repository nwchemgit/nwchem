      subroutine mp2_make_tot_wdens(g_w_ij_a,g_w_ij_b,g_w_ab_a,
     $     g_w_ab_b,g_w_ai_a,g_w_ai_b,g_w_tot_a,g_w_tot_b,
     $     noa,nva,nob,nvb,nbf,g_vecs_a,g_vecs_b,g_rhs,
     $     evals_a, evals_b, g_p_tot_a, g_p_tot_b)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf,noa,nva,nob,nvb
      integer g_vecs_a,g_vecs_b,g_rhs
      integer ierr
      integer g_w_ij_a,g_w_ij_b
      integer g_w_ab_a,g_w_ab_b
      integer g_w_ai_a,g_w_ai_b
      integer g_w_tot_a,g_w_tot_b
      integer g_w_ai_a_tr,g_w_ai_b_tr
      integer g_scratch
c
      double precision evals_a(*), evals_b(*)
      integer g_p_tot_a, g_p_tot_b
      double precision one,zero,two
      logical file_write_ga
      external file_write_ga
      character*256 w_a_ao
      character*256 w_b_ao
      integer w_a_ao_u,w_b_ao_u
      parameter(w_a_ao_u=21)
      parameter(w_b_ao_u=22)
c
c$$$      write(6,*) 'on entry to tot_wdens'
c$$$      call ga_print(g_w_ij_a)
c$$$      call ga_print(g_w_ab_a)
c$$$      call ga_print(g_w_ai_a)
c$$$      write(6,*) ' wai_a on input to tot_wdens'
c$$$      call ga_print(g_w_ai_a)
c$$$      write(6,*) 'wij on input to tot_wdens'
c$$$      call ga_print(g_w_ij_a)
c
      call util_file_name('w_a_ao',.true.,.true.,w_a_ao)
      call waopen(w_a_ao_u,w_a_ao,1,0,ierr)
      if(ierr.ne.0)call errquit(
     $     'could not open w_a_ao file',1)
      call util_file_name('w_b_ao',.true.,.true.,w_b_ao)
      call waopen(w_b_ao_u,w_b_ao,1,0,ierr)
      if(ierr.ne.0)call errquit(
     $     'could not open w_b_ao file',1)
      one=1.0d0
      zero=0.0d0
      two=2.0d0
      if(.not.ga_create(mt_dbl,nva,noa,'wai_a_tr',nva,1,
     $     g_w_ai_a_tr))call errquit(
     $     'could not alloc wai_a_tr',1)
      if(.not.ga_create(mt_dbl,nvb,nob,'wai_b_tr',nvb,1,
     $     g_w_ai_b_tr))call errquit(
     $     'could not alloc wai_b_tr',1)
c$$$      call ga_dscal(g_w_ai_a,two)
c$$$      call ga_dscal(g_w_ai_b,two)
c$$$      call ga_dscal(g_w_ai_a,-1.0d0)
c$$$      call ga_dscal(g_w_ai_b,-1.0d0)
      call ga_zero(g_w_tot_a)
      call ga_zero(g_w_tot_b)
      call ga_copy_patch('n',g_w_ai_a,1,noa,1,nva,g_w_tot_a,
     $     1, noa, noa+1, nbf)
c$$$      write(6,*) 'g_w_tot_a after first copy_patch'
c$$$      call ga_print(g_w_tot_a)
      call ga_transpose(g_w_ai_a,g_w_ai_a_tr)
      call ga_copy_patch('n',g_w_ai_a_tr,1,nva,1,noa,g_w_tot_a,
     $     noa+1,nbf,1,noa)
      call ga_copy_patch('n',g_w_ij_a,1,noa,1,noa,g_w_tot_a,
     $     1,noa,1,noa)
      call ga_copy_patch('n',g_w_ab_a,1,nva,1,nva,g_w_tot_a,
     $     noa+1,nbf,noa+1,nbf)
      call ga_copy_patch('n',g_w_ai_b,1,nob,1,nvb,g_w_tot_b,
     $     1, nob, nob+1,nbf)
      call ga_transpose(g_w_ai_b,g_w_ai_b_tr)
      call ga_copy_patch('n',g_w_ai_b_tr,1,nvb,1,nob,g_w_tot_b,
     $     nob+1,nbf,1,nob)
      call ga_copy_patch('n',g_w_ij_b,1,nob,1,nob,g_w_tot_b,
     $     1,nob,1,nob)
      call ga_copy_patch('n',g_w_ab_b,1,nvb,1,nvb,g_w_tot_b,
     $     nob+1,nbf,nob+1,nbf)
c
c     Wij += - Ei Pij
c     Wab += - Ea Pab
c     Wai += - Ei Pai
c
c     Brute force for now
c
c$$$      write(6,*) 'w before hack'
c$$$      call ga_print(g_w_tot_a)
      call mp2_w_hack(g_w_tot_a, g_p_tot_a, evals_a, noa, nva)
      call mp2_w_hack(g_w_tot_b, g_p_tot_b, evals_b, nob, nvb)
c
c...transform to AO basis and write
c
c$$$      write(6,*) ' TOTAL ALPHAS MO W'
c$$$      call ga_print(g_w_tot_a)
c     
c$$$      write(6,*) 'eigenvalues'
c$$$      call output(evals_a,1,nbf,1,1,nbf,1,1)
      if(.not.ga_create(mt_dbl,nbf,nbf,'scratch',1,1,g_scratch))
     $     call errquit('scratch',1)
      call ga_dgemm('n','t',nbf,nbf,nbf,one,g_w_tot_a,g_vecs_a,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_vecs_a,g_scratch,
     $     zero,g_w_tot_a)
      call ga_dgemm('n','t',nbf,nbf,nbf,one,g_w_tot_b,g_vecs_b,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_vecs_b,g_scratch,
     $     zero,g_w_tot_b)
c
      call ga_symmetrize(g_w_tot_a)
      call ga_symmetrize(g_w_tot_b)
c
c$$$      write(6,*) ' TOTAL ALPHAS AO W'
c$$$      call ga_print(g_w_tot_a)
c
c$$$      call ga_zero(g_w_tot_a)
c$$$      call ga_zero(g_w_tot_b)
      if(.not.file_write_ga(w_a_ao,g_w_tot_a))
     $     call errquit('could not write w_tot_a',1)
      if(.not.file_write_ga(w_b_ao,g_w_tot_b))
     $     call errquit('could not write w_tot_b',1)
      if(.not.ga_destroy(g_scratch))call errquit(
     $     'could not destroy g_scratch',1)
      if(.not.ga_destroy(g_w_ai_a_tr))call errquit(
     $     'could not destroy g_w_ai_a_tr handle',1)
      if(.not.ga_destroy(g_w_ai_b_tr))call errquit(
     $     'could not destroy g_w_ai_b_tr handle',1)
      end




      subroutine mp2_w_hack(g_w, g_p, e, no, nv)
      implicit none
#include "global.fh"
      integer g_w, g_p, no, nv
      double precision e(*)
c
      integer i, j, a, b
      double precision p
c
      do i = 1+ga_nodeid(), no,ga_nnodes()
         do j = 1, no
            call ga_get(g_p, i, i, j, j, p, 1)
            call ga_acc(g_w, i, i, j, j, -p*e(i), 1, 1.0d0)
         enddo
      enddo
c
      do a = no+1+ga_nodeid(), no+nv,ga_nnodes()
         do i = 1, no
            call ga_get(g_p, a, a, i, i, p, 1)
            call ga_acc(g_w, a, a, i, i, -p*e(i), 1, 1.0d0)
         enddo
      enddo
      do a = no+1+ga_nodeid(), no+nv,ga_nnodes()
         do i = 1, no
            call ga_get(g_p, i, i, a, a, p, 1)
            call ga_acc(g_w, i, i, a, a, -p*e(i), 1, 1.0d0)
         enddo
      enddo
c$$$c
      do a = no+1+ga_nodeid(), no+nv,ga_nnodes()
         do b = no+1, no+nv
            call ga_get(g_p, a, a, b, b, p, 1)
            call ga_acc(g_w, a, a, b, b, -p*e(a), 1, 1.0d0)
         enddo
      enddo
c$$$      write(6,*) 'w with hacks added'
c$$$      call ga_print(g_w)
c
      end

