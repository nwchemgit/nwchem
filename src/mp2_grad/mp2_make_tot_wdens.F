      subroutine mp2_make_tot_wdens(g_w_ij_a,g_w_ij_b,g_w_ab_a,
     $     g_w_ab_b,g_w_ai_a,g_w_ai_b,g_w_tot_a,g_w_tot_b,
     $     noa_lo, noa_hi,nva,nob_lo,nob_hi,
     $     nvb,nbf,nmo,g_vecs_a,g_vecs_b,
     $     evals_a, evals_b, g_p_tot_a, g_p_tot_b)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nbf,noa_hi,nva,nob_hi,nvb, nmo, noa_lo, nob_lo
      integer g_vecs_a,g_vecs_b
      integer g_w_ij_a,g_w_ij_b
      integer g_w_ab_a,g_w_ab_b
      integer g_w_ai_a,g_w_ai_b
      integer g_w_tot_a,g_w_tot_b
      integer g_w_ai_a_tr,g_w_ai_b_tr
      integer g_scratch, g_tmp
c     
      double precision evals_a(*), evals_b(*)
      integer g_p_tot_a, g_p_tot_b
      double precision one,zero,two
      logical file_write_ga
      external file_write_ga
      character*256 w_a_ao
      character*256 w_b_ao
c     
      call util_file_name('w_a_ao',.true.,.true.,w_a_ao)
      call util_file_name('w_b_ao',.true.,.true.,w_b_ao)
      one=1.0d0
      zero=0.0d0
      two=2.0d0
*ga:1:0
      if(.not.ga_create(mt_dbl,nva,noa_hi,'wai_a_tr',nva,0,
     $     g_w_ai_a_tr))call errquit(
     $     'could not alloc wai_a_tr',1, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nvb,nob_hi,'wai_b_tr',nvb,0,
     $     g_w_ai_b_tr))call errquit(
     $     'could not alloc wai_b_tr',1, GA_ERR)
c     
      call ga_zero(g_w_tot_a)
      call ga_zero(g_w_tot_b)
      call ga_copy_patch('n',g_w_ai_a,1,noa_hi,1,nva,g_w_tot_a,
     $     1, noa_hi, noa_hi+1, nmo)
c     
      call mp2_transpwai(g_w_ai_a,g_w_ai_a_tr)
      call ga_copy_patch('n',g_w_ai_a_tr,1,nva,1,noa_hi,g_w_tot_a,
     $     noa_hi+1,nmo,1,noa_hi)
      call ga_copy_patch('n',g_w_ij_a,1,noa_hi,1,noa_hi,g_w_tot_a,
     $     1,noa_hi,1,noa_hi)
      call ga_copy_patch('n',g_w_ab_a,1,nva,1,nva,g_w_tot_a,
     $     noa_hi+1,nmo,noa_hi+1,nmo)
      call ga_copy_patch('n',g_w_ai_b,1,nob_hi,1,nvb,g_w_tot_b,
     $     1, nob_hi, nob_hi+1,nmo)
      call mp2_transpwai(g_w_ai_b,g_w_ai_b_tr)
      call ga_copy_patch('n',g_w_ai_b_tr,1,nvb,1,nob_hi,g_w_tot_b,
     $     nob_hi+1,nmo,1,nob_hi)
      call ga_copy_patch('n',g_w_ij_b,1,nob_hi,1,nob_hi,g_w_tot_b,
     $     1,nob_hi,1,nob_hi)
      call ga_copy_patch('n',g_w_ab_b,1,nvb,1,nvb,g_w_tot_b,
     $     nob_hi+1,nmo,nob_hi+1,nmo)
c
      if(.not.ga_destroy(g_w_ai_a_tr))call errquit(
     $     'could not destroy g_w_ai_a_tr handle',1, GA_ERR)
      if(.not.ga_destroy(g_w_ai_b_tr))call errquit(
     $     'could not destroy g_w_ai_b_tr handle',1, GA_ERR)
c     
c     Wij += - Ei Pij
c     Wab += - Ea Pab
c     Wai += - Ei Pai
c     
c     Brute force for now
c     
      call mp2_w_hack(g_w_tot_a, g_p_tot_a, evals_a, noa_hi, nva)
      call mp2_w_hack(g_w_tot_b, g_p_tot_b, evals_b, nob_hi, nvb)
c     
c...  transform to AO basis and write
c     
*ga:1:0
      if(.not.ga_create(mt_dbl,nmo,nbf,'scratch',0,0,g_scratch))
     $     call errquit('mp2_grad: scratch',1, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nbf,nbf,'scratch',0,0,g_tmp))
     $     call errquit('mp2_grad: scratch',2, GA_ERR)
c
      call ga_dgemm('n','t',nmo,nbf,nmo,one,g_w_tot_a,g_vecs_a,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,one,g_vecs_a,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
      if(.not.file_write_ga(w_a_ao,g_tmp))
     $     call errquit('mp2_grad: could not write w_tot_a',1, DISK_ERR)
c     
      call ga_dgemm('n','t',nmo,nbf,nmo,one,g_w_tot_b,g_vecs_b,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,one,g_vecs_b,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
      if(.not.file_write_ga(w_b_ao,g_tmp))
     $     call errquit('mp2_grad: could not write w_tot_b',1, DISK_ERR)
c     
      if(.not.ga_destroy(g_tmp))call errquit(
     $     'mp2_grad: could not destroy g_scratch',1, GA_ERR)
      if(.not.ga_destroy(g_scratch))call errquit(
     $     'mp2_grad: could not destroy g_scratch',1, GA_ERR)
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
      call ga_sync()
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
c
      do a = no+1+ga_nodeid(), no+nv,ga_nnodes()
         do b = no+1, no+nv
            call ga_get(g_p, a, a, b, b, p, 1)
            call ga_acc(g_w, a, a, b, b, -p*e(a), 1, 1.0d0)
         enddo
      enddo
      call ga_sync()
c
      end

      subroutine mp2_transpwai(g,g_tr)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      integer g,g_tr
c
      integer ilo,ihi,jlo,jhi
      integer i,j
      integer k_buf,l_buf
      integer n_row,n_col,gtype
      integer dowork
c
      call ga_inquire(g_tr, gtype, n_row, n_col)
      call ga_distribution(g_tr, ga_nodeid(), 
     T     ilo, ihi, jlo, jhi)
      if ((ilo.eq.0).and.(ihi.eq.-1))  then
         dowork=0
      else
c     check if distrib by columns
         if(ilo.ne.1.and.ihi.ne.n_row) call
     R        errquit(' mp2transpwai bug ',0,0)
         dowork=1
      endif
      if(dowork.eq.1) then
         if (.not. ma_push_get(mt_dbl, n_row, 'mp2: wdens',
     $        l_buf, k_buf)) call errquit
     $        ('mp2: insufficient memory : bug ', n_row, MA_ERR)
      do j=jlo,jhi
         call ga_get(g,j,j,1,n_row,dbl_mb(k_buf),1)
          if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
          call ga_put(g_tr,1,n_row,j,j,dbl_mb(k_buf),n_row)
           if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
      enddo
      if (.not. ma_pop_stack(l_buf)) call errquit
     $     ('mp2transp pop_stack failed ', n_row, MA_ERR)
       
      endif
      call ga_sync()
      return
      end
