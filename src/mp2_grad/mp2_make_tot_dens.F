      subroutine mp2_make_tot_dens(g_p_ij_a,g_p_ij_b,g_p_ab_a,
     $     g_p_ab_b,g_rhs,g_p_tot_a,g_p_tot_b,noa,nva,nob,nvb,nbf,
     $     nmo,map_a, map_b,
     $     g_vecs_a,g_vecs_b,basis,scftype)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      character*(*) scftype
      integer map_a(*), map_b(*)
      integer nbf,nmo,noa,nva,nob,nvb
      integer g_p_ij_a,g_p_ij_b
      integer g_p_ab_a,g_p_ab_b
      integer g_rhs
      integer g_p_tot_a,g_p_tot_b
      integer g_scratch,g_vecs_a,g_vecs_b,g_tmp
      integer i,a, ia, ioff, ii, aa
      double precision xx
      integer basis
      double precision one,zero
      logical file_write_ga
      external file_write_ga
      character*256 p_a_ao
      character*256 p_b_ao
      integer base
      integer l_buf,k_buf
      integer ilo,ihi,jlo,jhi,dowork
c     
      call util_file_name('p_a_ao',.true.,.true.,p_a_ao)
      call util_file_name('p_b_ao',.true.,.true.,p_b_ao)
      one=1.0d0
      zero=0.0d0
c     
      call ga_dscal(g_rhs,2.0d0) ! Since our hessian is twice Frisch's
      call ga_zero(g_p_tot_a)
      call ga_zero(g_p_tot_b)
c     
c$$$  call ga_copy_patch('n',g_rhs,1,nva*noa,1,1,g_p_tot_a,
c$$$  $     noa+1,nmo,1,noa)
c$$$  c
c$$$  call ga_sync()
c$$$  do a=1+ga_nodeid(),nva,ga_nnodes()
c$$$  do i=1,noa
c$$$  indx=(i-1)*nva+a
c$$$  call ga_get(g_rhs,indx,indx,1,1,temp,1)
c$$$  call ga_put(g_p_tot_a,i,i,a+noa,a+noa,temp,1)
c$$$  end do
c$$$  end do
c     
      call ga_distribution(g_p_tot_a, ga_nodeid(), 
     T     ilo, ihi, jlo, jhi)
c     check if distrib by columns
      if ((ilo.eq.0).and.(ihi.eq.-1))  then
         dowork=0
      else
         if(ilo.ne.1.and.ihi.ne.nmo) call
     R        errquit(' mp2maketotdens bug ',0,0)
         dowork=1
      endif
      ioff = 0
      if(dowork.eq.1) then
      do i = 1,noa
         ii = map_a(i)
         if(ii.ge.jlo.and.ii.le.jhi) then
         do a = 1,nva
            aa = map_a(a+noa)-noa
            ia = ioff+a+(i-1)*nva
            call ga_get(g_rhs,ia,ia,1,1,xx,1)
            call ga_put(g_p_tot_a,aa+noa,aa+noa,ii,ii,xx,1)
            call ga_put(g_p_tot_a,ii,ii,aa+noa,aa+noa,xx,1)
         end do
         endif
      end do
      if (.not. ma_push_get(mt_dbl, nmo, 'mp2: maket',
     $     l_buf, k_buf)) call errquit
     $     ('mp2: insufficient memory : bug ', nmo, MA_ERR)

      do i=jlo,jhi
         if(i.le.noa) then
            call ga_get(g_p_ij_a,1,noa,i,i,dbl_mb(k_buf),noa)
            call ga_put(g_p_tot_a,1,noa,i,i,dbl_mb(k_buf),noa)
         else
            call ga_get(g_p_ab_a,1,nva,i-noa,i-noa,dbl_mb(k_buf),nva)
            call ga_put(g_p_tot_a,noa+1,nmo,i,i,dbl_mb(k_buf),nva)
         endif
      enddo
      if (.not. ma_pop_stack(l_buf)) call errquit
     $     ('mp2: pop_stack failed ', nmo, MA_ERR)
      endif
      call ga_mask_sync(.true.,.false.)
      call ga_sync()
c     
      if (scftype .eq. 'UHF') then
         base = (nva*noa)
      else
         base = 0               ! Alpha is the same as beta
      end if
      ioff = base
      if (dowork.eq.1) then
      do i = 1,nob
         ii = map_b(i)
         if(ii.ge.jlo.and.ii.le.jhi) then
         do a = 1,nvb
            aa = map_b(a+nob)-nob
            ia = ioff+a+(i-1)*nvb
            call ga_get(g_rhs,ia,ia,1,1,xx,1)
            call ga_put(g_p_tot_b,aa+nob,aa+nob,ii,ii,xx,1)
            call ga_put(g_p_tot_b,ii,ii,aa+nob,aa+nob,xx,1)
         end do
         endif
      end do
c$$$  call ga_copy_patch('n',g_rhs,base+1,base+(nvb*nob),
c$$$  $     1,1,g_p_tot_b,nob+1,nmo,1,nob)
c$$$  do a=1+ga_nodeid(),nvb,ga_nnodes()
c$$$  do i=1,nob
c$$$  indx=(i-1)*nvb+a+base
c$$$  call ga_get(g_rhs,indx,indx,1,1,temp,1)
c$$$  call ga_put(g_p_tot_b,i,i,a+nob,a+nob,temp,1)
c$$$  end do
c$$$  end do
c
      if (.not. ma_push_get(mt_dbl, nmo, 'mp2: maket',
     $     l_buf, k_buf)) call errquit
     $     ('mp2: insufficient memory : bug ', nmo, MA_ERR)

      do i=jlo,jhi
         if(i.le.nob) then
            call ga_get(g_p_ij_b,1,nob,i,i,dbl_mb(k_buf),nob)
            call ga_put(g_p_tot_b,1,nob,i,i,dbl_mb(k_buf),nob)
         else
            call ga_get(g_p_ab_b,1,nvb,i-nob,i-nob,dbl_mb(k_buf),nvb)
            call ga_put(g_p_tot_b,nob+1,nmo,i,i,dbl_mb(k_buf),nvb)
         endif
      enddo
      if (.not. ma_pop_stack(l_buf)) call errquit
     $     ('mp2: pop_stack failed ', nmo, MA_ERR)
      
      endif
      call ga_sync()
c     
*ga:1:0
      if(.not.ga_create(mt_dbl,nmo,nbf,'scratch',0,0,g_scratch))
     $     call errquit('mp2_grad: scratch',1, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nbf,nbf,'scratch',0,0,g_tmp))
     $     call errquit('mp2_grad: scratch',2, GA_ERR)
c     
      call ga_dgemm('n','t',nmo,nbf,nmo,one,g_p_tot_a,g_vecs_a,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,one,g_vecs_a,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
c     
      if(.not.file_write_ga(p_a_ao,g_tmp))
     $     call errquit('could not write p_tot_a',1, DISK_ERR)
c     
      call ga_dgemm('n','t',nmo,nbf,nmo,one,g_p_tot_b,g_vecs_b,
     $     zero,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,one,g_vecs_b,g_scratch,
     $     zero,g_tmp)
      call ga_symmetrize(g_tmp)
c     
      if(.not.file_write_ga(p_b_ao,g_tmp))
     $     call errquit('could not write p_tot_b',1, DISK_ERR)
c     
      if(.not.ga_destroy(g_scratch))call errquit(
     $     'could not destroy g_scratch',1, GA_ERR)
      if(.not.ga_destroy(g_tmp))call errquit(
     $     'could not destroy g_scratch',1, GA_ERR)
c     
      end


