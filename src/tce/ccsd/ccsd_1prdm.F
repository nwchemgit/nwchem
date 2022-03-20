      subroutine ccsd_1prdm(rtdb,d_t1,d_t2,d_lambda1,d_lambda2,
     1     k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "stdio.fh"
c
      integer rtdb
      integer d_t1,d_t2,d_lambda1,d_lambda2
      integer k_t1_offset,k_t2_offset
      integer k_l1_offset,k_l2_offset
c
      logical nodezero
      integer l_hh_offset,k_hh_offset
      integer d_hh, size_hh
      integer l_ph_offset,k_ph_offset
      integer d_ph, size_ph
      integer l_pp_offset,k_pp_offset
      integer d_pp, size_pp
      integer l_hp_offset,k_hp_offset
      integer d_hp, size_hp
      integer dim_rdm_ao, l_rdm_ao, k_rdm_ao
      integer nh,np
      integer dim_rdm_mo_hh,dim_rdm_mo_hp
      integer dim_rdm_mo_ph,dim_rdm_mo_pp
      integer l_rdm_mo_hh,k_rdm_mo_hh
      integer l_rdm_mo_hp,k_rdm_mo_hp
      integer l_rdm_mo_ph,k_rdm_mo_ph
      integer l_rdm_mo_pp,k_rdm_mo_pp
      integer l_mo_p,k_mo_p,l_mo_h,k_mo_h
      integer l_mo_p_tmp,k_mo_p_tmp
      integer l_mo_h_tmp,k_mo_h_tmp
      integer i,j
      integer dim_mu_h,dim_mu_p,dim_mo_h,dim_mo_p
      integer l_mu_h,k_mu_h,l_mu_p,k_mu_p
      integer hole,particle
      integer g_dens
      character*256 filename
      logical  ao_rdm_write
      external ao_rdm_write
c
c
      nodezero=(ga_nodeid().eq.0)      
      dim_rdm_ao = nbf*nbf
      if(.not.ma_push_get(mt_dbl,dim_rdm_ao,'rdm_ao',
     1     l_rdm_ao,k_rdm_ao)) 
     2     call errquit('ccsd_1prdm: ma problem',1,ma_err)
      call dfill(dim_rdm_ao,0d0,dbl_mb(k_rdm_ao),1)
c     
c     allocate memory for MOs
c     
      nh=nocc(1)+nocc(ipol)
      np=nmo(1)+nmo(ipol)-nh
      dim_mo_h=nh*nbf
      dim_mo_p=np*nbf
      if(.not.ma_push_get(mt_dbl,dim_mo_h,'mo_h',
     1     l_mo_h,k_mo_h))
     2     call errquit('ccsd_1prdm: ma problem', 2, ma_err)
      do i=1,dim_mo_h
         dbl_mb(k_mo_h+i-1)=0.d0
      enddo
      if(.not.ma_push_get(mt_dbl,dim_mo_p,'mo_p',
     1     l_mo_p,k_mo_p))
     2     call errquit('ccsd_1prdm: ma problem', 3, ma_err)
      do i=1,dim_mo_p
         dbl_mb(k_mo_p+i-1)=0.d0
      enddo
c     
c     get the MOs from GA and make them sorted (list them according column)
c     
      if(.not.ma_push_get(mt_dbl,dim_mo_h,'mo_h_tmp',
     1     l_mo_h_tmp,k_mo_h_tmp))
     2     call errquit('ccsd_1prdm: ma problem', 2, ma_err)
      do i=1,dim_mo_h
         dbl_mb(k_mo_h_tmp+i-1)=0.d0
      enddo           
c     
c     hole alpha
c     
      do hole=1,nocc(1)
         i=2*hole-1
         call ga_get(g_movecs(1),1,nbf,hole,hole,
     1        dbl_mb(k_mo_h_tmp+(i-1)*nbf),nbf)
      enddo
c     
c     hole beta
c     
      do hole=1,nocc(ipol)
         i=2*hole
         call ga_get(g_movecs(ipol),1,nbf,hole,hole,
     1        dbl_mb(k_mo_h_tmp+(i-1)*nbf),nbf)
      enddo
c     
c     make them sorted according to column index
c     
      do i=1,nh
         do j=1,nbf                 
            dbl_mb(k_mo_h+(j-1)*nh+(i-1))= 
     1           dbl_mb(k_mo_h_tmp+(i-1)*nbf+(j-1))
         enddo
      enddo
      if(.not.ma_pop_stack(l_mo_h_tmp))
     1     call errquit('ccsd_1prdm: ma problem',00,ma_err)
c     
c     particle alpha
c    
      if(.not.ma_push_get(mt_dbl,dim_mo_p,'mo_p_tmp',
     1     l_mo_p_tmp,k_mo_p_tmp))
     2     call errquit('ccsd_1prdm: ma problem', 2, ma_err)
      do i=1,dim_mo_p
         dbl_mb(k_mo_p_tmp+i-1)=0.d0
      enddo
c
      do particle=nocc(1)+1, nmo(1)
         i=2*particle-nh-1
         call ga_get(g_movecs(1),1,nbf,particle,particle,
     1        dbl_mb(k_mo_p_tmp+(i-1)*nbf),nbf)
      enddo
c     
c     particle beta
c
      do particle=nocc(ipol)+1,nmo(ipol)
         i=2*particle-nh
         call ga_get(g_movecs(ipol),1,nbf,particle,particle,
     1        dbl_mb(k_mo_p_tmp+(i-1)*nbf),nbf)
      enddo
c     
c     make them sorted according to column index
c     
      do i=1,np
         do j=1,nbf
            dbl_mb(k_mo_p+(j-1)*np+(i-1))=
     1           dbl_mb(k_mo_p_tmp+(i-1)*nbf+(j-1))
         enddo
      enddo
      if(.not.ma_pop_stack(l_mo_p_tmp))
     1     call errquit('ccsd_1prdm: ma problem',00,ma_err)
      
c->debug     
c     for debugging density on AO
c     
c      if(nodezero) then
c         open(1, status='unknown',file='mo.dat',
c     1        form='formatted')
c         do i=1,nh
c            do j=1,nbf
c               write(1,'(i5,i5,f20.16)') i, j, 
c     1           dbl_mb(k_mo_h+(j-1)*nh+(i-1))
c            enddo
c         enddo
c         do i=1,np
c            do j=1,nbf
c               write(1,'(i5,i5,f20.16)') i+nh,j,
c     1           dbl_mb(k_mo_p+(j-1)*np+(i-1))
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
c     
c     allocate memory for intermediates
c
      dim_mu_h = nbf*nh
      dim_mu_p = nbf*np
      if(.not.ma_push_get(mt_dbl,dim_mu_h,'mu_h',
     1     l_mu_h,k_mu_h))
     2     call errquit('ccsd_1prdm: ma problem', 2, ma_err)
      do i=1,dim_mu_h
         dbl_mb(k_mu_h+i-1)=0.d0
      enddo
      if(.not.ma_push_get(mt_dbl,dim_mu_p,'mu_p',
     1     l_mu_p,k_mu_p))
     2     call errquit('ccsd_1prdm: ma problem', 3, ma_err)
      do i=1, dim_mu_p
         dbl_mb(k_mu_p+i-1)=0.d0
      enddo
c     
c     ccsd_1prdm_hh
c     
c     allocate memory for hh block
c
      dim_rdm_mo_hh = nh*nh
      if(.not.ma_push_get(mt_dbl,dim_rdm_mo_hh,'rdm_mo_hh',
     1     l_rdm_mo_hh,k_rdm_mo_hh))
     2     call errquit('ccsd_1prdm: ma problem',4,ma_err)
      do i=1,dim_rdm_mo_hh
         dbl_mb(k_rdm_mo_hh+i-1)=0.d0
      enddo
      call tce_filename('hh',filename)
      call tce_dens_hh_offset(l_hh_offset,k_hh_offset,size_hh)
      call createfile(filename,d_hh,size_hh)
      call ccsd_1prdm_hh(d_hh,d_t1,d_t2,
     1     d_lambda1,d_lambda2,
     2     k_hh_offset,k_t1_offset,k_t2_offset,
     1     k_l1_offset,k_l2_offset)
      call reconcilefile(d_hh,size_hh)
      call get_mo_rdm_hh(d_hh,k_hh_offset,k_rdm_mo_hh,1.d0)
      call deletefile(d_hh)
      if (.not.ma_pop_stack(l_hh_offset))
     1     call errquit('ccsd_1prdm: ma problem',5,ma_err)
c->debug
c      if(nodezero) then
c         open(1,file='rdm_mo_hh.dat',status='unknown',
c     1        form='formatted')
c         do i=1,nh
c            do j=1,nh
c               write(1,'(i5,i5,f20.16)'),i,j,
c     1         dbl_mb(k_rdm_mo_hh+(j-1)*nh+(i-1))
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
c     
c     do the matrix multiplication C_{\mu}^{h} Gamma_{h}^{h}
c     
      call dgemm('t','n',nbf,nh,nh,1.d0,dbl_mb(k_mo_h),
     1     nh,dbl_mb(k_rdm_mo_hh),nh,
     2     0.d0,dbl_mb(k_mu_h),nbf)
c->debug
c      if(nodezero) then
c         open(1,file='mu_h_hh.dat',status='unknown',
c     1        form='formatted')
c         do i=1,nbf
c            do j=1,nh
c              write(1,'(i5,i5,f20.16)') i,j, 
c     1              dbl_mb(k_mu_h+(j-1)*nbf+(i-1))
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
c     
c     do the matrxi multiplication Gamma'_{\mu}^{h} C_{h}^{\nu}
c     
      call dgemm('n','n',nbf,nbf,nh,1.0d0,dbl_mb(k_mu_h),nbf,
     1     dbl_mb(k_mo_h),nh,0.d0, dbl_mb(k_rdm_ao),nbf)
c->debug
c      if(nodezero) then
c         open(1,file='rdm_ao_hh.dat',status='unknown',
c     1           form='formatted')
c         do i=1,nbf
c            do j=1,nbf
c               write(1,'(i5,i5,f20.16)') i,j, 
c     1          dbl_mb(k_rdm_ao+(j-1)*nbf +(i-1))
c            enddo
c         enddo
c         close(1)
c       endif
c<-debug
c     
c     release hh block memory
c     
c
      if(.not.ma_pop_stack(l_rdm_mo_hh))
     1     call errquit('ccsd_1prdm: ma problem',6,ma_err)
      do i=1,dim_mu_h
         dbl_mb(k_mu_h+i-1)=0.d0
      enddo
c     
c     ccsd_1prdm_hp
c     
c     
c     allocate memory for hp block
c     
c
      dim_rdm_mo_hp = nh*np
      if(.not.ma_push_get(mt_dbl,dim_rdm_mo_hp,'rdm_mo_hp',
     1     l_rdm_mo_hp,k_rdm_mo_hp))
     2     call errquit('ccsd_1prdm: ma problem',7,ma_err)
      do i=1,dim_rdm_mo_hp
         dbl_mb(k_rdm_mo_hp+i-1)=0.d0
      enddo
c     
      call tce_filename('hp',filename)
      call tce_dens_hp_offset(l_hp_offset,k_hp_offset,size_hp)
      call createfile(filename,d_hp,size_hp)
      call ccsd_1prdm_hp(d_hp,d_lambda1,
     1     k_hp_offset,k_l1_offset)
      call reconcilefile(d_hp,size_hp)
      call get_mo_rdm_hp(d_hp,k_hp_offset,k_rdm_mo_hp,1.d0)
      call deletefile(d_hp)
      if (.not.ma_pop_stack(l_hp_offset))
     1     call errquit('ccsd_1prdm: ma problem',8,MA_ERR)
c->debug
c      if(nodezero) then
c         open(1,file='rdm_mo_hp.dat',status='unknown',
c     1       form='formatted')
c         do i=1,nh
c            do j=1,np
c               write(1,'(i5,i5,f20.16)') i,j+nh,
c     1              dbl_mb(k_rdm_mo_hp+(i-1)+(j-1)*nh)
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
c     
c     do the matrix multiplication
c     
c     C_{\mu}^{h} Gamma_{h}^{p}
c     
      call dgemm('t','n',nbf,np,nh,1.0d0,dbl_mb(k_mo_h),
     1     nh,dbl_mb(k_rdm_mo_hp),nh,
     2     0.d0,dbl_mb(k_mu_p),nbf)
c->debug
c      if(nodezero) then
c         open(1,file='mu_p_hp.dat',status='unknown',
c     1        form='formatted')
c         do i=1,nbf
c            do j=1,np
c               write(1,'(i5,i5,f20.16)') i,j+nh,
c     1           dbl_mb(k_mu_p+(i-1)+(j-1)*nbf)
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
c     
c     Gamma_{\mu}^{p} C_{p}^{\nu}
c     
      call dgemm('n','n',nbf,nbf,np,1.0d0,dbl_mb(k_mu_p),nbf,
     1     dbl_mb(k_mo_p),np,1.0d0, dbl_mb(k_rdm_ao),nbf)
c->debug
c      if(nodezero) then
c         open(1,file='rdm_ao_hp.dat',status='unknown',
c     1       form='formatted')
c         do i=1,nbf
c            do j=1,nbf
c               write(1,'(i5,i5,f20.16)') i,j,
c     1          dbl_mb(k_rdm_ao+(i-1)+(j-1)*nbf)
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
c     
c     release hp block memory
c     
      if(.not.ma_pop_stack(l_rdm_mo_hp))
     1     call errquit('ccsd_1prdm: ma problem',9,ma_err)
      do i=1,dim_mu_p
         dbl_mb(k_mu_p+i-1)=0.d0
      enddo
c
c     ccsd_1prdm_ph
c     
c     allocate memory for ph block
c     
      dim_rdm_mo_ph=np*nh
      if(.not.ma_push_get(mt_dbl,dim_rdm_mo_ph,'rdm_mo_ph',
     1     l_rdm_mo_ph, k_rdm_mo_ph))
     2     call errquit('ccsd_1prdm: ma problem',91,ma_err)
      do i=1,dim_rdm_mo_ph
         dbl_mb(k_rdm_mo_ph+i-1)=0.d0
      enddo
      call tce_filename('ph',filename)
      call tce_dens_ph_offset(l_ph_offset,k_ph_offset,size_ph)
      call createfile(filename,d_ph,size_ph)
      call ccsd_1prdm_ph(d_ph,d_t1,d_t2,d_lambda1,d_lambda2,
     1     k_ph_offset,k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset)
      call reconcilefile(d_ph,size_ph)
      call get_mo_rdm_ph(d_ph,k_ph_offset,k_rdm_mo_ph,1.d0)
      call deletefile(d_ph)
      if (.not.ma_pop_stack(l_ph_offset))
     1     call errquit('ccsd_1prdm: ma problem',10,ma_err)
c
c->debug
c
c      if(nodezero) then
c         open(1,file='rdm_mo_ph.dat',status='unknown',
c     1        form='formatted')
c         do i=1,np
c            do j=1,nh
c               write(1,'(i5,i5,f20.16)') i+nh,j,
c     1              dbl_mb(k_rdm_mo_ph+(i-1)+(j-1)*np)
c            enddo
c         enddo
c         close(1)
c      endif
c
c<-debug
c
c     
c     do the matrix multiplication
c     
      call dgemm('t','n',nbf,nh,np,1.0d0,dbl_mb(k_mo_p),
     1     np,dbl_mb(k_rdm_mo_ph),np,
     2     0.d0,dbl_mb(k_mu_h),nbf)
c->debug
c      if(nodezero) then
c         open(1,file='mu_h_ph.dat',status='unknown',
c     1      form='formatted')
c         do i=1,nbf
c            do j=1,nh
c               write(1,'(i5,i5,f20.16)') i,j,
c     1         dbl_mb(k_mu_h+(i-1)+(j-1)*nbf)
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
      call dgemm('n','n',nbf,nbf,nh,1.0d0,dbl_mb(k_mu_h),nbf,
     1     dbl_mb(k_mo_h),nh,1.0d0, dbl_mb(k_rdm_ao),nbf)
c->debug
      if(nodezero.and.util_print('densmat',print_high)) then
        write(luout,*) '==================================='
        write(luout,*) 'Debug information of density matrix'
        write(luout,*) '==================================='
        do i=1,nbf
           do j=1,nbf
              if(abs(dbl_mb(k_rdm_ao+(i-1)+(j-1)*nbf)).gt.1.d-8)
     1           write(luout,'(i5,i5,f20.16)') i,j,
     1            (dbl_mb(k_rdm_ao+(i-1)+(j-1)*nbf)+
     1             dbl_mb(k_rdm_ao+(j-1)+(i-1)*nbf))/2.d0
           enddo
        enddo
        write(luout,*) '==================================='
        write(luout,*) '      End of debug information     '
        write(luout,*) '==================================='
      endif
c<-debug
c     
c     release memory for ph block
c     
      if(.not.ma_pop_stack(l_rdm_mo_ph))
     1     call errquit('ccsd_1prdm: ma problem',11,ma_err)
      do i=1,dim_mu_h
         dbl_mb(k_mu_h+i-1)=0.d0
      enddo
c     
c     -rdens_1e_pp
c     
c     allocate memory for pp block
c     
      dim_rdm_mo_pp=np*np
      if(.not.ma_push_get(mt_dbl,dim_rdm_mo_pp,'rdm_mo_pp',
     1     l_rdm_mo_pp,k_rdm_mo_pp))
     2     call errquit('ccsd_1prdm: ma problem',111,ma_err)
      do i=1,dim_rdm_mo_pp
         dbl_mb(k_rdm_mo_pp+i-1)=0.d0
      enddo
      call tce_filename('pp',filename)
      call tce_dens_pp_offset(l_pp_offset,k_pp_offset,size_pp)
      call createfile(filename,d_pp,size_pp)
      call ccsd_1prdm_pp(d_pp,d_t1,d_t2,
     1     d_lambda1,d_lambda2,
     2     k_pp_offset,k_t1_offset,k_t2_offset,
     3     k_l1_offset,k_l2_offset)
      call reconcilefile(d_pp,size_pp)
      call get_mo_rdm_pp(d_pp,k_pp_offset,k_rdm_mo_pp,1.d0)
      call deletefile(d_pp)
      if (.not.ma_pop_stack(l_pp_offset))
     1     call errquit('ccsd_1prdm: ma problem',12,ma_err)
c->debug
c      if(nodezero) then
c         open(1,file='rdm_mo_pp.dat',status='unknown',
c     1       form='formatted')
c         do i=1,np
c            do j=1,np
c               write(1,'(i5,i5,f20.16)') i+nh,j+nh,
c     1              dbl_mb(k_rdm_mo_pp+(i-1)+(j-1)*np)
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
c     
c     do the matrix multiplication
c     
      call dgemm('t','n',nbf,np,np,1.0d0,dbl_mb(k_mo_p),
     1     np,dbl_mb(k_rdm_mo_pp),np,
     2     0.d0,dbl_mb(k_mu_p),nbf)
c->debug
c      if(nodezero) then
c         open(1,file='mu_p_pp.dat',status='unknown',
c     1     form='formatted')
c         do i=1,nbf
c            do j=1,np
c               write(1,'(i5,i5,f20.16)') i,j+nh,
c     1          dbl_mb(k_mu_p+(i-1)+(j-1)*nbf)
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug
      call dgemm('n','n',nbf,nbf,np,1.0d0,dbl_mb(k_mu_p),nbf,
     1     dbl_mb(k_mo_p),np,1.0d0, dbl_mb(k_rdm_ao),nbf)
c->debug
c      if(nodezero) then
c         open(1,file='rdm_ao_pp.dat',status='unknown',
c     1        form='formatted')
c         do i=1,nbf
c            do j=1,nbf
c               write(1,'(i5,i5,f30.16)') i,j,
c     1          dbl_mb(k_rdm_ao+(i-1)+(j-1)*nbf)
c            enddo
c         enddo
c         close(1)
c      endif
c<-debug      
c     
c     release the memory for pp block
c     
      if(.not.ma_pop_stack(l_rdm_mo_pp))
     1     call errquit('ccsd_1prdm: ma problem',12,ma_err)
c     
c     release the memory for intermediates
c     
      if(.not.ma_pop_stack(l_mu_p))
     1     call errquit('ccsd_1prdm: ma problem',13,ma_err)
      if(.not.ma_pop_stack(l_mu_h))
     1     call errquit('ccsd_1prdm: ma problem',14,ma_err)
c     
c     release the memory for MOs
c     
      if(.not.ma_pop_stack(l_mo_p))
     1     call errquit('ccsd_1prdm: ma problem',15,ma_err)
      if(.not.ma_pop_stack(l_mo_h))
     1     call errquit('ccsd_1prdm: ma problem',16,ma_err)
c     
c     dump the ao rdm to a file
c     
      if (.not.rtdb_cget(rtdb,'tce:file_densmat',1,filename))
     1   call errquit('ccsd_1prdm: rtdb_cgett failed - file_densmat',0,
     1        RTDB_ERR)
c
#if 0
      if(.not.ao_rdm_write(filename,k_rdm_ao))
     1     call errquit('ccsd_1prdm: disk problem',1,disk_err)
#else
c     nodezero loads DM into g_dens
      if(.not.ga_duplicate(g_movecs(1),g_dens,'dm'))
     1     call errquit('ccsd_1prdm: ga_duplicate failed',0,GA_ERR)
      call ga_put(g_dens,1,nbf,1,nbf,dbl_mb(k_rdm_ao),nbf)
      call ao_1prdm_write(nbf,g_dens,filename)
      if (.not. ga_destroy(g_dens))
     1     call errquit('ccsd_1prdm: ga_destroy failed',0,GA_ERR)
#endif
      if(.not.ma_pop_stack(l_rdm_ao))
     1     call errquit('ccsd_1prdm: ma problem',03,ma_err)
c     
      end
c $Id$
