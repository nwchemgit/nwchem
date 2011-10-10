      subroutine eomccsd_1prdm(rtdb,omega,d_f1,d_v2,
     1                         k_f1_offset,k_v2_offset,
     1                         d_t1,d_t2,d_x1,d_x2,d_y1,d_y2,
     1                         k_t1_offset,k_t2_offset,
     2                         k_x1_offset,k_x2_offset,
     3                         k_y1_offset,k_y2_offset)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "stdio.fh"
      integer rtdb
      double precision omega
      integer d_f1,d_v2
      integer k_f1_offset,k_v2_offset
      integer d_t1,d_t2,d_x1,d_x2,d_y1,d_y2
      integer k_t1_offset,k_t2_offset
      integer k_x1_offset,k_x2_offset
      integer k_y1_offset,k_y2_offset
c
      integer d_e,l_e_offset,k_e_offset,size_e
      double precision r0
      integer d_d0,l_d0_offset,k_d0_offset,size_d0
      double precision denominator
      integer d_x0,l_x0_offset,k_x0_offset,size_x0
c
      character*256 filename
      logical nodezero
      integer d_hh,l_hh_offset,k_hh_offset,size_hh
      integer d_hp,l_hp_offset,k_hp_offset,size_hp
      integer d_ph,l_ph_offset,k_ph_offset,size_ph
      integer d_pp,l_pp_offset,k_pp_offset,size_pp
      integer dim_rdm_ao,l_rdm_ao,k_rdm_ao
      integer dim_mo_h,l_mo_h,k_mo_h,l_mo_h_tmp,k_mo_h_tmp
      integer dim_mo_p,l_mo_p,k_mo_p,l_mo_p_tmp,k_mo_p_tmp
      integer dim_mu_h,l_mu_h,k_mu_h
      integer dim_mu_p,l_mu_p,k_mu_p
      integer dim_rdm_mo_hh,l_rdm_mo_hh,k_rdm_mo_hh
      integer dim_rdm_mo_hp,l_rdm_mo_hp,k_rdm_mo_hp
      integer dim_rdm_mo_ph,l_rdm_mo_ph,k_rdm_mo_ph
      integer dim_rdm_mo_pp,l_rdm_mo_pp,k_rdm_mo_pp
      integer nh,np,i,j,particle,hole
      logical ao_rdm_write
      external ao_rdm_write
      
c
c
      nodezero=(ga_nodeid().eq.0)
c
c     calculate R0
c
      call tce_filename('e',filename)
      call tce_e_offset(l_e_offset,k_e_offset,size_e)
      call createfile(filename,d_e,size_e)
      call tce_zero(d_e,size_e)
      call  nr0(d_f1,d_e,d_t1,d_v2,d_x1,d_x2,k_f1_offset,
     1          k_e_offset,k_t1_offset,k_v2_offset,
     2          k_x1_offset,k_x2_offset)
      call reconcilefile(d_e,1)
      call get_block(d_e,r0,1,0)
      if(dabs(omega).gt.1.0d-7) then
         r0 = r0/omega
      else
         call errquit('eomccsd_1prdm: r0 is infinity',0,ma_err)
      end if
      call deletefile(d_e)
      if(.not.ma_pop_stack(l_e_offset))
     1 call errquit('eomccsd_1prdm: ma problem',1,ma_err)
c
c     calculate denominator
c
      call tce_filename('d_d0',filename)
      call tce_e_offset(l_d0_offset,k_d0_offset,size_d0)
      call createfile(filename,d_d0,size_d0)
      call tce_zero(d_d0,size_d0)
      call eomccsd_denominator(d_d0,d_x1,d_x2,d_y1,d_y2,
     1                k_d0_offset,k_x1_offset,
     2                k_x2_offset,k_y1_offset,k_y2_offset)
      call reconcilefile(d_d0,size_d0)
      call get_block(d_d0,denominator,1,0)
      call deletefile(d_d0)
      if(.not.ma_pop_stack(l_d0_offset))
     1 call errquit('eomccsd_1prdm: ma problem',2,ma_err)
c
c     put r0 value to d_x0
c
      call tce_filename('d_x0',filename)
      call tce_e_offset(l_x0_offset,k_x0_offset,size_x0)
      call createfile(filename,d_x0,size_x0)
      call tce_zero(d_x0,size_x0)
      call put_block(d_x0,r0,1,0)
c
c     allocate memory for 1PRDM with AO basis
c
      dim_rdm_ao=nbf*nbf
      if(.not.ma_push_get(mt_dbl,dim_rdm_ao,'rdm_ao',
     1     l_rdm_ao,k_rdm_ao))
     2     call errquit('eomccsd_1prdm: ma problem',1,ma_err)
      do i=1,dim_rdm_ao
         dbl_mb(k_rdm_ao+i-1)=0.d0
      enddo
c
c     allocate memory for MOs
c
      nh=nocc(1)+nocc(ipol)
      np=nmo(1)+nmo(ipol)-nh
      dim_mo_h=nh*nbf
      dim_mo_p=np*nbf
      if(.not.ma_push_get(mt_dbl,dim_mo_h,'mo_h',
     1 l_mo_h,k_mo_h))
     2 call errquit('eomccsd_1prdm: ma problem',2,ma_err)
      do i=1,dim_mo_h
         dbl_mb(k_mo_h+i-1)=0.d0
      enddo
      if(.not.ma_push_get(mt_dbl,dim_mo_p,'mo_p',
     1 l_mo_p,k_mo_p))
     2 call errquit('eomccsd_1prdm: ma problem',3,ma_err)
      do i=1,dim_mo_p
         dbl_mb(k_mo_p+i-1)=0.d0
      enddo
c
c     get the MOs from GA and make them sorted (list them according column)
c
      if(.not.ma_push_get(mt_dbl,dim_mo_h,'mo_h_tmp',
     1 l_mo_h_tmp,k_mo_h_tmp))
     2 call errquit('eomccsd_1prdm: ma problem',4,ma_err)
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
     1     call errquit('eomccsd_1prdm: ma problem',5,ma_err)
c
c     particle alpha
c
      if(.not.ma_push_get(mt_dbl,dim_mo_p,'mo_p_tmp',
     1     l_mo_p_tmp,k_mo_p_tmp))
     2     call errquit('eomccsd_1prdm: ma problem',6, ma_err)
      do i=1,dim_mo_p
         dbl_mb(k_mo_p_tmp+i-1)=0.d0
      enddo
c
      do particle=nocc(1)+1, nmo(1)
         i=2*particle-1-nh
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
     1     call errquit('eomccsd_1prdm: ma problem',7,ma_err)
c
c     allocate memory for intermediates
c
      dim_mu_h = nbf*nh
      dim_mu_p = nbf*np
      if(.not.ma_push_get(mt_dbl,dim_mu_h,'mu_h',
     1     l_mu_h,k_mu_h))
     2     call errquit('eomccsd_1prdm: ma problem',8, ma_err)
      do i=1,dim_mu_h
         dbl_mb(k_mu_h+i-1)=0.d0
      enddo
      if(.not.ma_push_get(mt_dbl,dim_mu_p,'mu_p',
     1     l_mu_p,k_mu_p))
     2     call errquit('eomccsd_1prdm: ma problem',9, ma_err)
      do i=1, dim_mu_p
         dbl_mb(k_mu_p+i-1)=0.d0
      enddo
c
c-eomccsd_1prdm_hh
c
c     allocate memory for hh block
      dim_rdm_mo_hh = nh*nh
      if(.not.ma_push_get(mt_dbl,dim_rdm_mo_hh,'rdm_mo_hh',
     1     l_rdm_mo_hh,k_rdm_mo_hh))
     2     call errquit('eomccsd_1prdm: ma problem',10,ma_err)
      do i=1,dim_rdm_mo_hh
         dbl_mb(k_rdm_mo_hh+i-1)=0.d0
      enddo
c
      call tce_filename('hh',filename)
      call tce_dens_hh_offset(l_hh_offset,k_hh_offset,size_hh)
      call createfile(filename,d_hh,size_hh)
      call dratoga(d_x1)
      call dratoga(d_x2)
      call dratoga(d_y1)
      call dratoga(d_y2)
      call eomccsd_1prdm_hh(d_hh,d_t1,d_t2,d_x0,d_x1,d_x2,d_y1,d_y2,
     1   k_hh_offset,k_t1_offset,k_t2_offset,k_x0_offset,
     2   k_x1_offset,k_x2_offset,
     3   k_y1_offset,k_y2_offset)
      call reconcilefile(d_hh,size_hh)
      call gatodra(d_y2)
      call gatodra(d_y1)
      call gatodra(d_x2)
      call gatodra(d_x1)
      call get_mo_rdm_hh(d_hh,k_hh_offset,k_rdm_mo_hh,denominator)
      call deletefile(d_hh)
      if (.not.ma_pop_stack(l_hh_offset))
     1 call errquit("eomccsd_1prdm: ma problem",11,ma_err)
c
c    do the matrix multiplication
c
      call dgemm('t','n',nbf,nh,nh,1.d0,dbl_mb(k_mo_h),
     1     nh,dbl_mb(k_rdm_mo_hh),nh,
     2     0.d0,dbl_mb(k_mu_h),nbf)
      call dgemm('n','n',nbf,nbf,nh,1.0d0,dbl_mb(k_mu_h),nbf,
     1     dbl_mb(k_mo_h),nh,0.d0, dbl_mb(k_rdm_ao),nbf)
c
c     release hh block memory
c
      if(.not.ma_pop_stack(l_rdm_mo_hh))
     1     call errquit('ccsd_1prdm: ma problem',12,ma_err)
      do i=1,dim_mu_h
         dbl_mb(k_mu_h+i-1)=0.d0
      enddo
c
c     eomccsd_1prdm_hp
c
c
c     allocate memory for hp block
c
      dim_rdm_mo_hp = nh*np
      if(.not.ma_push_get(mt_dbl,dim_rdm_mo_hp,'rdm_mo_hp',
     1     l_rdm_mo_hp,k_rdm_mo_hp))
     2     call errquit('eomccsd_1prdm: ma problem',13,ma_err)
      do i=1,dim_rdm_mo_hp
         dbl_mb(k_rdm_mo_hp+i-1)=0.d0
      enddo
     
c
c-eomccsd_1prdm_hp
c
      call tce_filename('hp',filename)
      call tce_dens_hp_offset(l_hp_offset,k_hp_offset,size_hp)
      call createfile(filename,d_hp,size_hp)
      call dratoga(d_x1)
      call dratoga(d_x2)
      call dratoga(d_y1)
      call dratoga(d_y2)
      call eomccsd_1prdm_hp(d_hp,d_x0,d_x1,d_y1,d_y2,
     1      k_hp_offset,k_x0_offset,k_x1_offset,
     2      k_y1_offset,k_y2_offset)
      call reconcilefile(d_hp,size_hp)
      call gatodra(d_y2)
      call gatodra(d_y1)
      call gatodra(d_x2)
      call gatodra(d_x1)
      call get_mo_rdm_hp(d_hp,k_hp_offset,k_rdm_mo_hp,denominator)
      call deletefile(d_hp)
      if (.not.ma_pop_stack(l_hp_offset))
     1 call errquit("eomccsd_1prdm: ma problem",14,ma_err)
c
c     do the matrices multiplication
c
      call dgemm('t','n',nbf,np,nh,1.0d0,dbl_mb(k_mo_h),
     1     nh,dbl_mb(k_rdm_mo_hp),nh,
     2     0.d0,dbl_mb(k_mu_p),nbf)
      call dgemm('n','n',nbf,nbf,np,1.0d0,dbl_mb(k_mu_p),nbf,
     1     dbl_mb(k_mo_p),np,1.0d0, dbl_mb(k_rdm_ao),nbf)
c
c     release hp block memory
c
      if(.not.ma_pop_stack(l_rdm_mo_hp))
     1     call errquit('eomccsd_1prdm: ma problem',15,ma_err)
      do i=1,dim_mu_p
         dbl_mb(k_mu_p+i-1)=0.d0
      enddo
c
c     eomccsd_1prdm_ph
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
      call tce_dens_ph_offset(l_ph_offset,k_ph_offset,
     1  size_ph)
      call createfile(filename,d_ph,size_ph)
      call dratoga(d_x1)
      call dratoga(d_x2)
      call dratoga(d_y1)
      call dratoga(d_y2)
      call eomccsd_1prdm_ph(d_ph,d_t1,d_t2,d_x0,d_x1,d_x2,d_y1,d_y2,
     1 k_ph_offset,k_t1_offset,k_t2_offset,
     2 k_x0_offset,k_x1_offset,k_x2_offset,
     2 k_y1_offset,k_y2_offset)
      call reconcilefile(d_ph,size_ph)
      call gatodra(d_y2)
      call gatodra(d_y1)
      call gatodra(d_x2)
      call gatodra(d_x1)
      call get_mo_rdm_ph(d_ph,k_ph_offset,k_rdm_mo_ph,denominator)
      call deletefile(d_ph)
      if(.not.ma_pop_stack(l_ph_offset))
     1 call errquit("eomccsd_1prdm: ma problem",16,ma_err)
c
c
c     do the matrix multiplication
c
      call dgemm('t','n',nbf,nh,np,1.0d0,dbl_mb(k_mo_p),
     1     np,dbl_mb(k_rdm_mo_ph),np,
     2     0.d0,dbl_mb(k_mu_h),nbf)
      call dgemm('n','n',nbf,nbf,nh,1.0d0,dbl_mb(k_mu_h),nbf,
     1     dbl_mb(k_mo_h),nh,1.0d0, dbl_mb(k_rdm_ao),nbf)
c
c     release memory for ph block
c
      if(.not.ma_pop_stack(l_rdm_mo_ph))
     1     call errquit('eomccsd_1prdm: ma problem',17,ma_err)
      do i=1,dim_mu_h
         dbl_mb(k_mu_h+i-1)=0.d0
      enddo
c
c-eomccsd_1prdm_pp
c
c     allocate memory for pp block
c
      dim_rdm_mo_pp=np*np
      if(.not.ma_push_get(mt_dbl,dim_rdm_mo_pp,'rdm_mo_pp',
     1     l_rdm_mo_pp,k_rdm_mo_pp))
     2     call errquit('eomccsd_1prdm: ma problem',18,ma_err)
      do i=1,dim_rdm_mo_pp
         dbl_mb(k_rdm_mo_pp+i-1)=0.d0
      enddo
      call tce_filename('pp',filename)
      call tce_dens_pp_offset(l_pp_offset,k_pp_offset,size_pp)
      call createfile(filename,d_pp,size_pp)
      call dratoga(d_x1)
      call dratoga(d_x2)
      call dratoga(d_y1)
      call dratoga(d_y2)
      call eomccsd_1prdm_pp(d_pp,d_t1,d_t2,d_x0,d_x1,d_x2,
     1  d_y1,d_y2,k_pp_offset,k_t1_offset,k_t2_offset,
     2  k_x0_offset,k_x1_offset,k_x2_offset,
     3  k_y1_offset,k_y2_offset)
      call reconcilefile(d_pp,size_pp)
      call gatodra(d_y2)
      call gatodra(d_y1)
      call gatodra(d_x2)
      call gatodra(d_x1)
      call get_mo_rdm_pp(d_pp,k_pp_offset,k_rdm_mo_pp,denominator)
      call deletefile(d_pp)
      if (.not.ma_pop_stack(l_pp_offset))
     1 call errquit("eomccsd_1prdm: ma problem",19,ma_err)
      call dgemm('t','n',nbf,np,np,1.0d0,dbl_mb(k_mo_p),
     1     np,dbl_mb(k_rdm_mo_pp),np,
     2     0.d0,dbl_mb(k_mu_p),nbf)
      call dgemm('n','n',nbf,nbf,np,1.0d0,dbl_mb(k_mu_p),nbf,
     1     dbl_mb(k_mo_p),np,1.0d0, dbl_mb(k_rdm_ao),nbf)
c
c
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
c
c     release the memory for pp block
c
      if(.not.ma_pop_stack(l_rdm_mo_pp))
     1     call errquit('eomccsd_1prdm: ma problem',20,ma_err)
c
c     release the memory for intermediates
c
      if(.not.ma_pop_stack(l_mu_p))
     1     call errquit('eomccsd_1prdm: ma problem',21,ma_err)
      if(.not.ma_pop_stack(l_mu_h))
     1     call errquit('eomccsd_1prdm: ma problem',22,ma_err)
c
c     release the memory for MOs
c
      if(.not.ma_pop_stack(l_mo_p))
     1     call errquit('eomccsd_1prdm: ma problem',23,ma_err)
      if(.not.ma_pop_stack(l_mo_h))
     1     call errquit('eomccsd_1prdm: ma problem',24,ma_err)
c
c     dump the ao rdm to a file
c
      if (.not.rtdb_cget(rtdb,'tce:file_densmat',1,filename))
     1   call errquit('eomccsd_1prdm: rtdb_cget failed - file_densmat',
     1        0,RTDB_ERR)
      if(.not.ao_rdm_write(filename,k_rdm_ao))
     1     call errquit('eomccsd_1prdm: disk problem',1,disk_err)
      if(.not.ma_pop_stack(l_rdm_ao))
     1     call errquit('eomccsd_1prdm: ma problem',03,ma_err)
c
c
c
      call deletefile(d_x0)
      if(.not.ma_pop_stack(l_x0_offset))
     1     call errquit('eomccsd_1prdm: ma problem',04,ma_err)

      end
c $Id$
