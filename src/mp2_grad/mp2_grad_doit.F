      subroutine mp2_grad_rhf(g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $     irs_a,sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,  
     $     num_oa,num_va,eval_a,basis,geom,tol2e,
     $     e2,nbf,eval_e2,eval_grad)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_vecs_a
      integer noa_lo,noa_hi
      integer irs_a(*)
      integer nva_lo,nva_hi
      integer sym_lo_oa(0:7),sym_hi_oa(0:7)
      integer sym_lo_va(0:7),sym_hi_va(0:7)
      integer num_oa(0:*),num_va(0:*)
      integer basis
      integer geom
      integer nbf
      integer num_i_perseg,npass,ierr
      integer k_amos,l_amos
      integer k_toff_a
      integer l_toff_a
      integer noa,nva
      integer maxireps, maxops
      parameter (maxireps=20, maxops=120)
      integer iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      character*8 group_name
      integer nop,nir
      logical sym_char_table
      external sym_char_table
      integer k_v,l_v
      integer k_temp,l_temp
      integer k_uv_ia,l_uv_ia
      integer k_jb_ia,l_jb_ia
      integer itmp,kunit
      integer tunita
      parameter (tunita=12)
      character*256 fnameta,tmpname
      integer nnbf,ioff
      character*256 fname
      integer nva_lo_local,nva_hi_local
      parameter (kunit=11)
      logical ocoul,oexch
      double precision tol2e,e2
      double precision eval_a(*)
      integer oseg_lo,oseg_hi,pass
      logical eval_e2,eval_grad
      integer k_t_a_in,l_t_a_in
c     
      if(eval_grad)then
         ocoul=.true.
      else
         ocoul=.false.
      endif
      oexch=.true.
c     
      noa=noa_hi-noa_lo+1
      nva=nva_hi-nva_lo+1
c     
      call sym_group_name(geom,group_name)
      if (.not. sym_char_table(group_name, nop, nir, iclass_dim,
     $     zir, zclass, chars))
     $     call errquit('sym_print_char_table: no char table',0)
c     
      num_i_perseg=noa_hi
      npass=0
      do oseg_lo=noa_lo,noa_hi,num_i_perseg
         npass=npass+1
         oseg_hi=min((oseg_lo+num_i_perseg-1),noa_hi)
         call moints_semi(basis,tol2e,ocoul,oexch,oseg_lo,oseg_hi,
     $        noa_lo,noa_hi,nva_lo,nva_hi,g_vecs_a)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'amos',l_amos,k_amos))
     $        call errquit('mp2_grad_doit: error with ma',1)
         if(.not.ma_push_get(mt_int,nva*noa*nir,'toff_a',
     $        l_toff_a,k_toff_a))
     $        call errquit('could not alloc toff_a',1)
         call ga_get(g_vecs_a,1,nbf,1,nbf,dbl_mb(k_amos),nbf)
         call util_file_prefix('khalf',tmpname)
         call util_pname(tmpname,fname)
         call waopen(kunit,fname,1,0,itmp)
         if(itmp.ne.0)call errquit('mp2_grad_doit waopen?',1)
         call moints_vrange(kunit,nva_lo_local,nva_hi_local,nnbf,ioff)
         if(.not.ma_push_get(mt_int,(nnbf+mod(nnbf,2)),'scatter',
     $        l_v,k_v))call errquit('could not alloc k_v',1)
         call moints_getscattv(kunit,nnbf,int_mb(k_v))
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'uv_ia',l_uv_ia,k_uv_ia))
     $        call errquit('could not alloc uv_ia',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'jb_ia',l_jb_ia,k_jb_ia))
     $        call errquit('could not alloc jb_ia',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'temp',l_temp,k_temp))
     $        call errquit('could not alloc temp',1)
         call util_file_prefix('t_alpha',tmpname)
         call util_pname(tmpname,fnameta)
         call waopen(tunita,fnameta,1,0,ierr)
         if(ierr.ne.0)call errquit(
     $        'could not open t_alpha file',1)
         pass=0
         call mp2_make_trhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $        oseg_lo,oseg_hi,basis,irs_a,dbl_mb(k_amos),
     $        dbl_mb(k_uv_ia),dbl_mb(k_jb_ia),eval_a,
     $        int_mb(k_toff_a),nva,noa,nir,
     $        dbl_mb(k_temp),ioff,nnbf,int_mb(k_v),
     $        nva_lo_local,nva_hi_local,
     $        e2,pass,eval_e2,eval_grad,tunita)
         pass=1
         call mp2_make_trhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $        oseg_lo,oseg_hi,basis,irs_a,dbl_mb(k_amos),
     $        dbl_mb(k_uv_ia),dbl_mb(k_jb_ia),eval_a,
     $        int_mb(k_toff_a),nva,noa,nir,
     $        dbl_mb(k_temp),ioff,nnbf,int_mb(k_v),
     $        nva_lo_local,nva_hi_local,e2,pass,
     $        eval_e2,eval_grad,tunita)
         if(.not.ma_chop_stack(l_amos))call errquit(
     $        'chop stack problem',1)
         call waclose(kunit,ierr)
         if(ierr.ne.0)call errquit('problem closing ia_uv file',1)
         if(eval_grad)then
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'t_a_in',
     $        l_t_a_in,k_t_a_in))call errquit(
     $        'could not alloc t_a_in in grad_rhf',1)
         call mp2_read_trhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,num_oa,
     $     oseg_lo,oseg_hi,irs_a,nva_lo_local,nva_hi_local,
     $     dbl_mb(k_t_a_in),tunita)
         if(.not.ma_chop_stack(l_t_a_in))call errquit(
     $        'could not chop stack at end of grad_rhf',1)
         call waclose(tunita,ierr)
         if(ierr.ne.0)call errquit('problem closing tunita',1)
         endif
      enddo
      end
c     
c     
      subroutine mp2_grad_uhf(g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $     irs_a,sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,  
     $     num_oa,num_va,
     $     g_vecs_b,nob_lo,nob_hi,nvb_lo,nvb_hi,
     $     irs_b,sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $     num_ob,num_vb,
     $     eval_a,eval_b,
     $     basis,geom,tol2e,e2,nbf,eval_e2,eval_grad)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_vecs_a,g_vecs_b
      integer noa_lo,noa_hi,nob_lo,nob_hi
      integer irs_a(*),irs_b(*)
      integer nva_lo,nva_hi,nvb_lo,nvb_hi
      integer sym_lo_oa(0:7),sym_hi_oa(0:7)
      integer sym_lo_va(0:7),sym_hi_va(0:7)
      integer sym_lo_ob(0:7),sym_hi_ob(0:7)
      integer sym_lo_vb(0:7),sym_hi_vb(0:7)
      integer num_oa(0:*),num_va(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer basis
      integer geom
      integer nbf
      integer num_i_perseg,npass,ierr
      integer k_amos,k_bmos,l_amos,l_bmos
      integer k_toff_a,k_toff_b
      integer l_toff_a,l_toff_b
      integer noa,nva,nob,nvb
      integer maxireps, maxops
      parameter (maxireps=20, maxops=120)
      integer iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      character*8 group_name
      integer nop,nir
      logical sym_char_table
      external sym_char_table
      integer k_v,l_v
      integer k_temp,l_temp
      integer k_uv_ia,l_uv_ia
      integer k_jb_ia,l_jb_ia
      integer itmp,kunit
      integer nnbf,ioff
      character*256 fname,tmpname
      character*256 fnameta,fnametb
      integer nva_lo_local,nva_hi_local
      parameter(kunit=11)
      integer tunita,tunitb
      parameter(tunita=13)
      parameter(tunitb=14)
      logical ocoul,oexch
      double precision tol2e,e2
      double precision eval_a(*),eval_b(*)
      integer oseg_lo,oseg_hi
      logical eval_e2,eval_grad
      integer k_t_a_in,k_t_b_in,l_t_a_in,l_t_b_in
c
      if(eval_grad)then
         ocoul=.true.
      else
         ocoul=.false.
      endif
      oexch=.true.
c     
      noa=noa_hi-noa_lo+1
      nva=nva_hi-nva_lo+1
      nob=nob_hi-nob_lo+1
      nvb=nvb_hi-nvb_lo+1
c     
      call sym_group_name(geom,group_name)
      if (.not. sym_char_table(group_name, nop, nir, iclass_dim,
     $     zir, zclass, chars))
     $     call errquit('sym_print_char_table: no char table',0)
c     
      num_i_perseg=noa_hi
      npass=0
      do oseg_lo=noa_lo,noa_hi,num_i_perseg
         npass=npass+1
         oseg_hi=min((oseg_lo+num_i_perseg-1),noa_hi)
         call moints_semi(basis,tol2e,ocoul,oexch,oseg_lo,oseg_hi,
     $        noa_lo,noa_hi,nva_lo,nva_hi,g_vecs_a)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'amos',l_amos,k_amos))
     $        call errquit('mp2_grad_doit: error with ma',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'bmos',l_bmos,k_bmos))
     $        call errquit('mp2_grad_doit: error with ma',1)
         if(.not.ma_push_get(mt_int,nva*noa*nir,'toff_a',
     $        l_toff_a,k_toff_a))
     $        call errquit('could not alloc toff_a',1)
         if(.not.ma_push_get(mt_int,nvb*nob*nir,'toff_b',
     $        l_toff_b,k_toff_b))
     $        call errquit('could not alloc toff_a',1)
         call ga_get(g_vecs_a,1,nbf,1,nbf,dbl_mb(k_amos),nbf)
         call ga_get(g_vecs_b,1,nbf,1,nbf,dbl_mb(k_bmos),nbf)
         call util_file_prefix('khalf',tmpname)
         call util_pname(tmpname,fname)
         call waopen(kunit,fname,1,0,itmp)
         if(itmp.ne.0)call errquit('mp2_grad_doit waopen?',1)
         call moints_vrange(kunit,nva_lo_local,nva_hi_local,nnbf,ioff)
         if(.not.ma_push_get(mt_int,(nnbf+mod(nnbf,2)),'scatter',
     $        l_v,k_v))call errquit('could not alloc k_v',1)
         call moints_getscattv(kunit,nnbf,int_mb(k_v))
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'uv_ia',l_uv_ia,k_uv_ia))
     $        call errquit('could not alloc uv_ia',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'jb_ia',l_jb_ia,k_jb_ia))
     $        call errquit('could not alloc jb_ia',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'temp',l_temp,k_temp))
     $        call errquit('could not alloc temp',1)
         call util_file_prefix('t_alpha',tmpname)
         call util_pname(tmpname,fnameta)
         call waopen(tunita,fnameta,1,0,ierr)
         if(ierr.ne.0)call errquit(
     $        'could not open t_alpha file',1)
         call util_file_prefix('t_beta',tmpname)
         call util_pname(tmpname,fnametb)
         call waopen(tunitb,fnametb,1,0,ierr)
         if(ierr.ne.0)call errquit(
     $        'could not open t_beta file',1)
         call mp2_make_tuhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $        nob_lo,nob_hi,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $        sym_lo_vb,sym_hi_vb,num_ob,num_vb,oseg_lo,oseg_hi,
     $        basis,irs_a,dbl_mb(k_amos),irs_b,dbl_mb(k_bmos),
     $        dbl_mb(k_uv_ia),dbl_mb(k_jb_ia),
     $        eval_a,eval_b,int_mb(k_toff_a),int_mb(k_toff_b),
     $        nva,noa,nvb,nob,nir,dbl_mb(k_temp),ioff,nnbf,
     $        int_mb(k_v),nva_lo_local,nva_hi_local,
     $        e2,eval_e2,eval_grad,tunita,tunitb)
         if(.not.ma_chop_stack(l_amos))call errquit(
     $        'chop stack problem',1)
         call waclose(kunit,ierr)
         if(ierr.ne.0)call errquit('problem closing ia_uv file',1)
         if(eval_grad)then
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'t_a_in',
     $        l_t_a_in,k_t_a_in))call errquit(
     $        'could not alloc t_a_in in grad_uhf',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'t_b_in',
     $        l_t_b_in,k_t_b_in))call errquit(
     $        'could not alloc t_b_in in grad_uhf',1)
         call mp2_read_tuhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,num_oa,
     $        nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,num_ob,oseg_lo,oseg_hi,
     $        irs_a,irs_b,nva_lo_local,nva_hi_local,dbl_mb(k_t_a_in),
     $        dbl_mb(k_t_b_in),tunita,tunitb)
         call waclose(tunitb,ierr)
         if(ierr.ne.0)call errquit('problem closing t_b file',ierr)
         call waclose(tunita,ierr)
         if(ierr.ne.0)call errquit('problem closing t_a file',ierr)
         if(.not.ma_chop_stack(l_t_a_in))call errquit(
     $        'chop stack problem',1)
         endif
      enddo
      end
         

