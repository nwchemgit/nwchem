      subroutine mp2_grad_rhf(g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $     irs_a,sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,  
     $     num_oa,num_va,eval_a,basis,geom,tol2e,
     $     e2,nbf,task)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_vecs_a
      character*8 task
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
      parameter(maxireps=20,maxops=120)
      integer iclass_dim(maxireps)
      character*8 zir(maxireps),zclass(maxireps)
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
      parameter(tunita=12)
      character*256 fnameta,tmpname
      integer nnbf,ioff
      character*256 fname
      integer nva_lo_local,nva_hi_local
      parameter(kunit=11)
      logical ocoul,oexch
      double precision tol2e,e2
      double precision eval_a(*)
      integer oseg_lo,oseg_hi,pass
      integer k_t_a_in,l_t_a_in
c     
      if(task.eq.'gradient')then
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
      if(.not.sym_char_table(group_name,nop,nir,iclass_dim,
     $     zir,zclass,chars))
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
         call util_file_name('khalf',.true.,.true.,fname)
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
c         if(task.eq.'gradient')then
c$$$  call util_file_prefix('t_alpha',tmpname)
c$$$  call util_pname(tmpname,fnameta)
c$$$  call waopen(tunita,fnameta,1,0,ierr)
c$$$  if(ierr.ne.0)call errquit(
c$$$  $           'could not open t_alpha file',1)
            call util_file_name('t_alpha',.true.,.true.,fnameta)
            call waopen(tunita,fnameta,1,0,ierr)
            if(ierr.ne.0)call errquit(
     $           'could not open t_alpha file',1)
c         endif
         pass=0
         call mp2_make_trhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $        oseg_lo,oseg_hi,basis,irs_a,dbl_mb(k_amos),
     $        dbl_mb(k_uv_ia),dbl_mb(k_jb_ia),eval_a,
     $        int_mb(k_toff_a),nva,noa,nir,
     $        dbl_mb(k_temp),ioff,nnbf,int_mb(k_v),
     $        nva_lo_local,nva_hi_local,
     $        e2,pass,task,tunita)
         pass=1
         call mp2_make_trhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $        oseg_lo,oseg_hi,basis,irs_a,dbl_mb(k_amos),
     $        dbl_mb(k_uv_ia),dbl_mb(k_jb_ia),eval_a,
     $        int_mb(k_toff_a),nva,noa,nir,
     $        dbl_mb(k_temp),ioff,nnbf,int_mb(k_v),
     $        nva_lo_local,nva_hi_local,e2,pass,
     $        task,tunita)
         if(.not.ma_chop_stack(l_amos))call errquit(
     $        'chop stack problem',1)
c         if(task.eq.'gradient')then
            if(.not.ma_push_get(mt_dbl,noa*nva,'t_a_in',
     $           l_t_a_in,k_t_a_in))call errquit(
     $           'could not alloc t_a_in in grad_rhf',1)
            call mp2_read_trhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $           num_oa,oseg_lo,oseg_hi,irs_a,nva_lo_local,
     $           nva_hi_local,dbl_mb(k_t_a_in),tunita)
            if(.not.ma_chop_stack(l_t_a_in))call errquit(
     $           'could not chop stack at end of grad_rhf',1)
            call waclose(kunit,ierr)
            if(ierr.ne.0)call errquit('problem closing ia_uv file',1)
            call waclose(tunita,ierr)
            if(ierr.ne.0)call errquit('problem closing tunita',1)
c         endif
      enddo
      end
c     
c     
      subroutine mp2_grad_uhf(g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $     irs_a,sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,  
     $     num_oa,num_va,g_vecs_b,nob_lo,nob_hi,nvb_lo,nvb_hi,
     $     irs_b,sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $     num_ob,num_vb,eval_a,eval_b,basis,geom,tol2e,e2,nbf,
     $     task,pass,e2back,oblk,p_ij_a,p_ij_b,
     $     p_ab_a,p_ab_b,g_p_ij_a,g_p_ij_b,g_p_ab_a,
     $     g_p_ab_b,g_lai_a,g_lai_b)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_vecs_a,g_vecs_b,pass
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
      integer k_ioff_a,k_ioff_b
      integer l_ioff_a,l_ioff_b
      integer k_ioff_vooo_a,l_ioff_vooo_a
      integer k_ioff_vooo_b,l_ioff_vooo_b
      integer noa,nva,nob,nvb
      integer maxireps, maxops
      parameter(maxireps=20,maxops=120)
      integer iclass_dim(maxireps)
      character*8 zir(maxireps),zclass(maxireps)
      double precision chars(maxireps*maxireps)
      character*8 group_name
      integer nop,nir
      logical sym_char_table
      external sym_char_table
      integer k_v,l_v
      integer k_temp,l_temp
      integer k_uv_ia,l_uv_ia
      integer k_jb_ia,l_jb_ia
      integer itmp,kunit,junit
      integer nnbf,ioff
      character*256 fname,tmpname
      character*256 fnameta1,fnametb1
      character*256 fnameta2,fnametb2
      integer nva_lo_local,nva_hi_local
      integer noa_lo_local,noa_hi_local
      parameter(kunit=11)
      parameter(junit=10)
      integer tunita1,tunitb1,tunita2,tunitb2
      parameter(tunita1=13)
      parameter(tunitb1=14)
      parameter(tunita2=15)
      parameter(tunitb2=16)
      logical ocoul,oexch
      double precision tol2e,e2,eaa,eab,e2back
      double precision eval_a(*),eval_b(*)
      integer oseg_lo,oseg_hi
      character*8 task
      integer k_t_a_in,k_t_b_in,l_t_a_in,l_t_b_in
      integer k_t_eps_a,k_t_eps_b,l_t_eps_a,l_t_eps_b
      double precision p_ij_a(*),p_ij_b(*),p_ab_a(*),p_ab_b(*)
c$$$      integer k_p_ij_a,k_p_ij_b
c$$$      integer l_p_ij_a,l_p_ij_b
c$$$      integer k_p_ab_a,k_p_ab_b
c$$$      integer l_p_ab_a,l_p_ab_b
      integer k_w_ij_a,k_w_ij_b
      integer l_w_ij_a,l_w_ij_b
      integer k_w_ab_a,k_w_ab_b
      integer l_w_ab_a,l_w_ab_b
      integer k_lai_a,k_lai_b
      integer l_lai_a,l_lai_b
      integer num_i,num_a
      integer g_vvvo_a,g_vvvo_b
      integer g_vvoo_a,g_vvoo_b
      integer g_iajk_a,g_iajk_b
      integer g_lai_a,g_lai_b
      integer k_vvvo,l_vvvo
      integer k_vvvo_off_a,l_vvvo_off_a
      integer k_vvvo_off_b,l_vvvo_off_b,xx
      integer k_p_ai_a,l_p_ai_a,k_p_ai_b,l_p_ai_b
      integer k_w_ai_a,l_w_ai_a,k_w_ai_b,l_w_ai_b
      integer g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b
c     
c...  for texas blocking in moints_semi
c     
      logical oblk
c     
c$$$      write(6,*) 'inside mp2_grad_doit'
c$$$      write(6,*) 'pass',pass
c$$$      call ga_print(g_p_ij_a)
c$$$      call ga_print(g_p_ij_b)
c$$$      call ga_print(g_p_ab_a)
c$$$      call ga_print(g_p_ab_b)
c$$$      call ga_print(g_lai_a)
c$$$      call ga_print(g_lai_b)
c     
      if(task.eq.'gradient')then
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
      if(.not.sym_char_table(group_name,nop,nir,iclass_dim,
     $     zir,zclass,chars))
     $     call errquit('sym_print_char_table: no char table',0)
c     
      num_i_perseg=noa_hi
      npass=0
      do oseg_lo=noa_lo,noa_hi,num_i_perseg
         npass=npass+1
         oseg_hi=min((oseg_lo+num_i_perseg-1),noa_hi)
         call moints_semi(basis,tol2e,ocoul,oexch,oseg_lo,oseg_hi,
     $        noa_lo,noa_hi,nva_lo,nva_hi,g_vecs_a,oblk)
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
         call util_file_name('khalf',.true.,.true.,fname)
         call waopen(kunit,fname,1,0,itmp)
         if(itmp.ne.0)call errquit('mp2_grad_doit waopen?',1)
         call moints_vrange(kunit,nva_lo_local,nva_hi_local,nnbf,ioff)
         if(task.eq.'gradient')then
            call util_file_prefix('jhalf',tmpname)
            call util_pname(tmpname,fname)
            call waopen(junit,fname,1,0,itmp)
            if(itmp.ne.0)call errquit('mp2_grad_doit waopen?',1)
            call moints_vrange(junit,noa_lo_local,
     $           noa_hi_local,nnbf,ioff)
         endif
         if(.not.ma_push_get(mt_int,(nnbf+mod(nnbf,2)),'scatter',
     $        l_v,k_v))call errquit('could not alloc k_v',1)
         call moints_getscattv(kunit,nnbf,int_mb(k_v))
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'jb_ia',l_jb_ia,k_jb_ia))
     $        call errquit('could not alloc jb_ia',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'temp',l_temp,k_temp))
     $        call errquit('could not alloc temp',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'uv_ia',l_uv_ia,k_uv_ia))
     $        call errquit('could not alloc uv_ia',1)
c         if(task.eq.'energy')then
            if(pass.eq.1)then
               call util_file_name('t_alpha1',.true.,.true.,fnameta1)
               call waopen(tunita1,fnameta1,1,0,ierr)
               if(ierr.ne.0)call errquit(
     $              'could not open t_alpha1 file',1)
               call util_file_name('t_beta1',.true.,.true.,fnametb1)
               call waopen(tunitb1,fnametb1,1,0,ierr)
               if(ierr.ne.0)call errquit(
     $              'could not open t_beta1 file',1)
            else
               call util_file_name('t_alpha2',.true.,.true.,fnameta2)
               call waopen(tunita2,fnameta2,1,0,ierr)
               if(ierr.ne.0)call errquit(
     $              'could not open t_alpha2 file',1)
               call util_file_name('t_beta2',.true.,.true.,fnametb2)
               call waopen(tunitb2,fnametb2,1,0,ierr)
               if(ierr.ne.0)call errquit(
     $              'could not open t_beta2 file',1)
            endif
c         endif
         num_i=oseg_hi-oseg_lo+1
         num_a=nva_hi_local-nva_lo_local+1
         if(pass.eq.1)then
            call mp2_make_tuhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $           sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $           nob_lo,nob_hi,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $           sym_lo_vb,sym_hi_vb,num_ob,num_vb,oseg_lo,oseg_hi,
     $           irs_a,dbl_mb(k_amos),irs_b,dbl_mb(k_bmos),
     $           dbl_mb(k_uv_ia),dbl_mb(k_jb_ia),
     $           eval_a,eval_b,int_mb(k_toff_a),int_mb(k_toff_b),
     $           nva,noa,nvb,nob,nir,dbl_mb(k_temp),ioff,nnbf,
     $           int_mb(k_v),nva_lo_local,nva_hi_local,
     $           e2,task,tunita1,tunitb1)
         else
            call mp2_make_tuhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $           sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $           nob_lo,nob_hi,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $           sym_lo_vb,sym_hi_vb,num_ob,num_vb,oseg_lo,oseg_hi,
     $           irs_a,dbl_mb(k_amos),irs_b,dbl_mb(k_bmos),
     $           dbl_mb(k_uv_ia),dbl_mb(k_jb_ia),
     $           eval_a,eval_b,int_mb(k_toff_a),int_mb(k_toff_b),
     $           nva,noa,nvb,nob,nir,dbl_mb(k_temp),ioff,nnbf,
     $           int_mb(k_v),nva_lo_local,nva_hi_local,
     $           e2,task,tunita2,tunitb2)
         endif
         if(task.eq.'gradient')then
c$$$            if(.not.ma_push_get(mt_dbl,noa*noa,'p_ij_a',
c$$$     $           l_p_ij_a,k_p_ij_a))call errquit(
c$$$     $           'could not alloc p_ij_a',1)
c$$$            if(.not.ma_push_get(mt_dbl,nva*nva,'p_ab_a',
c$$$     $           l_p_ab_a,k_p_ab_a))call errquit(
c$$$     $           'could not alloc p_ab_a',1)
c$$$            if(.not.ma_push_get(mt_dbl,nob*nob,'p_ij_b',
c$$$     $           l_p_ij_b,k_p_ij_b))call errquit(
c$$$     $           'could not alloc p_ij_b',1)
c$$$            if(.not.ma_push_get(mt_dbl,nvb*nvb,'p_ab_b',
c$$$     $           l_p_ab_b,k_p_ab_b))call errquit(
c$$$     $           'could not alloc p_ab_b',1)
            if(.not.ma_push_get(mt_dbl,noa*nva,
     $           't_a_in',l_t_a_in,k_t_a_in))call errquit(
     $           'could not alloc t_a_in in grad_uhf',1)
            if(.not.ma_push_get(mt_dbl,nob*nvb,
     $           't_b_in',l_t_b_in,k_t_b_in))call errquit(
     $           'could not alloc t_b_in in grad_uhf',1)
            if(.not.ma_push_get(mt_dbl,noa*nva,
     $           't_eps_a',l_t_eps_a,k_t_eps_a))call errquit(
     $           'could not alloc t_eps_a in grad_uhf',1)
            if(.not.ma_push_get(mt_dbl,nob*nvb,
     $           't_eps_b',l_t_eps_b,k_t_eps_b))call errquit(
     $           'could not alloc t_eps_b in grad_uhf',1)
            if(.not.ma_push_get(mt_dbl,noa*noa,'w_ij_a',
     $           l_w_ij_a,k_w_ij_a))call errquit(
     $           'could not alloc w_ij_a',1)
            if(.not.ma_push_get(mt_dbl,nva*nva,'w_ab_a',
     $           l_w_ab_a,k_w_ab_a))call errquit(
     $           'could not alloc w_ab_a',1)
            if(.not.ma_push_get(mt_dbl,nob*nob,'w_ij_b',
     $           l_w_ij_b,k_w_ij_b))call errquit(
     $           'could not alloc w_ij_b',1)
            if(.not.ma_push_get(mt_dbl,nvb*nvb,'w_ab_b',
     $           l_w_ab_b,k_w_ab_b))call errquit(
     $           'could not alloc w_ab_b',1)
            if(pass.eq.1)then
c$$$               call mp2_pijab_uhf(nva_lo,nva_hi,sym_lo_oa,
c$$$     $              sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
c$$$     $              sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
c$$$     $              nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
c$$$     $              nva_lo_local,nva_hi_local,noa_lo,noa_hi,
c$$$     $              nob_lo,nob_hi,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
c$$$     $              tunita1,tunitb1,p_ij_a,p_ij_b,p_ab_a,p_ab_b,
c$$$     $              eval_a,eval_b,e2back)
               call mp2_wijab_uhf(nir,nva_lo,nva_hi,
     $              sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,sym_lo_vb,
     $              sym_hi_vb,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $              num_oa,num_ob,num_va,num_vb,
     $              oseg_lo,oseg_hi,irs_a,irs_b,
     $              nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $              nob_lo,nob_hi,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
     $              dbl_mb(k_t_eps_a),dbl_mb(k_t_eps_b),
     $              tunita1,tunitb1,dbl_mb(k_w_ij_a),dbl_mb(k_w_ab_a),
     $              dbl_mb(k_w_ij_b),dbl_mb(k_w_ab_b),p_ij_a,
     $              p_ab_a,p_ij_b,p_ab_b,
     $              eval_a,eval_b)
            endif
            if(pass.eq.2)then
c$$$               call mp2_pijab_uhf(nva_lo,nva_hi,sym_lo_oa,
c$$$     $              sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
c$$$     $              sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
c$$$     $              nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
c$$$     $              nva_lo_local,nva_hi_local,noa_lo,noa_hi,
c$$$     $              nob_lo,nob_hi,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
c$$$     $              tunita2,tunitb2,p_ij_a,p_ij_b,p_ab_a,p_ab_b,
c$$$     $              eval_a,eval_b,e2back)
c$$$               if(ga_nodeid().eq.0)then
c$$$                  call ga_put(g_p_ij_b,1,nob,1,nob,
c$$$     $                 dbl_mb(k_p_ij_b),noa)
c$$$                  call ga_put(g_p_ij_a,1,noa,1,noa,
c$$$     $                 dbl_mb(k_p_ij_a),noa)
c$$$                  call ga_put(g_p_ab_b,1,nvb,1,nvb,
c$$$     $                 dbl_mb(k_p_ab_b),nvb)
c$$$                  call ga_put(g_p_ab_a,1,nva,1,nva,
c$$$     $                 dbl_mb(k_p_ab_a),nva)
c$$$               endif
               call mp2_wijab_uhf(nir,nva_lo,nva_hi,
     $              sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,sym_lo_vb,
     $              sym_hi_vb,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $              num_oa,num_ob,num_va,num_vb,
     $              oseg_lo,oseg_hi,irs_a,irs_b,
     $              nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $              nob_lo,nob_hi,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
     $              dbl_mb(k_t_eps_a),dbl_mb(k_t_eps_b),
     $              tunita2,tunitb2,dbl_mb(k_w_ij_a),dbl_mb(k_w_ab_a),
     $              dbl_mb(k_w_ij_b),dbl_mb(k_w_ab_b),p_ij_a,
     $              p_ab_a,p_ij_b,p_ab_b,
     $              eval_a,eval_b)
            endif
            return
            if(.not.ga_create(mt_dbl,nva*nva,noa,'vvoo_a',
     $           nva*nva,1,g_vvoo_a))call errquit(
     $           'could not alloc vvoo_a',1)
            if(.not.ga_create(mt_dbl,nvb*nvb,noa,'vvoo_b',
     $           nvb*nvb,1,g_vvoo_b))call errquit(
     $           'could not alloc vvoo_b',1)
            if(.not.ma_push_get(mt_int,nva*nva*nir,'ioff_a',
     $           l_ioff_a,k_ioff_a))
     $           call errquit('could not alloc ioff_a',1)
            if(.not.ma_push_get(mt_int,nvb*nvb*nir,'ioff_b',
     $           l_ioff_b,k_ioff_b))
     $           call errquit('could not alloc ioff_b',1)
            call mp2_vvoo_uhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $           sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,
     $           num_va,nob_lo,nob_hi,nvb_lo,nvb_hi,sym_lo_ob,
     $           sym_hi_ob,sym_lo_vb,sym_hi_vb,num_ob,num_vb,
     $           oseg_lo,oseg_hi,basis,irs_a,dbl_mb(k_amos),
     $           irs_b,dbl_mb(k_bmos),dbl_mb(k_uv_ia),
     $           dbl_mb(k_jb_ia),int_mb(k_ioff_a),int_mb(k_ioff_b),
     $           nva,noa,nvb,nob,nir,dbl_mb(k_temp),ioff,nnbf,
     $           int_mb(k_v),noa_lo_local,noa_hi_local,
     $           g_vvoo_a,g_vvoo_b)
            if(.not.ma_push_get(mt_int,noa*noa*nir,'ioff_vooo_a',
     $           l_ioff_vooo_a,k_ioff_vooo_a))
     $           call errquit('could not alloc ioff_vooo_a',1)
            if(.not.ma_push_get(mt_int,nob*nob*nir,'ioff_vooo_b',
     $           l_ioff_vooo_b,k_ioff_vooo_b))
     $           call errquit('could not alloc ioff_vooo_a',1)
            if(.not.ga_create(mt_dbl,noa*noa,nva,'iajk_a',
     $           (noa*noa),1,g_iajk_a))
     $           call errquit('could not alloc ga for iajk_a',1)
            if(.not.ga_create(mt_dbl,nob*nob,nva,'iajk_b',
     $           (nob*nob),1,g_iajk_b))
     $           call errquit('could not alloc ga for iajk_b',1)
            call mp2_vooo_uhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $           sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,
     $           num_va,nob_lo,nob_hi,nvb_lo,nvb_hi,sym_lo_ob,
     $           sym_hi_ob,sym_lo_vb,sym_hi_vb,num_ob,num_vb,
     $           oseg_lo,oseg_hi,basis,irs_a,dbl_mb(k_amos),
     $           irs_b,dbl_mb(k_bmos),dbl_mb(k_uv_ia),
     $           int_mb(k_ioff_vooo_a),int_mb(k_ioff_vooo_b),
     $           nva,noa,nvb,nob,nir,dbl_mb(k_temp),ioff,nnbf,
     $           int_mb(k_v),noa_lo_local,noa_hi_local,
     $           nva_lo_local,nva_hi_local,g_iajk_a,g_iajk_b)
            if(.not.ga_create(mt_dbl,nva*nva,nva,'vvvo_a',(nva*nva),
     $           1,g_vvvo_a))
     $           call errquit('could not alloc ga for vvvo_a',1)
            if(.not.ga_create(mt_dbl,nvb*nvb,nva,'vvvo_b',(nvb*nvb),
     $           1,g_vvvo_b))
     $           call errquit('could not alloc ga for vvvo_b',1) 
            if(.not.ma_push_get(mt_dbl,nbf*nbf,'vvvo',l_vvvo,k_vvvo))
     $           call errquit('could not alloc vvvo',1)
            if(.not.ma_push_get(mt_int,nva*nva*nir,'vvvo_off_a',
     $           l_vvvo_off_a,k_vvvo_off_a))
     $           call errquit('could not alloc vvvo_off_a',1)
            if(.not.ma_push_get(mt_int,nvb*nvb*nir,'vvvo_off_b',
     $           l_vvvo_off_b,k_vvvo_off_b))
     $           call errquit('could not alloc vvvo_off_b',1)
            call mp2_vvvo_uhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $           sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,
     $           num_va,nob_lo,nob_hi,nvb_lo,nvb_hi, sym_lo_ob,
     $           sym_hi_ob,sym_lo_vb,sym_hi_vb,num_ob,num_vb,
     $           oseg_lo,oseg_hi,basis,irs_a,dbl_mb(k_amos),irs_b,
     $           dbl_mb(k_bmos),dbl_mb(k_uv_ia),dbl_mb(k_vvvo),
     $           int_mb(k_vvvo_off_a),int_mb(k_vvvo_off_b),nva,
     $           noa,nvb,nob,nir,g_vvvo_a,g_vvvo_b,dbl_mb(k_temp),
     $           ioff,nnbf,int_mb(k_v),nva_lo_local,nva_hi_local)
            if(.not.ma_push_get(mt_dbl,noa*nva,'lai_a',
     $           l_lai_a,k_lai_a))call errquit(
     $           'could not alloc lai_a',1)
            if(.not.ma_push_get(mt_dbl,nob*nvb,'lai_b',
     $           l_lai_b,k_lai_b))call errquit(
     $           'could not alloc lai_b',1)
            if(pass.eq.1)then
c$$$               call mp2_lai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
c$$$     $              sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
c$$$     $              sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
c$$$     $              nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
c$$$     $              nva_lo_local,nva_hi_local,noa_lo_local,
c$$$     $              noa_hi_local,noa_lo,noa_hi,nob_lo,nob_hi,nva,nvb,
c$$$     $              noa,nob,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
c$$$     $              dbl_mb(k_vvvo),dbl_mb(k_jb_ia),
c$$$     $              int_mb(k_vvvo_off_a),int_mb(k_vvvo_off_b),
c$$$     $              g_vvvo_a,g_vvvo_b,g_iajk_a,g_iajk_b,
c$$$     $              int_mb(k_ioff_vooo_a),int_mb(k_ioff_vooo_b),
c$$$     $              tunita1,tunitb1,
c$$$     $              dbl_mb(k_p_ij_a),dbl_mb(k_p_ij_b),
c$$$     $              dbl_mb(k_p_ab_a),dbl_mb(k_p_ab_b),
c$$$     $              dbl_mb(k_lai_a),dbl_mb(k_lai_b),
c$$$     $              g_lai_a,g_lai_b,g_p_ij_a,g_p_ij_b,
c$$$     $              g_p_ab_a,g_p_ab_b)
c$$$            else
c$$$               call mp2_lai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
c$$$     $              sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
c$$$     $              sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
c$$$     $              nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
c$$$     $              nva_lo_local,nva_hi_local,noa_lo_local,
c$$$     $              noa_hi_local,noa_lo,noa_hi,nob_lo,nob_hi,nva,nvb,
c$$$     $              noa,nob,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
c$$$     $              dbl_mb(k_vvvo),dbl_mb(k_jb_ia),
c$$$     $              int_mb(k_vvvo_off_a),int_mb(k_vvvo_off_b),
c$$$     $              g_vvvo_a,g_vvvo_b,g_iajk_a,g_iajk_b,
c$$$     $              int_mb(k_ioff_vooo_a),int_mb(k_ioff_vooo_b),
c$$$     $              tunita2,tunitb2,
c$$$     $              dbl_mb(k_p_ij_a),dbl_mb(k_p_ij_b),
c$$$     $              dbl_mb(k_p_ab_a),dbl_mb(k_p_ab_b),
c$$$     $              dbl_mb(k_lai_a),dbl_mb(k_lai_b),
c$$$     $              g_lai_a,g_lai_b,g_p_ij_a,g_p_ij_b,
c$$$     $              g_p_ab_a,g_p_ab_b)
            endif
            if(.not.ma_push_get(mt_dbl,noa*nva,'p_ai_a',
     $           l_p_ai_a,k_p_ai_a))call errquit(
     $           'could not alloc p_ai_a',1)
            if(.not.ma_push_get(mt_dbl,nob*nvb,'p_ai_b',
     $           l_p_ai_b,k_p_ai_b))call errquit(
     $           'could not alloc p_ai_b',1)
            if(.not.ma_push_get(mt_dbl,noa*nva,'w_ai_a',
     $           l_w_ai_a,k_w_ai_a))call errquit(
     $           'could not alloc w_ai_a',1)
            if(.not.ma_push_get(mt_dbl,nob*nvb,'w_ai_b',
     $           l_w_ai_b,k_w_ai_b))call errquit(
     $           'could not alloc w_ai_b',1)
            if(pass.eq.1)then
               call mp2_wai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $              sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $              sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $              nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $              nva_lo_local,nva_hi_local,noa_lo,noa_hi,nob_lo,
     $              nob_hi,nva,nvb,noa,nob,dbl_mb(k_t_a_in),
     $              dbl_mb(k_t_b_in),dbl_mb(k_jb_ia),
     $              g_iajk_a,g_iajk_b,
     $              int_mb(k_ioff_vooo_a),int_mb(k_ioff_vooo_b),
     $              tunita1,tunitb1,
     $              dbl_mb(k_p_ai_a),dbl_mb(k_p_ai_b),
     $              dbl_mb(k_w_ai_a),dbl_mb(k_w_ai_b))
            else
               call mp2_wai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $              sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $              sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $              nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $              nva_lo_local,nva_hi_local,noa_lo,noa_hi,nob_lo,
     $              nob_hi,nva,nvb,noa,nob,dbl_mb(k_t_a_in),
     $              dbl_mb(k_t_b_in),dbl_mb(k_jb_ia),
     $              g_iajk_a,g_iajk_b,
     $              int_mb(k_ioff_vooo_a),int_mb(k_ioff_vooo_b),
     $              tunita2,tunitb2,
     $              dbl_mb(k_p_ai_a),dbl_mb(k_p_ai_b),
     $              dbl_mb(k_w_ai_a),dbl_mb(k_w_ai_b))
            endif
            if(.not.ga_destroy(g_vvvo_b))call errquit(
     $           'could not destroy vvvo_b handle',1)
            if(.not.ga_destroy(g_vvvo_a))call errquit(
     $           'could not destroy vvvo_a handle',1)
            if(.not.ga_destroy(g_iajk_b))call errquit(
     $           'could not destroy iajk_b handle',1)
            if(.not.ga_destroy(g_iajk_a))call errquit(
     $           'could not destroy iajk_a handle',1)
            if(.not.ga_destroy(g_vvoo_b))call errquit(
     $           'could not destroy vvoo_b handle',1)
            if(.not.ga_destroy(g_vvoo_a))call errquit(
     $           'could not destroy vvoo_a handle',1)
c     
            call waclose(kunit,ierr)
            if(ierr.ne.0)call errquit('problem closing ia_uv file',1)
            if(pass.eq.1)then
               call waclose(tunita1,ierr)
               if(ierr.ne.0)call errquit('problem closing t_a1 file',
     $              ierr)
               call waclose(tunitb1,ierr)
               if(ierr.ne.0)call errquit('problem closing t_b1 file',
     $              ierr)
            else
               call waclose(tunitb2,ierr)
               if(ierr.ne.0)call errquit('problem closing t_b2 file',
     $              ierr)
               call waclose(tunita2,ierr)
               if(ierr.ne.0)call errquit('problem closing t_a2 file',
     $              ierr)
            endif
            if(.not.ma_chop_stack(l_amos))call errquit(
     $           'could not chop stack at amos',1)
         endif
      enddo
      end
         





