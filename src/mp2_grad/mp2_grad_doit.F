      subroutine mp2_grad_uhf(rtdb,g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $     irs_a,sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,  
     $     num_oa,num_va,g_vecs_b,nob_lo,nob_hi,nvb_lo,nvb_hi,
     $     irs_b,sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $     num_ob,num_vb,eval_a,eval_b,basis,geom,tol2e,e2,nbf,
     $     task,pass,oblk,
     $     g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b,
     $     g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b,
     $     g_wai_a,g_wai_b,g_lai_a,g_lai_b,
     $     grad)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "cmp2ps.fh"
#include "eaf.fh"
#include "util.fh"
      integer rtdb, g_vecs_a,g_vecs_b,pass
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
      integer k_jb_ia,l_jb_ia
      integer itmp,kunit
      integer nnbf,ioff
      character*256 fnamek
      character*256 fnameta,fnametb
      integer nva_lo_local,nva_hi_local
      parameter(kunit=11)
      integer tunita,tunitb
      logical ocoul,oexch
      double precision tol2e,e2
      double precision eval_a(*),eval_b(*)
      integer oseg_lo,oseg_hi
      character*8 task
      integer k_t_a_in,k_t_b_in,l_t_a_in,l_t_b_in
      integer k_t_eps_a,k_t_eps_b,l_t_eps_a,l_t_eps_b
      integer k_lai3_a,k_lai3_b
      integer l_lai3_a,l_lai3_b
      integer k_lai1_a,k_lai1_b
      integer l_lai1_a,l_lai1_b
      integer g_lai_a,g_lai_b
      integer g_wai_a,g_wai_b
      integer k_vvvo,l_vvvo
      integer k_vooo,l_vooo
c
      integer g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b
      integer g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b
      double precision grad(3,*)
c     
c...  for texas blocking in moints_semi
c     
      logical oblk
c     
      if(task.eq.'gradient')then
         ocoul=.false.
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
      if(num_i_perseg.lt.1)call errquit(
     $     'segmentation for multipassing incorrect',1)
      npass=0
      do oseg_lo=noa_lo,noa_hi,num_i_perseg
         npass=npass+1
         if(ga_nodeid().eq.0)
     $        write(6,*) '# of passes through ints in grad_doit = ',
     $        npass
         oseg_hi=min((oseg_lo+num_i_perseg-1),noa_hi)
c     
c     Generate the part transformed integrals on disk (ia|uv)
c     
c     i in batch, a distributed, (u>=v) with zeroes removed.
c     
         call pstat_on(ps_moints)
         call int_init(rtdb,1,basis)
         call schwarz_init(geom,basis)
         call moints_semi(basis,tol2e,ocoul,oexch,oseg_lo,oseg_hi,
     $        noa_lo,noa_hi,nva_lo,nva_hi,g_vecs_a,oblk)
         call int_terminate()
         call schwarz_tidy()
         call pstat_off(ps_moints)
c     
c     make t-pure-spin in tunita, t-mixed-spin in tunitb
c
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'amos',l_amos,k_amos))
     $        call errquit('mp2_grad_doit: error with ma',1)
         if(.not.ma_push_get(mt_dbl,nbf*nbf,'bmos',l_bmos,k_bmos))
     $        call errquit('mp2_grad_doit: error with ma',1)
         if(.not.ma_push_get(mt_int,nva*nir,'toff_a',
     $        l_toff_a,k_toff_a))
     $        call errquit('could not alloc toff_a',1)
         if(.not.ma_push_get(mt_int,nvb*nir,'toff_b',
     $        l_toff_b,k_toff_b))
     $        call errquit('could not alloc toff_b',1)
         call ga_get(g_vecs_a,1,nbf,1,nbf,dbl_mb(k_amos),nbf)
         call ga_get(g_vecs_b,1,nbf,1,nbf,dbl_mb(k_bmos),nbf)
c     
         call util_file_name('khlf',.true.,.true.,fnamek)
         call waopen(kunit,fnamek,1,0,itmp)
         if(itmp.ne.0)call errquit('mp2_grad_doit waopen?',1)
         call moints_vrange(kunit,nva_lo_local,nva_hi_local,nnbf,ioff)
c     
         if(.not.ma_push_get(mt_int,(nnbf+mod(nnbf,2)),'scatter',
     $        l_v,k_v))call errquit('could not alloc k_v',1)
         call moints_getscattv(kunit,nnbf,int_mb(k_v))
c     
         if(task.eq.'gradient')then
            call util_file_name('ta',.true.,.true.,fnameta)
            call util_file_name('tb',.true.,.true.,fnametb)
            if (eaf_open(fnameta,eaf_rw,tunita) .ne. 0)
     $           call errquit('mp2_grad: eaf open of ta failed',0)
            if (eaf_open(fnametb,eaf_rw,tunitb) .ne. 0)
     $           call errquit('mp2_grad: eaf open of tb failed',0)
         endif
         call pstat_on(ps_maket)
         call mp2_make_tuhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $        nob_lo,nob_hi,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $        sym_lo_vb,sym_hi_vb,num_ob,num_vb,oseg_lo,oseg_hi,
     $        irs_a,dbl_mb(k_amos),irs_b,dbl_mb(k_bmos),
     $        eval_a,eval_b,int_mb(k_toff_a),int_mb(k_toff_b),
     $        nva,noa,nvb,nob,nir,ioff,nnbf,
     $        int_mb(k_v),nva_lo_local,nva_hi_local,
     $        e2,task,tunita,tunitb,kunit)
         call pstat_off(ps_maket)
c     
         if(task.eq.'gradient')then
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
c     
            call pstat_on(ps_pijab)
            call mp2_pijab_uhf(nva_lo,nva_hi,sym_lo_oa,
     $           sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $           sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $           nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $           nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $           nob_lo,nob_hi,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
     $           tunita,tunitb,
     $           eval_a,eval_b,g_p_ij_a,g_p_ij_b,
     $           g_p_ab_a,g_p_ab_b)
            call pstat_off(ps_pijab)
            call pstat_on(ps_wijab)
            call mp2_wijab_uhf(nir,nva_lo,nva_hi,
     $           sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,sym_lo_vb,
     $           sym_hi_vb,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $           num_oa,num_ob,num_va,num_vb,
     $           oseg_lo,oseg_hi,irs_a,irs_b,
     $           nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $           nob_lo,nob_hi,nbf,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
     $           dbl_mb(k_t_eps_a),dbl_mb(k_t_eps_b),
     $           tunita,tunitb,eval_a, eval_b,
     $           g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b)
            call pstat_off(ps_wijab)
c
c...  Form the terms in the Lagrangian which are not Fock builds:
c     
            if(.not.ma_push_get(mt_dbl,nbf*nbf,'vvvo',l_vvvo,k_vvvo))
     $           call errquit('could not alloc vvvo',1)
            if(.not.ma_push_get(mt_dbl,nbf*nbf,'vooo',l_vooo,k_vooo))
     $           call errquit('could not alloc vooo',1)
            if(.not.ma_push_get(mt_dbl,noa*nva,'lai3_a',
     $           l_lai3_a,k_lai3_a))call errquit(
     $           'could not alloc lai3_a',1)
            if(.not.ma_push_get(mt_dbl,nob*nvb,'lai3_b',
     $           l_lai3_b,k_lai3_b))call errquit(
     $           'could not alloc lai3_b',1)
            if(.not.ma_push_get(mt_dbl,noa*nva,'lai1_a',
     $           l_lai1_a,k_lai1_a))call errquit(
     $           'could not alloc lai1_a',1)
            if(.not.ma_push_get(mt_dbl,nob*nvb,'lai1_b',
     $           l_lai1_b,k_lai1_b))call errquit(
     $           'could not alloc lai1_b',1)
            if(.not.ma_push_get(mt_dbl,nbf*nbf,'jb_ia',l_jb_ia,k_jb_ia))
     $           call errquit('could not alloc jb_ia',1)
            if(.not.ma_push_get(mt_dbl,nbf*nbf,'temp',l_temp,k_temp))
     $           call errquit('could not alloc temp',1)
c
            call pstat_on(ps_lai)
            call mp2_lai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $           sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $           sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $           nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $           dbl_mb(k_amos),dbl_mb(k_bmos),int_mb(k_v),
     $           nbf,ioff,nnbf,nva_lo_local,nva_hi_local,
     $           noa_lo,noa_hi,nob_lo,nob_hi,nva,nvb,
     $           noa,nob,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
     $           dbl_mb(k_vvvo),dbl_mb(k_vooo),
     $           dbl_mb(k_jb_ia),dbl_mb(k_temp),
     $           tunita,tunitb,dbl_mb(k_lai3_a),dbl_mb(k_lai3_b),
     $           dbl_mb(k_lai1_a),dbl_mb(k_lai1_b),
     $           g_lai_a,g_lai_b,g_wai_a,g_wai_b,kunit)
            call pstat_off(ps_lai)
c
            if(.not. ma_chop_stack(l_vvvo))
     $           call errquit('mp2_grad_doit: failed chopping vvvo',0)
c     
            call waclose(kunit,ierr)
            if(ierr.ne.0)call errquit('problem closing kunit', ierr)
            call util_file_unlink(fnamek)
c     
            call mp2_nonsep_uhf( rtdb, geom,
     $           basis,
     $           nbf, 
     $           nir, 
     $           noa, nva, nva_lo, nva_hi, num_va, num_oa, 
     $           nob, nvb, nvb_lo, nvb_hi, num_vb, num_ob,
     $           sym_lo_oa, sym_hi_oa, sym_lo_va, sym_hi_va,
     $           sym_lo_ob, sym_hi_ob, sym_lo_vb, sym_hi_vb,
     $           oseg_lo, oseg_hi, 
     $           irs_a, irs_b,
     $           dbl_mb(k_amos), dbl_mb(k_bmos), 
     $           nva_lo_local, nva_hi_local, 
     $           tunita, tunitb, grad)
c     
            if (util_print('iostats', print_high) .and. 
     $           ga_nodeid().eq.0) then
               call eaf_print_stats(tunita)
               call eaf_print_stats(tunitb)
            endif
            if (eaf_close(tunita) .ne. 0) call errquit
     $           ('mp2_grad: error closing ta',0)
            if (eaf_close(tunitb) .ne. 0) call errquit
     $           ('mp2_grad: error closing tb',0)
            call util_file_unlink(fnameta)
            call util_file_unlink(fnametb)
         endif
         if(.not.ma_chop_stack(l_amos))call errquit(
     $        'could not chop stack at amos',1)
      enddo
      end
