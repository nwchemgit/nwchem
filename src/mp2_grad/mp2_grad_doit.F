      subroutine mp2_grad_uhf(rtdb,g_vecs_a,noa_lo,noa_hi,nva_lo,nva_hi,
     $     irs_a,sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,  
     $     num_oa,num_va,g_vecs_b,nob_lo,nob_hi,nvb_lo,nvb_hi,
     $     irs_b,sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $     num_ob,num_vb,

     $     sym_lo_foa,sym_hi_foa,num_foa,
     $     sym_lo_fob,sym_hi_fob,num_fob,

     $     eval_a,eval_b,basis,geom,tol2e,e2,esing,etrip,nbf,nmo,
     $     task,passab,oblk,oskel,osym,
     $     g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b,g_liy_a,g_liy_b,
     $     g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b,
     $     g_wai_a,g_wai_b,g_lai_a,g_lai_b,
     $     grad,scftype)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cmp2ps.fh"
#include "eaf.fh"
#include "util.fh"
#include "rtdb.fh"
      integer rtdb, g_vecs_a,g_vecs_b
      integer passab            ! 1 on first call (aa+ab) 2 on second (bb+ba)
      integer noa_lo,noa_hi,nob_lo,nob_hi
      integer irs_a(*),irs_b(*)
      integer nva_lo,nva_hi,nvb_lo,nvb_hi
      integer sym_lo_foa(0:7),sym_hi_foa(0:7)
      integer sym_lo_oa(0:7),sym_hi_oa(0:7)
      integer sym_lo_va(0:7),sym_hi_va(0:7)
      integer sym_lo_fob(0:7),sym_hi_fob(0:7)
      integer sym_lo_ob(0:7),sym_hi_ob(0:7)
      integer sym_lo_vb(0:7),sym_hi_vb(0:7)
      integer num_foa(0:*),num_oa(0:*),num_va(0:*)
      integer num_fob(0:*),num_ob(0:*),num_vb(0:*)
      integer basis
      integer geom
      integer nbf,nmo
      logical oskel,osym
      character*(*) scftype
      integer num_i_perseg,npass,npasses
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
      integer kunit
      integer nnbf,ioff
      character*256 fnamek
      character*256 fnameta,fnametb
      integer nva_lo_local,nva_hi_local
*     parameter(kunit=11)  ! it becomes output argument if EAF is used
      integer tunita,tunitb
      logical oprint_mem, oprint_where
      double precision tol2e,e2,esing,etrip
      double precision eval_a(*),eval_b(*)
      integer oseg_lo,oseg_hi, noa_frozen, nob_frozen
      character*8 task
      integer k_t_a_in,k_t_b_in,l_t_a_in,l_t_b_in
      integer k_t_eps_a,k_t_eps_b,l_t_eps_a,l_t_eps_b
      integer k_lai3_a,k_lai3_b,k_liy3_a,k_liy3_b
      integer l_lai3_a,l_lai3_b,l_liy3_a,l_liy3_b
      integer k_lai1_a,k_lai1_b
      integer l_lai1_a,l_lai1_b
      integer g_lai_a,g_lai_b
      integer g_wai_a,g_wai_b
      integer k_vvvo,l_vvvo
      integer k_vooo,l_vooo
c     
      integer g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b,g_liy_a
      integer g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b,g_liy_b
      double precision grad(3,*)
      double precision util_scratch_dir_avail_for_me
      external util_scratch_dir_avail_for_me
      integer max_ga,max_ma,avail_ga,avail_ma,iavail_file
      integer eaftype,eaf_size_in_mb,npass_first
      double precision max_file, avail_file
      double precision k_file_size,t_file_size,p_file_size
      integer inntsize
      logical status
      integer mp2_eaftype
      external mp2_eaftype
      character*20 info
c     
c...  for texas blocking in moints_semi
c     
      logical oblk
c     
      oprint_mem = util_print('multipass',print_default) .and. 
     $     ga_nodeid().eq.0
      oprint_where = util_print('where',print_default) .and. 
     $     ga_nodeid().eq.0
c     
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      if (scftype .eq. 'RHF') then
         info = ' for RHF alpha+beta'
      else if (passab .eq. 1) then
         info = ' for UHF alpha'
      else
         info = ' for UHF beta'
      end if
c     
      noa_frozen = noa_lo-1
      nob_frozen = nob_lo-1
      noa=noa_hi-noa_lo+1  !!!!!!! NOTE NOA/B HERE REFER TO ACTIVE OCCUPIED
      nva=nva_hi-nva_lo+1
      nob=nob_hi-nob_lo+1
      nvb=nvb_hi-nvb_lo+1
c     
      if (osym) then
         call sym_group_name(geom,group_name)
         if(.not.sym_char_table(group_name,nop,nir,iclass_dim,
     $        zir,zclass,chars))
     $        call errquit('mp2_grad: no char table',0, INPUT_ERR)
      else
         nir = 1
      endif
c     
c
      if (task .eq. 'gradient') then
        call intd_init(rtdb,1,basis) ! int_init allocates some memory
      else
        call int_init(rtdb,1,basis) ! int_init allocates some memory
      endif
      avail_ma = MA_inquire_avail(mt_dbl)
      call int_terminate()
      call util_file_name('kh',.true.,.true.,fnamek)
      if(rtdb_get(rtdb,'mp2:scratchdisk',mt_dbl,1,avail_file)) then
         avail_file=avail_file*1024d0*128d0
         write(0,*) ' #avail ',avail_file
      else
         avail_file=util_scratch_dir_avail_for_me()*128.0d0
      endif
cold         avail_file=0.99d0*((avail_ma*(9d0*3d0))/(16d0*4d0))
cold         avail_ma=(7d0*avail_ma)/16d0
#ifdef NOIO
cold      avail_ma=avail_ma/2
      avail_ma=avail_ma*7/8
      avail_file=0
#endif
c
      if (ga_uses_ma() .and. ga_memory_limited()) 
     $     avail_ma = avail_ma - ga_memory_avail()
      if(ga_memory_limited())then
         avail_ga=ga_memory_avail()
      else
         avail_ga=0
      end if
      call ga_igop(11,avail_ma,1,'min')
      call ga_igop(12,avail_ga,1,'min')
      avail_ga = avail_ga/8  ! convert to words
      call ga_dgop(13,avail_file,1,'min')
c     
      if (oprint_mem) then
         write(6,201) ga_uses_ma(), ga_memory_limited()
 201     format(/' GA uses MA = ',l1,'    GA memory limited = ',l1)
         write(6,21) dble(avail_ma), dble(avail_ga), avail_file
 21      format(/' Available: '/'  local mem=', 1p, 1d10.2,
     ,        /' global mem=', 1p, 1d10.2,
     $        /' local disk=', 1p, 1d10.2)
         call util_flush(6)
      end if
c     
      if(.not.rtdb_get(rtdb,'mp2:npasses',mt_int,1,npass_first)) 
     C   npass_first=1
      do npasses = npass_first, noa_hi-noa_lo+1
         num_i_perseg = (noa_hi-noa_lo)/npasses + 1
         call mp2_memory(basis,geom,rtdb,tol2e,oskel,nbf,noa,nva,nob,
     $        nvb,
     $        num_i_perseg,max_ma,max_ga,max_file,task,scftype,
     k        k_file_size,t_file_size,p_file_size)
         if(oprint_mem) then
            write(6,22) npasses, num_i_perseg, max_ma, max_ga, max_file
 22         format(1x,i3,' passes of',i4,':  ',
     $           ' ',i12,'            ',i15,
     $           '            ',f15.0)
            call util_flush(6)
         end if
#ifdef NOIO
cold         avail_file=(avail_ma-max_ma)/2
                  avail_file=(avail_ma-max_ma)
      if(avail_file.lt.0) call errquit(
     C     'grad_doit: increase MA memory by at least',int(-avail_file),
     M     MA_err)
         if(oprint_mem) then
            write(6,222) npasses, num_i_perseg, avail_ma, avail_ga, 
     C           dble(avail_file)
 222         format(1x,i3,' passes of',i4,': M',
     $           ' ',i12,'           G',i15,
     $           '           F',f15.0)
            call util_flush(6)
         end if
#endif
         if(ga_uses_ma().and..not.ga_memory_limited())then
            if(max_ma+max_ga.le.avail_ma.and.max_file.le.
     $           avail_file) goto 999
         elseif(ga_memory_limited())then
            if(max_ma.le.avail_ma.and.max_file.le.avail_file.
     $           and.max_ga.le.avail_ga) goto 999
         elseif(max_ma.le.avail_ma.and.max_file.le.avail_file)then
            if(ga_nodeid().eq.0) goto 999
         end if
      end do
      call errquit('mp2_grad:insufficient resources',1, INPUT_ERR)
 999  continue
      if (oprint_mem) then
         write(6,*)
         call util_flush(6)
      end if
c     
      if (npasses .gt. 200) call errquit
     $     ('mp2_grad: too many passes',npasses, CALC_ERR)
c     
      npass=0
      do oseg_lo=noa_lo,noa_hi,num_i_perseg
         npass=npass+1
         if(ga_nodeid().eq.0 .and. 
     $        util_print('passes',print_default)) then
            write(6,1101) npass, npasses,info,util_wallsec()
 1101       format(' Semi-direct pass number', i4, ' of',i4,1x,a,
     $           ' at ',f10.1,'s')
            call util_flush(6)
         end if
         oseg_hi=min((oseg_lo+num_i_perseg-1),noa_hi)
c     
c     Generate the part transformed integrals on disk (ia|uv)
c     
c     i in batch, a distributed, (u>=v) with zeroes removed.
c     
         call pstat_on(ps_moints)
         call int_init(rtdb,1,basis)
         call schwarz_init(geom,basis)
         call moints_semi(rtdb,basis,tol2e,oseg_lo,oseg_hi,
     $        noa_lo,noa_hi,nva_lo,nva_hi,g_vecs_a,oblk,
     K        k_file_size)
         call int_terminate()
         call schwarz_tidy()
         call pstat_off(ps_moints)
c     
         if (oprint_where) then
            write(6,111) 'moints_semi', util_wallsec()
 111        format(' Done ',a,' at ', f10.1,'s')
            call util_flush(6)
         end if
c     
c     make t-pure-spin in tunita, t-mixed-spin in tunitb
c
c     Note that for RHF we share memory between alpha and beta mos
c     
         if(.not.ma_push_get(mt_dbl,nbf*nmo,'amos',l_amos,k_amos))
     $        call errquit('mp2_grad_doit: error with ma',1, MA_ERR)
         if (scftype .eq. 'UHF') then
            if(.not.ma_push_get(mt_dbl,nbf*nmo,'bmos',l_bmos,k_bmos))
     $           call errquit('mp2_grad_doit: error with ma',1, MA_ERR)
         else
            k_bmos = k_amos
         endif
c
         if(.not.ma_push_get(mt_int,nva*nir,'toff_a',
     $        l_toff_a,k_toff_a))
     $        call errquit('mp2_grad_doit: could not alloc toff_a',1,
     &       MA_ERR)
         if(.not.ma_push_get(mt_int,nvb*nir,'toff_b',
     $        l_toff_b,k_toff_b))
     $        call errquit('mp2_grad_doit: could not alloc toff_b',1,
     &       MA_ERR)
c
#if 0
         if (ga_nodeid() .eq. 0) then
            call ga_get(g_vecs_a,1,nbf,1,nmo,dbl_mb(k_amos),nbf)
            if (scftype .eq. 'UHF') 
     $           call ga_get(g_vecs_b,1,nbf,1,nmo,dbl_mb(k_bmos),nbf)
         endif
         call ga_brdcst(134, dbl_mb(k_amos), nbf*nmo*8, 0)
         if (scftype .eq. 'UHF') 
     $        call ga_brdcst(135, dbl_mb(k_bmos), nbf*nmo*8, 0)
#else
         call util_mygabcast2(g_vecs_a,1,nbf,1,nmo,dbl_mb(k_amos),nbf)
         if (scftype .eq. 'UHF') call util_mygabcast2(
     G        g_vecs_b,1,nbf,1,nmo,dbl_mb(k_bmos),nbf)
#endif
c     
#ifdef NOIO
         eaftype=mp2_eaftype(rtdb,k_file_size)
         if(ga_nodeid().eq.0) write(6,*) ' mp2_eaf for kh ',eaftype
#else

         eaftype=EAF_RW
#endif
         call util_file_name('kh',.true.,.true.,fnamek)
         if (eaf_open(fnamek, eaftype, kunit).ne.0) 
     $        call errquit('tunitb: eaf_open failed', 0,
     &        DISK_ERR)
         call moints_vrange(kunit,nva_lo_local,nva_hi_local,nnbf,ioff)
c     
         if(.not.ma_push_get(mt_int,(nnbf+mod(nnbf,2)),'scatter',
     $        l_v,k_v))call errquit('mp2_grad: could not alloc k_v',1,
     &       MA_ERR)
         call moints_getscattv(kunit,nnbf,int_mb(k_v))
c     
         if(task.eq.'gradient')then
#ifdef NOIO
         eaftype=mp2_eaftype(rtdb,t_file_size)
         if(ga_nodeid().eq.0) write(6,*) ' mp2_eaf for ta tb ',eaftype
#else
         eaftype=EAF_RW
#endif
            call util_file_name('ta',.true.,.true.,fnameta)
            call util_file_name('tb',.true.,.true.,fnametb)
         if(ga_nodeid().eq.0) write(6,*) ' mp2_eaf for ta  ',eaftype
            if (eaf_open(fnameta,eaftype,tunita) .ne. 0)
     $           call errquit('mp2_grad: eaf open of ta failed',0,
     &       DISK_ERR)
         if(ga_nodeid().eq.0) write(6,*) ' mp2_eaf for tb  ',eaftype
            if (eaf_open(fnametb,eaftype,tunitb) .ne. 0)
     $           call errquit('mp2_grad: eaf open of tb failed',0,
     &       DISK_ERR)

         end if
         call pstat_on(ps_maket)
         call mp2_make_tuhf(rtdb,nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $        sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $        nob_lo,nob_hi,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,
     $        sym_lo_vb,sym_hi_vb,num_ob,num_vb,oseg_lo,oseg_hi,
     $        irs_a,dbl_mb(k_amos),irs_b,dbl_mb(k_bmos),
     $        eval_a,eval_b,int_mb(k_toff_a),int_mb(k_toff_b),
     $        nva,noa,nvb,nob,nir,ioff,nnbf,
     $        int_mb(k_v),nva_lo_local,nva_hi_local,
     $        e2,esing,etrip,task,tunita,tunitb,kunit)
         call pstat_off(ps_maket)
c     
         if (oprint_where) then
            write(6,111) 'maket',util_wallsec()
            call util_flush(6)
         end if
c     
         if(task.eq.'gradient')then
            if(.not.ma_push_get(mt_dbl,noa*nva,
     $           't_a_in',l_t_a_in,k_t_a_in))call errquit
     $           ('mp2_grad: could not alloc t_a',1, MA_ERR)
            if(.not.ma_push_get(mt_dbl,nob*nvb,
     $           't_b_in',l_t_b_in,k_t_b_in))call errquit
     $           ('mp2_grad: could not alloc t_b',1, MA_ERR)
            if(.not.ma_push_get(mt_dbl,noa*nva,
     $           't_eps_a',l_t_eps_a,k_t_eps_a))call errquit
     $           ('mp2_grad: could not alloc t_eps_a',1, MA_ERR)
            if(.not.ma_push_get(mt_dbl,nob*nvb,
     $           't_eps_b',l_t_eps_b,k_t_eps_b))call errquit
     $           ('mp2_grad: could not alloc t_eps_b',1, MA_ERR)
c     
            call pstat_on(ps_pijab)
            call mp2_pijab_uhf(rtdb,nva_lo,nva_hi,sym_lo_oa,
     $           sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $           sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $           nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $           nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $           nob_lo,nob_hi,dbl_mb(k_t_a_in),dbl_mb(k_t_b_in),
     $           tunita,tunitb,
     $           eval_a,eval_b,g_p_ij_a,g_p_ij_b,
     $           g_p_ab_a,g_p_ab_b)
            call pstat_off(ps_pijab)
c     
            if (oprint_where) then
               write(6,111) 'pijab', util_wallsec()
               call util_flush(6)
            end if
c     
            call pstat_on(ps_wijab)
            call mp2_wijab_uhf(rtdb,nir,nva_lo,nva_hi,
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
            if (oprint_where) then
               write(6,111) 'wijab', util_wallsec()
               call util_flush(6)
            end if
c     
            if (.not. ma_chop_stack(l_t_a_in)) call errquit
     $            ('mp2_grad: failed chopping at t_a_in', 0, MA_ERR)
*            if(.not.ma_pop_stack(l_t_eps_b))
*     $           call errquit('could not pop stack',1)
*            if(.not.ma_pop_stack(l_t_eps_a))
*     $           call errquit('could not pop stack',1)
c     
c...  Form the terms in the Lagrangian which are not Fock builds:
c
c     (logic assumes noa_frozen=nob_frozen, noa>=nob, nva<=nvb)
c
            if(.not.ma_push_get(mt_dbl, 
     $           max(max(nva,noa_frozen)*nva,max(nvb,nob_frozen)*nvb),
     $           'vvvo',l_vvvo,k_vvvo))
     $           call errquit('mp2_grad: could not alloc vvvo',1,
     &       MA_ERR)
            if(.not.ma_push_get(mt_dbl, 
     $           max(max(noa,noa_frozen)*noa,max(nob,nob_frozen)*nob),
     $           'vooo',l_vooo,k_vooo))
     $           call errquit('mp2_grad: could not alloc vooo',1,
     &       MA_ERR)
c
            if(.not.ma_push_get(mt_dbl,noa*nva,'lai3_a',
     $           l_lai3_a,k_lai3_a))call errquit
     $           ('mp2_grad: could not alloc lai3_a',1, MA_ERR)
            if(.not.ma_push_get(mt_dbl,nob*nvb,'lai3_b',
     $           l_lai3_b,k_lai3_b))call errquit
     $           ('mp2_grad: could not alloc lai3_b',1, MA_ERR)
c
            if(.not.ma_push_get(mt_dbl,max(noa*nva,nob*nvb),
     $           't_a_in',l_t_a_in,k_t_a_in))call errquit
     $           ('mp2_grad: could not alloc t_a',1, MA_ERR)
c
            if (noa_frozen .gt. 0) then
               if(.not.ma_push_get(mt_dbl,noa_frozen*noa,'lyi3_a',
     $              l_liy3_a,k_liy3_a))call errquit
     $              ('mp2_grad: could not alloc lyi3_a',1, MA_ERR)
               if(.not.ma_push_get(mt_dbl,nob_frozen*nob,'lyi3_b',
     $              l_liy3_b,k_liy3_b))call errquit
     $              ('mp2_grad: could not alloc lyi3_b',1, MA_ERR)
            else
               k_liy3_a = 1     ! To ensure valid pointer
               k_liy3_b = 1
            endif
c
            if(.not.ma_push_get(mt_dbl,(noa+noa_frozen)*nva,'lai1_a',
     $           l_lai1_a,k_lai1_a))call errquit
     $           ('mp2_grad: could not alloc lai1_a',1, MA_ERR)
            if(.not.ma_push_get(mt_dbl,(nob+nob_frozen)*nvb,'lai1_b',
     $           l_lai1_b,k_lai1_b)) call errquit(
     $           'could not alloc lai1_b',1,0)
c
            if(.not.ma_push_get(mt_dbl,nbf*nbf,'jb_ia',l_jb_ia,k_jb_ia))
     $           call errquit('mp2_grad: could not alloc jb_ia',1, 
     $             MA_ERR) ! ia_uv
            if(.not.ma_push_get(mt_dbl,nbf*nbf,'temp',l_temp,k_temp))
     $           call errquit('mp2_grad: could not alloc temp',1,
     &       MA_ERR)
c     
            call pstat_on(ps_lai)
            call mp2_lai_uhf(rtdb,nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $           sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $           sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $           sym_lo_foa,sym_hi_foa,num_foa,
     $           sym_lo_fob,sym_hi_fob,num_fob,
     $           nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $           dbl_mb(k_amos),dbl_mb(k_bmos),int_mb(k_v),
     $           nbf,ioff,nnbf,nva_lo_local,nva_hi_local,
     $           noa_lo,noa_hi,nob_lo,nob_hi,nva,nvb,
     $           noa,nob,dbl_mb(k_t_a_in),
     $           dbl_mb(k_vvvo),dbl_mb(k_vooo),
     $           dbl_mb(k_jb_ia),dbl_mb(k_temp),
     $           tunita,tunitb,dbl_mb(k_lai3_a),dbl_mb(k_lai3_b),
     $           dbl_mb(k_lai1_a),dbl_mb(k_lai1_b),
     $           dbl_mb(k_liy3_a), dbl_mb(k_liy3_b),
     $           g_lai_a,g_lai_b,g_wai_a,g_wai_b,
     $           g_liy_a, g_liy_b,
     $           kunit)
            call pstat_off(ps_lai)
c     
            if (oprint_where) then
               write(6,111) 'lai', util_wallsec()
               call util_flush(6)
            end if
c     
            if(.not. ma_chop_stack(l_toff_a)) call errquit(
     $           'mp2_grad_doit: failed chopping l_toff_a',0, MA_ERR)
c     
            if (eaf_close(kunit).ne.0)
     $        call errquit('mp2_grad_doit: close half int file',0,
     &       DISK_ERR)
            call util_file_unlink(fnamek)
c     
            call mp2_nonsep_uhf( rtdb, geom,
     $           basis, oskel,
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
     $           tunita, tunitb, grad,
     P           p_file_size)
c     
            if (oprint_where) then
               write(6,111) 'nonsep', util_wallsec()
               call util_flush(6)
            end if
c     
            if (util_print('iostats', print_high) .and. 
     $           ga_nodeid().eq.0) then
               call eaf_print_stats(tunita)
               call eaf_print_stats(tunitb)
            end if

            if (eaf_close(tunita) .ne. 0) call errquit
     $           ('mp2_grad: error closing ta',0, DISK_ERR)
            if (eaf_close(tunitb) .ne. 0) call errquit
     $           ('mp2_grad: error closing tb',0, DISK_ERR)
            call util_file_unlink(fnameta)
            call util_file_unlink(fnametb)
         end if
         if(.not.ma_chop_stack(l_amos))call errquit(
     $        'mp2_grad: could not chop stack at amos',1, MA_ERR)
      end do
c     
      if (task .ne. 'gradient') then
        if (eaf_close(kunit).ne.0)
     $    call errquit('mp2_grad: problem closing kunit',0, DISK_ERR)
         call util_file_unlink(fnamek)
      end if
c     
c     
      if (oprint_where) then
         write(6,111) 'multipass loop', util_wallsec()
         call util_flush(6)
      end if
c     
      end
      integer function mp2_eaftype(rtdb,
     D     file_size)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer rtdb
      double precision file_size
c
      integer eaf_size_in_mb,avail_ma_in_mb
      integer inntsize
c     
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      if (.not. rtdb_get(rtdb, 'mp2:eaf_size_in_mb',
     $     MT_INT, 1, mp2_eaftype)) mp2_eaftype=
     E     inntsize*file_size/1024/1024
c     make it at least 10
      mp2_eaftype=max(10,mp2_eaftype)
c
c     check if mem  is enough
c
      avail_ma_in_mb = MA_inquire_avail(mt_byte)/1024/1024
      if(mp2_eaftype.gt.avail_ma_in_mb) then
         write(6,*) ga_nodeid(),' mp2_eaftype ',mp2_eaftype,
     M ' avail_ma_in_mb ',avail_ma_in_mb
         call errquit(
     M        'mp2_eaftype gt avail_ma by',avail_ma_in_mb-mp2_eaftype,
     M        MA_ERR)
      endif
      if(ga_nodeid().eq.0) write(6,*) ' mp2_eaftype ',mp2_eaftype,
     M ' avail_ma_in_mb ',avail_ma_in_mb

      if(ga_nodeid().eq.0) then
         write(luout,'(A,I12,F12.0)')
     I        ' EAF : mp2_eaftype filesize',
     A        mp2_eaftype,
     C        file_size
         call util_flush(luout)
      endif
      return
      end
