      subroutine mp2_make_tuhf(rtdb,nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $     sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $     nob_lo,nob_hi,nvb_lo,nvb_hi, sym_lo_ob,sym_hi_ob,
     $     sym_lo_vb,sym_hi_vb,num_ob,num_vb,   
     $     oseg_lo,oseg_hi,irs_a,c_a,irs_b,c_b,
     $     eval_a,eval_b,toff_a,toff_b,nva,noa,nvb,nob,nir,
     $     ioff,nnbf,indx_v,nva_lo_local,nva_hi_local,
     $     e2,esing,etrip,task,tunita,tunitb,kunit)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "eaf.fh"
#include "util.fh"
      integer rtdb
      character*8 task
      integer nbf,i
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer sym_lo_oa(0:*),sym_hi_oa(0:*)
      integer sym_lo_va(0:*),sym_hi_va(0:*)
      integer num_oa(0:*),num_va(0:*)
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_ob(0:*),sym_hi_ob(0:*)
      integer sym_lo_vb(0:*),sym_hi_vb(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer noa,nva,nob,nvb,counter,g_t_a,g_t_b
      integer a_a,tunita,tunitb,counta,countb
      integer oseg_lo,oseg_hi
      integer nir
      integer toff_a(nva_lo:nva_hi,0:nir-1),
     $     toff_b(nvb_lo:nvb_hi,0:nir-1)
      double precision eaa,eab
c
      double precision tunitptra,tunitptrb
c
      double precision c_a(nbf,nbf),c_b(nbf,nbf)
      integer i_batch
      double precision e2,esing,etrip,denom,t,minusone
      integer syma,symi,symb,symj,symia
      integer irs_a(*),irs_b(*),b,j
      integer indx_v(*)
      integer nnbf,ioff,kunit
      integer nva_lo_local,nva_hi_local,jb
      double precision eval_a(*),eval_b(*)
      double precision eps
      parameter (eps=1d-14)
c
      logical otdebug
      logical do_mirra,do_mirrb,do_mirr
      integer g_ta_mirr,g_tb_mirr,g_ta_org,g_tb_org
      logical util_mirrmat,tuhfreplicated
      external util_mirrmat
      integer util_mirrmem,mem_mirr
      external util_mirrmem
c
      integer l_ia_uv, k_ia_uv, l_tmp, k_tmp, l_ia_jb, k_ia_jb
      double precision ddot
      external ddot
c
#include "bitops.fh"
c
      otdebug = util_print('tdebug',print_debug)
c     
      eaa=0.0d0
      eab=0.0d0
      minusone=-1.0d0
      noa=noa_hi-noa_lo+1
      nva=nva_hi-nva_lo+1
      nob=nob_hi-nob_lo+1
      nvb=nvb_hi-nvb_lo+1
c
      if (.not. ma_push_get(mt_dbl, nbf*nbf, 'tmp', l_tmp, k_tmp))
     $     call errquit('mp2_make_t: no room for tmp', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf, 'ia_uv', l_ia_uv, k_ia_uv))
     $     call errquit('mp2_make_t: no room for ia_uv', nbf*nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, max(noa*nva,nob*nvb), 'ia_jb', 
     $     l_ia_jb, k_ia_jb)) call errquit
     $     ('mp2_make_t: no room for ia_jb', max(noa*nva,nob*nvb),
     &       MA_ERR)
c
*ga:1:0
      if(.not.ga_create(mt_dbl,noa*nva,nva,'T_a',(noa*nva),0,g_t_a))
     $     call errquit('could not alloc ga for T_a',1, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nob*nvb,nva,'T_b',(nob*nvb),0,g_t_b))
     $     call errquit('could not alloc ga for T_b',1, GA_ERR)
      if (.not.rtdb_get(rtdb, 'mp2:tuhfreplicated', mt_log, 1, 
     T     tuhfreplicated)) tuhfreplicated=.false.
      do_mirr=.false.
      do_mirra=.false.
      do_mirrb=.false.
      if(tuhfreplicated) then
         mem_mirr=
     M        util_mirrmem(1,g_t_a)+
     M        util_mirrmem(1,g_t_b)
         if(ga_nodeid().eq.0) then
            write(6,14) ' needed by mirroring: ', mem_mirr
            write(6,14) ' ga_memory_avail()  : ', ga_memory_avail()
         endif
 14      format(a,i10,' bytes')
         if(ga_memory_avail().gt. mem_mirr) then
            do_mirra=util_mirrmat(1,g_t_a,g_ta_mirr,
     D           .false., .true.)
            do_mirrb=util_mirrmat(1,g_t_b,g_tb_mirr,
     D           .false., .true.)
            else
               if(ga_nodeid().eq.0)
     W              write(6,*) ' insuff. memory  for mirroring '
            endif
      endif
      if(do_mirra.and.do_mirrb) then
         if(ga_nodeid().eq.0)
     w        write(6,*) 'mp2_make_tuhf mirroring'
         do_mirr=.true.
         g_ta_org=g_t_a
         g_tb_org=g_t_b
         g_t_a=g_ta_mirr
         g_t_b=g_tb_mirr
      else
         if(do_mirra) then
            if(.not.ga_destroy(g_ta_mirr))call errquit(
     $           'could not destroy mirrt_a handle',1, GA_ERR)
            write(6,10) ga_nodeid(),'destroyed g_ta_mirr',g_ta_mirr
 10         format(i5,a,i5)
         endif
         g_ta_org=0
         g_tb_org=0
      endif
      tunitptra=1
      tunitptrb=1
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         call ga_zero(g_t_a)
         call ga_zero(g_t_b)
         do syma=0,nir-1
            symia=ieor(syma,symi)
            counter=0
            do symb=0,nir-1
               do b=sym_lo_va(symb),sym_hi_va(symb)
                  symj=ieor(symia,symb)
                  toff_a(b,syma)=counter
                  counter=counter+num_oa(symj)
               enddo
            enddo
         enddo
         do syma=0,nir-1
            symia=ieor(syma,symi)
            counter=0
            do symb=0,nir-1
               do b=sym_lo_vb(symb),sym_hi_vb(symb)
                  symj=ieor(symia,symb)
                  toff_b(b,syma)=counter
                  counter=counter+num_ob(symj)
               enddo
            enddo
         enddo
         i_batch=i-oseg_lo+1
         do a_a=nva_lo_local,nva_hi_local 
            syma=irs_a(a_a)
            symia=ieor(syma,symi)
            call dfill((nbf*nbf),0.0d0,dbl_mb(k_tmp),1)
            call moints_rdhfint(kunit,nva_lo_local,nva_hi_local,
     $           i_batch,a_a,nnbf,ioff,dbl_mb(k_ia_uv))
            call scatter(nnbf,dbl_mb(k_tmp),indx_v,dbl_mb(k_ia_uv))
            call dfill((nbf*nbf),0.0d0,dbl_mb(k_ia_uv),1)
            call upper2square(nbf,dbl_mb(k_tmp),dbl_mb(k_ia_uv))
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_oa(symj).gt.0.and.num_va(symb).gt.0)then
                  call dgemm('n','n',nbf,num_oa(symj),nbf,1.0d0,
     $                 dbl_mb(k_ia_uv),nbf,c_a(1,sym_lo_oa(symj)),nbf,
     $                 0.0d0,dbl_mb(k_tmp),nbf)
                  call dgemm('t','n',num_oa(symj),num_va(symb),
     $                 nbf,1.0d0,dbl_mb(k_tmp),nbf,
     $                 c_a(1,sym_lo_va(symb)),
     $                 nbf,0.0d0,dbl_mb(k_ia_jb),num_oa(symj))
#if 1
                  call ga_acc(g_t_a,toff_a(sym_lo_va(symb),syma)+1,
     $                    toff_a(sym_hi_va(symb),syma)+num_oa(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    dbl_mb(k_ia_jb),1,1.0d0)
                  jb=1
                  do b=sym_lo_va(symb),sym_hi_va(symb)
                     call ga_acc(g_t_a,toff_a(a_a,symb)+1,
     $                    toff_a(a_a,symb)+num_oa(symj),
     $                    b-nva_lo+1,b-nva_lo+1,
     $                    dbl_mb(k_ia_jb-1+jb),1,minusone)
                     jb=jb+num_oa(symj)
                  enddo
#else
                  jb=1
                  do b=sym_lo_va(symb),sym_hi_va(symb)
                     if(a_a.ne.b) then
                     call ga_acc(g_t_a,toff_a(b,syma)+1,
     $                    toff_a(b,syma)+num_oa(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    dbl_mb(k_ia_jb-1+jb),1,1.0d0)
                     call ga_acc(g_t_a,toff_a(a_a,symb)+1,
     $                    toff_a(a_a,symb)+num_oa(symj),
     $                    b-nva_lo+1,b-nva_lo+1,
     $                    dbl_mb(k_ia_jb-1+jb),1,minusone)
                     endif
                     jb=jb+num_oa(symj)
                  enddo
#endif
               endif
               if(num_ob(symj).gt.0.and.num_vb(symb).gt.0)then
                  call dgemm('n','n',nbf,num_ob(symj),nbf,1.0d0,
     $                 dbl_mb(k_ia_uv),nbf,c_b(1,sym_lo_ob(symj)),nbf,
     $                 0.0d0,dbl_mb(k_tmp),nbf)
                  call dgemm('t','n',num_ob(symj),num_vb(symb),
     $                 nbf,1.0d0,dbl_mb(k_tmp),nbf,
     $                 c_b(1,sym_lo_vb(symb)),
     $                 nbf,0.0d0,dbl_mb(k_ia_jb),num_ob(symj))
#if 1
                     call ga_put(g_t_b,toff_b(sym_lo_vb(symb),syma)+1,
     $                    toff_b(sym_hi_vb(symb),syma)+num_ob(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    dbl_mb(k_ia_jb),
     L                 toff_b(sym_hi_vb(symb),syma)+num_ob(symj)-
     M                 toff_b(sym_lo_vb(symb),syma))
#else
                  jb=1
                  do b=sym_lo_vb(symb),sym_hi_vb(symb)
                     call ga_put(g_t_b,toff_b(b,syma)+1,
     $                    toff_b(b,syma)+num_ob(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    dbl_mb(k_ia_jb-1+jb),num_ob(symj))
                     jb=jb+num_ob(symj)
                  enddo
#endif
               endif
            enddo
         enddo
         call ga_sync()
         if (do_mirr) then
            call ga_mask_sync(.true.,.false.)
            call ga_merge_mirrored(g_t_a)
            call ga_mask_sync(.true.,.false.)
            call ga_merge_mirrored(g_t_b)
         endif
         do a_a=nva_lo_local,nva_hi_local
            counta = 0
            do b=nva_lo,nva_hi
               symb=irs_a(b)
               syma=irs_a(a_a)
               symj=ieor(symb,ieor(symi,syma))
               if(num_oa(symj).ne.0)then
                  call ga_get(g_t_a,toff_a(b,syma)+1,
     $                 toff_a(b,syma)+num_oa(symj),
     $                 a_a-nva_lo+1,a_a-nva_lo+1,
     $                 dbl_mb(counta+k_ia_jb),num_oa(symj))
                  do j=sym_lo_oa(symj),sym_hi_oa(symj)
                     denom=eval_a(a_a)+eval_a(b)
     $                    -eval_a(i)-eval_a(j)
                     t=dbl_mb(counta+k_ia_jb+j-sym_lo_oa(symj))/denom
                     if(abs(t).lt.eps)t=0.0d0
c
                     eaa=eaa+t*t*denom
                     dbl_mb(counta+k_ia_jb+j-sym_lo_oa(symj))=t
                  enddo
                  counta=counta+(sym_hi_oa(symj)-sym_lo_oa(symj)+1)
               endif
            enddo
            if(task.eq.'gradient' .and. counta.gt.0)then
               if (eaf_write(tunita,tunitptra,dbl_mb(k_ia_jb),8*counta) 
     $              .ne. 0) call errquit('mp2_make_t:IO on ta?',0,
     &       DISK_ERR)
               tunitptra=tunitptra+8*counta
               if (otdebug) then
                  write(6,12) 'aa', i, a_a, ddot(counta,dbl_mb(k_ia_jb),
     $                 1,dbl_mb(k_ia_jb),1)
 12               format(1x,a,2i5,f20.8)
               endif
            endif
c
            countb = 0
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               syma=irs_a(a_a)
               symj=ieor(symb,ieor(symi,syma))
               if(num_ob(symj).ne.0)then
                  call ga_get(g_t_b,toff_b(b,syma)+1,
     $                 toff_b(b,syma)+num_ob(symj),
     $                 a_a-nva_lo+1,a_a-nva_lo+1,
     $                 dbl_mb(countb+k_ia_jb),num_ob(symj))
                  do j=sym_lo_ob(symj),sym_hi_ob(symj)
                     denom=eval_a(a_a)+eval_b(b)
     $                    -eval_a(i)-eval_b(j)
                     t=dbl_mb(countb+k_ia_jb+j-sym_lo_ob(symj))/denom
                     if(abs(t).lt.eps)t=0.0d0
                     eab=eab+t*t*denom
                     dbl_mb(countb+k_ia_jb+j-sym_lo_ob(symj))=t
                  enddo
                  countb=countb+(sym_hi_ob(symj)-sym_lo_ob(symj)+1)
               endif
            enddo
            if(task.eq.'gradient' .and. countb.gt.0)then
               if (eaf_write(tunitb,tunitptrb,dbl_mb(k_ia_jb),8*countb)
     $              .ne. 0) call errquit('mp2_make_t:IO on tb?',0,
     &       DISK_ERR)
               tunitptrb=tunitptrb+8*countb
               if (otdebug) then
                  write(6,12) 'ab', i, a_a, ddot(countb,dbl_mb(k_ia_jb),
     $                 1,dbl_mb(k_ia_jb),1)
               endif
            endif
         enddo
      enddo
c
      call ga_dgop(1,eaa,1,'+')
      call ga_dgop(1,eab,1,'+')
c
      e2=e2-0.25d0*eaa-0.5d0*eab
      esing = esing + (-4.0d0*eab + eaa)/4.0d0
      etrip = etrip + (-3.0d0*eaa)/4.0d0
c
      if(do_mirr) then
         if(.not.ga_destroy(g_t_b))call errquit(
     $        'could not destroy mirr t_b handle',1, GA_ERR)
         if(.not.ga_destroy(g_t_a))call errquit(
     $        'could not destroy mirr t_a handle',1, GA_ERR)
             g_t_a=g_ta_org
             g_t_b=g_tb_org
       endif

       if(.not.ga_destroy(g_t_b))then
          if(ga_nodeid().eq.0) call ga_summarize(1)
          call ga_sync()
        call errquit(
     $         'could not destroy t_b handle',1, GA_ERR)
        endif
      if(.not.ga_destroy(g_t_a))call errquit(
     $     'could not destroy t_a handle',1, GA_ERR)
c
      if (.not. ma_chop_stack(l_tmp)) call errquit
     $     ('mp2_make_t: failed chopping stack', 0, MA_ERR)
c
      call ga_sync()
c
      end
      
      
      
      
      
      


