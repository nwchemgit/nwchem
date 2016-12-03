      subroutine mp2_wijab_uhf(rtdb,nir,nva_lo,nva_hi,
     $     sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,sym_lo_vb,
     $     sym_hi_vb,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,num_oa,
     $     num_ob,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $     nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $     nob_lo,nob_hi,nbf,t_a,t_b,t_eps_a,t_eps_b,
     $     tunita,tunitb,
     $     eval_a,eval_b,
     $     g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "bitops_decls.fh"
#include "rtdb.fh"
      integer rtdb,i,a_a,nbf
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_oa(0:7),sym_hi_oa(0:7)
      integer sym_lo_va(0:7),sym_hi_va(0:7)
      integer sym_lo_ob(0:7),sym_hi_ob(0:7)
      integer sym_lo_vb(0:7),sym_hi_vb(0:7)
      integer noa,nva,nob,nvb
      integer num_oa(0:7),num_va(0:7)
      integer num_ob(0:7),num_vb(0:7)
      integer tunita,tunitb
      integer oseg_lo,oseg_hi
      integer nir
      integer g_w_ij_a,g_w_ij_b,g_w_ab_a,g_w_ab_b
      integer tcounta,tcountb
      double precision eval_a(*),eval_b(*)
      integer irs_a(*),irs_b(*)
      integer nva_lo_local,nva_hi_local
      double precision t_a(*),t_b(*),t_eps_a(*),t_eps_b(*)
      double precision zero,one,minusone,half,mhalf
      integer symi,syma,symia,symb,symj
      integer b,counta,countb
      double precision tunitptra,tunitptrb
      integer j
      double precision epsilon_a,epsilon_b
      integer msg_w_ij_a_sum,msg_w_ij_b_sum
      integer msg_w_ab_a_sum,msg_w_ab_b_sum
      parameter(msg_w_ij_a_sum=9996)
      parameter(msg_w_ij_b_sum=9995)
      parameter(msg_w_ab_a_sum=9994)
      parameter(msg_w_ab_b_sum=9993)
      integer wijacnt,wijbcnt,wabacnt,wabbcnt
c
      logical scs                    ! flag for SCS-MP2
      double precision fss,fos,hfss  ! scaling factors for SCS-MP2
c     
      integer k_w_ij_a,k_w_ij_b
      integer l_w_ij_a,l_w_ij_b
      integer k_w_ab_a,k_w_ab_b
      integer l_w_ab_a,l_w_ab_b
c     
c     inline functions to compute offset to diagonal symmetry
c     blocks in square arrays of dimension occupied and virtual
c     
      integer off_oa, off_va, off_ob, off_vb
      off_oa(symj) = sym_lo_oa(symj)-noa_lo+(sym_lo_oa(symj)-noa_lo)*noa
      off_ob(symj) = sym_lo_ob(symj)-nob_lo+(sym_lo_ob(symj)-noa_lo)*nob
      off_va(symb) = sym_lo_va(symb)-nva_lo+(sym_lo_va(symb)-nva_lo)*nva
      off_vb(symb) = sym_lo_vb(symb)-nvb_lo+(sym_lo_vb(symb)-nvb_lo)*nvb
c
#include "bitops_funcs.fh"
c     
      noa=noa_hi-noa_lo+1
      nva=nva_hi-nva_lo+1
      nob=nob_hi-nob_lo+1
      nvb=nvb_hi-nvb_lo+1
      if(.not.ma_push_get(mt_dbl,noa*noa,'w_ij_a',l_w_ij_a,k_w_ij_a))
     $     call errquit('could not alloc w_ij_a',noa*noa, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nva*nva,'w_ab_a',l_w_ab_a,k_w_ab_a))
     $     call errquit('could not alloc w_ab_a',nva*nva, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nob*nob,'w_ij_b',l_w_ij_b,k_w_ij_b))
     $     call errquit('could not alloc w_ij_b',nob*nob, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nvb*nvb,'w_ab_b',l_w_ab_b,k_w_ab_b))
     $     call errquit('could not alloc w_ab_b',nvb*nvb, MA_ERR)
c     
      zero=0.0d0
      one=1.0d0
      minusone=-1.0d0
      half=0.5d0
      mhalf=-0.5d0
      tunitptra=1
      tunitptrb=1
      call dfill((noa*noa),zero,dbl_mb(k_w_ij_a),1)
      call dfill((nva*nva),zero,dbl_mb(k_w_ab_a),1)
      call dfill((nob*nob),zero,dbl_mb(k_w_ij_b),1)
      call dfill((nvb*nvb),zero,dbl_mb(k_w_ab_b),1)

c     flag and scaling factors for SCS-MP2

      if (.not. rtdb_get(rtdb, 'mp2:scs', mt_log, 1, scs))
     $  scs = .false.       ! flag for SCS gradient
      if(scs)then
        if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $    fss = 1.0d0/3.0d0   ! same spin scaling factor for SCS-MP2
        if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $    fos = 1.2d0         ! opposite spin scaling factor for SCS-MP2
      else
        if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $    fss = one
        if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $    fos = one
        if ((fss.ne.one).or.(fos.ne.one)) then
          if (ga_nodeid().eq.0) then
            write(6,*) "You set fss or fos but not scs; setting scs=T"
          endif
          scs = .true.
        endif
      endif
      hfss=half*fss
c
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         do a_a=nva_lo_local,nva_hi_local
            syma=irs_a(a_a)
            symia=ieor(syma,symi)
c     
            call mp2_read_tijab(nva_lo, nva_hi, irs_a, symia,
     $           num_oa, sym_hi_oa, sym_lo_oa, tunita, tunitptra, t_a)
c     
            tcounta=1
            do b=nva_lo,nva_hi
               symb=irs_a(b)
               symj=ieor(symb,symia)
               if(num_oa(symj).ne.0)then
                  counta=sym_hi_oa(symj)-sym_lo_oa(symj)+1
                  do j=sym_lo_oa(symj),sym_hi_oa(symj)
                     epsilon_a=eval_a(a_a)+eval_a(b)
     $                    -eval_a(i)-eval_a(j)
                     t_eps_a(tcounta+j-sym_lo_oa(symj))
     $                    =t_a(tcounta+j-sym_lo_oa(symj))*epsilon_a
                  enddo
                  tcounta=tcounta+counta
               endif
            enddo
c     
            call mp2_read_tijab(nvb_lo, nvb_hi, irs_b, symia,
     $           num_ob, sym_hi_ob, sym_lo_ob, tunitb, tunitptrb, t_b)
c     
            tcountb=1
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               symj=ieor(symb,symia)
               if(num_ob(symj).ne.0)then
                  countb=sym_hi_ob(symj)-sym_lo_ob(symj)+1
                  do j=sym_lo_ob(symj),sym_hi_ob(symj)
                     epsilon_b=eval_a(a_a)+eval_b(b)
     $                    -eval_a(i)-eval_b(j)
                     t_eps_b(tcountb+j-sym_lo_ob(symj))
     $                    =t_b(tcountb+j-sym_lo_ob(symj))*epsilon_b
                  enddo
                  tcountb=tcountb+countb
               endif
            enddo
            tcounta=1
            tcountb=1
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_oa(symj).gt.0.and.num_va(symb).gt.0)then
                  wijacnt = off_oa(symj)+k_w_ij_a
                  call dgemm('n','t',num_oa(symj),
     $                 num_oa(symj),num_va(symb),hfss,
     $                 t_a(tcounta),num_oa(symj),t_eps_a(tcounta),
     $                 num_oa(symj),one,dbl_mb(wijacnt),noa)
                  wabacnt = off_va(symb)+k_w_ab_a
                  call dgemm('t','n',num_va(symb),
     $                 num_va(symb),num_oa(symj),hfss,
     $                 t_a(tcounta),num_oa(symj),t_eps_a(tcounta),
     $                 num_oa(symj),one,dbl_mb(wabacnt),nva)
                  tcounta=tcounta+num_oa(symj)*num_va(symb)
                  wabacnt=wabacnt+num_va(syma)*num_va(symb)
               endif
               if(num_ob(symj).gt.0.and.num_vb(symb).gt.0)then
                  wijbcnt = off_ob(symj)+k_w_ij_b
                  call dgemm('n','t',num_ob(symj),
     $                 num_ob(symj),num_vb(symb),fos,
     $                 t_b(tcountb),num_ob(symj),t_eps_b(tcountb),
     $                 num_ob(symj),one,dbl_mb(wijbcnt),nob)
                  wabbcnt = off_vb(symb)+k_w_ab_b
                  call dgemm('t','n',num_vb(symb),
     $                 num_vb(symb),num_ob(symj),fos,
     $                 t_b(tcountb),num_ob(symj),t_eps_b(tcountb),
     $                 num_ob(symj),one,dbl_mb(wabbcnt),nvb)
                  tcountb=tcountb+num_ob(symj)*num_vb(symb)
               endif
            enddo
         enddo
      enddo
      call ga_dgop(msg_w_ij_a_sum,dbl_mb(k_w_ij_a),noa*noa,'+')
      call ga_dgop(msg_w_ab_a_sum,dbl_mb(k_w_ab_a),nva*nva,'+')
      call ga_dgop(msg_w_ij_b_sum,dbl_mb(k_w_ij_b),nob*nob,'+')
      call ga_dgop(msg_w_ab_b_sum,dbl_mb(k_w_ab_b),nvb*nvb,'+')
      call ga_sync()
      if(ga_nodeid().eq.0)then
         call util_inplace_transpose(noa, dbl_mb(k_w_ij_a))
         call ga_acc(g_w_ij_a,noa_lo,noa_hi,noa_lo,noa_hi,
     $        dbl_mb(k_w_ij_a),noa,one)
         call util_inplace_transpose(nob, dbl_mb(k_w_ij_b))
         call ga_acc(g_w_ij_b,nob_lo,nob_hi,nob_lo,nob_hi,
     $        dbl_mb(k_w_ij_b),nob,one)
         call util_inplace_transpose(nva, dbl_mb(k_w_ab_a))
         call ga_acc(g_w_ab_a,1,nva,1,nva,dbl_mb(k_w_ab_a),nva,one)
         call util_inplace_transpose(nvb, dbl_mb(k_w_ab_b))
         call ga_acc(g_w_ab_b,1,nvb,1,nvb,dbl_mb(k_w_ab_b),nvb,one)
      endif
c     
      if (.not. ma_chop_stack(l_w_ij_a)) call errquit
     $     ('mp2_wijab_uhf: chop stack failed',0, MA_ERR)
c     
      end

