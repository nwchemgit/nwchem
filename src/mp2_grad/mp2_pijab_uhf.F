      subroutine mp2_pijab_uhf(rtdb,nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $     sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,sym_lo_ob,
     $     sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,nir,num_va,num_vb,
     $     oseg_lo,oseg_hi,irs_a,irs_b,nva_lo_local,nva_hi_local,
     $     noa_lo,noa_hi,nob_lo,nob_hi,t_a,t_b,
     $     tunita,tunitb,eval_a,eval_b,
     $     g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b)
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
      integer rtdb,i
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer sym_lo_oa(0:*),sym_hi_oa(0:*)
      integer sym_lo_va(0:*),sym_hi_va(0:*)
      integer num_oa(0:*),num_va(0:*)
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_ob(0:*),sym_hi_ob(0:*)
      integer sym_lo_vb(0:*),sym_hi_vb(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer noa,nva,nob,nvb
      integer a_a,tunita,tunitb
      integer oseg_lo,oseg_hi
      integer nir
      double precision zero,one,minusone,half,mhalf
      integer syma,symi,symb,symj,symia
      integer irs_a(*),irs_b(*)
      integer nva_lo_local,nva_hi_local
      double precision t_a(*),t_b(*)
      double precision tunitptra,tunitptrb
      integer g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b
      integer msg_p_ij_a_sum,msg_p_ij_b_sum
      integer msg_p_ab_a_sum,msg_p_ab_b_sum
      parameter(msg_p_ij_a_sum=10000)
      parameter(msg_p_ij_b_sum=9999)
      parameter(msg_p_ab_a_sum=9998)
      parameter(msg_p_ab_b_sum=9997)
      double precision eval_a(*),eval_b(*)
      integer tcounta,tcountb
      integer pijacnt,pijbcnt,pabacnt,pabbcnt
c     
      integer k_p_ij_a,k_p_ij_b
      integer l_p_ij_a,l_p_ij_b
      integer k_p_ab_a,k_p_ab_b
      integer l_p_ab_a,l_p_ab_b
c
      logical scs                               ! flag for SCS-MP2
      double precision fss,fos,hfss,mhfss,mfos  ! scaling factors for SCS-MP2
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
c     
c$$$  do symj = 0, nir-1
c$$$  write(6,*) symj, ' occa ', sym_lo_oa(symj), sym_hi_oa(symj)
c$$$  write(6,*) symj, ' occb ', sym_lo_ob(symj), sym_hi_ob(symj)
c$$$  write(6,*) symj, ' vira ', sym_lo_va(symj), sym_hi_va(symj)
c$$$  write(6,*) symj, ' virb ', sym_lo_vb(symj), sym_hi_vb(symj)
c$$$  write(6,*) symj, ' oa ', off_oa(symj)
c$$$  write(6,*) symj, ' ob ', off_ob(symj)
c$$$  write(6,*) symj, ' va ', off_va(symj)
c$$$  write(6,*) symj, ' vb ', off_vb(symj)
c$$$  end do
c     
      if(.not.ma_push_get(mt_dbl,noa*noa,'p_ij_a',l_p_ij_a,k_p_ij_a))
     $     call errquit('could not alloc p_ij_a',noa*noa, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nva*nva,'p_ab_a',l_p_ab_a,k_p_ab_a))
     $     call errquit('could not alloc p_ab_a',nva*nva, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nob*nob,'p_ij_b',l_p_ij_b,k_p_ij_b))
     $     call errquit('could not alloc p_ij_b',nob*nob, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nvb*nvb,'p_ab_b',l_p_ab_b,k_p_ab_b))
     $     call errquit('could not alloc p_ab_b',nvb*nvb, MA_ERR)
c     
      zero=0.0d0
      one=1.0d0
      half=0.5d0
      minusone=-1.0d0
      mhalf=-0.5d0
      tunitptra=1
      tunitptrb=1
      call dfill((noa*noa),zero,dbl_mb(k_p_ij_a),1)
      call dfill((nva*nva),zero,dbl_mb(k_p_ab_a),1)
      call dfill((nob*nob),zero,dbl_mb(k_p_ij_b),1)
      call dfill((nvb*nvb),zero,dbl_mb(k_p_ab_b),1)

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
      mhfss=mhalf*fss
      mfos=minusone*fos
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
            call mp2_read_tijab(nvb_lo, nvb_hi, irs_b, symia,
     $           num_ob, sym_hi_ob, sym_lo_ob, tunitb, tunitptrb, t_b)
c     
            tcounta=1
            tcountb=1
c     
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_oa(symj).gt.0 .and. num_va(symb).gt.0) then
                  pijacnt = off_oa(symj)+k_p_ij_a
                  call dgemm('n','t',num_oa(symj),
     $                 num_oa(symj),num_va(symb),mhfss,
     $                 t_a(tcounta),num_oa(symj),
     $                 t_a(tcounta),num_oa(symj),
     $                 one,dbl_mb(pijacnt),noa)
                  pabacnt = off_va(symb)+k_p_ab_a
                  call dgemm('t','n',num_va(symb),
     $                 num_va(symb),num_oa(symj),hfss,
     $                 t_a(tcounta),num_oa(symj),
     $                 t_a(tcounta),num_oa(symj),
     $                 one,dbl_mb(pabacnt),nva)
                  tcounta=tcounta+num_oa(symj)*num_va(symb)
               end if
               if(num_ob(symj).gt.0.and.num_vb(symb).gt.0)then
                  pijbcnt = off_ob(symj)+k_p_ij_b
                  call dgemm('n','t',num_ob(symj),
     $                 num_ob(symj),num_vb(symb),mfos,
     $                 t_b(tcountb),num_ob(symj),
     $                 t_b(tcountb),num_ob(symj),
     $                 one,dbl_mb(pijbcnt),nob)
                  pabbcnt = off_vb(symb)+k_p_ab_b
                  call dgemm('t','n',num_vb(symb),
     $                 num_vb(symb),num_ob(symj),fos,
     $                 t_b(tcountb),num_ob(symj),
     $                 t_b(tcountb),num_ob(symj),
     $                 one,dbl_mb(pabbcnt),nvb)
                  tcountb=tcountb+num_ob(symj)*num_vb(symb)
               end if
            end do
         end do
      end do
c     
      call ga_dgop(msg_p_ij_a_sum,dbl_mb(k_p_ij_a),noa*noa,'+')
      call ga_dgop(msg_p_ab_a_sum,dbl_mb(k_p_ab_a),nva*nva,'+')
      call ga_dgop(msg_p_ij_b_sum,dbl_mb(k_p_ij_b),nob*nob,'+')
      call ga_dgop(msg_p_ab_b_sum,dbl_mb(k_p_ab_b),nvb*nvb,'+')
      call ga_sync()
      if(ga_nodeid().eq.0)then
         call ga_acc(g_p_ij_a, noa_lo,noa_hi,noa_lo,noa_hi,
     $        dbl_mb(k_p_ij_a),noa,one)
         call ga_acc(g_p_ij_b,nob_lo,nob_hi,nob_lo,nob_hi,
     $        dbl_mb(k_p_ij_b),nob,one)
         call ga_acc(g_p_ab_a,1,nva,1,nva,dbl_mb(k_p_ab_a),nva,one)
         call ga_acc(g_p_ab_b,1,nvb,1,nvb,dbl_mb(k_p_ab_b),nvb,one)
      end if
c     
      if (.not. ma_chop_stack(l_p_ij_a)) call errquit
     $     ('mp2_pijab_uhf: cannot chop stack',0, MA_ERR)
c     
      end
      
      
      
      
      
      


