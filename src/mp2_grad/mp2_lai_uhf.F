      subroutine mp2_lai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $     sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $     sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $     nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $     c_a,c_b,indx_v,nbf,ioff,nnbf,
     $     nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $     nob_lo,nob_hi,nva,nvb,
     $     noa,nob,t_a,t_b,ia_cb,ia_jk,ia_uv,tmp,tunita,tunitb,
     $     l_ai3_a,l_ai3_b,l_ai1_a,l_ai1_b,g_lai_a,g_lai_b,
     $     g_wai_a,g_wai_b,kunit)
*
* $Id: mp2_lai_uhf.F,v 1.9 1997-03-27 08:49:52 d3e129 Exp $
*
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "eaf.fh"
      integer i
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer sym_lo_oa(0:*),sym_hi_oa(0:*)
      integer sym_lo_va(0:*),sym_hi_va(0:*)
      integer num_oa(0:*),num_va(0:*)
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_ob(0:*),sym_hi_ob(0:*)
      integer sym_lo_vb(0:*),sym_hi_vb(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer indx_v(*)
      integer nbf,nnbf,ioff,kunit
      integer noa,nva,nob,nvb
      integer nir
      integer a_a
      integer oseg_lo,oseg_hi
      double precision c_a(nbf,nbf),c_b(nbf,nbf)
      double precision zero,one,minusone,half
      integer syma,symi,symb,symj,symia,symc,symk
      integer irs_a(*),irs_b(*)
      integer nva_lo_local,nva_hi_local
      double precision l_ai1_a(noa_lo:noa_hi,nva_lo:nva_hi)
      double precision l_ai1_b(nob_lo:nob_hi,nvb_lo:nvb_hi)
      double precision l_ai3_a(noa_lo:noa_hi,nva_lo:nva_hi)
      double precision l_ai3_b(nob_lo:nob_hi,nvb_lo:nvb_hi)
      double precision t_a(*),t_b(*),ia_cb(*),ia_jk(*),tmp(*),ia_uv(*)
      integer g_lai_a,g_lai_b
      integer tunita,tunitb
      double precision tunitptra,tunitptrb
      integer tcounta,tcountb
      integer msg_l_ai3_a_sum,msg_l_ai3_b_sum
      parameter(msg_l_ai3_a_sum=9992)
      parameter(msg_l_ai3_b_sum=9991)
      integer msg_l_ai1_a_sum,msg_l_ai1_b_sum
      parameter(msg_l_ai1_a_sum=9990)
      parameter(msg_l_ai1_b_sum=9989)
      integer i_batch
      integer g_wai_a,g_wai_b
c
#include "bitops.fh"
c     
      zero=0.0d0
      one=1.0d0
      minusone=-1.0d0
      half=0.5d0
      tunitptra=1
      tunitptrb=1
      call dfill((noa*nva),zero,l_ai3_a,1)
      call dfill((nob*nvb),zero,l_ai3_b,1)
      call dfill((noa*nva),zero,l_ai1_a,1)
      call dfill((nob*nvb),zero,l_ai1_b,1)
c     
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         i_batch=i-oseg_lo+1
         do a_a=nva_lo_local,nva_hi_local
            syma=irs_a(a_a)
            symia=ieor(syma,symi)
            call dfill((nbf*nbf),zero,tmp,1)
            call moints_rdhfint(kunit,nva_lo_local,nva_hi_local,
     $           i_batch,a_a,nnbf,ioff,ia_uv)
            call scatter(nnbf,tmp,indx_v,ia_uv)
            call dfill((nbf*nbf),zero,ia_uv,1)
            call upper2square(nbf,tmp,ia_uv)
c     
c...  form pure spin t*(vv|vo) and t*(vo|oo) terms 
c     
            call mp2_read_tijab(nva_lo, nva_hi, irs_a, symia,
     $           num_oa, sym_hi_oa, sym_lo_oa, tunita, tunitptra, t_a)
            tcounta=1
            do symb=0,nir-1
               symc=ieor(symia,symb)
               symj=symc
               if(num_va(symc).gt.0.and.num_va(symb).gt.0.
     $              and.num_oa(symj).gt.0)then
                  call dgemm('n','n',nbf,num_va(symc),nbf,one,
     $                 ia_uv,nbf,c_a(1,sym_lo_va(symc)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_va(symc),num_va(symb),
     $                 nbf,one,tmp,nbf,c_a(1,sym_lo_va(symb)),
     $                 nbf,zero,ia_cb,num_va(symc))
c     
c     Ljc <- sum(iab) Tjb,ia*(cb|ia)
c     
                  call dgemm('n','t',num_oa(symj),
     $                 num_va(symc),num_va(symb),
     $                 one,t_a(tcounta),num_oa(symj),
     $                 ia_cb,num_va(symc),one,
     $                 l_ai3_a(sym_lo_oa(symj),sym_lo_va(symc)),noa)
               end if
               tcounta=tcounta+num_oa(symj)*num_va(symb)
            end do
            tcounta=1
            do symb=0,nir-1
               symk=symb
               symj=ieor(symia,symk)
               if(num_oa(symj).gt.0.and.num_oa(symk).gt.0.
     $              and.num_va(symb).gt.0)then
                  call dgemm('n','n',nbf,num_oa(symj),nbf,one,
     $                 ia_uv,nbf,c_a(1,sym_lo_oa(symj)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_oa(symj),num_oa(symk),
     $                 nbf,one,tmp,nbf,c_a(1,sym_lo_oa(symk)),
     $                 nbf,zero,ia_jk,num_oa(symj))
c     
c     Lkb = sum(iaj) (ia|jk)*Tjb,ia
c     
                  call dgemm('t','n',num_oa(symk),
     $                 num_va(symb),num_oa(symj),
     $                 minusone,ia_jk,num_oa(symj),
     $                 t_a(tcounta),num_oa(symj),one,
     $                 l_ai1_a(sym_lo_oa(symk),sym_lo_va(symb)),noa)
               end if
               tcounta=tcounta+num_oa(symj)*num_va(symb)
            end do
c     
c...  form mixed spin t*(vv|vo) and t*(vo|oo) terms
c     
            call mp2_read_tijab(nvb_lo, nvb_hi, irs_b, symia,
     $           num_ob, sym_hi_ob, sym_lo_ob, tunitb, tunitptrb, t_b)
            tcountb=1
            do symb=0,nir-1
               symc=ieor(symia,symb)
               symj=symc
               if(num_vb(symc).gt.0.and.num_vb(symb).gt.0.
     $              and.num_ob(symj).gt.0)then
                  call dgemm('n','n',nbf,num_vb(symc),nbf,one,
     $                 ia_uv,nbf,c_b(1,sym_lo_vb(symc)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_vb(symc),num_vb(symb),
     $                 nbf,one,tmp,nbf,c_b(1,sym_lo_vb(symb)),
     $                 nbf,zero,ia_cb,num_vb(symc))
c     
c     Ljc <- sum(iab) Tjb,ia*(cb|ia)
c     
                  call dgemm('n','t',num_ob(symj),
     $                 num_vb(symc),num_vb(symb),
     $                 one,t_b(tcountb),num_ob(symj),
     $                 ia_cb,num_vb(symc),one,
     $                 l_ai3_b(sym_lo_ob(symj),sym_lo_vb(symc)),nob)
               end if
               tcountb=tcountb+num_ob(symj)*num_vb(symb)
            end do
            tcountb=1
            do symb=0,nir-1
               symk=symb
               symj=ieor(symia,symk)
               if(num_ob(symj).gt.0.and.num_ob(symk).gt.0.
     $              and.num_vb(symb).gt.0)then
                  call dgemm('n','n',nbf,num_ob(symj),nbf,one,
     $                 ia_uv,nbf,c_b(1,sym_lo_ob(symj)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_ob(symj),num_ob(symk),
     $                 nbf,one,tmp,nbf,c_b(1,sym_lo_ob(symk)),
     $                 nbf,zero,ia_jk,num_ob(symj))
c     
c     Lkb = sum(iaj) (ia|kj)*Tjb,ia
c     
                  call dgemm('t','n',num_ob(symk),
     $                 num_vb(symb),num_ob(symj),
     $                 minusone,ia_jk,num_ob(symj),
     $                 t_b(tcountb),num_ob(symj),one,
     $                 l_ai1_b(sym_lo_ob(symk),sym_lo_vb(symb)),nob)
               end if
               tcountb=tcountb+num_ob(symj)*num_vb(symb)
            end do
         end do
      end do
c     
      call ga_dgop(msg_l_ai3_a_sum,l_ai3_a,noa*nva,'+')
      call ga_dgop(msg_l_ai3_b_sum,l_ai3_b,nob*nvb,'+')
      call ga_dgop(msg_l_ai1_a_sum,l_ai1_a,noa*nva,'+')
      call ga_dgop(msg_l_ai1_b_sum,l_ai1_b,nob*nvb,'+')
      call ga_sync()
      if(ga_nodeid().eq.0)then
         call ga_acc(g_lai_a,1,noa,1,nva,l_ai3_a,noa,one)
         call ga_acc(g_lai_b,1,nob,1,nvb,l_ai3_b,nob,one)
         call ga_acc(g_lai_a,1,noa,1,nva,l_ai1_a,noa,one)
         call ga_acc(g_lai_b,1,nob,1,nvb,l_ai1_b,nob,one)
         call ga_acc(g_wai_a,1,noa,1,nva,l_ai1_a,noa,minusone)
         call ga_acc(g_wai_b,1,nob,1,nvb,l_ai1_b,nob,minusone)
      end if
c     
      end


