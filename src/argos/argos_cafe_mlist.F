      subroutine argos_cafe_mlist(nums,xs,fs,isgan,isfrc,isdyn,isct1,
     + msb,ibnd,msh,iang,msd,idih,mso,iimp,ixmul,imul,xmul,fmul,ess)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nums
      real*8 xs(msa,3),fs(msa,3)
      integer isgan(msa),isfrc(msa),isdyn(msa),isct1(msa)
      integer msb,msh,msd,mso
      integer ibnd(msb,3),iang(msh,4),idih(msd,5),iimp(mso,5)
      integer ixmul(mmuli),imul(mmult,4)
      real*8 xmul(mmult,3),fmul(mmult,3),ess(msf,msf,mpe,2)
c
      integer i,j,ix
c
      call timer_start(14)
      do 1 i=1,nmul
      fmul(i,1)=zero
      fmul(i,2)=zero
      fmul(i,3)=zero
      xmul(i,1)=zero
      xmul(i,2)=zero
      xmul(i,3)=zero
      imul(i,2)=0
      imul(i,3)=0
      imul(i,4)=0
      do 2 j=1,nums
      if(isgan(j).eq.imul(i,1)) then
      imul(i,2)=isfrc(j)
      imul(i,3)=isdyn(j)
      imul(i,4)=isct1(j)
      xmul(i,1)=xs(j,1)
      xmul(i,2)=xs(j,2)
      xmul(i,3)=xs(j,3)
      endif
    2 continue
      if(me.gt.0) imul(i,1)=0
    1 continue
      call timer_stop(14)
c
      call timer_start(15)
      call ga_dgop(mcf_46,xmul,3*mmult,'+')
      call ga_igop(mcf_47,imul,4*mmult,'+')
      call timer_stop(15)
c
      call timer_start(16)
      if(me.eq.0) then
      if(nmult(1).gt.0) then
      ix=1
      call argos_cafe_fsb(nmult(1),ixmul(ix),mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),imul(1,4),dbl_mb(i_chg),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,me.eq.0,dbl_mb(i_sti),
     + int_mb(i_lseq))
      endif
c
      if(nmult(2).gt.0) then
      ix=1+nmult(1)
      call argos_cafe_fsh(nmult(2),ixmul(ix),mht(2),mhp(2),
     + int_mb(i_iang(2)),dbl_mb(i_ang(2)),dbl_mb(i_rang(2)),
     + dbl_mb(i_rub(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),imul(1,4),dbl_mb(i_chg),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,me.eq.0,dbl_mb(i_sti))
      endif
c
      if(nmult(3).gt.0) then
      ix=1+nmult(1)+nmult(2)
      call argos_cafe_fsd(nmult(3),ixmul(ix),mdt(2),mdp(2),
     + int_mb(i_idih(2)),dbl_mb(i_dih(2)),dbl_mb(i_rdih(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,me.eq.0,dbl_mb(i_sti))
      endif
c
      if(nmult(4).gt.0) then
      ix=1+nmult(1)+nmult(2)+nmult(3)
      call argos_cafe_fso(nmult(4),ixmul(ix),mit(2),mip(2),
     + int_mb(i_iimp(2)),dbl_mb(i_imp(2)),dbl_mb(i_rimp(2)),
     + nmul,mmult,imul(1,1),imul(1,2),imul(1,3),
     + xmul,fmul,ess,.true.,.true.,me.eq.0,me.eq.0,dbl_mb(i_sti))
      endif
      endif
      call timer_stop(16)
c
      call timer_start(15)
      call ga_brdcst(mcf_74,fmul,3*mmult*ma_sizeof(mt_dbl,1,mt_byte),0)
      call timer_stop(15)
c
      call timer_start(16)
      do 3 i=1,nmul
      do 4 j=1,nums
      if(isgan(j).eq.imul(i,1)) then
      fs(j,1)=fs(j,1)+fmul(i,1)
      fs(j,2)=fs(j,2)+fmul(i,2)
      fs(j,3)=fs(j,3)+fmul(i,3)
      endif
    4 continue
    3 continue
      call timer_stop(16)
c
      return
      end
c $Id$
