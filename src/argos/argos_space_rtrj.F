      logical function argos_space_rtrj(lfntri,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,iwl,iwlp,xw,vw,fw,xwcr,
     + isl,islp,xs,vs,fs,ipl,nb,ibw,bxw,bvw,bfw,brw,ibs,bxs,bvs,bfs)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfntri,nb
      logical lxw,lvw,lfw,lxs,lvs,lfs
      real*8 stime,pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),fs(msa,3)
      integer ipl(mbox,mip2),ibw(nb),ibs(nb)
      real*8 bxw(nb,3,mwa),bvw(nb,3,mwa),bfw(nb,3,mwa),brw(nb,3)
      real*8 bxs(nb,3),bvs(nb,3),bfs(nb,3)
c
      integer i,j,k,node,ncyc,icyc,numw,nums,number,nwmn,nsan
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ilw,ihw,jlw,jhw
      integer ils,ihs,jls,jhs
      character*10 rdate,rtime
      character*80 card
c
      logical lpw,lps
c
      lpw=.false.
      lps=.false.
c
      if(me.eq.0) then
c
  100 continue
      read(lfntri,1000,end=9999) card
 1000 format(a)
      if(card(1:5).ne.'frame') goto 100
c
      read(lfntri,1001) stime,temp,pres,rdate,rtime
 1001 format(2f12.6,1pe12.5,1x,2a10)
      tmp=temp
      prs=pres
      read(lfntri,1002) ((vlat(i,j),j=1,3),i=1,3)
 1002 format(3f12.6)
      read(lfntri,1000,end=9999) card
c
      if(card(8:8).eq.'T'.or.card(8:8).eq.'F') then
      read(card,1003) lxw,lvw,lfw,lpw,lxs,lvs,lfs,lps,nwm,nwa,nsa
 1003 format(8l1,3i10)
      elseif(card(6:6).eq.'T'.or.card(6:6).eq.'F') then
      read(card,1023) lxw,lvw,lfw,lxs,lvs,lfs,nwm,nwa,nsa
 1023 format(6l1,3i10)
      lpw=.false.
      lps=.false.
      else
      read(card,1033) lxw,lvw,lxs,lvs,nwm,nwa,nsa
 1033 format(4l1,3i10)
      lpw=.false.
      lps=.false.
      lfw=.false.
      lfs=.false.
      endif
c
      if((lxw.or.lvw.or.lfw).and.nwm.gt.0) then
      number=0
      ncyc=nwm/nb+1
      numw=nb
      do 1 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
c
      if(nprec.eq.0) then
      if(lfw) then
      if(lvw) then
      do 61 i=1,numw
      read(lfntri,1014) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),
     + (bfw(i,j,k),j=1,3),k=1,nwa)
 1014 format(6f8.3,3f8.1)
   61 continue
      else
      do 71 i=1,numw
      read(lfntri,1015) ((bxw(i,j,k),j=1,3),(bfw(i,j,k),j=1,3),k=1,nwa)
 1015 format(3f8.3,3f8.1)
   71 continue
      endif
      else
      if(lvw) then
      do 6 i=1,numw
      read(lfntri,1004) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
 1004 format(6f8.3)
    6 continue
      else
      do 7 i=1,numw
      read(lfntri,1005) ((bxw(i,j,k),j=1,3),k=1,nwa)
 1005 format(3f8.3)
    7 continue
      endif
      endif
      else
      if(lfw) then
      if(lvw) then
      do 261 i=1,numw
      read(lfntri,2014) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),
     + (bfw(i,j,k),j=1,3),k=1,nwa)
 2014 format(6e12.6,/,3e12.6)
  261 continue
      else
      do 271 i=1,numw
      read(lfntri,2015) ((bxw(i,j,k),j=1,3),(bfw(i,j,k),j=1,3),k=1,nwa)
 2015 format(3e12.6,/,3e12.6)
  271 continue
      endif
      else
      if(lvw) then
      do 26 i=1,numw
      read(lfntri,2004) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
 2004 format(6e12.6)
   26 continue
      else
      do 27 i=1,numw
      read(lfntri,2005) ((bxw(i,j,k),j=1,3),k=1,nwa)
 2005 format(3e12.6)
   27 continue
      endif
      endif
      endif
c
      do 2 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nwmn=ipl(1,2)
      if(nwmn.gt.0) then
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      if(npackw.eq.0) then
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+lwdyn-1,iwl,mwm)
      else
      call ga_get(ga_iw,ili,ili+nwmn-1,jli,jli+npackw-1,iwlp,mwm)
      call argos_space_unpackw(nwmn,iwl,iwlp)
      endif
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      if(lvw)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      if(lfw)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa+3,jlw+9*mwa+2,fw,mwm)
      do 3 i=1,nwmn
      j=iwl(i,lwgmn)-number
      if(j.gt.0.and.j.le.numw) then
      do 4 k=1,nwa
      xw(i,1,k)=bxw(j,1,k)
      xw(i,2,k)=bxw(j,2,k)
      xw(i,3,k)=bxw(j,3,k)
    4 continue
      if(lvw) then
      do 5 k=1,nwa
      vw(i,1,k)=bvw(j,1,k)
      vw(i,2,k)=bvw(j,2,k)
      vw(i,3,k)=bvw(j,3,k)
    5 continue
      endif
      if(lfw) then
      do 51 k=1,nwa
      fw(i,1,k)=bfw(j,1,k)
      fw(i,2,k)=bfw(j,2,k)
      fw(i,3,k)=bfw(j,3,k)
   51 continue
      endif
      endif
    3 continue
      call ga_put(ga_w,ilw,ilw+nwmn-1,jlw,jlw+3*mwa-1,xw,mwm)
      if(lvw)
     + call ga_put(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      if(lfw)
     + call ga_put(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa+3,jlw+9*mwa+2,fw,mwm)
      endif
    2 continue
      number=number+numw
    1 continue
      endif
c
      if((lxs.or.lvs.or.lfs).and.nsa.gt.0) then
      number=0
      ncyc=nsa/nb+1
      nums=nb
      do 8 icyc=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
c
      if(nprec.eq.0) then
      if(lfs) then
      if(lvs) then
      do 111 i=1,nums
      read(lfntri,1016) (bxs(i,j),j=1,3),(bvs(i,j),j=1,3),
     + (bfs(i,j),j=1,3)
 1016 format(6f8.3,3f8.1)
  111 continue
      else
      do 121 i=1,nums
      read(lfntri,1017) (bxs(i,j),j=1,3),(bfs(i,j),j=1,3)
 1017 format(3f8.3,3f8.1)
  121 continue
      endif
      else
      if(lvs) then
      do 11 i=1,nums
      read(lfntri,1006) (bxs(i,j),j=1,3),(bvs(i,j),j=1,3)
 1006 format(6f8.3)
   11 continue
      else
      do 12 i=1,nums
      read(lfntri,1007) (bxs(i,j),j=1,3)
 1007 format(3f8.3)
   12 continue
      endif
      endif
      else
      if(lfs) then
      if(lvs) then
      do 2111 i=1,nums
      read(lfntri,2016) (bxs(i,j),j=1,3),(bvs(i,j),j=1,3),
     + (bfs(i,j),j=1,3)
 2016 format(6e12.6,/,3e12.6)
 2111 continue
      else
      do 2121 i=1,nums
      read(lfntri,2017) (bxs(i,j),j=1,3),(bfs(i,j),j=1,3)
 2017 format(3e12.6,/,3e12.6)
 2121 continue
      endif
      else
      if(lvs) then
      do 211 i=1,nums
      read(lfntri,2006) (bxs(i,j),j=1,3),(bvs(i,j),j=1,3)
 2006 format(6e12.6)
  211 continue
      else
      do 212 i=1,nums
      read(lfntri,2007) (bxs(i,j),j=1,3)
 2007 format(3e12.6)
  212 continue
      endif
      endif
      endif
c
      do 9 node=np-1,0,-1
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ihp,jlp,jhp,ipl,mbox)
      nsan=ipl(2,2)
      if(nsan.gt.0) then
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      if(npack.eq.0) then
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+lsdyn-1,isl,msa)
      else
      call ga_get(ga_is,ili,ili+nsan-1,jli,jli+npack-1,islp,msa)
      call argos_space_unpack(nsan,isl,islp)
      endif
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      call ga_get(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      if(lvs) call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      if(lfs) call ga_get(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      do 10 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      xs(i,1)=bxs(j,1)
      xs(i,2)=bxs(j,2)
      xs(i,3)=bxs(j,3)
      if(lvs) then
      vs(i,1)=bvs(j,1)
      vs(i,2)=bvs(j,2)
      vs(i,3)=bvs(j,3)
      endif
      if(lfs) then
      fs(i,1)=bfs(j,1)
      fs(i,2)=bfs(j,2)
      fs(i,3)=bfs(j,3)
      endif
      endif
   10 continue
      call ga_put(ga_s,ils,ils+nsan-1,jls,jls+2,xs,msa)
      if(lvs) call ga_put(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      if(lfs) call ga_put(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      endif
    9 continue
      number=number+nums
    8 continue
      endif
c
      endif
c
      call ga_brdcst(msp_15,vlat,
     + 9*ma_sizeof(mt_dbl,1,mt_byte),0)
c
      box(1)=vlat(1,1)
      box(2)=vlat(2,2)
      box(3)=vlat(3,3)
      boxh(1)=5.0d-1*box(1)
      boxh(2)=5.0d-1*box(2)
      boxh(3)=5.0d-1*box(3)
c
      argos_space_rtrj=.true.
c
      return
c
 9999 continue
c
      argos_space_rtrj=.false.
c
      return
c
      end
c
c $Id$
