      subroutine argos_space_wtrst(lfnrst,filrst,lveloc,pres,
     + temp,tempw,temps,
     + iwl,iwlp,xw,vw,fw,xwcr,isl,islp,xs,vs,fs,xscr,
     + ipl,nb,ibw,bxw,bvw,bfw,brw,ibs,bxs,bvs,bfs,lseq)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer lfnrst,nb
      character*255 filrst
      logical lveloc
      real*8 pres,temp,tempw,temps
      integer iwl(mwm,miw2),isl(msa,mis2),lseq(mseq)
      integer iwlp(mwm,npackw),islp(msa,npack)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa),xwcr(mwm,3)
      real*8 xs(msa,3),vs(msa,3),fs(msa,3),xscr(msm,3)
      integer ipl(mbox,mip2),ibw(nb),ibs(nb,2)
      real*8 bxw(nb,3,mwa),bvw(nb,3,mwa),bfw(nb,3,mwa),brw(nb,3)
      real*8 bxs(nb,3),bvs(nb,3),bfs(nb,3)
c
      integer i,j,k,node,ncyc,icyc,numw,nums,number,nwmn,nsan
      integer ilp,ihp,jlp,jhp,ili,ihi,jli,jhi,ilw,ihw,jlw,jhw
      integer ils,ihs,jls,jhs
      character*10 rdate,rtime
      character*18 user
#ifdef USE_POSIXF
      integer ilen,ierror
#endif
      integer idyn,idynp,ihop
      logical lforces
c
      lforces=iguide.ne.0
c
      if(ga_nodeid().eq.0) then
c
      call swatch(rdate,rtime)
#ifdef USE_POSIXF
      call pxfgetlogin(user, ilen, ierror)
#else
      call getlog(user)
#endif
      if(user(18:18).ne.' ') user='                  '
c
      rewind(lfnrst)
      write(lfnrst,1000)
 1000 format('Restart file',/,' ',/,' ')
      write(lfnrst,1001) 4.2,rdate,rtime,nhist,lforces
 1001 format(f12.6,2a10,i5,4x,l1)
      hist(nhist)(1:18)=user
      hist(nhist)(19:28)=rdate
      hist(nhist)(29:48)=rtime
      hist(nhist)(49:108)=project(1:60)
      do 10 i=1,nhist
      write(lfnrst,1009) hist(i)
 1009 format(a)
   10 continue
      write(lfnrst,1002) npbtyp,nbxtyp,rsgm,((vlat(i,j),j=1,3),i=1,3)
 1002 format(2i5,f12.6,/,(3f12.6))
      write(lfnrst,1003) pres
 1003 format(1pe12.5)
      write(lfnrst,1004) temp,tempw,temps
 1004 format(3f12.6)
      write(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc,nsf,nseq,0,0
 1005 format(7i10,2i5)
c
      if(nwm.gt.0) then
      number=0
      ncyc=nwm/nb+1
      numw=nb
      do 1 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
c
c     begin test code 10/31/2001
c     initialize ibw to check that all atoms have been received
c
      do 1112 i=1,nb
      ibw(i)=-1
 1112 continue
c
c     end test code
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
      if(lveloc)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+3*mwa,jlw+6*mwa-1,vw,mwm)
      if(lforces)
     + call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa+3,jlw+9*mwa+2,fw,mwm)
      call ga_get(ga_w,ilw,ilw+nwmn-1,jlw+6*mwa,jlw+6*mwa+2,xwcr,mwm)
      do 3 i=1,nwmn
      j=iwl(i,lwgmn)-number
      if(j.gt.0.and.j.le.numw) then
      do 4 k=1,nwa
      bxw(j,1,k)=xw(i,1,k)
      bxw(j,2,k)=xw(i,2,k)
      bxw(j,3,k)=xw(i,3,k)
      bvw(j,1,k)=vw(i,1,k)
      bvw(j,2,k)=vw(i,2,k)
      bvw(j,3,k)=vw(i,3,k)
      if(lforces) then
      bfw(j,1,k)=fw(i,1,k)
      bfw(j,2,k)=fw(i,2,k)
      bfw(j,3,k)=fw(i,3,k)
      endif
    4 continue
      brw(j,1)=xwcr(i,1)
      brw(j,2)=xwcr(i,2)
      brw(j,3)=xwcr(i,3)
      ibw(j)=iwl(i,lwdyn)
      endif
    3 continue
      endif
    2 continue
      do 5 i=1,numw
      if(lveloc) then
      write(lfnrst,1006) ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
      else
      write(lfnrst,1006) ((bxw(i,j,k),j=1,3),(zero,j=1,3),k=1,nwa)
      endif
 1006 format(2x,6f13.8)
      if(lforces) write(lfnrst,1106) ((bfw(i,j,k),j=1,3),k=1,nwa)
 1106 format(2x,6e13.6)
      idyn=iand(ibw(i),12)/4
      idynp=iand(ibw(i),3)
      write(lfnrst,1007) idynp,idyn,(brw(i,k),k=1,3)
 1007 format(2i1,3f13.8)
c
c     begin test code 10/31/2001
c     check if al atoms have been received
c
      if(ibw(i).lt.0)
     + call md_abort('Missing solvent in wtrst',i)
c
c     end test code
c
    5 continue
      number=number+numw
    1 continue
      endif
c
      if(nsa.gt.0) then
      number=0
      ncyc=nsa/nb+1
      nums=nb
      do 6 icyc=1,ncyc
      if(nsa-number.lt.nums) nums=nsa-number
c
c     begin test code 10/31/2001
c     initialize ibw to check that all atoms have been received
c
      do 1117 i=1,nb
      ibs(i,1)=-1
      ibs(i,2)=0
 1117 continue
c
c     end test code
c
      do 7 node=np-1,0,-1
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
      if(lveloc) call ga_get(ga_s,ils,ils+nsan-1,jls+3,jls+5,vs,msa)
      if(lforces) call ga_get(ga_s,ils,ils+nsan-1,jls+6,jls+8,fs,msa)
      do 8 i=1,nsan
      j=isl(i,lsgan)-number
      if(j.gt.0.and.j.le.nums) then
      bxs(j,1)=xs(i,1)
      bxs(j,2)=xs(i,2)
      bxs(j,3)=xs(i,3)
      bvs(j,1)=vs(i,1)
      bvs(j,2)=vs(i,2)
      bvs(j,3)=vs(i,3)
      if(lforces) then
      bfs(j,1)=fs(i,1)
      bfs(j,2)=fs(i,2)
      bfs(j,3)=fs(i,3)
      endif
      ibs(j,1)=isl(i,lsdyn)
      ibs(j,2)=isl(i,lshop)
      endif
    8 continue
      endif
    7 continue
      do 9 i=1,nums
      idyn=iand(ibs(i,1),12)/4
      idynp=iand(ibs(i,1),3)
      ihop=ibs(i,2)
      if(iand(ihop,1).eq.1) then
      ihop=-(ihop/2)
      else
      ihop=ihop/2
      endif
      if(lveloc) then
      write(lfnrst,1008) idynp,idyn,(bxs(i,j),j=1,3),(bvs(i,j),j=1,3),
     + ihop
      else
      write(lfnrst,1008) idynp,idyn,(bxs(i,j),j=1,3),(zero,j=1,3),ihop
      endif
 1008 format(2i1,6f13.8,i5)
      if(lforces) write(lfnrst,1108) (bfs(i,j),j=1,3)
 1108 format(2x,3e13.6)
c
c     begin test code 10/31/2001
c     check if al atoms have been received
c
      if(ibs(i,1).lt.0)
     + call md_abort('Missing solute atom in wtrst',i)
c
c     end test code
c
    9 continue
      number=number+nums
    6 continue
      endif
c
      if(nsm.gt.0) then
      do 21 i=1,nsm
      write(lfnrst,1109) (xscr(i,j),j=1,3)
 1109 format(2x,3f13.8)
   21 continue
      endif
c
      if(nseq.gt.0) then
      write(lfnrst,1013) (lseq(i),i=1,nseq)
 1013 format(20i3)
      endif
c
      endif
c
      return
 9999 continue
      call md_abort('Failed to open restart file',me)
      return
      end
c $Id$
