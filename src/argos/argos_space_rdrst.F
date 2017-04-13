      subroutine argos_space_rdrst(lfnrst,filrst,lfntop,filtop,
     + temp,tempw,temps,ipl,xw,vw,fw,xwcr,iwl,xs,vs,fs,xscr,isl,
     + bxw,bvw,bfw,brw,ibw,nw,bxs,bvs,bfs,ibs,ns,ndx,
     + ibownr,boxsiz,lseq,isndx)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer argos_space_btop
      external argos_space_btop
c
      integer lfnrst,lfntop
      character*255 filrst,filtop
      integer nw,ns
      integer ibownr(maxbox,3)
      integer ipl(mbox,mip2),iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xs(msa,3),xwcr(mwm,3)
      real*8 vw(mwm,3,mwa),vs(msa,3),xscr(msm,3)
      real*8 fw(mwm,3,mwa),fs(msa,3)
      real*8 bxw(nw,3,mwa),bxs(ns,3),brw(nw,3)
      real*8 bvw(nw,3,mwa),bvs(ns,3)
      real*8 bfw(nw,3,mwa),bfs(ns,3)
      integer ibw(nw,2),ibs(ns,mis2)
      real*8 boxsiz(maxbox,3)
      integer ndx(nw),lseq(mseq),isndx(mseq)
      real*8 temp,tempw,temps
c
      character*1 cdum
      real*8 rdum,cgx,cgy,cgz
      integer i,j,k,idum,jdum,kdum,number,ncyc,numw
      integer icyc,ibx,iby,ibz,ipx,ipy,ipz,node,new,nold
      integer ilw,ihw,jlw,jhw,ils,ihs,jls,jhs
      integer ili,ihi,jli,jhi,ilp,ihp,jlp,jhp
c
      integer nat,nqt,naw,nbw,nhw,ndw,now,ntw,nnw,nsmr
      real*8 boxi(3)
      logical lforces
      character*80 card
c
      integer l,m,ib(3),nbox,nb(3),joff,npars,icount,lasts
      real*8 xtmin,xtmax,xtx,xt(3)
c
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY USE ENTERING ARGOS_SPACE_RDRST',0)
      endif
c
      lpbc9=.true.
      nsmr=0
      icount=0
      lasts=0
      ndums=0
c
      if(me.eq.0) then
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfnrst)
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + status='old',form='formatted',err=9899)
      rewind(lfntop)
c
      write(lfnout,6000)
 6000 format(/,' TOPOLOGY FILE INFORMATION',/)
c
      read(lfntop,2001,end=9897,err=9898) card
 2001 format(a80)
      write(lfnout,6001) card
 6001 format(' Title',t15,a)
      read(lfntop,2001,end=9897,err=9898) card
      write(lfnout,6002) card
 6002 format(t15,a)
      read(lfntop,2001,end=9897,err=9898) card
      write(lfnout,6002) card
c
      read(lfntop,2001,end=9897,err=9898) card
      write(lfnout,6003) card(1:12),card(13:32),card(33:42)
 6003 format(' Version',t15,a,/,' Date',t15,a,/,' Force field',t15,a)
c
      read(lfntop,2002,end=9897,err=9898) npars
      read(lfntop,2002,end=9897,err=9898) nat
      read(lfntop,2002,end=9897,err=9898) nqt
      read(lfntop,2002,end=9897,err=9898) nseq
 2002 format(i5)
      read(lfntop,2001,end=9897,err=9898)
      do 1103 k=1,npars
      do 102 i=1,nat
      read(lfntop,2001,end=9897,err=9898)
      do 103 j=i,nat
      read(lfntop,2001,end=9897,err=9898)
  103 continue
  102 continue
 1103 continue
      do 104 i=1,nqt*npars
      read(lfntop,2001,end=9897,err=9898)
  104 continue
      do 4104 i=1,nseq
      read(lfntop,2001,end=9897,err=9898)
 4104 continue
      read(lfntop,2003,end=9897,err=9898) naw,nbw,nhw,ndw,now,ntw,nnw
 2003 format(5i7,2i10)
      read(lfntop,2001,end=9897,err=9898)
      do 105 i=1,naw
      read(lfntop,2001,end=9897,err=9898)
  105 continue
      do 106 i=1,nbw*(npars+1)
      read(lfntop,2001,end=9897,err=9898)
  106 continue
      do 107 i=1,nhw*(npars+1)
      read(lfntop,2001,end=9897,err=9898)
  107 continue
      do 108 i=1,ndw*(npars+1)
      read(lfntop,2001,end=9897,err=9898)
  108 continue
      do 109 i=1,now*(npars+1)
      read(lfntop,2001,end=9897,err=9898)
  109 continue
      if(ntw.gt.0) then
      read(lfntop,2004,end=9897,err=9898)
      read(lfntop,2004,end=9897,err=9898)
 2004 format(11i7)
      endif
      if(nnw.gt.0) then
      read(lfntop,2005,end=9997,err=9998)
      read(lfntop,2005,end=9997,err=9998)
 2005 format(11i7)
      endif
      read(lfntop,2001,end=9897,err=9898)
      do 204 i=1,npars
      read(lfntop,2001,end=9897,err=9898)
  204 continue
c
      write(lfnout,6100)
 6100 format(/,' RESTART FILE INFORMATION',/)
c
      read(lfnrst,1001,end=9997,err=9998) card
 1001 format(a80)
      write(lfnout,6101) card
 6101 format(' Title',t15,a)
      read(lfnrst,1001,end=9997,err=9998) card
      write(lfnout,6102) card
 6102 format(t15,a)
      read(lfnrst,1001,end=9997,err=9998) card
      write(lfnout,6002) card
c
      read(lfnrst,1016) card(1:32),nhist,lforces
 1016 format(a32,i5,4x,l1)
      write(lfnout,6103) card(1:12),card(13:32)
 6103 format(' Version',t15,a,/,' Date',t15,a)
c
      if(nhist.gt.0) write(lfnout,6104)
 6104 format(/,' History',/)
c
      if(nhist.gt.0) then
      do 21 i=1,nhist
      read(lfnrst,1017) hist(i)
 1017 format(a)
      write(lfnout,6105) hist(i)
 6105 format(1x,a)
   21 continue
      endif
      if(nhist.lt.mxhist) then
      nhist=nhist+1
      else
      do 22 i=1,nhist-1
      hist(i)=hist(i+1)
   22 continue
      endif
      do 23 i=1,80
      hist(nhist)(i:i)=' '
   23 continue
      read(lfnrst,1002,end=9997,err=9998) npbtyp,nbxtyp
 1002 format(i5,i5)
      read(lfnrst,1003,end=9997,err=9998) ((vlat(i,j),j=1,3),i=1,3)
 1003 format(3f12.6)
      read(lfnrst,1004,end=9997,err=9998) jdum
 1004 format(40x,i5)
      read(lfnrst,1005,end=9997,err=9998) temp,tempw,temps
 1005 format(3f12.6)
      do 2 i=1,3
      box(i)=vlat(i,i)
      boxh(i)=half*box(i)
      boxi(i)=one/box(i)
      do 3 j=1,3
      vlati(i,j)=vlat(i,j)
    3 continue
    2 continue
      call matinv(vlati,3,3)
      if(jdum.ne.0) then
      read(lfnrst,1001,end=9997,err=9998) cdum
      endif
      read(lfnrst,1006,end=9997,err=9998) idum
 1006 format(70x,i5)
      if(idum.gt.0) then
      read(lfnrst,1007,end=9997,err=9998) idum,jdum,kdum
 1007 format(3i5)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,idum)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,jdum)
      read(lfnrst,1008,end=9997,err=9998) (rdum,i=1,kdum)
 1008 format(4e20.12)
      endif
c
      if(nwm.gt.0) then
      number=0
      ncyc=nwm/nw+1
      numw=nw
      do 4 icyc=1,ncyc
      if(nwm-number.lt.numw) numw=nwm-number
      do 44 i=1,numw
      read(lfnrst,1009,end=9997,err=9998)
     + ((bxw(i,j,k),j=1,3),(bvw(i,j,k),j=1,3),k=1,nwa)
 1009 format(2x,6f13.8)
      if(lforces) read(lfnrst,1109,end=9997,err=9998)
     + ((bfw(i,j,k),j=1,3),k=1,nwa)
 1109 format(2x,6e13.6)
      read(lfnrst,1010,end=9997,err=9998) ibw(i,2),(brw(i,k),k=1,3)
 1010 format(i1,1x,3f13.8)
   44 continue
      do 5 i=1,numw
      cgx=zero
      cgy=zero
      cgz=zero
      do 6 k=1,nwa
      cgx=cgx+bxw(i,1,k)
      cgy=cgy+bxw(i,2,k)
      cgz=cgz+bxw(i,3,k)
      if(.not.lforces) then
      bfw(i,1,k)=zero
      bfw(i,2,k)=zero
      bfw(i,3,k)=zero
      endif
    6 continue
      ibx=0
      iby=0
      ibz=0
      if(nbxtyp.ne.1) then
      xt(1)=cgx
      xt(2)=cgy
      xt(3)=cgz
      else
      xt(1)=box(1)*(vlati(1,1)*cgx+vlati(1,2)*cgy+vlati(1,3)*cgz)
      xt(2)=box(2)*(vlati(2,1)*cgx+vlati(2,2)*cgy+vlati(2,3)*cgz)
      xt(3)=box(3)*(vlati(3,1)*cgx+vlati(3,2)*cgy+vlati(3,3)*cgz)
      endif
      do 7 j=1,nbx-1
      if(xt(1)/nwa+boxh(1).gt.boxsiz(j,1)) ibx=j
    7 continue
      do 8 j=1,nby-1
      if(xt(2)/nwa+boxh(2).gt.boxsiz(j,2)) iby=j
    8 continue
      do 9 j=1,nbz-1
      if(xt(3)/nwa+boxh(3).gt.boxsiz(j,3)) ibz=j
    9 continue
c
      if(npbtyp.gt.0) then
      if(ibx.ge.nbx) ibx=ibx-nbx
      if(iby.ge.nby) iby=iby-nby
      if(ibx.lt.0) ibx=ibx+nbx
      if(iby.lt.0) iby=iby+nby
      if(npbtyp.eq.1) then
      if(ibz.ge.nbz) ibz=ibz-nbz
      if(ibz.lt.0) ibz=ibz+nbz
      else
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibz.lt.0) ibz=0
      endif
      else
      if(ibx.ge.nbx) ibx=nbx-1
      if(iby.ge.nby) iby=nby-1
      if(ibz.ge.nbz) ibz=nbz-1
      if(ibx.lt.0) ibx=0
      if(iby.lt.0) iby=0
      if(ibz.lt.0) ibz=0
      endif
      ipx=ibownr(ibx+1,1)
      ipy=ibownr(iby+1,2)
      ipz=ibownr(ibz+1,3)
c
      ndx(i)=(ipz*npy+ipy)*npx+ipx
      ibw(i,1)=(ibz*nby+iby)*nbx+ibx
      if(lnode0) then
      ndx(i)=0
      ibw(i,1)=0
      endif
    5 continue
      do 10 node=0,np-1
      new=0
      do 11 i=1,numw
      if(ndx(i).eq.node) then
      new=new+1
      iwl(new,lwgmn)=number+i
      iwl(new,lwnod)=node
      do 124 j=1,3
      do 12 k=1,nwa
      xw(new,j,k)=bxw(i,j,k)
      vw(new,j,k)=bvw(i,j,k)
      fw(new,j,k)=bfw(i,j,k)
   12 continue
      if(nserie.eq.0) then
      xwcr(new,j)=zero
      else
      xwcr(new,j)=brw(i,j)
      endif
  124 continue
      iwl(new,lwbox)=ibw(i,1)
      iwl(new,lwdyn)=5*ibw(i,2)
      endif
   11 continue
c
      if(new.gt.0) then
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      nold=ipl(1,2)
      if(nold+new.gt.mwm) call md_abort('Dimension mwm too small',0)
      ipl(1,2)=ipl(1,2)+new
      call ga_put(ga_ip,ilp,ilp,jlp,jhp,ipl,mbox)
      call ga_distribution(ga_iw,node,ili,ihi,jli,jhi)
      ili=ili+nold
      ihi=ili+new-1
      call ga_put(ga_iw,ili,ihi,jli,jhi,iwl,mwm)
      call ga_distribution(ga_w,node,ilw,ihw,jlw,jhw)
      ilw=ilw+nold
      ihw=ilw+new-1
      call ga_put(ga_w,ilw,ihw,jlw,jlw+3*mwa-1,xw(1,1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+3*mwa,jlw+6*mwa-1,vw(1,1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+6*mwa,jlw+6*mwa+2,xwcr(1,1),mwm)
      call ga_put(ga_w,ilw,ihw,jlw+6*mwa+3,jlw+9*mwa+2,fw(1,1,1),mwm)
      endif
c
   10 continue
      number=number+numw
    4 continue
      endif
c
      if(nsa.gt.0) then
      nb(1)=nbx
      nb(2)=nby
      nb(3)=nbz
      nbox=0
      k=1
      joff=0
      do 13 i=1,nsa
      if(k.gt.ns) call md_abort
     + ('Increase memory for argos_space_rdrst buffer',0)
      read(lfnrst,1011,end=9997,err=9998)
     + ibs(k,11),(bxs(k,j),j=1,3),(bvs(k,j),j=1,3),ibs(k,12)
 1011 format(i1,1x,6f13.8,i5)
      if(ibs(k,12).lt.0) ndums=ndums+1
      if(lforces) then
      read(lfnrst,1111,end=9997,err=9998) (bfs(k,j),j=1,3)
 1111 format(2x,6e13.6)
      else
      bfs(k,1)=zero
      bfs(k,2)=zero
      bfs(k,3)=zero
      endif
      read(lfntop,2009,end=9897,err=9898) (ibs(k,j),j=1,10)
 2009 format(16x,i3,4i7,5i5)
c 2009 format(16x,10i5)
      if(ibs(k,3).ne.lasts) then
      icount=icount+1
      lasts=ibs(k,3)
      isndx(icount)=lasts
      endif
      ibs(k,3)=icount
      if(nsmr.lt.ibs(k,2)) nsmr=ibs(k,2)
      if(nserie.eq.0) then
      if(nsf.lt.ibs(k,1)) nsf=ibs(k,1)
      else
      if(nsf.lt.ibs(k,1))
     + call md_abort('Error in number of solute fractions',nsf)
      endif
c
c     if segment of this atom is different from the segment of previous atom
c     then distribute all previous atoms in the list
c
      if(k.gt.1) then
      if(ibs(k,3).ne.ibs(k-1,3)) then
      new=k-1
      goto 14
      endif
      endif
c
c     if this is the last atom distribute
c
      if(i.eq.nsa) then
      new=k
      goto 14
      endif
c
c     read next atom
c
      k=k+1
      goto 13
c
c     distribute atoms 1 through new
c
   14 continue
c
c     determine the center of geometry
c
      do 15 l=1,3
      xtmax=bxs(1,l)
      xtmin=bxs(1,l)
      do 16 j=1,new
      xtmax=max(xtmax,bxs(j,l))
      xtmin=min(xtmin,bxs(j,l))
   16 continue
      xtx=0.5d0*(xtmax+xtmin)
      if(npbtyp.ne.0) then
      if(abs(xtx).gt.boxh(l)) then
      xtx=xtx-nint(xtx*boxi(l))*box(l)
      endif
      endif
      ib(l)=0
      do 17 m=1,nb(l)-1
      if(xtx+boxh(l).gt.boxsiz(m,l)) ib(l)=m
   17 continue
   15 continue
c
c     periodic boundaries
c
      if(npbtyp.gt.0) then
      m=2
      if(npbtyp.eq.1) m=3
      do 18 l=1,m
      if(ib(l).ge.nb(l)) ib(l)=ib(l)-nb(l)
      if(ib(l).lt.0) ib(l)=ib(l)+nb(l)
   18 continue
      if(npbtyp.gt.1) then
      if(ib(3).ge.nb(3)) ib(3)=nb(3)-1
      if(ib(3).lt.0) ib(3)=0
      endif
      else
      do 19 l=1,3
      if(ib(l).ge.nb(l)) ib(l)=nb(l)-1
      if(ib(l).lt.0) ib(l)=0
   19 continue
      endif
c
c     determine owning node
c
      if(.not.lnode0) nbox=(ib(3)*nb(2)+ib(2))*nb(1)+ib(1)
      node=argos_space_btop(nbox,ibownr)
c
      do 120 j=1,new
      isl(j,lsgan)=joff+j
      isl(j,lsfrc)=ibs(j,1)
      isl(j,lsmol)=ibs(j,2)
      isl(j,lssgm)=ibs(j,3)
      isl(j,lsgrp)=ibs(j,4)
      isl(j,lspgr)=ibs(j,5)
      isl(j,lsatt)=ibs(j,6)
      isl(j,lsct1)=ibs(j,7)
      isl(j,lsct2)=ibs(j,8)
      isl(j,lsct3)=ibs(j,9)
      isl(j,lssss)=ibs(j,10)
      isl(j,lsdyn)=5*ibs(j,11)
      isl(j,lsbox)=nbox
      isl(j,lsnod)=node
      isl(j,lshop)=0
      if(ibs(j,12).gt.0) isl(j,lshop)=ibs(j,12)*2
      if(ibs(j,12).lt.0) isl(j,lshop)=(-ibs(j,12))*2+1
      do 121 l=1,3
      xs(j,l)=bxs(j,l)
      vs(j,l)=bvs(j,l)
      fs(j,l)=bfs(j,l)
  121 continue
  120 continue
      joff=joff+new
c
c     communicate data to node
c
      call ga_distribution(ga_ip,node,ilp,ihp,jlp,jhp)
      call ga_get(ga_ip,ilp,ilp+1,jlp,jhp,ipl,mbox)
      nold=ipl(2,2)
      if(nold+new.gt.msa)
     + call md_abort('Dimension msa too small (1)',nold+new)
      ipl(2,2)=ipl(2,2)+new
      call ga_put(ga_ip,ilp,ilp+1,jlp,jhp,ipl,mbox)
      call ga_distribution(ga_is,node,ili,ihi,jli,jhi)
      ili=ili+nold
      ihi=ili+new-1
      call ga_put(ga_is,ili,ihi,jli,jhi,isl,msa)
      call ga_distribution(ga_s,node,ils,ihs,jls,jhs)
      ils=ils+nold
      ihs=ils+new-1
      call ga_put(ga_s,ils,ihs,jls,jls+2,xs(1,1),msa)
      call ga_put(ga_s,ils,ihs,jls+3,jls+5,vs(1,1),msa)
      call ga_put(ga_s,ils,ihs,jls+9,jls+11,fs(1,1),msa)
c
c     make first atom of next segment first in the list
c
      if(k.gt.new) then
      do 122 j=1,12
      ibs(1,j)=ibs(k,j)
  122 continue
      do 123 j=1,3
      bxs(1,j)=bxs(k,j)
      bvs(1,j)=bvs(k,j)
      bfs(1,j)=bfs(k,j)
  123 continue
      k=2
      if(i.eq.nsa) then
      new=1
      k=1
      goto 14
      endif
      endif
c
   13 continue
      endif
c
      if(nsm.gt.0) then
      do 31 i=1,nsm
      xscr(i,1)=zero
      xscr(i,2)=zero
      xscr(i,3)=zero
   31 continue
      do 32 i=1,nsm
      read(lfnrst,1012,end=99,err=99) (xscr(i,j),j=1,3)
 1012 format(2x,3f13.8)
      if(nserie.eq.0) then
      xscr(i,1)=zero
      xscr(i,2)=zero
      xscr(i,3)=zero
      endif
   32 continue
      endif
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY IN ARGOS_SPACE_RDRST',1)
      endif
c
      if(nseq.gt.0) then
      read(lfnrst,1013) (lseq(i),i=1,nseq)
 1013 format(20i3)
      endif
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY IN ARGOS_SPACE_RDRST',2)
      endif
c
   99 continue
c
      close(unit=lfnrst,status='keep')
      close(unit=lfntop,status='keep')
c
      endif
c
      if(nseq.gt.0) then
      call ga_brdcst(msp_22,lseq,nseq*ma_sizeof(mt_int,1,mt_byte),0)
      endif
      call ga_brdcst(msp_28,ndums,ma_sizeof(mt_int,1,mt_byte),0)
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY USE EXITING ARGOS_SPACE_RDRST',0)
      endif
c
      return
 9897 call md_abort('EOF encountered on topology file',0)
 9898 call md_abort('Error reading topology file',1)
 9899 call md_abort('Error opening topology file',2)
 9997 call md_abort('EOF encountered on restart file',0)
 9998 call md_abort('Error reading restart file',0)
 9999 call md_abort('Error opening restart file',0)
      return
      end
c $Id$
